"""
card_pipeline_final.py

End-to-end CARD -> tBLASTn -> GenBank extractor.
- Directly scrapes Protein FASTA from CARD.
- Uses tBLASTn to find ALL synonymous mutations with 100% protein identity in K. pneumoniae.
- Uses proper BLAST translation frames to handle genomic strand orientation natively.
- STRICT FILTER 1: Ensures exactly 1000bp upstream and downstream are available.
- STRICT FILTER 2: Ensures the extracted CDS begins with a valid start codon (ATG, GTG, TTG, CTG).
- Logs comprehensive details (and reasons for skipping) into an Excel summary.
"""

import time
import re
import os
import logging
from datetime import datetime
from urllib.parse import urljoin
import subprocess
import io
import sys

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
import pandas as pd
from bs4 import BeautifulSoup
from Bio import SeqIO, Entrez
from Bio.Blast import NCBIWWW, NCBIXML
from Bio.Seq import Seq

try:
    from tqdm import tqdm
    TQDM = True
except Exception:
    TQDM = False

# ---------- USER INPUT ----------
default_ontology = "36023"
ont_input = input(f"Enter CARD ontology number (digits after 'ontology/'), default {default_ontology}: ").strip() or default_ontology
EMAIL = input("Enter your Entrez email (required by NCBI): ").strip()
API_KEY = input("Optional NCBI API key (press Enter to skip): ").strip() or None

# ---------- CONFIG ----------
BASE_URL = "https://card.mcmaster.ca"
PARENT_ONTOLOGY = f"ontology/{ont_input}"
START_URL = urljoin(BASE_URL, PARENT_ONTOLOGY)

OUTPUT_PREFIX = f"card_{ont_input}_final"
OUT_FASTA_CDS = OUTPUT_PREFIX + "_cds_with_stop.fasta"
OUT_FASTA_UPSTREAM = OUTPUT_PREFIX + "_upstream_1000.fasta"
OUT_FASTA_PROMOTER = OUTPUT_PREFIX + "_promoter_300.fasta"
OUT_FASTA_DOWNSTREAM = OUTPUT_PREFIX + "_downstream_1000.fasta"
OUT_EXCEL = OUTPUT_PREFIX + "_summary.xlsx"

REQUEST_DELAY = 1.0
UPSTREAM_BP = 1000
PROMOTER_BP = 300
DOWNSTREAM_BP = 1000

# Valid start codons for bacteria (including the rare alternative CTG)
VALID_START_CODONS = {"ATG", "GTG", "TTG", "CTG"}

BLAST_HITLIST_SIZE = 500 
HTTP_TIMEOUT = 30
MAX_SAFEGET_RETRIES = 3
RETRY_BACKOFF_FACTOR = 1.5

BLAST_TIMEOUT = 300
BLAST_RETRY_ATTEMPTS = 3

Entrez.email = EMAIL
if API_KEY:
    Entrez.api_key = API_KEY

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)-7s %(message)s", datefmt="%H:%M:%S")
log = logging.getLogger("card-pipeline")

session = requests.Session()
session.headers.update({"User-Agent": "CARD-tblastn-FinalHits/6.0"})
retry_strategy = Retry(total=5, backoff_factor=0.8, status_forcelist=[429, 500, 502, 503, 504])
adapter = HTTPAdapter(max_retries=retry_strategy)
session.mount("https://", adapter)
session.mount("http://", adapter)

def now_ts(): return datetime.now().strftime("%Y-%m-%d %H:%M:%S")

def safe_get(url, timeout=HTTP_TIMEOUT, retries=MAX_SAFEGET_RETRIES):
    attempt, backoff = 0, 1.0
    while attempt < retries:
        try:
            resp = session.get(url, timeout=timeout)
            resp.raise_for_status()
            return resp
        except Exception as e:
            attempt += 1
            time.sleep(backoff)
            backoff *= RETRY_BACKOFF_FACTOR
    raise RuntimeError(f"GET permanently failed for {url}")

def extract_protein_from_soup(soup):
    div = soup.find(id="fastaProtein")
    if not div:
        return None, None
        
    text = div.get_text("\n").strip()
    match = re.search(r'(>[^\n\r]+)[\n\r]+([A-Z\s\-]+)', text, re.IGNORECASE)
    if match:
        header = match.group(1).strip()
        seq = re.sub(r'[^A-Z]', '', match.group(2).upper())
        return header, seq
    return None, None

def get_child_links(parent_soup):
    cont = parent_soup.find(id="show-children")
    tags = cont.find_all("a", href=True) if cont else parent_soup.find_all('a', href=re.compile(r"^/ontology/\d+$"))
    return [(urljoin(BASE_URL, a['href']), a.get_text(strip=True)) for a in tags]

def run_tblastn_remote(protein_seq, entrez_query="Klebsiella pneumoniae[Organism]"):
    log.info(f"  BLAST: starting remote tBLASTn (max hits: {BLAST_HITLIST_SIZE})...")
    attempt_params = [
        {"hitlist_size": BLAST_HITLIST_SIZE, "expect": 1e-10},
        {"hitlist_size": max(5, BLAST_HITLIST_SIZE // 2), "expect": 1e-5},
        {"hitlist_size": max(3, BLAST_HITLIST_SIZE // 5), "expect": 1e-2}
    ]

    for i, params in enumerate(attempt_params[:BLAST_RETRY_ATTEMPTS], start=1):
        py_snippet = f"""
import sys
from Bio.Blast import NCBIWWW
seq = sys.stdin.read()
handle = NCBIWWW.qblast(
    "tblastn", "nt", seq,
    entrez_query={repr(entrez_query)},
    hitlist_size={params["hitlist_size"]}, expect={params["expect"]}, format_type="XML"
)
sys.stdout.write(handle.read())
handle.close()
"""
        try:
            proc = subprocess.run([sys.executable, "-c", py_snippet], input=protein_seq.encode(), stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=BLAST_TIMEOUT)
            if proc.returncode == 0 and proc.stdout.strip():
                return NCBIXML.read(io.StringIO(proc.stdout.decode()))
        except Exception:
            pass
        time.sleep(1.0)
    raise Exception("tBLASTn failed after all attempts.")

def find_100pct_protein_hits(blast_record, prot_len):
    """Finds ALL hits that share 100% identity over the full length."""
    out = []
    for aln in blast_record.alignments:
        acc = getattr(aln, "accession", None) or (aln.hit_def.split()[0] if aln.hit_def else None)
        for hsp in aln.hsps:
            if getattr(hsp, "identities", 0) == prot_len and getattr(hsp, "align_length", 0) == prot_len:
                out.append((acc, hsp, aln.hit_def))
                break 
    return out

def fetch_genbank(accession):
    handle = Entrez.efetch(db="nucleotide", id=accession, rettype="gb", retmode="text")
    rec = SeqIO.read(handle, "genbank")
    handle.close()
    return rec

def extract_regions_with_stop(rec, hsp):
    s_start = hsp.sbjct_start
    s_end = hsp.sbjct_end
    seq_len = len(rec.seq)
    
    # CORRECT STRAND DETECTION: Use BLAST frame, not coordinates
    sbjct_frame = hsp.frame[1] 
    
    # Standardize low/high coordinates
    low_coord = min(s_start, s_end)
    high_coord = max(s_start, s_end)
    
    if sbjct_frame > 0:
        strand = '+'
        low = low_coord
        high = high_coord + 3  # Add 3bp to 3' end for stop codon
        
        up_start = max(1, low - UPSTREAM_BP)
        promoter_start = max(1, low - PROMOTER_BP)
        up_end = low - 1
        
        down_start = high + 1
        down_end = min(seq_len, high + DOWNSTREAM_BP)
    else:
        strand = '-'
        low = low_coord - 3  # Add 3bp to 3' end (lower coordinate) for stop codon
        high = high_coord
        
        up_start = high + 1
        up_end = min(seq_len, high + UPSTREAM_BP)
        promoter_end = min(seq_len, high + PROMOTER_BP)
        
        down_start = max(1, low - DOWNSTREAM_BP)
        down_end = low - 1

    # Broad bounds check (is it too close to the edge of the contig?)
    if up_start < 1 or down_end > seq_len:
        return {"error": "INCOMPLETE_FLANKS_OOB"}

    # Extract the raw genomic sequence
    cds_genomic = rec.seq[low - 1 : high]
    
    if strand == '+':
        up_seq = rec.seq[up_start - 1 : up_end]
        promoter_seq = rec.seq[promoter_start - 1 : up_end]
        down_seq = rec.seq[down_start - 1 : down_end]
        cds_seq = cds_genomic
    else:
        # Reverse complement everything for minus strand so output reads 5' -> 3'
        up_seq = Seq(rec.seq[high : up_end]).reverse_complement()
        promoter_seq = Seq(rec.seq[high : promoter_end]).reverse_complement()
        down_seq = Seq(rec.seq[down_start - 1 : down_end]).reverse_complement()
        cds_seq = Seq(cds_genomic).reverse_complement()

    # Strict length check for flanks
    if len(up_seq) != UPSTREAM_BP or len(down_seq) != DOWNSTREAM_BP:
        return {"error": f"INCOMPLETE_FLANKS_EXACT (Up: {len(up_seq)}, Down: {len(down_seq)})"}

    # Strict Start Codon Check (Now correctly evaluating reverse-complemented sequences)
    start_codon = str(cds_seq[:3]).upper()
    if start_codon not in VALID_START_CODONS:
        return {"error": f"INVALID_START_CODON ({start_codon})"}

    return {
        "strand": strand,
        "up_seq": up_seq,
        "promoter_seq": promoter_seq,
        "down_seq": down_seq,
        "cds_seq": cds_seq,
        "low": low, "high": high
    }

def append_fasta_file(path, header, seq, wrap=70):
    with open(path, "a") as fh:
        fh.write(header + "\n")
        s = str(seq)
        for i in range(0, len(s), wrap):
            fh.write(s[i:i+wrap] + "\n")

def append_and_save_row(row):
    df_row = pd.DataFrame([row])
    if not os.path.exists(OUT_EXCEL):
        df_row.to_excel(OUT_EXCEL, index=False)
    else:
        current = pd.read_excel(OUT_EXCEL)
        # Avoid the Future Warning by filtering out empty/NA dataframes
        if not current.empty:
            updated = pd.concat([current, df_row], ignore_index=True)
        else:
            updated = df_row
        updated.to_excel(OUT_EXCEL, index=False)

def main():
    # Initialize/Clear FASTA files at the start of the script
    open(OUT_FASTA_CDS, "w").close()
    open(OUT_FASTA_UPSTREAM, "w").close()
    open(OUT_FASTA_PROMOTER, "w").close()
    open(OUT_FASTA_DOWNSTREAM, "w").close()
    # Delete the Excel sheet if it exists so we start fresh
    if os.path.exists(OUT_EXCEL):
        os.remove(OUT_EXCEL)

    try: 
        parent_soup = BeautifulSoup(safe_get(START_URL).text, "lxml")
    except Exception as e: 
        return log.error(f"Failed parent: {e}")

    children = get_child_links(parent_soup)
    iterator = tqdm(children, desc="CARD genes") if TQDM else enumerate(children, 1)
    
    for item in iterator:
        url, gene = item if not TQDM else item
        process_child(url, gene)
        
    log.info(f"Pipeline complete! Details saved to {OUT_EXCEL}")

def process_child(url, gene):
    base_row = {
        "timestamp": now_ts(),
        "card_gene": gene,
        "protein_len": None,
        "hit_accession": None,
        "hit_definition": None,
        "strand": None,
        "genomic_low": None,
        "genomic_high": None,
        "cds_len_with_stop": None,
        "status": None
    }
    
    try:
        psoup = BeautifulSoup(safe_get(url).text, "lxml")
        header, prot_seq = extract_protein_from_soup(psoup)
        
        if not prot_seq:
            base_row["status"] = "ERROR: No Protein FASTA found"
            return append_and_save_row(base_row)

        prot_len = len(prot_seq)
        base_row["protein_len"] = prot_len
        log.info(f"Targeting {gene} (Protein Length: {prot_len} AA)")
        
        blast_record = run_tblastn_remote(prot_seq)
        hits = find_100pct_protein_hits(blast_record, prot_len)
        
        if not hits:
            base_row["status"] = "ERROR: No 100% protein hits found in BLAST"
            return append_and_save_row(base_row)

        log.info(f"  Found {len(hits)} matching 100% hits for {gene}. Processing all...")

        for acc, hsp, hit_def in hits:
            if not acc: continue
            
            hit_row = base_row.copy()
            hit_row["hit_accession"] = acc
            hit_row["hit_definition"] = hit_def[:150] 
            
            try:
                rec = fetch_genbank(acc)
                regions = extract_regions_with_stop(rec, hsp)
                
                # Handle our custom strict filtering errors
                if "error" in regions:
                    err_msg = regions["error"]
                    hit_row["status"] = f"SKIPPED: {err_msg}"
                    append_and_save_row(hit_row)
                    log.info(f"    Skipping {acc}: {err_msg}")
                    continue 
                    
                # Format Fasta Headers
                hdr_base = f">{gene}|{acc}|{regions['low']}-{regions['high']}|strand={regions['strand']}"
                
                # Append to respective files
                append_fasta_file(OUT_FASTA_CDS, hdr_base + "|part=cds_with_stop", regions['cds_seq'])
                append_fasta_file(OUT_FASTA_UPSTREAM, hdr_base + "|part=upstream_1000", regions['up_seq'])
                append_fasta_file(OUT_FASTA_PROMOTER, hdr_base + "|part=promoter_300", regions['promoter_seq'])
                append_fasta_file(OUT_FASTA_DOWNSTREAM, hdr_base + "|part=downstream_1000", regions['down_seq'])
                
                # Log successful extraction
                hit_row["strand"] = regions['strand']
                hit_row["genomic_low"] = regions['low']
                hit_row["genomic_high"] = regions['high']
                hit_row["cds_len_with_stop"] = len(regions['cds_seq'])
                hit_row["status"] = "SUCCESS: Extracted strictly verified regions"
                append_and_save_row(hit_row)
                log.info(f"    Saved: {acc} (CDS: {len(regions['cds_seq'])}bp, Valid Start/Flanks)")
                
            except Exception as e:
                hit_row["status"] = f"ERROR: Failed during GenBank fetch/parse ({str(e)})"
                append_and_save_row(hit_row)
                log.warning(f"    Error processing {acc}: {e}")
            
            time.sleep(REQUEST_DELAY)
                
    except Exception as e:
        base_row["status"] = f"ERROR: Pipeline failure ({str(e)})"
        append_and_save_row(base_row)

if __name__ == "__main__":
    main()
