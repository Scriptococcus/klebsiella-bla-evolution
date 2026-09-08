#!/usr/bin/env python3
"""
FASTA -> GenBank pipeline (row-preserving version with advanced reporting)

What this script does
---------------------
1. Reads every FASTA record as its own row.
   - If the same accession appears 2, 3, or 100 times, all rows stay.
   - Nothing is dropped from the FASTA input.
2. Extracts the accession from each FASTA header.
3. Fetches GenBank records for unique accessions only for efficiency.
4. Writes an Excel workbook with an executive summary and expanded rows.

Usage
-----
python fasta_to_genbank_pipeline_rowwise.py \
    --fasta kpccds.fasta \
    --output kpcgen_final.xlsx \
    --email your_email@example.com \
    --api-key YOUR_NCBI_API_KEY
"""

from __future__ import annotations

import argparse
import hashlib
import os
import re
import time
from collections import Counter
from io import StringIO
from pathlib import Path
from typing import Dict, List, Tuple

import pandas as pd
from Bio import Entrez, SeqIO
from Bio.SeqRecord import SeqRecord
from tqdm import tqdm


# -------------------------
# Accession detection
# -------------------------
ACC_PATTERNS = [
    re.compile(
        r"(?:[A-Z]{1,3}_[A-Z]{1,3}\d{5,9}(?:\.\d+)?|"
        r"[A-Z]{1,3}\d{5,9}(?:\.\d+)?)",
        re.IGNORECASE
    ),
]


def normalize_accession(x) -> str:
    if pd.isna(x):
        return ""
    return str(x).strip().upper()


def accession_base(acc: str) -> str:
    acc = normalize_accession(acc)
    return acc.split(".")[0]


def looks_like_accession(token: str) -> bool:
    token = normalize_accession(token)
    return any(p.fullmatch(token) for p in ACC_PATTERNS)


def extract_accession_from_fasta_header(record: SeqRecord) -> Tuple[str, str]:
    header = (record.description or record.id or "").strip()
    if not header:
        return "", ""

    parts = [p.strip() for p in header.split("|") if p.strip()]

    # 1) Look in pipe-separated fields
    for i, part in enumerate(parts):
        if looks_like_accession(part):
            variant = "|".join(parts[:i]).strip()
            return variant, normalize_accession(part)

    # 2) Look in whitespace-separated tokens
    for token in re.split(r"\s+", header):
        token = token.strip().strip(",;:()[]{}<>")
        if looks_like_accession(token):
            return "", normalize_accession(token)

    # 3) Search inside the whole header
    for pat in ACC_PATTERNS:
        m = pat.search(header)
        if m:
            return "", normalize_accession(m.group(0))

    return "", ""


def md5_of_sequence(seq: str) -> str:
    return hashlib.md5(seq.encode("utf-8")).hexdigest()


def chunked(items: List[str], size: int) -> List[List[str]]:
    return [items[i:i + size] for i in range(0, len(items), size)]


def join_qualifier(v) -> str:
    if v is None:
        return ""
    if isinstance(v, list):
        return " | ".join(str(i) for i in v)
    return str(v)


def get_source_feature(record: SeqRecord):
    for feat in record.features:
        if feat.type == "source":
            return feat
    return None


def classify_replicon(record: SeqRecord) -> Tuple[str, str]:
    source_feat = get_source_feature(record)
    evidence = []

    record_text = " ".join([
        record.description or "",
        str(record.annotations.get("source", "")),
        str(record.annotations.get("comment", "")),
    ]).lower()

    if source_feat is not None:
        for key in ["plasmid", "chromosome", "note", "db_xref"]:
            if key in source_feat.qualifiers:
                val = join_qualifier(source_feat.qualifiers.get(key, []))
                if val:
                    evidence.append(f"{key}={val}")

    if source_feat is not None and "plasmid" in source_feat.qualifiers:
        return "plasmid", "; ".join(evidence) if evidence else "source qualifier plasmid"
    if source_feat is not None and "chromosome" in source_feat.qualifiers:
        return "chromosome", "; ".join(evidence) if evidence else "source qualifier chromosome"

    if "plasmid" in record_text:
        return "plasmid", "; ".join(evidence) if evidence else "record text contains plasmid"
    if "chromosome" in record_text:
        return "chromosome", "; ".join(evidence) if evidence else "record text contains chromosome"

    return "unknown", "; ".join(evidence) if evidence else "no explicit chromosome/plasmid label"


# -------------------------
# Entrez fetch helpers
# -------------------------
def fetch_batch_records(accessions: List[str], max_retries: int = 3) -> Dict[str, SeqRecord]:
    records: Dict[str, SeqRecord] = {}
    if not accessions:
        return records

    ids = ",".join(accessions)
    last_err = None

    for attempt in range(1, max_retries + 1):
        try:
            with Entrez.efetch(
                db="nuccore",
                id=ids,
                rettype="gb",
                retmode="text",
            ) as handle:
                text = handle.read()

            parsed_records = list(SeqIO.parse(StringIO(text), "genbank"))

            for record in parsed_records:
                full = normalize_accession(record.id)
                base = accession_base(full)

                records[full] = record
                records[base] = record

                for acc in record.annotations.get("accessions", []):
                    accn = normalize_accession(acc)
                    records[accn] = record
                    records[accession_base(accn)] = record

            return records

        except Exception as e:
            last_err = e
            if attempt < max_retries:
                time.sleep(1.2 * attempt)
            else:
                raise RuntimeError(f"Batch fetch failed for {len(accessions)} accessions: {e}") from e

    raise RuntimeError(f"Batch fetch failed: {last_err}")


def fetch_single_record(accession: str, max_retries: int = 3) -> SeqRecord:
    last_err = None
    for attempt in range(1, max_retries + 1):
        try:
            with Entrez.efetch(
                db="nuccore",
                id=accession,
                rettype="gb",
                retmode="text",
            ) as handle:
                text = handle.read()
            return SeqIO.read(StringIO(text), "genbank")
        except Exception as e:
            last_err = e
            if attempt < max_retries:
                time.sleep(1.2 * attempt)
            else:
                raise RuntimeError(f"Failed to fetch {accession}: {e}") from e
    raise RuntimeError(f"Failed to fetch {accession}: {last_err}")


# -------------------------
# Record -> rows
# -------------------------
def summary_from_record(record: SeqRecord, input_accession: str) -> dict:
    source_feat = get_source_feature(record)
    replicon_type, replicon_evidence = classify_replicon(record)
    feat_counts = Counter(f.type for f in record.features)

    source_quals = {}
    if source_feat is not None:
        for k, v in source_feat.qualifiers.items():
            source_quals[k] = join_qualifier(v)

    row = {
        "input_accession": input_accession,
        "record_id": record.id,
        "record_name": record.name,
        "accession": record.annotations.get("accessions", [""])[0] if record.annotations.get("accessions") else "",
        "definition": record.description or "",
        "organism": record.annotations.get("organism", ""),
        "source": record.annotations.get("source", ""),
        "molecule_type": record.annotations.get("molecule_type", ""),
        "topology": record.annotations.get("topology", ""),
        "data_file_division": record.annotations.get("data_file_division", ""),
        "sequence_length_bp": len(record.seq),
        "references_count": len(record.annotations.get("references", [])),
        "feature_count_total": len(record.features),
        "feature_count_source": feat_counts.get("source", 0),
        "feature_count_gene": feat_counts.get("gene", 0),
        "feature_count_cds": feat_counts.get("CDS", 0),
        "feature_count_trna": feat_counts.get("tRNA", 0),
        "feature_count_rrna": feat_counts.get("rRNA", 0),
        "replicon_type": replicon_type,
        "replicon_evidence": replicon_evidence,
        "taxonomy": " ; ".join(record.annotations.get("taxonomy", [])) if record.annotations.get("taxonomy") else "",
        "source_qualifiers": str(source_quals),
    }

    for q in [
        "organism", "mol_type", "strain", "isolation_source", "host", "plasmid",
        "chromosome", "geo_loc_name", "collection_date", "country", "isolate",
        "note", "db_xref"
    ]:
        row[f"source_{q}"] = source_quals.get(q, "")

    return row


def metadata_from_record(record: SeqRecord, input_accession: str) -> dict:
    source_feat = get_source_feature(record)
    qualifiers = source_feat.qualifiers if source_feat is not None else {}

    def q(key: str) -> str:
        return join_qualifier(qualifiers.get(key, []))

    return {
        "input_accession": input_accession,
        "record_id": record.id,
        "record_name": record.name,
        "description": record.description,
        "sequence_length_bp": len(record.seq),
        "organism": q("organism"),
        "mol_type": q("mol_type"),
        "strain": q("strain"),
        "isolation_source": q("isolation_source"),
        "host": q("host"),
        "plasmid": q("plasmid"),
        "chromosome": q("chromosome"),
        "geo_loc_name": q("geo_loc_name"),
        "collection_date": q("collection_date"),
        "country": q("country"),
        "isolate": q("isolate"),
        "note": q("note"),
        "db_xref": q("db_xref"),
        "taxonomy": "; ".join(record.annotations.get("taxonomy", [])),
    }


def features_from_record(record: SeqRecord, input_accession: str) -> List[dict]:
    rows = []
    for i, feat in enumerate(record.features, start=1):
        q = feat.qualifiers or {}
        rows.append({
            "input_accession": input_accession,
            "record_id": record.id,
            "feature_index": i,
            "feature_type": feat.type,
            "location": str(feat.location),
            "strand": getattr(feat.location, "strand", None),
            "feature_length_bp": len(feat),
            "gene": join_qualifier(q.get("gene", [])),
            "locus_tag": join_qualifier(q.get("locus_tag", [])),
            "product": join_qualifier(q.get("product", [])),
            "note": join_qualifier(q.get("note", [])),
            "protein_id": join_qualifier(q.get("protein_id", [])),
            "translation": join_qualifier(q.get("translation", [])),
            "db_xref": join_qualifier(q.get("db_xref", [])),
            "codon_start": join_qualifier(q.get("codon_start", [])),
            "pseudo": "yes" if "pseudo" in q else "",
            "all_qualifiers": str({k: join_qualifier(v) for k, v in q.items()}),
        })
    return rows


# -------------------------
# Main
# -------------------------
def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Read FASTA, fetch GenBank data, and write merged Excel output.")
    p.add_argument("--fasta", required=True, help="Input FASTA file")
    p.add_argument("--output", default="final_output.xlsx", help="Output Excel workbook")
    p.add_argument("--email", required=True, help="NCBI Entrez email")
    p.add_argument("--api-key", default="", help="NCBI API key (optional)")
    p.add_argument("--batch-size", type=int, default=100, help="Batch size for Entrez fetch")
    p.add_argument("--sleep", type=float, default=0.34, help="Seconds to sleep between batches")
    p.add_argument("--max-retries", type=int, default=3, help="Retries per request")
    return p


def main() -> None:
    args = build_parser().parse_args()

    fasta_path = Path(args.fasta)
    if not fasta_path.exists():
        raise FileNotFoundError(f"FASTA file not found: {fasta_path}")

    Entrez.email = args.email
    Entrez.tool = "fasta_genbank_pipeline_rowwise"
    if args.api_key.strip():
        Entrez.api_key = args.api_key.strip()

    # -------------------------
    # Read FASTA as-is
    # -------------------------
    fasta_rows = []
    print(f"Reading FASTA: {fasta_path}")
    for idx, rec in enumerate(SeqIO.parse(str(fasta_path), "fasta"), start=1):
        variant, acc = extract_accession_from_fasta_header(rec)
        seq = str(rec.seq).upper()

        fasta_rows.append({
            "row_id": idx,
            "fasta_id": rec.id,
            "fasta_description": rec.description,
            "variant": variant,
            "accession_raw": acc,
            "accession_base": accession_base(acc),
            "sequence_length_bp": len(seq),
            "sequence": seq,
            "sequence_md5": md5_of_sequence(seq),
        })

    fasta_df = pd.DataFrame(fasta_rows)
    if fasta_df.empty:
        raise ValueError("No FASTA records found.")

    fasta_df["accession_count_in_fasta"] = fasta_df.groupby("accession_base")["accession_base"].transform("count")
    fasta_df["distinct_sequence_count_for_accession"] = fasta_df.groupby("accession_base")["sequence_md5"].transform("nunique")
    fasta_df["is_duplicate_accession"] = fasta_df["accession_count_in_fasta"] > 1
    fasta_df["same_sequence_as_other_same_accession"] = fasta_df["distinct_sequence_count_for_accession"] == 1

    accession_report = (
        fasta_df.groupby("accession_base")
        .agg(
            fasta_sequence_count=("accession_base", "size"),
            distinct_sequence_count=("sequence_md5", "nunique"),
            min_sequence_length_bp=("sequence_length_bp", "min"),
            max_sequence_length_bp=("sequence_length_bp", "max"),
        )
        .reset_index()
        .rename(columns={"accession_base": "accession"})
    )
    accession_report["duplicate_sequences"] = accession_report["fasta_sequence_count"] > 1
    accession_report["sequence_variants_present"] = accession_report["distinct_sequence_count"] > 1
    accession_report["duplicate_note"] = accession_report.apply(
        lambda r: (
            "same sequence repeated" if r["duplicate_sequences"] and not r["sequence_variants_present"]
            else "different sequences for same accession" if r["sequence_variants_present"]
            else ""
        ),
        axis=1,
    )

    unique_accessions = (
        fasta_df["accession_base"]
        .replace("", pd.NA)
        .dropna()
        .astype(str)
        .str.upper()
        .unique()
        .tolist()
    )

    print(f"Total FASTA sequences: {len(fasta_df)}")
    print(f"Unique accessions to fetch: {len(unique_accessions)}")
    print("Starting GenBank fetch...\n")

    # -------------------------
    # Fetch GenBank
    # -------------------------
    accession_to_record: Dict[str, SeqRecord] = {}
    failures_rows = []

    batches = chunked(unique_accessions, args.batch_size)
    for batch_no, batch in enumerate(tqdm(batches, desc="Fetching batches", unit="batch"), start=1):
        try:
            batch_records = fetch_batch_records(batch, max_retries=args.max_retries)
            accession_to_record.update(batch_records)
        except Exception as e:
            tqdm.write(f"Batch {batch_no} failed ({e}); falling back to single-accession fetches.")
            for acc in tqdm(batch, desc=f"Fallback batch {batch_no}", leave=False, unit="acc"):
                try:
                    rec = fetch_single_record(acc, max_retries=args.max_retries)
                    accession_to_record[normalize_accession(acc)] = rec
                    accession_to_record[accession_base(acc)] = rec
                    for a in rec.annotations.get("accessions", []):
                        an = normalize_accession(a)
                        accession_to_record[an] = rec
                        accession_to_record[accession_base(an)] = rec
                except Exception as ex:
                    failures_rows.append({
                        "input_accession": acc,
                        "error": str(ex),
                    })

        tqdm.write(
            f"Batch {batch_no}/{len(batches)} done | "
            f"cached={len(accession_to_record)} | "
            f"failures={len(failures_rows)}"
        )
        time.sleep(args.sleep)

    # -------------------------
    # Build data frames
    # -------------------------
    summary_rows: List[dict] = []
    metadata_rows: List[dict] = []
    feature_rows: List[dict] = []

    print("\nParsing fetched records...")
    for acc in tqdm(unique_accessions, desc="Parsing records", unit="acc"):
        rec = accession_to_record.get(accession_base(acc)) or accession_to_record.get(normalize_accession(acc))
        if rec is None:
            failures_rows.append({
                "input_accession": acc,
                "error": "No GenBank record found",
            })
            continue

        summary_rows.append(summary_from_record(rec, input_accession=acc))
        metadata_rows.append(metadata_from_record(rec, input_accession=acc))
        feature_rows.extend(features_from_record(rec, input_accession=acc))

    summary_df = pd.DataFrame(summary_rows)
    metadata_df = pd.DataFrame(metadata_rows)
    features_df = pd.DataFrame(feature_rows)
    failures_df = pd.DataFrame(failures_rows)

    # ==========================================================
    # Expand summary and metadata back to FASTA row level
    # ==========================================================
    if not summary_df.empty:
        summary_df["accession_base"] = summary_df["input_accession"].astype(str).str.upper().str.split(".").str[0]
        summary_expanded = fasta_df.merge(summary_df, on="accession_base", how="left", suffixes=("", "_gb"))
    else:
        summary_expanded = pd.DataFrame()

    if not metadata_df.empty:
        metadata_df["accession_base"] = metadata_df["input_accession"].astype(str).str.upper().str.split(".").str[0]
        metadata_expanded = fasta_df.merge(metadata_df, on="accession_base", how="left", suffixes=("", "_gb"))
    else:
        metadata_expanded = pd.DataFrame()

    # -------------------------
    # Final consolidated data
    # -------------------------
    if not summary_df.empty:
        final_df = fasta_df.merge(summary_df, how="left", on="accession_base", suffixes=("", "_gb"))
    else:
        final_df = fasta_df.copy()

    final_df["fetch_status"] = final_df["record_id"].apply(lambda x: "found" if pd.notna(x) and str(x).strip() else "missing")

    preferred_front = [
        "row_id", "fasta_id", "variant", "accession_raw", "accession_base",
        "sequence_length_bp", "accession_count_in_fasta",
        "distinct_sequence_count_for_accession", "is_duplicate_accession",
        "same_sequence_as_other_same_accession", "fetch_status",
        "record_id", "record_name", "organism", "replicon_type", "plasmid", "chromosome",
        "country", "host", "strain", "geo_loc_name", "collection_date", "note"
    ]
    preferred_front = [c for c in preferred_front if c in final_df.columns]
    remaining = [c for c in final_df.columns if c not in preferred_front]
    final_df = final_df[preferred_front + remaining]

    # ==========================================================
    # Advanced Statistical Summaries
    # ==========================================================
    # Count replicons based on final expanded data rows
    if "replicon_type" in summary_expanded.columns:
        replicon_counts = summary_expanded["replicon_type"].fillna("unknown").value_counts()
        plasmid_count = replicon_counts.get("plasmid", 0)
        chromosome_count = replicon_counts.get("chromosome", 0)
        unknown_count = replicon_counts.get("unknown", 0)
    else:
        plasmid_count = chromosome_count = unknown_count = 0

    # Calculate duplicate statistics
    unique_accessions_with_duplicates = (accession_report["fasta_sequence_count"] > 1).sum()
    total_duplicate_fasta_rows = fasta_df["is_duplicate_accession"].sum()

    # Build Run Summary metrics dataframe
    summary_metrics_data = [
        {"Metric Category": "FASTA Statistics", "Metric Name": "Total Sequences (Rows)", "Count Value": len(fasta_df)},
        {"Metric Category": "FASTA Statistics", "Metric Name": "Unique Accession Bases", "Count Value": len(unique_accessions)},
        {"Metric Category": "Duplicate Statistics", "Metric Name": "Unique Accessions with Duplicates", "Count Value": int(unique_accessions_with_duplicates)},
        {"Metric Category": "Duplicate Statistics", "Metric Name": "Total Redundant/Duplicate Rows", "Count Value": int(total_duplicate_fasta_rows)},
        {"Metric Category": "GenBank Replicon Distribution", "Metric Name": "Plasmid Replicons", "Count Value": int(plasmid_count)},
        {"Metric Category": "GenBank Replicon Distribution", "Metric Name": "Chromosome Replicons", "Count Value": int(chromosome_count)},
        {"Metric Category": "GenBank Replicon Distribution", "Metric Name": "Unknown/Unclassified Replicons", "Count Value": int(unknown_count)},
    ]
    summary_metrics_df = pd.DataFrame(summary_metrics_data)

    # -------------------------
    # Write workbook
    # -------------------------
    with pd.ExcelWriter(args.output, engine="openpyxl") as writer:
        summary_metrics_df.to_excel(writer, index=False, sheet_name="run_summary")
        fasta_df.to_excel(writer, index=False, sheet_name="fasta_input")
        accession_report.to_excel(writer, index=False, sheet_name="accession_report")
        summary_expanded.to_excel(writer, index=False, sheet_name="genbank_summary")
        metadata_expanded.to_excel(writer, index=False, sheet_name="metadata")
        features_df.to_excel(writer, index=False, sheet_name="genbank_features")
        final_df.to_excel(writer, index=False, sheet_name="final_data")
        if not failures_df.empty:
            failures_df.to_excel(writer, index=False, sheet_name="fetch_failures")

    # -------------------------
    # Comprehensive Console Report
    # -------------------------
    print("\n" + "="*50)
    print("                 PIPELINE RUN SUMMARY  ")
    print("="*50)
    print(f"FASTA Sequences Processed : {len(fasta_df)}")
    print(f"Unique Accessions Fetched : {len(unique_accessions)}")
    print(f"Total Summary Rows        : {len(summary_expanded)}")
    print(f"Total Metadata Rows       : {len(metadata_expanded)}")
    print(f"Feature Rows Parsed       : {len(features_df)}")
    print(f"Failures Encountered      : {len(failures_df)}")
    print("-"*50)
    print("REPLICON DISTRIBUTION (Expanded Row Level):")
    print(f"  - Chromosomes           : {chromosome_count}")
    print(f"  - Plasmids              : {plasmid_count}")
    print(f"  - Unknown               : {unknown_count}")
    print("-"*50)
    print("DUPLICATE TRACKING:")
    print(f"  - Unique Accessions featuring duplicates : {unique_accessions_with_duplicates}")
    print(f"  - Total Row Footprint of duplicate items : {total_duplicate_fasta_rows}")
    print("="*50)
    print(f"Output saved safely to: {Path(args.output).resolve()}\n")


if __name__ == "__main__":
    main()