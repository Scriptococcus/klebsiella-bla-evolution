# Comparative Genomic Architecture of *Klebsiella pneumoniae* β-Lactamase Loci

Computational workflows for comparative analysis of β-lactamase genes and their genomic neighbourhoods in *Klebsiella pneumoniae*, including phylogenetic, recombination, and mobile genetic element analyses.

This repository accompanies the study:

**β-Lactamases and Their Genomic Neighbourhoods in *Klebsiella pneumoniae* Exhibit a Partially Modular Evolutionary Architecture**

## Data availability

This repository contains the analysis scripts and software-environment specifications.

## Repository structure

```text
.
├── README.md
├── CITATION.cff
├── environment.yml
├── environment-r.yml
├── requirements.txt
├── .gitignore
└── scripts/
    ├── 01_card_tblastn_extract.py
    ├── 02_fasta_genbank_metadata.py
    ├── 03_variant_counts.py
    ├── 04_variant_distribution_plot.py
    ├── 05_shorten_fasta_headers.py
    ├── 06_generate_itol_annotations.py
    ├── 07_generate_shv_itol_annotations.py
    ├── 08_patristic_clustering_networks.R
    ├── 09_topological_discordance.R
    ├── 10_topological_discordance_plot.R
    ├── 11_rdp_breakpoint_maps.py
    ├── 12_genomic_feature_map.py
    ├── 13_compile_genomic_feature_maps.py
    └── 14_extract_rdp_detection_methods.py
```

## Workflow

The workflow covers sequence identification and extraction, metadata recovery, phylogenetic reconstruction, regional clustering, topological comparison, recombination analysis, and genomic-context visualisation.

```text
CARD protein targets
        │
        ▼
01_card_tblastn_extract.py
        │
        ├── CDS
        ├── 1-kb upstream region
        └── 1-kb downstream region
        │
        ▼
02_fasta_genbank_metadata.py
        │
        ├── GenBank metadata
        └── genomic feature information
        │
        ├── 03_variant_counts.py
        ├── 04_variant_distribution_plot.py
        └── 05_shorten_fasta_headers.py
        │
        ▼
     IQ-TREE
        │
        └── maximum-likelihood trees
            for upstream / CDS / downstream
                    │
                    ├── 06_generate_itol_annotations.py
                    ├── 07_generate_shv_itol_annotations.py
                    └── 08_patristic_clustering_networks.R
                            │
                            ├── patristic-distance clustering
                            └── regional cluster comparison
                    │
                    ├── 09_topological_discordance.R
                    ├── 10_topological_discordance_plot.R
                    ├── 11_rdp_breakpoint_maps.py
                    ├── 14_extract_rdp_detection_methods.py
                    └── 12_genomic_feature_map.py
                         13_compile_genomic_feature_maps.py
```

## Software requirements

### Python

A Conda environment is provided:

```bash
conda env create -f environment.yml
conda activate bla-regulatory-python
```

Alternatively:

```bash
pip install -r requirements.txt
```

### R

A separate Conda environment is provided:

```bash
conda env create -f environment-r.yml
conda activate bla-regulatory-r
```

If necessary:

```bash
conda config --add channels conda-forge
conda config --set channel_priority strict
```

### IQ-TREE

IQ-TREE is used to generate the maximum-likelihood phylogenies that serve as inputs to the downstream analyses. IQ-TREE should be installed separately when reproducing the tree-building step.

## Input data

Download the required datasets from the linked Zenodo record before running the workflows.

Depending on the analysis, inputs include:

- FASTA sequence datasets
- GenBank metadata
- maximum-likelihood tree files
- RDP4 CSV exports
- matched-taxon bootstrap data
- SnapGene `.dna` records

For the patristic-distance workflow, the phylogenetic trees should follow this naming convention:

```text
imp_upstream.treefile
imp_cds.treefile
imp_downstream.treefile

kpc_upstream.treefile
kpc_cds.treefile
kpc_downstream.treefile

ndm_upstream.treefile
ndm_cds.treefile
ndm_downstream.treefile

tem_upstream.treefile
tem_cds.treefile
tem_downstream.treefile

shv_upstream.treefile
shv_cds.treefile
shv_downstream.treefile
```

## Analysis scripts

### 1. CARD → tBLASTn → regional sequence extraction

`01_card_tblastn_extract.py` retrieves CARD protein targets, identifies 100%-identity tBLASTn matches in *K. pneumoniae*, and extracts the CDS together with the upstream, promoter, and downstream regions.

An NCBI Entrez email is required; an API key is optional.

```bash
python scripts/01_card_tblastn_extract.py
```

### 2. FASTA → GenBank metadata

```bash
python scripts/02_fasta_genbank_metadata.py \
  --fasta path/to/input.fasta \
  --output path/to/genbank_metadata.xlsx \
  --email your_email@example.com \
  --api-key YOUR_NCBI_API_KEY
```

### 3. Variant counts and distribution

```bash
python scripts/03_variant_counts.py \
  --input path/to/variant_sequences.fasta \
  --output path/to/variant_counts.csv
```

```bash
python scripts/04_variant_distribution_plot.py \
  --inputs imp.csv kpc.csv ndm.csv tem.csv shv.csv \
  --families IMP KPC NDM TEM SHV \
  --output path/to/variant_distribution.png
```

`05_shorten_fasta_headers.py` provides a FASTA-header utility for downstream identifier matching.

### 4. iTOL annotations

Generate standard iTOL annotations:

```bash
python scripts/06_generate_itol_annotations.py \
  --tree path/to/treefile \
  --metadata path/to/metadata.xlsx \
  --gene IMP \
  --id_col fasta_id
```

Generate SHV-specific annotations:

```bash
python scripts/07_generate_shv_itol_annotations.py \
  --tree path/to/shv.treefile \
  --metadata path/to/metadata.xlsx \
  --gene SHV \
  --id_col accession
```

### 5. Patristic-distance clustering and regional comparison

```bash
Rscript scripts/08_patristic_clustering_networks.R
```

Pairwise patristic distances are calculated directly from the maximum-likelihood trees, and clustering is performed on the full distance matrices.

PAM is used for datasets containing ≤300 sequences, while average-linkage hierarchical clustering is used for larger datasets. Candidate values of `k` range from 2 to `min(15, floor(n/15))`, and the solution with the highest average silhouette width is selected.

For datasets containing >1,000 sequences, silhouette evaluation uses a fixed random subset of 1,000 sequences (seed = 42), while final hierarchical clustering assignments are obtained from the full dataset.

The workflow generates regional cluster assignments, cluster-concordance metrics, and the cluster-level visualisations reported in the study.

### 6. Topological discordance

```bash
Rscript scripts/09_topological_discordance.R \
  --tree-dir path/to/trees \
  --output-dir path/to/output
```

Generate the bootstrap discordance plot:

```bash
Rscript scripts/10_topological_discordance_plot.R \
  --input path/to/matched_bootstrap_replicates.csv \
  --output-dir path/to/output
```

The analysis evaluates topological discordance using normalized clustering-information distance, normalized Robinson–Foulds distance, quartet distance, permutation testing, matched-taxon bootstrap comparisons, and tanglegram visualisation.

### 7. RDP4 breakpoint analysis

Generate upstream breakpoint/event maps:

```bash
python scripts/11_rdp_breakpoint_maps.py \
  --input-dir path/to/RDP4/CSV_FILES \
  --region upstream
```

Generate downstream breakpoint/event maps:

```bash
python scripts/11_rdp_breakpoint_maps.py \
  --input-dir path/to/RDP4/CSV_FILES \
  --region downstream
```

Extract recombination-detection method information from an RDP4 export:

```bash
python scripts/14_extract_rdp_detection_methods.py \
  --input path/to/rdp4_export.csv \
  --output path/to/rdp_detection_methods.xlsx
```

### 8. Genomic-context visualisation

Convert a SnapGene `.dna` record into a publication-style genomic feature map:

```bash
python scripts/12_genomic_feature_map.py \
  path/to/record.dna \
  path/to/output
```

Compile existing PNG maps:

```bash
python scripts/13_compile_genomic_feature_maps.py \
  path/to/generated_maps \
  path/to/output
```

## Reproducibility

The datasets required for the analyses are archived in Zenodo, while the analysis scripts and software-environment specifications are provided in this repository.

Use relative paths or command-line arguments rather than machine-specific paths.

Do not commit NCBI API keys, credentials, private datasets, or other sensitive configuration files.

## Citation

Please cite the associated publication when using this repository:

**β-Lactamases and Their Genomic Neighbourhoods in *Klebsiella pneumoniae* Exhibit a Partially Modular Evolutionary Architecture**

**Code:** https://github.com/Scriptococcus/klebsiella-bla-evolution

See `CITATION.cff` for citation metadata.
