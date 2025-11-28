# 🦪 Oyster Host & Microbiome Transcriptomics (for manuscript reproducibility)

## Overview: Reproducible raw-data processing + WGCNA analysis
This repository provides the materials necessary to reproduce the core computational workflows used in our study of Crassostrea virginica host and microbiome transcriptomes. Specifically, it contains (1) the complete bioinformatics pipeline used to generate raw count matrices and annotations for both the host and the microbiome, and (2) the full weighted gene co-expression network analysis (WGCNA) performed on filtered, log transformed expression matrices. These resources ensure transparency and reproducibility of the primary data processing and co-expression analyses associated with the manuscript.

Raw reads for this project are publicly available under NCBI BioProject accession  
[PRJNA1313129](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1313129).

---

## 📁 Repository contents

### Bioinformatics workflow
- `docs/pipeline.md`  
  End-to-end processing pipeline (QC → assembly/mapping → counting → annotation), including software versions and parameters used to generate the host and microbiome count matrices.

---

### Raw count matrices
- `data/counts/host_counts.tsv[.gz]` — raw host transcript counts  
- `data/counts/microbiome_counts.tsv[.gz]` — raw microbiome ORF counts  

These matrices allow users to apply alternative normalization or additional downstream analyses if desired.

---

### Supporting scripts
- `scripts/make_microbiome_counts.py` — generates ORF-level microbiome count tables  
- `scripts/emapper_to_tsv.py` — converts eggNOG-mapper annotation outputs to TSV format  

---

### WGCNA analysis
- `docs/WGCNA_markdown.Rmd`  
  R Markdown document containing the complete co-expression analysis for host and microbiome datasets.
- `docs/WGCNA_markdown.html`  
  Rendered HTML version for easy browser-based viewing.
- `data/counts/host_counts_WGCNA.csv`  
  Log-normalized and filtered host expression matrix used as input for WGCNA.
- `data/counts/microbiome_counts_WGCNA.csv`  
  Log-normalized and filtered microbiome expression matrix used as input for WGCNA.
- `data/counts/sample_traits_WGCNA.csv`  
  Sample trait metadata file used for module–trait relationships in the WGCNA workflow.

> All inputs needed to fully reproduce the WGCNA analysis are included in the repository.

---

## 🚀 Quick start: Load count matrices in R

```r
# Raw count matrices
host <- read.delim("data/counts/host_counts.tsv", check.names = FALSE)
micro <- read.delim("data/counts/microbiome_counts.tsv", check.names = FALSE)

# WGCNA inputs
host_wgcna <- read.csv("data/counts/host_counts_WGCNA.csv", check.names = FALSE, row.names = 1)
micro_wgcna <- read.csv("data/counts/microbiome_counts_WGCNA.csv", check.names = FALSE, row.names = 1)
traits_wgcna <- read.csv("data/counts/sample_traits_WGCNA.csv", check.names = FALSE, row.names = 1)
