# 🦪 Oyster Host & Microbiome Transcriptomics (for manuscript reproducibility)

## Overview: raw-data processing + WGCNA analysis
This repository provides the materials necessary to reproduce the core computational workflows used in our study of Crassostrea virginica host and microbiome transcriptomes. Specifically, it contains (1) the complete bioinformatics pipeline used to generate raw count matrices and annotations for both the host and the microbiome, and (2) the full weighted gene co-expression network analysis (WGCNA) performed on filtered, log transformed expression matrices. These resources ensure transparency and reproducibility of the primary data processing and co-expression analyses associated with the manuscript.

Raw reads for this project are publicly available under NCBI BioProject accession  
[PRJNA1313129](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1313129).

---

## 📁 Repository contents

### Bioinformatics workflow
- `1_raw_reads_to_raw_counts/pipeline.md`  
  End-to-end processing pipeline (QC → assembly/mapping → counting → annotation), including software versions and parameters used to generate the host and microbiome count matrices.

---

### Raw count matrices
- `1_raw_reads_to_raw_counts/host_counts.tsv[.gz]` — raw host transcript counts  
- `1_raw_reads_to_raw_counts/microbiome_counts.tsv[.gz]` — raw microbiome ORF counts  

These matrices allow users to apply alternative normalization or additional downstream analyses if desired.

---

### Supporting scripts
- `1_raw_reads_to_raw_counts/make_microbiome_counts.py` — generates ORF-level microbiome count tables  
- `1_raw_reads_to_raw_counts/emapper_to_tsv.py` — converts eggNOG-mapper annotation outputs to TSV format  

---

## Host transcriptomic analyses

### Differential expression and biomineralization toolkit genes
- `2_downstream_analyses/2.2_Host_analyses/  
  Scripts in this folder integrate host DESeq2 results across timepoints, joins transcript-level GO annotations, and identifies differentially expressed host genes that belong to a curated molluscan biomineralization toolkit. This analysis links statistically significant host responses to known biomineralization-related molecular functions described in the literature.

---

## Microbiome transcriptomic analyses

### Differential expression and biomineralization toolkit genes
- `2_downstream_analyses/2.3_Microbiome_Analyses/  
  Scripts in 2_downstream_analyses/2.3_Microbiome_Analyses/ process microbiome metatranscriptomic data from per-sample ORF count tables and eggNOG-mapper annotations through quality control, filtering, and functional aggregation. These analyses retain ORFs with both expression and functional annotations, remove eukaryotic and ribosomal transcripts, and construct phyloseq objects for downstream analyses. ORF-level counts are collapsed to KEGG orthologues (KOs), filtered for rarity and prevalence, and analyzed using DESeq2 to identify differentially abundant microbial functions across treatments and timepoints. Targeted downstream scripts then extract KOs involved in biochemical pathways relevant to calcifying fluid chemistry (e.g., nitrogen and sulfur cycling, carbonic anhydrase, urease) and generate tornado plots summarizing effect sizes across timepoints.
 
---


### WGCNA analysis
- `2.4_Host_Microbiome_WGCNA/WGCNA_markdown.Rmd`  
  R Markdown document containing the complete co-expression analysis for host and microbiome datasets.
- `2.4_Host_Microbiome_WGCNA/WGCNA_markdown.html`  
  Rendered HTML version for easy browser-based viewing.
- `2.4_Host_Microbiome_WGCNA/host_counts_WGCNA.csv`  
  Log-normalized and filtered host expression matrix used as input for WGCNA.
- `2.4_Host_Microbiome_WGCNA/microbiome_counts_WGCNA.csv`  
  Log-normalized and filtered microbiome expression matrix used as input for WGCNA.
- `2.4_Host_Microbiome_WGCNA/sample_traits_WGCNA.csv`  
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
