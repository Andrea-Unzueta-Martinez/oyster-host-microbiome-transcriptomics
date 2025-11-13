# Oyster host & microbiome transcript counts (for manuscript reproducibility)

## Overview
This repository provides the count matrices and scripts used to generate the oyster (*Crassostrea virginica*) host and microbiome transcriptomic datasets.  
These files allow others to reproduce downstream statistical analyses (e.g., differential expression, co-expression network analysis) without reprocessing the raw reads.

Raw reads for this project are publicly available under NCBI BioProject accession  
[PRJNA1313129](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1313129).

## What’s here
- `data/counts/host_counts.tsv[.gz]` — host transcript counts (rows = transcripts, columns = samples)  
- `data/counts/microbiome_counts.tsv[.gz]` — microbiome ORF counts (rows = ORFs, columns = samples)  
- `scripts/make_microbiome_counts.py` — Python script for generating per-sample microbiome ORF counts  
- `scripts/emapper_to_tsv.py` — script for converting eggNOG-mapper annotation files to tab-delimited TSV format  
- `docs/pipeline.md` — full step-by-step workflow (QC → mapping → counting → integration) including software versions and parameters  

## Data description
- Each count matrix contains **raw, unnormalized counts** produced after merging all samples.  
- For downstream analyses, counts were normalized in R to **RPKM (Reads Per Kilobase per Million mapped reads)**.  
- The uploaded raw matrices allow users to apply alternative normalization or transformation methods (e.g., TPM, DESeq2 size-factor normalization).  

## Quick start
Load count matrices in R:

```r
# R example
host <- read.delim("data/counts/host_counts.tsv", check.names = FALSE)
micro <- read.delim("data/counts/microbiome_counts.tsv", check.names = FALSE)

# or if gzipped:
host <- read.delim(gzfile("data/counts/host_counts.tsv.gz"), check.names = FALSE)
micro <- read.delim(gzfile("data/counts/microbiome_counts.tsv.gz"), check.names = FALSE)
