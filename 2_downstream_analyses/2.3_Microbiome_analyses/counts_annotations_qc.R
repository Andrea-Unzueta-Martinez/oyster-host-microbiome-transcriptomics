# ------------------------------------------------------------
# Script: counts_annotations_qc_phyloseq.R
#
# Purpose:
#   Process microbiome transcriptomic data from per-sample outputs
#   into analysis-ready objects. This script:
#     (1) loads ORF-level count tables and eggNOG-mapper annotations,
#     (2) performs quality control to retain ORFs with both counts
#         and functional annotations,
#     (3) merges all samples into master counts and annotation tables,
#     (4) aligns sample metadata and constructs phyloseq objects,
#     (5) applies taxonomic and functional filtering, and
#     (6) normalizes transcript abundances (RPKM).
#
# Outputs:
#   - transcripts_phy        : raw ORF-level phyloseq object
#   - transcripts_phy_filt2  : filtered (non-eukaryotic, non-ribosomal) phyloseq object
#   - RPKM_phy               : RPKM-normalized phyloseq object
# ------------------------------------------------------------


#load packages
library(readr)
library(readxl)
library(dplyr)
library(stringr)
library(tidyr)
library(purrr)
library(data.table) 
library(tibble)

# ---- paths ----
counts_and_annotations_dir <- "/Users/andreaunzueta/Dropbox/NSF_PRFB_GirguisLab/Bioinformatics/Metatranscriptomics_analysis/transcriptomics/counts_and_annotations"  # directory where counts tables are for each sample

# ---- 1) LOAD eggNOG annotation tables ----

#extract smaple IDs
extract_sample_id <- function(x) {
  # Capture SAMPLE, optionally followed by _S<number>, right before _non_oyster_transcripts_
  m <- str_match(x, "_([^_]+)(?:_S\\d+)?_non_oyster_transcripts_")
  if (!is.na(m[1,2])) return(m[1,2])
  
  stop("Could not extract sample_id from filename: ", x)
}

# Identify all eggNOG-mapper annotation files (one per sample)
anno_files <- list.files(
  counts_and_annotations_dir,
  pattern = "\\.emapper\\.annotations\\.tsv$",
  full.names = TRUE
)

#check sample names are correct
sample_names_df <- data.frame(
  file = basename(anno_files),
  sample_id = vapply(basename(anno_files), extract_sample_id, character(1)),
  row.names = NULL
)

head(sample_names_df)

#load annotations for each sample
walk(anno_files, function(f) {
  samp <- extract_sample_id(basename(f))
  obj  <- paste0("ORF_annotations_", samp)
  
  assign(
    x = obj,
    value = as.data.table(read_tsv(f, show_col_types = FALSE)),
    envir = .GlobalEnv
  )
})

#check they all loaded
ls(pattern = "^ORF_annotations_")


### ---- 2) LOAD counts tables ----
# This extractor works for your counts filenames: "62A_ORF_counts.csv", "37_ORF_counts.csv", "neg3_ORF_counts.csv"
extract_sample_id_counts <- function(x) {
  m <- str_match(x, "^([^_]+)_ORF_counts\\.csv$")
  if (!is.na(m[1,2])) return(m[1,2])
  stop("Could not extract sample_id from counts filename: ", x)
}

# Standardize a counts file to: orf_id + <sample>_counts
read_orf_counts <- function(filepath) {
  fname <- basename(filepath)
  samp  <- extract_sample_id_counts(fname)
  
  df <- read_csv(filepath, show_col_types = FALSE)
  
  if (ncol(df) < 2) stop("Counts file has <2 columns: ", filepath)
  
  # Find ORF id column (prefer something containing 'orf', else first column)
  orf_col <- names(df)[str_detect(names(df), regex("orf", ignore_case = TRUE))]
  orf_col <- if (length(orf_col) >= 1) orf_col[1] else names(df)[1]
  
  # Count column = first non-orf column (your files appear 2-col; this is safe + simple)
  count_col <- setdiff(names(df), orf_col)[1]
  
  out <- df %>%
    transmute(
      orf_id = .data[[orf_col]],
      counts = as.integer(.data[[count_col]])
    )
  
  # Rename "counts" to "<sample>_counts"
  setnames(as.data.table(out), "counts", paste0(samp, "_counts"))
}

# ---- load all counts files into individual objects ----
count_files <- list.files(
  counts_and_annotations_dir,
  pattern = "_ORF_counts\\.csv$",
  full.names = TRUE
)

walk(count_files, function(f) {
  samp <- extract_sample_id_counts(basename(f))
  obj  <- paste0("ORF_counts_", samp)
  
  assign(
    x = obj,
    value = read_orf_counts(f),
    envir = .GlobalEnv
  )
})

# quick check: see what loaded
ls(pattern = "^ORF_counts_")

##############################################################
# Quality control: match ORFs between counts and annotations
#
# For each sample, we retain only ORFs that:
#   (1) have non-zero expression counts AND
#   (2) have functional annotations from eggNOG-mapper
#
# This ensures a one-to-one correspondence between expression
# and functional data for all downstream analyses.
##############################################################

# Identify samples with both counts and annotation tables
count_objs <- ls(pattern = "^ORF_counts_")
anno_objs  <- ls(pattern = "^ORF_annotations_")

count_ids <- str_remove(count_objs, "^ORF_counts_")
anno_ids  <- str_remove(anno_objs,  "^ORF_annotations_")

sample_ids <- intersect(count_ids, anno_ids)

# QC function: subset to shared ORFs only
qc_match_orfs <- function(sample_id) {
  
  counts <- get(paste0("ORF_counts_", sample_id), envir = .GlobalEnv)
  anno   <- get(paste0("ORF_annotations_", sample_id), envir = .GlobalEnv)
  
  # Defensive checks
  stopifnot("orf_id" %in% names(counts))
  stopifnot("orf_id" %in% names(anno))
  
  shared_orfs <- intersect(counts$orf_id, anno$orf_id)
  
  counts_qc <- counts %>% filter(orf_id %in% shared_orfs)
  anno_qc   <- anno   %>% filter(orf_id %in% shared_orfs)
  
  # Save QC-passed objects
  assign(paste0("ORF_counts_qc_", sample_id), counts_qc, envir = .GlobalEnv)
  assign(paste0("ORF_annotations_qc_", sample_id), anno_qc, envir = .GlobalEnv)
  
  # Return QC metrics for reporting
  tibble(
    sample_id = sample_id,
    counts_total = nrow(counts),
    annotations_total = nrow(anno),
    orfs_retained = length(shared_orfs),
    counts_removed = nrow(counts) - length(shared_orfs),
    annotations_removed = nrow(anno) - length(shared_orfs)
  )
}

# Run QC across all samples
orf_qc_summary <- map_dfr(sample_ids, qc_match_orfs)

orf_qc_summary



############################################################
# Build master tables across ALL samples:
#   (1) ORF counts matrix (rows = ORFs, cols = samples)
#   (2) ORF annotations table (rows = ORFs, cols = annotation fields)
# Then keep only ORFs present in BOTH tables and align row order.
############################################################

# Collect the per-sample QC objects made earlier
anno_objs  <- ls(pattern = "^ORF_annotations_qc_")
count_objs <- ls(pattern = "^ORF_counts_qc_")

anno_list  <- setNames(lapply(anno_objs,  get, envir = .GlobalEnv),
                       str_remove(anno_objs,  "^ORF_annotations_qc_"))
count_list <- setNames(lapply(count_objs, get, envir = .GlobalEnv),
                       str_remove(count_objs, "^ORF_counts_qc_"))

# ---- 1) Master annotations table (ORF x annotation fields) ----
# Stack all annotation rows, then keep a single row per ORF.
# (eggNOG annotations are ORF-level, not sample-level)
allsamples_ORF_table_annotations <- bind_rows(anno_list) %>%
  distinct(orf_id, .keep_all = TRUE) %>%
  column_to_rownames("orf_id")

# ---- 2) Master counts table (ORF x samples) ----
# Full-join all sample count columns by orf_id; fill missing with 0.
allsamples_ORF_counts_table <- reduce(count_list, full_join, by = "orf_id") %>%
  mutate(across(where(is.numeric), ~ replace_na(.x, 0))) %>%
  column_to_rownames("orf_id")

# ---- Align ORFs between tables (intersection + same order) ----
shared_ids <- intersect(rownames(allsamples_ORF_table_annotations),
                        rownames(allsamples_ORF_counts_table))

allsamples_ORF_table_annotations <- allsamples_ORF_table_annotations[shared_ids, , drop = FALSE]
allsamples_ORF_counts_table      <- allsamples_ORF_counts_table[shared_ids, , drop = FALSE]

stopifnot(identical(rownames(allsamples_ORF_table_annotations),
                    rownames(allsamples_ORF_counts_table)))


library(dplyr)
library(tibble)
library(phyloseq)

############################################################
# Load sample metadata + build phyloseq object
############################################################

# ---- 1) Read metadata ----
meta_table <- read_csv("/Users/andreaunzueta/Dropbox/NSF_PRFB_GirguisLab/Bioinformatics/Metatranscriptomics_analysis/transcriptomics/Transcriptomics_metadata.csv", show_col_types = FALSE) %>%
  as.data.frame()

# ---- 2) Standardize sample IDs in the counts matrix ----
# Merge operations can add suffixes (.x/.y); counts files often include "_counts"
colnames(allsamples_ORF_counts_table) <-
  colnames(allsamples_ORF_counts_table) |>
  str_remove("\\.x$") |>
  str_remove("\\.y$") |>
  str_remove("_counts$")

head(colnames(allsamples_ORF_counts_table)) #check

# ---- 3) Standardize sample IDs in metadata to match counts ----
# Extract the sample token at the end of Sample_ID (e.g., e_6h_62A -> 62A)
meta_table_filt <- meta_table_filt %>%
  mutate(sample_id = str_extract(Sample_ID, "(neg\\d+|mock\\d+|\\d+[A-Z]?)$")) %>%
  as.data.frame()

rownames(meta_table_filt) <- meta_table_filt$sample_id

head(rownames(meta_table_filt))

# ---- 4) Keep only samples present in BOTH metadata and counts; align order ----
shared_samples <- intersect(
  colnames(allsamples_ORF_counts_table),
  rownames(meta_table_filt)
)

# subset and reorder explicitly
allsamples_ORF_counts_table <- allsamples_ORF_counts_table[, shared_samples, drop = FALSE]
meta_table_filt <- meta_table_filt[shared_samples, , drop = FALSE]

# QC: enforce perfect alignment (prevents silent sample mis-labeling)
stopifnot(
  identical(colnames(allsamples_ORF_counts_table),
            rownames(meta_table_filt))
)


# ---- 5) Build phyloseq object ----
# phyloseq expects matrices; tax_table works best with character data
library(phyloseq)

allsamples_ORF_table_annotations[] <-
  lapply(allsamples_ORF_table_annotations, as.character)

otu <- otu_table(as.matrix(allsamples_ORF_counts_table), taxa_are_rows = TRUE). # ORFs x samples
tax <- tax_table(as.matrix(allsamples_ORF_table_annotations))                   # ORFs x annotation fields
sam <- sample_data(meta_table_filt)                                             # samples x metadata fields

transcripts_phy <- phyloseq(otu, tax, sam)
transcripts_phy

############################################################
# QC filtering: remove eukaryotic transcripts and ribosomal ORFs
############################################################
library(phyloseq)
library(dplyr)

# Define taxonomic groups to exclude (based on max_annot_lvl)
exclude_taxa <- c("2759|Eukaryota", "33208|Metazoa", "4751|Fungi")

# Drop eukaryote-associated ORFs
transcripts_phy_filt <- subset_taxa(transcripts_phy, !(max_annot_lvl %in% exclude_taxa))

# QC: how many ORFs removed?
n_total <- ntaxa(transcripts_phy)
n_removed_tax <- ntaxa(subset_taxa(transcripts_phy, max_annot_lvl %in% exclude_taxa))
pct_removed_tax <- 100 * n_removed_tax / n_total

qc_tax_filter <- tibble(
  total_orfs = n_total,
  removed_euk_taxa = n_removed_tax,
  pct_removed_euk_taxa = pct_removed_tax,
  remaining_orfs = ntaxa(transcripts_phy_filt)
)
qc_tax_filter

# Optional QC: inspect remaining max_annot_lvl values
# unique(tax_table(transcripts_phy_filt)[, "max_annot_lvl"])

# Identify ribosomal modules to exclude
ribosomal_modules <- c("M00178", "M00179", "M00177")

n_ribo <- ntaxa(subset_taxa(transcripts_phy_filt, KEGG_Module %in% ribosomal_modules))
pct_ribo <- 100 * n_ribo / ntaxa(transcripts_phy_filt)

qc_ribosome <- tibble(
  remaining_orfs_post_tax_filter = ntaxa(transcripts_phy_filt),
  ribosomal_orfs = n_ribo,
  pct_ribosomal_orfs = pct_ribo
)
qc_ribosome

# Remove ribosomal ORFs
transcripts_phy_filt2 <- subset_taxa(transcripts_phy_filt, !(KEGG_Module %in% ribosomal_modules))
transcripts_phy_filt2


############################################################
# RPKM normalization (matrix-based; fast + reproducible)
############################################################
library(stringr)
library(phyloseq)

# Extract count matrix (ORFs x samples)
counts_mat <- as(otu_table(transcripts_phy_filt2), "matrix")
if (!taxa_are_rows(transcripts_phy_filt2)) counts_mat <- t(counts_mat)

# Parse gene length from ORF IDs like: NODE_15083_length_940_cov_...
orf_ids <- rownames(counts_mat)

gene_length <- str_match(orf_ids, "length_(\\d+)")[, 2]
if (any(is.na(gene_length))) {
  stop("Could not parse gene length for some ORF IDs. Example: ",
       paste(head(orf_ids[is.na(gene_length)]), collapse = ", "))
}
gene_length <- as.numeric(gene_length)

# Total reads per sample (library size)
total_reads <- colSums(counts_mat)

# Compute RPKM: counts / (length_kb * total_reads_million)
length_kb <- gene_length / 1000
reads_million <- total_reads / 1e6

rpkm_mat <- sweep(counts_mat, 1, length_kb, "/")
rpkm_mat <- sweep(rpkm_mat, 2, reads_million, "/")

# Build a new phyloseq object with RPKM abundances
RPKM_phy <- phyloseq(
  otu_table(rpkm_mat, taxa_are_rows = TRUE),
  tax_table(transcripts_phy_filt2),
  sample_data(transcripts_phy_filt2)
)

RPKM_phy
