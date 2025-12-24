
# ------------------------------------------------------------
# Microbiome differential abundance (DESeq2) at KEGG KO level
# Input:  transcripts_phy_filt2 (ORF-level phyloseq, raw counts)
# Output: KO-level DESeq2 object + normalized KO count matrix
# ------------------------------------------------------------

library(phyloseq)
library(dplyr)
library(stringr)
library(tidyr)
library(DESeq2)
library(tibble)
library(readr)

#-----------------------------#
# Helpers
#-----------------------------#

filter_phyloseq_samples <- function(ps,
                                    exclude_treatment = "procedural_control",
                                    min_library_size = 1000,
                                    min_taxa_sum = 1) {
  
  samp <- as(sample_data(ps), "data.frame")
  
  keep_samples <- rownames(samp)[
    samp$Treatment != exclude_treatment &
      sample_sums(ps) > min_library_size
  ]
  
  ps1 <- prune_samples(keep_samples, ps)
  ps2 <- prune_taxa(taxa_sums(ps1) > min_taxa_sum, ps1)
  
  ps2
}

collapse_orfs_to_kos <- function(ps, ko_col = "KEGG_ko", sep = ",") {
  
  stopifnot(ko_col %in% colnames(tax_table(ps)))
  
  counts <- as(otu_table(ps), "matrix")
  if (!taxa_are_rows(ps)) counts <- t(counts)
  
  tax <- as.data.frame(tax_table(ps))
  ko_raw <- tax[[ko_col]]
  names(ko_raw) <- rownames(tax)
  
  keep_orfs <- !is.na(ko_raw) & ko_raw != ""
  counts <- counts[keep_orfs, , drop = FALSE]
  ko_raw <- ko_raw[keep_orfs]
  
  map_df <- tibble(orf_id = names(ko_raw), KO = ko_raw) %>%
    separate_rows(KO, sep = sep) %>%
    mutate(KO = str_trim(KO)) %>%
    filter(KO != "")
  
  counts_df <- as.data.frame(counts) %>%
    rownames_to_column("orf_id") %>%
    inner_join(map_df, by = "orf_id")
  
  ko_counts <- counts_df %>%
    select(-orf_id) %>%
    group_by(KO) %>%
    summarise(across(where(is.numeric), ~ sum(.x, na.rm = TRUE)), .groups = "drop") %>%
    column_to_rownames("KO")
  
  ko_counts <- round(as.matrix(ko_counts))
  storage.mode(ko_counts) <- "integer"
  ko_counts
}

filter_ko_matrix <- function(ko_mat, min_total_counts = 5, min_samples_present = 2) {
  total_counts <- rowSums(ko_mat)
  n_present <- rowSums(ko_mat > 0)
  keep <- total_counts >= min_total_counts & n_present >= min_samples_present
  ko_mat[keep, , drop = FALSE]
}

make_coldata <- function(ps, vars = c("Treatment", "Timepoint")) {
  cd <- as(sample_data(ps), "data.frame")[, vars, drop = FALSE]
  cd$Treatment <- factor(cd$Treatment)
  cd$Timepoint <- factor(gsub("\\s+", "", cd$Timepoint))  # "6 hour" -> "6hour"
  cd$Interaction <- factor(paste(cd$Treatment, cd$Timepoint, sep = "_"))
  cd
}

#-----------------------------#
# Main
#-----------------------------#

stopifnot(exists("transcripts_phy_filt2"))
ps_in <- transcripts_phy_filt2

# 1) Filter samples + prune low-abundance ORFs
ps_filt <- filter_phyloseq_samples(
  ps_in,
  exclude_treatment = "procedural_control",
  min_library_size = 1000,
  min_taxa_sum = 1
)

# 2) Collapse ORFs -> KO
ko_counts <- collapse_orfs_to_kos(ps_filt, ko_col = "KEGG_ko", sep = ",")

# 3) Filter rare KOs
ko_counts_clean <- filter_ko_matrix(ko_counts, min_total_counts = 5, min_samples_present = 2)
nrow(ko_counts_clean)

# 4) Build colData aligned to KO matrix columns
coldata <- make_coldata(ps_filt)
coldata <- coldata[colnames(ko_counts_clean), , drop = FALSE]
stopifnot(identical(rownames(coldata), colnames(ko_counts_clean)))

# 5) DESeq2 (zero-tolerant size factors) 

# Build DESeq2 dataset
dds_clean <- DESeqDataSetFromMatrix(
  countData = ko_counts_clean,
  colData   = coldata,
  design    = ~ Interaction
)

# Set reference levels
dds_clean$Treatment <- relevel(dds_clean$Treatment, ref = "control")
dds_clean$Timepoint <- relevel(dds_clean$Timepoint, ref = "9hour")

# Drop KOs that are all-zero across all samples (safety)
dds_clean <- dds_clean[rowSums(counts(dds_clean)) > 0, ]

# Run DESeq2 using a size-factor method that handles zeros
dds_clean <- DESeq(dds_clean, sfType = "poscounts")

# Extract normalized counts
norm_counts <- counts(dds_clean, normalized = TRUE)

# 6) Global results summary 
res <- results(dds_clean)
res <- res[order(res$padj), ]

sig05 <- sum(res$padj < 0.05, na.rm = TRUE)
sig10 <- sum(res$padj < 0.10, na.rm = TRUE)

sig05
sig10

plotMA(res, ylim = c(-10, 10))

# Check the available coefficient names (important for contrasts)
resultsNames(dds_clean)

# 7) Treatment effect within each timepoint 
res_6h  <- results(dds_clean, contrast = c("Interaction", "exposure_6hour",  "control_6hour"))
res_9h  <- results(dds_clean, contrast = c("Interaction", "exposure_9hour",  "control_9hour"))
res_12h <- results(dds_clean, contrast = c("Interaction", "exposure_12hour", "control_12hour"))




