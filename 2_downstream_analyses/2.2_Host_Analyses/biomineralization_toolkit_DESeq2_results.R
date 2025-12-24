# ============================================================
# Host DESeq2 follow-up: biomineralization "toolkit" genes
#
# Goal:
#   Identify host DE transcripts (UP/DOWN) whose GO terms overlap
#   with a curated set of molluscan biomineralization toolkit GO IDs.
#
# Inputs assumed in your workspace:
#   - all_go_df: long GO annotation table with columns:
#       transcript_id, GOID  (one row per transcript x GO term)
#   - DESeq2 results (already LFC-shrunk + annotated with diffexpressed):
#       control_6hour_vs_exposure_6hour_resLFC
#       control_9hour_vs_exposure_9hour_resLFC
#       control_12hour_vs_exposure_12hour_resLFC
#     Each should contain transcript_id and a column "diffexpressed" in {UP, DOWN, ...}
#
# External file:
#   - Molluscan_Biomineralization_Toolkit_with_GO_Terms_and_References.csv
#     (contains biomineralization categories + Representative GO terms)
# ============================================================

library(dplyr)
library(stringr)
library(tidyr)
library(readr)

# -----------------------------
# 1) Define biomineralization GO terms of interest
# -----------------------------
toolkit_go <- unique(c(
  "GO:0006030", # chitin synthase activity
  "GO:0004568", # chitinase activity
  "GO:0005576", # extracellular region (often used in shell matrix / shematrin context)
  "GO:0051920","GO:0051213","GO:0004867","GO:0005515","GO:0004089","GO:0005509",
  "GO:0030154","GO:0003674","GO:0005615","GO:0008270","GO:0006955","GO:0007155",
  "GO:0030414","GO:0005388","GO:0005389","GO:0005245","GO:0015701","GO:0005391",
  "GO:0016021","GO:0015385"
))

# -----------------------------
# 2) Combine DE transcripts across timepoints
#    (keep only UP/DOWN calls)
# -----------------------------
deg_df <- bind_rows(
  as.data.frame(control_6hour_vs_exposure_6hour_resLFC)  %>%
    filter(diffexpressed %in% c("UP", "DOWN")) %>%
    mutate(Timepoint = "12:50"),
  as.data.frame(control_9hour_vs_exposure_9hour_resLFC)  %>%
    filter(diffexpressed %in% c("UP", "DOWN")) %>%
    mutate(Timepoint = "16:02"),
  as.data.frame(control_12hour_vs_exposure_12hour_resLFC) %>%
    filter(diffexpressed %in% c("UP", "DOWN")) %>%
    mutate(Timepoint = "19:14")
)

# -----------------------------
# 3) Attach GO term lists to each transcript
# -----------------------------
# Collapse GO annotations into a single comma-separated field per transcript
go_collapsed <- all_go_df %>%
  group_by(transcript_id) %>%
  summarise(
    GO_terms = paste(sort(unique(GOID)), collapse = ", "),
    .groups = "drop"
  )

deg_annot <- deg_df %>%
  left_join(go_collapsed, by = "transcript_id")

# -----------------------------
# 4) Keep only DE transcripts that overlap toolkit GO IDs
# -----------------------------
# Expand GO terms to long form and filter to toolkit list
deg_toolkit_long <- deg_annot %>%
  filter(!is.na(GO_terms) & GO_terms != "") %>%
  separate_rows(GO_terms, sep = ",\\s*") %>%
  filter(GO_terms %in% toolkit_go)

# Quick sanity checks
# n_distinct(deg_df$transcript_id)
# n_distinct(deg_toolkit_long$transcript_id)

# -----------------------------
# 5) Map toolkit GO IDs -> toolkit categories / gene names (from your CSV)
# -----------------------------
toolkit_ref <- read_csv(
  "Molluscan_Biomineralization_Toolkit_with_GO_Terms_and_References.csv",
  show_col_types = FALSE
)

# Extract GO IDs from the "Representative.GO.Term.s." column in the toolkit file
toolkit_ref_go <- toolkit_ref %>%
  mutate(GOID = str_extract_all(Representative.GO.Term.s., "GO:\\d+")) %>%
  unnest(GOID) %>%
  distinct(GOID, .keep_all = TRUE)

# Join DE transcripts -> GO term -> toolkit category/labels
deg_toolkit_annot <- deg_toolkit_long %>%
  left_join(toolkit_ref_go, by = c("GO_terms" = "GOID"))

# Optional: keep only records where the toolkit CSV provided a category
deg_toolkit_annot <- deg_toolkit_annot %>%
  filter(!is.na(Category))

# -----------------------------
# 6) Identify transcripts that map to a single biomineralization category
# -----------------------------
single_role_transcripts <- deg_toolkit_annot %>%
  group_by(transcript_id) %>%
  summarise(
    n_roles = n_distinct(Category),
    Categories = paste(sort(unique(Category)), collapse = ", "),
    .groups = "drop"
  ) %>%
  filter(n_roles == 1)

single_role_df <- deg_toolkit_annot %>%
  semi_join(single_role_transcripts, by = "transcript_id")
