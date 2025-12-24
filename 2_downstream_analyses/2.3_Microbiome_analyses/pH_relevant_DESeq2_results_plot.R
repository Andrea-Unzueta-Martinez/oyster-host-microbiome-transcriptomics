library(dplyr)
library(stringr)
library(tidyr)
library(tibble)
library(ggplot2)

# -----------------------------
# 1) Define KO set of interest
# -----------------------------
# Curated KOs relevant to carbonate chemistry / pH regulation / redox cycling
kegg_kos_ph <- c(
  # Carbonic anhydrase
  "ko:K01672","ko:K01673","ko:K01674","ko:K01743","ko:K18245","ko:K18246",
  # Urease
  "ko:K01427","ko:K01428","ko:K01429","ko:K01430","ko:K14048",
  # Denitrification (subset you listed)
  "ko:K00370","ko:K00371","ko:K00374","ko:K02567","ko:K02568","ko:K00368","ko:K15864","ko:K04561","ko:K02305","ko:K00376",
  # Assimilatory sulfate reduction
  "ko:K13811","ko:K00958","ko:K00860","ko:K00955","ko:K00957","ko:K00956","ko:K00390","ko:K05907","ko:K00380","ko:K00381","ko:K00392",
  # Dissimilatory sulfate reduction
  "ko:K00394","ko:K00395","ko:K11180","ko:K11181","ko:K27196","ko:K27187","ko:K27188","ko:K27189","ko:K27190","ko:K27191",
  # Assimilatory nitrate reduction
  "ko:K00367","ko:K10534","ko:K00372","ko:K00360","ko:K00366","ko:K17877","ko:K26139","ko:K26138","ko:K00361",
  # Dissimilatory nitrate reduction
  "ko:K00362","ko:K00363","ko:K03385","ko:K15876",
  # Methanogenesis-related (combined list from your code)
  "ko:K00925","ko:K00625","ko:K01895","ko:K00193","ko:K00197","ko:K00194",
  "ko:K00577","ko:K00578","ko:K00579","ko:K00580","ko:K00581","ko:K00582","ko:K00583","ko:K00584",
  "ko:K00399","ko:K00401","ko:K00402","ko:K22480","ko:K22481","ko:K22482",
  "ko:K03388","ko:K03389","ko:K03390","ko:K08264","ko:K08265",
  "ko:K14127","ko:K14126","ko:K14128","ko:K22516","ko:K00125",
  "ko:K00200","ko:K00201","ko:K00202","ko:K00203","ko:K11261","ko:K00205","ko:K11260","ko:K00204","ko:K00672","ko:K01499",
  "ko:K00319","ko:K13942","ko:K00320",
  # Sulfide oxidation (SOX)
  "ko:K17222","ko:K17223","ko:K17224","ko:K17225","ko:K22622","ko:K17226","ko:K17227",
  # Ammonium oxidation
  "ko:K20932","ko:K20933","ko:K20934","ko:K20935",
  "ko:K10944","ko:K10945","ko:K10946","ko:K10535",
  # Extra (in your methylotrophic block)
  "ko:K14080","ko:K04480","ko:K14081"
) %>% unique()

# -------------------------------------------------------
# 2) Helper: subset DESeq2 results to significant KOs only
# -------------------------------------------------------
# - Adds KEGG_KO from rownames
# - Computes diffexpressed label from padj + direction of log2FC
# - Filters to your KO list of interest
make_tornado_df <- function(res_obj, time_label, ko_keep, padj_cutoff = 0.05) {
  as.data.frame(res_obj) %>%
    rownames_to_column("KEGG_KO") %>%
    mutate(
      Timepoint = time_label,
      diffexpressed = if_else(!is.na(padj) & padj < padj_cutoff,
                              if_else(log2FoldChange > 0, "exposure", "control"),
                              "NO")
    ) %>%
    filter(diffexpressed != "NO") %>%
    filter(KEGG_KO %in% ko_keep)
}

# ------------------------------------------
# 3) Build one table across all timepoints
# ------------------------------------------
tornado_dat <- bind_rows(
  make_tornado_df(control_6hour_vs_exposure_6hour,   "12:50", kegg_kos_ph),
  make_tornado_df(control_9hour_vs_exposure_9hour,   "16:02", kegg_kos_ph),
  make_tornado_df(control_12hour_vs_exposure_12hour, "19:14", kegg_kos_ph)
) %>%
  mutate(Timepoint = factor(Timepoint, levels = c("12:50","16:02","19:14")))

# ------------------------------------------
# 4) Join KO definitions + module metadata
# ------------------------------------------
# Expectation:
#   KO_definitions:      KEGG_KO, KO_definition
#   KO_to_Module_key:    KEGG_KO, KEGG_Pathway, KEGG_Module
#   Module_definitions:  KEGG_Module, <module definition columns>

colnames(KO_definitions)   <- c("KEGG_KO", "KO_definition")
colnames(KO_to_Module_key) <- c("KEGG_KO", "KEGG_Pathway", "KEGG_Module")

tornado_dat_annot <- tornado_dat %>%
  left_join(KO_definitions, by = "KEGG_KO") %>%
  mutate(
    # KO symbol = first chunk before ";" in the KO definition field
    KO_Symbol = str_split(KO_definition, ";", simplify = TRUE)[, 1] %>% str_trim()
  ) %>%
  left_join(KO_to_Module_key, by = "KEGG_KO") %>%
  left_join(Module_definitions, by = "KEGG_Module") %>%
  # Keep one row per KO x timepoint (prevents duplicate joins exploding rows)
  distinct(KEGG_KO, Timepoint, .keep_all = TRUE)

# ----------------------------------------------------
# 5) Clean labels for module coloring (stable mapping)
# ----------------------------------------------------
# Map KEGG module IDs to the exact “clean” labels you want to show in the plot.
module_clean_map <- c(
  "M00529" = "M00529 Denitrification, nitrate => nitrogen",
  "M00530" = "M00530 Dissimilatory nitrate reduction, nitrate => ammonia",
  "M00531" = "M00531 Assimilatory nitrate reduction, nitrate => ammonia",
  "M00176" = "M00176 Assimilatory sulfate reduction, sulfate => H2S",
  "M00595" = "M00595 Thiosulfate oxidation by SOX complex, thiosulfate => sulfate"
)

tornado_dat_annot <- tornado_dat_annot %>%
  mutate(
    # Override to keep CA + urease colored as their own categories
    KEGG_Module_definitions_clean = case_when(
      str_detect(KO_Symbol, regex("^ure", ignore_case = TRUE)) ~ "urease subunit alpha",
      str_detect(KO_Symbol, regex("cah|cynT|can", ignore_case = TRUE)) ~ "carbonic anhydrase",
      TRUE ~ unname(module_clean_map[KEGG_Module])
    ),
    KEGG_Module_definitions_clean = factor(
      KEGG_Module_definitions_clean,
      levels = c(
        "M00529 Denitrification, nitrate => nitrogen",
        "M00530 Dissimilatory nitrate reduction, nitrate => ammonia",
        "M00531 Assimilatory nitrate reduction, nitrate => ammonia",
        "M00176 Assimilatory sulfate reduction, sulfate => H2S",
        "M00595 Thiosulfate oxidation by SOX complex, thiosulfate => sulfate",
        "urease subunit alpha",
        "carbonic anhydrase"
      )
    )
  )

# ------------------------------------------
# 6) Tornado plot (faceted by timepoint)
# ------------------------------------------
important_genes <- ggplot(
  tornado_dat_annot,
  aes(x = log2FoldChange, y = KO_Symbol, fill = KEGG_Module_definitions_clean)
) +
  geom_col() +
  facet_grid(. ~ Timepoint, scales = "free_y") +
  scale_y_discrete(limits = rev) +
  labs(
    x = "log2 fold change (exposure vs control)",
    y = "KEGG Orthologue",
    fill = "Biochemical reaction (KEGG module)"
  ) +
  theme_classic()

important_genes
