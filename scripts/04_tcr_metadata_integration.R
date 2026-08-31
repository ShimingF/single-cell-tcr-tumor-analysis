# ==============================================================================
# Script: 04_tcr_metadata_integration.R
# Purpose: Merge Strict Paired TCRs with Clean Metadata to build Master Dataset
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
})

cat("==> [04_tcr_metadata_integration.R] Starting Metadata & TCR Synergy Merge...\n")

# 1. Load Data
metadata     <- readRDS("data/processed/clean_metadata.rds")
strict_cells <- read_csv("data/processed/tcr_strict_paired_cells.csv", show_col_types = FALSE)

# 2. Inner Join via cell_key
master_df <- inner_join(
  metadata,
  strict_cells %>% select(-sample, -barcode),
  by = "cell_key"
) %>%
  mutate(
    # True Sequence Identifier (CDR3a + CDR3b)
    tcr_sequence_id = paste(CDR3a, CDR3b, sep = "_"),
    # Patient-level clonotype string
    patient_clonotype = ifelse(!is.na(patient_clonotype_legacy), patient_clonotype_legacy, paste(patient, "unassigned", sep = "."))
  )

cat(sprintf("   - Total Clean Metadata Cells: %d\n", nrow(metadata)))
cat(sprintf("   - Total Strict 1a1b TCR Cells: %d\n", nrow(strict_cells)))
cat(sprintf("   - Final High-Confidence Matched Cells: %d\n", nrow(master_df)))

# 3. Export Master Dataset
saveRDS(master_df, "data/processed/final_paired_tcell_metadata.rds")
dir.create("results/tables", showWarnings = FALSE, recursive = TRUE)
write_csv(master_df, "results/tables/final_paired_tcell_metadata.csv")

cat("==> [04_tcr_metadata_integration.R] Completed. Master Dataset saved to data/processed/final_paired_tcell_metadata.rds\n\n")