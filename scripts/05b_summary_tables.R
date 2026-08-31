# ==============================================================================
# Script: 05b_summary_tables.R
# Purpose: Generate comprehensive summary tables across Patients, Cancer Types, and Sources
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(ggplot2)
})


cat("==> [05b_summary_tables.R] Generating Multi-dimensional Summary Tables...\n")

# 1. Load Master Integrated Dataset
master_path <- "data/processed/final_paired_tcell_metadata.rds"

if (!file.exists(master_path)) {
  stop("Master dataset not found! Please run scripts/04 and 05 first.")
}

df <- readRDS(master_path)
dir.create("results/tables", showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------------------------
# Table 1: Patient-level Comprehensive Summary
# ------------------------------------------------------------------------------
table_patient <- df %>%
  group_by(cancer_type, patient) %>%
  summarise(
    total_strict_cells = n(),
    unique_clonotypes  = n_distinct(tcr_sequence_id),
    n_tumor_cells      = sum(source == "Tumor", na.rm = TRUE),
    n_nat_cells        = sum(source == "NAT", na.rm = TRUE),
    n_blood_cells      = sum(source == "Blood", na.rm = TRUE),
    hyper_expanded_clones = n_distinct(tcr_sequence_id[expansion_status == "Hyper-expanded"]),
    expanded_clones       = n_distinct(tcr_sequence_id[expansion_status %in% c("Expanded", "Highly expanded", "Hyper-expanded")]),
    singleton_clones      = n_distinct(tcr_sequence_id[expansion_status == "Singleton"]),
    .groups = "drop"
  ) %>%
  arrange(cancer_type, patient)

write_csv(table_patient, "results/tables/summary_by_patient.csv")
cat("   - Exported: results/tables/summary_by_patient.csv\n")

# ------------------------------------------------------------------------------
# Table 2: Cancer Type-level Summary
# ------------------------------------------------------------------------------
table_cancer <- df %>%
  group_by(cancer_type) %>%
  summarise(
    n_patients         = n_distinct(patient),
    total_strict_cells = n(),
    unique_clonotypes  = n_distinct(tcr_sequence_id),
    n_tumor_cells      = sum(source == "Tumor", na.rm = TRUE),
    n_nat_cells        = sum(source == "NAT", na.rm = TRUE),
    n_blood_cells      = sum(source == "Blood", na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    cell_pct = round(total_strict_cells / sum(total_strict_cells) * 100, 2)
  )

write_csv(table_cancer, "results/tables/summary_by_cancer_type.csv")
cat("   - Exported: results/tables/summary_by_cancer_type.csv\n")

# ------------------------------------------------------------------------------
# Table 3: Tissue Source (Sampling Location) Summary
# ------------------------------------------------------------------------------
table_source <- df %>%
  group_by(source) %>%
  summarise(
    total_strict_cells = n(),
    unique_clonotypes  = n_distinct(tcr_sequence_id),
    n_patients_covered = n_distinct(patient),
    .groups = "drop"
  ) %>%
  mutate(
    cell_pct = round(total_strict_cells / sum(total_strict_cells) * 100, 2)
  )

write_csv(table_source, "results/tables/summary_by_tissue_source.csv")
cat("   - Exported: results/tables/summary_by_tissue_source.csv\n")

# ------------------------------------------------------------------------------
# Table 4: Multi-dimensional Matrix (Cancer Type x Tissue Source)
# ------------------------------------------------------------------------------
table_cancer_x_source <- df %>%
  group_by(cancer_type, source) %>%
  summarise(cell_count = n(), .groups = "drop") %>%
  pivot_wider(names_from = source, values_from = cell_count, values_fill = 0)

write_csv(table_cancer_x_source, "results/tables/summary_cancer_x_tissue_matrix.csv")
cat("   - Exported: results/tables/summary_cancer_x_tissue_matrix.csv\n")

# ------------------------------------------------------------------------------
# Table 5: Expansion Tier Distribution per Tissue Source
# ------------------------------------------------------------------------------
table_expansion_x_source <- df %>%
  group_by(source, expansion_status) %>%
  summarise(
    cell_count = n(),
    n_clonotypes = n_distinct(tcr_sequence_id),
    .groups = "drop"
  ) %>%
  arrange(source, expansion_status)

write_csv(table_expansion_x_source, "results/tables/summary_expansion_by_tissue.csv")
cat("   - Exported: results/tables/summary_expansion_by_tissue.csv\n")

cat("==> [05b_summary_tables.R] All 5 summary tables exported successfully!\n\n")
