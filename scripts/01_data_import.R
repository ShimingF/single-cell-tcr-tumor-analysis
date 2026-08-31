# ==============================================================================
# Script: 01_data_import.R
# Purpose: Clean GSE139555 metadata with exact raw TSV row name handling
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(tibble)
})

cat("==> [01_data_import.R] Starting Metadata Import & Standardization...\n")

# 1. Paths
meta_path   <- "data/raw/GSE139555_tcell_metadata.txt"
output_path <- "data/processed/clean_metadata.rds"

if (!file.exists(meta_path)) {
  meta_path <- list.files("data/raw", pattern = "metadata", full.names = TRUE)[1]
}

# 2. Read Metadata
df_raw <- read_tsv(meta_path, show_col_types = FALSE)

# 3. Handle TSV Unnamed First Column (Row Names)
# If the first column name is automatically assigned as ...1 or missing, set as cell_id
if ("...1" %in% colnames(df_raw)) {
  df_raw <- df_raw %>% rename(cell_id = ...1)
} else if (colnames(df_raw)[1] != "cell_id") {
  colnames(df_raw)[1] <- "cell_id"
}

# 4. Standardize Fields
df_clean <- df_raw %>%
  filter(!is.na(cell_id), !is.na(sample)) %>%
  mutate(
    # cell_id itself is already formatted as "LT1_AAACCTGAGGATATAC-1", making it the perfect cell_key
    cell_key = cell_id,
    # Extract raw 10x barcode with -1 suffix (e.g. "AAACCTGAGGATATAC-1")
    barcode = str_extract(cell_id, "[ACGT]+-[0-9]+$"),
    # Parse cancer type (e.g., "Lung1" -> cancer="Lung", patient_num="1")
    cancer_type = str_extract(patient, "^[a-zA-Z]+"),
    patient_number = str_extract(patient, "[0-9]+$"),
    # Store legacy metadata clonotype
    patient_clonotype_legacy = ifelse("clonotype" %in% colnames(.), clonotype, NA_character_)
  ) %>%
  distinct(cell_key, .keep_all = TRUE)

cat(sprintf("   - Total Metadata Cells Processed: %d\n", nrow(df_clean)))
cat(sprintf("   - Patients: %d | Samples: %d\n", n_distinct(df_clean$patient), n_distinct(df_clean$sample)))

# 5. Export Clean Metadata
dir.create("data/processed", showWarnings = FALSE, recursive = TRUE)
saveRDS(df_clean, output_path)
cat("==> [01_data_import.R] Completed. Exported to:", output_path, "\n\n")