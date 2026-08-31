# ==============================================================================
# Script: 02_tcr_contig_qc.R
# Purpose: Import 32 GEO TCR files and parse exact sample IDs (e.g., -lt1 -> LT1)
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(purrr)
  library(stringr)
  library(fs)
})

cat("==> [02_tcr_contig_qc.R] Starting Contig-level Quality Control...\n")

# 1. Target Directory & File Search
tcr_dir <- Sys.getenv("TCR_DATA_DIR", unset = "data/raw/TCR")

if (!dir.exists(tcr_dir)) {
  stop(
    "TCR data directory not found: ", tcr_dir,
    "\nPlace the GEO TCR files under data/raw/TCR ",
    "or set environment variable TCR_DATA_DIR."
  )
}

tcr_files <- list.files(tcr_dir, pattern = "\\.filtered_contig_annotations\\.csv$", full.names = TRUE)

# Fallback pattern if named slightly differently
if (length(tcr_files) == 0) {
  tcr_files <- list.files(tcr_dir, pattern = "\\.csv$", full.names = TRUE)
}

cat(sprintf("   - Target Path: %s\n", tcr_dir))
cat(sprintf("   - Found %d TCR contig files.\n", length(tcr_files)))

if (length(tcr_files) == 0) {
  stop("Error: No CSV files found in target directory!")
}

# 2. Contig Processing Function with Precise Sample ID Parsing
process_geo_tcr_file <- function(file_path) {
  fname <- basename(file_path)
  
  # Extract sample_id after "-" and before "." (e.g. "lt1" from "...-lt1.filtered...")
  sample_raw <- str_extract(fname, "(?<=-)[a-zA-Z0-9]+(?=\\.)")
  
  # Standardize sample_id to uppercase (e.g., "lt1" -> "LT1") to match Metadata
  sample_id <- toupper(sample_raw)
  
  if (is.na(sample_id)) {
    warning(sprintf("Could not parse sample ID from %s, skipping.", fname))
    return(NULL)
  }
  
  df <- read_csv(file_path, show_col_types = FALSE, col_types = cols(.default = "c"))
  
  # Validate required columns
  if (!all(c("barcode", "chain", "cdr3") %in% colnames(df))) {
    warning(sprintf("Skipping %s: Missing core columns.", fname))
    return(NULL)
  }
  
  df_qc <- df %>%
    mutate(
      is_cell = if("is_cell" %in% colnames(.)) as.logical(is_cell) else TRUE,
      high_confidence = if("high_confidence" %in% colnames(.)) as.logical(high_confidence) else TRUE,
      full_length = if("full_length" %in% colnames(.)) as.logical(full_length) else TRUE,
      productive = if("productive" %in% colnames(.)) as.logical(productive) else TRUE
    ) %>%
    # Strict Contig Filters
    filter(
      is_cell == TRUE,
      high_confidence == TRUE,
      full_length == TRUE,
      productive == TRUE,
      chain %in% c("TRA", "TRB"),
      !is.na(cdr3), cdr3 != "None"
    ) %>%
    mutate(
      sample = sample_id,
      # Retain original barcode suffix (-1)
      raw_barcode = str_extract(barcode, "[ACGT]+-[0-9]+$"),
      raw_barcode = ifelse(is.na(raw_barcode), barcode, raw_barcode),
      # Construct cell_key matching Script 01 (e.g., "LT1_AAACCTGAGGATATAC-1")
      cell_key = paste(sample_id, raw_barcode, sep = "_")
    )
  
  return(df_qc)
}

# 3. Batch Process 32 Files
all_contigs_qc <- map_dfr(tcr_files, process_geo_tcr_file)

cat(sprintf("   - Passed Productive Contigs: %d\n", nrow(all_contigs_qc)))
cat(sprintf("   - Unique Cells represented in TCR data: %d\n", n_distinct(all_contigs_qc$cell_key)))

# 4. Export Result
dir.create("data/processed", showWarnings = FALSE, recursive = TRUE)
write_csv(all_contigs_qc, "data/processed/tcr_contigs_qc.csv")
cat("==> [02_tcr_contig_qc.R] Completed successfully! Saved to: data/processed/tcr_contigs_qc.csv\n\n")