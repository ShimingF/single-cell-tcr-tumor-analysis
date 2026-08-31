# ==============================================================================
# Script: 10_TCRex.R
# Purpose: Export strict TCRex-standard TSV format (CDR3_beta, TRBV_gene, TRBJ_gene)
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
})

cat("==> [10_TCRex.R] Preparing Data Export with Exact TCRex Headers...\n")

# 1. Load Master Integrated Dataset
master_path <- "data/processed/final_paired_tcell_metadata.rds"

if (!file.exists(master_path)) {
  stop("Master dataset not found! Please run scripts/04 and 05 first.")
}

master_df <- readRDS(master_path)

# 2. Extract Expanded CDR3b Sequences & Rename to Exact TCRex Standard Headers
tcrex_input <- master_df %>%
  filter(expansion_status %in% c("Expanded", "Highly expanded", "Hyper-expanded")) %>%
  distinct(patient, tcr_sequence_id, .keep_all = TRUE) %>%
  mutate(
    # Clean sequence and gene names
    CDR3_beta = str_trim(as.character(CDR3b)),
    TRBV_gene = str_trim(as.character(TRB_v)),
    TRBJ_gene = str_trim(as.character(TRB_j))
  ) %>%
  # Filter strictly for valid amino acid sequences starting with C
  filter(
    !is.na(CDR3_beta), CDR3_beta != "", CDR3_beta != "None",
    str_detect(CDR3_beta, "^C[A-Z]+"),
    !is.na(TRBV_gene), TRBV_gene != "", TRBV_gene != "None"
  ) %>%
  select(CDR3_beta, TRBV_gene, TRBJ_gene) %>%
  distinct()

cat(sprintf("   - Extracted %d high-quality expanded clonotypes.\n", nrow(tcrex_input)))

# 3. Export Plain TSV without Quotes
dir.create("results/tables", showWarnings = FALSE, recursive = TRUE)
out_tsv_path <- "results/tables/10_tcrex_input_expanded_clones.tsv"

write.table(
  tcrex_input,
  file = out_tsv_path,
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE,
  na = ""
)

cat(sprintf("   - Exported file to: %s\n", out_tsv_path))
cat("   - Preview of exported TCRex TSV:\n")
print(head(tcrex_input, 6))

cat("\n==> [10_TCRex.R] Completed successfully!\n\n")