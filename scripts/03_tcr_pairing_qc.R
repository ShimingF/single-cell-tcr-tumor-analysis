# ==============================================================================
# Script: 03_tcr_pairing_qc.R
# Purpose: Cell-level chain aggregation & Strict 1a1b Pairing Rule Enforcement
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
})

cat("==> [03_tcr_pairing_qc.R] Starting Strict 1a1b TCR Pairing...\n")

# 1. Load Contigs
contigs <- read_csv("data/processed/tcr_contigs_qc.csv", show_col_types = FALSE)

# 2. Aggregate Productive Chains per Cell Key
cell_summary <- contigs %>%
  group_by(cell_key, sample, barcode, chain) %>%
  summarise(
    n_productive = n(),
    cdr3_seq = paste(unique(cdr3), collapse = ";"),
    v_gene   = paste(unique(v_gene), collapse = ";"),
    j_gene   = paste(unique(j_gene), collapse = ";"),
    raw_clonotype_id = first(raw_clonotype_id),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = chain,
    values_from = c(n_productive, cdr3_seq, v_gene, j_gene),
    values_fill = list(n_productive = 0, cdr3_seq = NA_character_, v_gene = NA_character_, j_gene = NA_character_)
  )

# 3. Apply Strict Pairing Logic
cell_paired <- cell_summary %>%
  mutate(
    n_TRA = n_productive_TRA,
    n_TRB = n_productive_TRB,
    pairing_status = case_when(
      n_TRA == 1 & n_TRB == 1 ~ "Strict_Paired_1a1b",
      n_TRA > 1  & n_TRB == 1 ~ "Dual_TRA",
      n_TRA == 1 & n_TRB > 1 ~ "Multi_TRB",
      n_TRA > 1  & n_TRB > 1 ~ "Multi_Chain_Complex",
      n_TRA == 1 & n_TRB == 0 ~ "Orphan_TRA",
      n_TRA == 0 & n_TRB == 1 ~ "Orphan_TRB",
      TRUE                    ~ "Other"
    )
  )

# 4. Print Breakdown
cat("\n--- Pairing Status Breakdown ---\n")
print(count(cell_paired, pairing_status) %>% mutate(pct = round(n/sum(n)*100, 2)))

# 5. Extract ONLY Strict 1a1b Cells
strict_cells <- cell_paired %>%
  filter(pairing_status == "Strict_Paired_1a1b") %>%
  select(
    cell_key, sample, barcode,
    CDR3a = cdr3_seq_TRA, TRA_v = v_gene_TRA, TRA_j = j_gene_TRA,
    CDR3b = cdr3_seq_TRB, TRB_v = v_gene_TRB, TRB_j = j_gene_TRB,
    raw_clonotype_id, pairing_status
  )

cat(sprintf("\n   - Retained Strict Paired 1a1b Cells: %d\n", nrow(strict_cells)))

# 6. Export
write_csv(cell_paired, "data/processed/tcr_pairing_all_cells_qc.csv")
write_csv(strict_cells, "data/processed/tcr_strict_paired_cells.csv")
cat("==> [03_tcr_pairing_qc.R] Completed. Exported to: data/processed/tcr_strict_paired_cells.csv\n\n")