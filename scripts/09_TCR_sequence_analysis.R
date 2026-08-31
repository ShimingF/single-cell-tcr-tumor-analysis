# ==============================================================================
# Script: 09_TCR_sequence_analysis.R
# Purpose: Comprehensive V(D)J gene usage for both TRAV and TRBV paired chains
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(stringr)
  library(tidyr)
  library(patchwork)
})

cat("==> [09_TCR_sequence_analysis.R] Starting V(D)J Gene & CDR3 Sequence Analysis (TRAV + TRBV)...\n")

# 1. Load Master Integrated Dataset
master_path <- "data/processed/final_paired_tcell_metadata.rds"

if (!file.exists(master_path)) {
  stop("Master dataset not found! Please run scripts/04 and 05 first.")
}

master_df <- readRDS(master_path)
cat(sprintf("   - Loaded Master Dataset: %d strict 1a1b paired cells.\n", nrow(master_df)))

# 2. Extract Unique Clonotypes per Patient
clonotypes <- master_df %>%
  distinct(patient, tcr_sequence_id, .keep_all = TRUE) %>%
  mutate(
    cdr3a_len = ifelse(!is.na(CDR3a), nchar(CDR3a), NA_integer_),
    cdr3b_len = ifelse(!is.na(CDR3b), nchar(CDR3b), NA_integer_)
  )

# 3. Analyze Both TRAV & TRBV Gene Usage Frequency
v_gene_usage_trb <- clonotypes %>%
  filter(!is.na(TRB_v), TRB_v != "", TRB_v != "None") %>%
  count(TRB_v, sort = TRUE) %>%
  rename(v_gene = TRB_v) %>%
  mutate(percentage = round(n / sum(n) * 100, 2), chain = "TRBV")

v_gene_usage_tra <- clonotypes %>%
  filter(!is.na(TRA_v), TRA_v != "", TRA_v != "None") %>%
  count(TRA_v, sort = TRUE) %>%
  rename(v_gene = TRA_v) %>%
  mutate(percentage = round(n / sum(n) * 100, 2), chain = "TRAV")

# 4. Export Comprehensive Tables
dir.create("results/tables", showWarnings = FALSE, recursive = TRUE)
write_csv(v_gene_usage_trb, "results/tables/09_trb_vgene_usage.csv")
write_csv(v_gene_usage_tra, "results/tables/09_tra_vgene_usage.csv")

cat("\n--- Top 5 TRAV Genes ---\n")
print(head(v_gene_usage_tra, 5))

cat("\n--- Top 5 TRBV Genes ---\n")
print(head(v_gene_usage_trb, 5))

# 5. Generate Figures for both TRAV and TRBV
dir.create("results/figures", showWarnings = FALSE, recursive = TRUE)

# Plot A: Top 15 TRAV Barplot
p_trav <- ggplot(head(v_gene_usage_tra, 15), aes(x = reorder(v_gene, n), y = percentage, fill = percentage)) +
  geom_col(show.legend = FALSE) +
  coord_flip() +
  scale_fill_viridis_c(option = "viridis", direction = -1) +
  theme_bw() +
  labs(
    title = "Top 15 TRAV Gene Usage Preference",
    x = "TRA V Gene",
    y = "Frequency (%)"
  ) +
  theme(plot.title = element_text(face = "bold", size = 11))

# Plot B: Top 15 TRBV Barplot
p_trbv <- ggplot(head(v_gene_usage_trb, 15), aes(x = reorder(v_gene, n), y = percentage, fill = percentage)) +
  geom_col(show.legend = FALSE) +
  coord_flip() +
  scale_fill_viridis_c(option = "mako", direction = -1) +
  theme_bw() +
  labs(
    title = "Top 15 TRBV Gene Usage Preference",
    x = "TRB V Gene",
    y = "Frequency (%)"
  ) +
  theme(plot.title = element_text(face = "bold", size = 11))

# Save Single Chain Plots
ggsave("results/figures/09_top15_trav_vgene_usage.png", plot = p_trav, width = 7, height = 5, dpi = 300)
ggsave("results/figures/09_top15_trb_vgene_usage.png", plot = p_trbv, width = 7, height = 5, dpi = 300)

# Plot C: Side-by-Side Dual Chain Comparison (TRAV + TRBV)
p_combined <- p_trav + p_trbv + 
  plot_annotation(
    title = "GSE139555 T-Cell Receptor V-Gene Pairing Landscape (TRAV vs TRBV)",
    theme = theme(plot.title = element_text(face = "bold", size = 13))
  )

ggsave("results/figures/09_top15_trav_and_trbv_usage.png", plot = p_combined, width = 13, height = 5.5, dpi = 300)

# Plot D: CDR3 Length Distribution (?? & ?? chains)
p_len <- ggplot(clonotypes, aes(x = cdr3b_len)) +
  geom_histogram(binwidth = 1, fill = "#2b5c8f", color = "white", alpha = 0.85) +
  theme_bw() +
  labs(
    title = "TRB CDR3 Amino Acid Length Distribution",
    x = "CDR3b Length (aa)",
    y = "Number of Unique Clonotypes"
  )

ggsave("results/figures/09_cdr3b_length_distribution.png", plot = p_len, width = 7, height = 5, dpi = 300)

cat("\n==> [09_TCR_sequence_analysis.R] Completed successfully!")
cat("\n   - Exported Figure 1: results/figures/09_top15_trav_vgene_usage.png")
cat("\n   - Exported Figure 2: results/figures/09_top15_trb_vgene_usage.png")
cat("\n   - Exported Figure 3: results/figures/09_top15_trav_and_trbv_usage.png\n\n")