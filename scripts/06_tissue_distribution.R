# ==============================================================================
# Script: 06_tissue_distribution.R
# Purpose: Tumor / NAT / Blood cross-tissue clonal sharing & distribution analysis
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(tidyr)
})

cat("==> [06_tissue_distribution.R] Starting Tissue Distribution Analysis...\n")

# 1. Load Master Integrated Dataset (Correct File Path)
master_path <- "data/processed/final_paired_tcell_metadata.rds"

if (!file.exists(master_path)) {
  stop("Master dataset not found! Please ensure scripts/04 and 05 ran successfully.")
}

master_df <- readRDS(master_path)

cat(sprintf("   - Loaded Master Dataset: %d strict 1a1b paired cells.\n", nrow(master_df)))

# 2. Tissue Sharing Classification per Clonotype
tissue_sharing_summary <- master_df %>%
  group_by(patient, tcr_sequence_id, patient_clonotype_clean, expansion_status) %>%
  summarise(
    n_tumor = sum(source == "Tumor", na.rm = TRUE),
    n_nat   = sum(source == "NAT", na.rm = TRUE),
    n_blood = sum(source == "Blood", na.rm = TRUE),
    total_size = n(),
    .groups = "drop"
  ) %>%
  mutate(
    tissue_sharing_type = case_when(
      n_tumor > 0 & n_nat > 0 & n_blood > 0 ~ "Tumor + NAT + Blood",
      n_tumor > 0 & n_nat > 0 & n_blood == 0 ~ "Tumor + NAT",
      n_tumor > 0 & n_blood > 0 & n_nat == 0 ~ "Tumor + Blood",
      n_nat > 0 & n_blood > 0 & n_tumor == 0 ~ "NAT + Blood",
      n_tumor > 0 & n_nat == 0 & n_blood == 0 ~ "Tumor Only",
      n_nat > 0 & n_tumor == 0 & n_blood == 0 ~ "NAT Only",
      n_blood > 0 & n_tumor == 0 & n_nat == 0 ~ "Blood Only",
      TRUE ~ "Other"
    )
  )

# 3. Export Summary Tables
dir.create("results/tables", showWarnings = FALSE, recursive = TRUE)

# Overall Sharing Category Breakdown
sharing_stats <- tissue_sharing_summary %>%
  count(tissue_sharing_type) %>%
  mutate(pct = round(n / sum(n) * 100, 2))

write_csv(tissue_sharing_summary, "results/tables/06_tissue_sharing_per_clonotype.csv")
write_csv(sharing_stats, "results/tables/06_tissue_sharing_category_stats.csv")

cat("\n--- Tissue Sharing Breakdown ---\n")
print(sharing_stats)

# 4. Plot Tissue Distribution Visualization
dir.create("results/figures", showWarnings = FALSE, recursive = TRUE)

p_sharing <- ggplot(tissue_sharing_summary, aes(x = tissue_sharing_type, fill = tissue_sharing_type)) +
  geom_bar() +
  theme_minimal() +
  labs(
    title = "Cross-Tissue Clonal Distribution (Strict 1a1b Paired)",
    x = "Tissue Compartment Sharing",
    y = "Number of Unique Clonotypes"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

ggsave("results/figures/06_tissue_distribution_sharing.png", plot = p_sharing, width = 8, height = 5)

cat("\n==> [06_tissue_distribution.R] Completed successfully!\n")
cat("   - Exported Table:  results/tables/06_tissue_sharing_per_clonotype.csv\n")
cat("   - Exported Figure: results/figures/06_tissue_distribution_sharing.png\n\n")