# ==============================================================================
# Script: 08_patient_level.R
# Purpose: Patient-level TCR repertoire diversity analysis and automated figure exports
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(tidyr)
})

cat("==> [08_patient_level.R] Starting Patient-Level Diversity Analysis...\n")

# 1. Load Master Integrated Dataset
master_path <- "data/processed/final_paired_tcell_metadata.rds"

if (!file.exists(master_path)) {
  stop("Master dataset not found! Please run scripts/04 and 05 first.")
}

master_df <- readRDS(master_path)

# 2. Calculate Diversity Indices per Patient & Tissue Source
patient_diversity <- master_df %>%
  group_by(patient, cancer_type, source) %>%
  summarise(
    total_cells = n(),
    unique_clones = n_distinct(tcr_sequence_id),
    # Shannon Entropy Index
    shannon_index = -sum((table(tcr_sequence_id) / n()) * log(table(tcr_sequence_id) / n())),
    # Gini-Simpson Index
    simpson_index = 1 - sum((table(tcr_sequence_id) / n())^2),
    .groups = "drop"
  ) %>%
  mutate(
    shannon_index = ifelse(is.nan(shannon_index), 0, round(shannon_index, 3)),
    simpson_index = ifelse(is.nan(simpson_index), 0, round(simpson_index, 3))
  )

cat("\n--- Patient Diversity Index Summary ---\n")
print(head(patient_diversity, 6))

# 3. Export Data Tables
dir.create("results/tables", showWarnings = FALSE, recursive = TRUE)
write_csv(patient_diversity, "results/tables/08_patient_repertoire_diversity.csv")

# 4. Generate & Save Publication-Ready Figures
dir.create("results/figures", showWarnings = FALSE, recursive = TRUE)

# Figure A: Boxplot of Shannon Diversity grouped by Tissue Source & Cancer Type
p_box <- ggplot(patient_diversity, aes(x = source, y = shannon_index, fill = source)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6, width = 0.5) +
  geom_jitter(aes(color = cancer_type), width = 0.15, size = 2.5, alpha = 0.8) +
  theme_bw() +
  labs(
    title = "TCR Repertoire Shannon Diversity Across Tissue Compartments",
    subtitle = "Calculated on Strict 1a1b Paired Master Dataset",
    x = "Tissue Source Compartment",
    y = "Shannon Diversity Index (H')",
    fill = "Tissue Compartment",
    color = "Cancer Type"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    axis.text = element_text(size = 10),
    panel.grid.minor = element_blank()
  )

ggsave("results/figures/08_patient_shannon_diversity_by_tissue.png", plot = p_box, width = 8, height = 5, dpi = 300)

# Figure B: Tile/Heatmap of Shannon Diversity (Patient x Source)
p_heat <- ggplot(patient_diversity, aes(x = source, y = patient, fill = shannon_index)) +
  geom_tile(color = "white", linewidth = 0.5) +
  scale_fill_viridis_c(option = "magma", name = "Shannon Index") +
  theme_minimal() +
  labs(
    title = "Patient-Level TCR Diversity Landscape",
    x = "Tissue Compartment",
    y = "Patient ID"
  ) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    panel.grid = element_blank()
  )

ggsave("results/figures/08_patient_diversity_heatmap.png", plot = p_heat, width = 6, height = 8, dpi = 300)

cat("\n==> [08_patient_level.R] Completed successfully!")
cat("\n   - Exported Figure 1: results/figures/08_patient_shannon_diversity_by_tissue.png")
cat("\n   - Exported Figure 2: results/figures/08_patient_diversity_heatmap.png\n\n")