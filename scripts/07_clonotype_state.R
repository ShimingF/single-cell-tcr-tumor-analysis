# ==============================================================================
# Script: 07_clonotype_state.R
# Purpose: Map clone expansion status to single-cell functional states (ident)
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(tidyr)
})

cat("==> [07_clonotype_state.R] Starting Clonotype Cell-State Analysis...\n")

# 1. Load Master Integrated Dataset (Correct File Path)
master_path <- "data/processed/final_paired_tcell_metadata.rds"

if (!file.exists(master_path)) {
  stop("Master dataset not found! Please run scripts/04 and 05 first.")
}

master_df <- readRDS(master_path)

cat(sprintf("   - Loaded Master Dataset: %d strict 1a1b paired cells.\n", nrow(master_df)))

# 2. Cell State Composition per Expansion Status Category
state_expansion_matrix <- master_df %>%
  filter(!is.na(ident), !is.na(expansion_status)) %>%
  group_by(expansion_status, ident) %>%
  summarise(cell_count = n(), .groups = "drop") %>%
  group_by(expansion_status) %>%
  mutate(
    percentage = round(cell_count / sum(cell_count) * 100, 2)
  ) %>%
  ungroup()

# 3. Export Summary Tables
dir.create("results/tables", showWarnings = FALSE, recursive = TRUE)

# Wide-format matrix table (Expansion Tier x Cell State)
state_matrix_wide <- state_expansion_matrix %>%
  select(expansion_status, ident, cell_count) %>%
  pivot_wider(names_from = ident, values_from = cell_count, values_fill = 0)

write_csv(state_expansion_matrix, "results/tables/07_cell_state_by_expansion_long.csv")
write_csv(state_matrix_wide, "results/tables/07_cell_state_by_expansion_matrix.csv")

cat("\n--- Cell State Distribution across Expansion Tiers (Top States) ---\n")
print(head(state_matrix_wide, 6))

# 4. Plot Cell State Stacked Percentage Bar Plot
dir.create("results/figures", showWarnings = FALSE, recursive = TRUE)

p_state <- ggplot(state_expansion_matrix, aes(x = expansion_status, y = percentage, fill = ident)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  theme_minimal() +
  labs(
    title = "T-Cell State Composition Across Clonal Expansion Tiers",
    x = "Clonal Expansion Tier",
    y = "Cell Proportion",
    fill = "Cell State (ident)"
  ) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    legend.position = "right"
  )

ggsave("results/figures/07_clonotype_cell_state_composition.png", plot = p_state, width = 9, height = 6)

cat("\n==> [07_clonotype_state.R] Completed successfully!\n")
cat("   - Exported Table:  results/tables/07_cell_state_by_expansion_matrix.csv\n")
cat("   - Exported Figure: results/figures/07_clonotype_cell_state_composition.png\n\n")