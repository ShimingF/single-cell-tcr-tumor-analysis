# ==============================================================================
# Script: 05_clonal_architecture.R
# Purpose: Re-calculate true clonal expansion with updated 6-level status categories
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(tidyr)
})

cat("==> [05_clonal_architecture.R] Starting Clonal Architecture Re-analysis...\n")

# 1. Load Master Integrated Dataset
master_path <- "data/processed/final_paired_tcell_metadata.rds"

if (!file.exists(master_path)) {
  stop("Master dataset not found! Please run scripts/04_tcr_metadata_integration.R first.")
}

master_df <- readRDS(master_path)

cat(sprintf("   - Loaded Master Dataset: %d strict 1a1b paired cells.\n", nrow(master_df)))

# 2. Define Clones per Patient & Re-calculate True Clone Size with New Expansion Status
clonal_summary <- master_df %>%
  group_by(patient, tcr_sequence_id, CDR3a, CDR3b) %>%
  summarise(
    true_clone_size = n(),
    n_tumor = sum(source == "Tumor", na.rm = TRUE),
    n_nat   = sum(source == "NAT", na.rm = TRUE),
    n_blood = sum(source == "Blood", na.rm = TRUE),
    cell_states = paste(unique(ident), collapse = ";"),
    .groups = "drop"
  ) %>%
  arrange(patient, desc(true_clone_size)) %>%
  group_by(patient) %>%
  mutate(
    patient_clonotype_clean = paste(patient, "clone", row_number(), sep = "_")
  ) %>%
  ungroup() %>%
  mutate(
    # Updated 6-tier expansion status classification
    expansion_status = factor(
      case_when(
        true_clone_size == 1 ~ "Singleton",
        true_clone_size >= 2  & true_clone_size <= 4  ~ "Small",
        true_clone_size >= 5  & true_clone_size <= 9  ~ "Intermediate",
        true_clone_size >= 10 & true_clone_size <= 49 ~ "Expanded",
        true_clone_size >= 50 & true_clone_size <= 99 ~ "Highly expanded",
        true_clone_size >= 100                        ~ "Hyper-expanded"
      ),
      levels = c("Singleton", "Small", "Intermediate", "Expanded", "Highly expanded", "Hyper-expanded")
    )
  ) %>%
  arrange(desc(true_clone_size))

# 3. Print Breakdown & Top Expanded Clones
cat("\n--- Clonal Expansion Status Breakdown ---\n")
print(count(clonal_summary, expansion_status) %>% mutate(pct = round(n / sum(n) * 100, 2)))

cat("\n--- Top 10 Expanded Clones ---\n")
print(head(clonal_summary %>% select(patient, patient_clonotype_clean, true_clone_size, expansion_status, n_tumor, n_nat, n_blood), 10))

# 4. Merge Clean Clone Metrics Back to Master Dataset
master_df <- master_df %>%
  select(-any_of(c("true_clone_size", "expansion_status", "patient_clonotype_clean"))) %>%
  left_join(
    clonal_summary %>% select(patient, tcr_sequence_id, true_clone_size, expansion_status, patient_clonotype_clean),
    by = c("patient", "tcr_sequence_id")
  )

# 5. Export Summary Table & Updated Master RDS
dir.create("results/tables", showWarnings = FALSE, recursive = TRUE)
write_csv(clonal_summary, "results/tables/recomputed_clonal_summary.csv")
saveRDS(master_df, "data/processed/final_paired_tcell_metadata.rds")

# 6. Plot Expansion Architecture
dir.create("results/figures", showWarnings = FALSE, recursive = TRUE)

p_expansion <- ggplot(clonal_summary, aes(x = expansion_status, fill = expansion_status)) +
  geom_bar() +
  scale_fill_brewer(palette = "YlOrRd") +
  theme_minimal() +
  labs(
    title = "Clonal Expansion Architecture (Strict 1a1b Filtered)",
    x = "Expansion Category",
    y = "Number of Unique Clonotypes"
  ) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    legend.position = "none"
  )

ggsave("results/figures/05_clonal_expansion_distribution.png", plot = p_expansion, width = 7, height = 5)

cat("==> [05_clonal_architecture.R] Completed successfully!\n\n")