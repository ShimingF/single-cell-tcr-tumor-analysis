# ==============================================================================
# Script: 15_figure4_patient_reproducibility.R
# Purpose: Generate Figure 4: Patient-Level Heterogeneity and Clonal Architecture
#
# Terminology:
#   - "Clonally expanded" = global patient-specific clone size >= 2.
#   - "Expanded" remains reserved for the 10-49 cell expansion tier.
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(readr)
  library(scales)
})

# Locate shared helper script under either common project layout.
.common_candidates <- c(
  "scripts/figure/11_figure_common.R",
  "scripts/11_figure_common.R",
  "11_figure_common.R"
)
.common_path <- .common_candidates[file.exists(.common_candidates)][1]
if (is.na(.common_path)) {
  stop(
    "Cannot find 11_figure_common.R. Expected one of: ",
    paste(.common_candidates, collapse = ", ")
  )
}
source(.common_path)
rm(.common_candidates, .common_path)

cat("\n============================================================\n")
cat(" Figure 4: Patient-Level Heterogeneity (Clonally Expanded >= 2)\n")
cat("============================================================\n\n")

# 1. Input / Output Paths -------------------------------------------------------
master_path <- "data/processed/final_paired_tcell_metadata.rds"
figure_dir  <- "results/figures"
table_dir   <- "results/tables"

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(master_path)) {
  stop("Master dataset not found: ", master_path, "\nRun upstream preprocessing scripts first.")
}

master_df <- readRDS(master_path)

# 2. Data Cleaning & Canonical Expansion Tier Assignment -----------------------
required_columns <- c("patient", "source", "patient_clonotype_clean")
missing_cols <- setdiff(required_columns, colnames(master_df))
if (length(missing_cols) > 0) {
  stop("Missing columns: ", paste(missing_cols, collapse = ", "))
}

df <- master_df %>%
  # Remove any upstream derived columns that are recomputed canonically below.
  # This prevents dplyr joins from creating .x/.y suffixes and silently drifting definitions.
  select(-any_of(c("expansion_status", "global_expansion_status", "global_clone_size", "is_clonally_expanded"))) %>%
  mutate(
    patient = as.character(patient),
    source = as.character(source),
    patient_clonotype_clean = as.character(patient_clonotype_clean)
  ) %>%
  filter(
    !is.na(patient),
    !is.na(source),
    source %in% TISSUE_LEVELS,
    !is.na(patient_clonotype_clean),
    patient_clonotype_clean != "",
    patient_clonotype_clean != "NA"
  ) %>%
  mutate(source = factor(source, levels = TISSUE_LEVELS))

clonotype_global_sizes <- df %>%
  count(patient, patient_clonotype_clean, name = "global_clone_size") %>%
  mutate(global_expansion_status = assign_expansion_tier(global_clone_size))

df <- df %>%
  left_join(
    clonotype_global_sizes,
    by = c("patient", "patient_clonotype_clean")
  ) %>%
  mutate(is_clonally_expanded = global_clone_size >= 2)

# 3. Patient x Tissue Clonal-Expansion Summary ---------------------------------
patient_tissue_summary <- df %>%
  group_by(patient, source) %>%
  summarise(
    total_cells = n(),
    clonally_expanded_cells = sum(is_clonally_expanded),
    nonexpanded_cells = total_cells - clonally_expanded_cells,
    clonally_expanded_fraction = clonally_expanded_cells / total_cells * 100,
    .groups = "drop"
  ) %>%
  complete(
    patient = unique(df$patient),
    source = factor(TISSUE_LEVELS, levels = TISSUE_LEVELS),
    fill = list(
      total_cells = 0,
      clonally_expanded_cells = 0,
      nonexpanded_cells = 0,
      clonally_expanded_fraction = NA_real_
    )
  )

write_csv(
  patient_tissue_summary,
  file.path(table_dir, "Figure4_patient_tissue_expansion_summary.csv")
)

# 4. Panel A: Patient x Tissue Heatmap ------------------------------------------
patient_order <- patient_tissue_summary %>%
  filter(source == "Tumor") %>%
  arrange(desc(clonally_expanded_fraction)) %>%
  pull(patient)

patient_tissue_summary <- patient_tissue_summary %>%
  mutate(patient = factor(patient, levels = rev(patient_order)))

p4a <- ggplot(patient_tissue_summary, aes(x = source, y = patient, fill = clonally_expanded_fraction)) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(
    aes(label = ifelse(is.na(clonally_expanded_fraction), "N/A", sprintf("%.1f%%", clonally_expanded_fraction))),
    size = 3.2, fontface = "bold",
    color = ifelse(
      !is.na(patient_tissue_summary$clonally_expanded_fraction) &
        patient_tissue_summary$clonally_expanded_fraction > 50,
      "white", "black"
    )
  ) +
  scale_fill_distiller(
    palette = "YlOrRd", direction = 1, na.value = "#e0e0e0",
    name = "Cells in clonally\nexpanded clones (%)\n(Size >= 2)",
    labels = function(x) paste0(x, "%")
  ) +
  labs(
    title = "A. Patient-level clonal expansion heatmap",
    subtitle = "Fraction of T cells in patient-specific clones with global size >= 2",
    x = "Tissue compartment",
    y = "Patient ID"
  ) +
  project_theme() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(face = "bold", size = 9.5),
    axis.text.y = element_text(size = 9),
    legend.position = "right"
  )

# 5. Panel B: Tumor vs NAT Log2 Odds Ratio of Clonal Expansion -----------------
# This avoids unstable raw percentage ratios when a compartment has a small
# expanded fraction. Positive values indicate higher odds of a T cell belonging
# to a clonally expanded clone in Tumor than in matched NAT.
expansion_or_df <- patient_tissue_summary %>%
  filter(source %in% c("Tumor", "NAT"), total_cells > 0) %>%
  select(
    patient, source, total_cells,
    clonally_expanded_cells, nonexpanded_cells
  ) %>%
  pivot_wider(
    names_from = source,
    values_from = c(total_cells, clonally_expanded_cells, nonexpanded_cells)
  ) %>%
  filter(!is.na(total_cells_Tumor), !is.na(total_cells_NAT)) %>%
  mutate(
    log2_or = log2_or_ha(
      clonally_expanded_cells_Tumor,
      nonexpanded_cells_Tumor,
      clonally_expanded_cells_NAT,
      nonexpanded_cells_NAT
    ),
    patient = factor(as.character(patient), levels = patient_order)
  )

write_csv(expansion_or_df, file.path(table_dir, "Figure4_patient_tumor_vs_nat_log2_or.csv"))

p4b <- ggplot(expansion_or_df, aes(x = reorder(patient, log2_or), y = log2_or, fill = log2_or > 0)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", linewidth = 0.6) +
  geom_col(color = "black", linewidth = 0.25, width = 0.68, show.legend = FALSE) +
  geom_text(
    aes(label = sprintf("%.2f", log2_or), hjust = ifelse(log2_or >= 0, -0.15, 1.15)),
    size = 2.8, fontface = "bold"
  ) +
  scale_fill_manual(values = c("TRUE" = TISSUE_COLORS[["Tumor"]], "FALSE" = "#3182bd")) +
  coord_flip() +
  scale_y_continuous(expand = expansion(mult = c(0.18, 0.22))) +
  labs(
    title = "B. Tumor-vs-NAT odds of clonal expansion",
    subtitle = "Haldane-Anscombe log2 OR; >0 indicates higher tumor expansion",
    x = "Patient ID",
    y = "Log2 odds ratio (Tumor vs NAT)"
  ) +
  project_theme() +
  theme(
    axis.text.y = element_text(size = 8.5),
    panel.grid.major.y = element_blank()
  )

# 6. Panel C: Patient-Level Expansion Tier Breakdown ---------------------------
patient_tier_comp <- df %>%
  count(patient, source, global_expansion_status, name = "cell_count") %>%
  group_by(patient, source) %>%
  mutate(percentage = cell_count / sum(cell_count) * 100) %>%
  ungroup() %>%
  mutate(patient = factor(patient, levels = rev(patient_order)))

p4c <- ggplot(patient_tier_comp, aes(x = source, y = percentage, fill = global_expansion_status)) +
  geom_col(color = "black", linewidth = 0.15, width = 0.75) +
  scale_fill_manual(values = EXPANSION_COLORS, drop = FALSE, name = "Expansion Tier") +
  facet_wrap(~ patient, ncol = 7) +
  scale_y_continuous(labels = function(x) paste0(x, "%"), expand = c(0, 0)) +
  labs(
    title = "C. Patient-level expansion-tier composition across compartments",
    subtitle = "Descriptive cell composition across the six canonical global expansion tiers",
    x = "Tissue compartment",
    y = "Percentage of cells (%)",
    fill = "Expansion Tier"
  ) +
  project_theme() +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1, size = 8),
    axis.text.y = element_text(size = 8),
    strip.background = element_rect(fill = "grey92", color = "black"),
    strip.text = element_text(face = "bold", size = 8.5),
    legend.position = "bottom"
  )

# 7. Combine 3 Panels & Export --------------------------------------------------
figure4_3panel <- (p4a | p4b) / p4c +
  plot_layout(heights = c(1, 1.1)) +
  plot_annotation(
    title = "Figure 4. Cross-patient heterogeneity and compartment-specific clonal architecture",
    subtitle = "Patient-resolved evaluation of TCR clonal expansion across matched tissue compartments",
    theme = theme(
      plot.title = element_text(face = "bold", size = 15),
      plot.subtitle = element_text(size = 10, color = "grey30")
    )
  )

paths <- save_figure_pair(
  figure4_3panel,
  file.path(figure_dir, "Figure4_patient_reproducibility"),
  width = 13,
  height = 11
)

cat("============================================================\n")
cat(" Figure 4 completed successfully.\n")
cat(" PNG:", paths[["png"]], "\n")
cat(" SVG:", paths[["svg"]], "\n")
cat("============================================================\n\n")
