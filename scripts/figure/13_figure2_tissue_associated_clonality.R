# ==============================================================================
# Script: 13_figure2_tissue_associated_clonality.R
# Purpose: Generate Figure 2: Tissue-Associated TCR Clonal Architecture
#
# Design principle:
#   - Descriptive panels may pool cells/clonotypes across the dataset.
#   - Inferential tumor-enrichment analysis is patient-aware and compares
#     matched Tumor vs NAT samples within each patient.
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
cat(" Figure 2: Tissue-Associated Clonal Architecture\n")
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
cat("Input dataset loaded. Total cells:", nrow(master_df), "\n\n")

# 2. Basic Cleaning & Canonical Expansion Tier Definition ----------------------
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

clonotype_global_df <- df %>%
  count(patient, patient_clonotype_clean, source, name = "cell_count") %>%
  pivot_wider(names_from = source, values_from = cell_count, values_fill = 0) %>%
  mutate(
    global_clone_size = Blood + NAT + Tumor,
    global_expansion_status = assign_expansion_tier(global_clone_size)
  )

df <- df %>%
  left_join(
    clonotype_global_df %>%
      select(patient, patient_clonotype_clean, global_clone_size, global_expansion_status),
    by = c("patient", "patient_clonotype_clean")
  )

# 3. Panel A: Tissue-Specific Clone-Size Tail Distribution ---------------------
# Clone size is discrete and strongly right-skewed; a CCDF displays the tail
# directly and avoids violin-density distortion around the large singleton mass.
clone_tissue_df <- df %>%
  count(patient, source, patient_clonotype_clean, name = "clone_size")

ccdf_df <- clone_tissue_df %>%
  count(source, clone_size, name = "n_clones") %>%
  group_by(source) %>%
  arrange(desc(clone_size), .by_group = TRUE) %>%
  mutate(ccdf = cumsum(n_clones) / sum(n_clones)) %>%
  arrange(source, clone_size) %>%
  ungroup()

p2a <- ggplot(ccdf_df, aes(x = clone_size, y = ccdf, color = source)) +
  geom_step(linewidth = 0.9, direction = "hv") +
  scale_x_log10(labels = comma) +
  scale_y_log10(
    breaks = c(1, 0.1, 0.01, 0.001, 0.0001),
    labels = c("100%", "10%", "1%", "0.1%", "0.01%")
  ) +
  scale_color_manual(values = TISSUE_COLORS, drop = FALSE) +
  labs(
    title = "A. Tissue-specific clone-size tail distribution",
    subtitle = "Complementary cumulative distribution: P(Clone size >= x)",
    x = "Clone size (cells, log scale)",
    y = "Fraction of clonotypes >= x (log scale)",
    color = "Tissue"
  ) +
  project_theme() +
  theme(legend.position = "right")

# 4. Panel B: Patient-Level Repertoire Clonality -------------------------------
patient_clonality <- df %>%
  count(patient, source, patient_clonotype_clean, name = "clone_size") %>%
  group_by(patient, source) %>%
  summarise(
    total_cells = sum(clone_size),
    unique_clones = n(),
    shannon = -sum((clone_size / sum(clone_size)) * log(clone_size / sum(clone_size))),
    .groups = "drop"
  ) %>%
  mutate(
    pielou_evenness = ifelse(unique_clones > 1, shannon / log(unique_clones), 0),
    clonality = ifelse(unique_clones > 1, 1 - pielou_evenness, 1)
  )

write_csv(patient_clonality, file.path(table_dir, "Figure2_patient_clonality.csv"))

p2b <- ggplot(patient_clonality, aes(x = source, y = clonality)) +
  geom_boxplot(
    aes(group = source), width = 0.35, fill = "white", color = "black",
    linewidth = 0.4, outlier.shape = NA
  ) +
  geom_line(
    data = patient_clonality %>% group_by(patient) %>% filter(n() >= 2) %>% ungroup(),
    aes(group = patient), color = "grey70", alpha = 0.6, linewidth = 0.45
  ) +
  geom_jitter(aes(color = source), width = 0.06, size = 2, alpha = 0.9) +
  scale_color_manual(values = TISSUE_COLORS, drop = FALSE) +
  scale_y_continuous(limits = c(0, 1.05)) +
  labs(
    title = "B. Patient-level repertoire clonality",
    subtitle = "Descriptive comparison; each point represents one patient-tissue sample",
    x = "Tissue compartment",
    y = "TCR clonality (1 - Pielou)"
  ) +
  project_theme() +
  theme(legend.position = "none")

# 5. Panel C: Global Tissue Occupancy Across Expansion Classes -----------------
occupancy_df <- df %>%
  count(global_expansion_status, source, name = "n_cells", .drop = FALSE) %>%
  group_by(global_expansion_status) %>%
  mutate(percentage = n_cells / sum(n_cells) * 100) %>%
  ungroup()

p2c <- ggplot(occupancy_df, aes(x = global_expansion_status, y = percentage, fill = source)) +
  geom_col(color = "black", linewidth = 0.2, width = 0.7) +
  geom_text(
    aes(label = ifelse(percentage >= 5, sprintf("%.0f%%", percentage), "")),
    position = position_stack(vjust = 0.5), size = 2.8,
    color = "white", fontface = "bold"
  ) +
  scale_y_continuous(labels = function(x) paste0(x, "%"), limits = c(0, 100.1), expand = c(0, 0)) +
  scale_fill_manual(values = TISSUE_COLORS) +
  labs(
    title = "C. Tissue occupancy across global expansion classes",
    subtitle = "Descriptive pooled-cell composition; not an enrichment test",
    x = "Global clone expansion class",
    y = "Proportion of cells (%)",
    fill = "Tissue"
  ) +
  project_theme() +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    legend.position = "bottom"
  )

# 6. Panel D: Patient-Aware Tumor vs NAT Enrichment ----------------------------
# For every patient and expansion tier, construct:
#                       Tier       Other tiers
#   Tumor                 a             b
#   NAT                   c             d
# and report the Haldane-Anscombe corrected log2 odds ratio. Patient-level
# values are the inferential units; confidence intervals are across patients.
paired_totals <- df %>%
  filter(source %in% c("NAT", "Tumor")) %>%
  count(patient, source, name = "total_cells") %>%
  pivot_wider(names_from = source, values_from = total_cells, values_fill = 0) %>%
  filter(NAT > 0, Tumor > 0)

paired_patients <- paired_totals$patient

tier_counts_observed <- df %>%
  filter(patient %in% paired_patients, source %in% c("NAT", "Tumor")) %>%
  count(patient, source, global_expansion_status, name = "tier_cells") %>%
  mutate(
    source = as.character(source),
    global_expansion_status = as.character(global_expansion_status)
  )

tier_counts_patient <- tidyr::expand_grid(
  patient = paired_patients,
  source = c("NAT", "Tumor"),
  global_expansion_status = EXPANSION_LEVELS
) %>%
  left_join(
    tier_counts_observed,
    by = c("patient", "source", "global_expansion_status")
  ) %>%
  mutate(
    tier_cells = replace_na(tier_cells, 0),
    global_expansion_status = factor(global_expansion_status, levels = EXPANSION_LEVELS)
  ) %>%
  pivot_wider(names_from = source, values_from = tier_cells, values_fill = 0) %>%
  left_join(paired_totals, by = "patient", suffix = c("_tier", "_total")) %>%
  mutate(
    a = Tumor_tier,
    b = Tumor_total - Tumor_tier,
    c = NAT_tier,
    d = NAT_total - NAT_tier,
    informative = (a + c) > 0,
    log2_or_tumor_vs_nat = ifelse(informative, log2_or_ha(a, b, c, d), NA_real_)
  )

te_summary <- tier_counts_patient %>%
  group_by(global_expansion_status) %>%
  summarise(
    n_patients = sum(is.finite(log2_or_tumor_vs_nat)),
    mean_log2_or = mean(log2_or_tumor_vs_nat, na.rm = TRUE),
    median_log2_or = median(log2_or_tumor_vs_nat, na.rm = TRUE),
    sd_log2_or = sd(log2_or_tumor_vs_nat, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    se_log2_or = sd_log2_or / sqrt(n_patients),
    t_crit = ifelse(n_patients > 1, qt(0.975, df = n_patients - 1), NA_real_),
    ci_lower = mean_log2_or - t_crit * se_log2_or,
    ci_upper = mean_log2_or + t_crit * se_log2_or,
    n_label = paste0("n = ", n_patients)
  )

write_csv(tier_counts_patient, file.path(table_dir, "Figure2_patient_tier_tumor_vs_nat_enrichment.csv"))

# With only 8-14 patients per tier, raw patient points are more informative
# than kernel-density violins. Use a symmetric, data-adaptive y-axis.
or_extent <- max(
  abs(tier_counts_patient$log2_or_tumor_vs_nat),
  abs(te_summary$ci_lower),
  abs(te_summary$ci_upper),
  na.rm = TRUE
)
or_extent <- max(1, or_extent * 1.08)
te_summary <- te_summary %>% mutate(label_y = or_extent * 1.06)

p2d <- ggplot(tier_counts_patient, aes(x = global_expansion_status, y = log2_or_tumor_vs_nat)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", linewidth = 0.6) +
  geom_jitter(
    aes(color = global_expansion_status), width = 0.09,
    size = 2.1, alpha = 0.85, na.rm = TRUE
  ) +
  geom_errorbar(
    data = te_summary,
    aes(x = global_expansion_status, ymin = ci_lower, ymax = ci_upper),
    inherit.aes = FALSE, width = 0.14, linewidth = 0.75, color = "black"
  ) +
  geom_point(
    data = te_summary,
    aes(x = global_expansion_status, y = mean_log2_or),
    inherit.aes = FALSE, shape = 23, size = 3.2,
    fill = "white", color = "black", stroke = 0.8
  ) +
  geom_point(
    data = te_summary,
    aes(x = global_expansion_status, y = median_log2_or),
    inherit.aes = FALSE, shape = 21, size = 1.9,
    fill = "black", color = "white"
  ) +
  geom_text(
    data = te_summary,
    aes(x = global_expansion_status, y = label_y, label = n_label),
    inherit.aes = FALSE, size = 2.7, fontface = "bold", color = "grey20"
  ) +
  scale_color_manual(values = EXPANSION_COLORS, drop = FALSE) +
  scale_y_continuous(
    limits = c(-or_extent, or_extent * 1.13),
    breaks = pretty(c(-or_extent, or_extent), n = 7),
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  labs(
    title = "D. Patient-level tumor enrichment by expansion class",
    subtitle = "Tumor vs matched NAT; points = patients, diamond = mean with 95% t-CI, black dot = median",
    x = "Global clone expansion class",
    y = "Log2 odds ratio (Tumor vs NAT)"
  ) +
  project_theme() +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    legend.position = "none"
  )

# 7. Save Descriptive Tables ----------------------------------------------------
write_csv(
  clone_tissue_df %>%
    group_by(source) %>%
    summarise(
      n_clonotypes = n(),
      median_clone_size = median(clone_size),
      mean_clone_size = mean(clone_size),
      max_clone_size = max(clone_size),
      .groups = "drop"
    ),
  file.path(table_dir, "Figure2_clone_statistics.csv")
)

# 8. Combine 4 Panels & Export --------------------------------------------------
figure2_4panel <- (p2a | p2b) / (p2c | p2d) +
  plot_annotation(
    title = "Figure 2. Tissue-associated TCR clonal architecture and tumor enrichment",
    subtitle = "Pooled descriptive views are separated from patient-aware Tumor-vs-NAT inference",
    theme = theme(
      plot.title = element_text(face = "bold", size = 15),
      plot.subtitle = element_text(size = 10, color = "grey30")
    )
  )

paths <- save_figure_pair(
  figure2_4panel,
  file.path(figure_dir, "Figure2_tissue_associated_clonality"),
  width = 12,
  height = 9.5
)

cat("============================================================\n")
cat(" Figure 2 completed successfully.\n")
cat(" PNG:", paths[["png"]], "\n")
cat(" SVG:", paths[["svg"]], "\n")
cat("============================================================\n\n")
