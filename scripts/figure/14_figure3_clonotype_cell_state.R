# ==============================================================================
# Script: 14_figure3_clonotype_cell_state.R
# Purpose: Generate Figure 3: Integrated scRNA-seq + TCR Repertoire Analysis
#
# Design principle:
#   - UMAP and raw compositions are descriptive and may pool cells.
#   - Cell-state association is estimated within each patient and summarized
#     across patients, avoiding pooled-cell pseudoreplication for inference.
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
cat(" Figure 3: Clonotype x T-Cell Transcriptional State\n")
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
required_columns <- c("patient", "source", "patient_clonotype_clean", "ident")
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
    patient_clonotype_clean = as.character(patient_clonotype_clean),
    ident = as.character(ident)
  ) %>%
  filter(
    !is.na(patient),
    !is.na(source),
    source %in% TISSUE_LEVELS,
    !is.na(patient_clonotype_clean),
    patient_clonotype_clean != "",
    patient_clonotype_clean != "NA",
    !is.na(ident),
    ident != ""
  )

clonotype_global_sizes <- df %>%
  count(patient, patient_clonotype_clean, name = "global_clone_size") %>%
  mutate(global_expansion_status = assign_expansion_tier(global_clone_size))

df <- df %>%
  left_join(
    clonotype_global_sizes,
    by = c("patient", "patient_clonotype_clean")
  )

# 3. Panel A: Patient-Stratified Cell-State Association -------------------------
MIN_PATIENTS_FOR_STATE_TEST <- 5
# Bubble size remains a pooled descriptive quantity (% of cells within a tier).
# Bubble color is the median patient-level log2 OR. Statistical evidence is a
# one-sample Wilcoxon test across patient-level log2 ORs, BH-adjusted globally.
pooled_fraction_df <- df %>%
  count(global_expansion_status, ident, name = "cell_count") %>%
  complete(
    global_expansion_status = factor(EXPANSION_LEVELS, levels = EXPANSION_LEVELS),
    ident = sort(unique(df$ident)),
    fill = list(cell_count = 0)
  ) %>%
  group_by(global_expansion_status) %>%
  mutate(
    tier_total_cells = sum(cell_count),
    fraction_within_tier = ifelse(tier_total_cells > 0, cell_count / tier_total_cells * 100, 0)
  ) %>%
  ungroup()

patient_totals <- df %>% count(patient, name = "patient_total")
patient_tier_totals <- df %>%
  count(patient, global_expansion_status, name = "tier_total")
patient_state_totals <- df %>%
  count(patient, ident, name = "state_total")

patient_or_df <- df %>%
  count(patient, global_expansion_status, ident, name = "a") %>%
  complete(
    patient = unique(df$patient),
    global_expansion_status = factor(EXPANSION_LEVELS, levels = EXPANSION_LEVELS),
    ident = sort(unique(df$ident)),
    fill = list(a = 0)
  ) %>%
  left_join(patient_totals, by = "patient") %>%
  left_join(patient_tier_totals, by = c("patient", "global_expansion_status")) %>%
  left_join(patient_state_totals, by = c("patient", "ident")) %>%
  mutate(
    tier_total = replace_na(tier_total, 0),
    state_total = replace_na(state_total, 0),
    b = tier_total - a,
    c = state_total - a,
    d = patient_total - a - b - c,
    informative =
      tier_total > 0 & tier_total < patient_total &
      state_total > 0 & state_total < patient_total,
    log2_odds_ratio = ifelse(informative, log2_or_ha(a, b, c, d), NA_real_)
  )

effect_summary <- patient_or_df %>%
  group_by(global_expansion_status, ident) %>%
  summarise(
    n_patients = sum(is.finite(log2_odds_ratio)),
    median_log2_or = ifelse(
      n_patients > 0,
      median(log2_odds_ratio, na.rm = TRUE),
      NA_real_
    ),
    mean_log2_or = ifelse(
      n_patients > 0,
      mean(log2_odds_ratio, na.rm = TRUE),
      NA_real_
    ),
    p_value = ifelse(
      n_patients >= MIN_PATIENTS_FOR_STATE_TEST,
      safe_wilcox_p(log2_odds_ratio, min_n = MIN_PATIENTS_FOR_STATE_TEST),
      NA_real_
    ),
    .groups = "drop"
  ) %>%
  mutate(
    p_adj = p.adjust(p_value, method = "BH"),
    fdr_status = case_when(
      n_patients < MIN_PATIENTS_FOR_STATE_TEST ~ "Insufficient coverage",
      !is.na(p_adj) & p_adj < 0.05 ~ "FDR < 0.05",
      TRUE ~ "Not significant"
    )
  ) %>%
  left_join(
    pooled_fraction_df %>% select(global_expansion_status, ident, fraction_within_tier),
    by = c("global_expansion_status", "ident")
  )

write_csv(patient_or_df, file.path(table_dir, "Figure3_patient_level_state_odds_ratios.csv"))
write_csv(effect_summary, file.path(table_dir, "Figure3_cell_state_expansion_enrichment.csv"))

p3a <- ggplot(effect_summary, aes(x = ident, y = global_expansion_status)) +
  geom_point(
    aes(
      size = fraction_within_tier,
      fill = median_log2_or,
      color = fdr_status
    ),
    shape = 21, alpha = 0.94, stroke = 0.85
  ) +
  # A second black ring makes cross-patient statistical support immediately
  # visible without changing the effect-size color scale.
  geom_point(
    data = effect_summary %>% filter(fdr_status == "FDR < 0.05"),
    aes(size = fraction_within_tier, fill = median_log2_or),
    shape = 21, color = "black", alpha = 1, stroke = 1.45,
    show.legend = FALSE
  ) +
  scale_size_continuous(range = c(1.5, 8), name = "% within Tier") +
  scale_fill_gradient2(
    low = "#2166ac", mid = "#f7f7f7", high = "#b2182c",
    midpoint = 0, na.value = "grey90",
    name = "Median patient\nLog2(Odds Ratio)"
  ) +
  scale_color_manual(
    values = c(
      "FDR < 0.05" = "black",
      "Not significant" = "grey78",
      "Insufficient coverage" = "grey92"
    ),
    name = "Patient-level evidence",
    drop = FALSE
  ) +
  labs(
    title = "A. Patient-stratified cell-state associations",
    subtitle = paste0(
      "Fill = median patient log2 OR; thick black ring = BH-FDR < 0.05; ",
      "tests require >= ", MIN_PATIENTS_FOR_STATE_TEST, " informative patients"
    ),
    x = "T-cell transcriptional state (ident)",
    y = "Global clone expansion class"
  ) +
  project_theme() +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1, size = 9),
    panel.grid.major = element_line(color = "grey92"),
    legend.position = "right"
  )

# 4. Panel B: UMAP Expansion Space Overlay (Descriptive) -----------------------
umap_x_col <- intersect(c("UMAP_1", "umap_1", "Umap_1"), colnames(df))
umap_y_col <- intersect(c("UMAP_2", "umap_2", "Umap_2"), colnames(df))

if (length(umap_x_col) > 0 && length(umap_y_col) > 0) {
  df_umap <- df %>%
    mutate(
      UMAP_X = .data[[umap_x_col[1]]],
      UMAP_Y = .data[[umap_y_col[1]]]
    ) %>%
    filter(!is.na(UMAP_X), !is.na(UMAP_Y)) %>%
    arrange(global_expansion_status)
  
  p3b <- ggplot(df_umap, aes(x = UMAP_X, y = UMAP_Y, color = global_expansion_status)) +
    geom_point(size = 0.4, alpha = 0.58) +
    scale_color_manual(values = EXPANSION_COLORS, drop = FALSE, name = "Expansion Tier") +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
    labs(
      title = "B. Expansion tiers on the T-cell UMAP",
      subtitle = "Descriptive pooled-cell visualization",
      x = "UMAP 1",
      y = "UMAP 2"
    ) +
    project_theme() +
    theme(
      panel.grid = element_blank(),
      legend.position = "right"
    )
} else {
  p3b <- ggplot(df, aes(x = ident, fill = global_expansion_status)) +
    geom_bar(position = "fill", color = "black", linewidth = 0.2, width = 0.7) +
    scale_y_continuous(labels = percent_format()) +
    scale_fill_manual(values = EXPANSION_COLORS, drop = FALSE, name = "Expansion Tier") +
    labs(
      title = "B. Expansion tier composition per T-cell functional cluster",
      subtitle = "UMAP coordinates were unavailable; displaying pooled composition instead",
      x = "T-cell transcriptional state (ident)",
      y = "Cell proportion (%)"
    ) +
    project_theme() +
    theme(
      axis.text.x = element_text(angle = 35, hjust = 1, size = 9),
      legend.position = "right"
    )
}

# 5. Panel C: Raw Expansion-Tier Composition per Cell State --------------------
state_tier_df <- df %>%
  count(ident, global_expansion_status, name = "cell_count")

state_order <- state_tier_df %>%
  group_by(ident) %>%
  summarise(total_cells = sum(cell_count), .groups = "drop") %>%
  arrange(desc(total_cells)) %>%
  pull(ident)

state_tier_df <- state_tier_df %>%
  mutate(ident = factor(ident, levels = state_order))

p3c <- ggplot(state_tier_df, aes(x = ident, y = cell_count, fill = global_expansion_status)) +
  geom_col(color = "black", linewidth = 0.2, width = 0.72) +
  scale_fill_manual(values = EXPANSION_COLORS, drop = FALSE, name = "Expansion Tier") +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.05))) +
  labs(
    title = "C. Raw expansion-tier composition by T-cell state",
    subtitle = "Descriptive absolute single-cell counts; no inferential claim is made from this panel",
    x = "T-cell transcriptional state (ident)",
    y = "Number of single cells",
    fill = "Expansion Tier"
  ) +
  project_theme() +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1, size = 9),
    legend.position = "right"
  )

# 6. Combine 3 Panels & Export --------------------------------------------------
figure3_3panel <- (p3a | p3b) / p3c +
  plot_layout(heights = c(1.1, 0.9)) +
  plot_annotation(
    title = "Figure 3. Association between T-cell clonal expansion and transcriptional cell states",
    subtitle = "Patient-stratified inference integrated with descriptive scRNA-seq/TCR visualizations",
    theme = theme(
      plot.title = element_text(face = "bold", size = 15),
      plot.subtitle = element_text(size = 10, color = "grey30")
    )
  )

paths <- save_figure_pair(
  figure3_3panel,
  file.path(figure_dir, "Figure3_clonotype_cell_state"),
  width = 12.5,
  height = 10
)

cat("============================================================\n")
cat(" Figure 3 completed successfully.\n")
cat(" PNG:", paths[["png"]], "\n")
cat(" SVG:", paths[["svg"]], "\n")
cat("============================================================\n\n")
