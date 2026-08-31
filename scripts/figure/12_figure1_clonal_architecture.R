# ==============================================================================
# Script: 12_figure1_clonal_architecture.R
# Purpose: Generate Figure 1 (Panel A, B, C) for T-cell Clonal Architecture
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(patchwork)
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

cat("==> [12_figure1_clonal_architecture.R] Generating Figure 1 Panels...\n")

# 1. Load Master Integrated Dataset --------------------------------------------
master_path <- "data/processed/final_paired_tcell_metadata.rds"
figure_dir  <- "results/figures"
table_dir   <- "results/tables"

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(master_path)) {
  stop("Master dataset not found! Please run upstream preprocessing scripts first.")
}

master_df <- readRDS(master_path)

required_columns <- c("patient", "patient_clonotype_clean")
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
    patient_clonotype_clean = as.character(patient_clonotype_clean)
  ) %>%
  filter(
    !is.na(patient),
    !is.na(patient_clonotype_clean),
    patient_clonotype_clean != "",
    patient_clonotype_clean != "NA"
  )

# 2. Recompute Patient-Specific Global Clone Size ------------------------------
# All downstream figures use exactly the same patient x clonotype definition.
clone_summary <- df %>%
  count(patient, patient_clonotype_clean, name = "global_clone_size") %>%
  mutate(expansion_status = assign_expansion_tier(global_clone_size))

# Fail fast if an upstream clone-size field exists but disagrees with the
# canonical patient-specific definition used by this figure set.
if ("true_clone_size" %in% colnames(df)) {
  upstream_check <- df %>%
    group_by(patient, patient_clonotype_clean) %>%
    summarise(
      upstream_clone_size = first(true_clone_size),
      n_upstream_sizes = n_distinct(true_clone_size),
      .groups = "drop"
    ) %>%
    left_join(clone_summary, by = c("patient", "patient_clonotype_clean"))
  
  if (any(upstream_check$n_upstream_sizes > 1, na.rm = TRUE)) {
    stop("Upstream true_clone_size is not constant within patient-clonotype.")
  }
  
  n_mismatch <- upstream_check %>%
    filter(!is.na(upstream_clone_size), upstream_clone_size != global_clone_size) %>%
    nrow()
  
  if (n_mismatch > 0) {
    stop(
      "Detected ", n_mismatch,
      " patient-clonotypes whose upstream true_clone_size disagrees with the canonical definition."
    )
  }
}

write_csv(clone_summary, file.path(table_dir, "Figure1_clonotype_summary.csv"))

# 3. Panel A: Log-Log Rank-Abundance Plot --------------------------------------
rank_df <- clone_summary %>%
  arrange(desc(global_clone_size)) %>%
  mutate(rank = row_number())

p1a <- ggplot(rank_df, aes(x = rank, y = global_clone_size, color = expansion_status)) +
  geom_point(alpha = 0.82, size = 1.5) +
  scale_x_log10(labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(labels = trans_format("log10", math_format(10^.x))) +
  scale_color_manual(values = EXPANSION_COLORS, drop = FALSE) +
  labs(
    title = "A. Clone-size rank-abundance distribution",
    x = "Clonotype rank",
    y = "Clone size (cells)",
    color = "Expansion Tier"
  ) +
  project_theme() +
  theme(legend.position = "right")

# 4. Panel B: Expansion Class Distribution (Clonotype Level) -------------------
p1b_df <- rank_df %>% count(expansion_status, .drop = FALSE)

p1b <- ggplot(p1b_df, aes(x = expansion_status, y = n, fill = expansion_status)) +
  geom_col(show.legend = FALSE, color = "black", linewidth = 0.2) +
  geom_text(aes(label = comma(n)), vjust = -0.3, size = 3) +
  scale_fill_manual(values = EXPANSION_COLORS, drop = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(
    title = "B. Clonotype abundance by tier",
    x = "Expansion Tier",
    y = "Number of unique patient-clonotypes"
  ) +
  project_theme() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

# 5. Panel C: Cellular Mass Contribution ---------------------------------------
# Use an explicit canonical column name here so this panel is immune to any
# upstream column named expansion_status that may still exist in master_df.
clone_lookup <- clone_summary %>%
  transmute(
    patient,
    patient_clonotype_clean,
    canonical_expansion_status = expansion_status
  )

cell_contrib_df <- df %>%
  left_join(
    clone_lookup,
    by = c("patient", "patient_clonotype_clean")
  ) %>%
  filter(!is.na(canonical_expansion_status)) %>%
  count(canonical_expansion_status, .drop = FALSE) %>%
  rename(expansion_status = canonical_expansion_status) %>%
  mutate(
    expansion_status = factor(expansion_status, levels = EXPANSION_LEVELS),
    percentage = n / sum(n) * 100,
    dataset = "All Paired T-Cells"
  )

p1c <- ggplot(cell_contrib_df, aes(x = dataset, y = percentage, fill = expansion_status)) +
  geom_col(color = "black", linewidth = 0.3, width = 0.4) +
  geom_text(
    aes(label = sprintf("%.1f%%", percentage)),
    position = position_stack(vjust = 0.5),
    size = 3
  ) +
  scale_fill_manual(values = EXPANSION_COLORS, drop = FALSE) +
  labs(
    title = "C. Cellular mass by tier",
    x = "",
    y = "Percentage of total cells (%)",
    fill = "Expansion Tier"
  ) +
  project_theme() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none"
  )

# 6. Combine Panels & Export ----------------------------------------------------
figure_subtitle <- sprintf(
  "Analysis of %s strict paired T cells across %d patients; clonotypes are patient-specific",
  comma(nrow(df)), n_distinct(df$patient)
)

bottom_row <- (p1b | p1c) +
  plot_layout(widths = c(1.25, 0.75))

p_fig1_final <- (p1a / bottom_row) +
  plot_layout(heights = c(1.15, 1)) +
  plot_annotation(
    title = "Figure 1. T-cell clonal architecture and cellular mass contribution",
    subtitle = figure_subtitle,
    theme = theme(
      plot.title = element_text(face = "bold", size = 15),
      plot.subtitle = element_text(size = 10, color = "grey30")
    )
  )

paths <- save_figure_pair(
  p_fig1_final,
  file.path(figure_dir, "Figure1_clonal_architecture_panels"),
  width = 10.5,
  height = 8.5
)

cat("\n==> Figure 1 generated successfully.\n")
cat("    PNG:", paths[["png"]], "\n")
cat("    SVG:", paths[["svg"]], "\n\n")
