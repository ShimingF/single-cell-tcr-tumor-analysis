# ==============================================================================
# Script: 11_figure_common.R
# Purpose: Shared definitions and helper functions for Figures 1-5
# ==============================================================================

EXPANSION_LEVELS <- c(
  "Singleton", "Small", "Intermediate", "Expanded",
  "Highly expanded", "Hyper-expanded"
)

EXPANSION_COLORS <- c(
  "Singleton"       = "#ffffb2",
  "Small"           = "#fed976",
  "Intermediate"    = "#feb24c",
  "Expanded"        = "#fd8d3c",
  "Highly expanded" = "#f03b20",
  "Hyper-expanded"  = "#bd0026"
)

TISSUE_LEVELS <- c("Blood", "NAT", "Tumor")
TISSUE_COLORS <- c(
  "Blood" = "#1b9e77",
  "NAT"   = "#7570b3",
  "Tumor" = "#d95f02"
)

assign_expansion_tier <- function(clone_size) {
  tier <- dplyr::case_when(
    clone_size == 1 ~ "Singleton",
    clone_size >= 2  & clone_size <= 4  ~ "Small",
    clone_size >= 5  & clone_size <= 9  ~ "Intermediate",
    clone_size >= 10 & clone_size <= 49 ~ "Expanded",
    clone_size >= 50 & clone_size <= 99 ~ "Highly expanded",
    clone_size >= 100 ~ "Hyper-expanded",
    TRUE ~ NA_character_
  )

  factor(tier, levels = EXPANSION_LEVELS)
}

log2_or_ha <- function(a, b, c, d, correction = 0.5) {
  # Haldane-Anscombe correction. Here the 2x2 table is interpreted as:
  #             Group 1   Group 2
  # Feature       a         c
  # Other         b         d
  # so OR = (a*d)/(b*c).
  log2(
    ((a + correction) * (d + correction)) /
      ((b + correction) * (c + correction))
  )
}

safe_wilcox_p <- function(x, min_n = 3) {
  x <- x[is.finite(x) & !is.na(x)]
  if (length(x) < min_n || length(unique(x)) < 2) {
    return(NA_real_)
  }

  suppressWarnings(
    tryCatch(
      stats::wilcox.test(x, mu = 0, exact = FALSE)$p.value,
      error = function(e) NA_real_
    )
  )
}

project_theme <- function(base_size = 10) {
  ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 11),
      plot.subtitle = ggplot2::element_text(size = 8.5, color = "grey30"),
      panel.grid.minor = ggplot2::element_line(color = "grey95", linewidth = 0.25),
      panel.grid.major = ggplot2::element_line(color = "grey90", linewidth = 0.35)
    )
}

save_figure_pair <- function(plot, basename, width, height, dpi = 300) {
  png_path <- paste0(basename, ".png")
  svg_path <- paste0(basename, ".svg")

  ggplot2::ggsave(
    png_path, plot = plot, width = width, height = height,
    units = "in", dpi = dpi, bg = "white"
  )
  ggplot2::ggsave(
    svg_path, plot = plot, width = width, height = height,
    units = "in", bg = "white"
  )

  invisible(c(png = png_path, svg = svg_path))
}
