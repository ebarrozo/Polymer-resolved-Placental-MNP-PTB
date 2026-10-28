# ============================================================
# R/figure_theme_mnp.R
# Shared plotting theme, palettes, polymer labels, and export helpers
# for the PTB-MNP manuscript figure set.
# ============================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(scales)
  library(dplyr)
  library(grid)
})

# ----------------------------
# Canonical manuscript ordering / labels
# ----------------------------
mnp_polymer_order <- c(
  "Cumulative_MNPs",
  "PU", "PET", "PVC", "N66", "PC", "SBR",
  "PMMA", "N6", "PE", "PP", "ABS", "PS"
)

mnp_polymer_order_no_total <- setdiff(mnp_polymer_order, "Cumulative_MNPs")

mnp_polymer_label <- function(poly) {
  poly <- as.character(poly)
  vapply(poly, function(p) {
    if (p == "Cumulative_MNPs") return("Cumulative MNP burden")
    if (exists("POLYMER_FULL_NAMES", inherits = TRUE) && p %in% names(POLYMER_FULL_NAMES)) {
      return(paste0(POLYMER_FULL_NAMES[[p]], " (", p, ")"))
    }
    p
  }, character(1))
}

mnp_polymer_short_label <- function(poly) {
  poly <- as.character(poly)
  ifelse(poly == "Cumulative_MNPs", "Cumulative MNPs", poly)
}

mnp_factor_polymer_labels <- function(polymer_label, polymer_code = NULL, include_total = TRUE) {
  if (is.null(polymer_code)) {
    return(factor(polymer_label, levels = rev(unique(polymer_label))))
  }
  order_codes <- if (include_total) mnp_polymer_order else mnp_polymer_order_no_total
  ordered_labels <- mnp_polymer_label(order_codes[order_codes %in% polymer_code])
  factor(polymer_label, levels = rev(ordered_labels))
}

# ----------------------------
# Core colors
# ----------------------------
mnp_delivery_colors <- c(
  Term = "#8DB3D3",
  Preterm = "#F28E8E"
)

mnp_sample_type_colors <- c(
  `Procedural Control` = "#F8766D",
  `Procedural Controls` = "#F8766D",
  Control = "#F8766D",
  `Placental Sample` = "#00BFC4",
  `Placental Samples` = "#00BFC4",
  Placenta = "#00BFC4"
)

mnp_family_colors <- c(
  cumulative = "grey55",
  `Cumulative burden` = "grey55",
  polyesters = "#E69F00",
  `Oxygen-containing polymers` = "#E69F00",
  vinyls = "#CC79A7",
  `Halogenated polymers` = "#CC79A7",
  polyurethanes = "#D55E00",
  Polyurethane = "#D55E00",
  polyolefins = "#009E73",
  polycarbonates = "#7A3E9D",
  acrylics = "#999999",
  methacrylates = "#999999",
  polyamides = "#56B4E9",
  Polyamides = "#56B4E9",
  styrenics = "#F0E442",
  elastomers = "#F0E442",
  `Styrenic/rubber polymers` = "#F0E442",
  other = "#000000",
  Other = "#000000"
)

mnp_sensitivity_colors <- c(
  Crude = "#56B4E9",
  CRUDE = "#56B4E9",
  `Primary + HTN subtype` = "#E69F00",
  `Primary + race/ethnicity` = "#CC79A7",
  `Primary + SVI` = "#009E73",
  `Primary without SDI` = "#D55E00",
  `Primary adjusted, excluding CLR/PCA outliers` = "#7A3E9D",
  `Overall without GA` = "#56B4E9",
  `GA + hospital` = "#E69F00",
  `GA + SDI` = "#009E73",
  `GA + prior PTB` = "#CC79A7",
  `GA + hypertension` = "#7A3E9D",
  `Primary without SDI` = "#E69F00",
  `M1 + SVI` = "#009E73",
  `Primary + race/ethnicity` = "#CC79A7"
)

mnp_q_shapes <- c(`FALSE` = 16, `TRUE` = 17)

mnp_evidence_shapes <- c(
  `not nominal` = 16,
  `nominal p < 0.05 only` = 1,
  `BH q < 0.05, prespecified PVC/N66` = 23,
  `BH q < 0.05, all analytes` = 17,
  `global BH q < 0.05` = 15,
  `family BH q < 0.05` = 22,
  `BH q < 0.05` = 17
)

# ----------------------------
# Theme functions
# ----------------------------
theme_mnp <- function(base_size = 9.5, base_family = "") {
  ggplot2::theme_bw(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(linewidth = 0.25, color = "grey90"),
      plot.title = ggplot2::element_text(face = "bold", hjust = 0, size = base_size * 1.00),
      plot.subtitle = ggplot2::element_text(size = base_size * 0.82, hjust = 0),
      axis.title = ggplot2::element_text(size = base_size * 0.95),
      axis.text = ggplot2::element_text(size = base_size * 0.78),
      legend.title = ggplot2::element_text(face = "bold", size = base_size * 0.80),
      legend.text = ggplot2::element_text(size = base_size * 0.74),
      strip.background = ggplot2::element_rect(fill = "grey92", color = "grey60"),
      strip.text = ggplot2::element_text(size = base_size * 0.78, face = "bold"),
      plot.margin = ggplot2::margin(5, 6, 5, 5)
    )
}

theme_mnp_forest <- function(base_size = 11, base_family = "") {
  theme_mnp(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      axis.title.y = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(size = base_size * 0.78),
      legend.position = "right"
    )
}

theme_mnp_facet_forest <- function(base_size = 10.5, base_family = "") {
  theme_mnp_forest(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      legend.position = "bottom",
      legend.box = "vertical",
      legend.box.just = "left",
      legend.spacing.y = grid::unit(0.05, "cm"),
      axis.text.x = ggplot2::element_text(size = base_size * 0.70),
      axis.text.y = ggplot2::element_text(size = base_size * 0.70)
    )
}

theme_mnp_boxplot <- function(base_size = 9.5, base_family = "") {
  theme_mnp(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      legend.position = "none"
    )
}

# ----------------------------
# Scale helpers
# ----------------------------
scale_mnp_delivery_fill <- function(...) ggplot2::scale_fill_manual(values = mnp_delivery_colors, drop = FALSE, ...)
scale_mnp_delivery_color <- function(...) ggplot2::scale_color_manual(values = mnp_delivery_colors, drop = FALSE, ...)
scale_mnp_sample_type_fill <- function(...) ggplot2::scale_fill_manual(values = mnp_sample_type_colors, drop = FALSE, ...)
scale_mnp_sample_type_color <- function(...) ggplot2::scale_color_manual(values = mnp_sample_type_colors, drop = FALSE, ...)
scale_mnp_family_color <- function(...) ggplot2::scale_color_manual(values = mnp_family_colors, drop = FALSE, ...)
scale_mnp_family_fill <- function(...) ggplot2::scale_fill_manual(values = mnp_family_colors, drop = FALSE, ...)
scale_mnp_sensitivity_color <- function(...) ggplot2::scale_color_manual(values = mnp_sensitivity_colors, drop = FALSE, ...)
scale_mnp_q_shape <- function(name = "BH q < 0.05", ...) ggplot2::scale_shape_manual(values = mnp_q_shapes, name = name, ...)
scale_mnp_evidence_shape <- function(name = "Evidence", ...) ggplot2::scale_shape_manual(values = mnp_evidence_shapes, name = name, ...)

# ----------------------------
# Log-axis / p-value helpers
# ----------------------------
fmt_log_axis <- function(x) {
  dplyr::case_when(
    x < 0.01 ~ formatC(x, format = "f", digits = 3),
    x < 1 ~ formatC(x, format = "f", digits = 2),
    x < 10 ~ formatC(x, format = "f", digits = 1),
    TRUE ~ formatC(x, format = "f", digits = 0)
  )
}

common_log_breaks <- function(lims) {
  br <- c(0.001, 0.003, 0.01, 0.03, 0.05, 0.10, 0.20, 0.30, 0.50, 1, 2, 3, 4, 10, 20, 30, 50, 100)
  br[br >= lims[1] & br <= lims[2]]
}

scale_mnp_log10_or <- function(limits = c(0.03, 10), breaks = NULL, ...) {
  if (is.null(breaks)) breaks <- common_log_breaks(limits)
  ggplot2::scale_x_log10(breaks = breaks, limits = limits, labels = fmt_log_axis, ...)
}

fmt_mnp_p <- function(p) {
  ifelse(is.na(p), NA_character_, ifelse(p < 0.001, "<0.001", formatC(p, format = "f", digits = 3)))
}

# ----------------------------
# Patchwork panel tags and saving
# ----------------------------
mnp_panel_annotation <- function(tag_levels = "a", tag_size = 13) {
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("patchwork is required for mnp_panel_annotation().", call. = FALSE)
  }
  patchwork::plot_annotation(tag_levels = tag_levels) &
    ggplot2::theme(
      plot.tag = ggplot2::element_text(face = "bold", size = tag_size),
      plot.tag.position = c(0, 1)
    )
}

ggsave_mnp <- function(filename, plot, width, height, dpi = 300, units = "in", bg = "white", ...) {
  dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(filename = filename, plot = plot, width = width, height = height, dpi = dpi, units = units, bg = bg, ...)
  invisible(filename)
}
