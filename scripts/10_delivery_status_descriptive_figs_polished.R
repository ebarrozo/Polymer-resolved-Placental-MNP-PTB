#!/usr/bin/env Rscript
# POLISHED VERSION: sources R/figure_theme_mnp.R and writes *_polished outputs where applicable.
# ============================================================
# 10_delivery_status_descriptive_figs.R
# Delivery-status descriptive figures for MNP burden
# ============================================================
#
# Run after:
#   source("run_analysis_setup.R")
#
# Input:
#   outputs/tables/mnp_ptb_cohort__analysis_ready_model_table.csv
#
# Outputs:
#   outputs/figs/fig_delivery_status__cumulative_mnp_boxplot.png
#   outputs/figs/fig_delivery_status__polymer_effect_size_median_diff.png
#   outputs/figs/fig_delivery_status__combined_cumulative_and_effect_size.png
#
# Tables:
#   outputs/tables/10_delivery_status_descriptives__<timestamp>/
#     delivery_status_cumulative_summary.csv
#     delivery_status_polymer_effect_sizes.csv
#
# Figure logic:
#   Panel A:
#     Cumulative MNP burden by delivery class.
#     Boxplot + points + Wilcoxon p-value.
#
#   Panel B:
#     Polymer-specific effect size:
#       Preterm median log10(burden + pseudocount) -
#       Term median log10(burden + pseudocount)
#     95% CI from stratified bootstrap.
#     Point color = polymer family.
#     Point shape = BH q < 0.05 from Wilcoxon test across analytes.
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
  library(stringr)
  library(ggplot2)
  library(scales)
  library(patchwork)
})

source("00_config.R")
source(file.path("R", "utils_mnp.R"))
source(file.path("R", "figure_theme_mnp.R"))

# ----------------------------
# Settings
# ----------------------------
IN_MODEL_TABLE <- file.path(DIR_TABLES, "mnp_ptb_cohort__analysis_ready_model_table.csv")

STAMP <- format(Sys.time(), "%Y%m%d_%H%M%S")
OUT_DIR <- file.path(DIR_TABLES, paste0("10_delivery_status_descriptives__", STAMP))
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

OUT_CUM_SUMMARY <- file.path(OUT_DIR, "delivery_status_cumulative_summary_polished.csv")
OUT_EFFECTS <- file.path(OUT_DIR, "delivery_status_polymer_effect_sizes_polished.csv")
OUT_LOG <- file.path(DIR_LOGS, "10_delivery_status_descriptive_figs__run_log_polished.txt")

OUT_FIG_BOX <- file.path(DIR_FIGS, "fig_delivery_status__cumulative_mnp_boxplot_polished.png")
OUT_FIG_EFFECT <- file.path(DIR_FIGS, "fig_delivery_status__polymer_effect_size_median_diff_polished.png")
OUT_FIG_COMBINED <- file.path(DIR_FIGS, "fig_delivery_status__combined_cumulative_and_effect_size_polished.png")

BASE_SIZE <- 10.5
N_BOOT <- 5000
BOOT_SEED <- 20260519

# Consistent polymer order, matching prior figures.
POLYMER_ORDER_DELIVERY <- c(
  "PET", "PVC", "PU", "PP", "PS", "PC",
  "PE", "Cumulative_MNPs", "PMMA", "N66", "N6", "SBR", "ABS"
)

# Use available analytes only.
POLYMER_ORDER_DELIVERY <- POLYMER_ORDER_DELIVERY[
  POLYMER_ORDER_DELIVERY %in% c("Cumulative_MNPs", POLYMERS)
]

# ----------------------------
# Helpers
# ----------------------------
safe_num2 <- function(x) suppressWarnings(as.numeric(as.character(x)))

poly_family <- function(poly) {
  dplyr::case_when(
    poly == "Cumulative_MNPs" ~ "cumulative",
    poly %in% c("PE", "PP") ~ "polyolefins",
    poly %in% c("PET") ~ "polyesters",
    poly %in% c("PVC") ~ "vinyls",
    poly %in% c("PU") ~ "polyurethanes",
    poly %in% c("PC") ~ "polycarbonates",
    poly %in% c("PMMA") ~ "acrylics",
    poly %in% c("N6", "N66") ~ "polyamides",
    poly %in% c("SBR") ~ "styrenics",
    poly %in% c("ABS", "PS") ~ "styrenics",
    TRUE ~ "other"
  )
}

poly_full <- function(poly) {
  if (poly == "Cumulative_MNPs") return("Cumulative MNPs")
  if (poly %in% names(POLYMER_FULL_NAMES)) return(POLYMER_FULL_NAMES[[poly]])
  poly
}

poly_plot_label <- function(poly) {
  if (poly == "Cumulative_MNPs") return("Cumulative MNPs")
  paste0(poly_full(poly), " (", poly, ")")
}

get_ug_col <- function(poly) {
  if (poly == "Cumulative_MNPs") "Cumulative_MNPs_ug_per_g" else legacy_col_ug_g(poly)
}

get_log_col <- function(poly) {
  if (poly == "Cumulative_MNPs") "Cumulative_MNPs_log10ugpg" else legacy_col_log10(poly)
}

read_model_table <- function(path) {
  readr::read_csv(
    path,
    show_col_types = FALSE,
    guess_max = 10000,
    col_types = cols(
      SampleID = col_character(),
      PregnancyID = col_character(),
      `Pregnancy ID` = col_character(),
      home_zip_2nd_3rd_trim = col_character(),
      .default = col_guess()
    )
  ) %>%
    mutate(
      SampleID = as.character(SampleID),
      PregnancyID = as.character(PregnancyID),
      Delivery_class = factor(as.character(Delivery_class), levels = c("Term", "Preterm"))
    )
}

wilcox_p <- function(x, g) {
  d <- tibble(x = x, g = g) %>%
    filter(is.finite(x), !is.na(g))

  if (n_distinct(d$g) != 2) return(NA_real_)
  if (any(table(d$g) < 2)) return(NA_real_)
  if (length(unique(d$x)) < 2) return(NA_real_)

  suppressWarnings(
    tryCatch(
      wilcox.test(x ~ g, data = d, exact = FALSE)$p.value,
      error = function(e) NA_real_
    )
  )
}

median_diff <- function(x, g) {
  d <- tibble(x = x, g = g) %>%
    filter(is.finite(x), !is.na(g))

  if (!all(c("Term", "Preterm") %in% unique(as.character(d$g)))) return(NA_real_)

  med_preterm <- median(d$x[d$g == "Preterm"], na.rm = TRUE)
  med_term <- median(d$x[d$g == "Term"], na.rm = TRUE)
  med_preterm - med_term
}

bootstrap_median_diff <- function(x, g, n_boot = N_BOOT, seed = BOOT_SEED) {
  d <- tibble(x = x, g = g) %>%
    filter(is.finite(x), !is.na(g))

  term <- d$x[d$g == "Term"]
  preterm <- d$x[d$g == "Preterm"]

  if (length(term) < 3 || length(preterm) < 3) {
    return(tibble(diff = NA_real_, lo = NA_real_, hi = NA_real_))
  }

  set.seed(seed)
  boot <- replicate(n_boot, {
    median(sample(preterm, length(preterm), replace = TRUE), na.rm = TRUE) -
      median(sample(term, length(term), replace = TRUE), na.rm = TRUE)
  })

  tibble(
    diff = median_diff(x, g),
    lo = unname(quantile(boot, 0.025, na.rm = TRUE)),
    hi = unname(quantile(boot, 0.975, na.rm = TRUE))
  )
}

# ----------------------------
# Load data
# ----------------------------
df0 <- read_model_table(IN_MODEL_TABLE)

df <- df0 %>%
  filter(Delivery_class %in% c("Term", "Preterm")) %>%
  mutate(
    Delivery_class = droplevels(Delivery_class),
    Cumulative_MNPs_ug_per_g = safe_num2(Cumulative_MNPs_ug_per_g),
    Cumulative_MNPs_log10ugpg = safe_num2(Cumulative_MNPs_log10ugpg)
  )

# ----------------------------
# Cumulative MNP burden boxplot summary
# ----------------------------
cum_summary <- df %>%
  group_by(Delivery_class) %>%
  summarise(
    n = sum(is.finite(Cumulative_MNPs_ug_per_g)),
    mean = mean(Cumulative_MNPs_ug_per_g, na.rm = TRUE),
    sd = sd(Cumulative_MNPs_ug_per_g, na.rm = TRUE),
    median = median(Cumulative_MNPs_ug_per_g, na.rm = TRUE),
    q1 = quantile(Cumulative_MNPs_ug_per_g, 0.25, na.rm = TRUE),
    q3 = quantile(Cumulative_MNPs_ug_per_g, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

cum_p <- wilcox_p(df$Cumulative_MNPs_ug_per_g, df$Delivery_class)
cum_log_p <- wilcox_p(df$Cumulative_MNPs_log10ugpg, df$Delivery_class)

cum_summary <- cum_summary %>%
  mutate(
    wilcox_p_raw_ug_per_g = cum_p,
    wilcox_p_log10 = cum_log_p
  )

write_csv_clean(cum_summary, OUT_CUM_SUMMARY)

# ----------------------------
# Effect-size table
# ----------------------------
effect_tbl <- purrr::map_dfr(POLYMER_ORDER_DELIVERY, function(poly) {
  ug_col <- get_ug_col(poly)
  log_col <- get_log_col(poly)

  if (!ug_col %in% names(df)) {
    return(tibble(
      Polymer = poly,
      PolymerLabel = poly_plot_label(poly),
      family = poly_family(poly),
      n_term = NA_integer_,
      n_preterm = NA_integer_,
      median_term = NA_real_,
      median_preterm = NA_real_,
      median_diff_preterm_minus_term = NA_real_,
      ci_lo = NA_real_,
      ci_hi = NA_real_,
      p = NA_real_,
      skip_reason = "missing_ug_column"
    ))
  }

  x_ug <- safe_num2(df[[ug_col]])
  # Build log10(burden + c) directly from ug/g.
  # c is half of the smallest positive value for this analyte.
  pos_min <- suppressWarnings(min(x_ug[is.finite(x_ug) & x_ug > 0], na.rm = TRUE))
  pseudocount <- ifelse(is.finite(pos_min), pos_min / 2, 1e-6)
  x_log <- log10(pmax(x_ug, 0) + pseudocount)

  b <- bootstrap_median_diff(x_log, df$Delivery_class, n_boot = N_BOOT, seed = BOOT_SEED + match(poly, POLYMER_ORDER_DELIVERY))
  p <- wilcox_p(x_log, df$Delivery_class)

  tibble(
    Polymer = poly,
    PolymerLabel = poly_plot_label(poly),
    family = poly_family(poly),
    ug_col = ug_col,
    pseudocount_ug_per_g = pseudocount,
    n_term = sum(is.finite(x_log) & df$Delivery_class == "Term"),
    n_preterm = sum(is.finite(x_log) & df$Delivery_class == "Preterm"),
    median_term = median(x_log[df$Delivery_class == "Term"], na.rm = TRUE),
    median_preterm = median(x_log[df$Delivery_class == "Preterm"], na.rm = TRUE),
    median_diff_preterm_minus_term = b$diff,
    ci_lo = b$lo,
    ci_hi = b$hi,
    p = p,
    skip_reason = NA_character_
  )
}) %>%
  mutate(
    q_BH = p.adjust(p, method = "BH"),
    sig_q05 = q_BH < 0.05,
    family = factor(family, levels = c(
      "cumulative", "polyesters", "vinyls", "polyurethanes", "polyolefins",
      "polycarbonates", "acrylics", "polyamides", "styrenics", "other"
    )),
    PolymerLabel = factor(PolymerLabel, levels = rev(vapply(POLYMER_ORDER_DELIVERY, poly_plot_label, character(1))))
  )

write_csv_clean(effect_tbl, OUT_EFFECTS)

# ----------------------------
# Figure A: cumulative boxplot
# ----------------------------
p_label <- ifelse(is.na(cum_p), "p = NA", paste0("p = ", scales::pvalue(cum_p, accuracy = 0.001)))
y_max <- max(df$Cumulative_MNPs_ug_per_g, na.rm = TRUE)

p_box <- ggplot(df, aes(x = Delivery_class, y = Cumulative_MNPs_ug_per_g, fill = Delivery_class)) +
  geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.65, linewidth = 0.45) +
  geom_jitter(width = 0.12, height = 0, size = 1.1, alpha = 0.75, color = "black") +
  annotate("segment", x = 1, xend = 2, y = y_max * 1.08, yend = y_max * 1.08, linewidth = 0.35) +
  annotate("text", x = 1.5, y = y_max * 1.13, label = p_label, size = 3.4) +
  scale_fill_manual(values = c("Term" = "#8DB3D3", "Preterm" = "#F28E8E"), guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.18))) +
  labs(
    title = "Cumulative MNP burden by delivery status",
    x = NULL,
    y = "Cumulative MNP mass (\u00b5g/g tissue)"
  ) +
  theme_mnp(base_size = BASE_SIZE) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text.x = element_text(size = 9)
  )

ggsave_mnp(OUT_FIG_BOX, p_box, width = 4.1, height = 3.4, dpi = 300)

# ----------------------------
# Figure B: polymer-specific median differences
# ----------------------------
p_eff <- ggplot(effect_tbl, aes(y = PolymerLabel, x = median_diff_preterm_minus_term, color = family)) +
  geom_vline(xintercept = 0, linewidth = 0.4, color = "grey35") +
  geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi), height = 0.22, linewidth = 0.65, alpha = 0.95) +
  geom_point(aes(shape = sig_q05), size = 2.7, stroke = 0.85) +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17), name = "BH q < 0.05") +
  scale_color_manual(
    values = c(
      "cumulative" = "grey55",
      "polyesters" = "#E69F00",
      "vinyls" = "#CC79A7",
      "polyurethanes" = "#D55E00",
      "polyolefins" = "#009E73",
      "polycarbonates" = "#7A3E9D",
      "acrylics" = "#999999",
      "polyamides" = "#56B4E9",
      "styrenics" = "#F0E442",
      "other" = "#000000"
    ),
    drop = FALSE
  ) +
  labs(
    title = "Effect size of polymer-specific MNP burden by delivery status",
    subtitle = "Effect size = Preterm median \u2212 Term median on log10(burden + pseudocount) scale",
    x = "Preterm \u2212 Term median log10(burden + c)",
    y = NULL,
    color = "Polymer family"
  ) +
  theme_mnp(base_size = BASE_SIZE) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 8.8, hjust = 0.5),
    axis.text.y = element_text(size = 7.5),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 8.5),
    legend.text = element_text(size = 7.5)
  )

ggsave_mnp(OUT_FIG_EFFECT, p_eff, width = 8.4, height = 4.6, dpi = 300)

# ----------------------------
# Combined figure
# ----------------------------
p_combined <- p_box / p_eff +
  patchwork::plot_layout(heights = c(0.95, 1.35)) +
  patchwork::plot_annotation(tag_levels = "A")

ggsave_mnp(OUT_FIG_COMBINED, p_combined, width = 8.7, height = 7.6, dpi = 300)

# ----------------------------
# Log
# ----------------------------
run_log <- c(
  "10_delivery_status_descriptive_figs.R",
  paste0("PROJECT_ROOT: ", PROJECT_ROOT),
  paste0("Rows used: ", nrow(df)),
  "",
  "Delivery counts:",
  capture.output(print(df %>% count(Delivery_class), n = Inf)),
  "",
  "Cumulative summary:",
  capture.output(print(cum_summary, n = Inf)),
  "",
  "Effect table:",
  capture.output(print(effect_tbl, n = Inf)),
  "",
  "Wrote:",
  OUT_CUM_SUMMARY,
  OUT_EFFECTS,
  OUT_FIG_BOX,
  OUT_FIG_EFFECT,
  OUT_FIG_COMBINED
)

cat(paste(run_log, collapse = "\n"), file = OUT_LOG)
message(paste(run_log, collapse = "\n"))
