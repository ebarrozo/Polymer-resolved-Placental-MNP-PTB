#!/usr/bin/env Rscript
# POLISHED VERSION: sources R/figure_theme_mnp.R and writes *_polished outputs where applicable.
# ============================================================
# 03_qc_figures_abs_mass.R
# QC figures and tables using absolute polymer mass (µg/sample)
# ============================================================
# Inputs from scripts/01_build_analysis_ready_inputs.R:
#   outputs/tables/mnp_all_records__qc_flags_and_outlier_flags_wide.csv
#   outputs/tables/mnp_controls__lob_lod_loq_maxblank_sdmethod_by_polymer.csv
#   outputs/tables/mnp_controls__blank_summary_stats_by_polymer.csv
#
# Outputs:
#   outputs/figs/Fig1A_cumulative_mnp_absug_controls_vs_placenta.png
#   outputs/figs/Fig1B_thresholds_facets_absug_controls_vs_placenta.png
#   outputs/figs/SuppFig_S1A_qcflag_distribution_controls_vs_placenta.png
#   outputs/tables/Fig1A_cumulative_mnp_absug_summary.csv
#   outputs/tables/Fig1B_retained_counts_by_polymer.csv
#   outputs/tables/SuppFig_S1A_qcflag_percent_by_polymer.csv
#   outputs/tables/Fig1_controls_outlier_exclusion_audit.csv
#   outputs/logs/03_qc_figures_abs_mass__run_log.txt
#
# Notes:
#   - PE blanks are expected to be zero if PE was not detected in controls.
#   - Control outlier columns are optional. If absent, no controls are excluded.
#   - Figures are PNG only.
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
  library(stringr)
  library(ggplot2)
  library(scales)
  library(glue)
})

source("00_config.R")
source(file.path("R", "utils_mnp.R"))
source(file.path("R", "figure_theme_mnp.R"))

# ----------------------------
# Inputs / outputs
# ----------------------------
in_wide <- file.path(DIR_TABLES, "mnp_all_records__qc_flags_and_outlier_flags_wide.csv")
in_thr  <- file.path(DIR_TABLES, "mnp_controls__lob_lod_loq_maxblank_sdmethod_by_polymer.csv")

stopifnot(file.exists(in_wide), file.exists(in_thr))

out_fig1a <- file.path(DIR_FIGS, "Fig1A_cumulative_mnp_absug_controls_vs_placenta_polished.png")
out_fig1b <- file.path(DIR_FIGS, "Fig1B_thresholds_facets_absug_controls_vs_placenta_polished.png")
out_s1a   <- file.path(DIR_FIGS, "SuppFig_S1A_qcflag_distribution_controls_vs_placenta_polished.png")

out_tbl_fig1a <- file.path(DIR_TABLES, "Fig1A_cumulative_mnp_absug_summary_polished.csv")
out_tbl_fig1b <- file.path(DIR_TABLES, "Fig1B_retained_counts_by_polymer_polished.csv")
out_tbl_s1a   <- file.path(DIR_TABLES, "SuppFig_S1A_qcflag_percent_by_polymer_polished.csv")
out_tbl_audit <- file.path(DIR_TABLES, "Fig1_controls_outlier_exclusion_audit_polished.csv")
out_log       <- file.path(DIR_LOGS, "03_qc_figures_abs_mass__run_log_polished.txt")

BASE_SIZE <- 14

# ----------------------------
# Helpers
# ----------------------------
format_p_nc <- function(p) {
  if (is.na(p) || p >= 0.05) return(NA_character_)
  if (p < 0.001) return("p < 0.001")
  paste0("p = ", formatC(p, format = "f", digits = 3))
}

poly_label <- function(poly) {
  full <- POLYMER_FULL_NAMES[[poly]]
  if (is.null(full) || is.na(full)) return(poly)
  paste0(full, " (", poly, ")")
}

threshold_col <- function(thr, candidates) {
  hit <- candidates[candidates %in% names(thr)]
  if (length(hit) == 0) return(rep(NA_real_, nrow(thr)))
  safe_num(thr[[hit[[1]]]])
}

# ----------------------------
# Load
# ----------------------------
df <- readr::read_csv(
  in_wide,
  show_col_types = FALSE,
  col_types = cols(
    SampleID = col_character(),
    .default = col_guess()
  )
)

thr <- readr::read_csv(in_thr, show_col_types = FALSE) %>%
  mutate(Analyte = toupper(as.character(Analyte))) %>%
  filter(Analyte %in% POLYMERS) %>%
  mutate(
    LOB_ug = threshold_col(cur_data(), c("LOB_ug", "lob_ug", "LOB")),
    LOD_ug = threshold_col(cur_data(), c("LOD_ug", "lod_ug", "LOD")),
    LOQ_ug = threshold_col(cur_data(), c("LOQ_ug", "loq_ug", "LOQ")),
    MaxBlank_ug = threshold_col(cur_data(), c("MaxBlank_ug", "max_blank_ug", "MaxBlank", "maxblank"))
  )

missing_poly <- setdiff(POLYMERS, names(df))
if (length(missing_poly) > 0) {
  stop("Input QC file is missing polymer mass columns: ", paste(missing_poly, collapse = ", "), call. = FALSE)
}

qc_cols <- paste0("QC_flag_", POLYMERS)
missing_qc <- setdiff(qc_cols, names(df))
if (length(missing_qc) > 0) {
  stop("Input QC file is missing QC flag columns: ", paste(missing_qc, collapse = ", "), call. = FALSE)
}

stopifnot(all(c("Dataset", "RecordID") %in% names(df)))

# ----------------------------
# Controls: optional consensus outlier exclusion
# ----------------------------
cons_cols <- paste0("Consensus_Outlier_TRUE_", POLYMERS)
cons_cols <- cons_cols[cons_cols %in% names(df)]

df <- df %>%
  mutate(
    Dataset = as.character(Dataset),
    is_control  = Dataset == "Control",
    is_placenta = Dataset == "Placenta",
    control_has_consensus_outlier =
      if (length(cons_cols) > 0) {
        is_control & (rowSums(across(all_of(cons_cols), ~ .x %in% TRUE), na.rm = TRUE) > 0)
      } else {
        FALSE
      }
  )

df_plotbase <- df %>%
  filter(Dataset %in% c("Control", "Placenta")) %>%
  filter(!(Dataset == "Control" & control_has_consensus_outlier))

audit <- tibble(
  n_controls_total = sum(df$Dataset == "Control", na.rm = TRUE),
  n_controls_excluded_consensus_outlier = sum(df$Dataset == "Control" & df$control_has_consensus_outlier, na.rm = TRUE),
  n_controls_used = sum(df_plotbase$Dataset == "Control", na.rm = TRUE),
  n_placenta_used = sum(df_plotbase$Dataset == "Placenta", na.rm = TRUE),
  pe_control_blank_forced_zero = FORCE_PE_CONTROL_BLANK_ZERO
)

write_csv_clean(audit, out_tbl_audit)

# ============================================================
# Fig 1A: cumulative absolute mass, controls vs placenta
# ============================================================
df_fig1a <- df_plotbase %>%
  mutate(
    Sample_Type = if_else(Dataset == "Control", "Procedural Control", "Placental Samples"),
    Sample_Type = factor(Sample_Type, levels = c("Procedural Control", "Placental Samples")),
    Cumulative_MNPs_absug = rowSums(across(all_of(POLYMERS), safe_num), na.rm = TRUE)
  )

w <- suppressWarnings(wilcox.test(Cumulative_MNPs_absug ~ Sample_Type, data = df_fig1a, exact = FALSE))
p_lab <- format_p_nc(w$p.value)

tbl_fig1a <- df_fig1a %>%
  group_by(Sample_Type) %>%
  summarise(
    n = n(),
    n_detected_gt0 = sum(is.finite(Cumulative_MNPs_absug) & Cumulative_MNPs_absug > 0),
    median_absug = median(Cumulative_MNPs_absug, na.rm = TRUE),
    q1_absug = quantile(Cumulative_MNPs_absug, 0.25, na.rm = TRUE, names = FALSE),
    q3_absug = quantile(Cumulative_MNPs_absug, 0.75, na.rm = TRUE, names = FALSE),
    mean_absug = mean(Cumulative_MNPs_absug, na.rm = TRUE),
    sd_absug = sd(Cumulative_MNPs_absug, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(wilcox_p = w$p.value)

write_csv_clean(tbl_fig1a, out_tbl_fig1a)

subtitle_1a <- glue(
  "Procedural controls: n = {sum(df_fig1a$Sample_Type == 'Procedural Control')}. ",
  "Placental samples: n = {sum(df_fig1a$Sample_Type == 'Placental Samples')}."
)

p1a <- ggplot(df_fig1a, aes(x = Sample_Type, y = Cumulative_MNPs_absug)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.20, width = 0.55) +
  geom_jitter(aes(color = Sample_Type), width = 0.18, alpha = 0.55, size = 2) +
  labs(
    title = "Cumulative MNP Mass by Sample Type",
    subtitle = subtitle_1a,
    y = expression("Cumulative MNP mass ("*mu*"g per sample)"),
    x = NULL
  ) +
  theme_classic(base_size = BASE_SIZE) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = BASE_SIZE),
    plot.title = element_text(face = "bold")
  )

if (!is.na(p_lab)) {
  y_pos <- max(df_fig1a$Cumulative_MNPs_absug, na.rm = TRUE) * 1.05
  p1a <- p1a + annotate("text", x = 1.5, y = y_pos, label = p_lab, size = 4.2)
}

ggsave_mnp(out_fig1a, p1a, width = 7.2, height = 5.4, dpi = 300)

# ============================================================
# Fig 1B: per-polymer mass + thresholds
# ============================================================
df_long <- df_plotbase %>%
  select(RecordID, SampleID, Dataset, all_of(POLYMERS), all_of(qc_cols)) %>%
  pivot_longer(
    cols = all_of(POLYMERS),
    names_to = "Analyte",
    values_to = "mass_ug"
  ) %>%
  mutate(
    Analyte = factor(Analyte, levels = POLYMERS),
    Sample_Type = if_else(Dataset == "Control", "Procedural Control", "Placental Samples"),
    Sample_Type = factor(Sample_Type, levels = c("Procedural Control", "Placental Samples")),
    mass_ug = safe_num(mass_ug)
  ) %>%
  left_join(thr, by = c("Analyte" = "Analyte"))

# Count retained/above LOD/above LOQ by polymer among placentas.
tbl_fig1b <- df_long %>%
  group_by(Analyte, Sample_Type) %>%
  summarise(
    n = n(),
    n_nonmissing = sum(!is.na(mass_ug)),
    n_detected_gt0 = sum(is.finite(mass_ug) & mass_ug > 0),
    n_above_LOD = sum(is.finite(mass_ug) & is.finite(LOD_ug) & mass_ug >= LOD_ug),
    n_above_LOQ = sum(is.finite(mass_ug) & is.finite(LOQ_ug) & mass_ug >= LOQ_ug),
    median_ug = median(mass_ug, na.rm = TRUE),
    q1_ug = quantile(mass_ug, 0.25, na.rm = TRUE, names = FALSE),
    q3_ug = quantile(mass_ug, 0.75, na.rm = TRUE, names = FALSE),
    .groups = "drop"
  ) %>%
  mutate(
    Analyte = as.character(Analyte),
    Polymer = vapply(Analyte, poly_label, character(1))
  )

write_csv_clean(tbl_fig1b, out_tbl_fig1b)

# Dynamic ordering by placental median mass, but keep PE first only if desired globally?
poly_order_fig <- tbl_fig1b %>%
  filter(Sample_Type == "Placental Samples") %>%
  arrange(desc(median_ug)) %>%
  pull(Analyte) %>%
  unique()

if (length(poly_order_fig) == 0) poly_order_fig <- POLYMERS

df_long <- df_long %>%
  mutate(
    Analyte = factor(as.character(Analyte), levels = poly_order_fig),
    Polymer = factor(vapply(as.character(Analyte), poly_label, character(1)),
                     levels = vapply(poly_order_fig, poly_label, character(1)))
  )

thr_plot <- thr %>%
  filter(Analyte %in% poly_order_fig) %>%
  mutate(
    Analyte = factor(Analyte, levels = poly_order_fig),
    Polymer = factor(vapply(as.character(Analyte), poly_label, character(1)),
                     levels = vapply(poly_order_fig, poly_label, character(1)))
  )

# ----------------------------
# Fig 1B retention labels
# Retained = placental samples >= LOD
# ----------------------------
fig1b_retention_labels <- df_long %>%
  dplyr::filter(Sample_Type == "Placental Samples") %>%
  dplyr::group_by(Analyte, Polymer) %>%
  dplyr::summarise(
    n_total = dplyr::n(),
    n_retained = sum(is.finite(mass_ug) & is.finite(LOD_ug) & mass_ug >= LOD_ug, na.rm = TRUE),
    pct_retained = 100 * n_retained / n_total,
    y_lab = max(mass_ug, LOD_ug, LOQ_ug, na.rm = TRUE),
    label = sprintf("Retained: %d/%d (%.1f%%)", n_retained, n_total, pct_retained),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    y_lab = ifelse(is.finite(y_lab), y_lab, 1),
    y_lab = y_lab * 1.08
  )

write_csv_clean(
  fig1b_retention_labels,
  file.path(DIR_TABLES, "Fig1B_percent_retained_by_polymer_polished.csv")
)

# Use pseudo-log to show zeros without dropping them.
p1b <- ggplot(df_long, aes(x = Sample_Type, y = mass_ug)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.18, width = 0.50) +
  geom_jitter(aes(shape = Sample_Type), width = 0.15, alpha = 0.45, size = 1.1) +
  geom_hline(
    data = thr_plot,
    aes(yintercept = LOD_ug),
    linetype = "dashed",
    linewidth = 0.35
  ) +
  geom_hline(
    data = thr_plot,
    aes(yintercept = LOQ_ug),
    linetype = "dotted",
    linewidth = 0.35
  ) +
  geom_text(
    data = fig1b_retention_labels,
    aes(x = 1.5, y = y_lab, label = label),
    inherit.aes = FALSE,
    size = 2.7,
    fontface = "plain"
  ) +
  facet_wrap(~ Polymer, scales = "free_y", ncol = 4) +
  scale_y_continuous(
    trans = scales::pseudo_log_trans(base = 10),
    labels = scales::label_number(),
    expand = expansion(mult = c(0.04, 0.20))
  ) +
  labs(
    title = "Polymer-Specific MNP Mass and Control-Derived Thresholds",
    subtitle = "Dashed line = LOD; dotted line = LOQ; retained = placental samples >= LOD",
    x = NULL,
    y = expression("Polymer mass ("*mu*"g per sample)")
  ) +
  theme_mnp(base_size = 11) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 30, hjust = 1),
    plot.title = element_text(face = "bold")
  )

ggsave_mnp(out_fig1b, p1b, width = 12, height = 8.8, dpi = 300)

# ============================================================
# Supp Fig S1A: QC flag distribution
# ============================================================
qc_long <- df_plotbase %>%
  select(RecordID, SampleID, Dataset, all_of(qc_cols)) %>%
  pivot_longer(
    cols = all_of(qc_cols),
    names_to = "Analyte",
    values_to = "QC_flag"
  ) %>%
  mutate(
    Analyte = str_remove(Analyte, "^QC_flag_"),
    Analyte = factor(Analyte, levels = POLYMERS),
    Sample_Type = if_else(Dataset == "Control", "Procedural Control", "Placental Samples"),
    Sample_Type = factor(Sample_Type, levels = c("Procedural Control", "Placental Samples")),
    QC_flag = as.character(QC_flag),
    QC_flag = if_else(is.na(QC_flag), "missing", QC_flag)
  )

qc_pct <- qc_long %>%
  count(Sample_Type, Analyte, QC_flag, name = "n") %>%
  group_by(Sample_Type, Analyte) %>%
  mutate(
    total = sum(n),
    pct = 100 * n / total,
    Polymer = vapply(as.character(Analyte), poly_label, character(1))
  ) %>%
  ungroup()

write_csv_clean(qc_pct, out_tbl_s1a)

qc_order <- c("above_LOQ", "above_LOD_below_LOQ", "above_LOB_below_LOD", "below_LOB", "not_detected", "missing")
qc_pct <- qc_pct %>%
  mutate(
    QC_flag = factor(QC_flag, levels = qc_order),
    Analyte = factor(as.character(Analyte), levels = POLYMERS),
    Polymer = factor(Polymer, levels = vapply(POLYMERS, poly_label, character(1)))
  )

p_s1a <- ggplot(qc_pct, aes(x = Polymer, y = pct, fill = QC_flag)) +
  geom_col(width = 0.82) +
  facet_wrap(~ Sample_Type, ncol = 1) +
  coord_flip() +
  labs(
    title = "QC Flag Distribution by Polymer",
    x = NULL,
    y = "Percent of records",
    fill = "QC flag"
  ) +
  theme_mnp(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold")
  )

ggsave_mnp(out_s1a, p_s1a, width = 8.5, height = 9.5, dpi = 300)

# ----------------------------
# Log
# ----------------------------
run_log <- c(
  "03_qc_figures_abs_mass.R",
  paste0("PROJECT_ROOT: ", PROJECT_ROOT),
  paste0("Input rows: ", nrow(df)),
  paste0("Plot rows after control outlier exclusion: ", nrow(df_plotbase)),
  paste0("Controls used: ", audit$n_controls_used),
  paste0("Placenta used: ", audit$n_placenta_used),
  paste0("PE control blank forced zero setting: ", FORCE_PE_CONTROL_BLANK_ZERO),
  "",
  "Wrote:",
  out_fig1a,
  out_fig1b,
  out_s1a,
  out_tbl_fig1a,
  out_tbl_fig1b,
  file.path(DIR_TABLES, "Fig1B_percent_retained_by_polymer_polished.csv"),
  out_tbl_s1a,
  out_tbl_audit
)

cat(paste(run_log, collapse = "\n"), file = out_log)
message(paste(run_log, collapse = "\n"))
