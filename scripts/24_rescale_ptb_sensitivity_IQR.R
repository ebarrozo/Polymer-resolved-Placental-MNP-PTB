#!/usr/bin/env Rscript
# ============================================================
# 24_rescale_ptb_sensitivity_IQR.R
#
# Propagate reviewer Point 1 across EXISTING PTB sensitivity outputs
# that are all based on the same underlying PTB cohort/model exposure
# definition. This script:
#   - applies the polymer-specific IQR scale from script 22
#   - rewrites PTB model and leave-one-out tables with IQR effects
#   - produces two clean candidate supplementary panels:
#       A) model-specification forest
#       B) leave-one-out range plot
#
# It does NOT touch:
#   - binary threshold models (>=LOD / >LOQ), because those ORs are
#     already correctly interpreted and require no IQR rescaling
#   - ordinal gestational-age models
#   - labor augmentation models
#   - pathology models
# Those require analysis-population-specific IQRs and will be handled
# after the global code audit identifies their exact source scripts.
# ============================================================

source("00_config.R")

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
  library(ggplot2)
  library(stringr)
})

find_latest_dir <- function(parent, prefix) {
  dirs <- list.dirs(parent, recursive = FALSE, full.names = TRUE)
  dirs <- dirs[grepl(paste0("^", prefix), basename(dirs))]
  if (length(dirs) == 0) stop("No directory beginning with ", prefix, call. = FALSE)
  dirs[order(file.info(dirs)$mtime, decreasing = TRUE)][1]
}

AUDIT_DIR <- find_latest_dir(DIR_TABLES, "22_ptb_effect_scale_audit__")
MODEL_DIR <- find_latest_dir(DIR_TABLES, "06_primary_ptb_models__")

iqr_path <- file.path(AUDIT_DIR, "PTB_primary_IQR_reference_values.csv")
all_path <- file.path(MODEL_DIR, "PTB__ALL_models_long.csv")

stopifnot(file.exists(iqr_path), file.exists(all_path))

iqr_ref <- read_csv(iqr_path, show_col_types = FALSE)
all_models <- read_csv(all_path, show_col_types = FALSE)

OUT_DIR <- file.path(DIR_TABLES, "24_ptb_sensitivity_IQR")
OUT_FIG_DIR <- file.path(DIR_FIGS, "24_ptb_sensitivity_IQR")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_FIG_DIR, recursive = TRUE, showWarnings = FALSE)

delta_old <- log10(2)

rescale_from_old <- function(df) {
  df %>%
    left_join(iqr_ref %>% select(Polymer, IQR_log10_1px), by = "Polymer") %>%
    mutate(
      beta_recovered = log(OR_2x) / delta_old,
      beta_lo_recovered = log(OR_2x_lo) / delta_old,
      beta_hi_recovered = log(OR_2x_hi) / delta_old,
      se_beta_recovered = (beta_hi_recovered - beta_lo_recovered) / (2 * 1.96),
      logOR_IQR = beta_recovered * IQR_log10_1px,
      SE_logOR_IQR = se_beta_recovered * IQR_log10_1px,
      OR_IQR = exp(logOR_IQR),
      OR_IQR_lo = exp((beta_recovered - 1.96 * se_beta_recovered) * IQR_log10_1px),
      OR_IQR_hi = exp((beta_recovered + 1.96 * se_beta_recovered) * IQR_log10_1px)
    )
}

all_iqr <- rescale_from_old(all_models)
write_csv(all_iqr, file.path(OUT_DIR, "PTB_ALL_models_long_IQR.csv"))

# ------------------------------------------------------------
# Search for leave-one-out table in current model directory
# ------------------------------------------------------------

csvs <- list.files(MODEL_DIR, pattern = "\\.csv$", full.names = TRUE)
loo_candidates <- csvs[grepl("leave|loo", basename(csvs), ignore.case = TRUE)]

loo_out <- NULL

for (f in loo_candidates) {
  tmp <- tryCatch(read_csv(f, show_col_types = FALSE), error = function(e) NULL)
  if (is.null(tmp)) next
  if (all(c("Polymer", "OR_2x", "OR_2x_lo", "OR_2x_hi") %in% names(tmp))) {
    tmp$source_file <- basename(f)
    loo_out <- bind_rows(loo_out, tmp)
  }
}

if (!is.null(loo_out) && nrow(loo_out) > 0) {
  loo_iqr <- rescale_from_old(loo_out)
  write_csv(loo_iqr, file.path(OUT_DIR, "PTB_leave_one_out_IQR.csv"))
} else {
  loo_iqr <- NULL
  warning("No leave-one-out table with OR_2x columns found in ", MODEL_DIR)
}

# ------------------------------------------------------------
# Candidate SFig4 panel A: all model specifications
# ------------------------------------------------------------

model_order <- all_iqr %>%
  distinct(ModelTag, ModelLabel) %>%
  mutate(ModelLabel = if_else(is.na(ModelLabel) | ModelLabel == "", ModelTag, ModelLabel))

key_polys <- c("PU", "PET", "SBR", "ABS", "Cumulative_MNPs")

plot_a <- all_iqr %>%
  filter(Polymer %in% key_polys) %>%
  mutate(
    PolymerLabel = if_else(Polymer == "Cumulative_MNPs", "Total MNP burden", PolymerLabel)
  )

if (nrow(plot_a) > 0) {
  pA <- ggplot(plot_a, aes(x = OR_IQR, y = ModelLabel)) +
    geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.4) +
    geom_errorbarh(aes(xmin = OR_IQR_lo, xmax = OR_IQR_hi), height = 0.15) +
    geom_point(size = 2.1) +
    facet_wrap(~ PolymerLabel, scales = "free_y") +
    scale_x_log10() +
    labs(
      title = "PTB sensitivity models on the corrected IQR effect scale",
      x = "OR per polymer-specific IQR increase in log10(1 + burden)",
      y = NULL
    ) +
    theme_bw(base_size = 10) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold")
    )

  ggsave(
    file.path(OUT_FIG_DIR, "SFig4_candidate_A_model_specifications_IQR.png"),
    pA, width = 10, height = 7, dpi = 300
  )
  ggsave(
    file.path(OUT_FIG_DIR, "SFig4_candidate_A_model_specifications_IQR.pdf"),
    pA, width = 10, height = 7
  )
}

# ------------------------------------------------------------
# Candidate SFig4 panel B: leave-one-out estimate range
# ------------------------------------------------------------

if (!is.null(loo_iqr) && nrow(loo_iqr) > 0) {
  # Identify deleted-sample rows if column exists; otherwise summarize all.
  id_col <- intersect(c("leave_out_sample", "left_out", "excluded_sample", "PregnancyID"), names(loo_iqr))
  plot_b <- loo_iqr %>%
    filter(Polymer %in% c("PU", "PET", "SBR", "ABS")) %>%
    group_by(Polymer, PolymerLabel) %>%
    summarise(
      OR_min = min(OR_IQR, na.rm = TRUE),
      OR_max = max(OR_IQR, na.rm = TRUE),
      OR_median = median(OR_IQR, na.rm = TRUE),
      .groups = "drop"
    )

  pB <- ggplot(plot_b, aes(x = OR_median, y = reorder(PolymerLabel, OR_median))) +
    geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.4) +
    geom_segment(aes(x = OR_min, xend = OR_max, yend = reorder(PolymerLabel, OR_median)),
                 linewidth = 0.9) +
    geom_point(size = 2.5) +
    scale_x_log10() +
    labs(
      title = "Leave-one-out robustness on the corrected IQR effect scale",
      subtitle = "Horizontal ranges show the minimum and maximum one-participant-deleted OR estimates",
      x = "OR per polymer-specific IQR increase in log10(1 + burden)",
      y = NULL
    ) +
    theme_bw(base_size = 10) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      plot.title = element_text(face = "bold")
    )

  ggsave(
    file.path(OUT_FIG_DIR, "SFig4_candidate_B_leave_one_out_IQR.png"),
    pB, width = 7.2, height = 4.5, dpi = 300
  )
  ggsave(
    file.path(OUT_FIG_DIR, "SFig4_candidate_B_leave_one_out_IQR.pdf"),
    pB, width = 7.2, height = 4.5
  )
}

legend_note <- paste0(
  "Continuous PTB burden models should now be reported as odds ratios per polymer-specific ",
  "interquartile-range increase in log10(1 + placental polymer burden [µg/g]). ",
  "Binary analytical-threshold models (>=LOD or >LOQ) are unaffected by reviewer Point 1 ",
  "and should retain their original odds-ratio interpretation."
)
writeLines(legend_note, file.path(OUT_DIR, "supplementary_reporting_note_point1.txt"))

message("PTB sensitivity rescaling complete.")
message("Tables: ", OUT_DIR)
message("Figures: ", OUT_FIG_DIR)
