#!/usr/bin/env Rscript
# ============================================================
# 22_ptb_effect_scale_audit.R
#
# Purpose
#   Address reviewer Point 1: the existing PTB models use
#     log10(1 + burden_ug_per_g)
#   and the current OR_2x values were obtained by multiplying the
#   coefficient by log10(2). That quantity is NOT an exact OR for
#   x -> 2x burden; it is the OR for a log10(2) increment in the
#   transformed exposure, equivalently a 2-fold increase in (1 + x).
#
#   This script does NOT change/refit the existing primary PTB models.
#   It:
#     1) recovers beta and SE(beta) from the current OR_2x outputs;
#     2) rescales the same fitted coefficients to an exact OR per IQR
#        increase in log10(1 + ug/g);
#     3) computes exact x -> 2x contrasts at representative raw
#        burdens (P25, median, P75 among positive measurements);
#     4) generates a candidate Figure 3 using the IQR effect scale;
#     5) exports an audit table comparing all scales.
#
# Run from PROJECT_ROOT:
#   source("scripts/22_ptb_effect_scale_audit.R")
#
# PROJECT_ROOT (from 00_config.R):
#   /Users/enricobarrozo/Library/CloudStorage/Box-Box/BarrozoLab/
#   1_Mentoring and Publications/Manuscripts/MNP_NatureComms/
#   data_RawfromMarcus/
#
# This script intentionally does not overwrite any existing figures,
# tables, or model outputs.
# ============================================================

source("00_config.R")

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
  library(ggplot2)
  library(stringr)
  library(scales)
})

# ------------------------------------------------------------
# Settings
# ------------------------------------------------------------

PRIMARY_MODEL_TAG <- "M1_ADJ_PRIMARY"

DISPLAY_SENS_MODELS <- c(
  "M0_CRUDE",
  "M1_ADJ_PRIMARY_NO_CLR_OUTLIERS",
  "M1_PLUS_HTNsubtype"
)

AUDIT_DIR <- file.path(
  DIR_TABLES,
  paste0("22_ptb_effect_scale_audit__", format(Sys.time(), "%Y%m%d_%H%M%S"))
)
dir.create(AUDIT_DIR, recursive = TRUE, showWarnings = FALSE)

AUDIT_FIG_DIR <- file.path(DIR_FIGS, "22_ptb_effect_scale_audit")
dir.create(AUDIT_FIG_DIR, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------

find_latest_dir <- function(parent, prefix) {
  dirs <- list.dirs(parent, recursive = FALSE, full.names = TRUE)
  dirs <- dirs[grepl(paste0("^", prefix), basename(dirs))]
  if (length(dirs) == 0) {
    stop("No directory beginning with '", prefix, "' found under ", parent, call. = FALSE)
  }
  dirs[order(file.info(dirs)$mtime, decreasing = TRUE)][1]
}

first_existing <- function(paths) {
  hit <- paths[file.exists(paths)]
  if (length(hit) == 0) return(NA_character_)
  hit[[1]]
}

safe_quantile <- function(x, p) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  as.numeric(stats::quantile(x, probs = p, na.rm = TRUE, names = FALSE, type = 7))
}

# Exact contrast for raw x -> 2x under log10(1+x) model
true_double_delta <- function(x) {
  ifelse(
    is.na(x) | x < 0,
    NA_real_,
    log10(1 + 2 * x) - log10(1 + x)
  )
}

or_from_beta_delta <- function(beta, se_beta, delta) {
  tibble(
    OR = exp(beta * delta),
    OR_lo = exp((beta - 1.96 * se_beta) * delta),
    OR_hi = exp((beta + 1.96 * se_beta) * delta)
  )
}

# ------------------------------------------------------------
# Locate latest primary-model output
# ------------------------------------------------------------

MODEL_DIR <- find_latest_dir(DIR_TABLES, "06_primary_ptb_models__")
message("Using primary-model directory: ", MODEL_DIR)

all_models_path <- file.path(MODEL_DIR, "PTB__ALL_models_long.csv")
primary_path    <- file.path(MODEL_DIR, "table_primary_ptb_models.csv")

stopifnot(file.exists(all_models_path), file.exists(primary_path))

all_models <- readr::read_csv(all_models_path, show_col_types = FALSE)
primary_tbl <- readr::read_csv(primary_path, show_col_types = FALSE)

required_model_cols <- c(
  "ModelTag", "ModelLabel", "Polymer", "PolymerLabel", "ExposureCol",
  "N", "Events", "NonEvents", "Covariates",
  "OR_2x", "OR_2x_lo", "OR_2x_hi", "q_BH"
)
missing_model_cols <- setdiff(required_model_cols, names(all_models))
if (length(missing_model_cols) > 0) {
  stop(
    "Missing required columns in PTB__ALL_models_long.csv: ",
    paste(missing_model_cols, collapse = ", "),
    call. = FALSE
  )
}

# ------------------------------------------------------------
# Locate canonical analysis-ready subject-level model table
# ------------------------------------------------------------

candidate_model_tables <- c(
  file.path(DIR_TABLES, "mnp_ptb_cohort__analysis_ready_model_table.csv"),
  file.path(DIR_MODEL_READY, "placenta_model_ready_wide.csv"),
  file.path("outputs", "tables", "mnp_ptb_cohort__analysis_ready_model_table.csv")
)

MODEL_TABLE_PATH <- first_existing(candidate_model_tables)

if (is.na(MODEL_TABLE_PATH)) {
  stop(
    "Could not find an analysis-ready model table.\nTried:\n  ",
    paste(candidate_model_tables, collapse = "\n  "),
    call. = FALSE
  )
}

message("Using analysis-ready model table: ", MODEL_TABLE_PATH)
dat <- readr::read_csv(MODEL_TABLE_PATH, show_col_types = FALSE, guess_max = 10000)

if (!"PTB01" %in% names(dat)) {
  stop("PTB01 not found in model-ready table.", call. = FALSE)
}

# ------------------------------------------------------------
# Recover beta and SE(beta) from existing reported scale
# ------------------------------------------------------------
# Existing:
#   logOR_2x = beta * log10(2)
# therefore:
#   beta = log(OR_2x) / log10(2)
#
# This exactly recovers the fitted exposure coefficient as long as
# OR_2x was generated from exp(beta * log10(2)), which is the current
# pipeline definition.

delta_old <- log10(2)

models_recovered <- all_models %>%
  mutate(
    beta_recovered = log(OR_2x) / delta_old,
    beta_lo_recovered = log(OR_2x_lo) / delta_old,
    beta_hi_recovered = log(OR_2x_hi) / delta_old,
    se_beta_recovered = (beta_hi_recovered - beta_lo_recovered) / (2 * 1.96)
  )

# ------------------------------------------------------------
# Build polymer/exposure/raw-column map
# ------------------------------------------------------------

primary_rows <- models_recovered %>%
  filter(ModelTag == PRIMARY_MODEL_TAG)

if (nrow(primary_rows) == 0) {
  stop("No rows found for ModelTag == ", PRIMARY_MODEL_TAG, call. = FALSE)
}

# Prefer explicit raw ug/g names. Infer them from ExposureCol if needed.
exposure_map <- primary_rows %>%
  distinct(Polymer, PolymerLabel, ExposureCol, Covariates) %>%
  mutate(
    UgCol = case_when(
      Polymer == "Cumulative_MNPs" ~ "Cumulative_MNPs_ug_per_g",
      TRUE ~ paste0(Polymer, "_ug_per_g")
    )
  )

missing_exposure_cols <- setdiff(exposure_map$ExposureCol, names(dat))
if (length(missing_exposure_cols) > 0) {
  stop(
    "Model-ready table is missing transformed exposure columns: ",
    paste(missing_exposure_cols, collapse = ", "),
    call. = FALSE
  )
}

missing_raw_cols <- setdiff(exposure_map$UgCol, names(dat))
if (length(missing_raw_cols) > 0) {
  stop(
    "Model-ready table is missing raw ug/g columns required for true doubling audit: ",
    paste(missing_raw_cols, collapse = ", "),
    call. = FALSE
  )
}

# ------------------------------------------------------------
# Compute PRIMARY-model analysis-sample IQR for each exposure
# ------------------------------------------------------------
# We use the primary adjusted model's complete-case analysis sample
# separately for each polymer. This makes the displayed IQR scale
# correspond to the actual primary model sample. The same polymer-
# specific IQR is then applied to sensitivity-model coefficients so
# all model specifications are displayed on one common exposure unit.

parse_covariates <- function(x) {
  if (is.na(x) || !nzchar(x)) return(character(0))
  trimws(unlist(strsplit(x, "\\+")))
}

iqr_rows <- vector("list", nrow(exposure_map))

for (i in seq_len(nrow(exposure_map))) {
  row <- exposure_map[i, ]
  covs <- parse_covariates(row$Covariates)
  needed <- unique(c("PTB01", row$ExposureCol, row$UgCol, covs))

  missing_needed <- setdiff(needed, names(dat))
  if (length(missing_needed) > 0) {
    stop(
      "Missing variables needed to reconstruct primary analysis sample for ",
      row$Polymer, ": ", paste(missing_needed, collapse = ", "),
      call. = FALSE
    )
  }

  dsub <- dat %>%
    select(all_of(needed)) %>%
    filter(if_all(everything(), ~ !is.na(.x)))

  z <- dsub[[row$ExposureCol]]
  x <- dsub[[row$UgCol]]
  x_pos <- x[is.finite(x) & x > 0]

  q25_z <- safe_quantile(z, 0.25)
  q50_z <- safe_quantile(z, 0.50)
  q75_z <- safe_quantile(z, 0.75)

  q25_x_pos <- safe_quantile(x_pos, 0.25)
  q50_x_pos <- safe_quantile(x_pos, 0.50)
  q75_x_pos <- safe_quantile(x_pos, 0.75)

  iqr_rows[[i]] <- tibble(
    Polymer = row$Polymer,
    PolymerLabel = row$PolymerLabel,
    ExposureCol = row$ExposureCol,
    UgCol = row$UgCol,
    N_primary_reconstructed = nrow(dsub),
    Events_primary_reconstructed = sum(dsub$PTB01 == 1, na.rm = TRUE),
    NonEvents_primary_reconstructed = sum(dsub$PTB01 == 0, na.rm = TRUE),
    Q25_log10_1px = q25_z,
    Median_log10_1px = q50_z,
    Q75_log10_1px = q75_z,
    IQR_log10_1px = q75_z - q25_z,
    PositiveRaw_N = length(x_pos),
    PositiveRaw_P25_ugpg = q25_x_pos,
    PositiveRaw_Median_ugpg = q50_x_pos,
    PositiveRaw_P75_ugpg = q75_x_pos
  )
}

iqr_tbl <- bind_rows(iqr_rows)

# Verify reconstructed N against model output. Do not stop because a
# legacy script could contain an additional filter; flag mismatches.
n_check <- primary_rows %>%
  select(Polymer, N_reported = N, Events_reported = Events, NonEvents_reported = NonEvents) %>%
  left_join(
    iqr_tbl %>%
      select(
        Polymer,
        N_primary_reconstructed,
        Events_primary_reconstructed,
        NonEvents_primary_reconstructed
      ),
    by = "Polymer"
  ) %>%
  mutate(
    N_match = N_reported == N_primary_reconstructed,
    Events_match = Events_reported == Events_primary_reconstructed,
    NonEvents_match = NonEvents_reported == NonEvents_primary_reconstructed
  )

write_csv(n_check, file.path(AUDIT_DIR, "PTB_primary_analysis_sample_reconstruction_check.csv"))

if (any(!n_check$N_match, na.rm = TRUE)) {
  warning(
    "One or more reconstructed complete-case N values do not match the primary model output. ",
    "Inspect PTB_primary_analysis_sample_reconstruction_check.csv before treating IQR estimates as final."
  )
}

# ------------------------------------------------------------
# Rescale ALL model coefficients to primary-sample IQR units
# ------------------------------------------------------------

models_iqr <- models_recovered %>%
  left_join(
    iqr_tbl %>%
      select(
        Polymer,
        IQR_log10_1px,
        Q25_log10_1px,
        Median_log10_1px,
        Q75_log10_1px,
        PositiveRaw_N,
        PositiveRaw_P25_ugpg,
        PositiveRaw_Median_ugpg,
        PositiveRaw_P75_ugpg
      ),
    by = "Polymer"
  ) %>%
  mutate(
    logOR_IQR = beta_recovered * IQR_log10_1px,
    SE_logOR_IQR = se_beta_recovered * IQR_log10_1px,
    OR_IQR = exp(logOR_IQR),
    OR_IQR_lo = exp((beta_recovered - 1.96 * se_beta_recovered) * IQR_log10_1px),
    OR_IQR_hi = exp((beta_recovered + 1.96 * se_beta_recovered) * IQR_log10_1px),
    OR_IQR_CI = sprintf("%.2f (%.2f-%.2f)", OR_IQR, OR_IQR_lo, OR_IQR_hi)
  )

# ------------------------------------------------------------
# Exact raw x -> 2x contrasts for PRIMARY model
# ------------------------------------------------------------

primary_contrasts <- models_iqr %>%
  filter(ModelTag == PRIMARY_MODEL_TAG) %>%
  rowwise() %>%
  mutate(
    delta_true2x_P25 = true_double_delta(PositiveRaw_P25_ugpg),
    delta_true2x_Median = true_double_delta(PositiveRaw_Median_ugpg),
    delta_true2x_P75 = true_double_delta(PositiveRaw_P75_ugpg),

    OR_true2x_P25 = exp(beta_recovered * delta_true2x_P25),
    OR_true2x_P25_lo = exp((beta_recovered - 1.96 * se_beta_recovered) * delta_true2x_P25),
    OR_true2x_P25_hi = exp((beta_recovered + 1.96 * se_beta_recovered) * delta_true2x_P25),

    OR_true2x_Median = exp(beta_recovered * delta_true2x_Median),
    OR_true2x_Median_lo = exp((beta_recovered - 1.96 * se_beta_recovered) * delta_true2x_Median),
    OR_true2x_Median_hi = exp((beta_recovered + 1.96 * se_beta_recovered) * delta_true2x_Median),

    OR_true2x_P75 = exp(beta_recovered * delta_true2x_P75),
    OR_true2x_P75_lo = exp((beta_recovered - 1.96 * se_beta_recovered) * delta_true2x_P75),
    OR_true2x_P75_hi = exp((beta_recovered + 1.96 * se_beta_recovered) * delta_true2x_P75),

    old_delta_log10_2 = delta_old,
    delta_ratio_P25_vs_old = delta_true2x_P25 / delta_old,
    delta_ratio_Median_vs_old = delta_true2x_Median / delta_old,
    delta_ratio_P75_vs_old = delta_true2x_P75 / delta_old
  ) %>%
  ungroup()

# ------------------------------------------------------------
# Compact primary audit table
# ------------------------------------------------------------

primary_audit <- primary_contrasts %>%
  select(
    Polymer, PolymerLabel, ExposureCol,
    N, Events, NonEvents,
    beta_recovered, se_beta_recovered,
    OR_2x, OR_2x_lo, OR_2x_hi,
    IQR_log10_1px, OR_IQR, OR_IQR_lo, OR_IQR_hi,
    PositiveRaw_N,
    PositiveRaw_P25_ugpg, OR_true2x_P25, OR_true2x_P25_lo, OR_true2x_P25_hi,
    PositiveRaw_Median_ugpg, OR_true2x_Median, OR_true2x_Median_lo, OR_true2x_Median_hi,
    PositiveRaw_P75_ugpg, OR_true2x_P75, OR_true2x_P75_lo, OR_true2x_P75_hi,
    delta_ratio_P25_vs_old, delta_ratio_Median_vs_old, delta_ratio_P75_vs_old,
    p, q_BH
  )

write_csv(
  primary_audit,
  file.path(AUDIT_DIR, "PTB_effect_scale_comparison_primary.csv")
)

write_csv(
  models_iqr,
  file.path(AUDIT_DIR, "PTB_ALL_models_IQR_rescaled.csv")
)

write_csv(
  primary_contrasts,
  file.path(AUDIT_DIR, "PTB_true_doubling_reference_contrasts.csv")
)

write_csv(
  iqr_tbl,
  file.path(AUDIT_DIR, "PTB_primary_IQR_reference_values.csv")
)

# ------------------------------------------------------------
# Candidate Figure 3 using OR per IQR of log10(1+ug/g)
# ------------------------------------------------------------

fig_dat <- models_iqr %>%
  filter(ModelTag %in% c(PRIMARY_MODEL_TAG, DISPLAY_SENS_MODELS)) %>%
  mutate(
    sig_q05 = ifelse(is.na(q_BH), FALSE, q_BH < 0.05),
    PolymerLabel = case_when(
      Polymer == "Cumulative_MNPs" ~ "Total MNP burden",
      TRUE ~ PolymerLabel
    )
  )

# Keep the same ordering concept as current Figure 3, but order on
# the corrected candidate display scale.
poly_order_primary <- fig_dat %>%
  filter(ModelTag == PRIMARY_MODEL_TAG, Polymer != "Cumulative_MNPs") %>%
  arrange(desc(OR_IQR)) %>%
  pull(PolymerLabel)

if (length(poly_order_primary) >= 4) {
  poly_levels <- c(
    poly_order_primary[1:4],
    "Total MNP burden",
    poly_order_primary[5:length(poly_order_primary)]
  )
} else {
  poly_levels <- c(poly_order_primary, "Total MNP burden")
}
poly_levels <- unique(poly_levels[!is.na(poly_levels)])

base_positions <- tibble(
  PolymerLabel = poly_levels,
  y_base = rev(seq_along(poly_levels))
)

primary_df <- fig_dat %>%
  filter(ModelTag == PRIMARY_MODEL_TAG) %>%
  left_join(base_positions, by = "PolymerLabel")

offset_map <- tibble(
  ModelTag = DISPLAY_SENS_MODELS,
  y_offset = c(-0.18, 0.00, 0.18)
)

sens_df <- fig_dat %>%
  filter(ModelTag %in% DISPLAY_SENS_MODELS) %>%
  left_join(base_positions, by = "PolymerLabel") %>%
  left_join(offset_map, by = "ModelTag") %>%
  mutate(y = y_base + y_offset)

# Dynamic x-range based on displayed CIs.
all_lows <- primary_df$OR_IQR_lo[is.finite(primary_df$OR_IQR_lo) & primary_df$OR_IQR_lo > 0]
all_highs <- primary_df$OR_IQR_hi[is.finite(primary_df$OR_IQR_hi) & primary_df$OR_IQR_hi > 0]

X_LOWER <- if (length(all_lows) > 0) max(0.02, 10^(floor(log10(min(all_lows)) * 2) / 2)) else 0.05
X_UPPER <- if (length(all_highs) > 0) min(20, 10^(ceiling(log10(max(all_highs)) * 2) / 2)) else 4

# Preserve current Figure 3 colors if present in the polished workflow.
COLOR_MAP <- c(
  "Crude" = "#56B4E9",
  "Primary adjusted, excluding CLR/PCA outliers" = "#7A3E9D",
  "Primary + HTN subtype" = "#E69F00"
)

subtitle_txt <- paste(
  "Primary model: maternal age + BMI + SDI + twins + prior PTB + hypertension + hospital.",
  "Effects are rescaled from the existing log10(1 + µg/g) models to a polymer-specific IQR."
)

p <- ggplot() +
  geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.4, color = "grey30") +
  geom_segment(
    data = primary_df,
    aes(x = OR_IQR_lo, xend = OR_IQR_hi, y = y_base, yend = y_base),
    linewidth = 0.85,
    color = "black"
  ) +
  geom_point(
    data = sens_df,
    aes(x = OR_IQR, y = y, color = ModelLabel),
    size = 2.3,
    alpha = 0.95
  ) +
  geom_point(
    data = primary_df,
    aes(x = OR_IQR, y = y_base, shape = sig_q05),
    size = 3.0,
    color = "black",
    fill = "black"
  ) +
  scale_x_log10() +
  coord_cartesian(xlim = c(X_LOWER, X_UPPER)) +
  scale_y_continuous(
    breaks = base_positions$y_base,
    labels = base_positions$PolymerLabel,
    expand = expansion(mult = c(0.03, 0.03))
  ) +
  scale_color_manual(
    values = COLOR_MAP,
    breaks = c(
      "Crude",
      "Primary adjusted, excluding CLR/PCA outliers",
      "Primary + HTN subtype"
    )
  ) +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17)) +
  labs(
    title = "Polymer-specific associations with preterm birth",
    subtitle = stringr::str_wrap(subtitle_txt, width = 82),
    x = "OR per IQR increase in log10(1 + placental burden [µg/g]) (95% CI)",
    y = NULL,
    color = "Sensitivity models",
    shape = "Primary model\nBH q < 0.05"
  ) +
  guides(
    color = guide_legend(order = 1, override.aes = list(size = 3)),
    shape = guide_legend(order = 2)
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 10.8, lineheight = 1.05),
    legend.position = "right",
    legend.box = "vertical"
  )

FIG_PNG <- file.path(AUDIT_FIG_DIR, "fig_primary_ptb_models__figure3_IQR_audit.png")
FIG_PDF <- file.path(AUDIT_FIG_DIR, "fig_primary_ptb_models__figure3_IQR_audit.pdf")

ggsave(FIG_PNG, p, width = 9.0, height = 5.9, dpi = 300)
ggsave(FIG_PDF, p, width = 9.0, height = 5.9)

# ------------------------------------------------------------
# Suggested manuscript wording (NOT automatically applied)
# ------------------------------------------------------------

methods_text <- paste0(
  "Placental polymer burdens were modeled as log10(1 + µg/g). ",
  "For continuous regression models, effect estimates were expressed as odds ratios ",
  "per interquartile-range increase in the transformed polymer burden. ",
  "Because log10(1 + x) is nonlinear on the raw concentration scale, a fixed ",
  "log10(2) increment corresponds to a twofold increase in (1 + x), not an exact ",
  "doubling of x."
)

results_note <- paste0(
  "Reviewer-point-1 audit: the fitted coefficients, standard errors, p values, and ",
  "BH-adjusted q values are unchanged. Only the reporting scale is re-expressed. ",
  "Use PTB_effect_scale_comparison_primary.csv to compare the former log10(2)-scaled ",
  "ORs with ORs per IQR and exact x-to-2x contrasts at representative observed burdens."
)

writeLines(
  c(
    "PROPOSED METHODS WORDING",
    methods_text,
    "",
    "INTERNAL RESULTS NOTE",
    results_note
  ),
  file.path(AUDIT_DIR, "reviewer_point1_manuscript_wording.txt")
)

# ------------------------------------------------------------
# Console summary
# ------------------------------------------------------------

message("")
message("============================================================")
message("Reviewer Point 1 effect-scale audit complete")
message("============================================================")
message("Project root: ", PROJECT_ROOT)
message("Existing model directory: ", MODEL_DIR)
message("Audit tables: ", AUDIT_DIR)
message("Candidate figure directory: ", AUDIT_FIG_DIR)
message("")
message("Primary outputs:")
message("  - PTB_effect_scale_comparison_primary.csv")
message("  - PTB_ALL_models_IQR_rescaled.csv")
message("  - PTB_true_doubling_reference_contrasts.csv")
message("  - PTB_primary_IQR_reference_values.csv")
message("  - PTB_primary_analysis_sample_reconstruction_check.csv")
message("  - reviewer_point1_manuscript_wording.txt")
message("  - ", FIG_PNG)
message("  - ", FIG_PDF)
message("")
message("IMPORTANT: This script does not overwrite Figure 3 or any existing model output.")
