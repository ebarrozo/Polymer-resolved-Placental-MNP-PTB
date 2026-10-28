#!/usr/bin/env Rscript
# POLISHED VERSION: sources R/figure_theme_mnp.R and writes *_polished outputs where applicable.
# ============================================================
# 07_ptb_robustness_figure4.R
# Robustness / Figure 4 module for PTB-MNP models
# ============================================================
# Run after:
#   source("run_primary_models.R")
#
# Inputs:
#   outputs/tables/mnp_ptb_cohort__analysis_ready_model_table.csv
#   latest outputs/tables/06_primary_ptb_models__*/PTB__ALL_models_long.csv
#
# Outputs:
#   outputs/tables/07_ptb_robustness__<timestamp>/
#     robustness_primary_vs_sensitivity_summary.csv
#     robustness_leave_one_out_models.csv
#     robustness_leave_one_out_summary.csv
#     robustness_influence_sample_frequency.csv
#     robustness_model_key.csv
#
#   outputs/figs/fig4_ptb_robustness__model_sensitivity_heatmap.png
#   outputs/figs/fig4_ptb_robustness__leave_one_out_forest.png
#   outputs/figs/fig4_ptb_robustness__leave_one_out_range.png
#
# Notes:
#   - Main Figure 3 is the clean primary/sensitivity point figure.
#   - This module is designed as a robustness figure set.
#   - Leave-one-out is run for key polymers only by default:
#       PU, PET, SBR, ABS, plus cumulative burden.
#   - Effect scale is OR per doubling in µg/g tissue.
#   - Primary model retains all eligible samples.
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
  library(stringr)
  library(ggplot2)
  library(scales)
})

source("00_config.R")
source(file.path("R", "utils_mnp.R"))
source(file.path("R", "figure_theme_mnp.R"))

# ----------------------------
# Settings
# ----------------------------
IN_MODEL_TABLE <- file.path(DIR_TABLES, "mnp_ptb_cohort__analysis_ready_model_table.csv")

STAMP <- format(Sys.time(), "%Y%m%d_%H%M%S")
OUT_DIR <- file.path(DIR_TABLES, paste0("07_ptb_robustness__", STAMP))
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

OUT_SUMMARY <- file.path(OUT_DIR, "robustness_primary_vs_sensitivity_summary_polished.csv")
OUT_LOO <- file.path(OUT_DIR, "robustness_leave_one_out_models_polished.csv")
OUT_LOO_SUMMARY <- file.path(OUT_DIR, "robustness_leave_one_out_summary_polished.csv")
OUT_INFLUENCE <- file.path(OUT_DIR, "robustness_influence_sample_frequency_polished.csv")
OUT_MODEL_KEY <- file.path(OUT_DIR, "robustness_model_key_polished.csv")

OUT_HEATMAP <- file.path(DIR_FIGS, "fig4_ptb_robustness__model_sensitivity_heatmap_polished.png")
OUT_LOO_FOREST <- file.path(DIR_FIGS, "fig4_ptb_robustness__leave_one_out_forest_polished.png")
OUT_LOO_RANGE <- file.path(DIR_FIGS, "fig4_ptb_robustness__leave_one_out_range_polished.png")
OUT_LOG <- file.path(DIR_LOGS, "07_ptb_robustness_figure4__run_log_polished.txt")

PRIMARY_TAG <- "M1_ADJ_PRIMARY"
KEY_POLYMERS <- c("PU", "PET", "SBR", "ABS", "Cumulative_MNPs")
K2 <- log10(2)

PRIMARY_COVARIATES <- c("Maternal_age", "BMI", "SDI_score", "Twins01", "Prior_PTB", "Hypertension01", "Hospital")
BASE_SIZE <- 11

# ----------------------------
# Helpers
# ----------------------------
find_latest_primary_dir <- function() {
  dirs <- list.dirs(DIR_TABLES, recursive = FALSE, full.names = TRUE)
  dirs <- dirs[grepl("^06_primary_ptb_models__", basename(dirs))]
  if (length(dirs) == 0) {
    stop("No 06_primary_ptb_models__* directory found. Run run_primary_models.R first.", call. = FALSE)
  }
  dirs[order(file.info(dirs)$mtime, decreasing = TRUE)][1]
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
      Delivery_class = factor(as.character(Delivery_class), levels = c("Term", "Preterm")),
      PTB01 = case_when(
        Delivery_class == "Preterm" ~ 1L,
        Delivery_class == "Term" ~ 0L,
        TRUE ~ NA_integer_
      ),
      Maternal_age = safe_num(Maternal_age),
      BMI = safe_num(BMI),
      SDI_score = safe_num(SDI_score),
      Twins01 = as.integer(as_01(Twins01)),
      Prior_PTB = as.integer(as_01(Prior_PTB)),
      Hypertension01 = as.integer(as_01(Hypertension01)),
      Hospital = factor(safe_chr(Hospital))
    ) %>%
    filter(!is.na(PTB01))
}

fmt_p <- function(p) {
  ifelse(is.na(p), NA_character_,
         ifelse(p < 0.001, "<0.001", formatC(p, format = "f", digits = 3)))
}

poly_label <- function(poly) {
  if (poly == "Cumulative_MNPs") return("Cumulative MNP burden")
  if (poly %in% names(POLYMER_FULL_NAMES)) return(paste0(POLYMER_FULL_NAMES[[poly]], " (", poly, ")"))
  poly
}

usable_covariates <- function(dat, covs) {
  keep <- c()
  for (cv in covs) {
    if (!cv %in% names(dat)) next
    x <- dat[[cv]]
    nonmiss <- !is.na(x)
    if (sum(nonmiss) == 0) next
    if (is.numeric(x) || is.integer(x)) {
      if (length(unique(x[nonmiss])) < 2 || sd(as.numeric(x), na.rm = TRUE) == 0) next
    } else {
      if (length(unique(as.character(x[nonmiss]))) < 2) next
    }
    keep <- c(keep, cv)
  }
  keep
}

fit_ptb_model <- function(dat, exposure_col, polymer, covariates = PRIMARY_COVARIATES) {
  covs <- usable_covariates(dat, covariates)
  vars <- c("SampleID", "PTB01", exposure_col, covs)

  d <- dat %>%
    select(all_of(vars)) %>%
    mutate(
      exposure = safe_num(.data[[exposure_col]]),
      PTB01 = as.integer(PTB01)
    ) %>%
    select(-all_of(exposure_col)) %>%
    drop_na() %>%
    mutate(across(where(is.factor), droplevels))

  base <- tibble(
    Polymer = polymer,
    PolymerLabel = poly_label(polymer),
    ExposureCol = exposure_col,
    N = nrow(d),
    Events = sum(d$PTB01 == 1, na.rm = TRUE),
    NonEvents = sum(d$PTB01 == 0, na.rm = TRUE),
    Covariates = paste(covs, collapse = " + "),
    OR_2x = NA_real_,
    OR_2x_lo = NA_real_,
    OR_2x_hi = NA_real_,
    p = NA_real_,
    fit_status = NA_character_
  )

  if (nrow(d) < 25 || base$Events < 5 || base$NonEvents < 5 || sd(d$exposure, na.rm = TRUE) == 0) {
    return(base %>% mutate(fit_status = "skipped_insufficient_data"))
  }

  fml <- as.formula(paste("PTB01 ~", paste(c("exposure", covs), collapse = " + ")))

  fit <- try(glm(fml, data = d, family = binomial(), control = glm.control(maxit = 100)), silent = TRUE)

  if (inherits(fit, "try-error")) {
    return(base %>% mutate(fit_status = paste0("glm_error: ", as.character(fit)[1])))
  }

  b_est <- NA_real_
  se_est <- NA_real_
  p_val <- NA_real_
  status <- NA_character_

  if (requireNamespace("sandwich", quietly = TRUE) && requireNamespace("lmtest", quietly = TRUE)) {
    ct <- try(
      suppressWarnings(lmtest::coeftest(fit, vcov. = sandwich::vcovHC(fit, type = "HC1"))),
      silent = TRUE
    )

    if (!inherits(ct, "try-error") && "exposure" %in% rownames(ct)) {
      b_est <- unname(ct["exposure", "Estimate"])
      se_est <- unname(ct["exposure", "Std. Error"])
      p_val <- 2 * pnorm(abs(b_est / se_est), lower.tail = FALSE)
      status <- "ok_HC1"
    }
  }

  if (!is.finite(b_est) || !is.finite(se_est) || !is.finite(p_val)) {
    sm <- summary(fit)$coefficients
    if (!"exposure" %in% rownames(sm)) {
      return(base %>% mutate(fit_status = "exposure_coefficient_missing"))
    }
    b_est <- unname(sm["exposure", "Estimate"])
    se_est <- unname(sm["exposure", "Std. Error"])
    p_val <- unname(sm["exposure", "Pr(>|z|)"])
    status <- "ok_model_based_fallback"
  }

  logOR <- b_est * K2
  se_logOR <- se_est * K2

  base %>%
    mutate(
      OR_2x = exp(logOR),
      OR_2x_lo = exp(logOR - 1.96 * se_logOR),
      OR_2x_hi = exp(logOR + 1.96 * se_logOR),
      p = p_val,
      fit_status = status
    )
}

# ----------------------------
# Load latest primary outputs
# ----------------------------
latest_dir <- find_latest_primary_dir()
all_models_path <- file.path(latest_dir, "PTB__ALL_models_long.csv")
all_models <- readr::read_csv(all_models_path, show_col_types = FALSE)

df <- read_model_table(IN_MODEL_TABLE)

# ----------------------------
# Primary vs sensitivity summary
# ----------------------------
sens_summary <- all_models %>%
  filter(ModelTag %in% c(
    "M0_CRUDE",
    "M1_ADJ_PRIMARY",
    "M1_ADJ_PRIMARY_NO_CLR_OUTLIERS",
    "M1_PLUS_SVI",
    "M1_NO_SDI",
    "M1_PLUS_HTNsubtype",
    "M1_PLUS_RACE"
  )) %>%
  mutate(
    PolymerLabel = ifelse(Polymer == "Cumulative_MNPs", "Cumulative MNP burden", PolymerLabel),
    ModelLabel = factor(ModelLabel)
  ) %>%
  group_by(Polymer, PolymerLabel) %>%
  summarise(
    n_models = sum(!is.na(OR_2x)),
    primary_OR = OR_2x[match(TRUE, ModelTag == PRIMARY_TAG)],
    primary_lo = OR_2x_lo[match(TRUE, ModelTag == PRIMARY_TAG)],
    primary_hi = OR_2x_hi[match(TRUE, ModelTag == PRIMARY_TAG)],
    primary_p = p[match(TRUE, ModelTag == PRIMARY_TAG)],
    primary_q = q_BH[match(TRUE, ModelTag == PRIMARY_TAG)],
    OR_min = min(OR_2x, na.rm = TRUE),
    OR_max = max(OR_2x, na.rm = TRUE),
    all_same_direction_vs_null = all((OR_2x > 1) == (primary_OR > 1), na.rm = TRUE),
    any_CI_crosses_1 = any(OR_2x_lo <= 1 & OR_2x_hi >= 1, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    primary_q_sig = primary_q < 0.05,
    robustness_class = case_when(
      primary_q_sig & all_same_direction_vs_null ~ "Primary significant; direction stable",
      primary_q_sig & !all_same_direction_vs_null ~ "Primary significant; direction unstable",
      !primary_q_sig & all_same_direction_vs_null ~ "Primary not significant; direction stable",
      TRUE ~ "Primary not significant; direction unstable"
    )
  )

write_csv_clean(sens_summary, OUT_SUMMARY)

model_key <- all_models %>%
  filter(ModelTag %in% c(
    "M0_CRUDE",
    "M1_ADJ_PRIMARY",
    "M1_ADJ_PRIMARY_NO_CLR_OUTLIERS",
    "M1_PLUS_SVI",
    "M1_NO_SDI",
    "M1_PLUS_HTNsubtype",
    "M1_PLUS_RACE"
  )) %>%
  group_by(ModelTag, ModelLabel) %>%
  summarise(
    Covariates = first(Covariates),
    N_range = paste0(min(N, na.rm = TRUE), " to ", max(N, na.rm = TRUE)),
    .groups = "drop"
  )

write_csv_clean(model_key, OUT_MODEL_KEY)

# ----------------------------
# Leave-one-out robustness
# ----------------------------
exposure_map <- tibble(
  Polymer = KEY_POLYMERS,
  ExposureCol = c("PU_log10ugpg", "PET_log10ugpg", "SBR_log10ugpg", "ABS_log10ugpg", "Cumulative_MNPs_log10ugpg")
) %>%
  filter(ExposureCol %in% names(df))

loo_models <- purrr::pmap_dfr(exposure_map, function(Polymer, ExposureCol) {
  message("Running leave-one-out for ", Polymer)

  primary <- fit_ptb_model(df, ExposureCol, Polymer, PRIMARY_COVARIATES) %>%
    mutate(leave_out_sample = "NONE", analysis = "primary_all_samples")

  loo <- purrr::map_dfr(df$SampleID, function(id) {
    fit_ptb_model(df %>% filter(SampleID != id), ExposureCol, Polymer, PRIMARY_COVARIATES) %>%
      mutate(leave_out_sample = id, analysis = "leave_one_out")
  })

  bind_rows(primary, loo)
}) %>%
  group_by(Polymer) %>%
  mutate(q_BH_within_polymer = p.adjust(p, method = "BH")) %>%
  ungroup()

write_csv_clean(loo_models, OUT_LOO)

loo_summary <- loo_models %>%
  filter(analysis == "leave_one_out") %>%
  group_by(Polymer, PolymerLabel) %>%
  summarise(
    n_loo = sum(!is.na(OR_2x)),
    OR_median = median(OR_2x, na.rm = TRUE),
    OR_min = min(OR_2x, na.rm = TRUE),
    OR_max = max(OR_2x, na.rm = TRUE),
    OR_q1 = quantile(OR_2x, 0.25, na.rm = TRUE, names = FALSE),
    OR_q3 = quantile(OR_2x, 0.75, na.rm = TRUE, names = FALSE),
    p_median = median(p, na.rm = TRUE),
    p_min = min(p, na.rm = TRUE),
    p_max = max(p, na.rm = TRUE),
    pct_OR_same_direction_as_primary = {
      primary_or <- loo_models %>%
        filter(Polymer == first(Polymer), analysis == "primary_all_samples") %>%
        pull(OR_2x) %>%
        .[1]
      100 * mean((OR_2x > 1) == (primary_or > 1), na.rm = TRUE)
    },
    .groups = "drop"
  )

write_csv_clean(loo_summary, OUT_LOO_SUMMARY)

# Samples that create the largest absolute logOR shift for any key polymer.
influence <- loo_models %>%
  filter(analysis == "leave_one_out") %>%
  group_by(Polymer) %>%
  mutate(
    primary_OR = loo_models$OR_2x[match(Polymer, loo_models$Polymer[loo_models$analysis == "primary_all_samples"])],
    abs_logOR_shift = abs(log(OR_2x) - log(primary_OR))
  ) %>%
  ungroup() %>%
  filter(is.finite(abs_logOR_shift)) %>%
  group_by(leave_out_sample) %>%
  summarise(
    n_key_polymers_top10_shift = sum(abs_logOR_shift >= quantile(abs_logOR_shift, 0.90, na.rm = TRUE)),
    max_abs_logOR_shift = max(abs_logOR_shift, na.rm = TRUE),
    polymers_with_top_shift = paste(unique(Polymer[abs_logOR_shift >= quantile(abs_logOR_shift, 0.90, na.rm = TRUE)]), collapse = ", "),
    .groups = "drop"
  ) %>%
  arrange(desc(n_key_polymers_top10_shift), desc(max_abs_logOR_shift))

write_csv_clean(influence, OUT_INFLUENCE)

# ----------------------------
# Figures
# ----------------------------
# Figure 4A: sensitivity heatmap for all modeled polymers
heat_dat <- all_models %>%
  filter(ModelTag %in% c(
    "M0_CRUDE",
    "M1_ADJ_PRIMARY",
    "M1_ADJ_PRIMARY_NO_CLR_OUTLIERS",
    "M1_PLUS_SVI",
    "M1_NO_SDI",
    "M1_PLUS_HTNsubtype",
    "M1_PLUS_RACE"
  )) %>%
  mutate(
    PolymerLabel = ifelse(Polymer == "Cumulative_MNPs", "Cumulative MNP burden", PolymerLabel),
    ModelLabel = factor(ModelLabel, levels = c(
      "Crude",
      "Primary adjusted",
      "Primary adjusted, excluding CLR/PCA outliers",
      "Primary + SVI",
      "Primary without SDI",
      "Primary + HTN subtype",
      "Primary + race/ethnicity"
    )),
    log2OR = log2(OR_2x),
    sig = q_BH < 0.05 | p < 0.05
  )

poly_order <- heat_dat %>%
  filter(ModelTag == PRIMARY_TAG) %>%
  arrange(desc(log2OR)) %>%
  pull(PolymerLabel)

heat_dat <- heat_dat %>%
  mutate(PolymerLabel = factor(PolymerLabel, levels = rev(poly_order)))

p_heat <- ggplot(heat_dat, aes(x = ModelLabel, y = PolymerLabel, fill = log2OR)) +
  geom_tile(color = "white", linewidth = 0.35) +
  geom_text(aes(label = ifelse(sig, "*", "")), size = 4) +
  scale_fill_gradient2(
    midpoint = 0,
    labels = label_number(accuracy = 0.01),
    name = "log2(OR)"
  ) +
  labs(
    title = "PTB model robustness across covariate specifications",
    subtitle = "* indicates nominal p < 0.05 or BH q < 0.05 within each model",
    x = NULL,
    y = NULL
  ) +
  theme_mnp_forest(base_size = BASE_SIZE) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 35, hjust = 1),
    plot.title = element_text(face = "bold")
  )

ggsave_mnp(OUT_HEATMAP, p_heat, width = 9.5, height = 6.2, dpi = 300)

# Figure 4B: leave-one-out forest for key polymers
loo_plot <- loo_models %>%
  mutate(
    PolymerLabel = factor(PolymerLabel, levels = rev(unique(loo_summary$PolymerLabel))),
    is_primary = analysis == "primary_all_samples"
  )

p_loo <- ggplot() +
  geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.35) +
  geom_point(
    data = loo_plot %>% filter(analysis == "leave_one_out"),
    aes(x = OR_2x, y = PolymerLabel),
    alpha = 0.18,
    size = 1.5
  ) +
  geom_errorbarh(
    data = loo_plot %>% filter(is_primary),
    aes(xmin = OR_2x_lo, xmax = OR_2x_hi, y = PolymerLabel),
    height = 0.20,
    linewidth = 0.8
  ) +
  geom_point(
    data = loo_plot %>% filter(is_primary),
    aes(x = OR_2x, y = PolymerLabel),
    size = 3.0
  ) +
  scale_x_log10(
    breaks = c(0.05, 0.10, 0.20, 0.50, 1, 2, 4, 8),
    labels = c("0.05", "0.10", "0.20", "0.50", "1", "2", "4", "8")
  ) +
  labs(
    title = "Leave-one-out robustness for key PTB-associated polymers",
    subtitle = "Grey points show one-sample-deleted ORs; black point and CI show the primary model",
    x = "OR per doubling",
    y = NULL
  ) +
  theme_mnp_forest(base_size = BASE_SIZE) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"))

ggsave_mnp(OUT_LOO_FOREST, p_loo, width = 8.4, height = 4.8, dpi = 300)

# Figure 4C: leave-one-out OR range
range_plot <- loo_summary %>%
  mutate(
    PolymerLabel = factor(PolymerLabel, levels = PolymerLabel[order(OR_median)])
  )

p_range <- ggplot(range_plot, aes(y = PolymerLabel)) +
  geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.35) +
  geom_segment(aes(x = OR_min, xend = OR_max, yend = PolymerLabel), linewidth = 1.0, alpha = 0.55) +
  geom_point(aes(x = OR_median), size = 2.8) +
  scale_x_log10(
    breaks = c(0.05, 0.10, 0.20, 0.50, 1, 2, 4, 8),
    labels = c("0.05", "0.10", "0.20", "0.50", "1", "2", "4", "8")
  ) +
  labs(
    title = "Range of leave-one-out OR estimates",
    subtitle = "Segments show minimum to maximum one-sample-deleted OR; points show median",
    x = "Leave-one-out OR per doubling",
    y = NULL
  ) +
  theme_mnp_forest(base_size = BASE_SIZE) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"))

ggsave_mnp(OUT_LOO_RANGE, p_range, width = 7.8, height = 4.8, dpi = 300)

# ----------------------------
# Log
# ----------------------------
run_log <- c(
  "07_ptb_robustness_figure4.R",
  paste0("PROJECT_ROOT: ", PROJECT_ROOT),
  paste0("Latest primary model dir: ", latest_dir),
  paste0("Rows in model table: ", nrow(df)),
  paste0("Leave-one-out polymers: ", paste(exposure_map$Polymer, collapse = ", ")),
  "",
  "Sensitivity summary:",
  capture.output(print(sens_summary)),
  "",
  "LOO summary:",
  capture.output(print(loo_summary)),
  "",
  "Wrote:",
  OUT_SUMMARY,
  OUT_LOO,
  OUT_LOO_SUMMARY,
  OUT_INFLUENCE,
  OUT_MODEL_KEY,
  OUT_HEATMAP,
  OUT_LOO_FOREST,
  OUT_LOO_RANGE
)

cat(paste(run_log, collapse = "\n"), file = OUT_LOG)
message(paste(run_log, collapse = "\n"))
