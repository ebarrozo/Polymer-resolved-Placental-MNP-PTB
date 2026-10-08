#!/usr/bin/env Rscript
# INTERNAL FIT STAGE — NOT A MANUSCRIPT-FACING EFFECT SCALE
# -----------------------------------------------------------------------------
# This script is retained to reproduce the original fitted model coefficients and
# analysis samples. Some intermediate columns/labels use the historical fixed
# log10(2) scaling. Do not report those intermediate effect estimates. The
# manuscript-facing scripts in ../ re-express continuous effects on the final
# IQR scale without refitting the models.
# -----------------------------------------------------------------------------
# ============================================================
# 06_primary_ptb_models.R
# Primary PTB logistic models + sensitivity models
# ============================================================
# Input:
#   outputs/tables/mnp_ptb_cohort__analysis_ready_model_table.csv
#   outputs/tables/mnp_ptb_cohort__polymer_detection_summary.csv
#
# Outputs:
#   outputs/tables/06_primary_ptb_models__<stamp>/
#     table_primary_ptb_models.csv
#     table_primary_ptb_models.html
#     table_primary_ptb_models.png
#     PTB__ALL_models_long.csv
#     PTB__DELTA_vs_M0.csv
#     model_covariate_audit.csv
#     model_sparse_polymer_audit.csv
#
#   outputs/figs/fig_primary_ptb_models__forest_primary.png
#   outputs/figs/fig_primary_ptb_models__forest_all_sensitivity.png
#
# Notes:
#   - Outcome is PTB <37w, encoded as PTB01 where 1=Preterm, 0=Term.
#   - Models intentionally do NOT adjust for GA because GA defines the outcome.
#   - Exposure effect is reported as OR per doubling in µg/g tissue.
#   - Uses log10(ug/g + 1) exposure columns and scales beta by log10(2).
#   - Sparse polymers, especially PS, are retained in audit outputs but skipped
#     from continuous regression when detected n < threshold.
#   - Primary BH q-values are adjusted across eligible polymers within M1.
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
  library(stringr)
  library(ggplot2)
  library(scales)
  library(gt)
})

source("00_config.R")
source(file.path("R", "utils_mnp.R"))

# ----------------------------
# Settings
# ----------------------------
IN_CSV <- file.path(DIR_TABLES, "mnp_ptb_cohort__analysis_ready_model_table.csv")
IN_DETECT <- file.path(DIR_TABLES, "mnp_ptb_cohort__polymer_detection_summary.csv")
IN_CLUSTER_OUTLIERS <- file.path(DIR_TABLES, "clustering__clr_pca_outlier_audit.csv")

STAMP <- format(Sys.time(), "%Y%m%d_%H%M%S")
OUT_DIR <- file.path(DIR_TABLES, paste0("06_primary_ptb_models__", STAMP))
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

OUT_PRIMARY <- file.path(OUT_DIR, "table_primary_ptb_models.csv")
OUT_PRIMARY_HTML <- file.path(OUT_DIR, "table_primary_ptb_models.html")
OUT_PRIMARY_PNG <- file.path(OUT_DIR, "table_primary_ptb_models.png")
OUT_ALL <- file.path(OUT_DIR, "PTB__ALL_models_long.csv")
OUT_DELTA <- file.path(OUT_DIR, "PTB__DELTA_vs_M0.csv")
OUT_COV_AUDIT <- file.path(OUT_DIR, "model_covariate_audit.csv")
OUT_SPARSE <- file.path(OUT_DIR, "model_sparse_polymer_audit.csv")
OUT_CLUSTER_OUTLIER_AUDIT <- file.path(OUT_DIR, "model_cluster_outlier_sensitivity_audit.csv")

OUT_FOREST_PRIMARY <- file.path(DIR_FIGS, "fig_primary_ptb_models__forest_primary.png")
OUT_FOREST_ALL <- file.path(DIR_FIGS, "fig_primary_ptb_models__forest_all_sensitivity.png")
OUT_LOG <- file.path(DIR_LOGS, "06_primary_ptb_models__run_log.txt")

BASE_SIZE <- 11
K2 <- log10(2)

# Covariate settings.
# Keep primary model parsimonious and clinically defensible.
COV_M1_PRIMARY <- c("Maternal_age", "BMI", "SDI_score", "Twins01", "Prior_PTB", "Hypertension01", "Hospital")
COV_M1_PLUS_SVI <- c("Maternal_age", "BMI", "SDI_score", "SVI_score", "Twins01", "Prior_PTB", "Hypertension01", "Hospital")
COV_M1_NO_SDI <- c("Maternal_age", "BMI", "Twins01", "Prior_PTB", "Hypertension01", "Hospital")
COV_M1_PLUS_HTN_SUBTYPE <- c("Maternal_age", "BMI", "SDI_score", "Twins01", "Prior_PTB", "Hypertension_type", "Hospital")
COV_M1_PLUS_RACE <- c("Maternal_age", "BMI", "SDI_score", "Twins01", "Prior_PTB", "Hypertension01", "Hospital", "Ethnicity_Race")

MODEL_SPECS <- list(
  M0_CRUDE = character(0),
  M1_ADJ_PRIMARY = COV_M1_PRIMARY,
  M1_ADJ_PRIMARY_NO_CLR_OUTLIERS = COV_M1_PRIMARY,
  M1_PLUS_SVI = COV_M1_PLUS_SVI,
  M1_NO_SDI = COV_M1_NO_SDI,
  M1_PLUS_HTNsubtype = COV_M1_PLUS_HTN_SUBTYPE,
  M1_PLUS_RACE = COV_M1_PLUS_RACE
)

MODEL_EXCLUDE_CLUSTER_OUTLIERS <- c(
  M0_CRUDE = FALSE,
  M1_ADJ_PRIMARY = FALSE,
  M1_ADJ_PRIMARY_NO_CLR_OUTLIERS = TRUE,
  M1_PLUS_SVI = FALSE,
  M1_NO_SDI = FALSE,
  M1_PLUS_HTNsubtype = FALSE,
  M1_PLUS_RACE = FALSE
)

MODEL_LABELS <- c(
  M0_CRUDE = "Crude",
  M1_ADJ_PRIMARY = "Primary adjusted",
  M1_ADJ_PRIMARY_NO_CLR_OUTLIERS = "Primary adjusted, excluding CLR/PCA outliers",
  M1_PLUS_SVI = "Primary + SVI",
  M1_NO_SDI = "Primary without SDI",
  M1_PLUS_HTNsubtype = "Primary + HTN subtype",
  M1_PLUS_RACE = "Primary + race/ethnicity"
)

PRIMARY_MODEL_TAG <- "M1_ADJ_PRIMARY"

# ----------------------------
# Helpers
# ----------------------------
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
      Preterm = PTB01 == 1L,

      Maternal_age = safe_num(Maternal_age),
      BMI = safe_num(BMI),
      SDI_score = safe_num(SDI_score),
      SVI_score = safe_num(SVI_score),
      Twins01 = as.integer(as_01(Twins01)),
      Prior_PTB = as.integer(as_01(Prior_PTB)),
      Hypertension01 = as.integer(as_01(Hypertension01)),
      Hospital = factor(safe_chr(Hospital)),
      Ethnicity_Race = factor(safe_chr(Ethnicity_Race)),
      Hypertension_type = factor(safe_chr(Hypertension_type))
    )
}

detected_n <- function(x) sum(is.finite(x) & x > 0, na.rm = TRUE)

fmt_or_ci <- function(or, lo, hi) {
  ifelse(is.na(or), NA_character_, sprintf("%.2f (%.2f-%.2f)", or, lo, hi))
}

fmt_p <- function(p) {
  ifelse(is.na(p), NA_character_,
         ifelse(p < 0.001, "<0.001", formatC(p, format = "f", digits = 3)))
}

# Drop covariates that are absent, all missing, constant, or otherwise unusable.
usable_covariates <- function(dat, covs) {
  keep <- c()
  dropped <- tibble(covariate = character(), reason = character())

  for (cv in covs) {
    if (!cv %in% names(dat)) {
      dropped <- bind_rows(dropped, tibble(covariate = cv, reason = "not_in_data"))
      next
    }

    x <- dat[[cv]]
    nonmiss <- !is.na(x)

    if (sum(nonmiss) == 0) {
      dropped <- bind_rows(dropped, tibble(covariate = cv, reason = "all_missing"))
      next
    }

    if (is.numeric(x) || is.integer(x)) {
      if (sd(as.numeric(x), na.rm = TRUE) == 0 || length(unique(x[nonmiss])) < 2) {
        dropped <- bind_rows(dropped, tibble(covariate = cv, reason = "constant_numeric"))
        next
      }
    } else {
      if (length(unique(as.character(x[nonmiss]))) < 2) {
        dropped <- bind_rows(dropped, tibble(covariate = cv, reason = "constant_factor"))
        next
      }
    }

    keep <- c(keep, cv)
  }

  list(keep = keep, dropped = dropped)
}

# Build model frame and require enough observations/events.
make_model_data <- function(df, exposure_col, covariates) {
  vars <- c("PTB01", exposure_col, covariates)
  d <- df %>%
    select(all_of(vars)) %>%
    mutate(
      PTB01 = as.integer(PTB01),
      exposure = safe_num(.data[[exposure_col]])
    ) %>%
    select(-all_of(exposure_col)) %>%
    drop_na()

  # remove unused factor levels
  d <- d %>% mutate(across(where(is.factor), droplevels))
  d
}

fit_one_model <- function(df, exposure_col, polymer_label, model_tag, covariates) {
  d0 <- df %>%
    mutate(exposure_raw_ug_g = safe_num(.data[[str_replace(exposure_col, "_log10ugpg$", "_ug_per_g")]]))

  cov_use <- usable_covariates(d0, covariates)
  covariates_kept <- cov_use$keep

  d <- make_model_data(d0, exposure_col, covariates_kept)

  n <- nrow(d)
  events <- sum(d$PTB01 == 1, na.rm = TRUE)
  nonevents <- sum(d$PTB01 == 0, na.rm = TRUE)

  base_out <- tibble(
    ModelTag = model_tag,
    ModelLabel = MODEL_LABELS[[model_tag]],
    Polymer = polymer_label,
    PolymerLabel = ifelse(polymer_label == "Cumulative_MNPs", "Total MNP burden", polymer_label),
    ExposureCol = exposure_col,
    N = n,
    Events = events,
    NonEvents = nonevents,
    Covariates = paste(covariates_kept, collapse = " + "),
    OR_2x = NA_real_,
    OR_2x_lo = NA_real_,
    OR_2x_hi = NA_real_,
    logOR_2x = NA_real_,
    SE_logOR_2x = NA_real_,
    z = NA_real_,
    p = NA_real_,
    fit_status = NA_character_
  )

  if (n < 25 || events < 5 || nonevents < 5) {
    return(base_out %>% mutate(fit_status = "skipped_too_few_observations_or_events"))
  }

  if (sd(d$exposure, na.rm = TRUE) == 0) {
    return(base_out %>% mutate(fit_status = "skipped_zero_exposure_variance"))
  }

  # formula
  rhs <- c("exposure", covariates_kept)
  fml <- as.formula(paste("PTB01 ~", paste(rhs, collapse = " + ")))

  fit <- try(
    glm(fml, data = d, family = binomial(), control = glm.control(maxit = 100, epsilon = 1e-8)),
    silent = TRUE
  )

  if (inherits(fit, "try-error")) {
    return(base_out %>% mutate(fit_status = paste0("glm_error: ", as.character(fit)[1])))
  }

  # robust HC3 if available; fallback model-based.
  b_est <- NA_real_
  se_est <- NA_real_
  p_val <- NA_real_
  status_val <- NA_character_

  # Use HC1 by default. HC3 can become unstable with high-leverage samples
  # and produced warnings in this dataset. If HC1 fails or is non-finite,
  # fall back to model-based SEs.
  if (requireNamespace("sandwich", quietly = TRUE) && requireNamespace("lmtest", quietly = TRUE)) {
    ct <- try(
      suppressWarnings(lmtest::coeftest(fit, vcov. = sandwich::vcovHC(fit, type = "HC1"))),
      silent = TRUE
    )
    if (!inherits(ct, "try-error") && "exposure" %in% rownames(ct)) {
      b_est <- unname(ct["exposure", "Estimate"])
      se_est <- unname(ct["exposure", "Std. Error"])
      p_val <- 2 * pnorm(abs(b_est / se_est), lower.tail = FALSE)
      status_val <- "ok_HC1"
    }
  }

  if (!is.finite(b_est) || !is.finite(se_est) || !is.finite(p_val)) {
    sm <- summary(fit)$coefficients
    if (!"exposure" %in% rownames(sm)) {
      return(base_out %>% mutate(fit_status = "exposure_coefficient_missing"))
    }
    b_est <- unname(sm["exposure", "Estimate"])
    se_est <- unname(sm["exposure", "Std. Error"])
    p_val <- unname(sm["exposure", "Pr(>|z|)"])
    status_val <- "ok_model_based_fallback"
  }

  logOR_2x_val <- b_est * K2
  se_logOR_2x_val <- se_est * K2

  base_out %>%
    mutate(
      OR_2x = exp(logOR_2x_val),
      OR_2x_lo = exp(logOR_2x_val - 1.96 * se_logOR_2x_val),
      OR_2x_hi = exp(logOR_2x_val + 1.96 * se_logOR_2x_val),
      logOR_2x = logOR_2x_val,
      SE_logOR_2x = se_logOR_2x_val,
      z = logOR_2x_val / se_logOR_2x_val,
      p = p_val,
      fit_status = status_val
    )
}

make_gt_primary <- function(tbl) {
  tbl %>%
    transmute(
      Polymer = PolymerLabel,
      N,
      Events,
      `OR per doubling (95% CI)` = OR_CI,
      `p-value` = p_fmt,
      `BH q-value` = q_fmt,
      `Fit status` = fit_status
    ) %>%
    gt::gt() %>%
    gt::tab_header(title = gt::md("**Primary Adjusted PTB Models**")) %>%
    gt::tab_source_note(gt::md("Outcome: preterm birth <37 weeks. Effect estimates are ORs per doubling in µg/g tissue. BH q-values are adjusted across eligible polymers in the primary model.")) %>%
    gt::cols_align(align = "left", columns = Polymer) %>%
    gt::cols_align(align = "center", columns = -Polymer) %>%
    gt::tab_options(table.font.size = gt::px(11))
}

save_gt_html_png <- function(gt_obj, html_path, png_path) {
  gt::gtsave(gt_obj, filename = html_path)
  png_try <- try(gt::gtsave(gt_obj, filename = png_path), silent = TRUE)
  if (inherits(png_try, "try-error")) {
    warning(
      "Could not save PNG table: ", png_path, "\n",
      "HTML was saved. For PNG export, install webshot2 and ensure Chrome/Chromium is available.\n",
      "Error: ", as.character(png_try)
    )
  }
  invisible(TRUE)
}

# ----------------------------
# Load / audit
# ----------------------------
df0 <- read_model_table(IN_CSV)
df <- df0 %>% filter(!is.na(PTB01))

cluster_outlier_ids <- character(0)
if (file.exists(IN_CLUSTER_OUTLIERS)) {
  cluster_outliers <- readr::read_csv(
    IN_CLUSTER_OUTLIERS,
    show_col_types = FALSE,
    col_types = cols(
      SampleID = col_character(),
      PregnancyID = col_character(),
      .default = col_guess()
    )
  ) %>%
    mutate(SampleID = as.character(SampleID), PregnancyID = as.character(PregnancyID))

  cluster_outlier_ids <- unique(cluster_outliers$SampleID)

  write_csv_clean(
    cluster_outliers %>%
      mutate(used_for_sensitivity_exclusion = SampleID %in% df$SampleID),
    OUT_CLUSTER_OUTLIER_AUDIT
  )
} else {
  write_csv_clean(
    tibble(
      SampleID = character(),
      PregnancyID = character(),
      used_for_sensitivity_exclusion = logical(),
      note = "No clustering__clr_pca_outlier_audit.csv found; outlier-exclusion sensitivity model uses full data."
    ),
    OUT_CLUSTER_OUTLIER_AUDIT
  )
}

df <- df %>%
  mutate(cluster_clr_pca_outlier = SampleID %in% cluster_outlier_ids)

if (!all(c(0, 1) %in% unique(df$PTB01))) {
  stop("Need both term and preterm records for PTB models.", call. = FALSE)
}

detect_tbl <- if (file.exists(IN_DETECT)) {
  readr::read_csv(IN_DETECT, show_col_types = FALSE) %>%
    mutate(Polymer = toupper(as.character(Polymer)))
} else {
  tibble(
    Polymer = POLYMERS,
    n_detected_gt0 = vapply(POLYMERS, function(p) detected_n(df[[legacy_col_ug_g(p)]]), numeric(1)),
    continuous_model_ok = n_detected_gt0 >= MIN_DETECTED_N_FOR_CONTINUOUS_MODEL
  )
}

# Candidate exposures: total + polymers.
exposure_tbl <- tibble(
  Polymer = c("Cumulative_MNPs", POLYMERS),
  ExposureCol = c("Cumulative_MNPs_log10ugpg", legacy_col_log10(POLYMERS)),
  UgCol = c("Cumulative_MNPs_ug_per_g", legacy_col_ug_g(POLYMERS))
) %>%
  mutate(
    n_detected_gt0 = vapply(UgCol, function(cn) if (cn %in% names(df)) detected_n(df[[cn]]) else 0, numeric(1)),
    continuous_model_ok = Polymer == "Cumulative_MNPs" | n_detected_gt0 >= MIN_DETECTED_N_FOR_CONTINUOUS_MODEL,
    skip_reason = case_when(
      !ExposureCol %in% names(df) ~ "missing_exposure_column",
      !continuous_model_ok ~ "sparse_detection",
      TRUE ~ NA_character_
    )
  )

write_csv_clean(exposure_tbl, OUT_SPARSE)

eligible_exposures <- exposure_tbl %>%
  filter(is.na(skip_reason))

if (nrow(eligible_exposures) == 0) {
  stop("No eligible exposures for modeling.", call. = FALSE)
}

# Covariate audit on full analytic table.
covariate_audit <- purrr::imap_dfr(MODEL_SPECS, function(covs, tag) {
  cov_use <- usable_covariates(df, covs)

  if (length(covs) == 0) {
    return(tibble(
      ModelTag = as.character(tag),
      ModelLabel = as.character(MODEL_LABELS[[tag]]),
      covariate = NA_character_,
      kept = NA,
      reason = "crude_model_no_covariates"
    ))
  }

  tibble(
    ModelTag = as.character(tag),
    ModelLabel = as.character(MODEL_LABELS[[tag]]),
    covariate = as.character(covs)
  ) %>%
    mutate(
      kept = covariate %in% cov_use$keep,
      reason = dplyr::case_when(
        kept ~ "kept",
        covariate %in% cov_use$dropped$covariate ~ as.character(cov_use$dropped$reason[match(covariate, cov_use$dropped$covariate)]),
        TRUE ~ "not_evaluated"
      ),
      reason = as.character(reason)
    )
})

write_csv_clean(covariate_audit, OUT_COV_AUDIT)

# ----------------------------
# Fit models
# ----------------------------
all_models <- purrr::imap_dfr(MODEL_SPECS, function(covs, model_tag) {
  df_model <- df
  if (isTRUE(MODEL_EXCLUDE_CLUSTER_OUTLIERS[[model_tag]])) {
    df_model <- df_model %>% filter(!cluster_clr_pca_outlier)
  }

  purrr::pmap_dfr(
    eligible_exposures,
    function(Polymer, ExposureCol, UgCol, n_detected_gt0, continuous_model_ok, skip_reason) {
      fit_one_model(df_model, ExposureCol, Polymer, model_tag, covs)
    }
  )
}) %>%
  group_by(ModelTag) %>%
  mutate(
    q_BH = ifelse(
      fit_status %in% c("ok_HC1", "ok_HC3", "ok_model_based", "ok_model_based_fallback"),
      p.adjust(p, method = "BH"),
      NA_real_
    )
  ) %>%
  ungroup() %>%
  mutate(
    OR_CI = fmt_or_ci(OR_2x, OR_2x_lo, OR_2x_hi),
    p_fmt = fmt_p(p),
    q_fmt = fmt_p(q_BH),
    sig_q05 = q_BH < 0.05,
    PolymerLabel = case_when(
      Polymer == "Cumulative_MNPs" ~ "Total MNP burden",
      Polymer %in% names(POLYMER_FULL_NAMES) ~ paste0(POLYMER_FULL_NAMES[Polymer], " (", Polymer, ")"),
      TRUE ~ PolymerLabel
    )
  )

write_csv_clean(all_models, OUT_ALL)

primary_tbl <- all_models %>%
  filter(ModelTag == PRIMARY_MODEL_TAG) %>%
  arrange(desc(Polymer == "Cumulative_MNPs"), OR_2x) %>%
  mutate(
    ModelDescription = "Primary adjusted model"
  )

write_csv_clean(primary_tbl, OUT_PRIMARY)
save_gt_html_png(make_gt_primary(primary_tbl), OUT_PRIMARY_HTML, OUT_PRIMARY_PNG)

# Delta vs crude
delta_tbl <- all_models %>%
  filter(ModelTag %in% c("M0_CRUDE", PRIMARY_MODEL_TAG)) %>%
  select(ModelTag, Polymer, PolymerLabel, OR_2x, OR_2x_lo, OR_2x_hi, p, q_BH, N, fit_status) %>%
  pivot_wider(
    names_from = ModelTag,
    values_from = c(OR_2x, OR_2x_lo, OR_2x_hi, p, q_BH, N, fit_status),
    names_sep = "__"
  ) %>%
  mutate(
    OR_ratio_M1_vs_M0 = OR_2x__M1_ADJ_PRIMARY / OR_2x__M0_CRUDE,
    abs_logOR_delta_M1_vs_M0 = abs(log(OR_2x__M1_ADJ_PRIMARY) - log(OR_2x__M0_CRUDE))
  )

write_csv_clean(delta_tbl, OUT_DELTA)

# ----------------------------
# Figures
# ----------------------------
poly_order <- primary_tbl %>%
  filter(!is.na(OR_2x), Polymer != "Cumulative_MNPs") %>%
  arrange(OR_2x) %>%
  pull(PolymerLabel)

# Add total burden to top/bottom explicitly.
primary_plot <- primary_tbl %>%
  filter(!is.na(OR_2x)) %>%
  mutate(
    PolymerLabel = factor(PolymerLabel, levels = c(poly_order, "Total MNP burden")),
    clipped_lo = pmax(OR_2x_lo, 0.05),
    clipped_hi = pmin(OR_2x_hi, 10),
    clipped = OR_2x_lo < 0.05 | OR_2x_hi > 10
  )

p_primary <- ggplot(primary_plot, aes(y = PolymerLabel, x = OR_2x)) +
  geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.4) +
  geom_errorbarh(aes(xmin = clipped_lo, xmax = clipped_hi), height = 0.18, linewidth = 0.75) +
  geom_point(aes(shape = sig_q05), size = 2.8) +
  geom_text(
    aes(label = ifelse(clipped, "CI clipped", "")),
    nudge_y = 0.22,
    size = 2.6
  ) +
  scale_x_log10(breaks = c(0.05, 0.1, 0.25, 0.5, 1, 2, 4, 10)) +
  coord_cartesian(xlim = c(0.05, 10)) +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17), na.value = 16) +
  labs(
    title = "Primary Adjusted PTB Models",
    subtitle = "OR per doubling in placental MNP burden; BH q-values adjusted across modeled analytes",
    x = "OR per doubling (95% CI)",
    y = NULL,
    shape = "BH q < 0.05"
  ) +
  theme_bw(base_size = BASE_SIZE) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"))

ggsave(OUT_FOREST_PRIMARY, p_primary, width = 8.2, height = 6.2, dpi = 300)

# All sensitivity overlay
sens_keep <- c("M0_CRUDE", "M1_ADJ_PRIMARY", "M1_ADJ_PRIMARY_NO_CLR_OUTLIERS", "M1_PLUS_SVI", "M1_NO_SDI", "M1_PLUS_HTNsubtype", "M1_PLUS_RACE")
sens_plot <- all_models %>%
  filter(ModelTag %in% sens_keep, !is.na(OR_2x), Polymer != "Cumulative_MNPs") %>%
  mutate(
    ModelLabel = factor(ModelLabel, levels = MODEL_LABELS[sens_keep]),
    PolymerLabel = factor(PolymerLabel, levels = poly_order),
    clipped_lo = pmax(OR_2x_lo, 0.05),
    clipped_hi = pmin(OR_2x_hi, 10)
  )

p_all <- ggplot(sens_plot, aes(y = PolymerLabel, x = OR_2x, shape = ModelLabel)) +
  geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.35) +
  geom_errorbarh(aes(xmin = clipped_lo, xmax = clipped_hi), height = 0.16, alpha = 0.65, position = position_dodge(width = 0.55)) +
  geom_point(position = position_dodge(width = 0.55), size = 2.1, alpha = 0.80) +
  scale_x_log10(breaks = c(0.05, 0.1, 0.25, 0.5, 1, 2, 4, 10)) +
  coord_cartesian(xlim = c(0.05, 10)) +
  scale_shape_manual(values = c(16, 17, 15, 3, 7, 8, 4)) +
  labs(
    title = "PTB Model Sensitivity Analysis",
    subtitle = "OR per doubling in placental MNP burden",
    x = "OR per doubling (95% CI)",
    y = NULL,
    shape = "Model"
  ) +
  theme_bw(base_size = BASE_SIZE) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"), legend.position = "bottom")

ggsave(OUT_FOREST_ALL, p_all, width = 9.2, height = 7.2, dpi = 300)

# ----------------------------
# Log
# ----------------------------
run_log <- c(
  "06_primary_ptb_models.R",
  paste0("PROJECT_ROOT: ", PROJECT_ROOT),
  paste0("Rows in model table: ", nrow(df0)),
  paste0("Rows with PTB outcome: ", nrow(df)),
  paste0("Cluster CLR/PCA outliers available for sensitivity exclusion: ", length(cluster_outlier_ids)),
  paste0("Cluster CLR/PCA outliers overlapping model data: ", sum(df$cluster_clr_pca_outlier, na.rm = TRUE)),
  "",
  "Outcome counts:",
  capture.output(print(table(df$Delivery_class, useNA = "ifany"))),
  "",
  "Eligible exposures:",
  capture.output(print(eligible_exposures)),
  "",
  "Primary model:",
  capture.output(print(primary_tbl %>% select(Polymer, OR_2x, OR_2x_lo, OR_2x_hi, p, q_BH, N, fit_status), n = Inf)),
  "",
  "Wrote:",
  OUT_PRIMARY,
  OUT_PRIMARY_HTML,
  OUT_PRIMARY_PNG,
  OUT_ALL,
  OUT_DELTA,
  OUT_COV_AUDIT,
  OUT_SPARSE,
  OUT_CLUSTER_OUTLIER_AUDIT,
  OUT_FOREST_PRIMARY,
  OUT_FOREST_ALL
)

cat(paste(run_log, collapse = "\n"), file = OUT_LOG)
message(paste(run_log, collapse = "\n"))
