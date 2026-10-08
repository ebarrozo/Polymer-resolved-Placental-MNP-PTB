#!/usr/bin/env Rscript
# INTERNAL FIT STAGE — NOT A MANUSCRIPT-FACING EFFECT SCALE
# -----------------------------------------------------------------------------
# This script is retained to reproduce the original fitted model coefficients and
# analysis samples. Some intermediate columns/labels use the historical fixed
# log10(2) scaling. Do not report those intermediate effect estimates. The
# manuscript-facing scripts in ../ re-express continuous effects on the final
# IQR scale without refitting the models.
# -----------------------------------------------------------------------------
# POLISHED VERSION: sources R/figure_theme_mnp.R and writes *_polished outputs where applicable.
# ============================================================
# 08_labor_outcome_supplement.R
# Parsimonious labor outcome models:
# Spontaneous augmented labor vs spontaneous non-augmented labor
# v3.5
# ============================================================
#
# Run after:
#   source("run_analysis_setup.R")
#
# Input:
#   outputs/tables/mnp_ptb_cohort__analysis_ready_model_table.csv
#   outputs/tables/mnp_ptb_cohort__polymer_detection_summary.csv
#
# Outcome derivation:
#   Derived strictly from Type_of_Labor:
#     "Spontaneous: augmented" -> 1
#     "Spontaneous"            -> 0
#     "No labor"               -> excluded
#     missing/other            -> excluded
#
# Primary labor model:
#   Overall eligible laboring deliveries
#   Adjustment: gestational age at delivery only
#
# Sensitivity models:
#   Crude
#   GA + hospital
#   GA + SDI
#   GA + prior PTB
#   GA + hypertension
#
# Rationale:
#   The labor outcome has few augmented events. Fully adjusted models
#   are unstable and can produce quasi-separation, huge ORs, and huge CIs.
#   These models are therefore exploratory and intentionally parsimonious.
#
# Outputs:
#   outputs/tables/08_labor_models__<timestamp>/
#     labor_outcome_derivation_audit.csv
#     labor_model_all_results.csv
#     labor_model_primary_GA_only_results.csv
#     labor_model_multiple_testing_audit.csv
#     labor_model_key.csv
#     labor_model_sparse_polymer_audit.csv
#     labor_model_primary_GA_only_results.html
#     labor_model_primary_GA_only_results.png
#
#   outputs/figs/
#     supp_labor_models__type_of_labor_counts.png
#     supp_labor_models__GA_only_primary_forest.png
#     supp_labor_models__GA_only_primary_plus_parsimonious_sensitivities.png
#     supp_labor_models__GA_only_primary_forest_full_axis.png
#     supp_labor_models__GA_only_primary_plus_parsimonious_sensitivities_full_axis.png
#
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
source(file.path("R", "figure_theme_mnp.R"))

# ----------------------------
# Settings
# ----------------------------
IN_CSV <- file.path(DIR_TABLES, "mnp_ptb_cohort__analysis_ready_model_table.csv")
IN_DETECT <- file.path(DIR_TABLES, "mnp_ptb_cohort__polymer_detection_summary.csv")

STAMP <- format(Sys.time(), "%Y%m%d_%H%M%S")
OUT_DIR <- file.path(DIR_TABLES, paste0("08_labor_models__", STAMP))
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

OUT_DERIVE <- file.path(OUT_DIR, "labor_outcome_derivation_audit_polished.csv")
OUT_ALL <- file.path(OUT_DIR, "labor_model_all_results.csv")
OUT_PRIMARY <- file.path(OUT_DIR, "labor_model_primary_GA_only_results.csv")
OUT_Q_AUDIT <- file.path(OUT_DIR, "labor_model_multiple_testing_audit_polished.csv")
OUT_SPARSE <- file.path(OUT_DIR, "labor_model_sparse_polymer_audit_polished.csv")
OUT_KEY <- file.path(OUT_DIR, "labor_model_key_polished.csv")
OUT_PRIMARY_HTML <- file.path(OUT_DIR, "labor_model_primary_GA_only_results_polished.html")
OUT_PRIMARY_PNG <- file.path(OUT_DIR, "labor_model_primary_GA_only_results_polished.png")

OUT_FIG_COUNTS <- file.path(DIR_FIGS, "supp_labor_models__type_of_labor_counts_polished.png")
OUT_FIG_PRIMARY <- file.path(DIR_FIGS, "supp_labor_models__GA_only_primary_forest_polished.png")
OUT_FIG_SENS <- file.path(DIR_FIGS, "supp_labor_models__GA_only_primary_plus_parsimonious_sensitivities_polished.png")
OUT_FIG_PRIMARY_FULL_AXIS <- file.path(DIR_FIGS, "supp_labor_models__GA_only_primary_forest_full_axis_polished.png")
OUT_FIG_SENS_FULL_AXIS <- file.path(DIR_FIGS, "supp_labor_models__GA_only_primary_plus_parsimonious_sensitivities_full_axis_polished.png")
OUT_LEGEND <- file.path(DIR_TABLES, "supp_labor_models_figure_legend_polished.txt")
OUT_LOG <- file.path(DIR_LOGS, "08_labor_outcome_supplement__run_log_polished.txt")

K2 <- log10(2)
BASE_SIZE <- 11

# Use these fixed display limits for manuscript-style plots.
# Full-axis versions are also generated from actual finite CI ranges.
DISPLAY_XMIN <- 0.03
DISPLAY_XMAX <- 10
DISPLAY_BREAKS <- c(0.03, 0.05, 0.10, 0.20, 0.50, 1, 2, 4, 10)

# Pre-specified labor-hypothesis polymers from the manuscript framework.
# This is a secondary/sensitivity FDR family, not the main all-analyte FDR.
LABOR_PRESPEC_POLYMERS <- c("PVC", "N66")

MODEL_SPECS <- list(
  CRUDE = list(
    label = "Crude",
    covariates = character(0),
    note = "Unadjusted exploratory labor model."
  ),
  GA_ONLY = list(
    label = "GA-adjusted primary",
    covariates = c("GA"),
    note = "Primary parsimonious labor model adjusted only for gestational age at delivery."
  ),
  GA_PLUS_HOSPITAL = list(
    label = "GA + hospital",
    covariates = c("GA", "Hospital"),
    note = "Site/clinical-system sensitivity model."
  ),
  GA_PLUS_SDI = list(
    label = "GA + SDI",
    covariates = c("GA", "SDI_score"),
    note = "Social deprivation sensitivity model."
  ),
  GA_PLUS_PRIOR_PTB = list(
    label = "GA + prior PTB",
    covariates = c("GA", "Prior_PTB"),
    note = "Obstetric-history sensitivity model."
  ),
  GA_PLUS_HYPERTENSION = list(
    label = "GA + hypertension",
    covariates = c("GA", "Hypertension01"),
    note = "Hypertension sensitivity model."
  )
)

PRIMARY_MODEL_TAG <- "GA_ONLY"

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
      GA = safe_num(GA),
      Maternal_age = safe_num(Maternal_age),
      BMI = safe_num(BMI),
      SDI_score = safe_num(SDI_score),
      SVI_score = safe_num(SVI_score),
      Twins01 = as.integer(as_01(Twins01)),
      Prior_PTB = as.integer(as_01(Prior_PTB)),
      Hypertension01 = as.integer(as_01(Hypertension01)),
      Hospital = factor(safe_chr(Hospital)),
      Type_of_Labor = dplyr::coalesce(
        safe_chr(Type_of_Labor),
        safe_chr(`Type of labor`)
      ),
      Type_of_Labor_raw = Type_of_Labor
    )
}

fmt_or_ci <- function(or, lo, hi) {
  ifelse(is.na(or), NA_character_, sprintf("%.2f (%.2f-%.2f)", or, lo, hi))
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

poly_family_local <- function(poly) {
  case_when(
    poly == "Cumulative_MNPs" ~ "Cumulative burden",
    poly %in% c("PE", "PP") ~ "Polyolefins",
    poly %in% c("PET", "PC", "PMMA") ~ "Oxygen-containing polymers",
    poly %in% c("PVC") ~ "Halogenated polymers",
    poly %in% c("N6", "N66") ~ "Polyamides",
    poly %in% c("PU") ~ "Polyurethane",
    poly %in% c("SBR", "ABS", "PS") ~ "Styrenic/rubber polymers",
    TRUE ~ "Other"
  )
}

derive_labor_outcome <- function(df) {
  labor_raw <- safe_chr(df$Type_of_Labor_raw)
  labor_norm <- labor_raw %>%
    str_to_lower() %>%
    str_replace_all("[^a-z0-9]+", " ") %>%
    str_squish()

  augmented <- str_detect(labor_norm, "^spontaneou?s augmented$|^spontaneou?s augment|augmented|augmentation")
  spontaneous_no_aug <- str_detect(labor_norm, "^spontaneou?s$|^spontaneous labor$") & !augmented
  no_labor <- str_detect(labor_norm, "^no labor$|^none$|^no$")

  labor_augmented_binary <- case_when(
    augmented ~ 1L,
    spontaneous_no_aug ~ 0L,
    no_labor ~ NA_integer_,
    TRUE ~ NA_integer_
  )

  labor_outcome_label <- case_when(
    labor_augmented_binary == 1L ~ "Spontaneous labor with augmentation",
    labor_augmented_binary == 0L ~ "Spontaneous labor without augmentation",
    no_labor ~ "Excluded: no labor",
    is.na(labor_raw) | labor_raw == "" ~ "Excluded/unknown: missing Type_of_Labor",
    TRUE ~ paste0("Excluded/unknown: Type_of_Labor = ", labor_raw)
  )

  df %>%
    mutate(
      Type_of_Labor_clean = labor_norm,
      spontaneous_labor_augmented01 = labor_augmented_binary,
      labor_outcome_label = labor_outcome_label
    )
}

detected_n <- function(x) sum(is.finite(x) & x > 0, na.rm = TRUE)

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

make_model_subset <- function(df) {
  df %>%
    filter(!is.na(spontaneous_labor_augmented01))
}

fit_labor_model <- function(df, exposure_col, polymer, model_tag, model_label, covariates, note) {
  d0 <- make_model_subset(df)
  covs <- usable_covariates(d0, covariates)

  vars <- c("SampleID", "spontaneous_labor_augmented01", exposure_col, covs)

  d <- d0 %>%
    select(all_of(vars)) %>%
    mutate(
      outcome = as.integer(spontaneous_labor_augmented01),
      exposure = safe_num(.data[[exposure_col]])
    ) %>%
    select(-all_of(c(exposure_col, "spontaneous_labor_augmented01"))) %>%
    drop_na() %>%
    mutate(across(where(is.factor), droplevels))

  n <- nrow(d)
  events <- sum(d$outcome == 1, na.rm = TRUE)
  nonevents <- sum(d$outcome == 0, na.rm = TRUE)

  base <- tibble(
    ModelTag = model_tag,
    ModelLabel = model_label,
    ModelNote = note,
    Polymer = polymer,
    PolymerLabel = poly_label(polymer),
    ExposureCol = exposure_col,
    N = n,
    Events_Augmented = events,
    NonEvents_SpontaneousNoAug = nonevents,
    Covariates = paste(covs, collapse = " + "),
    OR_2x = NA_real_,
    OR_2x_lo = NA_real_,
    OR_2x_hi = NA_real_,
    p = NA_real_,
    fit_status = NA_character_
  )

  if (n < 20 || events < 5 || nonevents < 5) {
    return(base %>% mutate(fit_status = "skipped_too_few_observations_or_events"))
  }

  if (sd(d$exposure, na.rm = TRUE) == 0) {
    return(base %>% mutate(fit_status = "skipped_zero_exposure_variance"))
  }

  rhs <- c("exposure", covs)
  fml <- as.formula(paste("outcome ~", paste(rhs, collapse = " + ")))

  fit <- try(glm(fml, data = d, family = binomial(), control = glm.control(maxit = 100, epsilon = 1e-8)), silent = TRUE)

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

add_multiple_testing_columns <- function(res) {
  ok_status <- c("ok_HC1", "ok_model_based_fallback", "ok_model_based")

  res %>%
    mutate(
      polymer_family = poly_family_local(Polymer),
      is_cumulative = Polymer == "Cumulative_MNPs",
      is_polymer = !is_cumulative,
      is_labor_prespecified = Polymer %in% LABOR_PRESPEC_POLYMERS,
      p_for_adjust = ifelse(fit_status %in% ok_status, p, NA_real_)
    ) %>%
    group_by(ModelTag) %>%
    mutate(
      q_BH_all_analytes = {
        q <- rep(NA_real_, n())
        idx <- !is.na(p_for_adjust)
        q[idx] <- p.adjust(p_for_adjust[idx], method = "BH")
        q
      },
      q_BH_polymers_only = {
        q <- rep(NA_real_, n())
        idx <- !is.na(p_for_adjust) & is_polymer
        q[idx] <- p.adjust(p_for_adjust[idx], method = "BH")
        q
      },
      q_BH_labor_prespec = {
        q <- rep(NA_real_, n())
        idx <- !is.na(p_for_adjust) & is_labor_prespecified
        if (sum(idx) > 0) q[idx] <- p.adjust(p_for_adjust[idx], method = "BH")
        q
      }
    ) %>%
    ungroup() %>%
    group_by(ModelTag, polymer_family) %>%
    mutate(
      q_BH_within_family = {
        q <- rep(NA_real_, n())
        idx <- !is.na(p_for_adjust) & is_polymer
        if (sum(idx) > 0) q[idx] <- p.adjust(p_for_adjust[idx], method = "BH")
        q
      }
    ) %>%
    ungroup() %>%
    mutate(
      q_BH = q_BH_all_analytes,
      nominal_p05 = p < 0.05,
      sig_q05_all_analytes = q_BH_all_analytes < 0.05,
      sig_q05_polymers_only = q_BH_polymers_only < 0.05,
      sig_q05_within_family = q_BH_within_family < 0.05,
      sig_q05_labor_prespec = q_BH_labor_prespec < 0.05,
      sig_q05 = sig_q05_all_analytes,
      OR_CI = fmt_or_ci(OR_2x, OR_2x_lo, OR_2x_hi),
      p_fmt = fmt_p(p),
      q_fmt = fmt_p(q_BH_all_analytes),
      q_all_fmt = fmt_p(q_BH_all_analytes),
      q_poly_fmt = fmt_p(q_BH_polymers_only),
      q_family_fmt = fmt_p(q_BH_within_family),
      q_labor_prespec_fmt = fmt_p(q_BH_labor_prespec)
    )
}

make_gt_table <- function(tbl) {
  tbl %>%
    transmute(
      Polymer = PolymerLabel,
      N,
      `Augmented / spontaneous` = paste0(Events_Augmented, "/", NonEvents_SpontaneousNoAug),
      `OR per doubling (95% CI)` = fmt_or_ci(OR_2x, OR_2x_lo, OR_2x_hi),
      `p-value` = fmt_p(p),
      `BH q, all analytes` = fmt_p(q_BH_all_analytes),
      `BH q, polymers only` = fmt_p(q_BH_polymers_only),
      `BH q, within family` = fmt_p(q_BH_within_family),
      `BH q, labor prespec` = fmt_p(q_BH_labor_prespec),
      `Fit status` = fit_status
    ) %>%
    gt::gt() %>%
    gt::tab_header(title = gt::md("**Primary labor models: spontaneous augmented labor**")) %>%
    gt::tab_source_note(gt::md("Primary labor model compares spontaneous augmented labor with spontaneous labor without augmentation across eligible laboring deliveries, adjusted only for gestational age at delivery.")) %>%
    gt::cols_align(align = "left", columns = Polymer) %>%
    gt::cols_align(align = "center", columns = -Polymer) %>%
    gt::tab_options(table.font.size = gt::px(11))
}

save_gt_html_png <- function(gt_obj, html_path, png_path) {
  gt::gtsave(gt_obj, filename = html_path)
  png_try <- try(gt::gtsave(gt_obj, filename = png_path), silent = TRUE)
  if (inherits(png_try, "try-error")) {
    warning("Could not save PNG table. HTML was saved. Install webshot2/Chrome for PNG export.")
  }
  invisible(TRUE)
}

make_axis_limits <- function(dat, xmin_floor = 0.03, xmax_cap = 100) {
  vals <- c(dat$OR_2x, dat$OR_2x_lo, dat$OR_2x_hi)
  vals <- vals[is.finite(vals) & vals > 0]
  if (length(vals) == 0) return(c(xmin_floor, 10))

  lo <- 10 ^ floor(log10(min(vals, na.rm = TRUE)))
  hi <- 10 ^ ceiling(log10(max(vals, na.rm = TRUE)))

  lo <- max(lo, xmin_floor)
  hi <- min(hi, xmax_cap)

  if (!is.finite(lo) || !is.finite(hi) || lo >= hi) c(xmin_floor, 10) else c(lo, hi)
}

common_log_breaks <- function(lims) {
  br <- c(0.001, 0.003, 0.01, 0.03, 0.05, 0.10, 0.20, 0.50, 1, 2, 4, 10, 20, 50, 100)
  br[br >= lims[1] & br <= lims[2]]
}

prep_forest_data <- function(tbl, order_levels = NULL, xmin = DISPLAY_XMIN, xmax = DISPLAY_XMAX) {
  dat <- tbl %>%
    filter(!is.na(OR_2x)) %>%
    mutate(
      PolymerLabel = ifelse(Polymer == "Cumulative_MNPs", "Cumulative MNP burden", PolymerLabel),
      xmin = pmax(OR_2x_lo, xmin),
      xmax = pmin(OR_2x_hi, xmax),
      x_plot = pmin(pmax(OR_2x, xmin), xmax),
      left_trunc = OR_2x_lo < xmin | OR_2x < xmin,
      right_trunc = OR_2x_hi > xmax | OR_2x > xmax
    )

  if (is.null(order_levels)) {
    order_levels <- dat %>% arrange(desc(OR_2x)) %>% pull(PolymerLabel)
  }

  dat %>%
    mutate(PolymerLabel = factor(PolymerLabel, levels = rev(order_levels)))
}

plot_primary_forest <- function(primary_results, out_path, full_axis = FALSE) {
  lims <- if (full_axis) make_axis_limits(primary_results, xmin_floor = 0.001, xmax_cap = 100) else c(DISPLAY_XMIN, DISPLAY_XMAX)
  breaks <- if (full_axis) common_log_breaks(lims) else DISPLAY_BREAKS

  plot_dat <- prep_forest_data(primary_results, xmin = lims[1], xmax = lims[2]) %>%
    mutate(
      evidence = case_when(
        sig_q05_all_analytes ~ "BH q < 0.05, all analytes",
        sig_q05_labor_prespec ~ "BH q < 0.05, prespecified PVC/N66",
        nominal_p05 ~ "nominal p < 0.05 only",
        TRUE ~ "not nominal"
      ),
      evidence = factor(
        evidence,
        levels = c(
          "not nominal",
          "nominal p < 0.05 only",
          "BH q < 0.05, prespecified PVC/N66",
          "BH q < 0.05, all analytes"
        )
      )
    )

  subtitle_txt <- str_wrap(
    "Primary model: gestational age only. Black points and 95% CI bars = GA-adjusted primary model. Point shape denotes evidence level. Arrows indicate OR or CI truncated for display.",
    width = 95
  )

  p <- ggplot(plot_dat, aes(y = PolymerLabel)) +
    geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.35, color = "grey30") +
    geom_segment(aes(x = xmin, xend = xmax, yend = PolymerLabel), linewidth = 0.85, color = "black") +
    geom_point(aes(x = x_plot, shape = evidence), size = 3.4, color = "black", stroke = 0.9) +
    geom_text(data = plot_dat %>% filter(left_trunc), aes(x = lims[1], label = "\u2190"), hjust = 1.05, size = 3.6) +
    geom_text(data = plot_dat %>% filter(right_trunc), aes(x = lims[2], label = "\u2192"), hjust = -0.05, size = 3.6) +
    scale_x_log10(breaks = breaks, limits = lims, labels = label_number(accuracy = 0.001, trim = TRUE)) +
    scale_shape_manual(
      values = c(
        "not nominal" = 16,
        "nominal p < 0.05 only" = 1,
        "BH q < 0.05, prespecified PVC/N66" = 23,
        "BH q < 0.05, all analytes" = 17
      ),
      name = "Evidence"
    ) +
    labs(
      title = "Associations with spontaneous labor requiring augmentation",
      subtitle = subtitle_txt,
      x = "OR per doubling (95% CI)",
      y = NULL
    ) +
    theme_mnp_forest(base_size = BASE_SIZE) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(size = 10.5),
      legend.position = "right",
      legend.title = element_text(face = "bold")
    )

  ggsave_mnp(out_path, p, width = 8.8, height = 6.6, dpi = 300)
}

plot_sensitivity_forest <- function(all_results, primary_results, out_path, full_axis = FALSE) {
  sens_tags <- c("CRUDE", "GA_PLUS_HOSPITAL", "GA_PLUS_SDI", "GA_PLUS_PRIOR_PTB", "GA_PLUS_HYPERTENSION")
  sens_all <- all_results %>% filter(ModelTag %in% sens_tags, !is.na(OR_2x))

  axis_dat <- bind_rows(primary_results, sens_all)
  lims <- if (full_axis) make_axis_limits(axis_dat, xmin_floor = 0.001, xmax_cap = 100) else c(DISPLAY_XMIN, DISPLAY_XMAX)
  breaks <- if (full_axis) common_log_breaks(lims) else DISPLAY_BREAKS

  order_levels <- primary_results %>%
    filter(!is.na(OR_2x)) %>%
    arrange(desc(OR_2x)) %>%
    mutate(PolymerLabel = ifelse(Polymer == "Cumulative_MNPs", "Cumulative MNP burden", PolymerLabel)) %>%
    pull(PolymerLabel)

  evidence_levels <- c(
    "not nominal",
    "nominal p < 0.05 only",
    "BH q < 0.05, prespecified PVC/N66",
    "BH q < 0.05, all analytes"
  )

  primary_dat <- prep_forest_data(primary_results, order_levels, xmin = lims[1], xmax = lims[2]) %>%
    mutate(
      evidence = case_when(
        sig_q05_all_analytes ~ "BH q < 0.05, all analytes",
        sig_q05_labor_prespec ~ "BH q < 0.05, prespecified PVC/N66",
        nominal_p05 ~ "nominal p < 0.05 only",
        TRUE ~ "not nominal"
      ),
      evidence = factor(evidence, levels = evidence_levels)
    )

  sens_dat <- prep_forest_data(sens_all, order_levels, xmin = lims[1], xmax = lims[2]) %>%
    mutate(
      ModelLabel = factor(ModelLabel, levels = c(
        "Crude",
        "GA + hospital",
        "GA + SDI",
        "GA + prior PTB",
        "GA + hypertension"
      )),
      evidence = case_when(
        sig_q05_all_analytes ~ "BH q < 0.05, all analytes",
        sig_q05_labor_prespec ~ "BH q < 0.05, prespecified PVC/N66",
        nominal_p05 ~ "nominal p < 0.05 only",
        TRUE ~ "not nominal"
      ),
      evidence = factor(evidence, levels = evidence_levels)
    )

  subtitle_txt <- str_wrap(
    "Black points and 95% CI bars = GA-adjusted primary model. Colored points = parsimonious sensitivity models. Point shape denotes evidence level. Arrows indicate OR or CI truncated for display.",
    width = 95
  )

  p <- ggplot() +
    geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.35, color = "grey30") +
    geom_segment(data = primary_dat, aes(y = PolymerLabel, yend = PolymerLabel, x = xmin, xend = xmax), linewidth = 0.8, color = "black") +
    geom_point(data = primary_dat, aes(y = PolymerLabel, x = x_plot, shape = evidence), size = 3.0, color = "black", stroke = 0.9) +
    geom_point(data = sens_dat, aes(y = PolymerLabel, x = x_plot, color = ModelLabel, shape = evidence), size = 2.7, stroke = 0.9) +
    geom_text(data = primary_dat %>% filter(left_trunc), aes(y = PolymerLabel, x = lims[1], label = "\u2190"), hjust = 1.05, size = 3.5) +
    geom_text(data = primary_dat %>% filter(right_trunc), aes(y = PolymerLabel, x = lims[2], label = "\u2192"), hjust = -0.05, size = 3.5) +
    scale_x_log10(breaks = breaks, limits = lims, labels = label_number(accuracy = 0.001, trim = TRUE)) +
    scale_color_manual(values = c(
      "Crude" = "#56B4E9",
      "GA + hospital" = "#E69F00",
      "GA + SDI" = "#009E73",
      "GA + prior PTB" = "#CC79A7",
      "GA + hypertension" = "#7A3E9D"
    )) +
    scale_shape_manual(
      values = c(
        "not nominal" = 16,
        "nominal p < 0.05 only" = 1,
        "BH q < 0.05, prespecified PVC/N66" = 23,
        "BH q < 0.05, all analytes" = 17
      )
    ) +
    labs(
      title = "Labor outcome model sensitivity analysis",
      subtitle = subtitle_txt,
      x = "OR per doubling (95% CI)",
      y = NULL,
      color = "Sensitivity models",
      shape = "Evidence"
    ) +
    guides(
      color = guide_legend(order = 1, override.aes = list(shape = 16, size = 3.0)),
      shape = guide_legend(order = 2, override.aes = list(color = "black", size = 3.0))
    ) +
    theme_mnp_forest(base_size = BASE_SIZE) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(size = 10.5),
      legend.position = "right",
      legend.title = element_text(face = "bold")
    )

  ggsave_mnp(out_path, p, width = 9.2, height = 6.7, dpi = 300)
}

plot_counts <- function(df, out_path) {
  counts_plot <- df %>%
    mutate(
      labor_outcome_plot = case_when(
        spontaneous_labor_augmented01 == 1 ~ "Spontaneous: augmented",
        spontaneous_labor_augmented01 == 0 ~ "Spontaneous",
        str_detect(Type_of_Labor_clean, "^no labor$|^none$|^no$") ~ "No labor",
        TRUE ~ "Missing/other"
      ),
      labor_outcome_plot = factor(
        labor_outcome_plot,
        levels = c("Spontaneous", "Spontaneous: augmented", "No labor", "Missing/other")
      )
    ) %>%
    filter(!is.na(Delivery_class)) %>%
    count(Delivery_class, labor_outcome_plot, name = "n") %>%
    group_by(Delivery_class) %>%
    mutate(total = sum(n), Delivery_label = paste0(as.character(Delivery_class), " (n=", total, ")")) %>%
    ungroup()

  p <- ggplot(counts_plot, aes(x = Delivery_label, y = n, fill = labor_outcome_plot)) +
    geom_col(position = "stack", width = 0.75, color = "black", linewidth = 0.2) +
    geom_text(aes(label = ifelse(n > 0, n, "")), position = position_stack(vjust = 0.5), size = 3.2) +
    scale_fill_manual(values = c(
      "Spontaneous" = "#F8766D",
      "Spontaneous: augmented" = "#7CAE00",
      "No labor" = "#00BFC4",
      "Missing/other" = "#C77CFF"
    )) +
    labs(
      title = "Type of labor by delivery class",
      subtitle = "Labor model outcome: Spontaneous: augmented vs Spontaneous",
      x = NULL,
      y = "Samples",
      fill = "Type of labor"
    ) +
    theme_mnp_forest(base_size = BASE_SIZE) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold"),
      legend.position = "right",
      legend.title = element_text(face = "bold")
    )

  ggsave_mnp(out_path, p, width = 8.2, height = 5.9, dpi = 300)
}

# ----------------------------
# Load and derive outcome
# ----------------------------
df0 <- read_model_table(IN_CSV)
df <- derive_labor_outcome(df0)

derive_audit <- df %>%
  count(Delivery_class, Type_of_Labor_raw, Type_of_Labor_clean, labor_outcome_label, name = "n") %>%
  arrange(Delivery_class, desc(n))

write_csv_clean(derive_audit, OUT_DERIVE)

# Candidate exposures
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

# ----------------------------
# Fit models
# ----------------------------
all_results <- purrr::imap_dfr(MODEL_SPECS, function(spec, tag) {
  purrr::pmap_dfr(
    eligible_exposures,
    function(Polymer, ExposureCol, UgCol, n_detected_gt0, continuous_model_ok, skip_reason) {
      fit_labor_model(
        df = df,
        exposure_col = ExposureCol,
        polymer = Polymer,
        model_tag = tag,
        model_label = spec$label,
        covariates = spec$covariates,
        note = spec$note
      )
    }
  )
}) %>%
  add_multiple_testing_columns()

write_csv_clean(all_results, OUT_ALL)

primary_results <- all_results %>%
  filter(ModelTag == PRIMARY_MODEL_TAG) %>%
  arrange(desc(Polymer == "Cumulative_MNPs"), OR_2x)

write_csv_clean(primary_results, OUT_PRIMARY)

q_audit <- all_results %>%
  filter(ModelTag == PRIMARY_MODEL_TAG) %>%
  select(
    ModelTag, ModelLabel, Polymer, PolymerLabel, polymer_family,
    p, q_BH_all_analytes, q_BH_polymers_only, q_BH_within_family,
    q_BH_labor_prespec, nominal_p05, sig_q05_all_analytes,
    sig_q05_polymers_only, sig_q05_within_family, sig_q05_labor_prespec,
    is_labor_prespecified, fit_status
  ) %>%
  arrange(p)

write_csv_clean(q_audit, OUT_Q_AUDIT)
save_gt_html_png(make_gt_table(primary_results), OUT_PRIMARY_HTML, OUT_PRIMARY_PNG)

model_key <- purrr::imap_dfr(MODEL_SPECS, function(spec, tag) {
  tibble(
    ModelTag = tag,
    ModelLabel = spec$label,
    CovariatesRequested = paste(spec$covariates, collapse = " + "),
    Note = spec$note
  )
})

write_csv_clean(model_key, OUT_KEY)

# ----------------------------
# Figures
# ----------------------------
plot_counts(df, OUT_FIG_COUNTS)

plot_primary_forest(primary_results, OUT_FIG_PRIMARY, full_axis = FALSE)
plot_sensitivity_forest(all_results, primary_results, OUT_FIG_SENS, full_axis = FALSE)

plot_primary_forest(primary_results, OUT_FIG_PRIMARY_FULL_AXIS, full_axis = TRUE)
plot_sensitivity_forest(all_results, primary_results, OUT_FIG_SENS_FULL_AXIS, full_axis = TRUE)

# ----------------------------
# Legend text
# ----------------------------
legend_text <- paste0(
  "Supplementary Figure | Polymer-specific associations with spontaneous labor requiring augmentation.\n",
  "Polymer-specific associations with spontaneous labor requiring augmentation are shown as odds ratios (ORs) per doubling of placental polymer burden with 95% confidence intervals. ",
  "The outcome was derived strictly from the Type_of_Labor column, comparing 'Spontaneous: augmented' with 'Spontaneous' and excluding 'No labor' or missing/other labor categories. ",
  "Because augmented labor events were limited, labor analyses were considered exploratory and used parsimonious adjustment. ",
  "The primary labor model adjusted only for gestational age at delivery. Point shape denotes BH q < 0.05. ",
  "Sensitivity models included crude models and one-covariate-at-a-time models adding hospital, Social Deprivation Index (SDI), prior preterm birth, or hypertension to the gestational age-adjusted model. ",
  "Multiple-testing correction is reported using explicitly labeled families: all analytes, polymers only, polymer family, and a pre-specified labor-hypothesis family. ",
  "Point shape distinguishes nominal p < 0.05, BH q < 0.05 within the pre-specified PVC/N66 labor-hypothesis set, and BH q < 0.05 across all modeled analytes. The x-axis is shown on a log scale; arrows indicate odds ratios or confidence intervals truncated for display."
)

writeLines(legend_text, OUT_LEGEND)

# ----------------------------
# Log
# ----------------------------
run_log <- c(
  "08_labor_outcome_supplement.R",
  paste0("PROJECT_ROOT: ", PROJECT_ROOT),
  paste0("Rows in input table: ", nrow(df0)),
  paste0("Labor pre-specified polymers: ", paste(LABOR_PRESPEC_POLYMERS, collapse = ", ")),
  "",
  "Labor outcome derivation:",
  capture.output(print(derive_audit, n = Inf)),
  "",
  "Primary GA-only model results:",
  capture.output(print(primary_results %>% select(Polymer, OR_2x, OR_2x_lo, OR_2x_hi, p, q_BH_all_analytes, q_BH_polymers_only, q_BH_within_family, q_BH_labor_prespec, N, Events_Augmented, NonEvents_SpontaneousNoAug, fit_status), n = Inf)),
  "",
  "Multiple-testing audit:",
  capture.output(print(q_audit, n = Inf)),
  "",
  "Model key:",
  capture.output(print(model_key, n = Inf)),
  "",
  "Wrote:",
  OUT_DERIVE,
  OUT_ALL,
  OUT_PRIMARY,
  OUT_Q_AUDIT,
  OUT_SPARSE,
  OUT_KEY,
  OUT_FIG_COUNTS,
  OUT_FIG_PRIMARY,
  OUT_FIG_SENS,
  OUT_FIG_PRIMARY_FULL_AXIS,
  OUT_FIG_SENS_FULL_AXIS,
  OUT_LEGEND
)

cat(paste(run_log, collapse = "\n"), file = OUT_LOG)
message(paste(run_log, collapse = "\n"))
