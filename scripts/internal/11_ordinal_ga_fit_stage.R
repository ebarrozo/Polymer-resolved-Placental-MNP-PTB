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
# 11_ordinal_ga_category_models.R
# Ordinal gestational-age category models for MNP burden
# v4.4: self-contained helper + recursion fix
# ============================================================
#
# Run after run_analysis_setup.R has already created:
#   outputs/tables/mnp_ptb_cohort__analysis_ready_model_table.csv
#
# You do NOT need an in-memory df.
# This script reads the analysis-ready CSV from disk.
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

IN_MODEL_TABLE <- file.path(DIR_TABLES, "mnp_ptb_cohort__analysis_ready_model_table.csv")

STAMP <- format(Sys.time(), "%Y%m%d_%H%M%S")
OUT_DIR <- file.path(DIR_TABLES, paste0("11_ordinal_ga_category_models__", STAMP))
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

OUT_ANALYSIS_TABLE <- file.path(OUT_DIR, "ordinal_ga_category_analysis_table.csv")
OUT_CATEGORY_COUNTS <- file.path(OUT_DIR, "ordinal_ga_category_counts.csv")
OUT_ALL <- file.path(OUT_DIR, "ordinal_ga_category_all_models.csv")
OUT_PRIMARY <- file.path(OUT_DIR, "ordinal_ga_category_primary_models.csv")
OUT_MODEL_KEY <- file.path(OUT_DIR, "ordinal_ga_category_model_key.csv")
OUT_POLYMER_AUDIT <- file.path(OUT_DIR, "ordinal_ga_category_polymer_audit.csv")
OUT_ASSUMPTION_AUDIT <- file.path(OUT_DIR, "ordinal_ga_category_assumption_audit.csv")
OUT_LEGEND <- file.path(DIR_TABLES, "ordinal_ga_category_figure_legend.txt")
OUT_LOG <- file.path(DIR_LOGS, "11_ordinal_ga_category_models__run_log.txt")

OUT_FIG_COUNTS <- file.path(DIR_FIGS, "fig_ordinal_ga_category__counts.png")
OUT_FIG_PRIMARY <- file.path(DIR_FIGS, "fig_ordinal_ga_category__primary_forest_earlier_delivery.png")
OUT_FIG_SENS <- file.path(DIR_FIGS, "fig_ordinal_ga_category__sensitivity_forest_earlier_delivery.png")
OUT_FIG_PRIMARY_FULL <- file.path(DIR_FIGS, "fig_ordinal_ga_category__primary_forest_earlier_delivery_full_axis.png")
OUT_FIG_SENS_FULL <- file.path(DIR_FIGS, "fig_ordinal_ga_category__sensitivity_forest_earlier_delivery_full_axis.png")

BASE_SIZE <- 11
K2 <- log10(2)

DISPLAY_XMIN <- 0.03
DISPLAY_XMAX <- 10
DISPLAY_BREAKS <- c(0.03, 0.05, 0.10, 0.20, 0.50, 1, 2, 4, 10)

GA_LEVELS <- c(
  "Very/extremely preterm (<32w)",
  "Moderate preterm (32-<34w)",
  "Late preterm (34-<37w)",
  "Early term (37-<39w)",
  "Full term (39-41w)",
  "Late term (>41w)"
)

POLYMERS_LOCAL <- c("PE", "SBR", "PP", "PVC", "N66", "PET", "N6", "PMMA", "PC", "PU", "ABS", "PS")
ANALYTES <- c("Cumulative_MNPs", POLYMERS_LOCAL)

POLYMER_FULL_NAMES_LOCAL <- c(
  PE   = "Polyethylene",
  SBR  = "Styrene-butadiene rubber",
  PP   = "Polypropylene",
  PVC  = "Polyvinyl chloride",
  N66  = "Nylon-66",
  PET  = "Polyethylene terephthalate",
  N6   = "Nylon-6",
  PMMA = "Polymethyl methacrylate",
  PC   = "Polycarbonate",
  PU   = "Polyurethane",
  ABS  = "Acrylonitrile-butadiene-styrene",
  PS   = "Polystyrene"
)

MODEL_SPECS <- list(
  CRUDE = list(
    label = "Crude",
    covariates = character(0)
  ),
  M1_PRIMARY = list(
    label = "Primary adjusted",
    covariates = c("Maternal_age", "BMI", "SDI_score", "Twins01", "Prior_PTB", "Hypertension01", "Hospital")
  ),
  M1_NO_SDI = list(
    label = "Primary without SDI",
    covariates = c("Maternal_age", "BMI", "Twins01", "Prior_PTB", "Hypertension01", "Hospital")
  ),
  M1_PLUS_SVI = list(
    label = "Primary + SVI",
    covariates = c("Maternal_age", "BMI", "SDI_score", "SVI_score", "Twins01", "Prior_PTB", "Hypertension01", "Hospital")
  ),
  M1_PLUS_RACE = list(
    label = "Primary + race/ethnicity",
    covariates = c("Maternal_age", "BMI", "SDI_score", "Twins01", "Prior_PTB", "Hypertension01", "Hospital", "Ethnicity_Race")
  )
)

PRIMARY_MODEL_TAG <- "M1_PRIMARY"

# ----------------------------
# Self-contained helpers
# ----------------------------
safe_num2 <- function(x) suppressWarnings(as.numeric(as.character(x)))

safe_chr2 <- function(x) {
  if (missing(x) || is.null(x)) return(rep(NA_character_, 0))
  y <- as.character(x)
  y <- stringr::str_squish(y)
  y[y %in% c("", "NA", "NaN", "N/A", "NULL", "null", "Unknown", "unknown")] <- NA_character_
  y
}

as_01_local <- function(x) {
  if (missing(x) || is.null(x)) return(rep(NA_integer_, 0))
  if (is.logical(x)) return(as.integer(x))
  if (is.numeric(x) || is.integer(x)) {
    y <- as.numeric(x)
    return(dplyr::case_when(
      is.na(y) ~ NA_integer_,
      y == 0 ~ 0L,
      y == 1 ~ 1L,
      TRUE ~ as.integer(y > 0)
    ))
  }
  z <- stringr::str_to_lower(safe_chr2(x))
  dplyr::case_when(
    is.na(z) ~ NA_integer_,
    z %in% c("1", "yes", "y", "true", "t", "present", "positive", "pos", "current", "preterm") ~ 1L,
    z %in% c("0", "no", "n", "false", "f", "absent", "negative", "neg", "never", "term") ~ 0L,
    TRUE ~ NA_integer_
  )
}

write_csv_clean_local <- function(x, file, na = "") {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(x, file = file, na = na)
  invisible(file)
}

fmt_est_ci_local <- function(est, lo, hi, digits = 2) {
  ifelse(
    is.na(est),
    NA_character_,
    sprintf(
      paste0("%.", digits, "f (%.", digits, "f-%.", digits, "f)"),
      est, lo, hi
    )
  )
}

fmt_log_axis_local <- function(x) {
  dplyr::case_when(
    x < 0.01 ~ formatC(x, format = "f", digits = 3),
    x < 1 ~ formatC(x, format = "f", digits = 2),
    x < 10 ~ formatC(x, format = "f", digits = 1),
    TRUE ~ formatC(x, format = "f", digits = 0)
  )
}

legacy_col_ug_g_local <- function(polymer) {
  ifelse(polymer == "Cumulative_MNPs", "Cumulative_MNPs_ug_per_g", paste0(polymer, "_ug_per_g"))
}

legacy_col_log10_local <- function(polymer) {
  ifelse(polymer == "Cumulative_MNPs", "Cumulative_MNPs_log10ugpg", paste0(polymer, "_log10ugpg"))
}

poly_label_local <- function(poly) {
  poly <- as.character(poly)
  dplyr::case_when(
    poly == "Cumulative_MNPs" ~ "Total MNP burden",
    poly %in% names(POLYMER_FULL_NAMES_LOCAL) ~ paste0(POLYMER_FULL_NAMES_LOCAL[poly], " (", poly, ")"),
    TRUE ~ poly
  )
}

poly_family_local <- function(poly) {
  poly <- as.character(poly)
  dplyr::case_when(
    poly == "Cumulative_MNPs" ~ "cumulative",
    poly %in% c("PE", "PP") ~ "polyolefins",
    poly %in% c("PET") ~ "polyesters",
    poly %in% c("PVC") ~ "vinyls",
    poly %in% c("PU") ~ "polyurethanes",
    poly %in% c("PC") ~ "polycarbonates",
    poly %in% c("PMMA") ~ "acrylics",
    poly %in% c("N6", "N66") ~ "polyamides",
    poly %in% c("SBR", "ABS", "PS") ~ "styrenics",
    TRUE ~ "other"
  )
}

detected_n_local <- function(x) {
  x <- safe_num2(x)
  sum(is.finite(x) & x > 0, na.rm = TRUE)
}

usable_covariates_local <- function(dat, covs) {
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

derive_ga_category <- function(ga_weeks) {
  ga <- safe_num2(ga_weeks)
  out <- case_when(
    is.na(ga) ~ NA_character_,
    ga < 32 ~ "Very/extremely preterm (<32w)",
    ga >= 32 & ga < 34 ~ "Moderate preterm (32-<34w)",
    ga >= 34 & ga < 37 ~ "Late preterm (34-<37w)",
    ga >= 37 & ga < 39 ~ "Early term (37-<39w)",
    ga >= 39 & ga <= 41 ~ "Full term (39-41w)",
    ga > 41 ~ "Late term (>41w)",
    TRUE ~ NA_character_
  )
  factor(out, levels = GA_LEVELS, ordered = TRUE)
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
      GA = safe_num2(GA),
      Maternal_age = safe_num2(Maternal_age),
      BMI = safe_num2(BMI),
      SDI_score = safe_num2(SDI_score),
      SVI_score = safe_num2(SVI_score),
      Twins01 = as.integer(as_01_local(Twins01)),
      Prior_PTB = as.integer(as_01_local(Prior_PTB)),
      Hypertension01 = as.integer(as_01_local(Hypertension01)),
      Hospital = factor(safe_chr2(Hospital)),
      Ethnicity_Race = factor(safe_chr2(Ethnicity_Race))
    )
}

make_axis_limits_local <- function(dat, xmin_floor = 0.03, xmax_cap = 100) {
  vals <- c(dat$OR_earlier_2x, dat$OR_earlier_2x_lo, dat$OR_earlier_2x_hi)
  vals <- vals[is.finite(vals) & vals > 0]
  if (length(vals) == 0) return(c(xmin_floor, 10))
  lo <- 10 ^ floor(log10(min(vals, na.rm = TRUE)))
  hi <- 10 ^ ceiling(log10(max(vals, na.rm = TRUE)))
  lo <- max(lo, xmin_floor)
  hi <- min(hi, xmax_cap)
  if (!is.finite(lo) || !is.finite(hi) || lo >= hi) c(xmin_floor, 10) else c(lo, hi)
}

log_breaks_local <- function(lims) {
  br <- c(0.001, 0.003, 0.01, 0.03, 0.05, 0.10, 0.20, 0.50, 1, 2, 4, 10, 20, 50, 100)
  br[br >= lims[1] & br <= lims[2]]
}

# ----------------------------
# Model function
# ----------------------------
fit_ordinal_model <- function(df, exposure_col, polymer, model_tag, model_label, covariates, exclude_outliers = FALSE) {
  d0 <- df

  if (exclude_outliers && "cluster_clr_pca_outlier" %in% names(d0)) {
    d0 <- d0 %>% filter(is.na(cluster_clr_pca_outlier) | !as.logical(cluster_clr_pca_outlier))
  }

  covs <- usable_covariates_local(d0, covariates)
  vars <- c("SampleID", "PregnancyID", "GA_category_ordered", exposure_col, covs)

  d <- d0 %>%
    dplyr::select(all_of(vars)) %>%
    mutate(
      outcome = droplevels(GA_category_ordered),
      exposure = safe_num2(.data[[exposure_col]])
    ) %>%
    dplyr::select(-all_of(c("GA_category_ordered", exposure_col))) %>%
    drop_na() %>%
    mutate(across(where(is.factor), droplevels))

  n <- nrow(d)
  n_categories <- nlevels(d$outcome)
  category_counts <- paste(names(table(d$outcome)), as.integer(table(d$outcome)), sep = "=", collapse = "; ")

  base <- tibble(
    ModelTag = model_tag,
    ModelLabel = model_label,
    Polymer = polymer,
    PolymerLabel = poly_label_local(polymer),
    PolymerFamily = poly_family_local(polymer),
    ExposureCol = exposure_col,
    N = n,
    N_categories = n_categories,
    CategoryCounts = category_counts,
    Covariates = paste(covs, collapse = " + "),
    OR_later_2x = NA_real_,
    OR_later_2x_lo = NA_real_,
    OR_later_2x_hi = NA_real_,
    OR_earlier_2x = NA_real_,
    OR_earlier_2x_lo = NA_real_,
    OR_earlier_2x_hi = NA_real_,
    p = NA_real_,
    fit_engine = NA_character_,
    fit_status = NA_character_
  )

  if (n < 30 || n_categories < 3 || sd(d$exposure, na.rm = TRUE) == 0) {
    return(base %>% mutate(fit_status = "skipped_insufficient_n_categories_or_exposure_variance"))
  }

  rhs <- c("exposure", covs)
  fml <- as.formula(paste("outcome ~", paste(rhs, collapse = " + ")))

  mm <- try(model.matrix(as.formula(paste("~", paste(rhs, collapse = " + "))), data = d), silent = TRUE)
  if (!inherits(mm, "try-error") && qr(mm)$rank >= n - 3) {
    return(base %>% mutate(fit_status = "skipped_rank_too_high_for_n"))
  }

  if (requireNamespace("ordinal", quietly = TRUE)) {
    fit <- try(ordinal::clm(fml, data = d, link = "logit", Hess = TRUE), silent = TRUE)

    if (!inherits(fit, "try-error")) {
      sm <- try(summary(fit)$coefficients, silent = TRUE)
      if (!inherits(sm, "try-error") && "exposure" %in% rownames(sm)) {
        b <- unname(sm["exposure", "Estimate"])
        se <- unname(sm["exposure", "Std. Error"])
        pval <- unname(sm["exposure", "Pr(>|z|)"])

        log_or_later <- b * K2
        se_log_or_later <- se * K2

        or_later <- exp(log_or_later)
        or_later_lo <- exp(log_or_later - 1.96 * se_log_or_later)
        or_later_hi <- exp(log_or_later + 1.96 * se_log_or_later)

        return(base %>%
          mutate(
            OR_later_2x = or_later,
            OR_later_2x_lo = or_later_lo,
            OR_later_2x_hi = or_later_hi,
            OR_earlier_2x = 1 / or_later,
            OR_earlier_2x_lo = 1 / or_later_hi,
            OR_earlier_2x_hi = 1 / or_later_lo,
            p = pval,
            fit_engine = "ordinal::clm",
            fit_status = "ok"
          ))
      }
    }
  }

  if (requireNamespace("MASS", quietly = TRUE)) {
    fit2 <- try(MASS::polr(fml, data = d, method = "logistic", Hess = TRUE), silent = TRUE)

    if (!inherits(fit2, "try-error")) {
      cf <- try(coef(summary(fit2)), silent = TRUE)
      if (!inherits(cf, "try-error") && "exposure" %in% rownames(cf)) {
        b <- unname(cf["exposure", "Value"])
        se <- unname(cf["exposure", "Std. Error"])
        z <- b / se
        pval <- 2 * pnorm(abs(z), lower.tail = FALSE)

        log_or_later <- b * K2
        se_log_or_later <- se * K2

        or_later <- exp(log_or_later)
        or_later_lo <- exp(log_or_later - 1.96 * se_log_or_later)
        or_later_hi <- exp(log_or_later + 1.96 * se_log_or_later)

        return(base %>%
          mutate(
            OR_later_2x = or_later,
            OR_later_2x_lo = or_later_lo,
            OR_later_2x_hi = or_later_hi,
            OR_earlier_2x = 1 / or_later,
            OR_earlier_2x_lo = 1 / or_later_hi,
            OR_earlier_2x_hi = 1 / or_later_lo,
            p = pval,
            fit_engine = "MASS::polr",
            fit_status = "ok"
          ))
      }
    }

    return(base %>% mutate(fit_status = paste0("fit_error: ", as.character(fit2)[1])))
  }

  base %>% mutate(fit_status = "skipped_no_ordinal_or_MASS_package")
}

add_q_values <- function(res) {
  res %>%
    group_by(ModelTag) %>%
    mutate(
      q_BH_all_analytes = {
        q <- rep(NA_real_, n())
        idx <- fit_status == "ok" & !is.na(p)
        q[idx] <- p.adjust(p[idx], method = "BH")
        q
      },
      q_BH_polymers_only = {
        q <- rep(NA_real_, n())
        idx <- fit_status == "ok" & !is.na(p) & Polymer != "Cumulative_MNPs"
        if (sum(idx) > 0) q[idx] <- p.adjust(p[idx], method = "BH")
        q
      }
    ) %>%
    ungroup() %>%
    group_by(ModelTag, PolymerFamily) %>%
    mutate(
      q_BH_within_family = {
        q <- rep(NA_real_, n())
        idx <- fit_status == "ok" & !is.na(p) & Polymer != "Cumulative_MNPs"
        if (sum(idx) > 0) q[idx] <- p.adjust(p[idx], method = "BH")
        q
      }
    ) %>%
    ungroup() %>%
    mutate(
      q_BH = q_BH_all_analytes,
      sig_q05 = q_BH_all_analytes < 0.05,
      nominal_p05 = p < 0.05,
      OR_earlier_CI = fmt_est_ci_local(OR_earlier_2x, OR_earlier_2x_lo, OR_earlier_2x_hi)
    )
}

# ----------------------------
# Plot helpers
# ----------------------------
prep_forest_data <- function(tbl, order_levels = NULL, lims = c(DISPLAY_XMIN, DISPLAY_XMAX)) {
  dat <- tbl %>%
    filter(fit_status == "ok", !is.na(OR_earlier_2x)) %>%
    mutate(
      xmin = pmax(OR_earlier_2x_lo, lims[1]),
      xmax = pmin(OR_earlier_2x_hi, lims[2]),
      x_plot = pmin(pmax(OR_earlier_2x, lims[1]), lims[2]),
      left_trunc = OR_earlier_2x_lo < lims[1] | OR_earlier_2x < lims[1],
      right_trunc = OR_earlier_2x_hi > lims[2] | OR_earlier_2x > lims[2]
    )

  if (is.null(order_levels)) {
    order_levels <- dat %>% arrange(desc(OR_earlier_2x)) %>% pull(PolymerLabel)
  }

  dat %>%
    mutate(
      PolymerLabel = factor(PolymerLabel, levels = rev(order_levels)),
      sig_q05 = factor(ifelse(sig_q05, "TRUE", "FALSE"), levels = c("FALSE", "TRUE"))
    )
}

plot_primary <- function(primary_results, out_path, full_axis = FALSE) {
  lims <- if (full_axis) make_axis_limits_local(primary_results, xmin_floor = 0.001, xmax_cap = 100) else c(DISPLAY_XMIN, DISPLAY_XMAX)
  breaks <- if (full_axis) log_breaks_local(lims) else DISPLAY_BREAKS

  plot_dat <- prep_forest_data(primary_results, lims = lims)

  subtitle_txt <- str_wrap(
    "Cumulative-link ordinal regression across gestational-age categories. OR >1 indicates greater odds of an earlier gestational-age category per doubling of placental MNP burden.",
    width = 95
  )

  p <- ggplot(plot_dat, aes(y = PolymerLabel)) +
    geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.35, color = "grey30") +
    geom_segment(aes(x = xmin, xend = xmax, yend = PolymerLabel), color = "black", linewidth = 0.85) +
    geom_point(aes(x = x_plot, shape = sig_q05), color = "black", size = 3.2) +
    geom_text(data = plot_dat %>% filter(left_trunc), aes(x = lims[1], label = "\u2190"), hjust = 1.05, size = 3.5) +
    geom_text(data = plot_dat %>% filter(right_trunc), aes(x = lims[2], label = "\u2192"), hjust = -0.05, size = 3.5) +
    scale_x_log10(breaks = breaks, limits = lims, labels = fmt_log_axis_local) +
    scale_shape_manual(values = c("FALSE" = 16, "TRUE" = 17), name = "BH q < 0.05") +
    labs(
      title = "Ordinal gestational-age category models",
      subtitle = subtitle_txt,
      x = "OR for earlier gestational-age category per doubling (95% CI)",
      y = NULL
    ) +
    theme_bw(base_size = BASE_SIZE) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(size = 10.2),
      legend.position = "right",
      legend.title = element_text(face = "bold")
    )

  ggsave(out_path, p, width = 8.8, height = 6.4, dpi = 300)
}

plot_sensitivity <- function(all_results, primary_results, out_path, full_axis = FALSE) {
  sens <- all_results %>% filter(ModelTag != PRIMARY_MODEL_TAG)
  axis_dat <- bind_rows(primary_results, sens)

  lims <- if (full_axis) make_axis_limits_local(axis_dat, xmin_floor = 0.001, xmax_cap = 100) else c(DISPLAY_XMIN, DISPLAY_XMAX)
  breaks <- if (full_axis) log_breaks_local(lims) else DISPLAY_BREAKS

  order_levels <- primary_results %>%
    filter(fit_status == "ok") %>%
    arrange(desc(OR_earlier_2x)) %>%
    pull(PolymerLabel)

  primary_dat <- prep_forest_data(primary_results, order_levels, lims)
  sens_dat <- prep_forest_data(sens, order_levels, lims) %>%
    mutate(
      ModelLabel = factor(ModelLabel, levels = c(
        "Crude",
        "Primary without SDI",
        "Primary + SVI",
        "Primary + race/ethnicity",
        "Primary adjusted, excluding CLR/PCA outliers"
      ))
    )

  subtitle_txt <- str_wrap(
    "Black points and 95% CI bars = primary adjusted model. Colored points = sensitivity models. Point shape denotes BH q<0.05.",
    width = 95
  )

  p <- ggplot() +
    geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.35, color = "grey30") +
    geom_segment(data = primary_dat, aes(y = PolymerLabel, yend = PolymerLabel, x = xmin, xend = xmax), color = "black", linewidth = 0.85) +
    geom_point(data = primary_dat, aes(y = PolymerLabel, x = x_plot, shape = sig_q05), color = "black", size = 3.0) +
    geom_point(data = sens_dat, aes(y = PolymerLabel, x = x_plot, color = ModelLabel, shape = sig_q05), size = 2.5, stroke = 0.8) +
    geom_text(data = primary_dat %>% filter(left_trunc), aes(y = PolymerLabel, x = lims[1], label = "\u2190"), hjust = 1.05, size = 3.5) +
    geom_text(data = primary_dat %>% filter(right_trunc), aes(y = PolymerLabel, x = lims[2], label = "\u2192"), hjust = -0.05, size = 3.5) +
    scale_x_log10(breaks = breaks, limits = lims, labels = fmt_log_axis_local) +
    scale_shape_manual(values = c("FALSE" = 16, "TRUE" = 17), name = "BH q < 0.05") +
    scale_color_manual(
      values = c(
        "Crude" = "#56B4E9",
        "Primary without SDI" = "#7A3E9D",
        "Primary + SVI" = "#009E73",
        "Primary + race/ethnicity" = "#CC79A7",
        "Primary adjusted, excluding CLR/PCA outliers" = "#E69F00"
      ),
      drop = FALSE
    ) +
    labs(
      title = "Ordinal gestational-age category model sensitivity analysis",
      subtitle = subtitle_txt,
      x = "OR for earlier gestational-age category per doubling (95% CI)",
      y = NULL,
      color = "Sensitivity models"
    ) +
    guides(
      color = guide_legend(order = 1, override.aes = list(shape = 16, size = 3)),
      shape = guide_legend(order = 2, override.aes = list(color = "black", size = 3))
    ) +
    theme_bw(base_size = BASE_SIZE) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(size = 10.2),
      legend.position = "right",
      legend.title = element_text(face = "bold")
    )

  ggsave(out_path, p, width = 9.3, height = 6.6, dpi = 300)
}

plot_counts <- function(dat, out_path) {
  counts <- dat %>%
    count(GA_category_ordered, name = "n") %>%
    mutate(GA_category_ordered = factor(GA_category_ordered, levels = GA_LEVELS, ordered = TRUE))

  p <- ggplot(counts, aes(x = GA_category_ordered, y = n)) +
    geom_col(fill = "grey70", color = "black", linewidth = 0.25, width = 0.75) +
    geom_text(aes(label = n), vjust = -0.25, size = 3.3) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.14))) +
    labs(
      title = "Gestational-age category distribution",
      x = NULL,
      y = "Samples"
    ) +
    theme_bw(base_size = BASE_SIZE) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold"),
      axis.text.x = element_text(angle = 35, hjust = 1)
    )

  ggsave(out_path, p, width = 8.2, height = 4.6, dpi = 300)
}

# ----------------------------
# Load data from disk
# ----------------------------
if (!file.exists(IN_MODEL_TABLE)) {
  stop("Cannot find analysis-ready model table: ", IN_MODEL_TABLE, call. = FALSE)
}

df0 <- read_model_table(IN_MODEL_TABLE)

df <- df0 %>%
  mutate(GA_category_ordered = derive_ga_category(GA)) %>%
  filter(!is.na(GA_category_ordered)) %>%
  mutate(GA_category_ordered = droplevels(GA_category_ordered))

# Detect optional outlier flag if already merged into the table.
outlier_col <- names(df)[stringr::str_detect(names(df), regex("clr.*pca.*outlier|cluster.*outlier", ignore_case = TRUE))]
if (length(outlier_col) > 0 && !"cluster_clr_pca_outlier" %in% names(df)) {
  df <- df %>% mutate(cluster_clr_pca_outlier = as.logical(.data[[outlier_col[1]]]))
}

write_csv_clean_local(df, OUT_ANALYSIS_TABLE)

category_counts <- df %>%
  count(GA_category_ordered, name = "n") %>%
  mutate(prop = n / sum(n))

write_csv_clean_local(category_counts, OUT_CATEGORY_COUNTS)

# ----------------------------
# Exposure eligibility
# ----------------------------
exposure_tbl <- tibble(
  Polymer = ANALYTES,
  ExposureCol = vapply(ANALYTES, legacy_col_log10_local, character(1)),
  UgCol = vapply(ANALYTES, legacy_col_ug_g_local, character(1))
) %>%
  mutate(
    exposure_available = ExposureCol %in% names(df),
    n_detected_gt0 = vapply(UgCol, function(cn) if (cn %in% names(df)) detected_n_local(df[[cn]]) else 0, numeric(1)),
    continuous_model_ok = Polymer == "Cumulative_MNPs" | n_detected_gt0 >= MIN_DETECTED_N_FOR_CONTINUOUS_MODEL,
    skip_reason = case_when(
      !exposure_available ~ "missing_exposure_column",
      !continuous_model_ok ~ "sparse_detection",
      TRUE ~ NA_character_
    )
  )

write_csv_clean_local(exposure_tbl, OUT_POLYMER_AUDIT)
eligible_exposures <- exposure_tbl %>% filter(is.na(skip_reason))

model_specs <- MODEL_SPECS
if ("cluster_clr_pca_outlier" %in% names(df)) {
  model_specs$M1_NO_CLR_OUTLIERS <- list(
    label = "Primary adjusted, excluding CLR/PCA outliers",
    covariates = MODEL_SPECS$M1_PRIMARY$covariates,
    exclude_outliers = TRUE
  )
}

model_key <- purrr::imap_dfr(model_specs, function(spec, tag) {
  tibble(
    ModelTag = tag,
    ModelLabel = spec$label,
    CovariatesRequested = paste(spec$covariates, collapse = " + "),
    ExcludeCLRorPCAOutliers = isTRUE(spec$exclude_outliers)
  )
})
write_csv_clean_local(model_key, OUT_MODEL_KEY)

# ----------------------------
# Fit models
# ----------------------------
all_results <- purrr::imap_dfr(model_specs, function(spec, tag) {
  purrr::pmap_dfr(
    eligible_exposures,
    function(Polymer, ExposureCol, UgCol, exposure_available, n_detected_gt0, continuous_model_ok, skip_reason) {
      fit_ordinal_model(
        df = df,
        exposure_col = ExposureCol,
        polymer = Polymer,
        model_tag = tag,
        model_label = spec$label,
        covariates = spec$covariates,
        exclude_outliers = isTRUE(spec$exclude_outliers)
      )
    }
  )
}) %>%
  add_q_values()

write_csv_clean_local(all_results, OUT_ALL)

primary_results <- all_results %>%
  filter(ModelTag == PRIMARY_MODEL_TAG) %>%
  arrange(desc(OR_earlier_2x))

write_csv_clean_local(primary_results, OUT_PRIMARY)

assumption_audit <- tibble(
  engine_preferred = ifelse(requireNamespace("ordinal", quietly = TRUE), "ordinal::clm", "MASS::polr fallback"),
  proportional_odds_note = "Cumulative-link model assumes proportional odds. This script estimates ordinal associations but does not formally test all analyte-specific proportional-odds assumptions.",
  outcome_order = paste(GA_LEVELS, collapse = " < ")
)
write_csv_clean_local(assumption_audit, OUT_ASSUMPTION_AUDIT)

# ----------------------------
# Figures
# ----------------------------
plot_counts(df, OUT_FIG_COUNTS)
plot_primary(primary_results, OUT_FIG_PRIMARY, full_axis = FALSE)
plot_sensitivity(all_results, primary_results, OUT_FIG_SENS, full_axis = FALSE)
plot_primary(primary_results, OUT_FIG_PRIMARY_FULL, full_axis = TRUE)
plot_sensitivity(all_results, primary_results, OUT_FIG_SENS_FULL, full_axis = TRUE)

legend_text <- paste0(
  "Figure | Ordinal gestational-age category models. ",
  "Gestational age at delivery was modeled as an ordered categorical outcome: very/extremely preterm (<32 weeks), moderate preterm (32 to <34 weeks), late preterm (34 to <37 weeks), early term (37 to <39 weeks), full term (39 to 41 weeks), and late term (>41 weeks). ",
  "Cumulative-link ordinal regression was used to estimate odds ratios per doubling of placental MNP burden. ",
  "Estimates are plotted as odds ratios for being in an earlier gestational-age category, so values greater than 1 indicate earlier delivery. ",
  "The primary adjusted model included maternal age, BMI, Social Deprivation Index, twin gestation, prior preterm birth, hypertension, and hospital; gestational age was not included as a covariate because it defines the outcome. ",
  "Colored points represent sensitivity models. Point shape denotes BH q<0.05 across modeled analytes. ",
  "The x-axis is log-scaled; arrows indicate estimates or confidence intervals truncated for display."
)
writeLines(legend_text, OUT_LEGEND)

# ----------------------------
# Log
# ----------------------------
run_log <- c(
  "11_ordinal_ga_category_models.R v4.4",
  paste0("PROJECT_ROOT: ", PROJECT_ROOT),
  paste0("Rows in input table: ", nrow(df0)),
  paste0("Rows with GA category: ", nrow(df)),
  "",
  "Category counts:",
  capture.output(print(category_counts, n = Inf)),
  "",
  "Polymer audit:",
  capture.output(print(exposure_tbl, n = Inf)),
  "",
  "Model key:",
  capture.output(print(model_key, n = Inf)),
  "",
  "Primary model results:",
  capture.output(print(primary_results %>% dplyr::select(Polymer, OR_earlier_2x, OR_earlier_2x_lo, OR_earlier_2x_hi, p, q_BH_all_analytes, N, fit_engine, fit_status), n = Inf)),
  "",
  "Wrote:",
  OUT_ANALYSIS_TABLE,
  OUT_CATEGORY_COUNTS,
  OUT_ALL,
  OUT_PRIMARY,
  OUT_MODEL_KEY,
  OUT_POLYMER_AUDIT,
  OUT_ASSUMPTION_AUDIT,
  OUT_FIG_COUNTS,
  OUT_FIG_PRIMARY,
  OUT_FIG_SENS,
  OUT_FIG_PRIMARY_FULL,
  OUT_FIG_SENS_FULL,
  OUT_LEGEND
)

cat(paste(run_log, collapse = "\n"), file = OUT_LOG)
message(paste(run_log, collapse = "\n"))
