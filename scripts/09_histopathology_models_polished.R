#!/usr/bin/env Rscript
# POLISHED VERSION: sources R/figure_theme_mnp.R and writes *_polished outputs where applicable.
# ============================================================
# 09_histopathology_models.R
# v4.0: dplyr::select namespace fix plus histopath axis/legend polish
# Histopathology severity models joined to analysis-ready MNP table
# ============================================================
#
# Run after:
#   source("run_analysis_setup.R")
#
# Required input in PROJECT_ROOT or current working directory:
#   restricted pathology-score CSV
#
# Required existing pipeline input:
#   outputs/tables/mnp_ptb_cohort__analysis_ready_model_table.csv
#   outputs/tables/mnp_ptb_cohort__polymer_detection_summary.csv
#
# What this script does:
#   1. Reads pathology score file.
#   2. Drops empty "Unnamed" columns.
#   3. Standardizes PregnancyID.
#   4. Joins to the analysis-ready MNP model table by PregnancyID.
#   5. Subsets samples with path scores.
#   6. Runs stratified models by Delivery_class:
#        A. Microcalcification severity among positives:
#           NB model on microcalc_1 + microcalc_2 among microcalc-positive samples.
#           Effect = IRR per doubling in polymer burden.
#        B. Perivillous fibrin deposition severity among PVF-positive samples:
#           Gamma-log model on average_pvf among PVF-positive samples.
#           Effect = RR per doubling in polymer burden.
#   7. Generates PTB-style forest figures with:
#        Black CI bars = M1 clinical primary model
#        Colored points = sensitivity models
#        M1 point shape/fill = global BH q < 0.05
#        Thick M1 outline = family BH q < 0.05
#
# Notes:
#   - Models are exploratory and stratified by term/preterm.
#   - q_global is BH across modeled polymers within outcome/model/delivery class.
#   - q_family is BH within polymer family within outcome/model/delivery class.
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
  library(stringr)
  library(ggplot2)
  library(scales)
  library(MASS)
})

source("00_config.R")
source(file.path("R", "utils_mnp.R"))
source(file.path("R", "figure_theme_mnp.R"))

# ----------------------------
# File locations
# ----------------------------
PATH_FILE <- Sys.getenv(
  "MNP_PATHOLOGY_INPUT",
  unset = file.path(PROJECT_ROOT, "data", "restricted", "pathology_scores.csv")
)
if (!file.exists(PATH_FILE)) {
  stop(
    "Pathology input not found: ", PATH_FILE, "\n",
    "Set MNP_PATHOLOGY_INPUT to the restricted pathology-score CSV.",
    call. = FALSE
  )
}

IN_MODEL_TABLE <- file.path(DIR_TABLES, "mnp_ptb_cohort__analysis_ready_model_table.csv")
IN_DETECT <- file.path(DIR_TABLES, "mnp_ptb_cohort__polymer_detection_summary.csv")

STAMP <- format(Sys.time(), "%Y%m%d_%H%M%S")
OUT_DIR <- file.path(DIR_TABLES, paste0("09_histopathology_models__", STAMP))
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

OUT_JOINED <- file.path(OUT_DIR, "histopath_joined_analysis_table_polished.csv")
OUT_JOIN_AUDIT <- file.path(OUT_DIR, "histopath_join_audit_polished.csv")
OUT_MODEL_ALL <- file.path(OUT_DIR, "histopath_model_all_results_polished.csv")
OUT_MODEL_MICRO <- file.path(OUT_DIR, "histopath_microcalc_NB_positive_results_polished.csv")
OUT_MODEL_PVF <- file.path(OUT_DIR, "histopath_pvf_gamma_positive_results_polished.csv")
OUT_Q_AUDIT <- file.path(OUT_DIR, "histopath_multiple_testing_audit_polished.csv")
OUT_MODEL_KEY <- file.path(OUT_DIR, "histopath_model_key_polished.csv")
OUT_POLYMER_AUDIT <- file.path(OUT_DIR, "histopath_polymer_audit_polished.csv")
OUT_LEGEND <- file.path(DIR_TABLES, "supp_histopathology_models_figure_legend_polished.txt")
OUT_LOG <- file.path(DIR_LOGS, "09_histopathology_models__run_log_polished.txt")

OUT_FIG_MICRO <- file.path(DIR_FIGS, "supp_histopath__microcalcification_severity_NB_positive_polished.png")
OUT_FIG_PVF <- file.path(DIR_FIGS, "supp_histopath__pvf_severity_gamma_positive_polished.png")
OUT_FIG_MICRO_FULL <- file.path(DIR_FIGS, "supp_histopath__microcalcification_severity_NB_positive_full_axis_polished.png")
OUT_FIG_PVF_FULL <- file.path(DIR_FIGS, "supp_histopath__pvf_severity_gamma_positive_full_axis_polished.png")

K2 <- log10(2)
BASE_SIZE <- 10.5
DISPLAY_XMIN <- 0.03
DISPLAY_XMAX <- 20
DISPLAY_BREAKS <- c(0.03, 0.05, 0.10, 0.20, 0.50, 1, 2, 4, 10, 20)

# Clean log-axis labels for manuscript figures.
fmt_log_axis <- function(x) {
  dplyr::case_when(
    x < 0.01 ~ formatC(x, format = "f", digits = 3),
    x < 1 ~ formatC(x, format = "f", digits = 2),
    x < 10 ~ formatC(x, format = "f", digits = 1),
    TRUE ~ formatC(x, format = "f", digits = 0)
  )
}

# Use fewer breaks on full-axis plots so labels fit.
log_breaks <- function(lims) {
  br <- c(0.001, 0.003, 0.01, 0.03, 0.10, 0.30, 1, 3, 10, 30, 100)
  br[br >= lims[1] & br <= lims[2]]
}

# ----------------------------
# Polymer and family settings
# ----------------------------
POLYMER_ORDER_HISTO <- c("N6", "N66", "PET", "PU", "PE", "PP", "PVC", "PMMA", "PC", "SBR", "ABS")
POLYMER_ORDER_HISTO <- POLYMER_ORDER_HISTO[POLYMER_ORDER_HISTO %in% POLYMERS]

poly_family_local <- function(poly) {
  dplyr::case_when(
    poly %in% c("N6", "N66") ~ "polyamides",
    poly %in% c("PET") ~ "polyesters",
    poly %in% c("PU") ~ "polyurethanes",
    poly %in% c("PE", "PP") ~ "polyolefins",
    poly %in% c("PVC") ~ "vinyls",
    poly %in% c("PMMA") ~ "acrylics",
    poly %in% c("PC") ~ "polycarbonates",
    poly %in% c("SBR") ~ "styrenics",
    poly %in% c("ABS") ~ "other",
    TRUE ~ "other"
  )
}

poly_full <- function(poly) {
  if (poly %in% names(POLYMER_FULL_NAMES)) POLYMER_FULL_NAMES[[poly]] else poly
}

poly_plot_label <- function(poly) {
  if (poly %in% names(POLYMER_FULL_NAMES)) {
    paste0(POLYMER_FULL_NAMES[[poly]], " (", poly, ")")
  } else {
    poly
  }
}

# ----------------------------
# Model specifications
# ----------------------------
# M1_CLIN is the black primary model in figures.
# Covariates are automatically dropped if not available or non-informative.
MODEL_SPECS <- list(
  CRUDE = list(
    label = "CRUDE",
    covariates = character(0),
    color = "#56B4E9"
  ),
  M1_CLIN = list(
    label = "M1",
    covariates = c("GA", "Maternal_age", "BMI", "SDI_score", "Twins01", "Prior_PTB", "Hypertension01", "Hospital"),
    color = "black"
  ),
  M1_NO_SDI = list(
    label = "Primary without SDI",
    covariates = c("GA", "Maternal_age", "BMI", "Twins01", "Prior_PTB", "Hypertension01", "Hospital"),
    color = "#E69F00"
  ),
  M1_PLUS_SVI = list(
    label = "M1 + SVI",
    covariates = c("GA", "Maternal_age", "BMI", "SDI_score", "SVI_score", "Twins01", "Prior_PTB", "Hypertension01", "Hospital"),
    color = "#009E73"
  ),
  M1_PLUS_RACE = list(
    label = "Primary + race/ethnicity",
    covariates = c("GA", "Maternal_age", "BMI", "SDI_score", "Twins01", "Prior_PTB", "Hypertension01", "Hospital", "Ethnicity_Race"),
    color = "#CC79A7"
  )
)

PRIMARY_MODEL_TAG <- "M1_CLIN"

# ----------------------------
# Helpers
# ----------------------------
safe_num2 <- function(x) suppressWarnings(as.numeric(as.character(x)))

read_analysis_table <- function(path) {
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
      Twins01 = as.integer(as_01(Twins01)),
      Prior_PTB = as.integer(as_01(Prior_PTB)),
      Hypertension01 = as.integer(as_01(Hypertension01)),
      Hospital = factor(safe_chr(Hospital)),
      Ethnicity_Race = factor(safe_chr(Ethnicity_Race))
    )
}

read_path_scores <- function(path) {
  readr::read_csv(path, show_col_types = FALSE, guess_max = 10000) %>%
    dplyr::select(where(~ any(!is.na(.x)))) %>%
    dplyr::rename(PregnancyID = `Pregnancy ID`) %>%
    mutate(
      PregnancyID = as.character(PregnancyID),
      pvf_1 = safe_num2(pvf_1),
      pvf_2 = safe_num2(pvf_2),
      average_pvf = safe_num2(average_pvf),
      massive_pvf = as.integer(safe_num2(massive_pvf)),
      microcalc_1 = safe_num2(microcalc_1),
      microcalc_2 = safe_num2(microcalc_2),
      average_microcalcifications = safe_num2(average_microcalcifications),
      microcalc_presence = as.integer(safe_num2(microcalc_presence)),
      microcalc_total = microcalc_1 + microcalc_2,
      pvf_positive = is.finite(average_pvf) & average_pvf > 0,
      microcalc_positive = is.finite(microcalc_total) & microcalc_total > 0
    )
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

fit_one_model <- function(dat, outcome_col, exposure_col, polymer, outcome_name, model_type, model_tag, model_label, covariates, delivery_class) {
  covs <- usable_covariates(dat, covariates)

  vars <- c("PregnancyID", "SampleID", outcome_col, exposure_col, covs)
  d <- dat %>%
    dplyr::select(all_of(vars)) %>%
    mutate(
      outcome = safe_num2(.data[[outcome_col]]),
      exposure = safe_num2(.data[[exposure_col]])
    ) %>%
    dplyr::select(-all_of(c(outcome_col, exposure_col))) %>%
    drop_na() %>%
    mutate(across(where(is.factor), droplevels))

  n <- nrow(d)
  outcome_sum <- sum(d$outcome, na.rm = TRUE)

  base <- tibble(
    Outcome = outcome_name,
    ModelType = model_type,
    Delivery_class = as.character(delivery_class),
    ModelTag = model_tag,
    ModelLabel = model_label,
    Polymer = polymer,
    PolymerLabel = poly_plot_label(polymer),
    PolymerFull = poly_full(polymer),
    PolymerFamily = poly_family_local(polymer),
    ExposureCol = exposure_col,
    OutcomeCol = outcome_col,
    N = n,
    OutcomeSum = outcome_sum,
    Covariates = paste(covs, collapse = " + "),
    Ratio_2x = NA_real_,
    Ratio_2x_lo = NA_real_,
    Ratio_2x_hi = NA_real_,
    p = NA_real_,
    fit_status = NA_character_
  )

  if (n < 10 || length(unique(d$outcome)) < 2 || sd(d$exposure, na.rm = TRUE) == 0) {
    return(base %>% mutate(fit_status = "skipped_insufficient_variation_or_n"))
  }

  rhs <- c("exposure", covs)
  fml <- as.formula(paste("outcome ~", paste(rhs, collapse = " + ")))

  # Avoid heavily singular overfit models.
  mm <- try(model.matrix(fml, data = d), silent = TRUE)
  if (!inherits(mm, "try-error")) {
    if (qr(mm)$rank >= n - 2) {
      return(base %>% mutate(fit_status = "skipped_rank_too_high_for_n"))
    }
  }

  fit <- NULL
  if (model_type == "NB") {
    # NB is intended for integer-ish counts. microcalc_total is count-based.
    fit <- try(MASS::glm.nb(fml, data = d, control = glm.control(maxit = 100)), silent = TRUE)
  } else if (model_type == "Gamma") {
    d <- d %>% filter(outcome > 0)
    if (nrow(d) < 10 || length(unique(d$outcome)) < 2) {
      return(base %>% mutate(N = nrow(d), fit_status = "skipped_gamma_nonpositive_or_low_n"))
    }
    fit <- try(glm(fml, data = d, family = Gamma(link = "log"), control = glm.control(maxit = 100)), silent = TRUE)
  } else {
    stop("Unknown model_type: ", model_type)
  }

  if (inherits(fit, "try-error")) {
    return(base %>% mutate(fit_status = paste0("fit_error: ", as.character(fit)[1])))
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

  log_ratio <- b_est * K2
  se_log_ratio <- se_est * K2

  base %>%
    mutate(
      N = nrow(d),
      Ratio_2x = exp(log_ratio),
      Ratio_2x_lo = exp(log_ratio - 1.96 * se_log_ratio),
      Ratio_2x_hi = exp(log_ratio + 1.96 * se_log_ratio),
      p = p_val,
      fit_status = status
    )
}

add_q_values <- function(res) {
  ok_status <- c("ok_HC1", "ok_model_based_fallback", "ok_model_based")

  res %>%
    mutate(
      p_for_adjust = ifelse(fit_status %in% ok_status, p, NA_real_)
    ) %>%
    group_by(Outcome, Delivery_class, ModelTag) %>%
    mutate(
      q_global = {
        q <- rep(NA_real_, n())
        idx <- !is.na(p_for_adjust)
        q[idx] <- p.adjust(p_for_adjust[idx], method = "BH")
        q
      },
      sig_global_q05 = q_global < 0.05
    ) %>%
    ungroup() %>%
    group_by(Outcome, Delivery_class, ModelTag, PolymerFamily) %>%
    mutate(
      q_family = {
        q <- rep(NA_real_, n())
        idx <- !is.na(p_for_adjust)
        if (sum(idx) > 0) q[idx] <- p.adjust(p_for_adjust[idx], method = "BH")
        q
      },
      sig_family_q05 = q_family < 0.05
    ) %>%
    ungroup() %>%
    mutate(
      nominal_p05 = p < 0.05,
      Ratio_CI = ifelse(is.na(Ratio_2x), NA_character_, sprintf("%.2f (%.2f-%.2f)", Ratio_2x, Ratio_2x_lo, Ratio_2x_hi))
    )
}

make_axis_limits <- function(dat, xmin_floor = 0.001, xmax_cap = 100) {
  vals <- c(dat$Ratio_2x, dat$Ratio_2x_lo, dat$Ratio_2x_hi)
  vals <- vals[is.finite(vals) & vals > 0]
  if (length(vals) == 0) return(c(DISPLAY_XMIN, DISPLAY_XMAX))

  lo <- 10 ^ floor(log10(min(vals, na.rm = TRUE)))
  hi <- 10 ^ ceiling(log10(max(vals, na.rm = TRUE)))

  lo <- max(lo, xmin_floor)
  hi <- min(hi, xmax_cap)

  if (!is.finite(lo) || !is.finite(hi) || lo >= hi) c(DISPLAY_XMIN, DISPLAY_XMAX) else c(lo, hi)
}


prep_plot_dat <- function(tbl, lims) {
  tbl %>%
    filter(!is.na(Ratio_2x)) %>%
    mutate(
      Polymer = factor(Polymer, levels = rev(POLYMER_ORDER_HISTO)),
      PolymerLabel = factor(PolymerLabel, levels = rev(vapply(POLYMER_ORDER_HISTO, poly_plot_label, character(1)))),
      xmin = pmax(Ratio_2x_lo, lims[1]),
      xmax = pmin(Ratio_2x_hi, lims[2]),
      x_plot = pmin(pmax(Ratio_2x, lims[1]), lims[2]),
      left_trunc = Ratio_2x_lo < lims[1] | Ratio_2x < lims[1],
      right_trunc = Ratio_2x_hi > lims[2] | Ratio_2x > lims[2]
    )
}

make_histopath_plot <- function(res, outcome_name, model_type, out_path, full_axis = FALSE) {
  dat0 <- res %>%
    filter(Outcome == outcome_name)

  if (nrow(dat0) == 0) {
    warning("No rows available for outcome: ", outcome_name)
    return(invisible(NULL))
  }

  lims <- if (full_axis) make_axis_limits(dat0, xmin_floor = 0.001, xmax_cap = 100) else c(DISPLAY_XMIN, DISPLAY_XMAX)
  breaks <- if (full_axis) log_breaks(lims) else DISPLAY_BREAKS

  primary <- dat0 %>%
    filter(ModelTag == PRIMARY_MODEL_TAG) %>%
    prep_plot_dat(lims) %>%
    mutate(
      evidence = case_when(
        sig_global_q05 ~ "global BH q < 0.05",
        sig_family_q05 ~ "family BH q < 0.05",
        nominal_p05 ~ "nominal p < 0.05",
        TRUE ~ "not nominal"
      ),
      evidence = factor(evidence, levels = c(
        "not nominal",
        "nominal p < 0.05",
        "family BH q < 0.05",
        "global BH q < 0.05"
      )),
      stroke_val = ifelse(sig_family_q05, 1.5, 0.6)
    )

  sens <- dat0 %>%
    filter(ModelTag != PRIMARY_MODEL_TAG) %>%
    prep_plot_dat(lims) %>%
    mutate(
      ModelLabel = factor(ModelLabel, levels = vapply(MODEL_SPECS[names(MODEL_SPECS) != PRIMARY_MODEL_TAG], `[[`, character(1), "label"))
    )

  facet_counts <- primary %>%
    group_by(Delivery_class) %>%
    summarise(n_facet = max(N, na.rm = TRUE), .groups = "drop")

  primary <- primary %>%
    left_join(facet_counts, by = "Delivery_class") %>%
    mutate(FacetLabel = paste0(toupper(Delivery_class), "\nN=", n_facet))

  sens <- sens %>%
    left_join(facet_counts, by = "Delivery_class") %>%
    mutate(FacetLabel = paste0(toupper(Delivery_class), "\nN=", n_facet))

  model_colors <- vapply(MODEL_SPECS[names(MODEL_SPECS) != PRIMARY_MODEL_TAG], `[[`, character(1), "color")
  names(model_colors) <- vapply(MODEL_SPECS[names(MODEL_SPECS) != PRIMARY_MODEL_TAG], `[[`, character(1), "label")

  ratio_label <- ifelse(model_type == "NB", "IRR per doubling (log scale)", "RR per doubling (log scale)")
  title_txt <- ifelse(
    outcome_name == "microcalc_severity_positive",
    "Microcalcification severity among positives (negative binomial)",
    "Perivillous fibrin deposition severity among positives (gamma)"
  )

  subtitle_txt <- str_wrap(
    "Black intervals show the M1 model; colored points show crude and sensitivity models. Point shape denotes BH evidence.",
    width = 100
  )

  p <- ggplot() +
    geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.35, color = "grey45") +
    geom_segment(
      data = primary,
      aes(y = PolymerLabel, yend = PolymerLabel, x = xmin, xend = xmax),
      linewidth = 0.65,
      color = "black"
    ) +
    geom_point(
      data = sens,
      aes(y = PolymerLabel, x = x_plot, color = ModelLabel),
      size = 1.7,
      alpha = 0.9
    ) +
    geom_point(
      data = primary,
      aes(y = PolymerLabel, x = x_plot, shape = evidence),
      size = 2.5,
      color = "black",
      stroke = primary$stroke_val
    ) +
    geom_text(
      data = primary %>% filter(left_trunc),
      aes(y = PolymerLabel, x = lims[1], label = "\u2190"),
      hjust = 1.05,
      size = 3.0
    ) +
    geom_text(
      data = primary %>% filter(right_trunc),
      aes(y = PolymerLabel, x = lims[2], label = "\u2192"),
      hjust = -0.05,
      size = 3.0
    ) +
    facet_wrap(~FacetLabel, nrow = 1, scales = "free_y") +
    scale_x_log10(breaks = breaks, limits = lims, labels = fmt_log_axis) +
    scale_color_manual(values = model_colors, drop = FALSE) +
    scale_shape_manual(
      values = c(
        "not nominal" = 1,
        "nominal p < 0.05" = 16,
        "family BH q < 0.05" = 22,
        "global BH q < 0.05" = 15
      ),
      name = "M1 evidence"
    ) +
    labs(
      title = title_txt,
      subtitle = subtitle_txt,
      x = ratio_label,
      y = NULL,
      color = "Sensitivity model"
    ) +
    guides(
      shape = guide_legend(order = 1, nrow = 1, byrow = TRUE),
      color = guide_legend(order = 2, nrow = 2, byrow = TRUE, override.aes = list(size = 2.4))
    ) +
    theme_mnp_forest(base_size = BASE_SIZE) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold", size = 11),
      plot.subtitle = element_text(size = 8.5),
      strip.background = element_rect(fill = "grey90", color = "grey50"),
      strip.text = element_text(size = 8.5),
      legend.position = "bottom",
      legend.box = "vertical",
      legend.box.just = "left",
      legend.title = element_text(size = 8, face = "bold"),
      legend.text = element_text(size = 7),
      legend.spacing.y = grid::unit(0.05, "cm"),
      axis.text.y = element_text(size = 7),
      axis.text.x = element_text(size = 7)
    )

  ggsave_mnp(out_path, p, width = 10.8, height = 6.2, dpi = 300)
}

# ----------------------------
# Load and join
# ----------------------------
message("Reading pathology scores: ", PATH_FILE)
path_df <- read_path_scores(PATH_FILE)
model_df <- read_analysis_table(IN_MODEL_TABLE)

joined <- model_df %>%
  left_join(path_df, by = "PregnancyID") %>%
  mutate(
    has_path_scores = !is.na(average_pvf) | !is.na(average_microcalcifications)
  )

write_csv_clean(joined, OUT_JOINED)

join_audit <- tibble(
  path_file = PATH_FILE,
  model_rows = nrow(model_df),
  path_rows = nrow(path_df),
  joined_rows = nrow(joined),
  rows_with_any_path_score = sum(joined$has_path_scores, na.rm = TRUE),
  path_ids_not_in_model = sum(!path_df$PregnancyID %in% model_df$PregnancyID),
  model_ids_with_path = sum(model_df$PregnancyID %in% path_df$PregnancyID)
)

write_csv_clean(join_audit, OUT_JOIN_AUDIT)

# Polymer audit
exposure_tbl <- tibble(
  Polymer = POLYMER_ORDER_HISTO,
  ExposureCol = legacy_col_log10(POLYMER_ORDER_HISTO),
  UgCol = legacy_col_ug_g(POLYMER_ORDER_HISTO)
) %>%
  mutate(
    exposure_available = ExposureCol %in% names(joined),
    n_detected_gt0 = vapply(UgCol, function(cn) if (cn %in% names(joined)) detected_n(joined[[cn]]) else 0, numeric(1)),
    continuous_model_ok = exposure_available & n_detected_gt0 >= MIN_DETECTED_N_FOR_CONTINUOUS_MODEL,
    skip_reason = case_when(
      !exposure_available ~ "missing_exposure_column",
      !continuous_model_ok ~ "sparse_detection",
      TRUE ~ NA_character_
    )
  )

write_csv_clean(exposure_tbl, OUT_POLYMER_AUDIT)

eligible_exposures <- exposure_tbl %>%
  filter(is.na(skip_reason))

# ----------------------------
# Fit models
# ----------------------------
model_inputs <- list(
  microcalc_severity_positive = list(
    outcome_col = "microcalc_total",
    model_type = "NB",
    subset_expr = expr(microcalc_positive)
  ),
  pvf_severity_positive = list(
    outcome_col = "average_pvf",
    model_type = "Gamma",
    subset_expr = expr(pvf_positive)
  )
)

all_results <- purrr::imap_dfr(model_inputs, function(info, outcome_name) {
  dat_outcome <- joined %>%
    filter(has_path_scores) %>%
    filter(!!info$subset_expr)

  purrr::map_dfr(levels(joined$Delivery_class), function(dc) {
    dat_dc <- dat_outcome %>%
      filter(Delivery_class == dc)

    purrr::imap_dfr(MODEL_SPECS, function(spec, model_tag) {
      purrr::pmap_dfr(
        eligible_exposures,
        function(Polymer, ExposureCol, UgCol, exposure_available, n_detected_gt0, continuous_model_ok, skip_reason) {
          fit_one_model(
            dat = dat_dc,
            outcome_col = info$outcome_col,
            exposure_col = ExposureCol,
            polymer = Polymer,
            outcome_name = outcome_name,
            model_type = info$model_type,
            model_tag = model_tag,
            model_label = spec$label,
            covariates = spec$covariates,
            delivery_class = dc
          )
        }
      )
    })
  })
}) %>%
  add_q_values()

write_csv_clean(all_results, OUT_MODEL_ALL)

micro_res <- all_results %>% filter(Outcome == "microcalc_severity_positive")
pvf_res <- all_results %>% filter(Outcome == "pvf_severity_positive")

write_csv_clean(micro_res, OUT_MODEL_MICRO)
write_csv_clean(pvf_res, OUT_MODEL_PVF)

q_audit <- all_results %>%
  filter(ModelTag == PRIMARY_MODEL_TAG) %>%
  arrange(Outcome, Delivery_class, p) %>%
  dplyr::select(
    Outcome, Delivery_class, ModelTag, Polymer, PolymerLabel, PolymerFamily,
    N, Ratio_2x, Ratio_2x_lo, Ratio_2x_hi, p, q_global, q_family,
    nominal_p05, sig_global_q05, sig_family_q05, fit_status
  )

write_csv_clean(q_audit, OUT_Q_AUDIT)

model_key <- purrr::imap_dfr(MODEL_SPECS, function(spec, tag) {
  tibble(
    ModelTag = tag,
    ModelLabel = spec$label,
    CovariatesRequested = paste(spec$covariates, collapse = " + "),
    Color = spec$color
  )
})
write_csv_clean(model_key, OUT_MODEL_KEY)

# ----------------------------
# Figures
# ----------------------------
make_histopath_plot(
  res = all_results,
  outcome_name = "microcalc_severity_positive",
  model_type = "NB",
  out_path = OUT_FIG_MICRO,
  full_axis = FALSE
)
make_histopath_plot(
  res = all_results,
  outcome_name = "pvf_severity_positive",
  model_type = "Gamma",
  out_path = OUT_FIG_PVF,
  full_axis = FALSE
)
make_histopath_plot(
  res = all_results,
  outcome_name = "microcalc_severity_positive",
  model_type = "NB",
  out_path = OUT_FIG_MICRO_FULL,
  full_axis = TRUE
)
make_histopath_plot(
  res = all_results,
  outcome_name = "pvf_severity_positive",
  model_type = "Gamma",
  out_path = OUT_FIG_PVF_FULL,
  full_axis = TRUE
)

legend_text <- paste0(
  "Supplementary Figure | Polymer-specific associations with histopathology severity.\n",
  "Pathology scores were joined to the analysis-ready MNP table by PregnancyID. ",
  "Microcalcification severity was modeled among microcalcification-positive samples using negative binomial regression of the summed duplicate microcalcification score. ",
  "Perivillous fibrin deposition severity was modeled among PVF-positive samples using gamma regression with a log link applied to the average PVF score. ",
  "Models were stratified by delivery class. Estimates are shown as incidence rate ratios (IRRs; negative binomial models) or relative ratios (RRs; gamma models) per doubling of placental polymer burden. ",
  "Black points and confidence intervals represent the M1 clinical model. Colored points represent crude and sensitivity models. ",
  "Global BH correction was calculated across modeled polymers within each outcome, delivery class, and model. Family BH correction was calculated within polymer family. ",
  "The x-axis is shown on a log scale; arrows indicate estimates or confidence intervals truncated for display."
)
writeLines(legend_text, OUT_LEGEND)

# ----------------------------
# Log
# ----------------------------
run_log <- c(
  "09_histopathology_models.R",
  paste0("PROJECT_ROOT: ", PROJECT_ROOT),
  paste0("Pathology file: ", PATH_FILE),
  "",
  "Join audit:",
  capture.output(print(join_audit)),
  "",
  "Polymer audit:",
  capture.output(print(exposure_tbl, n = Inf)),
  "",
  "Path score counts:",
  capture.output(joined %>% filter(has_path_scores) %>% count(Delivery_class, microcalc_positive, pvf_positive) %>% print(n = Inf)),
  "",
  "Primary M1 q audit:",
  capture.output(print(q_audit, n = Inf)),
  "",
  "Wrote:",
  OUT_JOINED,
  OUT_JOIN_AUDIT,
  OUT_MODEL_ALL,
  OUT_MODEL_MICRO,
  OUT_MODEL_PVF,
  OUT_Q_AUDIT,
  OUT_MODEL_KEY,
  OUT_POLYMER_AUDIT,
  OUT_FIG_MICRO,
  OUT_FIG_PVF,
  OUT_FIG_MICRO_FULL,
  OUT_FIG_PVF_FULL,
  OUT_LEGEND
)

cat(paste(run_log, collapse = "\n"), file = OUT_LOG)
message(paste(run_log, collapse = "\n"))
