#!/usr/bin/env Rscript
# ============================================================
# 12_clinical_outcome_correlations.R
# v4.6: dplyr::select namespace fix for masked select()
# Spearman correlations between polymer-specific MNP burden and clinical outcomes
# ============================================================
#
# Run after run_analysis_setup.R has already created:
#   outputs/tables/mnp_ptb_cohort__analysis_ready_model_table.csv
#
# You do NOT need an in-memory df.
# This script reads the analysis-ready CSV from disk.
#
# Outputs:
#   Fig. 2a-style correlation heatmap
#   Supplementary Table 1c-style correlation table
#
# Analysis:
#   Spearman correlation between log10 polymer burden and clinical outcomes.
#   BH correction is performed across the full polymer/outcome correlation family.
#
# Default clinical outcomes:
#   GA
#   Birth_weight
#   Prior_PTB
#   SDI_score
#   SVI_score
#   BMI
#   Maternal_age
#   Hypertension01
#   Twins01
#
# Notes:
#   - Prior_PTB is binary, so Spearman rho is rank-biserial-like but acceptable
#     for an exploratory correlation heatmap.
#   - Outcome labels can be modified below.
#   - Figure includes cumulative MNP burden but exact wording should emphasize
#     polymer-specific signals if cumulative burden is weaker.
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
OUT_DIR <- file.path(DIR_TABLES, paste0("12_clinical_correlations__", STAMP))
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

OUT_ANALYSIS_TABLE <- file.path(OUT_DIR, "clinical_correlation_analysis_table.csv")
OUT_LONG <- file.path(OUT_DIR, "clinical_correlations_long.csv")
OUT_SUPP_TABLE <- file.path(OUT_DIR, "Supplementary_Table_1c_clinical_correlations.csv")
OUT_TOP <- file.path(OUT_DIR, "clinical_correlations_top_hits.csv")
OUT_OUTCOME_AUDIT <- file.path(OUT_DIR, "clinical_correlation_outcome_audit.csv")
OUT_POLYMER_AUDIT <- file.path(OUT_DIR, "clinical_correlation_polymer_audit.csv")
OUT_FIG_DATA <- file.path(OUT_DIR, "Fig2a_clinical_correlation_heatmap_data.csv")
OUT_LEGEND <- file.path(DIR_TABLES, "Fig2a_clinical_correlation_figure_legend.txt")
OUT_LOG <- file.path(DIR_LOGS, "12_clinical_outcome_correlations__run_log.txt")

OUT_FIG_HEATMAP <- file.path(DIR_FIGS, "Fig2a_clinical_correlation_heatmap.png")
OUT_FIG_HEATMAP_PDF <- file.path(DIR_FIGS, "Fig2a_clinical_correlation_heatmap.pdf")
OUT_FIG_DOT <- file.path(DIR_FIGS, "Fig2a_clinical_correlation_dotplot.png")
OUT_FIG_DOT_PDF <- file.path(DIR_FIGS, "Fig2a_clinical_correlation_dotplot.pdf")

BASE_SIZE <- 10.5

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

# Put clinically central outcomes first.
DEFAULT_OUTCOMES <- c(
  "GA",
  "Birth_weight",
  "Prior_PTB",
  "SDI_score",
  "SVI_score",
  "BMI",
  "Maternal_age",
  "Hypertension01",
  "Twins01"
)

OUTCOME_LABELS <- c(
  GA = "Gestational age",
  Birth_weight = "Birth weight",
  Prior_PTB = "Prior PTB",
  SDI_score = "SDI",
  SVI_score = "SVI",
  BMI = "BMI",
  Maternal_age = "Maternal age",
  Hypertension01 = "Hypertension",
  Twins01 = "Twin gestation"
)

# ----------------------------
# Helpers
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

fmt_p_local <- function(p, digits = 3) {
  ifelse(
    is.na(p),
    NA_character_,
    ifelse(p < 0.001, "<0.001", formatC(p, format = "f", digits = digits))
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

poly_short_label <- function(poly) {
  poly <- as.character(poly)
  dplyr::case_when(
    poly == "Cumulative_MNPs" ~ "Total MNP burden",
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
      Hypertension01 = as.integer(as_01_local(Hypertension01))
    )
}

find_birth_weight_col <- function(dat) {
  candidates <- c(
    "Birth_weight",
    "Birth_weight_g",
    "BirthWeight",
    "birth_weight",
    "birthweight",
    "Birth.weight",
    "birth_weight_g",
    "Birth weight",
    "Birth weight (g)",
    "Birthweight_g"
  )
  hit <- candidates[candidates %in% names(dat)]
  if (length(hit) > 0) return(hit[1])
  regex_hit <- names(dat)[stringr::str_detect(names(dat), regex("birth.*weight|weight.*birth", ignore_case = TRUE))]
  if (length(regex_hit) > 0) return(regex_hit[1])
  NA_character_
}

spearman_one <- function(dat, exposure_col, outcome_col, polymer, outcome_label) {
  d <- dat %>%
    dplyr::transmute(
      exposure = safe_num2(.data[[exposure_col]]),
      outcome = safe_num2(.data[[outcome_col]])
    ) %>%
    filter(is.finite(exposure), is.finite(outcome))

  n <- nrow(d)

  base <- tibble(
    Polymer = polymer,
    PolymerLabel = poly_label_local(polymer),
    PolymerShort = poly_short_label(polymer),
    PolymerFamily = poly_family_local(polymer),
    ExposureCol = exposure_col,
    Outcome = outcome_col,
    OutcomeLabel = outcome_label,
    N = n,
    rho = NA_real_,
    p = NA_real_,
    fit_status = NA_character_
  )

  if (n < 6) {
    return(base %>% mutate(fit_status = "skipped_n_lt_6"))
  }

  if (length(unique(d$exposure)) < 2 || length(unique(d$outcome)) < 2) {
    return(base %>% mutate(fit_status = "skipped_no_variation"))
  }

  ct <- suppressWarnings(
    tryCatch(
      cor.test(d$exposure, d$outcome, method = "spearman", exact = FALSE),
      error = function(e) e
    )
  )

  if (inherits(ct, "error")) {
    return(base %>% mutate(fit_status = paste0("cor_error: ", ct$message)))
  }

  base %>%
    mutate(
      rho = unname(ct$estimate),
      p = ct$p.value,
      fit_status = "ok"
    )
}

# ----------------------------
# Load and prepare data
# ----------------------------
if (!file.exists(IN_MODEL_TABLE)) {
  stop("Cannot find analysis-ready model table: ", IN_MODEL_TABLE, call. = FALSE)
}

df0 <- read_model_table(IN_MODEL_TABLE)

bw_col <- find_birth_weight_col(df0)
if (!is.na(bw_col) && bw_col != "Birth_weight") {
  df0 <- df0 %>% mutate(Birth_weight = safe_num2(.data[[bw_col]]))
}
if (is.na(bw_col) && !"Birth_weight" %in% names(df0)) {
  warning("No birth weight column found. Birth_weight will be omitted.")
}

# Keep available clinical outcomes only.
available_outcomes <- DEFAULT_OUTCOMES[DEFAULT_OUTCOMES %in% names(df0)]

df <- df0 %>%
  mutate(across(any_of(available_outcomes), safe_num2)) %>%
  mutate(
    Prior_PTB = if ("Prior_PTB" %in% names(.)) as.integer(as_01_local(Prior_PTB)) else NA_integer_,
    Hypertension01 = if ("Hypertension01" %in% names(.)) as.integer(as_01_local(Hypertension01)) else NA_integer_,
    Twins01 = if ("Twins01" %in% names(.)) as.integer(as_01_local(Twins01)) else NA_integer_
  )

write_csv_clean_local(df, OUT_ANALYSIS_TABLE)

outcome_audit <- tibble(
  OutcomeRequested = DEFAULT_OUTCOMES,
  OutcomeAvailable = DEFAULT_OUTCOMES %in% names(df),
  OutcomeLabel = ifelse(DEFAULT_OUTCOMES %in% names(OUTCOME_LABELS), OUTCOME_LABELS[DEFAULT_OUTCOMES], DEFAULT_OUTCOMES),
  N_nonmissing = vapply(DEFAULT_OUTCOMES, function(cn) if (cn %in% names(df)) sum(is.finite(safe_num2(df[[cn]]))) else 0, integer(1)),
  N_unique = vapply(DEFAULT_OUTCOMES, function(cn) if (cn %in% names(df)) length(unique(safe_num2(df[[cn]])[is.finite(safe_num2(df[[cn]]))])) else 0, integer(1))
)

write_csv_clean_local(outcome_audit, OUT_OUTCOME_AUDIT)

exposure_tbl <- tibble(
  Polymer = ANALYTES,
  ExposureCol = vapply(ANALYTES, legacy_col_log10_local, character(1)),
  UgCol = vapply(ANALYTES, legacy_col_ug_g_local, character(1))
) %>%
  mutate(
    PolymerLabel = vapply(Polymer, poly_label_local, character(1)),
    PolymerShort = vapply(Polymer, poly_short_label, character(1)),
    PolymerFamily = vapply(Polymer, poly_family_local, character(1)),
    exposure_available = ExposureCol %in% names(df),
    n_nonmissing = vapply(ExposureCol, function(cn) if (cn %in% names(df)) sum(is.finite(safe_num2(df[[cn]]))) else 0, integer(1)),
    n_detected_gt0 = vapply(UgCol, function(cn) if (cn %in% names(df)) detected_n_local(df[[cn]]) else 0, numeric(1)),
    correlation_ok = exposure_available & n_nonmissing >= 6 & (Polymer == "Cumulative_MNPs" | n_detected_gt0 >= MIN_DETECTED_N_FOR_CONTINUOUS_MODEL),
    skip_reason = case_when(
      !exposure_available ~ "missing_exposure_column",
      n_nonmissing < 6 ~ "too_few_nonmissing",
      Polymer != "Cumulative_MNPs" & n_detected_gt0 < MIN_DETECTED_N_FOR_CONTINUOUS_MODEL ~ "sparse_detection",
      TRUE ~ NA_character_
    )
  )

write_csv_clean_local(exposure_tbl, OUT_POLYMER_AUDIT)

eligible_exposures <- exposure_tbl %>% filter(correlation_ok)
eligible_outcomes <- outcome_audit %>% filter(OutcomeAvailable, N_nonmissing >= 6, N_unique >= 2)

# ----------------------------
# Spearman correlations
# ----------------------------
corr_long <- purrr::pmap_dfr(
  tidyr::crossing(
    eligible_exposures %>% select(Polymer, ExposureCol),
    eligible_outcomes %>% select(OutcomeRequested, OutcomeLabel)
  ),
  function(Polymer, ExposureCol, OutcomeRequested, OutcomeLabel) {
    spearman_one(
      dat = df,
      exposure_col = ExposureCol,
      outcome_col = OutcomeRequested,
      polymer = Polymer,
      outcome_label = OutcomeLabel
    )
  }
) %>%
  mutate(
    q_BH_global = {
      q <- rep(NA_real_, n())
      idx <- fit_status == "ok" & !is.na(p)
      q[idx] <- p.adjust(p[idx], method = "BH")
      q
    }
  ) %>%
  group_by(Outcome) %>%
  mutate(
    q_BH_by_outcome = {
      q <- rep(NA_real_, n())
      idx <- fit_status == "ok" & !is.na(p)
      if (sum(idx) > 0) q[idx] <- p.adjust(p[idx], method = "BH")
      q
    }
  ) %>%
  ungroup() %>%
  group_by(Polymer) %>%
  mutate(
    q_BH_by_polymer = {
      q <- rep(NA_real_, n())
      idx <- fit_status == "ok" & !is.na(p)
      if (sum(idx) > 0) q[idx] <- p.adjust(p[idx], method = "BH")
      q
    }
  ) %>%
  ungroup() %>%
  mutate(
    q_BH = q_BH_global,
    sig_q05 = q_BH_global < 0.05,
    nominal_p05 = p < 0.05,
    rho_label = ifelse(is.na(rho), "", sprintf("%.2f", rho)),
    p_label = fmt_p_local(p),
    q_label = fmt_p_local(q_BH_global),
    significance_label = case_when(
      q_BH_global < 0.05 ~ "BH q<0.05",
      p < 0.05 ~ "nominal p<0.05",
      TRUE ~ "not significant"
    )
  )

write_csv_clean_local(corr_long, OUT_LONG)
write_csv_clean_local(corr_long, OUT_SUPP_TABLE)

top_hits <- corr_long %>%
  filter(fit_status == "ok") %>%
  arrange(q_BH_global, desc(abs(rho))) %>%
  dplyr::select(
    Polymer, PolymerLabel, PolymerFamily,
    Outcome, OutcomeLabel, N, rho, p, q_BH_global, q_BH_by_outcome, q_BH_by_polymer,
    significance_label, fit_status
  )

write_csv_clean_local(top_hits, OUT_TOP)

# ----------------------------
# Figure data
# ----------------------------
poly_order <- c("Cumulative_MNPs", "PE", "SBR", "PP", "PVC", "N66", "PET", "N6", "PMMA", "PC", "PU", "ABS", "PS")
poly_order <- poly_order[poly_order %in% corr_long$Polymer]

outcome_order <- available_outcomes[available_outcomes %in% corr_long$Outcome]

fig_dat <- corr_long %>%
  mutate(
    PolymerShort = factor(PolymerShort, levels = rev(vapply(poly_order, poly_short_label, character(1)))),
    OutcomeLabel = factor(OutcomeLabel, levels = OUTCOME_LABELS[outcome_order]),
    sig_symbol = case_when(
      q_BH_global < 0.05 ~ "\u25B2",
      p < 0.05 ~ "\u25CF",
      TRUE ~ ""
    ),
    cell_label = ifelse(sig_symbol == "", rho_label, paste0(rho_label, "\n", sig_symbol))
  )

write_csv_clean_local(fig_dat, OUT_FIG_DATA)

# ----------------------------
# Heatmap
# ----------------------------
heat_title <- "Polymer-specific correlations with clinical outcomes"
heat_subtitle <- str_wrap(
  "Spearman rho between log10 polymer burden and clinical outcomes. Triangles indicate BH q<0.05 across the full matrix; circles indicate nominal p<0.05 only.",
  width = 110
)

p_heat <- ggplot(fig_dat, aes(x = OutcomeLabel, y = PolymerShort, fill = rho)) +
  geom_tile(color = "white", linewidth = 0.45) +
  geom_text(aes(label = cell_label), size = 2.7, lineheight = 0.85) +
  scale_fill_gradient2(
    low = "#2C7BB6",
    mid = "white",
    high = "#D7191C",
    midpoint = 0,
    limits = c(-1, 1),
    oob = scales::squish,
    name = "Spearman\nrho"
  ) +
  labs(
    title = heat_title,
    subtitle = heat_subtitle,
    x = NULL,
    y = NULL
  ) +
  theme_bw(base_size = BASE_SIZE) +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 9.2),
    axis.text.x = element_text(angle = 35, hjust = 1),
    axis.text.y = element_text(size = 8),
    legend.title = element_text(face = "bold")
  )

ggsave(OUT_FIG_HEATMAP, p_heat, width = 8.8, height = 6.2, dpi = 300)
ggsave(OUT_FIG_HEATMAP_PDF, p_heat, width = 8.8, height = 6.2)

# ----------------------------
# Dotplot alternative
# ----------------------------
p_dot <- ggplot(fig_dat, aes(x = OutcomeLabel, y = PolymerShort)) +
  geom_point(aes(size = abs(rho), color = rho, shape = sig_q05), alpha = 0.95, stroke = 0.7) +
  scale_color_gradient2(
    low = "#2C7BB6",
    mid = "grey85",
    high = "#D7191C",
    midpoint = 0,
    limits = c(-1, 1),
    oob = scales::squish,
    name = "Spearman\nrho"
  ) +
  scale_size_continuous(range = c(1.2, 6.0), name = "|rho|") +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17), name = "BH q < 0.05") +
  labs(
    title = heat_title,
    subtitle = "Dot size reflects absolute Spearman rho; color reflects direction.",
    x = NULL,
    y = NULL
  ) +
  theme_bw(base_size = BASE_SIZE) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 9.2),
    axis.text.x = element_text(angle = 35, hjust = 1),
    axis.text.y = element_text(size = 8),
    legend.title = element_text(face = "bold")
  )

ggsave(OUT_FIG_DOT, p_dot, width = 8.8, height = 6.2, dpi = 300)
ggsave(OUT_FIG_DOT_PDF, p_dot, width = 8.8, height = 6.2)

# ----------------------------
# Legend text
# ----------------------------
legend_text <- paste0(
  "Figure 2a | Polymer-specific correlations with clinical outcomes. ",
  "Spearman correlation analysis was used to evaluate relationships between log10-transformed polymer-specific placental burden and clinical outcomes. ",
  "Clinical outcomes included gestational age at delivery, birth weight, prior preterm birth, Social Deprivation Index, Social Vulnerability Index, BMI, maternal age, hypertension, and twin gestation when available. ",
  "Benjamini-Hochberg correction was performed across the full polymer-by-outcome correlation matrix. ",
  "Heatmap color indicates Spearman rho. Triangles indicate BH q<0.05 across the full matrix; circles indicate nominal p<0.05 only. ",
  "Exact rho, p values, and q values are reported in Supplementary Table 1c."
)
writeLines(legend_text, OUT_LEGEND)

# ----------------------------
# Log
# ----------------------------
run_log <- c(
  "12_clinical_outcome_correlations.R",
  paste0("PROJECT_ROOT: ", PROJECT_ROOT),
  paste0("Input table: ", IN_MODEL_TABLE),
  paste0("Rows in input table: ", nrow(df0)),
  paste0("Birth weight source column: ", ifelse(is.na(bw_col), "not found", bw_col)),
  "",
  "Available outcomes:",
  capture.output(print(outcome_audit, n = Inf)),
  "",
  "Eligible exposures:",
  capture.output(print(exposure_tbl, n = Inf)),
  "",
  "Top correlations:",
  capture.output(print(top_hits %>% select(Polymer, OutcomeLabel, N, rho, p, q_BH_global, significance_label) %>% head(30), n = Inf)),
  "",
  "Wrote:",
  OUT_ANALYSIS_TABLE,
  OUT_LONG,
  OUT_SUPP_TABLE,
  OUT_TOP,
  OUT_OUTCOME_AUDIT,
  OUT_POLYMER_AUDIT,
  OUT_FIG_DATA,
  OUT_FIG_HEATMAP,
  OUT_FIG_HEATMAP_PDF,
  OUT_FIG_DOT,
  OUT_FIG_DOT_PDF,
  OUT_LEGEND
)

cat(paste(run_log, collapse = "\n"), file = OUT_LOG)
message(paste(run_log, collapse = "\n"))
