#!/usr/bin/env Rscript
# POLISHED VERSION: sources R/figure_theme_mnp.R and writes *_polished outputs where applicable.
# ============================================================
# 04_descriptive_and_univariate.R
# Cohort descriptive tables + polymer landscape + univariate PTB/GA analyses
# ============================================================
# Input:
#   outputs/tables/mnp_ptb_cohort__analysis_ready_model_table.csv
#
# Outputs:
#   outputs/tables/table1__cohort_characteristics_term_vs_ptb.csv
#   outputs/tables/table2__detection_and_burden_by_ptb.csv
#   outputs/tables/univariate__ugpg_summary_term_vs_ptb.csv
#   outputs/tables/univariate__log10ugpg_summary_term_vs_ptb.csv
#   outputs/tables/univariate__GA_spearman__log10ugpg.csv
#   outputs/figs/fig2__polymer_landscape__detection_and_burden.png
#   outputs/figs/fig3__univariate_ptb__burden_by_polymer.png
#   outputs/figs/fig3__univariate_ptb__volcano.png
#   outputs/figs/univariate__GA_spearman__lollipop.png
#   outputs/logs/04_descriptive_and_univariate__run_log.txt
#
# Notes:
#   - PTB status is derived from Preterm/Delivery_class already produced by script 01.
#   - PS is retained in descriptive outputs but flagged as sparse when n detected is low.
#   - Cumulative MNP burden is listed first in Table 2.
#   - PNG only.
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
# Inputs / outputs
# ----------------------------
IN_CSV <- file.path(DIR_TABLES, "mnp_ptb_cohort__analysis_ready_model_table.csv")
stopifnot(file.exists(IN_CSV))

OUT_TABLE1 <- file.path(DIR_TABLES, "table1__cohort_characteristics_term_vs_ptb_polished.csv")
OUT_TABLE1_HTML <- file.path(DIR_TABLES, "table1__cohort_characteristics_term_vs_ptb_polished.html")
OUT_TABLE1_PNG <- file.path(DIR_TABLES, "table1__cohort_characteristics_term_vs_ptb_polished.png")

OUT_TABLE2 <- file.path(DIR_TABLES, "table2__detection_and_burden_by_ptb_polished.csv")
OUT_TABLE2_HTML <- file.path(DIR_TABLES, "table2__detection_and_burden_by_ptb_polished.html")
OUT_TABLE2_PNG <- file.path(DIR_TABLES, "table2__detection_and_burden_by_ptb_polished.png")

OUT_UNI_UG <- file.path(DIR_TABLES, "univariate__ugpg_summary_term_vs_ptb_polished.csv")
OUT_UNI_LOG <- file.path(DIR_TABLES, "univariate__log10ugpg_summary_term_vs_ptb_polished.csv")
OUT_GA <- file.path(DIR_TABLES, "univariate__GA_spearman__log10ugpg_polished.csv")

OUT_FIG2 <- file.path(DIR_FIGS, "fig2__polymer_landscape__detection_and_burden_polished.png")
OUT_FIG3_BURDEN <- file.path(DIR_FIGS, "fig3__univariate_ptb__burden_by_polymer_polished.png")
OUT_FIG3_VOLCANO <- file.path(DIR_FIGS, "fig3__univariate_ptb__volcano_polished.png")
OUT_FIG_GA <- file.path(DIR_FIGS, "univariate__GA_spearman__lollipop_polished.png")
OUT_LOG <- file.path(DIR_LOGS, "04_descriptive_and_univariate__run_log_polished.txt")

BASE_SIZE <- 11
set.seed(1234)

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
      `Pregnancy ID` = as.character(`Pregnancy ID`),
      home_zip_2nd_3rd_trim = if ("home_zip_2nd_3rd_trim" %in% names(.)) pad_zip(home_zip_2nd_3rd_trim) else NA_character_,
      Delivery_class = factor(as.character(Delivery_class), levels = c("Term", "Preterm")),
      PTB01 = case_when(
        Delivery_class == "Preterm" ~ 1L,
        Delivery_class == "Term" ~ 0L,
        TRUE ~ NA_integer_
      ),
      Preterm = PTB01 == 1L,
      GA = safe_num(GA),
      GA_weeks = safe_num(GA_weeks)
    )
}

p_cont <- function(x, g) {
  d <- tibble(x = safe_num(x), g = factor(g)) %>% filter(is.finite(x), !is.na(g))
  if (nrow(d) < 3 || n_distinct(d$g) < 2) return(NA_real_)
  suppressWarnings(stats::wilcox.test(x ~ g, data = d, exact = FALSE)$p.value)
}

p_cat <- function(x, g) {
  d <- tibble(x = as.character(x), g = factor(g)) %>%
    filter(!is.na(x), x != "", !is.na(g))
  if (nrow(d) < 3 || n_distinct(d$g) < 2 || n_distinct(d$x) < 2) return(NA_real_)
  tab <- table(d$g, d$x)
  suppressWarnings({
    out <- try(stats::fisher.test(tab)$p.value, silent = TRUE)
    if (inherits(out, "try-error")) {
      out <- try(stats::chisq.test(tab)$p.value, silent = TRUE)
      if (inherits(out, "try-error")) NA_real_ else out
    } else {
      out
    }
  })
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

make_table1_gt <- function(tbl) {
  tbl %>%
    transmute(
      Variable = variable,
      Overall,
      Term,
      Preterm,
      `p-value` = p_fmt,
      SMD
    ) %>%
    gt::gt() %>%
    gt::tab_header(
      title = gt::md("**Table 1. Cohort Characteristics by Delivery Class**")
    ) %>%
    gt::cols_align(align = "left", columns = c(Variable, Overall, Term, Preterm)) %>%
    gt::cols_align(align = "center", columns = c(`p-value`, SMD)) %>%
    gt::tab_options(table.font.size = gt::px(12))
}

make_table2_gt <- function(tbl) {
  tbl %>%
    mutate(
      `Term detection` = ifelse(!is.na(n_detected_Term), sprintf("%d/%d (%.1f%%)", n_detected_Term, N_Term, pct_detected_Term), NA_character_),
      `Preterm detection` = ifelse(!is.na(n_detected_Preterm), sprintf("%d/%d (%.1f%%)", n_detected_Preterm, N_Preterm, pct_detected_Preterm), NA_character_),
      `Term burden, mean ± SD` = ifelse(!is.na(mean_ug_g_Term), sprintf("%.2f ± %.2f", mean_ug_g_Term, sd_ug_g_Term), NA_character_),
      `Preterm burden, mean ± SD` = ifelse(!is.na(mean_ug_g_Preterm), sprintf("%.2f ± %.2f", mean_ug_g_Preterm, sd_ug_g_Preterm), NA_character_),
      `Term burden, median (IQR)` = ifelse(!is.na(median_ug_g_Term), sprintf("%.2f (%.2f, %.2f)", median_ug_g_Term, q1_ug_g_Term, q3_ug_g_Term), NA_character_),
      `Preterm burden, median (IQR)` = ifelse(!is.na(median_ug_g_Preterm), sprintf("%.2f (%.2f, %.2f)", median_ug_g_Preterm, q1_ug_g_Preterm, q3_ug_g_Preterm), NA_character_)
    ) %>%
    transmute(
      Analyte = PolymerLabel,
      `Term detection`,
      `Preterm detection`,
      `Term burden, mean ± SD`,
      `Preterm burden, mean ± SD`,
      `Term burden, median (IQR)`,
      `Preterm burden, median (IQR)`,
      `p-value` = p_fmt,
      `BH q-value` = q_fmt,
      `Sparse detection` = ifelse(sparse_detection, "Yes", "No")
    ) %>%
    gt::gt() %>%
    gt::tab_header(
      title = gt::md("**Table 2. Detection and Placental MNP Burden by Delivery Class**")
    ) %>%
    gt::tab_source_note(gt::md("Burden is reported as µg/g tissue. BH q-values are adjusted across polymers for the univariate comparison.")) %>%
    gt::cols_align(align = "left", columns = Analyte) %>%
    gt::cols_align(align = "center", columns = -Analyte) %>%
    gt::tab_options(table.font.size = gt::px(11))
}

smd_cont <- function(x, g) {
  x <- safe_num(x)
  g <- droplevels(factor(g))
  if (nlevels(g) != 2) return(NA_real_)
  lv <- levels(g)
  x1 <- x[g == lv[1]]
  x2 <- x[g == lv[2]]
  x1 <- x1[is.finite(x1)]
  x2 <- x2[is.finite(x2)]
  if (!length(x1) || !length(x2)) return(NA_real_)
  sp <- sqrt((sd(x1)^2 + sd(x2)^2) / 2)
  if (!is.finite(sp) || sp == 0) return(NA_real_)
  (mean(x2) - mean(x1)) / sp
}

cat_summary <- function(x) {
  x <- as.character(x)
  x <- x[!is.na(x) & x != ""]
  if (!length(x)) return(NA_character_)
  tab <- sort(table(x), decreasing = TRUE)
  paste(sprintf("%s: %d (%.1f%%)", names(tab), as.integer(tab), 100 * as.integer(tab) / sum(tab)), collapse = "; ")
}

cont_summary <- function(x) {
  x <- safe_num(x)
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_character_)
  paste0(median_iqr(x, digits = 2), "; mean ", mean_sd(x, digits = 2))
}

make_table1_cont <- function(df, var, label) {
  x <- df[[var]]
  g <- df$Delivery_class
  tibble(
    variable = label,
    variable_name = var,
    type = "continuous",
    level = NA_character_,
    Overall = cont_summary(x),
    Term = cont_summary(x[g == "Term"]),
    Preterm = cont_summary(x[g == "Preterm"]),
    p_value = p_cont(x, g),
    SMD = smd_cont(x, g),
    n_nonmissing = sum(is.finite(safe_num(x)))
  )
}

make_table1_cat <- function(df, var, label) {
  x <- df[[var]]
  g <- df$Delivery_class
  tibble(
    variable = label,
    variable_name = var,
    type = "categorical",
    level = NA_character_,
    Overall = cat_summary(x),
    Term = cat_summary(x[g == "Term"]),
    Preterm = cat_summary(x[g == "Preterm"]),
    p_value = p_cat(x, g),
    SMD = NA_real_,
    n_nonmissing = sum(!is.na(x) & as.character(x) != "")
  )
}

summarise_exposure_one <- function(df, col, label, scale_label) {
  x <- safe_num(df[[col]])
  g <- df$PTB01
  d <- tibble(value = x, PTB01 = g) %>% filter(is.finite(value), !is.na(PTB01))

  if (nrow(d) == 0) {
    return(tibble(
      PolymerLabel = label, ExposureCol = col, Scale = scale_label,
      N_all = 0L, N_term = 0L, N_preterm = 0L,
      median_term = NA_real_, median_preterm = NA_real_,
      mean_term = NA_real_, mean_preterm = NA_real_,
      delta_median_preterm_minus_term = NA_real_,
      test_used = NA_character_, p_value = NA_real_, q_BH = NA_real_
    ))
  }

  d0 <- d %>% filter(PTB01 == 0)
  d1 <- d %>% filter(PTB01 == 1)
  p <- if (nrow(d0) > 0 && nrow(d1) > 0) {
    suppressWarnings(wilcox.test(value ~ PTB01, data = d, exact = FALSE)$p.value)
  } else {
    NA_real_
  }

  tibble(
    PolymerLabel = label,
    ExposureCol = col,
    Scale = scale_label,
    N_all = nrow(d),
    N_term = nrow(d0),
    N_preterm = nrow(d1),
    median_all = median(d$value, na.rm = TRUE),
    q1_all = quantile(d$value, 0.25, na.rm = TRUE, names = FALSE),
    q3_all = quantile(d$value, 0.75, na.rm = TRUE, names = FALSE),
    mean_all = mean(d$value, na.rm = TRUE),
    sd_all = sd(d$value, na.rm = TRUE),
    median_term = median(d0$value, na.rm = TRUE),
    median_preterm = median(d1$value, na.rm = TRUE),
    mean_term = mean(d0$value, na.rm = TRUE),
    mean_preterm = mean(d1$value, na.rm = TRUE),
    delta_median_preterm_minus_term = median_preterm - median_term,
    test_used = "Wilcoxon rank-sum",
    p_value = p
  )
}

# ----------------------------
# Load / validate
# ----------------------------
df0 <- read_model_table(IN_CSV)

df <- df0 %>%
  filter(!is.na(Delivery_class))

if (!all(c("Term", "Preterm") %in% unique(as.character(df$Delivery_class)))) {
  stop("Need both Term and Preterm after excluding missing Delivery_class.", call. = FALSE)
}

ug_cols <- legacy_col_ug_g(POLYMERS)
log_cols <- legacy_col_log10(POLYMERS)
missing_ug <- setdiff(ug_cols, names(df))
missing_log <- setdiff(log_cols, names(df))
if (length(missing_ug) > 0) stop("Missing ug/g columns: ", paste(missing_ug, collapse = ", "), call. = FALSE)
if (length(missing_log) > 0) stop("Missing log10 columns: ", paste(missing_log, collapse = ", "), call. = FALSE)

# ----------------------------
# Table 1
# ----------------------------
cont_vars <- c(
  GA = "Gestational age at delivery, weeks",
  Maternal_age = "Maternal age, years",
  BMI = "BMI",
  Prepregnancy_BMI = "Prepregnancy BMI",
  SDI_score = "Social Deprivation Index",
  SVI_score = "Social Vulnerability Index",
  Birth_weight_g = "Birth weight, g",
  Gravida = "Gravida"
)
cont_vars <- cont_vars[names(cont_vars) %in% names(df)]

cat_vars <- c(
  Hospital = "Hospital",
  Ethnicity_Race = "Race and ethnicity",
  Infant_gender = "Infant sex",
  Twins = "Twin pregnancy",
  Prior_PTB = "Prior preterm birth",
  `Maternal smoking` = "Maternal smoking",
  Hypertension = "Hypertension",
  `Preeclampsia/gestational hypertension` = "Preeclampsia/gestational hypertension",
  Gestational_diabetes = "Gestational diabetes",
  Chronic_diabetes = "Chronic diabetes",
  PPROM = "PPROM/PROM",
  Delivery_route = "Delivery route",
  Type_of_Labor = "Type of labor",
  Delivery_classification = "Delivery classification",
  Chorioamnionitis = "Chorioamnionitis",
  IUGR = "IUGR",
  Oligohydramnios = "Oligohydramnios",
  Abruption = "Abruption",
  Substance_use = "Substance use"
)
cat_vars <- cat_vars[names(cat_vars) %in% names(df)]

table1 <- bind_rows(
  purrr::imap_dfr(cont_vars, ~ make_table1_cont(df, .y, .x)),
  purrr::imap_dfr(cat_vars, ~ make_table1_cat(df, .y, .x))
) %>%
  mutate(
    p_fmt = pretty_p(p_value),
    SMD = round(SMD, 3)
  )

write_csv_clean(table1, OUT_TABLE1)
save_gt_html_png(make_table1_gt(table1), OUT_TABLE1_HTML, OUT_TABLE1_PNG)

# ----------------------------
# Table 2: detection + burden by PTB
# ----------------------------
det_long <- df %>%
  select(SampleID, Delivery_class, all_of(ug_cols)) %>%
  pivot_longer(
    cols = all_of(ug_cols),
    names_to = "ExposureCol",
    values_to = "mnp_ug_g"
  ) %>%
  mutate(
    Polymer = str_remove(ExposureCol, "_ug_per_g$"),
    PolymerLabel = paste0(POLYMER_FULL_NAMES[Polymer], " (", Polymer, ")"),
    detected = is.finite(mnp_ug_g) & mnp_ug_g > 0
  )

det_summary <- det_long %>%
  group_by(Polymer, PolymerLabel, Delivery_class) %>%
  summarise(
    N = n(),
    n_nonmissing = sum(!is.na(mnp_ug_g)),
    n_detected = sum(detected, na.rm = TRUE),
    pct_detected = 100 * n_detected / N,
    mean_ug_g = mean(mnp_ug_g, na.rm = TRUE),
    sd_ug_g = sd(mnp_ug_g, na.rm = TRUE),
    median_ug_g = median(mnp_ug_g, na.rm = TRUE),
    q1_ug_g = quantile(mnp_ug_g, 0.25, na.rm = TRUE, names = FALSE),
    q3_ug_g = quantile(mnp_ug_g, 0.75, na.rm = TRUE, names = FALSE),
    .groups = "drop"
  )

burden_tests <- purrr::map_dfr(POLYMERS, function(poly) {
  col <- legacy_col_ug_g(poly)
  summarise_exposure_one(df, col, paste0(POLYMER_FULL_NAMES[[poly]], " (", poly, ")"), "ug/g")
}) %>%
  mutate(q_BH = p.adjust(p_value, method = "BH"))

table2_poly <- det_summary %>%
  pivot_wider(
    id_cols = c(Polymer, PolymerLabel),
    names_from = Delivery_class,
    values_from = c(N, n_detected, pct_detected, mean_ug_g, sd_ug_g, median_ug_g, q1_ug_g, q3_ug_g),
    names_sep = "_"
  ) %>%
  left_join(
    burden_tests %>% select(PolymerLabel, p_value, q_BH),
    by = "PolymerLabel"
  ) %>%
  mutate(
    row_type = "Polymer",
    sparse_detection = Polymer == "PS" | pmin(coalesce(n_detected_Term, 0), coalesce(n_detected_Preterm, 0)) < MIN_DETECTED_N_FOR_CONTINUOUS_MODEL
  )

cum_test <- summarise_exposure_one(df, "Cumulative_MNPs_ug_per_g", "Cumulative MNPs", "ug/g") %>%
  mutate(q_BH = p_value, row_type = "Cumulative", sparse_detection = FALSE)

cum_row <- df %>%
  select(SampleID, Delivery_class, Cumulative_MNPs_ug_per_g) %>%
  group_by(Delivery_class) %>%
  summarise(
    N = n(),
    n_detected = sum(is.finite(Cumulative_MNPs_ug_per_g) & Cumulative_MNPs_ug_per_g > 0),
    pct_detected = 100 * n_detected / N,
    mean_ug_g = mean(Cumulative_MNPs_ug_per_g, na.rm = TRUE),
    sd_ug_g = sd(Cumulative_MNPs_ug_per_g, na.rm = TRUE),
    median_ug_g = median(Cumulative_MNPs_ug_per_g, na.rm = TRUE),
    q1_ug_g = quantile(Cumulative_MNPs_ug_per_g, 0.25, na.rm = TRUE, names = FALSE),
    q3_ug_g = quantile(Cumulative_MNPs_ug_per_g, 0.75, na.rm = TRUE, names = FALSE),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = Delivery_class,
    values_from = c(N, n_detected, pct_detected, mean_ug_g, sd_ug_g, median_ug_g, q1_ug_g, q3_ug_g),
    names_sep = "_"
  ) %>%
  mutate(
    Polymer = "Cumulative_MNPs",
    PolymerLabel = "Cumulative MNPs"
  ) %>%
  left_join(cum_test %>% select(PolymerLabel, p_value, q_BH, row_type, sparse_detection), by = "PolymerLabel")

table2 <- bind_rows(
  cum_row,
  table2_poly
) %>%
  mutate(
    p_fmt = pretty_p(p_value),
    q_fmt = pretty_p(q_BH)
  ) %>%
  select(row_type, Polymer, PolymerLabel, everything())

write_csv_clean(table2, OUT_TABLE2)
save_gt_html_png(make_table2_gt(table2), OUT_TABLE2_HTML, OUT_TABLE2_PNG)

# ----------------------------
# Univariate summaries: ug/g and log10
# ----------------------------
uni_ug <- purrr::map_dfr(POLYMERS, function(poly) {
  summarise_exposure_one(df, legacy_col_ug_g(poly), poly, "ug/g")
}) %>%
  mutate(q_BH = p.adjust(p_value, method = "BH"),
         sparse_detection = PolymerLabel == "PS")

uni_log <- purrr::map_dfr(POLYMERS, function(poly) {
  summarise_exposure_one(df, legacy_col_log10(poly), poly, "log10(ug/g + 1)")
}) %>%
  mutate(
    q_BH = p.adjust(p_value, method = "BH"),
    polymer_family = vapply(PolymerLabel, poly_family, character(1)),
    q_BH_within_family = ave(p_value, polymer_family, FUN = function(z) p.adjust(z, method = "BH")),
    sparse_detection = PolymerLabel == "PS"
  )

write_csv_clean(uni_ug, OUT_UNI_UG)
write_csv_clean(uni_log, OUT_UNI_LOG)
write_csv_clean(
  uni_log %>% select(PolymerLabel, polymer_family, p_value, q_BH, q_BH_within_family, sparse_detection),
  file.path(DIR_TABLES, "univariate__log10ugpg_summary_term_vs_ptb_with_family_q_polished.csv")
)

# ----------------------------
# GA Spearman correlations
# ----------------------------
ga_tbl <- purrr::map_dfr(c("Cumulative_MNPs_log10ugpg", log_cols), function(col) {
  if (!col %in% names(df)) return(NULL)
  d <- df %>%
    transmute(
      GA = safe_num(GA),
      exposure = safe_num(.data[[col]])
    ) %>%
    filter(is.finite(GA), is.finite(exposure))

  if (nrow(d) < 5 || sd(d$exposure) == 0) {
    return(tibble(ExposureCol = col, PolymerLabel = str_remove(col, "_log10ugpg$"), N = nrow(d), rho = NA_real_, p_value = NA_real_))
  }

  ct <- suppressWarnings(cor.test(d$GA, d$exposure, method = "spearman", exact = FALSE))
  tibble(
    ExposureCol = col,
    PolymerLabel = ifelse(col == "Cumulative_MNPs_log10ugpg", "Cumulative MNP burden", str_remove(col, "_log10ugpg$")),
    N = nrow(d),
    rho = unname(ct$estimate),
    p_value = unname(ct$p.value)
  )
}) %>%
  mutate(q_BH = p.adjust(p_value, method = "BH")) %>%
  arrange(p_value)

write_csv_clean(ga_tbl, OUT_GA)

# ----------------------------
# Figures
# ----------------------------
# Fig 2: landscape detection + burden
landscape <- det_long %>%
  group_by(Polymer) %>%
  summarise(
    PolymerLabel = first(PolymerLabel),
    n = n(),
    n_detected = sum(detected, na.rm = TRUE),
    pct_detected = 100 * n_detected / n,
    median_ug_g = median(mnp_ug_g, na.rm = TRUE),
    mean_ug_g = mean(mnp_ug_g, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Polymer = factor(Polymer, levels = POLYMERS),
    PolymerLabel = factor(PolymerLabel, levels = paste0(POLYMER_FULL_NAMES[POLYMERS], " (", POLYMERS, ")"))
  )

p_landscape <- ggplot(landscape, aes(x = PolymerLabel)) +
  geom_col(aes(y = pct_detected), alpha = 0.35, width = 0.75) +
  geom_point(aes(y = rescale(log10(mean_ug_g + 1), to = c(0, 100))), size = 2.5) +
  coord_flip() +
  labs(
    title = "Placental MNP Detection and Mean Burden",
    subtitle = "Bars show detection percent; points show scaled log10 mean ug/g",
    x = NULL,
    y = "Detection percent"
  ) +
  theme_mnp(base_size = BASE_SIZE) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"))

ggsave_mnp(OUT_FIG2, p_landscape, width = 8, height = 6.5, dpi = 300)

# Fig 3 burden by PTB
ptb_stat_labels <- uni_log %>%
  mutate(
    Polymer = factor(PolymerLabel, levels = POLYMERS),
    PolymerLabelPretty = paste0(POLYMER_FULL_NAMES[PolymerLabel], " (", PolymerLabel, ")"),
    q_label = ifelse(
      is.na(q_BH),
      "q = NA",
      ifelse(q_BH < 0.001, "q < 0.001", paste0("q = ", formatC(q_BH, format = "f", digits = 3)))
    ),
    q_label = ifelse(sparse_detection, paste0(q_label, "; sparse"), q_label)
  ) %>%
  select(Polymer = PolymerLabel, q_label)

plot_long <- df %>%
  select(SampleID, Delivery_class, all_of(ug_cols)) %>%
  pivot_longer(all_of(ug_cols), names_to = "ExposureCol", values_to = "mnp_ug_g") %>%
  mutate(
    Polymer = str_remove(ExposureCol, "_ug_per_g$"),
    Polymer = factor(Polymer, levels = POLYMERS),
    PolymerLabel = factor(paste0(POLYMER_FULL_NAMES[as.character(Polymer)], " (", as.character(Polymer), ")"),
                          levels = paste0(POLYMER_FULL_NAMES[POLYMERS], " (", POLYMERS, ")"))
  )

ptb_stat_labels <- plot_long %>%
  group_by(Polymer, PolymerLabel) %>%
  summarise(
    y_lab = max(mnp_ug_g, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    y_lab = ifelse(is.finite(y_lab), y_lab, 1),
    y_lab = y_lab * 1.12,
    Polymer_chr = as.character(Polymer)
  ) %>%
  left_join(ptb_stat_labels %>% mutate(Polymer_chr = as.character(Polymer)), by = "Polymer_chr") %>%
  select(-Polymer_chr, -Polymer.y) %>%
  rename(Polymer = Polymer.x)

p_burden <- ggplot(plot_long, aes(x = Delivery_class, y = mnp_ug_g)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.25, width = 0.55) +
  geom_jitter(width = 0.15, alpha = 0.35, size = 0.8) +
  geom_text(
    data = ptb_stat_labels,
    aes(x = 1.5, y = y_lab, label = q_label),
    inherit.aes = FALSE,
    size = 2.8
  ) +
  facet_wrap(~ PolymerLabel, scales = "free_y", ncol = 4) +
  scale_y_continuous(
    trans = scales::pseudo_log_trans(base = 10),
    labels = label_number(),
    expand = expansion(mult = c(0.04, 0.22))
  ) +
  labs(
    title = "Polymer Burden by Delivery Class",
    x = NULL,
    y = expression("MNP burden ("*mu*"g/g tissue)")
  ) +
  theme_mnp(base_size = 10.5) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 30, hjust = 1),
    plot.title = element_text(face = "bold")
  )

ggsave_mnp(OUT_FIG3_BURDEN, p_burden, width = 12, height = 8.6, dpi = 300)

# Volcano
volc <- uni_log %>%
  mutate(
    Polymer = PolymerLabel,
    neglog10q = -log10(q_BH),
    sig = q_BH < 0.05,
    sparse_detection = Polymer == "PS"
  )

if (requireNamespace("ggrepel", quietly = TRUE)) {
  p_volc <- ggplot(volc, aes(x = delta_median_preterm_minus_term, y = neglog10q, label = Polymer)) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.35) +
    geom_vline(xintercept = 0, linetype = "dotted", linewidth = 0.35) +
    geom_point(aes(shape = sparse_detection), size = 2.8, alpha = 0.80) +
    ggrepel::geom_text_repel(size = 3, max.overlaps = Inf, seed = 1234) +
    labs(
      title = "Univariate PTB Differences in Polymer Burden",
      x = "Median difference: Preterm - Term, log10(ug/g + 1)",
      y = "-log10(BH q-value)",
      shape = "Sparse detection"
    ) +
    theme_mnp(base_size = BASE_SIZE) +
    theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"))
} else {
  p_volc <- ggplot(volc, aes(x = delta_median_preterm_minus_term, y = neglog10q)) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.35) +
    geom_vline(xintercept = 0, linetype = "dotted", linewidth = 0.35) +
    geom_point(aes(shape = sparse_detection), size = 2.8, alpha = 0.80) +
    geom_text(aes(label = Polymer), size = 3, vjust = -0.5, check_overlap = TRUE) +
    labs(
      title = "Univariate PTB Differences in Polymer Burden",
      x = "Median difference: Preterm - Term, log10(ug/g + 1)",
      y = "-log10(BH q-value)",
      shape = "Sparse detection"
    ) +
    theme_mnp(base_size = BASE_SIZE) +
    theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"))
}

ggsave_mnp(OUT_FIG3_VOLCANO, p_volc, width = 7.5, height = 5.5, dpi = 300)

# GA lollipop
ga_plot <- ga_tbl %>%
  filter(ExposureCol != "Cumulative_MNPs_log10ugpg") %>%
  mutate(
    PolymerLabel = factor(PolymerLabel, levels = PolymerLabel[order(rho, na.last = TRUE)]),
    sig = q_BH < 0.05
  )

p_ga <- ggplot(ga_plot, aes(x = rho, y = PolymerLabel)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.35) +
  geom_segment(aes(x = 0, xend = rho, yend = PolymerLabel), linewidth = 0.6) +
  geom_point(aes(shape = sig), size = 2.8) +
  labs(
    title = "Spearman Correlation Between GA and Polymer Burden",
    x = "Spearman rho",
    y = NULL,
    shape = "BH q < 0.05"
  ) +
  theme_mnp(base_size = BASE_SIZE) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"))

ggsave_mnp(OUT_FIG_GA, p_ga, width = 7.2, height = 5.5, dpi = 300)

# ----------------------------
# Log
# ----------------------------
run_log <- c(
  "04_descriptive_and_univariate.R",
  paste0("PROJECT_ROOT: ", PROJECT_ROOT),
  paste0("Input rows: ", nrow(df0)),
  paste0("Rows with Delivery_class: ", nrow(df)),
  "",
  "Delivery_class counts:",
  capture.output(print(table(df$Delivery_class, useNA = "ifany"))),
  "",
  "Wrote:",
  OUT_TABLE1,
  OUT_TABLE1_HTML,
  OUT_TABLE1_PNG,
  OUT_TABLE2,
  OUT_TABLE2_HTML,
  OUT_TABLE2_PNG,
  OUT_UNI_UG,
  OUT_UNI_LOG,
  file.path(DIR_TABLES, "univariate__log10ugpg_summary_term_vs_ptb_with_family_q_polished.csv"),
  OUT_GA,
  OUT_FIG2,
  OUT_FIG3_BURDEN,
  OUT_FIG3_VOLCANO,
  OUT_FIG_GA
)

cat(paste(run_log, collapse = "\n"), file = OUT_LOG)
message(paste(run_log, collapse = "\n"))
