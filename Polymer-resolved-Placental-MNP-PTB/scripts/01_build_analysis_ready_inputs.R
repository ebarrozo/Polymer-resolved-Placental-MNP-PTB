#!/usr/bin/env Rscript
# ============================================================
# 01_build_analysis_ready_inputs.R
# Build legacy-compatible analysis inputs from cleaned model-ready files
# ============================================================
# This is the bridge script. It consumes the cleaned files generated from:
#   placenta_model_ready_wide.csv
#   controls_blank_ready_wide.csv
#   controls_blank_ready_long.csv
#
# It writes the filenames expected by existing downstream scripts:
#   outputs/tables/mnp_ptb_cohort__analysis_ready_model_table.csv
#   outputs/tables/mnp_all_records__qc_flags_and_outlier_flags_wide.csv
#   outputs/tables/mnp_controls__lob_lod_loq_maxblank_sdmethod_by_polymer.csv
#   outputs/tables/mnp_controls__blank_summary_stats_by_polymer.csv
#   outputs/tables/mnp_ptb_cohort__polymer_detection_summary.csv
#   outputs/logs/01_build_analysis_ready_inputs__run_log.txt
#
# Design choices:
#   - PE blank mass is set to 0 when absent from control blank export.
#   - PS remains in descriptive outputs but is flagged by detection count.
#   - Legacy columns are generated without deleting cleaned canonical columns.
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
  library(stringr)
})

source("00_config.R")
source(file.path("R", "utils_mnp.R"))

# ----------------------------
# Inputs
# ----------------------------
placenta_wide_path <- file.path(DIR_MODEL_READY, "placenta_model_ready_wide.csv")
placenta_long_path <- file.path(DIR_MODEL_READY, "placenta_model_ready_long.csv")
controls_wide_path <- file.path(DIR_MODEL_READY, "controls_blank_ready_wide.csv")
controls_long_path <- file.path(DIR_MODEL_READY, "controls_blank_ready_long.csv")

stopifnot(file.exists(placenta_wide_path))
stopifnot(file.exists(controls_wide_path))
stopifnot(file.exists(controls_long_path))

plac0 <- readr::read_csv(placenta_wide_path, show_col_types = FALSE, guess_max = 5000)
ctrl0 <- readr::read_csv(controls_wide_path, show_col_types = FALSE, guess_max = 5000)
ctrl_long0 <- readr::read_csv(controls_long_path, show_col_types = FALSE, guess_max = 5000)

# ----------------------------
# Optional PeriBank metadata + SDI/SVI enrichment
# ----------------------------
peribank_aug <- NULL
if (exists("PERIBANK_XLSX") && file.exists(PERIBANK_XLSX)) {
  if (!requireNamespace("readxl", quietly = TRUE)) {
    stop("PERIBANK_XLSX exists but package 'readxl' is not installed. Run install.packages('readxl').", call. = FALSE)
  }

  message("Reading PeriBank metadata: ", PERIBANK_XLSX)

  peribank_raw <- readxl::read_excel(PERIBANK_XLSX, guess_max = 10000) %>%
    as_tibble(.name_repair = "unique") %>%
    mutate(PregnancyID = as.character(PregnancyID))

  if (!"PregnancyID" %in% names(peribank_raw)) {
    stop("PeriBank workbook is missing PregnancyID.", call. = FALSE)
  }

  # Deduplicate conservatively, keeping first row per PregnancyID.
  # The file name says wodups, but this keeps the pipeline safe.
  peribank_aug <- peribank_raw %>%
    group_by(PregnancyID) %>%
    slice(1) %>%
    ungroup() %>%
    mutate(
      home_zip_2nd_3rd_trim = pad_zip(pull_existing(., c(
        "Home zip code at 2nd/3rd trim",
        "Home zip code at 2nd/3rd trim...12",
        "Home_zip_code_at_2nd_3rd_trim"
      ), default = NA_character_)),

      PP_BMI = safe_num(pull_existing(., c("PP BMI", "PP_BMI", "Prepregnancy BMI"), default = NA_real_)),
      Parity_raw = safe_chr(pull_existing(., c("Parity"), default = NA_character_)),
      Prior_PTB = derive_prior_ptb_from_parity(Parity_raw),

      Hospital_from_peribank = collapse_hospital(pull_existing(., c("Hospital", "Hospital...11", "Hospital...68"), default = NA_character_)),
      Ethnicity_Race = derive_ethnicity_race(.),

      Maternal_smoking = safe_chr(pull_existing(., c("Maternal smoking"), default = NA_character_)),
      Substance_use = safe_chr(pull_existing(., c("Substance use (alcohol or other)"), default = NA_character_)),

      Gravida = safe_num(pull_existing(., c("Gravida"), default = NA_real_)),
      PP_weight = safe_num(pull_existing(., c("PP weight"), default = NA_real_)),
      PP_height_ft = safe_num(pull_existing(., c("PP height (ft)"), default = NA_real_)),
      PP_height_in = safe_num(pull_existing(., c("PP height (in)"), default = NA_real_)),

      Type_I_diabetes = as_01(pull_existing(., c("Type I diabetes"), default = NA_character_)),
      Type_II_diabetes = as_01(pull_existing(., c("Type II diabetes"), default = NA_character_)),
      Prior_GDMA1 = as_01(pull_existing(., c("Prior GDMA1"), default = NA_character_)),
      Prior_GDMA2 = as_01(pull_existing(., c("Prior GDMA2"), default = NA_character_)),

      H1N1_vaccine = as_01(pull_existing(., c("H1N1 vaccine"), default = NA_character_)),
      HPV_vaccine = as_01(pull_existing(., c("HPV vaccine"), default = NA_character_)),
      Influenza_vaccine = as_01(pull_existing(., c("Influenza (seasonal) vaccine"), default = NA_character_)),

      Chlamydia = as_01(pull_existing(., c("Chlamydia"), default = NA_character_)),
      HIV_1st_OB = as_01(pull_existing(., c("HIV at 1st OB"), default = NA_character_)),
      CMV = as_01(pull_existing(., c("CMV"), default = NA_character_)),
      HSV_1_2 = as_01(pull_existing(., c("HSV 1/2"), default = NA_character_)),
      HepB = as_01(pull_existing(., c("HepB"), default = NA_character_)),
      GBS = as_01(pull_existing(., c("GBS"), default = NA_character_)),
      Rubella = as_01(pull_existing(., c("Rubella"), default = NA_character_)),
      Hx_COVID19 = as_01(pull_existing(., c("Hx of COVID-19"), default = NA_character_))
    ) %>%
    select(
      PregnancyID,
      home_zip_2nd_3rd_trim,
      PP_BMI,
      Parity_raw,
      Prior_PTB,
      Hospital_from_peribank,
      Ethnicity_Race,
      Maternal_smoking,
      Substance_use,
      Gravida,
      PP_weight,
      PP_height_ft,
      PP_height_in,
      Type_I_diabetes,
      Type_II_diabetes,
      Prior_GDMA1,
      Prior_GDMA2,
      H1N1_vaccine,
      HPV_vaccine,
      Influenza_vaccine,
      Chlamydia,
      HIV_1st_OB,
      CMV,
      HSV_1_2,
      HepB,
      GBS,
      Rubella,
      Hx_COVID19
    )

  # Join SDI/SVI by home ZIP if files are present.
  if (exists("SVI_PATH") && file.exists(SVI_PATH)) {
    svi_raw <- readr::read_csv(SVI_PATH, guess_max = 1000000, show_col_types = FALSE) %>%
      as_tibble(.name_repair = "unique")

    svi_zip_col <- first_existing_col(svi_raw, c("Zip", "zip", "ZIP", "ZCTA", "zcta", "ZCTA5", "zcta5"))
    if (is.na(svi_zip_col)) {
      warning("SVI file found but no Zip-like column detected. Skipping SVI join.")
    } else {
      svi <- svi_raw %>%
        mutate(home_zip_2nd_3rd_trim = pad_zip(.data[[svi_zip_col]])) %>%
        group_by(home_zip_2nd_3rd_trim) %>%
        slice(1) %>%
        ungroup() %>%
        rename_with(~ paste0("svi_", .x), .cols = -home_zip_2nd_3rd_trim)

      peribank_aug <- peribank_aug %>% left_join(svi, by = "home_zip_2nd_3rd_trim")
    }
  } else {
    message("SVI file not found; skipping SVI join: ", SVI_PATH)
  }

  if (exists("SDI_PATH") && file.exists(SDI_PATH)) {
    sdi_raw <- readr::read_csv(SDI_PATH, guess_max = 1000000, show_col_types = FALSE) %>%
      as_tibble(.name_repair = "unique")

    sdi_zip_col <- first_existing_col(sdi_raw, c("Zip", "zip", "ZIP", "ZCTA", "zcta", "ZCTA5", "zcta5"))
    if (is.na(sdi_zip_col)) {
      warning("SDI file found but no Zip-like column detected. Skipping SDI join.")
    } else {
      sdi <- sdi_raw %>%
        mutate(home_zip_2nd_3rd_trim = pad_zip(.data[[sdi_zip_col]])) %>%
        group_by(home_zip_2nd_3rd_trim) %>%
        slice(1) %>%
        ungroup() %>%
        rename_with(~ paste0("sdi_", .x), .cols = -home_zip_2nd_3rd_trim)

      peribank_aug <- peribank_aug %>% left_join(sdi, by = "home_zip_2nd_3rd_trim")
    }
  } else {
    message("SDI file not found; skipping SDI join: ", SDI_PATH)
  }

  # Identify likely SDI/SVI score columns for model compatibility.
  sdi_cols <- names(peribank_aug)[str_detect(names(peribank_aug), regex("^sdi_.*(score|sdi|rank|index)$|^sdi_score$", ignore_case = TRUE))]
  svi_cols <- names(peribank_aug)[str_detect(names(peribank_aug), regex("^svi_.*(rpl_themes|rpl_theme|overall|svi|score|rank|index)$", ignore_case = TRUE))]

  peribank_aug <- peribank_aug %>%
    mutate(
      SDI_score_from_zip = if (length(sdi_cols) > 0) safe_num(.data[[sdi_cols[[1]]]]) else NA_real_,
      SVI_score_from_zip = if (length(svi_cols) > 0) safe_num(.data[[svi_cols[[1]]]]) else NA_real_
    )

  write_csv_clean(peribank_aug, file.path(DIR_DIAG, "peribank_metadata_with_sdi_svi_joined.csv"))

} else {
  message("PeriBank workbook not found; metadata enrichment skipped: ", PERIBANK_XLSX)
}

# ----------------------------
# Build placenta legacy-compatible model table
# ----------------------------
plac <- plac0 %>%
  mutate(
    PregnancyID = as.character(PregnancyID),
    SampleID = as.character(PregnancyID),
    `Pregnancy ID` = as.character(PregnancyID),
    Dataset = "Placenta",
    Sample_Type = "Placental Samples",

    # Core outcome/covariates expected by legacy scripts
    GA = safe_num(ga_delivery_weeks),
    GA_weeks = safe_num(ga_delivery_weeks),
    GA_days = safe_num(ga_delivery_days),

    Preterm = case_when(
      !is.na(preterm_birth) ~ as.logical(preterm_birth == 1),
      is.finite(GA) ~ GA < 37,
      TRUE ~ NA
    ),
    PTB01 = as.integer(Preterm),
    Delivery_class = factor(ifelse(Preterm, "Preterm", "Term"), levels = c("Term","Preterm")),

    Maternal_age = safe_num(maternal_age_years),
    Maternal_age_years = safe_num(maternal_age_years),

    BMI = safe_num(bmi),
    Prepregnancy_BMI = safe_num(bmi),

    Twins01 = as.integer(twin_pregnancy),
    Twins = factor(ifelse(is.na(Twins01), NA_character_, ifelse(Twins01 == 1, "Yes", "No")),
                   levels = c("No","Yes")),

    Infant_gender = safe_chr(fetal_sex),
    Fetal_sex = safe_chr(fetal_sex),

    Birth_weight_g = safe_num(birth_weight_g),

    Delivery_classification = safe_chr(delivery_classification),
    `Delivery classification` = safe_chr(delivery_classification),
    Delivery_route = safe_chr(delivery_route),
    `Delivery route` = safe_chr(delivery_route),
    Type_of_Labor = safe_chr(labor_type),
    `Type of labor` = safe_chr(labor_type),

    PROM = as_yesno_factor(prom),
    PPROM = as_yesno_factor(prom),

    Hypertension = as_yesno_factor(hypertensive_disorder),
    Hypertension01 = as.integer(hypertensive_disorder),

    `Preeclampsia/gestational hypertension` = safe_chr(preeclampsia_or_gest_htn_type),
    Hypertension_type = safe_chr(preeclampsia_or_gest_htn_type),

    Preeclampsia01 = as.integer(str_detect(str_to_lower(Hypertension_type %||% ""), "preeclamp|hellp")),
    Gestational_Hypertension01 = as.integer(str_detect(str_to_lower(Hypertension_type %||% ""), "gestational")),

    Gestational_diabetes = as_yesno_factor(gestational_diabetes),
    `Gestational diabetes` = as_yesno_factor(gestational_diabetes),
    Chronic_diabetes = as_yesno_factor(chronic_diabetes),
    Chorioamnionitis = as_yesno_factor(chorioamnionitis),
    Infection_any = as.integer(infection_any),

    Oligohydramnios = as_yesno_factor(oligohydramnios),
    IUGR = as_yesno_factor(iugr),
    Abruption = as_yesno_factor(abruption),

    # Placeholders retained for older models. Fill later if the metadata source is restored.
    Prior_PTB = NA_integer_,
    `Maternal smoking` = NA_character_,
    SDI_score = NA_real_,
    Hospital = NA_character_,
    Ethnicity_Race = NA_character_
  )

# Attach PeriBank/ZIP-derived metadata when available.
if (!is.null(peribank_aug)) {
  plac <- plac %>%
    mutate(PregnancyID = as.character(PregnancyID)) %>%
    left_join(peribank_aug %>% mutate(PregnancyID = as.character(PregnancyID)), by = "PregnancyID") %>%
    mutate(
      home_zip_2nd_3rd_trim = coalesce(home_zip_2nd_3rd_trim, NA_character_),

      # Fill placeholders / weak fields from PeriBank.
      Prior_PTB = coalesce(as.integer(Prior_PTB.y), as.integer(Prior_PTB.x)),
      `Maternal smoking` = coalesce(Maternal_smoking, `Maternal smoking`),
      SDI_score = coalesce(SDI_score_from_zip, safe_num(SDI_score)),
      SVI_score = SVI_score_from_zip,
      Hospital = coalesce(Hospital_from_peribank, Hospital),
      Ethnicity_Race = coalesce(Ethnicity_Race.y, Ethnicity_Race.x),
      BMI = coalesce(PP_BMI, BMI),
      Prepregnancy_BMI = coalesce(PP_BMI, Prepregnancy_BMI)
    ) %>%
    select(
      -any_of(c(
        "Prior_PTB.x", "Prior_PTB.y",
        "Ethnicity_Race.x", "Ethnicity_Race.y"
      ))
    )
} else {
  plac <- plac %>%
    mutate(
      home_zip_2nd_3rd_trim = NA_character_,
      SVI_score = NA_real_
    )
}

# Add legacy polymer concentration + log10 columns.
for (poly in POLYMERS) {
  src <- poly_col_ug_g(poly)
  out_ug <- legacy_col_ug_g(poly)
  out_log <- legacy_col_log10(poly)

  if (src %in% names(plac)) {
    plac[[out_ug]] <- safe_num(plac[[src]])
  } else {
    plac[[out_ug]] <- NA_real_
  }

  plac[[out_log]] <- safe_log10p1(plac[[out_ug]])
}

legacy_ug_cols <- legacy_col_ug_g(POLYMERS)
legacy_log_cols <- legacy_col_log10(POLYMERS)

plac <- plac %>%
  mutate(
    Cumulative_MNPs_ug_per_g = rowSums(across(all_of(legacy_ug_cols)), na.rm = TRUE),
    Cumulative_MNPs_log10ugpg = log10(Cumulative_MNPs_ug_per_g + 1),
    Cumulative_MNPs_ug = rowSums(across(all_of(poly_col_sum_ug(POLYMERS))), na.rm = TRUE),
    Cumulative_MNPs_log10ug_per_sample = log10(Cumulative_MNPs_ug + 1)
  )

# ----------------------------
# Detection and sparse-polymer audit
# ----------------------------
detection_summary <- purrr::map_dfr(POLYMERS, function(poly) {
  col <- legacy_col_ug_g(poly)
  x <- plac[[col]]
  tibble(
    Polymer = poly,
    column = col,
    n_total = nrow(plac),
    n_nonmissing = sum(!is.na(x)),
    n_detected_gt0 = sum(is.finite(x) & x > 0),
    detection_prop = n_detected_gt0 / n_total,
    continuous_model_ok = n_detected_gt0 >= MIN_DETECTED_N_FOR_CONTINUOUS_MODEL
  )
})

# ----------------------------
# Build control blank thresholds
# ----------------------------
# controls_long0 has no PE blank mass in current export because PE was not detected.
# Add PE = 0 explicitly so blank subtraction/QC is intentional and reproducible.
ctrl_long <- ctrl_long0 %>%
  mutate(
    polymer = toupper(as.character(polymer)),
    blank_ug = safe_num(blank_ug)
  ) %>%
  filter(polymer %in% POLYMERS)

if (FORCE_PE_CONTROL_BLANK_ZERO && !("PE" %in% unique(ctrl_long$polymer))) {
  pe_controls <- ctrl0 %>%
    transmute(
      sample_type,
      control_id,
      control_key,
      contents,
      source,
      tube,
      replicate,
      run_id,
      source_file,
      batch,
      batch_date,
      polymer = "PE",
      blank_ug = 0
    )
  ctrl_long <- bind_rows(ctrl_long, pe_controls)
}

blank_summary <- ctrl_long %>%
  group_by(Analyte = polymer) %>%
  summarise(
    n_controls = n(),
    n_detected_gt0 = sum(is.finite(blank_ug) & blank_ug > 0),
    mean_blank_ug = mean(blank_ug, na.rm = TRUE),
    sd_blank_ug = sd(blank_ug, na.rm = TRUE),
    median_blank_ug = median(blank_ug, na.rm = TRUE),
    q95_blank_ug = quantile(blank_ug, 0.95, na.rm = TRUE, names = FALSE),
    max_blank_ug = max(blank_ug, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    sd_blank_ug = ifelse(is.na(sd_blank_ug), 0, sd_blank_ug),
    LOB_ug = mean_blank_ug + 1.645 * sd_blank_ug,
    LOD_ug = LOB_ug + 1.645 * sd_blank_ug,
    LOQ_ug = mean_blank_ug + 10 * sd_blank_ug,
    threshold_rule = "overall controls; PE forced to 0 if absent from blank export"
  ) %>%
  arrange(match(Analyte, POLYMERS))

thresholds <- blank_summary %>%
  transmute(
    Analyte,
    n_controls,
    n_detected_gt0,
    MeanBlank_ug = mean_blank_ug,
    SDBlank_ug = sd_blank_ug,
    MedianBlank_ug = median_blank_ug,
    Q95Blank_ug = q95_blank_ug,
    MaxBlank_ug = max_blank_ug,
    LOB_ug,
    LOD_ug,
    LOQ_ug,
    SD_method = "sample SD",
    threshold_rule
  )

# ----------------------------
# Add QC flag columns to placenta model table
# ----------------------------
for (poly in POLYMERS) {
  sum_col <- poly_col_sum_ug(poly)
  qc_col <- paste0("QC_flag_", poly)
  thr_lod <- thresholds$LOD_ug[match(poly, thresholds$Analyte)]
  thr_lob <- thresholds$LOB_ug[match(poly, thresholds$Analyte)]
  thr_loq <- thresholds$LOQ_ug[match(poly, thresholds$Analyte)]

  x <- if (sum_col %in% names(plac)) safe_num(plac[[sum_col]]) else rep(NA_real_, nrow(plac))

  plac[[qc_col]] <- case_when(
    is.na(x) ~ "missing",
    x == 0 ~ "not_detected",
    is.finite(thr_loq) & x >= thr_loq ~ "above_LOQ",
    is.finite(thr_lod) & x >= thr_lod ~ "above_LOD_below_LOQ",
    is.finite(thr_lob) & x >= thr_lob ~ "above_LOB_below_LOD",
    TRUE ~ "below_LOB"
  )
}

# ----------------------------
# Build old-style all-records QC table for Fig 1 scripts
# ----------------------------
plac_records <- plac %>%
  transmute(
    RecordID = row_number(),
    SampleID,
    Dataset = "Placenta",
    Sample_Type = "Placental Samples",
    Sample_Weight_g = total_tissue_weight_g,
    Tube = NA_character_,
    Source = NA_character_
  )

# Add absolute polymer mass columns expected by legacy QC/Fig 1 scripts.
# These columns do not exist in plac as PE/SBR/etc, so they must be created
# after transmute rather than selected with across(all_of(POLYMERS)).
for (poly in POLYMERS) {
  mass_col <- poly_col_sum_ug(poly)
  plac_records[[poly]] <- if (mass_col %in% names(plac)) safe_num(plac[[mass_col]]) else NA_real_
}

ctrl_records <- ctrl0 %>%
  mutate(
    RecordID = row_number() + nrow(plac_records),
    SampleID = as.character(control_id),
    Dataset = "Control",
    Sample_Type = "Procedural Control",
    Sample_Weight_g = NA_real_,
    Tube = as.character(tube),
    Source = as.character(source)
  ) %>%
  select(RecordID, SampleID, Dataset, Sample_Type, Sample_Weight_g, Tube, Source, everything())

for (poly in POLYMERS) {
  blank_col <- paste0("blank_", tolower(poly), "_ug")
  if (blank_col %in% names(ctrl_records)) {
    ctrl_records[[poly]] <- safe_num(ctrl_records[[blank_col]])
  } else if (poly == "PE" && FORCE_PE_CONTROL_BLANK_ZERO) {
    ctrl_records[[poly]] <- 0
  } else {
    ctrl_records[[poly]] <- NA_real_
  }
}

ctrl_records <- ctrl_records %>%
  select(RecordID, SampleID, Dataset, Sample_Type, Sample_Weight_g, Tube, Source, all_of(POLYMERS), everything())

qc_all <- bind_rows(
  plac_records %>% select(RecordID, SampleID, Dataset, Sample_Type, Sample_Weight_g, Tube, Source, all_of(POLYMERS)),
  ctrl_records %>% select(RecordID, SampleID, Dataset, Sample_Type, Sample_Weight_g, Tube, Source, all_of(POLYMERS))
) %>%
  mutate(Cumulative_MNPs_ug_per_sample = rowSums(across(all_of(POLYMERS)), na.rm = TRUE))

# Attach QC flags to all records using absolute mass thresholds.
for (poly in POLYMERS) {
  qc_col <- paste0("QC_flag_", poly)
  thr_lod <- thresholds$LOD_ug[match(poly, thresholds$Analyte)]
  thr_lob <- thresholds$LOB_ug[match(poly, thresholds$Analyte)]
  thr_loq <- thresholds$LOQ_ug[match(poly, thresholds$Analyte)]

  x <- safe_num(qc_all[[poly]])
  qc_all[[qc_col]] <- case_when(
    is.na(x) ~ "missing",
    x == 0 ~ "not_detected",
    is.finite(thr_loq) & x >= thr_loq ~ "above_LOQ",
    is.finite(thr_lod) & x >= thr_lod ~ "above_LOD_below_LOQ",
    is.finite(thr_lob) & x >= thr_lob ~ "above_LOB_below_LOD",
    TRUE ~ "below_LOB"
  )

  # Placeholder consensus outlier columns expected by old Fig 1 QC code.
  # Current cleaned-control table has no control outlier flags, so default FALSE.
  qc_all[[paste0("Consensus_Outlier_TRUE_", poly)]] <- FALSE
}

# ----------------------------
# Model-ready CSV column ordering
# ----------------------------
front_cols <- c(
  "SampleID", "PregnancyID", "Pregnancy ID", "Dataset", "Sample_Type",
  "GA", "GA_weeks", "GA_days", "Preterm", "PTB01", "Delivery_class",
  "Maternal_age", "Maternal_age_years", "BMI", "Prepregnancy_BMI",
  "Twins", "Twins01", "Prior_PTB", "SDI_score", "SVI_score", "home_zip_2nd_3rd_trim", "Hospital", "Ethnicity_Race", "Parity_raw", "Gravida",
  "Infant_gender", "Fetal_sex", "Birth_weight_g",
  "Delivery_classification", "Delivery classification", "Delivery_route", "Delivery route",
  "Type_of_Labor", "Type of labor",
  "PROM", "PPROM", "Hypertension", "Hypertension01",
  "Preeclampsia/gestational hypertension", "Hypertension_type",
  "Preeclampsia01", "Gestational_Hypertension01",
  "Gestational_diabetes", "Gestational diabetes", "Chronic_diabetes", "Type_I_diabetes", "Type_II_diabetes", "Prior_GDMA1", "Prior_GDMA2",
  "Chorioamnionitis", "Infection_any", "Oligohydramnios", "IUGR", "Abruption", "Substance_use", "H1N1_vaccine", "HPV_vaccine", "Influenza_vaccine", "Chlamydia", "HIV_1st_OB", "CMV", "HSV_1_2", "HepB", "GBS", "Rubella", "Hx_COVID19",
  legacy_ug_cols, legacy_log_cols,
  "Cumulative_MNPs_ug_per_g", "Cumulative_MNPs_log10ugpg",
  "Cumulative_MNPs_ug", "Cumulative_MNPs_log10ug_per_sample"
)

plac_out <- plac %>%
  select(any_of(front_cols), everything())

# ----------------------------
# Write outputs
# ----------------------------
out_model <- file.path(DIR_TABLES, "mnp_ptb_cohort__analysis_ready_model_table.csv")
out_all_qc <- file.path(DIR_TABLES, "mnp_all_records__qc_flags_and_outlier_flags_wide.csv")
out_thr <- file.path(DIR_TABLES, "mnp_controls__lob_lod_loq_maxblank_sdmethod_by_polymer.csv")
out_blank <- file.path(DIR_TABLES, "mnp_controls__blank_summary_stats_by_polymer.csv")
out_detect <- file.path(DIR_TABLES, "mnp_ptb_cohort__polymer_detection_summary.csv")
out_log <- file.path(DIR_LOGS, "01_build_analysis_ready_inputs__run_log.txt")

write_csv_clean(plac_out, out_model)
write_csv_clean(qc_all, out_all_qc)
write_csv_clean(thresholds, out_thr)
write_csv_clean(blank_summary, out_blank)
write_csv_clean(detection_summary, out_detect)

run_log <- c(
  "01_build_analysis_ready_inputs.R",
  paste0("PROJECT_ROOT: ", PROJECT_ROOT),
  paste0("Input placenta rows: ", nrow(plac0)),
  paste0("Output placenta rows: ", nrow(plac_out)),
  paste0("Output placenta columns: ", ncol(plac_out)),
  paste0("Input controls rows: ", nrow(ctrl0)),
  paste0("QC all-records rows: ", nrow(qc_all)),
  paste0("PE control blank forced to zero: ", FORCE_PE_CONTROL_BLANK_ZERO),
  paste0("PeriBank metadata joined: ", !is.null(peribank_aug)),
  "",
  "Delivery_class counts:",
  capture.output(print(table(plac_out$Delivery_class, useNA = "ifany"))),
  "",
  "Detection summary:",
  capture.output(print(detection_summary, n = Inf)),
  "",
  "Wrote:",
  out_model,
  out_all_qc,
  out_thr,
  out_blank,
  out_detect
)

cat(paste(run_log, collapse = "\n"), file = out_log)
message(paste(run_log, collapse = "\n"))
