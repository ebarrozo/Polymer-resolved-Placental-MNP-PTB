#!/usr/bin/env Rscript
# =============================================================================
# 20_rebuild_supplementary_table_S1_v8_IQR.R
#
# Rebuilds Supplementary Table S1 after the reviewer Point 1 effect-scale
# correction and the current 171-subject cohort lock.
#
# IMPORTANT
# ---------
# This is intentionally a rebuild of the manuscript-facing analytical sheets,
# not a cosmetic search/replace of "per doubling".
#
# The previous v7 workbook contains a mixture of current validation content and
# older 173-row analytical tables. This script:
#
#   S1a  -> current analysis-ready 171-subject dataset
#   S1b  -> current expanded cohort characteristics (85 term / 86 preterm)
#   S1c  -> PRESERVED from the current v7 workbook
#   S1d  -> historical blank-derived thresholds PRESERVED + current 171-subject
#           polymer detection summary
#   S1e  -> current clinical-outcome Spearman correlation table
#   S1f  -> current PTB continuous models on the corrected IQR effect scale +
#           leave-one-out robustness on the same fixed primary-cohort IQR scale
#   S1g  -> exploratory labor-augmentation AND ordinal gestational-age models
#           on analysis-specific IQR scales (mirrors revised Supplementary Fig. 5)
#   S1h  -> current histopathology models on stratum-specific IQR scales
#   S1i  -> PRESERVED validation workbook content, except the continuous
#           max-blank sensitivity subsection is updated to the IQR effect scale
#
# Output:
#   04-TableS1_v8_IQR_final.xlsx
#   outputs/tables/31_SuppTableS1_v8_audit.csv
#
# Run from PROJECT_ROOT:
#   source("scripts/20_rebuild_supplementary_table_S1_v8_IQR.R")
# =============================================================================

source("00_config.R")

required_pkgs <- c("openxlsx", "readr", "dplyr", "stringr", "purrr", "tidyr")
missing_pkgs <- required_pkgs[
  !vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_pkgs)) {
  stop(
    "Missing required R packages: ", paste(missing_pkgs, collapse=", "),
    "\nInstall with: install.packages(c(",
    paste(sprintf('"%s"', missing_pkgs), collapse=", "), "))",
    call.=FALSE
  )
}

suppressPackageStartupMessages({
  library(openxlsx)
  library(readr)
  library(dplyr)
  library(stringr)
  library(purrr)
  library(tidyr)
})

# -----------------------------------------------------------------------------
# Paths and helpers
# -----------------------------------------------------------------------------

latest_dir <- function(parent, prefix, required_files = NULL) {
  dirs <- list.dirs(parent, recursive = FALSE, full.names = TRUE)
  dirs <- dirs[grepl(paste0("^", prefix), basename(dirs))]
  if (!length(dirs)) return(NA_character_)

  if (!is.null(required_files)) {
    keep <- vapply(
      dirs,
      function(d) all(file.exists(file.path(d, required_files))),
      logical(1)
    )
    dirs <- dirs[keep]
    if (!length(dirs)) return(NA_character_)
  }

  dirs[order(file.info(dirs)$mtime, decreasing = TRUE)][1]
}

first_existing <- function(paths) {
  hit <- paths[file.exists(paths)]
  if (!length(hit)) return(NA_character_)
  hit[[1]]
}

clean_names <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  x
}

TABLE_INPUT_ENV <- Sys.getenv("MNP_TABLE_S1_INPUT", unset = "")
TABLE_INPUT <- first_existing(c(
  if (nzchar(TABLE_INPUT_ENV)) TABLE_INPUT_ENV else NA_character_,
  file.path(PROJECT_ROOT, "data", "supplementary", "04-TableS1_v7_spike_validation.xlsx"),
  file.path(PROJECT_ROOT, "04-TableS1_v7_spike_validation.xlsx"),
  file.path(DIR_TABLES, "04-TableS1_v7_spike_validation.xlsx")
))
if (is.na(TABLE_INPUT)) {
  stop(
    "Could not locate 04-TableS1_v7_spike_validation.xlsx. ",
    "Place it in PROJECT_ROOT or set TABLE_INPUT manually.",
    call.=FALSE
  )
}

TABLE_OUTPUT <- file.path(PROJECT_ROOT, "04-TableS1_v8_IQR_final.xlsx")
AUDIT_OUTPUT <- file.path(DIR_TABLES, "31_SuppTableS1_v8_audit.csv")

MODEL_TABLE <- file.path(DIR_TABLES, "mnp_ptb_cohort__analysis_ready_model_table.csv")
if (!file.exists(MODEL_TABLE)) {
  stop("Missing current analysis-ready model table: ", MODEL_TABLE, call.=FALSE)
}

message("Base workbook: ", TABLE_INPUT)
message("Current model table: ", MODEL_TABLE)
message("Output workbook: ", TABLE_OUTPUT)

# -----------------------------------------------------------------------------
# Current cohort lock
# -----------------------------------------------------------------------------

cohort <- read_csv(MODEL_TABLE, show_col_types=FALSE, guess_max=10000) %>%
  filter(as.character(Delivery_class) %in% c("Term", "Preterm"))

n_all <- nrow(cohort)
n_term <- sum(as.character(cohort$Delivery_class) == "Term")
n_preterm <- sum(as.character(cohort$Delivery_class) == "Preterm")

if (!identical(c(n_all, n_term, n_preterm), c(171L, 85L, 86L))) {
  stop(
    "Cohort-lock check failed. Expected 171 total = 85 term + 86 preterm; found ",
    n_all, " total = ", n_term, " term + ", n_preterm, " preterm.\n",
    "Do not rebuild the supplement until the analysis-ready table matches the manuscript cohort.",
    call.=FALSE
  )
}

message("Cohort lock passed: 171 total; 85 term; 86 preterm.")

# -----------------------------------------------------------------------------
# Locate current analysis outputs
# -----------------------------------------------------------------------------

PTB_DIR <- file.path(DIR_TABLES, "24_ptb_sensitivity_IQR")
if (!file.exists(file.path(PTB_DIR, "PTB_ALL_models_long_IQR.csv"))) {
  stop("Missing Point 1 PTB IQR outputs in ", PTB_DIR, call.=FALSE)
}

SFIG4_DIR <- file.path(DIR_TABLES, "27_SuppFigure4_IQR_final")
if (!file.exists(file.path(SFIG4_DIR, "robustness_leave_one_out_models_IQR.csv"))) {
  stop("Missing finalized SFig4 IQR robustness outputs in ", SFIG4_DIR, call.=FALSE)
}

SECONDARY_DIR <- latest_dir(
  DIR_TABLES,
  "26_secondary_models_IQR__",
  c("labor_model_all_results_IQR.csv", "ordinal_ga_category_all_models_IQR.csv")
)
if (is.na(SECONDARY_DIR)) {
  stop(
    "No complete 26_secondary_models_IQR__* directory found containing both labor and ordinal IQR tables.",
    call.=FALSE
  )
}

HIST_DIR <- file.path(DIR_TABLES, "29_histopathology_IQR_final")
if (!file.exists(file.path(HIST_DIR, "histopath_model_all_results_IQR.csv"))) {
  stop("Missing finalized histopathology IQR outputs in ", HIST_DIR, call.=FALSE)
}

CORR_DIR <- latest_dir(
  DIR_TABLES,
  "12_clinical_correlations__",
  "Supplementary_Table_1c_clinical_correlations.csv"
)
if (is.na(CORR_DIR)) {
  stop(
    "No current clinical-correlation output found. Run scripts/12_clinical_outcome_correlations.R first.",
    call.=FALSE
  )
}

VALIDATION_DIR <- latest_dir(
  DIR_TABLES,
  "18_placenta_spike_recovery_validation__",
  c(
    "pu_pet_historical_maxblank_subtracted_continuous_models.csv",
    "pu_pet_historical_maxblank_subtraction_specimen_audit.csv"
  )
)
if (is.na(VALIDATION_DIR)) {
  stop("No complete script-18 validation output directory found.", call.=FALSE)
}

# Corrected continuous validation effects generated during Point 1.
CONT_IQR_CANDIDATES <- c(
  file.path(
    latest_dir(
      DIR_TABLES,
      "26_secondary_models_IQR__",
      "PU_PET_continuous_analytical_sensitivity_IQR.csv"
    ),
    "PU_PET_continuous_analytical_sensitivity_IQR.csv"
  )
)
CONT_IQR_PATH <- CONT_IQR_CANDIDATES[file.exists(CONT_IQR_CANDIDATES)][1]
if (is.na(CONT_IQR_PATH) || length(CONT_IQR_PATH)==0) {
  stop(
    "Could not locate PU_PET_continuous_analytical_sensitivity_IQR.csv. ",
    "Run script 26 first.",
    call.=FALSE
  )
}

# Current cohort descriptive outputs.
TABLE1_PATH <- file.path(DIR_TABLES, "table1__cohort_characteristics_term_vs_ptb.csv")
if (!file.exists(TABLE1_PATH)) {
  message("Current cohort-characteristics output is missing; running script 04 once.")
  source(file.path("scripts", "04_descriptive_and_univariate.R"))
}
if (!file.exists(TABLE1_PATH)) stop("Could not generate current cohort-characteristics table.", call.=FALSE)

DETECT_PATH <- file.path(DIR_TABLES, "mnp_ptb_cohort__polymer_detection_summary.csv")
if (!file.exists(DETECT_PATH)) {
  # This is generated upstream; use the current model table to generate a minimal current summary.
  warning(
    "Current polymer detection summary not found. A current 171-row detection summary will be ",
    "derived directly from the analysis-ready table for S1d."
  )
}

# -----------------------------------------------------------------------------
# Read current result tables
# -----------------------------------------------------------------------------

tbl1 <- read_csv(TABLE1_PATH, show_col_types=FALSE)
corr <- read_csv(
  file.path(CORR_DIR, "Supplementary_Table_1c_clinical_correlations.csv"),
  show_col_types=FALSE
)
ptb <- read_csv(file.path(PTB_DIR, "PTB_ALL_models_long_IQR.csv"), show_col_types=FALSE)
loo <- read_csv(file.path(SFIG4_DIR, "robustness_leave_one_out_models_IQR.csv"), show_col_types=FALSE)
loo_summary <- read_csv(file.path(SFIG4_DIR, "robustness_leave_one_out_summary_IQR.csv"), show_col_types=FALSE)
labor <- read_csv(file.path(SECONDARY_DIR, "labor_model_all_results_IQR.csv"), show_col_types=FALSE)
ordinal <- read_csv(file.path(SECONDARY_DIR, "ordinal_ga_category_all_models_IQR.csv"), show_col_types=FALSE)
hist <- read_csv(file.path(HIST_DIR, "histopath_model_all_results_IQR.csv"), show_col_types=FALSE)
cont_iqr <- read_csv(CONT_IQR_PATH, show_col_types=FALSE)

# -----------------------------------------------------------------------------
# Detect current 171-subject polymer summary for S1d
# -----------------------------------------------------------------------------

polymer_cols <- grep("_ug_per_g$", names(cohort), value=TRUE)
polymer_cols <- polymer_cols[!grepl("^Cumulative_", polymer_cols)]

detection_current <- map_dfr(polymer_cols, function(cn) {
  x <- suppressWarnings(as.numeric(cohort[[cn]]))
  tibble(
    Polymer = sub("_ug_per_g$", "", cn),
    column = cn,
    n_total = nrow(cohort),
    n_nonmissing = sum(is.finite(x)),
    n_detected_gt0 = sum(is.finite(x) & x > 0),
    detection_prop = n_detected_gt0 / n_total,
    continuous_model_ok = n_detected_gt0 >= 6
  )
})

# -----------------------------------------------------------------------------
# Load base workbook
# -----------------------------------------------------------------------------

wb <- loadWorkbook(TABLE_INPUT)
sheets <- names(wb)

if (length(sheets) < 9) {
  stop("Expected at least 9 worksheets (S1a-S1i); found ", length(sheets), call.=FALSE)
}

sheet_s1a <- sheets[1]
sheet_s1b <- sheets[2]
sheet_s1c <- sheets[3]  # preserved
sheet_s1d <- sheets[4]
sheet_s1e <- sheets[5]
sheet_s1f <- sheets[6]
sheet_s1g <- sheets[7]
sheet_s1h <- sheets[8]
sheet_s1i <- sheets[9]

# -----------------------------------------------------------------------------
# Workbook styles
# -----------------------------------------------------------------------------

style_title <- createStyle(
  fontSize=11, textDecoration="bold", wrapText=TRUE,
  valign="top", fgFill="#D9EAF7", border="Bottom"
)
style_section <- createStyle(
  fontSize=10, textDecoration="bold", fgFill="#EDEDED",
  border="TopBottom", wrapText=TRUE
)
style_header <- createStyle(
  fontSize=9, textDecoration="bold", fgFill="#F2F2F2",
  border="Bottom", halign="center", valign="center", wrapText=TRUE
)
style_body <- createStyle(fontSize=8.5, valign="top")
style_note <- createStyle(fontSize=8.5, textDecoration="italic", wrapText=TRUE, valign="top")

clear_sheet <- function(wb, sheet, max_rows=5000, max_cols=220) {
  deleteData(
    wb, sheet=sheet,
    cols=seq_len(max_cols), rows=seq_len(max_rows),
    gridExpand=TRUE
  )
}

write_title <- function(wb, sheet, text, width=40) {
  writeData(wb, sheet, text, startRow=1, startCol=1, colNames=FALSE)
  addStyle(wb, sheet, style_title, rows=1, cols=1, gridExpand=TRUE, stack=TRUE)
  setColWidths(wb, sheet, cols=1, widths=width)
  setRowHeights(wb, sheet, rows=1, heights=55)
}

write_section <- function(wb, sheet, title, dat, start_row, start_col=1) {
  writeData(wb, sheet, title, startRow=start_row, startCol=start_col, colNames=FALSE)
  addStyle(wb, sheet, style_section, rows=start_row, cols=start_col, gridExpand=TRUE, stack=TRUE)

  if (is.null(dat) || nrow(dat)==0) {
    return(start_row + 2)
  }

  writeData(
    wb, sheet, dat,
    startRow=start_row+1, startCol=start_col,
    colNames=TRUE, rowNames=FALSE, keepNA=TRUE
  )
  addStyle(
    wb, sheet, style_header,
    rows=start_row+1,
    cols=start_col:(start_col+ncol(dat)-1),
    gridExpand=TRUE, stack=TRUE
  )
  if (nrow(dat) > 0) {
    addStyle(
      wb, sheet, style_body,
      rows=(start_row+2):(start_row+1+nrow(dat)),
      cols=start_col:(start_col+ncol(dat)-1),
      gridExpand=TRUE, stack=TRUE
    )
  }
  setColWidths(wb, sheet, cols=start_col:(start_col+ncol(dat)-1), widths="auto")
  start_row + nrow(dat) + 4
}

# =============================================================================
# S1a — current 171-subject analysis-ready table
# =============================================================================

clear_sheet(wb, sheet_s1a)
s1a_desc <- paste0(
  "Supplementary Table S1a | Current analysis-ready clinical metadata and placental MNP dataset used for downstream modeling. ",
  "This table contains the locked manuscript cohort of 171 pregnancies (85 term and 86 preterm) and includes clinical covariates, ",
  "polymer-specific placental burdens, transformed burden variables, cumulative burden, and analytical QC flags used in the current analyses."
)
write_title(wb, sheet_s1a, s1a_desc, width=55)
writeData(wb, sheet_s1a, cohort, startRow=3, startCol=1, colNames=TRUE, rowNames=FALSE, keepNA=TRUE)
addStyle(wb, sheet_s1a, style_header, rows=3, cols=1:ncol(cohort), gridExpand=TRUE, stack=TRUE)
addStyle(wb, sheet_s1a, style_body, rows=4:(3+nrow(cohort)), cols=1:ncol(cohort), gridExpand=TRUE, stack=TRUE)
freezePane(wb, sheet_s1a, firstActiveRow=4, firstActiveCol=2)
setColWidths(wb, sheet_s1a, cols=1:ncol(cohort), widths="auto")

# =============================================================================
# S1b — current cohort characteristics
# =============================================================================

clear_sheet(wb, sheet_s1b)
s1b_desc <- paste0(
  "Supplementary Table S1b | Expanded maternal, pregnancy, delivery, and labor characteristics stratified by delivery timing. ",
  "Continuous variables are summarized as median (interquartile range) and mean ± standard deviation; categorical variables are summarized as counts and percentages. ",
  "P values compare term and preterm groups using Wilcoxon rank-sum tests for continuous variables and Fisher exact or chi-square tests for categorical variables as appropriate. ",
  "The current cohort includes 85 term and 86 preterm deliveries."
)
write_title(wb, sheet_s1b, s1b_desc, width=55)
writeData(wb, sheet_s1b, tbl1, startRow=3, startCol=1, colNames=TRUE, rowNames=FALSE, keepNA=TRUE)
addStyle(wb, sheet_s1b, style_header, rows=3, cols=1:ncol(tbl1), gridExpand=TRUE, stack=TRUE)
addStyle(wb, sheet_s1b, style_body, rows=4:(3+nrow(tbl1)), cols=1:ncol(tbl1), gridExpand=TRUE, stack=TRUE)
freezePane(wb, sheet_s1b, firstActiveRow=4)
setColWidths(wb, sheet_s1b, cols=1:ncol(tbl1), widths="auto")

# =============================================================================
# S1c — preserved exactly from v7
# =============================================================================
# No changes: analytical Py-GC/MS parameters / chromatogram-standard support.

# =============================================================================
# S1d — preserve historical threshold block, update cohort detection summary
# =============================================================================

# Read threshold values from base workbook before clearing.
s1d_old <- read.xlsx(
  TABLE_INPUT, sheet=sheet_s1d,
  colNames=FALSE, skipEmptyRows=FALSE, skipEmptyCols=FALSE
)

# Historical threshold header is on row 3/4 in the current workbook.
# Identify rows using known column names where possible.
header_row <- which(apply(s1d_old, 1, function(z) any(as.character(z)=="Analyte", na.rm=TRUE)))[1]
if (is.na(header_row)) header_row <- 3

threshold_headers <- clean_names(unlist(s1d_old[header_row, 1:13]))
threshold_dat <- s1d_old[(header_row+1):nrow(s1d_old), 1:13, drop=FALSE]
names(threshold_dat) <- make.unique(threshold_headers)

# dplyr's .data pronoun requires string column names; use the first
# column name explicitly rather than numeric indexing.
first_threshold_col <- names(threshold_dat)[1]
threshold_dat <- threshold_dat %>%
  filter(
    !is.na(.data[[first_threshold_col]]) &
      as.character(.data[[first_threshold_col]]) != ""
  )

clear_sheet(wb, sheet_s1d)
s1d_desc <- paste0(
  "Supplementary Table S1d | Historical blank-derived analytical thresholds and current polymer detection summary for placental MNP measurements. ",
  "Historical procedural controls used with the original clinical cohort define the blank-derived threshold values shown at left. ",
  "The detection summary at right has been regenerated from the locked 171-pregnancy analysis-ready dataset. ",
  "The 2026 spike/recovery validation controls are reported separately in S1i and were not used to redefine the historical cohort thresholds."
)
write_title(wb, sheet_s1d, s1d_desc, width=55)

writeData(wb, sheet_s1d, "Historical blank-derived analytical thresholds", startRow=3, startCol=1, colNames=FALSE)
addStyle(wb, sheet_s1d, style_section, rows=3, cols=1, stack=TRUE)
writeData(wb, sheet_s1d, threshold_dat, startRow=4, startCol=1, colNames=TRUE, rowNames=FALSE, keepNA=TRUE)
addStyle(wb, sheet_s1d, style_header, rows=4, cols=1:ncol(threshold_dat), gridExpand=TRUE, stack=TRUE)

start_det_col <- ncol(threshold_dat) + 3
writeData(wb, sheet_s1d, "Current 171-subject polymer detection summary", startRow=3, startCol=start_det_col, colNames=FALSE)
addStyle(wb, sheet_s1d, style_section, rows=3, cols=start_det_col, stack=TRUE)
writeData(wb, sheet_s1d, detection_current, startRow=4, startCol=start_det_col, colNames=TRUE, rowNames=FALSE, keepNA=TRUE)
addStyle(
  wb, sheet_s1d, style_header, rows=4,
  cols=start_det_col:(start_det_col+ncol(detection_current)-1),
  gridExpand=TRUE, stack=TRUE
)
setColWidths(wb, sheet_s1d, cols=1:(start_det_col+ncol(detection_current)-1), widths="auto")

# =============================================================================
# S1e — current clinical correlations
# =============================================================================

clear_sheet(wb, sheet_s1e)
s1e_desc <- paste0(
  "Supplementary Table S1e | Spearman correlations between placental MNP burden and clinical variables. ",
  "Spearman correlation analysis evaluated associations between log10-transformed cumulative or polymer-specific placental burden and clinical variables. ",
  "The table reports nominal P values and Benjamini–Hochberg-adjusted q values across the full polymer-by-outcome matrix, together with outcome- and polymer-specific q values. ",
  "Positive rho values indicate direct associations and negative values inverse associations."
)
write_title(wb, sheet_s1e, s1e_desc, width=55)
writeData(wb, sheet_s1e, corr, startRow=3, startCol=1, colNames=TRUE, rowNames=FALSE, keepNA=TRUE)
addStyle(wb, sheet_s1e, style_header, rows=3, cols=1:ncol(corr), gridExpand=TRUE, stack=TRUE)
addStyle(wb, sheet_s1e, style_body, rows=4:(3+nrow(corr)), cols=1:ncol(corr), gridExpand=TRUE, stack=TRUE)
freezePane(wb, sheet_s1e, firstActiveRow=4)
setColWidths(wb, sheet_s1e, cols=1:ncol(corr), widths="auto")

# =============================================================================
# S1f — PTB models on corrected IQR scale + robustness
# =============================================================================

clear_sheet(wb, sheet_s1f)
s1f_desc <- paste0(
  "Supplementary Table S1f | Primary preterm-birth models and robustness analyses. ",
  "Continuous exposure models use log10(1 + placental burden [µg/g]) and are reported as odds ratios per polymer-specific interquartile-range (IQR) increase in the transformed burden. ",
  "For PTB sensitivity and leave-one-out comparisons, each polymer uses the same IQR defined in the primary PTB analysis cohort so that variation across models reflects coefficient stability rather than changes in the reporting unit. ",
  "P values and multiple-testing conclusions are unchanged by this effect-scale re-expression."
)
write_title(wb, sheet_s1f, s1f_desc, width=60)

ptb_out <- ptb %>%
  select(
    any_of(c(
      "ModelTag","ModelLabel","Polymer","PolymerLabel","ExposureCol",
      "N","Events","NonEvents","Covariates","IQR_log10_1px",
      "OR_IQR","OR_IQR_lo","OR_IQR_hi","logOR_IQR","SE_logOR_IQR",
      "p","q_BH","fit_status"
    ))
  )

loo_out <- loo %>%
  select(
    any_of(c(
      "Polymer","PolymerLabel","leave_out_sample","analysis",
      "N","Events","NonEvents","IQR_log10_1px",
      "OR_IQR","OR_IQR_lo","OR_IQR_hi",
      "p","q_BH_within_polymer","fit_status"
    ))
  )

loo_sum_out <- loo_summary %>%
  select(any_of(c(
    "Polymer","PolymerLabel","n_iterations","OR_min","OR_median","OR_max",
    "CIlo_min","CIhi_max"
  )))

r <- 3
r <- write_section(wb, sheet_s1f, "S1f-a | PTB logistic-regression models on the corrected IQR scale", ptb_out, r)
r <- write_section(wb, sheet_s1f, "S1f-b | Leave-one-out PTB robustness on the fixed primary-cohort IQR scale", loo_out, r)
r <- write_section(wb, sheet_s1f, "S1f-c | Leave-one-out PTB summary", loo_sum_out, r)
freezePane(wb, sheet_s1f, firstActiveRow=4)

# =============================================================================
# S1g — revised SF5: labor augmentation + ordinal gestational-age models
# =============================================================================

clear_sheet(wb, sheet_s1g)
s1g_desc <- paste0(
  "Supplementary Table S1g | Exploratory analyses of spontaneous labor augmentation and delivery timing. ",
  "Labor-augmentation models compare augmented versus non-augmented spontaneous labor and report odds ratios per analysis-specific IQR increase in log10(1 + placental polymer burden [µg/g]). ",
  "Ordinal gestational-age models evaluate ordered delivery-timing categories; odds ratios >1 indicate greater odds of delivery in an earlier gestational-age category. ",
  "Because these models use different analysis populations, each analysis uses its own reconstructed IQR of the transformed exposure."
)
write_title(wb, sheet_s1g, s1g_desc, width=60)

labor_out <- labor %>%
  select(
    any_of(c(
      "ModelTag","ModelLabel","Polymer","PolymerLabel","ExposureCol",
      "N","Events_Augmented","NonEvents_SpontaneousNoAug",
      "Covariates","IQR_transformed","OR_IQR","OR_IQR_lo","OR_IQR_hi",
      "p","q_BH_all_analytes","q_BH_polymers_only","q_BH_within_family",
      "q_BH_labor_prespec","nominal_p05","sig_q05_all_analytes",
      "sig_q05_labor_prespec","fit_status"
    ))
  )

ordinal_out <- ordinal %>%
  select(
    any_of(c(
      "ModelTag","ModelLabel","Polymer","PolymerLabel","ExposureCol",
      "N","Covariates","IQR_transformed",
      "OR_earlier_IQR","OR_earlier_IQR_lo","OR_earlier_IQR_hi",
      "OR_later_IQR","OR_later_IQR_lo","OR_later_IQR_hi",
      "p","q_BH_all_analytes","q_BH_polymers_only","q_BH_within_family",
      "fit_engine","fit_status"
    ))
  )

r <- 3
r <- write_section(wb, sheet_s1g, "S1g-a | Spontaneous labor requiring augmentation", labor_out, r)
r <- write_section(wb, sheet_s1g, "S1g-b | Ordinal gestational-age category models", ordinal_out, r)
freezePane(wb, sheet_s1g, firstActiveRow=4)

# =============================================================================
# S1h — histopathology on corrected stratum-specific IQR scale
# =============================================================================

clear_sheet(wb, sheet_s1h)
s1h_desc <- paste0(
  "Supplementary Table S1h | Polymer-specific associations with placental histopathology severity. ",
  "Microcalcification severity among microcalcification-positive placentae was modeled using negative-binomial regression and is reported as incidence rate ratios (IRRs). ",
  "Perivillous fibrin deposition severity among positive placentae was modeled using gamma regression with a log link and is reported as relative ratios (RRs). ",
  "Models are stratified by delivery class and reported per outcome-, delivery-class-, and polymer-specific IQR increase in log10(1 + placental polymer burden [µg/g]). ",
  "Within each stratum, the primary-model IQR is held fixed across sensitivity models."
)
write_title(wb, sheet_s1h, s1h_desc, width=60)

hist_out <- hist %>%
  select(
    any_of(c(
      "Outcome","OutcomeCol","ModelType","Delivery_class",
      "ModelTag","ModelLabel","Polymer","PolymerLabel","PolymerFamily",
      "ExposureCol","N","Covariates","IQR_primary",
      "Ratio_IQR","Ratio_IQR_lo","Ratio_IQR_hi",
      "p","q_global","q_family","nominal_p05",
      "sig_global_q05","sig_family_q05","fit_status"
    ))
  )
writeData(wb, sheet_s1h, hist_out, startRow=3, startCol=1, colNames=TRUE, rowNames=FALSE, keepNA=TRUE)
addStyle(wb, sheet_s1h, style_header, rows=3, cols=1:ncol(hist_out), gridExpand=TRUE, stack=TRUE)
addStyle(wb, sheet_s1h, style_body, rows=4:(3+nrow(hist_out)), cols=1:ncol(hist_out), gridExpand=TRUE, stack=TRUE)
freezePane(wb, sheet_s1h, firstActiveRow=4)
setColWidths(wb, sheet_s1h, cols=1:ncol(hist_out), widths="auto")

# =============================================================================
# S1i — preserve validation content, update continuous max-blank subsection
# =============================================================================

s1i_matrix <- read.xlsx(
  TABLE_INPUT, sheet=sheet_s1i,
  colNames=FALSE, skipEmptyRows=FALSE, skipEmptyCols=FALSE
)

# Locate section labels by first column text.
first_col <- as.character(s1i_matrix[[1]])
row_q <- which(str_detect(first_col, fixed("S1i-q", ignore_case=TRUE)))[1]
row_r <- which(str_detect(first_col, fixed("S1i-r", ignore_case=TRUE)))[1]

if (is.na(row_q) || is.na(row_r) || row_r <= row_q) {
  stop(
    "Could not locate S1i-q / S1i-r section boundaries in the v7 workbook. ",
    "S1i was not modified.",
    call.=FALSE
  )
}

# Build corrected S1i-q table.
cont_out <- cont_iqr %>%
  select(
    any_of(c(
      "Polymer","Sensitivity","Historical_MaxBlank_ug_per_run",
      "N","Events","IQR_log2_1px",
      "OR_IQR","OR_IQR_lo","OR_IQR_hi",
      "p","se_method","n_zero_after_subtraction","Covariates"
    ))
  )

# Clear the current q-section only, leaving S1i-r and everything after it intact.
deleteData(
  wb, sheet=sheet_s1i,
  cols=1:40,
  rows=row_q:(row_r-1),
  gridExpand=TRUE
)

q_title <- paste0(
  "S1i-q | Conservative historical max-blank subtraction sensitivity: continuous PU/PET PTB models. ",
  "Continuous effects are reported as odds ratios per analysis-specific IQR increase in transformed burden; ",
  "the max-blank-subtracted exposure is a conservative sensitivity analysis and is not treated as a corrected primary exposure."
)
writeData(wb, sheet_s1i, q_title, startRow=row_q, startCol=1, colNames=FALSE)
addStyle(wb, sheet_s1i, style_section, rows=row_q, cols=1, stack=TRUE)
writeData(
  wb, sheet_s1i, cont_out,
  startRow=row_q+1, startCol=1,
  colNames=TRUE, rowNames=FALSE, keepNA=TRUE
)
addStyle(
  wb, sheet_s1i, style_header,
  rows=row_q+1, cols=1:ncol(cont_out),
  gridExpand=TRUE, stack=TRUE
)

# Update top S1i description if present.
s1i_top <- paste0(
  "Supplementary Table S1i | Whole-process spike-in recovery, validation controls, and PU/PET analytical sensitivity analyses. ",
  "This worksheet reports calibration-standard composition, paired placental recovery, spike-only process recovery, validation controls and blanks, ",
  "analytical-threshold sensitivity models, bias-reduced and leave-one-out analyses, and conservative historical max-blank subtraction. ",
  "Binary threshold-model odds ratios retain their original interpretation; continuous max-blank sensitivity models are reported per analysis-specific IQR."
)
writeData(wb, sheet_s1i, s1i_top, startRow=1, startCol=1, colNames=FALSE)
addStyle(wb, sheet_s1i, style_title, rows=1, cols=1, stack=TRUE)

# =============================================================================
# Workbook-wide formatting touches
# =============================================================================

for (sh in sheets) {
  freezePane(wb, sh, firstActiveRow=2)
}

# Prevent pathological auto-widths in text-heavy analysis sheets.
for (sh in c(sheet_s1a, sheet_s1b, sheet_s1d, sheet_s1e, sheet_s1f, sheet_s1g, sheet_s1h, sheet_s1i)) {
  setColWidths(wb, sh, cols=1:80, widths="auto")
}

# -----------------------------------------------------------------------------
# QA audit before save
# -----------------------------------------------------------------------------

audit <- tibble(
  check = c(
    "Locked cohort total",
    "Locked cohort term",
    "Locked cohort preterm",
    "S1a rows",
    "S1b source rows",
    "S1d detection summary n_total unique",
    "S1e correlation rows",
    "S1f PTB model rows",
    "S1f LOO model rows",
    "S1g labor rows",
    "S1g ordinal rows",
    "S1h histopathology rows",
    "S1i continuous maxblank rows"
  ),
  value = c(
    n_all,
    n_term,
    n_preterm,
    nrow(cohort),
    nrow(tbl1),
    paste(sort(unique(detection_current$n_total)), collapse=","),
    nrow(corr),
    nrow(ptb_out),
    nrow(loo_out),
    nrow(labor_out),
    nrow(ordinal_out),
    nrow(hist_out),
    nrow(cont_out)
  ),
  expected = c(
    "171","85","86","171",
    "current script-04 output",
    "171",
    "current script-12 output",
    "current Point-1 output",
    "current Point-1 output",
    "current script-26 output",
    "current script-26 output",
    "current script-29 output",
    "4 (PU/PET × original/maxblank), unless validation workflow changed"
  )
)

write_csv(audit, AUDIT_OUTPUT)

# Save, then reload and scan manuscript-facing sheets for known stale language.
saveWorkbook(wb, TABLE_OUTPUT, overwrite=TRUE)

check_wb <- loadWorkbook(TABLE_OUTPUT)
stale_hits <- list()

for (sh in c(sheet_s1f, sheet_s1g, sheet_s1h, sheet_s1i)) {
  x <- read.xlsx(
    TABLE_OUTPUT, sheet=sh,
    colNames=FALSE, skipEmptyRows=FALSE, skipEmptyCols=FALSE
  )
  txt <- as.character(unlist(x))
  bad <- txt[
    str_detect(
      txt,
      regex(
        "per doubling of placental polymer burden|OR_2x|Ratio_2x|OR_per_doubling_plus1",
        ignore_case=TRUE
      )
    )
  ]
  if (length(bad)) {
    stale_hits[[sh]] <- unique(bad)
  }
}

if (length(stale_hits)) {
  stale_df <- bind_rows(lapply(names(stale_hits), function(sh) {
    tibble(sheet=sh, residual_text=stale_hits[[sh]])
  }))
  stale_path <- file.path(DIR_TABLES, "31_SuppTableS1_v8_residual_scale_terms.csv")
  write_csv(stale_df, stale_path)
  warning(
    "Residual legacy effect-scale terms remain in the saved workbook. ",
    "Inspect: ", stale_path
  )
} else {
  message("Effect-scale audit passed: no legacy per-doubling/OR_2x terms in S1f-S1i.")
}

message("")
message("============================================================")
message("Supplementary Table S1 v8 rebuild complete")
message("============================================================")
message("Output workbook: ", TABLE_OUTPUT)
message("Audit: ", AUDIT_OUTPUT)
message("Cohort: 171 total; 85 term; 86 preterm")
message("S1c preserved from v7.")
message("S1i preserved except corrected continuous max-blank subsection and top description.")
