# ============================================================
# 20_update_supplementary_table_S1_spike_validation.R
#
# Updates Supplementary Table S1 for the PTB-MNP placenta manuscript.
#
# Starting workbook:
#   04-TableS1_v6.xlsx
#
# Changes:
#   - updates descriptive text in S1c, S1d, and S1f
#   - adds new S1i worksheet with:
#       * calibration-standard composition
#       * paired placenta whole-process recovery
#       * spike-only recovery by control block
#       * process-vs-placenta recovery comparison
#       * validation blank/control summaries
#       * PU/PET QC-category distributions
#       * PU/PET threshold models
#       * Fisher 2x2 analyses
#       * ordinal QC models
#       * Firth models
#       * >LOQ leave-one-out analyses
#       * conservative max-blank subtraction sensitivity
#
# The original workbook is NOT overwritten.
#
# Run from PROJECT_ROOT:
#   source("scripts/20_update_supplementary_table_S1_spike_validation_v2_confirmed_controls.R")
# ============================================================

source("00_config.R")

# -----------------------------
# Packages
# -----------------------------
req_pkgs <- c("openxlsx", "readr", "dplyr", "stringr", "purrr")
missing_pkgs <- req_pkgs[!vapply(req_pkgs, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing_pkgs) > 0) {
  stop(
    "Missing required packages: ",
    paste(missing_pkgs, collapse = ", "),
    "\nInstall with:\ninstall.packages(c(",
    paste(sprintf('"%s"', missing_pkgs), collapse = ", "),
    "))",
    call. = FALSE
  )
}

library(openxlsx)
library(readr)
library(dplyr)
library(stringr)
library(purrr)

# ============================================================
# USER-EDITABLE PATHS
# ============================================================

# If your Table S1 workbook lives elsewhere, edit this path.
TABLE_S1_INPUT <- file.path(PROJECT_ROOT, "04-TableS1_v6.xlsx")

# Output workbook; original remains unchanged.
TABLE_S1_OUTPUT <- file.path(PROJECT_ROOT, "04-TableS1_v7_spike_validation.xlsx")

# Validation output directory.
# Leave NULL to automatically use the newest:
#   outputs/tables/18_placenta_spike_recovery_validation__*
VALIDATION_DIR <- NULL

# ============================================================
# Locate files
# ============================================================

if (!file.exists(TABLE_S1_INPUT)) {
  # A few fallback locations commonly used in this project.
  candidates <- c(
    file.path(PROJECT_ROOT, "04-TableS1_v6.xlsx"),
    file.path(PROJECT_ROOT, "outputs", "04-TableS1_v6.xlsx"),
    file.path(PROJECT_ROOT, "outputs", "tables", "04-TableS1_v6.xlsx")
  )
  candidates <- candidates[file.exists(candidates)]

  if (length(candidates) == 1) {
    TABLE_S1_INPUT <- candidates
  } else {
    stop(
      "Could not locate 04-TableS1_v6.xlsx.\n",
      "Set TABLE_S1_INPUT near the top of this script to the correct path.",
      call. = FALSE
    )
  }
}

if (is.null(VALIDATION_DIR)) {
  validation_candidates <- list.dirs(
    file.path(PROJECT_ROOT, "outputs", "tables"),
    full.names = TRUE,
    recursive = FALSE
  )

  validation_candidates <- validation_candidates[
    grepl(
      "^18_placenta_spike_recovery_validation__",
      basename(validation_candidates)
    )
  ]

  if (length(validation_candidates) == 0) {
    stop(
      "No 18_placenta_spike_recovery_validation__* output folder found under outputs/tables.\n",
      "Run scripts/18_placenta_spike_recovery_validation.R first.",
      call. = FALSE
    )
  }

  VALIDATION_DIR <- validation_candidates[
    order(file.info(validation_candidates)$mtime, decreasing = TRUE)
  ][1]
}

message("Input workbook: ", TABLE_S1_INPUT)
message("Validation directory: ", VALIDATION_DIR)
message("Output workbook: ", TABLE_S1_OUTPUT)

# ============================================================
# Required validation files
# ============================================================

required_files <- c(
  "standard_composition_derived_from_tissue_workbook.csv",
  "placenta_overall_recovery_by_pair.csv",
  "placenta_polymer_recovery_summary.csv",
  "spike_only_process_recovery_summary_by_control_block.csv",
  "spike_only_process_recovery_summary_placenta.csv",
  "process_vs_placenta_recovery_comparison.csv",
  "new_control_totals_and_flags.csv",
  "new_blank_summary_by_polymer.csv",
  "new_controls_long.csv",
  "pu_pet_qc_category_counts_by_delivery_class.csv",
  "pu_pet_qc_2x2_counts.csv",
  "pu_pet_qc_2x2_fisher_stats.csv",
  "pu_pet_qc_binary_ptb_sensitivity_models.csv",
  "pu_pet_qc_ordinal_ptb_sensitivity_models.csv",
  "pu_pet_qc_bias_reduced_models.csv",
  "pu_pet_aboveLOQ_leave_one_out_summary.csv",
  "pu_pet_aboveLOQ_leave_one_out_models.csv",
  "pu_pet_historical_maxblank_subtracted_continuous_models.csv",
  "pu_pet_historical_maxblank_subtraction_specimen_audit.csv"
)

missing_files <- required_files[
  !file.exists(file.path(VALIDATION_DIR, required_files))
]

if (length(missing_files) > 0) {
  stop(
    "The following validation outputs are missing from:\n",
    VALIDATION_DIR, "\n\n  ",
    paste(missing_files, collapse = "\n  "),
    "\n\nRe-run the current scripts/18_placenta_spike_recovery_validation.R.",
    call. = FALSE
  )
}

# ============================================================
# Helpers
# ============================================================

read_validation <- function(filename) {
  readr::read_csv(
    file.path(VALIDATION_DIR, filename),
    show_col_types = FALSE,
    na = c("", "NA")
  )
}

fmt_numeric_df <- function(df) {
  # Keep raw numeric values in Excel; just return as dataframe.
  # Formatting is handled separately in the workbook.
  as.data.frame(df, stringsAsFactors = FALSE)
}

section_fill <- "#E7E6E6"
header_fill  <- "#D9EAF7"
note_fill    <- "#F7F7F7"

style_section <- createStyle(
  fontName = "Aptos Narrow",
  fontSize = 11,
  textDecoration = "bold",
  fgFill = section_fill,
  valign = "top"
)

style_header <- createStyle(
  fontName = "Aptos Narrow",
  fontSize = 10,
  textDecoration = "bold",
  fgFill = header_fill,
  border = "Bottom",
  borderStyle = "thin",
  wrapText = TRUE,
  valign = "top"
)

style_body <- createStyle(
  fontName = "Aptos Narrow",
  fontSize = 10,
  valign = "top"
)

style_note <- createStyle(
  fontName = "Aptos Narrow",
  fontSize = 10,
  textDecoration = "italic",
  fgFill = note_fill,
  wrapText = TRUE,
  valign = "top"
)

style_title <- createStyle(
  fontName = "Aptos Narrow",
  fontSize = 12,
  valign = "top"
)

style_pct <- createStyle(numFmt = "0.0")
style_num2 <- createStyle(numFmt = "0.00")
style_num3 <- createStyle(numFmt = "0.000")
style_p <- createStyle(numFmt = "0.0000")
style_sci <- createStyle(numFmt = "0.00E+00")

write_section <- function(
    wb, sheet, start_row, section_title, df,
    note = NULL,
    numeric_cols_2 = NULL,
    numeric_cols_3 = NULL,
    pct_cols = NULL,
    p_cols = NULL
) {
  r <- start_row

  writeData(
    wb, sheet,
    x = section_title,
    startCol = 1, startRow = r,
    colNames = FALSE
  )
  addStyle(wb, sheet, style_section, rows = r, cols = 1:30, gridExpand = TRUE)
  r <- r + 1

  if (!is.null(note)) {
    writeData(
      wb, sheet,
      x = note,
      startCol = 1, startRow = r,
      colNames = FALSE
    )
    addStyle(wb, sheet, style_note, rows = r, cols = 1:30, gridExpand = TRUE)
    setRowHeights(wb, sheet, rows = r, heights = 30)
    r <- r + 1
  }

  df <- fmt_numeric_df(df)

  writeData(
    wb, sheet,
    x = df,
    startCol = 1, startRow = r,
    colNames = TRUE,
    rowNames = FALSE,
    keepNA = TRUE,
    na.string = ""
  )

  addStyle(
    wb, sheet, style_header,
    rows = r,
    cols = seq_len(ncol(df)),
    gridExpand = TRUE
  )

  if (nrow(df) > 0) {
    addStyle(
      wb, sheet, style_body,
      rows = (r + 1):(r + nrow(df)),
      cols = seq_len(ncol(df)),
      gridExpand = TRUE
    )
  }

  # Optional number formatting by column names.
  apply_num_style <- function(cols, style) {
    if (is.null(cols)) return(invisible(NULL))
    idx <- which(names(df) %in% cols)
    if (length(idx) > 0 && nrow(df) > 0) {
      addStyle(
        wb, sheet, style,
        rows = (r + 1):(r + nrow(df)),
        cols = idx,
        gridExpand = TRUE,
        stack = TRUE
      )
    }
  }

  apply_num_style(numeric_cols_2, style_num2)
  apply_num_style(numeric_cols_3, style_num3)
  apply_num_style(pct_cols, style_pct)
  apply_num_style(p_cols, style_p)

  # Next section: two blank rows.
  r + nrow(df) + 3
}

# ============================================================
# Load original workbook
# ============================================================

wb <- loadWorkbook(TABLE_S1_INPUT)

sheet_names <- names(wb)

if (length(sheet_names) < 8) {
  stop(
    "Expected at least 8 worksheets corresponding to S1a-S1h; found ",
    length(sheet_names), ".",
    call. = FALSE
  )
}

# Existing workbook uses Sheet1 ... Sheet8.
sheet_s1c <- sheet_names[3]
sheet_s1d <- sheet_names[4]
sheet_s1f <- sheet_names[6]

# ============================================================
# Update S1c descriptive text
# ============================================================

s1c_text <- paste0(
  "Supplementary Table S1c | Pyrolysis-gas chromatography/mass spectrometry parameters, ",
  "extracted ion chromatograms, and spike-in standards for MNP polymer quantification. Raw extracted ",
  "ion chromatogram spectra and analytical parameters are provided for the 12 micro- and nanoplastic ",
  "polymers quantified by pyrolysis-gas chromatography/mass spectrometry. The table includes polymer-specific ",
  "extracted ion chromatogram spectra, diagnostic ions or fragments used for polymer identification and ",
  "quantification, retention windows, Py-GC/MS acquisition parameters, and calibration-standard composition. ",
  "The 12-polymer panel included PE, SBR, PP, PVC, N66, PET, N6, PMMA, PC, PU, ABS, and PS. For the new ",
  "whole-process validation experiment, approximately 4 mg of the 12-polymer calibration standard was added ",
  "at the beginning of digestion to the spiked half of fresh paired placental aliquots; matched unspiked aliquots ",
  "underwent otherwise identical processing, and spike-only controls underwent the same workflow without tissue. ",
  "These parameters define the analytical basis for polymer identification and the validation analyses reported in S1i."
)

writeData(
  wb, sheet_s1c,
  x = s1c_text,
  startCol = 1, startRow = 1,
  colNames = FALSE
)

# ============================================================
# Update S1d descriptive text
# ============================================================

s1d_text <- paste0(
  "Supplementary Table S1d | Blank-derived analytical thresholds and polymer detection summary for placental ",
  "MNP measurements. Analytical thresholding and detection summaries for polymer-specific placental micro- and ",
  "nanoplastic measurements are provided. Blank-derived limits were calculated from the historical procedural ",
  "control samples used with the original cohort and include limit of blank, limit of detection, limit of quantification, ",
  "maximum blank value, and SD-based threshold estimates. The detection summary reports the number of placental ",
  "samples evaluated for each polymer, number of non-missing measurements, number and proportion of samples detected ",
  "above zero after blank subtraction, and whether each polymer met the minimum detection criteria for continuous ",
  "downstream modeling. PE was handled according to the prespecified blank-subtraction workflow because it was not ",
  "detected in the historical procedural controls, whereas sparsely detected polymers were flagged for restricted or ",
  "excluded use in continuous models as indicated. The new 2026 spike-only, true, field, and procedural validation ",
  "controls are reported separately in S1i and were not used to redefine these historical cohort LOB/LOD/LOQ thresholds."
)

writeData(
  wb, sheet_s1d,
  x = s1d_text,
  startCol = 1, startRow = 1,
  colNames = FALSE
)

# ============================================================
# Update S1f descriptive text
# ============================================================

s1f_text <- paste0(
  "Supplementary Table S1f | Primary preterm birth models and robustness analyses. Logistic regression models ",
  "were used to evaluate associations between placental micro- and nanoplastic burden and preterm birth, defined ",
  "as delivery at gestational age <37 weeks. Models were fit separately for cumulative MNP burden and each ",
  "polymer-specific burden, with effect estimates reported as odds ratios per doubling of placental polymer burden ",
  "and 95% confidence intervals. The table includes the full set of primary and epidemiologic sensitivity PTB models, ",
  "leave-one-sample-out robustness models, and sample-influence frequency summaries. Primary adjusted models included ",
  "maternal age, BMI, Social Deprivation Index, twin gestation, prior preterm birth, hypertension, and hospital unless ",
  "otherwise specified in the model label. Sensitivity models evaluate robustness to alternative covariate structures, ",
  "exclusion of influential samples, and repeated leave-one-out perturbation. Benjamini-Hochberg q values were calculated ",
  "within the relevant model family. Analytical-threshold, Firth bias-reduced, spike-recovery, and conservative ",
  "blank-subtraction sensitivity analyses for PU and PET are reported separately in Supplementary Table S1i."
)

writeData(
  wb, sheet_s1f,
  x = s1f_text,
  startCol = 1, startRow = 1,
  colNames = FALSE
)

# Match existing workbook title-row style.
for (sh in c(sheet_s1c, sheet_s1d, sheet_s1f)) {
  addStyle(wb, sh, style_title, rows = 1, cols = 1, stack = TRUE)
}

# ============================================================
# Create/recreate S1i
# ============================================================

if ("Sheet9" %in% names(wb)) {
  removeWorksheet(wb, "Sheet9")
}

addWorksheet(wb, "Sheet9")

sheet <- "Sheet9"

s1i_title <- paste0(
  "Supplementary Table S1i | Whole-process spike-in recovery, validation controls, and PU/PET analytical ",
  "sensitivity analyses. Fresh placental tissue from a validation subset was split into matched spiked and unspiked ",
  "aliquots. Approximately 4 mg of a 12-polymer calibration standard was added to the spiked aliquot at the beginning ",
  "of digestion, and paired aliquots underwent identical digestion, centrifugation, secondary processing, transfer, and ",
  "Py-GC/MS analysis. Spike-only controls underwent the same workflow without biological tissue to estimate baseline ",
  "whole-process recovery. Because the standard was added before digestion, differences between spike-only and placental ",
  "recovery are interpreted as matrix-associated attenuation of whole-process recovery and do not isolate instrumental ",
  "matrix effects. New placenta-run true, field, and procedural blanks are reported for validation/background assessment but were not ",
  "used to redefine the historical cohort LOB/LOD/LOQ thresholds in S1d. Control block A corresponds to the placenta run; ",
  "control block B corresponds to the follicular-fluid run. Block B 'blank 3' was confirmed to be a spiked sample with ",
  "unknown spike amount and is retained for audit but excluded from blank and recovery summaries. PU and PET PTB sensitivity analyses include ",
  "analytical-confidence category counts, adjusted >=LOD and >LOQ models, bias-reduced Firth models, >LOQ leave-one-out ",
  "analyses, and a conservative historical max-blank-subtracted continuous-burden sensitivity analysis."
)

writeData(
  wb, sheet,
  x = s1i_title,
  startCol = 1, startRow = 1,
  colNames = FALSE
)
addStyle(wb, sheet, style_title, rows = 1, cols = 1, stack = TRUE)

# Freeze top two rows, consistent with workbook navigation.
freezePane(wb, sheet, firstActiveRow = 3)

# ============================================================
# Load all S1i datasets
# ============================================================

standard_comp <- read_validation(
  "standard_composition_derived_from_tissue_workbook.csv"
)

overall_rec <- read_validation(
  "placenta_overall_recovery_by_pair.csv"
)

placenta_rec <- read_validation(
  "placenta_polymer_recovery_summary.csv"
)

spike_by_block <- read_validation(
  "spike_only_process_recovery_summary_by_control_block.csv"
)

process_vs_tissue <- read_validation(
  "process_vs_placenta_recovery_comparison.csv"
)

control_totals <- read_validation(
  "new_control_totals_and_flags.csv"
)

blank_summary <- read_validation(
  "new_blank_summary_by_polymer.csv"
)

controls_long <- read_validation(
  "new_controls_long.csv"
)

qc_categories <- read_validation(
  "pu_pet_qc_category_counts_by_delivery_class.csv"
)

qc_2x2 <- read_validation(
  "pu_pet_qc_2x2_counts.csv"
)

qc_fisher <- read_validation(
  "pu_pet_qc_2x2_fisher_stats.csv"
)

qc_binary <- read_validation(
  "pu_pet_qc_binary_ptb_sensitivity_models.csv"
)

qc_ordinal <- read_validation(
  "pu_pet_qc_ordinal_ptb_sensitivity_models.csv"
)

qc_firth <- read_validation(
  "pu_pet_qc_bias_reduced_models.csv"
)

loo_summary <- read_validation(
  "pu_pet_aboveLOQ_leave_one_out_summary.csv"
)

loo_models <- read_validation(
  "pu_pet_aboveLOQ_leave_one_out_models.csv"
)

maxblank_models <- read_validation(
  "pu_pet_historical_maxblank_subtracted_continuous_models.csv"
)

maxblank_audit <- read_validation(
  "pu_pet_historical_maxblank_subtraction_specimen_audit.csv"
)

# ============================================================
# Write S1i sections
# ============================================================

r <- 3

r <- write_section(
  wb, sheet, r,
  "S1i-a | 12-polymer calibration-standard composition derived from the paired tissue workbook",
  standard_comp,
  note = paste0(
    "Expected polymer mass per mg of calibration standard was derived from internally consistent ",
    "spike rows in the validation workbook."
  ),
  numeric_cols_2 = c(
    "expected_ug_per_mg_standard",
    "min_ug_per_mg",
    "max_ug_per_mg",
    "cv_pct"
  )
)

r <- write_section(
  wb, sheet, r,
  "S1i-b | Overall apparent whole-process recovery by paired placental specimen",
  overall_rec,
  note = paste0(
    "Overall recovery is the summed paired spiked-minus-unspiked polymer signal divided by ",
    "the summed expected spike mass."
  ),
  numeric_cols_3 = c(
    "spike_mass_mg",
    "expected_total_ug",
    "delta_total_observed_ug"
  ),
  pct_cols = "overall_apparent_recovery_pct"
)

r <- write_section(
  wb, sheet, r,
  "S1i-c | Polymer-specific apparent recovery in paired placental aliquots",
  placenta_rec,
  note = "Summaries include the five complete paired placenta spike/unspike experiments.",
  pct_cols = names(placenta_rec)[grepl("recovery_pct", names(placenta_rec))]
)

r <- write_section(
  wb, sheet, r,
  "S1i-d | Spike-only whole-process recovery by control block",
  spike_by_block,
  note = paste0(
    "Control block A corresponds to the placenta validation run and control block B to the follicular-fluid run. ",
    "Only block A spike-only controls are used as the process comparator for placental recovery."
  ),
  pct_cols = names(spike_by_block)[grepl("recovery_pct", names(spike_by_block))]
)

r <- write_section(
  wb, sheet, r,
  "S1i-e | Placenta-run spike-only recovery versus placental recovery",
  process_vs_tissue,
  note = paste0(
    "The tissue-to-spike-only ratio summarizes additional attenuation observed in placental ",
    "matrix relative to baseline whole-process recovery measured in the matched placenta-run ",
    "spike-only controls (control block A)."
  ),
  pct_cols = names(process_vs_tissue)[
    grepl("recovery_pct|attenuation_pct", names(process_vs_tissue))
  ],
  numeric_cols_3 = "tissue_to_spike_only_ratio"
)

r <- write_section(
  wb, sheet, r,
  "S1i-f | Validation control totals and high-blank flag",
  control_totals,
  note = paste0(
    "Control block B 'blank 3' was confirmed by the analytical laboratory to be a spiked sample ",
    "with unknown spike amount rather than a procedural blank. It is retained in the audit table ",
    "and excluded from blank-background and recovery summaries on the basis of confirmed sample identity."
  ),
  numeric_cols_3 = "total_observed_ug_nd0"
)

r <- write_section(
  wb, sheet, r,
  "S1i-g | Polymer-level summary of new validation blanks",
  blank_summary,
  note = paste0(
    "Descriptive blank summaries use the placenta validation run (control block A) only. ",
    "The follicular-fluid block is not pooled into placental background estimates. These new controls ",
    "do not redefine the historical cohort thresholds in S1d."
  )
)

r <- write_section(
  wb, sheet, r,
  "S1i-h | Polymer-resolved validation control measurements",
  controls_long,
  note = paste0(
    "Reported non-detects are retained as reported_raw and represented separately from ",
    "observed_ug_nd0, where non-detects are treated as zero for descriptive summaries."
  ),
  numeric_cols_3 = c("spike_mass_mg", "observed_ug", "observed_ug_nd0")
)

r <- write_section(
  wb, sheet, r,
  "S1i-i | PU and PET analytical-confidence categories by delivery class",
  qc_categories,
  note = paste0(
    "QC categories are based on the historical polymer-specific thresholds in S1d. ",
    "Rows with missing PTB status are retained in the table but excluded from term/preterm comparisons."
  ),
  pct_cols = "pct_within_group"
)

r <- write_section(
  wb, sheet, r,
  "S1i-j | PU and PET threshold exposure counts by delivery class",
  qc_2x2,
  note = "Counts use the full delivery-class cohort (term n=85; preterm n=86)."
)

r <- write_section(
  wb, sheet, r,
  "S1i-k | Crude Fisher exact statistics for PU and PET analytical thresholds",
  qc_fisher,
  note = paste0(
    "These unadjusted 2x2 analyses are descriptive complements to the adjusted threshold models below."
  ),
  numeric_cols_2 = c("Crude_OR_Fisher", "Crude_CI_lo", "Crude_CI_hi"),
  p_cols = "Fisher_p"
)

r <- write_section(
  wb, sheet, r,
  "S1i-l | Adjusted PU and PET analytical-threshold PTB sensitivity models",
  qc_binary,
  note = paste0(
    "Adjusted complete-case models used n=158 (77 preterm births) and the primary covariate set: ",
    "maternal age, BMI, SDI, twin pregnancy, prior PTB, hypertension, and hospital."
  ),
  numeric_cols_2 = c("OR", "CI_lo", "CI_hi"),
  p_cols = c("p", "q_BH_within_4")
)

r <- write_section(
  wb, sheet, r,
  "S1i-m | Adjusted ordinal analytical-confidence PTB sensitivity models",
  qc_ordinal,
  note = "Odds ratios are reported per one-category increase in analytical confidence.",
  numeric_cols_2 = c("OR_per_category", "CI_lo", "CI_hi"),
  p_cols = c("p", "q_BH_within_2")
)

r <- write_section(
  wb, sheet, r,
  "S1i-n | Bias-reduced Firth logistic regression for PU and PET threshold models",
  qc_firth,
  note = paste0(
    "Firth models address potential small-cell bias/quasi-separation in the threshold-based sensitivity analyses."
  ),
  numeric_cols_2 = c("OR", "CI_lo", "CI_hi"),
  p_cols = c("p", "q_BH_within_4")
)

r <- write_section(
  wb, sheet, r,
  "S1i-o | Leave-one-out robustness summary for adjusted >LOQ models",
  loo_summary,
  note = paste0(
    "Each adjusted >LOQ model was refit after sequentially excluding one complete-case participant."
  ),
  numeric_cols_2 = names(loo_summary)[grepl("OR|CI", names(loo_summary))],
  p_cols = names(loo_summary)[grepl("^p|max_p|min_p", names(loo_summary))]
)

r <- write_section(
  wb, sheet, r,
  "S1i-p | Leave-one-out model estimates for adjusted >LOQ PU and PET analyses",
  loo_models,
  numeric_cols_2 = names(loo_models)[grepl("^OR$|CI_lo|CI_hi", names(loo_models))],
  p_cols = names(loo_models)[names(loo_models) %in% c("p", "q_BH")]
)

r <- write_section(
  wb, sheet, r,
  "S1i-q | Conservative historical max-blank-subtracted continuous PTB sensitivity models",
  maxblank_models,
  note = paste0(
    "The historical polymer-specific maximum blank mass was subtracted once per analytical run, ",
    "values were floored at zero, and tissue burden was recalculated before refitting the primary ",
    "adjusted continuous model. This is a conservative sensitivity analysis, not a replacement exposure metric."
  ),
  numeric_cols_3 = "Historical_MaxBlank_ug_per_run",
  numeric_cols_2 = c("OR_per_doubling_plus1", "CI_lo", "CI_hi"),
  p_cols = c("p", "q_BH_within_PU_PET")
)

r <- write_section(
  wb, sheet, r,
  "S1i-r | Specimen-level audit for the historical max-blank subtraction sensitivity analysis",
  maxblank_audit
)

# ============================================================
# S1i formatting
# ============================================================

# Useful widths. Long data blocks can still be filtered/sorted in Excel.
setColWidths(
  wb, sheet,
  cols = 1:26,
  widths = c(
    24, 24, 22, 22, 22, 22, 22, 18, 18, 18, 18, 18, 18,
    18, 18, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20
  )
)

# Keep the workbook usable rather than forcing all cells to wrap.
setRowHeights(wb, sheet, rows = 1, heights = 18)

# ============================================================
# Save workbook
# ============================================================

saveWorkbook(
  wb,
  file = TABLE_S1_OUTPUT,
  overwrite = TRUE
)

message("")
message("Supplementary Table S1 update complete.")
message("Created: ", TABLE_S1_OUTPUT)
message("")
message("Confirmed control provenance applied in S1i:")
message("  - control block A = placenta")
message("  - control block B = follicular fluid")
message("  - block-B 'blank 3' = spiked sample with unknown spike amount; excluded from blank/recovery summaries by identity")
message("  - new validation blanks do not redefine historical LOB/LOD/LOQ thresholds")
