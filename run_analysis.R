#!/usr/bin/env Rscript
# Canonical manuscript-analysis runner.
# Run from repository root after restricted inputs are configured.

source("scripts/01_build_analysis_ready_inputs.R")

# Manuscript cohort lock: 171 pregnancies = 85 term + 86 preterm.
.lock_path <- file.path("outputs", "tables", "mnp_ptb_cohort__analysis_ready_model_table.csv")
.lock <- readr::read_csv(.lock_path, show_col_types = FALSE)
.lock <- .lock[as.character(.lock$Delivery_class) %in% c("Term", "Preterm"), , drop = FALSE]
.lock_counts <- table(as.character(.lock$Delivery_class))
if (nrow(.lock) != 171L || unname(.lock_counts["Term"]) != 85L || unname(.lock_counts["Preterm"]) != 86L) {
  stop(
    "Manuscript cohort lock failed. Expected 171 total = 85 term + 86 preterm; found ",
    nrow(.lock), " total = ", unname(.lock_counts["Term"]), " term + ",
    unname(.lock_counts["Preterm"]), " preterm.", call. = FALSE
  )
}
message("Cohort lock passed: 171 total; 85 term; 86 preterm.")
rm(.lock, .lock_counts, .lock_path)
source("scripts/03_qc_figures_abs_mass_polished.R")
source("scripts/04_descriptive_and_univariate_polished.R")
source("scripts/05_clustering_polished.R")

# Internal coefficient-generating stages. Their historical fixed-log10(2)
# effect columns are intermediate only and are never manuscript-facing.
source("scripts/internal/06_ptb_fit_stage.R")
source("scripts/06a_ptb_effect_scale_reference.R")
source("scripts/06b_finalize_primary_ptb_IQR.R")
source("scripts/06c_finalize_ptb_sensitivity_IQR.R")

source("scripts/internal/07_ptb_robustness_fit_stage.R")
source("scripts/07_finalize_ptb_robustness_IQR.R")

source("scripts/internal/08_labor_fit_stage.R")
source("scripts/internal/09_histopathology_fit_stage.R")
source("scripts/10_delivery_status_descriptive_figs_polished.R")
source("scripts/internal/11_ordinal_ga_fit_stage.R")
source("scripts/12_clinical_outcome_correlations.R")
source("scripts/internal/18_validation_fit_stage.R")

# Final manuscript-facing effect scales and figures.
source("scripts/08_finalize_secondary_IQR.R")
source("scripts/08b_finalize_ordinal_labor_figures_IQR.R")
source("scripts/09_finalize_histopathology_IQR.R")

message("Canonical analysis complete. Final continuous-model effects are reported on IQR scales.")
message("Run scripts/20_rebuild_supplementary_table_S1_v8_IQR.R separately if the local Table S1 source workbook is available.")
