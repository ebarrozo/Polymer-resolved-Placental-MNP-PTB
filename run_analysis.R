#!/usr/bin/env Rscript
# PTB publication-code dispatcher. Run from repository root.
# Does not run any analysis without an explicit --run TARGET.
# Restricted inputs and previously frozen outputs are NOT shipped in this repository.

args <- commandArgs(trailingOnly = TRUE)
help <- function() {
  cat(paste(
    'Usage:',
    '  Rscript run_analysis.R --list',
    '  Rscript run_analysis.R --check TARGET',
    '  Rscript run_analysis.R --run TARGET',
    '',
    'Run from the repository root. --check checks known input paths,',
    'not packages, cohort content, or every historical dependency.',
    'No automatic --all mode: different scripts consume frozen upstream outputs.',
    sep='\n'), '\n')
}
if (identical(args, '--list')) args <- c('--list','all')
if (length(args) != 2L || !args[1] %in% c('--list','--check','--run') ||
    (args[1] == '--list' && args[2] != 'all')) {
  help()
  quit(status=if (length(args)==0L) 0L else 2L)
}
root <- normalizePath(getwd(), mustWork=TRUE)
if (!file.exists(file.path(root,'00_config.R')) || !dir.exists(file.path(root,'scripts'))) {
  stop('Run from the repository root (00_config.R and scripts/ must be present).')
}
# Only these named final outputs should be used in the manuscript. Other available
# source scripts are retained for upstream derivation, not publication plotting.
targets <- c(
  'prepare'       = '01_build_analysis_ready_inputs.R',
  'qc'            = '03_qc_figures_abs_mass_polished.R',
  'descriptive'   = '04_descriptive_and_univariate_polished.R',
  'clustering'    = '05_clustering_polished.R',
  'delivery'      = '10_delivery_status_descriptive_figs_polished.R',
  'primary-base'  = '06_primary_ptb_models.R',
  'primary-scale' = '22_ptb_effect_scale_audit.R',
  'primary-iqr'   = '23_finalize_primary_ptb_IQR.R',
  'robust-scale'  = '24_rescale_ptb_sensitivity_IQR.R',
  'robust-iqr'    = '27_finalize_suppfigure4_IQR.R',
  'labor-base'    = '08_labor_outcome_supplement_polished.R',
  'ordinal-base'  = '11_ordinal_ga_category_models.R',
  'labor-scale'   = '26_secondary_models_IQR_audit.R',
  'labor-iqr'     = '28_finalize_ordinal_labor_IQR_figures.R',
  'path-base'     = '09_histopathology_models_polished.R',
  'path-iqr'      = '29_finalize_histopathology_IQR.R',
  'correlations'  = '12_clinical_outcome_correlations.R',
  'heterogeneity' = '32_ptb_heterogeneity_sensitivity_LOCKED_PTB.R',
  'spline'        = '37_four_polymer_nonlinearity_robustness.R',
  'validation-base'= '18_placenta_spike_recovery_validation.R',
  'validation-s6'= '39_merged_placenta_analytical_suppfigure.R'
)
notes <- c(
  'prepare'='Restricted cleaned specimen/control input tables',
  'primary-base'='Legacy intermediate model output; NOT final IQR figure',
  'primary-scale'='Historical scale-audit intermediary (not final plot)',
  'primary-iqr'='Requires completed Script 22 audit folder',
  'robust-scale'='Requires Script 06 and Script 22 output folders',
  'robust-iqr'='Requires Script 24 and frozen leave-one-out/robustness outputs',
  'labor-scale'='Requires completed labor/ordinal model results',
  'labor-iqr'='Requires Script 26 output folder',
  'path-iqr'='Requires Script 09 pathology model outputs',
  'spline'='Requires outputs/tables/mnp_ptb_cohort__analysis_ready_model_table.csv',
  'validation-base'='Requires restricted paired tissue and control workbooks',
  'validation-s6'='Requires seven frozen Script 18 CSVs (including restricted specimen audit) in outputs/tables/39_merged_analytical_inputs'
)
if (args[1] == '--list') {
  for (nm in names(targets)) cat(sprintf('%-17s %s%s\n',nm,targets[[nm]],
        if(nm %in% names(notes)) paste0('  [',notes[[nm]],']') else ''))
  quit(status=0L)
}
target <- args[2]
if (!target %in% names(targets)) {help();stop('Unknown target: ',target)}
script <- file.path('scripts', targets[[target]])
if (!file.exists(script)) stop('Missing expected file: ',script)
# Only basic, deterministic checks, not a false claim of full dependency validation.
known <- switch(target,
  'prepare'=file.path('data','restricted','model_ready',c('placenta_model_ready_wide.csv','controls_blank_ready_wide.csv','controls_blank_ready_long.csv')),
  'spline'=file.path('outputs','tables','mnp_ptb_cohort__analysis_ready_model_table.csv'),
  'heterogeneity'=file.path('outputs','tables','mnp_ptb_cohort__analysis_ready_model_table.csv'),
  'validation-s6'=file.path('outputs','tables','39_merged_analytical_inputs',c(
    'spike_only_process_recovery.csv', 'placenta_paired_polymer_recovery.csv',
    'pu_pet_qc_category_counts_by_delivery_class.csv',
    'pu_pet_qc_binary_ptb_sensitivity_models.csv','pu_pet_qc_bias_reduced_models.csv',
    'pu_pet_historical_maxblank_subtracted_continuous_models.csv',
    'pu_pet_historical_maxblank_subtraction_specimen_audit.csv')),
  character(0))
if (identical(args[1],'--check')) {
  cat('Target:',target,'\nScript:',script,'\n')
  if(length(known)) {
    missing <- known[!file.exists(known)]
    if(length(missing)) {cat('Missing known prerequisites:\n',paste0('  ',missing,collapse='\n'),'\n',sep='');quit(status=1L)}
    cat('Known paths present; OTHER dependencies may still be required.\n')
  } else cat('Script found. Full input dependencies must be reviewed in script header.\n')
  quit(status=0L)
}
if (length(known) && any(!file.exists(known))) {
  stop('Missing known prerequisite(s) for ',target,':\n',paste(known[!file.exists(known)],collapse='\n'),
       '\nUse --check and consult README.md. No analysis run.')
}
cat('Running target:',target,'\nScript:',script,'\n')
source(script, local=.GlobalEnv, chdir=FALSE)
cat('Completed target:',target,'\n')
