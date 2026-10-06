# Polymer-resolved placental micro- and nanoplastics and preterm birth

Reproducibility code for the manuscript **“Polymer-resolved Placental Micro- and Nanoplastics Associate with Preterm Birth, Labor Dysfunction, and Placental Microcalcification.”**

This repository contains the analysis and figure-generation code used for the manuscript. Participant-level clinical data, specimen identifiers, pathology source files, and raw analytical validation workbooks are **not distributed in this repository**.

## Repository scope

The workflow includes:

- construction of analysis-ready placental MNP and control datasets;
- polymer-specific analytical QC and batch-matched procedural-control thresholds;
- descriptive and univariate analyses;
- compositional clustering and visualization;
- adjusted preterm-birth models and robustness analyses;
- term labor-augmentation analyses;
- placental histopathology analyses;
- gestational-age category and clinical correlation analyses;
- paired whole-process spike-recovery validation;
- PU/PET analytical-threshold, Firth, leave-one-out, and max-blank sensitivity analyses;
- generation of manuscript Figure 5, Supplementary Figures 4 and 7, and Supplementary Table S1 updates.

## Directory structure

```text
.
├── 00_config.R
├── R/
│   ├── figure_theme_mnp.R
│   └── utils_mnp.R
├── scripts/
├── data/
│   ├── README.md
│   ├── external/
│   ├── restricted/
│   └── supplementary/
└── outputs/                 # generated locally; ignored by git
```

## Data access and privacy

Human-subject and specimen-level inputs are intentionally excluded from the public repository. The scripts retain variable names such as `PregnancyID` and `SampleID` because those fields are required for joins and leave-one-out analyses, but **no identifier values are included here**.

See `data/README.md` for the expected local input structure. Data availability and access should follow the final manuscript Data Availability statement and applicable institutional/IRB restrictions.

## Software

Analyses were developed in **R 4.5.2**. Install required CRAN packages with:

```r
source("install_packages.R")
```

Exact package versions were not captured in the source archive; users seeking archival reproducibility should create a release-specific `renv.lock` from the analysis environment before publication.

## Configuration

Run analyses from the repository root. By default, `00_config.R` treats the current working directory as the project root.

For restricted local inputs, either use the default locations under `data/restricted/` or set environment variables before running the scripts. Common overrides include:

```r
Sys.setenv(
  MNP_PROJECT_ROOT = "/path/to/Polymer-resolved-Placental-MNP-PTB",
  MNP_PERIBANK_XLSX = "/secure/path/peribank_metadata.xlsx",
  MNP_PATHOLOGY_INPUT = "/secure/path/pathology_scores.csv",
  MNP_SPIKE_VALIDATION_DIR = "/secure/path/spike_validation",
  MNP_TABLE_S1_INPUT = "/path/to/04-TableS1_v6.xlsx"
)
```

For the spike-validation source workbooks, exact files may also be set with:

```r
Sys.setenv(
  MNP_TISSUE_SPIKE_XLSX = "/secure/path/paired_placenta_spike.xlsx",
  MNP_VALIDATION_CONTROLS_XLSX = "/secure/path/validation_controls.xlsx"
)
```

## Recommended run order

Run scripts from the repository root. Numbering reflects the manuscript-development workflow; deprecated/intermediate scripts have been removed.

```r
source("scripts/01_build_analysis_ready_inputs.R")
source("scripts/03_qc_figures_abs_mass_polished.R")
source("scripts/04_descriptive_and_univariate_polished.R")
source("scripts/05_clustering_polished.R")
source("scripts/06_primary_ptb_models.R")
source("scripts/06b_primary_ptb_figure3_clean_polished.R")
source("scripts/06c_primary_ptb_supp_sensitivity_figures_polished.R")
source("scripts/07_ptb_robustness_figure4_polished.R")
source("scripts/08_labor_outcome_supplement_polished.R")
source("scripts/08b_labor_figure_polish_polished.R")
source("scripts/09_histopathology_models_polished.R")
source("scripts/10_delivery_status_descriptive_figs_polished.R")
source("scripts/11_ordinal_ga_category_models.R")
source("scripts/12_clinical_outcome_correlations.R")
source("scripts/18_placenta_spike_recovery_validation.R")
source("scripts/19_figure5_analytical_validation_panels.R")
source("scripts/20_update_supplementary_table_S1_spike_validation.R")
source("scripts/21_update_supplementary_figures_PTB_robustness.R")
```

Scripts 18–21 require the validation inputs and/or Supplementary Table S1 source workbook described in `data/README.md`.

## Important analytical notes

- The primary PTB endpoint is **preterm birth <37 weeks overall**, not spontaneous PTB.
- Primary adjusted PTB models include maternal age, BMI, Social Deprivation Index, twin gestation, prior PTB, hypertension, and hospital.
- LOB/LOD/LOQ and maximum-blank sensitivity analyses use the **batch-matched procedural controls from the primary cohort**.
- The subsequent fresh-tissue spike-recovery experiment is used to characterize whole-process recovery and matrix-associated attenuation; it does **not** numerically recovery-correct the primary cohort.
- For the validation controls, control block A corresponds to the placenta run and control block B to the follicular-fluid run. The sample labeled block-B `blank 3` was confirmed to be a spike of unknown amount and is excluded from blank/recovery summaries by sample identity rather than statistical outlier exclusion.

## Generated outputs

Scripts write tables, figures, logs, and diagnostics under `outputs/`. This directory is ignored by git because some generated intermediate files may contain participant-level identifiers.

## Citation

Please cite the associated manuscript when using this code. 
