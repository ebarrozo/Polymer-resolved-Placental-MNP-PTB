# Polymer-resolved placental micro- and nanoplastics and preterm birth

Reproducibility code for the manuscript **“Polymer-resolved Placental Micro- and Nanoplastics Associate with Preterm Birth, Labor Augmentation, and Placental Microcalcification.”**

This repository contains analysis and figure-generation code used for the manuscript. Participant-level clinical data, specimen identifiers, pathology source files, and raw analytical-validation workbooks are **not distributed in this repository**.

## Repository scope

The workflow includes:

- construction of the locked analysis-ready placental MNP cohort;
- polymer-specific analytical QC and procedural-control summaries;
- descriptive and univariate analyses;
- compositional clustering and visualization;
- adjusted preterm-birth models and robustness analyses;
- exploratory spontaneous-labor augmentation analyses;
- ordered gestational-age sensitivity analyses;
- placental histopathology analyses;
- clinical-correlation analyses;
- paired whole-process spike-recovery validation;
- PU/PET analytical-threshold, Firth, leave-one-out, and max-blank sensitivity analyses;
- generation of the revised manuscript-facing IQR-scaled figures and Supplementary Table S1 support tables.

## Final effect-scale convention

Continuous regression models use `log10(1 + placental polymer burden [µg/g])` as the exposure. Manuscript-facing continuous-model effects are reported as effect estimates per **interquartile-range (IQR) increase in the transformed exposure**.

This replaces an earlier development-stage reporting convention that multiplied fitted coefficients by `log10(2)`. That earlier convention corresponds to a twofold increase in `(1 + x)`, not an exact doubling of raw burden `x`. The final IQR reporting scale changes only the presentation unit: fitted coefficients, standard errors, P values, and multiple-testing conclusions are unchanged.

For PTB sensitivity and leave-one-out analyses, each polymer uses the same IQR defined in the primary PTB analysis cohort. Labor, ordinal gestational-age, and histopathology analyses use IQRs appropriate to their own analysis populations or strata, as documented in the corresponding scripts.

## Cohort lock

The manuscript analytical cohort is locked at **171 pregnancies: 85 term and 86 preterm**. The public workflow expects the restricted analysis-ready input to reproduce this cohort. Downstream manuscript-facing scripts should not be run against an alternate cohort without deliberate reanalysis.

## Directory structure

```text
.
├── 00_config.R
├── run_analysis.R
├── R/
│   ├── figure_theme_mnp.R
│   └── utils_mnp.R
├── scripts/
│   ├── internal/             # coefficient-generating fit stages; not manuscript-facing effect scales
│   └── ...                   # finalization/figure scripts
├── data/
│   ├── README.md
│   ├── external/
│   ├── restricted/
│   └── supplementary/
└── outputs/                  # generated locally; ignored by git
```

## Data access and privacy

Human-subject and specimen-level inputs are intentionally excluded from the public repository. Scripts retain variable names such as `PregnancyID` and `SampleID` because those fields are required for joins and leave-one-out analyses, but **no identifier values are included here**.

See `data/README.md` for the expected local input structure. Data availability and access should follow the final manuscript Data Availability statement and applicable institutional/IRB restrictions.

## Software

Analyses were developed in **R 4.5.2**. Install required CRAN packages with:

```r
source("install_packages.R")
```

Exact package versions were not captured in the original source archive. For archival reproduction, create a release-specific `renv.lock` from the analysis environment before publication.

## Configuration

Run analyses from the repository root. By default, `00_config.R` treats the current working directory as the project root.

For restricted local inputs, use the default locations under `data/restricted/` or set environment variables before running. Common overrides include:

```r
Sys.setenv(
  MNP_PROJECT_ROOT = "/path/to/Polymer-resolved-Placental-MNP-PTB",
  MNP_PERIBANK_XLSX = "/secure/path/peribank_metadata.xlsx",
  MNP_PATHOLOGY_INPUT = "/secure/path/pathology_scores.csv",
  MNP_SPIKE_VALIDATION_DIR = "/secure/path/spike_validation",
  MNP_TABLE_S1_INPUT = "/path/to/04-TableS1_v7_spike_validation.xlsx"
)
```

For spike-validation workbooks, exact files may also be set with:

```r
Sys.setenv(
  MNP_TISSUE_SPIKE_XLSX = "/secure/path/paired_placenta_spike.xlsx",
  MNP_VALIDATION_CONTROLS_XLSX = "/secure/path/validation_controls.xlsx"
)
```

## Recommended run order

The canonical analysis can be run with:

```r
source("run_analysis.R")
```

The runner performs the following stages:

1. build and QC the analysis-ready dataset;
2. descriptive and clustering analyses;
3. generate original fitted PTB coefficients and analysis samples in `scripts/internal/`;
4. re-express the primary PTB and robustness estimates on the final fixed primary-cohort IQR scale;
5. generate labor, ordinal gestational-age, pathology, and analytical-validation fit stages;
6. re-express those secondary continuous-model effects on their final analysis-specific IQR scales;
7. generate the manuscript-facing revised figures and supporting tables.

The scripts under `scripts/internal/` are retained because they reproduce the fitted models and exact analysis samples. Some of those internal outputs preserve historical fixed-`log10(2)` columns for auditability. **Those columns are not manuscript-facing and must not be interpreted as raw-burden doubling effects.**

### Supplementary Table S1

If the local v7 Supplementary Table source workbook is available, rebuild the revised analytical supplement after `run_analysis.R`:

```r
source("scripts/20_rebuild_supplementary_table_S1_v8_IQR.R")
```

The v8 rebuild refreshes the locked 171-subject analytical sheets and replaces continuous-model effect reporting with the final IQR scale while preserving analytical-method and validation-support content.

## Important analytical notes

- The primary PTB endpoint is **preterm birth <37 weeks overall**, not spontaneous PTB.
- Primary adjusted PTB models include maternal age, BMI, Social Deprivation Index, twin gestation, prior PTB, hypertension, and hospital.
- PTB sensitivity/leave-one-out estimates use a fixed polymer-specific IQR from the primary PTB analysis cohort.
- Labor and ordinal gestational-age analyses use analysis-specific IQRs.
- Histopathology models use the primary-model IQR within each outcome × delivery-class × polymer stratum and hold that IQR fixed across sensitivity models in the same stratum.
- Binary analytical-threshold models (`≥LOD`, `>LOQ`) are unaffected by the continuous-effect-scale correction.
- Historical blank-derived thresholds and maximum-blank sensitivity analyses use the procedural controls associated with the primary cohort.
- The subsequent fresh-tissue spike-recovery experiment characterizes whole-process recovery and matrix-associated attenuation; it does **not** numerically recovery-correct the primary cohort.
- Validation control block A corresponds to the placenta run and block B to the follicular-fluid run. The sample labeled block-B `blank 3` was confirmed to be a spike of unknown amount and is excluded from blank/recovery summaries by sample identity rather than statistical outlier exclusion.

## Generated outputs

Scripts write tables, figures, logs, and diagnostics under `outputs/`. This directory is ignored by git because some generated intermediate files may contain participant-level identifiers.

## Revision note

The October 2026 repository revision updates manuscript-facing continuous-model effects from the historical fixed-`log10(2)` display unit to IQR-based reporting. This is an effect-scale correction, not a model refit. The revision also aligns Supplementary Figure 5 with the labor-augmentation and ordered gestational-age sensitivity analyses and updates the Supplementary Table S1 workflow to the locked 171-subject cohort.

## Citation

Please cite the associated manuscript when using this code. Add the final journal citation and DOI here after publication.

## License

No software license is assigned in this bundle. Select a license consistent with institutional policy before making the repository public.
