# Polymer-resolved placental micro- and nanoplastics and preterm birth

**Publication-facing R analysis code** for the placental MNP–preterm birth study. The repository includes scripts supporting **all analyses reported in the main manuscript and supplement**. It excludes unreported hypothetical recovery-correction simulations and superseded exploratory analyses. The manuscript and final supplementary tables remain the authorities for the reported numerical results.

> **Data availability and reproducibility:** This is a **code-only repository**. Restricted clinical/specimen-level data, validation-control workbooks, and some frozen intermediate analysis outputs are not distributed. The complete pipeline cannot be independently executed from this repository alone. This has not been validated by a fresh end-to-end run. See [data/README.md](data/README.md) and [REPRODUCIBILITY_NOTES.md](REPRODUCIBILITY_NOTES.md).

## Final figure and table crosswalk

| Manuscript exhibit | Analysis | Main publication script(s) |
|---|---|---|
| Table 1 | Cohort clinical characteristics | `01`; final Supplementary Table S1a–S1b |
| Figure 1 and Table 2 | Polymer burden, composition, detection | `03`, `04`, `10` |
| Figure 2 | Adjusted PTB association forest, transformed-IQR ORs | **`23`**, upstream `06`, `22` |
| Figure 3 | Microcalcification severity | **`29`**, upstream `09` |
| Supplementary Figure S1 | Analytical QC and detection | `03`, `20` |
| Supplementary Figure S2 | Descriptive/univariate burden and gestational age | `04`, `10` |
| Supplementary Figure S3 | CLR-PCA / composition clustering | `05` |
| Supplementary Figure S4 | PTB model-specification and influence robustness | **`24`, `27`**, upstream `06c`, `07` |
| Supplementary Figure S5 | ABS, SBR, PET, and PU nonlinearity | **`37`** |
| Supplementary Figure S6 | 12-polymer recovery and PU/PET analytical sensitivity | **`39`**, frozen `18` outputs |
| Supplementary Figure S7 | Labor augmentation and ordinal gestational age | **`28`**, upstream `08`, `11`, `26` |
| Supplementary Figure S8 | Perivillous fibrin severity | **`29`**, upstream `09` |
| Supplementary Table S1a–S1i | Cohort, assay, clinical models, pathology, spike validation | Final **`04-TableS1_v11_IQR_ANALYTICAL_ALIGNED.xlsx`**, separately held; see [ANALYSIS_MANIFEST.md](ANALYSIS_MANIFEST.md) |

The figure numbers above reflect the **manuscript display order**, not the older figure numbers embedded in script filenames or historical filenames. Analyses moved from main figures to supplementary figures are still represented in this public code.

## Exposure transformation and reporting scale

For the primary PTB logistic models, polymer-specific placental burden \(C\) (µg/g) is modeled as \(X=\log_{10}(1+C)\). The reported continuous exposure contrast is **one polymer-specific IQR of the transformed burden**, with \(\mathrm{OR}_{IQR}=\exp(\beta\,IQR(X))\). The source coefficients, robust HC1 uncertainty estimates, and model specifications define the inference; rescaling does not refit the underlying linear regression.

The continuous PU/PET sensitivity estimates in Supplementary Figure S6 and Supplementary Table S1i-q use a **fixed original-model transformed IQR** for both the original and maximum-blank-subtracted models. The threshold models (at/above LOD and above LOQ, including Firth) are **binary contrasts, not IQR contrasts**. Script `37` retains its declared model-specific IQR references for spline/influence diagnostics. Since some polymers show nonlinearity, a single linear IQR OR should not be interpreted as a constant exposure response throughout the observed range.

**Earlier upstream scripts may emit historical per-doubling-labelled intermediate results. They are retained for provenance/dependencies, not for manuscript reporting.** Use `23`, `27`, `28`, `29`, `37`, and `39` for the corresponding final displays. No historical per-doubling effect estimates should be inserted into the manuscript.

## Getting started

Requires R with the dependencies listed in `install_packages.R`. Start **at the repository root**:

```r
source('install_packages.R')
source('00_config.R')
```

The new `run_analysis.R` is a **targeted launcher**, not a one-command pipeline. It will not run anything unless a specific target is requested:

```sh
Rscript run_analysis.R --list
Rscript run_analysis.R --check spline
Rscript run_analysis.R --run spline
Rscript run_analysis.R --check validation-s6
Rscript run_analysis.R --run validation-s6
```

For interactive R, you can also source the chosen script directly from the repository root:

```r
source('scripts/37_four_polymer_nonlinearity_robustness.R')
source('scripts/39_merged_placenta_analytical_suppfigure.R')
```

**Prerequisite outputs:** the primary IQR sequence is `06` → `22` → `23`; PTB robustness uses `24` → `27` and previous `07`/leave-one-out outputs; ordinal/labor finalization uses `26` → `28`; pathology uses `09` → `29`. These are **dependency guides, not an automatically tested topological run order**. In particular, Script `39` needs seven frozen Script 18 validation CSV inputs in `outputs/tables/39_merged_analytical_inputs/`, including a restricted specimen-level audit. Script `37` requires `outputs/tables/mnp_ptb_cohort__analysis_ready_model_table.csv`. Use the script headers and `REPRODUCIBILITY_NOTES.md` to locate other necessary upstream outputs.

`--check` verifies only a small set of explicitly known file prerequisites; it is **not a complete validation** of packages, script-specific assumptions, or numerical replication. Each `--run` executes just one selected script and can overwrite its usual output paths: keep frozen manuscript output files backed up.

## Repository structure

- `00_config.R` — project-relative paths and shared constants; environment variables can override restricted input paths.
- `run_analysis.R` — explicit target listing, basic prerequisite checks, and single-target execution.
- `scripts/` — study scripts, including upstream/legacy intermediates required to generate publication results.
- `R/` — shared helpers and plotting theme.
- `data/README.md` — restricted input requirements; **do not commit** patient-level data.
- `ANALYSIS_MANIFEST.md` — authoritative versus historical/upstream scripts.
- `METHODS_SPECIFICATIONS.md` — analysis-specification summary.
- `RESULTS_CROSSWALK.csv` — exhibit-to-script provenance.
- `SCRIPT39_IQR_CONVERSION_REFERENCE.csv` — reporting-scale crosscheck.
- `REPRODUCIBILITY_NOTES.md` — limitations and release checklist.

## Release and data protection

The **v11 supplementary workbook** is held separately and is **not automatically public data**: S1a and validation appendices contain individual-level clinical or specimen data. Publication or repository distribution requires the study team's data-release approval. Do not commit local `outputs/`, private model-ready data, individual specimen audit files, or raw clinical/validation workbooks.

This README describes the curated code state as supplied for repository staging. It does not claim the contents have been committed to GitHub or that a remote GitHub Action has reproduced the figures. Complete the manuscript numerical crosscheck before assigning a final publication release tag.
