# Data inputs

The public repository contains **code only**. Human-subject, specimen-level, and internal analytical source files are intentionally excluded.

## Restricted model-ready inputs

By default, script 01 expects the following under `data/restricted/model_ready/`:

```text
placenta_model_ready_wide.csv
placenta_model_ready_long.csv
controls_blank_ready_wide.csv
controls_blank_ready_long.csv
```

The exact contents of these files should be distributed only as permitted by the study's IRB, consent, and data-use requirements.

## Clinical metadata

Optional clinical enrichment used by script 01 is configured with `MNP_PERIBANK_XLSX` and defaults to:

```text
data/restricted/peribank_metadata.xlsx
```

This file is not included publicly.

## Pathology data

Script 09 uses a restricted pathology-score CSV configured by `MNP_PATHOLOGY_INPUT`, defaulting to:

```text
data/restricted/pathology_scores.csv
```

## Spike-recovery validation

Script 18 uses the paired fresh-placenta spike/unspike workbook and the validation-control workbook. Place them under:

```text
data/restricted/spike_validation/
```

or set:

- `MNP_SPIKE_VALIDATION_DIR`
- `MNP_TISSUE_SPIKE_XLSX`
- `MNP_VALIDATION_CONTROLS_XLSX`

These source workbooks are not included because they contain internal specimen-level analytical records.

## Supplementary Table S1 source workbook

Script 20 updates the manuscript Supplementary Table S1 workbook. The default expected source is:

```text
data/supplementary/04-TableS1_v6.xlsx
```

Alternatively set `MNP_TABLE_S1_INPUT` to a local copy.

## External deprivation-index files

Script 01 optionally uses public SVI/SDI source files. Defaults are:

```text
data/external/SVI_2022_US_ZCTA.csv
data/external/rgcsdi-2015-2019-zcta.csv
```

These external source files are not redistributed here. Obtain them from their authoritative sources and set `MNP_SVI_PATH` / `MNP_SDI_PATH` if stored elsewhere.
