# Reproducibility and release limitations

1. This package is a code and manuscript-output **release candidate**, not a GitHub upload and not a complete rerun of all analysis scripts.
2. **Verified by study author:** updated Script 39 produced the uploaded final Figure S6. The IQR S1i values were reconciled in the separate authoritative v11 workbook.
3. **Inputs restricted:** anonymization of a specimen ID does not automatically make the underlying clinical data publicly redistributable. Source cohort, pathology, procedural controls, raw spike recovery, and individual blank-subtraction audits are not present.
4. **Intermediate outputs needed:** Scripts written sequentially throughout development use selected frozen output directories. Use `RESULTS_CROSSWALK.csv` to identify provenance, and retain frozen inputs locally for independent re-execution.
5. **Historical labeling:** some upstream scripts emit earlier per-doubling displays. Only designated final IQR scripts and the S1 v11 artifact should be used for manuscript numbers.
6. **Before pushing:** inspect `git status`, remove any locally tracked patient-level/clinical output files, ensure `.gitignore` excludes those files, compare SHA256 sums, and perform a dry run of critical figure-generation scripts in the secured study environment.
7. **Release tag suggestion:** `ptb-submission-candidate-v1` until final manuscript audit; replace with a final tag only after manuscript figure/table/legend consistency review.
