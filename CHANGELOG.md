# Changelog

## 2026-10-07 — manuscript effect-scale and cohort-lock update

- Updated manuscript-facing continuous-model effect reporting from the historical fixed `log10(2)` display increment to IQR-based reporting for `log10(1 + burden [µg/g])` exposures.
- Retained the original fitted model coefficients, standard errors, P values, and multiple-testing results; this is a reporting-scale correction rather than a model refit.
- Locked the manuscript analytical cohort at 171 pregnancies (85 term, 86 preterm).
- Updated PTB primary/sensitivity and leave-one-out outputs to a fixed polymer-specific IQR defined in the primary PTB analysis cohort.
- Updated labor and ordinal gestational-age models to analysis-specific IQR scales.
- Updated histopathology models to outcome × delivery-class × polymer-specific IQR scales, held fixed across sensitivity models within each stratum.
- Updated PU/PET continuous analytical-sensitivity reporting to IQR-scaled effects; binary LOD/LOQ models are unchanged.
- Revised Supplementary Figure 5 workflow to include labor augmentation and ordinal gestational-age analyses.
- Replaced the Supplementary Table S1 update workflow with the v8 IQR rebuild for the locked 171-subject cohort.
- Moved coefficient-generating scripts that retain historical intermediate effect columns to `scripts/internal/` and clearly marked those columns as non-manuscript-facing.
