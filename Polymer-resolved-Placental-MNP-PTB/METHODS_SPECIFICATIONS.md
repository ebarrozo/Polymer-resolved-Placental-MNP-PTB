# PTB Methods — locked analysis reference (draft for manuscript alignment)

- Primary exposure: polymer-specific placental burden in µg/g, transformed as `log10(1 + burden)`.
- Reported continuous odds ratios: per polymer-specific IQR increment **on the transformed scale**. Base choice is immaterial when increment and coefficient are transformed consistently.
- Primary PTB model: adjusted logistic regression with maternal age, BMI, Social Deprivation Index, twins, prior PTB, hypertension, and hospital; HC1 robust confidence intervals.
- Primary effect table/figure: Script 23 / Supplementary Table S1f; do not use older per-doubling display outputs.
- Nonlinearity: Script 37 four-polymer robust HC1 tests with BH correction across ABS/SBR/PET/PU; nonlinear associations should not be interpreted as constant ORs across the entire exposure range.
- Figure S6 continuous analytical sensitivity: Script 39 uses the *fixed original-burden transformed IQR* within each polymer for both original and max-blank-adjusted models; p-values and fitted regressions are inherited unchanged from Script 18. Threshold contrasts are not IQR effects.
- Recovery experiments do not identify outcome-dependent recovery or correct for it empirically. Do not describe hypothetical simulations as measured recovery correction.
- No additional range/scaling reference table is planned.

## Verify before textual lock
- Confirm Figure 2 printed OR/CI versus final S1f and Script 23.
- Confirm final source-of-truth heterogeneity Script 32 uses locked PTB and fixed reference IQR.
- Compare Figure S6 regenerated panel C and Table S1i continuous values after the same scale update.
- Verify exact figure/table numbering in manuscript and legends.
