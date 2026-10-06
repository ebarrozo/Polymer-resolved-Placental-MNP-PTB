# ============================================================
# 21_update_supplementary_figures_PTB_robustness.R
#
# Regenerates/organizes Supplementary Figure 4 and creates
# Supplementary Figure 7 for the PTB-MNP placenta manuscript.
#
# SUPPLEMENTARY FIGURE 4
#   Existing epidemiologic/model robustness:
#     a. all PTB sensitivity models
#     b. CLR/PCA outlier-exclusion sensitivity
#     c. leave-one-out OR range
#     d. leave-one-out detailed forest for key polymers
#     e. model-specification heatmap
#
# SUPPLEMENTARY FIGURE 7
#   New analytical robustness for PU/PET:
#     a. prevalence above analytical thresholds by delivery class
#     b. conventional adjusted vs Firth threshold models
#     c. >LOQ leave-one-out OR distributions
#     d. original vs conservative max-blank-subtracted burden
#
# Run from PROJECT_ROOT:
#   source("scripts/21_update_supplementary_figures_PTB_robustness.R")
#
# Prerequisites:
#   - current 06 primary PTB outputs
#   - current validation outputs from script 18
#   - existing scripts:
#       06c_primary_ptb_supp_sensitivity_figures_polished.R
#       07_ptb_robustness_figure4_polished.R
# ============================================================

source("00_config.R")

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
  library(ggplot2)
  library(scales)
})

source(file.path("R", "figure_theme_mnp.R"))

# ============================================================
# Settings
# ============================================================

RUN_EXISTING_SFIG4_SCRIPTS <- TRUE

SFIG4_DIR <- file.path(DIR_FIGS, "SupplementaryFigure4_updated")
SFIG7_DIR <- file.path(DIR_FIGS, "SupplementaryFigure7_analytical_robustness")

dir.create(SFIG4_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(SFIG7_DIR, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# Find latest validation output directory
# ------------------------------------------------------------

validation_dirs <- list.dirs(
  DIR_TABLES,
  recursive = FALSE,
  full.names = TRUE
)

validation_dirs <- validation_dirs[
  grepl("^18_placenta_spike_recovery_validation__", basename(validation_dirs))
]

if (length(validation_dirs) == 0) {
  stop(
    "No 18_placenta_spike_recovery_validation__* directory found under ",
    DIR_TABLES,
    ". Run script 18 first.",
    call. = FALSE
  )
}

VALIDATION_DIR <- validation_dirs[
  order(file.info(validation_dirs)$mtime, decreasing = TRUE)
][1]

message("Using validation directory: ", VALIDATION_DIR)

# ============================================================
# PART 1 — Regenerate and organize Supplementary Figure 4
# ============================================================

if (RUN_EXISTING_SFIG4_SCRIPTS) {

  script_06c <- file.path(
    "scripts",
    "06c_primary_ptb_supp_sensitivity_figures_polished.R"
  )

  script_07 <- file.path(
    "scripts",
    "07_ptb_robustness_figure4_polished.R"
  )

  if (!file.exists(script_06c)) {
    stop("Missing existing script: ", script_06c, call. = FALSE)
  }

  if (!file.exists(script_07)) {
    stop("Missing existing script: ", script_07, call. = FALSE)
  }

  message("")
  message("Regenerating existing Supplementary Figure 4 source panels...")
  source(script_06c)
  source(script_07)
}

# Existing polished panel outputs
sfig4_source <- c(
  a = file.path(
    DIR_FIGS,
    "supp_fig_primary_ptb_models__all_sensitivity_forest_polished.png"
  ),
  b = file.path(
    DIR_FIGS,
    "supp_fig_primary_ptb_models__outlier_sensitivity_focus_polished.png"
  ),
  c = file.path(
    DIR_FIGS,
    "fig4_ptb_robustness__leave_one_out_range_polished.png"
  ),
  d = file.path(
    DIR_FIGS,
    "fig4_ptb_robustness__leave_one_out_forest_polished.png"
  ),
  e = file.path(
    DIR_FIGS,
    "fig4_ptb_robustness__model_sensitivity_heatmap_polished.png"
  )
)

missing_sfig4 <- sfig4_source[!file.exists(sfig4_source)]

if (length(missing_sfig4) > 0) {
  stop(
    "Could not locate these expected Supplementary Fig. 4 panels:\n  ",
    paste(missing_sfig4, collapse = "\n  "),
    call. = FALSE
  )
}

sfig4_dest <- c(
  a = file.path(SFIG4_DIR, "SFig4a_all_PTB_sensitivity_models.png"),
  b = file.path(SFIG4_DIR, "SFig4b_CLR_PCA_outlier_sensitivity.png"),
  c = file.path(SFIG4_DIR, "SFig4c_leave_one_out_OR_range.png"),
  d = file.path(SFIG4_DIR, "SFig4d_leave_one_out_key_polymers.png"),
  e = file.path(SFIG4_DIR, "SFig4e_model_specification_heatmap.png")
)

for (nm in names(sfig4_source)) {
  file.copy(
    from = sfig4_source[[nm]],
    to = sfig4_dest[[nm]],
    overwrite = TRUE
  )
  message("Wrote: ", sfig4_dest[[nm]])
}

# ============================================================
# PART 2 — Supplementary Figure 7 analytical robustness
# ============================================================

read_val <- function(filename) {
  path <- file.path(VALIDATION_DIR, filename)
  if (!file.exists(path)) {
    stop(
      "Missing validation output: ", path,
      "\nRe-run the current script 18.",
      call. = FALSE
    )
  }
  read_csv(path, show_col_types = FALSE)
}

qc_2x2 <- read_val("pu_pet_qc_2x2_counts.csv")
qc_bin <- read_val("pu_pet_qc_binary_ptb_sensitivity_models.csv")
qc_firth <- read_val("pu_pet_qc_bias_reduced_models.csv")
loo_models <- read_val("pu_pet_aboveLOQ_leave_one_out_models.csv")
loo_summary <- read_val("pu_pet_aboveLOQ_leave_one_out_summary.csv")
maxblank_models <- read_val(
  "pu_pet_batchmatched_maxblank_subtracted_continuous_models.csv"
)
maxblank_audit <- read_val(
  "pu_pet_batchmatched_maxblank_subtraction_specimen_audit.csv"
)

poly_labs <- c(
  PU = "Polyurethane (PU)",
  PET = "Polyethylene terephthalate (PET)"
)

# Nature-style base theme, using project helper where practical.
supp_theme <- theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 11, hjust = 0),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 10),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 9),
    legend.title = element_blank(),
    legend.position = "bottom",
    plot.margin = margin(8, 8, 8, 8)
  )

save_png_pdf <- function(plot, stem, width, height) {
  png_path <- file.path(SFIG7_DIR, paste0(stem, ".png"))
  pdf_path <- file.path(SFIG7_DIR, paste0(stem, ".pdf"))

  ggsave(
    png_path, plot,
    width = width, height = height,
    units = "in", dpi = 600, bg = "white"
  )

  ggsave(
    pdf_path, plot,
    width = width, height = height,
    units = "in", bg = "white"
  )

  message("Wrote: ", png_path)
}

# ------------------------------------------------------------
# SFig7a — threshold prevalence by delivery class
# ------------------------------------------------------------

a_dat <- qc_2x2 %>%
  filter(
    Polymer %in% c("PU", "PET"),
    Delivery %in% c("Term", "Preterm"),
    Exposure == "Exposed"
  ) %>%
  mutate(
    Threshold = recode(
      Threshold,
      "at_or_above_LOD" = "≥LOD",
      "above_LOQ" = ">LOQ"
    ),
    Delivery = factor(Delivery, levels = c("Term", "Preterm")),
    Polymer = factor(
      Polymer,
      levels = c("PU", "PET"),
      labels = poly_labs[c("PU", "PET")]
    ),
    pct = case_when(
      Delivery == "Term" ~ 100 * n / 85,
      Delivery == "Preterm" ~ 100 * n / 86,
      TRUE ~ NA_real_
    )
  )

p7a <- ggplot(
  a_dat,
  aes(x = Threshold, y = pct, fill = Delivery)
) +
  geom_col(
    position = position_dodge(width = 0.75),
    width = 0.66
  ) +
  geom_text(
    aes(label = paste0(n, "/", ifelse(Delivery == "Term", 85, 86))),
    position = position_dodge(width = 0.75),
    vjust = -0.35,
    size = 3.0
  ) +
  facet_wrap(~ Polymer, nrow = 1) +
  scale_y_continuous(
    limits = c(0, 65),
    breaks = seq(0, 60, 10),
    labels = label_percent(scale = 1),
    expand = expansion(mult = c(0, 0.02))
  ) +
  scale_fill_manual(
    values = c(
      Term = "#A6A6A6",
      Preterm = "#4C78A8"
    )
  ) +
  labs(
    title = "PU and PET analytical-threshold prevalence by delivery class",
    x = NULL,
    y = "Placentae meeting threshold",
    fill = NULL
  ) +
  supp_theme

save_png_pdf(
  p7a,
  "SFig7a_PU_PET_threshold_prevalence",
  width = 7.2,
  height = 4.0
)

# ------------------------------------------------------------
# SFig7b — adjusted HC1 vs Firth threshold models
# ------------------------------------------------------------

b_hc1 <- qc_bin %>%
  filter(Polymer %in% c("PU", "PET")) %>%
  transmute(
    Polymer,
    Threshold = recode(
      Threshold,
      "at_or_above_LOD" = "≥LOD vs <LOD",
      "above_LOQ" = ">LOQ vs ≤LOQ"
    ),
    Method = "Adjusted logistic (HC1)",
    OR,
    CI_lo,
    CI_hi
  )

b_firth <- qc_firth %>%
  filter(Polymer %in% c("PU", "PET")) %>%
  transmute(
    Polymer,
    Threshold = recode(
      Threshold,
      "at_or_above_LOD" = "≥LOD vs <LOD",
      "above_LOQ" = ">LOQ vs ≤LOQ"
    ),
    Method = "Firth bias-reduced",
    OR,
    CI_lo,
    CI_hi
  )

b_dat <- bind_rows(b_hc1, b_firth) %>%
  mutate(
    Polymer = factor(
      Polymer,
      levels = c("PU", "PET"),
      labels = poly_labs[c("PU", "PET")]
    ),
    Threshold = factor(
      Threshold,
      levels = c("≥LOD vs <LOD", ">LOQ vs ≤LOQ")
    ),
    Method = factor(
      Method,
      levels = c("Adjusted logistic (HC1)", "Firth bias-reduced")
    )
  )

p7b <- ggplot(
  b_dat,
  aes(
    x = OR,
    y = Threshold,
    color = Method,
    shape = Method
  )
) +
  geom_vline(
    xintercept = 1,
    linetype = "dashed",
    linewidth = 0.45,
    color = "grey35"
  ) +
  geom_errorbarh(
    aes(xmin = CI_lo, xmax = CI_hi),
    position = position_dodge(width = 0.45),
    height = 0.18,
    linewidth = 0.75
  ) +
  geom_point(
    position = position_dodge(width = 0.45),
    size = 2.6
  ) +
  facet_wrap(~ Polymer, nrow = 1) +
  scale_x_log10(
    breaks = c(1, 2, 5, 10, 25),
    labels = c("1", "2", "5", "10", "25")
  ) +
  scale_color_manual(
    values = c(
      "Adjusted logistic (HC1)" = "black",
      "Firth bias-reduced" = "#D7301F"
    )
  ) +
  scale_shape_manual(
    values = c(
      "Adjusted logistic (HC1)" = 16,
      "Firth bias-reduced" = 17
    )
  ) +
  labs(
    title = "Threshold-based PU and PET PTB models are robust to bias reduction",
    x = "Odds ratio for preterm birth (95% CI)",
    y = NULL,
    color = NULL,
    shape = NULL
  ) +
  supp_theme +
  theme(
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.35)
  )

save_png_pdf(
  p7b,
  "SFig7b_PU_PET_HC1_vs_Firth",
  width = 7.8,
  height = 4.2
)

# ------------------------------------------------------------
# SFig7c — leave-one-out >LOQ OR distributions
# ------------------------------------------------------------

# Identify estimate / CI columns robustly.
or_col <- intersect(c("OR", "OR_estimate", "OR_threshold"), names(loo_models))[1]
lo_col <- intersect(c("CI_lo", "OR_lo", "CI_lower"), names(loo_models))[1]
hi_col <- intersect(c("CI_hi", "OR_hi", "CI_upper"), names(loo_models))[1]

if (is.na(or_col) || is.na(lo_col) || is.na(hi_col)) {
  stop(
    "Could not identify OR/CI columns in pu_pet_aboveLOQ_leave_one_out_models.csv.\n",
    "Columns found: ", paste(names(loo_models), collapse = ", "),
    call. = FALSE
  )
}

c_dat <- loo_models %>%
  filter(Polymer %in% c("PU", "PET")) %>%
  mutate(
    OR_plot = .data[[or_col]],
    CI_lo_plot = .data[[lo_col]],
    CI_hi_plot = .data[[hi_col]],
    Polymer = factor(
      Polymer,
      levels = c("PU", "PET"),
      labels = poly_labs[c("PU", "PET")]
    )
  )

# Full-model adjusted estimates, not Firth, for direct comparison.
c_full <- qc_bin %>%
  filter(
    Polymer %in% c("PU", "PET"),
    Threshold == "above_LOQ"
  ) %>%
  mutate(
    Polymer = factor(
      Polymer,
      levels = c("PU", "PET"),
      labels = poly_labs[c("PU", "PET")]
    )
  )

p7c <- ggplot(
  c_dat,
  aes(x = OR_plot, y = Polymer)
) +
  geom_vline(
    xintercept = 1,
    linetype = "dashed",
    linewidth = 0.45,
    color = "grey35"
  ) +
  geom_jitter(
    height = 0.14,
    width = 0,
    size = 1.45,
    alpha = 0.30
  ) +
  geom_errorbarh(
    data = c_full,
    aes(
      x = OR,
      xmin = CI_lo,
      xmax = CI_hi,
      y = Polymer
    ),
    inherit.aes = FALSE,
    height = 0.10,
    linewidth = 1.0,
    color = "#D7301F"
  ) +
  geom_point(
    data = c_full,
    aes(x = OR, y = Polymer),
    inherit.aes = FALSE,
    size = 3.0,
    color = "#D7301F"
  ) +
  scale_x_log10(
    breaks = c(1, 2, 5, 10, 15),
    labels = c("1", "2", "5", "10", "15")
  ) +
  labs(
    title = "Adjusted >LOQ associations are stable to leave-one-out analysis",
    subtitle = str_wrap(
      "Gray points show one-participant-deleted estimates; red points and bars show the full adjusted model and 95% CI",
      width = 82
    ),
    x = "Odds ratio for preterm birth",
    y = NULL
  ) +
  supp_theme +
  theme(
    legend.position = "none",
    panel.grid.major.y = element_blank(),
    plot.title = element_text(face = "bold", size = 10.5, hjust = 0),
    plot.subtitle = element_text(size = 9, hjust = 0, margin = margin(b = 5))
  )

save_png_pdf(
  p7c,
  "SFig7c_PU_PET_aboveLOQ_leave_one_out",
  width = 7.2,
  height = 3.4
)

# ------------------------------------------------------------
# SFig7d — specimen-level effect of conservative max-blank subtraction
# ------------------------------------------------------------

d_dat <- maxblank_audit %>%
  filter(Polymer %in% c("PU", "PET")) %>%
  mutate(
    Delivery = factor(
      PTB01,
      levels = c(0, 1),
      labels = c("Term", "Preterm")
    ),
    Polymer = factor(
      Polymer,
      levels = c("PU", "PET"),
      labels = poly_labs[c("PU", "PET")]
    )
  )

# Use log1p coordinates because zero is expected after subtraction.
p7d <- ggplot(
  d_dat,
  aes(
    x = log1p(original_burden_ug_g),
    y = log1p(blank_subtracted_ug_g),
    shape = Delivery
  )
) +
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dashed",
    linewidth = 0.5,
    color = "grey40"
  ) +
  geom_point(
    size = 1.8,
    alpha = 0.65
  ) +
  facet_wrap(~ Polymer, nrow = 1, scales = "free") +
  scale_shape_manual(values = c(Term = 16, Preterm = 17)) +
  labs(
    title = str_wrap(
      "Conservative batch-matched max-blank subtraction reduces burden without eliminating PU/PET exposure contrasts",
      width = 82
    ),
    subtitle = str_wrap(
      "Axes show log(1 + µg/g); dashed line indicates no change after subtraction",
      width = 82
    ),
    x = "Original placental burden, log(1 + µg/g)",
    y = "Max-blank-subtracted burden, log(1 + µg/g)",
    shape = NULL
  ) +
  supp_theme +
  theme(
    plot.title = element_text(face = "bold", size = 10.5, hjust = 0),
    plot.subtitle = element_text(size = 9, hjust = 0, margin = margin(b = 5))
  )

save_png_pdf(
  p7d,
  "SFig7d_PU_PET_maxblank_specimen_sensitivity",
  width = 7.5,
  height = 4.2
)

# ============================================================
# Write draft legends
# ============================================================

legend_sfig4 <- paste0(
  "Supplementary Figure 4. Robustness of polymer-specific associations with preterm birth across model specifications. ",
  "a, Sensitivity analysis of polymer-specific associations with preterm birth across crude and adjusted covariate specifications. ",
  "b, Outlier-exclusion sensitivity analysis comparing the primary adjusted model with the model excluding the two CLR-PCA composition outliers. ",
  "c, Leave-one-out analysis showing the range of odds-ratio estimates after sequential exclusion of individual subjects for selected polymers. ",
  "d, Leave-one-out robustness analysis for key preterm-birth-associated polymers. ",
  "e, Heatmap of preterm birth association estimates across covariate specifications. ",
  "Odds ratios for continuous burden models are reported per doubling of placental polymer burden."
)

legend_sfig7 <- paste0(
  "Supplementary Figure 7. Analytical-threshold and blank-background sensitivity analyses of PU and PET associations with preterm birth. ",
  "a, Proportions of term and preterm placentae with PU or PET measurements at or above the polymer-specific LOD and above the LOQ; labels show exposed samples over the total delivery-class denominator. ",
  "b, Adjusted analytical-threshold models comparing conventional logistic regression with HC1 robust standard errors and bias-reduced Firth logistic regression. ",
  "c, Leave-one-out robustness of adjusted >LOQ models. Gray points show odds-ratio estimates after sequential exclusion of one complete-case participant; red points and horizontal bars show the full adjusted model and 95% confidence interval. ",
  "d, Specimen-level PU and PET burden before and after conservative subtraction of the batch-matched polymer-specific maximum procedural blank contribution once per analytical run. ",
  "Threshold and sensitivity models used the same primary clinical covariate set where applicable."
)

writeLines(
  c(
    "SUPPLEMENTARY FIGURE 4",
    legend_sfig4,
    "",
    "SUPPLEMENTARY FIGURE 7",
    legend_sfig7
  ),
  file.path(SFIG7_DIR, "draft_SFig4_and_SFig7_legends.txt")
)

# ============================================================
# Run notes
# ============================================================

notes <- c(
  paste0("Validation directory: ", VALIDATION_DIR),
  paste0("SFig4 output directory: ", SFIG4_DIR),
  paste0("SFig7 output directory: ", SFIG7_DIR),
  "SFig4 is regenerated from the existing polished scripts to preserve the manuscript's prior visual style.",
  "SFig7 is new and focuses on analytical robustness rather than ordinary covariate/model sensitivity.",
  "The blank-background figure is intentionally not included here; it remains reserved for later review.",
  "SFig7d uses the batch-matched max-blank sensitivity analysis only; it should not be described as a corrected primary exposure."
)

writeLines(
  notes,
  file.path(SFIG7_DIR, "supplementary_figure_update_notes.txt")
)

message("")
message("Supplementary figure update complete.")
message("SFig4 panels: ", SFIG4_DIR)
message("SFig7 panels: ", SFIG7_DIR)
