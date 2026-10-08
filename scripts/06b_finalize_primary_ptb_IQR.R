#!/usr/bin/env Rscript
# ============================================================
# 06b_finalize_primary_ptb_IQR.R
#
# Finalize reviewer Point 1 for the PRIMARY PTB analysis:
#   - preserve the fitted log10(1 + ug/g) models
#   - express continuous effects per polymer-specific IQR increase
#     in log10(1 + ug/g)
#   - regenerate final Figure 3 and primary-model table/legend
#
# Requires:
#   scripts/06a_ptb_effect_scale_reference.R has been run successfully.
#
# Run from PROJECT_ROOT:
#   source("scripts/06b_finalize_primary_ptb_IQR.R")
#
# Does NOT overwrite legacy Figure 3/model outputs.
# ============================================================

source("00_config.R")

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
  library(ggplot2)
  library(stringr)
  library(scales)
})

PRIMARY_MODEL_TAG <- "M1_ADJ_PRIMARY"
DISPLAY_SENS_MODELS <- c(
  "M0_CRUDE",
  "M1_ADJ_PRIMARY_NO_CLR_OUTLIERS",
  "M1_PLUS_HTNsubtype"
)

COLOR_MAP <- c(
  "Crude" = "#56B4E9",
  "Primary adjusted, excluding CLR/PCA outliers" = "#7A3E9D",
  "Primary + HTN subtype" = "#E69F00"
)

find_latest_dir <- function(parent, prefix) {
  dirs <- list.dirs(parent, recursive = FALSE, full.names = TRUE)
  dirs <- dirs[grepl(paste0("^", prefix), basename(dirs))]
  if (length(dirs) == 0) stop("No directory beginning with ", prefix, call. = FALSE)
  dirs[order(file.info(dirs)$mtime, decreasing = TRUE)][1]
}

AUDIT_DIR <- find_latest_dir(DIR_TABLES, "22_ptb_effect_scale_audit__")
message("Using Point 1 audit directory: ", AUDIT_DIR)

all_iqr_path <- file.path(AUDIT_DIR, "PTB_ALL_models_IQR_rescaled.csv")
primary_iqr_path <- file.path(AUDIT_DIR, "PTB_effect_scale_comparison_primary.csv")

stopifnot(file.exists(all_iqr_path), file.exists(primary_iqr_path))

all_models <- read_csv(all_iqr_path, show_col_types = FALSE)
primary_tbl <- read_csv(primary_iqr_path, show_col_types = FALSE)

# ------------------------------------------------------------
# Output locations
# ------------------------------------------------------------

OUT_DIR <- file.path(DIR_TABLES, "23_primary_ptb_IQR_final")
OUT_FIG_DIR <- file.path(DIR_FIGS, "23_primary_ptb_IQR_final")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_FIG_DIR, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# Final primary model table
# ------------------------------------------------------------

primary_final <- all_models %>%
  filter(ModelTag == PRIMARY_MODEL_TAG) %>%
  transmute(
    ModelTag,
    ModelLabel,
    Polymer,
    PolymerLabel,
    ExposureCol,
    N,
    Events,
    NonEvents,
    Covariates,
    IQR_log10_1px,
    OR_IQR,
    OR_IQR_lo,
    OR_IQR_hi,
    logOR_IQR,
    SE_logOR_IQR,
    p,
    q_BH,
    OR_IQR_CI,
    sig_q05 = ifelse(is.na(q_BH), FALSE, q_BH < 0.05)
  )

write_csv(primary_final, file.path(OUT_DIR, "table_primary_ptb_models_IQR.csv"))

# ------------------------------------------------------------
# Figure 3
# ------------------------------------------------------------

fig_dat <- all_models %>%
  filter(ModelTag %in% c(PRIMARY_MODEL_TAG, DISPLAY_SENS_MODELS)) %>%
  mutate(
    sig_q05 = ifelse(is.na(q_BH), FALSE, q_BH < 0.05),
    PolymerLabel = if_else(
      Polymer == "Cumulative_MNPs",
      "Total MNP burden",
      PolymerLabel
    )
  )

poly_order_primary <- fig_dat %>%
  filter(ModelTag == PRIMARY_MODEL_TAG, Polymer != "Cumulative_MNPs") %>%
  arrange(desc(OR_IQR)) %>%
  pull(PolymerLabel)

if (length(poly_order_primary) >= 4) {
  poly_levels <- c(
    poly_order_primary[1:4],
    "Total MNP burden",
    poly_order_primary[5:length(poly_order_primary)]
  )
} else {
  poly_levels <- c(poly_order_primary, "Total MNP burden")
}
poly_levels <- unique(poly_levels[!is.na(poly_levels)])

base_positions <- tibble(
  PolymerLabel = poly_levels,
  y_base = rev(seq_along(poly_levels))
)

primary_df <- fig_dat %>%
  filter(ModelTag == PRIMARY_MODEL_TAG) %>%
  left_join(base_positions, by = "PolymerLabel")

offset_map <- tibble(
  ModelTag = DISPLAY_SENS_MODELS,
  y_offset = c(-0.18, 0.00, 0.18)
)

sens_df <- fig_dat %>%
  filter(ModelTag %in% DISPLAY_SENS_MODELS) %>%
  left_join(base_positions, by = "PolymerLabel") %>%
  left_join(offset_map, by = "ModelTag") %>%
  mutate(y = y_base + y_offset)

# Robust display range
all_lo <- primary_df$OR_IQR_lo[is.finite(primary_df$OR_IQR_lo) & primary_df$OR_IQR_lo > 0]
all_hi <- primary_df$OR_IQR_hi[is.finite(primary_df$OR_IQR_hi) & primary_df$OR_IQR_hi > 0]
X_LOWER <- if (length(all_lo)) max(0.02, 10^(floor(log10(min(all_lo))*2)/2)) else 0.05
X_UPPER <- if (length(all_hi)) min(20, 10^(ceiling(log10(max(all_hi))*2)/2)) else 4

subtitle_txt <- "Primary model: maternal age + BMI + SDI + twins + prior PTB + hypertension + hospital."

p <- ggplot() +
  geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.4, color = "grey30") +
  geom_segment(
    data = primary_df,
    aes(x = OR_IQR_lo, xend = OR_IQR_hi, y = y_base, yend = y_base),
    linewidth = 0.85,
    color = "black"
  ) +
  geom_point(
    data = sens_df,
    aes(x = OR_IQR, y = y, color = ModelLabel),
    size = 2.3,
    alpha = 0.95
  ) +
  geom_point(
    data = primary_df,
    aes(x = OR_IQR, y = y_base, shape = sig_q05),
    size = 3.0,
    color = "black",
    fill = "black"
  ) +
  scale_x_log10() +
  coord_cartesian(xlim = c(X_LOWER, X_UPPER)) +
  scale_y_continuous(
    breaks = base_positions$y_base,
    labels = base_positions$PolymerLabel,
    expand = expansion(mult = c(0.03, 0.03))
  ) +
  scale_color_manual(
    values = COLOR_MAP,
    breaks = c(
      "Crude",
      "Primary adjusted, excluding CLR/PCA outliers",
      "Primary + HTN subtype"
    )
  ) +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17)) +
  labs(
    title = "Polymer-specific associations with preterm birth",
    subtitle = subtitle_txt,
    x = "OR per IQR increase in log10(1 + burden) (95% CI)",
    y = NULL,
    color = "Sensitivity models",
    shape = "Primary model\nBH q < 0.05"
  ) +
  guides(
    color = guide_legend(order = 1, override.aes = list(size = 3)),
    shape = guide_legend(order = 2)
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 10.8, lineheight = 1.05),
    legend.position = "right",
    legend.box = "vertical"
  )

OUT_PNG <- file.path(OUT_FIG_DIR, "Figure3_primary_ptb_IQR_final.png")
OUT_PDF <- file.path(OUT_FIG_DIR, "Figure3_primary_ptb_IQR_final.pdf")
ggsave(OUT_PNG, p, width = 8.6, height = 5.7, dpi = 300)
ggsave(OUT_PDF, p, width = 8.6, height = 5.7)

# ------------------------------------------------------------
# Figure legend and model specification
# ------------------------------------------------------------

legend_text <- paste0(
  "Figure 3 | Polymer-specific associations with preterm birth and sensitivity to covariate adjustment. ",
  "Polymer-specific associations with preterm birth (PTB; delivery at gestational age <37.0 weeks) ",
  "are shown as odds ratios (ORs) per polymer-specific interquartile-range (IQR) increase in ",
  "log10(1 + placental polymer burden [µg/g]), with 95% confidence intervals. ",
  "Black points and horizontal bars represent the primary adjusted model, which included maternal age, ",
  "BMI, Social Deprivation Index (SDI), twin gestation, prior preterm birth, hypertension, and hospital. ",
  "Colored points indicate sensitivity models (crude model, the primary model excluding the two CLR/PCA ",
  "composition outliers identified during clustering, and the primary model using hypertension subtype ",
  "adjustment). Point shape denotes BH-FDR significance for the primary model (q<0.05). ",
  "The x-axis is shown on a log scale."
)
writeLines(legend_text, file.path(OUT_DIR, "Figure3_legend_IQR_final.txt"))

primary_covs <- all_models %>%
  filter(ModelTag == PRIMARY_MODEL_TAG) %>%
  pull(Covariates) %>%
  unique()
primary_covs <- primary_covs[!is.na(primary_covs) & nzchar(primary_covs)]

spec_tbl <- tibble(
  ModelTag = PRIMARY_MODEL_TAG,
  ModelLabel = "Primary adjusted",
  Outcome = "PTB01 (1 = Preterm, 0 = Term)",
  ExposureModel = "log10(1 + placental burden [µg/g])",
  DisplayScale = "OR per polymer-specific IQR increase in log10(1 + burden)",
  PrimaryCovariates = ifelse(length(primary_covs), primary_covs[[1]], NA_character_),
  Note = "Fitted model unchanged from prior primary analysis; only effect-size reporting scale was corrected."
)
write_csv(spec_tbl, file.path(OUT_DIR, "table_primary_model_spec_IQR.csv"))

# ------------------------------------------------------------
# Manuscript replacement text
# ------------------------------------------------------------

key <- primary_final %>%
  filter(Polymer %in% c("Cumulative_MNPs", "PU", "PET", "SBR", "ABS")) %>%
  mutate(
    p_fmt = if_else(p < 0.001, "<0.001", sprintf("%.4f", p)),
    q_fmt = if_else(q_BH < 0.001, "<0.001", sprintf("%.3f", q_BH))
  )

getrow <- function(poly) key %>% filter(Polymer == poly) %>% slice(1)

fmt_est <- function(r) {
  sprintf("%.2f, 95%% CI %.2f-%.2f, p = %s, q = %s",
          r$OR_IQR, r$OR_IQR_lo, r$OR_IQR_hi, r$p_fmt, r$q_fmt)
}

r_cum <- getrow("Cumulative_MNPs")
r_pu  <- getrow("PU")
r_pet <- getrow("PET")
r_sbr <- getrow("SBR")
r_abs <- getrow("ABS")

results_text <- paste0(
  "Cumulative MNP burden was not associated with PTB (OR per IQR = ",
  sprintf("%.2f", r_cum$OR_IQR), ", 95% CI ",
  sprintf("%.2f", r_cum$OR_IQR_lo), "-", sprintf("%.2f", r_cum$OR_IQR_hi),
  ", p = ", sprintf("%.3f", r_cum$p), "). ",
  "In contrast, higher PU and PET burdens were associated with increased odds of PTB ",
  "(PU: OR per IQR = ", fmt_est(r_pu), "; PET: OR per IQR = ", fmt_est(r_pet), "). ",
  "SBR and ABS were inversely associated with PTB ",
  "(SBR: OR per IQR = ", fmt_est(r_sbr), "; ABS: OR per IQR = ", fmt_est(r_abs), ")."
)

methods_text <- paste0(
  "PTB associations were modeled using logistic regression with log10(1 + µg/g) polymer burden as the exposure. ",
  "Continuous-model effect estimates were expressed as odds ratios per polymer-specific interquartile-range increase ",
  "in the transformed polymer burden, with 95% confidence intervals. HC1 heteroskedasticity-robust standard errors ",
  "were used for primary and sensitivity logistic models."
)

writeLines(
  c(
    "METHODS REPLACEMENT",
    methods_text,
    "",
    "RESULTS REPLACEMENT",
    results_text
  ),
  file.path(OUT_DIR, "manuscript_replacement_text_point1_primary.txt")
)

message("Finalized primary PTB IQR reporting.")
message("Tables: ", OUT_DIR)
message("Figures: ", OUT_FIG_DIR)
