#!/usr/bin/env Rscript
# POLISHED VERSION: sources R/figure_theme_mnp.R and writes *_polished outputs where applicable.
# ============================================================
# 06b_primary_ptb_figure3_clean.R
# Cleaner Figure 3 remake from the latest primary PTB model outputs
# v2.5: wrapped subtitle + only selected sensitivity models
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
  library(ggplot2)
  library(scales)
  library(stringr)
})

source("00_config.R")
source(file.path("R", "figure_theme_mnp.R"))

# ----------------------------
# Settings
# ----------------------------
PRIMARY_MODEL_TAG <- "M1_ADJ_PRIMARY"
OUTLIER_SENS_TAG <- "M1_ADJ_PRIMARY_NO_CLR_OUTLIERS"

# Main-figure sensitivity models only
DISPLAY_SENS_MODELS <- c(
  "M0_CRUDE",
  OUTLIER_SENS_TAG,
  "M1_PLUS_HTNsubtype"
)

X_LOWER <- 0.03
X_UPPER <- 4.0

# Keep primary model black. Give each sensitivity model a distinct color.
COLOR_MAP <- c(
  "Crude" = "#56B4E9",
  "Primary adjusted, excluding CLR/PCA outliers" = "#7A3E9D",
  "Primary + HTN subtype" = "#E69F00"
)

OUT_FIG <- file.path(DIR_FIGS, "fig_primary_ptb_models__figure3_clean_polished.png")
OUT_SPEC <- file.path(DIR_TABLES, "table_primary_model_spec_polished.csv")
OUT_LEGEND <- file.path(DIR_TABLES, "figure3_legend_revised_polished.txt")

# ----------------------------
# Find latest model directory
# ----------------------------
model_dirs <- list.dirs(DIR_TABLES, recursive = FALSE, full.names = TRUE)
model_dirs <- model_dirs[grepl("^06_primary_ptb_models__", basename(model_dirs))]

if (length(model_dirs) == 0) {
  stop("No 06_primary_ptb_models__* directory found in outputs/tables/. Run run_primary_models.R first.")
}

latest_dir <- model_dirs[order(file.info(model_dirs)$mtime, decreasing = TRUE)][1]
all_models_path <- file.path(latest_dir, "PTB__ALL_models_long.csv")
primary_path <- file.path(latest_dir, "table_primary_ptb_models.csv")

all_models <- readr::read_csv(all_models_path, show_col_types = FALSE)
primary_tbl <- readr::read_csv(primary_path, show_col_types = FALSE)

# ----------------------------
# Primary model specification table
# ----------------------------
primary_covs <- all_models %>%
  filter(ModelTag == PRIMARY_MODEL_TAG) %>%
  pull(Covariates) %>%
  unique()

primary_covs <- primary_covs[!is.na(primary_covs) & nzchar(primary_covs)]
primary_covs_chr <- if (length(primary_covs) > 0) primary_covs[1] else NA_character_

spec_tbl <- tibble(
  ModelTag = PRIMARY_MODEL_TAG,
  ModelLabel = "Primary adjusted",
  Outcome = "PTB01 (1 = Preterm, 0 = Term)",
  ExposureScale = "OR per doubling in placental burden from log10(µg/g + 1) exposure models",
  PrimaryCovariates = primary_covs_chr,
  RowsWithOutcome = 171L,
  CompleteCase_N_Range = paste0(min(primary_tbl$N, na.rm = TRUE), " to ", max(primary_tbl$N, na.rm = TRUE)),
  OutlierHandling = "Primary model retains all eligible samples. Sensitivity model excludes 2 CLR/PCA composition outliers identified during clustering.",
  Note = "Gestational age was not included as a covariate because it defines PTB."
)

write_csv(spec_tbl, OUT_SPEC)

# ----------------------------
# Data prep
# ----------------------------
fig_dat <- all_models %>%
  filter(ModelTag %in% c(PRIMARY_MODEL_TAG, DISPLAY_SENS_MODELS)) %>%
  mutate(
    ModelLabel = factor(ModelLabel),
    sig_q05 = ifelse(is.na(q_BH), FALSE, q_BH < 0.05),
    PolymerLabel = case_when(
      Polymer == "Cumulative_MNPs" ~ "Cumulative MNP burden",
      TRUE ~ PolymerLabel
    )
  )

poly_order_primary <- fig_dat %>%
  filter(ModelTag == PRIMARY_MODEL_TAG, Polymer != "Cumulative_MNPs") %>%
  arrange(desc(OR_2x)) %>%
  pull(PolymerLabel)

poly_levels <- c(poly_order_primary[1:4], "Cumulative MNP burden", poly_order_primary[5:length(poly_order_primary)])
poly_levels <- poly_levels[!is.na(poly_levels)]

base_positions <- tibble(
  PolymerLabel = poly_levels,
  y_base = rev(seq_along(poly_levels))
)

primary_df <- fig_dat %>%
  filter(ModelTag == PRIMARY_MODEL_TAG) %>%
  left_join(base_positions, by = "PolymerLabel") %>%
  mutate(
    xmin = pmax(OR_2x_lo, X_LOWER),
    xmax = pmin(OR_2x_hi, X_UPPER),
    left_trunc = OR_2x_lo < X_LOWER,
    right_trunc = OR_2x_hi > X_UPPER
  )

offset_map <- tibble(
  ModelTag = DISPLAY_SENS_MODELS,
  y_offset = c(-0.18, 0.00, 0.18)
)

sens_df <- fig_dat %>%
  filter(ModelTag %in% DISPLAY_SENS_MODELS) %>%
  left_join(base_positions, by = "PolymerLabel") %>%
  left_join(offset_map, by = "ModelTag") %>%
  mutate(y = y_base + y_offset)

# Subtitle wrapping
subtitle_txt <- paste(
  "Primary model: maternal age + BMI + SDI + twins + prior PTB + hypertension + hospital.",
  "Black points and 95% CI bars = primary adjusted model."
)

subtitle_wrapped <- stringr::str_wrap(subtitle_txt, width = 82)

# ----------------------------
# Plot
# ----------------------------
p <- ggplot() +
  geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.4, color = "grey30") +

  # Primary model CI bars
  geom_segment(
    data = primary_df,
    aes(x = xmin, xend = xmax, y = y_base, yend = y_base),
    linewidth = 0.85,
    color = "black"
  ) +

  # Truncation arrows for primary CIs
  geom_text(
    data = primary_df %>% filter(left_trunc),
    aes(x = X_LOWER, y = y_base, label = "\u2190"),
    hjust = 1.05,
    vjust = 0.35,
    size = 4.0
  ) +
  geom_text(
    data = primary_df %>% filter(right_trunc),
    aes(x = X_UPPER, y = y_base, label = "\u2192"),
    hjust = -0.10,
    vjust = 0.35,
    size = 4.0
  ) +

  # Sensitivity model points only
  geom_point(
    data = sens_df,
    aes(x = OR_2x, y = y, color = ModelLabel),
    size = 2.3,
    alpha = 0.95
  ) +

  # Primary model point, shaped by q < 0.05
  geom_point(
    data = primary_df,
    aes(x = OR_2x, y = y_base, shape = sig_q05),
    size = 3.0,
    color = "black",
    fill = "black"
  ) +

  scale_x_log10(
    breaks = c(0.03, 0.05, 0.10, 0.20, 0.50, 1.00, 2.00, 4.00),
    labels = c("0.03", "0.05", "0.10", "0.20", "0.50", "1.00", "2.00", "4.00")
  ) +
  coord_cartesian(xlim = c(X_LOWER, X_UPPER)) +
  scale_y_continuous(
    breaks = base_positions$y_base,
    labels = base_positions$PolymerLabel,
    expand = expansion(mult = c(0.03, 0.03))
  ) +
  scale_color_manual(
    values = COLOR_MAP,
    breaks = c("Crude", "Primary adjusted, excluding CLR/PCA outliers", "Primary + HTN subtype"),
    labels = c("Crude", "Primary adjusted, excluding CLR/PCA outliers", "Primary + HTN subtype")
  ) +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17)) +
  labs(
    title = "Polymer-specific associations with preterm birth",
    subtitle = subtitle_wrapped,
    x = "OR per doubling (95% CI)",
    y = NULL,
    color = "Sensitivity models",
    shape = "Primary model\nBH q < 0.05"
  ) +
  guides(
    color = guide_legend(order = 1, override.aes = list(size = 3)),
    shape = guide_legend(order = 2)
  ) +
  theme_mnp_forest(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 10.8, lineheight = 1.05),
    legend.position = "right",
    legend.box = "vertical"
  )

ggsave_mnp(OUT_FIG, p, width = 8.6, height = 5.7, dpi = 300)

# ----------------------------
# Revised legend text
# ----------------------------
legend_text <- paste0(
  "Figure 3 | Polymer-specific associations with preterm birth and sensitivity to covariate adjustment.\n",
  "Polymer-specific associations with preterm birth (PTB; delivery at gestational age <37.0 weeks) are shown as odds ratios (ORs) per doubling of placental polymer burden with 95% confidence intervals. ",
  "Black points and horizontal bars represent the primary adjusted model (M1), which included maternal age, BMI, Social Deprivation Index (SDI), twin gestation, prior preterm birth, hypertension, and hospital. ",
  "Colored points indicate sensitivity models (crude model, the primary model excluding the two CLR/PCA composition outliers identified during clustering, and the primary model using hypertension subtype adjustment). ",
  "Point shape denotes BH–FDR significance for the primary model (q < 0.05). ",
  "The x-axis is shown on a log scale. Arrows indicate confidence intervals truncated for display."
)

writeLines(legend_text, OUT_LEGEND)

message("Wrote:\n", OUT_FIG, "\n", OUT_SPEC, "\n", OUT_LEGEND)
