#!/usr/bin/env Rscript
# POLISHED VERSION: sources R/figure_theme_mnp.R and writes *_polished outputs where applicable.
# ============================================================
# 06c_primary_ptb_supp_sensitivity_figures.R
# Supplementary sensitivity figures for primary PTB models
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

PRIMARY_TAG <- "M1_ADJ_PRIMARY"

# Full sensitivity set for supplement
SENS_TAGS_ALL <- c(
  "M0_CRUDE",
  "M1_ADJ_PRIMARY_NO_CLR_OUTLIERS",
  "M1_PLUS_SVI",
  "M1_NO_SDI",
  "M1_PLUS_HTNsubtype",
  "M1_PLUS_RACE"
)

MODEL_ORDER <- c(PRIMARY_TAG, SENS_TAGS_ALL)

MODEL_LABELS_SHORT <- c(
  "M0_CRUDE" = "Crude",
  "M1_ADJ_PRIMARY" = "Primary adjusted",
  "M1_ADJ_PRIMARY_NO_CLR_OUTLIERS" = "Primary adjusted,\nexcluding CLR/PCA outliers",
  "M1_PLUS_SVI" = "Primary + SVI",
  "M1_NO_SDI" = "Primary without SDI",
  "M1_PLUS_HTNsubtype" = "Primary + HTN subtype",
  "M1_PLUS_RACE" = "Primary + race/ethnicity"
)

MODEL_COLORS <- c(
  "Crude" = "#56B4E9",
  "Primary adjusted,\nexcluding CLR/PCA outliers" = "#7A3E9D",
  "Primary + SVI" = "#009E73",
  "Primary without SDI" = "#D55E00",
  "Primary + HTN subtype" = "#E69F00",
  "Primary + race/ethnicity" = "#CC79A7"
)

OUT_FIG_ALL <- file.path(DIR_FIGS, "supp_fig_primary_ptb_models__all_sensitivity_forest_polished.png")
OUT_FIG_OUTLIERS <- file.path(DIR_FIGS, "supp_fig_primary_ptb_models__outlier_sensitivity_focus_polished.png")
OUT_MODEL_KEY <- file.path(DIR_TABLES, "supp_table_primary_ptb_model_key_polished.csv")
OUT_LEGENDS <- file.path(DIR_TABLES, "supp_legend_primary_ptb_sensitivity_figures_polished.txt")

X_LOWER <- 0.03
X_UPPER_ALL <- 6.0
X_UPPER_OUTLIERS <- 6.0

# ----------------------------
# Find latest primary model outputs
# ----------------------------
model_dirs <- list.dirs(DIR_TABLES, recursive = FALSE, full.names = TRUE)
model_dirs <- model_dirs[grepl("^06_primary_ptb_models__", basename(model_dirs))]

if (length(model_dirs) == 0) {
  stop("No 06_primary_ptb_models__* directory found in outputs/tables/. Run run_primary_models.R first.")
}

latest_dir <- model_dirs[order(file.info(model_dirs)$mtime, decreasing = TRUE)][1]
all_models_path <- file.path(latest_dir, "PTB__ALL_models_long.csv")

dat <- readr::read_csv(all_models_path, show_col_types = FALSE)

# ----------------------------
# Model key table
# ----------------------------
model_key <- dat %>%
  filter(ModelTag %in% MODEL_ORDER) %>%
  group_by(ModelTag, ModelLabel) %>%
  summarise(
    Covariates = first(Covariates),
    N_range = paste0(min(N, na.rm = TRUE), " to ", max(N, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(
    ModelLabelShort = recode(ModelTag, !!!MODEL_LABELS_SHORT),
    Description = case_when(
      ModelTag == "M0_CRUDE" ~ "Unadjusted model",
      ModelTag == "M1_ADJ_PRIMARY" ~ "Primary adjusted model",
      ModelTag == "M1_ADJ_PRIMARY_NO_CLR_OUTLIERS" ~ "Primary model excluding 2 CLR/PCA composition outliers",
      ModelTag == "M1_PLUS_SVI" ~ "Primary model plus SVI",
      ModelTag == "M1_NO_SDI" ~ "Primary model without SDI",
      ModelTag == "M1_PLUS_HTNsubtype" ~ "Primary model with hypertension subtype",
      ModelTag == "M1_PLUS_RACE" ~ "Primary model plus race/ethnicity",
      TRUE ~ NA_character_
    )
  ) %>%
  select(ModelTag, ModelLabelShort, Description, Covariates, N_range)

write_csv(model_key, OUT_MODEL_KEY)

# ----------------------------
# Common data prep
# ----------------------------
plot_dat <- dat %>%
  filter(ModelTag %in% MODEL_ORDER) %>%
  mutate(
    PolymerLabel = case_when(
      Polymer == "Cumulative_MNPs" ~ "Cumulative MNP burden",
      TRUE ~ PolymerLabel
    ),
    sig_q05 = ifelse(is.na(q_BH), FALSE, q_BH < 0.05),
    ModelLabelShort = recode(ModelTag, !!!MODEL_LABELS_SHORT)
  )

poly_order <- plot_dat %>%
  filter(ModelTag == PRIMARY_TAG, Polymer != "Cumulative_MNPs") %>%
  arrange(desc(OR_2x)) %>%
  pull(PolymerLabel)

poly_levels <- c(poly_order[1:4], "Cumulative MNP burden", poly_order[5:length(poly_order)])
poly_levels <- poly_levels[!is.na(poly_levels)]

base_positions <- tibble(
  PolymerLabel = poly_levels,
  y_base = rev(seq_along(poly_levels))
)

primary_df <- plot_dat %>%
  filter(ModelTag == PRIMARY_TAG) %>%
  left_join(base_positions, by = "PolymerLabel") %>%
  mutate(
    xmin = pmax(OR_2x_lo, X_LOWER),
    xmax = pmin(OR_2x_hi, X_UPPER_ALL),
    left_trunc = OR_2x_lo < X_LOWER,
    right_trunc = OR_2x_hi > X_UPPER_ALL
  )

# ----------------------------
# Supplement figure 1:
# all sensitivity models
# ----------------------------
offset_map_all <- tibble(
  ModelTag = SENS_TAGS_ALL,
  y_offset = c(-0.25, -0.15, -0.05, 0.05, 0.15, 0.25)
)

sens_all <- plot_dat %>%
  filter(ModelTag %in% SENS_TAGS_ALL) %>%
  left_join(base_positions, by = "PolymerLabel") %>%
  left_join(offset_map_all, by = "ModelTag") %>%
  mutate(y = y_base + y_offset)

subtitle_all <- str_wrap(
  "Primary adjusted model is shown as black points and 95% CI bars. Colored points indicate sensitivity models.",
  width = 92
)

p_all <- ggplot() +
  geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.4, color = "grey30") +
  geom_segment(
    data = primary_df,
    aes(x = xmin, xend = xmax, y = y_base, yend = y_base),
    linewidth = 0.8,
    color = "black"
  ) +
  geom_text(
    data = primary_df %>% filter(left_trunc),
    aes(x = X_LOWER, y = y_base, label = "\u2190"),
    hjust = 1.05, vjust = 0.35, size = 3.8
  ) +
  geom_text(
    data = primary_df %>% filter(right_trunc),
    aes(x = X_UPPER_ALL, y = y_base, label = "\u2192"),
    hjust = -0.05, vjust = 0.35, size = 3.8
  ) +
  geom_point(
    data = sens_all,
    aes(x = OR_2x, y = y, color = ModelLabelShort),
    size = 1.9, alpha = 0.95
  ) +
  geom_point(
    data = primary_df,
    aes(x = OR_2x, y = y_base, shape = sig_q05),
    size = 2.7, color = "black", fill = "black"
  ) +
  scale_x_log10(
    breaks = c(0.03, 0.05, 0.10, 0.20, 0.50, 1.00, 2.00, 4.00, 6.00),
    labels = c("0.03", "0.05", "0.10", "0.20", "0.50", "1.00", "2.00", "4.00", "6.00")
  ) +
  coord_cartesian(xlim = c(X_LOWER, X_UPPER_ALL)) +
  scale_y_continuous(
    breaks = base_positions$y_base,
    labels = base_positions$PolymerLabel,
    expand = expansion(mult = c(0.03, 0.03))
  ) +
  scale_color_manual(values = MODEL_COLORS) +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17)) +
  labs(
    title = "Sensitivity analysis of polymer-specific PTB associations",
    subtitle = subtitle_all,
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
    plot.subtitle = element_text(size = 10.5, lineheight = 1.05),
    legend.position = "right",
    legend.box = "vertical"
  )

ggsave_mnp(OUT_FIG_ALL, p_all, width = 9.2, height = 6.2, dpi = 300)

# ----------------------------
# Supplement figure 2:
# focused outlier sensitivity
# ----------------------------
focus_tags <- c(PRIMARY_TAG, "M1_ADJ_PRIMARY_NO_CLR_OUTLIERS")
focus_colors <- c(
  "Primary adjusted" = "#000000",
  "Primary adjusted,\nexcluding CLR/PCA outliers" = "#7A3E9D"
)

offset_map_focus <- tibble(
  ModelTag = focus_tags,
  y_offset = c(-0.10, 0.10)
)

focus_df <- plot_dat %>%
  filter(ModelTag %in% focus_tags) %>%
  left_join(base_positions, by = "PolymerLabel") %>%
  left_join(offset_map_focus, by = "ModelTag") %>%
  mutate(
    y = y_base + y_offset,
    xmin = pmax(OR_2x_lo, X_LOWER),
    xmax = pmin(OR_2x_hi, X_UPPER_OUTLIERS),
    left_trunc = OR_2x_lo < X_LOWER,
    right_trunc = OR_2x_hi > X_UPPER_OUTLIERS
  )

subtitle_outliers <- str_wrap(
  "Comparison of the primary adjusted model with the model excluding the 2 CLR/PCA composition outliers identified during clustering.",
  width = 92
)

p_out <- ggplot() +
  geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.4, color = "grey30") +
  geom_segment(
    data = focus_df,
    aes(x = xmin, xend = xmax, y = y, yend = y, color = ModelLabelShort),
    linewidth = 0.85
  ) +
  geom_text(
    data = focus_df %>% filter(left_trunc),
    aes(x = X_LOWER, y = y, label = "\u2190", color = ModelLabelShort),
    hjust = 1.05, vjust = 0.35, size = 3.8, show.legend = FALSE
  ) +
  geom_text(
    data = focus_df %>% filter(right_trunc),
    aes(x = X_UPPER_OUTLIERS, y = y, label = "\u2192", color = ModelLabelShort),
    hjust = -0.05, vjust = 0.35, size = 3.8, show.legend = FALSE
  ) +
  geom_point(
    data = focus_df,
    aes(x = OR_2x, y = y, color = ModelLabelShort),
    size = 2.5
  ) +
  scale_x_log10(
    breaks = c(0.03, 0.05, 0.10, 0.20, 0.50, 1.00, 2.00, 4.00, 6.00),
    labels = c("0.03", "0.05", "0.10", "0.20", "0.50", "1.00", "2.00", "4.00", "6.00")
  ) +
  coord_cartesian(xlim = c(X_LOWER, X_UPPER_OUTLIERS)) +
  scale_y_continuous(
    breaks = base_positions$y_base,
    labels = base_positions$PolymerLabel,
    expand = expansion(mult = c(0.03, 0.03))
  ) +
  scale_color_manual(values = focus_colors) +
  labs(
    title = "Outlier-exclusion sensitivity analysis",
    subtitle = subtitle_outliers,
    x = "OR per doubling (95% CI)",
    y = NULL,
    color = NULL
  ) +
  theme_mnp_forest(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 10.5, lineheight = 1.05),
    legend.position = "right"
  )

ggsave_mnp(OUT_FIG_OUTLIERS, p_out, width = 8.7, height = 5.8, dpi = 300)

# ----------------------------
# Supplement legends
# ----------------------------
legend_text <- paste0(
  "Supplementary Figure A | Sensitivity analysis of polymer-specific associations with preterm birth.\n",
  "Polymer-specific associations with preterm birth (PTB; delivery at gestational age <37.0 weeks) are shown as odds ratios (ORs) per doubling of placental polymer burden with 95% confidence intervals. ",
  "Black points and horizontal bars represent the primary adjusted model, which included maternal age, BMI, Social Deprivation Index (SDI), twin gestation, prior preterm birth, hypertension, and hospital. ",
  "Colored points represent sensitivity models including the crude model, the primary model excluding the two CLR/PCA composition outliers identified during clustering, the primary model plus SVI, the primary model without SDI, the primary model using hypertension subtype adjustment, and the primary model plus race/ethnicity. ",
  "Point shape denotes BH–FDR significance for the primary model (q < 0.05). The x-axis is shown on a log scale; arrows indicate confidence intervals truncated for display.\n\n",
  "Supplementary Figure B | Outlier-exclusion sensitivity analysis.\n",
  "Primary adjusted and outlier-exclusion sensitivity models are shown as odds ratios (ORs) per doubling of placental polymer burden with 95% confidence intervals. ",
  "The outlier-exclusion model removes the two CLR/PCA composition outliers identified during clustering. The x-axis is shown on a log scale; arrows indicate confidence intervals truncated for display."
)

writeLines(legend_text, OUT_LEGENDS)

message("Wrote:\n", OUT_FIG_ALL, "\n", OUT_FIG_OUTLIERS, "\n", OUT_MODEL_KEY, "\n", OUT_LEGENDS)
