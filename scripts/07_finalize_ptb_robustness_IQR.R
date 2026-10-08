#!/usr/bin/env Rscript
# ============================================================
# 07_finalize_ptb_robustness_IQR.R
#
# Finalizes Supplementary Figure 4 after reviewer Point 1.
#
# Uses:
#   outputs/tables/24_ptb_sensitivity_IQR/PTB_ALL_models_long_IQR.csv
#   latest outputs/tables/07_ptb_robustness__*/ robustness tables
#   latest outputs/tables/22_ptb_effect_scale_audit__*/ primary IQR refs
#
# Panels:
#   a. all PTB model specifications, all analytes
#   b. primary vs CLR/PCA outlier-exclusion sensitivity
#   c. leave-one-out OR range for key polymers
#   d. detailed leave-one-out estimates for key associated polymers
#   e. model-specification heatmap
#
# The same polymer-specific PRIMARY-cohort IQR is used for all PTB
# sensitivity and leave-one-out comparisons. This is intentional:
# it keeps the effect unit fixed while assessing coefficient stability.
#
# Does NOT overwrite legacy SFig4 files.
# ============================================================

source("00_config.R")

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
  library(stringr)
  library(ggplot2)
  library(scales)
})

find_latest_dir <- function(parent, prefix) {
  x <- list.dirs(parent, recursive = FALSE, full.names = TRUE)
  x <- x[grepl(paste0("^", prefix), basename(x))]
  if (!length(x)) return(NA_character_)
  x[order(file.info(x)$mtime, decreasing = TRUE)][1]
}

OUT_DIR <- file.path(DIR_TABLES, "27_SuppFigure4_IQR_final")
OUT_FIG_DIR <- file.path(DIR_FIGS, "27_SuppFigure4_IQR_final")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_FIG_DIR, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# Inputs
# ------------------------------------------------------------

all_iqr_path <- file.path(DIR_TABLES, "24_ptb_sensitivity_IQR", "PTB_ALL_models_long_IQR.csv")
if (!file.exists(all_iqr_path)) {
  stop("Missing ", all_iqr_path, ". Run Point 1 propagation first.", call. = FALSE)
}
all_iqr <- read_csv(all_iqr_path, show_col_types = FALSE)

audit_dir <- find_latest_dir(DIR_TABLES, "22_ptb_effect_scale_audit__")
if (is.na(audit_dir)) stop("No Point 1 audit directory found.", call. = FALSE)
iqr_ref <- read_csv(file.path(audit_dir, "PTB_primary_IQR_reference_values.csv"), show_col_types = FALSE)

find_latest_dir_with_files <- function(parent, prefix, required_files) {
  x <- list.dirs(parent, recursive = FALSE, full.names = TRUE)
  x <- x[grepl(paste0("^", prefix), basename(x))]
  if (!length(x)) return(NA_character_)

  ok <- vapply(
    x,
    function(d) all(file.exists(file.path(d, required_files))),
    logical(1)
  )
  x <- x[ok]

  if (!length(x)) return(NA_character_)
  x[order(file.info(x)$mtime, decreasing = TRUE)][1]
}

rob_required <- c(
  "robustness_leave_one_out_models.csv",
  "robustness_leave_one_out_summary.csv"
)

rob_dir <- find_latest_dir_with_files(
  DIR_TABLES,
  "07_ptb_robustness__",
  rob_required
)

if (is.na(rob_dir)) {
  # The newest 07 directory can be an incomplete/polish-only run.
  # If no complete robustness directory exists, regenerate the canonical
  # robustness tables once using the non-polish analysis script.
  canonical_07 <- file.path("scripts", "internal", "07_ptb_robustness_fit_stage.R")

  if (file.exists(canonical_07)) {
    message(
      "No complete 07 robustness directory containing leave-one-out tables was found.\n",
      "Running canonical robustness analysis once: ", canonical_07
    )
    source(canonical_07)

    rob_dir <- find_latest_dir_with_files(
      DIR_TABLES,
      "07_ptb_robustness__",
      rob_required
    )
  }
}

if (is.na(rob_dir)) {
  stop(
    "Could not locate or generate a complete outputs/tables/07_ptb_robustness__* ",
    "directory containing robustness_leave_one_out_models.csv and ",
    "robustness_leave_one_out_summary.csv.",
    call. = FALSE
  )
}

message("Using complete robustness directory: ", rob_dir)

loo_path <- file.path(rob_dir, "robustness_leave_one_out_models.csv")
loo_sum_path <- file.path(rob_dir, "robustness_leave_one_out_summary.csv")

loo <- read_csv(loo_path, show_col_types = FALSE)
loo_sum_old <- read_csv(loo_sum_path, show_col_types = FALSE)

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------

K2 <- log10(2)

poly_label <- function(poly, fallback = poly) {
  if (poly == "Cumulative_MNPs") return("Total MNP burden")
  fallback
}

rescale_old_2x <- function(df, or_col="OR_2x", lo_col="OR_2x_lo", hi_col="OR_2x_hi") {
  df %>%
    left_join(iqr_ref %>% select(Polymer, IQR_log10_1px), by="Polymer") %>%
    mutate(
      beta_recovered = log(.data[[or_col]]) / K2,
      beta_lo_recovered = log(.data[[lo_col]]) / K2,
      beta_hi_recovered = log(.data[[hi_col]]) / K2,
      se_beta_recovered = (beta_hi_recovered - beta_lo_recovered)/(2*1.96),
      logOR_IQR = beta_recovered * IQR_log10_1px,
      OR_IQR = exp(logOR_IQR),
      OR_IQR_lo = exp((beta_recovered - 1.96*se_beta_recovered)*IQR_log10_1px),
      OR_IQR_hi = exp((beta_recovered + 1.96*se_beta_recovered)*IQR_log10_1px)
    )
}

loo_iqr <- rescale_old_2x(loo)
write_csv(loo_iqr, file.path(OUT_DIR, "robustness_leave_one_out_models_IQR.csv"))

# Recompute LOO summary from rescaled estimates rather than trying to
# transform old min/max summaries.
loo_summary <- loo_iqr %>%
  group_by(Polymer, PolymerLabel) %>%
  summarise(
    n_iterations = sum(is.finite(OR_IQR)),
    OR_min = min(OR_IQR, na.rm=TRUE),
    OR_median = median(OR_IQR, na.rm=TRUE),
    OR_max = max(OR_IQR, na.rm=TRUE),
    CIlo_min = min(OR_IQR_lo, na.rm=TRUE),
    CIhi_max = max(OR_IQR_hi, na.rm=TRUE),
    .groups="drop"
  )
write_csv(loo_summary, file.path(OUT_DIR, "robustness_leave_one_out_summary_IQR.csv"))

# QA: primary-cohort IQR refs should be present for every LOO polymer.
qa <- loo_iqr %>%
  group_by(Polymer) %>%
  summarise(
    n=n(),
    n_missing_IQR=sum(!is.finite(IQR_log10_1px)),
    .groups="drop"
  )
write_csv(qa, file.path(OUT_DIR, "SFig4_IQR_QA.csv"))
if (any(qa$n_missing_IQR > 0)) stop("Missing IQR reference for one or more LOO polymers.", call. = FALSE)

# ------------------------------------------------------------
# Display constants
# ------------------------------------------------------------

MODEL_ORDER <- c(
  "Crude",
  "Primary adjusted",
  "Primary without SDI",
  "Primary + SVI",
  "Primary + race/ethnicity",
  "Primary + HTN subtype",
  "Primary adjusted, excluding CLR/PCA outliers"
)

MODEL_COLORS <- c(
  "Crude" = "#56B4E9",
  "Primary adjusted" = "black",
  "Primary without SDI" = "#999999",
  "Primary + SVI" = "#009E73",
  "Primary + race/ethnicity" = "#CC79A7",
  "Primary + HTN subtype" = "#E69F00",
  "Primary adjusted, excluding CLR/PCA outliers" = "#7A3E9D"
)

POLY_ORDER <- all_iqr %>%
  filter(ModelTag == "M1_ADJ_PRIMARY") %>%
  mutate(PolymerLabel = if_else(Polymer=="Cumulative_MNPs", "Total MNP burden", PolymerLabel)) %>%
  arrange(desc(OR_IQR)) %>%
  pull(PolymerLabel)

base_theme <- theme_bw(base_size=10.5) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face="bold", size=11),
    plot.subtitle = element_text(size=9.2),
    legend.title = element_text(face="bold"),
    strip.background = element_rect(fill="grey92", color="grey55"),
    strip.text = element_text(face="bold")
  )

save_both <- function(p, stem, w, h) {
  ggsave(file.path(OUT_FIG_DIR, paste0(stem, ".png")), p, width=w, height=h, dpi=300)
  ggsave(file.path(OUT_FIG_DIR, paste0(stem, ".pdf")), p, width=w, height=h)
}

# ============================================================
# PANEL A — all model specifications, all analytes
# ============================================================

a_dat <- all_iqr %>%
  mutate(
    PolymerLabel = if_else(Polymer=="Cumulative_MNPs", "Total MNP burden", PolymerLabel),
    PolymerLabel = factor(PolymerLabel, levels=rev(POLY_ORDER)),
    ModelLabel = factor(ModelLabel, levels=MODEL_ORDER)
  ) %>%
  filter(!is.na(ModelLabel))

# Draw primary CI and all model estimates. This is much less cramped than
# faceting each polymer into a narrow strip.
pA <- ggplot(a_dat, aes(x=OR_IQR, y=PolymerLabel, color=ModelLabel)) +
  geom_vline(xintercept=1, linetype="dashed", linewidth=0.35, color="grey40") +
  geom_errorbar(
    data=a_dat %>% filter(ModelTag=="M1_ADJ_PRIMARY"),
    aes(xmin=OR_IQR_lo, xmax=OR_IQR_hi),
    orientation="y", width=0.18, linewidth=0.7,
    color="black", inherit.aes=TRUE
  ) +
  geom_point(
    position=position_dodge(width=0.48),
    size=2.0, alpha=0.95
  ) +
  scale_x_log10() +
  scale_color_manual(values=MODEL_COLORS, drop=FALSE) +
  labs(
    title="a  PTB associations across model specifications",
    subtitle="Primary-model 95% CIs are shown in black; colored points show alternative specifications",
    x="OR per polymer-specific IQR increase in log10(1 + burden)",
    y=NULL,
    color="Model"
  ) +
  base_theme +
  theme(legend.position="right", panel.grid.major.y=element_blank())

save_both(pA, "SFig4a_all_PTB_sensitivity_models_IQR", 10.2, 6.9)

# ============================================================
# PANEL B — primary vs CLR/PCA outlier exclusion
# ============================================================

b_dat <- all_iqr %>%
  filter(ModelTag %in% c("M1_ADJ_PRIMARY", "M1_ADJ_PRIMARY_NO_CLR_OUTLIERS")) %>%
  mutate(
    PolymerLabel=if_else(Polymer=="Cumulative_MNPs","Total MNP burden",PolymerLabel),
    PolymerLabel=factor(PolymerLabel, levels=rev(POLY_ORDER)),
    ModelLabel=factor(
      ModelLabel,
      levels=c("Primary adjusted","Primary adjusted, excluding CLR/PCA outliers")
    )
  )

pB <- ggplot(b_dat, aes(x=OR_IQR, y=PolymerLabel, color=ModelLabel)) +
  geom_vline(xintercept=1, linetype="dashed", linewidth=0.35, color="grey40") +
  geom_errorbar(
    aes(xmin=OR_IQR_lo, xmax=OR_IQR_hi),
    orientation="y",
    position=position_dodge(width=0.5),
    width=0.18, linewidth=0.55
  ) +
  geom_point(position=position_dodge(width=0.5), size=2.4) +
  scale_x_log10() +
  scale_color_manual(values=c(
    "Primary adjusted"="black",
    "Primary adjusted, excluding CLR/PCA outliers"="#7A3E9D"
  )) +
  labs(
    title="b  CLR/PCA outlier-exclusion sensitivity",
    x="OR per polymer-specific IQR increase in log10(1 + burden)",
    y=NULL,
    color=NULL
  ) +
  base_theme +
  theme(legend.position="bottom", panel.grid.major.y=element_blank())

save_both(pB, "SFig4b_CLR_PCA_outlier_sensitivity_IQR", 8.8, 6.3)

# ============================================================
# PANEL C — LOO range
# ============================================================

c_dat <- loo_summary %>%
  mutate(
    PolymerLabel=if_else(Polymer=="Cumulative_MNPs","Total MNP burden",PolymerLabel),
    PolymerLabel=factor(PolymerLabel, levels=PolymerLabel[order(OR_median)])
  )

pC <- ggplot(c_dat, aes(y=PolymerLabel)) +
  geom_vline(xintercept=1, linetype="dashed", linewidth=0.35, color="grey40") +
  geom_segment(aes(x=OR_min, xend=OR_max, yend=PolymerLabel), linewidth=1.0, alpha=0.65) +
  geom_point(aes(x=OR_median), size=2.8) +
  scale_x_log10() +
  labs(
    title="c  Range of leave-one-out PTB estimates",
    subtitle="Segments show minimum to maximum one-participant-deleted estimate; points show median",
    x="OR per polymer-specific IQR increase in log10(1 + burden)",
    y=NULL
  ) +
  base_theme +
  theme(panel.grid.major.y=element_blank())

save_both(pC, "SFig4c_leave_one_out_OR_range_IQR", 7.4, 4.3)

# ============================================================
# PANEL D — detailed LOO distributions
# ============================================================

KEY_ASSOC <- c("PU","PET","SBR","ABS")

primary_key <- all_iqr %>%
  filter(ModelTag=="M1_ADJ_PRIMARY", Polymer %in% KEY_ASSOC) %>%
  transmute(
    Polymer,
    primary_OR=OR_IQR,
    primary_lo=OR_IQR_lo,
    primary_hi=OR_IQR_hi
  )

d_dat <- loo_iqr %>%
  filter(Polymer %in% KEY_ASSOC) %>%
  left_join(primary_key, by="Polymer") %>%
  mutate(
    PolymerLabel=case_when(
      Polymer=="PU" ~ "Polyurethane (PU)",
      Polymer=="PET" ~ "Polyethylene terephthalate (PET)",
      Polymer=="SBR" ~ "Styrene-butadiene rubber (SBR)",
      Polymer=="ABS" ~ "Acrylonitrile-butadiene-styrene (ABS)",
      TRUE ~ Polymer
    )
  )

pD <- ggplot(d_dat, aes(x=OR_IQR, y=1)) +
  geom_vline(xintercept=1, linetype="dashed", linewidth=0.35, color="grey45") +
  geom_jitter(height=0.07, width=0, alpha=0.32, size=1.15, color="grey35") +
  geom_errorbar(
    data=d_dat %>% distinct(PolymerLabel, primary_OR, primary_lo, primary_hi),
    aes(xmin=primary_lo, xmax=primary_hi, y=1),
    orientation="y", width=0.14, linewidth=0.85,
    color="black",
    inherit.aes=FALSE
  ) +
  geom_point(
    data=d_dat %>% distinct(PolymerLabel, primary_OR, primary_lo, primary_hi),
    aes(x=primary_OR, y=1),
    size=3.0, color="black",
    inherit.aes=FALSE
  ) +
  facet_wrap(~PolymerLabel, ncol=2, scales="free_x") +
  scale_x_log10() +
  scale_y_continuous(NULL, breaks=NULL) +
  labs(
    title="d  Leave-one-out robustness of PTB-associated polymers",
    subtitle="Gray points are one-participant-deleted estimates; black point and interval show the primary model",
    x="OR per polymer-specific IQR increase in log10(1 + burden)",
    y=NULL
  ) +
  base_theme

save_both(pD, "SFig4d_leave_one_out_key_polymers_IQR", 8.7, 5.2)

# ============================================================
# PANEL E — model-specification heatmap
# ============================================================

e_dat <- all_iqr %>%
  mutate(
    PolymerLabel=if_else(Polymer=="Cumulative_MNPs","Total MNP burden",PolymerLabel),
    PolymerLabel=factor(PolymerLabel, levels=rev(POLY_ORDER)),
    ModelLabel=factor(ModelLabel, levels=MODEL_ORDER),
    logOR_IQR=log(OR_IQR),
    evidence=case_when(
      !is.na(q_BH) & q_BH < 0.05 ~ "BH q < 0.05",
      !is.na(p) & p < 0.05 ~ "Nominal p < 0.05",
      TRUE ~ "Not nominal"
    )
  ) %>%
  filter(!is.na(ModelLabel))

lim <- max(abs(e_dat$logOR_IQR), na.rm=TRUE)

pE <- ggplot(e_dat, aes(x=ModelLabel, y=PolymerLabel, fill=logOR_IQR)) +
  geom_tile(color="white", linewidth=0.35) +
  geom_point(
    data=e_dat %>% filter(evidence != "Not nominal"),
    aes(shape=evidence),
    size=2.0,
    color="black"
  ) +
  scale_fill_gradient2(
    midpoint=0,
    limits=c(-lim,lim),
    name="log(OR)\nper IQR"
  ) +
  scale_shape_manual(
    values=c("Nominal p < 0.05"=1, "BH q < 0.05"=16),
    name="Evidence"
  ) +
  labs(
    title="e  PTB association stability across model specifications",
    x=NULL, y=NULL
  ) +
  base_theme +
  theme(
    axis.text.x=element_text(angle=38, hjust=1),
    legend.position="right"
  )

save_both(pE, "SFig4e_model_specification_heatmap_IQR", 9.6, 6.4)

# ============================================================
# Combined figure if patchwork is installed
# ============================================================

combined_created <- FALSE
if (requireNamespace("patchwork", quietly=TRUE)) {
  combined <- (pA / pB / (pC | pD) / pE) +
    patchwork::plot_layout(heights=c(1.25,1.1,0.95,1.15)) +
    patchwork::plot_annotation(
      title="Supplementary Figure 4 | Robustness of polymer-specific associations with preterm birth",
      theme=theme(plot.title=element_text(face="bold", size=14))
    )

  save_both(combined, "SupplementaryFigure4_IQR_FINAL", 12.0, 23.0)
  combined_created <- TRUE
}

# ============================================================
# Legend
# ============================================================

legend_txt <- paste0(
  "Supplementary Figure 4. Robustness of polymer-specific associations with preterm birth across model specifications. ",
  "a, Polymer-specific PTB estimates across crude and adjusted covariate specifications. Continuous effects are expressed as odds ratios (ORs) per polymer-specific interquartile-range (IQR) increase in log10(1 + placental polymer burden [µg/g]); black horizontal intervals show 95% confidence intervals for the primary adjusted model. ",
  "b, Outlier-exclusion sensitivity analysis comparing the primary adjusted model with the corresponding model after exclusion of the two CLR/PCA composition outliers. ",
  "c, Range of OR estimates obtained from leave-one-out analyses after sequential exclusion of individual participants for PU, PET, SBR, ABS, and cumulative MNP burden; segments show the minimum-to-maximum estimate and points show the median. ",
  "d, Detailed leave-one-out robustness analysis for PU, PET, SBR, and ABS. Gray points show one-participant-deleted estimates; black points and horizontal intervals show the primary-model estimate and 95% confidence interval. ",
  "e, Heatmap of polymer-specific PTB association estimates across alternative covariate specifications. Color indicates log(OR) per IQR and symbols denote nominal or Benjamini-Hochberg false-discovery-rate-adjusted statistical support. ",
  "For sensitivity and leave-one-out comparisons, each polymer is displayed using the same IQR defined in the primary PTB analysis cohort so that differences across models reflect coefficient stability rather than changes in the reporting unit."
)

writeLines(legend_txt, file.path(OUT_DIR, "SupplementaryFigure4_legend_IQR_FINAL.txt"))

notes <- c(
  paste0("Robustness source directory: ", rob_dir),
  paste0("Point 1 audit directory: ", audit_dir),
  paste0("Combined figure created: ", combined_created),
  "",
  "Decision:",
  "This script supersedes the earlier 24_ptb_sensitivity_IQR candidate figure.",
  "The earlier narrow-facet candidate should not be used in the manuscript.",
  "Panels C/D use the actual robustness_leave_one_out_models.csv output from script 07.",
  "All PTB sensitivity/LOO displays use a fixed primary-cohort polymer-specific IQR."
)
writeLines(notes, file.path(OUT_DIR, "SFig4_finalization_notes.txt"))

message("Supplementary Figure 4 finalization complete.")
message("Tables/legend: ", OUT_DIR)
message("Panels/combined figure: ", OUT_FIG_DIR)
