#!/usr/bin/env Rscript
# ============================================================
# 08b_finalize_ordinal_labor_figures_IQR.R
#
# Finalize reviewer Point 1 for:
#   A) ordinal gestational-age category models
#   B) spontaneous labor augmentation models
#
# Uses the latest outputs/tables/26_secondary_models_IQR__* tables.
# Does NOT overwrite legacy outputs.
# ============================================================

source("00_config.R")

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
  library(ggplot2)
  library(stringr)
})

find_latest_dir <- function(parent, prefix) {
  x <- list.dirs(parent, recursive = FALSE, full.names = TRUE)
  x <- x[grepl(paste0("^", prefix), basename(x))]
  if (!length(x)) return(NA_character_)
  x[order(file.info(x)$mtime, decreasing = TRUE)][1]
}

SRC_DIR <- find_latest_dir(DIR_TABLES, "26_secondary_models_IQR__")
if (is.na(SRC_DIR)) stop("No 26_secondary_models_IQR__* directory found.", call. = FALSE)

OUT_DIR <- file.path(DIR_TABLES, "28_ordinal_labor_IQR_final")
OUT_FIG_DIR <- file.path(DIR_FIGS, "28_ordinal_labor_IQR_final")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_FIG_DIR, recursive = TRUE, showWarnings = FALSE)

save_both <- function(p, stem, w, h) {
  ggsave(file.path(OUT_FIG_DIR, paste0(stem, ".png")), p, width=w, height=h, dpi=300)
  ggsave(file.path(OUT_FIG_DIR, paste0(stem, ".pdf")), p, width=w, height=h)
}

theme_nc <- theme_bw(base_size=10.5) +
  theme(
    panel.grid.minor=element_blank(),
    plot.title=element_text(face="bold"),
    strip.background=element_rect(fill="grey92", color="grey55"),
    strip.text=element_text(face="bold")
  )

# ============================================================
# A. ORDINAL GESTATIONAL-AGE MODELS
# ============================================================

ord_path <- file.path(SRC_DIR, "ordinal_ga_category_all_models_IQR.csv")
if (file.exists(ord_path)) {
  ord <- read_csv(ord_path, show_col_types=FALSE)

  # QA
  if ("N_match" %in% names(ord) && any(ord$N_match == FALSE, na.rm=TRUE)) {
    stop("Ordinal GA N reconstruction mismatch detected.", call. = FALSE)
  }

  model_levels <- c(
    "Crude",
    "Primary adjusted",
    "Primary without SDI",
    "Primary + SVI",
    "Primary + race/ethnicity"
  )

  ord <- ord %>%
    mutate(
      PolymerLabel = if_else(
        Polymer == "Cumulative_MNPs",
        "Total MNP burden",
        PolymerLabel
      ),
      ModelLabel = factor(ModelLabel, levels=model_levels)
    )

  primary <- ord %>% filter(ModelTag=="M1_PRIMARY", is.finite(OR_earlier_IQR))
  sens <- ord %>% filter(ModelTag!="M1_PRIMARY", is.finite(OR_earlier_IQR))

  poly_order <- primary %>% arrange(desc(OR_earlier_IQR)) %>% pull(PolymerLabel) %>% unique()
  primary <- primary %>% mutate(PolymerLabel=factor(PolymerLabel, levels=rev(poly_order)))
  sens <- sens %>% mutate(PolymerLabel=factor(PolymerLabel, levels=rev(poly_order)))

  p_ord <- ggplot() +
    geom_vline(xintercept=1, linetype="dashed", linewidth=0.35, color="grey40") +
    geom_errorbar(
      data=primary,
      aes(xmin=OR_earlier_IQR_lo, xmax=OR_earlier_IQR_hi, y=PolymerLabel),
      orientation="y", width=0.17, linewidth=0.75, color="black"
    ) +
    geom_point(
      data=primary,
      aes(x=OR_earlier_IQR, y=PolymerLabel),
      size=2.7, color="black"
    ) +
    geom_point(
      data=sens,
      aes(x=OR_earlier_IQR, y=PolymerLabel, color=ModelLabel),
      position=position_dodge(width=0.45),
      size=1.9
    ) +
    scale_x_log10() +
    scale_color_manual(values=c(
      "Crude"="#56B4E9",
      "Primary without SDI"="#999999",
      "Primary + SVI"="#009E73",
      "Primary + race/ethnicity"="#CC79A7"
    ), na.translate=FALSE) +
    labs(
      title="Ordinal gestational-age category models",
      subtitle="Odds ratios >1 indicate greater odds of delivery in an earlier gestational-age category",
      x="OR per analysis-specific IQR increase in log10(1 + burden) (95% CI)",
      y=NULL,
      color="Sensitivity model"
    ) +
    theme_nc +
    theme(legend.position="right", panel.grid.major.y=element_blank())

  save_both(p_ord, "ordinal_GA_category_models_IQR_FINAL", 9.2, 6.6)

  ord_primary <- primary %>%
    select(
      ModelTag, Polymer, PolymerLabel, N, IQR_transformed,
      OR_earlier_IQR, OR_earlier_IQR_lo, OR_earlier_IQR_hi,
      p, any_of(c("q_BH_all_analytes","fit_engine","fit_status"))
    )
  write_csv(ord_primary, file.path(OUT_DIR, "ordinal_GA_primary_IQR_FINAL.csv"))

  ord_legend <- paste0(
    "Supplementary Figure | Ordinal gestational-age category models. ",
    "Gestational age at delivery was modeled as an ordered categorical outcome: very/extremely preterm (<32 weeks), ",
    "moderate preterm (32 to <34 weeks), late preterm (34 to <37 weeks), early term (37 to <39 weeks), ",
    "full term (39 to 41 weeks), and late term (>41 weeks). Cumulative-link ordinal regression was used, and ",
    "estimates are shown as odds ratios (ORs) per analysis-specific interquartile-range (IQR) increase in ",
    "log10(1 + placental polymer burden [µg/g]). ORs >1 indicate greater odds of delivery in an earlier ",
    "gestational-age category. The primary adjusted model included maternal age, BMI, Social Deprivation Index, ",
    "twin gestation, prior preterm birth, hypertension, and hospital. Colored points show sensitivity models. ",
    "The x-axis is log scaled."
  )
  writeLines(ord_legend, file.path(OUT_DIR, "ordinal_GA_figure_legend_IQR_FINAL.txt"))
}

# ============================================================
# B. LABOR AUGMENTATION MODELS
# ============================================================

lab_path <- file.path(SRC_DIR, "labor_model_all_results_IQR.csv")
if (file.exists(lab_path)) {
  lab <- read_csv(lab_path, show_col_types=FALSE)

  if ("N_match" %in% names(lab) && any(lab$N_match == FALSE, na.rm=TRUE)) {
    stop("Labor N reconstruction mismatch detected.", call. = FALSE)
  }

  model_levels <- c(
    "Crude",
    "GA-adjusted primary",
    "GA + hospital",
    "GA + SDI",
    "GA + prior PTB",
    "GA + hypertension"
  )

  lab <- lab %>%
    mutate(
      PolymerLabel=if_else(Polymer=="Cumulative_MNPs","Total MNP burden",PolymerLabel),
      ModelLabel=factor(ModelLabel, levels=model_levels)
    )

  primary <- lab %>% filter(ModelTag=="GA_ONLY", is.finite(OR_IQR))
  sens <- lab %>% filter(ModelTag!="GA_ONLY", is.finite(OR_IQR))

  poly_order <- primary %>% arrange(desc(OR_IQR)) %>% pull(PolymerLabel) %>% unique()
  primary <- primary %>% mutate(PolymerLabel=factor(PolymerLabel, levels=rev(poly_order)))
  sens <- sens %>% mutate(PolymerLabel=factor(PolymerLabel, levels=rev(poly_order)))

  primary <- primary %>%
    mutate(
      evidence=case_when(
        !is.na(q_BH_all_analytes) & q_BH_all_analytes < 0.05 ~ "BH q<0.05, all analytes",
        !is.na(q_BH_labor_prespec) & q_BH_labor_prespec < 0.05 ~ "BH q<0.05, prespecified PVC/N66",
        !is.na(p) & p < 0.05 ~ "nominal p<0.05",
        TRUE ~ "not nominal"
      )
    )

  p_lab <- ggplot() +
    geom_vline(xintercept=1, linetype="dashed", linewidth=0.35, color="grey40") +
    geom_errorbar(
      data=primary,
      aes(xmin=OR_IQR_lo, xmax=OR_IQR_hi, y=PolymerLabel),
      orientation="y", width=0.17, linewidth=0.75, color="black"
    ) +
    geom_point(
      data=primary,
      aes(x=OR_IQR, y=PolymerLabel, shape=evidence),
      size=2.8, color="black"
    ) +
    geom_point(
      data=sens,
      aes(x=OR_IQR, y=PolymerLabel, color=ModelLabel),
      position=position_dodge(width=0.45),
      size=1.8
    ) +
    scale_x_log10() +
    scale_color_manual(values=c(
      "Crude"="#56B4E9",
      "GA + hospital"="#E69F00",
      "GA + SDI"="#009E73",
      "GA + prior PTB"="#CC79A7",
      "GA + hypertension"="#7A3E9D"
    ), na.translate=FALSE) +
    scale_shape_manual(values=c(
      "not nominal"=16,
      "nominal p<0.05"=1,
      "BH q<0.05, prespecified PVC/N66"=23,
      "BH q<0.05, all analytes"=17
    )) +
    labs(
      title="Polymer-specific associations with spontaneous labor requiring augmentation",
      subtitle="Exploratory primary model adjusted for gestational age at delivery",
      x="OR per analysis-specific IQR increase in log10(1 + burden) (95% CI)",
      y=NULL,
      color="Sensitivity model",
      shape="Evidence"
    ) +
    theme_nc +
    theme(legend.position="right", panel.grid.major.y=element_blank())

  save_both(p_lab, "labor_augmentation_models_IQR_FINAL", 9.3, 6.7)

  lab_primary <- primary %>%
    select(
      ModelTag, Polymer, PolymerLabel, N, Events_Augmented,
      NonEvents_SpontaneousNoAug, IQR_transformed,
      OR_IQR, OR_IQR_lo, OR_IQR_hi, p,
      any_of(c(
        "q_BH_all_analytes","q_BH_polymers_only",
        "q_BH_within_family","q_BH_labor_prespec","fit_status"
      ))
    )
  write_csv(lab_primary, file.path(OUT_DIR, "labor_primary_IQR_FINAL.csv"))

  lab_legend <- paste0(
    "Supplementary Figure | Polymer-specific associations with spontaneous labor requiring augmentation. ",
    "Associations are shown as odds ratios (ORs) per analysis-specific interquartile-range (IQR) increase in ",
    "log10(1 + placental polymer burden [µg/g]), with 95% confidence intervals. The outcome compared ",
    "'Spontaneous: augmented' with 'Spontaneous' deliveries and excluded no-labor and missing/other labor categories. ",
    "Because augmented-labor events were limited, analyses were exploratory and intentionally parsimonious. ",
    "The primary model adjusted only for gestational age at delivery. Colored points show crude and one-covariate-at-a-time ",
    "sensitivity models adding hospital, Social Deprivation Index, prior preterm birth, or hypertension. ",
    "Point shape denotes nominal or Benjamini-Hochberg-adjusted statistical support as indicated. The x-axis is log scaled."
  )
  writeLines(lab_legend, file.path(OUT_DIR, "labor_figure_legend_IQR_FINAL.txt"))
}

message("Ordinal/labor IQR finalization complete.")
message("Tables/legends: ", OUT_DIR)
message("Figures: ", OUT_FIG_DIR)
