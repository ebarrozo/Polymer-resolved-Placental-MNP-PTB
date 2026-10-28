#!/usr/bin/env Rscript
# ============================================================
# 29_finalize_histopathology_IQR.R
#
# Finalize reviewer Point 1 for histopathology models.
#
# Robust behavior:
#   - searches all 09_histopathology_models__* directories for a
#     complete set of canonical model + joined-analysis files;
#   - if none exists, runs scripts/09_histopathology_models.R once;
#   - reconstructs each model's complete-case population;
#   - calculates primary-model IQRs within each
#     Outcome x Delivery_class x Polymer stratum;
#   - uses that fixed primary-model IQR for all sensitivity models
#     within the same stratum so coefficient stability is compared
#     on a common reporting unit;
#   - regenerates microcalcification and PVF figures.
#
# Does NOT overwrite legacy outputs.
# ============================================================

source("00_config.R")

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
  library(stringr)
  library(ggplot2)
})

PRIMARY_MODEL_TAG <- "M1_CLIN"
K2 <- log10(2)

find_complete_hist_dir <- function() {
  dirs <- list.dirs(DIR_TABLES, recursive=FALSE, full.names=TRUE)
  dirs <- dirs[grepl("^09_histopathology_models__", basename(dirs))]
  req <- c("histopath_model_all_results.csv","histopath_joined_analysis_table.csv")
  ok <- vapply(dirs, function(d) all(file.exists(file.path(d, req))), logical(1))
  dirs <- dirs[ok]
  if (!length(dirs)) return(NA_character_)
  dirs[order(file.info(dirs)$mtime, decreasing=TRUE)][1]
}

HIST_DIR <- find_complete_hist_dir()
if (is.na(HIST_DIR)) {
  canonical <- file.path("scripts","09_histopathology_models.R")
  if (!file.exists(canonical)) stop("Missing canonical histopathology script.", call. = FALSE)
  message("No complete histopathology output directory found. Running ", canonical)
  source(canonical)
  HIST_DIR <- find_complete_hist_dir()
}
if (is.na(HIST_DIR)) stop("Unable to generate complete histopathology outputs.", call. = FALSE)

message("Using histopathology directory: ", HIST_DIR)

OUT_DIR <- file.path(DIR_TABLES, "29_histopathology_IQR_final")
OUT_FIG_DIR <- file.path(DIR_FIGS, "29_histopathology_IQR_final")
dir.create(OUT_DIR, recursive=TRUE, showWarnings=FALSE)
dir.create(OUT_FIG_DIR, recursive=TRUE, showWarnings=FALSE)

res <- read_csv(file.path(HIST_DIR,"histopath_model_all_results.csv"), show_col_types=FALSE)
dat <- read_csv(file.path(HIST_DIR,"histopath_joined_analysis_table.csv"), show_col_types=FALSE, guess_max=10000)

parse_covs <- function(x) {
  if (is.na(x) || !nzchar(trimws(x))) return(character())
  trimws(strsplit(x, "\\+")[[1]])
}

usable_covariates <- function(dat, covs) {
  keep <- character()
  for (cv in covs) {
    if (!cv %in% names(dat)) next
    x <- dat[[cv]]
    nonmiss <- !is.na(x)
    if (!sum(nonmiss)) next
    if (is.numeric(x) || is.integer(x)) {
      if (length(unique(x[nonmiss])) < 2 || sd(as.numeric(x), na.rm=TRUE)==0) next
    } else {
      if (length(unique(as.character(x[nonmiss]))) < 2) next
    }
    keep <- c(keep, cv)
  }
  keep
}

reconstruct_row <- function(r) {
  d <- dat
  if ("has_path_scores" %in% names(d)) d <- d %>% filter(as.logical(has_path_scores))
  if ("Delivery_class" %in% names(d)) d <- d %>% filter(as.character(Delivery_class)==as.character(r$Delivery_class))

  if (r$Outcome=="microcalc_severity_positive" && "microcalc_positive" %in% names(d)) {
    d <- d %>% filter(as.logical(microcalc_positive))
  }
  if (r$Outcome=="pvf_severity_positive" && "pvf_positive" %in% names(d)) {
    d <- d %>% filter(as.logical(pvf_positive))
  }

  covs <- usable_covariates(d, parse_covs(r$Covariates))
  vars <- unique(c(r$OutcomeCol, r$ExposureCol, covs))
  if (!all(vars %in% names(d))) {
    return(tibble(N_reconstructed=NA_integer_, IQR_row=NA_real_))
  }

  dd <- d %>%
    select(all_of(vars)) %>%
    mutate(
      outcome=suppressWarnings(as.numeric(.data[[r$OutcomeCol]])),
      exposure=suppressWarnings(as.numeric(.data[[r$ExposureCol]]))
    ) %>%
    drop_na() %>%
    filter(is.finite(outcome),is.finite(exposure))

  if (r$ModelType=="Gamma") dd <- dd %>% filter(outcome>0)

  tibble(
    N_reconstructed=nrow(dd),
    IQR_row=IQR(dd$exposure, na.rm=TRUE)
  )
}

recon <- map_dfr(seq_len(nrow(res)), function(i) {
  bind_cols(tibble(row_id=i), reconstruct_row(res[i,]))
})

res2 <- res %>%
  mutate(row_id=row_number()) %>%
  left_join(recon, by="row_id") %>%
  mutate(N_match = N == N_reconstructed)

write_csv(
  res2 %>% select(Outcome,Delivery_class,ModelTag,Polymer,N,N_reconstructed,N_match,IQR_row),
  file.path(OUT_DIR,"histopath_N_IQR_reconstruction_check.csv")
)

if (any(res2$N_match==FALSE, na.rm=TRUE)) {
  stop("Histopathology N reconstruction mismatch detected; inspect QA table.", call. = FALSE)
}

# Primary-model IQR fixed within Outcome x Delivery x Polymer.
iqr_primary <- res2 %>%
  filter(ModelTag==PRIMARY_MODEL_TAG) %>%
  select(Outcome,Delivery_class,Polymer,IQR_primary=IQR_row)

res_iqr <- res2 %>%
  left_join(iqr_primary, by=c("Outcome","Delivery_class","Polymer")) %>%
  mutate(
    beta=log(Ratio_2x)/K2,
    beta_lo=log(Ratio_2x_lo)/K2,
    beta_hi=log(Ratio_2x_hi)/K2,
    se=(beta_hi-beta_lo)/(2*1.96),
    Ratio_IQR=exp(beta*IQR_primary),
    Ratio_IQR_lo=exp((beta-1.96*se)*IQR_primary),
    Ratio_IQR_hi=exp((beta+1.96*se)*IQR_primary)
  )

write_csv(res_iqr, file.path(OUT_DIR,"histopath_model_all_results_IQR.csv"))
write_csv(
  res_iqr %>% filter(ModelTag==PRIMARY_MODEL_TAG),
  file.path(OUT_DIR,"histopath_primary_M1_IQR.csv")
)

theme_nc <- theme_bw(base_size=10.3) +
  theme(
    panel.grid.minor=element_blank(),
    plot.title=element_text(face="bold"),
    strip.background=element_rect(fill="grey92", color="grey55"),
    strip.text=element_text(face="bold"),
    legend.position="bottom"
  )

model_colors <- c(
  "CRUDE"="#56B4E9",
  "M1 no SDI"="#E69F00",
  "M1 + SVI"="#009E73",
  "M1 + race"="#CC79A7"
)

make_plot <- function(outcome_name, ratio_name, title, stem) {
  d0 <- res_iqr %>% filter(Outcome==outcome_name, is.finite(Ratio_IQR))
  if (!nrow(d0)) return(invisible(NULL))

  # nicer labels
  d0 <- d0 %>%
    mutate(
      PolymerLabel=case_when(
        Polymer=="PET" ~ "Polyethylene terephthalate (PET)",
        Polymer=="PU" ~ "Polyurethane (PU)",
        Polymer=="PE" ~ "Polyethylene (PE)",
        Polymer=="PP" ~ "Polypropylene (PP)",
        Polymer=="PVC" ~ "Polyvinyl chloride (PVC)",
        Polymer=="PMMA" ~ "Polymethyl methacrylate (PMMA)",
        Polymer=="PC" ~ "Polycarbonate (PC)",
        Polymer=="SBR" ~ "Styrene-butadiene rubber (SBR)",
        Polymer=="ABS" ~ "Acrylonitrile-butadiene-styrene (ABS)",
        Polymer=="N6" ~ "Nylon-6 (N6)",
        Polymer=="N66" ~ "Nylon-66 (N66)",
        TRUE ~ Polymer
      )
    )

  primary <- d0 %>% filter(ModelTag==PRIMARY_MODEL_TAG)
  sens <- d0 %>% filter(ModelTag!=PRIMARY_MODEL_TAG)

  # Order by preterm primary estimate if possible, otherwise all primary.
  ord <- primary %>%
    arrange(Delivery_class!="Preterm", desc(Ratio_IQR)) %>%
    pull(PolymerLabel) %>% unique()

  primary <- primary %>% mutate(PolymerLabel=factor(PolymerLabel, levels=rev(ord)))
  sens <- sens %>% mutate(PolymerLabel=factor(PolymerLabel, levels=rev(ord)))

  primary <- primary %>%
    mutate(
      evidence=case_when(
        !is.na(q_global) & q_global<0.05 ~ "global BH q<0.05",
        !is.na(q_family) & q_family<0.05 ~ "family BH q<0.05",
        !is.na(p) & p<0.05 ~ "nominal p<0.05",
        TRUE ~ "not nominal"
      )
    )

  p <- ggplot() +
    geom_vline(xintercept=1, linetype="dashed", linewidth=0.35, color="grey45") +
    geom_errorbar(
      data=primary,
      aes(xmin=Ratio_IQR_lo,xmax=Ratio_IQR_hi,y=PolymerLabel),
      orientation="y", width=0.16, linewidth=0.7, color="black"
    ) +
    geom_point(
      data=primary,
      aes(x=Ratio_IQR,y=PolymerLabel,shape=evidence),
      size=2.6, color="black"
    ) +
    geom_point(
      data=sens,
      aes(x=Ratio_IQR,y=PolymerLabel,color=ModelLabel),
      position=position_dodge(width=0.45),
      size=1.7
    ) +
    facet_wrap(~Delivery_class, nrow=1, scales="free_y") +
    scale_x_log10() +
    scale_color_manual(values=model_colors, na.translate=FALSE) +
    scale_shape_manual(values=c(
      "not nominal"=1,
      "nominal p<0.05"=16,
      "family BH q<0.05"=22,
      "global BH q<0.05"=15
    )) +
    labs(
      title=title,
      subtitle="Black points and 95% CIs show the primary clinical model; colored points show sensitivity models",
      x=paste0(ratio_name," per stratum-specific IQR increase in log10(1 + burden)"),
      y=NULL,
      color="Sensitivity model",
      shape="Primary-model evidence"
    ) +
    theme_nc

  ggsave(file.path(OUT_FIG_DIR,paste0(stem,".png")),p,width=10.2,height=6.3,dpi=300)
  ggsave(file.path(OUT_FIG_DIR,paste0(stem,".pdf")),p,width=10.2,height=6.3)
}

make_plot(
  "microcalc_severity_positive",
  "IRR",
  "Microcalcification severity among microcalcification-positive placentae",
  "histopath_microcalcification_IQR_FINAL"
)

make_plot(
  "pvf_severity_positive",
  "RR",
  "Perivillous fibrin deposition severity among positive placentae",
  "histopath_PVF_IQR_FINAL"
)

legend_txt <- paste0(
  "Supplementary Figure | Polymer-specific associations with placental histopathology severity. ",
  "Microcalcification severity was modeled among microcalcification-positive samples using negative binomial regression of the summed duplicate microcalcification score, ",
  "and perivillous fibrin deposition (PVF) severity was modeled among PVF-positive samples using gamma regression with a log link applied to the average PVF score. ",
  "Models were stratified by delivery class. Estimates are shown as incidence rate ratios (IRRs) for microcalcification or relative ratios (RRs) for PVF per ",
  "stratum-specific interquartile-range (IQR) increase in log10(1 + placental polymer burden [µg/g]). ",
  "Within each outcome, delivery class, and polymer, the IQR from the primary clinical model analysis population was used for all sensitivity-model estimates. ",
  "Black points and 95% confidence intervals represent the primary clinical model; colored points represent crude and sensitivity models. ",
  "Benjamini-Hochberg correction was calculated globally across modeled polymers and within polymer families for each outcome, delivery class, and model. ",
  "The x-axis is log scaled."
)
writeLines(legend_txt, file.path(OUT_DIR,"histopathology_figure_legend_IQR_FINAL.txt"))

message("Histopathology IQR finalization complete.")
message("Tables/legend: ",OUT_DIR)
message("Figures: ",OUT_FIG_DIR)
