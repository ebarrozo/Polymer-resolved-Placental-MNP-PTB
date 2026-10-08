#!/usr/bin/env Rscript
# ============================================================
# 08_finalize_secondary_IQR.R
#
# Reviewer Point 1 — secondary-model effect-scale correction.
#
# Post-processes the latest existing outputs WITHOUT refitting:
#   A) ordinal gestational-age models (script 11)
#   B) spontaneous labor augmentation models (script 08)
#   C) histopathology models (script 09)
#   D) PU/PET continuous analytical-sensitivity models (script 18/19)
#
# For A-C, the fitted model used log10(1 + ug/g). The old reported
# "per doubling" ratios are re-expressed as ratios per IQR increase
# in the transformed exposure, with the IQR calculated in the exact
# model-specific complete-case analysis population.
#
# For D, the validation model used log2(1 + ug/g). The IQR effect is
# invariant to log base; the script uses the actual log2 exposure and
# its own analysis population, then combines corrected continuous IQR
# estimates with the unchanged binary threshold estimates.
#
# IMPORTANT
#   - p values and BH q values do not change.
#   - binary >=LOD / >LOQ ORs do not change.
#   - legacy output files are not overwritten.
#
# Run from PROJECT_ROOT:
#   source("scripts/08_finalize_secondary_IQR.R")
# ============================================================

source("00_config.R")

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
  library(stringr)
  library(ggplot2)
})

STAMP <- format(Sys.time(), "%Y%m%d_%H%M%S")
OUT_DIR <- file.path(DIR_TABLES, paste0("26_secondary_models_IQR__", STAMP))
OUT_FIG_DIR <- file.path(DIR_FIGS, "26_secondary_models_IQR")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_FIG_DIR, recursive = TRUE, showWarnings = FALSE)

K2_LOG10 <- log10(2)

find_latest_dir <- function(parent, prefix) {
  x <- list.dirs(parent, recursive = FALSE, full.names = TRUE)
  x <- x[grepl(paste0("^", prefix), basename(x))]
  if (!length(x)) return(NA_character_)
  x[order(file.info(x)$mtime, decreasing = TRUE)][1]
}

parse_covs <- function(x) {
  if (is.na(x) || !nzchar(trimws(x))) return(character())
  trimws(strsplit(x, "\\+")[[1]])
}

usable_covariates_local <- function(dat, covs) {
  keep <- character()
  for (cv in covs) {
    if (!cv %in% names(dat)) next
    x <- dat[[cv]]
    nonmiss <- !is.na(x)
    if (!sum(nonmiss)) next
    if (is.numeric(x) || is.integer(x)) {
      if (length(unique(x[nonmiss])) < 2 || sd(as.numeric(x), na.rm = TRUE) == 0) next
    } else {
      if (length(unique(as.character(x[nonmiss]))) < 2) next
    }
    keep <- c(keep, cv)
  }
  keep
}

recover_and_rescale <- function(or, lo, hi, iqr, unit = K2_LOG10) {
  if (!all(is.finite(c(or, lo, hi, iqr))) || or <= 0 || lo <= 0 || hi <= 0) {
    return(tibble(beta = NA_real_, se = NA_real_,
                  Ratio_IQR = NA_real_, Ratio_IQR_lo = NA_real_, Ratio_IQR_hi = NA_real_))
  }
  beta <- log(or) / unit
  beta_lo <- log(lo) / unit
  beta_hi <- log(hi) / unit
  se <- (beta_hi - beta_lo) / (2 * 1.96)
  tibble(
    beta = beta,
    se = se,
    Ratio_IQR = exp(beta * iqr),
    Ratio_IQR_lo = exp((beta - 1.96 * se) * iqr),
    Ratio_IQR_hi = exp((beta + 1.96 * se) * iqr)
  )
}

# ============================================================
# A. ORDINAL GESTATIONAL-AGE MODELS
# ============================================================

ord_dir <- find_latest_dir(DIR_TABLES, "11_ordinal_ga_category_models__")
if (!is.na(ord_dir)) {
  ord_all_path <- file.path(ord_dir, "ordinal_ga_category_all_models.csv")
  ord_dat_path <- file.path(ord_dir, "ordinal_ga_category_analysis_table.csv")

  if (file.exists(ord_all_path) && file.exists(ord_dat_path)) {
    ord <- read_csv(ord_all_path, show_col_types = FALSE)
    ddat <- read_csv(ord_dat_path, show_col_types = FALSE, guess_max = 10000)

    ord_rows <- vector("list", nrow(ord))

    for (i in seq_len(nrow(ord))) {
      r <- ord[i,]
      covs <- usable_covariates_local(ddat, parse_covs(r$Covariates))
      d <- ddat

      # Match the outlier-exclusion model logic if present.
      if (grepl("OUTLIER|NO_CLR|EXCLUD", r$ModelTag, ignore.case = TRUE) &&
          "cluster_clr_pca_outlier" %in% names(d)) {
        d <- d %>% filter(is.na(cluster_clr_pca_outlier) | !as.logical(cluster_clr_pca_outlier))
      }

      vars <- unique(c("GA_category_ordered", r$ExposureCol, covs))
      if (!all(vars %in% names(d))) {
        ord_rows[[i]] <- tibble(row_id=i, IQR_transformed=NA_real_, N_reconstructed=NA_integer_)
        next
      }

      dd <- d %>%
        select(all_of(vars)) %>%
        mutate(exposure = suppressWarnings(as.numeric(.data[[r$ExposureCol]]))) %>%
        drop_na() %>%
        filter(is.finite(exposure))

      iq <- IQR(dd$exposure, na.rm = TRUE)

      # Earlier-delivery OR
      z_e <- recover_and_rescale(
        r$OR_earlier_2x, r$OR_earlier_2x_lo, r$OR_earlier_2x_hi, iq
      )
      # Later-delivery OR
      z_l <- recover_and_rescale(
        r$OR_later_2x, r$OR_later_2x_lo, r$OR_later_2x_hi, iq
      )

      ord_rows[[i]] <- bind_cols(
        tibble(
          row_id=i,
          IQR_transformed=iq,
          N_reconstructed=nrow(dd),
          N_reported=r$N,
          N_match=nrow(dd) == r$N
        ),
        z_e %>% rename(
          beta_earlier=beta, se_earlier=se,
          OR_earlier_IQR=Ratio_IQR,
          OR_earlier_IQR_lo=Ratio_IQR_lo,
          OR_earlier_IQR_hi=Ratio_IQR_hi
        ),
        z_l %>% rename(
          beta_later=beta, se_later=se,
          OR_later_IQR=Ratio_IQR,
          OR_later_IQR_lo=Ratio_IQR_lo,
          OR_later_IQR_hi=Ratio_IQR_hi
        )
      )
    }

    ord_iqr <- ord %>%
      mutate(row_id=row_number()) %>%
      left_join(bind_rows(ord_rows), by="row_id") %>%
      select(-row_id)

    write_csv(ord_iqr, file.path(OUT_DIR, "ordinal_ga_category_all_models_IQR.csv"))
    write_csv(
      ord_iqr %>% filter(ModelTag == "M1_PRIMARY"),
      file.path(OUT_DIR, "ordinal_ga_category_primary_models_IQR.csv")
    )
    write_csv(
      ord_iqr %>% select(ModelTag, Polymer, N_reported=N, N_reconstructed, N_match),
      file.path(OUT_DIR, "ordinal_ga_N_reconstruction_check.csv")
    )
  }
}

# ============================================================
# B. LABOR AUGMENTATION MODELS
# ============================================================

labor_dir <- find_latest_dir(DIR_TABLES, "08_labor_models__")
if (!is.na(labor_dir)) {
  labor_all_path <- file.path(labor_dir, "labor_model_all_results.csv")
  model_table_path <- file.path(DIR_TABLES, "mnp_ptb_cohort__analysis_ready_model_table.csv")

  if (file.exists(labor_all_path) && file.exists(model_table_path)) {
    labor <- read_csv(labor_all_path, show_col_types = FALSE)
    ldat <- read_csv(model_table_path, show_col_types = FALSE, guess_max = 10000)

    # Reproduce labor-outcome derivation exactly.
    if ("Type_of_Labor_raw" %in% names(ldat)) {
      labor_raw <- as.character(ldat$Type_of_Labor_raw)
    } else if ("Type_of_Labor" %in% names(ldat)) {
      labor_raw <- as.character(ldat$Type_of_Labor)
    } else if ("Type of labor" %in% names(ldat)) {
      labor_raw <- as.character(ldat[["Type of labor"]])
    } else {
      labor_raw <- rep(NA_character_, nrow(ldat))
    }

    labor_norm <- labor_raw %>%
      str_to_lower() %>%
      str_replace_all("[^a-z0-9]+", " ") %>%
      str_squish()

    augmented <- str_detect(labor_norm, "^spontaneou?s augmented$|^spontaneou?s augment|augmented|augmentation")
    spontaneous_no_aug <- str_detect(labor_norm, "^spontaneou?s$|^spontaneous labor$") & !augmented

    ldat$spontaneous_labor_augmented01 <- case_when(
      augmented ~ 1L,
      spontaneous_no_aug ~ 0L,
      TRUE ~ NA_integer_
    )

    labor_rows <- vector("list", nrow(labor))

    for (i in seq_len(nrow(labor))) {
      r <- labor[i,]
      d0 <- ldat %>% filter(!is.na(spontaneous_labor_augmented01))
      covs <- usable_covariates_local(d0, parse_covs(r$Covariates))
      vars <- unique(c("spontaneous_labor_augmented01", r$ExposureCol, covs))

      if (!all(vars %in% names(d0))) {
        labor_rows[[i]] <- tibble(row_id=i, IQR_transformed=NA_real_, N_reconstructed=NA_integer_)
        next
      }

      dd <- d0 %>%
        select(all_of(vars)) %>%
        mutate(exposure = suppressWarnings(as.numeric(.data[[r$ExposureCol]]))) %>%
        drop_na() %>%
        filter(is.finite(exposure))

      iq <- IQR(dd$exposure, na.rm = TRUE)
      z <- recover_and_rescale(r$OR_2x, r$OR_2x_lo, r$OR_2x_hi, iq)

      labor_rows[[i]] <- bind_cols(
        tibble(
          row_id=i,
          IQR_transformed=iq,
          N_reconstructed=nrow(dd),
          N_reported=r$N,
          N_match=nrow(dd) == r$N
        ),
        z %>% rename(
          beta=beta, se=se,
          OR_IQR=Ratio_IQR,
          OR_IQR_lo=Ratio_IQR_lo,
          OR_IQR_hi=Ratio_IQR_hi
        )
      )
    }

    labor_iqr <- labor %>%
      mutate(row_id=row_number()) %>%
      left_join(bind_rows(labor_rows), by="row_id") %>%
      select(-row_id)

    write_csv(labor_iqr, file.path(OUT_DIR, "labor_model_all_results_IQR.csv"))

    # Detect primary model tag from label rather than hard-code if possible.
    primary_labor <- labor_iqr %>%
      filter(str_detect(str_to_lower(ModelLabel), "primary|ga")) %>%
      arrange(Polymer)

    write_csv(primary_labor, file.path(OUT_DIR, "labor_primary_candidate_IQR.csv"))
    write_csv(
      labor_iqr %>% select(ModelTag, ModelLabel, Polymer, N_reported=N, N_reconstructed, N_match),
      file.path(OUT_DIR, "labor_N_reconstruction_check.csv")
    )
  }
}

# ============================================================
# C. HISTOPATHOLOGY MODELS
# ============================================================

hist_dir <- find_latest_dir(DIR_TABLES, "09_histopathology_models__")
if (!is.na(hist_dir)) {
  hist_all_path <- file.path(hist_dir, "histopath_model_all_results.csv")
  hist_dat_path <- file.path(hist_dir, "histopath_joined_analysis_table.csv")

  if (file.exists(hist_all_path) && file.exists(hist_dat_path)) {
    hist <- read_csv(hist_all_path, show_col_types = FALSE)
    hdat <- read_csv(hist_dat_path, show_col_types = FALSE, guess_max = 10000)

    hist_rows <- vector("list", nrow(hist))

    for (i in seq_len(nrow(hist))) {
      r <- hist[i,]

      d0 <- hdat
      if ("has_path_scores" %in% names(d0)) d0 <- d0 %>% filter(as.logical(has_path_scores))
      if ("Delivery_class" %in% names(d0)) d0 <- d0 %>% filter(as.character(Delivery_class) == as.character(r$Delivery_class))

      if (r$Outcome == "microcalc_severity_positive") {
        if ("microcalc_positive" %in% names(d0)) d0 <- d0 %>% filter(as.logical(microcalc_positive))
      } else if (r$Outcome == "pvf_severity_positive") {
        if ("pvf_positive" %in% names(d0)) d0 <- d0 %>% filter(as.logical(pvf_positive))
      }

      covs <- usable_covariates_local(d0, parse_covs(r$Covariates))
      vars <- unique(c(r$OutcomeCol, r$ExposureCol, covs))

      if (!all(vars %in% names(d0))) {
        hist_rows[[i]] <- tibble(row_id=i, IQR_transformed=NA_real_, N_reconstructed=NA_integer_)
        next
      }

      dd <- d0 %>%
        select(all_of(vars)) %>%
        mutate(
          outcome = suppressWarnings(as.numeric(.data[[r$OutcomeCol]])),
          exposure = suppressWarnings(as.numeric(.data[[r$ExposureCol]]))
        ) %>%
        drop_na() %>%
        filter(is.finite(exposure), is.finite(outcome))

      if (r$ModelType == "Gamma") dd <- dd %>% filter(outcome > 0)

      iq <- IQR(dd$exposure, na.rm = TRUE)
      z <- recover_and_rescale(r$Ratio_2x, r$Ratio_2x_lo, r$Ratio_2x_hi, iq)

      hist_rows[[i]] <- bind_cols(
        tibble(
          row_id=i,
          IQR_transformed=iq,
          N_reconstructed=nrow(dd),
          N_reported=r$N,
          N_match=nrow(dd) == r$N
        ),
        z %>% rename(
          beta=beta, se=se,
          Ratio_IQR=Ratio_IQR,
          Ratio_IQR_lo=Ratio_IQR_lo,
          Ratio_IQR_hi=Ratio_IQR_hi
        )
      )
    }

    hist_iqr <- hist %>%
      mutate(row_id=row_number()) %>%
      left_join(bind_rows(hist_rows), by="row_id") %>%
      select(-row_id)

    write_csv(hist_iqr, file.path(OUT_DIR, "histopath_model_all_results_IQR.csv"))
    write_csv(
      hist_iqr %>% filter(ModelTag == "M1"),
      file.path(OUT_DIR, "histopath_primary_M1_IQR.csv")
    )
    write_csv(
      hist_iqr %>% select(Outcome, Delivery_class, ModelTag, Polymer, N_reported=N, N_reconstructed, N_match),
      file.path(OUT_DIR, "histopath_N_reconstruction_check.csv")
    )
  }
}

# ============================================================
# D. FIGURE 5 CONTINUOUS PU/PET ANALYTICAL SENSITIVITY
# ============================================================

val_dir <- find_latest_dir(DIR_TABLES, "18_placenta_spike_recovery_validation__")
if (!is.na(val_dir)) {
  cont_path <- file.path(val_dir, "pu_pet_historical_maxblank_subtracted_continuous_models.csv")
  audit_path <- file.path(val_dir, "pu_pet_historical_maxblank_subtraction_specimen_audit.csv")
  bin_path <- file.path(val_dir, "pu_pet_qc_binary_ptb_sensitivity_models.csv")
  firth_path <- file.path(val_dir, "pu_pet_qc_bias_reduced_models.csv")

  if (all(file.exists(c(cont_path, audit_path)))) {
    cont <- read_csv(cont_path, show_col_types = FALSE)
    aud <- read_csv(audit_path, show_col_types = FALSE)

    cont_iqr_rows <- vector("list", nrow(cont))

    for (i in seq_len(nrow(cont))) {
      r <- cont[i,]
      aa <- aud %>% filter(Polymer == r$Polymer)

      x <- if (r$Sensitivity == "Original_recomputed_from_analysis_ready_table") {
        log2(aa$original_burden_ug_g + 1)
      } else {
        log2(aa$blank_subtracted_ug_g + 1)
      }

      iq <- IQR(x[is.finite(x)], na.rm = TRUE)

      # Here one model unit = 1 log2 unit, so unit = 1.
      z <- recover_and_rescale(
        r$OR_per_doubling_plus1, r$CI_lo, r$CI_hi, iq, unit = 1
      )

      cont_iqr_rows[[i]] <- bind_cols(
        tibble(row_id=i, IQR_log2_1px=iq),
        z %>% rename(
          beta_log2=beta, se_log2=se,
          OR_IQR=Ratio_IQR, OR_IQR_lo=Ratio_IQR_lo, OR_IQR_hi=Ratio_IQR_hi
        )
      )
    }

    cont_iqr <- cont %>%
      mutate(row_id=row_number()) %>%
      left_join(bind_rows(cont_iqr_rows), by="row_id") %>%
      select(-row_id)

    write_csv(cont_iqr, file.path(OUT_DIR, "PU_PET_continuous_analytical_sensitivity_IQR.csv"))

    # Build corrected Figure 5 continuous+threshold panel.
    if (file.exists(bin_path)) {
      bin <- read_csv(bin_path, show_col_types = FALSE)
    } else bin <- tibble()

    if (file.exists(firth_path)) {
      firth <- read_csv(firth_path, show_col_types = FALSE)
    } else firth <- tibble()

    d_cont <- cont_iqr %>%
      filter(Polymer %in% c("PU","PET")) %>%
      transmute(
        Polymer,
        Model = case_when(
          Sensitivity == "Original_recomputed_from_analysis_ready_table" ~ "Primary continuous\n(per IQR)",
          Sensitivity == "MaxBlank_subtracted_per_run" ~ "Max-blank subtracted\n(per IQR)",
          TRUE ~ Sensitivity
        ),
        OR=OR_IQR, CI_lo=OR_IQR_lo, CI_hi=OR_IQR_hi,
        Family="Continuous"
      )

    d_bin <- if (nrow(bin)) {
      bin %>%
        filter(Polymer %in% c("PU","PET"), Threshold %in% c("at_or_above_LOD","above_LOQ")) %>%
        transmute(
          Polymer,
          Model=case_when(
            Threshold=="at_or_above_LOD" ~ ">=LOD vs <LOD",
            Threshold=="above_LOQ" ~ ">LOQ vs <=LOQ"
          ),
          OR, CI_lo, CI_hi, Family="Threshold"
        )
    } else tibble()

    d_firth <- if (nrow(firth)) {
      firth %>%
        filter(Polymer %in% c("PU","PET"), Threshold=="above_LOQ", is.finite(OR)) %>%
        transmute(
          Polymer, Model=">LOQ vs <=LOQ\n(Firth)",
          OR, CI_lo, CI_hi, Family="Bias-reduced"
        )
    } else tibble()

    dplot <- bind_rows(d_cont, d_bin, d_firth) %>%
      mutate(
        Polymer=factor(
          Polymer, levels=c("PU","PET"),
          labels=c("Polyurethane (PU)","Polyethylene terephthalate (PET)")
        ),
        Model=factor(
          Model,
          levels=rev(c(
            "Primary continuous\n(per IQR)",
            "Max-blank subtracted\n(per IQR)",
            ">=LOD vs <LOD",
            ">LOQ vs <=LOQ",
            ">LOQ vs <=LOQ\n(Firth)"
          ))
        )
      )

    p <- ggplot(dplot, aes(x=OR, y=Model, color=Family)) +
      geom_vline(xintercept=1, linetype="dashed", linewidth=0.4, color="grey40") +
      geom_errorbar(aes(xmin=CI_lo, xmax=CI_hi), orientation="y", width=0.15, linewidth=0.65) +
      geom_point(size=2.5) +
      facet_wrap(~Polymer, nrow=1) +
      scale_x_log10() +
      labs(
        title="Analytical sensitivity of PU and PET associations with preterm birth",
        subtitle="Continuous models are shown per analysis-specific IQR; threshold-model ORs are unchanged",
        x="Odds ratio (95% CI)",
        y=NULL,
        color=NULL
      ) +
      theme_bw(base_size=10.5) +
      theme(
        panel.grid.minor=element_blank(),
        plot.title=element_text(face="bold"),
        legend.position="bottom"
      )

    ggsave(
      file.path(OUT_FIG_DIR, "Figure5_analytical_sensitivity_IQR_FINAL.png"),
      p, width=8.0, height=4.8, dpi=300
    )
    ggsave(
      file.path(OUT_FIG_DIR, "Figure5_analytical_sensitivity_IQR_FINAL.pdf"),
      p, width=8.0, height=4.8
    )
  }
}

# ============================================================
# SUMMARY / QA
# ============================================================

qa_files <- list.files(OUT_DIR, pattern="N_reconstruction_check\\.csv$", full.names=TRUE)
qa <- map_dfr(qa_files, function(f) {
  x <- read_csv(f, show_col_types=FALSE)
  tibble(
    file=basename(f),
    n_rows=nrow(x),
    n_mismatch=sum(!x$N_match, na.rm=TRUE),
    n_missing_match=sum(is.na(x$N_match))
  )
})
write_csv(qa, file.path(OUT_DIR, "QA_N_reconstruction_summary.csv"))

summary_lines <- c(
  "Reviewer Point 1 secondary-model IQR audit",
  "==========================================",
  paste0("Project root: ", PROJECT_ROOT),
  paste0("Output directory: ", OUT_DIR),
  "",
  "Core rule:",
  "Continuous log-transformed exposure models are now expressed per analysis-specific IQR.",
  "Binary >=LOD and >LOQ threshold models are unchanged.",
  "",
  "QA reconstruction summary:",
  capture.output(print(qa, n=Inf)),
  "",
  "Next decision:",
  "Only replace secondary figures after all N_match checks are TRUE (or any mismatch is explained)."
)
writeLines(summary_lines, file.path(OUT_DIR, "README_secondary_IQR_results.txt"))

message(paste(summary_lines, collapse="\n"))
message("Candidate Figure 5 panel: ", OUT_FIG_DIR)
