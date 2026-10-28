#!/usr/bin/env Rscript
# POLISHED VERSION: sources R/figure_theme_mnp.R and writes *_polished outputs where applicable.
# ============================================================
# 08b_labor_figure_polish.R
# v3.4: labor figure polish with data-derived full-CI log axes
# ============================================================
# Run after labor models have been generated:
#   source("run_labor_figure_polish.R")
#
# Key fix:
#   Axis limits are derived from the actual finite OR_2x_lo/OR_2x_hi
#   values in labor_model_all_results.csv / split model result tables.
#   This prevents misleading clipping from fixed 0.03-10 limits.
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
  library(stringr)
  library(ggplot2)
  library(scales)
})

source("00_config.R")
source(file.path("R", "utils_mnp.R"))
source(file.path("R", "figure_theme_mnp.R"))

BASE_SIZE <- 11
OVERWRITE_STANDARD_NAMES <- TRUE

# ----------------------------
# Helpers
# ----------------------------
latest_labor_dir <- function() {
  dirs <- list.dirs(file.path(DIR_TABLES), recursive = FALSE, full.names = TRUE)
  hits <- dirs[grepl("^08_labor_models__", basename(dirs))]
  if (length(hits) == 0) {
    stop("No 08_labor_models__ output directory found. Run run_labor_outcome_supplement.R first.", call. = FALSE)
  }
  hits[order(file.info(hits)$mtime, decreasing = TRUE)][1]
}

poly_label <- function(poly) {
  if (poly == "Cumulative_MNPs") return("Cumulative MNP burden")
  if (poly %in% names(POLYMER_FULL_NAMES)) return(paste0(POLYMER_FULL_NAMES[[poly]], " (", poly, ")"))
  poly
}

read_model_table_local <- function(path) {
  readr::read_csv(path, show_col_types = FALSE, guess_max = 10000) %>%
    mutate(
      Delivery_class = factor(as.character(Delivery_class), levels = c("Term", "Preterm")),
      Type_of_Labor = dplyr::coalesce(safe_chr(Type_of_Labor), safe_chr(`Type of labor`))
    )
}

derive_labor_counts <- function(df) {
  labor_raw <- safe_chr(df$Type_of_Labor)
  labor_norm <- labor_raw %>%
    str_to_lower() %>%
    str_replace_all("[^a-z0-9]+", " ") %>%
    str_squish()

  augmented <- str_detect(labor_norm, "^spontaneou?s augmented$|^spontaneou?s augment|augmented|augmentation")
  spontaneous_no_aug <- str_detect(labor_norm, "^spontaneou?s$|^spontaneous labor$") & !augmented
  no_labor <- str_detect(labor_norm, "^no labor$|^none$|^no$")

  df %>%
    mutate(
      labor_group = case_when(
        spontaneous_no_aug ~ "Spontaneous",
        augmented ~ "Spontaneous: augmented",
        no_labor ~ "No labor",
        TRUE ~ "Missing/other"
      ),
      labor_group = factor(
        labor_group,
        levels = c("Spontaneous", "Spontaneous: augmented", "No labor", "Missing/other")
      )
    ) %>%
    filter(!is.na(Delivery_class)) %>%
    count(Delivery_class, labor_group, name = "n")
}

# Derive full display limits from actual finite CIs/ORs.
# Rounds outward to powers of ten, with no arbitrary truncation.
compute_full_ci_axis <- function(..., min_floor = 1e-12, max_ceiling = 1e12) {
  dfs <- list(...)
  vals <- purrr::map(dfs, function(x) {
    if (is.null(x) || nrow(x) == 0) return(numeric(0))
    x %>%
      select(any_of(c("OR_2x_lo", "OR_2x", "OR_2x_hi"))) %>%
      unlist(use.names = FALSE) %>%
      as.numeric()
  }) %>% unlist(use.names = FALSE)

  vals <- vals[is.finite(vals) & vals > 0]
  if (length(vals) == 0) return(list(limits = c(0.03, 10), breaks = c(0.03, 0.1, 0.5, 1, 2, 10)))

  lo <- max(min_floor, 10^floor(log10(min(vals))))
  hi <- min(max_ceiling, 10^ceiling(log10(max(vals))))

  # Powers of ten are clearest when the range is very wide.
  pow_seq <- seq(floor(log10(lo)), ceiling(log10(hi)), by = 1)
  breaks <- 10^pow_seq

  list(limits = c(lo, hi), breaks = breaks)
}

label_log_axis <- function(x) {
  out <- ifelse(
    x >= 1000 | x < 0.01,
    scales::scientific(x, digits = 1),
    scales::number(x, accuracy = 0.01, trim = TRUE)
  )
  out
}

prep_plot_dat <- function(tbl, order_levels = NULL) {
  dat <- tbl %>%
    filter(!is.na(OR_2x)) %>%
    mutate(
      PolymerLabel = if ("PolymerLabel" %in% names(.)) PolymerLabel else vapply(Polymer, poly_label, character(1)),
      evidence = case_when(
        "q_BH_all_analytes" %in% names(.) & !is.na(q_BH_all_analytes) & q_BH_all_analytes < 0.05 ~ "BH q < 0.05",
        !is.na(p) & p < 0.05 ~ "nominal p < 0.05 only",
        TRUE ~ "not nominal"
      ),
      evidence = factor(evidence, levels = c("not nominal", "nominal p < 0.05 only", "BH q < 0.05"))
    )

  if (is.null(order_levels)) {
    order_levels <- dat %>% arrange(desc(OR_2x)) %>% pull(PolymerLabel)
  }

  dat %>% mutate(PolymerLabel = factor(PolymerLabel, levels = rev(order_levels)))
}

make_primary_fig <- function(primary_tbl, out_path) {
  axis <- compute_full_ci_axis(primary_tbl)
  plot_dat <- prep_plot_dat(primary_tbl)

  subtitle_txt <- str_wrap(
    paste0(
      "Primary model: gestational age + maternal age + BMI + SDI + twins + prior PTB + hypertension + hospital. ",
      "Black points and 95% CI bars = primary adjusted model. Axis spans the full finite CI range in the model table."
    ),
    width = 100
  )

  p <- ggplot(plot_dat, aes(y = PolymerLabel, x = OR_2x)) +
    geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.35, color = "grey30") +
    geom_segment(aes(x = OR_2x_lo, xend = OR_2x_hi, yend = PolymerLabel), linewidth = 0.85, color = "black") +
    geom_point(aes(shape = evidence, fill = evidence), size = 3.1, color = "black", stroke = 0.9) +
    scale_x_log10(
      limits = axis$limits,
      breaks = axis$breaks,
      labels = label_log_axis,
      expand = expansion(mult = c(0.01, 0.01))
    ) +
    scale_shape_manual(
      values = c("not nominal" = 21, "nominal p < 0.05 only" = 21, "BH q < 0.05" = 24),
      name = "Primary model\nevidence"
    ) +
    scale_fill_manual(
      values = c("not nominal" = "black", "nominal p < 0.05 only" = "white", "BH q < 0.05" = "black"),
      guide = "none"
    ) +
    labs(
      title = "Polymer-specific associations with spontaneous labor requiring augmentation",
      subtitle = subtitle_txt,
      x = "OR per doubling (95% CI)",
      y = NULL
    ) +
    theme_mnp_forest(base_size = BASE_SIZE) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(size = 10.2),
      legend.position = "right",
      legend.title = element_text(face = "bold")
    )

  ggsave_mnp(out_path, p, width = 10.5, height = 6.6, dpi = 300)
}

make_sensitivity_fig <- function(primary_tbl, noga_tbl, term_tbl, preterm_tbl, out_path) {
  axis <- compute_full_ci_axis(primary_tbl, noga_tbl, term_tbl, preterm_tbl)

  order_levels <- prep_plot_dat(primary_tbl) %>%
    arrange(desc(OR_2x)) %>%
    pull(PolymerLabel) %>%
    as.character()

  primary_dat <- prep_plot_dat(primary_tbl, order_levels)

  sens_dat <- bind_rows(
    noga_tbl %>% mutate(Model = "Overall without GA"),
    term_tbl %>% mutate(Model = "Term exploratory"),
    preterm_tbl %>% mutate(Model = "Preterm exploratory")
  ) %>%
    filter(!is.na(OR_2x))

  if (nrow(sens_dat) > 0) {
    sens_dat <- prep_plot_dat(sens_dat, order_levels) %>%
      mutate(Model = factor(Model, levels = c("Overall without GA", "Term exploratory", "Preterm exploratory")))
  }

  subtitle_txt <- str_wrap(
    "Black points and 95% CI bars = overall GA-adjusted primary model. Colored points and bars = sensitivity models. Axis spans the full finite CI range in the model table.",
    width = 100
  )

  p <- ggplot() +
    geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.35, color = "grey30") +
    geom_segment(
      data = primary_dat,
      aes(y = PolymerLabel, yend = PolymerLabel, x = OR_2x_lo, xend = OR_2x_hi),
      linewidth = 0.8,
      color = "black"
    ) +
    geom_point(data = primary_dat, aes(y = PolymerLabel, x = OR_2x), size = 3.0, shape = 16, color = "black")

  if (nrow(sens_dat) > 0) {
    p <- p +
      geom_segment(
        data = sens_dat,
        aes(y = PolymerLabel, yend = PolymerLabel, x = OR_2x_lo, xend = OR_2x_hi, color = Model),
        linewidth = 0.55,
        alpha = 0.55
      ) +
      geom_point(data = sens_dat, aes(y = PolymerLabel, x = OR_2x, color = Model, shape = Model), size = 2.6, stroke = 0.8)
  }

  p <- p +
    scale_x_log10(
      limits = axis$limits,
      breaks = axis$breaks,
      labels = label_log_axis,
      expand = expansion(mult = c(0.01, 0.01))
    ) +
    labs(
      title = "Labor outcome model sensitivity analysis",
      subtitle = subtitle_txt,
      x = "OR per doubling (95% CI)",
      y = NULL,
      color = "Sensitivity models",
      shape = "Sensitivity models"
    ) +
    theme_mnp_forest(base_size = BASE_SIZE) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(size = 10.2),
      legend.position = "right",
      legend.title = element_text(face = "bold")
    )

  if (nrow(sens_dat) > 0) {
    present_levels <- unique(as.character(na.omit(sens_dat$Model)))
    color_vals <- c("Overall without GA" = "#56B4E9", "Term exploratory" = "#E69F00", "Preterm exploratory" = "#009E73")
    shape_vals <- c("Overall without GA" = 16, "Term exploratory" = 17, "Preterm exploratory" = 15)
    p <- p +
      scale_color_manual(values = color_vals[present_levels], breaks = present_levels) +
      scale_shape_manual(values = shape_vals[present_levels], breaks = present_levels)
  }

  ggsave_mnp(out_path, p, width = 10.5, height = 6.6, dpi = 300)
}

make_counts_fig <- function(counts_tbl, out_path) {
  totals <- counts_tbl %>%
    group_by(Delivery_class) %>%
    summarise(n_total = sum(n), .groups = "drop") %>%
    mutate(lbl = paste0(as.character(Delivery_class), " (n=", n_total, ")"))

  counts_plot <- counts_tbl %>%
    left_join(totals, by = "Delivery_class") %>%
    mutate(
      Delivery_label = factor(lbl, levels = totals$lbl),
      labor_group = factor(labor_group, levels = c("Spontaneous", "Spontaneous: augmented", "No labor", "Missing/other"))
    )

  p <- ggplot(counts_plot, aes(x = Delivery_label, y = n, fill = labor_group)) +
    geom_col(width = 0.72, color = "black", linewidth = 0.25) +
    geom_text(aes(label = ifelse(n > 0, n, "")), position = position_stack(vjust = 0.5), size = 3.2) +
    scale_fill_manual(values = c(
      "Spontaneous" = "#F8766D",
      "Spontaneous: augmented" = "#7CAE00",
      "No labor" = "#00BFC4",
      "Missing/other" = "#C77CFF"
    )) +
    labs(
      title = "Type of labor by delivery class",
      subtitle = "Outcome definition for labor models: Spontaneous: augmented vs Spontaneous",
      x = NULL,
      y = "Samples",
      fill = "Type of labor"
    ) +
    theme_mnp_forest(base_size = BASE_SIZE) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold"),
      legend.position = "right",
      legend.title = element_text(face = "bold")
    )

  ggsave_mnp(out_path, p, width = 8.2, height = 5.9, dpi = 300)
}

# ----------------------------
# Run
# ----------------------------
latest <- latest_labor_dir()
message("Using labor output directory: ", latest)

primary_tbl <- readr::read_csv(file.path(latest, "labor_model_primary_overall_GA_adjusted_results_polished.csv"), show_col_types = FALSE)
term_tbl <- readr::read_csv(file.path(latest, "labor_model_sensitivity_term_results_polished.csv"), show_col_types = FALSE)
preterm_tbl <- readr::read_csv(file.path(latest, "labor_model_sensitivity_preterm_results_polished.csv"), show_col_types = FALSE)
noga_tbl <- readr::read_csv(file.path(latest, "labor_model_sensitivity_overall_no_GA_results_polished.csv"), show_col_types = FALSE)
model_tbl <- read_model_table_local(file.path(DIR_TABLES, "mnp_ptb_cohort__analysis_ready_model_table.csv"))
counts_tbl <- derive_labor_counts(model_tbl)

primary_full <- file.path(DIR_FIGS, "supp_labor_models__overall_GA_adjusted_primary_forest_full_ci_axis_polished.png")
sens_full <- file.path(DIR_FIGS, "supp_labor_models__overall_primary_plus_stratified_sensitivities_full_ci_axis_polished.png")
counts_out <- file.path(DIR_FIGS, "supp_labor_models__type_of_labor_counts_polished.png")

make_primary_fig(primary_tbl, primary_full)
make_sensitivity_fig(primary_tbl, noga_tbl, term_tbl, preterm_tbl, sens_full)
make_counts_fig(counts_tbl, counts_out)

if (isTRUE(OVERWRITE_STANDARD_NAMES)) {
  file.copy(primary_full, file.path(DIR_FIGS, "supp_labor_models__overall_GA_adjusted_primary_forest_polished.png"), overwrite = TRUE)
  file.copy(sens_full, file.path(DIR_FIGS, "supp_labor_models__overall_primary_plus_stratified_sensitivities_polished.png"), overwrite = TRUE)
}

message("Wrote full-CI-axis labor figures to outputs/figs/.")
