# =============================================================================
# 32_ptb_heterogeneity_sensitivity.R
#
# PTB heterogeneity sensitivity analyses.
#
# Public reproducibility version:
#   * Uses the locked PTB outcome for all inferential models.
#   * Uses continuous GA only as an audit of the locked PTB definition.
#   * Retains subjects with clear locked PTB status even if continuous GA cannot
#     be reconstructed.
#   * Uses fixed polymer-specific IQRs from the full locked primary cohort.
#   * Uses HC1 robust standard errors.
#   * Includes:
#       1) singleton-only PTB
#       2) HDP-excluded PTB
#       3) spontaneous PTB vs term
#       4) indicated PTB vs term
#
# Output:
#   outputs/tables/32_ptb_heterogeneity_sensitivity/
#
# Public-repository note:
#   Only aggregate model and classification-audit outputs are written.
#   No subject or specimen identifiers are exported by this script.
#
# Notes:
#   - Spontaneous PTB = spontaneous labor PTB + PPROM leading to spontaneous PTB.
#   - Indicated PTB = indicated PTB + PPROM leading to induction/cesarean.
#   - Locked PTBs carrying a "Term delivery" subtype label are excluded from
#     subtype-specific models because subtype is unresolved.
# =============================================================================

req_pkgs <- c(
  "dplyr", "readr", "tidyr", "stringr", "purrr", "tibble",
  "sandwich", "lmtest"
)

missing_pkgs <- req_pkgs[
  !vapply(req_pkgs, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_pkgs) > 0) {
  stop(
    "Missing required packages: ",
    paste(missing_pkgs, collapse = ", "),
    "\nInstall them before running this script.",
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(tibble)
  library(sandwich)
  library(lmtest)
})

# -----------------------------------------------------------------------------
# 0. Project configuration and paths
# -----------------------------------------------------------------------------

# Run from PROJECT_ROOT, consistent with the public analysis workflow.
source("00_config.R")

IN_MODEL_TABLE <- file.path(
  DIR_TABLES,
  "mnp_ptb_cohort__analysis_ready_model_table.csv"
)

OUT_DIR <- file.path(
  DIR_TABLES,
  "32_ptb_heterogeneity_sensitivity"
)

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(IN_MODEL_TABLE)) {
  stop(
    "Required analysis-ready model table not found: ",
    IN_MODEL_TABLE,
    "\nGenerate the analysis-ready inputs before running script 32.",
    call. = FALSE
  )
}

message("Reading: ", IN_MODEL_TABLE)
df <- readr::read_csv(IN_MODEL_TABLE, show_col_types = FALSE)

# -----------------------------------------------------------------------------
# 1. Helper functions
# -----------------------------------------------------------------------------

first_existing <- function(dat, candidates, required = TRUE) {
  hit <- candidates[candidates %in% names(dat)]
  if (length(hit) == 0) {
    if (required) {
      stop(
        "None of the expected columns were found: ",
        paste(candidates, collapse = ", ")
      )
    }
    return(NA_character_)
  }
  hit[[1]]
}

as_binary01 <- function(x) {
  if (is.logical(x)) return(ifelse(is.na(x), NA_integer_, as.integer(x)))

  x_chr <- trimws(tolower(as.character(x)))
  out <- rep(NA_integer_, length(x_chr))

  out[x_chr %in% c("1", "true", "yes", "y", "preterm", "ptb")] <- 1L
  out[x_chr %in% c("0", "false", "no", "n", "term")] <- 0L

  # Numeric-looking values not already mapped
  suppressWarnings(x_num <- as.numeric(x_chr))
  idx <- is.na(out) & !is.na(x_num) & x_num %in% c(0, 1)
  out[idx] <- as.integer(x_num[idx])

  out
}

safe_bh <- function(p) {
  out <- rep(NA_real_, length(p))
  ok <- is.finite(p)
  if (any(ok)) out[ok] <- p.adjust(p[ok], method = "BH")
  out
}

fmt_range <- function(x) {
  if (length(x) == 0 || all(is.na(x))) return(NA_character_)
  sprintf("%.3f to %.3f", min(x, na.rm = TRUE), max(x, na.rm = TRUE))
}

# -----------------------------------------------------------------------------
# 2. Resolve core variables
# -----------------------------------------------------------------------------

locked_ptb_col <- first_existing(
  df,
  c("Preterm", "preterm_birth", "PTB", "ptb")
)

weeks_col <- first_existing(
  df,
  c("ga_delivery_weeks", "GA_weeks"),
  required = FALSE
)

days_remainder_col <- first_existing(
  df,
  c("ga_delivery_days_remainder", "GA_days_remainder"),
  required = FALSE
)

total_days_col <- first_existing(
  df,
  c("ga_delivery_days", "GA_days"),
  required = FALSE
)

twins_col <- first_existing(
  df,
  c("Twins01", "twin_pregnancy", "Twins")
)

delivery_class_col <- first_existing(
  df,
  c("Delivery_classification", "delivery_classification", "Delivery classification")
)

# Primary covariates
age_col <- first_existing(df, c("Maternal_age", "maternal_age"))
bmi_col <- first_existing(df, c("BMI", "bmi"))
sdi_col <- first_existing(df, c("SDI_score", "sdi_score", "SDI"))
prior_ptb_col <- first_existing(df, c("Prior_PTB", "prior_ptb"))
htn01_col <- first_existing(df, c("Hypertension01", "hypertension01"))
hospital_col <- first_existing(df, c("Hospital", "hospital"))

# HDP variables
preE_col <- first_existing(
  df,
  c("Preeclampsia01", "preeclampsia01"),
  required = FALSE
)

gh_col <- first_existing(
  df,
  c("Gestational_Hypertension01", "gestational_hypertension01"),
  required = FALSE
)

htn_type_col <- first_existing(
  df,
  c("Hypertension_type", "Preeclampsia/gestational hypertension"),
  required = FALSE
)

# -----------------------------------------------------------------------------
# 3. Locked PTB outcome + continuous-GA audit
# -----------------------------------------------------------------------------

df <- df %>%
  mutate(
    PTB_locked = as_binary01(.data[[locked_ptb_col]]),
    Twins_locked = as_binary01(.data[[twins_col]])
  )

# Build continuous GA only for audit.
# Prefer completed weeks + day remainder/7.
df$GA_continuous_audit <- NA_real_
ga_source <- "unavailable"

if (!is.na(weeks_col) && !is.na(days_remainder_col)) {
  wk <- suppressWarnings(as.numeric(df[[weeks_col]]))
  dr <- suppressWarnings(as.numeric(df[[days_remainder_col]]))

  # Only treat day-remainder values 0-6 as valid remainder values.
  ok <- is.finite(wk) & is.finite(dr) & dr >= 0 & dr <= 6
  df$GA_continuous_audit[ok] <- wk[ok] + dr[ok] / 7
  ga_source <- paste0(weeks_col, " + ", days_remainder_col, "/7")
}

# Fallback to total gestational days / 7 only where still missing.
if (!is.na(total_days_col)) {
  td <- suppressWarnings(as.numeric(df[[total_days_col]]))
  ok <- is.na(df$GA_continuous_audit) & is.finite(td) & td > 100
  df$GA_continuous_audit[ok] <- td[ok] / 7

  if (any(ok)) {
    if (ga_source == "unavailable") {
      ga_source <- paste0(total_days_col, "/7")
    } else {
      ga_source <- paste0(ga_source, " with fallback ", total_days_col, "/7")
    }
  }
}

df <- df %>%
  mutate(
    PTB_from_GA_audit = case_when(
      is.na(GA_continuous_audit) ~ NA_integer_,
      GA_continuous_audit < 37 ~ 1L,
      TRUE ~ 0L
    ),
    GA_locked_disagree = case_when(
      is.na(PTB_from_GA_audit) | is.na(PTB_locked) ~ NA,
      TRUE ~ PTB_from_GA_audit != PTB_locked
    )
  )

ga_audit_summary <- tibble(
  metric = c(
    "Total rows",
    "Locked PTB nonmissing",
    "Continuous GA reconstructable",
    "GA/locked PTB comparable",
    "GA/locked PTB disagreements",
    "Locked PTB with missing continuous GA"
  ),
  value = c(
    nrow(df),
    sum(!is.na(df$PTB_locked)),
    sum(!is.na(df$GA_continuous_audit)),
    sum(!is.na(df$GA_locked_disagree)),
    sum(df$GA_locked_disagree %in% TRUE, na.rm = TRUE),
    sum(!is.na(df$PTB_locked) & is.na(df$GA_continuous_audit))
  )
)

# Public repository version intentionally writes aggregate audit information
# only; no subject/specimen identifiers are exported.

# -----------------------------------------------------------------------------
# 4. HDP exclusion variable
# -----------------------------------------------------------------------------

# Preferred definition: preeclampsia OR gestational hypertension OR HELLP.
df$HDP_exclude <- FALSE

if (!is.na(preE_col)) {
  pree <- as_binary01(df[[preE_col]])
  df$HDP_exclude <- df$HDP_exclude | (!is.na(pree) & pree == 1L)
}

if (!is.na(gh_col)) {
  gh <- as_binary01(df[[gh_col]])
  df$HDP_exclude <- df$HDP_exclude | (!is.na(gh) & gh == 1L)
}

if (!is.na(htn_type_col)) {
  ht <- tolower(trimws(as.character(df[[htn_type_col]])))
  df$HDP_exclude <- df$HDP_exclude |
    str_detect(ht, "preeclamps|gestational hypertension|hellp")
}

# -----------------------------------------------------------------------------
# 5. PTB subtype definition
# -----------------------------------------------------------------------------

df <- df %>%
  mutate(
    DeliveryClass_raw = as.character(.data[[delivery_class_col]]),
    PTB_subtype = case_when(
      PTB_locked != 1L ~ NA_character_,

      str_detect(
        str_to_lower(DeliveryClass_raw),
        "spontaneous labor with preterm delivery"
      ) ~ "Spontaneous",

      str_detect(
        str_to_lower(DeliveryClass_raw),
        "pprom leading to spontaneous preterm delivery"
      ) ~ "Spontaneous",

      str_detect(
        str_to_lower(DeliveryClass_raw),
        "indicated preterm delivery"
      ) ~ "Indicated",

      str_detect(
        str_to_lower(DeliveryClass_raw),
        "pprom leading to induction/cesarean"
      ) ~ "Indicated",

      TRUE ~ NA_character_
    ),
    subtype_usable = case_when(
      PTB_locked == 0L ~ TRUE,            # all locked term controls retained
      PTB_locked == 1L & !is.na(PTB_subtype) ~ TRUE,
      TRUE ~ FALSE
    )
  )

subtype_audit <- df %>%
  count(
    PTB_locked,
    DeliveryClass_raw,
    PTB_subtype,
    subtype_usable,
    name = "n",
    .drop = FALSE
  ) %>%
  arrange(PTB_locked, DeliveryClass_raw, PTB_subtype)

# -----------------------------------------------------------------------------
# 6. Exposure definitions + fixed IQRs from full locked cohort
# -----------------------------------------------------------------------------

polymer_map <- tribble(
  ~Polymer,          ~PolymerLabel,                              ~ExposureCol,
  "Cumulative_MNPs", "Total MNP burden",                         "Cumulative_MNPs_log10ugpg",
  "PE",              "Polyethylene (PE)",                        "PE_log10ugpg",
  "SBR",             "Styrene-butadiene rubber (SBR)",           "SBR_log10ugpg",
  "PP",              "Polypropylene (PP)",                       "PP_log10ugpg",
  "PVC",             "Polyvinyl chloride (PVC)",                 "PVC_log10ugpg",
  "N66",             "Nylon-66 (N66)",                           "N66_log10ugpg",
  "PET",             "Polyethylene terephthalate (PET)",         "PET_log10ugpg",
  "N6",              "Nylon-6 (N6)",                             "N6_log10ugpg",
  "PMMA",            "Polymethyl methacrylate (PMMA)",           "PMMA_log10ugpg",
  "PC",              "Polycarbonate (PC)",                       "PC_log10ugpg",
  "PU",              "Polyurethane (PU)",                        "PU_log10ugpg",
  "ABS",             "Acrylonitrile-butadiene-styrene (ABS)",    "ABS_log10ugpg",
  "PS",              "Polystyrene (PS)",                         "PS_log10ugpg"
)

missing_exposures <- setdiff(polymer_map$ExposureCol, names(df))
if (length(missing_exposures) > 0) {
  warning(
    "Exposure columns missing and will be skipped: ",
    paste(missing_exposures, collapse = ", ")
  )
  polymer_map <- polymer_map %>%
    filter(ExposureCol %in% names(df))
}

primary_locked <- df %>% filter(!is.na(PTB_locked))

iqr_table <- polymer_map %>%
  mutate(
    IQR_primary = map_dbl(
      ExposureCol,
      ~ IQR(primary_locked[[.x]], na.rm = TRUE)
    )
  )

if (any(!is.finite(iqr_table$IQR_primary) | iqr_table$IQR_primary <= 0)) {
  warning("One or more exposures have zero/invalid primary-cohort IQR.")
}

# -----------------------------------------------------------------------------
# 7. Analysis subsets
# -----------------------------------------------------------------------------

analysis_sets <- list(
  "Singleton-only" = df %>%
    filter(
      !is.na(PTB_locked),
      Twins_locked == 0L
    ) %>%
    mutate(outcome = PTB_locked),

  "Exclude preeclampsia/GH/HELLP" = df %>%
    filter(
      !is.na(PTB_locked),
      !HDP_exclude
    ) %>%
    mutate(outcome = PTB_locked),

  "Spontaneous PTB vs term" = df %>%
    filter(
      !is.na(PTB_locked),
      subtype_usable,
      PTB_locked == 0L | PTB_subtype == "Spontaneous"
    ) %>%
    mutate(
      outcome = case_when(
        PTB_locked == 0L ~ 0L,
        PTB_subtype == "Spontaneous" ~ 1L,
        TRUE ~ NA_integer_
      )
    ),

  "Indicated PTB vs term" = df %>%
    filter(
      !is.na(PTB_locked),
      subtype_usable,
      PTB_locked == 0L | PTB_subtype == "Indicated"
    ) %>%
    mutate(
      outcome = case_when(
        PTB_locked == 0L ~ 0L,
        PTB_subtype == "Indicated" ~ 1L,
        TRUE ~ NA_integer_
      )
    )
)

sample_size_summary <- bind_rows(
  tibble(
    Analysis = "Locked primary cohort",
    N_total = sum(!is.na(df$PTB_locked)),
    N_events = sum(df$PTB_locked == 1L, na.rm = TRUE),
    N_non_events = sum(df$PTB_locked == 0L, na.rm = TRUE)
  ),
  imap_dfr(
    analysis_sets,
    ~ tibble(
      Analysis = .y,
      N_total = nrow(.x),
      N_events = sum(.x$outcome == 1L, na.rm = TRUE),
      N_non_events = sum(.x$outcome == 0L, na.rm = TRUE)
    )
  )
)

# -----------------------------------------------------------------------------
# 8. Model fitting
# -----------------------------------------------------------------------------

fit_one <- function(dat, analysis_name, polymer, polymer_label, exposure_col, iqr_fixed) {

  # Singleton-only model: twin term omitted because constant.
  covars <- c(
    age_col,
    bmi_col,
    sdi_col,
    prior_ptb_col,
    htn01_col,
    hospital_col
  )

  if (analysis_name != "Singleton-only") {
    covars <- c(covars, twins_col)
  }

  # Keep only covariates that exist.
  covars <- unique(covars[covars %in% names(dat)])

  vars_needed <- c("outcome", exposure_col, covars)

  d <- dat %>%
    select(all_of(vars_needed)) %>%
    filter(if_all(everything(), ~ !is.na(.x)))

  base <- tibble(
    ModelLabel = analysis_name,
    Polymer = polymer,
    PolymerLabel = polymer_label,
    ExposureCol = exposure_col,
    N = nrow(d),
    Events = if (nrow(d) > 0) sum(d$outcome == 1L) else 0L,
    NonEvents = if (nrow(d) > 0) sum(d$outcome == 0L) else 0L,
    Covariates = paste(covars, collapse = " + "),
    IQR_primary = iqr_fixed,
    OR_IQR = NA_real_,
    OR_IQR_lo = NA_real_,
    OR_IQR_hi = NA_real_,
    logOR_IQR = NA_real_,
    SE_logOR_IQR = NA_real_,
    p = NA_real_,
    fit_status = NA_character_
  )

  if (
    nrow(d) < 10 ||
    length(unique(d$outcome)) < 2 ||
    !is.finite(iqr_fixed) ||
    iqr_fixed <= 0
  ) {
    return(base %>%
      mutate(fit_status = "insufficient_data_or_invalid_IQR"))
  }

  form <- reformulate(
    termlabels = c(exposure_col, covars),
    response = "outcome"
  )

  fit <- try(
    suppressWarnings(glm(form, data = d, family = binomial())),
    silent = TRUE
  )

  if (inherits(fit, "try-error")) {
    return(base %>% mutate(fit_status = "glm_fit_error"))
  }

  # HC1 robust SE
  ct <- try(
    suppressWarnings(
      lmtest::coeftest(
        fit,
        vcov. = sandwich::vcovHC(fit, type = "HC1")
      )
    ),
    silent = TRUE
  )

  if (
    inherits(ct, "try-error") ||
    !exposure_col %in% rownames(ct)
  ) {
    return(base %>% mutate(fit_status = "HC1_coefficient_error"))
  }

  beta <- unname(ct[exposure_col, "Estimate"])
  se_beta <- unname(ct[exposure_col, "Std. Error"])
  p_val <- unname(ct[exposure_col, "Pr(>|z|)"])

  logOR <- beta * iqr_fixed
  se_logOR <- se_beta * iqr_fixed

  base %>%
    mutate(
      OR_IQR = exp(logOR),
      OR_IQR_lo = exp(logOR - 1.96 * se_logOR),
      OR_IQR_hi = exp(logOR + 1.96 * se_logOR),
      logOR_IQR = logOR,
      SE_logOR_IQR = se_logOR,
      p = p_val,
      fit_status = "ok_HC1"
    )
}

all_results <- imap_dfr(
  analysis_sets,
  function(dat, analysis_name) {
    pmap_dfr(
      iqr_table,
      function(Polymer, PolymerLabel, ExposureCol, IQR_primary) {
        fit_one(
          dat = dat,
          analysis_name = analysis_name,
          polymer = Polymer,
          polymer_label = PolymerLabel,
          exposure_col = ExposureCol,
          iqr_fixed = IQR_primary
        )
      }
    )
  }
) %>%
  group_by(ModelLabel) %>%
  mutate(
    q_BH = safe_bh(p),
    nominal_p05 = is.finite(p) & p < 0.05,
    sig_q05 = is.finite(q_BH) & q_BH < 0.05
  ) %>%
  ungroup()

key_polymers <- all_results %>%
  filter(Polymer %in% c("PU", "PET", "SBR", "ABS")) %>%
  arrange(
    factor(
      ModelLabel,
      levels = c(
        "Singleton-only",
        "Exclude preeclampsia/GH/HELLP",
        "Spontaneous PTB vs term",
        "Indicated PTB vs term"
      )
    ),
    factor(Polymer, levels = c("SBR", "PET", "PU", "ABS"))
  )

# -----------------------------------------------------------------------------
# 9. Additional audits
# -----------------------------------------------------------------------------

# Confirm term controls are retained in subtype analyses even if GA is missing.
subtype_control_audit <- tibble(
  metric = c(
    "Locked term controls total",
    "Locked term controls with missing continuous GA",
    "Term controls in spontaneous-vs-term set",
    "Term controls in indicated-vs-term set",
    "Locked PTBs with unresolved subtype"
  ),
  value = c(
    sum(df$PTB_locked == 0L, na.rm = TRUE),
    sum(df$PTB_locked == 0L & is.na(df$GA_continuous_audit), na.rm = TRUE),
    sum(analysis_sets[["Spontaneous PTB vs term"]]$outcome == 0L),
    sum(analysis_sets[["Indicated PTB vs term"]]$outcome == 0L),
    sum(df$PTB_locked == 1L & is.na(df$PTB_subtype), na.rm = TRUE)
  )
)

# -----------------------------------------------------------------------------
# 10. Write outputs
# -----------------------------------------------------------------------------

readr::write_csv(
  all_results,
  file.path(OUT_DIR, "PTB_heterogeneity_all_models_IQR.csv"),
  na = ""
)

readr::write_csv(
  key_polymers,
  file.path(OUT_DIR, "PTB_heterogeneity_key_polymers_IQR.csv"),
  na = ""
)

readr::write_csv(
  sample_size_summary,
  file.path(OUT_DIR, "PTB_heterogeneity_sample_sizes.csv"),
  na = ""
)

readr::write_csv(
  subtype_audit,
  file.path(OUT_DIR, "PTB_subtype_classification_audit.csv"),
  na = ""
)

readr::write_csv(
  ga_audit_summary,
  file.path(OUT_DIR, "PTB_GA_definition_audit_summary.csv"),
  na = ""
)


readr::write_csv(
  subtype_control_audit,
  file.path(OUT_DIR, "PTB_subtype_term_control_audit.csv"),
  na = ""
)

readr::write_csv(
  iqr_table,
  file.path(OUT_DIR, "PTB_primary_fixed_IQRs_used.csv"),
  na = ""
)

# -----------------------------------------------------------------------------
# 11. Console summary
# -----------------------------------------------------------------------------

message("Point 4 PTB heterogeneity sensitivity analyses complete.")
message("Output directory:")
message("  ", OUT_DIR)
message("")
message("Locked PTB outcome column: ", locked_ptb_col)
message("Gestational-age audit source: ", ga_source)
message(
  "GA/locked-PTB disagreements: ",
  sum(df$GA_locked_disagree %in% TRUE, na.rm = TRUE)
)
message(
  "Locked PTB subjects missing continuous GA: ",
  sum(!is.na(df$PTB_locked) & is.na(df$GA_continuous_audit))
)
message("Effect scale: fixed polymer-specific IQR from full locked primary cohort")
message("")

cat("\nSample sizes:\n")
print(sample_size_summary)

cat("\nKey polymers (PU, PET, SBR, ABS):\n")
print(
  key_polymers %>%
    select(
      ModelLabel, Polymer, N, Events,
      OR_IQR, OR_IQR_lo, OR_IQR_hi, p, q_BH
    )
)

cat("\nSubtype classification audit:\n")
print(subtype_audit)

cat("\nSubtype term-control audit:\n")
print(subtype_control_audit)
