#!/usr/bin/env Rscript

# ============================================================
# 34_primary_complete_case_audit.R
# Complete-case and polymer-specific sample-size audit
# ============================================================

source("00_config.R")

req_pkgs <- c("dplyr","readr","tibble","purrr","stringr")
missing_pkgs <- req_pkgs[!vapply(req_pkgs, requireNamespace, logical(1), quietly=TRUE)]
if (length(missing_pkgs)>0) stop("Missing required packages: ", paste(missing_pkgs, collapse=", "), call.=FALSE)

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tibble); library(purrr); library(stringr)
})

INFILE <- file.path(DIR_TABLES, "mnp_ptb_cohort__analysis_ready_model_table.csv")
OUTDIR <- file.path(DIR_TABLES, "34_primary_complete_case_audit")
dir.create(OUTDIR, recursive=TRUE, showWarnings=FALSE)
if (!file.exists(INFILE)) stop("Input file not found: ", INFILE, call.=FALSE)

first_existing <- function(dat, candidates) {
  hit <- candidates[candidates %in% names(dat)]
  if (length(hit)==0) NA_character_ else hit[[1]]
}
safe_num <- function(x) suppressWarnings(as.numeric(x))
as_ptb01 <- function(x) {
  if (is.logical(x)) return(ifelse(is.na(x),NA_real_,as.numeric(x)))
  if (is.numeric(x) || is.integer(x)) {
    out <- suppressWarnings(as.numeric(x)); out[!out %in% c(0,1)] <- NA_real_; return(out)
  }
  z <- trimws(tolower(as.character(x)))
  case_when(z %in% c("1","true","yes","preterm","ptb") ~ 1,
            z %in% c("0","false","no","term") ~ 0,
            TRUE ~ NA_real_)
}
fmt_med_iqr <- function(x) {
  x <- safe_num(x); x <- x[is.finite(x)]
  if (length(x)==0) return(NA_character_)
  q <- quantile(x,c(.25,.50,.75),na.rm=TRUE,type=7)
  sprintf("%.3g [%.3g-%.3g]", q[2],q[1],q[3])
}
smd_cont <- function(x,g) {
  x <- safe_num(x); ok <- is.finite(x) & !is.na(g); x<-x[ok]; g<-g[ok]
  if (length(unique(g))!=2) return(NA_real_)
  x0<-x[g==0]; x1<-x[g==1]
  if (length(x0)<2 || length(x1)<2) return(NA_real_)
  sp <- sqrt((var(x0)+var(x1))/2)
  if (!is.finite(sp) || sp==0) return(NA_real_)
  (mean(x1)-mean(x0))/sp
}
smd_binary <- function(x,g) {
  x <- safe_num(x); ok <- x %in% c(0,1) & !is.na(g); x<-x[ok]; g<-g[ok]
  if (length(unique(g))!=2) return(NA_real_)
  p0<-mean(x[g==0]==1); p1<-mean(x[g==1]==1)
  den<-sqrt((p0*(1-p0)+p1*(1-p1))/2)
  if (!is.finite(den) || den==0) return(NA_real_)
  (p1-p0)/den
}

dat <- read_csv(INFILE, show_col_types=FALSE, guess_max=10000)

outcome_col <- first_existing(dat,c("Preterm","PTB_locked","PTB01","ptb01"))
age_col   <- first_existing(dat,c("Maternal_age","maternal_age"))
bmi_col   <- first_existing(dat,c("BMI","bmi"))
sdi_col   <- first_existing(dat,c("SDI_score","sdi_score"))
twins_col <- first_existing(dat,c("Twins01","twin_pregnancy"))
prior_col <- first_existing(dat,c("Prior_PTB","prior_ptb"))
htn_col   <- first_existing(dat,c("Hypertension01","hypertensive_disorder"))
hosp_col  <- first_existing(dat,c("Hospital","hospital"))

covars <- c(age_col,bmi_col,sdi_col,twins_col,prior_col,htn_col,hosp_col)
if (is.na(outcome_col) || any(is.na(covars))) stop("Could not resolve outcome/covariates.", call.=FALSE)

dat <- dat %>%
  mutate(
    PTB01_locked=as_ptb01(.data[[outcome_col]]),
    primary_covariate_complete=complete.cases(across(all_of(covars))),
    primary_analysis_eligible=!is.na(PTB01_locked),
    complete_case_group=case_when(
      !primary_analysis_eligible ~ NA_character_,
      primary_covariate_complete ~ "Included_complete_case",
      TRUE ~ "Excluded_missing_covariate"
    )
  )

audit_base <- dat %>%
  filter(primary_analysis_eligible, !is.na(complete_case_group)) %>%
  mutate(excluded01=if_else(complete_case_group=="Excluded_missing_covariate",1,0))

group_counts <- audit_base %>%
  count(complete_case_group, PTB01_locked, name="n") %>%
  mutate(Delivery=if_else(PTB01_locked==1,"Preterm","Term")) %>%
  select(complete_case_group,Delivery,n)
write_csv(group_counts, file.path(OUTDIR,"complete_case_group_counts_by_delivery.csv"), na="")

missingness <- tibble(
  Variable=covars,
  N_missing=vapply(covars, function(v) sum(is.na(audit_base[[v]])), numeric(1))
) %>% mutate(N_total=nrow(audit_base), pct_missing=100*N_missing/N_total)
write_csv(missingness, file.path(OUTDIR,"primary_covariate_missingness.csv"), na="")

continuous_vars <- c(
  age_col,bmi_col,sdi_col,
  first_existing(dat,c("GA","GA_weeks","ga_delivery_weeks")),
  first_existing(dat,c("Birth_weight_g","birth_weight_g")),
  "PU_ug_per_g","PET_ug_per_g","SBR_ug_per_g","ABS_ug_per_g"
)
continuous_vars <- unique(continuous_vars[!is.na(continuous_vars) & continuous_vars %in% names(audit_base)])

cont_rows <- map_dfr(continuous_vars, function(v) {
  x <- safe_num(audit_base[[v]]); g <- audit_base$excluded01
  p <- tryCatch(wilcox.test(x ~ g, exact=FALSE)$p.value, error=function(e) NA_real_)
  tibble(
    Variable=v, Type="continuous",
    Included_summary=fmt_med_iqr(x[g==0]),
    Excluded_summary=fmt_med_iqr(x[g==1]),
    Included_nonmissing=sum(is.finite(x[g==0])),
    Excluded_nonmissing=sum(is.finite(x[g==1])),
    SMD_excluded_minus_included=smd_cont(x,g),
    p_value=p
  )
})

binary_vars <- c(twins_col,prior_col,htn_col)
bin_rows <- map_dfr(binary_vars, function(v) {
  x <- safe_num(audit_base[[v]]); g <- audit_base$excluded01
  tab <- table(factor(g,levels=c(0,1)),factor(x,levels=c(0,1)))
  p <- tryCatch(fisher.test(tab)$p.value,error=function(e) NA_real_)
  n0<-sum(g==0 & x %in% c(0,1)); n1<-sum(g==1 & x %in% c(0,1))
  p0<-ifelse(n0>0,mean(x[g==0]==1,na.rm=TRUE),NA_real_)
  p1<-ifelse(n1>0,mean(x[g==1]==1,na.rm=TRUE),NA_real_)
  tibble(
    Variable=v, Type="binary",
    Included_summary=ifelse(is.finite(p0),sprintf("%d/%d (%.1f%%)",sum(x[g==0]==1,na.rm=TRUE),n0,100*p0),NA_character_),
    Excluded_summary=ifelse(is.finite(p1),sprintf("%d/%d (%.1f%%)",sum(x[g==1]==1,na.rm=TRUE),n1,100*p1),NA_character_),
    Included_nonmissing=n0, Excluded_nonmissing=n1,
    SMD_excluded_minus_included=smd_binary(x,g),
    p_value=p
  )
})

outcome_tab <- table(factor(audit_base$excluded01,levels=c(0,1)), factor(audit_base$PTB01_locked,levels=c(0,1)))
outcome_p <- fisher.test(outcome_tab)$p.value
outcome_row <- tibble(
  Variable="Locked_PTB_outcome", Type="binary",
  Included_summary=sprintf("%d/%d (%.1f%% PTB)",
                           sum(audit_base$PTB01_locked[audit_base$excluded01==0]==1),
                           sum(audit_base$excluded01==0),
                           100*mean(audit_base$PTB01_locked[audit_base$excluded01==0]==1)),
  Excluded_summary=sprintf("%d/%d (%.1f%% PTB)",
                           sum(audit_base$PTB01_locked[audit_base$excluded01==1]==1),
                           sum(audit_base$excluded01==1),
                           100*mean(audit_base$PTB01_locked[audit_base$excluded01==1]==1)),
  Included_nonmissing=sum(audit_base$excluded01==0),
  Excluded_nonmissing=sum(audit_base$excluded01==1),
  SMD_excluded_minus_included=smd_binary(audit_base$PTB01_locked,audit_base$excluded01),
  p_value=outcome_p
)

comparison <- bind_rows(outcome_row,cont_rows,bin_rows) %>%
  mutate(p_BH_descriptive=p.adjust(p_value,method="BH"))
write_csv(comparison,file.path(OUTDIR,"complete_case_included_vs_excluded_comparison.csv"),na="")

hospital_summary <- audit_base %>%
  count(complete_case_group,.data[[hosp_col]],name="n") %>%
  group_by(complete_case_group) %>% mutate(pct=100*n/sum(n)) %>% ungroup()
names(hospital_summary)[2] <- "Hospital"
write_csv(hospital_summary,file.path(OUTDIR,"complete_case_hospital_distribution.csv"),na="")

polymer_map <- tibble::tribble(
  ~Polymer,~ExposureCol,
  "Cumulative_MNPs","Cumulative_MNPs_log10ugpg",
  "PE","PE_log10ugpg","SBR","SBR_log10ugpg","PP","PP_log10ugpg",
  "PVC","PVC_log10ugpg","N66","N66_log10ugpg","PET","PET_log10ugpg",
  "N6","N6_log10ugpg","PMMA","PMMA_log10ugpg","PC","PC_log10ugpg",
  "PU","PU_log10ugpg","ABS","ABS_log10ugpg"
)

polymer_n <- polymer_map %>%
  rowwise() %>%
  mutate(
    N_primary_covariate_complete=sum(audit_base$primary_covariate_complete,na.rm=TRUE),
    N_exposure_nonmissing_among_complete=if (ExposureCol %in% names(audit_base))
      sum(audit_base$primary_covariate_complete & is.finite(safe_num(audit_base[[ExposureCol]])),na.rm=TRUE)
      else NA_integer_,
    N_missing_exposure_among_complete=if (is.finite(N_exposure_nonmissing_among_complete))
      N_primary_covariate_complete-N_exposure_nonmissing_among_complete
      else NA_integer_
  ) %>% ungroup()
write_csv(polymer_n,file.path(OUTDIR,"polymer_specific_primary_model_sample_sizes.csv"),na="")

writeLines(c(
  "34_primary_complete_case_audit.R",
  paste0("Input: ",INFILE),
  paste0("Outcome column: ",outcome_col),
  "",
  "Group counts:", capture.output(print(group_counts,n=Inf)),
  "",
  "Missingness:", capture.output(print(missingness,n=Inf)),
  "",
  "Included vs excluded:", capture.output(print(comparison,n=Inf)),
  "",
  "Polymer-specific model N:", capture.output(print(polymer_n,n=Inf))
), file.path(OUTDIR,"34_primary_complete_case_audit_run_log.txt"))

message("Done: ",OUTDIR)
