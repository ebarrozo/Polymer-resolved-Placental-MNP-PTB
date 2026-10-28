#!/usr/bin/env Rscript
# Script 35 — Placental PU/PET primary-model reconciliation (NO recovery corrections).
# Run from the placenta project root, containing 00_config.R:
# source('scripts/35_reconcile_PU_PET_primary_models.R') or source('35_reconcile_PU_PET_primary_models.R')
# No original inputs are modified; no participant identifiers are exported.
source('00_config.R')
required <- c('readr','dplyr','sandwich','lmtest')
missing <- required[!vapply(required,requireNamespace,logical(1),quietly=TRUE)]
if(length(missing)) stop('Install: ',paste(missing,collapse=', '))
MODEL_FILE <- file.path(DIR_TABLES,'mnp_ptb_cohort__analysis_ready_model_table.csv')
OUT <- file.path(DIR_TABLES,'35_PU_PET_primary_model_reconciliation')
dir.create(OUT,recursive=TRUE,showWarnings=FALSE)
if(!file.exists(MODEL_FILE))stop('Missing model-ready file: ',MODEL_FILE)
d <- as.data.frame(readr::read_csv(MODEL_FILE,show_col_types=FALSE,guess_max=10000))
first <- function(x,required=TRUE){ z<-x[x%in%names(d)]; if(!length(z)){if(required) stop('Missing columns: ',paste(x,collapse='/'));return(NA_character_)};z[1] }
col_ptb<-first(c('Preterm','PTB_locked','PTB01','ptb01'))
cv <- c(age=first(c('Maternal_age','maternal_age')),bmi=first(c('BMI','bmi')),
        sdi=first(c('SDI_score','sdi_score')),twins=first(c('Twins01','twin_pregnancy')),
        prior=first(c('Prior_PTB','prior_ptb')),hypertension=first(c('Hypertension01','hypertensive_disorder')),
        hospital=first(c('Hospital','hospital')))
parse_binary <- function(x) {z<-tolower(trimws(as.character(x)));ans<-rep(NA_real_,length(z));ans[z%in%c('1','true','yes','preterm','ptb')]<-1;ans[z%in%c('0','false','no','term')]<-0;ans}
d$.ptb<-parse_binary(d[[col_ptb]])
for(j in cv[names(cv)!='hospital']) d[[j]]<-suppressWarnings(as.numeric(d[[j]]))
d[[cv['hospital']]]<-as.factor(d[[cv['hospital']]])
poly_map<-c(PU='PU_ug_per_g',PET='PET_ug_per_g')
if(!all(poly_map%in%names(d)))stop('Missing raw PU/PET concentration columns')
for(p in names(poly_map)){
 x<-suppressWarnings(as.numeric(d[[poly_map[p]]]))
 if(any(x[is.finite(x)]<0))stop('Negative concentration for ',p)
 d[[paste0('.x_',p)]]<-log10(1+pmax(x,0))
}
# Define the complete-case samples used in Script 34, and audit attrition.
allvars <- unname(cv)
missingness<-data.frame(field=c('Locked PTB',names(cv),'PU concentration','PET concentration'),
  column=c(col_ptb,unname(cv),unname(poly_map)),
  n_missing=c(sum(is.na(d$.ptb)),vapply(d[unname(cv)],function(x)sum(is.na(x)),integer(1)),
              vapply(d[unname(poly_map)],function(x)sum(!is.finite(suppressWarnings(as.numeric(x)))),integer(1))))
readr::write_csv(missingness,file.path(OUT,'01_missingness_by_field.csv'))
# Sequential exclusions show the effect of each covariate ON THIS ORDER, not unique attribution.
seqvars<-c('.ptb',unname(cv));ok<-rep(TRUE,nrow(d));attrition<-data.frame(stage='All rows',n=nrow(d),n_ptb=sum(d$.ptb==1,na.rm=TRUE))
for(k in seq_along(seqvars)){
 ok<-ok & !is.na(d[[seqvars[k]]])
 attrition<-rbind(attrition,data.frame(stage=c('Locked outcome',paste0('Add ',names(cv)))[k],n=sum(ok),n_ptb=sum(d$.ptb[ok]==1)))
}
readr::write_csv(attrition,file.path(OUT,'02_sequential_complete_case_attrition.csv'))
# Primary candidate fits: covariate reduction is a diagnostic, not automatic model selection.
variants<-list(
 Script34_all7=unname(cv),
 Exclude_twins=unname(cv[names(cv)!='twins']),
 Exclude_hypertension=unname(cv[names(cv)!='hypertension']),
 Exclude_twins_and_hypertension=unname(cv[!names(cv)%in%c('twins','hypertension')]),
 Exclude_prior=unname(cv[names(cv)!='prior']),
 Exclude_SDI=unname(cv[names(cv)!='sdi']),
 Exclude_hospital=unname(cv[names(cv)!='hospital']),
 Minimal_age_BMI_SDI_prior_hospital=unname(cv[c('age','bmi','sdi','prior','hospital')])
)
fit_model<-function(poly,covars,label,subset_rule='candidate'){
 xc<-paste0('.x_',poly)
 columns<-c('.ptb',xc,covars)
 cc<-complete.cases(d[,columns,drop=FALSE]); dd<-d[cc,,drop=FALSE]
 ret<-data.frame(polymer=poly,model=label,N=nrow(dd),events=sum(dd$.ptb==1),non_events=sum(dd$.ptb==0),
   covariates=paste(covars,collapse=' + '),IQR_full_locked=NA_real_,IQR_complete_case=NA_real_,
   OR_fixed_IQR=NA_real_,CI_lo_fixed_IQR=NA_real_,CI_hi_fixed_IQR=NA_real_,
   OR_per_log10_2=NA_real_,CI_lo_per_log10_2=NA_real_,CI_hi_per_log10_2=NA_real_,
   beta_log10_1px=NA_real_,SE_HC1=NA_real_,p=NA_real_,fit_status='not_fit')
 allx<-d[[xc]][!is.na(d$.ptb)]
 ret$IQR_full_locked<-IQR(allx,na.rm=TRUE)
 if(nrow(dd))ret$IQR_complete_case<-IQR(dd[[xc]],na.rm=TRUE)
 if(nrow(dd)<30||length(unique(dd$.ptb))!=2||sd(dd[[xc]])==0){ret$fit_status<-'insufficient';return(ret)}
 formula<-stats::reformulate(c(xc,covars),response='.ptb')
 model<-tryCatch(suppressWarnings(glm(formula,data=dd,family=binomial())),error=function(e)NULL)
 if(is.null(model)){ret$fit_status<-'glm_error';return(ret)}
 ct<-tryCatch(lmtest::coeftest(model,vcov.=sandwich::vcovHC(model,type='HC1')),error=function(e)NULL)
 if(is.null(ct)||!xc%in%rownames(ct)){ret$fit_status<-'HC1_error';return(ret)}
 beta<-unname(ct[xc,1]);se<-unname(ct[xc,2]);pv<-unname(ct[xc,4]);ret$beta_log10_1px<-beta;ret$SE_HC1<-se;ret$p<-pv
 scale_effect<-function(z) c(exp(beta*z),exp((beta-1.96*se)*z),exp((beta+1.96*se)*z))
 v<-scale_effect(ret$IQR_full_locked);ret[c('OR_fixed_IQR','CI_lo_fixed_IQR','CI_hi_fixed_IQR')]<-as.list(v)
 v<-scale_effect(log10(2));ret[c('OR_per_log10_2','CI_lo_per_log10_2','CI_hi_per_log10_2')]<-as.list(v)
 ret$fit_status<-'ok';ret
}
fits<-do.call(rbind,lapply(names(poly_map),function(p) do.call(rbind,lapply(names(variants),function(v)fit_model(p,variants[[v]],v)))))
readr::write_csv(fits,file.path(OUT,'03_primary_model_reconciliation.csv'))
# Validate Script 34 (if present): exactly identical baseline n, OR, CI, P at machine precision.
f34<-file.path(DIR_TABLES,'34_PU_PET_recovery_sensitivity','01_all_model_fits.csv')
check<-data.frame(polymer=names(poly_map),Script34_exists=file.exists(f34),baseline_N_match=NA,baseline_OR_match=NA,baseline_p_match=NA,max_abs_difference=NA_real_)
if(file.exists(f34)){
 old<-as.data.frame(readr::read_csv(f34,show_col_types=FALSE))
 for(i in seq_len(nrow(check))){p<-check$polymer[i];a<-old[old$polymer==p & old$scenario=='Original measured',,drop=FALSE];b<-fits[fits$polymer==p & fits$model=='Script34_all7',,drop=FALSE]
   if(nrow(a)!=1||nrow(b)!=1)stop('Script 34 baseline does not contain exactly one row for ',p)
   vals<-abs(c(a$OR-b$OR_fixed_IQR,a$CI_low-b$CI_lo_fixed_IQR,a$CI_high-b$CI_hi_fixed_IQR,a$p-b$p))
   check$baseline_N_match[i]<-identical(as.integer(a$N),as.integer(b$N))
   check$baseline_OR_match[i]<-isTRUE(all.equal(a$OR,b$OR_fixed_IQR,tolerance=1e-8))
   check$baseline_p_match[i]<-isTRUE(all.equal(a$p,b$p,tolerance=1e-8))
   check$max_abs_difference[i]<-max(vals,na.rm=TRUE)
 }
 if(!all(check$baseline_N_match & check$baseline_OR_match & check$baseline_p_match))stop('Script 35 did not reproduce the Script 34 baseline. Inspect check file.')
}
readr::write_csv(check,file.path(OUT,'04_Script34_baseline_reproduction.csv'))
# All reference values below are discussion targets from prior manuscript versions,
# not authoritative machine-read primary model results. DO NOT declare primary match by OR alone.
reference<-data.frame(polymer=c('PU','PET'),manuscript_historical_OR_per_twofold=c(2.00,1.65),
 note='Historical stated estimate; verify exact analysis sample, covariates and exposure definition from manuscript final-model CSV')
readr::write_csv(reference,file.path(OUT,'05_historical_manuscript_OR_reference_NOT_VALIDATED.csv'))
report<-c('SCRIPT 35 — PU/PET primary-model reconciliation',
 paste('Input model table:',MODEL_FILE),paste('Rows in model table:',nrow(d)),
 paste('Locked PTB events:',sum(d$.ptb==1,na.rm=TRUE)),
 paste('Outcome column:',col_ptb),paste('Primary covariates:',paste(names(cv),'=',cv,collapse='; ')),
 'The per-log10(2) OR is for a twofold increase in (1 + burden), NOT an exact doubling of µg/g near zero.',
 'The baseline uses the same covariate set and full locked-cohort IQR as Script 34.',
 'Candidate covariate reductions are diagnostic comparisons; DO NOT select variants based on matching the published OR.',
 'HC1 robust standard errors used. No recoveries adjusted or simulated.',
 '', 'SEQUENTIAL ATTRITION:',capture.output(print(attrition,row.names=FALSE)),
 '', 'BASELINES + CANDIDATE FITS:',capture.output(print(fits[,c('polymer','model','N','events','OR_fixed_IQR','OR_per_log10_2','p','fit_status')],row.names=FALSE)),
 '', 'SCRIPT 34 REPRODUCTION:',capture.output(print(check,row.names=FALSE)),
 '', 'PUBLICATION BLOCKER: Please supply the final primary PU/PET model-results CSV or the script that generated final manuscript Table S1. The historical per-2x values alone cannot establish equivalence.')
writeLines(report,file.path(OUT,'QC_READ_FIRST.txt'))
message('Script 35 complete. Review: ',file.path(OUT,'QC_READ_FIRST.txt'))
