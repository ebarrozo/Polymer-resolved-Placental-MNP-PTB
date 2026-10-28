# Script 39 — Consolidated placental MNP analytical validation supplement
# Reporting-only: does not refit models; algebraically converts frozen continuous ORs and CIs to the original-burden IQR contrast.
# Run from placenta project root: source('scripts/39_merged_placenta_analytical_suppfigure.R')
# Required packages: ggplot2, dplyr, readr. Optional: patchwork.

if (file.exists('00_config.R')) {
  source('00_config.R')
}
suppressPackageStartupMessages({library(ggplot2);library(dplyr);library(readr)})

# EDIT if the input CSVs are elsewhere; bundled files can be copied to the following folder.
INPUT_DIR <- file.path('outputs','tables','39_merged_analytical_inputs')
# The source CSVs are the frozen Script 18 analytical-validation outputs (October 2026).
# The specimen-level audit is restricted: do NOT commit it to a public GitHub repo.
OUT_DIR <- if (exists('DIR_FIGS')) file.path(DIR_FIGS,'39_merged_analytical_suppfigure') else file.path('outputs','figures','39_merged_analytical_suppfigure')
TABLE_DIR <- if (exists('DIR_TABLES')) file.path(DIR_TABLES,'39_merged_analytical_suppfigure') else file.path('outputs','tables','39_merged_analytical_suppfigure')
dir.create(OUT_DIR,recursive=TRUE,showWarnings=FALSE)
dir.create(TABLE_DIR,recursive=TRUE,showWarnings=FALSE)

# True = create supplementary 4th panel illustrating observed effect of historical
# maximum-blank subtraction on concentrations. False = streamlined three-panel figure.
INCLUDE_BLANK_PANEL <- FALSE
# The validation run's placenta spike-only controls MUST be block A, not FF block B.
CONTROL_BLOCK <- 'A'

needed <- c('spike_only_process_recovery.csv','placenta_paired_polymer_recovery.csv',
 'pu_pet_qc_category_counts_by_delivery_class.csv',
 'pu_pet_qc_binary_ptb_sensitivity_models.csv','pu_pet_qc_bias_reduced_models.csv',
 'pu_pet_historical_maxblank_subtracted_continuous_models.csv',
 'pu_pet_historical_maxblank_subtraction_specimen_audit.csv')
missing <- needed[!file.exists(file.path(INPUT_DIR,needed))]
if(length(missing)) stop('Missing Script 18 validation files in ',INPUT_DIR,':\n',paste(missing,collapse='\n'),
 '\nCopy the bundled inputs/ CSVs to that folder, or edit INPUT_DIR.')
read_in <- function(f) suppressMessages(readr::read_csv(file.path(INPUT_DIR,f),show_col_types=FALSE))
colors <- c('Spike-only (placenta run)'='#d95f02','Placenta paired spike'='#1b9e77')
polymer_levels <- c('PE','PP','PS','ABS','SBR','PMMA','PC','PVC','PU','PET','N6','N66')
base_theme <- theme_bw(base_size=10)+theme(plot.title=element_text(face='bold',size=12),
 plot.subtitle=element_text(size=9),legend.position='bottom',
 strip.text=element_text(face='bold'),plot.margin=margin(7,8,7,8))
save_plot <- function(p,stem,w,h) {
 ggsave(file.path(OUT_DIR,paste0(stem,'.png')),p,width=w,height=h,units='in',dpi=450,bg='white',limitsize=FALSE)
 ggsave(file.path(OUT_DIR,paste0(stem,'.pdf')),p,width=w,height=h,units='in',bg='white',limitsize=FALSE)
}

# A. All-polymer experimental recovery; includes ABS/SBR/PET/PU naturally.
s <- read_in('spike_only_process_recovery.csv') %>% filter(ControlBlock==CONTROL_BLOCK)
p <- read_in('placenta_paired_polymer_recovery.csv') %>% filter(paired_unspiked_available %in% TRUE)
if(!nrow(s) || !nrow(p)) stop('Placenta spike-only control block A or paired tissue recovery is missing.')
rec <- bind_rows(s %>% transmute(Polymer,Sample=as.character(ControlLabel),Recovery=process_recovery_pct,Group='Spike-only (placenta run)'),
 p %>% transmute(Polymer,Sample=as.character(BaseSampleID),Recovery=apparent_recovery_pct,Group='Placenta paired spike')) %>%
 filter(Polymer %in% polymer_levels,is.finite(Recovery)) %>% mutate(Polymer=factor(Polymer,levels=polymer_levels),Group=factor(Group,levels=names(colors)))
summ <- rec %>% group_by(Polymer,Group) %>% summarise(n=n(),median=median(Recovery),q1=quantile(Recovery,.25),q3=quantile(Recovery,.75),.groups='drop')
write_csv(rec,file.path(TABLE_DIR,'panel_A_individual_recovery.csv'))
write_csv(summ,file.path(TABLE_DIR,'panel_A_recovery_summary.csv'))
pA <- ggplot(rec,aes(Polymer,Recovery,colour=Group,shape=Group)) + geom_hline(yintercept=100,linetype=2,colour='grey50')+
 geom_point(position=position_jitterdodge(jitter.width=.12,dodge.width=.55),alpha=.65,size=1.55)+
 geom_errorbar(data=summ,aes(x=Polymer,ymin=q1,ymax=q3,group=Group),inherit.aes=FALSE,position=position_dodge(width=.55),width=.18,colour='grey25')+
 geom_point(data=summ,aes(y=median,group=Group),position=position_dodge(width=.55),size=2.8,stroke=.9,fill='white')+
 scale_colour_manual(values=colors)+scale_shape_manual(values=c(16,21))+
 scale_y_continuous(expand=expansion(mult=c(.04,.08)))+
 labs(title='A  Experimental whole-process recovery (12 polymers)',subtitle='Five paired placentas; placenta-run matrix-free spike-only controls (block A)',x=NULL,y='Apparent recovery (%)')+
 base_theme+theme(axis.text.x=element_text(angle=45,hjust=1))
save_plot(pA,'Panel_A_12polymer_recovery',10,5.1)

# B. QC confidence: original PU/PET categories. No ABS/SBR thresholds invented.
qc <- read_in('pu_pet_qc_category_counts_by_delivery_class.csv') %>% filter(Polymer %in% c('PU','PET'),PTB %in% c('Term','Preterm'))
qc_levels <- c('not_detected','below_LOB','above_LOB_below_LOD','above_LOD_below_LOQ','above_LOQ')
qc_labels <- c('Not detected','Below LOB','LOB–LOD','LOD–LOQ','Above LOQ')
qc <- qc %>% mutate(QC=factor(QC,levels=qc_levels,labels=qc_labels),Polymer=factor(Polymer,levels=c('PU','PET')),
 PTB=factor(PTB,levels=c('Term','Preterm')))
if(any(is.na(qc$QC))) stop('Unexpected QC category; inspect Script 18 data before plotting')
write_csv(qc,file.path(TABLE_DIR,'panel_B_QC_counts.csv'))
pB <- ggplot(qc,aes(PTB,pct_within_group,fill=QC))+geom_col(width=.67,colour='white',linewidth=.2)+facet_wrap(~Polymer,nrow=1)+
 scale_fill_manual(values=c('#d9d9d9','#969696','#fdae6b','#f16913','#a63603'),drop=FALSE)+
 scale_y_continuous(breaks=seq(0,100,25),limits=c(0,101),expand=expansion(mult=c(0,0)))+
 labs(title='B  PU/PET analytical-confidence categories',subtitle='Descriptive composition by delivery class',x=NULL,y='Proportion of specimens (%)',fill=NULL)+base_theme
save_plot(pB,'Panel_B_PU_PET_QC',8.2,4.5)

# C. Actual threshold associations, plus optional continuous controls.
bi <- read_in('pu_pet_qc_binary_ptb_sensitivity_models.csv') %>% filter(Polymer %in% c('PU','PET'),Threshold %in% c('at_or_above_LOD','above_LOQ')) %>%
 transmute(Polymer,Model=if_else(Threshold=='at_or_above_LOD','At/above LOD vs below','Above LOQ vs at/below'),Method='HC1 logistic',OR,CI_lo,CI_hi,p)
fi <- read_in('pu_pet_qc_bias_reduced_models.csv') %>% filter(Polymer %in% c('PU','PET'),Threshold=='above_LOQ') %>%
 transmute(Polymer,Model='Above LOQ vs at/below (Firth)',Method='Firth logistic',OR,CI_lo,CI_hi,p)
# Convert the same frozen HC1 continuous models from their original log2(1+C) unit
# into a fixed, full-model IQR contrast. This does NOT refit any model.
# Script 18's specimen audit is a restricted local input containing sample identifiers.
co_raw <- read_in('pu_pet_historical_maxblank_subtracted_continuous_models.csv') %>%
  filter(Polymer %in% c('PU','PET'),
         Sensitivity %in% c('Original_recomputed_from_analysis_ready_table','MaxBlank_subtracted_per_run'))
aud_scale <- read_in('pu_pet_historical_maxblank_subtraction_specimen_audit.csv') %>%
  filter(Polymer %in% c('PU','PET'))
required_audit_cols <- c('Polymer','original_burden_ug_g')
if (!all(required_audit_cols %in% names(aud_scale))) stop('Missing original specimen burden in restricted Script 18 audit.')
reference_iqr <- aud_scale %>%
  group_by(Polymer) %>%
  summarise(IQR_log2_1p=IQR(log2(1 + original_burden_ug_g),na.rm=TRUE),
            n_reference=sum(is.finite(original_burden_ug_g)),.groups='drop')
if(any(!is.finite(reference_iqr$IQR_log2_1p)| reference_iqr$IQR_log2_1p <=0)) stop('Invalid reference IQR')
co <- co_raw %>% left_join(reference_iqr,by='Polymer') %>%
  mutate(Model=if_else(Sensitivity=='Original_recomputed_from_analysis_ready_table',
                'Continuous: original (IQR)','Continuous: max blank subtracted (IQR)'),
         Method='HC1 logistic',
         OR=OR_per_doubling_plus1^IQR_log2_1p,
         CI_lo=CI_lo^IQR_log2_1p,
         CI_hi=CI_hi^IQR_log2_1p) %>%
  select(Polymer,Model,Method,OR,CI_lo,CI_hi,p,IQR_log2_1p,n_reference)
# Identical N, event counts and p values must be preserved from the frozen source.
if(any(!is.finite(co$OR)|!is.finite(co$CI_lo)|!is.finite(co$CI_hi))) stop('IQR conversion yielded invalid estimates')
write_csv(co,file.path(TABLE_DIR,'panel_C_continuous_IQR_conversion_audit.csv'))
models <- bind_rows(co,bi,fi) %>% filter(is.finite(OR),is.finite(CI_lo),is.finite(CI_hi),OR>0,CI_lo>0,CI_hi>0)
model_names <- c('Continuous: original (IQR)','Continuous: max blank subtracted (IQR)',
 'At/above LOD vs below','Above LOQ vs at/below','Above LOQ vs at/below (Firth)')
models <- models %>% mutate(Model=factor(Model,levels=rev(model_names)),Polymer=factor(Polymer,levels=c('PU','PET')))
write_csv(models,file.path(TABLE_DIR,'panel_C_model_estimates.csv'))
# Important: different rows quantify DIFFERENT exposure contrasts; do not compare magnitudes directly.
pC <- ggplot(models,aes(OR,Model,colour=Method))+geom_vline(xintercept=1,linetype=2,colour='grey50')+
 geom_segment(aes(x=CI_lo,xend=CI_hi,yend=Model),linewidth=.7)+geom_point(size=2.2)+
 facet_wrap(~Polymer,nrow=1)+scale_x_log10()+scale_colour_manual(values=c('HC1 logistic'='#2166ac','Firth logistic'='#b2182b'))+
 labs(title='C  PU/PET analytical-threshold sensitivity',subtitle='Continuous and threshold-model ORs use distinct exposure contrasts',x='Adjusted odds ratio (95% CI; log scale)',y=NULL,colour=NULL)+base_theme
save_plot(pC,'Panel_C_PU_PET_threshold_models',10.5,5.2)

# D. Descriptive original vs maximum historical blank subtraction, optional.
aud <- read_in('pu_pet_historical_maxblank_subtraction_specimen_audit.csv') %>% filter(Polymer %in% c('PU','PET'))
write_csv(aud,file.path(TABLE_DIR,'panel_D_specimen_blank_audit.csv'))
pD <- ggplot(aud,aes(original_burden_ug_g,blank_subtracted_ug_g)) +
 geom_abline(intercept=0,slope=1,linetype=2,colour='grey50')+
 geom_point(alpha=.55,size=1.2)+facet_wrap(~Polymer,scales='free',nrow=1)+
 scale_x_continuous(trans='log1p')+scale_y_continuous(trans='log1p')+
 labs(title='D  Conservative historical max-blank subtraction',subtitle='Specimen-level comparison; not used as primary clinical exposure',x='Original burden (µg/g; log1p axis)',y='After maximum blank subtraction (µg/g; log1p axis)')+base_theme
if(INCLUDE_BLANK_PANEL) save_plot(pD,'Panel_D_historical_maxblank_comparison',9,4.3)

if(requireNamespace('patchwork',quietly=TRUE)) {
  combined <- if(INCLUDE_BLANK_PANEL) (pA/pB)/(pC/pD) else pA/pB/pC
  save_plot(combined,'SUPP_analytical_validation_COMBINED',13,if(INCLUDE_BLANK_PANEL) 19 else 16.5)
} else message('patchwork not installed: separate PDF/PNG panels generated; combined figure skipped.')

notes <- c('Script 39: plots existing Script 18 validation outputs; NO clinical refitting.',
 'Control block A = placenta run. Do NOT include block B (follicular fluid).',
 'Panel A includes all 12 polymers, including ABS, SBR, PET and PU.',
 'Panels B/C include ONLY existing PU/PET threshold analyses; no fabricated ABS or SBR models.',
 'Panel C continuous effects are per fixed polymer-specific IQR of log2(1+original burden), using the full Script 18 analysis sample; categorical threshold ORs are distinct contrasts.',
 'Panel D is an illustrative historical-maximum-blank sensitivity, NOT empirical outcome-dependent recovery validation.',
 'Five paired placenta spike-ins, predominantly PTB, do not identify term-vs-PTB recovery differences.',
 paste('Inputs:', normalizePath(INPUT_DIR,mustWork=TRUE)),
 paste('Generated:',Sys.time()))
writeLines(notes,file.path(TABLE_DIR,'QC_READ_FIRST.txt'))
message('Finished Script 39. Figure folder: ',OUT_DIR,'; QC: ',file.path(TABLE_DIR,'QC_READ_FIRST.txt'))
