#!/usr/bin/env Rscript
# Script 37: four-polymer PTB nonlinearity, robust tests, influence checks, figures and tables.
# Execute from placenta project root (where 00_config.R resides):
# source('scripts/37_four_polymer_nonlinearity_robustness.R')
# Primary observed concentrations are never altered. No multiple-model winner selection.

if (!file.exists('00_config.R')) stop('Set working directory to the placenta project root containing 00_config.R')
source('00_config.R')
req <- c('ggplot2','sandwich','lmtest','scales')
miss <- req[!vapply(req,requireNamespace,logical(1),quietly=TRUE)]
if(length(miss)) stop('Install required packages: ',paste(miss,collapse=', '))
if(!exists('DIR_TABLES') || !exists('DIR_FIGS')) stop('00_config.R must define DIR_TABLES and DIR_FIGS')
infile <- file.path(DIR_TABLES,'mnp_ptb_cohort__analysis_ready_model_table.csv')
if(!file.exists(infile)) stop('Missing input: ',infile)
out <- file.path(DIR_TABLES,'37_four_polymer_nonlinearity'); figs <- file.path(DIR_FIGS,'37_four_polymer_nonlinearity')
dir.create(out,recursive=TRUE,showWarnings=FALSE);dir.create(figs,recursive=TRUE,showWarnings=FALSE)
raw <- read.csv(infile,check.names=FALSE,stringsAsFactors=FALSE)
find_col <- function(cands) { x <- intersect(cands,names(raw)); if(!length(x))stop('Missing: ',paste(cands,collapse=' / '));x[1] }
cn <- c(PTB=find_col(c('Preterm','PTB_locked','PTB01','ptb01')),
        age=find_col(c('Maternal_age','maternal_age')),
        bmi=find_col(c('BMI','bmi')),
        sdi=find_col(c('SDI_score','sdi_score')),
        twins=find_col(c('Twins01','twin_pregnancy')),
        prior=find_col(c('Prior_PTB','prior_ptb')),
        htn=find_col(c('Hypertension01','hypertensive_disorder')),
        hosp=find_col(c('Hospital','hospital')))
polys <- c(ABS='ABS_ug_per_g',SBR='SBR_ug_per_g',PET='PET_ug_per_g',PU='PU_ug_per_g')
if(!all(polys %in% names(raw)))stop('Missing polymer burden columns: ',paste(setdiff(polys,names(raw)),collapse=', '))
ptb01 <- function(x){ z<-tolower(trimws(as.character(x)));ans<-rep(NA_integer_,length(z));ans[z%in%c('1','yes','true','ptb','preterm')]<-1L;ans[z%in%c('0','no','false','term')]<-0L;ans }
num <- function(x) suppressWarnings(as.numeric(as.character(x)))
COV <- c('age','bmi','sdi','twins','prior','htn','hosp')
base <- data.frame(PTB=ptb01(raw[[cn['PTB']]]),age=num(raw[[cn['age']]]),bmi=num(raw[[cn['bmi']]]),sdi=num(raw[[cn['sdi']]]),twins=num(raw[[cn['twins']]]),prior=num(raw[[cn['prior']]]),htn=num(raw[[cn['htn']]]),hosp=factor(raw[[cn['hosp']]]))
# Predeclared sensitivity sets: interior knot/boundaries 10/50/90, 5/50/95, 20/50/80;
# influence removal 1 maximum or upper 2.5% of exposure (no optimized cutoff).
settings <- list(primary=c(.10,.50,.90),wide=c(.05,.50,.95),narrow=c(.20,.50,.80))
summary_rows<-list(); curves<-list(); density_rows<-list(); audit<-list(); rnum<-0; cnum<-0
create_z <- function(x, k, which=2L, projection=NULL) {
  bas <- splines::ns(x,knots=k[2],Boundary.knots=k[c(1,3)])
  if(is.null(projection)) projection<-as.numeric(qr.solve(cbind(1,x),bas[,which]))
  list(z=as.numeric(bas[,which]-projection[1]-projection[2]*x),projection=projection)
}
fit_one <- function(d, poly, test, knot_prop, ref_knots=NULL, ref_projection=NULL, origin='full'){
  if(nrow(d)<45 || length(unique(d$PTB))<2) return(NULL)
  x<-d$x
  k<-if(is.null(ref_knots))as.numeric(quantile(x,knot_prop,names=FALSE)) else ref_knots
  if(length(unique(k))<3 || k[1]>=k[2] || k[2]>=k[3]) return(NULL)
  zdat<-tryCatch(create_z(x,k,projection=ref_projection),error=function(e)NULL)
  if(is.null(zdat))return(NULL)
  d$nonlinear_z<-zdat$z
  # A linear + nonlinear spline spanning the same subspace as a 1-knot natural cubic spline.
  if(qr(cbind(1,x,d$nonlinear_z))$rank<3) return(NULL)
  f0<-as.formula(paste('PTB ~',paste(COV,collapse='+')))
  f1<-as.formula(paste('PTB ~ x +',paste(COV,collapse='+')))
  f2<-as.formula(paste('PTB ~ x + nonlinear_z +',paste(COV,collapse='+')))
  a<-tryCatch(suppressWarnings(glm(f0,binomial(),data=d)),error=function(e)NULL)
  l<-tryCatch(suppressWarnings(glm(f1,binomial(),data=d)),error=function(e)NULL)
  s<-tryCatch(suppressWarnings(glm(f2,binomial(),data=d)),error=function(e)NULL)
  if(any(vapply(list(a,l,s),is.null,logical(1))))return(NULL)
  if(any(!is.finite(coef(s))))return(NULL)
  vc<-tryCatch(sandwich::vcovHC(s,type='HC1'),error=function(e)NULL)
  vl<-tryCatch(sandwich::vcovHC(l,type='HC1'),error=function(e)NULL)
  if(is.null(vc)||is.null(vl)||!all(is.finite(vc)))return(NULL)
  b<-coef(s)[c('x','nonlinear_z')]; V<-vc[c('x','nonlinear_z'),c('x','nonlinear_z'),drop=FALSE]
  wald_all<-tryCatch(as.numeric(t(b)%*%solve(V,b)),error=function(e)NA_real_)
  z_stat<-coef(s)['nonlinear_z']/sqrt(vc['nonlinear_z','nonlinear_z'])
  lr<-tryCatch(anova(l,s,test='Chisq'),error=function(e)NULL)
  nonlinear_lrt<-if(is.null(lr))NA_real_ else lr$`Pr(>Chi)`[2]
  iqr<-IQR(x)
  beta<-coef(l)['x']; se<-sqrt(vl['x','x'])
  record<-data.frame(polymer=poly,analysis=test,subset=origin,N=nrow(d),events=sum(d$PTB==1),
    minimum_ug_g=min(d$burden),maximum_ug_g=max(d$burden),IQR_log10_1px=iqr,
    fixed_knots=paste(signif(k,7),collapse=';'),
    linear_OR_IQR=exp(beta*iqr),linear_CI_lo=exp((beta-1.96*se)*iqr),linear_CI_hi=exp((beta+1.96*se)*iqr),
    p_linear_HC1=2*pnorm(-abs(beta/se)),
    p_nonlinearity_HC1=2*pnorm(-abs(z_stat)),p_nonlinearity_LRT=nonlinear_lrt,
    p_overall_spline_HC1=if(is.finite(wald_all))pchisq(wald_all,df=2,lower.tail=FALSE) else NA_real_,
    AIC_linear=AIC(l),AIC_spline=AIC(s),delta_AIC_linear_minus_spline=AIC(l)-AIC(s),
    max_std_pearson_resid=max(abs(residuals(s,type='pearson')),na.rm=TRUE),stringsAsFactors=FALSE)
  if(test=='primary' && origin=='full') {
    # Curve is relative to median transformed exposure, holding covariates identical at reference.
    qlo<-quantile(d$burden,.01); qhi<-quantile(d$burden,.99)
    grid<-seq(qlo,qhi,length.out=250); gx<-log10(1+pmax(grid,0))
    gspec<-create_z(gx,k,projection=zdat$projection)
    reference<-median(x); rz<-create_z(reference,k,projection=zdat$projection)$z
    vX<-cbind(gx-reference,gspec$z-rz); colnames(vX)<-c('x','nonlinear_z')
    pred<-as.numeric(vX%*%b); pse<-sqrt(pmax(0,rowSums((vX%*%V)*vX)))
    attr(record,'curve')<-data.frame(polymer=poly,burden_ug_g=grid,OR_vs_median=exp(pred),lo=exp(pred-1.96*pse),hi=exp(pred+1.96*pse),reference_ug_g=10^reference-1)
    attr(record,'knots')<-k;attr(record,'projection')<-zdat$projection
  }
  record
}
for(poly in names(polys)){
 d<-base;d$burden<-num(raw[[polys[[poly]]]]);d$x<-log10(1+pmax(d$burden,0))
 d<-d[complete.cases(d[,c('PTB','x',COV)]),,drop=FALSE]
 if(nrow(d)<45)stop('Too few complete cases for ',poly)
 if(any(d$burden<0,na.rm=TRUE))warning('Negative burden clipped to zero in transformed exposure for ',poly)
 audit[[poly]]<-data.frame(polymer=poly,N=nrow(d),PTB=sum(d$PTB),term=sum(d$PTB==0),zero=sum(d$burden==0),
        p95_ug_g=as.numeric(quantile(d$burden,.95)),max_ug_g=max(d$burden),stringsAsFactors=FALSE)
 for(spec in names(settings)){
  result<-fit_one(d,poly,spec,settings[[spec]])
  if(is.null(result))stop('Spline fit failed for ',poly,' / ',spec,' (inspect knot uniqueness / model separation)')
  rnum<-rnum+1;summary_rows[[rnum]]<-result
  if(spec=='primary'){
   cnum<-cnum+1;curves[[cnum]]<-attr(result,'curve')
   knots<-attr(result,'knots');proj<-attr(result,'projection')
   # Hold original knots and projection fixed when excluding high-exposure subjects.
   ord<-order(d$x,decreasing=TRUE);d1<-d[-ord[1],,drop=FALSE]
   threshold<-as.numeric(quantile(d$x,.975,names=FALSE));d2<-d[d$x<=threshold,,drop=FALSE]
   for(i in seq_along(list(d1,d2))){di<-list(d1,d2)[[i]];nm<-c('remove_max','remove_top_2.5pct')[i]
    rr<-fit_one(di,poly,'primary_knots_fixed',settings$primary,ref_knots=knots,ref_projection=proj,origin=nm)
    if(is.null(rr)) warning(poly,' ',nm,' could not be fit') else {rnum<-rnum+1;summary_rows[[rnum]]<-rr}
   }
  }
 }
 density_rows[[poly]]<-data.frame(polymer=poly,burden_ug_g=d$burden,PTB=d$PTB)
}
result<-do.call(rbind,summary_rows);curve<-do.call(rbind,curves);dens<-do.call(rbind,density_rows);aud<-do.call(rbind,audit)
write.csv(result,file.path(out,'SPLINE_MODEL_COMPARISON_AND_INFLUENCE.csv'),row.names=FALSE)
write.csv(curve,file.path(out,'SPLINE_CURVE_PLOTTED_VALUES.csv'),row.names=FALSE)
write.csv(aud,file.path(out,'COHORT_AND_TAIL_SUPPORT_AUDIT.csv'),row.names=FALSE)
# BH q-values applied only to predefined primary nonlinearity tests across four polymers.
prim<-result[result$analysis=='primary'&result$subset=='full',,drop=FALSE]
prim$q_nonlinearity_HC1_BH4<-p.adjust(prim$p_nonlinearity_HC1,method='BH',n=4)
prim$q_nonlinearity_LRT_BH4<-p.adjust(prim$p_nonlinearity_LRT,method='BH',n=4)
write.csv(prim,file.path(out,'PRIMARY_SPLINE_TESTS_4_POLYMERS.csv'),row.names=FALSE)

library(ggplot2)
poly_levels<-c('ABS','SBR','PET','PU')
curve$polymer<-factor(curve$polymer,levels=poly_levels)
dens$polymer<-factor(dens$polymer,levels=poly_levels)
p_curve<-ggplot(curve,aes(burden_ug_g,OR_vs_median))+
 geom_ribbon(aes(ymin=lo,ymax=hi),fill='#4575a1',alpha=.18)+
 geom_line(color='#345d87',linewidth=.85)+geom_hline(yintercept=1,linetype=2,color='grey45')+
 geom_rug(data=dens,aes(x=burden_ug_g),inherit.aes=FALSE,alpha=.3,sides='b')+
 facet_wrap(~polymer,scales='free_x',ncol=2)+
 scale_x_continuous(trans=scales::pseudo_log_trans(base=10))+
 scale_y_log10()+
 labs(title='A  Adjusted exposure–PTB associations by polymer',subtitle='Natural cubic spline (one median knot; 10th/90th percentile boundary knots) | reference: median burden',
      x='Placental polymer burden (µg/g; pseudo-log scale)',y='Adjusted OR relative to median burden (log scale)')+
 theme_bw(base_size=11)+theme(panel.grid.minor=element_blank())

plot_sens<-result[result$subset!='full',,drop=FALSE]
plot_sens$polymer<-factor(plot_sens$polymer,levels=poly_levels)
plot_sens$subset<-factor(plot_sens$subset,levels=c('remove_max','remove_top_2.5pct'))
# Two complementary summaries: predefined test p-values and stability of fitted shape.
p_tests<-ggplot(prim,aes(x=factor(polymer,levels=poly_levels),y=-log10(pmax(p_nonlinearity_HC1,1e-10))))+
 geom_col(fill='#4c729b',width=.6)+geom_hline(yintercept=-log10(.05),linetype=2)+
 geom_text(aes(label=sprintf('q=%.3f',q_nonlinearity_HC1_BH4)),vjust=-.5,size=3.3)+
 coord_cartesian(ylim=c(0,max(-log10(pmax(prim$p_nonlinearity_HC1,1e-10)))+0.6))+
 labs(title='B  HC1-robust nonlinearity test',subtitle='One-degree-of-freedom Wald test; BH correction across four polymers',x=NULL,y=expression(-log[10](italic(p))))+theme_bw(base_size=11)

p_influence<-ggplot(result[result$analysis%in%c('primary','primary_knots_fixed'),],
                    aes(x=subset,y=delta_AIC_linear_minus_spline,group=polymer,color=polymer))+
 geom_hline(yintercept=0,linetype=2)+geom_point(size=2)+geom_line() +
 facet_wrap(~polymer,ncol=2,scales='free_y')+
 scale_x_discrete(limits=c('full','remove_max','remove_top_2.5pct'),labels=c('Full','Exclude max','Exclude upper 2.5%'))+
 labs(title='C  Influence of upper exposure tail',subtitle='Linear-versus-spline AIC difference; original spline knots held fixed',x=NULL,y='AIC(linear) − AIC(spline)')+
 theme_bw(base_size=10)+theme(legend.position='none',axis.text.x=element_text(angle=20,hjust=1))

for(p in list(list(n='SPLINE_4_POLYMERS',o=p_curve,w=13,h=9),list(n='ROBUST_NONLINEARITY_TESTS',o=p_tests,w=9,h=5.5),list(n='INFLUENCE_AND_TAIL_SENSITIVITY',o=p_influence,w=11,h=7))){
  ggsave(file.path(figs,paste0(p$n,'.png')),p$o,width=p$w,height=p$h,dpi=400)
  ggsave(file.path(figs,paste0(p$n,'.pdf')),p$o,width=p$w,height=p$h)
}
if(requireNamespace('patchwork',quietly=TRUE)){
 combo<-p_curve / (p_tests | p_influence) + patchwork::plot_layout(heights=c(1.1,1))
 ggsave(file.path(figs,'FOUR_POLYMER_NONLINEARITY_COMBINED.png'),combo,width=17,height=15,dpi=400)
 ggsave(file.path(figs,'FOUR_POLYMER_NONLINEARITY_COMBINED.pdf'),combo,width=17,height=15)
} else message('patchwork not installed; individual panels saved. install.packages("patchwork") for combined figure.')
if(requireNamespace('writexl',quietly=TRUE))writexl::write_xlsx(list('Primary spline tests'=prim,'All sensitivity fits'=result,'Tail support'=aud),file.path(out,'NONLINEARITY_TABLES.xlsx'))
if(requireNamespace('officer',quietly=TRUE)&&requireNamespace('flextable',quietly=TRUE)){
  ft<-flextable::flextable(prim[,c('polymer','N','linear_OR_IQR','linear_CI_lo','linear_CI_hi','p_linear_HC1','p_nonlinearity_HC1','q_nonlinearity_HC1_BH4','p_overall_spline_HC1','delta_AIC_linear_minus_spline')]);ft<-flextable::autofit(ft)
  doc<-officer::read_docx();doc<-officer::body_add_par(doc,'Supplementary Table. Polymer-specific linear and nonlinear PTB associations',style='heading 1');doc<-flextable::body_add_flextable(doc,ft);print(doc,target=file.path(out,'NONLINEARITY_PRIMARY_TABLE.docx'))
}
writeLines(c('Script 37: 4-polymer nonlinearity and influence sensitivity',
  paste('Input:',infile),paste('Covariates:',paste(paste(names(cn),cn,sep='='),collapse='; ')),
  'Primary spline knots: p10, p50, p90 on log10(1+burden); alternatives p5/p50/p95 and p20/p50/p80.',
  'Nonlinearity: HC1 robust 1-df Wald test of nonlinear spline component; LRT is supplementary.',
  'Overall association: HC1 robust joint 2-df Wald test; BH nonlinearity q computed over four polymers.',
  'Influence: remove highest exposure or observations above 97.5th percentile; primary spline basis kept fixed.',
  'Curves are adjusted OR versus median exposure, using identical covariates at both exposure values.',
  'Caution: tail curves are imprecise; nominal shape tests and adjusted associations do not imply causality.',
  'This script does not select models according to significance or replace primary linear logistic models.',
  paste('Result rows:',nrow(result),'and primary:',nrow(prim))),file.path(out,'QC_READ_FIRST.txt'))
message('Script 37 complete. Tables: ',out,' | Figures: ',figs)
