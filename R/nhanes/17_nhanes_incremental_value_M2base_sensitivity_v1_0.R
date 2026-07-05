#!/usr/bin/env Rscript
# =============================================================================
# 17_nhanes_incremental_value_M2base_sensitivity_v1_0.R
# SENSITIVITY to script 16: repeat the incremental-value ladder on the M2
# socioeconomic covariate base (adds education + income-to-poverty ratio),
# which reduces the complete-case sample. Confirms the incremental findings are
# not an artifact of the M1 covariate choice.
# Outputs (AGGREGATE ONLY): results/tables/incremental_value_M2base_models_v1_0.csv
#                           results/tables/incremental_value_M2base_deltaC_v1_0.csv
#                           results/figures/incremental_value_M2base_cindex_v1_0.{png,pdf}
# =============================================================================
suppressPackageStartupMessages({library(survey);library(survival);library(dplyr);library(readr);library(tibble);library(ggplot2)})
set.seed(20260704); options(survey.lonely.psu="adjust"); VER<-"v1_0"
ROOT<-getwd(); if(!dir.exists(file.path(ROOT,"derived_sensitive"))){a<-commandArgs(FALSE);fa<-grep("^--file=",a,value=TRUE);if(length(fa))ROOT<-normalizePath(file.path(dirname(sub("^--file=","",fa[1])),"..",".."),mustWork=FALSE)}
IN<-file.path(ROOT,"derived_sensitive/nhanes/nhanes_mortality_analysis_ready_v30_0.rds")
DT<-file.path(ROOT,"results/tables"); DF<-file.path(ROOT,"results/figures")
for(d in c(DT,DF)) if(!dir.exists(d)) dir.create(d,recursive=TRUE)

d0<-readRDS(IN)
req<-c("resp_vulnerability_z","pef_l_min","gli_global2022_z_fev1","gli_global2022_z_fvc","gli_global2022_z_fev1_fvc","age_years","bmi","income_poverty_ratio")
dd<-d0 %>% mutate(sex_f=factor(sex),race_f=factor(race_ethnicity),
     smoking_f=factor(smoking_status,levels=c("never","former","current")),
     education_f=factor(education)) %>%
  filter(!is.na(all_cause_death),!is.na(followup_years_exam),!is.na(wtmec6yr),wtmec6yr>0,
         followup_years_exam>0,!is.na(sex_f),!is.na(race_f),!is.na(smoking_f),!is.na(education_f)) %>%
  filter(if_all(all_of(req),~!is.na(.x)))
dd$pef_l_min_z<-as.numeric(scale(dd$pef_l_min))
dd$w_scaled<-dd$wtmec6yr/mean(dd$wtmec6yr)
n_used<-nrow(dd); n_events<-sum(dd$all_cause_death==1L)
message(sprintf("M2-base common sample: N=%d, events=%d",n_used,n_events))

BASE<-"age_years + sex_f + race_f + bmi + smoking_f + education_f + income_poverty_ratio"
GLIP<-"gli_global2022_z_fev1 + gli_global2022_z_fvc + gli_global2022_z_fev1_fvc"
specs<-tibble::tribble(~id,~label,~rhs,
 "M0","Base (M2 socioeconomic)",BASE,
 "M1","Base + raw PEF (per SD)",paste(BASE,"+ pef_l_min_z"),
 "M2","Base + GLI FEV1 z",paste(BASE,"+ gli_global2022_z_fev1"),
 "M3","Base + residualized PEF marker",paste(BASE,"+ resp_vulnerability_z"),
 "M4","Base + GLI FEV1 z + PEF marker",paste(BASE,"+ gli_global2022_z_fev1 + resp_vulnerability_z"),
 "M5","Base + full GLI panel",paste(BASE,"+",GLIP),
 "M6","Base + full GLI panel + PEF marker",paste(BASE,"+",GLIP,"+ resp_vulnerability_z"))
des<-survey::svydesign(ids=~psu,strata=~strata,weights=~wtmec6yr,nest=TRUE,data=dd)
LHS<-"survival::Surv(followup_years_exam, all_cause_death)"
svy<-lapply(specs$rhs,function(r)survey::svycoxph(as.formula(paste(LHS,"~",r)),design=des)); names(svy)<-specs$id
wtd<-lapply(specs$rhs,function(r)survival::coxph(as.formula(paste(LHS,"~",r)),data=dd,weights=w_scaled,robust=TRUE,x=TRUE,y=TRUE)); names(wtd)<-specs$id
c_one<-function(f){ct<-survival::concordance(f);c(C=unname(coef(ct)),se=unname(sqrt(as.numeric(vcov(ct)))))}
cidx<-t(vapply(wtd,c_one,numeric(2)))|>as.data.frame()|>tibble::rownames_to_column("id")
grab<-function(f,term){s<-summary(f)$coefficients;if(!term%in%rownames(s))return(c(hr=NA,lo=NA,hi=NA,p=NA));b<-s[term,"coef"];se<-s[term,if("robust se"%in%colnames(s))"robust se" else "se(coef)"];c(hr=exp(b),lo=exp(b-1.96*se),hi=exp(b+1.96*se),p=s[term,ncol(s)])}
hrM3<-grab(svy[["M3"]],"resp_vulnerability_z");hrM4<-grab(svy[["M4"]],"resp_vulnerability_z");hrM6<-grab(svy[["M6"]],"resp_vulnerability_z")
models<-specs|>left_join(cidx,by="id")|>mutate(n=n_used,events=n_events,C_lo=C-1.96*se,C_hi=C+1.96*se,C_ci=sprintf("%.3f (%.3f-%.3f)",C,C_lo,C_hi))|>select(id,label,n,events,C,C_se=se,C_lo,C_hi,C_ci)
dC<-function(a,b,lab){ct<-survival::concordance(wtd[[a]],wtd[[b]]);cc<-coef(ct);V<-vcov(ct);w<-c(-1,1);d<-as.numeric(w%*%cc);se<-sqrt(as.numeric(t(w)%*%V%*%w));z<-d/se;tibble::tibble(contrast=lab,model_a=a,model_b=b,C_a=unname(cc[1]),C_b=unname(cc[2]),deltaC=d,se=se,lo=d-1.96*se,hi=d+1.96*se,z=z,p=2*pnorm(-abs(z)))}
delta<-dplyr::bind_rows(dC("M0","M3","PEF marker vs base"),dC("M0","M2","GLI FEV1 z vs base"),dC("M1","M3","Residualized marker vs raw PEF"),dC("M2","M4","PEF marker on top of GLI FEV1 z"),dC("M5","M6","PEF marker on top of FULL GLI panel [DECISIVE]"))
readr::write_csv(models,file.path(DT,paste0("incremental_value_M2base_models_",VER,".csv")))
readr::write_csv(delta,file.path(DT,paste0("incremental_value_M2base_deltaC_",VER,".csv")))
plot_df<-models|>mutate(id=factor(id,levels=rev(specs$id)))
p<-ggplot(plot_df,aes(x=C,y=id))+geom_errorbarh(aes(xmin=C_lo,xmax=C_hi),height=0.16,colour="grey40")+geom_point(size=2.6,colour="#7a3b1f")+geom_text(aes(label=label),hjust=0,nudge_y=0.30,size=2.9,colour="grey20")+scale_y_discrete(expand=expansion(add=c(0.6,0.9)))+labs(x="Weighted Harrell's C (95% CI)",y=NULL,title="Incremental discrimination (M2 socioeconomic base) \u2014 sensitivity",subtitle=sprintf("NHANES 2007-2012; complete-case M2 sample N=%d, deaths=%d",n_used,n_events))+theme_minimal(base_size=11)+theme(plot.title.position="plot",panel.grid.minor=element_blank())
ggsave(file.path(DF,paste0("incremental_value_M2base_cindex_",VER,".png")),plot=p,width=7.6,height=4.2,dpi=600,units="in")
ggsave(file.path(DF,paste0("incremental_value_M2base_cindex_",VER,".pdf")),plot=p,width=7.6,height=4.2,units="in",device=grDevices::pdf,useDingbats=FALSE)
cat("\n===== M2-BASE SENSITIVITY =====\n"); cat(sprintf("N=%d deaths=%d\n",n_used,n_events))
print(as.data.frame(models[,c("id","label","C_ci")]),row.names=FALSE)
cat("\n-- Delta-C --\n"); print(as.data.frame(delta[,c("contrast","deltaC","lo","hi","p")]),row.names=FALSE)
cat(sprintf("\nMarker HR: M3 %.2f(%.2f-%.2f)p=%.3g | M4|FEV1z %.2f(%.2f-%.2f)p=%.3g | M6|panel %.2f(%.2f-%.2f)p=%.3g\n",
 hrM3["hr"],hrM3["lo"],hrM3["hi"],hrM3["p"],hrM4["hr"],hrM4["lo"],hrM4["hi"],hrM4["p"],hrM6["hr"],hrM6["lo"],hrM6["hi"],hrM6["p"]))