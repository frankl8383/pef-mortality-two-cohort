#!/usr/bin/env Rscript
# =============================================================================
# 18_reviewer_response_frailty_rawPEF_designC_v1_0.R
# Reviewer-response analyses for the residualized-PEF mortality paper (V39).
# Addresses three points raised in external independent review:
#   PROBLEM 1  Does the residualized-PEF marker survive frailty/vitality
#              adjustment? (NHANES; CHARLS handled in R/charls/13 already)
#   PROBLEM 2  Raw PEF vs residualized marker, BOTH on the full GLI panel
#              (the "why residualize at all" decisive HR that was missing).
#   PROBLEM 3  Design-based variance for the decisive Delta-C (M5->M6),
#              replacing the design-naive concordance SE in script 16.
#
# DATA EXPORT RULE: aggregate statistics only (HR / C / Delta-C / SE / p).
#   No row-level data are written.
# ENV: conda 'nhanes-incr' (r-survey, r-survival, r-tidyverse, r-broom)
# =============================================================================
suppressPackageStartupMessages({library(survey); library(survival); library(dplyr); library(readr)})
set.seed(20260704); options(survey.lonely.psu = "adjust")
ROOT <- getwd()
IN   <- file.path(ROOT, "derived_sensitive/nhanes/nhanes_mortality_analysis_ready_v30_0.rds")
DT   <- file.path(ROOT, "results/tables"); if(!dir.exists(DT)) dir.create(DT, recursive=TRUE)

d0 <- readRDS(IN)
req <- c("resp_vulnerability_z","pef_l_min","gli_global2022_z_fev1",
         "gli_global2022_z_fvc","gli_global2022_z_fev1_fvc","age_years","bmi")
dd <- d0 %>%
  mutate(sex_f=factor(sex), race_f=factor(race_ethnicity),
         smoking_f=factor(smoking_status, levels=c("never","former","current"))) %>%
  filter(!is.na(all_cause_death), !is.na(followup_years_exam), !is.na(wtmec6yr),
         wtmec6yr>0, followup_years_exam>0, !is.na(sex_f), !is.na(race_f), !is.na(smoking_f)) %>%
  filter(if_all(all_of(req), ~ !is.na(.x)))
dd$pef_l_min_z <- as.numeric(scale(dd$pef_l_min))
dd$w_scaled    <- dd$wtmec6yr/mean(dd$wtmec6yr)
BASE <- "age_years + sex_f + race_f + bmi + smoking_f"
GLIP <- "gli_global2022_z_fev1 + gli_global2022_z_fvc + gli_global2022_z_fev1_fvc"
M6   <- paste(BASE,"+",GLIP,"+ resp_vulnerability_z")
LHS  <- "survival::Surv(followup_years_exam, all_cause_death)"
grab <- function(fit, term){s<-summary(fit)$coefficients; if(!term%in%rownames(s)) return(c(hr=NA,lo=NA,hi=NA,p=NA))
  b<-s[term,"coef"]; se<-s[term, if("robust se"%in%colnames(s))"robust se" else "se(coef)"]
  c(hr=exp(b),lo=exp(b-1.96*se),hi=exp(b+1.96*se),p=s[term,ncol(s)])}

des <- svydesign(ids=~psu, strata=~strata, weights=~wtmec6yr, nest=TRUE, data=dd)
des <- update(des, frailty_cc=!is.na(nhanes_frailty_proxy_count),
              frailty_count=nhanes_frailty_proxy_count, low_pa=low_physical_activity)

# ---- PROBLEM 1: frailty/vitality decomposition (marker = resp_vulnerability_z) ----
fa <- svycoxph(as.formula(paste(LHS,"~",M6)), design=des)                    # full sample
des_cc <- subset(des, frailty_cc)
fb <- svycoxph(as.formula(paste(LHS,"~",M6)), design=des_cc)                 # CC subsample, no frailty
fc <- svycoxph(as.formula(paste(LHS,"~",M6,"+ frailty_count")), design=des_cc)
des_pa <- subset(des, !is.na(low_pa))
fd <- svycoxph(as.formula(paste(LHS,"~",M6,"+ low_pa")), design=des_pa)
nof <- function(x) nrow(x$variables)
p1 <- rbind(
  data.frame(model="(a) M6 full sample [headline]",             n=nof(des),    t(grab(fa,"resp_vulnerability_z"))),
  data.frame(model="(b) M6 frailty-CC subsample, no frailty",   n=nof(des_cc), t(grab(fb,"resp_vulnerability_z"))),
  data.frame(model="(c) M6 + frailty proxy count",              n=nof(des_cc), t(grab(fc,"resp_vulnerability_z"))),
  data.frame(model="(d) M6 + low physical activity",            n=nof(des_pa), t(grab(fd,"resp_vulnerability_z"))))

# ---- PROBLEM 2: raw PEF on full GLI panel (the missing decisive cell) ----
f_raw <- svycoxph(as.formula(paste(LHS,"~",paste(BASE,"+",GLIP,"+ pef_l_min_z"))), design=des)
p2 <- rbind(
  data.frame(term="residualized marker (M6)", t(grab(fa,"resp_vulnerability_z"))),
  data.frame(term="raw PEF per SD (on full GLI panel)", t(grab(f_raw,"pef_l_min_z"))))

# ---- PROBLEM 3c: is raw PEF's larger |logHR| partly HEIGHT leakage? (reviewer round-2 Q3) ----
# BASE has no height; residualized marker is height-orthogonal by construction. Add height to the
# raw-PEF model and see how much of the raw-vs-residualized |logHR| gap it closes.
alh <- function(fit,term) unname(abs(coef(fit)[term]))
f_raw_h <- svycoxph(as.formula(paste(LHS,"~",BASE,"+ height_cm +",GLIP,"+ pef_l_min_z")), design=des)
a_raw  <- alh(f_raw,"pef_l_min_z"); a_rawh <- alh(f_raw_h,"pef_l_min_z"); a_mk <- alh(fa,"resp_vulnerability_z")
pct_height <- 100*(a_raw - a_rawh)/(a_raw - a_mk)
p3c <- data.frame(
  model=c("raw PEF, no height","raw PEF + height","residualized marker"),
  absloghr=c(a_raw,a_rawh,a_mk),
  height_gap_closed_pct=c(NA,round(pct_height,1),NA),
  height_term_p=c(NA, summary(f_raw_h)$coef["height_cm","Pr(>|z|)"], NA))
write_csv(p3c, file.path(DT,"reviewer_response_p3c_height_leakage_v1_0.csv"))

# ---- PROBLEM 3: design-based variance for decisive Delta-C (M5->M6) ----
M5f <- as.formula(paste(LHS,"~",paste(BASE,"+",GLIP)))
M6f <- as.formula(paste(LHS,"~",paste(BASE,"+",GLIP,"+ resp_vulnerability_z")))
c_of <- function(f) unname(coef(survival::concordance(f)))
f5 <- coxph(M5f,data=dd,weights=w_scaled,robust=TRUE,x=TRUE,y=TRUE)
f6 <- coxph(M6f,data=dd,weights=w_scaled,robust=TRUE,x=TRUE,y=TRUE)
dC <- c_of(f6)-c_of(f5)
ctp<-survival::concordance(f5,f6); w<-c(-1,1); se_naive<-sqrt(as.numeric(t(w)%*%vcov(ctp)%*%w))
key <- paste(dd$strata, dd$psu, sep="_")
sl  <- split(unique(dd[,c("strata","psu")])$psu, unique(dd[,c("strata","psu")])$strata)
ibc <- split(seq_len(nrow(dd)), key)
B<-400; dCb<-numeric(B)
for(b in 1:B){ rows<-list()
  for(s in names(sl)){ psus<-sl[[s]]
    samp<- if(length(psus)<2) psus else sample(psus,length(psus),replace=TRUE)
    for(p in samp) rows[[length(rows)+1]]<-ibc[[paste(s,p,sep="_")]] }
  db<-dd[unlist(rows),]
  f5b<-tryCatch(coxph(M5f,data=db,weights=w_scaled,x=TRUE,y=TRUE),error=function(e)NULL)
  f6b<-tryCatch(coxph(M6f,data=db,weights=w_scaled,x=TRUE,y=TRUE),error=function(e)NULL)
  dCb[b]<- if(is.null(f5b)||is.null(f6b)) NA else c_of(f6b)-c_of(f5b) }
se_boot<-sd(dCb,na.rm=TRUE)
# delete-one-PSU jackknife (survey-standard corroboration of the bootstrap SE)
psu_keys<-unique(key); th<-setNames(numeric(length(psu_keys)),psu_keys)
sok<-setNames(unique(data.frame(strata=dd$strata,key=key))$strata,
              unique(data.frame(strata=dd$strata,key=key))$key)
for(i in seq_along(psu_keys)){ db<-dd[key!=psu_keys[i],]
  f5j<-tryCatch(coxph(M5f,data=db,weights=w_scaled),error=function(e)NULL)
  f6j<-tryCatch(coxph(M6f,data=db,weights=w_scaled),error=function(e)NULL)
  th[i]<-if(is.null(f5j)||is.null(f6j)) NA else c_of(f6j)-c_of(f5j) }
th<-th[!is.na(th)]; var_jk<-0
for(s in unique(sok[names(th)])){ ks<-names(th)[sok[names(th)]==s]
  if(length(ks)<2) next; nh<-length(ks); var_jk<-var_jk+((nh-1)/nh)*sum((th[ks]-mean(th[ks]))^2) }
se_jk<-sqrt(var_jk)
# design df for the t reference (Problem 4): #PSUs - #strata
DFdes <- length(unique(interaction(dd$strata,dd$psu,drop=TRUE))) - length(unique(dd$strata))
p3 <- data.frame(
  method=c("design-naive (script 16)","PSU bootstrap B=400","PSU jackknife"),
  deltaC=dC, se=c(se_naive,se_boot,se_jk), design_df=DFdes,
  p_normal=c(2*pnorm(-abs(dC/se_naive)), 2*pnorm(-abs(dC/se_boot)), 2*pnorm(-abs(dC/se_jk))),
  p=c(2*pt(-abs(dC/se_naive),DFdes), 2*pt(-abs(dC/se_boot),DFdes), 2*pt(-abs(dC/se_jk),DFdes)))

# ---- PROBLEM 3b: does the marker improve MODEL FIT beyond the full GLI panel? ----
# Design-based Wald test (survey-correct LRT analogue) + working weighted LRT for reference.
M5s <- as.formula(paste(LHS,"~",paste(BASE,"+",GLIP)))
fit5s <- svycoxph(M5s, design=des); fit6s <- svycoxph(M6f, design=des)
wald  <- survey::regTermTest(fit6s, ~resp_vulnerability_z, method="Wald")
fit5n <- coxph(M5s, data=dd, weights=w_scaled, robust=FALSE)
fit6n <- coxph(M6f, data=dd, weights=w_scaled, robust=FALSE)
lrt   <- anova(fit5n, fit6n)
p3b <- data.frame(
  test  = c("design-based Wald (marker term)","working weighted LRT (reference)"),
  stat  = c(unname(wald$Ftest), unname(lrt$Chisq[2])),
  df    = c(paste(wald$df, wald$ddf, sep=","), as.character(lrt$Df[2])),
  p     = c(wald$p, lrt$`Pr(>|Chi|)`[2]))
write_csv(p3b, file.path(DT,"reviewer_response_p3b_LRT_v1_0.csv"))

# ---- PROBLEM 2b: M2-base decisive ΔC (M5->M6) with design-based variance (unify variance idiom) ----
# Reconstruct the M2 socioeconomic-base sample (adds education + income-to-poverty; N=6,389).
req2 <- c("resp_vulnerability_z","pef_l_min","gli_global2022_z_fev1","gli_global2022_z_fvc",
          "gli_global2022_z_fev1_fvc","age_years","bmi","income_poverty_ratio")
dd2 <- d0 %>% mutate(sex_f=factor(sex),race_f=factor(race_ethnicity),
        smoking_f=factor(smoking_status,levels=c("never","former","current")),
        education_f=factor(education)) %>%
  filter(!is.na(all_cause_death),!is.na(followup_years_exam),!is.na(wtmec6yr),wtmec6yr>0,
         followup_years_exam>0,!is.na(sex_f),!is.na(race_f),!is.na(smoking_f),!is.na(education_f)) %>%
  filter(if_all(all_of(req2),~!is.na(.x)))
dd2$w_scaled <- dd2$wtmec6yr/mean(dd2$wtmec6yr)
BASE2 <- "age_years + sex_f + race_f + bmi + smoking_f + education_f + income_poverty_ratio"
m5_2 <- as.formula(paste(LHS,"~",BASE2,"+",GLIP)); m6_2 <- as.formula(paste(LHS,"~",BASE2,"+",GLIP,"+ resp_vulnerability_z"))
# paired concordance ΔC (reproduces script-17 point value 0.004311)
# paired concordance needs the data reachable by NAME (concordance re-evaluates data=);
# assign the resampled frame to a stable global name so it is found in both single & paired calls.
dCpair2 <- function(d){ assign(".rs_d", d, envir=.GlobalEnv)
  a<-coxph(m5_2,data=.rs_d,weights=.rs_d$w_scaled,robust=TRUE,x=TRUE,y=TRUE)
  b<-coxph(m6_2,data=.rs_d,weights=.rs_d$w_scaled,robust=TRUE,x=TRUE,y=TRUE)
  cc<-coef(concordance(a,b)); as.numeric(cc[2]-cc[1]) }
dC2 <- dCpair2(dd2)
a2<-coxph(m5_2,data=dd2,weights=dd2$w_scaled,robust=TRUE,x=TRUE,y=TRUE); b2<-coxph(m6_2,data=dd2,weights=dd2$w_scaled,robust=TRUE,x=TRUE,y=TRUE)
ctp2<-concordance(a2,b2); se2_naive<-sqrt(as.numeric(t(c(-1,1))%*%vcov(ctp2)%*%c(-1,1)))
key2<-interaction(dd2$strata,dd2$psu,drop=TRUE)
sk2<-split(seq_len(nrow(dd2)),dd2$strata); pbs2<-lapply(sk2,function(ix) split(ix,key2[ix]))
set.seed(20260704); dCb2<-numeric(400)
for(bi in 1:400){ idx<-integer(0)
  for(s in names(pbs2)){ cl<-pbs2[[s]]; if(length(cl)<1)next; pk<-sample(names(cl),length(cl),replace=TRUE); idx<-c(idx,unlist(cl[pk])) }
  dCb2[bi]<-tryCatch(dCpair2(dd2[idx,]),error=function(e)NA) }
dCb2<-dCb2[!is.na(dCb2)]; se2_boot<-sd(dCb2)
pk2<-unique(key2); so2<-setNames(unique(data.frame(st=dd2$strata,k=key2))$st, unique(data.frame(st=dd2$strata,k=key2))$k)
th2<-setNames(rep(NA_real_,length(pk2)),pk2)
for(i in seq_along(pk2)) th2[i]<-tryCatch(dCpair2(dd2[key2!=pk2[i],]),error=function(e)NA)
th2<-th2[!is.na(th2)]; vjk2<-0
for(s in unique(so2[names(th2)])){ ks<-names(th2)[so2[names(th2)]==s]; if(length(ks)<2)next; nh<-length(ks); vjk2<-vjk2+((nh-1)/nh)*sum((th2[ks]-mean(th2[ks]))^2) }
se2_jk<-sqrt(vjk2)
DFdes2 <- length(unique(interaction(dd2$strata,dd2$psu,drop=TRUE))) - length(unique(dd2$strata))
p2b <- data.frame(base="M2 socioeconomic", contrast="M5->M6 decisive", n=nrow(dd2), events=sum(dd2$all_cause_death),
  deltaC=dC2, method=c("design-naive","PSU bootstrap B=400","PSU jackknife"),
  se=c(se2_naive,se2_boot,se2_jk), design_df=DFdes2,
  p_normal=c(2*pnorm(-abs(dC2/se2_naive)),2*pnorm(-abs(dC2/se2_boot)),2*pnorm(-abs(dC2/se2_jk))),
  p=c(2*pt(-abs(dC2/se2_naive),DFdes2),2*pt(-abs(dC2/se2_boot),DFdes2),2*pt(-abs(dC2/se2_jk),DFdes2)))
write_csv(p2b, file.path(DT,"reviewer_response_p2b_M2base_designC_v1_0.csv"))
if(exists(".rs_d", envir=.GlobalEnv)) rm(".rs_d", envir=.GlobalEnv)

# ---- write aggregate outputs ----
write_csv(p1, file.path(DT,"reviewer_response_p1_frailty_v1_0.csv"))
write_csv(p2, file.path(DT,"reviewer_response_p2_rawvsresid_v1_0.csv"))
write_csv(p3, file.path(DT,"reviewer_response_p3_designC_v1_0.csv"))
cat("PROBLEM 1 (frailty/vitality):\n"); print(p1, row.names=FALSE)
cat("\nPROBLEM 2 (raw vs residualized on full GLI panel):\n"); print(p2, row.names=FALSE)
cat("\nPROBLEM 3 (design-based Delta-C variance):\n"); print(p3, row.names=FALSE)
cat("\nPROBLEM 3b (model-fit / LRT beyond full GLI panel):\n"); print(p3b, row.names=FALSE)
cat("\nPROBLEM 2b (M2-base decisive ΔC, design-based variance):\n"); print(p2b, row.names=FALSE)
cat("\nPROBLEM 3c (height-leakage in raw-vs-residualized |logHR|):\n"); print(p3c, row.names=FALSE)
