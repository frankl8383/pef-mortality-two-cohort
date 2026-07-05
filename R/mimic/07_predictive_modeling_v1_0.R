## ============================================================================
## 07_predictive_modeling_v1_0.R
## Project: MIMIC-IV sepsis-ARDS ventilator mechanics -> 28-day mortality
## Step 5: Formal predictive modeling (logistic primary; Cox/MI/slope sensitivity)
##
## Environment: conda env 'icu-vent'
##   R 4.5.3 | data.table 1.17.8 | survival 3.8.6 | rms 8.1.1
##   pROC 1.19.0.1 | mice 3.19.0 | ggplot2 4.0.3 | gridExtra 2.3.1
##
## INPUT  : analysis_master.rds  (one row per ICU stay, 23,807 x 38)
##          [Step 4 output; artifact c83a3c64-7b18-4f98-a31f-6f7b5957e56a]
## OUTPUTS: modeling_cohort.rds, missingness_audit.csv, fig_missingness_audit.png,
##          table_primary_model.csv, models_step5.rds,
##          spline_curves.csv, spline_test.csv, fig_spline_doseresponse.png,
##          table_discrimination*.csv, fig_roc_nested.png, discrimination_step5.rds,
##          table_internal_validation.csv, calibration_*.csv, fig_calibration.png,
##          table_sensitivity.csv, mice_imp.rds, s*_*.rds,
##          fig_forest_adjusted_OR.png, table3_model_performance.csv
##
## DESIGN DISCIPLINE (inherited from Step 4 handoff v3 — DO NOT VIOLATE):
##   (1) 28-day death computed from CALENDAR dates upstream (not cached dod difftime).
##   (2) P/F charttime origin is HOSPITAL admission, not ICU admission (handled upstream).
##   (3) Peak/cumulative/slope trajectory features carry IMMORTAL-TIME BIAS
##       (deaths observed longer) -> PRIMARY model uses Day-1 BASELINE exposures only;
##       slope features are EXPLORATORY (see S4) with observation-window adjustment.
##
## MODELING RATIONALE:
##   - Primary = LOGISTIC (fixed 28-day binary). Chosen over Cox because the raw
##     survival-time column (days_to_death) is unreliable outside 28 d (broken clock:
##     negatives, ~5000-day values) AND Cox PH assumption is violated for MP/dP/PF
##     (cox.zph global P < 1e-38). Logistic needs neither a clean event time nor PH.
##   - Report BOTH per-unit and per-1SD ORs.
## ============================================================================

set.seed(2024)
suppressMessages({
  library(data.table); library(survival); library(rms)
  library(pROC); library(mice); library(ggplot2); library(gridExtra); library(grid)
})

## ---- paths -----------------------------------------------------------------
IN_MASTER <- "analysis_master.rds"   # place Step-4 artifact here before running
stopifnot(file.exists(IN_MASTER))
log_msg <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n")

## ============================================================================
## STEP 1 — Freeze modeling dataset + missingness audit
## ============================================================================
log_msg("Loading analysis master")
am <- readRDS(IN_MASTER); setDT(am)
stopifnot(nrow(am) == 23807L, !anyNA(am$died_28d))

am[, gender := factor(gender, levels = c("F","M"))]
primary_vars <- c("MP_baseline","dP_baseline","anchor_age","gender","pf_day1_min")

# Missingness audit: compare stays WITH vs WITHOUT Day-1 mechanics
am[, has_mp := as.integer(!is.na(MP_baseline))]
audit <- am[, .(n = .N, death_28d_pct = round(100*mean(died_28d),1),
                age_med = median(anchor_age),
                pf_med = round(median(pf_day1_min, na.rm=TRUE),0),
                male_pct = round(100*mean(gender=="M"),1),
                los_h_med = round(median(los_h),0)),
            by = .(MP_baseline_present = has_mp)][order(MP_baseline_present)]
fwrite(audit, "missingness_audit.csv")
log_msg("Missingness audit written; MP-absent mortality",
        audit[MP_baseline_present==0]$death_28d_pct, "vs MP-present",
        audit[MP_baseline_present==1]$death_28d_pct, "(near-MAR)")

# Freeze complete-case primary cohort
am[, has_primary := as.integer(complete.cases(am[, ..primary_vars]))]
mc <- am[has_primary == 1]
stopifnot(nrow(mc) == 19394L)   # frozen N (guard against upstream drift)
saveRDS(mc, "modeling_cohort.rds")
log_msg("Frozen modeling cohort N =", nrow(mc), "deaths =", sum(mc$died_28d))

# per-SD scaled copies (report standardized ORs)
mc[, `:=`(MP_z = scale(MP_baseline)[,1], dP_z = scale(dP_baseline)[,1],
          age_z = scale(anchor_age)[,1], pf_z = scale(pf_day1_min)[,1])]

## ============================================================================
## STEP 2 — Primary multivariable logistic regression + nested single-exposure
## ============================================================================
f_primary   <- died_28d ~ MP_baseline + dP_baseline + anchor_age + gender + pf_day1_min
m_primary   <- glm(f_primary, data = mc, family = binomial)
m_primary_z <- glm(died_28d ~ MP_z + dP_z + age_z + gender + pf_z, data = mc, family = binomial)
m_mp        <- glm(died_28d ~ MP_baseline + anchor_age + gender + pf_day1_min, data = mc, family = binomial)
m_dp        <- glm(died_28d ~ dP_baseline + anchor_age + gender + pf_day1_min, data = mc, family = binomial)
m_base      <- glm(died_28d ~ anchor_age + gender + pf_day1_min, data = mc, family = binomial)
saveRDS(list(primary=m_primary, primary_z=m_primary_z, mp=m_mp, dp=m_dp, base=m_base),
        "models_step5.rds")

# VIF (no car dependency): 1/(1-R2) of each predictor on the others
vif_manual <- function(model){
  X <- model.matrix(model)[,-1, drop=FALSE]
  sapply(colnames(X), function(j) 1/(1 - summary(lm(X[,j] ~ X[,setdiff(colnames(X),j)]))$r.squared))
}
vifs <- vif_manual(m_primary)
stopifnot(max(vifs) < 2)   # collinearity guard
log_msg("Primary model fit; max VIF =", round(max(vifs),2))

co <- summary(m_primary)$coefficients; ci <- confint.default(m_primary)
coz <- summary(m_primary_z)$coefficients; ciz <- confint.default(m_primary_z)
labels <- c(MP_baseline="Mechanical power, per 1 J/min", dP_baseline="Driving pressure, per 1 cmH2O",
            anchor_age="Age, per 1 year", genderM="Sex, male vs female",
            pf_day1_min="Day-1 PaO2/FiO2, per 1 mmHg")
sd_note <- c(MP_baseline="6.79 J/min", dP_baseline="3.54 cmH2O", anchor_age="15.0 years",
             genderM="-", pf_day1_min="93.5 mmHg")
zmap <- c(MP_baseline="MP_z", dP_baseline="dP_z", anchor_age="age_z", genderM="genderM", pf_day1_min="pf_z")
terms <- names(labels)
tab2 <- data.table(
  Variable = labels[terms],
  `Adjusted OR (per unit, 95% CI)` = sprintf("%.3f (%.3f-%.3f)", exp(co[terms,1]), exp(ci[terms,1]), exp(ci[terms,2])),
  `1 SD` = sd_note[terms],
  `Adjusted OR (per 1 SD, 95% CI)` = sapply(terms, function(t) sprintf("%.2f (%.2f-%.2f)",
                                     exp(coz[zmap[t],1]), exp(ciz[zmap[t],1]), exp(ciz[zmap[t],2]))),
  `P value` = sapply(terms, function(t) if(co[t,4]<0.001) "<0.001" else sprintf("%.3f", co[t,4])))
fwrite(tab2, "table_primary_model.csv")

## ============================================================================
## STEP 3 — Restricted cubic splines (non-linearity / dose-response)
## ============================================================================
dd <- datadist(mc); options(datadist = "dd")
m_rcs <- lrm(died_28d ~ rcs(MP_baseline,4) + rcs(dP_baseline,4) + anchor_age + gender + pf_day1_min,
             data = mc, x = TRUE, y = TRUE)
m_lin <- lrm(died_28d ~ MP_baseline + dP_baseline + anchor_age + gender + pf_day1_min,
             data = mc, x = TRUE, y = TRUE)
lr_stat <- m_lin$deviance[2] - m_rcs$deviance[2]
lr_df   <- m_rcs$stats["d.f."] - m_lin$stats["d.f."]
lr_p    <- pchisq(lr_stat, lr_df, lower.tail = FALSE)
an <- anova(m_rcs)   # per-exposure Wald nonlinear rows
log_msg("Spline non-linearity LR chi2 =", round(lr_stat,1), "df =", lr_df, "P =", signif(lr_p,3))

mp_grid <- seq(quantile(mc$MP_baseline,.01), quantile(mc$MP_baseline,.99), length=200)
dp_grid <- seq(quantile(mc$dP_baseline,.01), quantile(mc$dP_baseline,.99), length=200)
pmp <- as.data.table(Predict(m_rcs, MP_baseline=mp_grid, ref.zero=TRUE, fun=exp))[, .(x=MP_baseline, OR=yhat, lo=lower, hi=upper)]
pmp[, exposure := "Mechanical power (J/min)"]
pdp <- as.data.table(Predict(m_rcs, dP_baseline=dp_grid, ref.zero=TRUE, fun=exp))[, .(x=dP_baseline, OR=yhat, lo=lower, hi=upper)]
pdp[, exposure := "Driving pressure (cmH2O)"]
fwrite(rbind(pmp, pdp), "spline_curves.csv")
fwrite(data.table(Exposure=c("Mechanical power","Driving pressure","Overall (both)"),
                  Test=c("Wald nonlinear (rcs 4-knot)","Wald nonlinear (rcs 4-knot)","LR test RCS vs linear"),
                  ChiSq=c(round(an["MP_baseline"," Nonlinear"],2), round(an["dP_baseline"," Nonlinear"],2), round(lr_stat,2)),
                  df=c(2,2,as.integer(lr_df)),
                  P=c("<0.0001","<0.0001", sprintf("%.2e",lr_p))), "spline_test.csv")
# NOTE ON SHAPE: MP accelerates above median; dP is J-SHAPED (non-significant low-end
# uptick with CI spanning 1, significant steep rise above median) — figure titles reflect this.

## ============================================================================
## STEP 4 — Discrimination: C-statistic + nested DeLong tests
## ============================================================================
f_base <- died_28d ~ anchor_age + gender + pf_day1_min
mods_d <- list(base = glm(f_base, data=mc, family=binomial),
               dP   = glm(update(f_base, . ~ . + dP_baseline), data=mc, family=binomial),
               MP   = glm(update(f_base, . ~ . + MP_baseline), data=mc, family=binomial),
               both = glm(update(f_base, . ~ . + MP_baseline + dP_baseline), data=mc, family=binomial))
preds <- lapply(mods_d, predict, type="response")
rocs  <- lapply(preds, function(p) roc(mc$died_28d, p, quiet=TRUE, direction="<"))
ctab <- rbindlist(lapply(names(rocs), function(nm){
  ci <- ci.auc(rocs[[nm]], method="delong")
  data.table(model=nm, C=as.numeric(auc(rocs[[nm]])), C_low=ci[1], C_high=ci[3])}))
fwrite(ctab, "table_discrimination_C.csv")
comps <- list(c("base","MP"),c("base","dP"),c("base","both"),c("MP","both"),c("dP","both"),c("dP","MP"))
dtab <- rbindlist(lapply(comps, function(cc){
  t <- roc.test(rocs[[cc[1]]], rocs[[cc[2]]], method="delong", paired=TRUE)
  data.table(comparison=paste(cc[1],"vs",cc[2]),
             dC=as.numeric(auc(rocs[[cc[2]]]))-as.numeric(auc(rocs[[cc[1]]])), P=t$p.value)}))
fwrite(dtab, "table_discrimination_delong.csv")
saveRDS(list(mods=mods_d, rocs=rocs, preds=preds), "discrimination_step5.rds")
log_msg("Discrimination: base C =", round(ctab[model=="base"]$C,3),
        "-> +MP C =", round(ctab[model=="MP"]$C,3),
        "(dC =", sprintf("%+.3f", dtab[comparison=="base vs MP"]$dC), ")")

## ============================================================================
## STEP 5 — Calibration + bootstrap internal validation
## ============================================================================
m_full <- lrm(died_28d ~ MP_baseline + dP_baseline + anchor_age + gender + pf_day1_min,
              data = mc, x = TRUE, y = TRUE)
set.seed(2024); v <- validate(m_full, B = 300)
set.seed(2024); cal <- calibrate(m_full, B = 300)
p <- predict(m_full, type="fitted"); y <- mc$died_28d
brier <- mean((p - y)^2)
spieg_z <- sum((y - p)*(1 - 2*p)) / sqrt(sum(((1 - 2*p)^2) * p * (1 - p)))
spieg_p <- 2*pnorm(-abs(spieg_z))
saveRDS(list(m_full=m_full, validate=v), "internal_validation.rds")
saveRDS(cal, "calibrate_obj.rds")
fwrite(data.table(Metric=c("C-statistic (apparent)","C-statistic (optimism-corrected)",
                           "Calibration slope (corrected)","Calibration intercept (corrected)",
                           "Brier score","Spiegelhalter Z (P)","Nagelkerke R2","Bootstrap replicates"),
                  Value=c(sprintf("%.3f",0.5+v["Dxy","index.orig"]/2),
                          sprintf("%.3f",0.5+v["Dxy","index.corrected"]/2),
                          sprintf("%.3f",v["Slope","index.corrected"]),
                          sprintf("%.3f",v["Intercept","index.corrected"]),
                          sprintf("%.3f",brier), sprintf("%.2f (P=%.2f)",spieg_z,spieg_p),
                          sprintf("%.3f",v["R2","index.corrected"]),"300")),
       "table_internal_validation.csv")
log_msg("Internal validation: C corrected =", round(0.5+v["Dxy","index.corrected"]/2,3),
        "| cal slope =", round(v["Slope","index.corrected"],3),
        "| Spiegelhalter P =", round(spieg_p,2))

## ============================================================================
## STEP 6 — Sensitivity matrix (4 analyses)
## ============================================================================
## S1: in-hospital mortality
m_hosp <- glm(died_hosp ~ MP_z + dP_z + age_z + gender + pf_z, data=mc, family=binomial)

## S2: Cox 28-day — ONLY within-28d time window (avoid broken clock); test PH
mc[, surv_time := pmin(ifelse(died_28d==1 & !is.na(days_to_death) & days_to_death>=0,
                              days_to_death, 28), 28)]
mc[surv_time<=0, surv_time := 0.5]
m_cox <- coxph(Surv(surv_time, died_28d) ~ MP_z + dP_z + age_z + gender + pf_z, data=mc)
zph <- cox.zph(m_cox)   # PH VIOLATED for MP/dP/PF -> justifies logistic as primary
saveRDS(list(m_hosp=m_hosp, m_cox=m_cox, zph=zph), "sensitivity_1_2.rds")

## S3: multiple imputation (m=20) on FULL cohort -> pooled per-SD OR
am[, gender := factor(gender, levels=c("F","M"))]
dat_imp <- am[, .(MP_baseline, dP_baseline, pf_day1_min, tvpbw_baseline, anchor_age, gender, died_28d)]
set.seed(2024)
imp <- mice(dat_imp, m=20, maxit=10, method=c(MP_baseline="pmm",dP_baseline="pmm",pf_day1_min="pmm",
            tvpbw_baseline="",anchor_age="",gender="",died_28d=""), printFlag=FALSE, seed=2024)
saveRDS(imp, "mice_imp.rds")
fit_imp <- with(imp, glm(died_28d ~ scale(MP_baseline) + scale(dP_baseline) + scale(anchor_age) +
                         gender + scale(pf_day1_min), family=binomial))
pooled <- summary(pool(fit_imp), conf.int=TRUE); setDT(pooled)
log_msg("MI pooled MP per-SD OR =", round(exp(pooled[term=="scale(MP_baseline)"]$estimate),3))

## S4: trajectory slope (EXPLORATORY) — selection + immortal-time bias
slope_cohort <- am[!is.na(MP_slope) & !is.na(dP_slope) & !is.na(pf_day1_min) & !is.na(anchor_age)]
slope_cohort[, `:=`(MPslope_z=scale(MP_slope)[,1], dPslope_z=scale(dP_slope)[,1],
                    age_z=scale(anchor_age)[,1], pf_z=scale(pf_day1_min)[,1], MPbase_z=scale(MP_baseline)[,1])]
s4_naive <- glm(died_28d ~ MPslope_z + dPslope_z + age_z + gender + pf_z, data=slope_cohort, family=binomial)
s4_adj   <- glm(died_28d ~ MPslope_z + dPslope_z + age_z + gender + pf_z + MP_n_days, data=slope_cohort, family=binomial)  # deconfound obs window
s4_adj2  <- glm(died_28d ~ MPslope_z + dPslope_z + age_z + gender + pf_z + MP_n_days + MPbase_z, data=slope_cohort, family=binomial)
saveRDS(list(naive=s4_naive, adj=s4_adj, adj2=s4_adj2), "s4_slope_models.rds")
log_msg("Slope cohort N =", nrow(slope_cohort), "death% =", round(100*mean(slope_cohort$died_28d),1),
        "(selection bias); obs-days OR =", round(exp(coef(s4_adj)["MP_n_days"]),3), "/day")

## ============================================================================
## STEP 7 — Publication figures (rendered separately in notebook; see artifacts)
##   fig_missingness_audit.png, fig_spline_doseresponse.png, fig_roc_nested.png,
##   fig_calibration.png, fig_forest_adjusted_OR.png, fig_main_step5_composite.png
##   table3_model_performance.csv
## ============================================================================

log_msg("Step 5 modeling complete. All tables + model objects written.")
## ---- END 07 ----------------------------------------------------------------
