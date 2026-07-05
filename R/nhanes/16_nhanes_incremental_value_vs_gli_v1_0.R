#!/usr/bin/env Rscript
# =============================================================================
# 16_nhanes_incremental_value_vs_gli_v1_0.R
# Incremental prognostic value of the residualized PEF marker
# (resp_vulnerability_z) BEYOND externally-referenced GLI lung function.
# NHANES 2007-2012 linked mortality (all-cause death).
# -----------------------------------------------------------------------------
# WHY THIS SCRIPT EXISTS (novelty defense):
#   Reviewers will ask whether the residualized-PEF marker is anything more than
#   PEF %predicted / an FEV1 z-score renamed. This answers with a pre-specified
#   nested-model comparison on ONE common analytic sample:
#
#     M0  Base .................. age + sex + race + bmi + smoking
#     M1  Base + raw PEF (per SD)
#     M2  Base + GLI-2022 FEV1 z ................ single established reference
#     M3  Base + resp_vulnerability_z .......... the residualized PEF marker
#     M4  Base + GLI FEV1 z + PEF marker ........ marker on top of FEV1 z
#     M5  Base + GLI FEV1 z + FVC z + FEV1/FVC z  full external spirometry panel
#     M6  Base + full GLI panel + PEF marker .... DECISIVE: marker on top of the
#                                                 complete externally-referenced
#                                                 spirometry set
#
#   Primary (survey-valid) evidence : svycoxph Wald HR + p for the marker term
#                                     in M4 and M6.
#   Discrimination                  : weighted Harrell's C and paired Delta-C
#                                     (survival::concordance + coef()/vcov()).
#
# OUTPUTS (AGGREGATE STATISTICS ONLY — no row-level data written):
#   results/tables/incremental_value_models_v1_0.csv
#   results/tables/incremental_value_deltaC_v1_0.csv
#   results/figures/incremental_value_cindex_v1_0.{png,pdf}
#   results/logs/incremental_value_session_v1_0.txt
# =============================================================================

suppressPackageStartupMessages({
  library(survey); library(survival); library(dplyr); library(readr)
  library(tibble); library(ggplot2)
})
set.seed(20260704)
options(survey.lonely.psu = "adjust")
VER <- "v1_0"

# ---- paths ------------------------------------------------------------------
ROOT <- getwd()
if (!dir.exists(file.path(ROOT, "derived_sensitive"))) {
  a <- commandArgs(trailingOnly = FALSE); fa <- grep("^--file=", a, value = TRUE)
  if (length(fa)) ROOT <- normalizePath(file.path(dirname(sub("^--file=", "", fa[1])), "..", ".."), mustWork = FALSE)
}
stopifnot(dir.exists(file.path(ROOT, "derived_sensitive")))
IN_RDS  <- file.path(ROOT, "derived_sensitive/nhanes/nhanes_mortality_analysis_ready_v30_0.rds")
DIR_TAB <- file.path(ROOT, "results/tables")
DIR_FIG <- file.path(ROOT, "results/figures")
DIR_LOG <- file.path(ROOT, "results/logs")
for (d in c(DIR_TAB, DIR_FIG, DIR_LOG)) if (!dir.exists(d)) dir.create(d, recursive = TRUE)

# ---- load + common analytic sample ------------------------------------------
d0 <- readRDS(IN_RDS)
req <- c("resp_vulnerability_z","pef_l_min",
         "gli_global2022_z_fev1","gli_global2022_z_fvc","gli_global2022_z_fev1_fvc",
         "age_years","bmi")
dd <- d0 %>%
  mutate(
    sex_f = factor(sex),
    race_f = factor(race_ethnicity),
    smoking_f = factor(smoking_status, levels = c("never","former","current"))
  ) %>%
  filter(!is.na(all_cause_death), !is.na(followup_years_exam),
         !is.na(wtmec6yr), wtmec6yr > 0, followup_years_exam > 0,
         !is.na(sex_f), !is.na(race_f), !is.na(smoking_f)) %>%
  filter(if_all(all_of(req), ~ !is.na(.x)))
dd$pef_l_min_z <- as.numeric(scale(dd$pef_l_min))
# Weights rescaled to mean 1 (sum = n) for discrimination variance estimation.
# Point estimates of C are invariant to weight scaling; this only fixes the
# variance scale so paired Delta-C SEs are valid (population-scale weights, which
# sum to ~90 million, otherwise inflate concordance robust variance ~300-fold).
dd$w_scaled <- dd$wtmec6yr / mean(dd$wtmec6yr)

n_used <- nrow(dd); n_events <- sum(dd$all_cause_death == 1L)
message(sprintf("Common analytic sample: N=%d, events=%d, person-years=%.0f",
                n_used, n_events, sum(dd$followup_years_exam)))

# ---- model ladder -----------------------------------------------------------
BASE <- "age_years + sex_f + race_f + bmi + smoking_f"
GLIP <- "gli_global2022_z_fev1 + gli_global2022_z_fvc + gli_global2022_z_fev1_fvc"
specs <- tibble::tribble(
  ~id,  ~label,                                    ~rhs,
  "M0", "Base covariates",                          BASE,
  "M1", "Base + raw PEF (per SD)",                  paste(BASE, "+ pef_l_min_z"),
  "M2", "Base + GLI FEV1 z",                        paste(BASE, "+ gli_global2022_z_fev1"),
  "M3", "Base + residualized PEF marker",           paste(BASE, "+ resp_vulnerability_z"),
  "M4", "Base + GLI FEV1 z + PEF marker",           paste(BASE, "+ gli_global2022_z_fev1 + resp_vulnerability_z"),
  "M5", "Base + full GLI panel",                    paste(BASE, "+", GLIP),
  "M6", "Base + full GLI panel + PEF marker",       paste(BASE, "+", GLIP, "+ resp_vulnerability_z")
)

des <- survey::svydesign(ids = ~psu, strata = ~strata, weights = ~wtmec6yr, nest = TRUE, data = dd)
LHS <- "survival::Surv(followup_years_exam, all_cause_death)"
fit_svy <- function(rhs) survey::svycoxph(as.formula(paste(LHS, "~", rhs)), design = des)
fit_wtd <- function(rhs) survival::coxph(as.formula(paste(LHS, "~", rhs)),
                                         data = dd, weights = w_scaled, robust = TRUE, x = TRUE, y = TRUE)
svy_fits <- lapply(specs$rhs, fit_svy); names(svy_fits) <- specs$id
wtd_fits <- lapply(specs$rhs, fit_wtd); names(wtd_fits) <- specs$id

# ---- weighted Harrell's C per model (coef()/vcov() on concordance) ----------
c_one <- function(fit) {
  ct <- survival::concordance(fit)
  c(C = unname(coef(ct)), se = unname(sqrt(as.numeric(vcov(ct)))))
}
cidx <- t(vapply(wtd_fits, c_one, numeric(2))) |> as.data.frame() |>
  tibble::rownames_to_column("id")

# ---- survey-weighted marker/FEV1 HRs ----------------------------------------
grab <- function(fit, term) {
  s <- summary(fit)$coefficients
  if (!term %in% rownames(s)) return(c(hr=NA,lo=NA,hi=NA,p=NA))
  b <- s[term,"coef"]; se <- s[term, if ("robust se" %in% colnames(s)) "robust se" else "se(coef)"]
  c(hr=exp(b), lo=exp(b-1.96*se), hi=exp(b+1.96*se), p=s[term, ncol(s)])
}
hr_marker_M3 <- grab(svy_fits[["M3"]], "resp_vulnerability_z")
hr_marker_M4 <- grab(svy_fits[["M4"]], "resp_vulnerability_z")
hr_marker_M6 <- grab(svy_fits[["M6"]], "resp_vulnerability_z")
hr_fev1_M4   <- grab(svy_fits[["M4"]], "gli_global2022_z_fev1")

# ---- model-level table ------------------------------------------------------
models_tbl <- specs |>
  left_join(cidx, by = "id") |>
  mutate(n = n_used, events = n_events,
         C_lo = C - 1.96*se, C_hi = C + 1.96*se,
         C_ci = sprintf("%.3f (%.3f-%.3f)", C, C_lo, C_hi)) |>
  select(id, label, n, events, C, C_se = se, C_lo, C_hi, C_ci)

# ---- paired Delta-C (b - a) with covariance via vcov() ----------------------
deltaC <- function(a, b, lab) {
  ct <- survival::concordance(wtd_fits[[a]], wtd_fits[[b]])
  cc <- coef(ct); V <- vcov(ct); w <- c(-1, 1)
  d <- as.numeric(w %*% cc); se <- sqrt(as.numeric(t(w) %*% V %*% w)); z <- d/se
  tibble::tibble(contrast = lab, model_a = a, model_b = b,
                 C_a = unname(cc[1]), C_b = unname(cc[2]),
                 deltaC = d, se = se, lo = d-1.96*se, hi = d+1.96*se,
                 z = z, p = 2*pnorm(-abs(z)))
}
delta_tbl <- dplyr::bind_rows(
  deltaC("M0","M3","PEF marker vs base"),
  deltaC("M0","M2","GLI FEV1 z vs base"),
  deltaC("M1","M3","Residualized marker vs raw PEF"),
  deltaC("M2","M4","PEF marker on top of GLI FEV1 z"),
  deltaC("M5","M6","PEF marker on top of FULL GLI panel [DECISIVE]")
)

# ---- write aggregate outputs ------------------------------------------------
readr::write_csv(models_tbl, file.path(DIR_TAB, paste0("incremental_value_models_", VER, ".csv")))
readr::write_csv(delta_tbl,  file.path(DIR_TAB, paste0("incremental_value_deltaC_", VER, ".csv")))

# ---- figure -----------------------------------------------------------------
plot_df <- models_tbl |>
  mutate(id = factor(id, levels = rev(specs$id)))
p <- ggplot(plot_df, aes(x = C, y = id)) +
  geom_errorbarh(aes(xmin = C_lo, xmax = C_hi), height = 0.16, colour = "grey40") +
  geom_point(size = 2.6, colour = "#1f4e79") +
  geom_text(aes(label = label), hjust = 0, nudge_y = 0.30, size = 2.9, colour = "grey20") +
  scale_y_discrete(expand = expansion(add = c(0.6, 0.9))) +
  labs(x = "Weighted Harrell's C (95% CI)", y = NULL,
       title = "Incremental discrimination for all-cause mortality",
       subtitle = sprintf("NHANES 2007-2012 linked mortality; common sample N=%d, deaths=%d", n_used, n_events)) +
  theme_minimal(base_size = 11) +
  theme(plot.title.position = "plot", panel.grid.minor = element_blank())
ggsave(file.path(DIR_FIG, paste0("incremental_value_cindex_", VER, ".png")), plot = p,
       width = 7.6, height = 4.2, dpi = 600, units = "in")
ggsave(file.path(DIR_FIG, paste0("incremental_value_cindex_", VER, ".pdf")), plot = p,
       width = 7.6, height = 4.2, units = "in", device = grDevices::pdf, useDingbats = FALSE)

# ---- console summary --------------------------------------------------------
cat("\n================ INCREMENTAL VALUE — KEY RESULTS ================\n")
cat(sprintf("Common sample: N=%d, deaths=%d\n\n", n_used, n_events))
print(as.data.frame(models_tbl[, c("id","label","C_ci")]), row.names = FALSE)
cat("\n-- Paired Delta-C (model_b - model_a) --\n")
print(as.data.frame(delta_tbl[, c("contrast","C_a","C_b","deltaC","lo","hi","p")]), row.names = FALSE)
cat(sprintf("\nMarker HR (survey-weighted):\n  M3 alone : %.2f (%.2f-%.2f) p=%.3g\n  M4 |FEV1z: %.2f (%.2f-%.2f) p=%.3g\n  M6 |panel: %.2f (%.2f-%.2f) p=%.3g\n",
            hr_marker_M3["hr"],hr_marker_M3["lo"],hr_marker_M3["hi"],hr_marker_M3["p"],
            hr_marker_M4["hr"],hr_marker_M4["lo"],hr_marker_M4["hi"],hr_marker_M4["p"],
            hr_marker_M6["hr"],hr_marker_M6["lo"],hr_marker_M6["hi"],hr_marker_M6["p"]))
cat(sprintf("  (ref) GLI FEV1 z HR in M4: %.2f (%.2f-%.2f) p=%.3g\n",
            hr_fev1_M4["hr"],hr_fev1_M4["lo"],hr_fev1_M4["hi"],hr_fev1_M4["p"]))
cat("================================================================\n")
writeLines(capture.output(sessionInfo()), file.path(DIR_LOG, paste0("incremental_value_session_", VER, ".txt")))
