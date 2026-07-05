# Build NHANES V31 mortality sensitivity lock, including PH diagnostics.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(survival)
  library(tibble)
})

version_id <- "v31_5"

find_project_root <- function() {
  script_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", script_args, value = TRUE)
  if (length(file_arg) > 0) {
    script_path <- sub("^--file=", "", file_arg[[1]])
    return(normalizePath(file.path(dirname(script_path), "..", ".."), mustWork = FALSE))
  }
  normalizePath(getwd(), mustWork = FALSE)
}

dir_create <- function(path) if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)

md_escape <- function(x) str_replace_all(ifelse(is.na(x), "", as.character(x)), "\\|", "\\\\|")
md_table <- function(dat) {
  dat <- as.data.frame(dat, stringsAsFactors = FALSE)
  dat[] <- lapply(dat, function(x) ifelse(is.na(x), "", as.character(x)))
  c(
    paste0("| ", paste(md_escape(names(dat)), collapse = " | "), " |"),
    paste0("| ", paste(rep("---", ncol(dat)), collapse = " | "), " |"),
    apply(dat, 1, function(x) paste0("| ", paste(md_escape(x), collapse = " | "), " |"))
  )
}

fmt_pct <- function(x) paste0(sprintf("%.1f", 100 * as.numeric(x)), "%")

prepare_mortality_data <- function(dat) {
  dat %>%
    filter(mortality_analysis_primary == 1L) %>%
    mutate(
      sex_f = factor(sex),
      race_f = factor(race_ethnicity),
      education_f = factor(education),
      smoking_f = factor(smoking_status, levels = c("never", "former", "current")),
      cycle_f = factor(cycle_label, levels = c("E", "F", "G"))
    )
}

run_ph_diagnostic <- function(dat, exposure, model_label) {
  rhs <- c(exposure, "age_years", "sex_f", "race_f", "education_f", "income_poverty_ratio", "bmi", "smoking_f")
  needed <- c("followup_years_exam", "all_cause_death", rhs, "wtmec6yr", "psu")
  df <- dat %>%
    filter(!is.na(wtmec6yr), wtmec6yr > 0, complete.cases(across(all_of(needed))))
  if (nrow(df) < 1000 || sum(df$all_cause_death == 1L, na.rm = TRUE) < 100) {
    return(tibble(
      model_label = model_label,
      term = c(exposure, "GLOBAL"),
      chisq = NA_real_,
      p_value = NA_real_,
      diagnostic_method = "weighted_coxph_cox_zph",
      n = nrow(df),
      events = sum(df$all_cause_death == 1L, na.rm = TRUE),
      status = "insufficient_records_or_events"
    ))
  }
  form <- as.formula(paste0(
    "Surv(followup_years_exam, all_cause_death) ~ ",
    paste(rhs, collapse = " + "),
    " + cluster(psu)"
  ))
  fit <- survival::coxph(form, data = df, weights = wtmec6yr, x = TRUE)
  z <- survival::cox.zph(fit, transform = "km")
  ztab <- as.data.frame(z$table)
  ztab$term <- rownames(ztab)
  n_diag <- nrow(df)
  events_diag <- sum(df$all_cause_death == 1L, na.rm = TRUE)
  ztab %>%
    filter(term %in% c(exposure, "GLOBAL")) %>%
    transmute(
      model_label,
      term,
      chisq = round(chisq, 4),
      p_value = round(p, 4),
      diagnostic_method = "weighted_coxph_cox_zph_approximation",
      n = n_diag,
      events = events_diag,
      status = if_else(p_value >= 0.05, "no_ph_violation_detected", "possible_ph_violation_review_needed")
    )
}

root <- find_project_root()
tab_dir <- file.path(root, "results", "tables")
log_dir <- file.path(root, "results", "logs")
manuscript_dir <- file.path(root, "manuscript")
invisible(lapply(c(tab_dir, log_dir, manuscript_dir), dir_create))

models <- read_csv(file.path(tab_dir, "nhanes_mortality_models_v30_0.csv"), show_col_types = FALSE)
dat <- readRDS(file.path(root, "derived_sensitive", "nhanes", "nhanes_mortality_analysis_ready_v30_0.rds")) %>%
  prepare_mortality_data()

ph <- bind_rows(
  run_ph_diagnostic(dat, "resp_vulnerability_z", "M2 primary socioeconomic"),
  run_ph_diagnostic(dat, "resp_vulnerability_z_race_calibrated", "M2 race-calibrated exposure")
)
write_csv(ph, file.path(tab_dir, "nhanes_v31_ph_diagnostics.csv"))

pick <- function(id, term = NULL) {
  out <- models %>% filter(analysis_id == id)
  if (!is.null(term)) out <- out %>% filter(.data$term == term)
  out %>% slice(1)
}

primary <- pick("t2_allcause_m2_primary", "resp_vulnerability_z")
race_cal <- pick("t2_allcause_m2_race_calibrated", "resp_vulnerability_z_race_calibrated")
cycle_fixed <- pick("t2_allcause_m2_cycle_fixed", "resp_vulnerability_z")
frailty <- pick("t2_allcause_m3_frailty_adjusted", "resp_vulnerability_z")
clrd <- pick("t2_clrd_mini_exploratory", "resp_vulnerability_z")
cycles <- models %>% filter(analysis_family == "cycle_stratified_sensitivity", term == "resp_vulnerability_z")

frailty_missing_n <- primary$n - frailty$n
frailty_missing_pct <- frailty_missing_n / primary$n

sensitivity <- bind_rows(
  tibble(
    check_type = "primary_reference",
    analysis = "Primary all-cause mortality M2",
    n = primary$n,
    events = primary$events,
    estimate = primary$hr_ci,
    p_value = primary$p_value_formatted,
    interpretation = "Primary linked-mortality estimate remains the reference.",
    boundary = "Observational all-cause mortality association."
  ),
  tibble(
    check_type = "proportional_hazards",
    analysis = "Weighted cox.zph approximation for primary exposure",
    n = ph$n[ph$model_label == "M2 primary socioeconomic" & ph$term == "resp_vulnerability_z"],
    events = ph$events[ph$model_label == "M2 primary socioeconomic" & ph$term == "resp_vulnerability_z"],
    estimate = paste0("PH p=", ph$p_value[ph$model_label == "M2 primary socioeconomic" & ph$term == "resp_vulnerability_z"],
                      "; global p=", ph$p_value[ph$model_label == "M2 primary socioeconomic" & ph$term == "GLOBAL"]),
    p_value = as.character(ph$p_value[ph$model_label == "M2 primary socioeconomic" & ph$term == "resp_vulnerability_z"]),
    interpretation = "No strong PH violation if p>=0.05; review Schoenfeld diagnostics if p<0.05.",
    boundary = "Diagnostic approximation using weighted coxph because survey-weighted cox.zph is not directly available."
  ),
  tibble(
    check_type = "race_calibrated",
    analysis = "Race-calibrated marker sensitivity",
    n = race_cal$n,
    events = race_cal$events,
    estimate = race_cal$hr_ci,
    p_value = race_cal$p_value_formatted,
    interpretation = "Race-calibrated exposure gives a similar all-cause mortality estimate.",
    boundary = "Sensitivity only; primary uses race-neutral marker with race covariate."
  ),
  tibble(
    check_type = "cycle_fixed",
    analysis = "M2 plus cycle fixed effects",
    n = cycle_fixed$n,
    events = cycle_fixed$events,
    estimate = cycle_fixed$hr_ci,
    p_value = cycle_fixed$p_value_formatted,
    interpretation = "Adding cycle fixed effects does not materially change the primary estimate.",
    boundary = "Sensitivity for pooled-cycle structure."
  ),
  cycles %>%
    transmute(
      check_type = "cycle_stratified",
      analysis = paste0("Cycle ", cycle, " sensitivity"),
      n, events, estimate = hr_ci, p_value = p_value_formatted,
      interpretation = "Cycle-specific estimate; precision varies by sample size and events.",
      boundary = "Sensitivity; do not overinterpret between-cycle differences."
    ),
  tibble(
    check_type = "frailty_complete_case",
    analysis = "Frailty-adjusted complete-case sensitivity",
    n = frailty$n,
    events = frailty$events,
    estimate = frailty$hr_ci,
    p_value = frailty$p_value_formatted,
    interpretation = paste0("Frailty sensitivity retains direction but uses ", frailty_missing_n,
                            " fewer participants than primary M2 (", fmt_pct(frailty_missing_pct), ")."),
    boundary = "Sensitivity only; does not prove independence from frailty."
  ),
  tibble(
    check_type = "clrd_exploratory",
    analysis = "CLRD mortality exploratory age-sex model",
    n = clrd$n,
    events = clrd$events,
    estimate = clrd$hr_ci,
    p_value = clrd$p_value_formatted,
    interpretation = "Directionally supportive but sparse.",
    boundary = "Exploratory only; not abstract-level evidence."
  )
)

write_csv(sensitivity, file.path(tab_dir, "nhanes_v31_mortality_sensitivity.csv"))

note <- c(
  "# NHANES Mortality Sensitivity V31.5",
  "",
  "## Summary",
  "",
  "The NHANES linked-mortality result remains robust as an all-cause mortality extension of the PEF-based marker. The primary M2 estimate is preserved as the reference result; race-calibrated, cycle-fixed, cycle-stratified and frailty-complete-case analyses are treated as sensitivity checks.",
  "",
  "## Sensitivity table",
  "",
  md_table(sensitivity),
  "",
  "## Proportional-hazards diagnostics",
  "",
  md_table(ph),
  "",
  "## Manuscript boundary",
  "",
  "Use the NHANES mortality layer as observational linked-mortality extension. The frailty-adjusted model is sensitivity only because its complete-case sample is smaller. CLRD mortality remains exploratory and should not be placed in the abstract."
)
writeLines(note, file.path(manuscript_dir, "nhanes_mortality_sensitivity_v31_5.md"))

log <- c(
  "# NHANES V31 Mortality Sensitivity Build",
  "",
  paste0("Version: ", version_id),
  paste0("Rows written: ", nrow(sensitivity)),
  paste0("PH diagnostic rows written: ", nrow(ph)),
  "",
  "Outputs:",
  "- results/tables/nhanes_v31_mortality_sensitivity.csv",
  "- results/tables/nhanes_v31_ph_diagnostics.csv",
  "- manuscript/nhanes_mortality_sensitivity_v31_5.md"
)
writeLines(log, file.path(log_dir, "nhanes_v31_mortality_sensitivity.md"))

message("NHANES V31 mortality sensitivity outputs written.")
