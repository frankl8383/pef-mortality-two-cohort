# NHANES v0.3 survey-weighted replication models.
# Uses the local-only v0.2 NHANES analysis dataset and writes aggregate model
# summaries only.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(survey)
  library(tibble)
})

options(survey.lonely.psu = "adjust")

find_project_root <- function() {
  script_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", script_args, value = TRUE)
  if (length(file_arg) > 0) {
    script_path <- sub("^--file=", "", file_arg[[1]])
    return(normalizePath(file.path(dirname(script_path), "..", ".."), mustWork = FALSE))
  }
  normalizePath(getwd(), mustWork = FALSE)
}

format_p <- function(p) {
  ifelse(is.na(p), NA_character_, ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
}

markdown_table <- function(dat) {
  dat <- as.data.frame(dat, stringsAsFactors = FALSE)
  dat[] <- lapply(dat, function(x) ifelse(is.na(x), "", as.character(x)))
  header <- paste0("| ", paste(names(dat), collapse = " | "), " |")
  separator <- paste0("| ", paste(rep("---", ncol(dat)), collapse = " | "), " |")
  rows <- apply(dat, 1, function(x) paste0("| ", paste(x, collapse = " | "), " |"))
  c(header, separator, rows)
}

model_specs <- tibble::tribble(
  ~analysis_id, ~outcome, ~exposure, ~label, ~model_family, ~quality_scope,
  "nhanes_v0_3_obstruction_abc", "obstruction_fixed_ratio_abc", "resp_vulnerability_z", "Spirometric obstruction, A/B/C quality", "survey_quasibinomial", "main A/B/C",
  "nhanes_v0_3_obstruction_ab_strict", "obstruction_fixed_ratio_ab", "resp_vulnerability_z", "Spirometric obstruction, strict A/B quality", "survey_quasibinomial", "strict A/B",
  "nhanes_v0_3_current_asthma", "current_asthma", "resp_vulnerability_z", "Current asthma", "survey_quasibinomial", "PEF A/B/C exposure",
  "nhanes_v0_3_emphysema_or_bronchitis", "self_reported_emphysema_or_bronchitis", "resp_vulnerability_z", "Self-reported emphysema or chronic bronchitis", "survey_quasibinomial", "PEF A/B/C exposure",
  "nhanes_v0_3_frailty_proxy_ge2", "nhanes_frailty_proxy_ge2", "resp_vulnerability_z", "NHANES frailty proxy >=2 components", "survey_quasibinomial", "PEF A/B/C exposure",
  "nhanes_v0_3_obstruction_race_calibrated", "obstruction_fixed_ratio_abc", "resp_vulnerability_z_race_calibrated", "Spirometric obstruction, race-calibrated exposure", "survey_quasibinomial", "race-calibrated sensitivity"
)

prepare_model_data <- function(data, outcome, exposure) {
  covariates <- c("age_years", "sex", "race_ethnicity", "bmi", "smoking_status")
  required <- c(outcome, exposure, covariates, "wtmec6yr", "psu", "strata", "adult_45plus", "cycle")
  missing_cols <- setdiff(required, names(data))
  if (length(missing_cols) > 0) {
    stop("NHANES model data missing columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }
  data %>%
    filter(.data$adult_45plus == 1) %>%
    filter(!is.na(.data[[outcome]]), !is.na(.data[[exposure]])) %>%
    filter(!is.na(.data$age_years), !is.na(.data$sex), !is.na(.data$race_ethnicity), !is.na(.data$bmi), !is.na(.data$smoking_status)) %>%
    filter(!is.na(.data$wtmec6yr), .data$wtmec6yr > 0, !is.na(.data$psu), !is.na(.data$strata)) %>%
    mutate(
      sex = factor(.data$sex),
      race_ethnicity = factor(.data$race_ethnicity),
      smoking_status = factor(.data$smoking_status, levels = c("never", "former", "current")),
      cycle = factor(.data$cycle)
    )
}

fit_spec <- function(data, spec) {
  df <- prepare_model_data(data, spec$outcome, spec$exposure)
  if (nrow(df) < 100 || sum(df[[spec$outcome]] == 1, na.rm = TRUE) < 20) {
    stop("Insufficient records/events for model: ", spec$analysis_id, call. = FALSE)
  }
  design <- survey::svydesign(
    ids = ~psu,
    strata = ~strata,
    weights = ~wtmec6yr,
    nest = TRUE,
    data = df
  )
  formula <- stats::as.formula(paste(
    spec$outcome,
    "~",
    paste(c(spec$exposure, "age_years", "sex", "race_ethnicity", "bmi", "smoking_status"), collapse = " + ")
  ))
  warnings <- character()
  model <- withCallingHandlers(
    survey::svyglm(formula, design = design, family = quasibinomial()),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  coefficients <- summary(model)$coefficients
  if (!spec$exposure %in% rownames(coefficients)) {
    stop("Exposure coefficient not found for model: ", spec$analysis_id, call. = FALSE)
  }
  estimate <- coefficients[spec$exposure, "Estimate"]
  se <- coefficients[spec$exposure, "Std. Error"]
  p_value <- coefficients[spec$exposure, ncol(coefficients)]
  conf_low <- estimate - stats::qnorm(0.975) * se
  conf_high <- estimate + stats::qnorm(0.975) * se
  list(
    row = tibble(
      analysis_id = spec$analysis_id,
      outcome = spec$outcome,
      outcome_label = spec$label,
      exposure = spec$exposure,
      quality_scope = spec$quality_scope,
      model_type = spec$model_family,
      population = "NHANES 2007-2012 adults age >=45 with valid exposure/outcome/covariates",
      weight_variable = "wtmec6yr",
      psu_variable = "psu",
      strata_variable = "strata",
      n = nrow(df),
      events = sum(df[[spec$outcome]] == 1, na.rm = TRUE),
      estimate_log_or = estimate,
      se_log_or = se,
      or = exp(estimate),
      conf_low = exp(conf_low),
      conf_high = exp(conf_high),
      p_value = p_value,
      p_value_formatted = format_p(p_value),
      warnings_n = length(warnings)
    ),
    warnings = tibble(
      analysis_id = spec$analysis_id,
      warning = if (length(warnings) == 0) NA_character_ else warnings
    )
  )
}

main <- function() {
  root <- find_project_root()
  rds_path <- file.path(root, "derived_sensitive", "nhanes", "nhanes_replication_v0_2_analysis_ready.rds")
  if (!file.exists(rds_path)) {
    stop("Missing NHANES V0.2 analysis-ready RDS. Run R/nhanes/04_clean_nhanes_replication_v0_2.R first.", call. = FALSE)
  }
  table_dir <- file.path(root, "results", "tables")
  log_dir <- file.path(root, "results", "logs")
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

  data <- readRDS(rds_path)
  fits <- lapply(seq_len(nrow(model_specs)), function(i) fit_spec(data, model_specs[i, ]))
  model_table <- bind_rows(lapply(fits, `[[`, "row")) %>%
    mutate(
      or_ci = sprintf("%.2f (%.2f-%.2f)", .data$or, .data$conf_low, .data$conf_high)
    )
  warning_table <- bind_rows(lapply(fits, `[[`, "warnings")) %>%
    filter(!is.na(.data$warning))
  if (nrow(warning_table) == 0) {
    warning_table <- tibble(analysis_id = character(), warning = character())
  }

  key_terms <- model_table %>%
    transmute(
      analysis_id,
      outcome_label,
      exposure,
      n,
      events,
      or_per_1sd_higher_respiratory_vulnerability = round(.data$or, 3),
      ci = sprintf("%.3f-%.3f", .data$conf_low, .data$conf_high),
      p_value = .data$p_value_formatted
    )

  readr::write_csv(model_table, file.path(table_dir, "nhanes_v0_3_survey_model_table.csv"))
  readr::write_csv(key_terms, file.path(table_dir, "nhanes_v0_3_survey_key_terms.csv"))
  readr::write_csv(warning_table, file.path(table_dir, "nhanes_v0_3_survey_model_warnings.csv"))

  log_lines <- c(
    "# NHANES V0.3 Survey Replication Models",
    "",
    paste0("- Run date: ", Sys.Date()),
    "- Dataset: local-only `derived_sensitive/nhanes/nhanes_replication_v0_2_analysis_ready.rds`.",
    "- Design: `svydesign(ids = ~psu, strata = ~strata, weights = ~wtmec6yr, nest = TRUE)`.",
    "- Covariates: age, sex, race/ethnicity, BMI, smoking status.",
    "- Coefficient scale: odds ratio per 1-SD higher residualized respiratory vulnerability.",
    "",
    "## Key Results",
    "",
    markdown_table(key_terms),
    "",
    "## Warnings",
    "",
    if (nrow(warning_table) == 0) {
      "No survey model warnings."
    } else {
      markdown_table(warning_table)
    },
    "",
    "## Interpretation Boundary",
    "",
    "- These are first-pass cross-sectional replication models, not final manuscript models.",
    "- Spirometric obstruction is an objective NHANES outcome; self-reported respiratory diseases are secondary proxies.",
    "- Frailty is a proxy outcome/covariate layer, not a validated NHANES frailty index."
  )
  log_path <- file.path(log_dir, "nhanes_v0_3_survey_model_log.md")
  writeLines(log_lines, log_path)

  message("Wrote NHANES v0.3 survey model table: ", file.path(table_dir, "nhanes_v0_3_survey_model_table.csv"))
  message("Wrote log: ", log_path)
}

if (sys.nframe() == 0) {
  main()
}
