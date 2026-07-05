# Survey-weighted NHANES model templates.
# These functions expect an analysis-ready dataset with internal standardized names.

required_design_columns <- c(
  "participant_id",
  "cycle",
  "survey_weight",
  "psu",
  "strata"
)

assert_columns_present <- function(data, columns, label = "data") {
  missing_columns <- setdiff(columns, names(data))
  if (length(missing_columns) > 0) {
    stop(
      label, " is missing required columns: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

assert_no_codebook_placeholders <- function(data) {
  if (any(data == "TODO_CODEBOOK_CHECK", na.rm = TRUE)) {
    stop("Analysis-ready NHANES data still contains TODO_CODEBOOK_CHECK.", call. = FALSE)
  }
  invisible(TRUE)
}

require_survey <- function() {
  if (!requireNamespace("survey", quietly = TRUE)) {
    stop("Install the R package 'survey' before running NHANES survey-weighted models.", call. = FALSE)
  }
  invisible(TRUE)
}

build_nhanes_design <- function(data, weight_var = "survey_weight", psu_var = "psu", strata_var = "strata") {
  require_survey()
  assert_columns_present(data, c(required_design_columns, weight_var, psu_var, strata_var), "NHANES analysis data")
  assert_no_codebook_placeholders(data)

  data <- data[!is.na(data[[weight_var]]) & data[[weight_var]] > 0, , drop = FALSE]
  if (nrow(data) == 0) {
    stop("No records with positive survey weights remain.", call. = FALSE)
  }

  survey::svydesign(
    ids = stats::as.formula(paste0("~", psu_var)),
    strata = stats::as.formula(paste0("~", strata_var)),
    weights = stats::as.formula(paste0("~", weight_var)),
    nest = TRUE,
    data = data
  )
}

fit_survey_glm <- function(design, outcome, exposure, covariates = character(), family = stats::gaussian()) {
  require_survey()
  terms <- c(exposure, covariates)
  formula <- stats::as.formula(paste(outcome, "~", paste(terms, collapse = " + ")))
  survey::svyglm(formula, design = design, family = family)
}

fit_survey_logistic <- function(design, outcome, exposure, covariates = character()) {
  fit_survey_glm(
    design = design,
    outcome = outcome,
    exposure = exposure,
    covariates = covariates,
    family = stats::quasibinomial()
  )
}

extract_model_row <- function(model, analysis_id, outcome, exposure, model_type, population, cycles, weight_variable) {
  coefficients <- summary(model)$coefficients
  if (!exposure %in% rownames(coefficients)) {
    stop("Exposure coefficient not found in model summary: ", exposure, call. = FALSE)
  }

  estimate <- coefficients[exposure, "Estimate"]
  se <- coefficients[exposure, "Std. Error"]
  p_value <- coefficients[exposure, ncol(coefficients)]
  conf_low <- estimate - stats::qnorm(0.975) * se
  conf_high <- estimate + stats::qnorm(0.975) * se

  data.frame(
    analysis_id = analysis_id,
    outcome = outcome,
    exposure = exposure,
    model_type = model_type,
    population = population,
    cycles = cycles,
    weight_variable = weight_variable,
    psu_variable = "psu",
    strata_variable = "strata",
    n = stats::nobs(model),
    events = NA,
    estimate = estimate,
    conf_low = conf_low,
    conf_high = conf_high,
    p_value = p_value,
    notes = "Template output; verify outcome scale before interpretation.",
    stringsAsFactors = FALSE
  )
}

run_template_model <- function(data, outcome, exposure, covariates = c("age", "sex")) {
  design <- build_nhanes_design(data)
  model <- fit_survey_logistic(design, outcome = outcome, exposure = exposure, covariates = covariates)
  extract_model_row(
    model = model,
    analysis_id = "TODO_ANALYSIS_ID",
    outcome = outcome,
    exposure = exposure,
    model_type = "survey_weighted_quasibinomial",
    population = "TODO_CODEBOOK_CHECK",
    cycles = "TODO_CODEBOOK_CHECK",
    weight_variable = "survey_weight"
  )
}

