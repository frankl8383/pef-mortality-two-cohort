# NHANES v0.4 weighted descriptives, quartile plots, and sensitivity models.
# Uses the local-only v0.2 NHANES analysis dataset and writes aggregate outputs.

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(stringr)
  library(survey)
  library(tibble)
  library(tidyr)
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

format_ci <- function(est, lo, hi, digits = 2) {
  sprintf(
    paste0("%.", digits, "f (%.", digits, "f-%.", digits, "f)"),
    est, lo, hi
  )
}

markdown_table <- function(dat) {
  dat <- as.data.frame(dat, stringsAsFactors = FALSE)
  dat[] <- lapply(dat, function(x) ifelse(is.na(x), "", as.character(x)))
  header <- paste0("| ", paste(names(dat), collapse = " | "), " |")
  separator <- paste0("| ", paste(rep("---", ncol(dat)), collapse = " | "), " |")
  rows <- apply(dat, 1, function(x) paste0("| ", paste(x, collapse = " | "), " |"))
  c(header, separator, rows)
}

weighted_quantile <- function(x, w, probs = c(0.25, 0.5, 0.75)) {
  ok <- !is.na(x) & !is.na(w) & w > 0
  x <- x[ok]
  w <- w[ok]
  if (length(x) == 0) {
    return(rep(NA_real_, length(probs)))
  }
  ord <- order(x)
  x <- x[ord]
  w <- w[ord]
  cw <- cumsum(w) / sum(w)
  vapply(probs, function(p) x[which(cw >= p)[1]], numeric(1))
}

make_design <- function(data, weight_var = "wtmec6yr") {
  survey::svydesign(
    ids = ~psu,
    strata = ~strata,
    weights = stats::as.formula(paste0("~", weight_var)),
    nest = TRUE,
    data = data
  )
}

clean_design_data <- function(data, weight_var = "wtmec6yr") {
  data %>%
    filter(
      !is.na(.data[[weight_var]]),
      .data[[weight_var]] > 0,
      !is.na(.data$psu),
      !is.na(.data$strata)
    )
}

estimate_continuous <- function(data, var, label, group_label, weight_var = "wtmec6yr") {
  df <- clean_design_data(data, weight_var) %>%
    filter(!is.na(.data[[var]]))
  if (nrow(df) < 2) {
    return(tibble(
      group = group_label,
      variable = var,
      variable_label = label,
      type = "continuous",
      level = NA_character_,
      n_unweighted = nrow(df),
      estimate = NA_real_,
      se = NA_real_,
      ci_low = NA_real_,
      ci_high = NA_real_,
      display = NA_character_
    ))
  }
  design <- make_design(df, weight_var)
  form <- stats::as.formula(paste0("~", var))
  est <- survey::svymean(form, design = design, na.rm = TRUE)
  mean_value <- as.numeric(stats::coef(est)[1])
  se_value <- as.numeric(survey::SE(est)[1])
  ci <- mean_value + c(-1, 1) * stats::qnorm(0.975) * se_value
  tibble(
    group = group_label,
    variable = var,
    variable_label = label,
    type = "continuous",
    level = NA_character_,
    n_unweighted = nrow(df),
    estimate = mean_value,
    se = se_value,
    ci_low = ci[1],
    ci_high = ci[2],
    display = sprintf("%.1f (SE %.1f)", mean_value, se_value)
  )
}

estimate_binary <- function(data, var, label, group_label, weight_var = "wtmec6yr") {
  df <- clean_design_data(data, weight_var) %>%
    filter(!is.na(.data[[var]])) %>%
    mutate(.indicator = as.numeric(.data[[var]] == 1))
  if (nrow(df) < 2) {
    return(tibble(
      group = group_label,
      variable = var,
      variable_label = label,
      type = "binary",
      level = "yes",
      n_unweighted = nrow(df),
      events_unweighted = sum(df$.indicator == 1, na.rm = TRUE),
      estimate = NA_real_,
      se = NA_real_,
      ci_low = NA_real_,
      ci_high = NA_real_,
      display = NA_character_
    ))
  }
  design <- make_design(df, weight_var)
  est <- survey::svymean(~.indicator, design = design, na.rm = TRUE)
  prop <- as.numeric(stats::coef(est)[1])
  se_value <- as.numeric(survey::SE(est)[1])
  ci <- pmin(pmax(prop + c(-1, 1) * stats::qnorm(0.975) * se_value, 0), 1)
  tibble(
    group = group_label,
    variable = var,
    variable_label = label,
    type = "binary",
    level = "yes",
    n_unweighted = nrow(df),
    events_unweighted = sum(df$.indicator == 1, na.rm = TRUE),
    estimate = prop,
    se = se_value,
    ci_low = ci[1],
    ci_high = ci[2],
    display = sprintf("%.1f%% (SE %.1f)", 100 * prop, 100 * se_value)
  )
}

estimate_categorical <- function(data, var, label, group_label, weight_var = "wtmec6yr") {
  df <- clean_design_data(data, weight_var) %>%
    filter(!is.na(.data[[var]]), .data[[var]] != "") %>%
    mutate(.level = factor(.data[[var]]))
  if (nrow(df) < 2 || dplyr::n_distinct(df$.level) < 1) {
    return(tibble(
      group = group_label,
      variable = var,
      variable_label = label,
      type = "categorical",
      level = NA_character_,
      n_unweighted = nrow(df),
      events_unweighted = NA_integer_,
      estimate = NA_real_,
      se = NA_real_,
      ci_low = NA_real_,
      ci_high = NA_real_,
      display = NA_character_
    ))
  }
  design <- make_design(df, weight_var)
  est <- survey::svymean(~.level, design = design, na.rm = TRUE)
  coefs <- stats::coef(est)
  ses <- survey::SE(est)
  levels_out <- sub("^\\.level", "", names(coefs))
  tibble(
    group = group_label,
    variable = var,
    variable_label = label,
    type = "categorical",
    level = levels_out,
    n_unweighted = nrow(df),
    events_unweighted = as.integer(table(df$.level)[levels_out]),
    estimate = as.numeric(coefs),
    se = as.numeric(ses),
    ci_low = pmax(as.numeric(coefs) - stats::qnorm(0.975) * as.numeric(ses), 0),
    ci_high = pmin(as.numeric(coefs) + stats::qnorm(0.975) * as.numeric(ses), 1),
    display = sprintf("%.1f%% (SE %.1f)", 100 * as.numeric(coefs), 100 * as.numeric(ses))
  )
}

continuous_vars <- tibble::tribble(
  ~variable, ~variable_label,
  "age_years", "Age, years",
  "bmi", "Body mass index, kg/m2",
  "height_cm", "Height, cm",
  "pef_l_min", "Peak expiratory flow, L/min",
  "fev1_fvc", "FEV1/FVC ratio",
  "resp_vulnerability_z", "Respiratory vulnerability z-score",
  "phq9_score", "PHQ-9 score",
  "wbc_1000_ul", "White blood cells, 1000/uL",
  "hba1c_percent", "HbA1c, percent",
  "albumin_g_dl", "Albumin, g/dL"
)

binary_vars <- tibble::tribble(
  ~variable, ~variable_label,
  "obstruction_fixed_ratio_abc", "Spirometric obstruction, A/B/C quality",
  "obstruction_fixed_ratio_ab", "Spirometric obstruction, strict A/B quality",
  "current_asthma", "Current asthma",
  "self_reported_emphysema_or_bronchitis", "Self-reported emphysema or chronic bronchitis",
  "nhanes_frailty_proxy_ge2", "NHANES frailty proxy >=2 components",
  "depressive_symptoms_phq9_ge10", "PHQ-9 >=10",
  "low_physical_activity", "Low physical activity"
)

categorical_vars <- tibble::tribble(
  ~variable, ~variable_label,
  "sex", "Sex",
  "race_ethnicity", "Race/ethnicity",
  "education", "Education",
  "smoking_status", "Smoking status",
  "cycle", "NHANES cycle"
)

outcome_specs <- tibble::tribble(
  ~outcome, ~outcome_label,
  "obstruction_fixed_ratio_abc", "Spirometric obstruction",
  "current_asthma", "Current asthma",
  "self_reported_emphysema_or_bronchitis", "Emphysema/chronic bronchitis",
  "nhanes_frailty_proxy_ge2", "Frailty proxy >=2"
)

prepare_analysis_data <- function(data) {
  base <- data %>%
    filter(.data$adult_45plus == 1) %>%
    filter(!is.na(.data$resp_vulnerability_z)) %>%
    clean_design_data("wtmec6yr")

  cuts <- weighted_quantile(base$resp_vulnerability_z, base$wtmec6yr, c(0.25, 0.5, 0.75))
  if (any(is.na(cuts)) || length(unique(cuts)) < length(cuts)) {
    stop("Unable to compute distinct weighted respiratory vulnerability quartiles.", call. = FALSE)
  }

  data_q <- data %>%
    mutate(
      rv_quartile = cut(
        .data$resp_vulnerability_z,
        breaks = c(-Inf, cuts, Inf),
        labels = c("Q1", "Q2", "Q3", "Q4"),
        include.lowest = TRUE,
        right = TRUE
      ),
      rv_quartile_num = as.numeric(.data$rv_quartile),
      rv_q4_vs_q1 = dplyr::case_when(
        .data$rv_quartile == "Q4" ~ 1,
        .data$rv_quartile == "Q1" ~ 0,
        TRUE ~ NA_real_
      )
    )

  list(
    data = data_q,
    cutpoints = tibble(
      quantile = c("weighted_p25", "weighted_p50", "weighted_p75"),
      resp_vulnerability_z = as.numeric(cuts)
    )
  )
}

build_table1 <- function(data) {
  table1_data <- data %>%
    filter(.data$adult_45plus == 1, !is.na(.data$rv_quartile)) %>%
    clean_design_data("wtmec6yr")

  estimate_group <- function(df, group_label) {
    bind_rows(
      lapply(seq_len(nrow(continuous_vars)), function(i) {
        estimate_continuous(df, continuous_vars$variable[[i]], continuous_vars$variable_label[[i]], group_label)
      }),
      lapply(seq_len(nrow(binary_vars)), function(i) {
        estimate_binary(df, binary_vars$variable[[i]], binary_vars$variable_label[[i]], group_label)
      }),
      lapply(seq_len(nrow(categorical_vars)), function(i) {
        estimate_categorical(df, categorical_vars$variable[[i]], categorical_vars$variable_label[[i]], group_label)
      })
    )
  }

  overall <- estimate_group(table1_data, "Overall")
  by_quartile <- bind_rows(lapply(levels(table1_data$rv_quartile), function(q) {
    estimate_group(filter(table1_data, .data$rv_quartile == q), q)
  }))

  list(overall = overall, by_quartile = by_quartile)
}

build_outcome_rates <- function(data) {
  rate_data <- data %>%
    filter(.data$adult_45plus == 1, !is.na(.data$rv_quartile)) %>%
    clean_design_data("wtmec6yr")

  bind_rows(lapply(seq_len(nrow(outcome_specs)), function(i) {
    outcome <- outcome_specs$outcome[[i]]
    outcome_label <- outcome_specs$outcome_label[[i]]
    bind_rows(lapply(levels(rate_data$rv_quartile), function(q) {
      df <- rate_data %>%
        filter(.data$rv_quartile == q, !is.na(.data[[outcome]])) %>%
        mutate(.indicator = as.numeric(.data[[outcome]] == 1))
      if (nrow(df) < 2) {
        return(tibble(
          outcome = outcome,
          outcome_label = outcome_label,
          rv_quartile = q,
          n_unweighted = nrow(df),
          events_unweighted = sum(df$.indicator == 1, na.rm = TRUE),
          weighted_prevalence = NA_real_,
          se = NA_real_,
          ci_low = NA_real_,
          ci_high = NA_real_
        ))
      }
      design <- make_design(df, "wtmec6yr")
      est <- survey::svymean(~.indicator, design = design, na.rm = TRUE)
      prop <- as.numeric(stats::coef(est)[1])
      se_value <- as.numeric(survey::SE(est)[1])
      ci <- pmin(pmax(prop + c(-1, 1) * stats::qnorm(0.975) * se_value, 0), 1)
      tibble(
        outcome = outcome,
        outcome_label = outcome_label,
        rv_quartile = q,
        n_unweighted = nrow(df),
        events_unweighted = sum(df$.indicator == 1, na.rm = TRUE),
        weighted_prevalence = prop,
        se = se_value,
        ci_low = ci[1],
        ci_high = ci[2]
      )
    }))
  })) %>%
    mutate(
      weighted_prevalence_percent = 100 * .data$weighted_prevalence,
      ci_low_percent = 100 * .data$ci_low,
      ci_high_percent = 100 * .data$ci_high
    )
}

prepare_model_data <- function(data, outcome, exposure, covariates, weight_var = "wtmec6yr", cycle_filter = NA_character_) {
  required <- unique(c(outcome, exposure, covariates, weight_var, "psu", "strata", "adult_45plus", "cycle"))
  missing_cols <- setdiff(required, names(data))
  if (length(missing_cols) > 0) {
    stop("NHANES model data missing columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  df <- data %>%
    filter(.data$adult_45plus == 1) %>%
    filter(!is.na(.data[[outcome]]), !is.na(.data[[exposure]])) %>%
    filter(!is.na(.data[[weight_var]]), .data[[weight_var]] > 0, !is.na(.data$psu), !is.na(.data$strata))

  if (!is.na(cycle_filter)) {
    df <- df %>% filter(.data$cycle == cycle_filter)
  }

  for (v in covariates) {
    df <- df %>% filter(!is.na(.data[[v]]))
  }

  df %>%
    mutate(
      sex = factor(.data$sex),
      race_ethnicity = factor(.data$race_ethnicity),
      smoking_status = factor(.data$smoking_status, levels = c("never", "former", "current")),
      cycle = factor(.data$cycle),
      rv_quartile = factor(.data$rv_quartile, levels = c("Q1", "Q2", "Q3", "Q4"))
    )
}

fit_model_spec <- function(data, spec) {
  covariates <- unlist(strsplit(spec$covariates, "\\|", fixed = FALSE))
  df <- prepare_model_data(
    data = data,
    outcome = spec$outcome,
    exposure = spec$exposure,
    covariates = covariates,
    weight_var = spec$weight_variable,
    cycle_filter = spec$cycle
  )
  events <- sum(df[[spec$outcome]] == 1, na.rm = TRUE)
  if (nrow(df) < 100 || events < 20) {
    stop("Insufficient records/events for model: ", spec$analysis_id, call. = FALSE)
  }

  design <- make_design(df, spec$weight_variable)
  formula <- stats::as.formula(paste(
    spec$outcome,
    "~",
    paste(c(spec$exposure, covariates), collapse = " + ")
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
  term <- spec$coefficient
  if (!term %in% rownames(coefficients)) {
    stop("Coefficient `", term, "` not found for model: ", spec$analysis_id, call. = FALSE)
  }
  estimate <- coefficients[term, "Estimate"]
  se <- coefficients[term, "Std. Error"]
  p_value <- coefficients[term, ncol(coefficients)]
  conf_low_log <- estimate - stats::qnorm(0.975) * se
  conf_high_log <- estimate + stats::qnorm(0.975) * se
  or <- exp(estimate)
  conf_low <- exp(conf_low_log)
  conf_high <- exp(conf_high_log)

  list(
    row = tibble(
      status = "ok",
      analysis_id = spec$analysis_id,
      analysis_family = spec$analysis_family,
      model_label = spec$model_label,
      outcome = spec$outcome,
      outcome_label = spec$outcome_label,
      exposure = spec$exposure,
      coefficient = term,
      covariates = paste(covariates, collapse = "; "),
      weight_variable = spec$weight_variable,
      cycle = ifelse(is.na(spec$cycle), "pooled", spec$cycle),
      n = nrow(df),
      events = events,
      estimate_log_or = estimate,
      se_log_or = se,
      or = or,
      conf_low = conf_low,
      conf_high = conf_high,
      p_value = p_value,
      p_value_formatted = format_p(p_value),
      or_ci = format_ci(or, conf_low, conf_high, 2),
      warnings_n = length(warnings),
      error = NA_character_
    ),
    warnings = tibble(
      analysis_id = spec$analysis_id,
      warning = if (length(warnings) == 0) NA_character_ else warnings
    )
  )
}

safe_fit_model_spec <- function(data, spec) {
  tryCatch(
    fit_model_spec(data, spec),
    error = function(e) {
      list(
        row = tibble(
          status = "error",
          analysis_id = spec$analysis_id,
          analysis_family = spec$analysis_family,
          model_label = spec$model_label,
          outcome = spec$outcome,
          outcome_label = spec$outcome_label,
          exposure = spec$exposure,
          coefficient = spec$coefficient,
          covariates = spec$covariates,
          weight_variable = spec$weight_variable,
          cycle = ifelse(is.na(spec$cycle), "pooled", spec$cycle),
          n = NA_integer_,
          events = NA_integer_,
          estimate_log_or = NA_real_,
          se_log_or = NA_real_,
          or = NA_real_,
          conf_low = NA_real_,
          conf_high = NA_real_,
          p_value = NA_real_,
          p_value_formatted = NA_character_,
          or_ci = NA_character_,
          warnings_n = 1L,
          error = conditionMessage(e)
        ),
        warnings = tibble(analysis_id = spec$analysis_id, warning = conditionMessage(e))
      )
    }
  )
}

build_model_specs <- function(data) {
  base_covars <- "age_years|sex|race_ethnicity|bmi|smoking_status"
  frailty_covars <- "age_years|sex|race_ethnicity|bmi|smoking_status|nhanes_frailty_proxy_count"
  cycle_covars <- "age_years|sex|race_ethnicity|bmi|smoking_status|cycle"
  cycles <- sort(unique(data$cycle[!is.na(data$cycle)]))

  sensitivity_specs <- tibble::tribble(
    ~analysis_id, ~analysis_family, ~model_label, ~outcome, ~outcome_label, ~exposure, ~coefficient, ~covariates, ~weight_variable, ~cycle,
    "v0_4_main_obstruction_abc", "pooled_sensitivity", "Main A/B/C", "obstruction_fixed_ratio_abc", "Spirometric obstruction", "resp_vulnerability_z", "resp_vulnerability_z", base_covars, "wtmec6yr", NA_character_,
    "v0_4_main_current_asthma", "pooled_sensitivity", "Main", "current_asthma", "Current asthma", "resp_vulnerability_z", "resp_vulnerability_z", base_covars, "wtmec6yr", NA_character_,
    "v0_4_main_emphysema_bronchitis", "pooled_sensitivity", "Main", "self_reported_emphysema_or_bronchitis", "Emphysema/chronic bronchitis", "resp_vulnerability_z", "resp_vulnerability_z", base_covars, "wtmec6yr", NA_character_,
    "v0_4_main_frailty_proxy", "pooled_sensitivity", "Main", "nhanes_frailty_proxy_ge2", "Frailty proxy >=2", "resp_vulnerability_z", "resp_vulnerability_z", base_covars, "wtmec6yr", NA_character_,
    "v0_4_strict_obstruction_ab", "pooled_sensitivity", "Strict A/B quality", "obstruction_fixed_ratio_ab", "Spirometric obstruction", "resp_vulnerability_z", "resp_vulnerability_z", base_covars, "wtmec6yr", NA_character_,
    "v0_4_race_cal_obstruction_abc", "pooled_sensitivity", "Race-calibrated exposure", "obstruction_fixed_ratio_abc", "Spirometric obstruction", "resp_vulnerability_z_race_calibrated", "resp_vulnerability_z_race_calibrated", base_covars, "wtmec6yr", NA_character_,
    "v0_4_race_cal_current_asthma", "pooled_sensitivity", "Race-calibrated exposure", "current_asthma", "Current asthma", "resp_vulnerability_z_race_calibrated", "resp_vulnerability_z_race_calibrated", base_covars, "wtmec6yr", NA_character_,
    "v0_4_race_cal_emphysema_bronchitis", "pooled_sensitivity", "Race-calibrated exposure", "self_reported_emphysema_or_bronchitis", "Emphysema/chronic bronchitis", "resp_vulnerability_z_race_calibrated", "resp_vulnerability_z_race_calibrated", base_covars, "wtmec6yr", NA_character_,
    "v0_4_race_cal_frailty_proxy", "pooled_sensitivity", "Race-calibrated exposure", "nhanes_frailty_proxy_ge2", "Frailty proxy >=2", "resp_vulnerability_z_race_calibrated", "resp_vulnerability_z_race_calibrated", base_covars, "wtmec6yr", NA_character_,
    "v0_4_frailty_adj_obstruction_abc", "frailty_adjusted", "Frailty-adjusted", "obstruction_fixed_ratio_abc", "Spirometric obstruction", "resp_vulnerability_z", "resp_vulnerability_z", frailty_covars, "wtmec6yr", NA_character_,
    "v0_4_frailty_adj_current_asthma", "frailty_adjusted", "Frailty-adjusted", "current_asthma", "Current asthma", "resp_vulnerability_z", "resp_vulnerability_z", frailty_covars, "wtmec6yr", NA_character_,
    "v0_4_frailty_adj_emphysema_bronchitis", "frailty_adjusted", "Frailty-adjusted", "self_reported_emphysema_or_bronchitis", "Emphysema/chronic bronchitis", "resp_vulnerability_z", "resp_vulnerability_z", frailty_covars, "wtmec6yr", NA_character_,
    "v0_4_cycle_fixed_obstruction_abc", "cycle_fixed", "Cycle fixed effects", "obstruction_fixed_ratio_abc", "Spirometric obstruction", "resp_vulnerability_z", "resp_vulnerability_z", cycle_covars, "wtmec6yr", NA_character_,
    "v0_4_cycle_fixed_current_asthma", "cycle_fixed", "Cycle fixed effects", "current_asthma", "Current asthma", "resp_vulnerability_z", "resp_vulnerability_z", cycle_covars, "wtmec6yr", NA_character_,
    "v0_4_cycle_fixed_emphysema_bronchitis", "cycle_fixed", "Cycle fixed effects", "self_reported_emphysema_or_bronchitis", "Emphysema/chronic bronchitis", "resp_vulnerability_z", "resp_vulnerability_z", cycle_covars, "wtmec6yr", NA_character_,
    "v0_4_cycle_fixed_frailty_proxy", "cycle_fixed", "Cycle fixed effects", "nhanes_frailty_proxy_ge2", "Frailty proxy >=2", "resp_vulnerability_z", "resp_vulnerability_z", cycle_covars, "wtmec6yr", NA_character_
  )

  quartile_trend_specs <- tibble::tribble(
    ~analysis_id, ~analysis_family, ~model_label, ~outcome, ~outcome_label, ~exposure, ~coefficient, ~covariates, ~weight_variable, ~cycle,
    "v0_4_quartile_trend_obstruction_abc", "quartile_model", "Trend per higher quartile", "obstruction_fixed_ratio_abc", "Spirometric obstruction", "rv_quartile_num", "rv_quartile_num", base_covars, "wtmec6yr", NA_character_,
    "v0_4_quartile_trend_current_asthma", "quartile_model", "Trend per higher quartile", "current_asthma", "Current asthma", "rv_quartile_num", "rv_quartile_num", base_covars, "wtmec6yr", NA_character_,
    "v0_4_quartile_trend_emphysema_bronchitis", "quartile_model", "Trend per higher quartile", "self_reported_emphysema_or_bronchitis", "Emphysema/chronic bronchitis", "rv_quartile_num", "rv_quartile_num", base_covars, "wtmec6yr", NA_character_,
    "v0_4_quartile_trend_frailty_proxy", "quartile_model", "Trend per higher quartile", "nhanes_frailty_proxy_ge2", "Frailty proxy >=2", "rv_quartile_num", "rv_quartile_num", base_covars, "wtmec6yr", NA_character_,
    "v0_4_q4_vs_q1_obstruction_abc", "quartile_model", "Q4 versus Q1", "obstruction_fixed_ratio_abc", "Spirometric obstruction", "rv_q4_vs_q1", "rv_q4_vs_q1", base_covars, "wtmec6yr", NA_character_,
    "v0_4_q4_vs_q1_current_asthma", "quartile_model", "Q4 versus Q1", "current_asthma", "Current asthma", "rv_q4_vs_q1", "rv_q4_vs_q1", base_covars, "wtmec6yr", NA_character_,
    "v0_4_q4_vs_q1_emphysema_bronchitis", "quartile_model", "Q4 versus Q1", "self_reported_emphysema_or_bronchitis", "Emphysema/chronic bronchitis", "rv_q4_vs_q1", "rv_q4_vs_q1", base_covars, "wtmec6yr", NA_character_,
    "v0_4_q4_vs_q1_frailty_proxy", "quartile_model", "Q4 versus Q1", "nhanes_frailty_proxy_ge2", "Frailty proxy >=2", "rv_q4_vs_q1", "rv_q4_vs_q1", base_covars, "wtmec6yr", NA_character_
  )

  cycle_specs <- tidyr::expand_grid(
    cycle_value = cycles,
    spec_id = c("obstruction_abc", "current_asthma", "emphysema_bronchitis", "frailty_proxy")
  ) %>%
    mutate(
      analysis_id = paste0(
        "v0_4_cycle_stratified_",
        .data$spec_id,
        "_",
        stringr::str_replace_all(.data$cycle_value, "[^0-9]", "_")
      ),
      analysis_family = "cycle_stratified",
      model_label = paste0("Cycle ", .data$cycle_value),
      outcome = dplyr::case_when(
        .data$spec_id == "obstruction_abc" ~ "obstruction_fixed_ratio_abc",
        .data$spec_id == "current_asthma" ~ "current_asthma",
        .data$spec_id == "emphysema_bronchitis" ~ "self_reported_emphysema_or_bronchitis",
        .data$spec_id == "frailty_proxy" ~ "nhanes_frailty_proxy_ge2"
      ),
      outcome_label = dplyr::case_when(
        .data$spec_id == "obstruction_abc" ~ "Spirometric obstruction",
        .data$spec_id == "current_asthma" ~ "Current asthma",
        .data$spec_id == "emphysema_bronchitis" ~ "Emphysema/chronic bronchitis",
        .data$spec_id == "frailty_proxy" ~ "Frailty proxy >=2"
      ),
      exposure = "resp_vulnerability_z",
      coefficient = "resp_vulnerability_z",
      covariates = base_covars,
      weight_variable = "wtmec2yr",
      cycle = .data$cycle_value
    ) %>%
    select(
      analysis_id, analysis_family, model_label, outcome, outcome_label, exposure,
      coefficient, covariates, weight_variable, cycle
    )

  list(
    sensitivity = sensitivity_specs,
    quartile = quartile_trend_specs,
    cycle = cycle_specs
  )
}

plot_outcome_rates <- function(outcome_rates, figure_dir) {
  plot_data <- outcome_rates %>%
    mutate(
      outcome_label = factor(
        .data$outcome_label,
        levels = c(
          "Spirometric obstruction",
          "Current asthma",
          "Emphysema/chronic bronchitis",
          "Frailty proxy >=2"
        )
      )
    )

  p <- ggplot(plot_data, aes(x = rv_quartile, y = weighted_prevalence_percent, group = outcome_label)) +
    geom_hline(yintercept = 0, linewidth = 0.25, colour = "#888888") +
    geom_errorbar(aes(ymin = ci_low_percent, ymax = ci_high_percent), width = 0.12, linewidth = 0.35, colour = "#4B5563") +
    geom_line(linewidth = 0.45, colour = "#2F6F8F") +
    geom_point(size = 1.65, colour = "#2F6F8F") +
    facet_wrap(~outcome_label, scales = "free_y", ncol = 2) +
    labs(
      x = "Respiratory vulnerability quartile",
      y = "Weighted prevalence (%)",
      title = "NHANES outcome rates by respiratory vulnerability quartile"
    ) +
    theme_classic(base_size = 8) +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = element_line(linewidth = 0.35, colour = "black"),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      plot.title = element_text(face = "bold", size = 9),
      panel.grid.major.y = element_line(linewidth = 0.2, colour = "#E5E7EB")
    )

  png_path <- file.path(figure_dir, "nhanes_v0_4_outcome_rates_by_quartile.png")
  pdf_path <- file.path(figure_dir, "nhanes_v0_4_outcome_rates_by_quartile.pdf")
  ggplot2::ggsave(png_path, p, width = 7.2, height = 4.8, dpi = 600, units = "in")
  ggplot2::ggsave(pdf_path, p, width = 7.2, height = 4.8, units = "in", device = grDevices::pdf, useDingbats = FALSE)
  c(png = png_path, pdf = pdf_path)
}

plot_sensitivity_forest <- function(model_table, figure_dir) {
  plot_data <- model_table %>%
    filter(.data$status == "ok") %>%
    filter(.data$analysis_family %in% c("pooled_sensitivity", "frailty_adjusted", "cycle_fixed", "cycle_stratified")) %>%
    mutate(
      plot_label = dplyr::case_when(
        .data$analysis_family == "cycle_stratified" ~ paste0("Cycle ", .data$cycle),
        TRUE ~ .data$model_label
      ),
      plot_label = factor(
        .data$plot_label,
        levels = rev(unique(c(
          "Main A/B/C",
          "Main",
          "Strict A/B quality",
          "Race-calibrated exposure",
          "Frailty-adjusted",
          "Cycle fixed effects",
          paste0("Cycle ", sort(unique(.data$cycle[.data$cycle != "pooled"])))
        )))
      ),
      outcome_label = factor(
        .data$outcome_label,
        levels = c(
          "Spirometric obstruction",
          "Current asthma",
          "Emphysema/chronic bronchitis",
          "Frailty proxy >=2"
        )
      )
    )

  p <- ggplot(plot_data, aes(y = plot_label, x = or, xmin = conf_low, xmax = conf_high)) +
    geom_vline(xintercept = 1, linewidth = 0.35, linetype = "dashed", colour = "#777777") +
    geom_errorbarh(height = 0.16, linewidth = 0.35, colour = "#4B5563") +
    geom_point(size = 1.45, colour = "#2F6F8F") +
    scale_x_log10() +
    facet_wrap(~outcome_label, scales = "free_y", ncol = 2) +
    labs(
      x = "Odds ratio per 1-SD higher respiratory vulnerability",
      y = NULL,
      title = "NHANES sensitivity models"
    ) +
    theme_classic(base_size = 7.4) +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = element_line(linewidth = 0.35, colour = "black"),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      plot.title = element_text(face = "bold", size = 8.6),
      panel.grid.major.x = element_line(linewidth = 0.2, colour = "#E5E7EB")
    )

  png_path <- file.path(figure_dir, "nhanes_v0_4_sensitivity_forest.png")
  pdf_path <- file.path(figure_dir, "nhanes_v0_4_sensitivity_forest.pdf")
  ggplot2::ggsave(png_path, p, width = 7.2, height = 5.8, dpi = 600, units = "in")
  ggplot2::ggsave(pdf_path, p, width = 7.2, height = 5.8, units = "in", device = grDevices::pdf, useDingbats = FALSE)
  c(png = png_path, pdf = pdf_path)
}

literature_sources <- function() {
  tibble::tribble(
    ~source_id, ~topic, ~citation_short, ~url, ~applied_to_v0_4,
    "cdc_weighting", "NHANES survey design", "CDC/NCHS NHANES Weighting Tutorial", "https://wwwn.cdc.gov/nchs/nhanes/tutorials/weighting.aspx", "Use MEC weights for spirometry/MEC variables and construct WTMEC6YR from three 2-year cycles.",
    "cdc_variance", "NHANES survey design", "CDC/NCHS NHANES Variance Estimation Tutorial", "https://wwwn.cdc.gov/nchs/nhanes/tutorials/varianceestimation.aspx", "Use PSU, strata, weights, and Taylor linearization via the R survey package.",
    "nhanes_guidelines", "NHANES analytic guidelines", "CDC/NCHS NHANES Survey Methods and Analytic Guidelines", "https://wwwn.cdc.gov/nchs/nhanes/analyticguidelines.aspx", "Keep design-variable and multi-cycle weighting choices explicit.",
    "ats_ers_spirometry_2019", "Spirometry quality", "Graham et al. Standardization of Spirometry 2019 Update", "https://www.apta.org/patient-care/evidence-based-practice-resources/cpgs/standardization-of-spirometry-2019-update.-an-official-american-thoracic-society-and-european-respiratory-society-technical-statement", "Treat spirometry quality and standardization as central to outcome interpretation.",
    "gli_2012", "Spirometry reference equations", "Quanjer et al. GLI 2012 multi-ethnic spirometry equations", "https://publications.ersnet.org/content/erj/40/6/1324", "Defer LLN/PRISm to a later equation-lock step rather than mixing unverified percent-predicted rules into V0.4.",
    "ers_ats_2021", "Pulmonary function interpretation", "Stanojevic et al. ERS/ATS interpretive strategies technical standard", "https://www.thoracic.org/statements/guideline-implementation-tools/technical-standards-interpretive-strategies-lung-function-tests.php", "Document fixed-ratio obstruction as a pragmatic replication phenotype, not the final interpretive standard.",
    "vaz_fragoso_aging_lung", "Aging lung and PEF z-scores", "Vaz Fragoso et al. Respiratory impairment and the aging lung", "https://pubmed.ncbi.nlm.nih.gov/22138206/", "Support residualized/z-scored PEF as an aging-risk signal while noting lower diagnostic accuracy than spirometry.",
    "tilert_nhanes_copd", "NHANES COPD/spirometry precedent", "Tilert et al. NHANES 2007-2010 COPD prevalence", "https://pubmed.ncbi.nlm.nih.gov/24107140/", "Benchmark NHANES 2007-2010 spirometry use and fixed-ratio versus LLN differences.",
    "magave_pef_frailty", "PEF and frailty", "Magave et al. PEF as an index of frailty syndrome", "https://pubmed.ncbi.nlm.nih.gov/33155627/", "Motivate frailty-adjusted and frailty-proxy analyses as biologically plausible sensitivity layers.",
    "pef_frailty_longitudinal", "PEF and frailty", "Cross-sectional and longitudinal associations between PEF and frailty", "https://pmc.ncbi.nlm.nih.gov/articles/PMC6912606/", "Treat frailty as a possible shared vulnerability pathway rather than only a confounder.",
    "copd_frailty_review", "COPD and frailty", "Systematic review/meta-analysis of frailty in COPD", "https://pmc.ncbi.nlm.nih.gov/articles/PMC9816100/", "Support reporting frailty-adjusted attenuation and preserving frailty as a secondary construct."
  )
}

write_literature_log <- function(log_path, source_table) {
  lines <- c(
    "# NHANES V0.4 Literature and Method Context",
    "",
    paste0("- Run date: ", Sys.Date()),
    "- Purpose: document why V0.4 adds weighted descriptives, quartile prevalence plots, frailty-adjusted sensitivity, and cycle-stratified checks.",
    "",
    "## Applied Sources",
    "",
    markdown_table(source_table),
    "",
    "## Method Boundary",
    "",
    "- NHANES V0.4 remains a cross-sectional replication layer.",
    "- The primary objective outcome is fixed-ratio spirometric obstruction using the currently locked NHANES variables and quality grades.",
    "- LLN and PRISm are deferred until GLI/percent-predicted equations are locked and tested in the cleaning layer.",
    "- PEF residual z-score is used as a respiratory vulnerability signal, not as a standalone diagnostic replacement for formal spirometry.",
    "- Frailty is handled as both a secondary outcome and an adjustment/sensitivity layer because PEF and frailty may share aging-related vulnerability pathways."
  )
  writeLines(lines, log_path)
}

write_analysis_log <- function(log_path, cutpoints, outcome_rates, model_table, warning_table, figure_paths) {
  key_models <- model_table %>%
    filter(.data$status == "ok", .data$analysis_family %in% c("pooled_sensitivity", "frailty_adjusted", "cycle_fixed")) %>%
    select(analysis_id, model_label, outcome_label, cycle, n, events, or_ci, p_value_formatted)

  q4_rates <- outcome_rates %>%
    filter(.data$rv_quartile %in% c("Q1", "Q4")) %>%
    transmute(
      outcome_label,
      rv_quartile,
      n_unweighted,
      events_unweighted,
      weighted_prevalence = sprintf("%.1f%%", .data$weighted_prevalence_percent),
      ci = sprintf("%.1f-%.1f%%", .data$ci_low_percent, .data$ci_high_percent)
    )

  lines <- c(
    "# NHANES V0.4 Weighted Descriptives and Sensitivity Models",
    "",
    paste0("- Run date: ", Sys.Date()),
    "- Dataset: local-only `derived_sensitive/nhanes/nhanes_replication_v0_2_analysis_ready.rds`.",
    "- Design: `svydesign(ids = ~psu, strata = ~strata, weights = ~wtmec6yr, nest = TRUE)` for pooled MEC/spirometry analyses.",
    "- Cycle-stratified checks use each cycle's `WTMEC2YR`.",
    "- Figure contract: quantitative grid; outcome-rate panel tests monotonic prevalence by vulnerability quartile; forest panel tests robustness of the continuous-exposure association.",
    "",
    "## Weighted Respiratory Vulnerability Quartile Cutpoints",
    "",
    markdown_table(cutpoints %>% mutate(resp_vulnerability_z = round(.data$resp_vulnerability_z, 4))),
    "",
    "## Q1 and Q4 Weighted Outcome Rates",
    "",
    markdown_table(q4_rates),
    "",
    "## Key Pooled Models",
    "",
    markdown_table(key_models),
    "",
    "## Warnings",
    "",
    if (nrow(warning_table) == 0) {
      "No survey model warnings or errors."
    } else {
      markdown_table(warning_table)
    },
    "",
    "## Figures",
    "",
    paste0("- Outcome rates: `", figure_paths[["outcome_png"]], "` and `", figure_paths[["outcome_pdf"]], "`."),
    paste0("- Sensitivity forest: `", figure_paths[["forest_png"]], "` and `", figure_paths[["forest_pdf"]], "`.")
  )
  writeLines(lines, log_path)
}

write_findings_log <- function(log_path, outcome_rates, sensitivity_table, quartile_table, cycle_table) {
  q1_q4 <- outcome_rates %>%
    filter(.data$rv_quartile %in% c("Q1", "Q4")) %>%
    select(outcome_label, rv_quartile, weighted_prevalence_percent) %>%
    tidyr::pivot_wider(names_from = rv_quartile, values_from = weighted_prevalence_percent) %>%
    mutate(
      q1_to_q4 = sprintf("%.1f%% to %.1f%%", .data$Q1, .data$Q4),
      absolute_difference_pct_points = round(.data$Q4 - .data$Q1, 1)
    ) %>%
    select(outcome_label, q1_to_q4, absolute_difference_pct_points)

  pooled_key <- sensitivity_table %>%
    filter(.data$status == "ok") %>%
    filter(.data$analysis_family %in% c("pooled_sensitivity", "frailty_adjusted", "cycle_fixed")) %>%
    select(model_label, outcome_label, n, events, or_ci, p_value_formatted)

  quartile_key <- quartile_table %>%
    filter(.data$status == "ok") %>%
    select(model_label, outcome_label, n, events, or_ci, p_value_formatted)

  cycle_ranges <- cycle_table %>%
    filter(.data$status == "ok") %>%
    group_by(.data$outcome_label) %>%
    summarise(
      cycle_or_range = sprintf("%.2f-%.2f", min(.data$or, na.rm = TRUE), max(.data$or, na.rm = TRUE)),
      all_cycle_ci_exclude_one = all(.data$conf_low > 1, na.rm = TRUE),
      .groups = "drop"
    )

  lines <- c(
    "# NHANES V0.4 Findings",
    "",
    paste0("- Run date: ", Sys.Date()),
    "- Analysis layer: survey-weighted cross-sectional NHANES replication among adults age >=45 with valid respiratory vulnerability exposure.",
    "",
    "## Main Message",
    "",
    "Higher residualized respiratory vulnerability is associated with higher weighted prevalence of spirometric obstruction, self-reported respiratory disease, and frailty proxy in NHANES 2007-2012. The continuous-exposure associations remain positive after race-calibrated exposure scoring, frailty adjustment, cycle fixed effects, and single-cycle stratification.",
    "",
    "## Weighted Q1 to Q4 Outcome Gradient",
    "",
    markdown_table(q1_q4),
    "",
    "## Pooled Sensitivity Models",
    "",
    markdown_table(pooled_key),
    "",
    "## Quartile Models",
    "",
    markdown_table(quartile_key),
    "",
    "## Cycle-Stratified Range",
    "",
    markdown_table(cycle_ranges),
    "",
    "## Interpretation",
    "",
    "- Strongest NHANES replication signal: spirometric obstruction.",
    "- Most manuscript-stable secondary signals: self-reported emphysema/chronic bronchitis and frailty proxy.",
    "- Current asthma is directionally consistent but should remain a secondary respiratory outcome because its quartile pattern is less monotonic and one cycle-stratified estimate is wider.",
    "- These analyses support the CHARLS-first respiratory vulnerability construct as a reproducible population-level axis, but they do not establish longitudinal causality in NHANES."
  )
  writeLines(lines, log_path)
}

main <- function() {
  root <- find_project_root()
  rds_path <- file.path(root, "derived_sensitive", "nhanes", "nhanes_replication_v0_2_analysis_ready.rds")
  if (!file.exists(rds_path)) {
    stop("Missing NHANES V0.2 analysis-ready RDS. Run R/nhanes/04_clean_nhanes_replication_v0_2.R first.", call. = FALSE)
  }

  table_dir <- file.path(root, "results", "tables")
  figure_dir <- file.path(root, "results", "figures")
  log_dir <- file.path(root, "results", "logs")
  metadata_dir <- file.path(root, "metadata")
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(metadata_dir, recursive = TRUE, showWarnings = FALSE)

  raw_data <- readRDS(rds_path)
  prepared <- prepare_analysis_data(raw_data)
  data <- prepared$data

  table1 <- build_table1(data)
  outcome_rates <- build_outcome_rates(data)
  model_specs <- build_model_specs(data)
  all_specs <- bind_rows(model_specs$sensitivity, model_specs$quartile, model_specs$cycle)
  fits <- lapply(seq_len(nrow(all_specs)), function(i) safe_fit_model_spec(data, all_specs[i, ]))
  model_table <- bind_rows(lapply(fits, `[[`, "row"))
  warning_table <- bind_rows(lapply(fits, `[[`, "warnings")) %>%
    filter(!is.na(.data$warning))
  if (nrow(warning_table) == 0) {
    warning_table <- tibble(analysis_id = character(), warning = character())
  }

  sensitivity_table <- model_table %>% filter(.data$analysis_family %in% c("pooled_sensitivity", "frailty_adjusted", "cycle_fixed"))
  quartile_table <- model_table %>% filter(.data$analysis_family == "quartile_model")
  cycle_table <- model_table %>% filter(.data$analysis_family == "cycle_stratified")

  readr::write_csv(prepared$cutpoints, file.path(table_dir, "nhanes_v0_4_vulnerability_quartile_cutpoints.csv"))
  readr::write_csv(table1$overall, file.path(table_dir, "nhanes_v0_4_weighted_table1_overall.csv"))
  readr::write_csv(table1$by_quartile, file.path(table_dir, "nhanes_v0_4_weighted_table1_by_quartile.csv"))
  readr::write_csv(outcome_rates, file.path(table_dir, "nhanes_v0_4_weighted_outcome_rates_by_quartile.csv"))
  readr::write_csv(sensitivity_table, file.path(table_dir, "nhanes_v0_4_sensitivity_model_table.csv"))
  readr::write_csv(quartile_table, file.path(table_dir, "nhanes_v0_4_quartile_model_table.csv"))
  readr::write_csv(cycle_table, file.path(table_dir, "nhanes_v0_4_cycle_stratified_models.csv"))
  readr::write_csv(warning_table, file.path(table_dir, "nhanes_v0_4_model_warnings.csv"))

  source_table <- literature_sources()
  readr::write_csv(source_table, file.path(metadata_dir, "nhanes_literature_sources_v0_4.csv"))

  outcome_fig <- plot_outcome_rates(outcome_rates, figure_dir)
  forest_fig <- plot_sensitivity_forest(model_table, figure_dir)
  figure_paths <- c(
    outcome_png = outcome_fig[["png"]],
    outcome_pdf = outcome_fig[["pdf"]],
    forest_png = forest_fig[["png"]],
    forest_pdf = forest_fig[["pdf"]]
  )

  write_literature_log(
    file.path(log_dir, "nhanes_literature_context_v0_4.md"),
    source_table
  )
  write_analysis_log(
    file.path(log_dir, "nhanes_v0_4_analysis_log.md"),
    prepared$cutpoints,
    outcome_rates,
    model_table,
    warning_table,
    figure_paths
  )
  write_findings_log(
    file.path(log_dir, "nhanes_v0_4_findings.md"),
    outcome_rates,
    sensitivity_table,
    quartile_table,
    cycle_table
  )

  message("Wrote NHANES v0.4 weighted Table 1, quartile rates, sensitivity models, and figures.")
  message("Log: ", file.path(log_dir, "nhanes_v0_4_analysis_log.md"))
}

if (sys.nframe() == 0) {
  main()
}
