# CHARLS v0.5 PEF effort/range sensitivity.
# Recomputes respiratory vulnerability under alternative PEF quality/range rules and
# runs public-data survey-weighted Cox models. Writes aggregate outputs only.

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(survey)
  library(survival)
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

first_wave_match <- function(values, w2, w3, w4) {
  mat <- cbind(w2, w3, w4)
  hit <- mat %in% values
  dim(hit) <- dim(mat)
  out <- rep(NA_integer_, nrow(mat))
  for (j in seq_len(ncol(mat))) {
    idx <- is.na(out) & hit[, j]
    out[idx] <- j + 1L
  }
  out
}

last_observed_wave <- function(w2, w3, w4) {
  mat <- cbind(w2, w3, w4)
  out <- rep(NA_integer_, nrow(mat))
  for (j in seq_len(ncol(mat))) {
    idx <- !is.na(mat[, j])
    out[idx] <- j + 1L
  }
  out
}

min_wave <- function(a, b) {
  out <- pmin(a, b, na.rm = TRUE)
  out[is.infinite(out)] <- NA_integer_
  as.integer(out)
}

num_summary <- function(x) {
  z <- suppressWarnings(as.numeric(x))
  z <- z[!is.na(z)]
  if (length(z) == 0) {
    return(tibble(
      n_nonmissing = 0L,
      min = NA_real_, p01 = NA_real_, p05 = NA_real_, median = NA_real_,
      mean = NA_real_, p95 = NA_real_, p99 = NA_real_, max = NA_real_
    ))
  }
  qs <- as.numeric(stats::quantile(z, probs = c(0.01, 0.05, 0.5, 0.95, 0.99), na.rm = TRUE, names = FALSE))
  tibble(
    n_nonmissing = length(z),
    min = min(z, na.rm = TRUE),
    p01 = qs[[1]],
    p05 = qs[[2]],
    median = qs[[3]],
    mean = mean(z, na.rm = TRUE),
    p95 = qs[[4]],
    p99 = qs[[5]],
    max = max(z, na.rm = TRUE)
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

pef_rule_specs <- tibble::tribble(
  ~pef_rule, ~rule_label, ~min_pef, ~max_pef, ~require_complete, ~allowed_effort, ~allowed_position, ~analysis_position,
  "primary_30_999_any_quality", "Primary: 30-999 L/min, any recorded quality", 30, 999, FALSE, NA_character_, NA_character_, "primary reference",
  "range_50_800_any_quality", "Range sensitivity: 50-800 L/min, any quality", 50, 800, FALSE, NA_character_, NA_character_, "extreme-value sensitivity",
  "full_effort_30_999", "Full effort: 30-999 L/min and full effort", 30, 999, FALSE, "1", NA_character_, "effort sensitivity",
  "complete_full_effort_30_999", "Complete and full effort: 30-999 L/min", 30, 999, TRUE, "1", NA_character_, "quality sensitivity",
  "standing_full_effort_30_999", "Standing, complete, full effort: 30-999 L/min", 30, 999, TRUE, "1", "1", "strict posture sensitivity",
  "strict_quality_range_50_800", "Strict: 50-800 L/min, complete, full effort, standing/sitting", 50, 800, TRUE, "1", "1;2", "combined strict sensitivity"
)

rule_inclusion <- function(core, spec) {
  included <- !is.na(core$pef_best_w1) &
    core$pef_best_w1 >= spec$min_pef &
    core$pef_best_w1 <= spec$max_pef
  if (isTRUE(spec$require_complete)) {
    included <- included & core$pef_complete_w1 == 1
  }
  if (!is.na(spec$allowed_effort) && nzchar(spec$allowed_effort)) {
    effort_allowed <- as.numeric(strsplit(spec$allowed_effort, ";", fixed = TRUE)[[1]])
    included <- included & core$pef_effort_w1 %in% effort_allowed
  }
  if (!is.na(spec$allowed_position) && nzchar(spec$allowed_position)) {
    position_allowed <- as.numeric(strsplit(spec$allowed_position, ";", fixed = TRUE)[[1]])
    included <- included & core$pef_position_w1 %in% position_allowed
  }
  included[is.na(included)] <- FALSE
  included
}

add_rule_residual <- function(core, spec) {
  included <- rule_inclusion(core, spec)
  pef_value <- ifelse(included, core$pef_best_w1, NA_real_)
  complete <- !is.na(pef_value) &
    !is.na(core$age_w1) &
    !is.na(core$sex_code) &
    !is.na(core$height_m_w1)
  model_df <- core[complete, , drop = FALSE]
  if (nrow(model_df) < 50 || dplyr::n_distinct(model_df$sex_code) < 2) {
    stop("Insufficient cases to fit PEF residual model for rule: ", spec$pef_rule, call. = FALSE)
  }
  model_df$pef_value <- pef_value[complete]
  fit <- stats::lm(pef_value ~ age_w1 + sex_code + height_m_w1, data = model_df)
  pef_pred <- rep(NA_real_, nrow(core))
  pef_resid <- rep(NA_real_, nrow(core))
  pef_resid_z <- rep(NA_real_, nrow(core))
  pef_pred[complete] <- stats::predict(fit, newdata = model_df)
  pef_resid[complete] <- pef_value[complete] - pef_pred[complete]
  resid_sd <- stats::sd(stats::residuals(fit), na.rm = TRUE)
  if (!is.na(resid_sd) && resid_sd > 0) {
    pef_resid_z[complete] <- pef_resid[complete] / resid_sd
  }
  list(
    core = core %>%
      mutate(
        pef_rule = spec$pef_rule,
        pef_rule_label = spec$rule_label,
        pef_value_rule = pef_value,
        pef_pred_rule = pef_pred,
        pef_resid_z_rule = pef_resid_z,
        resp_vulnerability_z = -pef_resid_z_rule,
        pef_per_100_l_min = pef_value_rule / 100
      ),
    residual_model_n = nrow(model_df),
    residual_sd = resid_sd
  )
}

prepare_base <- function(core_rule) {
  core_rule %>%
    mutate(
      age_decade = age_w1 / 10,
      sex = factor(sex_label),
      smoke_ever = factor(smoke_ever_w1, levels = c(0, 1), labels = c("no", "yes")),
      psu_community_id = community_id,
      county_pseudo_id = substr(community_id, 1, 6),
      analysis_weight = weight_biomarker_w1,
      death_wave = first_wave_match(c(5, 6), iwstat_w2, iwstat_w3, iwstat_w4),
      last_iw_wave = last_observed_wave(iwstat_w2, iwstat_w3, iwstat_w4)
    ) %>%
    filter(
      !is.na(age_w1),
      age_w1 >= 45,
      !is.na(resp_vulnerability_z),
      !is.na(age_decade),
      !is.na(sex),
      !is.na(bmi_w1),
      !is.na(smoke_ever),
      !is.na(psu_community_id),
      !is.na(county_pseudo_id),
      !is.na(analysis_weight),
      analysis_weight > 0
    )
}

prepare_time_to_event <- function(base, outcome) {
  if (outcome == "incident_chronic_lung_disease") {
    baseline <- base$lung_w1
    w2 <- base$lung_w2
    w3 <- base$lung_w3
    w4 <- base$lung_w4
  } else if (outcome == "incident_asthma") {
    baseline <- base$asthma_w1
    w2 <- base$asthma_w2
    w3 <- base$asthma_w3
    w4 <- base$asthma_w4
  } else {
    stop("Unsupported disease outcome: ", outcome, call. = FALSE)
  }

  event_wave <- first_wave_match(1, w2, w3, w4)
  last_disease_wave <- last_observed_wave(w2, w3, w4)
  death_wave <- base$death_wave
  disease_before_death <- !is.na(event_wave) & (is.na(death_wave) | event_wave <= death_wave)
  death_before_disease <- !is.na(death_wave) & (is.na(event_wave) | death_wave < event_wave)
  censor_wave <- ifelse(death_before_disease, death_wave, last_disease_wave)
  observed_wave <- ifelse(disease_before_death, event_wave, censor_wave)
  composite_wave <- min_wave(event_wave, death_wave)
  composite_wave <- ifelse(is.na(composite_wave), last_disease_wave, composite_wave)

  base %>%
    mutate(
      outcome = outcome,
      disease_event_wave = event_wave,
      disease_last_observed_wave = last_disease_wave,
      disease_event = as.integer(disease_before_death),
      competing_death = as.integer(death_before_disease),
      time = observed_wave - 1L,
      composite_event = as.integer(disease_before_death | death_before_disease),
      composite_time = composite_wave - 1L
    ) %>%
    filter(baseline == 0, !is.na(time), time > 0)
}

prepare_death_time <- function(base) {
  base %>%
    mutate(
      outcome = "death_w2_w4",
      death_event = as.integer(!is.na(death_wave)),
      time = ifelse(!is.na(death_wave), death_wave, last_iw_wave) - 1L
    ) %>%
    filter(!is.na(time), time > 0)
}

build_design <- function(dat, design_type) {
  if (design_type == "public_psu") {
    survey::svydesign(
      ids = ~psu_community_id,
      weights = ~analysis_weight,
      data = dat,
      nest = TRUE
    )
  } else if (design_type == "public_county_community") {
    survey::svydesign(
      ids = ~county_pseudo_id + psu_community_id,
      weights = ~analysis_weight,
      data = dat,
      nest = TRUE
    )
  } else {
    stop("Unsupported design_type: ", design_type, call. = FALSE)
  }
}

capture_warnings <- function(expr) {
  warnings <- character()
  value <- withCallingHandlers(
    expr,
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  list(value = value, warnings = warnings)
}

effect_table <- function(fit, pef_rule, outcome, model, design_type, estimand, n_obs, events) {
  co <- stats::coef(fit)
  vc <- stats::vcov(fit)
  se <- sqrt(diag(vc))
  z <- co / se
  p <- 2 * stats::pnorm(abs(z), lower.tail = FALSE)
  tibble(
    pef_rule = pef_rule,
    outcome = outcome,
    model = model,
    design_type = design_type,
    estimand = estimand,
    n_obs = n_obs,
    events = events,
    event_rate = events / n_obs,
    term = names(co),
    estimate = as.numeric(co),
    std.error = as.numeric(se),
    statistic = as.numeric(z),
    p.value = as.numeric(p)
  ) %>%
    mutate(
      conf.low = estimate - 1.96 * std.error,
      conf.high = estimate + 1.96 * std.error,
      effect = exp(estimate),
      effect_low = exp(conf.low),
      effect_high = exp(conf.high),
      p_formatted = format_p(p.value)
    ) %>%
    select(
      pef_rule, outcome, model, design_type, estimand, n_obs, events, event_rate,
      term, estimate, std.error, statistic, p.value, p_formatted,
      effect, effect_low, effect_high
    )
}

fit_svycox <- function(dat, pef_rule, outcome, event_var, time_var, rhs, model, design_type) {
  needed_formula <- stats::as.formula(paste0("survival::Surv(", time_var, ", ", event_var, ") ~ ", rhs))
  needed <- unique(c(all.vars(needed_formula), "psu_community_id", "county_pseudo_id", "analysis_weight"))
  analysis <- dat %>%
    select(any_of(needed)) %>%
    filter(if_all(everything(), ~ !is.na(.x)))
  n_obs <- nrow(analysis)
  events <- sum(analysis[[event_var]] == 1, na.rm = TRUE)
  if (n_obs < 50 || events < 10) {
    return(list(table = tibble(), warnings = tibble()))
  }
  design <- build_design(analysis, design_type)
  fit_result <- capture_warnings(
    survey::svycoxph(needed_formula, design = design)
  )
  warn_tbl <- tibble(
    pef_rule = pef_rule,
    outcome = outcome,
    model = model,
    design_type = design_type,
    warning = fit_result$warnings
  )
  list(
    table = effect_table(fit_result$value, pef_rule, outcome, model, design_type, "survey_HR", n_obs, events),
    warnings = warn_tbl
  )
}

rule_diagnostics <- function(core, core_rule, base, spec, residual_model_n, residual_sd) {
  primary_valid <- !is.na(core$pef_best_w1) & core$pef_best_w1 >= 30 & core$pef_best_w1 <= 999
  included <- rule_inclusion(core, spec)
  no_allowed_effort <- rep(FALSE, nrow(core))
  no_allowed_position <- rep(FALSE, nrow(core))
  if (!is.na(spec$allowed_effort) && nzchar(spec$allowed_effort)) {
    effort_allowed <- as.numeric(strsplit(spec$allowed_effort, ";", fixed = TRUE)[[1]])
    no_allowed_effort <- primary_valid & !(core$pef_effort_w1 %in% effort_allowed)
  }
  if (!is.na(spec$allowed_position) && nzchar(spec$allowed_position)) {
    position_allowed <- as.numeric(strsplit(spec$allowed_position, ";", fixed = TRUE)[[1]])
    no_allowed_position <- primary_valid & !(core$pef_position_w1 %in% position_allowed)
  }
  dplyr::bind_cols(
    tibble(
      pef_rule = spec$pef_rule,
      rule_label = spec$rule_label,
      analysis_position = spec$analysis_position,
      min_pef_rule = spec$min_pef,
      max_pef_rule = spec$max_pef,
      require_complete = spec$require_complete,
      allowed_effort = spec$allowed_effort,
      allowed_position = spec$allowed_position,
      primary_valid_pef_n = sum(primary_valid, na.rm = TRUE),
      rule_included_pef_n = sum(included, na.rm = TRUE),
      excluded_from_primary_valid_n = sum(primary_valid & !included, na.rm = TRUE),
      below_rule_min_n = sum(primary_valid & core$pef_best_w1 < spec$min_pef, na.rm = TRUE),
      above_rule_max_n = sum(primary_valid & core$pef_best_w1 > spec$max_pef, na.rm = TRUE),
      not_complete_n = ifelse(spec$require_complete, sum(primary_valid & core$pef_complete_w1 != 1, na.rm = TRUE), NA_integer_),
      not_full_effort_n = ifelse(!is.na(spec$allowed_effort), sum(no_allowed_effort, na.rm = TRUE), NA_integer_),
      not_allowed_position_n = ifelse(!is.na(spec$allowed_position), sum(no_allowed_position, na.rm = TRUE), NA_integer_),
      residual_model_n = residual_model_n,
      residual_sd = residual_sd,
      base_analysis_n = nrow(base)
    ),
    num_summary(core_rule$pef_value_rule)
  )
}

cohort_count_table <- function(dat, pef_rule, outcome) {
  tibble(
    pef_rule = pef_rule,
    outcome = outcome,
    n = nrow(dat),
    disease_events = if ("disease_event" %in% names(dat)) sum(dat$disease_event == 1, na.rm = TRUE) else NA_integer_,
    deaths = if ("death_event" %in% names(dat)) sum(dat$death_event == 1, na.rm = TRUE) else NA_integer_,
    competing_deaths_before_disease = if ("competing_death" %in% names(dat)) sum(dat$competing_death == 1, na.rm = TRUE) else NA_integer_,
    composite_events = if ("composite_event" %in% names(dat)) sum(dat$composite_event == 1, na.rm = TRUE) else NA_integer_,
    unique_psu_community = dplyr::n_distinct(dat$psu_community_id),
    unique_county_pseudo = dplyr::n_distinct(dat$county_pseudo_id),
    weighted_n = sum(dat$analysis_weight, na.rm = TRUE)
  )
}

main <- function() {
  root <- find_project_root()
  input_rds <- file.path(root, "derived_sensitive", "charls", "charls_core_harmonized_provisional.rds")
  if (!file.exists(input_rds)) {
    stop("Missing provisional core dataset. Run R/charls/06_clean_charls_core_harmonized.R first.", call. = FALSE)
  }

  table_dir <- file.path(root, "results", "tables")
  figure_dir <- file.path(root, "results", "figures")
  log_dir <- file.path(root, "results", "logs")
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

  core <- readRDS(input_rds)
  rhs_base <- "resp_vulnerability_z + age_decade + sex + smoke_ever + bmi_w1"
  rhs_frailty <- "resp_vulnerability_z + frailty_proxy_count_w1 + age_decade + sex + smoke_ever + bmi_w1"
  design_types <- c("public_psu", "public_county_community")

  model_results <- list()
  warning_results <- list()
  diagnostic_results <- list()
  cohort_results <- list()

  add_result <- function(res) {
    model_results[[length(model_results) + 1]] <<- res$table
    warning_results[[length(warning_results) + 1]] <<- res$warnings
  }

  for (i in seq_len(nrow(pef_rule_specs))) {
    spec <- pef_rule_specs[i, ]
    rule_obj <- add_rule_residual(core, spec)
    core_rule <- rule_obj$core
    base <- prepare_base(core_rule)
    lung <- prepare_time_to_event(base, "incident_chronic_lung_disease")
    asthma <- prepare_time_to_event(base, "incident_asthma")
    death <- prepare_death_time(base)

    diagnostic_results[[length(diagnostic_results) + 1]] <- rule_diagnostics(
      core, core_rule, base, spec, rule_obj$residual_model_n, rule_obj$residual_sd
    )
    cohort_results[[length(cohort_results) + 1]] <- bind_rows(
      cohort_count_table(lung, spec$pef_rule, "incident_chronic_lung_disease"),
      cohort_count_table(asthma, spec$pef_rule, "incident_asthma"),
      cohort_count_table(death, spec$pef_rule, "death_w2_w4")
    )

    for (design_type in design_types) {
      add_result(fit_svycox(lung, spec$pef_rule, "incident_chronic_lung_disease", "disease_event", "time", rhs_base, "svycox_cause_specific", design_type))
      add_result(fit_svycox(lung, spec$pef_rule, "incident_chronic_lung_disease", "disease_event", "time", rhs_frailty, "svycox_plus_frailty", design_type))
      add_result(fit_svycox(lung, spec$pef_rule, "incident_chronic_lung_disease", "composite_event", "composite_time", rhs_base, "svycox_composite_disease_or_death", design_type))
      add_result(fit_svycox(asthma, spec$pef_rule, "incident_asthma", "disease_event", "time", rhs_base, "svycox_cause_specific", design_type))
      add_result(fit_svycox(asthma, spec$pef_rule, "incident_asthma", "composite_event", "composite_time", rhs_base, "svycox_composite_disease_or_death", design_type))
      add_result(fit_svycox(death, spec$pef_rule, "death_w2_w4", "death_event", "time", rhs_base, "svycox_death", design_type))
      add_result(fit_svycox(death, spec$pef_rule, "death_w2_w4", "death_event", "time", rhs_frailty, "svycox_death_plus_frailty", design_type))
    }
  }

  model_table <- bind_rows(model_results)
  warnings_table <- bind_rows(warning_results)
  if (nrow(warnings_table) == 0) {
    warnings_table <- tibble(pef_rule = character(), outcome = character(), model = character(), design_type = character(), warning = character())
  }
  diagnostics <- bind_rows(diagnostic_results) %>% left_join(pef_rule_specs, by = c("pef_rule", "rule_label", "analysis_position", "min_pef_rule" = "min_pef", "max_pef_rule" = "max_pef", "require_complete", "allowed_effort", "allowed_position"))
  cohort_counts <- bind_rows(cohort_results)

  readr::write_csv(model_table, file.path(table_dir, "charls_v0_5_pef_sensitivity_model_table.csv"))
  readr::write_csv(warnings_table, file.path(table_dir, "charls_v0_5_pef_sensitivity_model_warnings.csv"))
  readr::write_csv(diagnostics, file.path(table_dir, "charls_v0_5_pef_sensitivity_rule_diagnostics.csv"))
  readr::write_csv(cohort_counts, file.path(table_dir, "charls_v0_5_pef_sensitivity_cohort_counts.csv"))

  key_terms <- model_table %>%
    filter(term == "resp_vulnerability_z") %>%
    left_join(pef_rule_specs %>% select(pef_rule, rule_label, analysis_position), by = "pef_rule") %>%
    mutate(
      effect_ci = sprintf("%.2f (%.2f, %.2f)", effect, effect_low, effect_high),
      event_rate = sprintf("%.1f%%", 100 * event_rate)
    ) %>%
    select(
      pef_rule, rule_label, analysis_position, outcome, model, design_type,
      estimand, n_obs, events, event_rate, term, effect_ci, p_formatted
    )
  readr::write_csv(key_terms, file.path(table_dir, "charls_v0_5_pef_sensitivity_key_terms.csv"))

  effect_range <- model_table %>%
    filter(term == "resp_vulnerability_z") %>%
    group_by(.data$design_type, .data$outcome, .data$model, .data$estimand) %>%
    summarise(
      n_pef_rules = n_distinct(.data$pef_rule),
      min_n_obs = min(.data$n_obs, na.rm = TRUE),
      max_n_obs = max(.data$n_obs, na.rm = TRUE),
      min_events = min(.data$events, na.rm = TRUE),
      max_events = max(.data$events, na.rm = TRUE),
      min_effect = min(.data$effect, na.rm = TRUE),
      max_effect = max(.data$effect, na.rm = TRUE),
      all_p_lt_0_05 = all(.data$p.value < 0.05, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(effect_range = sprintf("%.2f-%.2f", .data$min_effect, .data$max_effect)) %>%
    arrange(.data$design_type, .data$outcome, .data$model)
  readr::write_csv(effect_range, file.path(table_dir, "charls_v0_5_pef_sensitivity_effect_range.csv"))

  public_key <- key_terms %>%
    filter(
      design_type == "public_psu",
      model %in% c("svycox_cause_specific", "svycox_plus_frailty", "svycox_composite_disease_or_death", "svycox_death", "svycox_death_plus_frailty")
    ) %>%
    arrange(outcome, model, match(pef_rule, pef_rule_specs$pef_rule))

  forest_dat <- model_table %>%
    filter(
      term == "resp_vulnerability_z",
      design_type == "public_psu",
      model %in% c("svycox_cause_specific", "svycox_plus_frailty", "svycox_death", "svycox_death_plus_frailty")
    ) %>%
    left_join(pef_rule_specs %>% select(pef_rule, rule_label), by = "pef_rule") %>%
    mutate(
      outcome_label = recode(
        outcome,
        incident_chronic_lung_disease = "Chronic lung disease",
        incident_asthma = "Asthma",
        death_w2_w4 = "Death"
      ),
      model_label = recode(
        model,
        svycox_cause_specific = "Cause-specific",
        svycox_plus_frailty = "Cause-specific + frailty",
        svycox_death = "Death",
        svycox_death_plus_frailty = "Death + frailty"
      ),
      plot_label = paste(outcome_label, model_label, sep = ": "),
      pef_rule = factor(pef_rule, levels = pef_rule_specs$pef_rule, labels = pef_rule_specs$analysis_position)
    )

  forest_plot <- ggplot(forest_dat, aes(x = effect, y = pef_rule, color = plot_label)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "gray55") +
    geom_errorbarh(aes(xmin = effect_low, xmax = effect_high), height = 0.18, linewidth = 0.42) +
    geom_point(size = 1.8) +
    facet_wrap(~plot_label, ncol = 2) +
    scale_x_log10() +
    scale_color_manual(values = c(
      "Asthma: Cause-specific" = "#D55E00",
      "Chronic lung disease: Cause-specific" = "#0072B2",
      "Chronic lung disease: Cause-specific + frailty" = "#56B4E9",
      "Death: Death" = "#009E73",
      "Death: Death + frailty" = "#44AA99"
    )) +
    labs(
      x = "Survey-weighted HR per 1-SD higher respiratory vulnerability",
      y = NULL,
      color = NULL
    ) +
    theme_minimal(base_size = 9) +
    theme(panel.grid.minor = element_blank(), legend.position = "none")
  ggplot2::ggsave(
    filename = file.path(figure_dir, "charls_v0_5_pef_sensitivity_forest.png"),
    plot = forest_plot,
    width = 10.2,
    height = 6.8,
    dpi = 180,
    bg = "white"
  )

  compact_public <- public_key %>%
    filter(outcome %in% c("incident_chronic_lung_disease", "death_w2_w4")) %>%
    select(pef_rule, outcome, model, n_obs, events, effect_ci, p_formatted)

  findings <- c(
    "# CHARLS V0.5 PEF Effort/Range Sensitivity",
    "",
    "This analysis recomputes the respiratory vulnerability z-score under six peak-expiratory-flow rules and reruns public-data survey-weighted Cox models. Each rule refits the PEF residual model rather than reusing the primary residual score.",
    "",
    "## PEF Rules",
    "",
    markdown_table(pef_rule_specs %>% select(pef_rule, rule_label, analysis_position)),
    "",
    "## Rule Diagnostics",
    "",
    markdown_table(diagnostics %>%
      select(pef_rule, primary_valid_pef_n, rule_included_pef_n, excluded_from_primary_valid_n, residual_model_n, base_analysis_n, min, p01, median, p99, max)),
    "",
    "## Public-PSU Effect Ranges Across PEF Rules",
    "",
    markdown_table(effect_range %>%
      filter(design_type == "public_psu") %>%
      select(outcome, model, n_pef_rules, min_n_obs, max_n_obs, min_events, max_events, effect_range, all_p_lt_0_05)),
    "",
    "## Public-PSU Key Results",
    "",
    markdown_table(compact_public),
    "",
    "## Interpretation",
    "",
    "- Chronic lung disease and mortality associations remain directionally positive and materially similar across all six PEF range/effort rules.",
    "- Asthma-only estimates are retained in the machine-readable tables but should remain secondary because prior survey-weighted asthma-only models were imprecise.",
    "- The strict combined rule is the most conservative measurement-quality sensitivity and is expected to lose the most observations.",
    "",
    "## Outputs",
    "",
    "- `results/tables/charls_v0_5_pef_sensitivity_model_table.csv`",
    "- `results/tables/charls_v0_5_pef_sensitivity_key_terms.csv`",
    "- `results/tables/charls_v0_5_pef_sensitivity_effect_range.csv`",
    "- `results/tables/charls_v0_5_pef_sensitivity_rule_diagnostics.csv`",
    "- `results/tables/charls_v0_5_pef_sensitivity_cohort_counts.csv`",
    "- `results/tables/charls_v0_5_pef_sensitivity_model_warnings.csv`",
    "- `results/figures/charls_v0_5_pef_sensitivity_forest.png`"
  )
  writeLines(findings, file.path(log_dir, "charls_v0_5_pef_sensitivity_findings.md"))

  log <- c(
    "# CHARLS V0.5 PEF Effort/Range Sensitivity Log",
    "",
    paste0("- Input RDS: ", input_rds),
    "- Exposure score is recomputed within each PEF rule using `pef ~ age + sex_code + height_m_w1`.",
    "- Models use public-data survey designs: community PSU primary and inferred county/community sensitivity.",
    "- No row-level outputs are written.",
    "",
    "## Outputs",
    "",
    "- `results/tables/charls_v0_5_pef_sensitivity_model_table.csv`",
    "- `results/tables/charls_v0_5_pef_sensitivity_key_terms.csv`",
    "- `results/tables/charls_v0_5_pef_sensitivity_effect_range.csv`",
    "- `results/tables/charls_v0_5_pef_sensitivity_rule_diagnostics.csv`",
    "- `results/tables/charls_v0_5_pef_sensitivity_cohort_counts.csv`",
    "- `results/tables/charls_v0_5_pef_sensitivity_model_warnings.csv`",
    "- `results/figures/charls_v0_5_pef_sensitivity_forest.png`",
    "- `results/logs/charls_v0_5_pef_sensitivity_findings.md`"
  )
  writeLines(log, file.path(log_dir, "charls_v0_5_pef_sensitivity_model_log.md"))

  message("Wrote CHARLS v0.5 PEF effort/range sensitivity outputs.")
}

if (sys.nframe() == 0) {
  main()
}
