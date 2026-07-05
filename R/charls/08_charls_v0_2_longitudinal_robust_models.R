# CHARLS v0.2 longitudinal robust respiratory vulnerability models.
# Reads local-only provisional RDS and writes aggregate outputs only.

suppressPackageStartupMessages({
  library(broom)
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(survival)
  library(tidyr)
})

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

normalize_weight <- function(w) {
  mean_w <- mean(w[w > 0], na.rm = TRUE)
  if (!is.finite(mean_w) || mean_w <= 0) {
    return(rep(NA_real_, length(w)))
  }
  ifelse(!is.na(w) & w > 0, w / mean_w, NA_real_)
}

prepare_base <- function(core) {
  core %>%
    mutate(
      resp_vulnerability_z = -pef_resid_z_w1,
      pef_per_100_l_min = pef_best_w1_valid_provisional / 100,
      age_decade = age_w1 / 10,
      sex = factor(sex_label),
      smoke_ever = factor(smoke_ever_w1, levels = c(0, 1), labels = c("no", "yes")),
      frailty_binary = factor(frailty_proxy_ge3_w1, levels = c(0, 1), labels = c("non_frail_proxy", "frail_proxy")),
      cluster_id = ifelse(!is.na(household_id) & household_id != "", household_id, participant_id),
      analysis_weight = normalize_weight(weight_biomarker_w1),
      death_wave = first_wave_match(c(5, 6), iwstat_w2, iwstat_w3, iwstat_w4),
      death_time = death_wave - 1L,
      last_iw_wave = last_observed_wave(iwstat_w2, iwstat_w3, iwstat_w4)
    ) %>%
    filter(
      baseline_respiratory_axis_eligible,
      !is.na(resp_vulnerability_z),
      !is.na(age_decade),
      !is.na(sex),
      !is.na(bmi_w1),
      !is.na(smoke_ever)
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

cox_table <- function(fit, outcome, model, n_obs, events, estimand = "HR") {
  s <- summary(fit)
  coef_df <- as.data.frame(s$coefficients)
  if (nrow(coef_df) == 0) {
    return(tibble())
  }
  coef_col <- "coef"
  se_col <- if ("robust se" %in% names(coef_df)) "robust se" else "se(coef)"
  z_col <- if ("robust z" %in% names(coef_df)) "robust z" else "z"
  p_col <- grep("^Pr", names(coef_df), value = TRUE)[1]
  tibble(
    outcome = outcome,
    model = model,
    estimand = estimand,
    n_obs = n_obs,
    events = events,
    event_rate = events / n_obs,
    term = rownames(coef_df),
    estimate = coef_df[[coef_col]],
    std.error = coef_df[[se_col]],
    statistic = coef_df[[z_col]],
    p.value = coef_df[[p_col]]
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
      outcome, model, estimand, n_obs, events, event_rate,
      term, estimate, std.error, statistic, p.value, p_formatted,
      effect, effect_low, effect_high
    )
}

fit_cox <- function(dat, outcome, event_var, time_var, rhs, model, weighted = FALSE) {
  formula <- stats::as.formula(paste0("survival::Surv(", time_var, ", ", event_var, ") ~ ", rhs, " + cluster(cluster_id)"))
  needed <- unique(c(all.vars(formula), "cluster_id", if (weighted) "analysis_weight" else character()))
  analysis <- dat %>%
    select(any_of(needed)) %>%
    filter(if_all(everything(), ~ !is.na(.x)))
  if (weighted) {
    analysis <- analysis %>% filter(analysis_weight > 0)
  }
  n_obs <- nrow(analysis)
  events <- sum(analysis[[event_var]] == 1, na.rm = TRUE)
  if (n_obs < 50 || events < 10) {
    return(tibble())
  }
  if (weighted) {
    fit <- survival::coxph(formula, data = analysis, weights = analysis_weight, robust = TRUE)
  } else {
    fit <- survival::coxph(formula, data = analysis, robust = TRUE)
  }
  cox_table(fit, outcome, model, n_obs, events)
}

glm_or_table <- function(fit, outcome, model, n_obs, events, estimand = "OR") {
  broom::tidy(fit) %>%
    mutate(
      outcome = outcome,
      model = model,
      estimand = estimand,
      n_obs = n_obs,
      events = events,
      event_rate = events / n_obs,
      conf.low = estimate - 1.96 * std.error,
      conf.high = estimate + 1.96 * std.error,
      effect = exp(estimate),
      effect_low = exp(conf.low),
      effect_high = exp(conf.high),
      p_formatted = format_p(p.value)
    ) %>%
    select(
      outcome, model, estimand, n_obs, events, event_rate,
      term, estimate, std.error, statistic, p.value, p_formatted,
      effect, effect_low, effect_high
    )
}

fit_discrete_time <- function(period_dat, outcome, rhs, model, weighted = FALSE) {
  formula <- stats::as.formula(paste0("period_event ~ ", rhs, " + factor(interval)"))
  needed <- unique(c(all.vars(formula), if (weighted) "analysis_weight" else character()))
  analysis <- period_dat %>%
    select(any_of(needed)) %>%
    filter(if_all(everything(), ~ !is.na(.x)))
  if (weighted) {
    analysis <- analysis %>% filter(analysis_weight > 0)
  }
  n_obs <- nrow(analysis)
  events <- sum(analysis$period_event == 1, na.rm = TRUE)
  if (n_obs < 50 || events < 10) {
    return(tibble())
  }
  if (weighted) {
    fit <- stats::glm(formula, data = analysis, weights = analysis_weight, family = stats::quasibinomial())
  } else {
    fit <- stats::glm(formula, data = analysis, family = stats::binomial())
  }
  glm_or_table(fit, outcome, model, n_obs, events)
}

make_period_data <- function(dat, outcome) {
  periods <- lapply(1:3, function(interval) {
    target_wave <- interval + 1L
    dat %>%
      mutate(
        interval = interval,
        target_wave = target_wave,
        period_event = as.integer(!is.na(disease_event_wave) & disease_event_wave == target_wave),
        at_risk = is.na(disease_event_wave) | disease_event_wave >= target_wave,
        before_death = is.na(death_wave) | death_wave > target_wave | (!is.na(disease_event_wave) & disease_event_wave == target_wave)
      ) %>%
      filter(at_risk, before_death)
  })
  bind_rows(periods) %>% mutate(outcome = outcome)
}

wave_count_table <- function(dat, outcome) {
  tibble(
    outcome = outcome,
    disease_event_wave_2 = sum(dat$disease_event_wave == 2, na.rm = TRUE),
    disease_event_wave_3 = sum(dat$disease_event_wave == 3, na.rm = TRUE),
    disease_event_wave_4 = sum(dat$disease_event_wave == 4, na.rm = TRUE),
    death_before_disease_wave_2 = sum(dat$competing_death == 1 & dat$death_wave == 2, na.rm = TRUE),
    death_before_disease_wave_3 = sum(dat$competing_death == 1 & dat$death_wave == 3, na.rm = TRUE),
    death_before_disease_wave_4 = sum(dat$competing_death == 1 & dat$death_wave == 4, na.rm = TRUE)
  )
}

cohort_count_table <- function(dat, outcome) {
  tibble(
    outcome = outcome,
    n = nrow(dat),
    disease_events = if ("disease_event" %in% names(dat)) sum(dat$disease_event == 1, na.rm = TRUE) else NA_integer_,
    deaths = if ("death_event" %in% names(dat)) sum(dat$death_event == 1, na.rm = TRUE) else NA_integer_,
    competing_deaths_before_disease = if ("competing_death" %in% names(dat)) sum(dat$competing_death == 1, na.rm = TRUE) else NA_integer_,
    composite_events = if ("composite_event" %in% names(dat)) sum(dat$composite_event == 1, na.rm = TRUE) else NA_integer_,
    median_time = stats::median(dat$time, na.rm = TRUE),
    weighted_n = sum(dat$analysis_weight, na.rm = TRUE)
  )
}

make_phenotype_lock_candidates <- function() {
  tibble::tribble(
    ~concept, ~selected_source, ~selected_variables, ~derived_variable, ~role, ~v0_2_rule, ~lock_status,
    "baseline respiratory reserve", "H_CHARLS_D_Data.dta", "r1puff; r1puff1-r1puff3; r1puffcomp; r1puffpos; r1puffeff", "pef_best_w1_valid_provisional", "primary exposure input", "Use Harmonized wave-1 maximum PEF; retain 30-999; flag >=900; confirm with official codebook.", "candidate_needs_codebook_lock",
    "respiratory vulnerability score", "derived", "r1puff; r1agey; ragender; r1mheight", "resp_vulnerability_z = -pef_resid_z_w1", "primary exposure", "PEF residual z adjusted for age, sex code, and measured height; higher means lower-than-expected reserve.", "candidate_needs_model_lock",
    "incident chronic lung disease", "H_CHARLS_D_Data.dta", "r1lunge-r4lunge", "disease_event; disease_event_wave", "primary outcome", "Baseline r1lunge == 0; first follow-up lung status == 1 in waves 2-4; death before event treated as competing censoring in cause-specific models.", "candidate_needs_codebook_lock",
    "incident asthma", "H_CHARLS_D_Data.dta", "r1asthmae-r4asthmae", "disease_event; disease_event_wave", "secondary outcome", "Baseline r1asthmae == 0; first follow-up asthma status == 1 in waves 2-4; death before event treated as competing censoring in cause-specific models.", "candidate_needs_codebook_lock",
    "mortality", "H_CHARLS_D_Data.dta", "r2iwstat-r4iwstat", "death_wave; death_event", "negative-control / severity outcome", "First wave with iwstat 5 or 6 defines death; time approximated by wave interval.", "candidate_needs_codebook_lock",
    "baseline frailty proxy", "H_CHARLS_D_Data.dta", "r1gripsum; r1walk1kma; r1walk100a; r1adlab_c; r1cesd10; r1vgact_c; r1mdact_c; r1ltact_c", "frailty_proxy_count_w1; frailty_proxy_ge3_w1", "effect modifier / covariate", "Five-component proxy: low grip, mobility difficulty, ADL difficulty, CESD-10 >=10, low activity.", "candidate_needs_clinical_lock",
    "core covariates", "H_CHARLS_D_Data.dta", "r1agey; ragender; r1mheight; r1mbmi; r1smokev", "age_decade; sex; bmi_w1; smoke_ever", "adjustment set", "Adjust for baseline age, sex, smoking ever, and BMI; height included in raw PEF sensitivity.", "candidate_needs_codebook_lock",
    "weights", "H_CHARLS_D_Data.dta", "r1wtrespbiob; r1wtrespbioa", "analysis_weight", "sensitivity analysis", "Use normalized positive wave-1 biomarker weight; formal survey design awaits PSU/strata lock.", "candidate_needs_design_lock"
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
  metadata_dir <- file.path(root, "metadata")
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(metadata_dir, recursive = TRUE, showWarnings = FALSE)

  core <- readRDS(input_rds)
  base <- prepare_base(core)
  lung <- prepare_time_to_event(base, "incident_chronic_lung_disease")
  asthma <- prepare_time_to_event(base, "incident_asthma")
  death <- prepare_death_time(base)

  rhs_base <- "resp_vulnerability_z + age_decade + sex + smoke_ever + bmi_w1"
  rhs_frailty <- "resp_vulnerability_z + frailty_proxy_count_w1 + age_decade + sex + smoke_ever + bmi_w1"
  rhs_raw_pef <- "pef_per_100_l_min + age_decade + sex + smoke_ever + bmi_w1 + height_m_w1"

  cox_results <- bind_rows(
    fit_cox(lung, "incident_chronic_lung_disease", "disease_event", "time", rhs_base, "cox_cause_specific_death_censored"),
    fit_cox(lung, "incident_chronic_lung_disease", "disease_event", "time", rhs_frailty, "cox_cause_specific_plus_frailty"),
    fit_cox(lung, "incident_chronic_lung_disease", "composite_event", "composite_time", rhs_base, "cox_composite_disease_or_death"),
    fit_cox(lung, "incident_chronic_lung_disease", "disease_event", "time", rhs_raw_pef, "cox_raw_pef_sensitivity"),
    fit_cox(lung, "incident_chronic_lung_disease", "disease_event", "time", rhs_base, "cox_weighted_cause_specific", weighted = TRUE),
    fit_cox(asthma, "incident_asthma", "disease_event", "time", rhs_base, "cox_cause_specific_death_censored"),
    fit_cox(asthma, "incident_asthma", "disease_event", "time", rhs_frailty, "cox_cause_specific_plus_frailty"),
    fit_cox(asthma, "incident_asthma", "composite_event", "composite_time", rhs_base, "cox_composite_disease_or_death"),
    fit_cox(asthma, "incident_asthma", "disease_event", "time", rhs_raw_pef, "cox_raw_pef_sensitivity"),
    fit_cox(asthma, "incident_asthma", "disease_event", "time", rhs_base, "cox_weighted_cause_specific", weighted = TRUE),
    fit_cox(death, "death_w2_w4", "death_event", "time", rhs_base, "cox_death"),
    fit_cox(death, "death_w2_w4", "death_event", "time", rhs_frailty, "cox_death_plus_frailty"),
    fit_cox(death, "death_w2_w4", "death_event", "time", rhs_raw_pef, "cox_death_raw_pef_sensitivity"),
    fit_cox(death, "death_w2_w4", "death_event", "time", rhs_base, "cox_weighted_death", weighted = TRUE)
  )
  readr::write_csv(cox_results, file.path(table_dir, "charls_v0_2_cox_model_table.csv"))

  lung_period <- make_period_data(lung, "incident_chronic_lung_disease")
  asthma_period <- make_period_data(asthma, "incident_asthma")
  discrete_results <- bind_rows(
    fit_discrete_time(lung_period, "incident_chronic_lung_disease", rhs_base, "discrete_time_death_censored"),
    fit_discrete_time(lung_period, "incident_chronic_lung_disease", rhs_frailty, "discrete_time_plus_frailty"),
    fit_discrete_time(lung_period, "incident_chronic_lung_disease", rhs_base, "discrete_time_weighted", weighted = TRUE),
    fit_discrete_time(asthma_period, "incident_asthma", rhs_base, "discrete_time_death_censored"),
    fit_discrete_time(asthma_period, "incident_asthma", rhs_frailty, "discrete_time_plus_frailty"),
    fit_discrete_time(asthma_period, "incident_asthma", rhs_base, "discrete_time_weighted", weighted = TRUE)
  )
  readr::write_csv(discrete_results, file.path(table_dir, "charls_v0_2_discrete_time_model_table.csv"))

  cohort_counts <- bind_rows(
    cohort_count_table(lung, "incident_chronic_lung_disease"),
    cohort_count_table(asthma, "incident_asthma"),
    cohort_count_table(death, "death_w2_w4")
  )
  readr::write_csv(cohort_counts, file.path(table_dir, "charls_v0_2_model_cohort_counts.csv"))

  wave_counts <- bind_rows(
    wave_count_table(lung, "incident_chronic_lung_disease"),
    wave_count_table(asthma, "incident_asthma")
  )
  readr::write_csv(wave_counts, file.path(table_dir, "charls_v0_2_wave_event_counts.csv"))

  competing_summary <- bind_rows(
    lung %>%
      summarise(
        outcome = "incident_chronic_lung_disease",
        n = n(),
        disease_events = sum(disease_event == 1, na.rm = TRUE),
        competing_deaths_before_disease = sum(competing_death == 1, na.rm = TRUE),
        composite_events = sum(composite_event == 1, na.rm = TRUE),
        censored_without_event_or_death = sum(disease_event == 0 & competing_death == 0, na.rm = TRUE)
      ),
    asthma %>%
      summarise(
        outcome = "incident_asthma",
        n = n(),
        disease_events = sum(disease_event == 1, na.rm = TRUE),
        competing_deaths_before_disease = sum(competing_death == 1, na.rm = TRUE),
        composite_events = sum(composite_event == 1, na.rm = TRUE),
        censored_without_event_or_death = sum(disease_event == 0 & competing_death == 0, na.rm = TRUE)
      )
  )
  readr::write_csv(competing_summary, file.path(table_dir, "charls_v0_2_competing_death_summary.csv"))

  key_terms <- bind_rows(cox_results, discrete_results) %>%
    filter(term %in% c("resp_vulnerability_z", "frailty_proxy_count_w1", "pef_per_100_l_min")) %>%
    mutate(
      effect_ci = sprintf("%.2f (%.2f, %.2f)", effect, effect_low, effect_high),
      event_rate = sprintf("%.1f%%", 100 * event_rate)
    ) %>%
    select(outcome, model, estimand, n_obs, events, event_rate, term, effect_ci, p_formatted)
  readr::write_csv(key_terms, file.path(table_dir, "charls_v0_2_key_terms.csv"))

  phenotype_candidates <- make_phenotype_lock_candidates()
  readr::write_csv(phenotype_candidates, file.path(metadata_dir, "charls_v0_2_phenotype_lock_candidates.csv"))

  forest_dat <- bind_rows(cox_results, discrete_results) %>%
    filter(
      term == "resp_vulnerability_z",
      model %in% c(
        "cox_cause_specific_death_censored",
        "cox_cause_specific_plus_frailty",
        "cox_composite_disease_or_death",
        "cox_weighted_cause_specific",
        "cox_death",
        "cox_death_plus_frailty",
        "cox_weighted_death",
        "discrete_time_death_censored",
        "discrete_time_plus_frailty",
        "discrete_time_weighted"
      )
    ) %>%
    mutate(
      outcome_label = recode(
        outcome,
        incident_chronic_lung_disease = "Incident chronic lung disease",
        incident_asthma = "Incident asthma",
        death_w2_w4 = "Death waves 2-4"
      ),
      model_label = recode(
        model,
        cox_cause_specific_death_censored = "Cox death-censored",
        cox_cause_specific_plus_frailty = "Cox + frailty",
        cox_composite_disease_or_death = "Cox disease/death composite",
        cox_weighted_cause_specific = "Cox weighted",
        cox_death = "Cox death",
        cox_death_plus_frailty = "Cox death + frailty",
        cox_weighted_death = "Cox death weighted",
        discrete_time_death_censored = "Discrete-time",
        discrete_time_plus_frailty = "Discrete-time + frailty",
        discrete_time_weighted = "Discrete-time weighted"
      ),
      plot_label = paste(outcome_label, model_label, sep = " - ")
    )

  forest_plot <- ggplot(forest_dat, aes(x = effect, y = reorder(plot_label, effect), color = outcome_label)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "gray55") +
    geom_errorbarh(aes(xmin = effect_low, xmax = effect_high), height = 0.16, linewidth = 0.45) +
    geom_point(size = 2.1) +
    scale_x_log10() +
    scale_color_manual(values = c(
      "Incident chronic lung disease" = "#0072B2",
      "Incident asthma" = "#D55E00",
      "Death waves 2-4" = "#009E73"
    )) +
    labs(
      x = "Effect per 1-SD higher respiratory vulnerability (log scale)",
      y = NULL,
      color = NULL
    ) +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "bottom"
    )
  ggplot2::ggsave(
    filename = file.path(figure_dir, "charls_v0_2_resp_vulnerability_forest.png"),
    plot = forest_plot,
    width = 8.4,
    height = 6.2,
    dpi = 180,
    bg = "white"
  )

  disease_curve <- bind_rows(
    lung %>%
      mutate(vulnerability_quartile = factor(dplyr::ntile(resp_vulnerability_z, 4), levels = 1:4, labels = c("Q1 least vulnerable", "Q2", "Q3", "Q4 most vulnerable"))) %>%
      group_by(vulnerability_quartile) %>%
      summarise(outcome = "Incident chronic lung disease", wave2 = mean(!is.na(disease_event_wave) & disease_event_wave <= 2), wave3 = mean(!is.na(disease_event_wave) & disease_event_wave <= 3), wave4 = mean(!is.na(disease_event_wave) & disease_event_wave <= 4), .groups = "drop"),
    asthma %>%
      mutate(vulnerability_quartile = factor(dplyr::ntile(resp_vulnerability_z, 4), levels = 1:4, labels = c("Q1 least vulnerable", "Q2", "Q3", "Q4 most vulnerable"))) %>%
      group_by(vulnerability_quartile) %>%
      summarise(outcome = "Incident asthma", wave2 = mean(!is.na(disease_event_wave) & disease_event_wave <= 2), wave3 = mean(!is.na(disease_event_wave) & disease_event_wave <= 3), wave4 = mean(!is.na(disease_event_wave) & disease_event_wave <= 4), .groups = "drop")
  ) %>%
    tidyr::pivot_longer(cols = starts_with("wave"), names_to = "wave", values_to = "cumulative_rate") %>%
    mutate(wave = recode(wave, wave2 = "Wave 2", wave3 = "Wave 3", wave4 = "Wave 4"))
  readr::write_csv(disease_curve, file.path(table_dir, "charls_v0_2_cumulative_event_by_wave.csv"))

  curve_plot <- ggplot(disease_curve, aes(x = wave, y = cumulative_rate, group = vulnerability_quartile, color = vulnerability_quartile)) +
    geom_line(linewidth = 0.75) +
    geom_point(size = 1.8) +
    facet_wrap(~outcome, scales = "free_y") +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    scale_color_manual(values = c("#56B4E9", "#009E73", "#E69F00", "#D55E00")) +
    labs(x = NULL, y = "Cumulative observed event rate", color = NULL) +
    theme_minimal(base_size = 10) +
    theme(panel.grid.minor = element_blank(), legend.position = "bottom")
  ggplot2::ggsave(
    filename = file.path(figure_dir, "charls_v0_2_cumulative_event_by_wave.png"),
    plot = curve_plot,
    width = 8.2,
    height = 4.8,
    dpi = 180,
    bg = "white"
  )

  key_resp <- key_terms %>%
    filter(term == "resp_vulnerability_z") %>%
    filter(model %in% c("cox_cause_specific_death_censored", "cox_cause_specific_plus_frailty", "cox_weighted_cause_specific", "cox_death", "cox_death_plus_frailty", "cox_weighted_death", "discrete_time_death_censored", "discrete_time_plus_frailty"))

  log <- c(
    "# CHARLS V0.2 Longitudinal Robust Model Log",
    "",
    paste0("- Input RDS: ", input_rds),
    paste0("- Base analysis rows after core covariate filters: ", nrow(base)),
    paste0("- Lung time-to-event rows: ", nrow(lung)),
    paste0("- Asthma time-to-event rows: ", nrow(asthma)),
    paste0("- Death time-to-event rows: ", nrow(death)),
    "- Formal `survey` package design was not available in this environment; weighted sensitivity uses normalized positive biomarker weights.",
    "- Time is approximated by follow-up wave interval because exact interview/death dates are not locked in the v0.2 dataset.",
    "- Death before disease is treated as censoring in cause-specific models and as an event in composite disease/death sensitivity models.",
    "",
    "## Key Respiratory Vulnerability Terms",
    "",
    paste(capture.output(print(key_resp, n = Inf, width = Inf)), collapse = "\n"),
    "",
    "## Outputs",
    "",
    "- `results/tables/charls_v0_2_cox_model_table.csv`",
    "- `results/tables/charls_v0_2_discrete_time_model_table.csv`",
    "- `results/tables/charls_v0_2_key_terms.csv`",
    "- `results/tables/charls_v0_2_model_cohort_counts.csv`",
    "- `results/tables/charls_v0_2_wave_event_counts.csv`",
    "- `results/tables/charls_v0_2_competing_death_summary.csv`",
    "- `results/tables/charls_v0_2_cumulative_event_by_wave.csv`",
    "- `metadata/charls_v0_2_phenotype_lock_candidates.csv`",
    "- `results/figures/charls_v0_2_resp_vulnerability_forest.png`",
    "- `results/figures/charls_v0_2_cumulative_event_by_wave.png`"
  )
  writeLines(log, file.path(log_dir, "charls_v0_2_longitudinal_model_log.md"))

  findings <- c(
    "# CHARLS V0.2 Findings",
    "",
    "This v0.2 analysis upgrades the first exploratory logistic models to wave-time Cox models, discrete-time hazard models, competing-death summaries, and phenotype-lock candidates. It remains provisional.",
    "",
    "## Model Frame",
    "",
    paste0("- Incident chronic lung disease cohort: N = ", nrow(lung), "; disease events = ", sum(lung$disease_event == 1, na.rm = TRUE), "; competing deaths before disease = ", sum(lung$competing_death == 1, na.rm = TRUE), "."),
    paste0("- Incident asthma cohort: N = ", nrow(asthma), "; disease events = ", sum(asthma$disease_event == 1, na.rm = TRUE), "; competing deaths before disease = ", sum(asthma$competing_death == 1, na.rm = TRUE), "."),
    paste0("- Death cohort: N = ", nrow(death), "; deaths = ", sum(death$death_event == 1, na.rm = TRUE), "."),
    "",
    "## Main Readout",
    "",
    "The key v0.2 tables are `charls_v0_2_key_terms.csv`, `charls_v0_2_cox_model_table.csv`, and `charls_v0_2_discrete_time_model_table.csv`. Interpret `resp_vulnerability_z` as higher lower-than-expected baseline PEF.",
    "",
    paste(capture.output(print(key_resp, n = Inf, width = Inf)), collapse = "\n"),
    "",
    "## Bottom Line",
    "",
    "The respiratory vulnerability signal persists when approximating event timing by follow-up wave and treating death before disease as competing censoring. The next methodological step is to lock official codebook rules, add exact interview/death timing if available, and move survey design from normalized-weight sensitivity to a formal PSU/strata design."
  )
  writeLines(findings, file.path(log_dir, "charls_v0_2_findings.md"))

  message("Wrote CHARLS v0.2 longitudinal robust model outputs.")
}

if (sys.nframe() == 0) {
  main()
}
