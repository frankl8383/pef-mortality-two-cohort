# CHARLS exploratory respiratory vulnerability models.
# Reads local-only provisional RDS and writes aggregate tables/figures only.

suppressPackageStartupMessages({
  library(broom)
  library(dplyr)
  library(ggplot2)
  library(readr)
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

wald_or_table <- function(fit, outcome, model, n_obs, events) {
  broom::tidy(fit) %>%
    mutate(
      outcome = outcome,
      model = model,
      n_obs = n_obs,
      events = events,
      event_rate = events / n_obs,
      conf.low = estimate - 1.96 * std.error,
      conf.high = estimate + 1.96 * std.error,
      or = exp(estimate),
      or_low = exp(conf.low),
      or_high = exp(conf.high),
      p_formatted = format_p(p.value)
    ) %>%
    select(
      outcome, model, n_obs, events, event_rate,
      term, estimate, std.error, statistic, p.value, p_formatted,
      or, or_low, or_high
    )
}

fit_logistic <- function(dat, outcome, event_var, formula, model, weighted = FALSE) {
  model_vars <- all.vars(formula)
  needed <- unique(c(event_var, model_vars, if (weighted) "analysis_weight" else character()))
  analysis <- dat %>%
    filter(!is.na(.data[[event_var]])) %>%
    select(any_of(needed)) %>%
    filter(if_all(everything(), ~ !is.na(.x)))

  if (weighted) {
    analysis <- analysis %>% filter(analysis_weight > 0)
  }
  n_obs <- nrow(analysis)
  events <- sum(analysis[[event_var]] == 1, na.rm = TRUE)
  if (n_obs < 50 || events < 10 || events >= n_obs - 10) {
    return(tibble(
      outcome = outcome,
      model = model,
      n_obs = n_obs,
      events = events,
      event_rate = ifelse(n_obs > 0, events / n_obs, NA_real_),
      term = NA_character_,
      estimate = NA_real_,
      std.error = NA_real_,
      statistic = NA_real_,
      p.value = NA_real_,
      p_formatted = NA_character_,
      or = NA_real_,
      or_low = NA_real_,
      or_high = NA_real_
    ))
  }

  if (weighted) {
    fit <- stats::glm(
      formula,
      data = analysis,
      family = stats::quasibinomial(),
      weights = analysis_weight
    )
  } else {
    fit <- stats::glm(
      formula,
      data = analysis,
      family = stats::binomial()
    )
  }
  wald_or_table(fit, outcome, model, n_obs, events)
}

event_rate_by_group <- function(dat, outcome, event_var, group_vars) {
  dat %>%
    filter(!is.na(.data[[event_var]])) %>%
    group_by(across(all_of(group_vars))) %>%
    summarise(
      outcome = outcome,
      n = n(),
      events = sum(.data[[event_var]] == 1, na.rm = TRUE),
      event_rate = events / n,
      .groups = "drop"
    ) %>%
    mutate(
      se = sqrt(event_rate * (1 - event_rate) / n),
      rate_low = pmax(0, event_rate - 1.96 * se),
      rate_high = pmin(1, event_rate + 1.96 * se)
    )
}

prep_outcome_data <- function(core, outcome) {
  if (outcome == "incident_chronic_lung_disease") {
    event_var <- "incident_lung_w2_w4"
    baseline_var <- "lung_w1"
  } else if (outcome == "incident_asthma") {
    event_var <- "incident_asthma_w2_w4"
    baseline_var <- "asthma_w1"
  } else if (outcome == "death_w2_w4") {
    event_var <- "died_w2_w4"
    baseline_var <- NA_character_
  } else {
    stop("Unsupported outcome: ", outcome, call. = FALSE)
  }

  dat <- core %>%
    mutate(
      resp_vulnerability_z = -pef_resid_z_w1,
      pef_per_100_l_min = pef_best_w1_valid_provisional / 100,
      age_decade = age_w1 / 10,
      sex = factor(sex_label),
      smoke_ever = factor(smoke_ever_w1, levels = c(0, 1), labels = c("no", "yes")),
      frailty_binary = factor(frailty_proxy_ge3_w1, levels = c(0, 1), labels = c("non_frail_proxy", "frail_proxy")),
      analysis_weight_raw = weight_biomarker_w1
    ) %>%
    filter(
      baseline_respiratory_axis_eligible,
      !is.na(resp_vulnerability_z),
      !is.na(age_decade),
      !is.na(sex),
      !is.na(bmi_w1),
      !is.na(smoke_ever)
    )

  if (!is.na(baseline_var)) {
    dat <- dat %>% filter(.data[[baseline_var]] == 0)
  }
  dat <- dat %>% filter(!is.na(.data[[event_var]]))
  mean_positive_weight <- mean(dat$analysis_weight_raw[dat$analysis_weight_raw > 0], na.rm = TRUE)
  dat <- dat %>%
    mutate(
      analysis_weight = ifelse(
        !is.na(analysis_weight_raw) & analysis_weight_raw > 0 & is.finite(mean_positive_weight),
        analysis_weight_raw / mean_positive_weight,
        NA_real_
      )
    )

  dat <- dat %>%
    mutate(
      vulnerability_quartile = dplyr::ntile(resp_vulnerability_z, 4),
      vulnerability_quartile = factor(
        vulnerability_quartile,
        levels = 1:4,
        labels = c("Q1 least vulnerable", "Q2", "Q3", "Q4 most vulnerable")
      )
    )

  list(data = dat, event_var = event_var)
}

fit_outcome_models <- function(core, outcome) {
  prepared <- prep_outcome_data(core, outcome)
  dat <- prepared$data
  event_var <- prepared$event_var
  event_symbol <- stats::as.formula(paste0(event_var, " ~ resp_vulnerability_z + age_decade + sex + smoke_ever + bmi_w1"))
  frailty_formula <- stats::as.formula(paste0(event_var, " ~ resp_vulnerability_z + frailty_proxy_count_w1 + age_decade + sex + smoke_ever + bmi_w1"))
  interaction_formula <- stats::as.formula(paste0(event_var, " ~ resp_vulnerability_z * frailty_binary + age_decade + sex + smoke_ever + bmi_w1"))
  raw_pef_formula <- stats::as.formula(paste0(event_var, " ~ pef_per_100_l_min + age_decade + sex + height_m_w1 + smoke_ever + bmi_w1"))

  dplyr::bind_rows(
    fit_logistic(dat, outcome, event_var, event_symbol, "unweighted_residual_pef_base"),
    fit_logistic(dat, outcome, event_var, frailty_formula, "unweighted_residual_pef_plus_frailty_count"),
    fit_logistic(dat, outcome, event_var, interaction_formula, "unweighted_residual_pef_x_frailty_binary"),
    fit_logistic(dat, outcome, event_var, raw_pef_formula, "unweighted_raw_pef_per_100_sensitivity"),
    fit_logistic(dat, outcome, event_var, event_symbol, "biomarker_weighted_residual_pef_base", weighted = TRUE)
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
  outcomes <- c("incident_chronic_lung_disease", "incident_asthma", "death_w2_w4")

  model_table <- dplyr::bind_rows(lapply(outcomes, function(outcome) {
    fit_outcome_models(core, outcome)
  }))
  readr::write_csv(model_table, file.path(table_dir, "charls_exploratory_model_table.csv"))

  rate_tables <- list()
  frailty_rate_tables <- list()
  cohorts <- list()
  for (outcome in outcomes) {
    prepared <- prep_outcome_data(core, outcome)
    dat <- prepared$data
    event_var <- prepared$event_var
    cohorts[[outcome]] <- tibble::tibble(
      outcome = outcome,
      n = nrow(dat),
      events = sum(dat[[event_var]] == 1, na.rm = TRUE),
      event_rate = events / n,
      pef_resid_z_mean = mean(dat$pef_resid_z_w1, na.rm = TRUE),
      resp_vulnerability_z_mean = mean(dat$resp_vulnerability_z, na.rm = TRUE),
      frailty_proxy_observed = sum(!is.na(dat$frailty_proxy_ge3_w1)),
      frailty_proxy_positive = sum(dat$frailty_proxy_ge3_w1 == 1, na.rm = TRUE)
    )
    rate_tables[[outcome]] <- event_rate_by_group(dat, outcome, event_var, "vulnerability_quartile")
    frailty_rate_tables[[outcome]] <- event_rate_by_group(
      dat %>% filter(!is.na(frailty_binary)),
      outcome,
      event_var,
      c("vulnerability_quartile", "frailty_binary")
    )
  }
  cohort_table <- dplyr::bind_rows(cohorts)
  event_rate_table <- dplyr::bind_rows(rate_tables)
  frailty_event_rate_table <- dplyr::bind_rows(frailty_rate_tables)

  readr::write_csv(cohort_table, file.path(table_dir, "charls_exploratory_cohort_counts.csv"))
  readr::write_csv(event_rate_table, file.path(table_dir, "charls_exploratory_event_rates_by_pef_quartile.csv"))
  readr::write_csv(frailty_event_rate_table, file.path(table_dir, "charls_exploratory_event_rates_by_pef_quartile_frailty.csv"))

  frailty_model_dat <- core %>%
    mutate(
      resp_vulnerability_z = -pef_resid_z_w1,
      age_decade = age_w1 / 10,
      sex = factor(sex_label),
      smoke_ever = factor(smoke_ever_w1, levels = c(0, 1), labels = c("no", "yes"))
    ) %>%
    filter(
      baseline_respiratory_axis_eligible,
      !is.na(frailty_proxy_ge3_w1),
      !is.na(resp_vulnerability_z),
      !is.na(age_decade),
      !is.na(sex),
      !is.na(smoke_ever),
      !is.na(bmi_w1)
    )
  frailty_fit <- stats::glm(
    frailty_proxy_ge3_w1 ~ resp_vulnerability_z + age_decade + sex + smoke_ever + bmi_w1,
    data = frailty_model_dat,
    family = stats::binomial()
  )
  frailty_model <- wald_or_table(
    frailty_fit,
    "baseline_frailty_proxy_ge3",
    "unweighted_residual_pef_base",
    nrow(frailty_model_dat),
    sum(frailty_model_dat$frailty_proxy_ge3_w1 == 1)
  )
  readr::write_csv(frailty_model, file.path(table_dir, "charls_exploratory_frailty_model_table.csv"))

  fig_dat <- event_rate_table %>%
    filter(outcome %in% c("incident_chronic_lung_disease", "incident_asthma"))
  p <- ggplot(fig_dat, aes(x = vulnerability_quartile, y = event_rate, group = outcome, color = outcome)) +
    geom_point(size = 2.2) +
    geom_line(linewidth = 0.8) +
    geom_errorbar(aes(ymin = rate_low, ymax = rate_high), width = 0.12, linewidth = 0.5) +
    scale_color_manual(
      values = c(
        incident_chronic_lung_disease = "#0072B2",
        incident_asthma = "#D55E00"
      ),
      labels = c(
        incident_chronic_lung_disease = "Incident chronic lung disease",
        incident_asthma = "Incident asthma"
      )
    ) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    labs(
      x = "Baseline respiratory vulnerability quartile (higher = lower-than-expected PEF)",
      y = "Incident event rate across follow-up waves 2-4",
      color = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "bottom"
    )
  ggplot2::ggsave(
    filename = file.path(figure_dir, "charls_exploratory_event_rate_by_pef_quartile.png"),
    plot = p,
    width = 7.2,
    height = 4.6,
    dpi = 180,
    bg = "white"
  )

  key_terms <- model_table %>%
    filter(term %in% c("resp_vulnerability_z", "pef_per_100_l_min", "frailty_proxy_count_w1")) %>%
    mutate(
      or_ci = sprintf("%.2f (%.2f, %.2f)", or, or_low, or_high),
      event_rate = sprintf("%.1f%%", 100 * event_rate)
    ) %>%
    select(outcome, model, n_obs, events, event_rate, term, or_ci, p_formatted)
  readr::write_csv(key_terms, file.path(table_dir, "charls_exploratory_key_terms.csv"))

  log <- c(
    "# CHARLS Exploratory Respiratory Vulnerability Models",
    "",
    paste0("- Input RDS: ", input_rds),
    paste0("- Rows in provisional core dataset: ", nrow(core)),
    paste0("- Outcomes modeled: ", paste(outcomes, collapse = ", ")),
    "",
    "## Interpretation Guardrails",
    "",
    "- These are provisional feasibility models, not publication-locked causal estimates.",
    "- `resp_vulnerability_z` is `-pef_resid_z_w1`, so higher values mean lower-than-expected baseline PEF after age, sex, and height adjustment.",
    "- Incident chronic lung disease/asthma require baseline no disease and any follow-up status equal to 1 in waves 2-4.",
    "- Weighted sensitivity uses the Harmonized wave-1 biomarker weight when available, with quasibinomial GLM.",
    "- All outputs are aggregate tables or figures; no row-level data are written outside `derived_sensitive/`.",
    "",
    "## Output Tables",
    "",
    "- `results/tables/charls_exploratory_model_table.csv`",
    "- `results/tables/charls_exploratory_key_terms.csv`",
    "- `results/tables/charls_exploratory_cohort_counts.csv`",
    "- `results/tables/charls_exploratory_event_rates_by_pef_quartile.csv`",
    "- `results/tables/charls_exploratory_event_rates_by_pef_quartile_frailty.csv`",
    "- `results/tables/charls_exploratory_frailty_model_table.csv`",
    "",
    "## Output Figure",
    "",
    "- `results/figures/charls_exploratory_event_rate_by_pef_quartile.png`"
  )
  writeLines(log, file.path(log_dir, "charls_exploratory_model_log.md"))

  message("Wrote CHARLS exploratory respiratory vulnerability model outputs.")
}

if (sys.nframe() == 0) {
  main()
}
