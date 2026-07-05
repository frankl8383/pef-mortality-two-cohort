# CHARLS v0.3 public-data survey-weighted models.
# Uses CHARLS public weights and communityID PSU clustering; writes aggregate outputs only.

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

prepare_base <- function(core) {
  core %>%
    mutate(
      resp_vulnerability_z = -pef_resid_z_w1,
      pef_per_100_l_min = pef_best_w1_valid_provisional / 100,
      age_decade = age_w1 / 10,
      sex = factor(sex_label),
      smoke_ever = factor(smoke_ever_w1, levels = c(0, 1), labels = c("no", "yes")),
      frailty_binary = factor(frailty_proxy_ge3_w1, levels = c(0, 1), labels = c("non_frail_proxy", "frail_proxy")),
      psu_community_id = community_id,
      county_pseudo_id = substr(community_id, 1, 6),
      analysis_weight = weight_biomarker_w1,
      death_wave = first_wave_match(c(5, 6), iwstat_w2, iwstat_w3, iwstat_w4),
      last_iw_wave = last_observed_wave(iwstat_w2, iwstat_w3, iwstat_w4)
    ) %>%
    filter(
      baseline_respiratory_axis_eligible,
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

effect_table <- function(fit, outcome, model, design_type, estimand, n_obs, events) {
  co <- stats::coef(fit)
  vc <- stats::vcov(fit)
  se <- sqrt(diag(vc))
  z <- co / se
  p <- 2 * stats::pnorm(abs(z), lower.tail = FALSE)
  tibble(
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
      outcome, model, design_type, estimand, n_obs, events, event_rate,
      term, estimate, std.error, statistic, p.value, p_formatted,
      effect, effect_low, effect_high
    )
}

fit_svycox <- function(dat, outcome, event_var, time_var, rhs, model, design_type) {
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
    outcome = outcome,
    model = model,
    design_type = design_type,
    warning = fit_result$warnings
  )
  list(
    table = effect_table(fit_result$value, outcome, model, design_type, "survey_HR", n_obs, events),
    warnings = warn_tbl
  )
}

fit_svyglm_discrete <- function(period_dat, outcome, rhs, model, design_type) {
  formula <- stats::as.formula(paste0("period_event ~ ", rhs, " + factor(interval)"))
  needed <- unique(c(all.vars(formula), "psu_community_id", "county_pseudo_id", "analysis_weight"))
  analysis <- period_dat %>%
    select(any_of(needed)) %>%
    filter(if_all(everything(), ~ !is.na(.x)))
  n_obs <- nrow(analysis)
  events <- sum(analysis$period_event == 1, na.rm = TRUE)
  if (n_obs < 50 || events < 10) {
    return(list(table = tibble(), warnings = tibble()))
  }
  design <- build_design(analysis, design_type)
  fit_result <- capture_warnings(
    survey::svyglm(formula, design = design, family = stats::quasibinomial())
  )
  warn_tbl <- tibble(
    outcome = outcome,
    model = model,
    design_type = design_type,
    warning = fit_result$warnings
  )
  list(
    table = effect_table(fit_result$value, outcome, model, design_type, "survey_OR", n_obs, events),
    warnings = warn_tbl
  )
}

design_diagnostics <- function(dat, label) {
  weight <- dat$analysis_weight
  tibble(
    analysis_set = label,
    n_rows = nrow(dat),
    unique_psu_community = dplyr::n_distinct(dat$psu_community_id),
    unique_county_pseudo = dplyr::n_distinct(dat$county_pseudo_id),
    min_weight = min(weight, na.rm = TRUE),
    median_weight = stats::median(weight, na.rm = TRUE),
    mean_weight = mean(weight, na.rm = TRUE),
    p99_weight = as.numeric(stats::quantile(weight, 0.99, na.rm = TRUE, names = FALSE)),
    max_weight = max(weight, na.rm = TRUE),
    kish_effective_n = (sum(weight, na.rm = TRUE)^2) / sum(weight^2, na.rm = TRUE)
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
    unique_psu_community = dplyr::n_distinct(dat$psu_community_id),
    unique_county_pseudo = dplyr::n_distinct(dat$county_pseudo_id),
    weighted_n = sum(dat$analysis_weight, na.rm = TRUE)
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
  base <- prepare_base(core)
  lung <- prepare_time_to_event(base, "incident_chronic_lung_disease")
  asthma <- prepare_time_to_event(base, "incident_asthma")
  death <- prepare_death_time(base)
  lung_period <- make_period_data(lung, "incident_chronic_lung_disease")
  asthma_period <- make_period_data(asthma, "incident_asthma")

  rhs_base <- "resp_vulnerability_z + age_decade + sex + smoke_ever + bmi_w1"
  rhs_frailty <- "resp_vulnerability_z + frailty_proxy_count_w1 + age_decade + sex + smoke_ever + bmi_w1"
  rhs_raw_pef <- "pef_per_100_l_min + age_decade + sex + smoke_ever + bmi_w1 + height_m_w1"
  design_types <- c("public_psu", "public_county_community")

  model_results <- list()
  warning_results <- list()
  add_result <- function(res) {
    model_results[[length(model_results) + 1]] <<- res$table
    warning_results[[length(warning_results) + 1]] <<- res$warnings
  }

  for (design_type in design_types) {
    add_result(fit_svycox(lung, "incident_chronic_lung_disease", "disease_event", "time", rhs_base, "svycox_cause_specific", design_type))
    add_result(fit_svycox(lung, "incident_chronic_lung_disease", "disease_event", "time", rhs_frailty, "svycox_plus_frailty", design_type))
    add_result(fit_svycox(lung, "incident_chronic_lung_disease", "composite_event", "composite_time", rhs_base, "svycox_composite_disease_or_death", design_type))
    add_result(fit_svycox(lung, "incident_chronic_lung_disease", "disease_event", "time", rhs_raw_pef, "svycox_raw_pef_sensitivity", design_type))
    add_result(fit_svycox(asthma, "incident_asthma", "disease_event", "time", rhs_base, "svycox_cause_specific", design_type))
    add_result(fit_svycox(asthma, "incident_asthma", "disease_event", "time", rhs_frailty, "svycox_plus_frailty", design_type))
    add_result(fit_svycox(asthma, "incident_asthma", "composite_event", "composite_time", rhs_base, "svycox_composite_disease_or_death", design_type))
    add_result(fit_svycox(asthma, "incident_asthma", "disease_event", "time", rhs_raw_pef, "svycox_raw_pef_sensitivity", design_type))
    add_result(fit_svycox(death, "death_w2_w4", "death_event", "time", rhs_base, "svycox_death", design_type))
    add_result(fit_svycox(death, "death_w2_w4", "death_event", "time", rhs_frailty, "svycox_death_plus_frailty", design_type))
    add_result(fit_svycox(death, "death_w2_w4", "death_event", "time", rhs_raw_pef, "svycox_death_raw_pef_sensitivity", design_type))

    add_result(fit_svyglm_discrete(lung_period, "incident_chronic_lung_disease", rhs_base, "svyglm_discrete_time", design_type))
    add_result(fit_svyglm_discrete(lung_period, "incident_chronic_lung_disease", rhs_frailty, "svyglm_discrete_time_plus_frailty", design_type))
    add_result(fit_svyglm_discrete(asthma_period, "incident_asthma", rhs_base, "svyglm_discrete_time", design_type))
    add_result(fit_svyglm_discrete(asthma_period, "incident_asthma", rhs_frailty, "svyglm_discrete_time_plus_frailty", design_type))
  }

  survey_table <- bind_rows(model_results)
  warnings_table <- bind_rows(warning_results)
  if (nrow(warnings_table) == 0) {
    warnings_table <- tibble(outcome = character(), model = character(), design_type = character(), warning = character())
  }

  readr::write_csv(survey_table, file.path(table_dir, "charls_v0_3_survey_model_table.csv"))
  readr::write_csv(warnings_table, file.path(table_dir, "charls_v0_3_survey_model_warnings.csv"))

  key_terms <- survey_table %>%
    filter(term %in% c("resp_vulnerability_z", "frailty_proxy_count_w1", "pef_per_100_l_min")) %>%
    mutate(
      effect_ci = sprintf("%.2f (%.2f, %.2f)", effect, effect_low, effect_high),
      event_rate = sprintf("%.1f%%", 100 * event_rate)
    ) %>%
    select(outcome, model, design_type, estimand, n_obs, events, event_rate, term, effect_ci, p_formatted)
  readr::write_csv(key_terms, file.path(table_dir, "charls_v0_3_survey_key_terms.csv"))

  cohort_counts <- bind_rows(
    cohort_count_table(lung, "incident_chronic_lung_disease"),
    cohort_count_table(asthma, "incident_asthma"),
    cohort_count_table(death, "death_w2_w4")
  )
  readr::write_csv(cohort_counts, file.path(table_dir, "charls_v0_3_survey_cohort_counts.csv"))

  diagnostics <- bind_rows(
    design_diagnostics(base, "base_after_core_filters"),
    design_diagnostics(lung, "lung_time_to_event"),
    design_diagnostics(asthma, "asthma_time_to_event"),
    design_diagnostics(death, "death_time_to_event"),
    design_diagnostics(lung_period, "lung_person_period"),
    design_diagnostics(asthma_period, "asthma_person_period")
  )
  readr::write_csv(diagnostics, file.path(table_dir, "charls_v0_3_survey_design_diagnostics.csv"))

  forest_dat <- survey_table %>%
    filter(
      term == "resp_vulnerability_z",
      model %in% c("svycox_cause_specific", "svycox_plus_frailty", "svycox_composite_disease_or_death", "svycox_death", "svycox_death_plus_frailty", "svyglm_discrete_time", "svyglm_discrete_time_plus_frailty")
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
        svycox_cause_specific = "svycox cause-specific",
        svycox_plus_frailty = "svycox + frailty",
        svycox_composite_disease_or_death = "svycox disease/death composite",
        svycox_death = "svycox death",
        svycox_death_plus_frailty = "svycox death + frailty",
        svyglm_discrete_time = "svyglm discrete-time",
        svyglm_discrete_time_plus_frailty = "svyglm discrete-time + frailty"
      ),
      design_label = recode(
        design_type,
        public_psu = "community PSU",
        public_county_community = "county/community"
      ),
      plot_label = paste(outcome_label, model_label, design_label, sep = " - ")
    )

  forest_plot <- ggplot(forest_dat, aes(x = effect, y = reorder(plot_label, effect), color = outcome_label)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "gray55") +
    geom_errorbarh(aes(xmin = effect_low, xmax = effect_high), height = 0.14, linewidth = 0.42) +
    geom_point(size = 1.9) +
    scale_x_log10() +
    scale_color_manual(values = c(
      "Incident chronic lung disease" = "#0072B2",
      "Incident asthma" = "#D55E00",
      "Death waves 2-4" = "#009E73"
    )) +
    labs(
      x = "Survey-weighted effect per 1-SD higher respiratory vulnerability (log scale)",
      y = NULL,
      color = NULL
    ) +
    theme_minimal(base_size = 9.5) +
    theme(panel.grid.minor = element_blank(), legend.position = "bottom")
  ggplot2::ggsave(
    filename = file.path(figure_dir, "charls_v0_3_survey_resp_vulnerability_forest.png"),
    plot = forest_plot,
    width = 9.2,
    height = 7.8,
    dpi = 180,
    bg = "white"
  )

  key_resp <- key_terms %>%
    filter(term == "resp_vulnerability_z") %>%
    filter(model %in% c(
      "svycox_cause_specific",
      "svycox_plus_frailty",
      "svycox_composite_disease_or_death",
      "svycox_death",
      "svycox_death_plus_frailty",
      "svyglm_discrete_time",
      "svyglm_discrete_time_plus_frailty"
    )) %>%
    arrange(design_type, outcome, model)

  format_key <- function(outcome_id, model_id, design_id = "public_psu") {
    row <- key_terms %>%
      filter(
        .data$outcome == .env$outcome_id,
        .data$model == .env$model_id,
        .data$design_type == .env$design_id,
        .data$term == "resp_vulnerability_z"
      ) %>%
      slice(1)
    if (nrow(row) == 0) {
      return("not estimated")
    }
    paste0(row$effect_ci, ", p=", row$p_formatted)
  }

  findings <- c(
    "# CHARLS V0.3 Survey-Weighted Findings",
    "",
    "This v0.3 analysis implements public-data CHARLS survey designs using `survey` 4.5. The primary design uses `community_id` as the public PSU and CHARLS biomarker weights. The sensitivity design uses the inferred `county_pseudo_id + community_id` hierarchy. No explicit strata are specified because no public strata variable was identified.",
    "",
    "## Design Frame",
    "",
    paste0("- Base analysis rows after complete-case and positive-weight filters: ", nrow(base), "."),
    paste0("- Lung time-to-event rows: ", nrow(lung), "; disease events: ", sum(lung$disease_event == 1, na.rm = TRUE), "."),
    paste0("- Asthma time-to-event rows: ", nrow(asthma), "; disease events: ", sum(asthma$disease_event == 1, na.rm = TRUE), "."),
    paste0("- Death rows: ", nrow(death), "; deaths: ", sum(death$death_event == 1, na.rm = TRUE), "."),
    paste0("- Unique community PSUs in base frame: ", dplyr::n_distinct(base$psu_community_id), "."),
    paste0("- Unique inferred county IDs in base frame: ", dplyr::n_distinct(base$county_pseudo_id), "."),
    "",
    "## Primary Public-PSU Results",
    "",
    paste0("- Chronic lung disease cause-specific `svycoxph`: ", format_key("incident_chronic_lung_disease", "svycox_cause_specific"), "."),
    paste0("- Chronic lung disease `svycoxph` plus frailty proxy: ", format_key("incident_chronic_lung_disease", "svycox_plus_frailty"), "."),
    paste0("- Chronic lung disease/death composite `svycoxph`: ", format_key("incident_chronic_lung_disease", "svycox_composite_disease_or_death"), "."),
    paste0("- Asthma cause-specific `svycoxph`: ", format_key("incident_asthma", "svycox_cause_specific"), "."),
    paste0("- Asthma/death composite `svycoxph`: ", format_key("incident_asthma", "svycox_composite_disease_or_death"), "."),
    paste0("- Death waves 2-4 `svycoxph`: ", format_key("death_w2_w4", "svycox_death"), "."),
    paste0("- Death waves 2-4 `svycoxph` plus frailty proxy: ", format_key("death_w2_w4", "svycox_death_plus_frailty"), "."),
    "",
    "## Key Respiratory Vulnerability Terms",
    "",
    markdown_table(key_resp),
    "",
    "## Interpretation",
    "",
    "- `resp_vulnerability_z` is higher when baseline PEF is lower than expected for age, sex, and height.",
    "- `survey_HR` rows are from `svycoxph`; `survey_OR` rows are from discrete-time `svyglm`.",
    "- These are public-data survey designs. Explicit GDP/region strata are not reconstructed."
  )
  writeLines(findings, file.path(log_dir, "charls_v0_3_survey_findings.md"))

  log <- c(
    "# CHARLS V0.3 Survey-Weighted Model Log",
    "",
    paste0("- Input RDS: ", input_rds),
    "- Primary design: ids = ~psu_community_id, weights = ~analysis_weight.",
    "- Sensitivity design: ids = ~county_pseudo_id + psu_community_id, weights = ~analysis_weight.",
    "- `survey.lonely.psu` set to `adjust`.",
    "- No explicit strata variable used; public files do not expose a clear official strata field.",
    "",
    "## Outputs",
    "",
    "- `results/tables/charls_v0_3_survey_model_table.csv`",
    "- `results/tables/charls_v0_3_survey_key_terms.csv`",
    "- `results/tables/charls_v0_3_survey_cohort_counts.csv`",
    "- `results/tables/charls_v0_3_survey_design_diagnostics.csv`",
    "- `results/tables/charls_v0_3_survey_model_warnings.csv`",
    "- `results/figures/charls_v0_3_survey_resp_vulnerability_forest.png`",
    "- `results/logs/charls_v0_3_survey_findings.md`"
  )
  writeLines(log, file.path(log_dir, "charls_v0_3_survey_model_log.md"))

  message("Wrote CHARLS v0.3 survey-weighted model outputs.")
}

if (sys.nframe() == 0) {
  main()
}
