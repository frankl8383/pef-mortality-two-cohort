# NHANES V30 linked mortality models.
# Runs survey-weighted Cox models after the T1 linked mortality ETL.

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(splines)
  library(survey)
  library(survival)
  library(tibble)
})

options(survey.lonely.psu = "adjust")
version_id <- "v30_0"

find_project_root <- function() {
  script_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", script_args, value = TRUE)
  if (length(file_arg) > 0) {
    script_path <- sub("^--file=", "", file_arg[[1]])
    return(normalizePath(file.path(dirname(script_path), "..", ".."), mustWork = FALSE))
  }
  normalizePath(getwd(), mustWork = FALSE)
}

dir_create <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
}

format_p <- function(p) {
  ifelse(is.na(p), NA_character_, ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
}

format_ci <- function(est, lo, hi, digits = 2) {
  sprintf(paste0("%.", digits, "f (%.", digits, "f-%.", digits, "f)"), est, lo, hi)
}

markdown_table <- function(dat) {
  dat <- as.data.frame(dat, stringsAsFactors = FALSE)
  dat[] <- lapply(dat, function(x) ifelse(is.na(x), "", as.character(x)))
  header <- paste0("| ", paste(names(dat), collapse = " | "), " |")
  separator <- paste0("| ", paste(rep("---", ncol(dat)), collapse = " | "), " |")
  rows <- apply(dat, 1, function(x) paste0("| ", paste(x, collapse = " | "), " |"))
  c(header, separator, rows)
}

make_design <- function(data, weight_var) {
  survey::svydesign(
    ids = ~psu,
    strata = ~strata,
    weights = stats::as.formula(paste0("~", weight_var)),
    nest = TRUE,
    data = data
  )
}

prepare_base_data <- function(data) {
  data %>%
    filter(mortality_analysis_primary == 1L) %>%
    mutate(
      sex_f = factor(sex),
      race_f = factor(race_ethnicity),
      education_f = factor(education),
      smoking_f = factor(smoking_status, levels = c("never", "former", "current")),
      cycle_f = factor(cycle_label, levels = c("E", "F", "G")),
      rv_quartile_f = factor(paste0("Q", rv_quartile_num), levels = c("Q1", "Q2", "Q3", "Q4"))
    )
}

complete_model_data <- function(data, rhs_terms, event_var, weight_var, cycle_filter = NA_character_) {
  required <- unique(c(
    "followup_years_exam", event_var, all.vars(stats::as.formula(paste0("~", rhs_terms))),
    weight_var, "psu", "strata"
  ))
  missing_cols <- setdiff(required, names(data))
  if (length(missing_cols) > 0) {
    stop("Missing model variables: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }
  out <- data
  if (!is.na(cycle_filter)) {
    out <- out %>% filter(cycle_label == cycle_filter)
  }
  out <- out %>%
    filter(
      !is.na(.data[[weight_var]]),
      .data[[weight_var]] > 0,
      !is.na(psu),
      !is.na(strata),
      complete.cases(across(all_of(required)))
    )
  out
}

term_label <- function(term) {
  dplyr::case_when(
    term == "resp_vulnerability_z" ~ "Per 1-SD higher respiratory vulnerability marker",
    term == "resp_vulnerability_z_race_calibrated" ~ "Per 1-SD higher race-calibrated marker",
    term == "rv_quartile_num" ~ "Per higher vulnerability quartile",
    term == "rv_quartile_fQ2" ~ "Q2 versus Q1",
    term == "rv_quartile_fQ3" ~ "Q3 versus Q1",
    term == "rv_quartile_fQ4" ~ "Q4 versus Q1",
    TRUE ~ term
  )
}

extract_cox_rows <- function(fit, spec, df, terms_to_keep, warnings) {
  beta <- stats::coef(fit)
  vc <- stats::vcov(fit)
  se <- sqrt(diag(vc))
  terms_to_keep <- intersect(terms_to_keep, names(beta))
  if (length(terms_to_keep) == 0) {
    stop("No requested coefficient found for ", spec$analysis_id, call. = FALSE)
  }
  bind_rows(lapply(terms_to_keep, function(term) {
    est <- beta[[term]]
    se_term <- se[[term]]
    z <- est / se_term
    p <- 2 * stats::pnorm(abs(z), lower.tail = FALSE)
    ci <- est + c(-1, 1) * stats::qnorm(0.975) * se_term
    tibble(
      status = "ok",
      analysis_id = spec$analysis_id,
      analysis_family = spec$analysis_family,
      outcome = spec$outcome,
      outcome_label = spec$outcome_label,
      model_label = spec$model_label,
      exposure = spec$exposure,
      term = term,
      term_label = term_label(term),
      covariates = spec$covariates,
      weight_variable = spec$weight_variable,
      cycle = ifelse(is.na(spec$cycle), "pooled", spec$cycle),
      n = nrow(df),
      events = sum(df[[spec$outcome]] == 1L, na.rm = TRUE),
      person_years = sum(df$followup_years_exam, na.rm = TRUE),
      log_hr = est,
      se_log_hr = se_term,
      hr = exp(est),
      ci_low = exp(ci[[1]]),
      ci_high = exp(ci[[2]]),
      p_value = p,
      p_value_formatted = format_p(p),
      hr_ci = format_ci(exp(est), exp(ci[[1]]), exp(ci[[2]]), 2),
      warnings_n = length(warnings),
      warning = ifelse(length(warnings) == 0, NA_character_, paste(unique(warnings), collapse = " | ")),
      error = NA_character_
    )
  }))
}

fit_model_spec <- function(data, spec) {
  df <- complete_model_data(
    data = data,
    rhs_terms = paste(c(spec$exposure, spec$covariates), collapse = " + "),
    event_var = spec$outcome,
    weight_var = spec$weight_variable,
    cycle_filter = spec$cycle
  )
  events <- sum(df[[spec$outcome]] == 1L, na.rm = TRUE)
  if (nrow(df) < spec$minimum_n || events < spec$minimum_events) {
    stop(
      "Insufficient records/events: n=", nrow(df), ", events=", events,
      ", required n>=", spec$minimum_n, " and events>=", spec$minimum_events,
      call. = FALSE
    )
  }
  design <- make_design(df, spec$weight_variable)
  formula <- stats::as.formula(paste0(
    "survival::Surv(followup_years_exam, ", spec$outcome, ") ~ ",
    paste(c(spec$exposure, spec$covariates), collapse = " + ")
  ))
  warnings <- character()
  fit <- withCallingHandlers(
    survey::svycoxph(formula, design = design),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  terms_to_keep <- if (nzchar(spec$extract_regex)) {
    grep(spec$extract_regex, names(stats::coef(fit)), value = TRUE)
  } else {
    spec$exposure
  }
  rows <- extract_cox_rows(fit, spec, df, terms_to_keep, warnings)
  list(rows = rows, fit = fit, data = df)
}

safe_fit_model_spec <- function(data, spec) {
  tryCatch(
    fit_model_spec(data, spec),
    error = function(e) {
      list(
        rows = tibble(
          status = "error",
          analysis_id = spec$analysis_id,
          analysis_family = spec$analysis_family,
          outcome = spec$outcome,
          outcome_label = spec$outcome_label,
          model_label = spec$model_label,
          exposure = spec$exposure,
          term = ifelse(nzchar(spec$extract_regex), spec$extract_regex, spec$exposure),
          term_label = ifelse(nzchar(spec$extract_regex), spec$extract_regex, term_label(spec$exposure)),
          covariates = spec$covariates,
          weight_variable = spec$weight_variable,
          cycle = ifelse(is.na(spec$cycle), "pooled", spec$cycle),
          n = NA_integer_,
          events = NA_integer_,
          person_years = NA_real_,
          log_hr = NA_real_,
          se_log_hr = NA_real_,
          hr = NA_real_,
          ci_low = NA_real_,
          ci_high = NA_real_,
          p_value = NA_real_,
          p_value_formatted = NA_character_,
          hr_ci = NA_character_,
          warnings_n = 1L,
          warning = NA_character_,
          error = conditionMessage(e)
        ),
        fit = NULL,
        data = NULL
      )
    }
  )
}

rate_by_quartile <- function(data, outcome, outcome_label) {
  bind_rows(lapply(levels(data$rv_quartile_f), function(q) {
    df <- data %>%
      filter(rv_quartile_f == q, !is.na(.data[[outcome]]), !is.na(wtmec6yr), wtmec6yr > 0)
    if (nrow(df) < 2) {
      return(tibble(
        outcome = outcome,
        outcome_label = outcome_label,
        rv_quartile = q,
        n = nrow(df),
        events = sum(df[[outcome]] == 1L, na.rm = TRUE),
        weighted_rate = NA_real_,
        se = NA_real_,
        ci_low = NA_real_,
        ci_high = NA_real_
      ))
    }
    df$.event <- as.numeric(df[[outcome]] == 1L)
    design <- make_design(df, "wtmec6yr")
    est <- survey::svymean(~.event, design = design, na.rm = TRUE)
    rate <- as.numeric(stats::coef(est)[1])
    se <- as.numeric(survey::SE(est)[1])
    ci <- pmin(pmax(rate + c(-1, 1) * stats::qnorm(0.975) * se, 0), 1)
    tibble(
      outcome = outcome,
      outcome_label = outcome_label,
      rv_quartile = q,
      n = nrow(df),
      events = sum(df[[outcome]] == 1L, na.rm = TRUE),
      weighted_rate = rate,
      se = se,
      ci_low = ci[[1]],
      ci_high = ci[[2]]
    )
  }))
}

build_spline <- function(data) {
  rhs <- paste(
    "splines::ns(resp_vulnerability_z, df = 3)",
    "age_years",
    "sex_f",
    "race_f",
    "education_f",
    "income_poverty_ratio",
    "bmi",
    "smoking_f",
    sep = " + "
  )
  df <- complete_model_data(
    data,
    rhs_terms = rhs,
    event_var = "all_cause_death",
    weight_var = "wtmec6yr"
  )
  design <- make_design(df, "wtmec6yr")
  formula <- stats::as.formula(paste0("survival::Surv(followup_years_exam, all_cause_death) ~ ", rhs))
  fit <- survey::svycoxph(formula, design = design)

  x_grid <- seq(
    stats::quantile(df$resp_vulnerability_z, 0.01, na.rm = TRUE),
    stats::quantile(df$resp_vulnerability_z, 0.99, na.rm = TRUE),
    length.out = 120
  )
  mode_value <- function(x) {
    ux <- unique(x[!is.na(x)])
    ux[which.max(tabulate(match(x, ux)))]
  }
  ref <- tibble(
    resp_vulnerability_z = 0,
    age_years = stats::median(df$age_years, na.rm = TRUE),
    sex_f = factor(mode_value(df$sex_f), levels = levels(df$sex_f)),
    race_f = factor(mode_value(df$race_f), levels = levels(df$race_f)),
    education_f = factor(mode_value(df$education_f), levels = levels(df$education_f)),
    income_poverty_ratio = stats::median(df$income_poverty_ratio, na.rm = TRUE),
    bmi = stats::median(df$bmi, na.rm = TRUE),
    smoking_f = factor(mode_value(df$smoking_f), levels = levels(df$smoking_f))
  )
  grid <- ref[rep(1, length(x_grid)), , drop = FALSE]
  grid$resp_vulnerability_z <- x_grid

  terms_obj <- stats::delete.response(stats::terms(fit))
  mm <- stats::model.matrix(terms_obj, grid)
  mm_ref <- stats::model.matrix(terms_obj, ref)
  beta <- stats::coef(fit)
  vc <- stats::vcov(fit)
  mm <- mm[, names(beta), drop = FALSE]
  mm_ref <- mm_ref[, names(beta), drop = FALSE]
  diff <- sweep(mm, 2, mm_ref[1, ], "-")
  log_hr <- as.numeric(diff %*% beta)
  se <- sqrt(diag(diff %*% vc %*% t(diff)))
  out <- tibble(
    resp_vulnerability_z = x_grid,
    log_hr = log_hr,
    se_log_hr = se,
    hr = exp(log_hr),
    ci_low = exp(log_hr - stats::qnorm(0.975) * se),
    ci_high = exp(log_hr + stats::qnorm(0.975) * se),
    reference_resp_vulnerability_z = 0,
    n = nrow(df),
    events = sum(df$all_cause_death == 1L, na.rm = TRUE)
  )
  list(source = out, fit = fit, data = df)
}

plot_forest <- function(model_results, figure_dir) {
  plot_data <- model_results %>%
    filter(
      status == "ok",
      outcome == "all_cause_death",
      term %in% c("resp_vulnerability_z", "resp_vulnerability_z_race_calibrated", "rv_quartile_num", "rv_quartile_fQ4")
    ) %>%
    mutate(
      display_label = paste0(model_label, " | ", term_label),
      display_label = factor(display_label, levels = rev(unique(display_label)))
    )

  p <- ggplot(plot_data, aes(x = hr, y = display_label)) +
    geom_vline(xintercept = 1, linewidth = 0.35, colour = "#6B7280") +
    geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0.14, linewidth = 0.35, colour = "#374151") +
    geom_point(size = 1.7, colour = "#1F6F8B") +
    scale_x_log10() +
    labs(x = "Hazard ratio for all-cause mortality", y = NULL) +
    theme_classic(base_size = 8) +
    theme(
      axis.line = element_line(linewidth = 0.35),
      axis.ticks = element_line(linewidth = 0.35),
      panel.grid.major.x = element_line(linewidth = 0.2, colour = "#E5E7EB")
    )

  ggplot2::ggsave(file.path(figure_dir, paste0("nhanes_mortality_model_forest_", version_id, ".png")), p, width = 7.2, height = 4.8, dpi = 600, units = "in")
  ggplot2::ggsave(file.path(figure_dir, paste0("nhanes_mortality_model_forest_", version_id, ".pdf")), p, width = 7.2, height = 4.8, units = "in", device = grDevices::pdf, useDingbats = FALSE)
}

plot_spline <- function(spline_source, figure_dir) {
  p <- ggplot(spline_source, aes(x = resp_vulnerability_z, y = hr)) +
    geom_hline(yintercept = 1, linewidth = 0.35, colour = "#6B7280") +
    geom_ribbon(aes(ymin = ci_low, ymax = ci_high), fill = "#BBD7E5", alpha = 0.7) +
    geom_line(linewidth = 0.55, colour = "#1F6F8B") +
    labs(
      x = "Respiratory vulnerability marker (SD)",
      y = "Hazard ratio for all-cause mortality"
    ) +
    theme_classic(base_size = 8) +
    theme(
      axis.line = element_line(linewidth = 0.35),
      axis.ticks = element_line(linewidth = 0.35),
      panel.grid.major.y = element_line(linewidth = 0.2, colour = "#E5E7EB")
    )
  ggplot2::ggsave(file.path(figure_dir, paste0("nhanes_mortality_spline_", version_id, ".png")), p, width = 4.8, height = 3.6, dpi = 600, units = "in")
  ggplot2::ggsave(file.path(figure_dir, paste0("nhanes_mortality_spline_", version_id, ".pdf")), p, width = 4.8, height = 3.6, units = "in", device = grDevices::pdf, useDingbats = FALSE)
}

plot_rates <- function(rate_table, figure_dir) {
  plot_data <- rate_table %>%
    mutate(
      weighted_rate_percent = 100 * weighted_rate,
      ci_low_percent = 100 * ci_low,
      ci_high_percent = 100 * ci_high,
      outcome_label = factor(outcome_label, levels = c("All-cause mortality", "CLRD mortality (exploratory)"))
    )
  p <- ggplot(plot_data, aes(x = rv_quartile, y = weighted_rate_percent, group = outcome_label)) +
    geom_errorbar(aes(ymin = ci_low_percent, ymax = ci_high_percent), width = 0.12, linewidth = 0.35, colour = "#4B5563") +
    geom_line(linewidth = 0.45, colour = "#1F6F8B") +
    geom_point(size = 1.65, colour = "#1F6F8B") +
    facet_wrap(~outcome_label, scales = "free_y", ncol = 2) +
    labs(x = "Respiratory vulnerability quartile", y = "Weighted mortality rate (%)") +
    theme_classic(base_size = 8) +
    theme(
      axis.line = element_line(linewidth = 0.35),
      axis.ticks = element_line(linewidth = 0.35),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      panel.grid.major.y = element_line(linewidth = 0.2, colour = "#E5E7EB")
    )
  ggplot2::ggsave(file.path(figure_dir, paste0("nhanes_mortality_rates_by_quartile_", version_id, ".png")), p, width = 6.4, height = 3.4, dpi = 600, units = "in")
  ggplot2::ggsave(file.path(figure_dir, paste0("nhanes_mortality_rates_by_quartile_", version_id, ".pdf")), p, width = 6.4, height = 3.4, units = "in", device = grDevices::pdf, useDingbats = FALSE)
}

root <- find_project_root()
derived_dir <- file.path(root, "derived_sensitive", "nhanes")
results_dir <- file.path(root, "results", "tables")
logs_dir <- file.path(root, "results", "logs")
figures_dir <- file.path(root, "results", "figures")
manuscript_dir <- file.path(root, "manuscript")
dir_create(results_dir)
dir_create(logs_dir)
dir_create(figures_dir)
dir_create(manuscript_dir)

input_rds <- file.path(derived_dir, "nhanes_mortality_analysis_ready_v30_0.rds")
if (!file.exists(input_rds)) {
  stop("Missing T1 NHANES mortality analysis-ready RDS. Run 08_nhanes_v30_linked_mortality_etl.R first.", call. = FALSE)
}
analysis_ready <- readRDS(input_rds)
base_data <- prepare_base_data(analysis_ready)

model_specs <- tibble::tribble(
  ~analysis_id, ~analysis_family, ~outcome, ~outcome_label, ~model_label, ~exposure, ~extract_regex, ~covariates, ~weight_variable, ~cycle, ~minimum_n, ~minimum_events,
  "t2_allcause_m0_demographic", "primary_continuous", "all_cause_death", "All-cause mortality", "M0 demographic", "resp_vulnerability_z", "", "age_years + sex_f + race_f", "wtmec6yr", NA_character_, 1000L, 100L,
  "t2_allcause_m1_clinical", "primary_continuous", "all_cause_death", "All-cause mortality", "M1 clinical", "resp_vulnerability_z", "", "age_years + sex_f + race_f + bmi + smoking_f", "wtmec6yr", NA_character_, 1000L, 100L,
  "t2_allcause_m2_primary", "primary_continuous", "all_cause_death", "All-cause mortality", "M2 primary socioeconomic", "resp_vulnerability_z", "", "age_years + sex_f + race_f + education_f + income_poverty_ratio + bmi + smoking_f", "wtmec6yr", NA_character_, 1000L, 100L,
  "t2_allcause_m3_frailty_adjusted", "frailty_sensitivity", "all_cause_death", "All-cause mortality", "M3 frailty-adjusted sensitivity", "resp_vulnerability_z", "", "age_years + sex_f + race_f + education_f + income_poverty_ratio + bmi + smoking_f + nhanes_frailty_proxy_count", "wtmec6yr", NA_character_, 1000L, 100L,
  "t2_allcause_m2_race_calibrated", "race_calibrated_sensitivity", "all_cause_death", "All-cause mortality", "M2 race-calibrated exposure", "resp_vulnerability_z_race_calibrated", "", "age_years + sex_f + race_f + education_f + income_poverty_ratio + bmi + smoking_f", "wtmec6yr", NA_character_, 1000L, 100L,
  "t2_allcause_m2_cycle_fixed", "cycle_fixed_sensitivity", "all_cause_death", "All-cause mortality", "M2 plus cycle fixed effects", "resp_vulnerability_z", "", "age_years + sex_f + race_f + education_f + income_poverty_ratio + bmi + smoking_f + cycle_f", "wtmec6yr", NA_character_, 1000L, 100L,
  "t2_allcause_quartile_trend", "quartile_model", "all_cause_death", "All-cause mortality", "M2 quartile trend", "rv_quartile_num", "", "age_years + sex_f + race_f + education_f + income_poverty_ratio + bmi + smoking_f", "wtmec6yr", NA_character_, 1000L, 100L,
  "t2_allcause_quartile_categories", "quartile_model", "all_cause_death", "All-cause mortality", "M2 quartile categories", "rv_quartile_f", "^rv_quartile_f", "age_years + sex_f + race_f + education_f + income_poverty_ratio + bmi + smoking_f", "wtmec6yr", NA_character_, 1000L, 100L,
  "t2_clrd_mini_exploratory", "respiratory_mortality_exploratory", "clrd_death", "CLRD mortality (exploratory)", "Exploratory age-sex model", "resp_vulnerability_z", "", "age_years + sex_f", "wtmec6yr", NA_character_, 1000L, 30L
)

cycle_specs <- tibble(
  cycle = c("E", "F", "G"),
  analysis_id = paste0("t2_allcause_cycle_", c("E", "F", "G")),
  analysis_family = "cycle_stratified_sensitivity",
  outcome = "all_cause_death",
  outcome_label = "All-cause mortality",
  model_label = paste0("Cycle ", c("E", "F", "G"), " sensitivity"),
  exposure = "resp_vulnerability_z",
  extract_regex = "",
  covariates = "age_years + sex_f + race_f + bmi + smoking_f",
  weight_variable = "wtmec2yr",
  minimum_n = 500L,
  minimum_events = 80L
) %>%
  transmute(
    analysis_id, analysis_family, outcome, outcome_label, model_label, exposure,
    extract_regex, covariates, weight_variable, cycle, minimum_n, minimum_events
  )
model_specs <- bind_rows(model_specs, cycle_specs)

fits <- lapply(seq_len(nrow(model_specs)), function(i) safe_fit_model_spec(base_data, model_specs[i, ]))
model_results <- bind_rows(lapply(fits, `[[`, "rows")) %>%
  mutate(
    hr = round(hr, 4),
    ci_low = round(ci_low, 4),
    ci_high = round(ci_high, 4),
    p_value = signif(p_value, 4),
    person_years = round(person_years, 1)
  )

model_cohorts <- model_results %>%
  group_by(analysis_id, analysis_family, outcome_label, model_label, weight_variable, cycle) %>%
  summarise(
    status = dplyr::first(status),
    n = dplyr::first(n),
    events = dplyr::first(events),
    person_years = dplyr::first(person_years),
    error = dplyr::first(error),
    .groups = "drop"
  )

rate_table <- bind_rows(
  rate_by_quartile(base_data, "all_cause_death", "All-cause mortality"),
  rate_by_quartile(base_data, "clrd_death", "CLRD mortality (exploratory)")
) %>%
  mutate(
    weighted_rate = round(weighted_rate, 5),
    se = round(se, 5),
    ci_low = round(ci_low, 5),
    ci_high = round(ci_high, 5)
  )

spline <- build_spline(base_data)
spline_source <- spline$source %>%
  mutate(
    log_hr = round(log_hr, 6),
    se_log_hr = round(se_log_hr, 6),
    hr = round(hr, 5),
    ci_low = round(ci_low, 5),
    ci_high = round(ci_high, 5)
  )

models_path <- file.path(results_dir, paste0("nhanes_mortality_models_", version_id, ".csv"))
cohorts_path <- file.path(results_dir, paste0("nhanes_mortality_model_cohort_counts_", version_id, ".csv"))
rates_path <- file.path(results_dir, paste0("nhanes_mortality_rates_by_quartile_", version_id, ".csv"))
spline_path <- file.path(results_dir, paste0("nhanes_mortality_spline_source_", version_id, ".csv"))
readr::write_csv(model_results, models_path)
readr::write_csv(model_cohorts, cohorts_path)
readr::write_csv(rate_table, rates_path)
readr::write_csv(spline_source, spline_path)

plot_forest(model_results, figures_dir)
plot_spline(spline_source, figures_dir)
plot_rates(rate_table, figures_dir)

limitations_path <- file.path(manuscript_dir, paste0("nhanes_mortality_limitations_", version_id, ".md"))
main_row <- model_results %>%
  filter(analysis_id == "t2_allcause_m2_primary", term == "resp_vulnerability_z") %>%
  slice(1)
q4_row <- model_results %>%
  filter(analysis_id == "t2_allcause_quartile_categories", term == "rv_quartile_fQ4") %>%
  slice(1)
clrd_row <- model_results %>%
  filter(analysis_id == "t2_clrd_mini_exploratory", term == "resp_vulnerability_z") %>%
  slice(1)

limitations_lines <- c(
  "# NHANES Linked Mortality Limitations",
  "",
  paste0("Version: ", version_id),
  "",
  "The NHANES linked mortality analysis strengthens the manuscript by adding a hard outcome, but it should remain observational and marker-focused. The primary estimand is an association between a lower-than-expected PEF marker and all-cause mortality, not a causal effect of respiratory vulnerability.",
  "",
  paste0("In the primary survey-weighted Cox model, each 1-SD higher marker was associated with all-cause mortality HR ", main_row$hr_ci, ", p=", main_row$p_value_formatted, ". Q4 versus Q1 in the quartile model was HR ", q4_row$hr_ci, ", p=", q4_row$p_value_formatted, "."),
  "",
  paste0("Chronic lower respiratory disease mortality is sparse in the primary cohort (41 events). The age-sex exploratory model gave HR ", clrd_row$hr_ci, ", p=", clrd_row$p_value_formatted, ", but this must be presented only as hypothesis-generating support, not as a definitive respiratory-specific mortality endpoint."),
  "",
  "Important limitations for manuscript wording:",
  "",
  "- The public-use NCHS linked mortality file follows participants through 2019 and includes public-use handling of selected follow-up and cause fields; all-cause vital status is the safest hard-outcome anchor.",
  "- The fully adjusted primary model uses complete cases for education, income-to-poverty ratio, BMI, and smoking status, so the modeled sample is smaller than the linked mortality cohort.",
  "- Frailty-adjusted models are sensitivity analyses because the NHANES frailty proxy is missing for many participants and may partly lie downstream of reduced respiratory reserve.",
  "- CLRD mortality has too few events for fully adjusted cause-specific inference.",
  "- The analysis supports a PEF-based marker of aging respiratory reserve and mortality risk; it does not establish a new causal axis or a clinical decision rule."
)
writeLines(limitations_lines, limitations_path)

log_path <- file.path(logs_dir, paste0("nhanes_mortality_models_", version_id, ".md"))
summary_rows <- model_results %>%
  filter(
    analysis_id %in% c(
      "t2_allcause_m2_primary",
      "t2_allcause_m3_frailty_adjusted",
      "t2_allcause_m2_race_calibrated",
      "t2_allcause_quartile_trend",
      "t2_allcause_quartile_categories",
      "t2_clrd_mini_exploratory"
    ),
    term %in% c("resp_vulnerability_z", "resp_vulnerability_z_race_calibrated", "rv_quartile_num", "rv_quartile_fQ4")
  ) %>%
  transmute(
    model = model_label,
    outcome = outcome_label,
    term = term_label,
    n,
    events,
    hr_ci,
    p = p_value_formatted
  )

log_lines <- c(
  "# NHANES V30 Mortality Models",
  "",
  paste0("Version: ", version_id),
  "",
  "## Status",
  "",
  "T2 survey-weighted mortality models completed. All-cause mortality is primary; CLRD mortality is exploratory because events are sparse.",
  "",
  "## Selected Results",
  "",
  markdown_table(summary_rows),
  "",
  "## Output Files",
  "",
  "- `results/tables/nhanes_mortality_models_v30_0.csv`",
  "- `results/tables/nhanes_mortality_model_cohort_counts_v30_0.csv`",
  "- `results/tables/nhanes_mortality_rates_by_quartile_v30_0.csv`",
  "- `results/tables/nhanes_mortality_spline_source_v30_0.csv`",
  "- `results/figures/nhanes_mortality_model_forest_v30_0.png` and `.pdf`",
  "- `results/figures/nhanes_mortality_spline_v30_0.png` and `.pdf`",
  "- `results/figures/nhanes_mortality_rates_by_quartile_v30_0.png` and `.pdf`",
  "- `manuscript/nhanes_mortality_limitations_v30_0.md`",
  "",
  "## Boundary",
  "",
  "Use bounded language: association, marker, mortality risk, supportive hard-outcome anchoring. Do not claim causality, validation of a new axis, or respiratory-specific mortality as definitive."
)
writeLines(log_lines, log_path)

message("NHANES mortality models complete.")
message("Primary M2 HR: ", main_row$hr_ci, ", p=", main_row$p_value_formatted)
message("Q4 vs Q1 HR: ", q4_row$hr_ci, ", p=", q4_row$p_value_formatted)
message("CLRD exploratory HR: ", clrd_row$hr_ci, ", p=", clrd_row$p_value_formatted)
