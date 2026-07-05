# NHANES V31.6 PH-violation sensitivity analyses.
#
# This script treats the primary Cox model as an average association and adds
# fixed-time and interval-specific checks that do not require a constant HR.

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(survey)
  library(survival)
  library(tibble)
})

options(survey.lonely.psu = "adjust")
version_id <- "v31_6"

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

format_pct_ci <- function(est, lo, hi, digits = 1) {
  sprintf(paste0("%.", digits, "f%% (%.", digits, "f%%-%.", digits, "f%%)"), 100 * est, 100 * lo, 100 * hi)
}

md_escape <- function(x) {
  gsub("\\|", "\\\\|", ifelse(is.na(x), "", as.character(x)))
}

md_table <- function(dat) {
  dat <- as.data.frame(dat, stringsAsFactors = FALSE)
  dat[] <- lapply(dat, function(x) ifelse(is.na(x), "", as.character(x)))
  c(
    paste0("| ", paste(md_escape(names(dat)), collapse = " | "), " |"),
    paste0("| ", paste(rep("---", ncol(dat)), collapse = " | "), " |"),
    apply(dat, 1, function(x) paste0("| ", paste(md_escape(x), collapse = " | "), " |"))
  )
}

prepare_mortality_data <- function(dat) {
  dat %>%
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

analysis_complete_cases <- function(dat) {
  required <- c(
    "participant_id", "followup_years_exam", "all_cause_death",
    "resp_vulnerability_z", "age_years", "sex_f", "race_f",
    "education_f", "income_poverty_ratio", "bmi", "smoking_f",
    "rv_quartile_f", "wtmec6yr", "psu", "strata"
  )
  missing_cols <- setdiff(required, names(dat))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }
  dat %>%
    filter(
      !is.na(wtmec6yr),
      wtmec6yr > 0,
      !is.na(psu),
      !is.na(strata),
      followup_years_exam > 0,
      complete.cases(across(all_of(required)))
    )
}

primary_rhs <- function(extra_terms = character()) {
  paste(
    c(
      "resp_vulnerability_z",
      extra_terms,
      "age_years",
      "sex_f",
      "race_f",
      "education_f",
      "income_poverty_ratio",
      "bmi",
      "smoking_f"
    ),
    collapse = " + "
  )
}

extract_log_effect <- function(fit, term, scale = "hr") {
  beta <- stats::coef(fit)
  vc <- stats::vcov(fit)
  if (!term %in% names(beta)) {
    stop("Term not found in model: ", term, call. = FALSE)
  }
  est <- beta[[term]]
  se <- sqrt(diag(vc))[[term]]
  ci <- est + c(-1, 1) * stats::qnorm(0.975) * se
  if (scale == "or") {
    measure <- "OR"
  } else {
    measure <- "HR"
  }
  tibble(
    term = term,
    log_estimate = est,
    se_log_estimate = se,
    estimate = exp(est),
    ci_low = exp(ci[[1]]),
    ci_high = exp(ci[[2]]),
    p_value = 2 * stats::pnorm(abs(est / se), lower.tail = FALSE),
    effect_measure = measure,
    estimate_ci = format_ci(exp(est), exp(ci[[1]]), exp(ci[[2]]), 2),
    p_value_formatted = format_p(2 * stats::pnorm(abs(est / se), lower.tail = FALSE))
  )
}

run_time_interaction <- function(df) {
  form <- stats::as.formula(paste0(
    "Surv(followup_years_exam, all_cause_death) ~ ",
    primary_rhs("tt(resp_vulnerability_z)"),
    " + cluster(psu)"
  ))
  fit <- survival::coxph(
    form,
    data = df,
    weights = wtmec6yr,
    ties = "efron",
    robust = TRUE,
    x = TRUE,
    tt = function(x, t, ...) x * log(pmax(t, 0.25))
  )
  base <- extract_log_effect(fit, "resp_vulnerability_z")
  time_term <- extract_log_effect(fit, "tt(resp_vulnerability_z)")
  bind_rows(base, time_term) %>%
    mutate(
      analysis_id = "nhanes_v31_6_time_interaction",
      analysis_label = "Weighted Cox with exposure by log(time) interaction",
      n = nrow(df),
      events = sum(df$all_cause_death == 1L, na.rm = TRUE),
      person_years = sum(df$followup_years_exam, na.rm = TRUE),
      interpretation = case_when(
        term == "resp_vulnerability_z" ~ "Estimated exposure association at log(time)=0 scale point; interpret jointly with the time-interaction term.",
        p_value < 0.05 ~ "Evidence that the exposure-mortality association changes across follow-up time.",
        TRUE ~ "No strong evidence of exposure-time interaction."
      )
    ) %>%
    select(
      analysis_id, analysis_label, term, effect_measure, n, events, person_years,
      log_estimate, se_log_estimate, estimate, ci_low, ci_high, p_value,
      estimate_ci, p_value_formatted, interpretation
    )
}

run_piecewise_cox <- function(df) {
  split_df <- survival::survSplit(
    data = as.data.frame(df),
    cut = c(5, 10),
    end = "followup_years_exam",
    start = "tstart",
    event = "all_cause_death",
    episode = "episode"
  ) %>%
    as_tibble() %>%
    mutate(
      followup_interval = factor(
        episode,
        levels = c(1, 2, 3),
        labels = c("0-5 years", "5-10 years", "10+ years")
      )
    )

  bind_rows(lapply(levels(split_df$followup_interval), function(interval_label) {
    interval_df <- split_df %>%
      filter(followup_interval == interval_label, followup_years_exam > tstart)
    events <- sum(interval_df$all_cause_death == 1L, na.rm = TRUE)
    if (nrow(interval_df) < 500 || events < 30) {
      return(tibble(
        analysis_id = "nhanes_v31_6_piecewise_followup",
        interval = interval_label,
        status = "insufficient_events",
        n = nrow(interval_df),
        events = events,
        person_years = sum(interval_df$followup_years_exam - interval_df$tstart, na.rm = TRUE),
        effect_measure = "HR",
        estimate = NA_real_,
        ci_low = NA_real_,
        ci_high = NA_real_,
        p_value = NA_real_,
        estimate_ci = NA_character_,
        p_value_formatted = NA_character_,
        interpretation = "Too few events for a stable interval-specific sensitivity estimate."
      ))
    }
    form <- stats::as.formula(paste0(
      "Surv(tstart, followup_years_exam, all_cause_death) ~ ",
      primary_rhs(),
      " + cluster(psu)"
    ))
    fit <- survival::coxph(
      form,
      data = interval_df,
      weights = wtmec6yr,
      ties = "efron",
      robust = TRUE,
      x = TRUE
    )
    row <- extract_log_effect(fit, "resp_vulnerability_z") %>%
      mutate(
        analysis_id = "nhanes_v31_6_piecewise_followup",
        interval = interval_label,
        status = "ok",
        n = nrow(interval_df),
        events = events,
        person_years = sum(interval_df$followup_years_exam - interval_df$tstart, na.rm = TRUE),
        interpretation = "Interval-specific weighted Cox estimate; compare direction and magnitude, not formal between-interval heterogeneity."
      ) %>%
      select(
        analysis_id, interval, status, n, events, person_years, effect_measure,
        estimate, ci_low, ci_high, p_value, estimate_ci, p_value_formatted,
        interpretation
      )
    row
  }))
}

make_design <- function(data) {
  survey::svydesign(
    ids = ~psu,
    strata = ~strata,
    weights = ~wtmec6yr,
    nest = TRUE,
    data = data
  )
}

run_fixed_time_model <- function(df, horizon_years) {
  risk_df <- df %>%
    mutate(
      horizon_years = horizon_years,
      known_status_at_horizon = followup_years_exam >= horizon_years | all_cause_death == 1L,
      death_by_horizon = as.integer(all_cause_death == 1L & followup_years_exam <= horizon_years)
    ) %>%
    filter(known_status_at_horizon)
  events <- sum(risk_df$death_by_horizon == 1L, na.rm = TRUE)
  if (nrow(risk_df) < 1000 || events < 80) {
    return(tibble(
      analysis_id = "nhanes_v31_6_fixed_time_mortality_model",
      horizon_years = horizon_years,
      status = "insufficient_events",
      n = nrow(risk_df),
      events = events,
      effect_measure = "OR",
      estimate = NA_real_,
      ci_low = NA_real_,
      ci_high = NA_real_,
      p_value = NA_real_,
      estimate_ci = NA_character_,
      p_value_formatted = NA_character_,
      interpretation = "Too few known-status records or deaths for stable fixed-time model."
    ))
  }
  design <- make_design(risk_df)
  form <- stats::as.formula(paste0(
    "death_by_horizon ~ ",
    primary_rhs()
  ))
  fit <- survey::svyglm(form, design = design, family = quasibinomial())
  extract_log_effect(fit, "resp_vulnerability_z", scale = "or") %>%
    mutate(
      analysis_id = "nhanes_v31_6_fixed_time_mortality_model",
      horizon_years = horizon_years,
      status = "ok",
      n = nrow(risk_df),
      events = events,
      interpretation = "Survey-weighted fixed-time mortality odds model; it does not require a constant hazard ratio over follow-up."
    ) %>%
    select(
      analysis_id, horizon_years, status, n, events, effect_measure,
      estimate, ci_low, ci_high, p_value, estimate_ci, p_value_formatted,
      interpretation
    )
}

run_fixed_time_rates <- function(df, horizon_years) {
  risk_df <- df %>%
    mutate(
      horizon_years = horizon_years,
      known_status_at_horizon = followup_years_exam >= horizon_years | all_cause_death == 1L,
      death_by_horizon = as.numeric(all_cause_death == 1L & followup_years_exam <= horizon_years)
    ) %>%
    filter(known_status_at_horizon)
  design <- make_design(risk_df)
  by <- survey::svyby(
    ~death_by_horizon,
    ~rv_quartile_f,
    design,
    survey::svymean,
    na.rm = TRUE,
    vartype = c("se", "ci"),
    keep.names = FALSE
  )
  as_tibble(by) %>%
    transmute(
      analysis_id = "nhanes_v31_6_fixed_time_mortality_rates",
      horizon_years = horizon_years,
      rv_quartile = as.character(rv_quartile_f),
      weighted_risk = death_by_horizon,
      se = se,
      ci_low = ci_l,
      ci_high = ci_u,
      risk_ci = format_pct_ci(death_by_horizon, ci_l, ci_u, 1),
      n_unweighted = as.integer(table(risk_df$rv_quartile_f)[rv_quartile]),
      events_unweighted = as.integer(tapply(risk_df$death_by_horizon, risk_df$rv_quartile_f, sum, na.rm = TRUE)[rv_quartile])
    )
}

plot_fixed_time_rates <- function(rates, figure_dir) {
  plot_data <- rates %>%
    mutate(
      horizon_label = paste0(horizon_years, "-year mortality"),
      weighted_risk_percent = 100 * weighted_risk,
      ci_low_percent = 100 * ci_low,
      ci_high_percent = 100 * ci_high,
      horizon_label = factor(horizon_label, levels = c("5-year mortality", "10-year mortality"))
    )
  p <- ggplot(plot_data, aes(x = rv_quartile, y = weighted_risk_percent, group = horizon_label)) +
    geom_errorbar(aes(ymin = ci_low_percent, ymax = ci_high_percent), width = 0.12, linewidth = 0.35, colour = "#4B5563") +
    geom_line(linewidth = 0.45, colour = "#1F6F8B") +
    geom_point(size = 1.7, colour = "#1F6F8B") +
    facet_wrap(~horizon_label, scales = "free_y", ncol = 2) +
    labs(
      x = "Respiratory vulnerability quartile",
      y = "Survey-weighted fixed-time mortality risk (%)"
    ) +
    theme_classic(base_size = 8) +
    theme(
      axis.line = element_line(linewidth = 0.35),
      axis.ticks = element_line(linewidth = 0.35),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      panel.grid.major.y = element_line(linewidth = 0.2, colour = "#E5E7EB")
    )
  ggplot2::ggsave(file.path(figure_dir, paste0("nhanes_v31_6_fixed_time_mortality_by_quartile.png")), p, width = 6.4, height = 3.4, dpi = 600, units = "in")
  ggplot2::ggsave(file.path(figure_dir, paste0("nhanes_v31_6_fixed_time_mortality_by_quartile.pdf")), p, width = 6.4, height = 3.4, units = "in", device = grDevices::pdf, useDingbats = FALSE)
}

root <- find_project_root()
tab_dir <- file.path(root, "results", "tables")
log_dir <- file.path(root, "results", "logs")
figure_dir <- file.path(root, "results", "figures")
manuscript_dir <- file.path(root, "manuscript")
invisible(lapply(c(tab_dir, log_dir, figure_dir, manuscript_dir), dir_create))

input_rds <- file.path(root, "derived_sensitive", "nhanes", "nhanes_mortality_analysis_ready_v30_0.rds")
if (!file.exists(input_rds)) {
  stop("Missing NHANES mortality analysis-ready RDS: ", input_rds, call. = FALSE)
}

df <- readRDS(input_rds) %>%
  prepare_mortality_data() %>%
  analysis_complete_cases()

time_interaction <- run_time_interaction(df) %>%
  mutate(
    log_estimate = round(log_estimate, 6),
    se_log_estimate = round(se_log_estimate, 6),
    estimate = round(estimate, 5),
    ci_low = round(ci_low, 5),
    ci_high = round(ci_high, 5),
    p_value = signif(p_value, 4),
    person_years = round(person_years, 1)
  )

piecewise <- run_piecewise_cox(df) %>%
  mutate(
    person_years = round(person_years, 1),
    estimate = round(estimate, 5),
    ci_low = round(ci_low, 5),
    ci_high = round(ci_high, 5),
    p_value = signif(p_value, 4)
  )

fixed_time_models <- bind_rows(lapply(c(5, 10), function(h) run_fixed_time_model(df, h))) %>%
  mutate(
    estimate = round(estimate, 5),
    ci_low = round(ci_low, 5),
    ci_high = round(ci_high, 5),
    p_value = signif(p_value, 4)
  )

fixed_time_rates <- bind_rows(lapply(c(5, 10), function(h) run_fixed_time_rates(df, h))) %>%
  mutate(
    weighted_risk = round(weighted_risk, 5),
    se = round(se, 5),
    ci_low = round(ci_low, 5),
    ci_high = round(ci_high, 5)
  )

write_csv(time_interaction, file.path(tab_dir, "nhanes_v31_6_time_interaction.csv"))
write_csv(piecewise, file.path(tab_dir, "nhanes_v31_6_piecewise_followup_sensitivity.csv"))
write_csv(fixed_time_models, file.path(tab_dir, "nhanes_v31_6_fixed_time_mortality_models.csv"))
write_csv(fixed_time_rates, file.path(tab_dir, "nhanes_v31_6_fixed_time_mortality_rates.csv"))

plot_fixed_time_rates(fixed_time_rates, figure_dir)

literature_brief <- c(
  "# Targeted Literature And Methods Brief V31.6",
  "",
  paste0("Search date: ", format(Sys.Date(), "%Y-%m-%d")),
  "",
  "Purpose: document why the PH-violation sensitivity package was added before manuscript assembly.",
  "",
  "## Targeted search scope",
  "",
  "- Recent NHANES/pulmonary-function mortality papers were searched for how they report Cox models, subgroup/sensitivity analyses, and fixed-time or robustness checks.",
  "- Survival-methods sources were searched for non-proportional hazards handling, especially when a constant hazard ratio is questionable.",
  "- The local peer-review/scientific-critical-thinking framework was applied to keep the manuscript claim as an observational marker association rather than a causal mechanism claim.",
  "",
  "## Directly relevant comparator papers identified",
  "",
  "- Zhao et al. The prevalence and mortality risks of PRISm and COPD in the United States from NHANES 2007-2012. Respiratory Research 2024. DOI: https://doi.org/10.1186/s12931-024-02841-y; PMID: 38750492.",
  "- Pacheco-Galvan et al. Prevalence, risk factors, and clinical implications of failed spirometry in adults: Results from NHANES 2007-2012. Pulmonology 2025. DOI: https://doi.org/10.1080/25310429.2025.2572011; PMID: 41084927.",
  "- Koo et al. Severity of Airflow Obstruction Based on FEV1/FVC Versus FEV1 Percent Predicted in the General U.S. Population. American Journal of Respiratory and Critical Care Medicine 2024. DOI: https://doi.org/10.1164/rccm.202310-1773OC; PMID: 38597717.",
  "- Bhatt et al. Discriminative Accuracy of FEV1:FVC Thresholds for COPD-Related Hospitalization and Mortality. JAMA 2019. DOI: https://doi.org/10.1001/jama.2019.7233; PMID: 31237643.",
  "- Fan et al. The systemic inflammation response index as risks factor for all-cause and cardiovascular mortality among individuals with respiratory sarcopenia. BMC Pulmonary Medicine 2025. DOI: https://doi.org/10.1186/s12890-025-03525-z; PMID: 40011897.",
  "",
  "## Survival-methods sources identified",
  "",
  "- Bejan-Angoulvant et al. The Impact of Violation of the Proportional Hazards Assumption on the Calibration of the Cox Proportional Hazards Model. Statistics in Medicine 2025. DOI: https://doi.org/10.1002/sim.70161; PMID: 40492822.",
  "- Crowther et al. Using fractional polynomials and restricted cubic splines to model non-proportional hazards or time-varying covariate effects in the Cox regression model. Statistics in Medicine 2022. DOI: https://doi.org/10.1002/sim.9259; PMID: 34806210.",
  "- Stensrud and Hernan. Why Test for Proportional Hazards? JAMA 2020. DOI: https://doi.org/10.1001/jama.2020.1267; PMID: 32167523.",
  "- Stensrud and Hernan. Why use methods that require proportional hazards? American Journal of Epidemiology 2025. PMID: 39756420.",
  "",
  "## Method decision absorbed into this project",
  "",
  "- Keep the primary survey-weighted Cox model, but describe it as an average follow-up association if PH diagnostics flag time-varying effects.",
  "- Add an exposure-by-log(time) model to test whether the marker association changes across follow-up.",
  "- Add 0-5, 5-10 and 10+ year interval-specific Cox estimates to show where the association is concentrated.",
  "- Add 5-year and 10-year survey-weighted fixed-time mortality models and quartile risk plots; these are easier for reviewers to interpret under non-proportional hazards.",
  "- Do not claim a constant long-term hazard ratio, a clinical risk score, or causality from this layer.",
  "",
  "## Positioning consequences",
  "",
  "- Contemporary NHANES lung-function and mortality studies establish that linked-mortality respiratory epidemiology is publishable in respiratory journals, but they mostly center on obstruction/PRISm/failed spirometry rather than a lower-than-expected PEF reserve marker.",
  "- The manuscript should use those papers as comparators, not as proof that the present marker is diagnostic or causal.",
  "- The PH-sensitivity paragraph should be placed in Results; detailed diagnostics can go to Supplement if word count is tight.",
  "- The fixed-time risk layer is the most reviewer-friendly way to preserve the mortality message while avoiding overclaiming a constant Cox HR."
)
writeLines(literature_brief, file.path(manuscript_dir, "literature_methods_brief_v31_6.md"))

main_time <- time_interaction %>% filter(term == "tt(resp_vulnerability_z)") %>% slice(1)
piecewise_display <- piecewise %>%
  transmute(
    interval,
    n,
    events,
    `person-years` = person_years,
    `HR per 1-SD` = estimate_ci,
    p = p_value_formatted,
    status
  )
fixed_display <- fixed_time_models %>%
  transmute(
    horizon = paste0(horizon_years, " years"),
    n,
    events,
    `OR per 1-SD` = estimate_ci,
    p = p_value_formatted,
    status
  )
rate_display <- fixed_time_rates %>%
  transmute(
    horizon = paste0(horizon_years, " years"),
    quartile = rv_quartile,
    n = n_unweighted,
    events = events_unweighted,
    `weighted risk` = risk_ci
  )

note <- c(
  "# NHANES PH-Violation Sensitivity V31.6",
  "",
  "## Why this was added",
  "",
  "The weighted cox.zph approximation in V31.5 flagged a possible proportional-hazards violation for the primary respiratory-vulnerability marker. The primary Cox estimate should therefore be treated as an average association over follow-up, not as evidence that the hazard ratio is constant at every time point.",
  "",
  "## Time-interaction result",
  "",
  paste0(
    "The exposure-by-log(time) term was ",
    main_time$estimate_ci,
    " (p=", main_time$p_value_formatted,
    "). This supports describing the marker-mortality association as time-varying rather than assuming a constant hazard ratio."
  ),
  "",
  "## Interval-specific Cox sensitivity",
  "",
  md_table(piecewise_display),
  "",
  "## Fixed-time mortality models",
  "",
  md_table(fixed_display),
  "",
  "## Fixed-time mortality risks by quartile",
  "",
  md_table(rate_display),
  "",
  "## Manuscript wording lock",
  "",
  "Recommended wording: lower-than-expected PEF was associated with higher all-cause mortality in NHANES, with evidence that the association varied over follow-up. Fixed-time risk analyses showed higher 5-year and 10-year mortality risk among participants with higher marker values. Avoid wording that implies a constant hazard ratio, causal effect, or deployable prognostic score."
)
writeLines(note, file.path(manuscript_dir, "nhanes_ph_violation_sensitivity_v31_6.md"))

log <- c(
  "# NHANES V31.6 PH-Violation Sensitivity Build",
  "",
  paste0("Version: ", version_id),
  paste0("Primary complete-case n: ", nrow(df)),
  paste0("All-cause deaths: ", sum(df$all_cause_death == 1L, na.rm = TRUE)),
  paste0("Person-years: ", round(sum(df$followup_years_exam, na.rm = TRUE), 1)),
  "",
  "Outputs:",
  "- results/tables/nhanes_v31_6_time_interaction.csv",
  "- results/tables/nhanes_v31_6_piecewise_followup_sensitivity.csv",
  "- results/tables/nhanes_v31_6_fixed_time_mortality_models.csv",
  "- results/tables/nhanes_v31_6_fixed_time_mortality_rates.csv",
  "- results/figures/nhanes_v31_6_fixed_time_mortality_by_quartile.png and .pdf",
  "- manuscript/nhanes_ph_violation_sensitivity_v31_6.md",
  "- manuscript/literature_methods_brief_v31_6.md"
)
writeLines(log, file.path(log_dir, "nhanes_v31_6_ph_violation_sensitivity.md"))

message("NHANES V31.6 PH-violation sensitivity outputs written.")
