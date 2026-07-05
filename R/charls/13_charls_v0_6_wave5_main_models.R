# CHARLS v0.6 Wave 5-extended main models.
# Adds raw CHARLS 2020 disease/death candidates to the existing Harmonized
# CHARLS D Wave 1-4 core. Writes aggregate outputs only.

suppressPackageStartupMessages({
  library(dplyr)
  library(haven)
  library(readr)
  library(survey)
  library(survival)
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

to_num <- function(x) suppressWarnings(as.numeric(haven::zap_labels(haven::zap_missing(x))))
to_chr <- function(x) as.character(haven::zap_labels(haven::zap_missing(x)))

read_zip_member <- function(zip_path, member, cols) {
  td <- tempfile("charls_v06_w5_")
  dir.create(td, recursive = TRUE)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)
  utils::unzip(zip_path, files = member, exdir = td, overwrite = TRUE)
  haven::read_dta(file.path(td, member), col_select = tidyselect::any_of(cols))
}

label_of <- function(dat, var) {
  lab <- attr(dat[[var]], "label")
  if (is.null(lab)) NA_character_ else as.character(lab)
}

value_label_text <- function(dat, var) {
  labs <- attr(dat[[var]], "labels")
  if (is.null(labs)) return(NA_character_)
  paste(paste0(unname(labs), "=", names(labs)), collapse = "; ")
}

first_nonmissing <- function(x) {
  idx <- which(!is.na(x))
  if (length(idx) == 0) return(x[NA_integer_][1])
  x[idx[[1]]]
}

first_wave_match <- function(values, ...) {
  mat <- cbind(...)
  hit <- mat %in% values
  dim(hit) <- dim(mat)
  out <- rep(NA_integer_, nrow(mat))
  for (j in seq_len(ncol(mat))) {
    idx <- is.na(out) & hit[, j]
    out[idx] <- j + 1L
  }
  out
}

last_observed_wave <- function(...) {
  mat <- cbind(...)
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

load_wave5 <- function(zip_path) {
  health <- read_zip_member(
    zip_path,
    "Health_Status_and_Functioning.dta",
    c("ID", "householdID", "communityID", "da002_5_", "da003_5_", "da004_5_")
  )
  sample_info <- read_zip_member(
    zip_path,
    "Sample_Infor.dta",
    c("ID", "householdID", "communityID", "died")
  )
  exit_mod <- read_zip_member(
    zip_path,
    "Exit_Module.dta",
    c("ID", "householdID", "communityID", "exb001_1", "exb001_2", "exb001_3", "xezdisease_2_")
  )
  weights <- read_zip_member(
    zip_path,
    "Weights.dta",
    c("ID", "householdID", "communityID", "INDV_weight", "INDV_weight_ad2")
  )

  lock <- bind_rows(
    tibble(
      concept = "wave5_chronic_lung_disease_status",
      module = "Health_Status_and_Functioning.dta",
      variable = "da003_5_",
      label = label_of(health, "da003_5_"),
      value_labels = value_label_text(health, "da003_5_"),
      rule = "Wave 5 chronic lung disease status is coded as 1 if da003_5_ == 1, 0 if da003_5_ == 2, and missing otherwise.",
      lock_status = "locked_v0_6_raw_2020_codebook"
    ),
    tibble(
      concept = "wave5_chronic_lung_disease_change",
      module = "Health_Status_and_Functioning.dta",
      variable = "da002_5_",
      label = label_of(health, "da002_5_"),
      value_labels = value_label_text(health, "da002_5_"),
      rule = "Used for context only; primary Wave 5 status uses da003_5_.",
      lock_status = "context_locked_v0_6"
    ),
    tibble(
      concept = "wave5_death_status",
      module = "Sample_Infor.dta",
      variable = "died",
      label = label_of(sample_info, "died"),
      value_labels = value_label_text(sample_info, "died"),
      rule = "Wave 5 death is coded as 1 if died == 1, 0 if died == 0, and missing otherwise.",
      lock_status = "locked_v0_6_raw_2020_codebook"
    ),
    tibble(
      concept = "wave5_death_date",
      module = "Exit_Module.dta",
      variable = "exb001_1/exb001_2/exb001_3",
      label = paste(label_of(exit_mod, "exb001_1"), label_of(exit_mod, "exb001_2"), label_of(exit_mod, "exb001_3"), sep = "; "),
      value_labels = NA_character_,
      rule = "Available for future exact-date mortality sensitivity, but main v0.6 keeps wave-interval timing for consistency.",
      lock_status = "candidate_future_exact_timing"
    ),
    tibble(
      concept = "wave5_exit_chronic_lung_disease",
      module = "Exit_Module.dta",
      variable = "xezdisease_2_",
      label = label_of(exit_mod, "xezdisease_2_"),
      value_labels = value_label_text(exit_mod, "xezdisease_2_"),
      rule = "Decedent-proxy context only; not used as the primary incident CLD endpoint.",
      lock_status = "context_only"
    ),
    tibble(
      concept = "wave5_individual_weight",
      module = "Weights.dta",
      variable = "INDV_weight/INDV_weight_ad2",
      label = paste(label_of(weights, "INDV_weight"), label_of(weights, "INDV_weight_ad2"), sep = "; "),
      value_labels = NA_character_,
      rule = "Documented for Wave 5; baseline biomarker weight remains the primary longitudinal analysis weight for comparability with prior CHARLS models.",
      lock_status = "documented_not_primary_weight"
    )
  )

  wave5 <- sample_info %>%
    transmute(
      participant_id = to_chr(.data$ID),
      household_id_w5 = to_chr(.data$householdID),
      community_id_w5 = to_chr(.data$communityID),
      died_w5_raw = to_num(.data$died),
      died_w5 = case_when(
        .data$died_w5_raw == 1 ~ 1L,
        .data$died_w5_raw == 0 ~ 0L,
        TRUE ~ NA_integer_
      )
    ) %>%
    full_join(
      health %>%
        transmute(
          participant_id = to_chr(.data$ID),
          lung_w5_raw = to_num(.data$da003_5_),
          lung_change_w5_raw = to_num(.data$da002_5_),
          lung_self_known_w5_raw = to_num(.data$da004_5_),
          lung_w5 = case_when(
            .data$lung_w5_raw == 1 ~ 1L,
            .data$lung_w5_raw == 2 ~ 0L,
            TRUE ~ NA_integer_
          )
        ),
      by = "participant_id"
    ) %>%
    full_join(
      exit_mod %>%
        transmute(
          participant_id = to_chr(.data$ID),
          death_year_w5 = to_num(.data$exb001_1),
          death_month_w5 = to_num(.data$exb001_2),
          death_day_w5 = to_num(.data$exb001_3),
          exit_lung_w5_raw = to_num(.data$xezdisease_2_),
          exit_lung_w5 = case_when(.data$exit_lung_w5_raw == 1 ~ 1L, TRUE ~ NA_integer_)
        ),
      by = "participant_id"
    ) %>%
    full_join(
      weights %>%
        transmute(
          participant_id = to_chr(.data$ID),
          indv_weight_w5 = to_num(.data$INDV_weight),
          indv_weight_ad2_w5 = to_num(.data$INDV_weight_ad2)
        ),
      by = "participant_id"
    ) %>%
    group_by(.data$participant_id) %>%
    summarise(across(everything(), first_nonmissing), .groups = "drop")

  list(wave5 = wave5, lock = lock, health = health, sample_info = sample_info, exit_mod = exit_mod, weights = weights)
}

prepare_base_w5 <- function(core_w5) {
  core_w5 %>%
    mutate(
      resp_vulnerability_z = -pef_resid_z_w1,
      pef_per_100_l_min = pef_best_w1_valid_provisional / 100,
      age_decade = age_w1 / 10,
      sex = factor(sex_label),
      smoke_ever = factor(smoke_ever_w1, levels = c(0, 1), labels = c("no", "yes")),
      psu_community_id = community_id,
      county_pseudo_id = substr(community_id, 1, 6),
      analysis_weight = weight_biomarker_w1,
      death_wave_w1_w5 = first_wave_match(c(5, 6), iwstat_w2, iwstat_w3, iwstat_w4, ifelse(died_w5 == 1, 5, NA_real_)),
      last_iw_wave_w1_w5 = last_observed_wave(iwstat_w2, iwstat_w3, iwstat_w4, ifelse(!is.na(died_w5), 0, NA_real_)),
      last_lung_wave_w1_w5 = last_observed_wave(lung_w2, lung_w3, lung_w4, lung_w5)
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

prepare_lung_time_w5 <- function(base) {
  event_wave <- first_wave_match(1, base$lung_w2, base$lung_w3, base$lung_w4, base$lung_w5)
  last_disease_wave <- last_observed_wave(base$lung_w2, base$lung_w3, base$lung_w4, base$lung_w5)
  death_wave <- base$death_wave_w1_w5
  disease_before_death <- !is.na(event_wave) & (is.na(death_wave) | event_wave <= death_wave)
  death_before_disease <- !is.na(death_wave) & (is.na(event_wave) | death_wave < event_wave)
  censor_wave <- ifelse(death_before_disease, death_wave, last_disease_wave)
  observed_wave <- ifelse(disease_before_death, event_wave, censor_wave)
  composite_wave <- min_wave(event_wave, death_wave)
  composite_wave <- ifelse(is.na(composite_wave), last_disease_wave, composite_wave)

  base %>%
    mutate(
      outcome = "incident_chronic_lung_disease_w2_w5",
      disease_event_wave = event_wave,
      disease_last_observed_wave = last_disease_wave,
      disease_event = as.integer(disease_before_death),
      competing_death = as.integer(death_before_disease),
      time = observed_wave - 1L,
      composite_event = as.integer(disease_before_death | death_before_disease),
      composite_time = composite_wave - 1L
    ) %>%
    filter(lung_w1 == 0, !is.na(time), time > 0)
}

prepare_death_time_w5 <- function(base) {
  base %>%
    mutate(
      outcome = "death_w2_w5",
      death_event = as.integer(!is.na(death_wave_w1_w5)),
      time = ifelse(!is.na(death_wave_w1_w5), death_wave_w1_w5, last_iw_wave_w1_w5) - 1L
    ) %>%
    filter(!is.na(time), time > 0)
}

make_period_data_w5 <- function(dat) {
  periods <- lapply(1:4, function(interval) {
    target_wave <- interval + 1L
    dat %>%
      mutate(
        interval = interval,
        target_wave = target_wave,
        period_event = as.integer(!is.na(disease_event_wave) & disease_event_wave == target_wave),
        at_risk = is.na(disease_event_wave) | disease_event_wave >= target_wave,
        before_death = is.na(death_wave_w1_w5) | death_wave_w1_w5 > target_wave | (!is.na(disease_event_wave) & disease_event_wave == target_wave)
      ) %>%
      filter(at_risk, before_death)
  })
  bind_rows(periods) %>% mutate(outcome = "incident_chronic_lung_disease_w2_w5")
}

build_design <- function(dat, design_type) {
  if (design_type == "public_psu") {
    survey::svydesign(ids = ~psu_community_id, weights = ~analysis_weight, data = dat, nest = TRUE)
  } else if (design_type == "public_county_community") {
    survey::svydesign(ids = ~county_pseudo_id + psu_community_id, weights = ~analysis_weight, data = dat, nest = TRUE)
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
  formula <- stats::as.formula(paste0("survival::Surv(", time_var, ", ", event_var, ") ~ ", rhs))
  needed <- unique(c(all.vars(formula), "psu_community_id", "county_pseudo_id", "analysis_weight"))
  analysis <- dat %>%
    select(any_of(needed)) %>%
    filter(if_all(everything(), ~ !is.na(.x)))
  n_obs <- nrow(analysis)
  events <- sum(analysis[[event_var]] == 1, na.rm = TRUE)
  if (n_obs < 50 || events < 10) {
    return(list(table = tibble(), warnings = tibble()))
  }
  fit_result <- capture_warnings(survey::svycoxph(formula, design = build_design(analysis, design_type)))
  list(
    table = effect_table(fit_result$value, outcome, model, design_type, "survey_HR", n_obs, events),
    warnings = tibble(outcome = outcome, model = model, design_type = design_type, warning = fit_result$warnings)
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
  fit_result <- capture_warnings(survey::svyglm(formula, design = build_design(analysis, design_type), family = stats::quasibinomial()))
  list(
    table = effect_table(fit_result$value, outcome, model, design_type, "survey_OR", n_obs, events),
    warnings = tibble(outcome = outcome, model = model, design_type = design_type, warning = fit_result$warnings)
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

wave_event_counts <- function(lung, death) {
  bind_rows(
    lung %>%
      count(disease_event_wave, name = "n") %>%
      transmute(endpoint = "incident_chronic_lung_disease", wave = disease_event_wave, n),
    death %>%
      filter(death_event == 1) %>%
      count(time, name = "n") %>%
      transmute(endpoint = "death", wave = time + 1L, n)
  ) %>%
    arrange(endpoint, wave)
}

format_key <- function(key_terms, outcome_id, model_id, design_id = "public_psu") {
  row <- key_terms %>%
    filter(
      .data$outcome == .env$outcome_id,
      .data$model == .env$model_id,
      .data$design_type == .env$design_id,
      .data$term == "resp_vulnerability_z"
    ) %>%
    slice(1)
  if (nrow(row) == 0) return("not estimated")
  paste0(row$effect_ci, ", p=", row$p_formatted)
}

main <- function() {
  root <- find_project_root()
  input_rds <- file.path(root, "derived_sensitive", "charls", "charls_core_harmonized_provisional.rds")
  zip_path <- "${CHARLS_RAW_ROOT}/2020年全国追踪调查/数据下载/CHARLS2020r.zip"
  if (!file.exists(input_rds)) stop("Missing CHARLS core RDS.", call. = FALSE)
  if (!file.exists(zip_path)) stop("Missing CHARLS 2020 zip: ", zip_path, call. = FALSE)

  table_dir <- file.path(root, "results", "tables")
  log_dir <- file.path(root, "results", "logs")
  metadata_dir <- file.path(root, "metadata")
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(metadata_dir, recursive = TRUE, showWarnings = FALSE)

  core <- readRDS(input_rds)
  w5_obj <- load_wave5(zip_path)
  wave5 <- w5_obj$wave5
  lock <- w5_obj$lock

  core_w5 <- core %>%
    left_join(wave5, by = "participant_id") %>%
    mutate(
      wave5_join_status = case_when(
        !is.na(lung_w5) | !is.na(died_w5) ~ "linked_wave5_status",
        TRUE ~ "no_wave5_status"
      )
    )

  base <- prepare_base_w5(core_w5)
  lung <- prepare_lung_time_w5(base)
  death <- prepare_death_time_w5(base)
  lung_period <- make_period_data_w5(lung)

  rhs_base <- "resp_vulnerability_z + age_decade + sex + smoke_ever + bmi_w1"
  rhs_frailty <- "resp_vulnerability_z + frailty_proxy_count_w1 + age_decade + sex + smoke_ever + bmi_w1"
  design_types <- c("public_psu", "public_county_community")

  model_results <- list()
  warning_results <- list()
  add_result <- function(res) {
    model_results[[length(model_results) + 1]] <<- res$table
    warning_results[[length(warning_results) + 1]] <<- res$warnings
  }

  for (design_type in design_types) {
    add_result(fit_svycox(lung, "incident_chronic_lung_disease_w2_w5", "disease_event", "time", rhs_base, "svycox_cause_specific", design_type))
    add_result(fit_svycox(lung, "incident_chronic_lung_disease_w2_w5", "disease_event", "time", rhs_frailty, "svycox_plus_frailty", design_type))
    add_result(fit_svycox(lung, "incident_chronic_lung_disease_w2_w5", "composite_event", "composite_time", rhs_base, "svycox_composite_disease_or_death", design_type))
    add_result(fit_svyglm_discrete(lung_period, "incident_chronic_lung_disease_w2_w5", rhs_base, "svyglm_discrete_time", design_type))
    add_result(fit_svycox(death, "death_w2_w5", "death_event", "time", rhs_base, "svycox_death", design_type))
    add_result(fit_svycox(death, "death_w2_w5", "death_event", "time", rhs_frailty, "svycox_death_plus_frailty", design_type))
  }

  model_table <- bind_rows(model_results)
  warnings_table <- bind_rows(warning_results)
  if (nrow(warnings_table) == 0) {
    warnings_table <- tibble(outcome = character(), model = character(), design_type = character(), warning = character())
  }

  key_terms <- model_table %>%
    filter(.data$term %in% c("resp_vulnerability_z", "frailty_proxy_count_w1")) %>%
    mutate(
      effect_ci = sprintf("%.2f (%.2f, %.2f)", effect, effect_low, effect_high),
      event_rate_display = sprintf("%.1f%%", 100 * event_rate)
    ) %>%
    select(outcome, model, design_type, estimand, n_obs, events, event_rate_display, term, effect_ci, p_formatted)

  cohort_counts <- bind_rows(
    cohort_count_table(lung, "incident_chronic_lung_disease_w2_w5"),
    cohort_count_table(death, "death_w2_w5")
  )
  event_counts <- wave_event_counts(lung, death)
  join_counts <- tibble(
    metric = c(
      "core_rows",
      "wave5_rows",
      "core_with_wave5_lung_status",
      "core_with_wave5_death_status",
      "wave5_lung_yes_linked_core",
      "wave5_deaths_linked_core",
      "baseline_lung_at_risk_with_wave5_or_prior_followup"
    ),
    value = c(
      nrow(core),
      nrow(wave5),
      sum(!is.na(core_w5$lung_w5)),
      sum(!is.na(core_w5$died_w5)),
      sum(core_w5$lung_w5 == 1, na.rm = TRUE),
      sum(core_w5$died_w5 == 1, na.rm = TRUE),
      nrow(lung)
    )
  )

  readr::write_csv(lock, file.path(metadata_dir, "charls_wave5_codebook_lock_v0_6.csv"))
  readr::write_csv(model_table, file.path(table_dir, "charls_v0_6_wave5_model_table.csv"))
  readr::write_csv(key_terms, file.path(table_dir, "charls_v0_6_wave5_key_terms.csv"))
  readr::write_csv(cohort_counts, file.path(table_dir, "charls_v0_6_wave5_cohort_counts.csv"))
  readr::write_csv(event_counts, file.path(table_dir, "charls_v0_6_wave5_event_counts.csv"))
  readr::write_csv(join_counts, file.path(table_dir, "charls_v0_6_wave5_join_counts.csv"))
  readr::write_csv(warnings_table, file.path(table_dir, "charls_v0_6_wave5_model_warnings.csv"))

  primary_lung <- format_key(key_terms, "incident_chronic_lung_disease_w2_w5", "svycox_cause_specific")
  frailty_lung <- format_key(key_terms, "incident_chronic_lung_disease_w2_w5", "svycox_plus_frailty")
  discrete_lung <- format_key(key_terms, "incident_chronic_lung_disease_w2_w5", "svyglm_discrete_time")
  death_key <- format_key(key_terms, "death_w2_w5", "svycox_death")

  findings <- c(
    "# CHARLS V0.6 Wave 5-Extended Findings",
    "",
    paste0("- Run date: ", Sys.Date()),
    "- Scope: Wave 1-5 CHARLS longitudinal extension using Harmonized CHARLS D Wave 1-4 plus raw 2020 module candidates.",
    "- Primary Wave 5 CLD variable: `Health_Status_and_Functioning.dta::da003_5_`, coded 1=Yes, 2=No.",
    "- Primary Wave 5 death variable: `Sample_Infor.dta::died`, coded 0=Alive, 1=Died.",
    "- Baseline biomarker weights and community-level PSU clustering remain the primary survey design for comparability with prior CHARLS models.",
    "",
    "## Primary Public-PSU Results",
    "",
    paste0("- Incident self-reported physician-diagnosed chronic lung disease through Wave 5 / 2020: ", primary_lung, "."),
    paste0("- Frailty-proxy adjusted incident chronic lung disease: ", frailty_lung, "."),
    paste0("- Discrete-time pooled model: ", discrete_lung, "."),
    paste0("- All-cause death through Wave 5 / 2020: ", death_key, "."),
    "",
    "## Boundary",
    "",
    "- Wave 5 is now included in the CHARLS main analysis, but exact death dates are still reserved for future sensitivity because current models use wave-interval timing.",
    "- Exit-module chronic lung disease is documented but not used as the primary incident CLD endpoint because it is decedent-proxy context."
  )
  writeLines(findings, file.path(log_dir, "charls_v0_6_wave5_findings.md"))

  log <- c(
    "# CHARLS V0.6 Wave 5 Model Log",
    "",
    paste0("- Input Harmonized core RDS: ", input_rds),
    paste0("- Input CHARLS 2020 archive: ", zip_path),
    "- Output model table: `results/tables/charls_v0_6_wave5_model_table.csv`.",
    "- Output key terms: `results/tables/charls_v0_6_wave5_key_terms.csv`.",
    "- Output codebook lock: `metadata/charls_wave5_codebook_lock_v0_6.csv`."
  )
  writeLines(log, file.path(log_dir, "charls_v0_6_wave5_model_log.md"))

  message("Wrote CHARLS v0.6 Wave 5-extended model outputs.")
}

if (sys.nframe() == 0) {
  main()
}
