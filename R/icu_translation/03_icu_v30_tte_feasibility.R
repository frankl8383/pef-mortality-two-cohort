# ICU V30 target-trial emulation feasibility only.
# This script audits whether HFNC-vs-NIV-type causal emulation is plausible.
# It does not estimate any treatment effect.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
})

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

markdown_table <- function(dat) {
  dat <- as.data.frame(dat, stringsAsFactors = FALSE)
  dat[] <- lapply(dat, function(x) ifelse(is.na(x), "", as.character(x)))
  header <- paste0("| ", paste(names(dat), collapse = " | "), " |")
  separator <- paste0("| ", paste(rep("---", ncol(dat)), collapse = " | "), " |")
  rows <- apply(dat, 1, function(x) paste0("| ", paste(x, collapse = " | "), " |"))
  c(header, separator, rows)
}

weighted_mean <- function(x, w) {
  ok <- !is.na(x) & !is.na(w) & w > 0
  if (!any(ok)) return(NA_real_)
  sum(w[ok] * x[ok]) / sum(w[ok])
}

weighted_var <- function(x, w) {
  ok <- !is.na(x) & !is.na(w) & w > 0
  if (sum(ok) < 2) return(NA_real_)
  mu <- weighted_mean(x[ok], w[ok])
  sum(w[ok] * (x[ok] - mu)^2) / sum(w[ok])
}

smd_one <- function(data, var, weight = NULL) {
  a <- data$support_hfnc
  x <- data[[var]]
  if (is.null(weight)) {
    mt <- mean(x[a == 1], na.rm = TRUE)
    mc <- mean(x[a == 0], na.rm = TRUE)
    vt <- stats::var(x[a == 1], na.rm = TRUE)
    vc <- stats::var(x[a == 0], na.rm = TRUE)
  } else {
    w <- data[[weight]]
    mt <- weighted_mean(x[a == 1], w[a == 1])
    mc <- weighted_mean(x[a == 0], w[a == 0])
    vt <- weighted_var(x[a == 1], w[a == 1])
    vc <- weighted_var(x[a == 0], w[a == 0])
  }
  denom <- sqrt((vt + vc) / 2)
  if (is.na(denom) || denom == 0) return(NA_real_)
  (mt - mc) / denom
}

prepare_ps_data <- function(data, dataset_label, covariates) {
  df <- data %>%
    filter(support_class %in% c("hfnc", "niv"), !is.na(support_hfnc), !is.na(outcome_failure_48h_l2)) %>%
    mutate(dataset_label = dataset_label)
  for (v in covariates) {
    miss <- is.na(df[[v]])
    med <- stats::median(df[[v]], na.rm = TRUE)
    if (is.na(med)) med <- 0
    df[[paste0(v, "_missing")]] <- as.integer(miss)
    df[[paste0(v, "_imp")]] <- df[[v]]
    df[[paste0(v, "_imp")]][miss] <- med
  }
  df
}

fit_ps <- function(df, covariates) {
  imp_vars <- paste0(covariates, "_imp")
  miss_vars <- paste0(covariates, "_missing")
  rhs <- paste(c(imp_vars, miss_vars), collapse = " + ")
  formula <- stats::as.formula(paste0("support_hfnc ~ ", rhs))
  warnings <- character()
  fit <- withCallingHandlers(
    stats::glm(formula, data = df, family = stats::binomial()),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  ps <- as.numeric(stats::predict(fit, type = "response"))
  ps <- pmin(pmax(ps, 1e-6), 1 - 1e-6)
  p_treat <- mean(df$support_hfnc == 1)
  df$ps_hfnc <- ps
  df$iptw_stabilized <- ifelse(df$support_hfnc == 1, p_treat / ps, (1 - p_treat) / (1 - ps))
  cap <- stats::quantile(df$iptw_stabilized, 0.99, na.rm = TRUE, names = FALSE)
  df$iptw_stabilized_trunc99 <- pmin(df$iptw_stabilized, cap)
  list(data = df, warnings = warnings, formula = deparse(formula))
}

ess <- function(w) {
  if (length(w) == 0 || all(is.na(w))) return(NA_real_)
  (sum(w, na.rm = TRUE)^2) / sum(w^2, na.rm = TRUE)
}

counts_for_dataset <- function(df, dataset_label) {
  bind_rows(
    df %>%
      summarise(
        dataset = .env$dataset_label,
        group = "overall",
        n = n(),
        events = sum(outcome_failure_48h_l2 == 1L, na.rm = TRUE),
        event_rate = events / n,
        hfnc_n = sum(support_hfnc == 1L, na.rm = TRUE),
        niv_n = sum(support_hfnc == 0L, na.rm = TRUE)
      ),
    df %>%
      group_by(support_class) %>%
      summarise(
        dataset = .env$dataset_label,
        group = support_class[1],
        n = n(),
        events = sum(outcome_failure_48h_l2 == 1L, na.rm = TRUE),
        event_rate = events / n,
        hfnc_n = sum(support_hfnc == 1L, na.rm = TRUE),
        niv_n = sum(support_hfnc == 0L, na.rm = TRUE),
        .groups = "drop"
      )
  )
}

weight_diagnostics <- function(df, dataset_label) {
  ps_t <- df$ps_hfnc[df$support_hfnc == 1L]
  ps_c <- df$ps_hfnc[df$support_hfnc == 0L]
  common_low <- max(min(ps_t, na.rm = TRUE), min(ps_c, na.rm = TRUE))
  common_high <- min(max(ps_t, na.rm = TRUE), max(ps_c, na.rm = TRUE))
  in_common <- df$ps_hfnc >= common_low & df$ps_hfnc <= common_high
  tibble(
    dataset = dataset_label,
    n = nrow(df),
    hfnc_n = sum(df$support_hfnc == 1L),
    niv_n = sum(df$support_hfnc == 0L),
    ps_min_hfnc = min(ps_t, na.rm = TRUE),
    ps_p05_hfnc = stats::quantile(ps_t, 0.05, na.rm = TRUE, names = FALSE),
    ps_median_hfnc = stats::median(ps_t, na.rm = TRUE),
    ps_p95_hfnc = stats::quantile(ps_t, 0.95, na.rm = TRUE, names = FALSE),
    ps_max_hfnc = max(ps_t, na.rm = TRUE),
    ps_min_niv = min(ps_c, na.rm = TRUE),
    ps_p05_niv = stats::quantile(ps_c, 0.05, na.rm = TRUE, names = FALSE),
    ps_median_niv = stats::median(ps_c, na.rm = TRUE),
    ps_p95_niv = stats::quantile(ps_c, 0.95, na.rm = TRUE, names = FALSE),
    ps_max_niv = max(ps_c, na.rm = TRUE),
    common_support_low = common_low,
    common_support_high = common_high,
    proportion_in_common_support = mean(in_common, na.rm = TRUE),
    iptw_mean = mean(df$iptw_stabilized, na.rm = TRUE),
    iptw_p99 = stats::quantile(df$iptw_stabilized, 0.99, na.rm = TRUE, names = FALSE),
    iptw_max = max(df$iptw_stabilized, na.rm = TRUE),
    ess_overall = ess(df$iptw_stabilized),
    ess_hfnc = ess(df$iptw_stabilized[df$support_hfnc == 1L]),
    ess_niv = ess(df$iptw_stabilized[df$support_hfnc == 0L]),
    iptw_trunc99_max = max(df$iptw_stabilized_trunc99, na.rm = TRUE),
    ess_trunc99_overall = ess(df$iptw_stabilized_trunc99),
    go_no_go = dplyr::case_when(
      hfnc_n < 500 ~ "no_go_treated_sample_sparse",
      proportion_in_common_support < 0.90 ~ "no_go_poor_common_support",
      iptw_p99 > 10 ~ "no_go_unstable_weights",
      TRUE ~ "feasibility_possible_but_no_effect_estimated"
    )
  )
}

balance_table <- function(df, dataset_label, covariates) {
  balance_vars <- c(paste0(covariates, "_imp"), paste0(covariates, "_missing"))
  bind_rows(lapply(balance_vars, function(v) {
    tibble(
      dataset = dataset_label,
      variable = v,
      smd_unweighted = smd_one(df, v, NULL),
      smd_iptw = smd_one(df, v, "iptw_stabilized"),
      smd_iptw_trunc99 = smd_one(df, v, "iptw_stabilized_trunc99")
    )
  })) %>%
    mutate(
      abs_smd_unweighted = abs(smd_unweighted),
      abs_smd_iptw = abs(smd_iptw),
      abs_smd_iptw_trunc99 = abs(smd_iptw_trunc99)
    )
}

root <- find_project_root()
derived_dir <- file.path(root, "derived_sensitive", "icu_translation")
results_dir <- file.path(root, "results", "tables")
logs_dir <- file.path(root, "results", "logs")
manuscript_dir <- file.path(root, "manuscript")
dir_create(results_dir)
dir_create(logs_dir)
dir_create(manuscript_dir)

mimic <- readr::read_csv(file.path(derived_dir, "mimiciv_icu_translation_modeling_v22_0.csv.gz"), show_col_types = FALSE)
eicu <- readr::read_csv(file.path(derived_dir, "eicu_icu_translation_modeling_v22_0.csv.gz"), show_col_types = FALSE)
covariates <- c("age", "male", "spo2_last", "fio2_last", "respiratory_rate_last", "heart_rate_last", "oxygen_flow_max", "peep_max", "rox_l2")

mimic_ps <- fit_ps(prepare_ps_data(mimic, "MIMIC-IV", covariates), covariates)
eicu_ps <- fit_ps(prepare_ps_data(eicu, "eICU-CRD", covariates), covariates)

counts <- bind_rows(
  counts_for_dataset(mimic_ps$data, "MIMIC-IV"),
  counts_for_dataset(eicu_ps$data, "eICU-CRD")
) %>%
  mutate(across(where(is.numeric), ~ round(.x, 6)))

weights <- bind_rows(
  weight_diagnostics(mimic_ps$data, "MIMIC-IV"),
  weight_diagnostics(eicu_ps$data, "eICU-CRD")
) %>%
  mutate(across(where(is.numeric), ~ round(.x, 6)))

balance <- bind_rows(
  balance_table(mimic_ps$data, "MIMIC-IV", covariates),
  balance_table(eicu_ps$data, "eICU-CRD", covariates)
) %>%
  mutate(across(where(is.numeric), ~ round(.x, 6)))

balance_summary <- balance %>%
  group_by(dataset) %>%
  summarise(
    max_abs_smd_unweighted = max(abs_smd_unweighted, na.rm = TRUE),
    max_abs_smd_iptw = max(abs_smd_iptw, na.rm = TRUE),
    max_abs_smd_iptw_trunc99 = max(abs_smd_iptw_trunc99, na.rm = TRUE),
    variables_abs_smd_iptw_gt_0_1 = sum(abs_smd_iptw > 0.1, na.rm = TRUE),
    variables_abs_smd_trunc99_gt_0_1 = sum(abs_smd_iptw_trunc99 > 0.1, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(across(where(is.numeric), ~ round(.x, 6)))

go_no_go <- weights %>%
  select(dataset, hfnc_n, niv_n, proportion_in_common_support, iptw_p99, iptw_max, ess_hfnc, ess_niv) %>%
  left_join(balance_summary, by = "dataset") %>%
  mutate(
    sample_gate = hfnc_n >= 500 & niv_n >= 500,
    positivity_gate = proportion_in_common_support >= 0.90,
    weight_gate = iptw_p99 <= 10 & ess_hfnc >= 50 & ess_niv >= 50,
    balance_gate = max_abs_smd_iptw_trunc99 <= 0.10 & variables_abs_smd_trunc99_gt_0_1 == 0,
    final_decision = dplyr::case_when(
      !sample_gate ~ "no_go_effect_estimation_sample_sparse",
      !positivity_gate ~ "no_go_effect_estimation_poor_common_support",
      !weight_gate ~ "no_go_effect_estimation_unstable_weights_or_low_ess",
      !balance_gate ~ "no_go_effect_estimation_balance_not_acceptable",
      TRUE ~ "protocol_review_required_before_any_effect_estimation"
    ),
    interpretation = dplyr::case_when(
      final_decision == "no_go_effect_estimation_balance_not_acceptable" ~ "Feasibility counts are adequate, but preliminary weighting does not achieve acceptable balance; do not estimate causal effects.",
      final_decision == "protocol_review_required_before_any_effect_estimation" ~ "Diagnostics pass screening, but treatment timing and confounder validity still require clinical protocol review.",
      TRUE ~ "Do not estimate causal effects from current feasibility layer."
    )
  )

protocol <- tibble::tribble(
  ~component, ~candidate_specification, ~feasibility_judgement,
  "Eligibility", "Adult first eligible ICU stays/patient-units receiving HFNC or NIV-type support at the locked L2 landmark.", "Feasible for counting in MIMIC and eICU; exact eICU HFNC mapping remains sparse.",
  "Treatment strategies", "HFNC-type support versus NIV-type support at support-class lock.", "MIMIC has both groups; eICU HFNC group is much smaller than NIV.",
  "Time zero", "Locked support initiation / L2 landmark from ICU translation pipeline.", "Acceptable for feasibility; timeline needs clinical review before causal emulation.",
  "Outcome", "Invasive ventilation or death within 48h from L2.", "Feasible as already used in ICU translation, but mapping differs across MIMIC/eICU.",
  "Confounding adjustment", "Age, sex, oxygenation/vital signs, oxygen flow, PEEP, ROX-like marker, missing indicators.", "Only preliminary; some variables are landmark/dynamic and may be post-treatment, so causal interpretation is unsafe.",
  "Positivity", "Propensity-score overlap and stabilized IPTW diagnostics.", "Use diagnostics only; stop before effect estimation if common support or weights are poor.",
  "Estimand", "Not estimated in V30 T5.", "No causal effect should be reported from this feasibility step."
)

counts_path <- file.path(results_dir, paste0("tte_feasibility_counts_", version_id, ".csv"))
weights_path <- file.path(results_dir, paste0("tte_weight_diagnostics_", version_id, ".csv"))
balance_path <- file.path(results_dir, paste0("tte_balance_preliminary_", version_id, ".csv"))
balance_summary_path <- file.path(results_dir, paste0("tte_balance_summary_", version_id, ".csv"))
protocol_path <- file.path(results_dir, paste0("tte_protocol_table_", version_id, ".csv"))
go_no_go_path <- file.path(results_dir, paste0("tte_go_no_go_", version_id, ".csv"))
readr::write_csv(counts, counts_path)
readr::write_csv(weights, weights_path)
readr::write_csv(balance, balance_path)
readr::write_csv(balance_summary, balance_summary_path)
readr::write_csv(protocol, protocol_path)
readr::write_csv(go_no_go, go_no_go_path)

note_path <- file.path(manuscript_dir, paste0("tte_feasibility_note_", version_id, ".md"))
note_lines <- c(
  "# ICU Target-Trial Emulation Feasibility Note",
  "",
  paste0("Version: ", version_id),
  "",
  "This T5 output is a feasibility audit only. It does not estimate the effect of HFNC versus NIV and should not be cited as causal evidence.",
  "",
  "## Counts",
  "",
  markdown_table(counts %>% select(dataset, group, n, events, event_rate, hfnc_n, niv_n)),
  "",
  "## Positivity and Weights",
  "",
  markdown_table(weights %>% select(dataset, hfnc_n, niv_n, proportion_in_common_support, iptw_p99, iptw_max, ess_overall, ess_hfnc, ess_niv, go_no_go)),
  "",
  "## Balance Summary",
  "",
  markdown_table(balance_summary),
  "",
  "## Go/No-Go Decision",
  "",
  markdown_table(go_no_go %>% select(dataset, sample_gate, positivity_gate, weight_gate, balance_gate, final_decision, interpretation)),
  "",
  "## Interpretation",
  "",
  "- MIMIC has enough HFNC and NIV rows for a feasibility exercise, but confounding and treatment-timing review remain necessary.",
  "- eICU HFNC rows are sparse relative to NIV, so exact HFNC-vs-NIV causal emulation is fragile.",
  "- Several candidate adjustment variables are L2 dynamic features and may be affected by treatment; they are useful for balance diagnostics but unsafe as a finalized causal adjustment set without protocol review.",
  "- Preliminary balance remains unacceptable after the screened IPTW approach, so no treatment-effect estimate should be run or reported from T5."
)
writeLines(note_lines, note_path)

log_path <- file.path(logs_dir, paste0("tte_feasibility_", version_id, ".md"))
log_lines <- c(
  "# ICU TTE Feasibility Log",
  "",
  paste0("Version: ", version_id),
  "",
  "T5 feasibility-only diagnostics completed. No causal effect was estimated.",
  "",
  "Outputs:",
  "",
  "- `results/tables/tte_feasibility_counts_v30_0.csv`",
  "- `results/tables/tte_weight_diagnostics_v30_0.csv`",
  "- `results/tables/tte_balance_preliminary_v30_0.csv`",
  "- `results/tables/tte_balance_summary_v30_0.csv`",
  "- `results/tables/tte_protocol_table_v30_0.csv`",
  "- `results/tables/tte_go_no_go_v30_0.csv`",
  "- `manuscript/tte_feasibility_note_v30_0.md`"
)
writeLines(log_lines, log_path)

message("ICU TTE feasibility outputs complete.")
