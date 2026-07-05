# Validate NHANES V30 linked mortality model outputs.

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

check_file <- function(path, label) {
  exists <- file.exists(path)
  size <- if (exists) file.info(path)$size else NA_real_
  tibble(
    check = paste0("file_exists_nonempty:", label),
    status = ifelse(exists && size > 0, "pass", "fail"),
    observed = ifelse(exists, as.character(size), "missing"),
    expected = ">0 bytes"
  )
}

check_scalar <- function(check, observed, expected, passed) {
  tibble(
    check = check,
    status = ifelse(passed, "pass", "fail"),
    observed = as.character(observed),
    expected = expected
  )
}

scan_for_absolute_user_paths <- function(paths) {
  hits <- character()
  for (p in paths[file.exists(paths)]) {
    txt <- paste(readLines(p, warn = FALSE), collapse = "\n")
    if (grepl("${HOME}", txt, fixed = TRUE)) {
      hits <- c(hits, basename(p))
    }
  }
  hits
}

root <- find_project_root()
results_dir <- file.path(root, "results", "tables")
logs_dir <- file.path(root, "results", "logs")
figures_dir <- file.path(root, "results", "figures")
manuscript_dir <- file.path(root, "manuscript")

models_path <- file.path(results_dir, paste0("nhanes_mortality_models_", version_id, ".csv"))
cohorts_path <- file.path(results_dir, paste0("nhanes_mortality_model_cohort_counts_", version_id, ".csv"))
rates_path <- file.path(results_dir, paste0("nhanes_mortality_rates_by_quartile_", version_id, ".csv"))
spline_path <- file.path(results_dir, paste0("nhanes_mortality_spline_source_", version_id, ".csv"))
limitations_path <- file.path(manuscript_dir, paste0("nhanes_mortality_limitations_", version_id, ".md"))
model_log_path <- file.path(logs_dir, paste0("nhanes_mortality_models_", version_id, ".md"))
validation_csv_path <- file.path(results_dir, paste0("nhanes_mortality_models_validation_", version_id, ".csv"))
validation_log_path <- file.path(logs_dir, paste0("nhanes_mortality_models_validation_", version_id, ".md"))

figure_paths <- c(
  file.path(figures_dir, paste0("nhanes_mortality_model_forest_", version_id, ".png")),
  file.path(figures_dir, paste0("nhanes_mortality_model_forest_", version_id, ".pdf")),
  file.path(figures_dir, paste0("nhanes_mortality_spline_", version_id, ".png")),
  file.path(figures_dir, paste0("nhanes_mortality_spline_", version_id, ".pdf")),
  file.path(figures_dir, paste0("nhanes_mortality_rates_by_quartile_", version_id, ".png")),
  file.path(figures_dir, paste0("nhanes_mortality_rates_by_quartile_", version_id, ".pdf"))
)

checks <- bind_rows(
  check_file(models_path, "models_csv"),
  check_file(cohorts_path, "cohorts_csv"),
  check_file(rates_path, "rates_csv"),
  check_file(spline_path, "spline_source_csv"),
  check_file(limitations_path, "limitations_md"),
  check_file(model_log_path, "model_log_md"),
  bind_rows(lapply(seq_along(figure_paths), function(i) check_file(figure_paths[[i]], paste0("figure_", i))))
)

if (!file.exists(models_path)) {
  stop("Cannot validate without model table.", call. = FALSE)
}
models <- readr::read_csv(models_path, show_col_types = FALSE)
cohorts <- readr::read_csv(cohorts_path, show_col_types = FALSE)
rates <- readr::read_csv(rates_path, show_col_types = FALSE)
spline <- readr::read_csv(spline_path, show_col_types = FALSE)

primary <- models %>% filter(analysis_id == "t2_allcause_m2_primary", term == "resp_vulnerability_z")
q4 <- models %>% filter(analysis_id == "t2_allcause_quartile_categories", term == "rv_quartile_fQ4")
frailty <- models %>% filter(analysis_id == "t2_allcause_m3_frailty_adjusted", term == "resp_vulnerability_z")
clrd <- models %>% filter(analysis_id == "t2_clrd_mini_exploratory", term == "resp_vulnerability_z")
errors <- models %>% filter(status != "ok")
absolute_path_hits <- scan_for_absolute_user_paths(c(models_path, cohorts_path, rates_path, spline_path, limitations_path, model_log_path))

checks <- bind_rows(
  checks,
  check_scalar("no_model_errors", nrow(errors), "0", nrow(errors) == 0),
  check_scalar("primary_model_present", nrow(primary), "1", nrow(primary) == 1),
  check_scalar("primary_model_n", ifelse(nrow(primary) == 1, primary$n, NA), ">=6000", nrow(primary) == 1 && primary$n >= 6000),
  check_scalar("primary_model_events", ifelse(nrow(primary) == 1, primary$events, NA), ">=800", nrow(primary) == 1 && primary$events >= 800),
  check_scalar("primary_hr_direction", ifelse(nrow(primary) == 1, primary$hr, NA), ">1", nrow(primary) == 1 && primary$hr > 1),
  check_scalar("primary_hr_significant", ifelse(nrow(primary) == 1, primary$p_value, NA), "<0.05", nrow(primary) == 1 && primary$p_value < 0.05),
  check_scalar("q4_model_present", nrow(q4), "1", nrow(q4) == 1),
  check_scalar("q4_hr_direction", ifelse(nrow(q4) == 1, q4$hr, NA), ">1", nrow(q4) == 1 && q4$hr > 1),
  check_scalar("frailty_sensitivity_present", nrow(frailty), "1", nrow(frailty) == 1),
  check_scalar("clrd_exploratory_present", nrow(clrd), "1", nrow(clrd) == 1),
  check_scalar("rate_table_has_two_outcomes", paste(unique(rates$outcome_label), collapse = ", "), "All-cause and CLRD", length(unique(rates$outcome_label)) == 2),
  check_scalar("rate_table_has_four_quartiles", paste(sort(unique(rates$rv_quartile)), collapse = ","), "Q1,Q2,Q3,Q4", identical(paste(sort(unique(rates$rv_quartile)), collapse = ","), "Q1,Q2,Q3,Q4")),
  check_scalar("spline_source_dense", nrow(spline), ">=100", nrow(spline) >= 100),
  check_scalar("cohort_table_nonempty", nrow(cohorts), ">0", nrow(cohorts) > 0),
  check_scalar("no_absolute_user_paths_in_public_outputs", paste(absolute_path_hits, collapse = ","), "none", length(absolute_path_hits) == 0)
)

all_passed <- all(checks$status == "pass")
readr::write_csv(checks, validation_csv_path)

if (all_passed) {
  ticket_path <- file.path(results_dir, "v30_high_impact_ticket_audit.csv")
  if (file.exists(ticket_path)) {
    ticket_audit <- readr::read_csv(ticket_path, show_col_types = FALSE)
    ticket_audit <- ticket_audit %>%
      mutate(status = ifelse(ticket == "T2_NHANES_mortality_models", "completed", status))
    readr::write_csv(ticket_audit, ticket_path)
  }
}

log_lines <- c(
  "# NHANES V30 Mortality Models Validation",
  "",
  paste0("Version: ", version_id),
  "",
  paste0("Overall status: ", ifelse(all_passed, "PASS", "FAIL")),
  "",
  "## Key Model Checks",
  "",
  paste0("- Primary M2 HR: ", ifelse(nrow(primary) == 1, primary$hr_ci, "missing"), ", p=", ifelse(nrow(primary) == 1, primary$p_value_formatted, "missing")),
  paste0("- Q4 versus Q1 HR: ", ifelse(nrow(q4) == 1, q4$hr_ci, "missing"), ", p=", ifelse(nrow(q4) == 1, q4$p_value_formatted, "missing")),
  paste0("- Frailty-adjusted sensitivity HR: ", ifelse(nrow(frailty) == 1, frailty$hr_ci, "missing"), ", p=", ifelse(nrow(frailty) == 1, frailty$p_value_formatted, "missing")),
  paste0("- CLRD exploratory HR: ", ifelse(nrow(clrd) == 1, clrd$hr_ci, "missing"), ", p=", ifelse(nrow(clrd) == 1, clrd$p_value_formatted, "missing")),
  "",
  "## Checks",
  "",
  paste0("| ", paste(names(checks), collapse = " | "), " |"),
  paste0("| ", paste(rep("---", ncol(checks)), collapse = " | "), " |"),
  apply(checks, 1, function(x) paste0("| ", paste(x, collapse = " | "), " |"))
)
writeLines(log_lines, validation_log_path)

message("NHANES mortality model validation status: ", ifelse(all_passed, "PASS", "FAIL"))
if (!all_passed) {
  print(checks %>% filter(status == "fail"))
  stop("Validation failed.", call. = FALSE)
}
