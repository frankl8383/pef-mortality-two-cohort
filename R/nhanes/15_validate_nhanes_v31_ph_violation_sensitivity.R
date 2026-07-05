# Validate NHANES V31.6 PH-violation sensitivity outputs.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(tibble)
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

dir_create <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
}

check_scalar <- function(check, observed, expected, passed) {
  tibble(
    check = check,
    status = ifelse(isTRUE(passed), "pass", "fail"),
    observed = as.character(observed),
    expected = expected
  )
}

md_table <- function(dat) {
  dat <- as.data.frame(dat, stringsAsFactors = FALSE)
  dat[] <- lapply(dat, function(x) ifelse(is.na(x), "", as.character(x)))
  c(
    paste0("| ", paste(names(dat), collapse = " | "), " |"),
    paste0("| ", paste(rep("---", ncol(dat)), collapse = " | "), " |"),
    apply(dat, 1, function(x) paste0("| ", paste(x, collapse = " | "), " |"))
  )
}

root <- find_project_root()
tab_dir <- file.path(root, "results", "tables")
log_dir <- file.path(root, "results", "logs")
figure_dir <- file.path(root, "results", "figures")
manuscript_dir <- file.path(root, "manuscript")
dir_create(log_dir)

time_path <- file.path(tab_dir, "nhanes_v31_6_time_interaction.csv")
piecewise_path <- file.path(tab_dir, "nhanes_v31_6_piecewise_followup_sensitivity.csv")
fixed_model_path <- file.path(tab_dir, "nhanes_v31_6_fixed_time_mortality_models.csv")
fixed_rate_path <- file.path(tab_dir, "nhanes_v31_6_fixed_time_mortality_rates.csv")
note_path <- file.path(manuscript_dir, "nhanes_ph_violation_sensitivity_v31_6.md")
brief_path <- file.path(manuscript_dir, "literature_methods_brief_v31_6.md")
fig_png <- file.path(figure_dir, "nhanes_v31_6_fixed_time_mortality_by_quartile.png")
fig_pdf <- file.path(figure_dir, "nhanes_v31_6_fixed_time_mortality_by_quartile.pdf")

time_interaction <- read_csv(time_path, show_col_types = FALSE)
piecewise <- read_csv(piecewise_path, show_col_types = FALSE)
fixed_models <- read_csv(fixed_model_path, show_col_types = FALSE)
fixed_rates <- read_csv(fixed_rate_path, show_col_types = FALSE)
note <- paste(readLines(note_path, warn = FALSE), collapse = "\n")
brief <- paste(readLines(brief_path, warn = FALSE), collapse = "\n")

checks <- bind_rows(
  check_scalar("time_interaction_has_base_and_tt_terms", nrow(time_interaction), ">=2", nrow(time_interaction) >= 2),
  check_scalar("time_interaction_has_tt_term", any(time_interaction$term == "tt(resp_vulnerability_z)"), "TRUE", any(time_interaction$term == "tt(resp_vulnerability_z)")),
  check_scalar("piecewise_three_intervals", paste(piecewise$interval, collapse = ", "), "0-5, 5-10, 10+", all(c("0-5 years", "5-10 years", "10+ years") %in% piecewise$interval)),
  check_scalar("piecewise_no_error_status", paste(unique(piecewise$status), collapse = ", "), "no error", !any(str_detect(piecewise$status, "error"))),
  check_scalar("fixed_models_two_horizons", paste(fixed_models$horizon_years, collapse = ", "), "5 and 10", all(c(5, 10) %in% fixed_models$horizon_years)),
  check_scalar("fixed_models_ok", paste(unique(fixed_models$status), collapse = ", "), "ok", all(fixed_models$status == "ok")),
  check_scalar("fixed_rates_two_horizons_four_quartiles", nrow(fixed_rates), ">=8", nrow(fixed_rates) >= 8),
  check_scalar("fixed_rates_all_bounded", paste(range(fixed_rates$weighted_risk, na.rm = TRUE), collapse = "-"), "0 to 1", all(fixed_rates$weighted_risk >= 0 & fixed_rates$weighted_risk <= 1, na.rm = TRUE)),
  check_scalar("manuscript_notes_average_association", str_detect(note, fixed("average association")), "TRUE", str_detect(note, fixed("average association"))),
  check_scalar("manuscript_warns_no_constant_hr", str_detect(note, fixed("constant hazard ratio")), "TRUE", str_detect(note, fixed("constant hazard ratio"))),
  check_scalar("literature_brief_present", nchar(brief), ">500 chars", nchar(brief) > 500),
  check_scalar("figure_png_exists", file.exists(fig_png), "TRUE", file.exists(fig_png)),
  check_scalar("figure_pdf_exists", file.exists(fig_pdf), "TRUE", file.exists(fig_pdf))
)

all_passed <- all(checks$status == "pass")
write_csv(checks, file.path(tab_dir, "nhanes_v31_6_ph_violation_sensitivity_validation.csv"))
writeLines(c(
  "# NHANES V31.6 PH-Violation Sensitivity Validation",
  "",
  paste0("Overall status: ", ifelse(all_passed, "PASS", "FAIL")),
  "",
  md_table(checks)
), file.path(log_dir, "nhanes_v31_6_ph_violation_sensitivity_validation.md"))

message("NHANES V31.6 PH-violation sensitivity validation status: ", ifelse(all_passed, "PASS", "FAIL"))
if (!all_passed) {
  print(checks %>% filter(status == "fail"))
  stop("Validation failed.", call. = FALSE)
}
