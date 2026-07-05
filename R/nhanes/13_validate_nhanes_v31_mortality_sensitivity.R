# Validate NHANES V31 mortality sensitivity outputs.

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

check_scalar <- function(check, observed, expected, passed) {
  tibble(check = check, status = ifelse(isTRUE(passed), "pass", "fail"), observed = as.character(observed), expected = expected)
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
manuscript_dir <- file.path(root, "manuscript")

sens_path <- file.path(tab_dir, "nhanes_v31_mortality_sensitivity.csv")
ph_path <- file.path(tab_dir, "nhanes_v31_ph_diagnostics.csv")
note_path <- file.path(manuscript_dir, "nhanes_mortality_sensitivity_v31_5.md")

sens <- read_csv(sens_path, show_col_types = FALSE)
ph <- read_csv(ph_path, show_col_types = FALSE)
note <- paste(readLines(note_path, warn = FALSE), collapse = "\n")

checks <- bind_rows(
  check_scalar("sensitivity_file_nonempty", nrow(sens), ">=8", nrow(sens) >= 8),
  check_scalar("ph_diagnostics_present", nrow(ph), ">=4", nrow(ph) >= 4),
  check_scalar("has_primary_reference", any(sens$check_type == "primary_reference"), "TRUE", any(sens$check_type == "primary_reference")),
  check_scalar("has_ph_check", any(sens$check_type == "proportional_hazards"), "TRUE", any(sens$check_type == "proportional_hazards")),
  check_scalar("has_race_calibrated", any(sens$check_type == "race_calibrated"), "TRUE", any(sens$check_type == "race_calibrated")),
  check_scalar("has_cycle_stratified", sum(sens$check_type == "cycle_stratified"), ">=3", sum(sens$check_type == "cycle_stratified") >= 3),
  check_scalar("has_frailty_complete_case", any(sens$check_type == "frailty_complete_case"), "TRUE", any(sens$check_type == "frailty_complete_case")),
  check_scalar("clrd_marked_exploratory", any(sens$check_type == "clrd_exploratory" & str_detect(sens$boundary, regex("Exploratory", ignore_case = TRUE))), "TRUE", any(sens$check_type == "clrd_exploratory" & str_detect(sens$boundary, regex("Exploratory", ignore_case = TRUE)))),
  check_scalar("note_keeps_clrd_out_of_abstract", str_detect(note, fixed("should not be placed in the abstract")), "TRUE", str_detect(note, fixed("should not be placed in the abstract"))),
  check_scalar("note_frailty_sensitivity_boundary", str_detect(note, fixed("sensitivity only")), "TRUE", str_detect(note, fixed("sensitivity only")))
)

all_passed <- all(checks$status == "pass")
write_csv(checks, file.path(tab_dir, "nhanes_v31_mortality_sensitivity_validation.csv"))
writeLines(c(
  "# NHANES V31 Mortality Sensitivity Validation",
  "",
  paste0("Overall status: ", ifelse(all_passed, "PASS", "FAIL")),
  "",
  md_table(checks)
), file.path(log_dir, "nhanes_v31_mortality_sensitivity_validation.md"))

message("NHANES V31 mortality sensitivity validation status: ", ifelse(all_passed, "PASS", "FAIL"))
if (!all_passed) {
  print(checks %>% filter(status == "fail"))
  stop("Validation failed.", call. = FALSE)
}
