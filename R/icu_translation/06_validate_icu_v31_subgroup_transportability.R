# Validate ICU V31 subgroup transportability outputs.

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

tab <- read_csv(file.path(tab_dir, "icu_v31_subgroup_transportability.csv"), show_col_types = FALSE)
note <- paste(readLines(file.path(manuscript_dir, "icu_subgroup_transportability_note_v31_5.md"), warn = FALSE), collapse = "\n")

checks <- bind_rows(
  check_scalar("has_hfnc_niv", any(tab$requested_domain == "HFNC/NIV support class"), "TRUE", any(tab$requested_domain == "HFNC/NIV support class")),
  check_scalar("has_age", any(tab$requested_domain == "Age group"), "TRUE", any(tab$requested_domain == "Age group")),
  check_scalar("has_sex", any(tab$requested_domain == "Sex"), "TRUE", any(tab$requested_domain == "Sex")),
  check_scalar("records_chronic_pulmonary_blocker", any(tab$requested_domain == "Chronic pulmonary disease" & tab$availability_status == "requires_row_level_rerun"), "TRUE", any(tab$requested_domain == "Chronic pulmonary disease" & tab$availability_status == "requires_row_level_rerun")),
  check_scalar("records_hypoxemic_hypercapnic_blocker", any(tab$requested_domain == "Hypoxemic/hypercapnic phenotype" & tab$availability_status == "requires_row_level_rerun"), "TRUE", any(tab$requested_domain == "Hypoxemic/hypercapnic phenotype" & tab$availability_status == "requires_row_level_rerun")),
  check_scalar("hfnc_sparse_boundary_present", str_detect(note, fixed("sparse")), "TRUE", str_detect(note, fixed("sparse"))),
  check_scalar("no_hfnc_specific_validation_claim", str_detect(note, fixed("HFNC-specific external validation is not claimed")), "TRUE", str_detect(note, fixed("HFNC-specific external validation is not claimed"))),
  check_scalar("tte_no_go_boundary_present", str_detect(note, fixed("TTE feasibility layer remains no-go")), "TRUE", str_detect(note, fixed("TTE feasibility layer remains no-go")))
)

all_passed <- all(checks$status == "pass")
write_csv(checks, file.path(tab_dir, "icu_v31_subgroup_transportability_validation.csv"))
writeLines(c(
  "# ICU V31 Subgroup Transportability Validation",
  "",
  paste0("Overall status: ", ifelse(all_passed, "PASS", "FAIL")),
  "",
  md_table(checks)
), file.path(log_dir, "icu_v31_subgroup_transportability_validation.md"))

message("ICU V31 subgroup transportability validation status: ", ifelse(all_passed, "PASS", "FAIL"))
if (!all_passed) {
  print(checks %>% filter(status == "fail"))
  stop("Validation failed.", call. = FALSE)
}
