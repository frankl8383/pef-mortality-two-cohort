# Validate ICU V31.6 row-level subgroup transportability outputs.

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
manuscript_dir <- file.path(root, "manuscript")
derived_dir <- file.path(root, "derived_sensitive", "icu_translation")
dir_create(log_dir)

perf_path <- file.path(tab_dir, "icu_v31_6_row_level_subgroup_performance.csv")
counts_path <- file.path(tab_dir, "icu_v31_6_row_level_subgroup_flag_counts.csv")
ledger_path <- file.path(tab_dir, "icu_v31_6_row_level_subgroup_extraction_ledger.csv")
note_path <- file.path(manuscript_dir, "icu_row_level_subgroup_transportability_note_v31_6.md")
aug_path <- file.path(derived_dir, "icu_v31_6_row_level_subgroup_augmented.csv.gz")

perf <- read_csv(perf_path, show_col_types = FALSE)
counts <- read_csv(counts_path, show_col_types = FALSE)
ledger <- read_csv(ledger_path, show_col_types = FALSE)
note <- paste(readLines(note_path, warn = FALSE), collapse = "\n")

checks <- bind_rows(
  check_scalar("augmented_row_level_file_exists", file.exists(aug_path), "TRUE", file.exists(aug_path)),
  check_scalar("has_chronic_pulmonary_domain", any(perf$subgroup_domain == "coded_chronic_pulmonary_disease"), "TRUE", any(perf$subgroup_domain == "coded_chronic_pulmonary_disease")),
  check_scalar("chronic_pulmonary_has_present_absent", paste(unique(perf$subgroup_level[perf$subgroup_domain == "coded_chronic_pulmonary_disease"]), collapse = ", "), "present and absent", all(c("chronic_pulmonary_present", "chronic_pulmonary_absent") %in% perf$subgroup_level[perf$subgroup_domain == "coded_chronic_pulmonary_disease"])),
  check_scalar("has_oxygenation_stress_domain", any(perf$subgroup_domain == "common_feature_oxygenation_stress"), "TRUE", any(perf$subgroup_domain == "common_feature_oxygenation_stress")),
  check_scalar("has_hypercapnia_domain", any(perf$subgroup_domain == "coarse_blood_gas_hypercapnia"), "TRUE", any(perf$subgroup_domain == "coarse_blood_gas_hypercapnia")),
  check_scalar("has_low_po2_domain", any(perf$subgroup_domain == "coarse_blood_gas_low_po2"), "TRUE", any(perf$subgroup_domain == "coarse_blood_gas_low_po2")),
  check_scalar("both_evaluations_present", paste(unique(perf$evaluation), collapse = ", "), "mimic_internal_test and eicu_external", all(c("mimic_internal_test", "eicu_external") %in% perf$evaluation)),
  check_scalar("counts_have_two_evaluations", nrow(counts), ">=2", nrow(counts) >= 2),
  check_scalar("chronic_blocker_resolved_in_ledger", ledger$status[ledger$domain == "Chronic pulmonary disease"], "estimated", any(ledger$domain == "Chronic pulmonary disease" & ledger$status == "estimated")),
  check_scalar("hypercapnia_boundary_present", any(str_detect(ledger$boundary, fixed("not baseline-window precise"))), "TRUE", any(str_detect(ledger$boundary, fixed("not baseline-window precise")))),
  check_scalar("note_blocks_adjudicated_copd", str_detect(note, fixed("not adjudicated COPD")), "TRUE", str_detect(note, fixed("not adjudicated COPD"))),
  check_scalar("note_blocks_treatment_effect", str_detect(note, fixed("treatment-effect")), "TRUE", str_detect(note, fixed("treatment-effect"))),
  check_scalar("note_blocks_direct_pef_validation", str_detect(note, fixed("direct validation of the population PEF marker")), "TRUE", str_detect(note, fixed("direct validation of the population PEF marker"))),
  check_scalar("metric_bounds_ok", "checked", "0<=rates/metrics<=1 where nonmissing", all(perf$event_rate >= 0 & perf$event_rate <= 1, na.rm = TRUE) && all(perf$mean_predicted_risk >= 0 & perf$mean_predicted_risk <= 1, na.rm = TRUE) && all(perf$auroc >= 0 & perf$auroc <= 1, na.rm = TRUE) && all(perf$auprc >= 0 & perf$auprc <= 1, na.rm = TRUE) && all(perf$brier >= 0 & perf$brier <= 1, na.rm = TRUE))
)

all_passed <- all(checks$status == "pass")
write_csv(checks, file.path(tab_dir, "icu_v31_6_row_level_subgroup_validation.csv"))
writeLines(c(
  "# ICU V31.6 Row-Level Subgroup Transportability Validation",
  "",
  paste0("Overall status: ", ifelse(all_passed, "PASS", "FAIL")),
  "",
  md_table(checks)
), file.path(log_dir, "icu_v31_6_row_level_subgroup_validation.md"))

message("ICU V31.6 row-level subgroup validation status: ", ifelse(all_passed, "PASS", "FAIL"))
if (!all_passed) {
  print(checks %>% filter(status == "fail"))
  stop("Validation failed.", call. = FALSE)
}
