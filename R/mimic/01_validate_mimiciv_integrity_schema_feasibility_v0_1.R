# Validate MIMIC-IV v0.1 integrity/schema feasibility outputs.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
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

require_file <- function(path, min_size = 1) {
  if (!file.exists(path)) stop("Missing required output: ", path, call. = FALSE)
  size <- file.info(path)$size
  if (is.na(size) || size < min_size) stop("Output is empty or too small: ", path, call. = FALSE)
  path
}

assert_true <- function(x, message) {
  if (!isTRUE(x)) stop(message, call. = FALSE)
}

main <- function() {
  root <- find_project_root()
  table_dir <- file.path(root, "results", "tables")
  log_dir <- file.path(root, "results", "logs")
  metadata_dir <- file.path(root, "metadata")

  required <- c(
    file.path(metadata_dir, "mimiciv_schema_lock_targets_v0_1.csv"),
    file.path(table_dir, "mimiciv_v0_1_download_manifest.csv"),
    file.path(table_dir, "mimiciv_v0_1_schema_columns.csv"),
    file.path(table_dir, "mimiciv_v0_1_schema_expectations.csv"),
    file.path(table_dir, "mimiciv_v0_1_d_items_respiratory_candidates.csv"),
    file.path(table_dir, "mimiciv_v0_1_endpoint_feasibility.csv"),
    file.path(log_dir, "mimiciv_integrity_schema_feasibility_v0_1.md")
  )
  invisible(lapply(required, require_file, min_size = 1))

  targets <- readr::read_csv(file.path(metadata_dir, "mimiciv_schema_lock_targets_v0_1.csv"), show_col_types = FALSE)
  manifest <- readr::read_csv(file.path(table_dir, "mimiciv_v0_1_download_manifest.csv"), show_col_types = FALSE)
  columns <- readr::read_csv(file.path(table_dir, "mimiciv_v0_1_schema_columns.csv"), show_col_types = FALSE)
  schema <- readr::read_csv(file.path(table_dir, "mimiciv_v0_1_schema_expectations.csv"), show_col_types = FALSE)
  candidates <- readr::read_csv(file.path(table_dir, "mimiciv_v0_1_d_items_respiratory_candidates.csv"), show_col_types = FALSE)
  endpoint <- readr::read_csv(file.path(table_dir, "mimiciv_v0_1_endpoint_feasibility.csv"), show_col_types = FALSE)
  log_text <- paste(readLines(file.path(log_dir, "mimiciv_integrity_schema_feasibility_v0_1.md"), warn = FALSE), collapse = "\n")

  expected_files <- c(
    "hosp/patients.csv.gz",
    "hosp/admissions.csv.gz",
    "hosp/transfers.csv.gz",
    "icu/icustays.csv.gz",
    "icu/d_items.csv.gz",
    "icu/chartevents.csv.gz",
    "icu/procedureevents.csv.gz"
  )
  missing_manifest <- setdiff(expected_files, manifest$file)
  assert_true(length(missing_manifest) == 0, paste("Manifest missing required MIMIC files:", paste(missing_manifest, collapse = "; ")))
  assert_true(nrow(targets) >= 12, "Schema lock target table is too small.")
  assert_true(nrow(schema) == nrow(targets), "Schema expectation rows do not match targets.")
  assert_true(nrow(endpoint) >= 6, "Endpoint feasibility table is too small.")
  assert_true(any(endpoint$feasibility_component == "no_modeling_gate" & endpoint$status == "enforced"), "No-modeling gate is missing.")
  assert_true(all(!str_detect(manifest$file, "^/")), "Manifest leaked absolute file paths.")
  assert_true(all(!str_detect(schema$table_file, "^/")), "Schema expectations leaked absolute file paths.")
  assert_true(str_detect(log_text, fixed("no row-level clinical exports")), "Log lacks row-level data boundary.")

  pending <- any(manifest$download_status != "present_not_active")
  if (pending) {
    assert_true(any(endpoint$status == "pending"), "Pending download state should produce pending feasibility rows.")
    assert_true(str_detect(log_text, fixed("STATUS: pending_download")), "Log should mark pending download.")
  } else {
    assert_true(all(manifest$gzip_status %in% c("pass", "not_applicable")), "Completed state has non-pass gzip statuses.")
    assert_true(all(manifest$sha256_status == "pass"), "Completed state has non-pass SHA256 statuses.")
    assert_true(str_detect(log_text, fixed("ready_for_full_gzip_sha_validation_or_schema_review")), "Log should mark ready state.")
  }

  validation <- c(
    "# MIMIC-IV Integrity And Schema Feasibility V0.1 Validation",
    "",
    paste0("- Run date: ", Sys.Date()),
    paste0("- Checked manifest rows: ", nrow(manifest), "."),
    paste0("- Checked schema target rows: ", nrow(targets), "."),
    paste0("- Checked schema column rows: ", nrow(columns), "."),
    paste0("- Checked d_items respiratory candidate rows: ", nrow(candidates), "."),
    paste0("- Checked endpoint feasibility rows: ", nrow(endpoint), "."),
    paste0("- Current state: ", ifelse(pending, "pending_download", "download_present_and_validated"), "."),
    "- Checked no absolute file paths and no row-level export boundary.",
    "",
    "PASS: MIMIC-IV integrity/schema feasibility V0.1 validation complete."
  )
  writeLines(validation, file.path(log_dir, "mimiciv_integrity_schema_feasibility_v0_1_validation.md"))
  message("PASS: MIMIC-IV integrity/schema feasibility V0.1 validation complete.")
}

if (sys.nframe() == 0) {
  main()
}
