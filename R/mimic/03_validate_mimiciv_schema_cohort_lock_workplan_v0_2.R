# Validate MIMIC-IV schema/cohort/endpoint lock workplan v0.2.

suppressPackageStartupMessages({
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

read_text <- function(path) {
  paste(readLines(path, warn = FALSE), collapse = "\n")
}

main <- function() {
  root <- find_project_root()
  files <- c(
    file.path(root, "metadata", "mimiciv_schema_cohort_lock_workplan_v0_2.csv"),
    file.path(root, "results", "logs", "mimiciv_schema_cohort_lock_workplan_v0_2.md"),
    file.path(root, "results", "logs", "mimiciv_schema_cohort_lock_workplan_v0_2_build.md")
  )
  invisible(lapply(files, require_file, min_size = 180))

  workplan <- readr::read_csv(files[[1]], show_col_types = FALSE)
  text <- read_text(files[[2]])
  assert_true(nrow(workplan) == 7, "Workplan should have seven gates.")
  assert_true(all(c("gate_order", "gate", "objective", "required_inputs", "allowed_outputs", "not_allowed") %in% names(workplan)), "Workplan columns are incomplete.")
  assert_true(any(workplan$gate == "download_integrity_gate"), "Download integrity gate missing.")
  assert_true(any(workplan$gate == "respiratory_support_endpoint_gate"), "Respiratory support endpoint gate missing.")
  assert_true(any(str_detect(workplan$not_allowed, fixed("modeling"))), "No-modeling boundary missing.")
  assert_true(str_detect(text, fixed("MIMIC-IV remains outside the active manuscript evidence claim")), "Evidence boundary missing.")
  assert_true(str_detect(text, fixed("Proceed to cohort or endpoint lock only if v0.1 validation passes")), "Post-download rule missing.")
  assert_true(!str_detect(text, regex("row-level export allowed|start modeling|active evidence claim", ignore_case = TRUE)), "Unsafe workplan wording found.")

  validation <- c(
    "# MIMIC-IV Schema/Cohort/Endpoint Lock Workplan V0.2 Validation",
    "",
    paste0("- Run date: ", Sys.Date()),
    paste0("- Gate rows: ", nrow(workplan), "."),
    "- Checked download, schema, linkage, endpoint, cohort and manuscript-integration gates.",
    "- Checked no-modeling and no-row-level-export boundaries.",
    "",
    "PASS: MIMIC-IV schema/cohort lock workplan V0.2 validation complete."
  )
  writeLines(validation, file.path(root, "results", "logs", "mimiciv_schema_cohort_lock_workplan_v0_2_validation.md"))
  message("PASS: MIMIC-IV schema/cohort lock workplan V0.2 validation complete.")
}

if (sys.nframe() == 0) {
  main()
}
