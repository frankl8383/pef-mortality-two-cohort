# CHARLS codebook audit skeleton.
# Run this before any CHARLS cleaning or analysis.

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || is.na(x)) y else x
}

find_project_root <- function() {
  script_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", script_args, value = TRUE)
  if (length(file_arg) > 0) {
    script_path <- sub("^--file=", "", file_arg[[1]])
    return(normalizePath(file.path(dirname(script_path), "..", ".."), mustWork = FALSE))
  }
  normalizePath(getwd(), mustWork = FALSE)
}

root <- find_project_root()
dict_path <- file.path(root, "data_dict", "charls_variable_dictionary_source.csv")
manifest_path <- file.path(root, "data_manifest", "charls_manifest.csv")
report_path <- file.path(root, "results", "logs", "charls_codebook_audit.md")

read_required_csv <- function(path) {
  if (!file.exists(path)) {
    stop("Missing required file: ", path, call. = FALSE)
  }
  readr::read_csv(path, show_col_types = FALSE, progress = FALSE)
}

dictionary <- read_required_csv(dict_path)
manifest <- read_required_csv(manifest_path)

todo_dictionary <- dictionary %>%
  filter(if_any(everything(), ~ .x == "TODO_CODEBOOK_CHECK"))

todo_manifest <- manifest %>%
  filter(if_any(everything(), ~ .x == "TODO_CODEBOOK_CHECK"))

codebook_dir <- Sys.getenv("CHARLS_CODEBOOK_DIR", unset = "")
raw_dir <- Sys.getenv("CHARLS_RAW_DIR", unset = "")
derived_dir <- Sys.getenv("CHARLS_DERIVED_DIR", unset = "")

dir_status <- tibble::tibble(
  env_var = c("CHARLS_CODEBOOK_DIR", "CHARLS_RAW_DIR", "CHARLS_DERIVED_DIR"),
  value = c(codebook_dir, raw_dir, derived_dir),
  exists = c(dir.exists(codebook_dir), dir.exists(raw_dir), dir.exists(derived_dir))
)

report <- c(
  "# CHARLS Codebook Audit",
  "",
  paste0("- Dictionary rows: ", nrow(dictionary)),
  paste0("- Dictionary rows with TODO_CODEBOOK_CHECK: ", nrow(todo_dictionary)),
  paste0("- Manifest rows: ", nrow(manifest)),
  paste0("- Manifest rows with TODO_CODEBOOK_CHECK: ", nrow(todo_manifest)),
  "",
  "## Local Directory Status",
  "",
  paste(capture.output(print(dir_status, n = Inf)), collapse = "\n"),
  "",
  "## Interpretation",
  "",
  if (nrow(todo_dictionary) > 0 || nrow(todo_manifest) > 0) {
    "CHARLS codebook mapping is not complete. Do not run cleaning or models."
  } else {
    "No TODO_CODEBOOK_CHECK values remain in the current manifest and dictionary. Proceed to source-file existence checks before cleaning."
  }
)

dir.create(dirname(report_path), recursive = TRUE, showWarnings = FALSE)
writeLines(report, report_path)
message("Wrote audit report: ", report_path)

if (nrow(todo_dictionary) > 0 || nrow(todo_manifest) > 0) {
  stop("CHARLS dictionary or manifest still contains TODO_CODEBOOK_CHECK.", call. = FALSE)
}
