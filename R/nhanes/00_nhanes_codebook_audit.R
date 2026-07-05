# NHANES codebook audit skeleton.
# Run this before any NHANES cleaning or survey-weighted modeling.

find_project_root <- function() {
  script_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", script_args, value = TRUE)
  if (length(file_arg) > 0) {
    script_path <- sub("^--file=", "", file_arg[[1]])
    return(normalizePath(file.path(dirname(script_path), "..", ".."), mustWork = FALSE))
  }
  normalizePath(getwd(), mustWork = FALSE)
}

read_required_csv <- function(path) {
  if (!file.exists(path)) {
    stop("Missing required file: ", path, call. = FALSE)
  }
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

count_todos <- function(data) {
  sum(data == "TODO_CODEBOOK_CHECK", na.rm = TRUE)
}

root <- find_project_root()
paths <- c(
  manifest = file.path(root, "data_manifest", "nhanes_manifest.csv"),
  dictionary = file.path(root, "data_dict", "nhanes_variable_dictionary_source.csv"),
  survey_design = file.path(root, "data_dict", "nhanes_survey_design_map.csv")
)

tables <- lapply(paths, read_required_csv)
todo_counts <- vapply(tables, count_todos, integer(1))

dir_status <- data.frame(
  env_var = c("NHANES_CODEBOOK_DIR", "NHANES_RAW_DIR", "NHANES_DERIVED_DIR", "NHANES_MORTALITY_DIR"),
  value = c(
    Sys.getenv("NHANES_CODEBOOK_DIR", unset = ""),
    Sys.getenv("NHANES_RAW_DIR", unset = ""),
    Sys.getenv("NHANES_DERIVED_DIR", unset = ""),
    Sys.getenv("NHANES_MORTALITY_DIR", unset = "")
  ),
  stringsAsFactors = FALSE
)
dir_status$exists <- dir.exists(dir_status$value)

report_path <- file.path(root, "results", "logs", "nhanes_codebook_audit.md")
dir.create(dirname(report_path), recursive = TRUE, showWarnings = FALSE)

report <- c(
  "# NHANES Codebook Audit",
  "",
  paste0("- Manifest rows: ", nrow(tables$manifest)),
  paste0("- Dictionary rows: ", nrow(tables$dictionary)),
  paste0("- Survey design map rows: ", nrow(tables$survey_design)),
  "",
  "## TODO_CODEBOOK_CHECK Counts",
  "",
  paste(capture.output(print(data.frame(table = names(todo_counts), todo_count = as.integer(todo_counts)))), collapse = "\n"),
  "",
  "## Local Directory Status",
  "",
  paste(capture.output(print(dir_status, row.names = FALSE)), collapse = "\n"),
  "",
  "## Interpretation",
  "",
  if (sum(todo_counts) > 0) {
    "NHANES codebook mapping is not complete. Do not run cleaning or survey-weighted models."
  } else {
    "No TODO_CODEBOOK_CHECK values remain in the current NHANES manifest, dictionary, and survey design map."
  }
)

writeLines(report, report_path)
message("Wrote audit report: ", report_path)

if (sum(todo_counts) > 0) {
  stop("NHANES manifest, dictionary, or survey design map still contains TODO_CODEBOOK_CHECK.", call. = FALSE)
}

