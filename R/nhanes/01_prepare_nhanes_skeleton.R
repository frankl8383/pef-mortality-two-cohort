# NHANES preparation skeleton.
# This script checks mapping and local raw/derived directories, but does not yet
# implement row-level NHANES cleaning.

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

assert_no_todos <- function(data, label) {
  todo_count <- sum(data == "TODO_CODEBOOK_CHECK", na.rm = TRUE)
  if (todo_count > 0) {
    stop(label, " contains ", todo_count, " TODO_CODEBOOK_CHECK values.", call. = FALSE)
  }
  invisible(TRUE)
}

assert_nhanes_mapping_verified <- function(root = find_project_root()) {
  manifest <- read_required_csv(file.path(root, "data_manifest", "nhanes_manifest.csv"))
  dictionary <- read_required_csv(file.path(root, "data_dict", "nhanes_variable_dictionary_source.csv"))
  design_map <- read_required_csv(file.path(root, "data_dict", "nhanes_survey_design_map.csv"))

  assert_no_todos(manifest, "NHANES manifest")
  assert_no_todos(dictionary, "NHANES dictionary")
  assert_no_todos(design_map, "NHANES survey design map")
  invisible(TRUE)
}

main <- function() {
  root <- find_project_root()
  assert_nhanes_mapping_verified(root)

  raw_dir <- Sys.getenv("NHANES_RAW_DIR", unset = "")
  derived_dir <- Sys.getenv("NHANES_DERIVED_DIR", unset = "")

  if (!dir.exists(raw_dir)) {
    stop("Set NHANES_RAW_DIR to a local folder containing raw XPT files.", call. = FALSE)
  }
  if (!dir.exists(derived_dir)) {
    stop("Set NHANES_DERIVED_DIR to a local-only derived-data folder.", call. = FALSE)
  }

  stop(
    "NHANES codebook mapping is verified, but this skeleton does not yet implement ",
    "row-level cleaning. Implement or run the V0.2 cleaner after downloading local ",
    "raw XPT files using metadata/nhanes_replication_file_manifest_v0_1.csv.",
    call. = FALSE
  )
}

if (sys.nframe() == 0) {
  main()
}
