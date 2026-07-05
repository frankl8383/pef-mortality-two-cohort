# CHARLS cleaning skeleton.
# This file defines the cleaning boundary but refuses to run until codebook mapping is complete.

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
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

read_dictionary <- function(root = find_project_root()) {
  path <- file.path(root, "data_dict", "charls_variable_dictionary_source.csv")
  if (!file.exists(path)) {
    stop("Missing CHARLS dictionary: ", path, call. = FALSE)
  }
  readr::read_csv(path, show_col_types = FALSE, progress = FALSE)
}

assert_dictionary_verified <- function(dictionary) {
  unresolved <- dictionary %>%
    filter(if_any(c(variable_name, source_table_or_file, codebook_status), ~ .x == "TODO_CODEBOOK_CHECK"))

  if (nrow(unresolved) > 0) {
    stop(
      "CHARLS dictionary has unresolved TODO_CODEBOOK_CHECK entries. ",
      "Update codebook mappings before cleaning.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

standardize_charls_wave <- function(raw_wave, dictionary, wave_label) {
  assert_dictionary_verified(dictionary)

  required_columns <- dictionary %>%
    filter(wave_or_cycle %in% c("all", wave_label, "baseline_and_followup")) %>%
    pull(variable_name) %>%
    unique()

  missing_columns <- setdiff(required_columns, names(raw_wave))
  if (length(missing_columns) > 0) {
    stop(
      "Raw CHARLS wave is missing verified columns: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }

  stop(
    "Implement wave-specific renaming and harmonization only after official codebook review. ",
    "Do not add source variables by memory.",
    call. = FALSE
  )
}

main <- function() {
  root <- find_project_root()
  dictionary <- read_dictionary(root)
  assert_dictionary_verified(dictionary)

  raw_dir <- Sys.getenv("CHARLS_RAW_DIR", unset = "")
  derived_dir <- Sys.getenv("CHARLS_DERIVED_DIR", unset = "")

  if (!dir.exists(raw_dir)) {
    stop("Set CHARLS_RAW_DIR to a local restricted-data folder.", call. = FALSE)
  }
  if (!dir.exists(derived_dir)) {
    stop("Set CHARLS_DERIVED_DIR to a local-only derived-data folder.", call. = FALSE)
  }

  stop("Cleaning implementation is intentionally blocked until codebook mappings are verified.", call. = FALSE)
}

if (sys.nframe() == 0) {
  main()
}

