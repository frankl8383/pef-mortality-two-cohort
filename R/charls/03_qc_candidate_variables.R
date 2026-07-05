# CHARLS P0 candidate-variable QC.
# Produces aggregate row counts, duplicate ID summaries, and missingness for candidate variables.

suppressPackageStartupMessages({
  library(dplyr)
  library(haven)
  library(readr)
  library(tidyselect)
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

extract_member <- function(archive_local_path, archive_type, inner_path) {
  td <- tempfile("charls_qc_")
  dir.create(td, recursive = TRUE)
  out <- file.path(td, basename(inner_path))
  if (archive_type == "zip") {
    utils::unzip(archive_local_path, files = inner_path, exdir = td, overwrite = TRUE)
    extracted <- file.path(td, inner_path)
    if (!file.exists(extracted)) {
      extracted <- file.path(td, basename(inner_path))
    }
    return(list(path = extracted, tmpdir = td))
  }
  if (archive_type == "rar") {
    err <- tempfile("bsdtar_err_")
    status <- system2("bsdtar", c("-xOf", archive_local_path, inner_path), stdout = out, stderr = err)
    err_text <- if (file.exists(err)) paste(readLines(err, warn = FALSE), collapse = " | ") else ""
    unlink(err)
    if (!identical(as.integer(status), 0L) || !file.exists(out) || file.info(out)$size == 0) {
      unlink(td, recursive = TRUE)
      stop("Failed to extract RAR member: ", archive_local_path, " :: ", inner_path, " ", err_text, call. = FALSE)
    }
    return(list(path = out, tmpdir = td))
  }
  stop("Unsupported archive type: ", archive_type, call. = FALSE)
}

read_selected <- function(index_row, selected_cols) {
  member <- extract_member(index_row$archive_local_path, index_row$archive_type, index_row$inner_path)
  on.exit(unlink(member$tmpdir, recursive = TRUE), add = TRUE)
  selected_cols <- unique(selected_cols[!is.na(selected_cols) & selected_cols != ""])
  haven::read_dta(member$path, col_select = tidyselect::any_of(selected_cols))
}

main <- function() {
  root <- find_project_root()
  index_path <- file.path(root, "derived", "charls_wave_file_index.csv")
  dictionary_path <- file.path(root, "metadata", "charls_variable_dictionary_draft.csv")
  key_path <- file.path(root, "metadata", "charls_key_variable_map.csv")
  if (!file.exists(index_path) || !file.exists(dictionary_path) || !file.exists(key_path)) {
    stop("Run 00_index_charls_archives.R and 02_build_variable_dictionary.R first.", call. = FALSE)
  }

  table_dir <- file.path(root, "results", "tables")
  log_dir <- file.path(root, "results", "logs")
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

  index <- read_csv(index_path, show_col_types = FALSE, progress = FALSE) %>%
    filter(preferred)
  dictionary <- read_csv(dictionary_path, show_col_types = FALSE, progress = FALSE)
  key_vars <- read_csv(key_path, show_col_types = FALSE, progress = FALSE)

  row_counts <- list()
  duplicate_reports <- list()
  missingness <- list()
  failures <- list()

  for (i in seq_len(nrow(index))) {
    row <- index[i, ]
    file_dict <- dictionary %>%
      filter(source_file == row$archive_path, inner_path == row$inner_path)
    file_keys <- key_vars %>%
      filter(source_file == row$archive_path, inner_path == row$inner_path)

    id_cols <- file_dict %>%
      filter(construct_domain == "id_linkage") %>%
      pull(variable_name)
    if (length(id_cols) == 0 && "ID" %in% file_dict$variable_name) {
      id_cols <- "ID"
    }
    first_col <- file_dict$variable_name[[1]]
    selected <- unique(c(id_cols, file_keys$variable_name, first_col))

    message("[", i, "/", nrow(index), "] QC: ", row$wave, " :: ", row$module_name)
    dat <- tryCatch(
      read_selected(row, selected),
      error = function(e) {
        failures[[length(failures) + 1]] <<- tibble::tibble(
          wave = row$wave,
          source_file = row$archive_path,
          inner_path = row$inner_path,
          error = conditionMessage(e)
        )
        NULL
      }
    )
    if (is.null(dat)) next

    row_counts[[length(row_counts) + 1]] <- tibble::tibble(
      wave = row$wave,
      source_file = row$archive_path,
      inner_path = row$inner_path,
      module_name = row$module_name,
      row_count = nrow(dat),
      selected_column_count = ncol(dat),
      total_dictionary_variable_count = nrow(file_dict),
      candidate_variable_count = nrow(file_keys)
    )

    if (length(id_cols) > 0 && id_cols[[1]] %in% names(dat)) {
      id_col <- id_cols[[1]]
      duplicate_reports[[length(duplicate_reports) + 1]] <- tibble::tibble(
        wave = row$wave,
        source_file = row$archive_path,
        inner_path = row$inner_path,
        module_name = row$module_name,
        id_variable = id_col,
        row_count = nrow(dat),
        nonmissing_id_count = sum(!is.na(dat[[id_col]])),
        unique_id_count = dplyr::n_distinct(dat[[id_col]], na.rm = TRUE),
        duplicated_id_rows = sum(duplicated(dat[[id_col]]) & !is.na(dat[[id_col]]))
      )
    } else {
      duplicate_reports[[length(duplicate_reports) + 1]] <- tibble::tibble(
        wave = row$wave,
        source_file = row$archive_path,
        inner_path = row$inner_path,
        module_name = row$module_name,
        id_variable = NA_character_,
        row_count = nrow(dat),
        nonmissing_id_count = NA_integer_,
        unique_id_count = NA_integer_,
        duplicated_id_rows = NA_integer_
      )
    }

    if (nrow(file_keys) > 0) {
      present_key_vars <- intersect(file_keys$variable_name, names(dat))
      if (length(present_key_vars) > 0) {
        miss <- lapply(present_key_vars, function(v) {
          meta <- file_keys %>% filter(variable_name == v) %>% slice(1)
          tibble::tibble(
            wave = row$wave,
            source_file = row$archive_path,
            inner_path = row$inner_path,
            module_name = row$module_name,
            variable_name = v,
            variable_label = meta$variable_label,
            construct_domain = meta$construct_domain,
            candidate_harmonized_name = meta$candidate_harmonized_name,
            confidence = meta$confidence,
            row_count = nrow(dat),
            missing_count = sum(is.na(dat[[v]])),
            nonmissing_count = sum(!is.na(dat[[v]])),
            missing_pct = ifelse(nrow(dat) > 0, sum(is.na(dat[[v]])) / nrow(dat), NA_real_)
          )
        })
        missingness[[length(missingness) + 1]] <- bind_rows(miss)
      }
    }
  }

  row_counts_df <- bind_rows(row_counts)
  duplicate_df <- bind_rows(duplicate_reports)
  missingness_df <- bind_rows(missingness)
  failures_df <- bind_rows(failures)

  write_csv(row_counts_df, file.path(table_dir, "charls_file_row_counts.csv"))
  write_csv(duplicate_df, file.path(table_dir, "charls_duplicate_id_report.csv"))
  write_csv(missingness_df, file.path(table_dir, "charls_candidate_missingness_table.csv"))
  write_csv(
    missingness_df %>% filter(construct_domain == "pef_breathing_test") %>% arrange(wave, source_file, variable_name),
    file.path(table_dir, "charls_pef_candidate_availability.csv")
  )
  write_csv(
    missingness_df %>% filter(construct_domain %in% c("chronic_lung_disease", "asthma")) %>% arrange(wave, source_file, variable_name),
    file.path(table_dir, "charls_lung_disease_candidate_availability.csv")
  )
  failure_path <- file.path(log_dir, "charls_qc_read_failures.csv")
  if (nrow(failures_df) > 0) {
    write_csv(failures_df, failure_path)
  } else if (file.exists(failure_path)) {
    unlink(failure_path)
  }

  pef_summary <- missingness_df %>%
    filter(construct_domain == "pef_breathing_test") %>%
    count(wave, module_name, name = "pef_candidate_variables")
  lung_summary <- missingness_df %>%
    filter(construct_domain %in% c("chronic_lung_disease", "asthma")) %>%
    count(wave, module_name, construct_domain, name = "lung_candidate_variables")

  log <- c(
    "# CHARLS P0 QC Log",
    "",
    paste0("- Preferred DTA files scanned: ", nrow(index)),
    paste0("- Successful QC reads: ", nrow(row_counts_df)),
    paste0("- Failed QC reads: ", nrow(failures_df)),
    paste0("- Candidate missingness rows: ", nrow(missingness_df)),
    "",
    "## Row Count Summary",
    "",
    paste(capture.output(print(row_counts_df %>% count(wave, name = "preferred_files") %>% left_join(row_counts_df %>% group_by(wave) %>% summarise(total_rows_across_files = sum(row_count), .groups = "drop"), by = "wave"), n = Inf)), collapse = "\n"),
    "",
    "## PEF Candidate Summary",
    "",
    paste(capture.output(print(pef_summary, n = Inf)), collapse = "\n"),
    "",
    "## Lung Disease Candidate Summary",
    "",
    paste(capture.output(print(lung_summary, n = Inf)), collapse = "\n"),
    "",
    "## Interpretation Boundary",
    "",
    "This is an automated P0 QC pass. Variable meanings, coding direction, and phenotype definitions still require manual codebook confirmation before modeling."
  )
  writeLines(log, file.path(log_dir, "charls_qc_log.md"))
  message("Wrote CHARLS P0 QC outputs.")
}

if (sys.nframe() == 0) {
  main()
}
