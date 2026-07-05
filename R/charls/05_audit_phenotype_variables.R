# CHARLS P1 phenotype-variable distribution audit.
# Reads only selected columns and writes aggregate summaries, never row-level outputs.

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
  td <- tempfile("charls_pheno_audit_")
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

num_summary <- function(x) {
  z <- suppressWarnings(as.numeric(x))
  z <- z[!is.na(z)]
  if (length(z) == 0) {
    return(tibble::tibble(
      min = NA_real_, p01 = NA_real_, p05 = NA_real_, median = NA_real_,
      mean = NA_real_, p95 = NA_real_, p99 = NA_real_, max = NA_real_
    ))
  }
  qs <- as.numeric(stats::quantile(z, probs = c(0.01, 0.05, 0.5, 0.95, 0.99), na.rm = TRUE, names = FALSE))
  tibble::tibble(
    min = min(z, na.rm = TRUE),
    p01 = qs[[1]],
    p05 = qs[[2]],
    median = qs[[3]],
    mean = mean(z, na.rm = TRUE),
    p95 = qs[[4]],
    p99 = qs[[5]],
    max = max(z, na.rm = TRUE)
  )
}

value_count_summary <- function(x, max_levels = 40) {
  z <- x[!is.na(x)]
  if (length(z) == 0) {
    return(tibble::tibble(value = character(), n = integer(), pct_nonmissing = numeric()))
  }
  if (dplyr::n_distinct(z) > max_levels) {
    return(tibble::tibble(value = character(), n = integer(), pct_nonmissing = numeric()))
  }
  tab <- as.data.frame(table(as.character(z), useNA = "no"), stringsAsFactors = FALSE)
  names(tab) <- c("value", "n")
  tibble::as_tibble(tab) %>%
    mutate(
      n = as.integer(n),
      pct_nonmissing = n / sum(n)
    ) %>%
    arrange(desc(n), value)
}

main <- function() {
  root <- find_project_root()
  index_path <- file.path(root, "derived", "charls_wave_file_index.csv")
  draft_path <- file.path(root, "metadata", "charls_phenotype_definition_draft.csv")
  if (!file.exists(index_path) || !file.exists(draft_path)) {
    stop("Run 00_index_charls_archives.R and 04_build_phenotype_definition_draft.R first.", call. = FALSE)
  }

  table_dir <- file.path(root, "results", "tables")
  log_dir <- file.path(root, "results", "logs")
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

  index <- read_csv(index_path, show_col_types = FALSE, progress = FALSE) %>%
    filter(preferred)
  draft <- read_csv(draft_path, show_col_types = FALSE, progress = FALSE) %>%
    filter(priority_tier %in% c("core_review", "standard_review", "supporting_review"))

  missingness <- list()
  numeric_summaries <- list()
  value_counts <- list()
  failures <- list()

  for (i in seq_len(nrow(index))) {
    row <- index[i, ]
    file_defs <- draft %>%
      filter(source_file == row$archive_path, inner_path == row$inner_path)
    if (nrow(file_defs) == 0) next

    selected <- unique(c("ID", file_defs$variable_name))
    message("[", i, "/", nrow(index), "] phenotype audit: ", row$wave, " :: ", row$module_name)
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

    present_vars <- intersect(file_defs$variable_name, names(dat))
    for (v in present_vars) {
      meta <- file_defs %>% filter(variable_name == v) %>% slice(1)
      x <- dat[[v]]
      row_count <- nrow(dat)
      missing_count <- sum(is.na(x))
      nonmissing_count <- sum(!is.na(x))
      distinct_nonmissing <- dplyr::n_distinct(x, na.rm = TRUE)
      base <- tibble::tibble(
        phenotype_id = meta$phenotype_id,
        priority_tier = meta$priority_tier,
        analysis_role = meta$analysis_role,
        proposed_concept = meta$proposed_concept,
        wave = row$wave,
        module_name = row$module_name,
        source_file = row$archive_path,
        inner_path = row$inner_path,
        variable_name = v,
        variable_label = meta$variable_label,
        value_labels = meta$value_labels,
        row_count = row_count,
        missing_count = missing_count,
        nonmissing_count = nonmissing_count,
        missing_pct = ifelse(row_count > 0, missing_count / row_count, NA_real_),
        distinct_nonmissing = distinct_nonmissing
      )
      missingness[[length(missingness) + 1]] <- base

      if (is.numeric(x) || inherits(x, "haven_labelled")) {
        numeric_summaries[[length(numeric_summaries) + 1]] <- bind_cols(base, num_summary(x))
      }

      vc <- value_count_summary(x)
      if (nrow(vc) > 0) {
        value_counts[[length(value_counts) + 1]] <- bind_cols(
          base %>% select(-row_count, -missing_count, -nonmissing_count, -missing_pct, -distinct_nonmissing),
          vc
        )
      }
    }
  }

  missingness_df <- bind_rows(missingness)
  numeric_df <- bind_rows(numeric_summaries)
  value_counts_df <- bind_rows(value_counts)
  failures_df <- bind_rows(failures)

  write_csv(missingness_df, file.path(table_dir, "charls_phenotype_variable_missingness.csv"))
  write_csv(numeric_df, file.path(table_dir, "charls_phenotype_variable_numeric_summary.csv"))
  write_csv(value_counts_df, file.path(table_dir, "charls_phenotype_variable_value_counts.csv"))

  failure_path <- file.path(log_dir, "charls_phenotype_audit_failures.csv")
  if (nrow(failures_df) > 0) {
    write_csv(failures_df, failure_path)
  } else if (file.exists(failure_path)) {
    unlink(failure_path)
  }

  pef_numeric <- numeric_df %>%
    filter(proposed_concept %in% c("pef_trial_value", "pef_harmonized_value_or_quality")) %>%
    select(wave, module_name, proposed_concept, variable_name, row_count, nonmissing_count, missing_pct, min, p01, median, p99, max)
  lung_counts <- value_counts_df %>%
    filter(proposed_concept %in% c("chronic_lung_disease_status", "asthma_status")) %>%
    select(wave, module_name, proposed_concept, variable_name, value, n, pct_nonmissing, value_labels)

  write_csv(pef_numeric, file.path(table_dir, "charls_pef_numeric_distribution_audit.csv"))
  write_csv(lung_counts, file.path(table_dir, "charls_lung_asthma_value_count_audit.csv"))

  role_counts <- missingness_df %>%
    count(priority_tier, analysis_role, proposed_concept, name = "audited_variables") %>%
    arrange(priority_tier, analysis_role, desc(audited_variables))

  log <- c(
    "# CHARLS P1 Phenotype Variable Audit Log",
    "",
    paste0("- Phenotype definitions considered: ", nrow(draft)),
    paste0("- Variables audited: ", nrow(missingness_df)),
    paste0("- Numeric summary rows: ", nrow(numeric_df)),
    paste0("- Value-count rows: ", nrow(value_counts_df)),
    paste0("- Failed DTA reads: ", nrow(failures_df)),
    "",
    "## Audited Variable Counts",
    "",
    paste(capture.output(print(role_counts, n = Inf)), collapse = "\n"),
    "",
    "## Boundary",
    "",
    "This audit writes only aggregate summaries. It is intended to support codebook locking and does not create an analysis-ready cohort."
  )
  writeLines(log, file.path(log_dir, "charls_phenotype_variable_audit_log.md"))
  message("Wrote CHARLS P1 phenotype-variable audit outputs.")
}

if (sys.nframe() == 0) {
  main()
}
