# Validate NHANES V30 linked mortality ETL outputs.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
})

version_id <- "v30_0"

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

check_file <- function(path, label, required = TRUE) {
  exists <- file.exists(path)
  size <- if (exists) file.info(path)$size else NA_real_
  passed <- exists && isTRUE(size > 0)
  tibble(
    check = paste0("file_exists_nonempty:", label),
    passed = ifelse(required, passed, TRUE),
    observed = ifelse(exists, as.character(size), "missing"),
    expected = ifelse(required, ">0 bytes", "optional")
  )
}

check_scalar <- function(check, observed, expected, passed) {
  tibble(
    check = check,
    passed = passed,
    observed = as.character(observed),
    expected = expected
  )
}

scan_for_absolute_user_paths <- function(paths) {
  hits <- character()
  for (p in paths[file.exists(paths)]) {
    txt <- paste(readLines(p, warn = FALSE), collapse = "\n")
    if (grepl("${HOME}", txt, fixed = TRUE)) {
      hits <- c(hits, basename(p))
    }
  }
  hits
}

root <- find_project_root()
derived_dir <- file.path(root, "derived_sensitive", "nhanes")
results_dir <- file.path(root, "results", "tables")
logs_dir <- file.path(root, "results", "logs")
metadata_dir <- file.path(root, "metadata")
dir_create(results_dir)
dir_create(logs_dir)

row_rds <- file.path(derived_dir, paste0("nhanes_mortality_analysis_ready_", version_id, ".rds"))
row_rds_latest <- file.path(derived_dir, "nhanes_mortality_analysis_ready.rds")
row_parquet <- file.path(derived_dir, paste0("nhanes_mortality_analysis_ready_", version_id, ".parquet"))
row_parquet_latest <- file.path(derived_dir, "nhanes_mortality_analysis_ready.parquet")
flow_path <- file.path(results_dir, paste0("nhanes_mortality_flow_", version_id, ".csv"))
counts_path <- file.path(results_dir, paste0("nhanes_mortality_counts_", version_id, ".csv"))
manifest_path <- file.path(results_dir, paste0("nhanes_mortality_download_manifest_", version_id, ".csv"))
codebook_path <- file.path(metadata_dir, paste0("nhanes_mortality_codebook_", version_id, ".md"))
etl_log_path <- file.path(logs_dir, paste0("nhanes_mortality_etl_", version_id, ".md"))
validation_csv_path <- file.path(results_dir, paste0("nhanes_mortality_validation_", version_id, ".csv"))
validation_log_path <- file.path(logs_dir, paste0("nhanes_mortality_validation_", version_id, ".md"))

checks <- bind_rows(
  check_file(row_rds, "row_rds"),
  check_file(row_rds_latest, "row_rds_latest"),
  check_file(row_parquet, "row_parquet", required = FALSE),
  check_file(row_parquet_latest, "row_parquet_latest", required = FALSE),
  check_file(flow_path, "flow_csv"),
  check_file(counts_path, "counts_csv"),
  check_file(manifest_path, "download_manifest_csv"),
  check_file(codebook_path, "codebook_md"),
  check_file(etl_log_path, "etl_log_md")
)

if (!file.exists(row_rds)) {
  stop("Cannot validate without row-level RDS output.", call. = FALSE)
}

analysis_ready <- readRDS(row_rds)
flow <- readr::read_csv(flow_path, show_col_types = FALSE)
counts <- readr::read_csv(counts_path, show_col_types = FALSE)
manifest <- readr::read_csv(manifest_path, show_col_types = FALSE)
primary <- analysis_ready %>% filter(mortality_analysis_primary == 1L)
adult45 <- analysis_ready %>% filter(!is.na(adult_45plus) & (adult_45plus == TRUE | adult_45plus == 1))
matched_n <- sum(analysis_ready$mortality_merge_status == "matched", na.rm = TRUE)
primary_deaths <- sum(primary$all_cause_death == 1L, na.rm = TRUE)
primary_clrd_deaths <- sum(primary$clrd_death == 1L, na.rm = TRUE)
positive_followup <- all(primary$followup_years_exam > 0, na.rm = TRUE)
cycle_levels <- paste(sort(unique(primary$cycle_label)), collapse = ",")
quartile_levels <- paste(sort(unique(primary$rv_quartile_num)), collapse = ",")
absolute_path_hits <- scan_for_absolute_user_paths(c(flow_path, counts_path, manifest_path, codebook_path, etl_log_path))

checks <- bind_rows(
  checks,
  check_scalar("expected_total_rows", nrow(analysis_ready), "30442", nrow(analysis_ready) == 30442),
  check_scalar("mortality_merge_complete", matched_n, "30442", matched_n == nrow(analysis_ready)),
  check_scalar("adult45_rows_plausible", nrow(adult45), ">=10000", nrow(adult45) >= 10000),
  check_scalar("primary_analysis_rows_plausible", nrow(primary), ">=6000", nrow(primary) >= 6000),
  check_scalar("primary_all_cause_deaths_plausible", primary_deaths, ">=500", primary_deaths >= 500),
  check_scalar("primary_clrd_deaths_exploratory_count", primary_clrd_deaths, ">=20 but treat as exploratory", primary_clrd_deaths >= 20),
  check_scalar("positive_followup_in_primary", positive_followup, "TRUE", isTRUE(positive_followup)),
  check_scalar("cycle_coverage", cycle_levels, "E,F,G", identical(cycle_levels, "E,F,G")),
  check_scalar("quartile_coverage", quartile_levels, "1,2,3,4", identical(quartile_levels, "1,2,3,4")),
  check_scalar("aggregate_tables_nonempty", paste(nrow(flow), nrow(counts), sep = "/"), "flow and counts rows > 0", nrow(flow) > 0 && nrow(counts) > 0),
  check_scalar("source_manifest_complete", paste(manifest$download_status, collapse = ","), "all present/downloaded", all(manifest$download_status %in% c("present", "downloaded"))),
  check_scalar("no_absolute_user_paths_in_public_outputs", paste(absolute_path_hits, collapse = ","), "none", length(absolute_path_hits) == 0)
)

all_passed <- all(checks$passed)
checks <- checks %>%
  mutate(status = ifelse(passed, "pass", "fail")) %>%
  select(check, status, observed, expected)
readr::write_csv(checks, validation_csv_path)

if (all_passed) {
  ticket_path <- file.path(results_dir, "v30_high_impact_ticket_audit.csv")
  if (file.exists(ticket_path)) {
    ticket_audit <- readr::read_csv(ticket_path, show_col_types = FALSE)
    ticket_audit <- ticket_audit %>%
      mutate(status = ifelse(ticket == "T1_NHANES_linked_mortality_ETL", "completed", status))
    readr::write_csv(ticket_audit, ticket_path)
  }
}

log_lines <- c(
  "# NHANES V30 Linked Mortality Validation",
  "",
  paste0("Version: ", version_id),
  "",
  paste0("Overall status: ", ifelse(all_passed, "PASS", "FAIL")),
  "",
  "## Key Counts",
  "",
  paste0("- Total linked rows: ", nrow(analysis_ready)),
  paste0("- Matched rows: ", matched_n),
  paste0("- Adult 45+ rows: ", nrow(adult45)),
  paste0("- Primary mortality analysis rows: ", nrow(primary)),
  paste0("- Primary all-cause deaths: ", primary_deaths),
  paste0("- Primary chronic lower respiratory disease deaths: ", primary_clrd_deaths, " (exploratory)"),
  "",
  "## Checks",
  "",
  paste0("| ", paste(names(checks), collapse = " | "), " |"),
  paste0("| ", paste(rep("---", ncol(checks)), collapse = " | "), " |"),
  apply(checks, 1, function(x) paste0("| ", paste(x, collapse = " | "), " |"))
)
writeLines(log_lines, validation_log_path)

message("NHANES linked mortality validation status: ", ifelse(all_passed, "PASS", "FAIL"))
if (!all_passed) {
  failed <- checks %>% filter(status == "fail")
  print(failed)
  stop("Validation failed.", call. = FALSE)
}
