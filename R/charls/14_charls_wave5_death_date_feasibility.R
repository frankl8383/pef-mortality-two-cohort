# CHARLS Wave 5 death-date feasibility audit.
# Writes aggregate feasibility summaries only. No row-level data are exported.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
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

require_input <- function(path) {
  if (!file.exists(path)) stop("Missing required input: ", path, call. = FALSE)
  path
}

safe_date <- function(year, month, day) {
  year <- suppressWarnings(as.integer(year))
  month <- suppressWarnings(as.integer(month))
  day <- suppressWarnings(as.integer(day))
  ok <- !is.na(year) & !is.na(month) & !is.na(day) &
    year >= 1900 & year <= 2030 & month >= 1 & month <= 12 & day >= 1 & day <= 31
  out <- rep(as.Date(NA), length(year))
  out[ok] <- suppressWarnings(as.Date(sprintf("%04d-%02d-%02d", year[ok], month[ok], day[ok])))
  out
}

pct <- function(num, den) {
  ifelse(den > 0, 100 * num / den, NA_real_)
}

main <- function() {
  root <- find_project_root()
  source(file.path(root, "R", "charls", "13_charls_v0_6_wave5_main_models.R"))

  input_rds <- require_input(file.path(root, "derived_sensitive", "charls", "charls_core_harmonized_provisional.rds"))
  zip_path <- require_input("${CHARLS_RAW_ROOT}/2020年全国追踪调查/数据下载/CHARLS2020r.zip")
  table_dir <- file.path(root, "results", "tables")
  log_dir <- file.path(root, "results", "logs")
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

  core <- readRDS(input_rds)
  wave5 <- load_wave5(zip_path)$wave5
  core_w5 <- core %>% left_join(wave5, by = "participant_id")
  base <- prepare_base_w5(core_w5)
  death <- prepare_death_time_w5(base) %>%
    mutate(
      wave5_death_date = safe_date(.data$death_year_w5, .data$death_month_w5, .data$death_day_w5),
      has_any_death_date_part = !is.na(.data$death_year_w5) | !is.na(.data$death_month_w5) | !is.na(.data$death_day_w5),
      has_valid_death_date = !is.na(.data$wave5_death_date),
      plausible_wave5_death_date = .data$has_valid_death_date &
        .data$wave5_death_date >= as.Date("2011-01-01") &
        .data$wave5_death_date <= as.Date("2021-12-31"),
      wave5_death_event = .data$death_wave_w1_w5 == 5
    )

  total_n <- nrow(death)
  total_deaths <- sum(death$death_event == 1, na.rm = TRUE)
  wave5_deaths <- sum(death$wave5_death_event, na.rm = TRUE)
  valid_exact_wave5 <- sum(death$wave5_death_event & death$plausible_wave5_death_date, na.rm = TRUE)
  any_part_wave5 <- sum(death$wave5_death_event & death$has_any_death_date_part, na.rm = TRUE)
  all_valid_dates <- sum(death$plausible_wave5_death_date, na.rm = TRUE)
  exact_status <- if (wave5_deaths > 0 && valid_exact_wave5 / wave5_deaths >= 0.80) {
    "feasible_for_primary_exact_date_sensitivity"
  } else if (wave5_deaths > 0 && valid_exact_wave5 / wave5_deaths >= 0.50) {
    "borderline_partial_sensitivity_only"
  } else {
    "not_feasible_for_primary_exact_date_sensitivity"
  }

  feasibility <- tibble(
    metric = c(
      "death_model_rows",
      "death_events_w2_w5",
      "wave5_death_events",
      "wave5_deaths_with_any_date_component",
      "wave5_deaths_with_valid_plausible_exact_date",
      "all_analysis_rows_with_valid_plausible_death_date",
      "exact_date_coverage_among_wave5_deaths_percent",
      "exact_date_coverage_among_all_deaths_percent",
      "feasibility_status"
    ),
    value = c(
      as.character(total_n),
      as.character(total_deaths),
      as.character(wave5_deaths),
      as.character(any_part_wave5),
      as.character(valid_exact_wave5),
      as.character(all_valid_dates),
      sprintf("%.1f", pct(valid_exact_wave5, wave5_deaths)),
      sprintf("%.1f", pct(valid_exact_wave5, total_deaths)),
      exact_status
    )
  )

  date_distribution <- death %>%
    filter(.data$plausible_wave5_death_date) %>%
    mutate(death_year = as.integer(format(.data$wave5_death_date, "%Y"))) %>%
    count(death_year, name = "n") %>%
    arrange(.data$death_year)
  if (nrow(date_distribution) == 0) {
    date_distribution <- tibble(death_year = integer(), n = integer())
  }

  year_month_distribution <- death %>%
    filter(.data$plausible_wave5_death_date) %>%
    mutate(
      death_year = as.integer(format(.data$wave5_death_date, "%Y")),
      death_month = as.integer(format(.data$wave5_death_date, "%m"))
    ) %>%
    count(death_year, death_month, name = "n") %>%
    arrange(.data$death_year, .data$death_month)
  if (nrow(year_month_distribution) == 0) {
    year_month_distribution <- tibble(death_year = integer(), death_month = integer(), n = integer())
  }

  readr::write_csv(feasibility, file.path(table_dir, "charls_wave5_death_date_feasibility_v0_4.csv"))
  readr::write_csv(date_distribution, file.path(table_dir, "charls_wave5_death_date_year_distribution_v0_4.csv"))
  readr::write_csv(year_month_distribution, file.path(table_dir, "charls_wave5_death_date_year_month_distribution_v0_4.csv"))

  lines <- c(
    "# CHARLS Wave 5 Death-Date Feasibility V0.4",
    "",
    paste0("- Run date: ", Sys.Date()),
    "- Scope: aggregate audit of whether CHARLS 2020 Exit Module death dates can support an exact-date mortality sensitivity.",
    paste0("- Death model rows: ", total_n, "."),
    paste0("- Total death events through Wave 5 / 2020: ", total_deaths, "."),
    paste0("- Wave 5 death events in the analysis cohort: ", wave5_deaths, "."),
    paste0("- Wave 5 deaths with plausible exact date: ", valid_exact_wave5, " (", sprintf("%.1f", pct(valid_exact_wave5, wave5_deaths)), "% of Wave 5 deaths)."),
    paste0("- Feasibility status: ", exact_status, "."),
    "",
    "## Interpretation",
    "",
    if (exact_status == "not_feasible_for_primary_exact_date_sensitivity") {
      "- Exact-date mortality sensitivity is not promoted to the main manuscript because date coverage among Wave 5 death events is too low. The current manuscript should retain wave-interval timing and mention exact-date mortality as a future sensitivity requiring fuller date capture."
    } else {
      "- Exact-date coverage appears sufficient for a date-informed sensitivity; a subsequent model script can be added if needed."
    }
  )
  writeLines(lines, file.path(log_dir, "charls_wave5_death_date_feasibility_v0_4.md"))
  message("Wrote CHARLS Wave 5 death-date feasibility outputs.")
}

if (sys.nframe() == 0) {
  main()
}
