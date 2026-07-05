# CHARLS Wave 5 endpoint feasibility audit.
# Reads local 2020 CHARLS modules only to produce aggregate feasibility outputs.

suppressPackageStartupMessages({
  library(dplyr)
  library(haven)
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

to_num <- function(x) suppressWarnings(as.numeric(haven::zap_labels(haven::zap_missing(x))))
to_chr <- function(x) as.character(haven::zap_labels(haven::zap_missing(x)))

read_zip_member <- function(zip_path, member, cols = NULL) {
  td <- tempfile("charls_wave5_")
  dir.create(td, recursive = TRUE)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)
  utils::unzip(zip_path, files = member, exdir = td, overwrite = TRUE)
  path <- file.path(td, member)
  if (is.null(cols)) {
    haven::read_dta(path)
  } else {
    haven::read_dta(path, col_select = tidyselect::any_of(cols))
  }
}

label_of <- function(dat, var) {
  if (!var %in% names(dat)) return(NA_character_)
  lab <- attr(dat[[var]], "label")
  if (is.null(lab)) NA_character_ else as.character(lab)
}

count_values <- function(dat, var) {
  if (!var %in% names(dat)) {
    return(tibble(variable = var, value = NA_real_, n = NA_integer_))
  }
  x <- to_num(dat[[var]])
  as.data.frame(table(x, useNA = "ifany"), stringsAsFactors = FALSE) %>%
    as_tibble() %>%
    transmute(
      variable = var,
      value = suppressWarnings(as.numeric(as.character(.data$x))),
      n = as.integer(.data$Freq)
    )
}

candidate_row <- function(module, variable, label, role, status, note) {
  tibble(
    module = module,
    variable = variable,
    label = label,
    endpoint_role = role,
    feasibility_status = status,
    note = note
  )
}

main <- function() {
  root <- find_project_root()
  zip_path <- "${CHARLS_RAW_ROOT}/2020年全国追踪调查/数据下载/CHARLS2020r.zip"
  if (!file.exists(zip_path)) {
    stop("Missing local CHARLS 2020 zip: ", zip_path, call. = FALSE)
  }

  metadata_dir <- file.path(root, "metadata")
  table_dir <- file.path(root, "results", "tables")
  log_dir <- file.path(root, "results", "logs")
  dir.create(metadata_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

  core_path <- file.path(root, "derived_sensitive", "charls", "charls_core_harmonized_provisional.rds")
  core <- readRDS(core_path)

  health <- read_zip_member(
    zip_path,
    "Health_Status_and_Functioning.dta",
    c("ID", "householdID", "communityID", "da002_5_", "da003_5_", "da004_5_", "xchrodistype_5_", "zdisease_5_")
  )
  sample_info <- read_zip_member(zip_path, "Sample_Infor.dta", c("ID", "householdID", "communityID", "died"))
  exit_mod <- read_zip_member(
    zip_path,
    "Exit_Module.dta",
    c("ID", "householdID", "communityID", "exb001_1", "exb001_2", "exb001_3", "xezdisease_2_")
  )
  weights <- read_zip_member(zip_path, "Weights.dta", c("ID", "householdID", "communityID", "INDV_weight", "INDV_weight_ad2"))

  health_ids <- unique(to_chr(health$ID))
  sample_ids <- unique(to_chr(sample_info$ID))
  exit_ids <- unique(to_chr(exit_mod$ID))
  weight_ids <- unique(to_chr(weights$ID))
  core_ids <- unique(as.character(core$participant_id))

  feasibility <- bind_rows(
    candidate_row(
      "Health_Status_and_Functioning.dta",
      "da003_5_",
      label_of(health, "da003_5_"),
      "candidate_wave5_prevalent_chronic_lung_disease",
      "candidate_not_locked",
      "2020 Codebook identifies this as Diagnosed Disease [5]; Harmonized D labels map RwLUNGE to chronic lung disease through Wave 4. Requires Wave 5 raw codebook/value-label lock before modeling."
    ),
    candidate_row(
      "Health_Status_and_Functioning.dta",
      "da002_5_",
      label_of(health, "da002_5_"),
      "candidate_wave5_incidence_since_last_interview",
      "candidate_not_locked",
      "May identify doctor-diagnosed disease [5] timing compared with 2020 interview time. Needs coding-rule confirmation before use as incident Wave 5 endpoint."
    ),
    candidate_row(
      "Health_Status_and_Functioning.dta",
      "da004_5_",
      label_of(health, "da004_5_"),
      "candidate_wave5_self_known_chronic_lung_disease",
      "candidate_not_primary",
      "Self-known disease [5] is not equivalent to physician diagnosis; may be sensitivity only after lock."
    ),
    candidate_row(
      "Sample_Infor.dta",
      "died",
      label_of(sample_info, "died"),
      "candidate_wave5_death_indicator",
      "candidate_feasible",
      "Can support Wave 5 mortality extension after ID linkage and denominator checks."
    ),
    candidate_row(
      "Exit_Module.dta",
      "exb001_1/exb001_2/exb001_3",
      paste(label_of(exit_mod, "exb001_1"), label_of(exit_mod, "exb001_2"), label_of(exit_mod, "exb001_3"), sep = "; "),
      "candidate_wave5_death_date",
      "candidate_feasible",
      "Can improve death timing for decedents after exit-module linkage."
    ),
    candidate_row(
      "Exit_Module.dta",
      "xezdisease_2_",
      label_of(exit_mod, "xezdisease_2_"),
      "candidate_exit_chronic_lung_disease",
      "candidate_context_only",
      "Exit-module chronic lung disease is decedent-proxy context, not a substitute for main interview incident CLD."
    ),
    candidate_row(
      "Weights.dta",
      "INDV_weight/INDV_weight_ad2",
      paste(label_of(weights, "INDV_weight"), label_of(weights, "INDV_weight_ad2"), sep = "; "),
      "candidate_wave5_weight",
      "candidate_feasible",
      "Wave 5 weights exist, but need harmonization with baseline biomarker weights before extending the longitudinal model."
    )
  )

  counts <- bind_rows(
    tibble(metric = "core_harmonized_rows", value = nrow(core)),
    tibble(metric = "core_harmonized_has_r5_variables", value = as.integer(any(grepl("^r5", names(core))))),
    tibble(metric = "health_2020_rows", value = nrow(health)),
    tibble(metric = "sample_info_2020_rows", value = nrow(sample_info)),
    tibble(metric = "exit_module_2020_rows", value = nrow(exit_mod)),
    tibble(metric = "weights_2020_rows", value = nrow(weights)),
    tibble(metric = "health_2020_ids_overlapping_core", value = length(intersect(health_ids, core_ids))),
    tibble(metric = "sample_info_2020_ids_overlapping_core", value = length(intersect(sample_ids, core_ids))),
    tibble(metric = "exit_module_2020_ids_overlapping_core", value = length(intersect(exit_ids, core_ids))),
    tibble(metric = "weights_2020_ids_overlapping_core", value = length(intersect(weight_ids, core_ids))),
    tibble(metric = "sample_info_died_yes", value = sum(to_num(sample_info$died) == 1, na.rm = TRUE))
  )

  value_counts <- bind_rows(
    count_values(health, "da002_5_"),
    count_values(health, "da003_5_"),
    count_values(health, "da004_5_"),
    count_values(sample_info, "died"),
    count_values(exit_mod, "xezdisease_2_")
  )

  feasibility_path <- file.path(metadata_dir, "charls_wave5_endpoint_feasibility_v0_2.csv")
  counts_path <- file.path(table_dir, "charls_wave5_candidate_counts_v0_2.csv")
  value_counts_path <- file.path(table_dir, "charls_wave5_candidate_value_counts_v0_2.csv")
  log_path <- file.path(log_dir, "charls_wave5_endpoint_feasibility.md")

  readr::write_csv(feasibility, feasibility_path)
  readr::write_csv(counts, counts_path)
  readr::write_csv(value_counts, value_counts_path)

  log <- c(
    "# CHARLS Wave 5 Endpoint Feasibility",
    "",
    paste0("- Run date: ", Sys.Date()),
    "- Local 2020 module archive was found and inspected at aggregate level only.",
    "- Current harmonized core RDS has no `r5*` variables.",
    "- Harmonized CHARLS D documentation is version D (2011-2018), through Wave 4.",
    "- Decision for current manuscript: keep the primary CHARLS analysis as 2011-2018 / through Wave 4.",
    "- Wave 5 can be added later only after a separate raw-2020 codebook lock, value-label lock, linkage audit, and weight harmonization.",
    "",
    "## Candidate Evidence",
    "",
    paste0("- Health 2020 rows: ", nrow(health), "; IDs overlapping current core: ", length(intersect(health_ids, core_ids)), "."),
    paste0("- Sample information 2020 rows: ", nrow(sample_info), "; IDs overlapping current core: ", length(intersect(sample_ids, core_ids)), "."),
    paste0("- Exit module 2020 rows: ", nrow(exit_mod), "; IDs overlapping current core: ", length(intersect(exit_ids, core_ids)), "."),
    paste0("- Wave 5 death indicator (`died == 1`) count in Sample_Infor: ", sum(to_num(sample_info$died) == 1, na.rm = TRUE), "."),
    "",
    "## Boundary",
    "",
    "The 2020 disease-array candidates are feasible but not locked. They must not be mixed into the current Wave 1-4 Harmonized D endpoint until the raw 2020 disease index and value-label rules are fully documented."
  )
  writeLines(log, log_path)

  message("Wrote CHARLS Wave 5 feasibility outputs.")
  message("Feasibility: ", feasibility_path)
  message("Counts: ", counts_path)
  message("Log: ", log_path)
}

if (sys.nframe() == 0) {
  main()
}

