# Build ICU V31 subgroup transportability summary from existing aggregate outputs.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
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

dir_create <- function(path) if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)

md_escape <- function(x) str_replace_all(ifelse(is.na(x), "", as.character(x)), "\\|", "\\\\|")
md_table <- function(dat) {
  dat <- as.data.frame(dat, stringsAsFactors = FALSE)
  dat[] <- lapply(dat, function(x) ifelse(is.na(x), "", as.character(x)))
  c(
    paste0("| ", paste(md_escape(names(dat)), collapse = " | "), " |"),
    paste0("| ", paste(rep("---", ncol(dat)), collapse = " | "), " |"),
    apply(dat, 1, function(x) paste0("| ", paste(md_escape(x), collapse = " | "), " |"))
  )
}

root <- find_project_root()
tab_dir <- file.path(root, "results", "tables")
log_dir <- file.path(root, "results", "logs")
manuscript_dir <- file.path(root, "manuscript")
invisible(lapply(c(tab_dir, log_dir, manuscript_dir), dir_create))

src <- read_csv(file.path(tab_dir, "icu_transportability_summary_v30_0.csv"), show_col_types = FALSE)

available <- src %>%
  mutate(
    requested_domain = case_when(
      group_family == "support_class" ~ "HFNC/NIV support class",
      group_family == "age_group" ~ "Age group",
      group_family == "sex" ~ "Sex",
      TRUE ~ group_family
    ),
    model_label = recode(model, dynamic_last_no_support = "Dynamic model", rox_only = "ROX-like comparator"),
    subgroup_label = recode(subgroup,
      hfnc = "HFNC-coded",
      niv = "NIV-coded",
      age_lt65 = "<65 years",
      age_65_79 = "65-79 years",
      age_ge80 = ">=80 years",
      female = "Female",
      male = "Male",
      .default = subgroup
    ),
    availability_status = "available_existing_v30_aggregate",
    manuscript_boundary = if_else(
      group_family == "support_class" & subgroup == "hfnc" & evaluation == "eicu_external",
      "Sparse eICU HFNC-coded stratum; descriptive only, not HFNC-specific validation.",
      "Descriptive transportability only; supportive ICU translation."
    )
  ) %>%
  transmute(
    requested_domain,
    evaluation,
    model,
    model_label,
    subgroup = subgroup_label,
    n,
    events,
    event_rate,
    auroc,
    auprc,
    brier,
    calibration_intercept_approx,
    calibration_slope,
    interpretation_flag,
    availability_status,
    manuscript_boundary
  )

missing_requested <- tibble::tribble(
  ~requested_domain, ~evaluation, ~model, ~model_label, ~subgroup, ~n, ~events, ~event_rate, ~auroc, ~auprc, ~brier, ~calibration_intercept_approx, ~calibration_slope, ~interpretation_flag, ~availability_status, ~manuscript_boundary,
  "Chronic pulmonary disease", NA_character_, NA_character_, NA_character_, "present/absent", NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, "not_estimable_from_current_aggregate", "requires_row_level_rerun", "P1 target not available in current v30 aggregate subgroup table; do not imply completed.",
  "Hypoxemic/hypercapnic phenotype", NA_character_, NA_character_, NA_character_, "hypoxemic/hypercapnic", NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, "not_estimable_from_current_aggregate", "requires_row_level_rerun", "P1 target not available in current v30 aggregate subgroup table; do not imply completed."
)

out <- bind_rows(available, missing_requested)
write_csv(out, file.path(tab_dir, "icu_v31_subgroup_transportability.csv"))

note <- c(
  "# ICU Subgroup Transportability V31.5",
  "",
  "## Summary",
  "",
  "Existing aggregate ICU transportability outputs support descriptive subgroup reporting for support class, age and sex. Chronic pulmonary disease and hypoxemic/hypercapnic subgroup transportability were requested by GPT5.5pro but are not present in the current v30 aggregate subgroup table; they require a row-level rerun before any result can be claimed.",
  "",
  "## Available subgroup summary",
  "",
  md_table(out %>%
    mutate(across(c(event_rate, auroc, auprc, brier, calibration_intercept_approx, calibration_slope), ~ ifelse(is.na(.x), "", sprintf("%.3f", .x)))) %>%
    select(requested_domain, evaluation, model_label, subgroup, n, events, event_rate, auroc, auprc, brier, calibration_slope, availability_status, manuscript_boundary)),
  "",
  "## Manuscript boundary",
  "",
  "Use this table as descriptive ICU transportability only. The eICU HFNC-coded dynamic-model stratum has n=566 and should be described as sparse; HFNC-specific external validation is not claimed. No ICU treatment-effect estimate is reported because the TTE feasibility layer remains no-go."
)
writeLines(note, file.path(manuscript_dir, "icu_subgroup_transportability_note_v31_5.md"))

log <- c(
  "# ICU V31 Subgroup Transportability Build",
  "",
  paste0("Available aggregate rows: ", nrow(available)),
  paste0("Requested-but-missing rows recorded: ", nrow(missing_requested)),
  "",
  "Outputs:",
  "- results/tables/icu_v31_subgroup_transportability.csv",
  "- manuscript/icu_subgroup_transportability_note_v31_5.md"
)
writeLines(log, file.path(log_dir, "icu_v31_subgroup_transportability.md"))

message("ICU V31 subgroup transportability outputs written.")
