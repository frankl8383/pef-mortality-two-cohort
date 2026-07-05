# NHANES V30 linked mortality ETL.
# Builds local-only row-level NHANES 2007-2012 + NCHS 2019 public-use
# linked mortality analysis data, plus aggregate audit tables.

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

is_true <- function(x) {
  !is.na(x) & (x == TRUE | x == 1)
}

fmt_num <- function(x, digits = 3) {
  ifelse(is.na(x), "", formatC(x, digits = digits, format = "f"))
}

markdown_table <- function(dat) {
  dat <- as.data.frame(dat, stringsAsFactors = FALSE)
  dat[] <- lapply(dat, function(x) ifelse(is.na(x), "", as.character(x)))
  header <- paste0("| ", paste(names(dat), collapse = " | "), " |")
  separator <- paste0("| ", paste(rep("---", ncol(dat)), collapse = " | "), " |")
  rows <- apply(dat, 1, function(x) paste0("| ", paste(x, collapse = " | "), " |"))
  c(header, separator, rows)
}

download_if_missing <- function(url, dest) {
  if (file.exists(dest) && file.info(dest)$size > 0) {
    return("present")
  }
  tmp <- paste0(dest, ".tmp")
  if (file.exists(tmp)) {
    unlink(tmp)
  }
  status <- tryCatch({
    utils::download.file(url, tmp, mode = "wb", quiet = TRUE)
    file.rename(tmp, dest)
    "downloaded"
  }, error = function(e) {
    if (file.exists(tmp)) {
      unlink(tmp)
    }
    paste0("failed: ", conditionMessage(e))
  })
  status
}

file_md5 <- function(path) {
  if (!file.exists(path)) {
    return(NA_character_)
  }
  unname(tools::md5sum(path))
}

read_nhanes_mortality <- function(path) {
  out <- readr::read_fwf(
    file = path,
    col_types = "iiiiiiii",
    col_positions = readr::fwf_cols(
      seqn = c(1, 6),
      eligstat = c(15, 15),
      mortstat = c(16, 16),
      ucod_leading = c(17, 19),
      diabetes = c(20, 20),
      hyperten = c(21, 21),
      permth_int = c(43, 45),
      permth_exm = c(46, 48)
    ),
    na = c("", ".")
  )
  out$mortality_source_file <- basename(path)
  out
}

source_file_to_cycle <- function(file_name) {
  dplyr::case_when(
    grepl("2007_2008", file_name) ~ "2007-2008",
    grepl("2009_2010", file_name) ~ "2009-2010",
    grepl("2011_2012", file_name) ~ "2011-2012",
    TRUE ~ NA_character_
  )
}

source_file_to_cycle_label <- function(file_name) {
  dplyr::case_when(
    grepl("2007_2008", file_name) ~ "E",
    grepl("2009_2010", file_name) ~ "F",
    grepl("2011_2012", file_name) ~ "G",
    TRUE ~ NA_character_
  )
}

ucod_label <- function(x) {
  dplyr::case_when(
    x == 1L ~ "Diseases of heart",
    x == 2L ~ "Malignant neoplasms",
    x == 3L ~ "Chronic lower respiratory diseases",
    x == 4L ~ "Accidents/unintentional injuries",
    x == 5L ~ "Cerebrovascular diseases",
    x == 6L ~ "Alzheimer's disease",
    x == 7L ~ "Diabetes mellitus",
    x == 8L ~ "Influenza and pneumonia",
    x == 9L ~ "Nephritis/nephrotic syndrome/nephrosis",
    x == 10L ~ "All other causes",
    TRUE ~ NA_character_
  )
}

summarise_mortality <- function(data, stratum_type, stratum) {
  tibble(
    stratum_type = stratum_type,
    stratum = stratum,
    n = nrow(data),
    all_cause_deaths = sum(data$all_cause_death == 1L, na.rm = TRUE),
    all_cause_death_rate = ifelse(nrow(data) > 0, all_cause_deaths / nrow(data), NA_real_),
    clrd_deaths = sum(data$clrd_death == 1L, na.rm = TRUE),
    clrd_death_rate = ifelse(nrow(data) > 0, clrd_deaths / nrow(data), NA_real_),
    followup_years_median = suppressWarnings(stats::median(data$followup_years_exam, na.rm = TRUE)),
    followup_years_p25 = suppressWarnings(stats::quantile(data$followup_years_exam, 0.25, na.rm = TRUE, names = FALSE)),
    followup_years_p75 = suppressWarnings(stats::quantile(data$followup_years_exam, 0.75, na.rm = TRUE, names = FALSE)),
    person_years = sum(data$followup_years_exam, na.rm = TRUE)
  )
}

root <- find_project_root()
derived_dir <- file.path(root, "derived_sensitive", "nhanes")
mortality_dir <- file.path(derived_dir, "linked_mortality_2019")
results_dir <- file.path(root, "results", "tables")
logs_dir <- file.path(root, "results", "logs")
metadata_dir <- file.path(root, "metadata")
dir_create(derived_dir)
dir_create(mortality_dir)
dir_create(results_dir)
dir_create(logs_dir)
dir_create(metadata_dir)

ftp_base <- "https://ftp.cdc.gov/pub/Health_Statistics/NCHS/datalinkage/linked_mortality"
source_manifest <- tibble(
  file_name = c(
    "NHANES_2007_2008_MORT_2019_PUBLIC.dat",
    "NHANES_2009_2010_MORT_2019_PUBLIC.dat",
    "NHANES_2011_2012_MORT_2019_PUBLIC.dat",
    "R_ReadInProgramAllSurveys.R",
    "public-use-linked-mortality-file-description.pdf",
    "public-use-linked-mortality-files-data-dictionary.pdf"
  ),
  source_url = c(
    paste0(ftp_base, "/NHANES_2007_2008_MORT_2019_PUBLIC.dat"),
    paste0(ftp_base, "/NHANES_2009_2010_MORT_2019_PUBLIC.dat"),
    paste0(ftp_base, "/NHANES_2011_2012_MORT_2019_PUBLIC.dat"),
    paste0(ftp_base, "/R_ReadInProgramAllSurveys.R"),
    "https://www.cdc.gov/nchs/data/datalinkage/public-use-linked-mortality-file-description.pdf",
    "https://www.cdc.gov/nchs/data/datalinkage/public-use-linked-mortality-files-data-dictionary.pdf"
  ),
  source_role = c(
    "mortality_data",
    "mortality_data",
    "mortality_data",
    "official_read_program",
    "official_file_description",
    "official_data_dictionary"
  )
)

source_manifest <- source_manifest %>%
  mutate(
    local_file = file.path(mortality_dir, file_name),
    download_status = vapply(seq_along(source_url), function(i) {
      download_if_missing(source_url[[i]], local_file[[i]])
    }, character(1)),
    size_bytes = ifelse(file.exists(local_file), file.info(local_file)$size, NA_real_),
    md5 = vapply(local_file, file_md5, character(1))
  )

manifest_public <- source_manifest %>%
  transmute(
    file_name,
    source_url,
    source_role,
    download_status,
    size_bytes,
    md5
  )
readr::write_csv(
  manifest_public,
  file.path(results_dir, paste0("nhanes_mortality_download_manifest_", version_id, ".csv"))
)

analytic_path <- file.path(derived_dir, "nhanes_replication_v0_5_gli_phenotypes.rds")
if (!file.exists(analytic_path)) {
  stop("Missing NHANES v0.5 GLI phenotype dataset.", call. = FALSE)
}
nhanes <- readRDS(analytic_path)

mortality_files <- source_manifest %>%
  filter(source_role == "mortality_data") %>%
  pull(local_file)
if (!all(file.exists(mortality_files))) {
  stop("One or more NCHS linked mortality data files are missing.", call. = FALSE)
}

mortality <- bind_rows(lapply(mortality_files, read_nhanes_mortality)) %>%
  mutate(
    linked_mortality_cycle = source_file_to_cycle(mortality_source_file),
    linked_mortality_cycle_label = source_file_to_cycle_label(mortality_source_file),
    eligstat_label = dplyr::case_when(
      eligstat == 1L ~ "Eligible",
      eligstat == 2L ~ "Under age 18, not available for public release",
      eligstat == 3L ~ "Ineligible",
      TRUE ~ NA_character_
    ),
    mortstat_label = dplyr::case_when(
      mortstat == 0L ~ "Assumed alive",
      mortstat == 1L ~ "Assumed deceased",
      TRUE ~ NA_character_
    ),
    ucod_leading_label = ucod_label(ucod_leading),
    mortality_eligible = as.integer(eligstat == 1L),
    all_cause_death = dplyr::case_when(
      eligstat == 1L & mortstat == 1L ~ 1L,
      eligstat == 1L & mortstat == 0L ~ 0L,
      TRUE ~ NA_integer_
    ),
    clrd_death = dplyr::case_when(
      eligstat == 1L & mortstat == 1L & ucod_leading == 3L ~ 1L,
      eligstat == 1L ~ 0L,
      TRUE ~ NA_integer_
    ),
    diabetes_mcod_death = dplyr::case_when(
      eligstat == 1L & mortstat == 1L & diabetes == 1L ~ 1L,
      eligstat == 1L & mortstat == 1L & diabetes == 0L ~ 0L,
      eligstat == 1L & mortstat == 0L ~ 0L,
      TRUE ~ NA_integer_
    ),
    hypertension_mcod_death = dplyr::case_when(
      eligstat == 1L & mortstat == 1L & hyperten == 1L ~ 1L,
      eligstat == 1L & mortstat == 1L & hyperten == 0L ~ 0L,
      eligstat == 1L & mortstat == 0L ~ 0L,
      TRUE ~ NA_integer_
    ),
    followup_months_interview = permth_int,
    followup_months_exam = permth_exm,
    followup_years_interview = permth_int / 12,
    followup_years_exam = permth_exm / 12
  )

analysis_ready <- nhanes %>%
  left_join(mortality, by = c("participant_id" = "seqn")) %>%
  mutate(
    mortality_merge_status = ifelse(!is.na(eligstat), "matched", "not_matched"),
    mortality_analysis_primary = as.integer(
      is_true(adult_45plus) &
        pef_quality_abc == 1L &
        !is.na(resp_vulnerability_z) &
        !is.na(wtmec6yr) & wtmec6yr > 0 &
        mortality_eligible == 1L &
        !is.na(all_cause_death) &
        !is.na(followup_years_exam) & followup_years_exam > 0
    ),
    mortality_time_origin = ifelse(
      mortality_analysis_primary == 1L,
      "NHANES MEC examination date",
      NA_character_
    ),
    mortality_followup_end = ifelse(
      mortality_analysis_primary == 1L,
      "Death or 2019-12-31 public-use linkage end",
      NA_character_
    )
  )

primary <- analysis_ready %>% filter(mortality_analysis_primary == 1L)
adult45 <- analysis_ready %>% filter(is_true(adult_45plus))
adult45_eligible <- adult45 %>% filter(mortality_eligible == 1L)
adult45_vulnerability <- adult45 %>% filter(!is.na(resp_vulnerability_z))
adult45_vulnerability_eligible <- adult45_vulnerability %>% filter(mortality_eligible == 1L)

flow_n <- c(
  nrow(analysis_ready),
  nrow(mortality),
  sum(analysis_ready$mortality_merge_status == "matched", na.rm = TRUE),
  nrow(adult45),
  nrow(adult45_eligible),
  nrow(adult45_vulnerability),
  nrow(adult45_vulnerability_eligible),
  nrow(primary)
)
flow_deaths <- c(
  sum(analysis_ready$all_cause_death == 1L, na.rm = TRUE),
  sum(mortality$all_cause_death == 1L, na.rm = TRUE),
  sum(analysis_ready$all_cause_death == 1L, na.rm = TRUE),
  sum(adult45$all_cause_death == 1L, na.rm = TRUE),
  sum(adult45_eligible$all_cause_death == 1L, na.rm = TRUE),
  sum(adult45_vulnerability$all_cause_death == 1L, na.rm = TRUE),
  sum(adult45_vulnerability_eligible$all_cause_death == 1L, na.rm = TRUE),
  sum(primary$all_cause_death == 1L, na.rm = TRUE)
)
flow_clrd <- c(
  sum(analysis_ready$clrd_death == 1L, na.rm = TRUE),
  sum(mortality$clrd_death == 1L, na.rm = TRUE),
  sum(analysis_ready$clrd_death == 1L, na.rm = TRUE),
  sum(adult45$clrd_death == 1L, na.rm = TRUE),
  sum(adult45_eligible$clrd_death == 1L, na.rm = TRUE),
  sum(adult45_vulnerability$clrd_death == 1L, na.rm = TRUE),
  sum(adult45_vulnerability_eligible$clrd_death == 1L, na.rm = TRUE),
  sum(primary$clrd_death == 1L, na.rm = TRUE)
)

flow <- tibble(
  step_order = seq_along(flow_n),
  step = c(
    "NHANES v0.5 analytic rows",
    "NCHS 2019 public-use mortality rows",
    "Rows matched by SEQN",
    "Age 45 years or older",
    "Age 45 years or older and mortality eligible",
    "Age 45 years or older with respiratory vulnerability marker available",
    "Age 45 years or older with marker available and mortality eligible",
    "Primary mortality analysis cohort"
  ),
  n = flow_n,
  excluded_from_previous = c(NA_integer_, NA_integer_, flow_n[1] - flow_n[3], diff(flow_n[3:8]) * -1L),
  all_cause_deaths = flow_deaths,
  clrd_deaths = flow_clrd,
  notes = c(
    "Local-only row-level NHANES GLI phenotype table.",
    "Three cycle-specific public-use linked mortality files: 2007-2008, 2009-2010, 2011-2012.",
    "Merge uses NHANES SEQN / participant_id.",
    "Age restriction follows the current respiratory vulnerability manuscript analysis.",
    "ELIGSTAT == 1.",
    "Primary PEF residual marker available.",
    "Eligible for mortality follow-up.",
    "Adult 45+, PEF marker available, MEC weight positive, mortality eligible, positive MEC follow-up."
  )
)

overall_counts <- bind_rows(
  summarise_mortality(analysis_ready, "cohort", "all_nhanes_rows"),
  summarise_mortality(adult45, "cohort", "adult_45plus"),
  summarise_mortality(adult45_eligible, "cohort", "adult_45plus_mortality_eligible"),
  summarise_mortality(primary, "cohort", "primary_mortality_analysis")
)

cycle_counts <- primary %>%
  group_by(cycle_label, cycle) %>%
  group_modify(~ summarise_mortality(.x, "cycle", paste0(.y$cycle_label, " / ", .y$cycle))) %>%
  ungroup() %>%
  select(-cycle_label, -cycle)

quartile_counts <- primary %>%
  group_by(rv_quartile_num, rv_quartile) %>%
  group_modify(~ summarise_mortality(.x, "respiratory_vulnerability_quartile", paste0("Q", .y$rv_quartile_num, ": ", .y$rv_quartile))) %>%
  ungroup() %>%
  select(-rv_quartile_num, -rv_quartile)

cause_counts <- primary %>%
  filter(all_cause_death == 1L) %>%
  count(ucod_leading, ucod_leading_label, name = "deaths") %>%
  arrange(ucod_leading) %>%
  mutate(
    stratum_type = "underlying_cause_among_primary_deaths",
    stratum = paste0(sprintf("%03d", ucod_leading), ": ", ucod_leading_label),
    n = nrow(primary),
    all_cause_deaths = deaths,
    all_cause_death_rate = deaths / nrow(primary),
    clrd_deaths = ifelse(ucod_leading == 3L, deaths, 0L),
    clrd_death_rate = clrd_deaths / nrow(primary),
    followup_years_median = NA_real_,
    followup_years_p25 = NA_real_,
    followup_years_p75 = NA_real_,
    person_years = NA_real_
  ) %>%
  select(stratum_type, stratum, n, all_cause_deaths, all_cause_death_rate, clrd_deaths,
         clrd_death_rate, followup_years_median, followup_years_p25, followup_years_p75,
         person_years)

counts <- bind_rows(overall_counts, cycle_counts, quartile_counts, cause_counts) %>%
  mutate(
    all_cause_death_rate = round(all_cause_death_rate, 5),
    clrd_death_rate = round(clrd_death_rate, 5),
    followup_years_median = round(followup_years_median, 3),
    followup_years_p25 = round(followup_years_p25, 3),
    followup_years_p75 = round(followup_years_p75, 3),
    person_years = round(person_years, 1)
  )

row_output_rds <- file.path(derived_dir, paste0("nhanes_mortality_analysis_ready_", version_id, ".rds"))
row_output_rds_latest <- file.path(derived_dir, "nhanes_mortality_analysis_ready.rds")
row_output_csv <- file.path(derived_dir, paste0("nhanes_mortality_analysis_ready_", version_id, ".csv.gz"))
row_output_csv_latest <- file.path(derived_dir, "nhanes_mortality_analysis_ready.csv.gz")
row_output_parquet <- file.path(derived_dir, paste0("nhanes_mortality_analysis_ready_", version_id, ".parquet"))
row_output_parquet_latest <- file.path(derived_dir, "nhanes_mortality_analysis_ready.parquet")

saveRDS(analysis_ready, row_output_rds)
file.copy(row_output_rds, row_output_rds_latest, overwrite = TRUE)
readr::write_csv(analysis_ready, row_output_csv)
file.copy(row_output_csv, row_output_csv_latest, overwrite = TRUE)

parquet_status <- "not_attempted"
python <- Sys.getenv("NHANES_PARQUET_PYTHON", unset = Sys.which("python3"))
if (nzchar(python)) {
  py_script <- tempfile(fileext = ".py")
  writeLines(c(
    "import sys",
    "import pandas as pd",
    "src, dst = sys.argv[1], sys.argv[2]",
    "df = pd.read_csv(src, low_memory=False)",
    "df.to_parquet(dst, index=False)"
  ), py_script)
  status <- suppressWarnings(system2(python, c(py_script, row_output_csv, row_output_parquet), stdout = TRUE, stderr = TRUE))
  exit_status <- attr(status, "status")
  if (is.null(exit_status) && file.exists(row_output_parquet) && file.info(row_output_parquet)$size > 0) {
    file.copy(row_output_parquet, row_output_parquet_latest, overwrite = TRUE)
    parquet_status <- "written"
  } else {
    parquet_status <- paste0("failed: ", paste(status, collapse = " | "))
  }
}

flow_path <- file.path(results_dir, paste0("nhanes_mortality_flow_", version_id, ".csv"))
counts_path <- file.path(results_dir, paste0("nhanes_mortality_counts_", version_id, ".csv"))
readr::write_csv(flow, flow_path)
readr::write_csv(counts, counts_path)

codebook_path <- file.path(metadata_dir, paste0("nhanes_mortality_codebook_", version_id, ".md"))
log_path <- file.path(logs_dir, paste0("nhanes_mortality_etl_", version_id, ".md"))

variable_codebook <- tibble(
  variable = c(
    "mortality_merge_status",
    "eligstat",
    "mortality_eligible",
    "mortstat",
    "all_cause_death",
    "ucod_leading",
    "ucod_leading_label",
    "clrd_death",
    "diabetes_mcod_death",
    "hypertension_mcod_death",
    "followup_months_exam",
    "followup_years_exam",
    "followup_months_interview",
    "followup_years_interview",
    "mortality_analysis_primary"
  ),
  definition = c(
    "Matched/not matched status after joining NCHS public-use linked mortality by SEQN.",
    "NCHS eligibility status for mortality follow-up: 1 eligible, 2 under age 18 not public-release mortality follow-up, 3 ineligible.",
    "Indicator for ELIGSTAT == 1.",
    "NCHS final vital status among eligible participants: 0 assumed alive, 1 assumed deceased.",
    "Primary event indicator for all-cause death among mortality-eligible participants.",
    "NCHS public-use leading underlying cause of death recode.",
    "Text label for UCOD_LEADING.",
    "Exploratory chronic lower respiratory disease death indicator; 1 only when eligible, deceased, and UCOD_LEADING == 3.",
    "Public-use diabetes multiple-cause-of-death flag mapped to an event indicator for deceased eligible participants.",
    "Public-use hypertension multiple-cause-of-death flag mapped to an event indicator for deceased eligible participants.",
    "NCHS person-months of follow-up from NHANES MEC examination date.",
    "followup_months_exam / 12.",
    "NCHS person-months of follow-up from NHANES interview date.",
    "followup_months_interview / 12.",
    "Adult 45+, PEF vulnerability marker available, MEC weight positive, mortality eligible, nonmissing all-cause death, and positive MEC follow-up."
  ),
  analysis_role = c(
    "merge audit",
    "eligibility audit",
    "analysis eligibility",
    "vital status",
    "T2 primary outcome",
    "exploratory cause-specific outcome source",
    "interpretation",
    "T2 exploratory respiratory mortality outcome if event count is adequate",
    "descriptive/sensitivity candidate only",
    "descriptive/sensitivity candidate only",
    "T2 primary survival time scale",
    "T2 primary survival time scale",
    "sensitivity/backup time scale",
    "sensitivity/backup time scale",
    "primary T2 analysis flag"
  )
)

codebook_lines <- c(
  "# NHANES V30 Linked Mortality Codebook",
  "",
  paste0("Version: ", version_id),
  "",
  "## Provenance",
  "",
  "This ETL links the existing local NHANES 2007-2012 respiratory vulnerability analysis table to the NCHS 2019 public-use linked mortality files. The row-level outputs are local-only and are not for manuscript upload. Public manuscript-facing outputs are aggregate counts and flow tables.",
  "",
  "Official source files are recorded in `results/tables/nhanes_mortality_download_manifest_v30_0.csv`. The official NCHS R read program defines SEQN and fixed-width fields for ELIGSTAT, MORTSTAT, UCOD_LEADING, DIABETES, HYPERTEN, PERMTH_INT, and PERMTH_EXM.",
  "",
  "## Row-Level Local Outputs",
  "",
  "- `derived_sensitive/nhanes/nhanes_mortality_analysis_ready_v30_0.rds`",
  "- `derived_sensitive/nhanes/nhanes_mortality_analysis_ready_v30_0.csv.gz`",
  paste0("- `derived_sensitive/nhanes/nhanes_mortality_analysis_ready_v30_0.parquet`: ", parquet_status),
  "",
  "## Aggregate Outputs",
  "",
  "- `results/tables/nhanes_mortality_flow_v30_0.csv`",
  "- `results/tables/nhanes_mortality_counts_v30_0.csv`",
  "- `results/tables/nhanes_mortality_download_manifest_v30_0.csv`",
  "",
  "## Variable Definitions",
  "",
  markdown_table(variable_codebook),
  "",
  "## Analysis Boundary",
  "",
  "All-cause mortality is the primary hard outcome for T2. Chronic lower respiratory disease mortality is retained as an exploratory cause-specific outcome only because the current primary mortality cohort has a modest event count.",
  "",
  "The primary survival time scale for T2 should use `followup_years_exam`, because PEF/spirometry measurements come from the NHANES MEC examination."
)
writeLines(codebook_lines, codebook_path)

primary_summary <- tibble(
  metric = c(
    "nhanes_rows",
    "mortality_rows",
    "matched_rows",
    "adult_45plus_rows",
    "primary_mortality_analysis_rows",
    "primary_all_cause_deaths",
    "primary_clrd_deaths",
    "primary_median_followup_years",
    "primary_person_years",
    "parquet_status"
  ),
  value = c(
    as.character(nrow(analysis_ready)),
    as.character(nrow(mortality)),
    as.character(sum(analysis_ready$mortality_merge_status == "matched", na.rm = TRUE)),
    as.character(nrow(adult45)),
    as.character(nrow(primary)),
    as.character(sum(primary$all_cause_death == 1L, na.rm = TRUE)),
    as.character(sum(primary$clrd_death == 1L, na.rm = TRUE)),
    fmt_num(stats::median(primary$followup_years_exam, na.rm = TRUE), 2),
    fmt_num(sum(primary$followup_years_exam, na.rm = TRUE), 1),
    parquet_status
  )
)

log_lines <- c(
  "# NHANES V30 Linked Mortality ETL Log",
  "",
  paste0("Version: ", version_id),
  "",
  "## Status",
  "",
  "T1 linked mortality ETL completed. This log is aggregate-only; row-level outputs remain in the local-only derived data directory.",
  "",
  "## Primary Summary",
  "",
  markdown_table(primary_summary),
  "",
  "## Cohort Flow",
  "",
  markdown_table(flow),
  "",
  "## Notes",
  "",
  "- Merge key: NHANES SEQN / `participant_id`.",
  "- Primary outcome for T2: all-cause death (`all_cause_death`).",
  "- Primary time scale for T2: follow-up years from MEC examination (`followup_years_exam`).",
  "- Exploratory respiratory mortality flag: chronic lower respiratory diseases (`ucod_leading == 3`).",
  "- Public-use linked mortality files have known public-use perturbation/synthetic handling for some follow-up and cause fields; vital status is the key hard outcome anchor."
)
writeLines(log_lines, log_path)

message("NHANES linked mortality ETL complete.")
message("Primary mortality analysis rows: ", nrow(primary))
message("Primary all-cause deaths: ", sum(primary$all_cause_death == 1L, na.rm = TRUE))
message("Primary CLRD deaths: ", sum(primary$clrd_death == 1L, na.rm = TRUE))
message("Parquet status: ", parquet_status)
