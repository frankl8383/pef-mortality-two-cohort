# MIMIC-IV v3.1 integrity and schema/cohort feasibility scaffold.
# This script outputs file-level and schema-level artifacts only. It does not
# export row-level clinical data.

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

mimic_root <- function() {
  root <- Sys.getenv("MIMICIV_ROOT", unset = "")
  if (identical(root, "")) root <- file.path(path.expand("~"), "secure_data", "mimiciv", "3.1")
  normalizePath(root, mustWork = FALSE)
}

require_input <- function(path) {
  if (!file.exists(path)) stop("Missing required input: ", path, call. = FALSE)
  path
}

rel_path <- function(path, root) {
  sub(paste0("^", gsub("([.|()\\^{}+$*?]|\\[|\\])", "\\\\\\1", root), "/?"), "", path)
}

latest_problem_files <- function(root) {
  files <- list.files(file.path(root, ".download"), pattern = "^problem_files_.*\\.txt$", full.names = TRUE)
  if (length(files) == 0) return(character())
  latest <- files[order(file.info(files)$mtime, decreasing = TRUE)][[1]]
  readLines(latest, warn = FALSE)
}

active_download <- function() {
  out <- tryCatch(
    system2("ps", c("-axo", "command"), stdout = TRUE, stderr = FALSE),
    error = function(e) character()
  )
  any(str_detect(out, "wget .*physionet\\.org/files/mimiciv/3\\.1|curl .*physionet\\.org/files/mimiciv/3\\.1"))
}

parse_sha256 <- function(root) {
  sha_path <- require_input(file.path(root, "SHA256SUMS.txt"))
  lines <- readLines(sha_path, warn = FALSE)
  dat <- tibble(raw = lines) %>%
    filter(str_detect(.data$raw, "\\S+\\s+\\S+")) %>%
    transmute(
      expected_sha256 = str_extract(.data$raw, "^[0-9a-fA-F]+"),
      file = str_trim(str_remove(.data$raw, "^[0-9a-fA-F]+\\s+"))
    )
  dat
}

gzip_test <- function(path) {
  status <- suppressWarnings(system2("gzip", c("-t", path), stdout = FALSE, stderr = FALSE))
  if (identical(status, 0L)) "pass" else "fail"
}

sha256_file <- function(path) {
  unname(tools::sha256sum(path))
}

read_header <- function(path) {
  con <- gzfile(path, open = "rt")
  on.exit(close(con), add = TRUE)
  header <- readLines(con, n = 1, warn = FALSE)
  if (length(header) == 0) return(character())
  str_split(header, ",", simplify = FALSE)[[1]]
}

schema_expectations <- function() {
  tribble(
    ~table_file, ~role, ~required_columns,
    "hosp/patients.csv.gz", "patient demographics and death date", "subject_id;gender;anchor_age;anchor_year;anchor_year_group;dod",
    "hosp/admissions.csv.gz", "hospital admissions and admission/discharge timing", "subject_id;hadm_id;admittime;dischtime;deathtime;admission_type;admission_location;discharge_location;insurance;race",
    "hosp/transfers.csv.gz", "ward and ICU location timing", "subject_id;hadm_id;transfer_id;eventtype;careunit;intime;outtime",
    "icu/icustays.csv.gz", "ICU stay denominator and timing", "subject_id;hadm_id;stay_id;first_careunit;last_careunit;intime;outtime;los",
    "icu/d_items.csv.gz", "ICU item dictionary for respiratory-support item lock", "itemid;label;abbreviation;linksto;category;unitname;param_type",
    "icu/chartevents.csv.gz", "bedside respiratory measurements and ventilator/oxygen charting", "subject_id;hadm_id;stay_id;charttime;storetime;itemid;value;valuenum;valueuom",
    "icu/inputevents.csv.gz", "ICU inputs and device-related event context", "subject_id;hadm_id;stay_id;starttime;endtime;itemid;amount;amountuom;rate;rateuom",
    "icu/outputevents.csv.gz", "ICU outputs and support context", "subject_id;hadm_id;stay_id;charttime;itemid;value;valueuom",
    "icu/procedureevents.csv.gz", "ICU procedure timing and respiratory-support procedure candidates", "subject_id;hadm_id;stay_id;starttime;endtime;itemid;value;valueuom",
    "hosp/diagnoses_icd.csv.gz", "hospital diagnosis code phenotyping", "subject_id;hadm_id;seq_num;icd_code;icd_version",
    "hosp/procedures_icd.csv.gz", "hospital procedure code phenotyping", "subject_id;hadm_id;seq_num;chartdate;icd_code;icd_version",
    "hosp/d_icd_diagnoses.csv.gz", "diagnosis code dictionary", "icd_code;icd_version;long_title",
    "hosp/d_icd_procedures.csv.gz", "procedure code dictionary", "icd_code;icd_version;long_title",
    "hosp/labevents.csv.gz", "laboratory respiratory physiology support variables", "labevent_id;subject_id;hadm_id;itemid;charttime;valuenum;valueuom",
    "hosp/d_labitems.csv.gz", "laboratory item dictionary", "itemid;label;fluid;category"
  )
}

candidate_patterns <- function() {
  c(
    "ventilator",
    "ventilation",
    "vent mode",
    "oxygen",
    "o2",
    "fio2",
    "peep",
    "cpap",
    "bipap",
    "high flow",
    "flow rate",
    "intub",
    "extub",
    "respiratory rate",
    "tidal volume",
    "minute ventilation",
    "plateau pressure",
    "driving pressure"
  )
}

read_d_items_candidates <- function(root, manifest) {
  row <- manifest %>% filter(.data$file == "icu/d_items.csv.gz")
  if (nrow(row) != 1 || row$gzip_status != "pass") {
    return(tibble(
      itemid = integer(),
      label = character(),
      abbreviation = character(),
      linksto = character(),
      category = character(),
      unitname = character(),
      matched_pattern = character()
    ))
  }
  d_items <- readr::read_csv(file.path(root, "icu", "d_items.csv.gz"), show_col_types = FALSE)
  cols <- intersect(c("itemid", "label", "abbreviation", "linksto", "category", "unitname"), names(d_items))
  pat <- candidate_patterns()
  d_items %>%
    mutate(search_text = str_to_lower(paste(!!!rlang::syms(cols), sep = " "))) %>%
    rowwise() %>%
    mutate(matched_pattern = paste(pat[str_detect(.data$search_text, fixed(pat))], collapse = ";")) %>%
    ungroup() %>%
    filter(.data$matched_pattern != "") %>%
    select(any_of(c("itemid", "label", "abbreviation", "linksto", "category", "unitname")), matched_pattern) %>%
    arrange(.data$linksto, .data$category, .data$label)
}

build_manifest <- function(root, run_sha = "auto") {
  sha <- parse_sha256(root)
  problem_files <- latest_problem_files(root)
  active <- active_download()
  pending_session <- active && length(problem_files) > 0
  pending_active_set <- if (active) problem_files else character()

  out <- sha %>%
    mutate(
      exists = file.exists(file.path(root, .data$file)),
      size_bytes = ifelse(.data$exists, file.info(file.path(root, .data$file))$size, NA_real_),
      is_gzip = str_ends(.data$file, "\\.gz"),
      download_status = case_when(
        !.data$exists ~ "missing",
        .data$file %in% pending_active_set ~ "pending_active_download",
        TRUE ~ "present_not_active"
      ),
      gzip_status = ifelse(.data$is_gzip, "not_tested_pending", "not_applicable"),
      sha256_status = "not_tested_pending"
    )

  if (!pending_session) {
    gzip_idx <- which(out$is_gzip & out$download_status == "present_not_active")
    for (idx in gzip_idx) {
      out$gzip_status[[idx]] <- gzip_test(file.path(root, out$file[[idx]]))
    }

    if (run_sha == "never") {
      out$sha256_status[out$download_status == "present_not_active"] <- "not_tested_disabled"
    } else {
      sha_idx <- which(out$download_status == "present_not_active")
      for (idx in sha_idx) {
        observed <- sha256_file(file.path(root, out$file[[idx]]))
        out$sha256_status[[idx]] <- ifelse(identical(observed, out$expected_sha256[[idx]]), "pass", "fail")
      }
    }
  }

  out %>%
    select(file, exists, size_bytes, download_status, gzip_status, sha256_status, expected_sha256)
}

build_schema_columns <- function(root, manifest) {
  idx <- which(manifest$exists & manifest$gzip_status == "pass" & str_detect(manifest$file, "\\.csv\\.gz$"))
  if (length(idx) == 0) {
    return(tibble(file = character(), columns = character()))
  }
  rows <- lapply(idx, function(i) {
    tibble(
      file = manifest$file[[i]],
      columns = paste(read_header(file.path(root, manifest$file[[i]])), collapse = ";")
    )
  })
  bind_rows(rows)
}

build_schema_expectation_check <- function(expectations, columns) {
  expectations %>%
    left_join(columns, by = c("table_file" = "file")) %>%
    rowwise() %>%
    mutate(
      available = !is.na(.data$columns),
      missing_columns = ifelse(
        .data$available,
        paste(setdiff(str_split(.data$required_columns, ";", simplify = FALSE)[[1]], str_split(.data$columns, ";", simplify = FALSE)[[1]]), collapse = ";"),
        "table_not_available_or_not_validated"
      ),
      schema_status = case_when(
        !.data$available ~ "pending",
        .data$missing_columns == "" ~ "pass",
        TRUE ~ "missing_expected_columns"
      )
    ) %>%
    ungroup()
}

build_endpoint_feasibility <- function(manifest, schema_check, candidates) {
  has_pass <- function(file) {
    any(manifest$file == file & manifest$gzip_status == "pass")
  }
  schema_pass <- function(file) {
    any(schema_check$table_file == file & schema_check$schema_status == "pass")
  }
  integrity_failures <- manifest %>%
    filter(
      .data$download_status != "present_not_active" |
        !.data$gzip_status %in% c("pass", "not_applicable") |
        .data$sha256_status != "pass"
    ) %>%
    pull(.data$file)
  integrity_status <- case_when(
    any(manifest$download_status != "present_not_active") ~ "pending_download",
    length(integrity_failures) == 0 ~ "validated",
    TRUE ~ "integrity_failed"
  )
  integrity_note <- ifelse(length(integrity_failures) == 0, "", paste(integrity_failures, collapse = ";"))
  tribble(
    ~feasibility_component, ~required_tables, ~status, ~notes,
    "download_integrity", "all SHA256SUMS files", integrity_status, integrity_note,
    "patient_admission_linkage", "hosp/patients.csv.gz;hosp/admissions.csv.gz", ifelse(schema_pass("hosp/patients.csv.gz") && schema_pass("hosp/admissions.csv.gz"), "ready", "pending"), "Needed for subject_id/hadm_id anchor and mortality context.",
    "icu_stay_denominator", "icu/icustays.csv.gz", ifelse(schema_pass("icu/icustays.csv.gz"), "ready", "pending"), "Needed for stay_id, ICU intime/outtime and ICU LOS.",
    "ward_icu_timing", "hosp/transfers.csv.gz;icu/icustays.csv.gz", ifelse(schema_pass("hosp/transfers.csv.gz") && schema_pass("icu/icustays.csv.gz"), "ready", "pending"), "Needed to check ICU location and transfer timing.",
    "respiratory_charting_items", "icu/d_items.csv.gz;icu/chartevents.csv.gz", ifelse(schema_pass("icu/d_items.csv.gz") && schema_pass("icu/chartevents.csv.gz") && nrow(candidates) > 0, "candidate_items_ready", "pending"), "D_items respiratory candidate labels are only metadata-level item candidates, not an endpoint definition.",
    "respiratory_support_procedures", "icu/procedureevents.csv.gz;hosp/procedures_icd.csv.gz;hosp/d_icd_procedures.csv.gz", ifelse(schema_pass("icu/procedureevents.csv.gz") && schema_pass("hosp/procedures_icd.csv.gz") && schema_pass("hosp/d_icd_procedures.csv.gz"), "ready_for_code_lock", "pending"), "Needed for ventilatory support procedure-code and procedure-event lock.",
    "oxygen_ventilation_endpoint", "icu/d_items.csv.gz;icu/chartevents.csv.gz;icu/procedureevents.csv.gz", ifelse(has_pass("icu/d_items.csv.gz") && has_pass("icu/chartevents.csv.gz") && has_pass("icu/procedureevents.csv.gz"), "ready_for_endpoint_lock", "pending"), "Endpoint lock must happen before any modeling.",
    "no_modeling_gate", "schema lock outputs", "enforced", "This script performs no cohort extraction and no predictive/association modeling."
  )
}

write_log <- function(root_out, manifest, schema_check, candidates, endpoint) {
  pending <- manifest %>% filter(.data$download_status != "present_not_active") %>% pull(.data$file)
  lines <- c(
    "# MIMIC-IV Integrity And Schema Feasibility V0.1",
    "",
    paste0("- Run date: ", Sys.Date()),
    "- Scope: file-level integrity, schema header checks, and first-pass respiratory-support feasibility only.",
    "- Data governance: no row-level clinical exports; outputs are aggregate/file/schema metadata only.",
    paste0("- Files listed in SHA256SUMS: ", nrow(manifest), "."),
    paste0("- Pending or active files: ", ifelse(length(pending) == 0, "none", paste(pending, collapse = "; ")), "."),
    paste0("- Schema expectations checked: ", nrow(schema_check), "."),
    paste0("- Respiratory d_items candidates exported: ", nrow(candidates), "."),
    paste0("- Endpoint feasibility rows: ", nrow(endpoint), "."),
    "- If pending files remain, do not start MIMIC modeling.",
    "",
    ifelse(length(pending) == 0, "STATUS: ready_for_full_gzip_sha_validation_or_schema_review.", "STATUS: pending_download.")
  )
  writeLines(lines, file.path(root_out, "results", "logs", "mimiciv_integrity_schema_feasibility_v0_1.md"))
}

main <- function() {
  project_root <- find_project_root()
  root <- mimic_root()
  dir.create(file.path(project_root, "results", "tables"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(project_root, "results", "logs"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(project_root, "metadata"), recursive = TRUE, showWarnings = FALSE)

  run_sha <- Sys.getenv("MIMIC_SHA_MODE", unset = "auto")
  if (!run_sha %in% c("auto", "always", "never")) stop("MIMIC_SHA_MODE must be auto, always, or never.", call. = FALSE)

  expectations <- schema_expectations()
  readr::write_csv(expectations, file.path(project_root, "metadata", "mimiciv_schema_lock_targets_v0_1.csv"))

  manifest <- build_manifest(root, run_sha = run_sha)
  columns <- build_schema_columns(root, manifest)
  schema_check <- build_schema_expectation_check(expectations, columns)
  candidates <- read_d_items_candidates(root, manifest)
  endpoint <- build_endpoint_feasibility(manifest, schema_check, candidates)

  readr::write_csv(manifest, file.path(project_root, "results", "tables", "mimiciv_v0_1_download_manifest.csv"))
  readr::write_csv(columns, file.path(project_root, "results", "tables", "mimiciv_v0_1_schema_columns.csv"))
  readr::write_csv(schema_check, file.path(project_root, "results", "tables", "mimiciv_v0_1_schema_expectations.csv"))
  readr::write_csv(candidates, file.path(project_root, "results", "tables", "mimiciv_v0_1_d_items_respiratory_candidates.csv"))
  readr::write_csv(endpoint, file.path(project_root, "results", "tables", "mimiciv_v0_1_endpoint_feasibility.csv"))
  write_log(project_root, manifest, schema_check, candidates, endpoint)
  message("Wrote MIMIC-IV integrity/schema feasibility V0.1 outputs.")
}

if (sys.nframe() == 0) {
  main()
}
