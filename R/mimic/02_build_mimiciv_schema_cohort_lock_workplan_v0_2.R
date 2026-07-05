# Build a MIMIC-IV schema/cohort/endpoint lock workplan for the first
# post-download step. This is a planning artifact only, not data extraction.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
  library(stringr)
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

write_doc <- function(path, lines) {
  writeLines(lines, path)
}

main <- function() {
  root <- find_project_root()
  dirs <- file.path(root, c("metadata", "results/tables", "results/logs"))
  invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))

  workplan <- tribble(
    ~gate_order, ~gate, ~objective, ~required_inputs, ~allowed_outputs, ~not_allowed,
    0L, "download_integrity_gate", "Confirm all files listed in SHA256SUMS are present and pass gzip/SHA256 checks.", "MIMIC v0.1 download manifest; SHA256SUMS", "aggregate file manifest; pass/fail log", "modeling; row-level export; manual deletion of partial files",
    1L, "schema_header_gate", "Confirm required columns for patients, admissions, icustays, d_items, chartevents, procedureevents, inputevents, outputevents and code dictionaries.", "metadata/mimiciv_schema_lock_targets_v0_1.csv; header-only reads", "schema column table; missing-column table", "reading full event tables before integrity pass",
    2L, "linkage_key_gate", "Check availability of subject_id, hadm_id and stay_id linkage keys and timing fields.", "patients; admissions; icustays; transfers headers and aggregate counts only", "aggregate linkage-feasibility table", "patient-level listing; row-level joins exported to disk",
    3L, "respiratory_item_gate", "Lock candidate d_items rows for oxygen, ventilation, airway, respiratory rate and device charting.", "icu/d_items.csv.gz after integrity pass", "candidate item dictionary with itemid/label/category only", "using item candidates as final endpoints without manual review",
    4L, "respiratory_support_endpoint_gate", "Separate candidate endpoints into invasive ventilation, non-invasive ventilation, high-flow oxygen, conventional oxygen and airway/procedure context.", "d_items candidates; procedureevents headers; ICD procedure dictionaries", "endpoint candidate ledger; ambiguity flags", "claiming acute respiratory failure outcome without endpoint lock",
    5L, "cohort_denominator_gate", "Define candidate ICU denominator before analysis: adult ICU stays, first ICU stay policy, admission timing and minimum data availability.", "patients/admissions/icustays aggregate counts; no row-level export", "aggregate cohort-count ladder", "predictive/association modeling",
    6L, "manuscript_integration_gate", "Decide whether MIMIC can become an active evidence layer after schema/cohort/endpoint lock.", "all prior gate logs", "go/no-go decision note", "adding MIMIC claims to manuscript before gates pass"
  )

  readr::write_csv(workplan, file.path(root, "metadata", "mimiciv_schema_cohort_lock_workplan_v0_2.csv"))

  lines <- c(
    "# MIMIC-IV Schema/Cohort/Endpoint Lock Workplan V0.2",
    "",
    "This workplan defines the first post-download MIMIC step. It is deliberately limited to integrity, schema, linkage, cohort-denominator and endpoint-lock feasibility. It does not permit MIMIC modeling.",
    "",
    "| Gate | Objective | Required inputs | Allowed outputs | Not allowed |",
    "| --- | --- | --- | --- | --- |",
    apply(workplan, 1, function(x) paste0("| ", x[["gate"]], " | ", str_replace_all(x[["objective"]], "\\|", "\\\\|"), " | ", str_replace_all(x[["required_inputs"]], "\\|", "\\\\|"), " | ", str_replace_all(x[["allowed_outputs"]], "\\|", "\\\\|"), " | ", str_replace_all(x[["not_allowed"]], "\\|", "\\\\|"), " |")),
    "",
    "## Immediate Rule After Download",
    "",
    "Run only:",
    "",
    "```bash",
    "Rscript R/mimic/00_mimiciv_integrity_schema_feasibility_v0_1.R",
    "Rscript R/mimic/01_validate_mimiciv_integrity_schema_feasibility_v0_1.R",
    "```",
    "",
    "Proceed to cohort or endpoint lock only if v0.1 validation passes outside pending-download mode.",
    "",
    "## Evidence Boundary",
    "",
    "MIMIC-IV remains outside the active manuscript evidence claim until integrity, schema, linkage, cohort-denominator and endpoint-candidate gates are complete."
  )
  write_doc(file.path(root, "results", "logs", "mimiciv_schema_cohort_lock_workplan_v0_2.md"), lines)

  log <- c(
    "# MIMIC-IV Schema/Cohort/Endpoint Lock Workplan V0.2",
    "",
    paste0("- Run date: ", Sys.Date()),
    paste0("- Gate rows: ", nrow(workplan), "."),
    "- This is a planning artifact only.",
    "- No MIMIC row-level data were read or exported.",
    "- No modeling is permitted by this workplan.",
    "",
    "PASS: MIMIC-IV schema/cohort lock workplan V0.2 build complete."
  )
  write_doc(file.path(root, "results", "logs", "mimiciv_schema_cohort_lock_workplan_v0_2_build.md"), log)
  message("PASS: MIMIC-IV schema/cohort lock workplan V0.2 built.")
}

if (sys.nframe() == 0) {
  main()
}
