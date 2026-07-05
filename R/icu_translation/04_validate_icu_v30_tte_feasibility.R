# Validate ICU V30 TTE feasibility-only outputs.

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

check_file <- function(path, label) {
  exists <- file.exists(path)
  size <- if (exists) file.info(path)$size else NA_real_
  tibble(check = paste0("file_exists_nonempty:", label), status = ifelse(exists && size > 0, "pass", "fail"), observed = ifelse(exists, as.character(size), "missing"), expected = ">0 bytes")
}

check_scalar <- function(check, observed, expected, passed) {
  tibble(check = check, status = ifelse(passed, "pass", "fail"), observed = as.character(observed), expected = expected)
}

root <- find_project_root()
results_dir <- file.path(root, "results", "tables")
logs_dir <- file.path(root, "results", "logs")
manuscript_dir <- file.path(root, "manuscript")

counts_path <- file.path(results_dir, paste0("tte_feasibility_counts_", version_id, ".csv"))
weights_path <- file.path(results_dir, paste0("tte_weight_diagnostics_", version_id, ".csv"))
balance_path <- file.path(results_dir, paste0("tte_balance_preliminary_", version_id, ".csv"))
balance_summary_path <- file.path(results_dir, paste0("tte_balance_summary_", version_id, ".csv"))
protocol_path <- file.path(results_dir, paste0("tte_protocol_table_", version_id, ".csv"))
go_no_go_path <- file.path(results_dir, paste0("tte_go_no_go_", version_id, ".csv"))
note_path <- file.path(manuscript_dir, paste0("tte_feasibility_note_", version_id, ".md"))
log_path <- file.path(logs_dir, paste0("tte_feasibility_", version_id, ".md"))
validation_csv <- file.path(results_dir, paste0("tte_feasibility_validation_", version_id, ".csv"))
validation_log <- file.path(logs_dir, paste0("tte_feasibility_validation_", version_id, ".md"))

checks <- bind_rows(
  check_file(counts_path, "counts"),
  check_file(weights_path, "weights"),
  check_file(balance_path, "balance"),
  check_file(balance_summary_path, "balance_summary"),
  check_file(protocol_path, "protocol"),
  check_file(go_no_go_path, "go_no_go"),
  check_file(note_path, "note"),
  check_file(log_path, "log")
)

counts <- readr::read_csv(counts_path, show_col_types = FALSE)
weights <- readr::read_csv(weights_path, show_col_types = FALSE)
balance_summary <- readr::read_csv(balance_summary_path, show_col_types = FALSE)
go_no_go <- readr::read_csv(go_no_go_path, show_col_types = FALSE)
note <- paste(readLines(note_path, warn = FALSE), collapse = "\n")

mimic_overall <- counts %>% filter(dataset == "MIMIC-IV", group == "overall")
eicu_overall <- counts %>% filter(dataset == "eICU-CRD", group == "overall")
effect_language <- grepl("\\b(HR|OR)\\s*[=:]?\\s*[0-9]|risk difference\\s*[=:]?\\s*[-0-9]|\\bATE\\b|causal estimate\\s*[=:]?\\s*[-0-9]", note, ignore.case = TRUE)

checks <- bind_rows(
  checks,
  check_scalar("mimic_hfnc_sample", mimic_overall$hfnc_n, ">=500", nrow(mimic_overall) == 1 && mimic_overall$hfnc_n >= 500),
  check_scalar("eicu_hfnc_sample", eicu_overall$hfnc_n, ">=100 but sparse flag needed", nrow(eicu_overall) == 1 && eicu_overall$hfnc_n >= 100),
  check_scalar("weight_diagnostics_two_datasets", nrow(weights), "2", nrow(weights) == 2),
  check_scalar("balance_summary_two_datasets", nrow(balance_summary), "2", nrow(balance_summary) == 2),
  check_scalar("go_no_go_two_datasets", nrow(go_no_go), "2", nrow(go_no_go) == 2),
  check_scalar("no_dataset_approved_for_effect_estimation", paste(go_no_go$final_decision, collapse = "; "), "all no_go or protocol_review", all(grepl("no_go|protocol_review", go_no_go$final_decision))),
  check_scalar("note_has_no_effect_estimate_language", effect_language, "FALSE", !effect_language),
  check_scalar("note_states_feasibility_only", grepl("feasibility audit only", note, fixed = TRUE), "TRUE", grepl("feasibility audit only", note, fixed = TRUE))
)

all_passed <- all(checks$status == "pass")
readr::write_csv(checks, validation_csv)

if (all_passed) {
  ticket_path <- file.path(results_dir, "v30_high_impact_ticket_audit.csv")
  if (file.exists(ticket_path)) {
    ticket_audit <- readr::read_csv(ticket_path, show_col_types = FALSE)
    ticket_audit <- ticket_audit %>%
      mutate(status = ifelse(ticket == "T5_TTE_feasibility_only", "completed", status))
    readr::write_csv(ticket_audit, ticket_path)
  }
}

log_lines <- c(
  "# ICU TTE Feasibility Validation",
  "",
  paste0("Version: ", version_id),
  "",
  paste0("Overall status: ", ifelse(all_passed, "PASS", "FAIL")),
  "",
  "## Checks",
  "",
  paste0("| ", paste(names(checks), collapse = " | "), " |"),
  paste0("| ", paste(rep("---", ncol(checks)), collapse = " | "), " |"),
  apply(checks, 1, function(x) paste0("| ", paste(x, collapse = " | "), " |"))
)
writeLines(log_lines, validation_log)

message("ICU TTE feasibility validation status: ", ifelse(all_passed, "PASS", "FAIL"))
if (!all_passed) {
  print(checks %>% filter(status == "fail"))
  stop("Validation failed.", call. = FALSE)
}
