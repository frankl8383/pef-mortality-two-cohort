# Validate ICU translation V30 decision-threshold and recalibration outputs.

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
  tibble(
    check = paste0("file_exists_nonempty:", label),
    status = ifelse(exists && size > 0, "pass", "fail"),
    observed = ifelse(exists, as.character(size), "missing"),
    expected = ">0 bytes"
  )
}

check_scalar <- function(check, observed, expected, passed) {
  tibble(check = check, status = ifelse(passed, "pass", "fail"), observed = as.character(observed), expected = expected)
}

scan_for_absolute_user_paths <- function(paths) {
  hits <- character()
  for (p in paths[file.exists(paths)]) {
    if (!grepl("\\.(csv|md|txt)$", p, ignore.case = TRUE)) {
      next
    }
    txt <- paste(readLines(p, warn = FALSE), collapse = "\n")
    if (grepl("${HOME}", txt, fixed = TRUE)) {
      hits <- c(hits, basename(p))
    }
  }
  hits
}

root <- find_project_root()
derived_dir <- file.path(root, "derived_sensitive", "icu_translation")
results_dir <- file.path(root, "results", "tables")
logs_dir <- file.path(root, "results", "logs")
figures_dir <- file.path(root, "results", "figures")
manuscript_dir <- file.path(root, "manuscript")

paths <- list(
  predictions = file.path(derived_dir, paste0("icu_translation_predictions_", version_id, ".rds")),
  preprocess = file.path(results_dir, paste0("icu_dynamic_prediction_preprocessing_", version_id, ".csv")),
  reconstruction = file.path(results_dir, paste0("icu_prediction_reconstruction_audit_", version_id, ".csv")),
  threshold = file.path(results_dir, paste0("icu_decision_threshold_table_", version_id, ".csv")),
  delta = file.path(results_dir, paste0("icu_decision_threshold_delta_", version_id, ".csv")),
  dca = file.path(results_dir, paste0("icu_decision_curve_reproduction_", version_id, ".csv")),
  recalibration = file.path(results_dir, paste0("icu_recalibration_transportability_", version_id, ".csv")),
  transport = file.path(results_dir, paste0("icu_transportability_summary_", version_id, ".csv")),
  note = file.path(manuscript_dir, paste0("icu_decision_recalibration_interpretation_", version_id, ".md")),
  log = file.path(logs_dir, paste0("icu_decision_recalibration_", version_id, ".md")),
  fig_dca_png = file.path(figures_dir, paste0("icu_decision_threshold_net_benefit_", version_id, ".png")),
  fig_dca_pdf = file.path(figures_dir, paste0("icu_decision_threshold_net_benefit_", version_id, ".pdf")),
  fig_cal_png = file.path(figures_dir, paste0("icu_recalibration_plot_", version_id, ".png")),
  fig_cal_pdf = file.path(figures_dir, paste0("icu_recalibration_plot_", version_id, ".pdf")),
  fig_transport_png = file.path(figures_dir, paste0("icu_transportability_subgroups_", version_id, ".png")),
  fig_transport_pdf = file.path(figures_dir, paste0("icu_transportability_subgroups_", version_id, ".pdf"))
)

checks <- bind_rows(lapply(names(paths), function(nm) check_file(paths[[nm]], nm)))

reconstruction <- readr::read_csv(paths$reconstruction, show_col_types = FALSE)
threshold <- readr::read_csv(paths$threshold, show_col_types = FALSE)
delta <- readr::read_csv(paths$delta, show_col_types = FALSE)
dca <- readr::read_csv(paths$dca, show_col_types = FALSE)
recalibration <- readr::read_csv(paths$recalibration, show_col_types = FALSE)
transport <- readr::read_csv(paths$transport, show_col_types = FALSE)

dynamic_recon <- reconstruction %>% filter(model == "dynamic_last_no_support", evaluation %in% c("mimic_internal_test", "eicu_external"))
max_auc_delta <- max(abs(dynamic_recon$auroc_delta_vs_v24), na.rm = TRUE)
max_brier_delta <- max(abs(dynamic_recon$brier_delta_vs_v24), na.rm = TRUE)
max_dca_delta <- max(abs(dca$net_benefit_delta_vs_v24), na.rm = TRUE)
eicu_delta_020 <- delta %>% filter(evaluation == "eicu_external", threshold == 0.2)
eicu_recal <- recalibration %>% filter(evaluation == "eicu_external", model == "dynamic_last_no_support", recalibration_method == "original_mimic_trained")
hf_eicu <- transport %>% filter(evaluation == "eicu_external", model == "dynamic_last_no_support", group_family == "support_class", subgroup == "hfnc")
public_paths <- unlist(paths[setdiff(names(paths), "predictions")])
absolute_path_hits <- scan_for_absolute_user_paths(public_paths)

checks <- bind_rows(
  checks,
  check_scalar("dynamic_reconstruction_auroc_matches_v24", max_auc_delta, "<1e-4", max_auc_delta < 1e-4),
  check_scalar("dynamic_reconstruction_brier_matches_v24", max_brier_delta, "<1e-4", max_brier_delta < 1e-4),
  check_scalar("dca_reproduction_matches_v24", max_dca_delta, "<2e-4", max_dca_delta < 2e-4),
  check_scalar("threshold_grid_complete", paste(sort(unique(threshold$threshold)), collapse = ","), "0.05 to 0.50", length(unique(threshold$threshold)) == 10),
  check_scalar("eicu_delta_020_present", nrow(eicu_delta_020), "1", nrow(eicu_delta_020) == 1),
  check_scalar("eicu_dynamic_net_benefit_above_rox_at_020", ifelse(nrow(eicu_delta_020) == 1, eicu_delta_020$net_benefit_delta_dynamic_minus_rox, NA), ">0", nrow(eicu_delta_020) == 1 && eicu_delta_020$net_benefit_delta_dynamic_minus_rox > 0),
  check_scalar("eicu_dynamic_calibration_gap_recorded", ifelse(nrow(eicu_recal) == 1, eicu_recal$calibration_gap_observed_minus_predicted, NA), "not missing", nrow(eicu_recal) == 1 && !is.na(eicu_recal$calibration_gap_observed_minus_predicted)),
  check_scalar("eicu_hfnc_boundary_recorded", ifelse(nrow(hf_eicu) == 1, hf_eicu$safe_interpretation, NA), "descriptive only", nrow(hf_eicu) == 1 && grepl("Descriptive only", hf_eicu$safe_interpretation)),
  check_scalar("no_absolute_user_paths_in_public_outputs", paste(absolute_path_hits, collapse = ","), "none", length(absolute_path_hits) == 0)
)

all_passed <- all(checks$status == "pass")
validation_csv <- file.path(results_dir, paste0("icu_decision_recalibration_validation_", version_id, ".csv"))
readr::write_csv(checks, validation_csv)

if (all_passed) {
  ticket_path <- file.path(results_dir, "v30_high_impact_ticket_audit.csv")
  if (file.exists(ticket_path)) {
    ticket_audit <- readr::read_csv(ticket_path, show_col_types = FALSE)
    ticket_audit <- ticket_audit %>%
      mutate(status = dplyr::case_when(
        ticket == "T3_ICU_decision_threshold_analysis" ~ "completed",
        ticket == "T4_ICU_recalibration_transportability" ~ "completed",
        TRUE ~ status
      ))
    readr::write_csv(ticket_audit, ticket_path)
  }
}

validation_log <- file.path(logs_dir, paste0("icu_decision_recalibration_validation_", version_id, ".md"))
log_lines <- c(
  "# ICU V30 Decision/Recalibration Validation",
  "",
  paste0("Version: ", version_id),
  "",
  paste0("Overall status: ", ifelse(all_passed, "PASS", "FAIL")),
  "",
  "## Key Checks",
  "",
  paste0("- Max dynamic AUROC reconstruction delta versus V24: ", signif(max_auc_delta, 4)),
  paste0("- Max dynamic Brier reconstruction delta versus V24: ", signif(max_brier_delta, 4)),
  paste0("- Max DCA net-benefit reproduction delta versus V24: ", signif(max_dca_delta, 4)),
  "",
  "## Checks",
  "",
  paste0("| ", paste(names(checks), collapse = " | "), " |"),
  paste0("| ", paste(rep("---", ncol(checks)), collapse = " | "), " |"),
  apply(checks, 1, function(x) paste0("| ", paste(x, collapse = " | "), " |"))
)
writeLines(log_lines, validation_log)

message("ICU decision/recalibration validation status: ", ifelse(all_passed, "PASS", "FAIL"))
if (!all_passed) {
  print(checks %>% filter(status == "fail"))
  stop("Validation failed.", call. = FALSE)
}
