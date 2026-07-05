# ICU translation V30 decision-threshold and recalibration outputs.
# Uses existing local-only V22 modeling rows plus validated V23/V24 aggregate
# coefficients/diagnostics. This script reconstructs the V23 best model
# predictions exactly from saved coefficients and MIMIC train preprocessing.

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(tibble)
  library(tidyr)
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

format_ci <- function(est, lo, hi, digits = 3) {
  sprintf(paste0("%.", digits, "f (%.", digits, "f-%.", digits, "f)"), est, lo, hi)
}

markdown_table <- function(dat) {
  dat <- as.data.frame(dat, stringsAsFactors = FALSE)
  dat[] <- lapply(dat, function(x) ifelse(is.na(x), "", as.character(x)))
  header <- paste0("| ", paste(names(dat), collapse = " | "), " |")
  separator <- paste0("| ", paste(rep("---", ncol(dat)), collapse = " | "), " |")
  rows <- apply(dat, 1, function(x) paste0("| ", paste(x, collapse = " | "), " |"))
  c(header, separator, rows)
}

auroc_rank <- function(y, pred) {
  ok <- !is.na(y) & !is.na(pred)
  y <- y[ok]
  pred <- pred[ok]
  n1 <- sum(y == 1L)
  n0 <- sum(y == 0L)
  if (n1 == 0 || n0 == 0) {
    return(NA_real_)
  }
  ranks <- rank(pred, ties.method = "average")
  (sum(ranks[y == 1L]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

reconstruct_dynamic_predictions <- function(mimic, eicu, coef_table) {
  coef <- coef_table %>% filter(model == "dynamic_last_no_support")
  value_terms <- coef %>%
    filter(term != "intercept", !grepl("_missing$", term)) %>%
    pull(term)
  train <- mimic %>% filter(split == "train")
  medians <- vapply(value_terms, function(v) stats::median(train[[v]], na.rm = TRUE), numeric(1))
  train_imp <- train
  for (v in value_terms) {
    train_imp[[v]][is.na(train_imp[[v]])] <- medians[[v]]
  }
  centers <- vapply(value_terms, function(v) mean(train_imp[[v]], na.rm = TRUE), numeric(1))
  scales <- vapply(value_terms, function(v) stats::sd(train_imp[[v]], na.rm = TRUE), numeric(1))
  scales[is.na(scales) | scales == 0] <- 1

  predict_one <- function(df) {
    lp <- rep(coef$coefficient[coef$term == "intercept"], nrow(df))
    for (i in seq_len(nrow(coef))) {
      term <- coef$term[[i]]
      beta <- coef$coefficient[[i]]
      if (term == "intercept") {
        next
      }
      if (grepl("_missing$", term)) {
        base <- sub("_missing$", "", term)
        lp <- lp + beta * as.numeric(is.na(df[[base]]))
      } else {
        x <- df[[term]]
        x[is.na(x)] <- medians[[term]]
        if (term != "male") {
          x <- (x - centers[[term]]) / scales[[term]]
        }
        lp <- lp + beta * x
      }
    }
    stats::plogis(lp)
  }

  list(
    mimic_pred = predict_one(mimic),
    eicu_pred = predict_one(eicu),
    preprocessing = tibble(
      feature = value_terms,
      imputation_median_mimic_train = as.numeric(medians[value_terms]),
      center_after_imputation_mimic_train = as.numeric(centers[value_terms]),
      scale_after_imputation_mimic_train = as.numeric(scales[value_terms]),
      standardized = value_terms != "male"
    )
  )
}

performance_row <- function(data, pred_col, model, evaluation, source) {
  df <- data %>% filter(!is.na(.data[[pred_col]]))
  y <- df$outcome_failure_48h_l2
  pred <- df[[pred_col]]
  tibble(
    model = model,
    evaluation = evaluation,
    source = source,
    n = nrow(df),
    events = sum(y == 1L, na.rm = TRUE),
    event_rate = mean(y == 1L, na.rm = TRUE),
    mean_predicted_risk = mean(pred, na.rm = TRUE),
    auroc = auroc_rank(y, pred),
    brier = mean((y - pred)^2, na.rm = TRUE)
  )
}

decision_metrics <- function(data, pred_col, model, evaluation, analysis_population, threshold) {
  df <- data %>% filter(!is.na(.data[[pred_col]]))
  y <- df$outcome_failure_48h_l2
  pred_pos <- df[[pred_col]] >= threshold
  tp <- sum(pred_pos & y == 1L, na.rm = TRUE)
  fp <- sum(pred_pos & y == 0L, na.rm = TRUE)
  fn <- sum(!pred_pos & y == 1L, na.rm = TRUE)
  tn <- sum(!pred_pos & y == 0L, na.rm = TRUE)
  n <- nrow(df)
  events <- sum(y == 1L, na.rm = TRUE)
  nonevents <- sum(y == 0L, na.rm = TRUE)
  tibble(
    evaluation = evaluation,
    analysis_population = analysis_population,
    model = model,
    threshold = threshold,
    n = n,
    events = events,
    event_rate = events / n,
    predicted_positive = tp + fp,
    predicted_positive_rate = (tp + fp) / n,
    events_captured = tp,
    event_capture_rate = ifelse(events > 0, tp / events, NA_real_),
    false_positives = fp,
    false_positive_rate = ifelse(nonevents > 0, fp / nonevents, NA_real_),
    ppv = ifelse(tp + fp > 0, tp / (tp + fp), NA_real_),
    npv = ifelse(tn + fn > 0, tn / (tn + fn), NA_real_),
    sensitivity = ifelse(events > 0, tp / events, NA_real_),
    specificity = ifelse(nonevents > 0, tn / nonevents, NA_real_),
    net_benefit = tp / n - fp / n * (threshold / (1 - threshold))
  )
}

root <- find_project_root()
derived_dir <- file.path(root, "derived_sensitive", "icu_translation")
results_dir <- file.path(root, "results", "tables")
logs_dir <- file.path(root, "results", "logs")
figures_dir <- file.path(root, "results", "figures")
manuscript_dir <- file.path(root, "manuscript")
dir_create(derived_dir)
dir_create(results_dir)
dir_create(logs_dir)
dir_create(figures_dir)
dir_create(manuscript_dir)

mimic <- readr::read_csv(file.path(derived_dir, "mimiciv_icu_translation_modeling_v22_0.csv.gz"), show_col_types = FALSE)
eicu <- readr::read_csv(file.path(derived_dir, "eicu_icu_translation_modeling_v22_0.csv.gz"), show_col_types = FALSE)
coef_v23 <- readr::read_csv(file.path(results_dir, "icu_translation_model_variant_coefficients_v23_0.csv"), show_col_types = FALSE)
perf_v24 <- readr::read_csv(file.path(results_dir, "icu_translation_model_performance_v24_0.csv"), show_col_types = FALSE)
dca_v24 <- readr::read_csv(file.path(results_dir, "icu_translation_dca_v24_0.csv"), show_col_types = FALSE)
recal_v24 <- readr::read_csv(file.path(results_dir, "icu_translation_recalibration_v24_0.csv"), show_col_types = FALSE)
bins_v24 <- readr::read_csv(file.path(results_dir, "icu_translation_calibration_bins_v24_0.csv"), show_col_types = FALSE)
subgroups_v24 <- readr::read_csv(file.path(results_dir, "icu_translation_subgroup_transport_v24_0.csv"), show_col_types = FALSE)

preds <- reconstruct_dynamic_predictions(mimic, eicu, coef_v23)
mimic$pred_dynamic_last_no_support_v30 <- preds$mimic_pred
eicu$pred_dynamic_last_no_support_v30 <- preds$eicu_pred

icu_predictions <- bind_rows(
  mimic %>% mutate(evaluation = ifelse(split == "test", "mimic_internal_test", "mimic_train_apparent")),
  eicu %>% mutate(evaluation = "eicu_external")
)
saveRDS(icu_predictions, file.path(derived_dir, paste0("icu_translation_predictions_", version_id, ".rds")))

preprocess_path <- file.path(results_dir, paste0("icu_dynamic_prediction_preprocessing_", version_id, ".csv"))
readr::write_csv(preds$preprocessing, preprocess_path)

reconstruction <- bind_rows(
  performance_row(mimic %>% filter(split == "test"), "pred_dynamic_last_no_support_v30", "dynamic_last_no_support", "mimic_internal_test", "reconstructed_v30"),
  performance_row(eicu, "pred_dynamic_last_no_support_v30", "dynamic_last_no_support", "eicu_external", "reconstructed_v30"),
  performance_row(mimic %>% filter(split == "test", !is.na(pred_rox_only)), "pred_rox_only", "rox_only", "mimic_internal_test", "existing_row_level"),
  performance_row(eicu %>% filter(!is.na(pred_rox_only)), "pred_rox_only", "rox_only", "eicu_external", "existing_row_level")
) %>%
  left_join(
    perf_v24 %>%
      filter(evaluation %in% c("mimic_internal_test", "eicu_external"), model %in% c("dynamic_last_no_support", "rox_only")) %>%
      select(model, evaluation, v24_n = n, v24_events = events, v24_mean_predicted_risk = mean_predicted_risk, v24_auroc = auroc, v24_brier = brier),
    by = c("model", "evaluation")
  ) %>%
  mutate(
    n_delta_vs_v24 = n - v24_n,
    events_delta_vs_v24 = events - v24_events,
    mean_predicted_risk_delta_vs_v24 = mean_predicted_risk - v24_mean_predicted_risk,
    auroc_delta_vs_v24 = auroc - v24_auroc,
    brier_delta_vs_v24 = brier - v24_brier
  )
readr::write_csv(reconstruction, file.path(results_dir, paste0("icu_prediction_reconstruction_audit_", version_id, ".csv")))

thresholds <- seq(0.05, 0.50, by = 0.05)
mimic_test <- mimic %>% filter(split == "test")
eicu_all <- eicu
mimic_common <- mimic_test %>% filter(!is.na(pred_rox_only))
eicu_common <- eicu_all %>% filter(!is.na(pred_rox_only))

threshold_table_raw <- bind_rows(lapply(thresholds, function(th) {
  bind_rows(
    decision_metrics(mimic_test, "pred_dynamic_last_no_support_v30", "dynamic_last_no_support", "mimic_internal_test", "full_dynamic_available", th),
    decision_metrics(eicu_all, "pred_dynamic_last_no_support_v30", "dynamic_last_no_support", "eicu_external", "full_dynamic_available", th),
    decision_metrics(mimic_common, "pred_dynamic_last_no_support_v30", "dynamic_last_no_support", "mimic_internal_test", "common_rox_complete_case", th),
    decision_metrics(eicu_common, "pred_dynamic_last_no_support_v30", "dynamic_last_no_support", "eicu_external", "common_rox_complete_case", th),
    decision_metrics(mimic_common, "pred_rox_only", "rox_only", "mimic_internal_test", "common_rox_complete_case", th),
    decision_metrics(eicu_common, "pred_rox_only", "rox_only", "eicu_external", "common_rox_complete_case", th)
  )
}))

threshold_table <- threshold_table_raw %>%
  mutate(across(where(is.numeric), ~ round(.x, 6)))

delta_table_raw <- threshold_table_raw %>%
  filter(analysis_population == "common_rox_complete_case", model %in% c("dynamic_last_no_support", "rox_only")) %>%
  select(evaluation, threshold, model, n, events, predicted_positive, predicted_positive_rate,
         events_captured, event_capture_rate, false_positives, ppv, sensitivity, specificity, net_benefit) %>%
  pivot_wider(
    names_from = model,
    values_from = c(predicted_positive, predicted_positive_rate, events_captured, event_capture_rate,
                    false_positives, ppv, sensitivity, specificity, net_benefit),
    names_sep = "__"
  ) %>%
  mutate(
    predicted_positive_delta_dynamic_minus_rox = predicted_positive__dynamic_last_no_support - predicted_positive__rox_only,
    events_captured_delta_dynamic_minus_rox = events_captured__dynamic_last_no_support - events_captured__rox_only,
    false_positives_delta_dynamic_minus_rox = false_positives__dynamic_last_no_support - false_positives__rox_only,
    net_benefit_delta_dynamic_minus_rox = net_benefit__dynamic_last_no_support - net_benefit__rox_only
  )

delta_table <- delta_table_raw %>%
  mutate(across(where(is.numeric), ~ round(.x, 6)))

dca_reproduction <- threshold_table_raw %>%
  filter(analysis_population == "common_rox_complete_case", model %in% c("dynamic_last_no_support", "rox_only")) %>%
  select(evaluation, analysis_population, model, threshold, computed_net_benefit = net_benefit) %>%
  left_join(
    dca_v24 %>%
      filter(strategy == "model", model %in% c("dynamic_last_no_support", "rox_only")) %>%
      select(evaluation, analysis_population, model, threshold, v24_net_benefit = net_benefit),
    by = c("evaluation", "analysis_population", "model", "threshold")
  ) %>%
  mutate(net_benefit_delta_vs_v24 = computed_net_benefit - v24_net_benefit)

recalibration <- recal_v24 %>%
  group_by(evaluation, model) %>%
  mutate(
    original_brier = brier[recalibration_method == "original_mimic_trained"][1],
    original_mean_predicted_risk = mean_predicted_risk[recalibration_method == "original_mimic_trained"][1],
    calibration_gap_observed_minus_predicted = event_rate - mean_predicted_risk,
    brier_delta_vs_original = brier - original_brier,
    diagnostic_only = recalibration_method != "original_mimic_trained"
  ) %>%
  ungroup() %>%
  mutate(across(where(is.numeric), ~ round(.x, 6)))

transport <- subgroups_v24 %>%
  mutate(
    safe_interpretation = dplyr::case_when(
      evaluation == "eicu_external" & group_family == "support_class" & subgroup == "hfnc" ~ "Descriptive only: eICU HFNC-coded rows are sparse relative to total eICU validation.",
      interpretation_flag == "interpret_cautiously_sparse" ~ "Interpret cautiously because subgroup is sparse.",
      TRUE ~ "Adequate for descriptive transportability check."
    ),
    claim_boundary = "Supportive ICU translation only; do not call this an actionable clinical decision tool."
  ) %>%
  mutate(across(where(is.numeric), ~ round(.x, 6)))

threshold_path <- file.path(results_dir, paste0("icu_decision_threshold_table_", version_id, ".csv"))
delta_path <- file.path(results_dir, paste0("icu_decision_threshold_delta_", version_id, ".csv"))
dca_repro_path <- file.path(results_dir, paste0("icu_decision_curve_reproduction_", version_id, ".csv"))
recalibration_path <- file.path(results_dir, paste0("icu_recalibration_transportability_", version_id, ".csv"))
transport_path <- file.path(results_dir, paste0("icu_transportability_summary_", version_id, ".csv"))
readr::write_csv(threshold_table, threshold_path)
readr::write_csv(delta_table, delta_path)
readr::write_csv(dca_reproduction, dca_repro_path)
readr::write_csv(recalibration, recalibration_path)
readr::write_csv(transport, transport_path)

plot_dca <- delta_table_raw %>%
  select(evaluation, threshold, net_benefit__dynamic_last_no_support, net_benefit__rox_only) %>%
  pivot_longer(cols = starts_with("net_benefit"), names_to = "model", values_to = "net_benefit") %>%
  mutate(
    model = recode(model, net_benefit__dynamic_last_no_support = "Dynamic model", net_benefit__rox_only = "ROX-like comparator"),
    evaluation = recode(evaluation, mimic_internal_test = "MIMIC holdout", eicu_external = "eICU external")
  )
p_dca <- ggplot(plot_dca, aes(x = threshold, y = net_benefit, colour = model)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "#6B7280") +
  geom_line(linewidth = 0.5) +
  geom_point(size = 1.2) +
  facet_wrap(~evaluation, ncol = 2) +
  scale_colour_manual(values = c("Dynamic model" = "#1F6F8B", "ROX-like comparator" = "#7C2D12")) +
  labs(x = "Risk threshold", y = "Net benefit", colour = NULL) +
  theme_classic(base_size = 8) +
  theme(legend.position = "bottom", axis.line = element_line(linewidth = 0.35), axis.ticks = element_line(linewidth = 0.35))
ggplot2::ggsave(file.path(figures_dir, paste0("icu_decision_threshold_net_benefit_", version_id, ".png")), p_dca, width = 6.6, height = 3.6, dpi = 600, units = "in")
ggplot2::ggsave(file.path(figures_dir, paste0("icu_decision_threshold_net_benefit_", version_id, ".pdf")), p_dca, width = 6.6, height = 3.6, units = "in", device = grDevices::pdf, useDingbats = FALSE)

plot_bins <- bins_v24 %>%
  filter(model %in% c("dynamic_last_no_support", "rox_only")) %>%
  mutate(
    model_label = recode(model, dynamic_last_no_support = "Dynamic model", rox_only = "ROX-like comparator"),
    evaluation_label = recode(evaluation, mimic_internal_test = "MIMIC holdout", eicu_external = "eICU external")
  )
p_cal <- ggplot(plot_bins, aes(x = mean_predicted_risk, y = observed_event_rate, colour = model_label)) +
  geom_abline(slope = 1, intercept = 0, linewidth = 0.3, colour = "#6B7280") +
  geom_line(linewidth = 0.45) +
  geom_point(aes(size = n), alpha = 0.85) +
  facet_wrap(~evaluation_label, ncol = 2) +
  scale_colour_manual(values = c("Dynamic model" = "#1F6F8B", "ROX-like comparator" = "#7C2D12")) +
  scale_size_continuous(range = c(1.0, 2.2), guide = "none") +
  labs(x = "Mean predicted risk", y = "Observed event rate", colour = NULL) +
  theme_classic(base_size = 8) +
  theme(legend.position = "bottom", axis.line = element_line(linewidth = 0.35), axis.ticks = element_line(linewidth = 0.35))
ggplot2::ggsave(file.path(figures_dir, paste0("icu_recalibration_plot_", version_id, ".png")), p_cal, width = 6.6, height = 3.6, dpi = 600, units = "in")
ggplot2::ggsave(file.path(figures_dir, paste0("icu_recalibration_plot_", version_id, ".pdf")), p_cal, width = 6.6, height = 3.6, units = "in", device = grDevices::pdf, useDingbats = FALSE)

plot_transport <- transport %>%
  filter(model == "dynamic_last_no_support", evaluation %in% c("mimic_internal_test", "eicu_external")) %>%
  mutate(
    label = paste(group_family, subgroup, sep = ": "),
    evaluation_label = recode(evaluation, mimic_internal_test = "MIMIC holdout", eicu_external = "eICU external")
  )
p_transport <- ggplot(plot_transport, aes(x = auroc, y = reorder(label, auroc))) +
  geom_vline(xintercept = 0.70, linewidth = 0.25, colour = "#9CA3AF") +
  geom_point(aes(size = n, colour = evaluation_label), alpha = 0.9) +
  facet_wrap(~evaluation_label, scales = "free_y", ncol = 2) +
  scale_colour_manual(values = c("MIMIC holdout" = "#1F6F8B", "eICU external" = "#7C2D12")) +
  scale_size_continuous(range = c(1.0, 3.0), guide = "none") +
  labs(x = "AUROC", y = NULL, colour = NULL) +
  theme_classic(base_size = 8) +
  theme(legend.position = "none", axis.line = element_line(linewidth = 0.35), axis.ticks = element_line(linewidth = 0.35))
ggplot2::ggsave(file.path(figures_dir, paste0("icu_transportability_subgroups_", version_id, ".png")), p_transport, width = 7.2, height = 4.8, dpi = 600, units = "in")
ggplot2::ggsave(file.path(figures_dir, paste0("icu_transportability_subgroups_", version_id, ".pdf")), p_transport, width = 7.2, height = 4.8, units = "in", device = grDevices::pdf, useDingbats = FALSE)

selected_delta <- delta_table %>%
  filter(evaluation == "eicu_external", threshold %in% c(0.10, 0.20, 0.30)) %>%
  transmute(
    threshold,
    dynamic_flagged = predicted_positive__dynamic_last_no_support,
    rox_flagged = predicted_positive__rox_only,
    dynamic_events_captured = events_captured__dynamic_last_no_support,
    rox_events_captured = events_captured__rox_only,
    delta_events = events_captured_delta_dynamic_minus_rox,
    delta_false_positives = false_positives_delta_dynamic_minus_rox,
    delta_net_benefit = net_benefit_delta_dynamic_minus_rox
  )
selected_recal <- recalibration %>%
  filter(evaluation == "eicu_external", model == "dynamic_last_no_support", recalibration_method %in% c("original_mimic_trained", "posthoc_intercept_only")) %>%
  transmute(method = recalibration_method, n, events, mean_predicted_risk, brier, calibration_slope, note)

note_path <- file.path(manuscript_dir, paste0("icu_decision_recalibration_interpretation_", version_id, ".md"))
note_lines <- c(
  "# ICU Decision-Threshold and Recalibration Interpretation",
  "",
  paste0("Version: ", version_id),
  "",
  "This V30 module extends the ICU translation layer with threshold-impact and recalibration summaries. It should be framed as supportive external translation of a dynamic respiratory-support risk signal, not as direct validation of the community PEF marker and not as a clinical decision tool.",
  "",
  "## eICU Common-Population Threshold Impact",
  "",
  markdown_table(selected_delta),
  "",
  "## eICU Dynamic-Model Recalibration Diagnostic",
  "",
  markdown_table(selected_recal),
  "",
  "## Wording Boundary",
  "",
  "- Allowed: the dynamic ICU model showed better decision-curve net benefit than a ROX-like comparator across examined thresholds in common complete-case rows.",
  "- Allowed: eICU external validation supports broader respiratory-support transportability.",
  "- Not allowed: the model is ready for clinical decision-making.",
  "- Not allowed: eICU validates exact HFNC-specific use; HFNC-coded eICU rows are sparse and mapping differs from MIMIC.",
  "- Not allowed: this ICU model validates the same construct as the community PEF marker."
)
writeLines(note_lines, note_path)

log_path <- file.path(logs_dir, paste0("icu_decision_recalibration_", version_id, ".md"))
recon_brief <- reconstruction %>%
  filter(model == "dynamic_last_no_support", evaluation %in% c("mimic_internal_test", "eicu_external")) %>%
  transmute(
    evaluation,
    n,
    events,
    auroc = round(auroc, 3),
    v24_auroc = round(v24_auroc, 3),
    brier = round(brier, 3),
    v24_brier = round(v24_brier, 3)
  )
log_lines <- c(
  "# ICU V30 Decision Threshold and Recalibration Log",
  "",
  paste0("Version: ", version_id),
  "",
  "## Status",
  "",
  "T3/T4 ICU decision-threshold, DCA reproduction, recalibration, and transportability outputs built from V23/V24 materials.",
  "",
  "## Prediction Reconstruction Audit",
  "",
  markdown_table(recon_brief),
  "",
  "## Key Files",
  "",
  "- `results/tables/icu_decision_threshold_table_v30_0.csv`",
  "- `results/tables/icu_decision_threshold_delta_v30_0.csv`",
  "- `results/tables/icu_decision_curve_reproduction_v30_0.csv`",
  "- `results/tables/icu_recalibration_transportability_v30_0.csv`",
  "- `results/tables/icu_transportability_summary_v30_0.csv`",
  "- `manuscript/icu_decision_recalibration_interpretation_v30_0.md`",
  "",
  "## Boundary",
  "",
  "Supportive ICU translation only. No claim of clinical deployment, exact HFNC external validation, or direct validation of the community PEF marker."
)
writeLines(log_lines, log_path)

message("ICU V30 decision/recalibration outputs complete.")
