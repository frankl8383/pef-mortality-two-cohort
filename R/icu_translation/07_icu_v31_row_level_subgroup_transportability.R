# ICU V31.6 row-level subgroup transportability.
#
# Adds row-level chronic pulmonary disease, oxygenation-stress and coarse
# blood-gas hypercapnia strata to the existing ICU translation predictions.

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(readr)
  library(stringr)
  library(tibble)
})

version_id <- "v31_6"

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

format_p <- function(p) {
  ifelse(is.na(p), NA_character_, ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
}

md_escape <- function(x) {
  str_replace_all(ifelse(is.na(x), "", as.character(x)), "\\|", "\\\\|")
}

md_table <- function(dat) {
  dat <- as.data.frame(dat, stringsAsFactors = FALSE)
  dat[] <- lapply(dat, function(x) ifelse(is.na(x), "", as.character(x)))
  c(
    paste0("| ", paste(md_escape(names(dat)), collapse = " | "), " |"),
    paste0("| ", paste(rep("---", ncol(dat)), collapse = " | "), " |"),
    apply(dat, 1, function(x) paste0("| ", paste(md_escape(x), collapse = " | "), " |"))
  )
}

clip_prob <- function(p) {
  pmin(pmax(as.numeric(p), 1e-6), 1 - 1e-6)
}

auroc_rank <- function(y, pred) {
  ok <- !is.na(y) & !is.na(pred)
  y <- as.integer(y[ok])
  pred <- as.numeric(pred[ok])
  n1 <- sum(y == 1L)
  n0 <- sum(y == 0L)
  if (n1 == 0L || n0 == 0L) return(NA_real_)
  ranks <- rank(pred, ties.method = "average")
  (sum(ranks[y == 1L]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

average_precision <- function(y, pred) {
  ok <- !is.na(y) & !is.na(pred)
  y <- as.integer(y[ok])
  pred <- as.numeric(pred[ok])
  positives <- sum(y == 1L)
  if (positives == 0L) return(NA_real_)
  ord <- order(pred, decreasing = TRUE)
  y <- y[ord]
  tp <- cumsum(y == 1L)
  fp <- cumsum(y == 0L)
  precision <- tp / pmax(tp + fp, 1)
  recall <- tp / positives
  recall_prev <- c(0, head(recall, -1))
  sum((recall - recall_prev) * precision)
}

calibration_fit <- function(y, pred) {
  ok <- !is.na(y) & !is.na(pred)
  y <- as.integer(y[ok])
  pred <- clip_prob(pred[ok])
  if (length(unique(y)) < 2L || length(y) < 30L) {
    return(tibble(calibration_intercept = NA_real_, calibration_slope = NA_real_))
  }
  lp <- qlogis(pred)
  fit <- tryCatch(
    suppressWarnings(stats::glm(y ~ lp, family = stats::binomial())),
    error = function(e) NULL
  )
  if (is.null(fit)) {
    return(tibble(calibration_intercept = NA_real_, calibration_slope = NA_real_))
  }
  co <- stats::coef(fit)
  tibble(
    calibration_intercept = unname(co[[1]]),
    calibration_slope = unname(co[[2]])
  )
}

performance_row <- function(data, subgroup_domain, subgroup_level, source_definition) {
  df <- data %>% filter(!is.na(pred_dynamic_last_no_support_v30), !is.na(outcome_failure_48h_l2))
  y <- as.integer(df$outcome_failure_48h_l2)
  pred <- as.numeric(df$pred_dynamic_last_no_support_v30)
  events <- sum(y == 1L, na.rm = TRUE)
  n <- nrow(df)
  status <- ifelse(n >= 100 && events >= 10, "ok", "sparse_descriptive_only")
  cal <- calibration_fit(y, pred)
  tibble(
    subgroup_domain = subgroup_domain,
    subgroup_level = subgroup_level,
    evaluation = unique(df$evaluation)[1],
    n = n,
    events = events,
    event_rate = mean(y == 1L, na.rm = TRUE),
    mean_predicted_risk = mean(pred, na.rm = TRUE),
    auroc = auroc_rank(y, pred),
    auprc = average_precision(y, pred),
    brier = mean((y - pred)^2, na.rm = TRUE),
    calibration_intercept = cal$calibration_intercept,
    calibration_slope = cal$calibration_slope,
    status = status,
    source_definition = source_definition,
    manuscript_boundary = case_when(
      status != "ok" ~ "Sparse subgroup; descriptive only.",
      str_detect(subgroup_domain, "hypercapnia") ~ "Coarse blood-gas subgroup; not baseline-window precise.",
      str_detect(subgroup_domain, "chronic") ~ "Diagnosis-history subgroup; coded-history sensitivity, not adjudicated COPD.",
      TRUE ~ "Row-level descriptive transportability subgroup."
    )
  )
}

icd9_chronic_pulmonary <- function(code) {
  code <- str_replace_all(as.character(code), "[^0-9A-Za-z]", "")
  first3 <- suppressWarnings(as.integer(substr(code, 1, 3)))
  (first3 >= 490 & first3 <= 496) |
    (first3 >= 500 & first3 <= 505) |
    str_starts(code, "5064")
}

icd10_chronic_pulmonary <- function(code) {
  code <- str_to_upper(str_replace_all(as.character(code), "[^0-9A-Za-z]", ""))
  str_detect(code, "^J4[0-7]|^J6[0-7]|^J684|^J701|^J703")
}

eicu_chronic_pulmonary_text <- function(text) {
  str_detect(
    str_to_lower(ifelse(is.na(text), "", as.character(text))),
    "copd|chronic obstruct|emphysema|chronic bronchitis|bronchiectasis|asthma|interstitial lung|pulmonary fibrosis"
  )
}

root <- find_project_root()
secure_root <- "${SECURE_DATA_ROOT}"
mimic_root <- file.path(secure_root, "mimiciv", "3.1")
eicu_root <- file.path(secure_root, "eicu-crd", "2.0")
derived_dir <- file.path(root, "derived_sensitive", "icu_translation")
tab_dir <- file.path(root, "results", "tables")
log_dir <- file.path(root, "results", "logs")
manuscript_dir <- file.path(root, "manuscript")
invisible(lapply(c(derived_dir, tab_dir, log_dir, manuscript_dir), dir_create))

pred_path <- file.path(derived_dir, "icu_translation_predictions_v30_0.rds")
if (!file.exists(pred_path)) {
  stop("Missing row-level ICU predictions: ", pred_path, call. = FALSE)
}
icu <- readRDS(pred_path) %>%
  as_tibble() %>%
  mutate(
    oxygenation_stress = case_when(
      !is.na(spo2_min) & spo2_min < 90 ~ TRUE,
      !is.na(fio2_max) & !is.na(spo2_min) & fio2_max >= 0.40 & spo2_min <= 94 ~ TRUE,
      !is.na(rox_l2) & rox_l2 < 5 ~ TRUE,
      TRUE ~ FALSE
    ),
    oxygenation_stress_label = if_else(oxygenation_stress, "oxygenation_stress_present", "oxygenation_stress_absent")
  )

mimic_stays <- data.table::fread(
  file.path(mimic_root, "icu", "icustays.csv.gz"),
  select = c("subject_id", "hadm_id", "stay_id")
) %>%
  as_tibble() %>%
  mutate(across(c(subject_id, hadm_id, stay_id), as.character))

mimic_cohort_hadm <- icu %>%
  filter(dataset == "MIMIC-IV") %>%
  mutate(stay_id = as.character(stay_id), subject_id = as.character(subject_id)) %>%
  inner_join(mimic_stays, by = c("subject_id", "stay_id")) %>%
  filter(!is.na(hadm_id)) %>%
  distinct(hadm_id) %>%
  pull(hadm_id)

mimic_hadm_key <- tempfile("mimic_hadm_key_")
writeLines(mimic_cohort_hadm, mimic_hadm_key)

mimic_diag <- data.table::fread(
  file.path(mimic_root, "hosp", "diagnoses_icd.csv.gz"),
  select = c("subject_id", "hadm_id", "icd_code", "icd_version")
) %>%
  as_tibble() %>%
  mutate(
    subject_id = as.character(subject_id),
    hadm_id = as.character(hadm_id),
    chronic_pulmonary_flag = case_when(
      as.integer(icd_version) == 9L ~ icd9_chronic_pulmonary(icd_code),
      as.integer(icd_version) == 10L ~ icd10_chronic_pulmonary(icd_code),
      TRUE ~ FALSE
    )
  ) %>%
  group_by(subject_id, hadm_id) %>%
  summarise(chronic_pulmonary_disease = any(chronic_pulmonary_flag, na.rm = TRUE), .groups = "drop")

mimic_labevents_path <- file.path(mimic_root, "hosp", "labevents.csv.gz")
mimic_blood_gas_cache <- file.path(derived_dir, "icu_v31_6_mimic_blood_gas_by_admission.csv")
mimic_lab_cmd <- paste(
  "gzip -dc",
  shQuote(mimic_labevents_path),
  "| awk -F,",
  paste0("-v key=", shQuote(mimic_hadm_key)),
  shQuote(paste0(
    "BEGIN {while ((getline line < key) > 0) keep[line]=1} ",
    "NR==1 {print \"subject_id,hadm_id,itemid,valuenum\"} ",
    "NR>1 && (($5==50818 || $5==50821) && ($3 in keep)) {print $2 \",\" $3 \",\" $5 \",\" $10}"
  ))
)
if (file.exists(mimic_blood_gas_cache)) {
  mimic_blood_gas <- read_csv(mimic_blood_gas_cache, show_col_types = FALSE) %>%
    mutate(subject_id = as.character(subject_id), hadm_id = as.character(hadm_id))
} else {
  mimic_blood_gas <- data.table::fread(cmd = mimic_lab_cmd, showProgress = FALSE) %>%
    as_tibble() %>%
    mutate(
      subject_id = as.character(subject_id),
      hadm_id = as.character(hadm_id),
      itemid = as.integer(itemid),
      valuenum = suppressWarnings(as.numeric(valuenum))
    ) %>%
    filter(!is.na(hadm_id), !is.na(valuenum)) %>%
    group_by(subject_id, hadm_id) %>%
    summarise(
      mimic_pco2_max = suppressWarnings(max(valuenum[itemid == 50818], na.rm = TRUE)),
      mimic_po2_min = suppressWarnings(min(valuenum[itemid == 50821], na.rm = TRUE)),
      mimic_pco2_count = sum(itemid == 50818 & !is.na(valuenum)),
      mimic_po2_count = sum(itemid == 50821 & !is.na(valuenum)),
      .groups = "drop"
    ) %>%
    mutate(
      mimic_pco2_max = ifelse(is.infinite(mimic_pco2_max), NA_real_, mimic_pco2_max),
      mimic_po2_min = ifelse(is.infinite(mimic_po2_min), NA_real_, mimic_po2_min)
    )
  write_csv(mimic_blood_gas, mimic_blood_gas_cache)
}

eicu_diag <- data.table::fread(
  file.path(eicu_root, "diagnosis.csv.gz"),
  select = c("patientunitstayid", "diagnosisstring", "icd9code")
) %>%
  as_tibble() %>%
  mutate(
    patientunitstayid = as.character(patientunitstayid),
    chronic_pulmonary_flag = eicu_chronic_pulmonary_text(paste(diagnosisstring, icd9code))
  ) %>%
  group_by(patientunitstayid) %>%
  summarise(chronic_pulmonary_disease_diagnosis = any(chronic_pulmonary_flag, na.rm = TRUE), .groups = "drop")

eicu_history <- data.table::fread(
  file.path(eicu_root, "pastHistory.csv.gz"),
  select = c("patientunitstayid", "pasthistorypath", "pasthistoryvalue", "pasthistoryvaluetext")
) %>%
  as_tibble() %>%
  mutate(
    patientunitstayid = as.character(patientunitstayid),
    chronic_pulmonary_flag = eicu_chronic_pulmonary_text(paste(pasthistorypath, pasthistoryvalue, pasthistoryvaluetext))
  ) %>%
  group_by(patientunitstayid) %>%
  summarise(chronic_pulmonary_disease_history = any(chronic_pulmonary_flag, na.rm = TRUE), .groups = "drop")

eicu_cohort_ids <- icu %>%
  filter(dataset == "eICU-CRD") %>%
  mutate(patientunitstayid = as.character(patientunitstayid)) %>%
  filter(!is.na(patientunitstayid)) %>%
  distinct(patientunitstayid) %>%
  pull(patientunitstayid)
eicu_patient_key <- tempfile("eicu_patient_key_")
writeLines(eicu_cohort_ids, eicu_patient_key)
eicu_lab_path <- file.path(eicu_root, "lab.csv.gz")
eicu_blood_gas_cache <- file.path(derived_dir, "icu_v31_6_eicu_blood_gas_by_stay.csv")
eicu_lab_cmd <- paste(
  "gzip -dc",
  shQuote(eicu_lab_path),
  "| awk -F,",
  paste0("-v key=", shQuote(eicu_patient_key)),
  shQuote(paste0(
    "BEGIN {while ((getline line < key) > 0) keep[line]=1} ",
    "NR==1 {print \"patientunitstayid,labname,labresult\"} ",
    "NR>1 && (($2 in keep) && ($5==\"paCO2\" || $5==\"paO2\")) {print $2 \",\" $5 \",\" $6}"
  ))
)
if (file.exists(eicu_blood_gas_cache)) {
  eicu_lab <- read_csv(eicu_blood_gas_cache, show_col_types = FALSE) %>%
    mutate(patientunitstayid = as.character(patientunitstayid))
} else {
  eicu_lab <- data.table::fread(cmd = eicu_lab_cmd, showProgress = FALSE) %>%
    as_tibble() %>%
    mutate(
      patientunitstayid = as.character(patientunitstayid),
      labname_l = str_to_lower(labname),
      labresult = suppressWarnings(as.numeric(labresult))
    ) %>%
    filter(labname_l %in% c("paco2", "pao2"), !is.na(labresult)) %>%
    group_by(patientunitstayid) %>%
    summarise(
      eicu_pco2_max = suppressWarnings(max(labresult[labname_l == "paco2"], na.rm = TRUE)),
      eicu_po2_min = suppressWarnings(min(labresult[labname_l == "pao2"], na.rm = TRUE)),
      eicu_pco2_count = sum(labname_l == "paco2" & !is.na(labresult)),
      eicu_po2_count = sum(labname_l == "pao2" & !is.na(labresult)),
      .groups = "drop"
    ) %>%
    mutate(
      eicu_pco2_max = ifelse(is.infinite(eicu_pco2_max), NA_real_, eicu_pco2_max),
      eicu_po2_min = ifelse(is.infinite(eicu_po2_min), NA_real_, eicu_po2_min)
    )
  write_csv(eicu_lab, eicu_blood_gas_cache)
}

mimic_aug <- icu %>%
  filter(dataset == "MIMIC-IV") %>%
  mutate(stay_id = as.character(stay_id), subject_id = as.character(subject_id)) %>%
  left_join(mimic_stays, by = c("subject_id", "stay_id")) %>%
  left_join(mimic_diag, by = c("subject_id", "hadm_id")) %>%
  left_join(mimic_blood_gas, by = c("subject_id", "hadm_id")) %>%
  mutate(
    subject_id = as.character(subject_id),
    stay_id = as.character(stay_id),
    patientunitstayid = as.character(patientunitstayid),
    chronic_pulmonary_disease = coalesce(chronic_pulmonary_disease, FALSE),
    pco2_max = mimic_pco2_max,
    po2_min = mimic_po2_min,
    pco2_count = mimic_pco2_count,
    po2_count = mimic_po2_count
  )

eicu_aug <- icu %>%
  filter(dataset == "eICU-CRD") %>%
  mutate(patientunitstayid = as.character(patientunitstayid)) %>%
  left_join(eicu_diag, by = "patientunitstayid") %>%
  left_join(eicu_history, by = "patientunitstayid") %>%
  left_join(eicu_lab, by = "patientunitstayid") %>%
  mutate(
    subject_id = as.character(subject_id),
    stay_id = as.character(stay_id),
    patientunitstayid = as.character(patientunitstayid),
    chronic_pulmonary_disease = coalesce(chronic_pulmonary_disease_diagnosis, FALSE) | coalesce(chronic_pulmonary_disease_history, FALSE),
    pco2_max = eicu_pco2_max,
    po2_min = eicu_po2_min,
    pco2_count = eicu_pco2_count,
    po2_count = eicu_po2_count
  )

aug <- bind_rows(mimic_aug, eicu_aug) %>%
  mutate(
    chronic_pulmonary_label = if_else(chronic_pulmonary_disease, "chronic_pulmonary_present", "chronic_pulmonary_absent"),
    hypercapnia_label = case_when(
      is.na(pco2_count) | pco2_count == 0 ~ "pco2_not_measured",
      pco2_max >= 50 ~ "hypercapnia_present",
      TRUE ~ "hypercapnia_absent_measured"
    ),
    low_po2_label = case_when(
      is.na(po2_count) | po2_count == 0 ~ "po2_not_measured",
      po2_min < 60 ~ "po2_low_present",
      TRUE ~ "po2_low_absent_measured"
    )
  )

aug_out <- aug %>%
  select(
    dataset, subject_id, stay_id, patientunitstayid, evaluation, split,
    outcome_failure_48h_l2, pred_dynamic_last_no_support_v30,
    age, male, support_class, spo2_min, fio2_max, rox_l2,
    chronic_pulmonary_disease, chronic_pulmonary_label,
    oxygenation_stress, oxygenation_stress_label,
    pco2_max, po2_min, pco2_count, po2_count, hypercapnia_label, low_po2_label
  )
readr::write_csv(aug_out, file.path(derived_dir, "icu_v31_6_row_level_subgroup_augmented.csv.gz"))

eval_data <- aug %>%
  filter(evaluation %in% c("mimic_internal_test", "eicu_external"))

subgroup_specs <- tibble::tribble(
  ~subgroup_domain, ~column, ~source_definition,
  "coded_chronic_pulmonary_disease", "chronic_pulmonary_label", "MIMIC ICD-9/10 chronic pulmonary codes; eICU diagnosis/past-history text/code screen.",
  "common_feature_oxygenation_stress", "oxygenation_stress_label", "SpO2_min <90, or FiO2_max >=0.40 with SpO2_min <=94, or ROX-like index <5 in the common feature window.",
  "coarse_blood_gas_hypercapnia", "hypercapnia_label", "Any same-admission/stay pCO2 >=50 mmHg when pCO2 was measured; not baseline-window precise.",
  "coarse_blood_gas_low_po2", "low_po2_label", "Any same-admission/stay pO2 <60 mmHg when pO2 was measured; not baseline-window precise."
)

perf <- bind_rows(lapply(seq_len(nrow(subgroup_specs)), function(i) {
  spec <- subgroup_specs[i, ]
  bind_rows(lapply(sort(unique(eval_data[[spec$column]])), function(level) {
    bind_rows(lapply(c("mimic_internal_test", "eicu_external"), function(eval_name) {
      performance_row(
        eval_data %>% filter(evaluation == eval_name, .data[[spec$column]] == level),
        subgroup_domain = spec$subgroup_domain,
        subgroup_level = as.character(level),
        source_definition = spec$source_definition
      )
    }))
  }))
})) %>%
  mutate(
    status = if_else(str_detect(subgroup_level, "not_measured"), "missingness_stratum_not_phenotype", status),
    manuscript_boundary = if_else(
      str_detect(subgroup_level, "not_measured"),
      "Missingness stratum; not a physiologic phenotype.",
      manuscript_boundary
    )
  ) %>%
  mutate(across(c(event_rate, mean_predicted_risk, auroc, auprc, brier, calibration_intercept, calibration_slope), ~ round(.x, 4)))

flag_counts <- aug %>%
  filter(evaluation %in% c("mimic_internal_test", "eicu_external")) %>%
  group_by(evaluation) %>%
  summarise(
    n = n(),
    events = sum(outcome_failure_48h_l2 == 1L, na.rm = TRUE),
    chronic_pulmonary_present = sum(chronic_pulmonary_disease, na.rm = TRUE),
    oxygenation_stress_present = sum(oxygenation_stress, na.rm = TRUE),
    pco2_measured = sum(!is.na(pco2_count) & pco2_count > 0, na.rm = TRUE),
    hypercapnia_present = sum(hypercapnia_label == "hypercapnia_present", na.rm = TRUE),
    po2_measured = sum(!is.na(po2_count) & po2_count > 0, na.rm = TRUE),
    po2_low_present = sum(low_po2_label == "po2_low_present", na.rm = TRUE),
    .groups = "drop"
  )

ledger <- tibble::tribble(
  ~domain, ~status, ~definition, ~boundary,
  "Chronic pulmonary disease", "estimated", "MIMIC ICD-9/10 chronic pulmonary disease codes; eICU diagnosis and past-history text/code screen.", "Coded-history sensitivity; not adjudicated COPD.",
  "Hypoxemic/oxygenation-stress phenotype", "estimated", "Common feature window using SpO2, FiO2 and ROX-like index.", "Comparable across MIMIC/eICU because it uses the V22 common model features.",
  "Hypercapnia phenotype", "estimated_with_boundary", "Any same-admission/stay pCO2 >=50 mmHg when pCO2 measured.", "Coarse blood-gas subgroup; not baseline-window precise because current V22 prediction rows do not retain T0 timestamps.",
  "Low pO2 phenotype", "estimated_with_boundary", "Any same-admission/stay pO2 <60 mmHg when pO2 measured.", "Coarse blood-gas subgroup; not baseline-window precise."
)

write_csv(perf, file.path(tab_dir, "icu_v31_6_row_level_subgroup_performance.csv"))
write_csv(flag_counts, file.path(tab_dir, "icu_v31_6_row_level_subgroup_flag_counts.csv"))
write_csv(ledger, file.path(tab_dir, "icu_v31_6_row_level_subgroup_extraction_ledger.csv"))

display_perf <- perf %>%
  mutate(
    event_rate = sprintf("%.1f%%", 100 * event_rate),
    mean_predicted_risk = sprintf("%.1f%%", 100 * mean_predicted_risk),
    auroc = ifelse(is.na(auroc), "", sprintf("%.3f", auroc)),
    auprc = ifelse(is.na(auprc), "", sprintf("%.3f", auprc)),
    brier = ifelse(is.na(brier), "", sprintf("%.3f", brier)),
    calibration_slope = ifelse(is.na(calibration_slope), "", sprintf("%.2f", calibration_slope))
  ) %>%
  select(subgroup_domain, subgroup_level, evaluation, n, events, event_rate, mean_predicted_risk, auroc, auprc, brier, calibration_slope, status, manuscript_boundary)

note <- c(
  "# ICU Row-Level Subgroup Transportability V31.6",
  "",
  "## Summary",
  "",
  "The V31.6 row-level rerun resolves the previous v31.5 aggregate-table blocker for chronic pulmonary disease and oxygenation-related subgroups. Chronic pulmonary disease is now available as a coded-history sensitivity; oxygenation stress is available from common SpO2/FiO2/ROX-like features. Hypercapnia and low pO2 are estimated from same-admission/stay blood-gas records and must be labeled coarse rather than baseline-window precise.",
  "",
  "## Extraction ledger",
  "",
  md_table(ledger),
  "",
  "## Flag counts",
  "",
  md_table(flag_counts),
  "",
  "## Subgroup performance",
  "",
  md_table(display_perf),
  "",
  "## Manuscript boundary",
  "",
  "Use the chronic pulmonary disease and oxygenation-stress rows as descriptive ICU transportability checks. Hypercapnia/low-pO2 rows can be included only as sensitivity or supplement because the current V22 row-level prediction table does not retain T0 timestamps for exact baseline-window blood-gas alignment. Do not frame these rows as treatment-effect estimates, direct validation of the population PEF marker, or adjudicated respiratory phenotypes."
)
writeLines(note, file.path(manuscript_dir, "icu_row_level_subgroup_transportability_note_v31_6.md"))

log <- c(
  "# ICU V31.6 Row-Level Subgroup Transportability Build",
  "",
  paste0("MIMIC evaluation rows: ", sum(eval_data$evaluation == "mimic_internal_test")),
  paste0("eICU evaluation rows: ", sum(eval_data$evaluation == "eicu_external")),
  paste0("Performance rows written: ", nrow(perf)),
  "",
  "Outputs:",
  "- derived_sensitive/icu_translation/icu_v31_6_row_level_subgroup_augmented.csv.gz",
  "- results/tables/icu_v31_6_row_level_subgroup_performance.csv",
  "- results/tables/icu_v31_6_row_level_subgroup_flag_counts.csv",
  "- results/tables/icu_v31_6_row_level_subgroup_extraction_ledger.csv",
  "- manuscript/icu_row_level_subgroup_transportability_note_v31_6.md"
)
writeLines(log, file.path(log_dir, "icu_v31_6_row_level_subgroup_transportability.md"))

message("ICU V31.6 row-level subgroup transportability outputs written.")
