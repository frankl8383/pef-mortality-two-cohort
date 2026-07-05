# CHARLS provisional core harmonized cleaning.
# Creates a local-only row-level analysis dataset plus aggregate QC outputs.

suppressPackageStartupMessages({
  library(dplyr)
  library(haven)
  library(readr)
  library(tidyselect)
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

extract_member <- function(archive_local_path, archive_type, inner_path) {
  td <- tempfile("charls_core_clean_")
  dir.create(td, recursive = TRUE)
  out <- file.path(td, basename(inner_path))
  if (archive_type == "zip") {
    utils::unzip(archive_local_path, files = inner_path, exdir = td, overwrite = TRUE)
    extracted <- file.path(td, inner_path)
    if (!file.exists(extracted)) {
      extracted <- out
    }
    return(list(path = extracted, tmpdir = td))
  }
  if (archive_type == "rar") {
    err <- tempfile("bsdtar_err_")
    status <- system2("bsdtar", c("-xOf", archive_local_path, inner_path), stdout = out, stderr = err)
    err_text <- if (file.exists(err)) paste(readLines(err, warn = FALSE), collapse = " | ") else ""
    unlink(err)
    if (!identical(as.integer(status), 0L) || !file.exists(out) || file.info(out)$size == 0) {
      unlink(td, recursive = TRUE)
      stop("Failed to extract RAR member: ", archive_local_path, " :: ", inner_path, " ", err_text, call. = FALSE)
    }
    return(list(path = out, tmpdir = td))
  }
  stop("Unsupported archive type: ", archive_type, call. = FALSE)
}

read_selected <- function(index_row, selected_cols) {
  member <- extract_member(index_row$archive_local_path, index_row$archive_type, index_row$inner_path)
  on.exit(unlink(member$tmpdir, recursive = TRUE), add = TRUE)
  haven::read_dta(member$path, col_select = tidyselect::any_of(unique(selected_cols)))
}

to_num <- function(x) {
  suppressWarnings(as.numeric(haven::zap_labels(haven::zap_missing(x))))
}

to_chr <- function(x) {
  as.character(haven::zap_labels(haven::zap_missing(x)))
}

col_num <- function(dat, name) {
  if (name %in% names(dat)) dat[[name]] else rep(NA_real_, nrow(dat))
}

any_yes <- function(...) {
  mat <- cbind(...)
  yes <- rowSums(mat == 1, na.rm = TRUE) > 0
  observed <- rowSums(!is.na(mat)) > 0
  ifelse(yes, 1L, ifelse(observed, 0L, NA_integer_))
}

any_value <- function(values, ...) {
  mat <- cbind(...)
  hit <- mat %in% values
  dim(hit) <- dim(mat)
  observed <- rowSums(!is.na(mat)) > 0
  ifelse(rowSums(hit, na.rm = TRUE) > 0, 1L, ifelse(observed, 0L, NA_integer_))
}

followup_incident <- function(baseline, ...) {
  mat <- cbind(...)
  event <- rowSums(mat == 1, na.rm = TRUE) > 0
  observed <- rowSums(!is.na(mat)) > 0
  ifelse(
    baseline == 0 & event,
    1L,
    ifelse(baseline == 0 & observed, 0L, NA_integer_)
  )
}

first_event_wave <- function(...) {
  mat <- cbind(...)
  out <- rep(NA_integer_, nrow(mat))
  for (j in seq_len(ncol(mat))) {
    hit <- is.na(out) & mat[, j] == 1
    out[hit] <- j + 1L
  }
  out
}

all_no_activity <- function(vigorous, moderate, light) {
  mat <- cbind(vigorous, moderate, light)
  any_active <- rowSums(mat == 1, na.rm = TRUE) > 0
  all_observed <- rowSums(!is.na(mat)) == ncol(mat)
  all_no <- rowSums(mat == 0, na.rm = TRUE) == ncol(mat)
  ifelse(any_active, 0L, ifelse(all_observed & all_no, 1L, NA_integer_))
}

num_summary <- function(x) {
  z <- suppressWarnings(as.numeric(x))
  z <- z[!is.na(z)]
  if (length(z) == 0) {
    return(tibble::tibble(
      n_nonmissing = 0L,
      min = NA_real_, p01 = NA_real_, p05 = NA_real_, median = NA_real_,
      mean = NA_real_, p95 = NA_real_, p99 = NA_real_, max = NA_real_
    ))
  }
  qs <- as.numeric(stats::quantile(z, probs = c(0.01, 0.05, 0.5, 0.95, 0.99), na.rm = TRUE, names = FALSE))
  tibble::tibble(
    n_nonmissing = length(z),
    min = min(z, na.rm = TRUE),
    p01 = qs[[1]],
    p05 = qs[[2]],
    median = qs[[3]],
    mean = mean(z, na.rm = TRUE),
    p95 = qs[[4]],
    p99 = qs[[5]],
    max = max(z, na.rm = TRUE)
  )
}

event_summary <- function(core, outcome, baseline, incident, followup_observed) {
  at_risk <- baseline == 0
  at_risk_followed <- at_risk & followup_observed
  events <- incident == 1
  tibble::tibble(
    outcome = outcome,
    total_rows = nrow(core),
    baseline_yes = sum(baseline == 1, na.rm = TRUE),
    baseline_no = sum(at_risk, na.rm = TRUE),
    baseline_missing = sum(is.na(baseline)),
    at_risk_with_followup = sum(at_risk_followed, na.rm = TRUE),
    incident_events = sum(events, na.rm = TRUE),
    incident_risk_among_followed = ifelse(sum(at_risk_followed, na.rm = TRUE) > 0,
      sum(events, na.rm = TRUE) / sum(at_risk_followed, na.rm = TRUE),
      NA_real_
    )
  )
}

write_missingness <- function(core, out_path) {
  miss <- lapply(names(core), function(v) {
    x <- core[[v]]
    n_missing <- sum(is.na(x))
    tibble::tibble(
      variable = v,
      n = length(x),
      n_missing = n_missing,
      n_nonmissing = length(x) - n_missing,
      pct_missing = n_missing / length(x)
    )
  })
  readr::write_csv(dplyr::bind_rows(miss), out_path)
}

add_pef_residual <- function(core) {
  core$pef_pred_w1 <- NA_real_
  core$pef_resid_w1 <- NA_real_
  core$pef_resid_z_w1 <- NA_real_
  complete <- !is.na(core$pef_best_w1_valid_provisional) &
    !is.na(core$age_w1) &
    !is.na(core$sex_code) &
    !is.na(core$height_m_w1)
  model_df <- core[complete, , drop = FALSE]
  if (nrow(model_df) < 50 || dplyr::n_distinct(model_df$sex_code) < 2) {
    attr(core, "pef_model_status") <- "not_fitted_insufficient_complete_cases"
    return(core)
  }
  fit <- stats::lm(pef_best_w1_valid_provisional ~ age_w1 + sex_code + height_m_w1, data = model_df)
  core$pef_pred_w1[complete] <- stats::predict(fit, newdata = model_df)
  core$pef_resid_w1[complete] <- core$pef_best_w1_valid_provisional[complete] - core$pef_pred_w1[complete]
  resid_sd <- stats::sd(stats::residuals(fit), na.rm = TRUE)
  if (!is.na(resid_sd) && resid_sd > 0) {
    core$pef_resid_z_w1[complete] <- core$pef_resid_w1[complete] / resid_sd
  }
  attr(core, "pef_model_status") <- paste0(
    "fitted_n=", nrow(model_df),
    "; formula=pef_best_w1_valid_provisional ~ age_w1 + sex_code + height_m_w1",
    "; residual_sd=", signif(resid_sd, 4)
  )
  core
}

main <- function() {
  root <- find_project_root()
  index_path <- file.path(root, "derived", "charls_wave_file_index.csv")
  if (!file.exists(index_path)) {
    stop("Run 00_index_charls_archives.R before core cleaning.", call. = FALSE)
  }

  table_dir <- file.path(root, "results", "tables")
  log_dir <- file.path(root, "results", "logs")
  sensitive_dir <- file.path(root, "derived_sensitive", "charls")
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(sensitive_dir, recursive = TRUE, showWarnings = FALSE)

  index <- readr::read_csv(index_path, show_col_types = FALSE, progress = FALSE)
  h_row <- index %>%
    filter(preferred, inner_path == "H_CHARLS_D_Data.dta") %>%
    slice(1)
  if (nrow(h_row) != 1) {
    stop("Could not locate preferred H_CHARLS_D_Data.dta in charls_wave_file_index.csv.", call. = FALSE)
  }

  id_vars <- c("ID", "householdID", "communityID", "hhid")
  fixed_vars <- c("rabyear", "ragender")
  wave_suffixes <- c(
    "iwstat", "wtresp", "wtrespbioa", "wtrespbiob", "agey",
    "lunge", "asthmae",
    "walk1kma", "walk100a", "adla_c", "adlab_c",
    "vgact_c", "mdact_c", "ltact_c",
    "smokev", "smoken", "cesd10"
  )
  biomarker_suffixes <- c(
    "mheight", "mweight", "mbmi", "gripsum",
    "puff1", "puff2", "puff3", "puff", "puffcomp", "puffpos", "puffeff"
  )
  selected_cols <- unique(c(
    id_vars,
    fixed_vars,
    as.vector(outer(paste0("r", 1:4), wave_suffixes, paste0)),
    as.vector(outer(paste0("r", 1:3), biomarker_suffixes, paste0))
  ))

  message("Reading selected columns from H_CHARLS_D_Data.dta...")
  raw <- read_selected(h_row, selected_cols)
  core <- tibble::as_tibble(raw)
  present_id_vars <- intersect(id_vars, names(core))
  for (v in present_id_vars) {
    core[[v]] <- to_chr(core[[v]])
  }
  numeric_vars <- setdiff(names(core), present_id_vars)
  for (v in numeric_vars) {
    core[[v]] <- to_num(core[[v]])
  }
  selected_present_n <- ncol(core)

  core <- core %>%
    mutate(
      provisional_rule_version = "harmonized_core_v0_1_2026_06_21",
      participant_id = ID,
      household_id = householdID,
      community_id = communityID,
      sex_code = ragender,
      sex_label = case_when(
        sex_code == 1 ~ "male",
        sex_code == 2 ~ "female",
        TRUE ~ NA_character_
      ),
      birth_year = rabyear,
      age_w1 = r1agey,
      age_w2 = col_num(core, "r2agey"),
      age_w3 = col_num(core, "r3agey"),
      age_w4 = col_num(core, "r4agey"),
      height_m_w1 = col_num(core, "r1mheight"),
      weight_kg_w1 = col_num(core, "r1mweight"),
      bmi_w1 = col_num(core, "r1mbmi"),
      pef_trial1_w1 = col_num(core, "r1puff1"),
      pef_trial2_w1 = col_num(core, "r1puff2"),
      pef_trial3_w1 = col_num(core, "r1puff3"),
      pef_best_w1 = col_num(core, "r1puff"),
      pef_best_w2 = col_num(core, "r2puff"),
      pef_best_w3 = col_num(core, "r3puff"),
      pef_complete_w1 = col_num(core, "r1puffcomp"),
      pef_position_w1 = col_num(core, "r1puffpos"),
      pef_effort_w1 = col_num(core, "r1puffeff"),
      pef_best_w1_high_tail_flag = ifelse(!is.na(pef_best_w1) & pef_best_w1 >= 900, 1L, 0L),
      pef_best_w1_valid_provisional = ifelse(
        !is.na(pef_best_w1) & pef_best_w1 >= 30 & pef_best_w1 <= 999,
        pef_best_w1,
        NA_real_
      ),
      lung_w1 = col_num(core, "r1lunge"),
      lung_w2 = col_num(core, "r2lunge"),
      lung_w3 = col_num(core, "r3lunge"),
      lung_w4 = col_num(core, "r4lunge"),
      asthma_w1 = col_num(core, "r1asthmae"),
      asthma_w2 = col_num(core, "r2asthmae"),
      asthma_w3 = col_num(core, "r3asthmae"),
      asthma_w4 = col_num(core, "r4asthmae"),
      followup_lung_observed_w2_w4 = rowSums(!is.na(cbind(lung_w2, lung_w3, lung_w4))) > 0,
      followup_asthma_observed_w2_w4 = rowSums(!is.na(cbind(asthma_w2, asthma_w3, asthma_w4))) > 0,
      incident_lung_w2_w4 = followup_incident(lung_w1, lung_w2, lung_w3, lung_w4),
      incident_asthma_w2_w4 = followup_incident(asthma_w1, asthma_w2, asthma_w3, asthma_w4),
      lung_event_wave = ifelse(incident_lung_w2_w4 == 1, first_event_wave(lung_w2, lung_w3, lung_w4), NA_integer_),
      asthma_event_wave = ifelse(incident_asthma_w2_w4 == 1, first_event_wave(asthma_w2, asthma_w3, asthma_w4), NA_integer_),
      iwstat_w1 = col_num(core, "r1iwstat"),
      iwstat_w2 = col_num(core, "r2iwstat"),
      iwstat_w3 = col_num(core, "r3iwstat"),
      iwstat_w4 = col_num(core, "r4iwstat"),
      died_w2_w4 = any_value(c(5, 6), iwstat_w2, iwstat_w3, iwstat_w4),
      smoke_ever_w1 = col_num(core, "r1smokev"),
      smoke_now_w1 = col_num(core, "r1smoken"),
      weight_resp_w1 = col_num(core, "r1wtresp"),
      weight_biomarker_w1 = ifelse(
        !is.na(col_num(core, "r1wtrespbiob")),
        col_num(core, "r1wtrespbiob"),
        col_num(core, "r1wtrespbioa")
      ),
      grip_kg_w1 = col_num(core, "r1gripsum"),
      walk_1km_diff_w1 = col_num(core, "r1walk1kma"),
      walk_100m_diff_w1 = col_num(core, "r1walk100a"),
      adl4_score_w1 = col_num(core, "r1adla_c"),
      adl6_score_w1 = col_num(core, "r1adlab_c"),
      cesd10_w1 = col_num(core, "r1cesd10"),
      vigorous_activity_w1 = col_num(core, "r1vgact_c"),
      moderate_activity_w1 = col_num(core, "r1mdact_c"),
      light_activity_w1 = col_num(core, "r1ltact_c"),
      mobility_difficulty_w1 = any_yes(walk_1km_diff_w1, walk_100m_diff_w1),
      adl_difficulty_w1 = ifelse(!is.na(adl6_score_w1), as.integer(adl6_score_w1 > 0), NA_integer_),
      depressive_symptoms_w1 = ifelse(!is.na(cesd10_w1), as.integer(cesd10_w1 >= 10), NA_integer_),
      low_activity_w1 = all_no_activity(vigorous_activity_w1, moderate_activity_w1, light_activity_w1),
      baseline_age45plus = !is.na(age_w1) & age_w1 >= 45,
      baseline_pef_available = !is.na(pef_best_w1_valid_provisional),
      baseline_core_covariates_available = !is.na(age_w1) & !is.na(sex_code) & !is.na(height_m_w1),
      baseline_respiratory_axis_eligible = baseline_age45plus & baseline_pef_available & baseline_core_covariates_available
    )

  low_grip <- rep(NA_integer_, nrow(core))
  grip_cutoffs <- list()
  for (sex in c(1, 2)) {
    idx <- which(core$sex_code == sex & !is.na(core$grip_kg_w1) & core$grip_kg_w1 > 0)
    if (length(idx) >= 50) {
      cutoff <- as.numeric(stats::quantile(core$grip_kg_w1[idx], probs = 0.2, na.rm = TRUE, names = FALSE))
      low_grip[idx] <- as.integer(core$grip_kg_w1[idx] <= cutoff)
      grip_cutoffs[[as.character(sex)]] <- cutoff
    }
  }
  core$low_grip_w1 <- low_grip
  frailty_components <- cbind(
    core$low_grip_w1,
    core$mobility_difficulty_w1,
    core$adl_difficulty_w1,
    core$depressive_symptoms_w1,
    core$low_activity_w1
  )
  core$frailty_proxy_components_observed_w1 <- rowSums(!is.na(frailty_components))
  core$frailty_proxy_count_w1 <- ifelse(
    core$frailty_proxy_components_observed_w1 >= 3,
    rowSums(frailty_components == 1, na.rm = TRUE),
    NA_integer_
  )
  core$frailty_proxy_ge3_w1 <- ifelse(
    !is.na(core$frailty_proxy_count_w1),
    as.integer(core$frailty_proxy_count_w1 >= 3),
    NA_integer_
  )

  core <- add_pef_residual(core)
  pef_model_status <- attr(core, "pef_model_status")

  out_rds <- file.path(sensitive_dir, "charls_core_harmonized_provisional.rds")
  saveRDS(core, out_rds)

  row_counts <- tibble::tibble(
    metric = c(
      "total_harmonized_person_rows",
      "baseline_age45plus",
      "baseline_pef_available",
      "baseline_core_covariates_available",
      "baseline_respiratory_axis_eligible",
      "baseline_lung_no_pef_covariate_eligible",
      "baseline_asthma_no_pef_covariate_eligible",
      "incident_lung_events_w2_w4",
      "incident_asthma_events_w2_w4",
      "died_w2_w4"
    ),
    n = c(
      nrow(core),
      sum(core$baseline_age45plus, na.rm = TRUE),
      sum(core$baseline_pef_available, na.rm = TRUE),
      sum(core$baseline_core_covariates_available, na.rm = TRUE),
      sum(core$baseline_respiratory_axis_eligible, na.rm = TRUE),
      sum(core$baseline_respiratory_axis_eligible & core$lung_w1 == 0, na.rm = TRUE),
      sum(core$baseline_respiratory_axis_eligible & core$asthma_w1 == 0, na.rm = TRUE),
      sum(core$incident_lung_w2_w4 == 1, na.rm = TRUE),
      sum(core$incident_asthma_w2_w4 == 1, na.rm = TRUE),
      sum(core$died_w2_w4 == 1, na.rm = TRUE)
    )
  )
  readr::write_csv(row_counts, file.path(table_dir, "charls_core_cleaning_row_counts.csv"))

  pef_summary <- dplyr::bind_rows(
    dplyr::bind_cols(tibble::tibble(variable = "pef_best_w1"), num_summary(core$pef_best_w1)),
    dplyr::bind_cols(tibble::tibble(variable = "pef_best_w1_valid_provisional"), num_summary(core$pef_best_w1_valid_provisional)),
    dplyr::bind_cols(tibble::tibble(variable = "pef_resid_z_w1"), num_summary(core$pef_resid_z_w1))
  ) %>%
    mutate(
      high_tail_n = ifelse(variable == "pef_best_w1", sum(core$pef_best_w1_high_tail_flag == 1, na.rm = TRUE), NA_integer_),
      high_tail_rule = ifelse(variable == "pef_best_w1", "flagged_when_pef_best_w1_ge_900_not_removed", NA_character_)
    )
  readr::write_csv(pef_summary, file.path(table_dir, "charls_core_cleaning_pef_summary.csv"))

  events <- dplyr::bind_rows(
    event_summary(core, "chronic_lung_disease", core$lung_w1, core$incident_lung_w2_w4, core$followup_lung_observed_w2_w4),
    event_summary(core, "asthma", core$asthma_w1, core$incident_asthma_w2_w4, core$followup_asthma_observed_w2_w4)
  )
  readr::write_csv(events, file.path(table_dir, "charls_core_cleaning_event_counts.csv"))

  frailty_summary <- tibble::tibble(
    component = c("low_grip", "mobility_difficulty", "adl_difficulty", "depressive_symptoms", "low_activity", "frailty_proxy_ge3"),
    n_observed = c(
      sum(!is.na(core$low_grip_w1)),
      sum(!is.na(core$mobility_difficulty_w1)),
      sum(!is.na(core$adl_difficulty_w1)),
      sum(!is.na(core$depressive_symptoms_w1)),
      sum(!is.na(core$low_activity_w1)),
      sum(!is.na(core$frailty_proxy_ge3_w1))
    ),
    n_positive = c(
      sum(core$low_grip_w1 == 1, na.rm = TRUE),
      sum(core$mobility_difficulty_w1 == 1, na.rm = TRUE),
      sum(core$adl_difficulty_w1 == 1, na.rm = TRUE),
      sum(core$depressive_symptoms_w1 == 1, na.rm = TRUE),
      sum(core$low_activity_w1 == 1, na.rm = TRUE),
      sum(core$frailty_proxy_ge3_w1 == 1, na.rm = TRUE)
    )
  ) %>%
    mutate(pct_positive_observed = ifelse(n_observed > 0, n_positive / n_observed, NA_real_))
  readr::write_csv(frailty_summary, file.path(table_dir, "charls_core_cleaning_frailty_summary.csv"))

  key_missingness_vars <- c(
    "participant_id", "age_w1", "sex_code", "height_m_w1", "weight_kg_w1", "bmi_w1",
    "pef_best_w1", "pef_best_w1_valid_provisional", "pef_resid_z_w1",
    "lung_w1", "lung_w2", "lung_w3", "lung_w4",
    "asthma_w1", "asthma_w2", "asthma_w3", "asthma_w4",
    "smoke_ever_w1", "smoke_now_w1", "grip_kg_w1", "cesd10_w1",
    "frailty_proxy_count_w1", "incident_lung_w2_w4", "incident_asthma_w2_w4",
    "weight_resp_w1", "weight_biomarker_w1"
  )
  write_missingness(core[, intersect(key_missingness_vars, names(core)), drop = FALSE],
    file.path(table_dir, "charls_core_cleaning_missingness.csv")
  )

  grip_cutoff_text <- paste(
    names(grip_cutoffs),
    unlist(grip_cutoffs),
    sep = "=",
    collapse = "; "
  )
  log <- c(
    "# CHARLS Core Harmonized Cleaning Log",
    "",
    paste0("- Source file: ", h_row$archive_path, " :: ", h_row$inner_path),
    paste0("- Source local path: ", h_row$archive_local_path),
    paste0("- Selected columns requested: ", length(selected_cols)),
    paste0("- Selected source columns present: ", selected_present_n),
    paste0("- Final columns after derived variables: ", ncol(core)),
    paste0("- Row-level output: ", out_rds),
    paste0("- Rows: ", nrow(core)),
    paste0("- Rule version: ", unique(core$provisional_rule_version)),
    "",
    "## Provisional Decisions",
    "",
    "- This is a compromise analysis-ready dataset, not a final CHARLS-wide master database.",
    "- Harmonized CHARLS is used as the primary source for age, sex, PEF, chronic lung disease, asthma, frailty proxy variables, interview status, and weights.",
    "- Baseline PEF is kept when 30 <= r1puff <= 999; values >= 900 are flagged but not removed pending official codebook confirmation.",
    "- Baseline PEF residual z-score is adjusted for age, sex code, and measured height.",
    paste0("- PEF residual model status: ", pef_model_status),
    "- Incident chronic lung disease/asthma requires baseline no disease and any follow-up ever-diagnosed status equal to 1 in waves 2-4.",
    "- Baseline frailty proxy uses five provisional components: sex-specific lowest grip quintile, walking difficulty, ADL difficulty, CESD-10 >= 10, and no vigorous/moderate/light activity.",
    paste0("- Grip low-quintile cutoffs by sex code: ", ifelse(nchar(grip_cutoff_text) > 0, grip_cutoff_text, "not available")),
    "",
    "## Outputs",
    "",
    "- Aggregate row counts: results/tables/charls_core_cleaning_row_counts.csv",
    "- Aggregate missingness: results/tables/charls_core_cleaning_missingness.csv",
    "- PEF summary: results/tables/charls_core_cleaning_pef_summary.csv",
    "- Event counts: results/tables/charls_core_cleaning_event_counts.csv",
    "- Frailty summary: results/tables/charls_core_cleaning_frailty_summary.csv",
    "",
    "## Governance",
    "",
    "The RDS file is row-level restricted data under derived_sensitive/ and must not be committed or shared. Only aggregate QC outputs are written under results/."
  )
  writeLines(log, file.path(log_dir, "charls_core_cleaning_log.md"))

  message("Wrote provisional CHARLS core harmonized dataset and aggregate QC outputs.")
  message("Sensitive row-level RDS: ", out_rds)
}

if (sys.nframe() == 0) {
  main()
}
