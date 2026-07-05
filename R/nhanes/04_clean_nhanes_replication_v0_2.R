# NHANES v0.2 replication cleaner.
# Merges public NHANES 2007-2012 XPT files and creates a local-only
# analysis-ready dataset for respiratory vulnerability replication.

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

local_raw_dir <- function(root) {
  env_dir <- Sys.getenv("NHANES_RAW_DIR", unset = "")
  if (nzchar(env_dir)) {
    return(normalizePath(env_dir, mustWork = FALSE))
  }
  file.path(root, "derived_sensitive", "nhanes", "raw_xpt")
}

markdown_table <- function(dat) {
  dat <- as.data.frame(dat, stringsAsFactors = FALSE)
  dat[] <- lapply(dat, function(x) ifelse(is.na(x), "", as.character(x)))
  header <- paste0("| ", paste(names(dat), collapse = " | "), " |")
  separator <- paste0("| ", paste(rep("---", ncol(dat)), collapse = " | "), " |")
  rows <- apply(dat, 1, function(x) paste0("| ", paste(x, collapse = " | "), " |"))
  c(header, separator, rows)
}

num <- function(x) suppressWarnings(as.numeric(x))

char_quality <- function(x) {
  out <- toupper(trimws(as.character(x)))
  out[!nzchar(out) | is.na(out)] <- NA_character_
  out
}

yes_no <- function(x) {
  z <- num(x)
  dplyr::case_when(
    z == 1 ~ 1L,
    z == 2 ~ 0L,
    TRUE ~ NA_integer_
  )
}

difficulty_binary <- function(x) {
  z <- num(x)
  dplyr::case_when(
    z == 1 ~ 0L,
    z %in% 2:5 ~ 1L,
    TRUE ~ NA_integer_
  )
}

row_any_positive <- function(mat) {
  mat <- as.matrix(mat)
  observed <- rowSums(!is.na(mat))
  positive <- rowSums(mat == 1, na.rm = TRUE)
  out <- rep(NA_integer_, nrow(mat))
  out[observed > 0] <- as.integer(positive[observed > 0] > 0)
  out
}

row_count_min_obs <- function(mat, min_observed = 1L) {
  mat <- as.matrix(mat)
  observed <- rowSums(!is.na(mat))
  positive <- rowSums(mat == 1, na.rm = TRUE)
  out <- rep(NA_integer_, nrow(mat))
  out[observed >= min_observed] <- positive[observed >= min_observed]
  out
}

num_summary <- function(x) {
  z <- num(x)
  z <- z[!is.na(z)]
  if (length(z) == 0) {
    return(tibble(
      n_nonmissing = 0L,
      min = NA_real_, p01 = NA_real_, p05 = NA_real_, median = NA_real_,
      mean = NA_real_, p95 = NA_real_, p99 = NA_real_, max = NA_real_
    ))
  }
  qs <- as.numeric(stats::quantile(z, probs = c(0.01, 0.05, 0.5, 0.95, 0.99), na.rm = TRUE, names = FALSE))
  tibble(
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

component_vars <- list(
  DEMO = c("SEQN", "SDDSRVYR", "RIAGENDR", "RIDAGEYR", "RIDRETH1", "DMDEDUC2", "INDFMPIR", "WTINT2YR", "WTMEC2YR", "SDMVPSU", "SDMVSTRA"),
  BMX = c("SEQN", "BMXHT", "BMXWT", "BMXBMI", "BMXWAIST"),
  SPX = c("SEQN", "SPXNPEF", "SPXNFEV1", "SPXNFVC", "SPXNQEFF", "SPXNQFVC", "SPXNQFV1", "SPDNACC", "SPXNSTAT", "SPDBRONC", "SPXBSTAT"),
  MCQ = c("SEQN", "MCQ010", "MCQ035", "MCQ160G", "MCQ160K", "MCQ170K"),
  SMQ = c("SEQN", "SMQ020", "SMQ040", "SMD030", "SMD055", "SMD057", "SMD650"),
  PFQ = c("SEQN", "PFQ061A", "PFQ061B", "PFQ061C", "PFQ061F", "PFQ061G", "PFQ061H", "PFQ061I", "PFQ061J", "PFQ061K", "PFQ061L"),
  DPQ = c("SEQN", paste0("DPQ0", seq(10, 90, 10))),
  PAQ = c("SEQN", "PAQ605", "PAQ620", "PAQ635", "PAQ650", "PAQ665"),
  CBC = c("SEQN", "LBXWBCSI", "LBXNEPCT", "LBXLYPCT", "LBXPLTSI"),
  GHB = c("SEQN", "LBXGH"),
  TCHOL = c("SEQN", "LBXTC"),
  HDL = c("SEQN", "LBDHDD"),
  BIOPRO = c("SEQN", "LBXSAL")
)

cycle_specs <- tibble(
  cycle = c("2007-2008", "2009-2010", "2011-2012"),
  suffix = c("E", "F", "G"),
  sddsrvyr_expected = c(5, 6, 7)
)

read_component <- function(raw_dir, component, suffix) {
  file_path <- file.path(raw_dir, paste0(component, "_", suffix, ".XPT"))
  if (!file.exists(file_path)) {
    stop("Missing NHANES XPT file: ", file_path, call. = FALSE)
  }
  data <- as.data.frame(haven::zap_labels(haven::read_xpt(file_path)), stringsAsFactors = FALSE)
  names(data) <- toupper(names(data))
  vars <- component_vars[[component]]
  present <- intersect(vars, names(data))
  out <- data[, present, drop = FALSE]
  for (v in setdiff(vars, names(out))) {
    out[[v]] <- NA
  }
  out <- out[, vars, drop = FALSE]
  out
}

merge_cycle <- function(raw_dir, suffix) {
  tables <- lapply(names(component_vars), function(component) read_component(raw_dir, component, suffix))
  names(tables) <- names(component_vars)
  Reduce(function(x, y) dplyr::full_join(x, y, by = "SEQN"), tables)
}

add_respiratory_residuals <- function(data) {
  primary_idx <- with(
    data,
    adult_45plus == 1 &
      pef_quality_abc == 1 &
      !is.na(pef_l_min) &
      !is.na(age_years) &
      !is.na(sex_code) &
      !is.na(height_cm)
  )
  model_data <- data[primary_idx, , drop = FALSE]
  if (nrow(model_data) < 100 || dplyr::n_distinct(model_data$sex_code) < 2) {
    stop("Insufficient NHANES records for primary PEF residual model.", call. = FALSE)
  }
  fit <- stats::lm(pef_l_min ~ age_years + factor(sex_code) + height_cm, data = model_data)
  pred <- rep(NA_real_, nrow(data))
  resid <- rep(NA_real_, nrow(data))
  resid_z <- rep(NA_real_, nrow(data))
  pred[primary_idx] <- stats::predict(fit, newdata = model_data)
  resid[primary_idx] <- data$pef_l_min[primary_idx] - pred[primary_idx]
  resid_sd <- stats::sd(stats::residuals(fit), na.rm = TRUE)
  resid_z[primary_idx] <- resid[primary_idx] / resid_sd

  calibration_idx <- primary_idx & !is.na(data$race_ethnicity_code)
  calibration_data <- data[calibration_idx, , drop = FALSE]
  race_pred <- rep(NA_real_, nrow(data))
  race_resid <- rep(NA_real_, nrow(data))
  race_resid_z <- rep(NA_real_, nrow(data))
  race_resid_sd <- NA_real_
  if (nrow(calibration_data) >= 100 && dplyr::n_distinct(calibration_data$race_ethnicity_code) > 1) {
    race_fit <- stats::lm(
      pef_l_min ~ age_years + factor(sex_code) + height_cm + factor(race_ethnicity_code),
      data = calibration_data
    )
    race_pred[calibration_idx] <- stats::predict(race_fit, newdata = calibration_data)
    race_resid[calibration_idx] <- data$pef_l_min[calibration_idx] - race_pred[calibration_idx]
    race_resid_sd <- stats::sd(stats::residuals(race_fit), na.rm = TRUE)
    race_resid_z[calibration_idx] <- race_resid[calibration_idx] / race_resid_sd
  }

  data$pef_pred_l_min <- pred
  data$pef_resid_l_min <- resid
  data$pef_resid_z <- resid_z
  data$resp_vulnerability_z <- -resid_z
  data$pef_pred_l_min_race_calibrated <- race_pred
  data$pef_resid_z_race_calibrated <- race_resid_z
  data$resp_vulnerability_z_race_calibrated <- -race_resid_z

  attr(data, "pef_residual_model") <- list(
    formula = "pef_l_min ~ age_years + factor(sex_code) + height_cm",
    n = nrow(model_data),
    residual_sd = resid_sd,
    coefficients = stats::coef(fit),
    race_calibrated_formula = "pef_l_min ~ age_years + factor(sex_code) + height_cm + factor(race_ethnicity_code)",
    race_calibrated_n = nrow(calibration_data),
    race_calibrated_residual_sd = race_resid_sd
  )
  data
}

build_clean_data <- function(raw_dir) {
  merged <- bind_rows(lapply(seq_len(nrow(cycle_specs)), function(i) {
    cycle_data <- merge_cycle(raw_dir, cycle_specs$suffix[[i]])
    cycle_data$cycle_label <- cycle_specs$suffix[[i]]
    cycle_data$cycle <- cycle_specs$cycle[[i]]
    cycle_data
  }))

  dpq_vars <- paste0("DPQ0", seq(10, 90, 10))
  pfq_mobility_vars <- c("PFQ061B", "PFQ061C")
  pfq_adl_iadl_vars <- c("PFQ061A", "PFQ061F", "PFQ061G", "PFQ061H", "PFQ061I", "PFQ061J", "PFQ061K", "PFQ061L")
  paq_vars <- c("PAQ605", "PAQ620", "PAQ635", "PAQ650", "PAQ665")

  clean <- merged %>%
    mutate(
      participant_id = num(.data$SEQN),
      sddsrvyr = num(.data$SDDSRVYR),
      cycle = dplyr::case_when(
        .data$sddsrvyr == 5 ~ "2007-2008",
        .data$sddsrvyr == 6 ~ "2009-2010",
        .data$sddsrvyr == 7 ~ "2011-2012",
        TRUE ~ .data$cycle
      ),
      sex_code = num(.data$RIAGENDR),
      sex = dplyr::case_when(.data$sex_code == 1 ~ "male", .data$sex_code == 2 ~ "female", TRUE ~ NA_character_),
      age_years = num(.data$RIDAGEYR),
      adult_45plus = as.integer(!is.na(.data$age_years) & .data$age_years >= 45),
      race_ethnicity_code = num(.data$RIDRETH1),
      race_ethnicity = dplyr::case_when(
        .data$race_ethnicity_code == 1 ~ "Mexican American",
        .data$race_ethnicity_code == 2 ~ "Other Hispanic",
        .data$race_ethnicity_code == 3 ~ "Non-Hispanic White",
        .data$race_ethnicity_code == 4 ~ "Non-Hispanic Black",
        .data$race_ethnicity_code == 5 ~ "Other/Multi-Racial",
        TRUE ~ NA_character_
      ),
      education_code = num(.data$DMDEDUC2),
      education = dplyr::case_when(
        .data$education_code == 1 ~ "<9th grade",
        .data$education_code == 2 ~ "9-11th grade",
        .data$education_code == 3 ~ "High school/GED",
        .data$education_code == 4 ~ "Some college/AA",
        .data$education_code == 5 ~ "College graduate+",
        TRUE ~ NA_character_
      ),
      income_poverty_ratio = num(.data$INDFMPIR),
      wtint2yr = num(.data$WTINT2YR),
      wtmec2yr = num(.data$WTMEC2YR),
      wtint6yr = .data$wtint2yr / 3,
      wtmec6yr = .data$wtmec2yr / 3,
      psu = num(.data$SDMVPSU),
      strata = num(.data$SDMVSTRA),
      height_cm = num(.data$BMXHT),
      weight_kg = num(.data$BMXWT),
      bmi = num(.data$BMXBMI),
      waist_cm = num(.data$BMXWAIST),
      pef_ml_sec = num(.data$SPXNPEF),
      pef_l_min = .data$pef_ml_sec * 0.06,
      fev1_ml = num(.data$SPXNFEV1),
      fvc_ml = num(.data$SPXNFVC),
      fev1_l = .data$fev1_ml / 1000,
      fvc_l = .data$fvc_ml / 1000,
      fev1_fvc = dplyr::if_else(!is.na(.data$fev1_ml) & !is.na(.data$fvc_ml) & .data$fvc_ml > 0, .data$fev1_ml / .data$fvc_ml, NA_real_),
      pef_quality = char_quality(.data$SPXNQEFF),
      fvc_quality = char_quality(.data$SPXNQFVC),
      fev1_quality = char_quality(.data$SPXNQFV1),
      acceptable_curves = num(.data$SPDNACC),
      spirometry_status = num(.data$SPXNSTAT),
      pef_quality_abc = as.integer(!is.na(.data$pef_l_min) & .data$pef_quality %in% c("A", "B", "C")),
      pef_quality_ab = as.integer(!is.na(.data$pef_l_min) & .data$pef_quality %in% c("A", "B")),
      spirometry_quality_abc = as.integer(
        !is.na(.data$fev1_fvc) &
          .data$fev1_quality %in% c("A", "B", "C") &
          .data$fvc_quality %in% c("A", "B", "C")
      ),
      spirometry_quality_ab = as.integer(
        !is.na(.data$fev1_fvc) &
          .data$fev1_quality %in% c("A", "B") &
          .data$fvc_quality %in% c("A", "B")
      ),
      obstruction_fixed_ratio_abc = dplyr::if_else(.data$spirometry_quality_abc == 1, as.integer(.data$fev1_fvc < 0.70), NA_integer_),
      obstruction_fixed_ratio_ab = dplyr::if_else(.data$spirometry_quality_ab == 1, as.integer(.data$fev1_fvc < 0.70), NA_integer_),
      ever_asthma = yes_no(.data$MCQ010),
      current_asthma = dplyr::case_when(
        .data$ever_asthma == 0 ~ 0L,
        .data$ever_asthma == 1 & num(.data$MCQ035) == 1 ~ 1L,
        .data$ever_asthma == 1 & num(.data$MCQ035) == 2 ~ 0L,
        TRUE ~ NA_integer_
      ),
      ever_emphysema = yes_no(.data$MCQ160G),
      ever_chronic_bronchitis = yes_no(.data$MCQ160K),
      current_chronic_bronchitis = dplyr::case_when(
        .data$ever_chronic_bronchitis == 0 ~ 0L,
        .data$ever_chronic_bronchitis == 1 & num(.data$MCQ170K) == 1 ~ 1L,
        .data$ever_chronic_bronchitis == 1 & num(.data$MCQ170K) == 2 ~ 0L,
        TRUE ~ NA_integer_
      ),
      self_reported_emphysema_or_bronchitis = dplyr::case_when(
        .data$ever_emphysema == 1 | .data$ever_chronic_bronchitis == 1 ~ 1L,
        .data$ever_emphysema == 0 & .data$ever_chronic_bronchitis == 0 ~ 0L,
        TRUE ~ NA_integer_
      ),
      smoke_ever = yes_no(.data$SMQ020),
      smoke_current = dplyr::case_when(
        .data$smoke_ever == 1 & num(.data$SMQ040) %in% c(1, 2) ~ 1L,
        .data$smoke_ever == 1 & num(.data$SMQ040) == 3 ~ 0L,
        .data$smoke_ever == 0 ~ 0L,
        TRUE ~ NA_integer_
      ),
      smoke_former = dplyr::case_when(
        .data$smoke_ever == 1 & num(.data$SMQ040) == 3 ~ 1L,
        .data$smoke_ever == 1 & num(.data$SMQ040) %in% c(1, 2) ~ 0L,
        .data$smoke_ever == 0 ~ 0L,
        TRUE ~ NA_integer_
      ),
      smoking_status = dplyr::case_when(
        .data$smoke_ever == 0 ~ "never",
        .data$smoke_current == 1 ~ "current",
        .data$smoke_former == 1 ~ "former",
        TRUE ~ NA_character_
      ),
      cigarettes_per_day_past_30 = dplyr::if_else(num(.data$SMD650) %in% 777:999, NA_real_, num(.data$SMD650)),
      wbc_1000_ul = num(.data$LBXWBCSI),
      neutrophil_percent = num(.data$LBXNEPCT),
      lymphocyte_percent = num(.data$LBXLYPCT),
      platelet_1000_ul = num(.data$LBXPLTSI),
      neutrophil_1000_ul = .data$wbc_1000_ul * .data$neutrophil_percent / 100,
      lymphocyte_1000_ul = .data$wbc_1000_ul * .data$lymphocyte_percent / 100,
      nlr = .data$neutrophil_1000_ul / .data$lymphocyte_1000_ul,
      sii = .data$platelet_1000_ul * .data$neutrophil_1000_ul / .data$lymphocyte_1000_ul,
      hba1c_percent = num(.data$LBXGH),
      total_cholesterol_mg_dl = num(.data$LBXTC),
      hdl_cholesterol_mg_dl = num(.data$LBDHDD),
      albumin_g_dl = num(.data$LBXSAL)
    )

  clean[dpq_vars] <- lapply(clean[dpq_vars], function(x) dplyr::if_else(num(x) %in% 0:3, num(x), NA_real_))
  clean[paste0("difficulty_", pfq_mobility_vars)] <- lapply(clean[pfq_mobility_vars], difficulty_binary)
  clean[paste0("difficulty_", pfq_adl_iadl_vars)] <- lapply(clean[pfq_adl_iadl_vars], difficulty_binary)
  clean[paste0("active_", paq_vars)] <- lapply(clean[paq_vars], yes_no)

  clean$mobility_difficulty <- row_any_positive(clean[paste0("difficulty_", pfq_mobility_vars)])
  clean$adl_iadl_difficulty <- row_any_positive(clean[paste0("difficulty_", pfq_adl_iadl_vars)])
  clean$phq9_score <- rowSums(clean[dpq_vars], na.rm = FALSE)
  clean$depressive_symptoms_phq9_ge10 <- dplyr::if_else(!is.na(clean$phq9_score), as.integer(clean$phq9_score >= 10), NA_integer_)
  active_mat <- as.matrix(clean[paste0("active_", paq_vars)])
  clean$physical_activity_positive_count <- rowSums(active_mat == 1, na.rm = TRUE)
  clean$physical_activity_observed_count <- rowSums(!is.na(active_mat))
  clean$low_physical_activity <- dplyr::if_else(
    clean$physical_activity_observed_count == length(paq_vars),
    as.integer(clean$physical_activity_positive_count == 0),
    NA_integer_
  )
  frailty_mat <- clean[c("mobility_difficulty", "adl_iadl_difficulty", "depressive_symptoms_phq9_ge10", "low_physical_activity")]
  clean$nhanes_frailty_proxy_observed_count <- rowSums(!is.na(frailty_mat))
  clean$nhanes_frailty_proxy_count <- row_count_min_obs(frailty_mat, min_observed = 3L)
  clean$nhanes_frailty_proxy_ge2 <- dplyr::if_else(!is.na(clean$nhanes_frailty_proxy_count), as.integer(clean$nhanes_frailty_proxy_count >= 2), NA_integer_)

  clean <- add_respiratory_residuals(clean)

  selected <- clean %>%
    select(
      participant_id, cycle, cycle_label, sddsrvyr,
      age_years, adult_45plus, sex_code, sex, race_ethnicity_code, race_ethnicity,
      education_code, education, income_poverty_ratio,
      wtint2yr, wtmec2yr, wtint6yr, wtmec6yr, psu, strata,
      height_cm, weight_kg, bmi, waist_cm,
      pef_ml_sec, pef_l_min, pef_quality, pef_quality_abc, pef_quality_ab,
      fev1_ml, fvc_ml, fev1_l, fvc_l, fev1_fvc, fev1_quality, fvc_quality,
      acceptable_curves, spirometry_status, spirometry_quality_abc, spirometry_quality_ab,
      obstruction_fixed_ratio_abc, obstruction_fixed_ratio_ab,
      pef_pred_l_min, pef_resid_l_min, pef_resid_z, resp_vulnerability_z,
      pef_pred_l_min_race_calibrated, pef_resid_z_race_calibrated, resp_vulnerability_z_race_calibrated,
      ever_asthma, current_asthma, ever_emphysema, ever_chronic_bronchitis,
      current_chronic_bronchitis, self_reported_emphysema_or_bronchitis,
      smoke_ever, smoke_current, smoke_former, smoking_status, cigarettes_per_day_past_30,
      mobility_difficulty, adl_iadl_difficulty, phq9_score, depressive_symptoms_phq9_ge10,
      physical_activity_positive_count, low_physical_activity,
      nhanes_frailty_proxy_observed_count, nhanes_frailty_proxy_count, nhanes_frailty_proxy_ge2,
      wbc_1000_ul, neutrophil_percent, lymphocyte_percent, platelet_1000_ul,
      neutrophil_1000_ul, lymphocyte_1000_ul, nlr, sii,
      hba1c_percent, total_cholesterol_mg_dl, hdl_cholesterol_mg_dl, albumin_g_dl
    )
  attr(selected, "pef_residual_model") <- attr(clean, "pef_residual_model")
  selected
}

main <- function() {
  root <- find_project_root()
  raw_dir <- local_raw_dir(root)
  if (!dir.exists(raw_dir)) {
    stop("NHANES raw XPT directory does not exist: ", raw_dir, call. = FALSE)
  }

  local_dir <- file.path(root, "derived_sensitive", "nhanes")
  table_dir <- file.path(root, "results", "tables")
  log_dir <- file.path(root, "results", "logs")
  dir.create(local_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

  message("Building NHANES V0.2 analysis-ready dataset from: ", raw_dir)
  data <- build_clean_data(raw_dir)
  model_info <- attr(data, "pef_residual_model")

  rds_path <- file.path(local_dir, "nhanes_replication_v0_2_analysis_ready.rds")
  saveRDS(data, rds_path)

  component_row_counts <- bind_rows(lapply(names(component_vars), function(component) {
    bind_rows(lapply(cycle_specs$suffix, function(suffix) {
      file_path <- file.path(raw_dir, paste0(component, "_", suffix, ".XPT"))
      tibble(
        component = component,
        cycle_label = suffix,
        rows = nrow(haven::read_xpt(file_path, col_select = "SEQN")),
        cols_selected = length(component_vars[[component]])
      )
    }))
  }))

  cohort_counts <- tibble(
    metric = c(
      "merged_person_rows",
      "adult_45plus_rows",
      "adult_45plus_positive_mec_weight",
      "adult_45plus_pef_available",
      "adult_45plus_pef_quality_abc",
      "primary_pef_residual_model_n",
      "adult_45plus_obstruction_abc_observed",
      "adult_45plus_obstruction_abc_yes",
      "adult_45plus_current_asthma_observed",
      "adult_45plus_current_asthma_yes",
      "adult_45plus_emphysema_or_bronchitis_observed",
      "adult_45plus_emphysema_or_bronchitis_yes",
      "adult_45plus_frailty_proxy_count_observed",
      "adult_45plus_frailty_proxy_ge2_yes"
    ),
    value = c(
      nrow(data),
      sum(data$adult_45plus == 1, na.rm = TRUE),
      sum(data$adult_45plus == 1 & !is.na(data$wtmec6yr) & data$wtmec6yr > 0, na.rm = TRUE),
      sum(data$adult_45plus == 1 & !is.na(data$pef_l_min), na.rm = TRUE),
      sum(data$adult_45plus == 1 & data$pef_quality_abc == 1, na.rm = TRUE),
      model_info$n,
      sum(data$adult_45plus == 1 & !is.na(data$obstruction_fixed_ratio_abc), na.rm = TRUE),
      sum(data$adult_45plus == 1 & data$obstruction_fixed_ratio_abc == 1, na.rm = TRUE),
      sum(data$adult_45plus == 1 & !is.na(data$current_asthma), na.rm = TRUE),
      sum(data$adult_45plus == 1 & data$current_asthma == 1, na.rm = TRUE),
      sum(data$adult_45plus == 1 & !is.na(data$self_reported_emphysema_or_bronchitis), na.rm = TRUE),
      sum(data$adult_45plus == 1 & data$self_reported_emphysema_or_bronchitis == 1, na.rm = TRUE),
      sum(data$adult_45plus == 1 & !is.na(data$nhanes_frailty_proxy_count), na.rm = TRUE),
      sum(data$adult_45plus == 1 & data$nhanes_frailty_proxy_ge2 == 1, na.rm = TRUE)
    )
  )

  adult <- data %>% filter(.data$adult_45plus == 1)
  selected_missing <- c(
    "wtmec6yr", "psu", "strata", "age_years", "sex", "race_ethnicity",
    "height_cm", "bmi", "smoke_ever", "pef_l_min", "resp_vulnerability_z",
    "fev1_fvc", "obstruction_fixed_ratio_abc", "current_asthma",
    "self_reported_emphysema_or_bronchitis", "nhanes_frailty_proxy_count",
    "hba1c_percent", "wbc_1000_ul", "albumin_g_dl"
  )
  missingness <- tibble(variable = selected_missing) %>%
    rowwise() %>%
    mutate(
      population = "adult_45plus",
      n = nrow(adult),
      n_missing = sum(is.na(adult[[.data$variable]])),
      pct_missing = round(100 * .data$n_missing / .data$n, 2)
    ) %>%
    ungroup()

  numeric_vars <- c(
    "age_years", "height_cm", "bmi", "pef_l_min", "resp_vulnerability_z",
    "fev1_l", "fvc_l", "fev1_fvc", "phq9_score", "nhanes_frailty_proxy_count",
    "hba1c_percent", "wbc_1000_ul", "nlr", "sii", "albumin_g_dl"
  )
  numeric_summary <- bind_rows(lapply(numeric_vars, function(v) {
    num_summary(adult[[v]]) %>% mutate(variable = v, .before = 1)
  }))

  quality_counts <- bind_rows(lapply(c("pef_quality", "fev1_quality", "fvc_quality"), function(v) {
    adult %>%
      count(value = .data[[v]], name = "n") %>%
      mutate(variable = v, value = ifelse(is.na(.data$value), "missing", .data$value), .before = 1)
  }))

  outcome_counts <- adult %>%
    group_by(.data$cycle) %>%
    summarise(
      n = n(),
      pef_quality_abc_n = sum(.data$pef_quality_abc == 1, na.rm = TRUE),
      residual_model_n = sum(!is.na(.data$resp_vulnerability_z)),
      obstruction_abc_observed = sum(!is.na(.data$obstruction_fixed_ratio_abc)),
      obstruction_abc_yes = sum(.data$obstruction_fixed_ratio_abc == 1, na.rm = TRUE),
      current_asthma_observed = sum(!is.na(.data$current_asthma)),
      current_asthma_yes = sum(.data$current_asthma == 1, na.rm = TRUE),
      emphysema_or_bronchitis_observed = sum(!is.na(.data$self_reported_emphysema_or_bronchitis)),
      emphysema_or_bronchitis_yes = sum(.data$self_reported_emphysema_or_bronchitis == 1, na.rm = TRUE),
      frailty_count_observed = sum(!is.na(.data$nhanes_frailty_proxy_count)),
      frailty_ge2_yes = sum(.data$nhanes_frailty_proxy_ge2 == 1, na.rm = TRUE),
      .groups = "drop"
    )

  readr::write_csv(component_row_counts, file.path(table_dir, "nhanes_v0_2_component_row_counts.csv"))
  readr::write_csv(cohort_counts, file.path(table_dir, "nhanes_v0_2_cohort_counts.csv"))
  readr::write_csv(missingness, file.path(table_dir, "nhanes_v0_2_missingness.csv"))
  readr::write_csv(numeric_summary, file.path(table_dir, "nhanes_v0_2_numeric_summary.csv"))
  readr::write_csv(quality_counts, file.path(table_dir, "nhanes_v0_2_resp_quality_counts.csv"))
  readr::write_csv(outcome_counts, file.path(table_dir, "nhanes_v0_2_outcome_counts_by_cycle.csv"))

  log_lines <- c(
    "# NHANES V0.2 Cleaning Log",
    "",
    paste0("- Run date: ", Sys.Date()),
    "- Source: local public NHANES XPT files under `derived_sensitive/nhanes/raw_xpt`.",
    "- Row-level output remains local-only under `derived_sensitive/nhanes/`.",
    "",
    "## Cohort Counts",
    "",
    markdown_table(cohort_counts),
    "",
    "## PEF Residual Model",
    "",
    paste0("- Formula: `", model_info$formula, "`"),
    paste0("- Model N: ", model_info$n),
    paste0("- Residual SD: ", round(model_info$residual_sd, 4)),
    paste0("- Race-calibrated sensitivity N: ", model_info$race_calibrated_n),
    "",
    "## Primary Analysis Rules",
    "",
    "- Primary exposure: `SPXNPEF * 0.06` converted to L/min and residualized against age, sex, and height.",
    "- Higher `resp_vulnerability_z` means lower-than-expected PEF reserve.",
    "- Primary objective outcome: `SPXNFEV1 / SPXNFVC < 0.70` among A/B/C FEV1 and FVC quality grades.",
    "- Survey design variables are carried as `wtmec6yr`, `psu`, and `strata`; modeling is a later script.",
    "- Frailty remains a proxy from mobility/ADL-IADL difficulty, PHQ-9 >=10, and low physical activity.",
    "",
    "## Files Written",
    "",
    "- `derived_sensitive/nhanes/nhanes_replication_v0_2_analysis_ready.rds`",
    "- `results/tables/nhanes_v0_2_component_row_counts.csv`",
    "- `results/tables/nhanes_v0_2_cohort_counts.csv`",
    "- `results/tables/nhanes_v0_2_missingness.csv`",
    "- `results/tables/nhanes_v0_2_numeric_summary.csv`",
    "- `results/tables/nhanes_v0_2_resp_quality_counts.csv`",
    "- `results/tables/nhanes_v0_2_outcome_counts_by_cycle.csv`"
  )
  log_path <- file.path(log_dir, "nhanes_v0_2_cleaning_log.md")
  writeLines(log_lines, log_path)

  message("Wrote local-only NHANES analysis RDS: ", rds_path)
  message("Wrote aggregate cleaning log: ", log_path)
}

if (sys.nframe() == 0) {
  main()
}
