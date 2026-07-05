# NHANES v0.5 GLI/LLN/PRISm phenotype lock.
# Adds GLI Global 2022 and GLI-2012 spirometry reference phenotypes to the
# local-only NHANES analysis dataset and writes aggregate lock/QC outputs.

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(rspiro)
  library(stringr)
  library(survey)
  library(tibble)
  library(tidyr)
})

options(survey.lonely.psu = "adjust")

find_project_root <- function() {
  script_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", script_args, value = TRUE)
  if (length(file_arg) > 0) {
    script_path <- sub("^--file=", "", file_arg[[1]])
    return(normalizePath(file.path(dirname(script_path), "..", ".."), mustWork = FALSE))
  }
  normalizePath(getwd(), mustWork = FALSE)
}

format_p <- function(p) {
  ifelse(is.na(p), NA_character_, ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
}

format_ci <- function(est, lo, hi, digits = 2) {
  sprintf(
    paste0("%.", digits, "f (%.", digits, "f-%.", digits, "f)"),
    est, lo, hi
  )
}

markdown_table <- function(dat) {
  dat <- as.data.frame(dat, stringsAsFactors = FALSE)
  dat[] <- lapply(dat, function(x) ifelse(is.na(x), "", as.character(x)))
  header <- paste0("| ", paste(names(dat), collapse = " | "), " |")
  separator <- paste0("| ", paste(rep("---", ncol(dat)), collapse = " | "), " |")
  rows <- apply(dat, 1, function(x) paste0("| ", paste(x, collapse = " | "), " |"))
  c(header, separator, rows)
}

weighted_quantile <- function(x, w, probs = c(0.25, 0.5, 0.75)) {
  ok <- !is.na(x) & !is.na(w) & w > 0
  x <- x[ok]
  w <- w[ok]
  if (length(x) == 0) {
    return(rep(NA_real_, length(probs)))
  }
  ord <- order(x)
  x <- x[ord]
  w <- w[ord]
  cw <- cumsum(w) / sum(w)
  vapply(probs, function(p) x[which(cw >= p)[1]], numeric(1))
}

make_design <- function(data, weight_var = "wtmec6yr") {
  survey::svydesign(
    ids = ~psu,
    strata = ~strata,
    weights = stats::as.formula(paste0("~", weight_var)),
    nest = TRUE,
    data = data
  )
}

clean_design_data <- function(data, weight_var = "wtmec6yr") {
  data %>%
    filter(
      !is.na(.data[[weight_var]]),
      .data[[weight_var]] > 0,
      !is.na(.data$psu),
      !is.na(.data$strata)
    )
}

calc_prevalence <- function(data, variable, label, weight_var = "wtmec6yr") {
  df <- data %>%
    filter(.data$adult_45plus == 1, !is.na(.data[[variable]])) %>%
    clean_design_data(weight_var) %>%
    mutate(.indicator = as.numeric(.data[[variable]] == 1))

  if (nrow(df) < 2) {
    return(tibble(
      variable = variable,
      phenotype_label = label,
      n_unweighted = nrow(df),
      events_unweighted = sum(df$.indicator == 1, na.rm = TRUE),
      weighted_prevalence = NA_real_,
      se = NA_real_,
      ci_low = NA_real_,
      ci_high = NA_real_
    ))
  }

  design <- make_design(df, weight_var)
  est <- survey::svymean(~.indicator, design = design, na.rm = TRUE)
  prop <- as.numeric(stats::coef(est)[1])
  se_value <- as.numeric(survey::SE(est)[1])
  ci <- pmin(pmax(prop + c(-1, 1) * stats::qnorm(0.975) * se_value, 0), 1)
  tibble(
    variable = variable,
    phenotype_label = label,
    n_unweighted = nrow(df),
    events_unweighted = sum(df$.indicator == 1, na.rm = TRUE),
    weighted_prevalence = prop,
    se = se_value,
    ci_low = ci[1],
    ci_high = ci[2]
  )
}

empty_numeric <- function(n) rep(NA_real_, n)
empty_integer <- function(n) rep(NA_integer_, n)

fill_gli_global <- function(data) {
  n <- nrow(data)
  idx_ref <- with(
    data,
    adult_45plus == 1 &
      !is.na(age_years) & age_years >= 3 & age_years <= 95 &
      !is.na(height_cm) & height_cm > 0 &
      sex_code %in% c(1, 2)
  )
  age <- data$age_years[idx_ref]
  height_m <- data$height_cm[idx_ref] / 100
  gender <- data$sex_code[idx_ref]

  pred <- rspiro::pred_GLIgl(age, height_m, gender, param = c("FEV1", "FVC", "FEV1FVC"))
  lln <- rspiro::LLN_GLIgl(age, height_m, gender, param = c("FEV1", "FVC", "FEV1FVC"))

  data$gli_global2022_pred_fev1_l <- empty_numeric(n)
  data$gli_global2022_pred_fvc_l <- empty_numeric(n)
  data$gli_global2022_pred_fev1_fvc <- empty_numeric(n)
  data$gli_global2022_lln_fev1_l <- empty_numeric(n)
  data$gli_global2022_lln_fvc_l <- empty_numeric(n)
  data$gli_global2022_lln_fev1_fvc <- empty_numeric(n)
  data$gli_global2022_pred_fev1_l[idx_ref] <- pred$pred.FEV1
  data$gli_global2022_pred_fvc_l[idx_ref] <- pred$pred.FVC
  data$gli_global2022_pred_fev1_fvc[idx_ref] <- pred$pred.FEV1FVC
  data$gli_global2022_lln_fev1_l[idx_ref] <- lln$LLN.FEV1
  data$gli_global2022_lln_fvc_l[idx_ref] <- lln$LLN.FVC
  data$gli_global2022_lln_fev1_fvc[idx_ref] <- lln$LLN.FEV1FVC

  idx_obs <- idx_ref & !is.na(data$fev1_l) & !is.na(data$fvc_l) & !is.na(data$fev1_fvc)
  z <- rspiro::zscore_GLIgl(
    age = data$age_years[idx_obs],
    height = data$height_cm[idx_obs] / 100,
    gender = data$sex_code[idx_obs],
    FEV1 = data$fev1_l[idx_obs],
    FVC = data$fvc_l[idx_obs],
    FEV1FVC = data$fev1_fvc[idx_obs]
  )
  pct <- rspiro::pctpred_GLIgl(
    age = data$age_years[idx_obs],
    height = data$height_cm[idx_obs] / 100,
    gender = data$sex_code[idx_obs],
    FEV1 = data$fev1_l[idx_obs],
    FVC = data$fvc_l[idx_obs],
    FEV1FVC = data$fev1_fvc[idx_obs]
  )

  data$gli_global2022_z_fev1 <- empty_numeric(n)
  data$gli_global2022_z_fvc <- empty_numeric(n)
  data$gli_global2022_z_fev1_fvc <- empty_numeric(n)
  data$gli_global2022_ppred_fev1 <- empty_numeric(n)
  data$gli_global2022_ppred_fvc <- empty_numeric(n)
  data$gli_global2022_ppred_fev1_fvc <- empty_numeric(n)
  data$gli_global2022_z_fev1[idx_obs] <- z$z.score.FEV1
  data$gli_global2022_z_fvc[idx_obs] <- z$z.score.FVC
  data$gli_global2022_z_fev1_fvc[idx_obs] <- z$z.score.FEV1FVC
  data$gli_global2022_ppred_fev1[idx_obs] <- pct$pctpred.FEV1
  data$gli_global2022_ppred_fvc[idx_obs] <- pct$pctpred.FVC
  data$gli_global2022_ppred_fev1_fvc[idx_obs] <- pct$pctpred.FEV1FVC
  data
}

fill_gli2012 <- function(data) {
  n <- nrow(data)
  data$gli2012_ethnicity_code <- dplyr::case_when(
    data$race_ethnicity_code == 3 ~ 1, # Caucasian
    data$race_ethnicity_code == 4 ~ 2, # African-American
    data$race_ethnicity_code %in% c(1, 2, 5) ~ 5, # Other/mixed
    TRUE ~ NA_real_
  )
  data$gli2012_ethnicity_label <- dplyr::case_when(
    data$gli2012_ethnicity_code == 1 ~ "Caucasian",
    data$gli2012_ethnicity_code == 2 ~ "African-American",
    data$gli2012_ethnicity_code == 5 ~ "Other/mixed",
    TRUE ~ NA_character_
  )

  idx_ref <- with(
    data,
    adult_45plus == 1 &
      !is.na(age_years) & age_years >= 3 & age_years <= 95 &
      !is.na(height_cm) & height_cm > 0 &
      sex_code %in% c(1, 2) &
      !is.na(gli2012_ethnicity_code)
  )
  age <- data$age_years[idx_ref]
  height_m <- data$height_cm[idx_ref] / 100
  gender <- data$sex_code[idx_ref]
  ethnicity <- data$gli2012_ethnicity_code[idx_ref]

  pred <- rspiro::pred_GLI(age, height_m, gender, ethnicity, param = c("FEV1", "FVC", "FEV1FVC"))
  lln <- rspiro::LLN_GLI(age, height_m, gender, ethnicity, param = c("FEV1", "FVC", "FEV1FVC"))

  data$gli2012_pred_fev1_l <- empty_numeric(n)
  data$gli2012_pred_fvc_l <- empty_numeric(n)
  data$gli2012_pred_fev1_fvc <- empty_numeric(n)
  data$gli2012_lln_fev1_l <- empty_numeric(n)
  data$gli2012_lln_fvc_l <- empty_numeric(n)
  data$gli2012_lln_fev1_fvc <- empty_numeric(n)
  data$gli2012_pred_fev1_l[idx_ref] <- pred$pred.FEV1
  data$gli2012_pred_fvc_l[idx_ref] <- pred$pred.FVC
  data$gli2012_pred_fev1_fvc[idx_ref] <- pred$pred.FEV1FVC
  data$gli2012_lln_fev1_l[idx_ref] <- lln$LLN.FEV1
  data$gli2012_lln_fvc_l[idx_ref] <- lln$LLN.FVC
  data$gli2012_lln_fev1_fvc[idx_ref] <- lln$LLN.FEV1FVC

  idx_obs <- idx_ref & !is.na(data$fev1_l) & !is.na(data$fvc_l) & !is.na(data$fev1_fvc)
  z <- rspiro::zscore_GLI(
    age = data$age_years[idx_obs],
    height = data$height_cm[idx_obs] / 100,
    gender = data$sex_code[idx_obs],
    ethnicity = data$gli2012_ethnicity_code[idx_obs],
    FEV1 = data$fev1_l[idx_obs],
    FVC = data$fvc_l[idx_obs],
    FEV1FVC = data$fev1_fvc[idx_obs]
  )
  pct <- rspiro::pctpred_GLI(
    age = data$age_years[idx_obs],
    height = data$height_cm[idx_obs] / 100,
    gender = data$sex_code[idx_obs],
    ethnicity = data$gli2012_ethnicity_code[idx_obs],
    FEV1 = data$fev1_l[idx_obs],
    FVC = data$fvc_l[idx_obs],
    FEV1FVC = data$fev1_fvc[idx_obs]
  )

  data$gli2012_z_fev1 <- empty_numeric(n)
  data$gli2012_z_fvc <- empty_numeric(n)
  data$gli2012_z_fev1_fvc <- empty_numeric(n)
  data$gli2012_ppred_fev1 <- empty_numeric(n)
  data$gli2012_ppred_fvc <- empty_numeric(n)
  data$gli2012_ppred_fev1_fvc <- empty_numeric(n)
  data$gli2012_z_fev1[idx_obs] <- z$z.score.FEV1
  data$gli2012_z_fvc[idx_obs] <- z$z.score.FVC
  data$gli2012_z_fev1_fvc[idx_obs] <- z$z.score.FEV1FVC
  data$gli2012_ppred_fev1[idx_obs] <- pct$pctpred.FEV1
  data$gli2012_ppred_fvc[idx_obs] <- pct$pctpred.FVC
  data$gli2012_ppred_fev1_fvc[idx_obs] <- pct$pctpred.FEV1FVC
  data
}

add_gli_phenotypes <- function(data) {
  data <- fill_gli_global(data)
  data <- fill_gli2012(data)
  z_cut <- -1.64485362695147

  data %>%
    mutate(
      obstruction_gli_global2022_abc = dplyr::if_else(
        .data$spirometry_quality_abc == 1 & !is.na(.data$gli_global2022_z_fev1_fvc),
        as.integer(.data$gli_global2022_z_fev1_fvc < z_cut),
        NA_integer_
      ),
      obstruction_gli_global2022_ab = dplyr::if_else(
        .data$spirometry_quality_ab == 1 & !is.na(.data$gli_global2022_z_fev1_fvc),
        as.integer(.data$gli_global2022_z_fev1_fvc < z_cut),
        NA_integer_
      ),
      low_fev1_gli_global2022_abc = dplyr::if_else(
        .data$spirometry_quality_abc == 1 & !is.na(.data$gli_global2022_z_fev1),
        as.integer(.data$gli_global2022_z_fev1 < z_cut),
        NA_integer_
      ),
      low_fvc_gli_global2022_abc = dplyr::if_else(
        .data$spirometry_quality_abc == 1 & !is.na(.data$gli_global2022_z_fvc),
        as.integer(.data$gli_global2022_z_fvc < z_cut),
        NA_integer_
      ),
      prism_gli_global2022_z_abc = dplyr::if_else(
        .data$spirometry_quality_abc == 1 &
          !is.na(.data$gli_global2022_z_fev1_fvc) &
          !is.na(.data$gli_global2022_z_fev1),
        as.integer(.data$gli_global2022_z_fev1_fvc >= z_cut & .data$gli_global2022_z_fev1 < z_cut),
        NA_integer_
      ),
      prism_gli_global2022_80_abc = dplyr::if_else(
        .data$spirometry_quality_abc == 1 &
          !is.na(.data$fev1_fvc) &
          !is.na(.data$gli_global2022_ppred_fev1),
        as.integer(.data$fev1_fvc >= 0.70 & .data$gli_global2022_ppred_fev1 < 80),
        NA_integer_
      ),
      obstruction_gli2012_abc = dplyr::if_else(
        .data$spirometry_quality_abc == 1 & !is.na(.data$gli2012_z_fev1_fvc),
        as.integer(.data$gli2012_z_fev1_fvc < z_cut),
        NA_integer_
      ),
      obstruction_gli2012_ab = dplyr::if_else(
        .data$spirometry_quality_ab == 1 & !is.na(.data$gli2012_z_fev1_fvc),
        as.integer(.data$gli2012_z_fev1_fvc < z_cut),
        NA_integer_
      ),
      low_fev1_gli2012_abc = dplyr::if_else(
        .data$spirometry_quality_abc == 1 & !is.na(.data$gli2012_z_fev1),
        as.integer(.data$gli2012_z_fev1 < z_cut),
        NA_integer_
      ),
      low_fvc_gli2012_abc = dplyr::if_else(
        .data$spirometry_quality_abc == 1 & !is.na(.data$gli2012_z_fvc),
        as.integer(.data$gli2012_z_fvc < z_cut),
        NA_integer_
      ),
      prism_gli2012_z_abc = dplyr::if_else(
        .data$spirometry_quality_abc == 1 &
          !is.na(.data$gli2012_z_fev1_fvc) &
          !is.na(.data$gli2012_z_fev1),
        as.integer(.data$gli2012_z_fev1_fvc >= z_cut & .data$gli2012_z_fev1 < z_cut),
        NA_integer_
      ),
      prism_gli2012_80_abc = dplyr::if_else(
        .data$spirometry_quality_abc == 1 &
          !is.na(.data$fev1_fvc) &
          !is.na(.data$gli2012_ppred_fev1),
        as.integer(.data$fev1_fvc >= 0.70 & .data$gli2012_ppred_fev1 < 80),
        NA_integer_
      )
    )
}

prepare_analysis_data <- function(data) {
  base <- data %>%
    filter(.data$adult_45plus == 1) %>%
    filter(!is.na(.data$resp_vulnerability_z)) %>%
    clean_design_data("wtmec6yr")

  cuts <- weighted_quantile(base$resp_vulnerability_z, base$wtmec6yr, c(0.25, 0.5, 0.75))
  if (any(is.na(cuts)) || length(unique(cuts)) < length(cuts)) {
    stop("Unable to compute distinct weighted respiratory vulnerability quartiles.", call. = FALSE)
  }

  data_q <- data %>%
    mutate(
      rv_quartile = cut(
        .data$resp_vulnerability_z,
        breaks = c(-Inf, cuts, Inf),
        labels = c("Q1", "Q2", "Q3", "Q4"),
        include.lowest = TRUE,
        right = TRUE
      ),
      rv_quartile_num = as.numeric(.data$rv_quartile),
      rv_q4_vs_q1 = dplyr::case_when(
        .data$rv_quartile == "Q4" ~ 1,
        .data$rv_quartile == "Q1" ~ 0,
        TRUE ~ NA_real_
      )
    )

  list(
    data = data_q,
    cutpoints = tibble(
      quantile = c("weighted_p25", "weighted_p50", "weighted_p75"),
      resp_vulnerability_z = as.numeric(cuts)
    )
  )
}

phenotype_specs <- tibble::tribble(
  ~variable, ~phenotype_label, ~phenotype_family, ~equation, ~quality_scope, ~lock_role,
  "obstruction_fixed_ratio_abc", "Fixed-ratio obstruction", "obstruction", "Fixed FEV1/FVC <0.70", "A/B/C", "legacy_anchor",
  "obstruction_fixed_ratio_ab", "Fixed-ratio obstruction, strict A/B", "obstruction", "Fixed FEV1/FVC <0.70", "A/B", "quality_sensitivity",
  "obstruction_gli_global2022_abc", "GLI Global 2022 obstruction", "obstruction", "GLI Global 2022 race-neutral LLN", "A/B/C", "primary_lln",
  "obstruction_gli_global2022_ab", "GLI Global 2022 obstruction, strict A/B", "obstruction", "GLI Global 2022 race-neutral LLN", "A/B", "quality_sensitivity",
  "obstruction_gli2012_abc", "GLI-2012 obstruction", "obstruction", "GLI-2012 race/ethnicity-specific LLN", "A/B/C", "reference_sensitivity",
  "obstruction_gli2012_ab", "GLI-2012 obstruction, strict A/B", "obstruction", "GLI-2012 race/ethnicity-specific LLN", "A/B", "quality_sensitivity",
  "low_fev1_gli_global2022_abc", "GLI Global 2022 low FEV1", "low_fev1", "GLI Global 2022 race-neutral z<-1.645", "A/B/C", "prism_component",
  "low_fvc_gli_global2022_abc", "GLI Global 2022 low FVC", "low_fvc", "GLI Global 2022 race-neutral z<-1.645", "A/B/C", "restriction_screen",
  "prism_gli_global2022_z_abc", "GLI Global 2022 PRISm-z", "prism", "Preserved FEV1/FVC z and low FEV1 z", "A/B/C", "primary_prism",
  "prism_gli_global2022_80_abc", "GLI Global 2022 PRISm-80", "prism", "FEV1/FVC >=0.70 and FEV1 <80% predicted", "A/B/C", "conventional_prism",
  "low_fev1_gli2012_abc", "GLI-2012 low FEV1", "low_fev1", "GLI-2012 race/ethnicity-specific z<-1.645", "A/B/C", "reference_sensitivity",
  "low_fvc_gli2012_abc", "GLI-2012 low FVC", "low_fvc", "GLI-2012 race/ethnicity-specific z<-1.645", "A/B/C", "reference_sensitivity",
  "prism_gli2012_z_abc", "GLI-2012 PRISm-z", "prism", "Preserved FEV1/FVC z and low FEV1 z", "A/B/C", "reference_sensitivity",
  "prism_gli2012_80_abc", "GLI-2012 PRISm-80", "prism", "FEV1/FVC >=0.70 and FEV1 <80% predicted", "A/B/C", "reference_sensitivity"
)

build_prevalence_table <- function(data) {
  bind_rows(lapply(seq_len(nrow(phenotype_specs)), function(i) {
    calc_prevalence(
      data = data,
      variable = phenotype_specs$variable[[i]],
      label = phenotype_specs$phenotype_label[[i]]
    )
  })) %>%
    left_join(phenotype_specs, by = c("variable", "phenotype_label")) %>%
    mutate(
      weighted_prevalence_percent = 100 * .data$weighted_prevalence,
      ci_low_percent = 100 * .data$ci_low,
      ci_high_percent = 100 * .data$ci_high,
      display = sprintf("%.1f%% (%.1f-%.1f)", .data$weighted_prevalence_percent, .data$ci_low_percent, .data$ci_high_percent)
    ) %>%
    select(
      variable, phenotype_label, phenotype_family, equation, quality_scope, lock_role,
      n_unweighted, events_unweighted, weighted_prevalence_percent, ci_low_percent,
      ci_high_percent, display
    )
}

cohen_kappa <- function(tab) {
  n <- sum(tab)
  if (n == 0) return(NA_real_)
  po <- sum(diag(tab)) / n
  pe <- sum(rowSums(tab) * colSums(tab)) / (n * n)
  if (isTRUE(all.equal(1, pe))) return(NA_real_)
  (po - pe) / (1 - pe)
}

agreement_pair <- function(data, var_a, label_a, var_b, label_b) {
  df <- data %>%
    filter(.data$adult_45plus == 1) %>%
    filter(!is.na(.data[[var_a]]), !is.na(.data[[var_b]])) %>%
    clean_design_data("wtmec6yr")
  if (nrow(df) == 0) {
    return(tibble(
      variable_a = var_a, label_a = label_a, variable_b = var_b, label_b = label_b,
      n = 0L, both_negative = NA_integer_, a_only = NA_integer_, b_only = NA_integer_,
      both_positive = NA_integer_, percent_agreement = NA_real_,
      weighted_percent_agreement = NA_real_, kappa = NA_real_
    ))
  }
  a <- as.integer(df[[var_a]])
  b <- as.integer(df[[var_b]])
  tab <- table(factor(a, levels = 0:1), factor(b, levels = 0:1))
  tibble(
    variable_a = var_a,
    label_a = label_a,
    variable_b = var_b,
    label_b = label_b,
    n = nrow(df),
    both_negative = as.integer(tab["0", "0"]),
    a_only = as.integer(tab["1", "0"]),
    b_only = as.integer(tab["0", "1"]),
    both_positive = as.integer(tab["1", "1"]),
    percent_agreement = 100 * mean(a == b, na.rm = TRUE),
    weighted_percent_agreement = 100 * sum(df$wtmec6yr * (a == b), na.rm = TRUE) / sum(df$wtmec6yr, na.rm = TRUE),
    kappa = cohen_kappa(tab)
  )
}

build_agreement_table <- function(data) {
  bind_rows(
    agreement_pair(data, "obstruction_fixed_ratio_abc", "Fixed ratio", "obstruction_gli_global2022_abc", "GLI Global 2022 LLN"),
    agreement_pair(data, "obstruction_fixed_ratio_abc", "Fixed ratio", "obstruction_gli2012_abc", "GLI-2012 LLN"),
    agreement_pair(data, "obstruction_gli_global2022_abc", "GLI Global 2022 LLN", "obstruction_gli2012_abc", "GLI-2012 LLN"),
    agreement_pair(data, "prism_gli_global2022_z_abc", "GLI Global 2022 PRISm-z", "prism_gli_global2022_80_abc", "GLI Global 2022 PRISm-80"),
    agreement_pair(data, "prism_gli_global2022_z_abc", "GLI Global 2022 PRISm-z", "prism_gli2012_z_abc", "GLI-2012 PRISm-z")
  )
}

build_quartile_rates <- function(data) {
  rate_data <- data %>%
    filter(.data$adult_45plus == 1, !is.na(.data$rv_quartile)) %>%
    clean_design_data("wtmec6yr")

  plot_specs <- phenotype_specs %>%
    filter(.data$variable %in% c(
      "obstruction_fixed_ratio_abc",
      "obstruction_gli_global2022_abc",
      "obstruction_gli2012_abc",
      "low_fvc_gli_global2022_abc",
      "prism_gli_global2022_z_abc",
      "prism_gli_global2022_80_abc"
    ))

  bind_rows(lapply(seq_len(nrow(plot_specs)), function(i) {
    outcome <- plot_specs$variable[[i]]
    outcome_label <- plot_specs$phenotype_label[[i]]
    bind_rows(lapply(levels(rate_data$rv_quartile), function(q) {
      df <- rate_data %>%
        filter(.data$rv_quartile == q, !is.na(.data[[outcome]])) %>%
        mutate(.indicator = as.numeric(.data[[outcome]] == 1))
      if (nrow(df) < 2) {
        return(tibble(
          variable = outcome,
          phenotype_label = outcome_label,
          rv_quartile = q,
          n_unweighted = nrow(df),
          events_unweighted = sum(df$.indicator == 1, na.rm = TRUE),
          weighted_prevalence = NA_real_,
          se = NA_real_,
          ci_low = NA_real_,
          ci_high = NA_real_
        ))
      }
      design <- make_design(df, "wtmec6yr")
      est <- survey::svymean(~.indicator, design = design, na.rm = TRUE)
      prop <- as.numeric(stats::coef(est)[1])
      se_value <- as.numeric(survey::SE(est)[1])
      ci <- pmin(pmax(prop + c(-1, 1) * stats::qnorm(0.975) * se_value, 0), 1)
      tibble(
        variable = outcome,
        phenotype_label = outcome_label,
        rv_quartile = q,
        n_unweighted = nrow(df),
        events_unweighted = sum(df$.indicator == 1, na.rm = TRUE),
        weighted_prevalence = prop,
        se = se_value,
        ci_low = ci[1],
        ci_high = ci[2]
      )
    }))
  })) %>%
    left_join(plot_specs, by = c("variable", "phenotype_label")) %>%
    mutate(
      weighted_prevalence_percent = 100 * .data$weighted_prevalence,
      ci_low_percent = 100 * .data$ci_low,
      ci_high_percent = 100 * .data$ci_high
    )
}

prepare_model_data <- function(data, outcome, exposure, covariates) {
  required <- unique(c(outcome, exposure, covariates, "wtmec6yr", "psu", "strata", "adult_45plus"))
  missing_cols <- setdiff(required, names(data))
  if (length(missing_cols) > 0) {
    stop("NHANES GLI model data missing columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  df <- data %>%
    filter(.data$adult_45plus == 1) %>%
    filter(!is.na(.data[[outcome]]), !is.na(.data[[exposure]])) %>%
    filter(!is.na(.data$wtmec6yr), .data$wtmec6yr > 0, !is.na(.data$psu), !is.na(.data$strata))

  for (v in covariates) {
    df <- df %>% filter(!is.na(.data[[v]]))
  }

  df %>%
    mutate(
      sex = factor(.data$sex),
      race_ethnicity = factor(.data$race_ethnicity),
      smoking_status = factor(.data$smoking_status, levels = c("never", "former", "current"))
    )
}

fit_model_spec <- function(data, spec) {
  covariates <- unlist(strsplit(spec$covariates, "\\|", fixed = FALSE))
  df <- prepare_model_data(data, spec$outcome, spec$exposure, covariates)
  events <- sum(df[[spec$outcome]] == 1, na.rm = TRUE)
  if (nrow(df) < 100 || events < 20) {
    stop("Insufficient records/events for model: ", spec$analysis_id, call. = FALSE)
  }

  design <- make_design(df, "wtmec6yr")
  formula <- stats::as.formula(paste(
    spec$outcome,
    "~",
    paste(c(spec$exposure, covariates), collapse = " + ")
  ))

  warnings <- character()
  model <- withCallingHandlers(
    survey::svyglm(formula, design = design, family = quasibinomial()),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  coefficients <- summary(model)$coefficients
  term <- spec$exposure
  if (!term %in% rownames(coefficients)) {
    stop("Coefficient `", term, "` not found for model: ", spec$analysis_id, call. = FALSE)
  }

  estimate <- coefficients[term, "Estimate"]
  se <- coefficients[term, "Std. Error"]
  p_value <- coefficients[term, ncol(coefficients)]
  conf_low_log <- estimate - stats::qnorm(0.975) * se
  conf_high_log <- estimate + stats::qnorm(0.975) * se
  or <- exp(estimate)
  conf_low <- exp(conf_low_log)
  conf_high <- exp(conf_high_log)

  list(
    row = tibble(
      status = "ok",
      analysis_id = spec$analysis_id,
      model_label = spec$model_label,
      outcome = spec$outcome,
      outcome_label = spec$outcome_label,
      exposure = spec$exposure,
      covariates = paste(covariates, collapse = "; "),
      n = nrow(df),
      events = events,
      estimate_log_or = estimate,
      se_log_or = se,
      or = or,
      conf_low = conf_low,
      conf_high = conf_high,
      p_value = p_value,
      p_value_formatted = format_p(p_value),
      or_ci = format_ci(or, conf_low, conf_high, 2),
      warnings_n = length(warnings),
      error = NA_character_
    ),
    warnings = tibble(
      analysis_id = spec$analysis_id,
      warning = if (length(warnings) == 0) NA_character_ else warnings
    )
  )
}

safe_fit_model_spec <- function(data, spec) {
  tryCatch(
    fit_model_spec(data, spec),
    error = function(e) {
      list(
        row = tibble(
          status = "error",
          analysis_id = spec$analysis_id,
          model_label = spec$model_label,
          outcome = spec$outcome,
          outcome_label = spec$outcome_label,
          exposure = spec$exposure,
          covariates = spec$covariates,
          n = NA_integer_,
          events = NA_integer_,
          estimate_log_or = NA_real_,
          se_log_or = NA_real_,
          or = NA_real_,
          conf_low = NA_real_,
          conf_high = NA_real_,
          p_value = NA_real_,
          p_value_formatted = NA_character_,
          or_ci = NA_character_,
          warnings_n = 1L,
          error = conditionMessage(e)
        ),
        warnings = tibble(analysis_id = spec$analysis_id, warning = conditionMessage(e))
      )
    }
  )
}

build_model_specs <- function() {
  base_covars <- "age_years|sex|race_ethnicity|bmi|smoking_status"
  tibble::tribble(
    ~analysis_id, ~model_label, ~outcome, ~outcome_label, ~exposure, ~covariates,
    "v0_5_fixed_ratio_obstruction", "Fixed ratio anchor", "obstruction_fixed_ratio_abc", "Fixed-ratio obstruction", "resp_vulnerability_z", base_covars,
    "v0_5_gli_global_obstruction", "GLI Global 2022 LLN", "obstruction_gli_global2022_abc", "GLI Global 2022 obstruction", "resp_vulnerability_z", base_covars,
    "v0_5_gli2012_obstruction", "GLI-2012 LLN", "obstruction_gli2012_abc", "GLI-2012 obstruction", "resp_vulnerability_z", base_covars,
    "v0_5_gli_global_low_fev1", "GLI Global 2022 low FEV1", "low_fev1_gli_global2022_abc", "GLI Global 2022 low FEV1", "resp_vulnerability_z", base_covars,
    "v0_5_gli_global_low_fvc", "GLI Global 2022 low FVC", "low_fvc_gli_global2022_abc", "GLI Global 2022 low FVC", "resp_vulnerability_z", base_covars,
    "v0_5_gli_global_prism_z", "GLI Global 2022 PRISm-z", "prism_gli_global2022_z_abc", "GLI Global 2022 PRISm-z", "resp_vulnerability_z", base_covars,
    "v0_5_gli_global_prism_80", "GLI Global 2022 PRISm-80", "prism_gli_global2022_80_abc", "GLI Global 2022 PRISm-80", "resp_vulnerability_z", base_covars,
    "v0_5_gli2012_prism_z", "GLI-2012 PRISm-z", "prism_gli2012_z_abc", "GLI-2012 PRISm-z", "resp_vulnerability_z", base_covars
  )
}

build_source_table <- function() {
  tibble::tribble(
    ~source_id, ~topic, ~citation_short, ~url, ~applied_rule,
    "rspiro_0_5", "R implementation", "rspiro 0.5 CRAN package", "https://cran.r-project.org/package=rspiro", "Used for GLI Global 2022, GLI-2012 predicted values, LLN, percent predicted, and z-scores.",
    "rspiro_gli_global", "GLI Global 2022 functions", "rspiro GLI Global 2022 documentation", "https://rdrr.io/cran/rspiro/man/zscore_GLIgl.html", "Race-neutral primary LLN/z-score implementation.",
    "rspiro_gli2012", "GLI-2012 functions", "rspiro GLI-2012 documentation", "https://rdrr.io/cran/rspiro/man/LLN_GLI.html", "Race/ethnicity-specific reference sensitivity implementation.",
    "gli_calculator", "Official GLI calculator", "Global Lung Function Initiative calculator", "https://gli-calculator.ersnet.org/", "Confirms GLI calculator output dimensions: predicted, z-score, LLN, and percent predicted.",
    "gli_global_2022", "Race-neutral equations", "Bowerman et al. GLI Global race-neutral equations", "https://pubmed.ncbi.nlm.nih.gov/36383197/", "Primary NHANES LLN phenotype family.",
    "gli_2012", "Multi-ethnic equations", "Quanjer et al. GLI-2012 spirometry equations", "https://publications.ersnet.org/content/erj/40/6/1324", "Secondary race/ethnicity-specific reference family.",
    "ers_ats_interpretation", "Interpretive strategy", "ERS/ATS interpretive strategies technical standard", "https://www.thoracic.org/statements/guideline-implementation-tools/technical-standards-interpretive-strategies-lung-function-tests.php", "Use LLN/z-score interpretation and keep fixed-ratio as pragmatic anchor.",
    "prism_conventional", "PRISm convention", "PRISm commonly defined as FEV1/FVC >=0.70 with FEV1 <80% predicted", "https://journal.copdfoundation.org/jcopdf/id/1371/Journal-Club-Respiratory-Impairment-With-A-Preserved-Spirometric-Ratio", "Conventional PRISm-80 sensitivity."
  )
}

build_lock_table <- function() {
  tibble::tribble(
    ~concept, ~status, ~source_variables, ~primary_rule, ~sensitivity_rule, ~notes,
    "GLI Global 2022 predicted and LLN", "locked_v0_5", "RIDAGEYR; RIAGENDR; BMXHT; SPXNFEV1; SPXNFVC", "rspiro::pred_GLIgl/LLN_GLIgl/zscore_GLIgl/pctpred_GLIgl", "None; race-neutral by design", "Primary equation layer because it avoids race/ethnicity coefficients.",
    "GLI-2012 predicted and LLN", "locked_sensitivity_v0_5", "RIDAGEYR; RIAGENDR; RIDRETH1; BMXHT; SPXNFEV1; SPXNFVC", "rspiro::pred_GLI/LLN_GLI/zscore_GLI/pctpred_GLI", "NHANES Hispanic/other categories mapped to GLI Other/mixed", "Sensitivity/reference layer only.",
    "LLN obstruction", "locked_v0_5", "SPXNFEV1; SPXNFVC; SPXNQFV1; SPXNQFVC", "FEV1/FVC z-score < -1.64485 using A/B/C quality", "Strict A/B quality sensitivity", "Fixed ratio remains a legacy anchor, not the final interpretive standard.",
    "Low FEV1", "locked_v0_5", "SPXNFEV1", "FEV1 z-score < -1.64485 using A/B/C quality", "FEV1 <80% predicted for conventional PRISm", "Used as vulnerability/PRISm component.",
    "Low FVC", "locked_v0_5", "SPXNFVC", "FVC z-score < -1.64485 using A/B/C quality", "None", "Spirometry-only screen; not a definitive restrictive diagnosis without lung volumes.",
    "PRISm", "locked_v0_5", "SPXNFEV1; SPXNFVC", "Preserved FEV1/FVC z-score >= -1.64485 and FEV1 z-score < -1.64485", "Conventional FEV1/FVC >=0.70 and FEV1 <80% predicted", "Use PRISm-z as the GLI-consistent primary phenotype; report PRISm-80 as conventional sensitivity."
  )
}

plot_quartile_rates <- function(quartile_rates, figure_dir) {
  plot_data <- quartile_rates %>%
    filter(.data$variable %in% c(
      "obstruction_fixed_ratio_abc",
      "obstruction_gli_global2022_abc",
      "obstruction_gli2012_abc",
      "low_fvc_gli_global2022_abc",
      "prism_gli_global2022_z_abc",
      "prism_gli_global2022_80_abc"
    )) %>%
    mutate(
      phenotype_label = factor(.data$phenotype_label, levels = c(
        "Fixed-ratio obstruction",
        "GLI Global 2022 obstruction",
        "GLI-2012 obstruction",
        "GLI Global 2022 low FVC",
        "GLI Global 2022 PRISm-z",
        "GLI Global 2022 PRISm-80"
      ))
    )

  p <- ggplot(plot_data, aes(x = rv_quartile, y = weighted_prevalence_percent, group = phenotype_label)) +
    geom_errorbar(aes(ymin = ci_low_percent, ymax = ci_high_percent), width = 0.10, linewidth = 0.32, colour = "#4B5563") +
    geom_line(linewidth = 0.48, colour = "#2F6F8F") +
    geom_point(size = 1.6, colour = "#2F6F8F") +
    facet_wrap(~phenotype_label, scales = "free_y", ncol = 3) +
    labs(
      x = "Respiratory vulnerability quartile",
      y = "Weighted prevalence (%)",
      title = "NHANES GLI/LLN phenotypes by respiratory vulnerability quartile"
    ) +
    theme_classic(base_size = 7.6) +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = element_line(linewidth = 0.35, colour = "black"),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      plot.title = element_text(face = "bold", size = 8.8),
      panel.grid.major.y = element_line(linewidth = 0.2, colour = "#E5E7EB")
    )

  png_path <- file.path(figure_dir, "nhanes_v0_5_gli_phenotype_rates_by_quartile.png")
  pdf_path <- file.path(figure_dir, "nhanes_v0_5_gli_phenotype_rates_by_quartile.pdf")
  ggplot2::ggsave(png_path, p, width = 7.4, height = 4.6, dpi = 600, units = "in")
  ggplot2::ggsave(pdf_path, p, width = 7.4, height = 4.6, units = "in", device = grDevices::pdf, useDingbats = FALSE)
  c(png = png_path, pdf = pdf_path)
}

plot_model_forest <- function(model_table, figure_dir) {
  plot_data <- model_table %>%
    filter(.data$status == "ok") %>%
    mutate(
      outcome_label = factor(.data$outcome_label, levels = rev(.data$outcome_label)),
      label_text = paste0(
        .data$or_ci,
        ", p",
        ifelse(stringr::str_detect(.data$p_value_formatted, "^<"), .data$p_value_formatted, paste0("=", .data$p_value_formatted))
      ),
      label_x = pmin(.data$conf_high * 1.16, 22)
    )

  p <- ggplot(plot_data, aes(y = outcome_label, x = or, xmin = conf_low, xmax = conf_high)) +
    geom_vline(xintercept = 1, linewidth = 0.35, linetype = "dashed", colour = "#777777") +
    geom_errorbarh(height = 0.16, linewidth = 0.4, colour = "#4B5563") +
    geom_point(size = 1.7, colour = "#2F6F8F") +
    geom_text(aes(x = label_x, label = label_text), hjust = 0, size = 2.15, colour = "#374151") +
    scale_x_log10(limits = c(0.75, 35), breaks = c(0.8, 1, 1.25, 1.5, 2, 3, 4, 6, 8, 12, 18, 30)) +
    labs(
      x = "Odds ratio per 1-SD higher respiratory vulnerability",
      y = NULL,
      title = "NHANES GLI/LLN phenotype models"
    ) +
    theme_classic(base_size = 7.6) +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = element_line(linewidth = 0.35, colour = "black"),
      plot.title = element_text(face = "bold", size = 8.8),
      panel.grid.major.x = element_line(linewidth = 0.2, colour = "#E5E7EB")
    )

  png_path <- file.path(figure_dir, "nhanes_v0_5_gli_model_forest.png")
  pdf_path <- file.path(figure_dir, "nhanes_v0_5_gli_model_forest.pdf")
  ggplot2::ggsave(png_path, p, width = 8.2, height = 4.6, dpi = 600, units = "in")
  ggplot2::ggsave(pdf_path, p, width = 8.2, height = 4.6, units = "in", device = grDevices::pdf, useDingbats = FALSE)
  c(png = png_path, pdf = pdf_path)
}

write_log <- function(path, prevalence_table, agreement_table, model_table, lock_table, source_table, figure_paths) {
  compact_prev <- prevalence_table %>%
    filter(.data$lock_role %in% c("legacy_anchor", "primary_lln", "primary_prism", "conventional_prism", "restriction_screen")) %>%
    select(phenotype_label, equation, n_unweighted, events_unweighted, display)

  compact_models <- model_table %>%
    select(model_label, outcome_label, n, events, or_ci, p_value_formatted, status, error)

  lines <- c(
    "# NHANES V0.5 GLI/LLN/PRISm Phenotype Lock",
    "",
    paste0("- Run date: ", Sys.Date()),
    paste0("- `rspiro` version: ", as.character(utils::packageVersion("rspiro")), "."),
    "- Dataset: local-only `derived_sensitive/nhanes/nhanes_replication_v0_2_analysis_ready.rds`.",
    "- New local row-level output: `derived_sensitive/nhanes/nhanes_replication_v0_5_gli_phenotypes.rds`.",
    "- Primary equation layer: GLI Global 2022 race-neutral.",
    "- Sensitivity/reference layer: GLI-2012 race/ethnicity-specific with NHANES Hispanic/other categories mapped to GLI Other/mixed.",
    "",
    "## Lock Table",
    "",
    markdown_table(lock_table),
    "",
    "## Key Weighted Prevalence",
    "",
    markdown_table(compact_prev),
    "",
    "## Agreement Checks",
    "",
    markdown_table(agreement_table %>% mutate(across(where(is.numeric), ~round(.x, 3)))),
    "",
    "## Respiratory Vulnerability Models",
    "",
    markdown_table(compact_models),
    "",
    "## Source Ledger",
    "",
    markdown_table(source_table),
    "",
    "## Figures",
    "",
    paste0("- Quartile rates: `", figure_paths[["rates_png"]], "` / `", figure_paths[["rates_pdf"]], "`."),
    paste0("- Model forest: `", figure_paths[["forest_png"]], "` / `", figure_paths[["forest_pdf"]], "`."),
    "",
    "## Boundary",
    "",
    "- PRISm-z is a GLI-consistent spirometry phenotype, but conventional PRISm-80 remains reported for comparability with much of the PRISm literature.",
    "- Low FVC from spirometry alone is a screen for possible restriction, not confirmed restriction without lung volumes.",
    "- NHANES remains cross-sectional; these models are replication/phenotype validation, not causal inference."
  )
  writeLines(lines, path)
}

main <- function() {
  root <- find_project_root()
  rds_path <- file.path(root, "derived_sensitive", "nhanes", "nhanes_replication_v0_2_analysis_ready.rds")
  if (!file.exists(rds_path)) {
    stop("Missing NHANES V0.2 analysis-ready RDS. Run R/nhanes/04_clean_nhanes_replication_v0_2.R first.", call. = FALSE)
  }
  if (!requireNamespace("rspiro", quietly = TRUE)) {
    stop("Package `rspiro` is required. Install with install.packages('rspiro').", call. = FALSE)
  }

  local_dir <- file.path(root, "derived_sensitive", "nhanes")
  table_dir <- file.path(root, "results", "tables")
  figure_dir <- file.path(root, "results", "figures")
  log_dir <- file.path(root, "results", "logs")
  metadata_dir <- file.path(root, "metadata")
  dir.create(local_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(metadata_dir, recursive = TRUE, showWarnings = FALSE)

  message("Reading NHANES V0.2 dataset: ", rds_path)
  data <- readRDS(rds_path)
  message("Computing GLI Global 2022 and GLI-2012 phenotypes with rspiro ", as.character(utils::packageVersion("rspiro")), ".")
  data <- add_gli_phenotypes(data)
  prepared <- prepare_analysis_data(data)
  data <- prepared$data

  local_out <- file.path(local_dir, "nhanes_replication_v0_5_gli_phenotypes.rds")
  saveRDS(data, local_out)

  prevalence_table <- build_prevalence_table(data)
  agreement_table <- build_agreement_table(data)
  quartile_rates <- build_quartile_rates(data)
  model_specs <- build_model_specs()
  fits <- lapply(seq_len(nrow(model_specs)), function(i) safe_fit_model_spec(data, model_specs[i, ]))
  model_table <- bind_rows(lapply(fits, `[[`, "row"))
  warning_table <- bind_rows(lapply(fits, `[[`, "warnings")) %>%
    filter(!is.na(.data$warning))
  if (nrow(warning_table) == 0) {
    warning_table <- tibble(analysis_id = character(), warning = character())
  }
  lock_table <- build_lock_table()
  source_table <- build_source_table()

  readr::write_csv(prepared$cutpoints, file.path(table_dir, "nhanes_v0_5_gli_vulnerability_quartile_cutpoints.csv"))
  readr::write_csv(prevalence_table, file.path(table_dir, "nhanes_v0_5_gli_phenotype_prevalence.csv"))
  readr::write_csv(agreement_table, file.path(table_dir, "nhanes_v0_5_gli_phenotype_agreement.csv"))
  readr::write_csv(quartile_rates, file.path(table_dir, "nhanes_v0_5_gli_phenotype_rates_by_quartile.csv"))
  readr::write_csv(model_table, file.path(table_dir, "nhanes_v0_5_gli_model_table.csv"))
  readr::write_csv(warning_table, file.path(table_dir, "nhanes_v0_5_gli_model_warnings.csv"))
  readr::write_csv(lock_table, file.path(metadata_dir, "nhanes_gli_lln_prism_lock_v0_5.csv"))
  readr::write_csv(source_table, file.path(metadata_dir, "nhanes_gli_reference_sources_v0_5.csv"))

  rates_paths <- plot_quartile_rates(quartile_rates, figure_dir)
  forest_paths <- plot_model_forest(model_table, figure_dir)
  figure_paths <- c(
    rates_png = rates_paths[["png"]],
    rates_pdf = rates_paths[["pdf"]],
    forest_png = forest_paths[["png"]],
    forest_pdf = forest_paths[["pdf"]]
  )

  write_log(
    path = file.path(log_dir, "nhanes_v0_5_gli_lln_prism_lock.md"),
    prevalence_table = prevalence_table,
    agreement_table = agreement_table,
    model_table = model_table,
    lock_table = lock_table,
    source_table = source_table,
    figure_paths = figure_paths
  )

  message("Wrote local GLI phenotype RDS: ", local_out)
  message("Wrote aggregate NHANES V0.5 GLI/LLN/PRISm outputs.")
}

if (sys.nframe() == 0) {
  main()
}
