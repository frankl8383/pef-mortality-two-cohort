# CHARLS P1 phenotype-definition draft builder.
# Converts automated P0 variable candidates into a reviewable phenotype ledger.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
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

lower <- function(x) tolower(ifelse(is.na(x), "", x))

is_original <- function(wave) wave != "Harmonized CHARLS"
is_harmonized <- function(wave) wave == "Harmonized CHARLS"

respondent_family <- function(wave, variable_name) {
  v <- lower(variable_name)
  if (is_original(wave)) return("respondent_original")
  if (grepl("^r[0-9]", v) || grepl("^ra", v)) return("respondent_harmonized")
  if (grepl("^s[0-9]", v)) return("spouse_harmonized")
  "harmonized_other"
}

infer_wave_slot <- function(wave, variable_name) {
  v <- lower(variable_name)
  if (grepl("2011", wave)) return("w1_2011")
  if (grepl("2013", wave)) return("w2_2013")
  if (grepl("2015", wave)) return("w3_2015")
  if (grepl("2018", wave)) return("w4_2018")
  if (grepl("2020", wave)) return("w5_2020")
  m <- regmatches(v, regexpr("[rs][1-5]", v))
  if (length(m) > 0 && nchar(m) > 0) return(paste0("harmonized_", substr(m, 2, 2)))
  "not_wave_specific"
}

proposed_concept <- function(wave, module_name, variable_name, variable_label, domain) {
  v <- lower(variable_name)
  label <- lower(variable_label)
  module <- lower(module_name)

  if (v == "id") return("participant_id")
  if (v == "householdid") return("household_id")
  if (v == "communityid") return("community_id")

  if (module %in% c("demographic_background", "demographic_backgrounds")) {
    if (v == "ba004") return("age_at_interview")
    if (v %in% c("ba002_1", "ba002_2", "ba002_3")) return("birth_date_component")
    if (v == "rgender") return("sex_or_gender")
    if (v == "bd001") return("education")
    if (v %in% c("bc001", "bc002", "bc005")) return("hukou_or_residence")
  }
  if (is_harmonized(wave)) {
    if (v == "ragender") return("sex_or_gender")
    if (v == "rabyear") return("birth_year")
    if (grepl("^r[1-5]agey$", v)) return("age_at_interview")
    if (grepl("^r[1-5]iwstat$", v)) return("interview_or_exit_status")
  }

  if (v %in% c("qb002", "qb003", "qb004")) return("pef_trial_value")
  if (grepl("^qb001s", v)) return("pef_noncompletion_reason")
  if (grepl("^r[1-5]puff[0-9]?$|^r[1-5]puff(comp|pos)$", v)) return("pef_harmonized_value_or_quality")

  if (v %in% c("da007_5_", "zda007_5_", "xezdisease_2_") || grepl("^r[1-5]lung(e|f)$", v)) {
    return("chronic_lung_disease_status")
  }
  if (domain == "chronic_lung_disease") return("chronic_lung_disease_supporting_item")

  if (v %in% c("da007_14_", "zda007_14_") || grepl("^r[1-5]asthma(e|f)$", v)) {
    return("asthma_status")
  }
  if (domain == "asthma") return("asthma_supporting_item")

  if (v %in% c("qh006", "qi002", "qh006_1") || grepl("^r[1-5]mheight$", v)) return("height")
  if (v == "ql002" || grepl("^r[1-5]mweight$", v)) return("weight")
  if (grepl("^r[1-5]mbmi$", v)) return("bmi")
  if (grepl("^r[1-5](ht|wt)comp", v)) return("anthropometry_completion_flag")

  if (v %in% c("qc003", "qc004", "qc005", "qc006") || grepl("^r[1-5][lr]grip[0-9]?$", v)) {
    return("grip_strength_trial")
  }
  if (grepl("^r[1-5]grip(comp|eff|pos)$", v)) return("grip_quality_or_position")

  if (v %in% c("qg002", "qg003") || grepl("^r[1-5]wspeed[0-9]?$", v)) return("walking_speed_trial")
  if (grepl("^r[1-5]walk(1km|100|comp)", v)) return("walking_limitation_or_completion")

  if (domain == "adl_iadl") return("adl_iadl_item")
  if (domain == "depression_cesd") return("depression_cesd_item")
  if (domain == "cognition_memory") return("cognition_memory_item")

  if (v %in% c("da059", "zda059", "da061", "da061_w3", "da061_w4", "da063", "da064", "zsmoke") ||
      grepl("^r[1-5]smoke(v|n)$", v)) {
    return("smoking_status_or_intensity")
  }
  if (domain == "smoking") return("smoking_supporting_item")
  if (domain == "alcohol") return("alcohol_use")
  if (domain == "physical_activity") return("physical_activity_item")

  if (domain == "survey_design_weight") return("survey_design_or_weight")
  if (grepl("^r[1-5]iwstat$", v)) return("interview_or_exit_status")
  if (domain == "death_exit" &&
      module %in% c("exit_interview", "exit_module", "h_charls_eol_a") &&
      (grepl("^exb00[1-9]|^exb01[0-5]$", v, perl = TRUE) ||
       grepl("^xe.*(death|exit|date|died|deceased)", v, perl = TRUE) ||
       grepl("date.*death|death place|death certificate|death was|r's death|deceased|died", label, perl = TRUE))) {
    return("death_or_exit_status")
  }
  if (domain == "blood_biomarker") return("blood_biomarker")
  ""
}

analysis_role <- function(concept, domain, family) {
  if (concept %in% c("participant_id", "household_id", "community_id")) return("key")
  if (concept %in% c("age_at_interview", "birth_date_component", "birth_year", "sex_or_gender", "education", "hukou_or_residence")) return("core_covariate")
  if (concept %in% c("pef_trial_value", "pef_harmonized_value_or_quality")) return("primary_exposure_candidate")
  if (concept == "pef_noncompletion_reason") return("exposure_qc")
  if (concept %in% c("chronic_lung_disease_status", "asthma_status")) return("outcome_candidate")
  if (concept %in% c("chronic_lung_disease_supporting_item", "asthma_supporting_item")) return("outcome_supporting_item")
  if (concept %in% c("height", "weight", "bmi")) return("pef_standardization_covariate")
  if (concept %in% c("grip_strength_trial", "walking_speed_trial", "adl_iadl_item", "depression_cesd_item", "physical_activity_item")) return("frailty_component_candidate")
  if (concept %in% c("grip_quality_or_position", "walking_limitation_or_completion", "anthropometry_completion_flag")) return("frailty_or_measurement_qc")
  if (concept %in% c("smoking_status_or_intensity", "alcohol_use", "survey_design_or_weight", "death_or_exit_status", "interview_or_exit_status", "blood_biomarker")) return("covariate_or_sensitivity")
  if (domain %in% c("cognition_memory", "comorbidity_other", "healthcare_hospitalization")) return("secondary_component_candidate")
  if (family == "spouse_harmonized") return("spouse_supporting_item")
  ""
}

priority_tier <- function(concept, role, family) {
  if (role %in% c("key", "primary_exposure_candidate", "outcome_candidate", "pef_standardization_covariate") &&
      family != "spouse_harmonized") {
    return("core_review")
  }
  if (role %in% c("core_covariate", "frailty_component_candidate", "covariate_or_sensitivity", "exposure_qc") &&
      family != "spouse_harmonized") {
    return("standard_review")
  }
  if (role %in% c("outcome_supporting_item", "frailty_or_measurement_qc", "secondary_component_candidate")) {
    return("supporting_review")
  }
  "defer_or_spouse"
}

proposed_transform <- function(concept) {
  switch(
    concept,
    participant_id = "Use as stable longitudinal respondent key after duplicate and wave-link checks.",
    household_id = "Use for household linkage only when needed.",
    community_id = "Use for cluster/geographic context only under disclosure rules.",
    age_at_interview = "Use numeric age in years; verify wave-specific derivation.",
    birth_date_component = "Use only to derive age if direct age is unavailable or inconsistent.",
    birth_year = "Use only to derive age if wave age is unavailable.",
    sex_or_gender = "Harmonize to analysis categories after confirming coding.",
    pef_trial_value = "Verify units and invalid codes; derive best-of-three or prespecified summary; standardize by age, sex, and height.",
    pef_harmonized_value_or_quality = "Use as harmonized sensitivity or cross-check after confirming RAND/Harmonized definitions.",
    pef_noncompletion_reason = "Use to flag missing-not-at-random breathing-test noncompletion.",
    chronic_lung_disease_status = "Harmonize physician-diagnosed/self-reported chronic lung disease to no/yes; exclude prevalent baseline disease for incidence analyses.",
    asthma_status = "Harmonize physician-diagnosed/self-reported asthma to no/yes; decide whether to combine or separate from chronic lung disease.",
    height = "Verify cm vs m and measured vs imputed; use for PEF standardization and BMI derivation.",
    weight = "Verify kg and measured vs imputed; use for BMI/frailty covariates.",
    bmi = "Use only if definition and unit are confirmed, otherwise derive from height/weight.",
    grip_strength_trial = "Verify kg and trial quality; derive max or mean by prespecified rule.",
    walking_speed_trial = "Verify distance/time units; derive walking speed or slow-walk indicator by prespecified rule.",
    adl_iadl_item = "Recode difficulty levels into deficit indicators after confirming direction.",
    depression_cesd_item = "Score CESD-like items after confirming positive/reverse coding.",
    smoking_status_or_intensity = "Build never/former/current and intensity variables after confirming skip patterns.",
    survey_design_or_weight = "Select analysis-specific weight only after sample eligibility is fixed.",
    interview_or_exit_status = "Use to define observed interview status, attrition, death, or exit eligibility after harmonized documentation review.",
    death_or_exit_status = "Use as absorbing or competing event only after verifying timing source.",
    "Manual codebook confirmation required before use."
  )
}

main <- function() {
  root <- find_project_root()
  dict_path <- file.path(root, "metadata", "charls_variable_dictionary_draft.csv")
  if (!file.exists(dict_path)) {
    stop("Missing metadata/charls_variable_dictionary_draft.csv. Run 02_build_variable_dictionary.R first.", call. = FALSE)
  }

  metadata_dir <- file.path(root, "metadata")
  log_dir <- file.path(root, "results", "logs")
  table_dir <- file.path(root, "results", "tables")
  dir.create(metadata_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

  dictionary <- read_csv(dict_path, show_col_types = FALSE, progress = FALSE)

  draft <- dictionary %>%
    mutate(
      respondent_family = mapply(respondent_family, wave, variable_name),
      wave_slot = mapply(infer_wave_slot, wave, variable_name),
      proposed_concept = mapply(proposed_concept, wave, module_name, variable_name, variable_label, construct_domain),
      analysis_role = mapply(analysis_role, proposed_concept, construct_domain, respondent_family),
      priority_tier = mapply(priority_tier, proposed_concept, analysis_role, respondent_family),
      proposed_transform = vapply(proposed_concept, proposed_transform, character(1)),
      codebook_status = if_else(
        proposed_concept == "",
        "not_selected_for_phenotype_review",
        "P1_DRAFT_NEEDS_CODEBOOK_LOCK"
      ),
      lock_status = if_else(
        proposed_concept == "",
        "not_in_scope",
        "not_locked"
      ),
      phenotype_notes = case_when(
        proposed_concept == "" ~ "Not selected for this P1 phenotype-definition draft.",
        respondent_family == "spouse_harmonized" ~ "Spouse variable; not part of respondent primary cohort unless explicitly needed.",
        TRUE ~ "Candidate selected from Stata labels and P0 rules; verify against official CHARLS documentation before constructing data."
      )
    ) %>%
    filter(proposed_concept != "") %>%
    mutate(
      phenotype_id = sprintf("CHARLS_P1_%04d", row_number())
    ) %>%
    select(
      phenotype_id,
      priority_tier,
      analysis_role,
      proposed_concept,
      respondent_family,
      wave_slot,
      wave,
      source_file,
      archive_type,
      inner_path,
      module_name,
      variable_name,
      variable_label,
      value_labels,
      storage_class,
      construct_domain,
      candidate_harmonized_name,
      proposed_transform,
      codebook_status,
      lock_status,
      phenotype_notes
    ) %>%
    arrange(priority_tier, analysis_role, proposed_concept, wave, module_name, variable_name)

  core <- draft %>% filter(priority_tier %in% c("core_review", "standard_review"))
  counts <- draft %>%
    count(priority_tier, analysis_role, proposed_concept, name = "n") %>%
    arrange(priority_tier, analysis_role, desc(n))

  write_csv(draft, file.path(metadata_dir, "charls_phenotype_definition_draft.csv"))
  write_csv(core, file.path(metadata_dir, "charls_phenotype_definition_core_review.csv"))
  write_csv(counts, file.path(table_dir, "charls_phenotype_definition_counts.csv"))

  log <- c(
    "# CHARLS P1 Phenotype Definition Draft Log",
    "",
    paste0("- Input dictionary rows: ", nrow(dictionary)),
    paste0("- P1 phenotype draft rows: ", nrow(draft)),
    paste0("- Core/standard review rows: ", nrow(core)),
    "",
    "## Counts",
    "",
    paste(capture.output(print(counts, n = Inf)), collapse = "\n"),
    "",
    "## Boundary",
    "",
    "This file is a phenotype-definition draft, not a locked analysis specification. All selected variables keep `P1_DRAFT_NEEDS_CODEBOOK_LOCK` until official codebook confirmation."
  )
  writeLines(log, file.path(log_dir, "charls_phenotype_definition_log.md"))
  message("Wrote CHARLS P1 phenotype-definition draft.")
}

if (sys.nframe() == 0) {
  main()
}
