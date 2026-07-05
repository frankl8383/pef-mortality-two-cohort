# NHANES v0.1 replication mapping.
# Locks codebook-verified NHANES variables for a CHARLS respiratory-vulnerability
# replication layer. Writes aggregate metadata only; does not download row-level XPTs.

suppressPackageStartupMessages({
  library(dplyr)
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

collapse_nonempty <- function(x, sep = "; ") {
  x <- x[!is.na(x) & nzchar(x)]
  if (length(x) == 0) {
    return(NA_character_)
  }
  paste(unique(x), collapse = sep)
}

split_vars <- function(x) {
  if (is.na(x) || !nzchar(x)) {
    return(character())
  }
  x <- gsub("\\s+", "", x)
  x <- unlist(strsplit(x, ";", fixed = TRUE), use.names = FALSE)
  x[nzchar(x) & !grepl("^NOT_", x)]
}

markdown_table <- function(dat) {
  dat <- as.data.frame(dat, stringsAsFactors = FALSE)
  dat[] <- lapply(dat, function(x) ifelse(is.na(x), "", as.character(x)))
  header <- paste0("| ", paste(names(dat), collapse = " | "), " |")
  separator <- paste0("| ", paste(rep("---", ncol(dat)), collapse = " | "), " |")
  rows <- apply(dat, 1, function(x) paste0("| ", paste(x, collapse = " | "), " |"))
  c(header, separator, rows)
}

fetch_url_text <- function(url) {
  if (!nzchar(Sys.which("curl"))) {
    warning("curl is not available; cannot verify URL: ", url, call. = FALSE)
    return(character())
  }
  out <- suppressWarnings(system2("curl", c("-L", "-s", "--fail", url), stdout = TRUE, stderr = TRUE))
  status <- attr(out, "status")
  if (!is.null(status) && status != 0) {
    warning("Could not fetch URL: ", url, call. = FALSE)
    return(character())
  }
  out
}

variable_present <- function(lines, variable_name) {
  if (length(lines) == 0 || is.na(variable_name) || !nzchar(variable_name)) {
    return(FALSE)
  }
  text <- paste(lines, collapse = "\n")
  grepl(paste0("\\b", variable_name, "\\b"), text, ignore.case = TRUE)
}

first_snippet <- function(lines, variable_name, context = 1) {
  if (length(lines) == 0 || is.na(variable_name) || !nzchar(variable_name)) {
    return(NA_character_)
  }
  idx <- grep(paste0("\\b", variable_name, "\\b"), lines, ignore.case = TRUE)
  if (length(idx) == 0) {
    return(NA_character_)
  }
  window <- seq.int(max(1L, idx[[1]] - context), min(length(lines), idx[[1]] + context))
  snippet <- gsub("<[^>]+>", " ", lines[window])
  snippet <- gsub("\\s+", " ", paste(trimws(snippet), collapse = " | "))
  trimws(snippet)
}

cycle_specs <- tribble(
  ~cycle, ~cycle_label, ~sddsrvyr, ~suffix, ~year_path,
  "2007-2008", "E", 5L, "E", "2007",
  "2009-2010", "F", 6L, "F", "2009",
  "2011-2012", "G", 7L, "G", "2011"
)

component_specs <- tribble(
  ~component, ~component_label, ~source_domain, ~core_role,
  "DEMO", "Demographics", "demographics", "identity, covariates, weights, design",
  "BMX", "Body Measures", "examination", "height, weight, BMI, waist",
  "SPX", "Spirometry Pre and Post-Bronchodilator", "examination", "respiratory reserve and spirometric outcomes",
  "MCQ", "Medical Conditions", "questionnaire", "self-reported asthma, emphysema, chronic bronchitis",
  "SMQ", "Smoking - Cigarette Use", "questionnaire", "smoking covariates",
  "PFQ", "Physical Functioning", "questionnaire", "functional limitation proxy",
  "DPQ", "Depression Screener", "questionnaire", "PHQ-9 depressive symptoms proxy",
  "PAQ", "Physical Activity", "questionnaire", "physical activity proxy",
  "CBC", "Complete Blood Count", "laboratory", "CBC inflammation indices",
  "GHB", "Glycohemoglobin", "laboratory", "HbA1c metabolic marker",
  "TCHOL", "Total Cholesterol", "laboratory", "total cholesterol metabolic marker",
  "HDL", "HDL Cholesterol", "laboratory", "HDL cholesterol metabolic marker",
  "BIOPRO", "Standard Biochemistry Profile", "laboratory", "albumin nutrition marker"
)

file_manifest <- merge(cycle_specs, component_specs, by = NULL) %>%
  as_tibble() %>%
  mutate(
    file_stub = paste0(.data$component, "_", .data$suffix),
    doc_url = paste0(
      "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/",
      .data$year_path,
      "/DataFiles/",
      .data$file_stub,
      ".htm"
    ),
    xpt_url = paste0(
      "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/",
      .data$year_path,
      "/DataFiles/",
      .data$file_stub,
      ".XPT"
    ),
    mapping_status = "official_codebook_url_locked_v0_1",
    row_level_allowed_in_repo = "no"
  ) %>%
  select(
    "cycle", "cycle_label", "sddsrvyr", "component", "component_label",
    "source_domain", "core_role", "file_stub", "doc_url", "xpt_url",
    "mapping_status", "row_level_allowed_in_repo"
  )

variable_targets <- tribble(
  ~concept_group, ~component, ~variable_name, ~analysis_role, ~notes,
  "identity_design", "DEMO", "SEQN", "participant_id", "Respondent sequence number.",
  "identity_design", "DEMO", "SDDSRVYR", "cycle", "Cycle number: 5, 6, 7 for 2007-2012.",
  "identity_design", "DEMO", "WTINT2YR", "interview_weight", "Two-year interview weight.",
  "identity_design", "DEMO", "WTMEC2YR", "exam_weight", "Two-year MEC exam weight.",
  "identity_design", "DEMO", "SDMVPSU", "psu", "Masked variance pseudo-PSU.",
  "identity_design", "DEMO", "SDMVSTRA", "strata", "Masked variance pseudo-stratum.",
  "demographics", "DEMO", "RIAGENDR", "sex", "Sex.",
  "demographics", "DEMO", "RIDAGEYR", "age", "Age in years; older ages may be top-coded in NHANES.",
  "demographics", "DEMO", "RIDRETH1", "race_ethnicity", "Race/ethnicity categories for 2007-2012.",
  "demographics", "DEMO", "DMDEDUC2", "education", "Adult education category.",
  "demographics", "DEMO", "INDFMPIR", "income_poverty_ratio", "Family income to poverty ratio.",
  "body_measures", "BMX", "BMXHT", "height", "Standing height in cm.",
  "body_measures", "BMX", "BMXWT", "weight", "Weight in kg.",
  "body_measures", "BMX", "BMXBMI", "bmi", "Body mass index in kg/m2.",
  "body_measures", "BMX", "BMXWAIST", "waist", "Waist circumference in cm.",
  "spirometry", "SPX", "SPXNPEF", "pef", "Baseline best-test peak expiratory flow in mL/s.",
  "spirometry", "SPX", "SPXNFEV1", "fev1", "Baseline best-test FEV1 in mL.",
  "spirometry", "SPX", "SPXNFVC", "fvc", "Baseline best-test FVC in mL.",
  "spirometry", "SPX", "SPXNQEFF", "effort_quality", "Baseline effort quality grade.",
  "spirometry", "SPX", "SPXNQFVC", "fvc_quality", "Baseline FVC quality grade.",
  "spirometry", "SPX", "SPXNQFV1", "fev1_quality", "Baseline FEV1 quality grade.",
  "spirometry", "SPX", "SPDNACC", "acceptable_curves", "Number of acceptable baseline curves.",
  "spirometry", "SPX", "SPXNSTAT", "spirometry_status", "Baseline spirometry exam status.",
  "spirometry", "SPX", "SPDBRONC", "bronchodilator_selected", "Selected for bronchodilator testing.",
  "spirometry", "SPX", "SPXBSTAT", "bronchodilator_status", "Second-test spirometry exam status.",
  "resp_questionnaire", "MCQ", "MCQ010", "ever_asthma", "Ever told by health professional that participant had asthma.",
  "resp_questionnaire", "MCQ", "MCQ035", "current_asthma", "Still has asthma.",
  "resp_questionnaire", "MCQ", "MCQ160G", "ever_emphysema", "Ever told had emphysema.",
  "resp_questionnaire", "MCQ", "MCQ160K", "ever_chronic_bronchitis", "Ever told had chronic bronchitis.",
  "resp_questionnaire", "MCQ", "MCQ170K", "current_chronic_bronchitis", "Still has chronic bronchitis.",
  "smoking", "SMQ", "SMQ020", "ever_100_cigarettes", "Smoked at least 100 cigarettes in life.",
  "smoking", "SMQ", "SMQ040", "current_smoking", "Current cigarette smoking frequency.",
  "smoking", "SMQ", "SMD030", "age_started_regular_smoking", "Age started smoking cigarettes regularly.",
  "smoking", "SMQ", "SMD055", "age_last_regular_smoking", "Age last smoked cigarettes regularly.",
  "smoking", "SMQ", "SMD057", "cigarettes_per_day_when_quit", "Cigarettes per day when quit.",
  "smoking", "SMQ", "SMD650", "cigarettes_per_day_past_30_days", "Average cigarettes per day on smoking days in past 30 days.",
  "function", "PFQ", "PFQ061A", "money_difficulty", "Managing money difficulty.",
  "function", "PFQ", "PFQ061B", "walk_quarter_mile_difficulty", "Walking quarter mile difficulty.",
  "function", "PFQ", "PFQ061C", "walk_ten_steps_difficulty", "Walking up ten steps difficulty.",
  "function", "PFQ", "PFQ061F", "house_chore_difficulty", "House chore difficulty.",
  "function", "PFQ", "PFQ061G", "prepare_meals_difficulty", "Preparing meals difficulty.",
  "function", "PFQ", "PFQ061H", "walk_rooms_difficulty", "Walking between rooms difficulty.",
  "function", "PFQ", "PFQ061I", "chair_stand_difficulty", "Standing up from armless chair difficulty.",
  "function", "PFQ", "PFQ061J", "bed_transfer_difficulty", "Getting in and out of bed difficulty.",
  "function", "PFQ", "PFQ061K", "eating_difficulty", "Using utensils or drinking from cup difficulty.",
  "function", "PFQ", "PFQ061L", "dressing_difficulty", "Dressing difficulty.",
  "depression", "DPQ", "DPQ010", "phq9_item_1", "Little interest in doing things.",
  "depression", "DPQ", "DPQ020", "phq9_item_2", "Feeling down, depressed, or hopeless.",
  "depression", "DPQ", "DPQ030", "phq9_item_3", "Trouble sleeping or sleeping too much.",
  "depression", "DPQ", "DPQ040", "phq9_item_4", "Feeling tired or having little energy.",
  "depression", "DPQ", "DPQ050", "phq9_item_5", "Poor appetite or overeating.",
  "depression", "DPQ", "DPQ060", "phq9_item_6", "Feeling bad about yourself.",
  "depression", "DPQ", "DPQ070", "phq9_item_7", "Trouble concentrating.",
  "depression", "DPQ", "DPQ080", "phq9_item_8", "Moving or speaking slowly or too fast.",
  "depression", "DPQ", "DPQ090", "phq9_item_9", "Thoughts of death or self-harm item.",
  "physical_activity", "PAQ", "PAQ605", "vigorous_work", "Vigorous work activity.",
  "physical_activity", "PAQ", "PAQ620", "moderate_work", "Moderate work activity.",
  "physical_activity", "PAQ", "PAQ635", "walk_or_bicycle", "Walk or bicycle for transportation.",
  "physical_activity", "PAQ", "PAQ650", "vigorous_recreation", "Vigorous recreational activity.",
  "physical_activity", "PAQ", "PAQ665", "moderate_recreation", "Moderate recreational activity.",
  "cbc", "CBC", "LBXWBCSI", "white_blood_cell_count", "White blood cell count.",
  "cbc", "CBC", "LBXNEPCT", "neutrophil_percent", "Segmented neutrophil percent.",
  "cbc", "CBC", "LBXLYPCT", "lymphocyte_percent", "Lymphocyte percent.",
  "cbc", "CBC", "LBXPLTSI", "platelet_count", "Platelet count.",
  "metabolic", "GHB", "LBXGH", "hba1c", "Glycohemoglobin percentage.",
  "metabolic", "TCHOL", "LBXTC", "total_cholesterol", "Total cholesterol mg/dL.",
  "metabolic", "HDL", "LBDHDD", "hdl_cholesterol", "Direct HDL cholesterol mg/dL.",
  "nutrition", "BIOPRO", "LBXSAL", "albumin", "Albumin g/dL."
)

concept_map <- tribble(
  ~mapping_id, ~charls_anchor, ~nhanes_layer, ~concept, ~analysis_role, ~nhanes_variable, ~source_component, ~cycles, ~units, ~derivation_rule, ~quality_rule, ~survey_weight, ~psu_variable, ~strata_variable, ~mapping_status, ~notes,
  "NHANES_REPL_001", "ID", "design", "participant_id", "key", "SEQN", "DEMO", "2007-2008;2009-2010;2011-2012", "not_applicable", "Direct respondent sequence number.", "Must be non-missing.", "not_applicable", "SDMVPSU", "SDMVSTRA", "locked_codebook_v0_1", "Merges files within each cycle.",
  "NHANES_REPL_002", "survey wave", "design", "cycle", "key", "SDDSRVYR", "DEMO", "2007-2008;2009-2010;2011-2012", "cycle_code", "Use SDDSRVYR values 5, 6, and 7; retain text cycle label from file suffix.", "Must be in selected cycles.", "not_applicable", "SDMVPSU", "SDMVSTRA", "locked_codebook_v0_1", "Used for pooling and sensitivity by cycle.",
  "NHANES_REPL_003", "CHARLS public PSU design caution", "design", "survey_design", "design", "WTMEC2YR;WTINT2YR;SDMVPSU;SDMVSTRA", "DEMO", "2007-2008;2009-2010;2011-2012", "weight/design", "For spirometry and MEC/lab analyses use WTMEC6YR = WTMEC2YR / 3. Use interview weights only for questionnaire-only models.", "Positive selected weight; non-missing pseudo-PSU and pseudo-stratum.", "WTMEC2YR/3 primary; WTINT2YR/3 questionnaire-only", "SDMVPSU", "SDMVSTRA", "locked_official_design_v0_1", "NHANES public files expose masked variance pseudo-PSU and pseudo-stratum.",
  "NHANES_REPL_004", "age/sex core covariates", "covariates", "demographics", "adjustment", "RIDAGEYR;RIAGENDR;RIDRETH1;DMDEDUC2;INDFMPIR", "DEMO", "2007-2008;2009-2010;2011-2012", "years/category/index", "Adult replication primary restriction age >=45; adjust age, sex, race/ethnicity, education, income poverty ratio as available.", "Set refused/don't know values to missing by item-specific codebook.", "WTMEC2YR/3 when joined to spirometry", "SDMVPSU", "SDMVSTRA", "locked_codebook_v0_1", "Race/ethnicity is a NHANES calibration covariate and subgroup descriptor.",
  "NHANES_REPL_005", "height and body size", "covariates", "body_measures", "adjustment/metabolic", "BMXHT;BMXWT;BMXBMI;BMXWAIST", "BMX", "2007-2008;2009-2010;2011-2012", "cm;kg;kg/m2;cm", "Use BMXHT in respiratory-reserve residualization; BMI and waist as body/metabolic covariates.", "Use MEC-examined records with plausible non-missing measures.", "WTMEC2YR/3", "SDMVPSU", "SDMVSTRA", "locked_codebook_v0_1", "Height is essential for CHARLS-like PEF residual modeling.",
  "NHANES_REPL_006", "r1puff peak expiratory flow", "primary_exposure_input", "baseline_pef", "exposure_component", "SPXNPEF", "SPX", "2007-2008;2009-2010;2011-2012", "mL/s; convert to L/min by multiplying by 0.06", "Primary respiratory reserve input is baseline best-test PEF. Convert to L/min for CHARLS scale comparability.", "Main rule: complete spirometry and quality grades A/B/C for FEV1/FVC context; strict sensitivity A/B only. Exclude D/F for spirometric models.", "WTMEC2YR/3", "SDMVPSU", "SDMVSTRA", "locked_codebook_v0_1", "This is the closest objective NHANES analog to CHARLS peak flow.",
  "NHANES_REPL_007", "resp_vulnerability_z", "primary_exposure", "respiratory_vulnerability_score", "derived_exposure", "SPXNPEF;RIDAGEYR;RIAGENDR;BMXHT;RIDRETH1", "SPX;DEMO;BMX", "2007-2008;2009-2010;2011-2012", "z_score", "Fit PEF_L_min ~ age + sex + height in primary CHARLS-compatible model; add race/ethnicity in NHANES calibration sensitivity. Higher vulnerability = -standardized residual.", "Require valid PEF and core covariates; run A/B/C and A/B quality sensitivity.", "WTMEC2YR/3", "SDMVPSU", "SDMVSTRA", "derived_rule_locked_v0_1", "Primary replication of CHARLS residualized PEF vulnerability.",
  "NHANES_REPL_008", "incident chronic lung disease not available cross-sectionally", "respiratory_outcome", "spirometric_obstruction", "primary_replication_outcome", "SPXNFEV1;SPXNFVC;SPXNQFV1;SPXNQFVC;SPDNACC;SPXNSTAT", "SPX", "2007-2008;2009-2010;2011-2012", "ratio", "Define FEV1/FVC = SPXNFEV1 / SPXNFVC; main obstruction threshold <0.70, with LLN sensitivity requiring a prespecified prediction equation.", "Main A/B/C quality grades for both FEV1 and FVC; strict A/B sensitivity.", "WTMEC2YR/3", "SDMVPSU", "SDMVSTRA", "derived_rule_locked_v0_1", "Cross-sectional objective outcome, not incident disease.",
  "NHANES_REPL_009", "incident chronic lung disease subtype proxy", "respiratory_outcome", "self_reported_emphysema_or_chronic_bronchitis", "secondary_outcome", "MCQ160G;MCQ160K;MCQ170K", "MCQ", "2007-2008;2009-2010;2011-2012", "binary", "Analyze emphysema and chronic bronchitis separately; composite if MCQ160G == 1 or MCQ160K == 1, with current bronchitis by MCQ170K sensitivity.", "Set refused/don't know to missing.", "WTMEC2YR/3 when joined to spirometry; WTINT2YR/3 for questionnaire-only", "SDMVPSU", "SDMVSTRA", "locked_codebook_v0_1", "Closest self-reported chronic respiratory disease proxy in NHANES.",
  "NHANES_REPL_010", "incident asthma secondary outcome", "respiratory_outcome", "self_reported_asthma", "secondary_outcome", "MCQ010;MCQ035", "MCQ", "2007-2008;2009-2010;2011-2012", "binary", "Ever asthma by MCQ010; current asthma by MCQ035 among ever-asthma respondents.", "Set refused/don't know to missing.", "WTMEC2YR/3 when joined to spirometry; WTINT2YR/3 for questionnaire-only", "SDMVPSU", "SDMVSTRA", "locked_codebook_v0_1", "Comparable to CHARLS asthma as secondary respiratory disease layer, not longitudinal incidence.",
  "NHANES_REPL_011", "smoke_ever_w1", "covariates", "smoking_status", "adjustment", "SMQ020;SMQ040;SMD030;SMD055;SMD057;SMD650", "SMQ", "2007-2008;2009-2010;2011-2012", "category;years;cigarettes/day", "Ever smoking from SMQ020; current smoking from SMQ040; intensity sensitivity from SMD650 and former-smoker intensity variables.", "Set refused/don't know to missing; avoid pack-years unless duration and intensity rules are prespecified.", "WTMEC2YR/3 when joined to spirometry", "SDMVPSU", "SDMVSTRA", "locked_codebook_v0_1", "Richer than CHARLS ever-smoking core; keep primary adjustment simple for comparability.",
  "NHANES_REPL_012", "frailty proxy: walking/ADL", "frailty_proxy", "functional_limitation", "covariate_or_outcome", "PFQ061A;PFQ061B;PFQ061C;PFQ061F;PFQ061G;PFQ061H;PFQ061I;PFQ061J;PFQ061K;PFQ061L", "PFQ", "2007-2008;2009-2010;2011-2012", "ordinal difficulty", "Build mobility difficulty from quarter-mile/ten-steps items; ADL/IADL difficulty from house chores, meals, room walking, chair stand, bed transfer, eating, dressing, money.", "Difficulty codes 2-5 can indicate limitation; code 1 no difficulty; refused/don't know missing.", "WTMEC2YR/3 when joined to spirometry; WTINT2YR/3 for questionnaire-only", "SDMVPSU", "SDMVSTRA", "partial_proxy_locked_v0_1", "NHANES does not provide a direct CHARLS-equivalent grip component for these cycles.",
  "NHANES_REPL_013", "frailty proxy: depressive symptoms", "frailty_proxy", "depressive_symptoms", "covariate", "DPQ010;DPQ020;DPQ030;DPQ040;DPQ050;DPQ060;DPQ070;DPQ080;DPQ090", "DPQ", "2007-2008;2009-2010;2011-2012", "PHQ-9 score 0-27", "Sum DPQ010-DPQ090 when valid item responses are 0-3; depressive symptom flag PHQ-9 >=10.", "Respect item nonresponse; refused/don't know to missing.", "WTMEC2YR/3 when joined to spirometry", "SDMVPSU", "SDMVSTRA", "partial_proxy_locked_v0_1", "Aligns to CHARLS CESD-10 depression component as a proxy, not identical scale.",
  "NHANES_REPL_014", "frailty proxy: low physical activity", "frailty_proxy", "physical_activity", "covariate", "PAQ605;PAQ620;PAQ635;PAQ650;PAQ665", "PAQ", "2007-2008;2009-2010;2011-2012", "binary items", "Low-activity proxy if no vigorous/moderate work, no walking/bicycling transport, and no vigorous/moderate recreation after harmonizing yes/no codes.", "Set refused/don't know to missing.", "WTMEC2YR/3 when joined to spirometry", "SDMVPSU", "SDMVSTRA", "partial_proxy_locked_v0_1", "Comparable to CHARLS low-activity component only as a proxy.",
  "NHANES_REPL_015", "inflammation/metabolic optional domain", "mechanistic_covariates", "cbc_inflammation_indices", "optional_covariate", "LBXWBCSI;LBXNEPCT;LBXLYPCT;LBXPLTSI", "CBC", "2007-2008;2009-2010;2011-2012", "1000 cells/uL; percent", "Derive neutrophil and lymphocyte counts from WBC times percent / 100; compute NLR/SII only after missingness and units checks.", "Use component eligibility and non-missing CBC values.", "WTMEC2YR/3 unless component notes require more restrictive weight", "SDMVPSU", "SDMVSTRA", "locked_codebook_v0_1_optional", "Mechanistic support layer; not needed for core replication.",
  "NHANES_REPL_016", "metabolic optional domain", "mechanistic_covariates", "metabolic_labs", "optional_covariate", "LBXGH;LBXTC;LBDHDD;LBXSAL", "GHB;TCHOL;HDL;BIOPRO", "2007-2008;2009-2010;2011-2012", "percent;mg/dL;g/dL", "Use HbA1c, total cholesterol, HDL cholesterol, and albumin as optional metabolic/nutrition covariates.", "Use non-missing lab values; avoid fasting triglycerides in core because fasting subsample weights would narrow the population.", "WTMEC2YR/3 for listed full-sample exam labs unless component notes indicate otherwise", "SDMVPSU", "SDMVSTRA", "locked_codebook_v0_1_optional", "Keep optional to protect primary model comparability and sample size."
)

dictionary_updates <- tribble(
  ~dictionary_id, ~variable_name, ~source_component, ~cycle_or_release, ~definition, ~allowed_values, ~units, ~derivation_rule, ~missingness_rule, ~analysis_role, ~codebook_status, ~qa_status, ~notes,
  "NHANES_ID_001", "SEQN", "DEMO", "2007-2012 E/F/G", "Unique respondent sequence number.", "numeric sequence", "not_applicable", "Direct from DEMO.", "Must be non-missing.", "key", "locked_codebook_v0_1", "official_html_verified", "Used to merge files within cycle.",
  "NHANES_ID_002", "SDDSRVYR", "DEMO", "2007-2012 E/F/G", "NHANES survey cycle identifier.", "5=2007-2008;6=2009-2010;7=2011-2012", "cycle_code", "Direct from DEMO; also retain file suffix labels.", "Must be one of selected cycles.", "key", "locked_codebook_v0_1", "official_html_verified", "Needed for cycle pooling.",
  "NHANES_DESIGN_001", "SDMVPSU", "DEMO", "2007-2012 E/F/G", "Masked variance pseudo-PSU.", "cycle-specific masked values", "not_applicable", "Direct from DEMO.", "Must be non-missing for survey design.", "design", "locked_official_design_v0_1", "official_html_verified", "Public masked variance unit, not true PSU.",
  "NHANES_DESIGN_002", "SDMVSTRA", "DEMO", "2007-2012 E/F/G", "Masked variance pseudo-stratum.", "cycle-specific masked values", "not_applicable", "Direct from DEMO.", "Must be non-missing for survey design.", "design", "locked_official_design_v0_1", "official_html_verified", "Public masked stratum.",
  "NHANES_DESIGN_003", "WTINT2YR", "DEMO", "2007-2012 E/F/G", "Full-sample two-year interview weight.", "positive numeric weight", "weight", "Use WTINT6YR = WTINT2YR / 3 for questionnaire-only 2007-2012 pooled analyses.", "Must be positive for included records.", "design_weight", "locked_official_design_v0_1", "official_html_verified", "Not primary when joined to spirometry.",
  "NHANES_DESIGN_004", "WTMEC2YR", "DEMO", "2007-2012 E/F/G", "Full-sample two-year MEC exam weight.", "positive numeric weight", "weight", "Use WTMEC6YR = WTMEC2YR / 3 for spirometry, body measure, and full-sample exam lab pooled analyses.", "Must be positive for included records.", "design_weight", "locked_official_design_v0_1", "official_html_verified", "Primary weight for NHANES spirometry replication.",
  "NHANES_DESIGN_005", "NOT_IN_CORE_V0_1", "lab", "deferred", "Subsample lab weight concept.", "not_applicable", "weight", "Not used in the V0.1 core mapping; avoid fasting-subsample labs in primary replication.", "not_applicable", "design_weight_optional", "deferred_not_core_v0_1", "not_applicable", "Use most restrictive subsample weight if fasting labs are added later.",
  "NHANES_DEMO_001", "RIDAGEYR", "DEMO", "2007-2012 E/F/G", "Age at screening in years.", "numeric years", "years", "Direct from DEMO; primary replication age >=45.", "Check top-coding rules; missing excluded.", "covariate", "locked_codebook_v0_1", "official_html_verified", "Core adjustment and restriction.",
  "NHANES_DEMO_002", "RIAGENDR", "DEMO", "2007-2012 E/F/G", "Gender/sex variable.", "1=male;2=female", "category", "Direct from DEMO.", "Unknown/refused to missing if present.", "covariate", "locked_codebook_v0_1", "official_html_verified", "Core adjustment.",
  "NHANES_DEMO_003", "RIDRETH1", "DEMO", "2007-2012 E/F/G", "Race/ethnicity category.", "cycle codebook categories", "category", "Direct from DEMO; harmonize labels across selected cycles.", "Missing/refused to missing.", "covariate_or_subgroup", "locked_codebook_v0_1", "official_html_verified", "Use for NHANES calibration and subgroup context.",
  "NHANES_DEMO_004", "DMDEDUC2", "DEMO", "2007-2012 E/F/G", "Education level for adults 20+.", "adult education categories", "category", "Direct from DEMO for adult models.", "Refused/don't know to missing.", "covariate", "locked_codebook_v0_1", "official_html_verified", "Adult replication uses age >=45.",
  "NHANES_DEMO_005", "INDFMPIR", "DEMO", "2007-2012 E/F/G", "Family income to poverty ratio.", "numeric ratio", "index", "Direct from DEMO.", "Document missingness.", "covariate", "locked_codebook_v0_1", "official_html_verified", "Socioeconomic covariate.",
  "NHANES_RESP_001", "SPXNFEV1", "SPX", "2007-2012 E/F/G", "Baseline best-test forced expiratory volume in 1 second.", "numeric", "mL", "Direct from SPX.", "Require spirometry quality rule for outcome models.", "exposure_or_outcome", "locked_codebook_v0_1", "official_html_verified", "Use with SPXNFVC to define obstruction.",
  "NHANES_RESP_002", "SPXNFVC", "SPX", "2007-2012 E/F/G", "Baseline best-test forced vital capacity.", "numeric", "mL", "Direct from SPX.", "Require spirometry quality rule for outcome models.", "exposure_or_outcome", "locked_codebook_v0_1", "official_html_verified", "Use with SPXNFEV1 to define obstruction.",
  "NHANES_RESP_003", "SPXNFEV1;SPXNFVC", "derived", "2007-2012 E/F/G", "Baseline FEV1/FVC ratio.", "numeric ratio", "ratio", "SPXNFEV1 / SPXNFVC after unit and quality checks.", "Missing if either component invalid.", "outcome", "derived_rule_locked_v0_1", "official_html_verified", "Primary objective cross-sectional respiratory outcome.",
  "NHANES_RESP_004", "SPXNPEF", "SPX", "2007-2012 E/F/G", "Baseline peak expiratory flow.", "numeric", "mL/s", "Direct from SPX; convert to L/min by multiplying by 0.06.", "Require valid spirometry and quality sensitivity rules.", "exposure_component", "locked_codebook_v0_1", "official_html_verified", "Closest NHANES analog to CHARLS PEF.",
  "NHANES_RESP_005", "SPXNFEV1;SPXNFVC;RIDAGEYR;RIAGENDR;BMXHT;RIDRETH1", "derived", "2007-2012 E/F/G", "PRISm phenotype candidate.", "binary", "not_applicable", "Requires preserved FEV1/FVC plus low FEV1 percent predicted; prediction equation must be prespecified.", "Missing if spirometry/covariates invalid.", "outcome_secondary", "mapped_requires_prediction_equation_v0_1", "official_html_verified", "Do not run until percent-predicted method is locked.",
  "NHANES_RESP_006", "SPXNFEV1;SPXNFVC;SPXNQFV1;SPXNQFVC", "derived", "2007-2012 E/F/G", "COPD-like spirometric obstruction phenotype.", "binary", "not_applicable", "FEV1/FVC <0.70 primary; LLN sensitivity after prediction equation lock.", "Use A/B/C main and A/B strict quality sensitivity.", "outcome", "derived_rule_locked_v0_1", "official_html_verified", "Objective cross-sectional replication outcome.",
  "NHANES_RESP_007", "MCQ010;MCQ035", "MCQ", "2007-2012 E/F/G", "Self-reported asthma history/current asthma.", "yes/no/refused/don't know", "binary", "Ever asthma MCQ010; current asthma MCQ035.", "Refused/don't know to missing.", "outcome_or_covariate", "locked_codebook_v0_1", "official_html_verified", "Secondary respiratory disease layer.",
  "NHANES_RESP_008", "MCQ160G;MCQ160K;MCQ170K", "MCQ", "2007-2012 E/F/G", "Self-reported emphysema or chronic bronchitis.", "yes/no/refused/don't know", "binary", "Separate emphysema and chronic bronchitis; optional composite.", "Refused/don't know to missing.", "outcome_or_covariate", "locked_codebook_v0_1", "official_html_verified", "Self-reported chronic respiratory disease proxy.",
  "NHANES_BODY_001", "BMXHT", "BMX", "2007-2012 E/F/G", "Standing height.", "numeric", "cm", "Direct from BMX.", "Flag implausible values; missing excluded from respiratory residualization.", "covariate", "locked_codebook_v0_1", "official_html_verified", "Needed for PEF residualization.",
  "NHANES_BODY_002", "BMXWT", "BMX", "2007-2012 E/F/G", "Body weight.", "numeric", "kg", "Direct from BMX.", "Flag implausible values.", "covariate", "locked_codebook_v0_1", "official_html_verified", "Body-size covariate.",
  "NHANES_BODY_003", "BMXBMI", "BMX", "2007-2012 E/F/G", "Body mass index.", "numeric", "kg/m2", "Use official BMI.", "Missing if invalid exam measure.", "covariate_or_component", "locked_codebook_v0_1", "official_html_verified", "Metabolic/body-size covariate.",
  "NHANES_BODY_004", "BMXWAIST", "BMX", "2007-2012 E/F/G", "Waist circumference.", "numeric", "cm", "Direct from BMX.", "Missing if invalid exam measure.", "metabolic_component", "locked_codebook_v0_1", "official_html_verified", "Optional metabolic domain.",
  "NHANES_SMOKE_001", "SMQ020;SMQ040", "SMQ", "2007-2012 E/F/G", "Current smoking status.", "yes/no/current frequency", "category", "Current smoker if SMQ020 == 1 and SMQ040 indicates every day or some days.", "Refused/don't know to missing.", "covariate", "locked_codebook_v0_1", "official_html_verified", "Core confounder.",
  "NHANES_SMOKE_002", "SMQ020;SMQ040", "SMQ", "2007-2012 E/F/G", "Former smoking status.", "yes/no/current frequency", "category", "Former smoker if SMQ020 == 1 and SMQ040 indicates not at all.", "Refused/don't know to missing.", "covariate", "locked_codebook_v0_1", "official_html_verified", "Core confounder.",
  "NHANES_SMOKE_003", "SMD030;SMD055;SMD057;SMD650", "SMQ", "2007-2012 E/F/G", "Smoking duration/intensity candidates.", "numeric age and cigarettes/day", "years;cigarettes/day", "Use intensity sensitivity only after pack-year derivation rule is prespecified.", "Refused/don't know to missing.", "covariate_or_sensitivity", "locked_codebook_v0_1_optional", "official_html_verified", "Avoid invented pack-years in primary model.",
  "NHANES_LAB_001", "LBXWBCSI", "CBC", "2007-2012 E/F/G", "White blood cell count.", "numeric", "1000 cells/uL", "Direct from CBC.", "Document missingness.", "component_optional", "locked_codebook_v0_1_optional", "official_html_verified", "Optional inflammation domain.",
  "NHANES_LAB_002", "LBXNEPCT;LBXWBCSI", "CBC", "2007-2012 E/F/G", "Neutrophil percent/count candidate.", "numeric", "percent;1000 cells/uL", "Neutrophil count = LBXWBCSI * LBXNEPCT / 100 if needed.", "Document missingness.", "component_optional", "locked_codebook_v0_1_optional", "official_html_verified", "Optional NLR/SII component.",
  "NHANES_LAB_003", "LBXLYPCT;LBXWBCSI", "CBC", "2007-2012 E/F/G", "Lymphocyte percent/count candidate.", "numeric", "percent;1000 cells/uL", "Lymphocyte count = LBXWBCSI * LBXLYPCT / 100 if needed.", "Document missingness.", "component_optional", "locked_codebook_v0_1_optional", "official_html_verified", "Optional NLR/SII component.",
  "NHANES_LAB_004", "LBXPLTSI", "CBC", "2007-2012 E/F/G", "Platelet count.", "numeric", "1000 cells/uL", "Direct from CBC.", "Document missingness.", "component_optional", "locked_codebook_v0_1_optional", "official_html_verified", "Optional SII component.",
  "NHANES_LAB_005", "NOT_IN_CORE_V0_1", "CRP_or_inflammation", "deferred", "C-reactive protein concept.", "not_applicable", "not_applicable", "Not included in V0.1 core mapping because comparable full-cycle availability across 2007-2012 is not locked here.", "not_applicable", "component_optional_deferred", "deferred_not_core_v0_1", "not_applicable", "Can be revisited as a separate lab-specific mapping.",
  "NHANES_LAB_006", "LBXGH", "GHB", "2007-2012 E/F/G", "Glycohemoglobin.", "numeric", "percent", "Direct from GHB.", "Document assay notes and missingness.", "component_optional", "locked_codebook_v0_1_optional", "official_html_verified", "Optional metabolic marker.",
  "NHANES_LAB_007", "LBXTC;LBDHDD", "TCHOL;HDL", "2007-2012 E/F/G", "Total cholesterol and HDL cholesterol.", "numeric", "mg/dL", "Direct from TCHOL and HDL.", "Document missingness and assay notes.", "component_optional", "locked_codebook_v0_1_optional", "official_html_verified", "Optional metabolic markers; fasting triglycerides excluded from core.",
  "NHANES_LAB_008", "LBXSAL", "BIOPRO", "2007-2012 E/F/G", "Albumin.", "numeric", "g/dL", "Direct from BIOPRO.", "Document missingness.", "component_optional", "locked_codebook_v0_1_optional", "official_html_verified", "Optional nutrition/inflammation marker.",
  "NHANES_FUNC_001", "PFQ061A;PFQ061B;PFQ061C;PFQ061F;PFQ061G;PFQ061H;PFQ061I;PFQ061J;PFQ061K;PFQ061L;DPQ010;DPQ020;DPQ030;DPQ040;DPQ050;DPQ060;DPQ070;DPQ080;DPQ090;PAQ605;PAQ620;PAQ635;PAQ650;PAQ665", "PFQ;DPQ;PAQ", "2007-2012 E/F/G", "Functional limitation, depression, and low-activity frailty proxy candidates.", "ordinal/binary/score", "mixed", "Build NHANES frailty proxy from mobility/ADL difficulty, PHQ-9 >=10, and low physical activity; no grip component available in this mapping.", "Item-specific refused/don't know to missing.", "outcome_or_covariate", "partial_proxy_locked_v0_1", "official_html_verified", "Proxy construct, not identical to CHARLS frailty proxy.",
  "NHANES_MORT_001", "NOT_IN_CORE_V0_1", "mortality_linkage", "deferred", "Death status from linked mortality.", "not_applicable", "binary", "Not included in V0.1 NHANES cross-sectional replication mapping.", "not_applicable", "outcome_optional", "deferred_not_core_v0_1", "not_applicable", "Can be added after public-use linked mortality file terms and variables are locked.",
  "NHANES_MORT_002", "NOT_IN_CORE_V0_1", "mortality_linkage", "deferred", "Follow-up time from linked mortality.", "not_applicable", "time", "Not included in V0.1 NHANES cross-sectional replication mapping.", "not_applicable", "outcome_time_optional", "deferred_not_core_v0_1", "not_applicable", "Can be added after mortality linkage lock."
)

survey_design_updates <- tribble(
  ~analysis_id, ~analysis_population, ~required_components, ~weight_concept, ~psu_concept, ~strata_concept, ~cycle_pooling_rule, ~eligible_outcomes, ~codebook_status, ~notes,
  "NHANES_SURVEY_001", "questionnaire_only_adults_age_45plus", "DEMO;MCQ;SMQ;PFQ;DPQ;PAQ", "WTINT2YR", "SDMVPSU", "SDMVSTRA", "WTINT6YR = WTINT2YR / 3 for pooled 2007-2012 questionnaire-only models.", "self_reported_asthma;self_reported_emphysema_or_chronic_bronchitis;functional_limitation;depressive_symptoms;physical_activity", "locked_design_v0_1", "Use only when no MEC/exam/lab/spirometry variables are included.",
  "NHANES_SURVEY_002", "exam_or_body_measure_adults_age_45plus", "DEMO;BMX;MCQ;SMQ", "WTMEC2YR", "SDMVPSU", "SDMVSTRA", "WTMEC6YR = WTMEC2YR / 3 for pooled 2007-2012 MEC models.", "bmi;waist_circumference;self_reported_respiratory_outcomes", "locked_design_v0_1", "MEC weight applies when body measures are included.",
  "NHANES_SURVEY_003", "spirometry_eligible_adults_age_45plus", "DEMO;SPX;BMX;MCQ;SMQ", "WTMEC2YR", "SDMVPSU", "SDMVSTRA", "WTMEC6YR = WTMEC2YR / 3 for pooled 2007-2012 spirometry models.", "baseline_pef;respiratory_vulnerability_score;fev1_fvc_ratio;spirometric_obstruction;self_reported_resp_outcomes", "locked_design_v0_1", "CDC SPX analytic notes specify full-sample 2-year MEC exam weight for spirometry.",
  "NHANES_SURVEY_004", "full_sample_exam_lab_adults_age_45plus_optional", "DEMO;BMX;CBC;GHB;TCHOL;HDL;BIOPRO", "WTMEC2YR", "SDMVPSU", "SDMVSTRA", "WTMEC6YR = WTMEC2YR / 3 for listed full-sample exam/lab optional markers unless component notes require otherwise.", "cbc_indices;hba1c;cholesterol;albumin", "locked_optional_design_v0_1", "Avoid fasting triglycerides in core because fasting subsample weights would narrow the sample.",
  "NHANES_SURVEY_005", "mortality_linked_adults", "DEMO;mortality_linkage;optional_baseline_components", "DEFERRED_MORTALITY_WEIGHT_LOCK", "SDMVPSU", "SDMVSTRA", "Deferred; not part of V0.1 cross-sectional NHANES replication mapping.", "death_status;followup_time", "deferred_not_core_v0_1", "Mortality extension requires a separate linked-mortality variable and eligibility lock."
)

manifest_updates <- tribble(
  ~dataset, ~cycle_or_release, ~file_role, ~expected_file_type, ~source_component, ~source_url_or_provider, ~local_path_env_var, ~access_status, ~codebook_status, ~row_level_allowed_in_repo, ~notes,
  "NHANES", "2007-2012 E/F/G", "official_codebooks", "HTML documentation", "documentation", "CDC NHANES official data documentation pages", "NHANES_CODEBOOK_DIR", "verified_online_v0_1", "official_codebook_url_locked_v0_1", "no", "Online official HTML pages verified by R/nhanes/02_build_nhanes_replication_mapping.R.",
  "NHANES", "2007-2012 E/F/G", "raw_survey_data", "XPT", "all", "CDC NHANES official XPT files", "NHANES_RAW_DIR", "not_downloaded_local", "official_xpt_url_locked_v0_1", "no", "Do not commit row-level NHANES files; download locally only when cleaning begins.",
  "NHANES", "2007-2012 E/F/G", "demographics", "XPT", "DEMO", "CDC NHANES DEMO_E/F/G", "NHANES_RAW_DIR", "not_downloaded_local", "official_codebook_verified_v0_1", "no", "Contains SEQN, demographics, weights, SDMVPSU, SDMVSTRA.",
  "NHANES", "2007-2012 E/F/G", "spirometry_or_respiratory_measurements", "XPT", "SPX", "CDC NHANES SPX_E/F/G", "NHANES_RAW_DIR", "not_downloaded_local", "official_codebook_verified_v0_1", "no", "Primary NHANES replication file for PEF, FEV1, FVC, and quality flags.",
  "NHANES", "2007-2012 E/F/G", "respiratory_questionnaire", "XPT", "MCQ", "CDC NHANES MCQ_E/F/G", "NHANES_RAW_DIR", "not_downloaded_local", "official_codebook_verified_v0_1", "no", "Asthma, emphysema, chronic bronchitis self-report items.",
  "NHANES", "2007-2012 E/F/G", "body_measures", "XPT", "BMX", "CDC NHANES BMX_E/F/G", "NHANES_RAW_DIR", "not_downloaded_local", "official_codebook_verified_v0_1", "no", "Height, weight, BMI, waist.",
  "NHANES", "2007-2012 E/F/G", "smoking_questionnaire", "XPT", "SMQ", "CDC NHANES SMQ_E/F/G", "NHANES_RAW_DIR", "not_downloaded_local", "official_codebook_verified_v0_1", "no", "Ever/current smoking and intensity candidates.",
  "NHANES", "2007-2012 E/F/G", "complete_blood_count", "XPT", "CBC", "CDC NHANES CBC_E/F/G", "NHANES_RAW_DIR", "not_downloaded_local", "official_codebook_verified_optional_v0_1", "no", "Optional CBC inflammation indices.",
  "NHANES", "deferred", "crp_or_inflammation_labs", "XPT", "CRP_or_inflammation", "CDC NHANES lab files", "NHANES_RAW_DIR", "deferred_not_core", "deferred_not_core_v0_1", "no", "CRP not locked in V0.1 core mapping; revisit separately if needed.",
  "NHANES", "2007-2012 E/F/G", "metabolic_labs", "XPT", "GHB;TCHOL;HDL;BIOPRO", "CDC NHANES GHB/TCHOL/HDL/BIOPRO_E/F/G", "NHANES_RAW_DIR", "not_downloaded_local", "official_codebook_verified_optional_v0_1", "no", "Optional HbA1c, cholesterol, HDL, albumin; fasting labs excluded from core.",
  "NHANES", "2007-2012 E/F/G", "function_limitations", "XPT", "PFQ;DPQ;PAQ", "CDC NHANES PFQ/DPQ/PAQ_E/F/G", "NHANES_RAW_DIR", "not_downloaded_local", "official_codebook_verified_v0_1", "no", "Functional limitation, depression, and physical activity proxy components.",
  "NHANES", "deferred", "mortality_linkage", "public_use_linkage_file", "mortality", "NCHS linked mortality files", "NHANES_MORTALITY_DIR", "deferred_not_core", "deferred_not_core_v0_1", "no", "Not needed for V0.1 cross-sectional replication; lock separately if survival extension is used.",
  "NHANES", "derived", "analysis_ready_output", "parquet/local-only", "local pipeline", "local pipeline", "NHANES_DERIVED_DIR", "not_created", "not_applicable", "no", "May be created only after raw XPT download and cleaning implementation; must remain local."
)

update_existing_csv <- function(path, updates, key_col) {
  if (!file.exists(path)) {
    stop("Missing file to update: ", path, call. = FALSE)
  }
  current <- readr::read_csv(path, show_col_types = FALSE, progress = FALSE)
  missing_ids <- setdiff(current[[key_col]], updates[[key_col]])
  if (length(missing_ids) > 0) {
    stop("Missing update rows for ", path, ": ", paste(missing_ids, collapse = ", "), call. = FALSE)
  }
  extra_ids <- setdiff(updates[[key_col]], current[[key_col]])
  if (length(extra_ids) > 0) {
    stop("Update contains unknown rows for ", path, ": ", paste(extra_ids, collapse = ", "), call. = FALSE)
  }
  unchanged_cols <- setdiff(names(current), names(updates))
  current %>%
    select(all_of(c(key_col, unchanged_cols))) %>%
    left_join(updates, by = key_col) %>%
    select(all_of(names(current)))
}

main <- function() {
  root <- find_project_root()
  metadata_dir <- file.path(root, "metadata")
  table_dir <- file.path(root, "results", "tables")
  log_dir <- file.path(root, "results", "logs")
  dir.create(metadata_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

  urls <- unique(file_manifest$doc_url)
  message("Fetching ", length(urls), " NHANES official codebook pages...")
  fetched <- setNames(lapply(urls, fetch_url_text), urls)

  evidence <- variable_targets %>%
    tidyr::crossing(cycle_specs %>% select("cycle", "suffix", "year_path")) %>%
    mutate(
      file_stub = paste0(.data$component, "_", .data$suffix),
      doc_url = paste0(
        "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/",
        .data$year_path,
        "/DataFiles/",
        .data$file_stub,
        ".htm"
      ),
      page_fetched = vapply(.data$doc_url, function(url) length(fetched[[url]]) > 0, logical(1)),
      variable_present = mapply(
        function(url, var) variable_present(fetched[[url]], var),
        .data$doc_url,
        .data$variable_name
      ),
      first_snippet = mapply(
        function(url, var) first_snippet(fetched[[url]], var),
        .data$doc_url,
        .data$variable_name,
        USE.NAMES = FALSE
      ),
      evidence_status = case_when(
        .data$page_fetched & .data$variable_present ~ "official_html_variable_present",
        .data$page_fetched ~ "official_html_variable_absent",
        TRUE ~ "official_html_fetch_failed"
      )
    ) %>%
    arrange(.data$component, .data$variable_name, .data$cycle)

  variable_summary <- evidence %>%
    group_by(.data$concept_group, .data$component, .data$variable_name, .data$analysis_role, .data$notes) %>%
    summarise(
      cycles_checked = paste(.data$cycle, collapse = ";"),
      cycles_present = paste(.data$cycle[.data$variable_present], collapse = ";"),
      n_cycles_present = sum(.data$variable_present),
      n_cycles_checked = n(),
      all_selected_cycles_present = all(.data$variable_present),
      evidence_status = case_when(
        all_selected_cycles_present ~ "present_all_selected_cycles",
        n_cycles_present > 0 ~ "present_partial_selected_cycles",
        TRUE ~ "not_found_selected_cycles"
      ),
      .groups = "drop"
    ) %>%
    arrange(.data$concept_group, .data$component, .data$variable_name)

  missing_required <- variable_summary %>%
    filter(!.data$all_selected_cycles_present)
  if (nrow(missing_required) > 0) {
    warning(
      "Some mapped variables were not verified in all selected cycles. See evidence table.",
      call. = FALSE
    )
  }

  concept_variable_status <- concept_map %>%
    rowwise() %>%
    mutate(
      variables_checked = paste(split_vars(.data$nhanes_variable), collapse = ";"),
      n_variables_checked = length(split_vars(.data$nhanes_variable)),
      n_variables_present_all_cycles = sum(split_vars(.data$nhanes_variable) %in% variable_summary$variable_name[variable_summary$all_selected_cycles_present]),
      evidence_status = case_when(
        n_variables_checked == 0 ~ "not_variable_based",
        n_variables_checked == n_variables_present_all_cycles ~ "all_variables_present_all_selected_cycles",
        n_variables_present_all_cycles > 0 ~ "partial_variables_present_all_selected_cycles",
        TRUE ~ "no_variables_present_all_selected_cycles"
      )
    ) %>%
    ungroup()

  mapping_coverage <- concept_variable_status %>%
    count(.data$mapping_status, .data$evidence_status, name = "n_mapping_rows") %>%
    arrange(.data$mapping_status, .data$evidence_status)

  cycle_coverage <- evidence %>%
    group_by(.data$cycle, .data$component) %>%
    summarise(
      n_variables_checked = n(),
      n_variables_present = sum(.data$variable_present),
      n_pages_fetched = sum(.data$page_fetched),
      all_variables_present = all(.data$variable_present),
      .groups = "drop"
    ) %>%
    arrange(.data$cycle, .data$component)

  dictionary_path <- file.path(root, "data_dict", "nhanes_variable_dictionary_source.csv")
  design_path <- file.path(root, "data_dict", "nhanes_survey_design_map.csv")
  manifest_path <- file.path(root, "data_manifest", "nhanes_manifest.csv")

  dictionary_out <- update_existing_csv(dictionary_path, dictionary_updates, "dictionary_id")
  design_out <- update_existing_csv(design_path, survey_design_updates, "analysis_id")

  manifest_current <- readr::read_csv(manifest_path, show_col_types = FALSE, progress = FALSE)
  if (nrow(manifest_current) != nrow(manifest_updates)) {
    stop("Manifest update row count does not match current manifest.", call. = FALSE)
  }
  manifest_out <- manifest_updates %>% select(all_of(names(manifest_current)))

  readr::write_csv(concept_variable_status, file.path(metadata_dir, "nhanes_replication_mapping_v0_1.csv"))
  readr::write_csv(variable_summary, file.path(metadata_dir, "nhanes_replication_variable_summary_v0_1.csv"))
  readr::write_csv(evidence, file.path(metadata_dir, "nhanes_replication_variable_evidence_v0_1.csv"))
  readr::write_csv(file_manifest, file.path(metadata_dir, "nhanes_replication_file_manifest_v0_1.csv"))
  readr::write_csv(mapping_coverage, file.path(table_dir, "nhanes_replication_mapping_coverage_v0_1.csv"))
  readr::write_csv(cycle_coverage, file.path(table_dir, "nhanes_replication_cycle_coverage_v0_1.csv"))
  readr::write_csv(dictionary_out, dictionary_path)
  readr::write_csv(design_out, design_path)
  readr::write_csv(manifest_out, manifest_path)

  core_rows <- concept_variable_status %>%
    filter(.data$mapping_id %in% c(
      "NHANES_REPL_006",
      "NHANES_REPL_007",
      "NHANES_REPL_008",
      "NHANES_REPL_009",
      "NHANES_REPL_010",
      "NHANES_REPL_012",
      "NHANES_REPL_013",
      "NHANES_REPL_014"
    )) %>%
    transmute(
      mapping_id,
      concept,
      analysis_role,
      nhanes_variable,
      mapping_status,
      evidence_status
    )

  log_lines <- c(
    "# NHANES Replication Mapping V0.1",
    "",
    paste0("- Run date: ", Sys.Date()),
    "- Scope: official-codebook variable mapping for NHANES 2007-2008, 2009-2010, and 2011-2012.",
    "- Primary replication layer: baseline spirometry PEF/FEV1/FVC with public NHANES survey design.",
    "- This step writes metadata and dictionary files only; it does not download or create row-level data.",
    "",
    "## Core Mapping Rows",
    "",
    markdown_table(core_rows),
    "",
    "## Survey Design Lock",
    "",
    "- Primary spirometry and MEC analyses: `WTMEC6YR = WTMEC2YR / 3`, `SDMVPSU`, `SDMVSTRA`.",
    "- Questionnaire-only sensitivity: `WTINT6YR = WTINT2YR / 3`, `SDMVPSU`, `SDMVSTRA`.",
    "- NHANES public design variables are masked variance pseudo-PSU and pseudo-strata.",
    "",
    "## Coverage",
    "",
    markdown_table(mapping_coverage),
    "",
    "## Important Boundaries",
    "",
    "- NHANES replication is cross-sectional for respiratory disease outcomes; it does not reproduce CHARLS incident disease follow-up.",
    "- PEF is mapped directly through `SPXNPEF` and converted from mL/s to L/min by multiplying by 0.06.",
    "- The primary objective outcome is spirometric obstruction from `SPXNFEV1 / SPXNFVC`; PRISm remains a secondary method task until percent-predicted equations are locked.",
    "- The frailty layer is a proxy using PFQ/DPQ/PAQ; no direct grip-strength analog is locked for these cycles.",
    "- Optional CBC/metabolic markers are mapped, but not required for the core respiratory replication.",
    "",
    "## Files Written",
    "",
    "- `metadata/nhanes_replication_mapping_v0_1.csv`",
    "- `metadata/nhanes_replication_variable_summary_v0_1.csv`",
    "- `metadata/nhanes_replication_variable_evidence_v0_1.csv`",
    "- `metadata/nhanes_replication_file_manifest_v0_1.csv`",
    "- `results/tables/nhanes_replication_mapping_coverage_v0_1.csv`",
    "- `results/tables/nhanes_replication_cycle_coverage_v0_1.csv`",
    "- `data_dict/nhanes_variable_dictionary_source.csv`",
    "- `data_dict/nhanes_survey_design_map.csv`",
    "- `data_manifest/nhanes_manifest.csv`"
  )
  log_path <- file.path(log_dir, "nhanes_replication_mapping_v0_1.md")
  writeLines(log_lines, log_path)

  todo_counts <- c(
    manifest = sum(manifest_out == "TODO_CODEBOOK_CHECK", na.rm = TRUE),
    dictionary = sum(dictionary_out == "TODO_CODEBOOK_CHECK", na.rm = TRUE),
    survey_design = sum(design_out == "TODO_CODEBOOK_CHECK", na.rm = TRUE)
  )
  if (sum(todo_counts) > 0) {
    stop("TODO placeholders remain after NHANES mapping: ", paste(names(todo_counts), todo_counts, sep = "=", collapse = ", "), call. = FALSE)
  }

  message("Wrote NHANES v0.1 replication mapping and updated scaffold dictionaries.")
  message("Log: ", log_path)
}

if (sys.nframe() == 0) {
  main()
}
