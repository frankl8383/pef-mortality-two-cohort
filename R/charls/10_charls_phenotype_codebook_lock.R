# CHARLS v0.4 phenotype/codebook lock.
# Cross-checks current CHARLS analysis concepts against Harmonized CHARLS labels/PDF.
# Writes metadata and aggregate-only logs; does not read or write row-level data.

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
  x <- gsub("\\s+", "", x)
  x <- unlist(strsplit(x, ";", fixed = TRUE), use.names = FALSE)
  x[nzchar(x)]
}

markdown_table <- function(dat) {
  dat <- as.data.frame(dat, stringsAsFactors = FALSE)
  dat[] <- lapply(dat, function(x) ifelse(is.na(x), "", as.character(x)))
  header <- paste0("| ", paste(names(dat), collapse = " | "), " |")
  separator <- paste0("| ", paste(rep("---", ncol(dat)), collapse = " | "), " |")
  rows <- apply(dat, 1, function(x) paste0("| ", paste(x, collapse = " | "), " |"))
  c(header, separator, rows)
}

extract_pdf_text <- function(pdf_path) {
  if (!file.exists(pdf_path) || Sys.which("pdftotext") == "") {
    return(character())
  }
  out <- tempfile("charls_harmonized_pdf_", fileext = ".txt")
  on.exit(unlink(out), add = TRUE)
  status <- suppressWarnings(system2("pdftotext", c(shQuote(pdf_path), shQuote(out)), stdout = TRUE, stderr = TRUE))
  if (!file.exists(out)) {
    warning("pdftotext did not produce text output: ", paste(status, collapse = " | "), call. = FALSE)
    return(character())
  }
  readLines(out, warn = FALSE)
}

pdf_hits_for_var <- function(pdf_lines, variable_name, context = 1) {
  if (length(pdf_lines) == 0 || is.na(variable_name) || !nzchar(variable_name)) {
    return(list(n = 0L, snippet = NA_character_))
  }
  hit_idx <- grep(paste0("\\b", variable_name, "\\b"), pdf_lines, ignore.case = TRUE)
  if (length(hit_idx) == 0) {
    return(list(n = 0L, snippet = NA_character_))
  }
  first <- hit_idx[[1]]
  idx <- seq.int(max(1L, first - context), min(length(pdf_lines), first + context))
  snippet <- paste(trimws(pdf_lines[idx]), collapse = " | ")
  list(n = length(hit_idx), snippet = snippet)
}

get_dictionary_rows <- function(dictionary, vars) {
  dictionary %>%
    filter(.data$inner_path == "H_CHARLS_D_Data.dta", .data$variable_name %in% vars) %>%
    group_by(.data$variable_name) %>%
    slice(1) %>%
    ungroup()
}

concept_var_table <- tribble(
  ~concept_id, ~concept_label, ~analysis_role, ~source_variables, ~derived_variables, ~current_rule,
  "id_linkage", "Participant, household, and community identifiers", "linkage/design", "ID;householdID;communityID", "participant_id;household_id;community_id", "Use `ID` as respondent key; retain household/community IDs for linkage and public-data PSU design only.",
  "core_covariates", "Baseline demographic and anthropometric covariates", "adjustment set", "rabyear;ragender;r1agey;r1mheight;r1mbmi;r1smokev", "birth_year;sex_code;age_w1;height_m_w1;bmi_w1;smoke_ever_w1;age_decade;sex;smoke_ever", "Adjust models for baseline age, sex, smoking ever, and BMI; height is included in PEF residualization and raw PEF sensitivity.",
  "respiratory_reserve_pef", "Baseline peak expiratory flow reserve", "primary exposure input", "r1puff;r1puff1;r1puff2;r1puff3;r1puffcomp;r1puffpos;r1puffeff", "pef_best_w1_valid_provisional", "Use Harmonized wave-1 maximum peak flow (`r1puff`); retain 30-999 L/min and flag values >=900 pending final clinical range decision.",
  "respiratory_vulnerability_score", "Respiratory vulnerability score", "primary exposure", "r1puff;r1agey;ragender;r1mheight", "resp_vulnerability_z = -pef_resid_z_w1", "Fit baseline PEF residual model adjusted for age, sex code, and measured height; multiply residual z-score by -1 so higher means lower-than-expected reserve.",
  "incident_chronic_lung_disease", "Incident chronic lung disease", "primary outcome", "r1lunge;r2lunge;r3lunge;r4lunge", "disease_event;disease_event_wave;incident_lung_w2_w4", "At risk if baseline `r1lunge == 0`; incident event is first follow-up wave with `r*lunge == 1`; death before disease is cause-specific censoring and composite-event sensitivity.",
  "incident_asthma", "Incident asthma", "secondary outcome", "r1asthmae;r2asthmae;r3asthmae;r4asthmae", "disease_event;disease_event_wave;incident_asthma_w2_w4", "At risk if baseline `r1asthmae == 0`; incident event is first follow-up wave with `r*asthmae == 1`; death before disease is cause-specific censoring and composite-event sensitivity.",
  "mortality", "Mortality waves 2-4", "severity/competing outcome", "r1iwstat;r2iwstat;r3iwstat;r4iwstat", "death_wave;death_event;died_w2_w4", "First follow-up wave with interview status 5 or 6 defines death; time is approximated by wave interval.",
  "frailty_proxy", "Baseline frailty proxy", "covariate/effect modifier", "r1gripsum;r1walk1kma;r1walk100a;r1adlab_c;r1cesd10;r1vgact_c;r1mdact_c;r1ltact_c", "low_grip_w1;mobility_difficulty_w1;adl_difficulty_w1;depressive_symptoms_w1;low_activity_w1;frailty_proxy_count_w1;frailty_proxy_ge3_w1", "Five-component proxy: sex-specific lowest grip quintile, walking difficulty, ADL difficulty, CESD-10 >=10, and no vigorous/moderate/light activity.",
  "survey_weights_design", "Survey weights and public PSU design", "survey design", "communityID;r1wtrespbiob;r1wtrespbioa", "analysis_weight;psu_community_id;county_pseudo_id", "Primary public design uses `communityID` as PSU and wave-1 biomarker weight preferring `r1wtrespbiob` over `r1wtrespbioa`; inferred county/community design is sensitivity only."
)

status_rules <- tribble(
  ~concept_id, ~lock_status, ~lock_scope, ~remaining_issue, ~publication_position,
  "id_linkage", "locked_public_harmonized", "Variable identity and analysis use are locked from Harmonized D labels plus local PDF evidence.", "Do not disclose small-area identifiers beyond allowed aggregate design use.", "Can be described as Harmonized respondent/household/community identifiers.",
  "core_covariates", "locked_public_harmonized", "Variable identity and model transform are locked for V0.4.", "Smoking is ever-smoking only in current models; richer smoking intensity can be a later sensitivity.", "Can be used in main models.",
  "respiratory_reserve_pef", "partially_locked_public_harmonized", "PEF variable identity is locked; QC variables are documented.", "Final clinical range/effort filter remains a sensitivity decision because Harmonized labels do not impose a numeric upper bound.", "Use as main exposure input with explicit range and effort-sensitivity note.",
  "respiratory_vulnerability_score", "locked_analysis_rule_v0_4", "Derived score formula is reproducible and locked for current CHARLS stage.", "Residual formula may be expanded in sensitivity models if NHANES replication suggests a different normative model.", "Can be reported as residualized peak-flow vulnerability, not as a validated external clinical index.",
  "incident_chronic_lung_disease", "locked_public_harmonized", "Variable identity, value labels, at-risk rule, and wave-time event rule are locked.", "Self-reported ever-diagnosis cannot distinguish COPD/chronic bronchitis/emphysema subtypes.", "Suitable primary CHARLS outcome.",
  "incident_asthma", "locked_public_harmonized", "Variable identity, value labels, at-risk rule, and wave-time event rule are locked.", "Lower event count and competing mortality weaken asthma-only inference.", "Suitable secondary/exploratory outcome; avoid overclaiming specificity.",
  "mortality", "locked_public_harmonized", "Interview-status death codes and wave-time rule are locked.", "Exact death dates are unavailable in this Harmonized regular-wave core.", "Suitable severity/competing outcome.",
  "frailty_proxy", "partially_locked_proxy_construct", "Source variables and component coding are locked; the composite is an analysis proxy.", "This is not an official CHARLS frailty index; label all models as frailty-proxy adjusted.", "Use as adjustment/sensitivity, not as a headline validated frailty phenotype.",
  "survey_weights_design", "partially_locked_public_design", "Weights and public PSU variable are locked; inferred county/community design is sensitivity.", "Explicit official strata are not available in public files and should not be invented.", "Report as public-data survey-weighted estimates with PSU clustering and no explicit strata."
)

doc_sources <- tribble(
  ~source_id, ~source_type, ~relative_path, ~available,
  "harmonized_regular_waves_pdf", "Harmonized codebook PDF", "Harmonized CHARLS/Regular Waves/Harmonized_CHARLS_D.pdf", NA,
  "harmonized_regular_waves_do", "Harmonized construction do file", "Harmonized CHARLS/Regular Waves/H_CHARLS_D_do_file.zip", NA,
  "wave1_2011_codebook", "Original wave codebook", "2011年全国基线调查/与数据使用相关的文档/CHARLS_codebook.rar", NA,
  "wave2_2013_codebook", "Original wave codebook", "2013年全国追踪调查/与数据使用相关的文档/CHARLS_Wave2_CodeBook.pdf", NA,
  "wave3_2015_codebook", "Original wave codebook", "2015年全国追踪调查/与数据使用相关的文档/CHARLS_2015_Codebook.pdf", NA,
  "wave4_2018_codebook", "Original wave codebook", "2018年全国追踪调查/与数据使用相关的文档/CHARLS_2018_Codebook.pdf", NA,
  "wave5_2020_codebook", "Original wave codebook", "2020年全国追踪调查/与数据使用相关的文档/CHARLS_2020_Codebook.pdf", NA
)

main <- function() {
  root <- find_project_root()
  dictionary_path <- file.path(root, "metadata", "charls_variable_dictionary_draft.csv")
  if (!file.exists(dictionary_path)) {
    stop("Missing variable dictionary. Run R/charls/02_build_variable_dictionary.R first.", call. = FALSE)
  }

  metadata_dir <- file.path(root, "metadata")
  table_dir <- file.path(root, "results", "tables")
  log_dir <- file.path(root, "results", "logs")
  dir.create(metadata_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

  charls_root <- Sys.getenv("CHARLS_RAW_DIR", unset = "${CHARLS_RAW_ROOT}")
  harmonized_pdf <- file.path(charls_root, "Harmonized CHARLS", "Regular Waves", "Harmonized_CHARLS_D.pdf")
  pdf_lines <- extract_pdf_text(harmonized_pdf)
  dictionary <- readr::read_csv(dictionary_path, show_col_types = FALSE, progress = FALSE)

  concept_vars <- concept_var_table %>%
    mutate(source_variable = strsplit(.data$source_variables, ";", fixed = TRUE)) %>%
    tidyr::unnest(cols = "source_variable") %>%
    mutate(source_variable = trimws(.data$source_variable))

  all_vars <- unique(concept_vars$source_variable)
  dict_rows <- get_dictionary_rows(dictionary, all_vars)

  evidence <- concept_vars %>%
    left_join(dict_rows, by = c("source_variable" = "variable_name")) %>%
    rowwise() %>%
    mutate(
      dictionary_present = !is.na(.data$variable_label),
      pdf_hit_count = pdf_hits_for_var(pdf_lines, .data$source_variable)$n,
      pdf_first_snippet = pdf_hits_for_var(pdf_lines, .data$source_variable)$snippet,
      value_label_present = !is.na(.data$value_labels) & nzchar(.data$value_labels),
      evidence_status = case_when(
        .data$dictionary_present & .data$pdf_hit_count > 0 ~ "dictionary_and_pdf",
        .data$dictionary_present ~ "dictionary_only",
        .data$pdf_hit_count > 0 ~ "pdf_only",
        TRUE ~ "missing_evidence"
      )
    ) %>%
    ungroup() %>%
    select(
      "concept_id", "concept_label", "source_variable", "dictionary_present",
      "variable_label", "value_labels", "storage_class", "value_label_present",
      "pdf_hit_count", "pdf_first_snippet", "evidence_status"
    )

  concept_evidence <- evidence %>%
    group_by(.data$concept_id) %>%
    summarise(
      n_source_variables = n(),
      n_dictionary_present = sum(.data$dictionary_present),
      n_pdf_present = sum(.data$pdf_hit_count > 0),
      n_with_value_labels = sum(.data$value_label_present),
      missing_source_variables = collapse_nonempty(.data$source_variable[!.data$dictionary_present]),
      evidence_status = case_when(
        n_dictionary_present == n_source_variables & n_pdf_present == n_source_variables ~ "all_dictionary_and_pdf",
        n_dictionary_present == n_source_variables & n_pdf_present > 0 ~ "all_dictionary_partial_pdf",
        n_dictionary_present == n_source_variables ~ "all_dictionary_only",
        TRUE ~ "incomplete_dictionary"
      ),
      .groups = "drop"
    )

  lock_ledger <- concept_var_table %>%
    left_join(status_rules, by = "concept_id") %>%
    left_join(concept_evidence, by = "concept_id") %>%
    mutate(
      all_source_variables_present = .data$n_dictionary_present == .data$n_source_variables,
      lock_version = "charls_phenotype_codebook_lock_v0_4_2026_06_21",
      lock_date = "2026-06-21"
    ) %>%
    select(
      "lock_version", "lock_date", "concept_id", "concept_label", "analysis_role",
      "lock_status", "lock_scope", "publication_position", "current_rule",
      "source_variables", "derived_variables", "n_source_variables",
      "n_dictionary_present", "n_pdf_present", "n_with_value_labels",
      "evidence_status", "all_source_variables_present", "missing_source_variables",
      "remaining_issue"
    )

  doc_sources_out <- doc_sources %>%
    mutate(available = file.exists(file.path(charls_root, .data$relative_path)))

  lock_summary <- lock_ledger %>%
    count(.data$lock_status, name = "n_concepts") %>%
    arrange(.data$lock_status)

  evidence_summary <- lock_ledger %>%
    count(.data$evidence_status, name = "n_concepts") %>%
    arrange(.data$evidence_status)

  readr::write_csv(lock_ledger, file.path(metadata_dir, "charls_phenotype_codebook_lock_v0_4.csv"))
  readr::write_csv(evidence, file.path(metadata_dir, "charls_phenotype_codebook_variable_evidence_v0_4.csv"))
  readr::write_csv(doc_sources_out, file.path(metadata_dir, "charls_codebook_source_inventory_v0_4.csv"))
  readr::write_csv(lock_summary, file.path(table_dir, "charls_phenotype_codebook_lock_summary_v0_4.csv"))
  readr::write_csv(evidence_summary, file.path(table_dir, "charls_phenotype_codebook_evidence_summary_v0_4.csv"))

  ready_main <- lock_ledger %>%
    filter(.data$concept_id %in% c(
      "core_covariates",
      "respiratory_vulnerability_score",
      "incident_chronic_lung_disease",
      "mortality"
    )) %>%
    select("concept_label", "lock_status", "publication_position")

  caution <- lock_ledger %>%
    filter(grepl("^partially", .data$lock_status)) %>%
    select("concept_label", "lock_status", "remaining_issue")

  log <- c(
    "# CHARLS V0.4 Phenotype/Codebook Lock",
    "",
    "This lock file cross-checks the active CHARLS analysis concepts against the local Harmonized CHARLS D dictionary, Stata labels/value labels, and the local Harmonized CHARLS D PDF. It does not read row-level data.",
    "",
    "## Source Documents",
    "",
    markdown_table(doc_sources_out %>% select("source_id", "source_type", "relative_path", "available")),
    "",
    "## Lock Summary",
    "",
    markdown_table(lock_summary),
    "",
    "## Evidence Summary",
    "",
    markdown_table(evidence_summary),
    "",
    "## Main-Analysis Ready Concepts",
    "",
    markdown_table(ready_main),
    "",
    "## Concepts Requiring Explicit Caution",
    "",
    markdown_table(caution),
    "",
    "## Interpretation",
    "",
    "- Chronic lung disease, mortality, core covariates, and the current respiratory vulnerability derivation are ready for CHARLS main modeling under the V0.4 lock.",
    "- Asthma is codebook-locked but should remain secondary because survey-weighted asthma-only estimates are imprecise.",
    "- Peak flow is variable-locked, but the exact clinical range/effort-filter decision should be handled as sensitivity analysis.",
    "- Frailty is a proxy construct assembled from locked source variables; do not present it as an official CHARLS frailty index.",
    "- Survey weights and community PSU are locked for public-data analysis; exact official strata remain unavailable in public files.",
    "",
    "## Outputs",
    "",
    "- `metadata/charls_phenotype_codebook_lock_v0_4.csv`",
    "- `metadata/charls_phenotype_codebook_variable_evidence_v0_4.csv`",
    "- `metadata/charls_codebook_source_inventory_v0_4.csv`",
    "- `results/tables/charls_phenotype_codebook_lock_summary_v0_4.csv`",
    "- `results/tables/charls_phenotype_codebook_evidence_summary_v0_4.csv`"
  )
  writeLines(log, file.path(log_dir, "charls_phenotype_codebook_lock_v0_4.md"))

  message("Wrote CHARLS v0.4 phenotype/codebook lock outputs.")
}

if (sys.nframe() == 0) {
  main()
}
