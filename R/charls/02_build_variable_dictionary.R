# CHARLS P0 variable dictionary builder.
# Extracts Stata variable labels and value labels from preferred CHARLS DTA files.

suppressPackageStartupMessages({
  library(dplyr)
  library(haven)
  library(readr)
})

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || is.na(x)) y else x
}

find_project_root <- function() {
  script_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", script_args, value = TRUE)
  if (length(file_arg) > 0) {
    script_path <- sub("^--file=", "", file_arg[[1]])
    return(normalizePath(file.path(dirname(script_path), "..", ".."), mustWork = FALSE))
  }
  normalizePath(getwd(), mustWork = FALSE)
}

strip_charls_prefix <- function(path) {
  sub("^CHARLS/", "", gsub("\\\\", "/", path))
}

extract_member <- function(archive_local_path, archive_type, inner_path) {
  td <- tempfile("charls_dta_")
  dir.create(td, recursive = TRUE)
  out <- file.path(td, basename(inner_path))
  if (archive_type == "zip") {
    utils::unzip(archive_local_path, files = inner_path, exdir = td, overwrite = TRUE)
    extracted <- file.path(td, inner_path)
    if (!file.exists(extracted)) {
      extracted <- file.path(td, basename(inner_path))
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

label_to_string <- function(x) {
  label <- attr(x, "label")
  if (is.null(label)) "" else as.character(label)
}

value_labels_to_string <- function(x, limit = 3000) {
  labels <- attr(x, "labels")
  if (is.null(labels) || length(labels) == 0) {
    return("")
  }
  text <- paste0(names(labels), "=", unname(labels), collapse = "; ")
  if (nchar(text) > limit) {
    paste0(substr(text, 1, limit), " ...")
  } else {
    text
  }
}

classify_domain <- function(variable_name, variable_label, module_name, wave) {
  text <- tolower(paste(variable_name, variable_label, module_name, wave))
  var <- tolower(variable_name)
  label_text <- tolower(variable_label)

  if (var %in% c("id", "householdid", "communityid", "hhid", "pn", "prim_key")) {
    return("id_linkage")
  }
  if (var %in% c("qb002", "qb003", "qb004") || grepl("breathing test|peak expir|expiratory|peak flow", text)) {
    return("pef_breathing_test")
  }
  if (grepl("chronic lung|lung disease|lung diseases|pulmonary|bronchitis|emphysema|da007_5|zda007_5", text)) {
    return("chronic_lung_disease")
  }
  if (grepl("asthma|da007_14|zda007_14", text)) {
    return("asthma")
  }
  if (grepl("psu|strata|sample|longitudinal weight|respondent weight|response adjustment|blood weight|biomarker weight|weight with", text)) {
    return("survey_design_weight")
  }
  if (grepl("\\b(height|bmi)\\b|body mass|body weight|measured weight|weight\\s*\\(kg", label_text, perl = TRUE) ||
      var %in% c("qh006", "qi002", "ql002")) {
    return("anthropometry")
  }
  if (grepl("grip|hand strength|left hand|right hand", text) || var %in% c("qc003", "qc004", "qc005", "qc006")) {
    return("grip_strength")
  }
  if (grepl("walking speed|walk 1km|walk 100|walking 100|walking 1 km", text) || var %in% c("qg002", "qg003")) {
    return("walking_speed_mobility")
  }
  if (grepl("\\b(adl|iadl)\\b", text, perl = TRUE) ||
      grepl("difficulty with (dressing|bathing|eating|using the toilet|getting out of bed|doing household chores|shopping|managing your money|taking medications)", label_text, perl = TRUE) ||
      grepl("help with (dressing|bathing|eating|using toilet|shopping|managing money)", label_text, perl = TRUE) ||
      grepl("^db0(10|11|12|13|14|15|16|17|18|19|20)|^exdb0(10|11|12|13|14|15|16|17|18|19|20)", var, perl = TRUE)) {
    return("adl_iadl")
  }
  if (grepl("cesd|depress|hopeful|felt lonely|effort|sleep was restless", text)) {
    return("depression_cesd")
  }
  if (grepl("cognition|memory|recall|word recognition|mental|tics|draw", text)) {
    return("cognition_memory")
  }
  if (grepl("smok|cigarette|tobacco", text)) {
    return("smoking")
  }
  if (grepl("alcohol|drink|liquor|beer|wine", text)) {
    return("alcohol")
  }
  if (grepl("death|dead|exit|verbal autopsy|end of life|dyear|dmonth", text)) {
    return("death_exit")
  }
  if (grepl("blood|crp|hba1c|glucose|cholesterol|hdl|ldl|triglyceride|hemoglobin|cystatin|creatinine", text)) {
    return("blood_biomarker")
  }
  if (grepl("vigorous|moderate|physical activ|exercise|walking at least 10 minutes", text)) {
    return("physical_activity")
  }
  if (grepl("hypertension|diabetes|cancer|heart|stroke|kidney|liver|arthritis|dyslipidemia", text)) {
    return("comorbidity_other")
  }
  if (grepl("hospital|inpatient|outpatient|health care|insurance|doctor visit|clinic|medical care", text)) {
    return("healthcare_hospitalization")
  }
  "other"
}

candidate_name <- function(variable_name, domain) {
  var <- tolower(variable_name)
  if (domain == "id_linkage") {
    if (var == "id") return("participant_id")
    if (var == "householdid") return("household_id")
    if (var == "communityid") return("community_id")
  }
  if (domain == "pef_breathing_test") {
    if (var == "qb002") return("pef_reading_1")
    if (var == "qb003") return("pef_reading_2")
    if (var == "qb004") return("pef_reading_3")
  }
  if (domain == "chronic_lung_disease") return("chronic_lung_disease")
  if (domain == "asthma") return("asthma")
  if (domain == "anthropometry" && var %in% c("qh006", "qi002")) return("height")
  if (domain == "anthropometry" && var == "ql002") return("weight")
  if (domain == "grip_strength") return("grip_strength_component")
  if (domain == "walking_speed_mobility") return("walking_speed_or_mobility")
  if (domain == "adl_iadl") return("adl_iadl_component")
  if (domain == "depression_cesd") return("depression_cesd_component")
  if (domain == "cognition_memory") return("cognition_memory_component")
  if (domain == "smoking") return("smoking_status")
  if (domain == "alcohol") return("alcohol_use")
  if (domain == "healthcare_hospitalization") return("healthcare_hospitalization")
  if (domain == "death_exit") return("death_exit")
  if (domain == "survey_design_weight") return("survey_design_weight")
  if (domain == "blood_biomarker") return("blood_biomarker")
  if (domain == "physical_activity") return("physical_activity")
  if (domain == "comorbidity_other") return("comorbidity_other")
  ""
}

confidence_for <- function(variable_name, variable_label, domain) {
  text <- tolower(paste(variable_name, variable_label))
  var <- tolower(variable_name)
  if (domain == "other") return("low")
  if (domain == "id_linkage") return("high")
  if (domain == "pef_breathing_test" && var %in% c("qb002", "qb003", "qb004")) return("medium")
  if (domain %in% c("chronic_lung_disease", "asthma") && grepl("chronic lung|asthma|doctor", text)) return("high")
  if (domain %in% c("anthropometry", "grip_strength", "walking_speed_mobility") &&
      grepl("height|weight|left hand|right hand|walking speed", text)) return("high")
  "medium"
}

read_metadata <- function(index_row) {
  member <- extract_member(index_row$archive_local_path, index_row$archive_type, index_row$inner_path)
  on.exit(unlink(member$tmpdir, recursive = TRUE), add = TRUE)
  x <- haven::read_dta(member$path, n_max = 0)
  tibble::tibble(
    wave = index_row$wave,
    source_file = index_row$archive_path,
    archive_type = index_row$archive_type,
    inner_path = index_row$inner_path,
    module_name = index_row$module_name,
    variable_name = names(x),
    variable_label = vapply(x, label_to_string, character(1)),
    value_labels = vapply(x, value_labels_to_string, character(1)),
    storage_class = vapply(x, function(v) paste(class(v), collapse = "|"), character(1))
  )
}

main <- function() {
  root <- find_project_root()
  file_index_path <- file.path(root, "derived", "charls_wave_file_index.csv")
  if (!file.exists(file_index_path)) {
    stop("Missing file index. Run R/charls/00_index_charls_archives.R first.", call. = FALSE)
  }

  metadata_dir <- file.path(root, "metadata")
  table_dir <- file.path(root, "results", "tables")
  log_dir <- file.path(root, "results", "logs")
  dir.create(metadata_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

  index <- read_csv(file_index_path, show_col_types = FALSE, progress = FALSE) %>%
    filter(preferred)

  pieces <- vector("list", nrow(index))
  failures <- list()
  for (i in seq_len(nrow(index))) {
    row <- index[i, ]
    message("[", i, "/", nrow(index), "] metadata: ", row$wave, " :: ", row$module_name)
    pieces[[i]] <- tryCatch(
      read_metadata(row),
      error = function(e) {
        failures[[length(failures) + 1]] <<- tibble::tibble(
          wave = row$wave,
          source_file = row$archive_path,
          inner_path = row$inner_path,
          error = conditionMessage(e)
        )
        NULL
      }
    )
  }

  dictionary <- bind_rows(pieces) %>%
    mutate(
      construct_domain = mapply(classify_domain, variable_name, variable_label, module_name, wave),
      candidate_harmonized_name = mapply(candidate_name, variable_name, construct_domain),
      confidence = mapply(confidence_for, variable_name, variable_label, construct_domain),
      notes = if_else(
        construct_domain == "other",
        "Not selected by automated P0 keyword/domain rules.",
        "Automated P0 candidate. Verify against official codebook before analysis."
      )
    ) %>%
    select(
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
      notes,
      confidence
    )

  key_variable_map <- dictionary %>%
    filter(construct_domain != "other") %>%
    arrange(wave, construct_domain, source_file, inner_path, variable_name)

  key_file_map <- key_variable_map %>%
    count(wave, source_file, archive_type, inner_path, module_name, construct_domain, name = "candidate_variable_count") %>%
    arrange(wave, module_name, desc(candidate_variable_count))

  write_csv(dictionary, file.path(metadata_dir, "charls_variable_dictionary_draft.csv"))
  write_csv(key_variable_map, file.path(metadata_dir, "charls_key_variable_map.csv"))
  write_csv(key_file_map, file.path(metadata_dir, "charls_key_file_map.csv"))

  failure_path <- file.path(log_dir, "charls_variable_dictionary_failures.csv")
  if (length(failures) > 0) {
    write_csv(bind_rows(failures), failure_path)
  } else if (file.exists(failure_path)) {
    unlink(failure_path)
  }

  domain_counts <- dictionary %>%
    count(construct_domain, sort = TRUE)
  write_csv(domain_counts, file.path(table_dir, "charls_variable_domain_counts.csv"))

  log <- c(
    "# CHARLS P0 Variable Dictionary Log",
    "",
    paste0("- Preferred DTA files attempted: ", nrow(index)),
    paste0("- Variable dictionary rows: ", nrow(dictionary)),
    paste0("- Key variable candidate rows: ", nrow(key_variable_map)),
    paste0("- Failed DTA metadata reads: ", length(failures)),
    "",
    "## Domain Counts",
    "",
    paste(capture.output(print(domain_counts, n = Inf)), collapse = "\n"),
    "",
    "## Important Boundary",
    "",
    "These mappings are automated candidates from Stata labels and codebook-informed keyword rules. They are not final phenotype definitions."
  )
  writeLines(log, file.path(log_dir, "charls_variable_dictionary_log.md"))
  message("Wrote CHARLS variable dictionary outputs.")
}

if (sys.nframe() == 0) {
  main()
}
