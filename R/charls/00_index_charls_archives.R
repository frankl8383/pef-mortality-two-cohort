# CHARLS P0 archive inventory.
# Indexes the local CHARLS folder and archive members without copying raw data.

suppressPackageStartupMessages({
  library(dplyr)
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

normalize_rel <- function(path) {
  gsub("\\\\", "/", path)
}

is_noise_path <- function(rel_path) {
  parts <- strsplit(normalize_rel(rel_path), "/", fixed = TRUE)[[1]]
  any(parts %in% c("__MACOSX", ".DS_Store")) ||
    any(grepl("^\\._", parts)) ||
    any(grepl("^\\.DS_Store$", parts))
}

path_ext_lower <- function(path) {
  ext <- tools::file_ext(path)
  ifelse(ext == "", "", tolower(ext))
}

file_stem <- function(path) {
  tools::file_path_sans_ext(basename(path))
}

list_zip_members <- function(path) {
  out <- tryCatch(
    utils::unzip(path, list = TRUE),
    error = function(e) data.frame(Name = character(), Length = numeric(), error = conditionMessage(e))
  )
  if (!"error" %in% names(out)) {
    out$error <- NA_character_
  }
  tibble::tibble(
    inner_path = normalize_rel(out$Name),
    inner_size_bytes = suppressWarnings(as.numeric(out$Length)),
    archive_status = ifelse(is.na(out$error), "ok", out$error)
  )
}

list_rar_members <- function(path) {
  if (Sys.which("bsdtar") == "") {
    return(tibble::tibble(
      inner_path = character(),
      inner_size_bytes = numeric(),
      archive_status = "bsdtar_not_found"
    ))
  }
  res <- system2("bsdtar", c("-tf", path), stdout = TRUE, stderr = TRUE)
  status <- attr(res, "status") %||% 0
  if (!identical(status, 0)) {
    return(tibble::tibble(
      inner_path = character(),
      inner_size_bytes = numeric(),
      archive_status = paste(res, collapse = " | ")
    ))
  }
  tibble::tibble(
    inner_path = normalize_rel(res[res != ""]),
    inner_size_bytes = NA_real_,
    archive_status = "ok"
  )
}

score_preferred_archive <- function(archive_path, inner_path) {
  archive_base <- tolower(file_stem(archive_path))
  module_base <- tolower(file_stem(inner_path))
  exact_score <- ifelse(archive_base == module_base, 0, 5)
  bundle_score <- ifelse(
    grepl("charls20|charls201|dataset|data$|household_and_community", archive_base),
    10,
    0
  )
  exact_score + bundle_score + nchar(archive_base) / 1000
}

main <- function() {
  root <- find_project_root()
  raw_dir <- Sys.getenv("CHARLS_RAW_DIR", unset = "${CHARLS_RAW_ROOT}")
  raw_dir <- normalizePath(raw_dir, mustWork = TRUE)

  data_manifest_dir <- file.path(root, "data_manifest")
  derived_dir <- file.path(root, "derived")
  log_dir <- file.path(root, "results", "logs")
  dir.create(data_manifest_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(derived_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

  all_files_abs <- list.files(raw_dir, recursive = TRUE, full.names = TRUE, all.files = TRUE, no.. = TRUE)
  all_files_abs <- all_files_abs[file.info(all_files_abs)$isdir == FALSE]
  rel <- normalize_rel(sub(paste0("^", normalizePath(raw_dir, winslash = "/"), "/?"), "", normalizePath(all_files_abs, winslash = "/")))
  keep <- !vapply(rel, is_noise_path, logical(1))
  all_files_abs <- all_files_abs[keep]
  rel <- rel[keep]

  parts <- strsplit(rel, "/", fixed = TRUE)
  inventory <- tibble::tibble(
    wave = vapply(parts, function(x) x[[1]] %||% NA_character_, character(1)),
    section = vapply(parts, function(x) x[[2]] %||% NA_character_, character(1)),
    filename = basename(rel),
    size_bytes = file.info(all_files_abs)$size,
    zip_path = file.path("CHARLS", rel) |> normalize_rel(),
    local_path = normalizePath(all_files_abs, winslash = "/", mustWork = FALSE),
    extension = path_ext_lower(rel),
    archive_type = dplyr::case_when(
      extension == "zip" ~ "zip",
      extension == "rar" ~ "rar",
      TRUE ~ ""
    )
  )

  archive_inputs <- inventory %>%
    filter(archive_type %in% c("zip", "rar"))

  archive_rows <- bind_rows(lapply(seq_len(nrow(archive_inputs)), function(i) {
    row <- archive_inputs[i, ]
    members <- if (row$archive_type[[1]] == "zip") {
      list_zip_members(row$local_path[[1]])
    } else {
      list_rar_members(row$local_path[[1]])
    }
    members %>%
      mutate(
        wave = row$wave[[1]],
        section = row$section[[1]],
        outer_zip = row$filename[[1]],
        archive_path = row$zip_path[[1]],
        archive_local_path = row$local_path[[1]],
        archive_type = row$archive_type[[1]],
        archive_size_bytes = row$size_bytes[[1]]
      )
  })) %>%
    select(
      wave,
      section,
      outer_zip,
      archive_type,
      archive_path,
      archive_local_path,
      archive_size_bytes,
      inner_path,
      inner_size_bytes,
      archive_status
    ) %>%
    filter(!is.na(inner_path), inner_path != "", !grepl("/$", inner_path))

  dta_rows <- archive_rows %>%
    filter(tolower(tools::file_ext(inner_path)) == "dta") %>%
    mutate(
      module_name = file_stem(inner_path),
      preferred_score = score_preferred_archive(archive_path, inner_path)
    ) %>%
    arrange(wave, tolower(module_name), preferred_score, archive_path) %>%
    group_by(wave, module_name_lower = tolower(module_name)) %>%
    mutate(preferred = row_number() == 1) %>%
    ungroup() %>%
    select(-module_name_lower)

  write_csv(inventory %>% select(-local_path), file.path(data_manifest_dir, "charls_inventory.csv"))
  write_csv(archive_rows, file.path(data_manifest_dir, "charls_inner_zip_inventory.csv"))
  write_csv(dta_rows, file.path(derived_dir, "charls_wave_file_index.csv"))

  log <- c(
    "# CHARLS P0 Inventory Log",
    "",
    paste0("- Raw directory: `", raw_dir, "`"),
    paste0("- Inventoried non-hidden files: ", nrow(inventory)),
    paste0("- Archive files: ", sum(inventory$archive_type %in% c("zip", "rar"))),
    paste0("- Archive member rows: ", nrow(archive_rows)),
    paste0("- DTA member rows: ", nrow(dta_rows)),
    paste0("- Preferred DTA rows: ", sum(dta_rows$preferred)),
    "",
    "## Archive Status Counts",
    "",
    paste(capture.output(print(table(archive_rows$archive_status), quote = FALSE)), collapse = "\n"),
    "",
    "## Outputs",
    "",
    "- `data_manifest/charls_inventory.csv`",
    "- `data_manifest/charls_inner_zip_inventory.csv`",
    "- `derived/charls_wave_file_index.csv`"
  )
  writeLines(log, file.path(log_dir, "charls_p0_inventory_log.md"))
  message("Wrote CHARLS P0 inventory outputs under: ", root)
}

if (sys.nframe() == 0) {
  main()
}
