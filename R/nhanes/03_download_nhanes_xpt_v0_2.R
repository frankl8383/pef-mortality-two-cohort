# NHANES v0.2 public XPT downloader.
# Downloads official public NHANES XPT files listed in the v0.1 mapping manifest.
# Row-level XPT files are written only to a local ignored directory.

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

local_raw_dir <- function(root) {
  env_dir <- Sys.getenv("NHANES_RAW_DIR", unset = "")
  if (nzchar(env_dir)) {
    return(normalizePath(env_dir, mustWork = FALSE))
  }
  file.path(root, "derived_sensitive", "nhanes", "raw_xpt")
}

safe_file_size <- function(path) {
  if (!file.exists(path)) {
    return(NA_real_)
  }
  as.numeric(file.info(path)$size)
}

download_one <- function(url, dest) {
  dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
  if (!nzchar(Sys.which("curl"))) {
    stop("curl is required for resumable downloads.", call. = FALSE)
  }
  args <- c(
    "--fail",
    "--location",
    "--continue-at", "-",
    "--retry", "5",
    "--retry-delay", "2",
    "--speed-time", "30",
    "--speed-limit", "1024",
    "--output", dest,
    url
  )
  out <- suppressWarnings(system2("curl", args, stdout = TRUE, stderr = TRUE))
  status <- attr(out, "status")
  list(
    status = ifelse(is.null(status), 0L, as.integer(status)),
    output = paste(out, collapse = " | ")
  )
}

main <- function() {
  root <- find_project_root()
  manifest_path <- file.path(root, "metadata", "nhanes_replication_file_manifest_v0_1.csv")
  if (!file.exists(manifest_path)) {
    stop("Missing v0.1 NHANES file manifest: ", manifest_path, call. = FALSE)
  }

  raw_dir <- local_raw_dir(root)
  dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
  local_log_dir <- file.path(root, "derived_sensitive", "nhanes")
  dir.create(local_log_dir, recursive = TRUE, showWarnings = FALSE)
  table_dir <- file.path(root, "results", "tables")
  log_dir <- file.path(root, "results", "logs")
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

  manifest <- readr::read_csv(manifest_path, show_col_types = FALSE, progress = FALSE)
  components_env <- Sys.getenv("NHANES_COMPONENTS", unset = "")
  if (nzchar(components_env)) {
    requested_components <- trimws(unlist(strsplit(components_env, ",", fixed = TRUE), use.names = FALSE))
    manifest <- manifest %>% filter(.data$component %in% requested_components)
  }
  if (nrow(manifest) == 0) {
    stop("No NHANES files selected for download.", call. = FALSE)
  }

  rows <- vector("list", nrow(manifest))
  for (i in seq_len(nrow(manifest))) {
    row <- manifest[i, ]
    dest <- file.path(raw_dir, paste0(row$file_stub, ".XPT"))
    before_size <- safe_file_size(dest)
    message("[", i, "/", nrow(manifest), "] ", row$file_stub, " -> ", dest)
    result <- download_one(row$xpt_url, dest)
    after_size <- safe_file_size(dest)
    status_label <- if (result$status == 0 && !is.na(after_size) && after_size > 0) {
      if (!is.na(before_size) && before_size == after_size) "present_or_resumed_unchanged" else "downloaded_or_resumed"
    } else {
      "download_failed"
    }
    rows[[i]] <- tibble(
      cycle = row$cycle,
      cycle_label = row$cycle_label,
      component = row$component,
      file_stub = row$file_stub,
      xpt_url = row$xpt_url,
      local_file = dest,
      size_bytes_before = before_size,
      size_bytes_after = after_size,
      download_status = status_label,
      curl_status = result$status,
      curl_output_tail = substr(result$output, max(1, nchar(result$output) - 500), nchar(result$output))
    )
  }

  download_log <- bind_rows(rows)
  local_manifest_path <- file.path(local_log_dir, "nhanes_xpt_download_manifest_v0_2.csv")
  readr::write_csv(download_log, local_manifest_path)

  summary <- download_log %>%
    group_by(.data$cycle, .data$component, .data$download_status) %>%
    summarise(
      n_files = n(),
      total_size_mb = round(sum(.data$size_bytes_after, na.rm = TRUE) / 1024^2, 3),
      .groups = "drop"
    ) %>%
    arrange(.data$cycle, .data$component, .data$download_status)
  readr::write_csv(summary, file.path(table_dir, "nhanes_xpt_download_summary_v0_2.csv"))

  failed <- download_log %>% filter(.data$download_status == "download_failed")
  log_lines <- c(
    "# NHANES XPT Download V0.2",
    "",
    paste0("- Run date: ", Sys.Date()),
    paste0("- Files selected: ", nrow(download_log)),
    paste0("- Files downloaded/present: ", sum(download_log$download_status != "download_failed")),
    paste0("- Files failed: ", nrow(failed)),
    "- Full local path manifest is stored under `derived_sensitive/nhanes/` and should remain local-only.",
    "",
    "## Components",
    "",
    paste(sort(unique(download_log$component)), collapse = ", "),
    "",
    "## Failure Summary",
    "",
    if (nrow(failed) == 0) {
      "No failed downloads."
    } else {
      paste(failed$file_stub, failed$curl_status, sep = ": ", collapse = "\n")
    }
  )
  log_path <- file.path(log_dir, "nhanes_xpt_download_v0_2.md")
  writeLines(log_lines, log_path)

  if (nrow(failed) > 0) {
    stop("Some NHANES XPT downloads failed. See: ", log_path, call. = FALSE)
  }
  message("Downloaded or verified ", nrow(download_log), " NHANES XPT files.")
  message("Local manifest: ", local_manifest_path)
  message("Summary: ", file.path(table_dir, "nhanes_xpt_download_summary_v0_2.csv"))
}

if (sys.nframe() == 0) {
  main()
}
