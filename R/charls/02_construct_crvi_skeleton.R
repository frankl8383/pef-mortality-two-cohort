# CRVI construction skeleton for analysis-ready CHARLS data.
# Inputs must already be cleaned and mapped to internal standardized names.

suppressPackageStartupMessages({
  library(dplyr)
})

required_crvi_columns <- c(
  "participant_id",
  "wave",
  "age",
  "sex",
  "pef_std_resid",
  "frailty_score"
)

optional_crvi_columns <- c(
  "inflammation_metabolic_score",
  "current_smoking",
  "education",
  "income_or_wealth",
  "urban_rural_or_region"
)

assert_analysis_ready <- function(data, required = required_crvi_columns) {
  missing_columns <- setdiff(required, names(data))
  if (length(missing_columns) > 0) {
    stop(
      "Analysis-ready CHARLS data is missing required internal columns: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

z_score <- function(x) {
  if (all(is.na(x))) {
    return(rep(NA_real_, length(x)))
  }
  scale_value <- stats::sd(x, na.rm = TRUE)
  if (is.na(scale_value) || scale_value == 0) {
    return(rep(NA_real_, length(x)))
  }
  (x - mean(x, na.rm = TRUE)) / scale_value
}

orient_component <- function(x, direction) {
  if (direction == "higher_is_more_vulnerable") {
    return(x)
  }
  if (direction == "lower_is_more_vulnerable") {
    return(-x)
  }
  stop("Unknown component direction: ", direction, call. = FALSE)
}

construct_crvi <- function(data, include_inflammation_metabolic = TRUE) {
  assert_analysis_ready(data)

  components <- data %>%
    mutate(
      crvi_low_respiratory_reserve = z_score(orient_component(pef_std_resid, "lower_is_more_vulnerable")),
      crvi_frailty_deficit = z_score(orient_component(frailty_score, "higher_is_more_vulnerable"))
    )

  if (include_inflammation_metabolic && "inflammation_metabolic_score" %in% names(components)) {
    components <- components %>%
      mutate(
        crvi_inflammation_metabolic = z_score(
          orient_component(inflammation_metabolic_score, "higher_is_more_vulnerable")
        ),
        crvi_available_domains = rowSums(!is.na(pick(
          crvi_low_respiratory_reserve,
          crvi_frailty_deficit,
          crvi_inflammation_metabolic
        ))),
        crvi = rowMeans(pick(
          crvi_low_respiratory_reserve,
          crvi_frailty_deficit,
          crvi_inflammation_metabolic
        ), na.rm = TRUE)
      )
  } else {
    components <- components %>%
      mutate(
        crvi_available_domains = rowSums(!is.na(pick(
          crvi_low_respiratory_reserve,
          crvi_frailty_deficit
        ))),
        crvi = rowMeans(pick(
          crvi_low_respiratory_reserve,
          crvi_frailty_deficit
        ), na.rm = TRUE)
      )
  }

  components %>%
    mutate(
      crvi = if_else(crvi_available_domains >= 2, crvi, NA_real_)
    )
}

