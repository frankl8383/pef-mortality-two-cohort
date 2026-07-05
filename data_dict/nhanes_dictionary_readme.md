# NHANES Dictionary Notes

Ticket 03 creates NHANES manifest and survey design scaffolding.

## Files

- `data_manifest/nhanes_manifest.csv`: file-role manifest for NHANES codebooks, raw XPT files, and optional mortality linkage.
- `data_dict/nhanes_variable_dictionary_source.csv`: concept-level dictionary. All unresolved variables use `TODO_CODEBOOK_CHECK`.
- `data_dict/nhanes_survey_design_map.csv`: survey design and weight-selection planning table.
- `R/nhanes/survey_models.R`: survey-weighted model template.
- `results/tables/model_table_nhanes.csv`: empty model-result table shell.

## Rule

Do not fill NHANES variable names from memory. Confirm variable names, cycles, eligibility, weights, PSU, strata, subsample weights, and documentation before cleaning or modeling.

## Required Before Analysis

1. Add official cycle-specific NHANES documentation locally.
2. Update `nhanes_manifest.csv`.
3. Replace `TODO_CODEBOOK_CHECK` only after codebook review.
4. Decide which cycles are eligible for spirometry and each lab component.
5. Select the correct survey weight for each analysis.
6. Run the NHANES audit script before cleaning or modeling.

