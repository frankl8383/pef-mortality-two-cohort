# CHARLS Variable Dictionary Notes

This folder contains the Ticket 02 CHARLS dictionary scaffold.

## Files

- `charls_variable_dictionary_source.csv`: source table for the formatted workbook.
- `charls_crvi_domain_map.csv`: construct-level mapping from CHARLS concepts to CRVI domains and outcomes.
- `charls_variable_dictionary.xlsx`: formatted workbook generated from the CSV sources.

## Rule

No actual CHARLS variable names have been filled because no official CHARLS codebook was found in the local workspace. Every unresolved source variable remains `TODO_CODEBOOK_CHECK`.

Before cleaning or analysis:

1. Add official CHARLS codebooks locally.
2. Update `data_manifest/charls_manifest.csv`.
3. Replace each `TODO_CODEBOOK_CHECK` variable only after checking the official codebook.
4. Expand parent concepts such as ADL, IADL, frailty deficits, and biomarkers into item-level rows.
5. Run the CHARLS audit script before any cleaning script.

