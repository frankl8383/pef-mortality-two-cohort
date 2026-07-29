# PEF and mortality in CHARLS and NHANES

This repository contains the analytic code for the accompanying manuscript.
There is one final code archive, `analysis_code_final.zip`; its checksum is
recorded in `analysis_code_final.sha256`.

```bash
shasum -a 256 -c analysis_code_final.sha256
unzip analysis_code_final.zip
cd pef_mortality_analysis_final
Rscript --vanilla run_analysis.R --prepare-only
```

The archive is organized by analytic function rather than development
version: cohort construction, anthropometry, NHANES analyses, CHARLS
analyses, diagnostics, sensitivity analyses, and validation. It also includes
the final machine-readable specifications, aggregate validation targets, run
instructions, and pinned R environment.

Participant-level data, completed imputations, fitted objects, credentials,
and author-specific paths are not included. CHARLS data must be obtained under
the [CHARLS data-use terms](https://charls.pku.edu.cn/); NHANES public-use
files are available from [NCHS](https://www.cdc.gov/nchs/nhanes/).
