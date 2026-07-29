# PEF and mortality in CHARLS and NHANES

This repository contains the analytic code for the accompanying manuscript.
The complete versioned archive is `analysis_code_v45_8.zip`; its checksum is
recorded in `analysis_code_v45_8.sha256`.

```bash
shasum -a 256 -c analysis_code_v45_8.sha256
unzip analysis_code_v45_8.zip
```

The archive includes the cohort-construction, multiple-imputation,
survey-weighted, survival-model, and validation scripts, together with their
machine-readable specifications and aggregate validation targets. Run
instructions and the pinned R environment are inside the archive.

Participant-level data, completed imputations, fitted objects, credentials,
and author-specific paths are not included. CHARLS data must be obtained under
the [CHARLS data-use terms](https://charls.pku.edu.cn/); NHANES public-use
files are available from [NCHS](https://www.cdc.gov/nchs/nhanes/).
