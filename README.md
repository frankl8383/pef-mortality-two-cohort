# Peak expiratory flow and all-cause mortality in CHARLS and NHANES

This repository accompanies the manuscript:

> Peak expiratory flow and all-cause mortality among middle-aged and older
> adults in China and the United States: parallel analyses of CHARLS and NHANES

It is a disclosure-safe reproducibility release. It contains no
participant-level records, completed imputations, direct identifiers,
credentials, local paths, or restricted CHARLS files.

## What can be reproduced publicly

The public workflow starts from six frozen aggregate inputs and rebuilds:

- Figures 1–3 in TIFF, PNG, SVG, and PDF formats;
- the source-data CSVs underlying Figures 1–3 and Tables 1–3;
- the three manuscript tables in LaTeX form; and
- a build manifest with input and output SHA-256 hashes.

The six source-data CSVs and three table files reproduce byte-for-byte. The
TIFF files are rebuilt as RGB, LZW-compressed, 600-dpi images. Font rendering
can vary across operating systems, so the public validator checks figure
dimensions and technical properties instead of requiring identical raster
bytes on every platform.

## Quick start

Python 3.11 or later is recommended.

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python scripts/validate_repository.py
```

Successful validation ends with:

```text
REPOSITORY VALIDATION PASS
```

Generated files are written to `rebuild/`. The directory is excluded from
version control because every file in it is reproducible from the tracked
aggregate inputs.

## Repository layout

```text
results/                 disclosure-safe aggregate inputs
scripts/                 display builder and repository validator
analysis_contract/       model-specification summary and prepared-data schema
validation/              frozen manuscript source-data references and log
```

## Row-level analyses and data access

The public rebuild above does not claim to recreate the cohort models from
redistributed row-level data.

- CHARLS files must be obtained by registered users from the China Health and
  Retirement Longitudinal Study under its data-use terms.
- NHANES public files and linked mortality files are available from the US
  National Center for Health Statistics.
- The model formulas, survey structures, exposure scaling, mortality-time
  handling, high-level imputation settings, and required prepared variables
  are summarized in `analysis_contract/`.

These materials support audit and independent reimplementation planning
without exposing restricted data or the authors' local authorization
machinery. They do not contain the full predictor matrix, every factor
reference level, every spline boundary and knot, or a one-command
raw-to-results pipeline.

## Software

The submitted analyses used R 4.5.1 with `survey` 4.5, `survival` 3.8-3,
`mice` 3.19.0, and `mitools` 2.4. The public display rebuild was validated
with Python 3.13.2 and the package versions in `requirements.txt`.

## License and citation

The MIT License applies to repository code. The aggregate CSV files are not
relicensed and remain subject to the terms of their underlying CHARLS and
NHANES sources; see `DATA_USE_NOTICE.md`. Please cite the associated article
and the fixed GitHub release described in `CITATION.cff`.
