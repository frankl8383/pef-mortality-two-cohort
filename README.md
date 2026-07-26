# PEF mortality study: display reproduction

This is the minimal public display-reproduction package for the accompanying
CHARLS and NHANES manuscript. It contains one Python entry point and seven
aggregate, disclosure-safe CSV files.

## Run

Python 3.11 or newer is required by the pinned environment.

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python reproduce.py
```

A successful run prints `REPRODUCTION PASS` and creates:

- three LaTeX tables in `output/tables/`;
- four figures in TIFF, PNG, PDF, and SVG formats in `output/figures/`; and
- `output/manifest.csv`, containing SHA-256 checksums for all inputs and
  generated files.

The TIFF files are RGB, LZW-compressed, and tagged at 600 dpi.

## Scope

The script verifies input hashes, schemas, row counts, cross-file consistency,
and prespecified numerical and interpretation anchors before rendering the
manuscript displays. It does **not** reconstruct either cohort, perform
multiple imputation, or refit survival models.

No participant-level data, identifiers, completed imputations, credentials,
or local file paths are included. Access to CHARLS microdata remains governed
by [CHARLS](https://charls.pku.edu.cn/). NHANES public-use files are available
from [NCHS](https://wwwn.cdc.gov/nchs/nhanes/).

The two cohorts were modeled separately. PEF and GLI z-score effects use
different scales and should not be ranked from their hazard-ratio magnitudes.
The paired coefficient result is a ratio of hazard ratios, not an exposure
hazard ratio. The GLI reference-range analysis should not be described as a
normal- or healthy-lung subgroup, and observation-IPW analyses address only
selection related to measured covariates.
