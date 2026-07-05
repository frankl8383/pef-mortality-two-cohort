# Lower-than-expected peak expiratory flow and all-cause mortality — analysis code

Reproducibility code and aggregate summary data for the two-cohort study of
lower-than-expected peak expiratory flow (PEF) and mortality in CHARLS (China)
and NHANES (United States), with an exploratory critical-care translation in
MIMIC-IV and eICU.

> **Manuscript:** *Lower-than-expected peak expiratory flow and all-cause
> mortality in older adults: a two-cohort study in China and the United States.*
> Submitted to BMC Pulmonary Medicine.

## What is (and is not) in this repository

This repository contains **analysis code** and the **aggregate summary data**
(model coefficients, hazard ratios, concordance statistics, weighted
prevalences) that underlie the published figures and tables.

It does **not** contain any individual-level / row-level records. The primary
datasets are governed by their providers' data-use agreements and must be
obtained directly from source:

- **NHANES** and its linked mortality files — publicly available from the US
  National Center for Health Statistics.
- **CHARLS** — available on registration through the CHARLS national
  data-access system.
- **MIMIC-IV** and **eICU** — available from PhysioNet to credentialed users
  who have completed the required training and data-use agreements.

## Repository layout

```
R/charls/            CHARLS harmonization, survival & sensitivity models (17 scripts)
R/nhanes/            NHANES cleaning, survey-weighted mortality models,
                     GLI-2022 phenotypes, incremental-value analyses (20 scripts)
R/icu_translation/   ICU decision/recalibration & transportability (R)
R/mimic/             MIMIC-IV cohort / endpoint construction (R)
python/icu_translation/  ICU modeling, calibration & DCA (Python)
python/mimic/            MIMIC-IV schema/cohort ETL (Python)
figures/             Published figures (PNG) + plotting_code/ (matplotlib
                     scripts that render Figs 1-4 from summary_data/)
summary_data/        Aggregate statistics behind figures/tables (51 CSVs)
data_dict/           Variable dictionaries and harmonization maps
protocol/            Prespecified analysis plan
```

Total: 51 R scripts, 10 Python scripts, 51 aggregate CSV tables.

## Data paths

Scripts reference external data through placeholders rather than absolute local
paths. Set these to your local copies before running:

- `${SECURE_DATA_ROOT}` — root holding `mimiciv/` and `eicu-crd/` (PhysioNet).
- `${CHARLS_RAW_ROOT}` — root holding CHARLS wave files.
- NHANES `.XPT` files are downloaded by `R/nhanes/03_download_nhanes_xpt_v0_2.R`.

## Environment

- **R** 4.5.x with `survey`, `survival`, `tidyverse`, `broom`.
- **Python** 3.11 with `numpy`, `pandas`, `scikit-learn`, `matplotlib`.

## Reproducing the figures

The four main figures were rendered with matplotlib; each script in
`figures/plotting_code/` reads the corresponding aggregate CSV in
`summary_data/` and needs no controlled data. The ICU supplementary figure
(Fig S1) is produced by the ICU-translation pipeline.

## Citation

Please cite the associated manuscript (details to be added on publication).

## License

Code is released under the MIT License (see `LICENSE`).
