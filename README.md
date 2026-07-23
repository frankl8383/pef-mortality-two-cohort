# PEF mortality study

Minimal display-reproduction files for the accompanying CHARLS and NHANES
manuscript. The repository contains one Python script and six aggregate CSVs.

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python reproduce.py
```

Successful execution prints `REPRODUCTION PASS` and writes Tables 1–3,
Figures 1–3, and a hash manifest to `output/`.

This release contains no participant-level data or completed imputations and
does not reproduce cohort construction, imputation, or model fitting. CHARLS
data require authorization from [CHARLS](https://charls.pku.edu.cn/); NHANES
public files are available from [NCHS](https://wwwn.cdc.gov/nchs/nhanes/).
