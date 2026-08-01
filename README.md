# PEF and mortality in CHARLS and NHANES

This repository contains one final analytic-code archive for the accompanying
manuscript. It includes the V46.0 CHARLS adjacent-window primary analysis and
the NHANES strict A/B spirometry-quality sensitivity.

```bash
shasum -a 256 -c analysis_code_final.sha256
unzip analysis_code_final.zip
cd pef_mortality_analysis_final
Rscript --vanilla run_analysis.R --prepare-only
```

The archive contains no participant-level data, completed imputations, fitted
objects, credentials, or restricted CHARLS files. CHARLS data must be obtained
under the [CHARLS data-use terms](https://charls.pku.edu.cn/); NHANES public-use
files are available from [NCHS](https://www.cdc.gov/nchs/nhanes/).
