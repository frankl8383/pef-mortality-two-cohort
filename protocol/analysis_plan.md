# Analysis Plan V2.1

## Project Goal

Build a multi-database respiratory vulnerability axis study connecting community-level chronic respiratory vulnerability with ICU-level non-invasive respiratory support failure.

## Primary Design Principle

Do not force all databases into one artificial table. Each database should contribute to a defined evidence layer, with harmonization at the construct level rather than by inventing shared variables.

## Construct Definitions

### CRVI

Community Respiratory Vulnerability Index. Domain-level construct for CHARLS and NHANES:

- Low respiratory reserve.
- Frailty or deficit accumulation.
- Inflammation-metabolic imbalance.
- Social determinants and behavioral risk, used carefully as covariates or subgroup/equity factors.

### ARVI

Acute Respiratory Vulnerability Index. ICU construct for MIMIC-IV and eICU:

- Acute oxygenation and ventilatory load.
- Chronic lung disease proxies and support intensity.
- Deficit proxies and comorbidity burden.
- Inflammation-metabolic and organ dysfunction markers.

## Core Aims

### Aim 0: Burden Layer

Describe chronic respiratory disease burden and geographic disparities using GBD, CDC WONDER, BRFSS, and PLACES. Do not make individual-level causal claims from ecological data.

### Aim 1: Community Discovery

Evaluate whether CRVI predicts chronic respiratory disease incidence, lung function decline, frailty progression, and transitions toward death or disability in CHARLS.

### Aim 2: Cross-National Replication

Replicate CRVI associations with objective respiratory and functional outcomes in NHANES using survey weights, strata, and PSU.

### Aim 3: Genetic Causality

Use MR, MVMR, bidirectional MR, and selected colocalization to test whether respiratory reserve, inflammation, metabolic traits, smoking, BMI, frailty, and functional limitation show directionally coherent genetic evidence.

### Aim 4: Mechanism and Targets

Use processed GEO, SRA, Expression Atlas, and optional TCGA/GDC open data to evaluate tissue and cell plausibility. Use Open Targets and ChEMBL for target prioritization.

### Aim 5: ICU Translation

Develop MIMIC-IV models and externally validate in eICU for 48-hour NIRS failure after HFNC, NIV, or high-concentration oxygen support.

### Aim 6: ICU Target Trial Emulation

Emulate early invasive mechanical ventilation versus continued non-invasive respiratory support among eligible patients at the 6-hour landmark.

## ICU Time Anchors

- T0: first HFNC, NIV, or high-concentration oxygen support.
- L2: T0 plus 2 hours.
- L6: T0 plus 6 hours.
- Prediction window: T0 to L6.
- Main prediction outcome: 48-hour NIRS failure, defined as invasive mechanical ventilation or death.
- TTE time zero: L6.
- Strategy A: invasive mechanical ventilation between L6 and L12.
- Strategy B: continue non-invasive respiratory support between L6 and L12, with rescue intubation allowed after L12.

## Ticket 01 Deliverables

- Repository skeleton.
- `README.md`.
- `protocol/analysis_plan.md`.
- `data_manifest/data_access_log.csv`.
- `data_dict/variable_dictionary_template.csv`.
- `data_dict/harmonization_map_template.csv`.
- `environment.yml`.
- `renv.lock` placeholder.
- QA and data-governance notes.

## QA Gates

Before any model is run:

1. Every dataset must have a manifest entry.
2. Every analysis variable must have a dictionary entry.
3. Every uncertain variable must be marked `TODO_CODEBOOK_CHECK`.
4. Every sensitive data path must be represented by an environment variable or local config ignored by git.
5. MIMIC/eICU work must begin with schema checks only.

## Go/No-Go Checks

- CHARLS: sample size, event count, PEF/FI reproducibility.
- NHANES: objective respiratory outcome availability and complete survey design fields.
- Genetics: sufficient instruments and non-conflicting sensitivity results.
- Transcriptomics: at least two respiratory datasets and one immune, aging, or muscle dataset with coherent direction.
- ICU: stable identification of T0, respiratory support, invasive ventilation, and death in both MIMIC-IV and eICU.
- TTE: clinical equipoise, adequate events in both strategies, and stable weight distributions.

