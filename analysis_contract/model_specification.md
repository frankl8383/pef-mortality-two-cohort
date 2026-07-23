# Row-level model-specification summary

This file summarizes the model stage used to create the disclosure-safe
aggregate inputs in `results/`. It supports methodological audit and
independent reimplementation planning; it is not a complete executable
raw-data pipeline or a redistribution of CHARLS or NHANES data.

Natural-spline basis variables were created before fitting. Exact spline
boundaries and knots, factor reference levels, the complete multiple-imputation
predictor matrix, and source-file harmonization code are not asserted here as
fully reconstructive metadata. The manuscript and supplementary methods remain
the authoritative description of the submitted analysis.

## Common exposure definitions

Within each cohort, peak expiratory flow (PEF) was standardized separately for
men and women using the design-weighted mean and standard deviation:

```text
E0 = (sex-specific weighted mean PEF - individual PEF) /
     sex-specific weighted SD PEF
E0b = -PEF / 100
```

Larger values therefore indicate lower PEF. Constants frozen in the primary
samples were:

| Cohort | Sex | Weighted mean, L/min | Weighted SD, L/min |
|---|---|---:|---:|
| CHARLS | Men | 352.0554437 | 139.6633763 |
| CHARLS | Women | 250.4158738 | 99.8964517 |
| NHANES | Men | 544.8197628 | 123.5911609 |
| NHANES | Women | 390.0234624 | 87.4123993 |

The same constants were retained in the corresponding measurement-quality
sensitivities.

## CHARLS

### Population and time structure

- Baseline: Wave 1, age at least 45 years.
- Required: positive biomarker weight, known sex and community, valid Wave 1
  maximum PEF from 30 to 890 L/min, and at least one explicit vital-status
  observation in Waves 2–5.
- Outcome data: one row per observed survey interval.
- Event: first newly reported death.
- Baseline hazard: indicator terms for transition end wave.
- Time offset: `log(nominal_years)`.
- Survey design: baseline biomarker weight and participant identifiers nested
  within community; no invented public stratum.
- Model: survey-weighted binomial generalized linear model with
  complementary-log-log link.

### Formula tiers

```r
# A0
period_event ~ E0 +
  age_ns1 + age_ns2 + age_ns3 + factor(sex_code) +
  height_ns1 + height_ns2 + height_ns3 + factor(end_wave) +
  offset(log(nominal_years))

# A1 (primary)
period_event ~ E0 +
  age_ns1 + age_ns2 + age_ns3 + factor(sex_code) +
  height_ns1 + height_ns2 + height_ns3 + factor(end_wave) +
  offset(log(nominal_years)) +
  factor(smoking_status_3) +
  bmi_ns1 + bmi_ns2 + bmi_ns3 +
  factor(raeducl) + factor(h1rural) + log_hh1cperc

# A2
update(A1, . ~ . +
  lung_w1 + asthma_w1 + nonpulmonary_comorbidity_count +
  frailty_proxy_count_w1)
```

The same tiers were fitted after replacing `E0` with `E0b`.

## NHANES

### Population and time structure

- Cycles: 2007–2008, 2009–2010, and 2011–2012.
- Age: 45–79 years.
- Required: Mobile Examination Center examination, mortality linkage,
  linkage eligibility, no official spirometry safety exclusion, positive PEF,
  positive examination weight, positive follow-up, and grade-A spirometry.
- Follow-up: examination to the public-use linked mortality time through 2019.
- Survey design: six-year MEC weight (two-year weight divided by three), 45
  design strata, and 94 primary sampling units, created before restricting to
  the grade-A analytic domain.
- Model: survey-weighted Cox proportional-hazards regression.

### Formula tiers

```r
# A0
survival::Surv(followup_years, death) ~ E0 +
  age_ns1 + age_ns2 + age_ns3 + factor(RIAGENDR) +
  height_ns1 + height_ns2 + height_ns3 +
  factor(race_ethnicity) + factor(cycle_label)

# A1 (primary)
survival::Surv(followup_years, death) ~ E0 +
  age_ns1 + age_ns2 + age_ns3 + factor(RIAGENDR) +
  height_ns1 + height_ns2 + height_ns3 +
  factor(race_ethnicity) + factor(cycle_label) +
  factor(smoking_status) +
  bmi_ns1 + bmi_ns2 + bmi_ns3 +
  factor(education) + income_poverty_ratio

# A2
update(A1, . ~ . +
  self_reported_emphysema_or_bronchitis + ever_asthma +
  nhanes_function_health_proxy_count)
```

The same tiers were fitted after replacing `E0` with `E0b`.

## Missing data and combination of estimates

- Fifty imputations were created separately within each cohort.
- PEF, survey weights, outcome status, and follow-up information were included
  in the imputation models.
- Variables constructed from completed parent components were recalculated
  within each imputation.
- The survey model was fitted in every completed dataset.
- Log-hazard estimates and covariance matrices were combined with Rubin's
  rules; hazard ratios and confidence limits were exponentiated afterward.

## Sensitivity definitions

The registered aggregate rows cover complete cases, conditional survivors,
alternative PEF-quality rules, CHARLS anthropometric-quality rules, NHANES
values above 900 L/min, and leave-one-cycle-out analyses. Exact labels,
participant counts, deaths, estimates, and confidence intervals are in the
aggregate registries and the rebuilt Table 3 source data.
