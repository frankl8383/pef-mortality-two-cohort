#!/usr/bin/env python3
"""Build a design-only minimal MIMIC-IV model plan.

This does not run modeling. It records the smallest scientifically defensible
MIMIC analysis that could be executed after endpoint/time-window lock.
"""

from __future__ import annotations

import csv
from datetime import date
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[2]
TABLE_DIR = PROJECT_ROOT / "results" / "tables"
LOG_DIR = PROJECT_ROOT / "results" / "logs"
MANUSCRIPT_DIR = PROJECT_ROOT / "manuscript"


def write_csv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    rows = [
        {
            "component": "analysis_status",
            "proposed_definition": "Design-only; no model execution in V19.1.",
            "source_or_method": "V19.0 schema/cohort/endpoint lock and decision memo",
            "status": "not_executed",
            "reason": "An ICU charted peak-flow candidate exists but is not directly equivalent to standardized field PEF; final endpoint/time-window lock is still required.",
        },
        {
            "component": "scientific_role",
            "proposed_definition": "Optional ICU translation layer.",
            "source_or_method": "MIMIC captures acute ICU physiology/support; an ICU peak-flow item needs density/context review before any PEF-adjacent use.",
            "status": "allowed_with_boundary",
            "reason": "Use only as translational stress-test, not as equivalent third-database validation.",
        },
        {
            "component": "primary_cohort",
            "proposed_definition": "Adult ICU stays linked to patients and admissions; primary analysis should use first ICU stay per hospital admission or first ICU stay per patient after a documented choice.",
            "source_or_method": "icu/icustays.csv.gz, hosp/patients.csv.gz, hosp/admissions.csv.gz",
            "status": "feasible",
            "reason": "V19 aggregate linkage found all ICU stays linked to both patient and admission tables.",
        },
        {
            "component": "time_zero",
            "proposed_definition": "ICU admission time.",
            "source_or_method": "icu/icustays.csv.gz intime",
            "status": "feasible",
            "reason": "ICU intime/outtime variables passed key-variable lock.",
        },
        {
            "component": "baseline_window",
            "proposed_definition": "First 6 hours after ICU admission as primary; first 24 hours as sensitivity.",
            "source_or_method": "chartevents/procedureevents timestamps anchored to icu intime",
            "status": "requires_endpoint_lock",
            "reason": "Needed to prevent immortal-time leakage and separate baseline markers from subsequent outcomes.",
        },
        {
            "component": "exposure_option_A",
            "proposed_definition": "Baseline coded chronic lung disease history stratum.",
            "source_or_method": "hosp/diagnoses_icd.csv.gz plus locked chronic lung disease ICD code set",
            "status": "feasible_after_code_lock",
            "reason": "Most stable chronic-vulnerability proxy, but not equivalent to lower-than-expected PEF.",
        },
        {
            "component": "exposure_option_B",
            "proposed_definition": "Early ICU respiratory physiology/support signature using FiO2, SpO2, respiratory rate, PEEP, oxygen delivery and ventilation-related item groups.",
            "source_or_method": "icu/d_items.csv.gz item lock plus chartevents/procedureevents in baseline window",
            "status": "feasible_but_construct_risk",
            "reason": "Clinically respiratory, but may reflect acute severity rather than pre-illness vulnerability.",
        },
        {
            "component": "primary_outcome_option_A",
            "proposed_definition": "New invasive ventilation/intubation after the baseline window among stays without invasive ventilation during baseline.",
            "source_or_method": "locked ventilation/intubation d_items and/or ICD procedure code set",
            "status": "preferred_if_endpoint_lock_passes",
            "reason": "Closest acute respiratory support failure endpoint, but requires careful exclusion of baseline ventilated stays.",
        },
        {
            "component": "primary_outcome_option_B",
            "proposed_definition": "Oxygen-support escalation after baseline, using a prespecified oxygen/NIV/IMV hierarchy.",
            "source_or_method": "locked oxygen delivery and noninvasive/invasive support item IDs",
            "status": "secondary_or_sensitivity",
            "reason": "Potentially informative but needs hierarchy and documentation-density checks.",
        },
        {
            "component": "secondary_outcome_mortality",
            "proposed_definition": "Hospital mortality.",
            "source_or_method": "hosp/admissions.csv.gz hospital_expire_flag and deathtime",
            "status": "feasible",
            "reason": "Available aggregate outcome, but not respiratory-specific mortality.",
        },
        {
            "component": "secondary_outcome_los",
            "proposed_definition": "ICU LOS and hospital LOS.",
            "source_or_method": "icu/icustays.csv.gz los; admissions admit/discharge times",
            "status": "feasible",
            "reason": "Useful severity/resource outcomes; interpret cautiously because discharge practices confound LOS.",
        },
        {
            "component": "covariates_minimal",
            "proposed_definition": "Age, sex, race, admission type, first ICU care unit, coded chronic lung disease history; add calendar/anchor year only if needed for version/time drift.",
            "source_or_method": "patients, admissions, icustays, diagnoses_icd",
            "status": "feasible_after_code_lock",
            "reason": "Avoid overly rich adjustment that turns a feasibility analysis into a fragile model.",
        },
        {
            "component": "model_family",
            "proposed_definition": "Logistic regression for binary outcomes; linear or quantile-style summaries for LOS, with sensitivity by ICU unit and baseline window.",
            "source_or_method": "aggregate-safe derived cohort within restricted local environment",
            "status": "design_only",
            "reason": "Prefer interpretable models over black-box prediction for manuscript alignment.",
        },
        {
            "component": "reporting_boundary",
            "proposed_definition": "Report only aggregate model tables and figure-ready summaries; no row-level derived cohort export.",
            "source_or_method": "local restricted MIMIC environment",
            "status": "required",
            "reason": "Preserves data governance and manuscript evidentiary boundaries.",
        },
        {
            "component": "go_no_go_rule",
            "proposed_definition": "Run MIMIC model only if endpoint lock passes and exposure is named as ICU respiratory physiology/support signature or chronic respiratory history, not lower-than-expected PEF.",
            "source_or_method": "V19.0 decision memo",
            "status": "required_before_execution",
            "reason": "Prevents a misleading three-database claim.",
        },
    ]
    write_csv(
        TABLE_DIR / "mimiciv_minimal_model_design_v19_1.csv",
        rows,
        ["component", "proposed_definition", "source_or_method", "status", "reason"],
    )

    lines = [
        "# MIMIC-IV Minimal Model Design V19.1",
        "",
        f"- Run date: {date.today().isoformat()}.",
        "- Status: design-only; no MIMIC model has been executed.",
        "- One-sentence argument: a minimal MIMIC analysis is feasible only as an ICU translation layer, because the MIMIC peak-flow candidate is not established as equivalent to lower-than-expected field PEF.",
        "- Boundary: this design must not be presented as completed three-database validation.",
        "",
        "## Go/No-Go Position",
        "",
        "Do not run MIMIC modeling until the final endpoint lock table specifies item IDs, ICD code sets, baseline window, outcome window and exclusion rules. If the analysis proceeds, the exposure should be framed as either coded chronic respiratory history or an early ICU respiratory physiology/support signature. It should not be called lower-than-expected PEF.",
        "",
        "## Design Table",
        "",
        "| Component | Status | Proposed definition |",
        "| --- | --- | --- |",
    ]
    for row in rows:
        lines.append(f"| {row['component']} | {row['status']} | {row['proposed_definition']} |")
    lines.extend(["", "PASS: MIMIC-IV minimal model design V19.1 built."])
    (MANUSCRIPT_DIR / "mimiciv_minimal_model_design_v19_1.md").write_text("\n".join(lines) + "\n")
    print("Wrote MIMIC-IV minimal model design V19.1 outputs.")


if __name__ == "__main__":
    main()
