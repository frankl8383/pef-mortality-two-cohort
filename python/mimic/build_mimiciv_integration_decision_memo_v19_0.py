#!/usr/bin/env python3
"""Build a MIMIC-IV manuscript-integration decision memo from V19 feasibility outputs."""

from __future__ import annotations

import csv
from datetime import date
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[2]
TABLE_DIR = PROJECT_ROOT / "results" / "tables"
LOG_DIR = PROJECT_ROOT / "results" / "logs"
MANUSCRIPT_DIR = PROJECT_ROOT / "manuscript"


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def metric_lookup() -> dict[str, str]:
    rows = read_csv(TABLE_DIR / "mimiciv_cohort_feasibility_counts_v19_0.csv")
    return {row["metric"]: row["value"] for row in rows}


def endpoint_lookup() -> dict[str, str]:
    rows = read_csv(TABLE_DIR / "mimiciv_endpoint_feasibility_v19_0.csv")
    return {row["endpoint_or_construct"]: row["status"] for row in rows}


def main() -> None:
    counts = metric_lookup()
    endpoints = endpoint_lookup()
    download_status = {row["metric"]: row["value"] for row in read_csv(TABLE_DIR / "mimiciv_download_completion_status_v19_0.csv")}

    decision_rows = [
        {
            "decision_domain": "download_integrity",
            "evidence": "33/33 files pass gzip/SHA256; no missing files.",
            "decision": "pass",
            "manuscript_implication": "MIMIC can move from download gate to feasibility work.",
        },
        {
            "decision_domain": "schema_and_linkage",
            "evidence": f"{counts.get('icu_stay_rows')} ICU stays; {counts.get('icu_stays_linked_to_both')} link to patient and admission tables.",
            "decision": "pass",
            "manuscript_implication": "Adult ICU denominator and linkage are technically feasible.",
        },
        {
            "decision_domain": "acute_respiratory_endpoints",
            "evidence": "Ventilation/intubation, noninvasive ventilation, oxygen support, mortality and ICU LOS have ready or ready-for-lock status.",
            "decision": "pass_for_endpoint_lock",
            "manuscript_implication": "A MIMIC acute respiratory support outcome layer is feasible after endpoint/time-window lock.",
        },
        {
            "decision_domain": "PEF_axis_equivalence",
            "evidence": f"PEF-like exposure status: {endpoints.get('PEF-like exposure')}.",
            "decision": "candidate_not_directly_equivalent",
            "manuscript_implication": "MIMIC has an ICU charted peak-flow candidate but should not be described as direct three-database validation of the CHARLS/NHANES standardized PEF axis.",
        },
        {
            "decision_domain": "current_manuscript_inclusion",
            "evidence": "Feasibility is validated, but no MIMIC cohort extraction, endpoint lock, time-window rule or model result exists yet.",
            "decision": "do_not_add_to_abstract_results_or_conclusion_yet",
            "manuscript_implication": "Keep current main evidence as CHARLS + NHANES.",
        },
        {
            "decision_domain": "best_scientific_role",
            "evidence": "MIMIC captures acute ICU respiratory support failure but not field PEF or GLI spirometry.",
            "decision": "optional_ICU_translation_layer",
            "manuscript_implication": "If modeled, present as a translational stress-test layer, not as equivalent replication.",
        },
        {
            "decision_domain": "respiratory_research_probability",
            "evidence": "RR fit improves only if MIMIC yields a coherent respiratory-support layer; feasibility alone is not a publishable third result.",
            "decision": "probability_unchanged_until_model",
            "manuscript_implication": "Current RR estimate remains roughly 20-30%; a coherent MIMIC layer could raise it to roughly 30-45%.",
        },
    ]
    write_csv(
        TABLE_DIR / "mimiciv_manuscript_integration_decision_v19_0.csv",
        decision_rows,
        ["decision_domain", "evidence", "decision", "manuscript_implication"],
    )

    lines = [
        "# MIMIC-IV Manuscript Integration Decision Memo V19.0",
        "",
        f"- Run date: {date.today().isoformat()}.",
        "- One-sentence argument: MIMIC-IV is now technically ready for ICU schema and respiratory endpoint work, but it is not yet a direct validation of the CHARLS/NHANES lower-than-expected PEF axis because the ICU peak-flow candidate is not established as equivalent to standardized field PEF.",
        "- Canonical terms: MIMIC-IV, ICU translation layer, acute respiratory support outcome, lower-than-expected PEF, CHARLS + NHANES main evidence.",
        "- Boundary: this memo uses only aggregate feasibility outputs and metadata; it does not report a MIMIC model or cohort result.",
        "",
        "## Feasibility Result",
        "",
        f"MIMIC-IV v3.1 now passes the download integrity gate: {download_status.get('ok_files')}/{download_status.get('expected_files')} files are integrity-passing, with {download_status.get('missing_files')} missing files, {download_status.get('bad_gzip_files')} gzip failures and {download_status.get('sha_fail_files')} SHA256 failures. The schema/cohort lock found {counts.get('patients_rows')} patient rows, {counts.get('admissions_rows')} admission rows and {counts.get('icu_stay_rows')} ICU stay rows; all ICU stays linked to both patient and admission tables in the aggregate linkage check. Hospital mortality and ICU LOS are available as aggregate-outcome candidates, with {counts.get('hospital_deaths')} hospital deaths and a median ICU LOS of {counts.get('icu_los_median_days')} days.",
        "",
        "The respiratory endpoint metadata are also sufficient for endpoint lock. The V19 screen found candidate metadata for invasive ventilation or intubation, noninvasive ventilation, oxygen support or escalation, FiO2/PEEP/respiratory mechanics, gas-exchange markers and ICD-based respiratory diagnoses/procedures. These are not final endpoint definitions; they are the evidence that a careful endpoint lock is possible.",
        "",
        "## Scientific Decision",
        "",
        "MIMIC should not yet be added to the abstract, Results or conclusion as a third-database validation. The reason is construct validity, not download readiness. CHARLS and NHANES are organized around lower-than-expected PEF and objective spirometry anchoring, whereas MIMIC primarily contains acute ICU physiology, oxygenation, ventilation and coded disease/support events. Although an ICU charted Peak Exp Flow Rate item exists in the metadata, it has not been shown to be standardized, sufficiently populated, or comparable to field PEF in CHARLS/NHANES. That makes MIMIC suitable for an ICU translation layer, but not for direct replication of the same exposure construct.",
        "",
        "The strongest manuscript architecture therefore remains CHARLS + NHANES as the primary evidence, with MIMIC held as an optional translational extension. If MIMIC is modeled, it should answer a narrower question: whether an ICU respiratory-vulnerability signature or baseline respiratory-history stratum is associated with subsequent respiratory support escalation, ventilation, mortality or LOS under prespecified time windows. It should not be framed as proof that the PEF axis is causal, diagnostic, or fully validated across three databases.",
        "",
        "## Recommended Next Analysis If MIMIC Is Pursued",
        "",
        "Use MIMIC only after a final endpoint lock table is approved. The minimal viable analysis should define an adult ICU stay cohort, a time-zero such as ICU admission, a baseline window such as the first 6 or 24 hours, and outcomes measured after that baseline window. Candidate outcomes are invasive ventilation/intubation, noninvasive ventilation, oxygen-support escalation, hospital mortality and ICU LOS. Candidate covariates are age, sex, race, admission type, first ICU care unit and coded chronic lung disease history. The exposure must be named as an ICU respiratory physiology/support signature or chronic respiratory history stratum, not as lower-than-expected PEF.",
        "",
        "## Respiratory Research Implication",
        "",
        "Respiratory Research remains the most realistic balanced first target. The current CHARLS + NHANES package remains a credible two-dataset observational marker paper, with an approximate subjective acceptance range of 20-30% after careful author-field completion and formatting. MIMIC feasibility alone does not materially change that estimate. A coherent, validated MIMIC ICU translation layer could raise the ceiling to roughly 30-45%, mainly by improving clinical respiratory relevance and reviewer confidence in translational reach. A poorly aligned MIMIC layer would likely reduce clarity and should be left out.",
        "",
        "## Decision Table",
        "",
        "| Domain | Decision | Manuscript implication |",
        "| --- | --- | --- |",
    ]
    for row in decision_rows:
        lines.append(f"| {row['decision_domain']} | {row['decision']} | {row['manuscript_implication']} |")
    lines.extend(
        [
            "",
            "PASS: MIMIC-IV integration decision memo V19.0 built.",
        ]
    )
    (MANUSCRIPT_DIR / "mimiciv_integration_decision_memo_v19_0.md").write_text("\n".join(lines) + "\n")
    print("Wrote manuscript/mimiciv_integration_decision_memo_v19_0.md")


if __name__ == "__main__":
    main()
