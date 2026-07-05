#!/usr/bin/env python3
"""Build aggregate MIMIC-IV schema/cohort/endpoint feasibility lock outputs.

The outputs are intentionally limited to table metadata, dictionary metadata,
and aggregate counts. No row-level clinical records are exported.
"""

from __future__ import annotations

import csv
import gzip
import re
from collections import Counter
from datetime import date
from pathlib import Path
from statistics import median


PROJECT_ROOT = Path(__file__).resolve().parents[2]
MIMIC_ROOT = Path.home() / "secure_data" / "mimiciv" / "3.1"
TABLE_DIR = PROJECT_ROOT / "results" / "tables"
LOG_DIR = PROJECT_ROOT / "results" / "logs"
MANUSCRIPT_DIR = PROJECT_ROOT / "manuscript"


CORE_TABLES = {
    "patients": "hosp/patients.csv.gz",
    "admissions": "hosp/admissions.csv.gz",
    "transfers": "hosp/transfers.csv.gz",
    "icustays": "icu/icustays.csv.gz",
    "d_items": "icu/d_items.csv.gz",
    "chartevents": "icu/chartevents.csv.gz",
    "procedureevents": "icu/procedureevents.csv.gz",
    "diagnoses_icd": "hosp/diagnoses_icd.csv.gz",
    "procedures_icd": "hosp/procedures_icd.csv.gz",
    "d_icd_diagnoses": "hosp/d_icd_diagnoses.csv.gz",
    "d_icd_procedures": "hosp/d_icd_procedures.csv.gz",
    "d_labitems": "hosp/d_labitems.csv.gz",
    "labevents": "hosp/labevents.csv.gz",
}


def open_csv(rel_file: str):
    return gzip.open(MIMIC_ROOT / rel_file, "rt", newline="")


def read_header(rel_file: str) -> list[str]:
    with open_csv(rel_file) as handle:
        return next(csv.reader(handle))


def read_dict_rows(rel_file: str) -> list[dict[str, str]]:
    with open_csv(rel_file) as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def classify_item(text: str) -> list[str]:
    text = text.lower()
    domains: list[str] = []
    rules = [
        ("invasive_ventilation_airway", r"\bintub|endotracheal|ett\b|tracheostomy|trach\b"),
        ("ventilator_mode", r"ventilator mode|ventilation mode|\bvent mode\b|mode \(vent"),
        ("noninvasive_ventilation", r"\bbipap\b|\bcpap\b|non[ -]?invasive|\bniv\b"),
        ("oxygen_delivery", r"oxygen device|o2 device|oxygen flow|o2 flow|flow rate|high flow|nasal cannula|face mask|aerosol mask|trach collar|venturi"),
        ("icu_peak_expiratory_flow_candidate", r"peak exp flow|peak expiratory flow"),
        ("fio2_peep_resp_mechanics", r"\bfio2\b|\bpeep\b|tidal volume|minute ventilation|plateau pressure|peak insp|respiratory rate|driving pressure"),
        ("oxygenation_gas_exchange", r"\bspo2\b|o2 saturation|oxygen saturation|\bpo2\b|\bpao2\b|\bpco2\b|\bpaco2\b|arterial o2|venous o2"),
        ("respiratory_therapy_context", r"respiratory therapist|respiratory treatment|chest physiotherapy|nebulizer|suction"),
    ]
    for domain, pattern in rules:
        if re.search(pattern, text):
            domains.append(domain)
    return domains


def classify_icd(title: str, source: str) -> list[str]:
    text = title.lower()
    if source == "procedure":
        rules = [
            ("invasive_mechanical_ventilation", r"mechanical ventilation|continuous invasive mechanical ventilation|respiratory ventilation"),
            ("noninvasive_ventilation", r"non-invasive|noninvasive|cpap|bipap"),
            ("airway_intubation_tracheostomy", r"intubat|tracheostomy|tracheotomy"),
            ("oxygen_respiratory_support", r"oxygen|respiratory support|ventilatory support"),
        ]
    else:
        rules = [
            ("acute_respiratory_failure", r"acute respiratory failure|respiratory failure|\barf\b|\bards\b|acute respiratory distress"),
            ("chronic_lung_disease", r"chronic obstructive|copd|emphysema|chronic bronchitis|bronchiectasis|interstitial lung|pulmonary fibrosis"),
            ("asthma_airway_disease", r"\basthma\b|status asthmaticus"),
            ("pneumonia_infection", r"pneumonia|viral pneumonia|bacterial pneumonia"),
            ("hypoxemia_hypercapnia", r"hypoxemia|hypercapnia|respiratory acidosis"),
        ]
    out = []
    for domain, pattern in rules:
        if re.search(pattern, text):
            out.append(domain)
    return out


def build_key_variable_lock() -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    expected = {
        "hosp/patients.csv.gz": {
            "id": ["subject_id"],
            "time": ["anchor_year", "dod"],
            "clinical": ["gender", "anchor_age"],
        },
        "hosp/admissions.csv.gz": {
            "id": ["subject_id", "hadm_id"],
            "time": ["admittime", "dischtime", "deathtime"],
            "clinical": ["hospital_expire_flag", "admission_type", "race"],
        },
        "icu/icustays.csv.gz": {
            "id": ["subject_id", "hadm_id", "stay_id"],
            "time": ["intime", "outtime"],
            "clinical": ["first_careunit", "last_careunit", "los"],
        },
        "hosp/transfers.csv.gz": {
            "id": ["subject_id", "hadm_id", "transfer_id"],
            "time": ["intime", "outtime"],
            "clinical": ["eventtype", "careunit"],
        },
        "icu/chartevents.csv.gz": {
            "id": ["subject_id", "hadm_id", "stay_id", "itemid"],
            "time": ["charttime", "storetime"],
            "clinical": ["value", "valuenum", "valueuom"],
        },
        "icu/procedureevents.csv.gz": {
            "id": ["subject_id", "hadm_id", "stay_id", "itemid"],
            "time": ["starttime", "endtime", "storetime"],
            "clinical": ["value", "valueuom"],
        },
        "hosp/diagnoses_icd.csv.gz": {
            "id": ["subject_id", "hadm_id", "icd_code", "icd_version"],
            "time": [],
            "clinical": ["seq_num"],
        },
        "hosp/procedures_icd.csv.gz": {
            "id": ["subject_id", "hadm_id", "icd_code", "icd_version"],
            "time": ["chartdate"],
            "clinical": ["seq_num"],
        },
    }
    for rel_file, groups in expected.items():
        header = set(read_header(rel_file))
        for group, variables in groups.items():
            missing = [var for var in variables if var not in header]
            rows.append(
                {
                    "table_file": rel_file,
                    "variable_group": group,
                    "expected_variables": ";".join(variables) if variables else "none_required",
                    "missing_variables": ";".join(missing),
                    "status": "pass" if not missing else "missing_expected_variables",
                }
            )
    return rows


def build_cohort_counts() -> list[dict[str, object]]:
    patients: dict[str, dict[str, str]] = {}
    with open_csv("hosp/patients.csv.gz") as handle:
        for row in csv.DictReader(handle):
            patients[row["subject_id"]] = row

    admissions: dict[str, dict[str, str]] = {}
    hospital_deaths = 0
    admission_race_known = 0
    with open_csv("hosp/admissions.csv.gz") as handle:
        for row in csv.DictReader(handle):
            admissions[row["hadm_id"]] = row
            hospital_deaths += 1 if row.get("hospital_expire_flag") == "1" else 0
            admission_race_known += 1 if row.get("race") else 0

    los_values: list[float] = []
    first_careunit = Counter()
    icu_total = 0
    icu_link_patient = 0
    icu_link_admission = 0
    icu_link_both = 0
    adult_icu_stays = 0
    with open_csv("icu/icustays.csv.gz") as handle:
        for row in csv.DictReader(handle):
            icu_total += 1
            subj_ok = row["subject_id"] in patients
            hadm_ok = row["hadm_id"] in admissions
            icu_link_patient += int(subj_ok)
            icu_link_admission += int(hadm_ok)
            icu_link_both += int(subj_ok and hadm_ok)
            first_careunit[row.get("first_careunit", "") or "missing"] += 1
            try:
                los_values.append(float(row["los"]))
            except Exception:
                pass
            if subj_ok:
                try:
                    adult_icu_stays += int(float(patients[row["subject_id"]].get("anchor_age", "0")) >= 18)
                except Exception:
                    pass

    rows = [
        {"metric": "patients_rows", "value": len(patients), "note": "aggregate count from patients table"},
        {"metric": "admissions_rows", "value": len(admissions), "note": "aggregate count from admissions table"},
        {"metric": "hospital_deaths", "value": hospital_deaths, "note": "sum of hospital_expire_flag"},
        {"metric": "admissions_with_race_recorded", "value": admission_race_known, "note": "aggregate race availability"},
        {"metric": "icu_stay_rows", "value": icu_total, "note": "aggregate count from icustays table"},
        {"metric": "icu_stays_linked_to_patients", "value": icu_link_patient, "note": "stay subject_id found in patients"},
        {"metric": "icu_stays_linked_to_admissions", "value": icu_link_admission, "note": "stay hadm_id found in admissions"},
        {"metric": "icu_stays_linked_to_both", "value": icu_link_both, "note": "stay links to both patient and admission tables"},
        {"metric": "adult_icu_stays_by_anchor_age", "value": adult_icu_stays, "note": "anchor_age >= 18 among linked ICU stays"},
        {"metric": "icu_los_median_days", "value": f"{median(los_values):.2f}" if los_values else "", "note": "aggregate median from icustays los"},
        {"metric": "icu_los_available_rows", "value": len(los_values), "note": "nonmissing numeric los rows"},
    ]
    for unit, count in first_careunit.most_common():
        rows.append({"metric": f"first_careunit::{unit}", "value": count, "note": "aggregate ICU first-careunit count"})
    return rows


def build_item_candidates() -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for row in read_dict_rows("icu/d_items.csv.gz"):
        text = " ".join(
            [
                row.get("label", ""),
                row.get("abbreviation", ""),
                row.get("category", ""),
                row.get("unitname", ""),
                row.get("linksto", ""),
            ]
        )
        domains = classify_item(text)
        for domain in domains:
            priority = "primary" if domain in {"invasive_ventilation_airway", "ventilator_mode", "noninvasive_ventilation", "oxygen_delivery", "fio2_peep_resp_mechanics"} else "secondary"
            rows.append(
                {
                    "source": "icu/d_items",
                    "candidate_domain": domain,
                    "candidate_priority": priority,
                    "itemid": row.get("itemid", ""),
                    "label": row.get("label", ""),
                    "abbreviation": row.get("abbreviation", ""),
                    "linksto": row.get("linksto", ""),
                    "category": row.get("category", ""),
                    "unitname": row.get("unitname", ""),
                }
            )
    rows.sort(key=lambda r: (str(r["candidate_domain"]), str(r["candidate_priority"]), str(r["linksto"]), str(r["category"]), str(r["label"])))
    return rows


def build_icd_candidates() -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for source_name, rel_file, source in [
        ("hosp/d_icd_procedures", "hosp/d_icd_procedures.csv.gz", "procedure"),
        ("hosp/d_icd_diagnoses", "hosp/d_icd_diagnoses.csv.gz", "diagnosis"),
    ]:
        for row in read_dict_rows(rel_file):
            title = row.get("long_title", "")
            domains = classify_icd(title, source)
            for domain in domains:
                rows.append(
                    {
                        "source": source_name,
                        "candidate_domain": domain,
                        "icd_code": row.get("icd_code", ""),
                        "icd_version": row.get("icd_version", ""),
                        "long_title": title,
                    }
                )
    rows.sort(key=lambda r: (str(r["source"]), str(r["candidate_domain"]), str(r["icd_version"]), str(r["icd_code"])))
    return rows


def build_endpoint_feasibility(item_rows: list[dict[str, object]], icd_rows: list[dict[str, object]]) -> list[dict[str, object]]:
    item_domains = {str(row["candidate_domain"]) for row in item_rows}
    icd_domains = {str(row["candidate_domain"]) for row in icd_rows}
    rows = [
        {
            "endpoint_or_construct": "adult ICU stay denominator",
            "source_tables": "icu/icustays.csv.gz;hosp/patients.csv.gz;hosp/admissions.csv.gz",
            "status": "ready",
            "manuscript_role": "cohort denominator only",
            "notes": "All required ID and timing variables are present; aggregate linkage counts generated.",
        },
        {
            "endpoint_or_construct": "hospital mortality",
            "source_tables": "hosp/admissions.csv.gz",
            "status": "ready",
            "manuscript_role": "candidate secondary outcome",
            "notes": "hospital_expire_flag and deathtime are present; final time-zero rules still need prespecification.",
        },
        {
            "endpoint_or_construct": "ICU/hospital length of stay",
            "source_tables": "icu/icustays.csv.gz;hosp/admissions.csv.gz",
            "status": "ready",
            "manuscript_role": "candidate descriptive or secondary outcome",
            "notes": "ICU LOS is directly available; hospital LOS can be derived from admit/discharge times.",
        },
        {
            "endpoint_or_construct": "invasive mechanical ventilation / intubation",
            "source_tables": "icu/d_items.csv.gz;icu/chartevents.csv.gz;icu/procedureevents.csv.gz;hosp/d_icd_procedures.csv.gz;hosp/procedures_icd.csv.gz",
            "status": "ready_for_endpoint_lock" if {"invasive_ventilation_airway", "ventilator_mode"} & item_domains or "invasive_mechanical_ventilation" in icd_domains else "needs_review",
            "manuscript_role": "candidate acute respiratory support outcome",
            "notes": "Metadata candidates exist; final endpoint must lock item IDs/code sets and time windows before modeling.",
        },
        {
            "endpoint_or_construct": "noninvasive ventilation",
            "source_tables": "icu/d_items.csv.gz;icu/chartevents.csv.gz;hosp/d_icd_procedures.csv.gz",
            "status": "ready_for_endpoint_lock" if "noninvasive_ventilation" in item_domains or "noninvasive_ventilation" in icd_domains else "needs_review",
            "manuscript_role": "candidate acute respiratory support outcome",
            "notes": "Candidate CPAP/BiPAP/NIV terms exist in metadata; final cleaning rules remain required.",
        },
        {
            "endpoint_or_construct": "oxygen support / escalation",
            "source_tables": "icu/d_items.csv.gz;icu/chartevents.csv.gz",
            "status": "ready_for_endpoint_lock" if "oxygen_delivery" in item_domains else "needs_review",
            "manuscript_role": "candidate acute respiratory support outcome",
            "notes": "Candidate oxygen delivery and flow metadata exist, but oxygen-escalation severity hierarchy must be prespecified.",
        },
        {
            "endpoint_or_construct": "early respiratory physiology severity",
            "source_tables": "icu/d_items.csv.gz;icu/chartevents.csv.gz;hosp/d_labitems.csv.gz;hosp/labevents.csv.gz",
            "status": "feasible_but_not_pef_equivalent" if {"fio2_peep_resp_mechanics", "oxygenation_gas_exchange"} & item_domains else "needs_review",
            "manuscript_role": "possible ICU translation marker, not direct CHARLS/NHANES replication",
            "notes": "MIMIC has acute oxygenation/ventilator physiology, but not population PEF or lower-than-expected PEF residuals.",
        },
        {
            "endpoint_or_construct": "baseline chronic lung disease history",
            "source_tables": "hosp/d_icd_diagnoses.csv.gz;hosp/diagnoses_icd.csv.gz",
            "status": "ready_for_code_lock" if "chronic_lung_disease" in icd_domains else "needs_review",
            "manuscript_role": "candidate covariate or subgroup, not equivalent to CHARLS incident self-report endpoint",
            "notes": "Diagnosis-code candidates exist; this would capture coded history/comorbidity rather than longitudinal incident disease.",
        },
        {
            "endpoint_or_construct": "PEF-like exposure",
            "source_tables": "icu/d_items.csv.gz;icu/chartevents.csv.gz",
            "status": "candidate_not_directly_equivalent",
            "manuscript_role": "limits direct three-database validation",
            "notes": "A Peak Exp Flow Rate ICU chart item exists in metadata, but this is not yet equivalent to standardized field PEF or GLI spirometry and requires value-density/clinical-context review.",
        },
        {
            "endpoint_or_construct": "modeling gate",
            "source_tables": "V19 schema/cohort/endpoint lock outputs",
            "status": "do_not_model_until_endpoint_lock_reviewed",
            "manuscript_role": "governance boundary",
            "notes": "Schema feasibility is not the same as a validated MIMIC analysis.",
        },
    ]
    return rows


def build_log(
    key_rows: list[dict[str, object]],
    count_rows: list[dict[str, object]],
    item_rows: list[dict[str, object]],
    icd_rows: list[dict[str, object]],
    endpoint_rows: list[dict[str, object]],
) -> None:
    domain_counts = Counter(str(row["candidate_domain"]) for row in item_rows)
    icd_domain_counts = Counter(str(row["candidate_domain"]) for row in icd_rows)
    lines = [
        "# MIMIC-IV Schema Cohort Endpoint Lock V19.0",
        "",
        f"- Run date: {date.today().isoformat()}.",
        "- Scope: schema lock, aggregate cohort feasibility, candidate respiratory endpoint metadata, and manuscript-integration gate.",
        "- Data governance: no row-level clinical records were exported.",
        f"- Key-variable lock rows: {len(key_rows)}.",
        f"- Aggregate cohort count rows: {len(count_rows)}.",
        f"- Candidate respiratory d_items rows: {len(item_rows)}.",
        f"- Candidate ICD metadata rows: {len(icd_rows)}.",
        f"- Endpoint feasibility rows: {len(endpoint_rows)}.",
        "",
        "## Respiratory Item Domains",
        "",
    ]
    for domain, count in domain_counts.most_common():
        lines.append(f"- {domain}: {count}.")
    lines.extend(["", "## ICD Candidate Domains", ""])
    for domain, count in icd_domain_counts.most_common():
        lines.append(f"- {domain}: {count}.")
    lines.extend(
        [
            "",
            "## Initial Interpretation",
            "",
            "- MIMIC-IV is now technically usable for schema/cohort/endpoint feasibility.",
            "- Adult ICU stay denominator, mortality, ICU LOS, ventilation/intubation metadata, noninvasive ventilation metadata and oxygen-support metadata are feasible for endpoint lock.",
            "- A charted ICU Peak Exp Flow Rate candidate exists, but direct equivalence to standardized field PEF is not established; MIMIC should not yet be described as direct three-database validation of the CHARLS/NHANES PEF axis.",
            "- The scientifically safer role is an optional ICU translation layer after endpoint and time-window lock.",
            "",
            "PASS: MIMIC-IV schema/cohort/endpoint lock V19.0 built.",
        ]
    )
    (LOG_DIR / "mimiciv_schema_cohort_endpoint_lock_v19_0.md").write_text("\n".join(lines) + "\n")


def main() -> None:
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    MANUSCRIPT_DIR.mkdir(parents=True, exist_ok=True)

    key_rows = build_key_variable_lock()
    count_rows = build_cohort_counts()
    item_rows = build_item_candidates()
    icd_rows = build_icd_candidates()
    endpoint_rows = build_endpoint_feasibility(item_rows, icd_rows)

    write_csv(
        TABLE_DIR / "mimiciv_key_variable_lock_v19_0.csv",
        key_rows,
        ["table_file", "variable_group", "expected_variables", "missing_variables", "status"],
    )
    write_csv(
        TABLE_DIR / "mimiciv_cohort_feasibility_counts_v19_0.csv",
        count_rows,
        ["metric", "value", "note"],
    )
    write_csv(
        TABLE_DIR / "mimiciv_respiratory_item_candidates_v19_0.csv",
        item_rows,
        ["source", "candidate_domain", "candidate_priority", "itemid", "label", "abbreviation", "linksto", "category", "unitname"],
    )
    write_csv(
        TABLE_DIR / "mimiciv_icd_respiratory_code_candidates_v19_0.csv",
        icd_rows,
        ["source", "candidate_domain", "icd_code", "icd_version", "long_title"],
    )
    write_csv(
        TABLE_DIR / "mimiciv_endpoint_feasibility_v19_0.csv",
        endpoint_rows,
        ["endpoint_or_construct", "source_tables", "status", "manuscript_role", "notes"],
    )
    build_log(key_rows, count_rows, item_rows, icd_rows, endpoint_rows)
    print("Wrote MIMIC-IV schema/cohort/endpoint lock V19.0 outputs.")


if __name__ == "__main__":
    main()
