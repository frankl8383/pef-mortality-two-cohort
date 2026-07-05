#!/usr/bin/env python3
"""Build aggregate ICU translation-layer feasibility outputs.

This is a feasibility and variable-lock stage, not a prediction model. It
keeps row-level MIMIC/eICU data in memory only and exports aggregate CSV/MD
outputs suitable for manuscript/go-no-go review.
"""

from __future__ import annotations

import csv
import gzip
import os
import re
from collections import Counter, defaultdict
from datetime import date
from pathlib import Path

import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[2]
MIMIC_ROOT = Path.home() / "secure_data" / "mimiciv" / "3.1"
EICU_ROOT = Path.home() / "secure_data" / "eicu-crd" / "2.0"
TABLE_DIR = PROJECT_ROOT / "results" / "tables"
LOG_DIR = PROJECT_ROOT / "results" / "logs"
METADATA_DIR = PROJECT_ROOT / "metadata"
MANUSCRIPT_DIR = PROJECT_ROOT / "manuscript"
VERSION = "v21_0"

CHUNK_ROWS = int(os.environ.get("ICU_TRANSLATION_CHUNK_ROWS", "750000"))

MIMIC_PROC_NIV = {225794: "Non-invasive Ventilation"}
MIMIC_PROC_IMV = {
    225792: "Invasive Ventilation",
    224385: "Intubation",
    226237: "Open Tracheostomy",
    225448: "Percutaneous Tracheostomy",
}
MIMIC_O2_DEVICE_ITEM = 226732
MIMIC_PHYS_ITEM_MAP = {
    220277: "spo2",
    223835: "fio2",
    220210: "respiratory_rate",
    224690: "respiratory_rate",
    224689: "respiratory_rate",
    220045: "heart_rate",
    223834: "oxygen_flow",
    224691: "oxygen_flow",
    227287: "oxygen_flow",
    227582: "oxygen_flow",
    220339: "peep",
    224700: "peep",
    220181: "mean_bp",
    220179: "systolic_bp",
    223761: "temperature",
    223762: "temperature",
    220739: "gcs_any",
    223900: "gcs_any",
    223901: "gcs_any",
}

MIMIC_REQUIRED_PHYS = ["spo2", "fio2", "respiratory_rate", "heart_rate"]
MIMIC_EXTENDED_PHYS = MIMIC_REQUIRED_PHYS + ["oxygen_flow", "peep", "mean_bp", "gcs_any"]

EICU_REQUIRED_PHYS = ["sao2", "respiration", "heartrate", "fio2"]
EICU_EXTENDED_PHYS = EICU_REQUIRED_PHYS + ["peep", "lpm_o2", "vent_rate"]


def ensure_dirs() -> None:
    for directory in [TABLE_DIR, LOG_DIR, METADATA_DIR, MANUSCRIPT_DIR]:
        directory.mkdir(parents=True, exist_ok=True)


def write_csv(path: Path, rows: list[dict[str, object]], fieldnames: list[str] | None = None) -> None:
    if fieldnames is None:
        fieldnames = list(rows[0].keys()) if rows else ["status"]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def read_gzip_header(path: Path) -> list[str]:
    with gzip.open(path, "rt", newline="") as handle:
        return next(csv.reader(handle))


def norm_text(value: object) -> str:
    return str(value or "").strip().lower()


def classify_support_text(text: object) -> str:
    value = norm_text(text)
    if not value:
        return ""
    if re.search(r"\bhfnc\b|high[- ]?flow|hi[- ]?flow|heated high", value):
        return "hfnc"
    if re.search(r"\bbipap\b|\bbi-pap\b|\bcpap\b|non[- ]?invasive|\bniv\b|autoset", value):
        return "niv"
    if re.search(r"tracheal suction|trach suction|trach collar|trach mask|suctioning", value):
        return ""
    if re.search(r"mechanical ventilation|invasive ventilation|\bintubat|endotracheal|\bett\b|tracheostomy|tracheotomy", value):
        return "imv"
    return ""


def first_nonnull_min(series: pd.Series) -> pd.Timestamp | float:
    series = series.dropna()
    return series.min() if len(series) else np.nan


def parse_eicu_age(value: object) -> float:
    text = str(value or "").strip()
    if not text:
        return np.nan
    if text.startswith(">"):
        nums = re.findall(r"\d+", text)
        return float(nums[0]) if nums else np.nan
    try:
        return float(text)
    except Exception:
        return np.nan


def load_mimic_core() -> pd.DataFrame:
    patients = pd.read_csv(
        MIMIC_ROOT / "hosp" / "patients.csv.gz",
        usecols=["subject_id", "gender", "anchor_age", "anchor_year_group"],
    )
    admissions = pd.read_csv(
        MIMIC_ROOT / "hosp" / "admissions.csv.gz",
        usecols=[
            "subject_id",
            "hadm_id",
            "admittime",
            "dischtime",
            "deathtime",
            "admission_type",
            "race",
            "hospital_expire_flag",
        ],
    )
    icu = pd.read_csv(
        MIMIC_ROOT / "icu" / "icustays.csv.gz",
        usecols=["subject_id", "hadm_id", "stay_id", "first_careunit", "last_careunit", "intime", "outtime", "los"],
    )
    for col in ["admittime", "dischtime", "deathtime"]:
        admissions[col] = pd.to_datetime(admissions[col], errors="coerce")
    for col in ["intime", "outtime"]:
        icu[col] = pd.to_datetime(icu[col], errors="coerce")
    core = icu.merge(patients, on="subject_id", how="left").merge(
        admissions.drop(columns=["subject_id"]),
        on="hadm_id",
        how="left",
    )
    core["adult"] = pd.to_numeric(core["anchor_age"], errors="coerce") >= 18
    return core


def extract_mimic_procedure_events(adult_stays: set[int]) -> tuple[pd.DataFrame, pd.DataFrame]:
    itemids = set(MIMIC_PROC_NIV) | set(MIMIC_PROC_IMV)
    path = MIMIC_ROOT / "icu" / "procedureevents.csv.gz"
    cols = ["stay_id", "itemid", "starttime", "endtime", "statusdescription"]
    events = pd.read_csv(path, usecols=cols)
    events = events[events["itemid"].isin(itemids) & events["stay_id"].isin(adult_stays)].copy()
    events["starttime"] = pd.to_datetime(events["starttime"], errors="coerce")
    events["endtime"] = pd.to_datetime(events["endtime"], errors="coerce")
    events = events.dropna(subset=["starttime"])
    niv = (
        events[events["itemid"].isin(MIMIC_PROC_NIV)]
        .groupby("stay_id", as_index=False)["starttime"]
        .min()
        .rename(columns={"starttime": "mimic_proc_niv_t0"})
    )
    imv = (
        events[events["itemid"].isin(MIMIC_PROC_IMV)]
        .groupby("stay_id", as_index=False)["starttime"]
        .min()
        .rename(columns={"starttime": "mimic_proc_imv_time"})
    )
    return niv, imv


def scan_mimic_o2_device(adult_stays: set[int]) -> tuple[pd.DataFrame, list[dict[str, object]], list[dict[str, object]]]:
    path = MIMIC_ROOT / "icu" / "chartevents.csv.gz"
    usecols = ["stay_id", "charttime", "itemid", "value"]
    t0_by_stay: dict[int, dict[str, object]] = {}
    value_counts: Counter[str] = Counter()
    support_counts: Counter[str] = Counter()
    scanned_rows = 0
    matched_rows = 0
    for chunk in pd.read_csv(path, usecols=usecols, chunksize=CHUNK_ROWS):
        scanned_rows += len(chunk)
        chunk = chunk[(chunk["itemid"] == MIMIC_O2_DEVICE_ITEM) & chunk["stay_id"].isin(adult_stays)]
        if chunk.empty:
            continue
        matched_rows += len(chunk)
        chunk["charttime"] = pd.to_datetime(chunk["charttime"], errors="coerce")
        chunk["support_class"] = chunk["value"].map(classify_support_text)
        for value, count in chunk["value"].fillna("").astype(str).value_counts().items():
            value_counts[value] += int(count)
        for klass, count in chunk["support_class"].value_counts().items():
            if klass:
                support_counts[klass] += int(count)
        elig = chunk[chunk["support_class"].isin(["hfnc", "niv"])].dropna(subset=["charttime"])
        for row in elig.itertuples(index=False):
            stay_id = int(row.stay_id)
            current = t0_by_stay.get(stay_id)
            if current is None or row.charttime < current["mimic_o2_device_t0"]:
                t0_by_stay[stay_id] = {
                    "stay_id": stay_id,
                    "mimic_o2_device_t0": row.charttime,
                    "mimic_o2_device_support": row.support_class,
                    "mimic_o2_device_value": row.value,
                }
    device_rows = list(t0_by_stay.values())
    device_t0 = pd.DataFrame(device_rows)
    if device_t0.empty:
        device_t0 = pd.DataFrame(columns=["stay_id", "mimic_o2_device_t0", "mimic_o2_device_support", "mimic_o2_device_value"])
    top_values = [
        {
            "dataset": "MIMIC-IV",
            "source": "chartevents.itemid_226732",
            "value": value,
            "n_events": count,
            "support_class": classify_support_text(value) or "other_or_unclear",
        }
        for value, count in value_counts.most_common(80)
    ]
    scan_rows = [
        {"dataset": "MIMIC-IV", "metric": "chartevents_rows_scanned_for_o2_device", "value": scanned_rows, "note": "all rows streamed in chunks"},
        {"dataset": "MIMIC-IV", "metric": "o2_device_rows_in_adult_stays", "value": matched_rows, "note": "itemid 226732 among adult ICU stays"},
    ]
    for support, count in support_counts.items():
        scan_rows.append({"dataset": "MIMIC-IV", "metric": f"o2_device_support_events::{support}", "value": count, "note": "classified by device value text"})
    return device_t0, top_values, scan_rows


def build_mimic_t0(core: pd.DataFrame, proc_niv: pd.DataFrame, device_t0: pd.DataFrame) -> pd.DataFrame:
    base = core[core["adult"]].copy()
    t0 = base[["subject_id", "hadm_id", "stay_id", "intime", "outtime", "deathtime", "hospital_expire_flag", "anchor_age", "gender", "race", "first_careunit", "los"]]
    t0 = t0.merge(proc_niv, on="stay_id", how="left").merge(device_t0, on="stay_id", how="left")

    def choose_t0(row: pd.Series) -> pd.Series:
        candidates: list[tuple[pd.Timestamp, str, str]] = []
        if pd.notna(row.get("mimic_proc_niv_t0")):
            candidates.append((row["mimic_proc_niv_t0"], "procedureevents", "niv"))
        if pd.notna(row.get("mimic_o2_device_t0")):
            candidates.append((row["mimic_o2_device_t0"], "chartevents_o2_device", str(row.get("mimic_o2_device_support") or "")))
        if not candidates:
            return pd.Series({"t0": pd.NaT, "t0_source": "", "t0_support_class": ""})
        candidates.sort(key=lambda x: x[0])
        return pd.Series({"t0": candidates[0][0], "t0_source": candidates[0][1], "t0_support_class": candidates[0][2]})

    t0 = pd.concat([t0, t0.apply(choose_t0, axis=1)], axis=1)
    t0 = t0.dropna(subset=["t0"])
    t0 = t0[(t0["t0"] >= t0["intime"]) & (t0["t0"] <= t0["outtime"])]
    t0 = t0.sort_values(["subject_id", "t0", "stay_id"]).groupby("subject_id", as_index=False).head(1)
    return t0.reset_index(drop=True)


def add_mimic_outcomes(t0: pd.DataFrame, proc_imv: pd.DataFrame) -> pd.DataFrame:
    out = t0.merge(proc_imv, on="stay_id", how="left")
    out["t0_plus_2h"] = out["t0"] + pd.Timedelta(hours=2)
    out["t0_plus_6h"] = out["t0"] + pd.Timedelta(hours=6)
    out["t0_plus_24h"] = out["t0"] + pd.Timedelta(hours=24)
    out["t0_plus_48h"] = out["t0"] + pd.Timedelta(hours=48)
    out["t0_plus_72h"] = out["t0"] + pd.Timedelta(hours=72)
    out["imv_before_or_at_t0"] = pd.notna(out["mimic_proc_imv_time"]) & (out["mimic_proc_imv_time"] <= out["t0"])
    out["imv_by_l2"] = pd.notna(out["mimic_proc_imv_time"]) & (out["mimic_proc_imv_time"] > out["t0"]) & (out["mimic_proc_imv_time"] <= out["t0_plus_2h"])
    out["imv_by_l6"] = pd.notna(out["mimic_proc_imv_time"]) & (out["mimic_proc_imv_time"] > out["t0"]) & (out["mimic_proc_imv_time"] <= out["t0_plus_6h"])
    out["death_by_l2"] = pd.notna(out["deathtime"]) & (out["deathtime"] > out["t0"]) & (out["deathtime"] <= out["t0_plus_2h"])
    out["death_by_l6"] = pd.notna(out["deathtime"]) & (out["deathtime"] > out["t0"]) & (out["deathtime"] <= out["t0_plus_6h"])
    out["l2_eligible"] = ~(out["imv_before_or_at_t0"] | out["imv_by_l2"] | out["death_by_l2"])
    out["l6_eligible"] = ~(out["imv_before_or_at_t0"] | out["imv_by_l6"] | out["death_by_l6"])
    out["failure_48h_l2"] = out["l2_eligible"] & (
        (pd.notna(out["mimic_proc_imv_time"]) & (out["mimic_proc_imv_time"] > out["t0_plus_2h"]) & (out["mimic_proc_imv_time"] <= out["t0_plus_48h"]))
        | (pd.notna(out["deathtime"]) & (out["deathtime"] > out["t0_plus_2h"]) & (out["deathtime"] <= out["t0_plus_48h"]))
    )
    out["failure_24h_l2"] = out["l2_eligible"] & (
        (pd.notna(out["mimic_proc_imv_time"]) & (out["mimic_proc_imv_time"] > out["t0_plus_2h"]) & (out["mimic_proc_imv_time"] <= out["t0_plus_24h"]))
        | (pd.notna(out["deathtime"]) & (out["deathtime"] > out["t0_plus_2h"]) & (out["deathtime"] <= out["t0_plus_24h"]))
    )
    out["failure_72h_l2"] = out["l2_eligible"] & (
        (pd.notna(out["mimic_proc_imv_time"]) & (out["mimic_proc_imv_time"] > out["t0_plus_2h"]) & (out["mimic_proc_imv_time"] <= out["t0_plus_72h"]))
        | (pd.notna(out["deathtime"]) & (out["deathtime"] > out["t0_plus_2h"]) & (out["deathtime"] <= out["t0_plus_72h"]))
    )
    out["failure_48h_l6"] = out["l6_eligible"] & (
        (pd.notna(out["mimic_proc_imv_time"]) & (out["mimic_proc_imv_time"] > out["t0_plus_6h"]) & (out["mimic_proc_imv_time"] <= out["t0_plus_48h"]))
        | (pd.notna(out["deathtime"]) & (out["deathtime"] > out["t0_plus_6h"]) & (out["deathtime"] <= out["t0_plus_48h"]))
    )
    return out


def scan_mimic_physiology(t0: pd.DataFrame) -> tuple[pd.DataFrame, list[dict[str, object]]]:
    if t0.empty:
        return t0, []
    item_to_var = MIMIC_PHYS_ITEM_MAP
    itemids = set(item_to_var)
    t0_index = t0.set_index("stay_id")["t0"]
    candidate_stays = set(int(x) for x in t0["stay_id"])
    var_bits = {var: 1 << i for i, var in enumerate(sorted(set(item_to_var.values())))}
    mask_by_stay: dict[int, int] = defaultdict(int)
    events_by_var: Counter[str] = Counter()
    path = MIMIC_ROOT / "icu" / "chartevents.csv.gz"
    usecols = ["stay_id", "charttime", "itemid", "valuenum", "value"]
    scanned_rows = 0
    matched_rows = 0
    for chunk in pd.read_csv(path, usecols=usecols, chunksize=CHUNK_ROWS):
        scanned_rows += len(chunk)
        chunk = chunk[chunk["itemid"].isin(itemids) & chunk["stay_id"].isin(candidate_stays)]
        if chunk.empty:
            continue
        matched_rows += len(chunk)
        chunk["charttime"] = pd.to_datetime(chunk["charttime"], errors="coerce")
        chunk["t0"] = chunk["stay_id"].map(t0_index)
        window = (chunk["charttime"] >= chunk["t0"] - pd.Timedelta(hours=6)) & (chunk["charttime"] <= chunk["t0"] + pd.Timedelta(hours=2))
        chunk = chunk[window]
        if chunk.empty:
            continue
        chunk["var"] = chunk["itemid"].map(item_to_var)
        chunk = chunk[pd.notna(chunk["var"])]
        for row in chunk.itertuples(index=False):
            stay_id = int(row.stay_id)
            var = str(row.var)
            mask_by_stay[stay_id] |= var_bits[var]
            events_by_var[var] += 1
    out = t0.copy()
    out["phys_mask_l2"] = out["stay_id"].map(lambda x: mask_by_stay.get(int(x), 0))
    for var, bit in var_bits.items():
        out[f"has_{var}_l2_window"] = (out["phys_mask_l2"] & bit) > 0
    rows = [
        {"dataset": "MIMIC-IV", "metric": "chartevents_rows_scanned_for_physiology", "value": scanned_rows, "note": "all rows streamed in chunks"},
        {"dataset": "MIMIC-IV", "metric": "candidate_physiology_rows_by_itemid", "value": matched_rows, "note": "candidate itemids among T0 cohort stays before time-window filtering"},
    ]
    for var, count in events_by_var.items():
        rows.append({"dataset": "MIMIC-IV", "metric": f"baseline_dynamic_events::{var}", "value": count, "note": "T0-6h through L2 window"})
    return out, rows


def summarize_bool_count(df: pd.DataFrame, col: str) -> int:
    if col not in df:
        return 0
    return int(df[col].fillna(False).astype(bool).sum())


def mimic_counts(core: pd.DataFrame, t0: pd.DataFrame, extra_rows: list[dict[str, object]]) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = [
        {"dataset": "MIMIC-IV", "metric": "icu_stays_total", "value": len(core), "note": "icu/icustays rows"},
        {"dataset": "MIMIC-IV", "metric": "adult_icu_stays", "value": int(core["adult"].sum()), "note": "anchor_age >= 18"},
        {"dataset": "MIMIC-IV", "metric": "first_eligible_hfnc_or_niv_stays", "value": len(t0), "note": "first eligible support stay per subject_id"},
    ]
    for source, count in t0["t0_source"].value_counts(dropna=False).items():
        rows.append({"dataset": "MIMIC-IV", "metric": f"t0_source::{source}", "value": int(count), "note": "earliest support source"})
    for klass, count in t0["t0_support_class"].value_counts(dropna=False).items():
        rows.append({"dataset": "MIMIC-IV", "metric": f"t0_support_class::{klass}", "value": int(count), "note": "earliest support class"})
    rows.extend(
        [
            {"dataset": "MIMIC-IV", "metric": "excluded_imv_before_or_at_t0", "value": summarize_bool_count(t0, "imv_before_or_at_t0"), "note": "likely already invasive-ventilated before eligible support start"},
            {"dataset": "MIMIC-IV", "metric": "immediate_failure_by_l2", "value": summarize_bool_count(t0, "imv_by_l2") + summarize_bool_count(t0, "death_by_l2"), "note": "IMV or death from T0 to L2"},
            {"dataset": "MIMIC-IV", "metric": "l2_eligible", "value": summarize_bool_count(t0, "l2_eligible"), "note": "primary landmark denominator"},
            {"dataset": "MIMIC-IV", "metric": "failure_24h_l2", "value": summarize_bool_count(t0, "failure_24h_l2"), "note": "IMV or death from L2 to T0+24h"},
            {"dataset": "MIMIC-IV", "metric": "failure_48h_l2", "value": summarize_bool_count(t0, "failure_48h_l2"), "note": "primary event count: IMV or death from L2 to T0+48h"},
            {"dataset": "MIMIC-IV", "metric": "failure_72h_l2", "value": summarize_bool_count(t0, "failure_72h_l2"), "note": "sensitivity event count"},
            {"dataset": "MIMIC-IV", "metric": "l6_eligible", "value": summarize_bool_count(t0, "l6_eligible"), "note": "L6 sensitivity denominator"},
            {"dataset": "MIMIC-IV", "metric": "failure_48h_l6", "value": summarize_bool_count(t0, "failure_48h_l6"), "note": "L6 sensitivity event count"},
        ]
    )
    if "phys_mask_l2" in t0:
        for var in sorted(set(MIMIC_PHYS_ITEM_MAP.values())):
            col = f"has_{var}_l2_window"
            rows.append({"dataset": "MIMIC-IV", "metric": f"l2_window_available::{var}", "value": summarize_bool_count(t0, col), "note": "at least one value in T0-6h to L2 window"})
        for name, required in [("required_rox_like", MIMIC_REQUIRED_PHYS), ("extended_hacor_like", MIMIC_EXTENDED_PHYS)]:
            cols = [f"has_{var}_l2_window" for var in required if f"has_{var}_l2_window" in t0]
            denom = t0[t0["l2_eligible"]].copy()
            complete = int(denom[cols].all(axis=1).sum()) if cols and len(denom) else 0
            rows.append({"dataset": "MIMIC-IV", "metric": f"l2_complete::{name}", "value": complete, "note": ";".join(required)})
    rows.extend(extra_rows)
    return rows


def load_eicu_patient() -> pd.DataFrame:
    patient = pd.read_csv(EICU_ROOT / "patient.csv.gz", dtype=str)
    patient["patientunitstayid"] = pd.to_numeric(patient["patientunitstayid"], errors="coerce")
    patient = patient[pd.notna(patient["patientunitstayid"])].copy()
    patient["patientunitstayid"] = patient["patientunitstayid"].astype("int64")
    patient["age_num"] = patient["age"].map(parse_eicu_age)
    for col in ["unitdischargeoffset", "hospitaldischargeoffset", "hospitaladmitoffset"]:
        patient[col] = pd.to_numeric(patient[col], errors="coerce")
    patient["adult"] = patient["age_num"] >= 18
    patient["hospital_death"] = patient["hospitaldischargestatus"].str.lower().eq("expired")
    patient["icu_death"] = patient["unitdischargestatus"].str.lower().eq("expired")
    return patient


def normalize_eicu_unit_id(df: pd.DataFrame) -> pd.DataFrame:
    if "patientunitstayid" not in df.columns:
        return df
    out = df.copy()
    out["patientunitstayid"] = pd.to_numeric(out["patientunitstayid"], errors="coerce")
    out = out[pd.notna(out["patientunitstayid"])].copy()
    out["patientunitstayid"] = out["patientunitstayid"].astype("int64")
    return out


def scan_eicu_treatment(adult_units: set[int]) -> tuple[pd.DataFrame, pd.DataFrame, list[dict[str, object]]]:
    path = EICU_ROOT / "treatment.csv.gz"
    rows = []
    label_counts: Counter[str] = Counter()
    usecols = ["patientunitstayid", "treatmentoffset", "treatmentstring"]
    for chunk in pd.read_csv(path, usecols=usecols, dtype=str, chunksize=CHUNK_ROWS):
        chunk["patientunitstayid"] = pd.to_numeric(chunk["patientunitstayid"], errors="coerce")
        chunk = chunk[chunk["patientunitstayid"].isin(adult_units)].copy()
        if chunk.empty:
            continue
        for value, count in chunk["treatmentstring"].fillna("").value_counts().items():
            if re.search(r"ventilat|cpap|bipap|peep|oxygen|high[- ]?flow|hfnc", value, flags=re.I):
                label_counts[value] += int(count)
        chunk["offset"] = pd.to_numeric(chunk["treatmentoffset"], errors="coerce")
        chunk["support_class"] = chunk["treatmentstring"].map(classify_support_text)
        elig = chunk[chunk["support_class"].isin(["hfnc", "niv", "imv"]) & pd.notna(chunk["offset"])]
        if not elig.empty:
            rows.append(elig[["patientunitstayid", "offset", "support_class", "treatmentstring"]])
    if rows:
        events = pd.concat(rows, ignore_index=True)
    else:
        events = pd.DataFrame(columns=["patientunitstayid", "offset", "support_class", "treatmentstring"])
    t0 = (
        events[events["support_class"].isin(["hfnc", "niv"])]
        .sort_values(["patientunitstayid", "offset"])
        .groupby("patientunitstayid", as_index=False)
        .first()
        .rename(columns={"offset": "eicu_treatment_t0", "support_class": "eicu_treatment_support", "treatmentstring": "eicu_treatment_t0_label"})
    )
    imv = (
        events[events["support_class"].eq("imv")]
        .sort_values(["patientunitstayid", "offset"])
        .groupby("patientunitstayid", as_index=False)
        .first()
        .rename(columns={"offset": "eicu_treatment_imv_time", "treatmentstring": "eicu_treatment_imv_label"})
    )
    top = [
        {
            "dataset": "eICU-CRD",
            "source": "treatment.treatmentstring",
            "value": value,
            "n_events": count,
            "support_class": classify_support_text(value) or "other_or_unclear",
        }
        for value, count in label_counts.most_common(100)
    ]
    return t0, imv[["patientunitstayid", "eicu_treatment_imv_time", "eicu_treatment_imv_label"]], top


def scan_eicu_respiratory_care(adult_units: set[int]) -> pd.DataFrame:
    path = EICU_ROOT / "respiratoryCare.csv.gz"
    cols = ["patientunitstayid", "respcarestatusoffset", "airwaytype", "ventstartoffset", "ventendoffset", "setapneafio2"]
    df = pd.read_csv(path, usecols=cols, dtype=str)
    df["patientunitstayid"] = pd.to_numeric(df["patientunitstayid"], errors="coerce")
    df = df[df["patientunitstayid"].isin(adult_units)].copy()
    df["offset"] = pd.to_numeric(df["ventstartoffset"], errors="coerce")
    fallback = pd.to_numeric(df["respcarestatusoffset"], errors="coerce")
    df["offset"] = df["offset"].where(pd.notna(df["offset"]) & (df["offset"] != 0), fallback)
    df["airway_lower"] = df["airwaytype"].fillna("").str.lower()
    df = df[df["airway_lower"].str.contains(r"ett|tracheostomy|artificial airway", regex=True, na=False)]
    if df.empty:
        return pd.DataFrame(columns=["patientunitstayid", "eicu_respcare_imv_time", "eicu_respcare_airwaytype"])
    out = (
        df[pd.notna(df["offset"])]
        .sort_values(["patientunitstayid", "offset"])
        .groupby("patientunitstayid", as_index=False)
        .first()
        .rename(columns={"offset": "eicu_respcare_imv_time", "airwaytype": "eicu_respcare_airwaytype"})
    )
    return out[["patientunitstayid", "eicu_respcare_imv_time", "eicu_respcare_airwaytype"]]


def scan_eicu_respchart(adult_units: set[int]) -> tuple[pd.DataFrame, pd.DataFrame, list[dict[str, object]]]:
    path = EICU_ROOT / "respiratoryCharting.csv.gz"
    cols = ["patientunitstayid", "respchartoffset", "respchartvaluelabel", "respchartvalue"]
    label_counts: Counter[str] = Counter()
    value_counts: Counter[str] = Counter()
    support_rows = []
    usecols = cols
    for chunk in pd.read_csv(path, usecols=usecols, dtype=str, chunksize=CHUNK_ROWS):
        chunk["patientunitstayid"] = pd.to_numeric(chunk["patientunitstayid"], errors="coerce")
        chunk = chunk[chunk["patientunitstayid"].isin(adult_units)].copy()
        if chunk.empty:
            continue
        for value, count in chunk["respchartvaluelabel"].fillna("").value_counts().items():
            label_counts[value] += int(count)
        for value, count in chunk["respchartvalue"].fillna("").value_counts().items():
            if re.search(r"vent|cpap|bipap|niv|high[- ]?flow|hfnc|continued|on|off", value, flags=re.I):
                value_counts[value] += int(count)
        joined = chunk["respchartvaluelabel"].fillna("") + " " + chunk["respchartvalue"].fillna("")
        chunk["support_class"] = joined.map(classify_support_text)
        chunk["offset"] = pd.to_numeric(chunk["respchartoffset"], errors="coerce")
        elig = chunk[chunk["support_class"].isin(["hfnc", "niv", "imv"]) & pd.notna(chunk["offset"])]
        if not elig.empty:
            support_rows.append(elig[["patientunitstayid", "offset", "support_class", "respchartvaluelabel", "respchartvalue"]])
    events = pd.concat(support_rows, ignore_index=True) if support_rows else pd.DataFrame(columns=["patientunitstayid", "offset", "support_class", "respchartvaluelabel", "respchartvalue"])
    t0 = (
        events[events["support_class"].isin(["hfnc", "niv"])]
        .sort_values(["patientunitstayid", "offset"])
        .groupby("patientunitstayid", as_index=False)
        .first()
        .rename(columns={"offset": "eicu_respchart_t0", "support_class": "eicu_respchart_support"})
    )
    imv = (
        events[events["support_class"].eq("imv")]
        .sort_values(["patientunitstayid", "offset"])
        .groupby("patientunitstayid", as_index=False)
        .first()
        .rename(columns={"offset": "eicu_respchart_imv_time"})
    )
    labels = [
        {
            "dataset": "eICU-CRD",
            "source": "respiratoryCharting.respchartvaluelabel",
            "value": value,
            "n_events": count,
            "support_class": classify_support_text(value) or "measurement_or_unclear",
        }
        for value, count in label_counts.most_common(120)
    ]
    labels.extend(
        {
            "dataset": "eICU-CRD",
            "source": "respiratoryCharting.respchartvalue_keyword_subset",
            "value": value,
            "n_events": count,
            "support_class": classify_support_text(value) or "other_or_unclear",
        }
        for value, count in value_counts.most_common(60)
    )
    return t0, imv[["patientunitstayid", "eicu_respchart_imv_time"]], labels


def build_eicu_t0(
    patient: pd.DataFrame,
    treatment_t0: pd.DataFrame,
    respchart_t0: pd.DataFrame,
    treatment_imv: pd.DataFrame,
    respcare_imv: pd.DataFrame,
    respchart_imv: pd.DataFrame,
) -> pd.DataFrame:
    base = patient[patient["adult"]].copy()
    keep = [
        "patientunitstayid",
        "uniquepid",
        "age",
        "age_num",
        "gender",
        "ethnicity",
        "hospitalid",
        "unittype",
        "unitstaytype",
        "unitvisitnumber",
        "hospitaldischargestatus",
        "unitdischargestatus",
        "hospitaldischargeoffset",
        "unitdischargeoffset",
        "hospital_death",
        "icu_death",
    ]
    base = base[keep]
    base = normalize_eicu_unit_id(base)
    treatment_t0 = normalize_eicu_unit_id(treatment_t0)
    respchart_t0 = normalize_eicu_unit_id(respchart_t0)
    treatment_imv = normalize_eicu_unit_id(treatment_imv)
    respcare_imv = normalize_eicu_unit_id(respcare_imv)
    respchart_imv = normalize_eicu_unit_id(respchart_imv)
    out = base.merge(treatment_t0, on="patientunitstayid", how="left").merge(respchart_t0, on="patientunitstayid", how="left")
    out = out.merge(treatment_imv, on="patientunitstayid", how="left").merge(respcare_imv, on="patientunitstayid", how="left").merge(respchart_imv, on="patientunitstayid", how="left")

    def choose_t0(row: pd.Series) -> pd.Series:
        candidates: list[tuple[float, str, str]] = []
        if pd.notna(row.get("eicu_treatment_t0")):
            candidates.append((float(row["eicu_treatment_t0"]), "treatment", str(row.get("eicu_treatment_support") or "")))
        if pd.notna(row.get("eicu_respchart_t0")):
            candidates.append((float(row["eicu_respchart_t0"]), "respiratoryCharting", str(row.get("eicu_respchart_support") or "")))
        if not candidates:
            return pd.Series({"t0_offset": np.nan, "t0_source": "", "t0_support_class": ""})
        candidates.sort(key=lambda x: x[0])
        return pd.Series({"t0_offset": candidates[0][0], "t0_source": candidates[0][1], "t0_support_class": candidates[0][2]})

    out = pd.concat([out, out.apply(choose_t0, axis=1)], axis=1)
    out = out[pd.notna(out["t0_offset"])].copy()
    # eICU unit encounters are the native unit of analysis. For V21 feasibility,
    # retain all eligible patient-unit stays and document this boundary.
    return out.reset_index(drop=True)


def add_eicu_outcomes(t0: pd.DataFrame) -> pd.DataFrame:
    out = t0.copy()
    imv_cols = ["eicu_treatment_imv_time", "eicu_respcare_imv_time", "eicu_respchart_imv_time"]
    for col in imv_cols:
        out[col] = pd.to_numeric(out[col], errors="coerce")
    out["first_imv_time"] = out[imv_cols].min(axis=1, skipna=True)
    out["t0_plus_2h"] = out["t0_offset"] + 120
    out["t0_plus_6h"] = out["t0_offset"] + 360
    out["t0_plus_24h"] = out["t0_offset"] + 1440
    out["t0_plus_48h"] = out["t0_offset"] + 2880
    out["t0_plus_72h"] = out["t0_offset"] + 4320
    out["death_time"] = np.where(out["hospital_death"], out["hospitaldischargeoffset"], np.nan)
    out["imv_before_or_at_t0"] = pd.notna(out["first_imv_time"]) & (out["first_imv_time"] <= out["t0_offset"])
    out["imv_by_l2"] = pd.notna(out["first_imv_time"]) & (out["first_imv_time"] > out["t0_offset"]) & (out["first_imv_time"] <= out["t0_plus_2h"])
    out["imv_by_l6"] = pd.notna(out["first_imv_time"]) & (out["first_imv_time"] > out["t0_offset"]) & (out["first_imv_time"] <= out["t0_plus_6h"])
    out["death_by_l2"] = pd.notna(out["death_time"]) & (out["death_time"] > out["t0_offset"]) & (out["death_time"] <= out["t0_plus_2h"])
    out["death_by_l6"] = pd.notna(out["death_time"]) & (out["death_time"] > out["t0_offset"]) & (out["death_time"] <= out["t0_plus_6h"])
    out["l2_eligible"] = ~(out["imv_before_or_at_t0"] | out["imv_by_l2"] | out["death_by_l2"])
    out["l6_eligible"] = ~(out["imv_before_or_at_t0"] | out["imv_by_l6"] | out["death_by_l6"])
    out["failure_48h_l2"] = out["l2_eligible"] & (
        (pd.notna(out["first_imv_time"]) & (out["first_imv_time"] > out["t0_plus_2h"]) & (out["first_imv_time"] <= out["t0_plus_48h"]))
        | (pd.notna(out["death_time"]) & (out["death_time"] > out["t0_plus_2h"]) & (out["death_time"] <= out["t0_plus_48h"]))
    )
    out["failure_24h_l2"] = out["l2_eligible"] & (
        (pd.notna(out["first_imv_time"]) & (out["first_imv_time"] > out["t0_plus_2h"]) & (out["first_imv_time"] <= out["t0_plus_24h"]))
        | (pd.notna(out["death_time"]) & (out["death_time"] > out["t0_plus_2h"]) & (out["death_time"] <= out["t0_plus_24h"]))
    )
    out["failure_72h_l2"] = out["l2_eligible"] & (
        (pd.notna(out["first_imv_time"]) & (out["first_imv_time"] > out["t0_plus_2h"]) & (out["first_imv_time"] <= out["t0_plus_72h"]))
        | (pd.notna(out["death_time"]) & (out["death_time"] > out["t0_plus_2h"]) & (out["death_time"] <= out["t0_plus_72h"]))
    )
    out["failure_48h_l6"] = out["l6_eligible"] & (
        (pd.notna(out["first_imv_time"]) & (out["first_imv_time"] > out["t0_plus_6h"]) & (out["first_imv_time"] <= out["t0_plus_48h"]))
        | (pd.notna(out["death_time"]) & (out["death_time"] > out["t0_plus_6h"]) & (out["death_time"] <= out["t0_plus_48h"]))
    )
    return out


def scan_eicu_physiology(t0: pd.DataFrame) -> tuple[pd.DataFrame, list[dict[str, object]]]:
    if t0.empty:
        return t0, []
    candidate_units = set(int(x) for x in t0["patientunitstayid"])
    t0_index = t0.set_index("patientunitstayid")["t0_offset"]
    var_bits = {var: 1 << i for i, var in enumerate(sorted(set(EICU_EXTENDED_PHYS)))}
    mask_by_unit: dict[int, int] = defaultdict(int)
    events_by_var: Counter[str] = Counter()

    vital_path = EICU_ROOT / "vitalPeriodic.csv.gz"
    vital_cols = ["patientunitstayid", "observationoffset", "sao2", "heartrate", "respiration"]
    scanned = 0
    matched = 0
    for chunk in pd.read_csv(vital_path, usecols=vital_cols, chunksize=CHUNK_ROWS):
        scanned += len(chunk)
        chunk = chunk[chunk["patientunitstayid"].isin(candidate_units)].copy()
        if chunk.empty:
            continue
        matched += len(chunk)
        chunk["t0"] = chunk["patientunitstayid"].map(t0_index)
        chunk["observationoffset"] = pd.to_numeric(chunk["observationoffset"], errors="coerce")
        chunk = chunk[(chunk["observationoffset"] >= chunk["t0"] - 360) & (chunk["observationoffset"] <= chunk["t0"] + 120)]
        for source_col, var in [("sao2", "sao2"), ("heartrate", "heartrate"), ("respiration", "respiration")]:
            ok = pd.to_numeric(chunk[source_col], errors="coerce").notna()
            for unit in chunk.loc[ok, "patientunitstayid"]:
                mask_by_unit[int(unit)] |= var_bits[var]
                events_by_var[var] += 1

    chart_path = EICU_ROOT / "respiratoryCharting.csv.gz"
    chart_cols = ["patientunitstayid", "respchartoffset", "respchartvaluelabel", "respchartvalue"]
    for chunk in pd.read_csv(chart_path, usecols=chart_cols, dtype=str, chunksize=CHUNK_ROWS):
        chunk["patientunitstayid"] = pd.to_numeric(chunk["patientunitstayid"], errors="coerce")
        chunk = chunk[chunk["patientunitstayid"].isin(candidate_units)].copy()
        if chunk.empty:
            continue
        chunk["offset"] = pd.to_numeric(chunk["respchartoffset"], errors="coerce")
        chunk["t0"] = chunk["patientunitstayid"].map(t0_index)
        chunk = chunk[(chunk["offset"] >= chunk["t0"] - 360) & (chunk["offset"] <= chunk["t0"] + 120)]
        if chunk.empty:
            continue
        labels = chunk["respchartvaluelabel"].fillna("").str.lower()
        mapping = [
            ("fio2", labels.str.contains(r"\bfio2\b", regex=True, na=False)),
            ("peep", labels.str.contains(r"\bpeep\b", regex=True, na=False)),
            ("lpm_o2", labels.str.contains(r"lpm o2|o2 flow|flow", regex=True, na=False)),
            ("vent_rate", labels.str.contains(r"vent rate|total rr|respiratory rate", regex=True, na=False)),
        ]
        for var, mask in mapping:
            valid = pd.to_numeric(chunk.loc[mask, "respchartvalue"], errors="coerce").notna()
            for unit in chunk.loc[mask].loc[valid, "patientunitstayid"]:
                mask_by_unit[int(unit)] |= var_bits[var]
                events_by_var[var] += 1

    out = t0.copy()
    out["phys_mask_l2"] = out["patientunitstayid"].map(lambda x: mask_by_unit.get(int(x), 0))
    for var, bit in var_bits.items():
        out[f"has_{var}_l2_window"] = (out["phys_mask_l2"] & bit) > 0
    rows = [
        {"dataset": "eICU-CRD", "metric": "vitalPeriodic_rows_scanned_for_physiology", "value": scanned, "note": "all rows streamed in chunks"},
        {"dataset": "eICU-CRD", "metric": "candidate_vitalPeriodic_rows", "value": matched, "note": "rows in candidate support unit stays before time-window filtering"},
    ]
    for var, count in events_by_var.items():
        rows.append({"dataset": "eICU-CRD", "metric": f"baseline_dynamic_events::{var}", "value": count, "note": "T0-6h through L2 window"})
    return out, rows


def eicu_counts(patient: pd.DataFrame, t0: pd.DataFrame, extra_rows: list[dict[str, object]]) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = [
        {"dataset": "eICU-CRD", "metric": "patient_unit_stays_total", "value": len(patient), "note": "patient.csv rows"},
        {"dataset": "eICU-CRD", "metric": "adult_patient_unit_stays", "value": int(patient["adult"].sum()), "note": "age >=18; age >89 treated as 89"},
        {"dataset": "eICU-CRD", "metric": "eligible_hfnc_or_niv_unit_stays", "value": len(t0), "note": "patient-unit stays; not yet deduplicated by uniquepid"},
    ]
    for source, count in t0["t0_source"].value_counts(dropna=False).items():
        rows.append({"dataset": "eICU-CRD", "metric": f"t0_source::{source}", "value": int(count), "note": "earliest support source"})
    for klass, count in t0["t0_support_class"].value_counts(dropna=False).items():
        rows.append({"dataset": "eICU-CRD", "metric": f"t0_support_class::{klass}", "value": int(count), "note": "earliest support class"})
    rows.extend(
        [
            {"dataset": "eICU-CRD", "metric": "excluded_imv_before_or_at_t0", "value": summarize_bool_count(t0, "imv_before_or_at_t0"), "note": "likely already invasive-ventilated before eligible support start"},
            {"dataset": "eICU-CRD", "metric": "immediate_failure_by_l2", "value": summarize_bool_count(t0, "imv_by_l2") + summarize_bool_count(t0, "death_by_l2"), "note": "IMV or death from T0 to L2"},
            {"dataset": "eICU-CRD", "metric": "l2_eligible", "value": summarize_bool_count(t0, "l2_eligible"), "note": "primary landmark denominator"},
            {"dataset": "eICU-CRD", "metric": "failure_24h_l2", "value": summarize_bool_count(t0, "failure_24h_l2"), "note": "IMV or death from L2 to T0+24h"},
            {"dataset": "eICU-CRD", "metric": "failure_48h_l2", "value": summarize_bool_count(t0, "failure_48h_l2"), "note": "primary event count: IMV or death from L2 to T0+48h"},
            {"dataset": "eICU-CRD", "metric": "failure_72h_l2", "value": summarize_bool_count(t0, "failure_72h_l2"), "note": "sensitivity event count"},
            {"dataset": "eICU-CRD", "metric": "l6_eligible", "value": summarize_bool_count(t0, "l6_eligible"), "note": "L6 sensitivity denominator"},
            {"dataset": "eICU-CRD", "metric": "failure_48h_l6", "value": summarize_bool_count(t0, "failure_48h_l6"), "note": "L6 sensitivity event count"},
        ]
    )
    if "phys_mask_l2" in t0:
        for var in sorted(set(EICU_EXTENDED_PHYS)):
            col = f"has_{var}_l2_window"
            rows.append({"dataset": "eICU-CRD", "metric": f"l2_window_available::{var}", "value": summarize_bool_count(t0, col), "note": "at least one value in T0-6h to L2 window"})
        for name, required in [("required_rox_like", EICU_REQUIRED_PHYS), ("extended_hacor_like", EICU_EXTENDED_PHYS)]:
            cols = [f"has_{var}_l2_window" for var in required if f"has_{var}_l2_window" in t0]
            denom = t0[t0["l2_eligible"]].copy()
            complete = int(denom[cols].all(axis=1).sum()) if cols and len(denom) else 0
            rows.append({"dataset": "eICU-CRD", "metric": f"l2_complete::{name}", "value": complete, "note": ";".join(required)})
    rows.extend(extra_rows)
    return rows


def build_variable_lock(mimic_counts_rows: list[dict[str, object]], eicu_counts_rows: list[dict[str, object]]) -> list[dict[str, object]]:
    rows = [
        {
            "dataset": "MIMIC-IV",
            "domain": "cohort_denominator",
            "source_table": "icu/icustays.csv.gz;hosp/patients.csv.gz;hosp/admissions.csv.gz",
            "variable_or_item": "stay_id;subject_id;anchor_age;intime;outtime",
            "role": "adult ICU denominator",
            "lock_status": "locked_for_feasibility",
            "notes": "Primary V21 unit is first eligible HFNC/NIV stay per subject, not all ICU stays.",
        },
        {
            "dataset": "MIMIC-IV",
            "domain": "time_zero",
            "source_table": "icu/procedureevents.csv.gz;icu/chartevents.csv.gz",
            "variable_or_item": "225794 Non-invasive Ventilation; 226732 O2 Delivery Device(s)",
            "role": "first HFNC/NIV support start",
            "lock_status": "provisionally_locked",
            "notes": "HFNC from O2 device text; NIV from procedureevents and device text. Requires clinical review of device values before final modeling.",
        },
        {
            "dataset": "MIMIC-IV",
            "domain": "primary_outcome",
            "source_table": "icu/procedureevents.csv.gz;hosp/admissions.csv.gz",
            "variable_or_item": "225792 Invasive Ventilation; 224385 Intubation; tracheostomy itemids; deathtime",
            "role": "IMV or death from L2 to T0+48h",
            "lock_status": "provisionally_locked",
            "notes": "Procedureevents gives clean invasive-support timing; death uses admission deathtime.",
        },
        {
            "dataset": "MIMIC-IV",
            "domain": "dynamic_predictors",
            "source_table": "icu/chartevents.csv.gz",
            "variable_or_item": "SpO2 220277; FiO2 223835; RR 220210/224690/224689; HR 220045; O2 flow; PEEP; BP; GCS",
            "role": "ROX/HACOR-like feature availability",
            "lock_status": "feasibility_scanned",
            "notes": "V21 reports window-level availability, not cleaned feature values.",
        },
        {
            "dataset": "eICU-CRD",
            "domain": "cohort_denominator",
            "source_table": "patient.csv.gz",
            "variable_or_item": "patientunitstayid;uniquepid;age;unit offsets; discharge status",
            "role": "adult patient-unit denominator",
            "lock_status": "locked_for_feasibility",
            "notes": "V21 retains patient-unit stays; patient-level deduplication remains a modeling-stage sensitivity.",
        },
        {
            "dataset": "eICU-CRD",
            "domain": "time_zero",
            "source_table": "treatment.csv.gz;respiratoryCharting.csv.gz",
            "variable_or_item": "treatmentstring CPAP/BiPAP/NIV/HFNC keywords; respiratory chart labels/values",
            "role": "candidate first HFNC/NIV support start",
            "lock_status": "mapping_feasibility_only",
            "notes": "HFNC availability is uncertain; eICU may support broader NIV/ventilation external validation rather than exact HFNC/NIV replication.",
        },
        {
            "dataset": "eICU-CRD",
            "domain": "primary_outcome",
            "source_table": "treatment.csv.gz;respiratoryCare.csv.gz;patient.csv.gz",
            "variable_or_item": "mechanical ventilation treatment; airwaytype ETT/tracheostomy; hospital death offset",
            "role": "IMV or death from L2 to T0+48h",
            "lock_status": "mapping_feasibility_only",
            "notes": "Outcome mapping is broader than MIMIC and must be reviewed before external validation claims.",
        },
        {
            "dataset": "eICU-CRD",
            "domain": "dynamic_predictors",
            "source_table": "vitalPeriodic.csv.gz;respiratoryCharting.csv.gz;apacheApsVar.csv.gz",
            "variable_or_item": "SaO2; respiration; heart rate; FiO2; PEEP; LPM O2; vent rate",
            "role": "ROX/HACOR-like feature availability",
            "lock_status": "feasibility_scanned",
            "notes": "V21 reports window-level availability only.",
        },
    ]
    # Add compact go/no-go metrics into notes-free aggregate form.
    for dataset, counts in [("MIMIC-IV", mimic_counts_rows), ("eICU-CRD", eicu_counts_rows)]:
        metric_map = {str(row["metric"]): row["value"] for row in counts if row.get("dataset") == dataset}
        rows.append(
            {
                "dataset": dataset,
                "domain": "go_no_go_snapshot",
                "source_table": "V21 aggregate feasibility outputs",
                "variable_or_item": "l2_eligible;failure_48h_l2",
                "role": "minimum sample/event gate",
                "lock_status": "pass_if_counts_meet_thresholds" if int(metric_map.get("l2_eligible", 0)) > 0 else "needs_rebuild",
                "notes": f"L2 eligible={metric_map.get('l2_eligible', 0)}; 48h events={metric_map.get('failure_48h_l2', 0)}.",
            }
        )
    return rows


def metric_value(rows: list[dict[str, object]], dataset: str, metric: str) -> int:
    for row in rows:
        if row.get("dataset") == dataset and row.get("metric") == metric:
            try:
                return int(float(str(row.get("value"))))
            except Exception:
                return 0
    return 0


def build_go_no_go(mimic_rows: list[dict[str, object]], eicu_rows: list[dict[str, object]]) -> list[dict[str, object]]:
    mimic_n = metric_value(mimic_rows, "MIMIC-IV", "l2_eligible")
    mimic_events = metric_value(mimic_rows, "MIMIC-IV", "failure_48h_l2")
    mimic_rox = metric_value(mimic_rows, "MIMIC-IV", "l2_complete::required_rox_like")
    eicu_n = metric_value(eicu_rows, "eICU-CRD", "l2_eligible")
    eicu_events = metric_value(eicu_rows, "eICU-CRD", "failure_48h_l2")
    eicu_rox = metric_value(eicu_rows, "eICU-CRD", "l2_complete::required_rox_like")
    return [
        {
            "criterion": "MIMIC primary cohort size",
            "threshold": "L2 eligible >= 1000",
            "observed": mimic_n,
            "status": "pass" if mimic_n >= 1000 else "fail_or_borderline",
            "interpretation": "Needed for internal derivation feasibility.",
        },
        {
            "criterion": "MIMIC primary events",
            "threshold": "48h failures >= 200",
            "observed": mimic_events,
            "status": "pass" if mimic_events >= 200 else "fail_or_borderline",
            "interpretation": "Needed before model development.",
        },
        {
            "criterion": "MIMIC ROX-like availability",
            "threshold": "complete rows >= 500",
            "observed": mimic_rox,
            "status": "pass" if mimic_rox >= 500 else "needs_feature_review",
            "interpretation": "SpO2/FiO2/RR/HR in T0-6h to L2 window.",
        },
        {
            "criterion": "eICU external cohort size",
            "threshold": "L2 eligible >= 500",
            "observed": eicu_n,
            "status": "pass" if eicu_n >= 500 else "fail_or_borderline",
            "interpretation": "Needed for external validation feasibility.",
        },
        {
            "criterion": "eICU external events",
            "threshold": "48h failures >= 100",
            "observed": eicu_events,
            "status": "pass" if eicu_events >= 100 else "fail_or_borderline",
            "interpretation": "Needed for stable external performance estimates.",
        },
        {
            "criterion": "eICU ROX-like availability",
            "threshold": "complete rows >= 250",
            "observed": eicu_rox,
            "status": "pass" if eicu_rox >= 250 else "needs_feature_review",
            "interpretation": "SaO2/FiO2/respiration/heart rate in T0-6h to L2 window.",
        },
        {
            "criterion": "Construct boundary",
            "threshold": "No ICU PEF-equivalence claim",
            "observed": "ICU translation layer only",
            "status": "pass",
            "interpretation": "Do not market as direct PEF third-database validation.",
        },
    ]


def build_log(mimic_rows: list[dict[str, object]], eicu_rows: list[dict[str, object]], go_rows: list[dict[str, object]]) -> None:
    def find(rows: list[dict[str, object]], dataset: str, metric: str) -> object:
        for row in rows:
            if row.get("dataset") == dataset and row.get("metric") == metric:
                return row.get("value")
        return ""

    lines = [
        "# ICU Translation Feasibility V21.0",
        "",
        f"- Run date: {date.today().isoformat()}.",
        "- Scope: aggregate-only MIMIC-IV/eICU feasibility for the ICU translation layer.",
        "- Design: first eligible HFNC/NIV-type support start as T0; L2 primary landmark; outcome is invasive ventilation or death from L2 to T0+48h.",
        "- Boundary: this is not direct validation of lower-than-expected field PEF and does not use ICU Peak Exp Flow as a primary exposure.",
        "- Data governance: no row-level derived cohort is exported.",
        "",
        "## MIMIC-IV Snapshot",
        "",
        f"- Adult ICU stays: {find(mimic_rows, 'MIMIC-IV', 'adult_icu_stays')}.",
        f"- First eligible HFNC/NIV stays: {find(mimic_rows, 'MIMIC-IV', 'first_eligible_hfnc_or_niv_stays')}.",
        f"- L2 eligible: {find(mimic_rows, 'MIMIC-IV', 'l2_eligible')}.",
        f"- 48h failures after L2: {find(mimic_rows, 'MIMIC-IV', 'failure_48h_l2')}.",
        f"- ROX-like complete rows: {find(mimic_rows, 'MIMIC-IV', 'l2_complete::required_rox_like')}.",
        "",
        "## eICU-CRD Snapshot",
        "",
        f"- Adult patient-unit stays: {find(eicu_rows, 'eICU-CRD', 'adult_patient_unit_stays')}.",
        f"- Eligible HFNC/NIV-type unit stays: {find(eicu_rows, 'eICU-CRD', 'eligible_hfnc_or_niv_unit_stays')}.",
        f"- L2 eligible: {find(eicu_rows, 'eICU-CRD', 'l2_eligible')}.",
        f"- 48h failures after L2: {find(eicu_rows, 'eICU-CRD', 'failure_48h_l2')}.",
        f"- ROX-like complete rows: {find(eicu_rows, 'eICU-CRD', 'l2_complete::required_rox_like')}.",
        "",
        "## Go/No-Go",
        "",
    ]
    for row in go_rows:
        lines.append(f"- {row['criterion']}: {row['status']} (observed: {row['observed']}; threshold: {row['threshold']}).")
    lines.extend(
        [
            "",
            "## Immediate Interpretation",
            "",
            "- If both MIMIC and eICU meet size/event gates, proceed to feature cleaning and transparent baseline models.",
            "- If eICU T0 mapping is dominated by broad CPAP/PEEP treatment strings and has few HFNC labels, use eICU as broader respiratory-support external validation, not exact HFNC replication.",
            "- If ROX-like completeness is low, the next step is targeted feature extraction/imputation rules before any model comparison.",
            "",
            "PASS: ICU translation feasibility V21.0 built.",
        ]
    )
    (LOG_DIR / f"icu_translation_feasibility_{VERSION}.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def build_manuscript_note(go_rows: list[dict[str, object]]) -> None:
    lines = [
        "# ICU Translation Layer Status Note V21.0",
        "",
        "The MIMIC/eICU extension is framed as an ICU translation layer rather than a third-database validation of the population peak-flow marker. The intended question is whether early physiologic and support-response signatures after HFNC/NIV-type support initiation identify short-term respiratory support failure risk.",
        "",
        "The V21.0 feasibility build uses T0 as first eligible HFNC/NIV-type support start, L2 as the primary landmark, and invasive ventilation or death from L2 to T0+48 hours as the primary acute-care endpoint. Any manuscript language should remain conditional until feature cleaning, model development, and external validation are completed.",
        "",
        "## Current Gate",
        "",
    ]
    for row in go_rows:
        lines.append(f"- {row['criterion']}: {row['status']}; observed {row['observed']} against {row['threshold']}.")
    lines.append("")
    lines.append("Do not describe this module as direct PEF replication or as evidence that ICU peak-flow measurements are equivalent to CHARLS/NHANES field PEF.")
    (MANUSCRIPT_DIR / f"icu_translation_layer_status_note_{VERSION}.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    ensure_dirs()

    print("Loading MIMIC core tables...")
    mimic_core = load_mimic_core()
    adult_stays = set(int(x) for x in mimic_core.loc[mimic_core["adult"], "stay_id"])
    print("Scanning MIMIC procedureevents...")
    mimic_proc_niv, mimic_proc_imv = extract_mimic_procedure_events(adult_stays)
    print("Scanning MIMIC O2 device chartevents for HFNC/NIV T0...")
    mimic_device_t0, mimic_o2_values, mimic_o2_scan_rows = scan_mimic_o2_device(adult_stays)
    mimic_t0 = build_mimic_t0(mimic_core, mimic_proc_niv, mimic_device_t0)
    mimic_t0 = add_mimic_outcomes(mimic_t0, mimic_proc_imv)
    print("Scanning MIMIC physiology windows...")
    mimic_t0, mimic_phys_rows = scan_mimic_physiology(mimic_t0)
    mimic_rows = mimic_counts(mimic_core, mimic_t0, mimic_o2_scan_rows + mimic_phys_rows)

    print("Loading eICU patient table...")
    eicu_patient = load_eicu_patient()
    adult_units = set(int(x) for x in eicu_patient.loc[eicu_patient["adult"], "patientunitstayid"])
    print("Scanning eICU treatment support strings...")
    eicu_treat_t0, eicu_treat_imv, eicu_treat_values = scan_eicu_treatment(adult_units)
    print("Scanning eICU respiratoryCare...")
    eicu_respcare_imv = scan_eicu_respiratory_care(adult_units)
    print("Scanning eICU respiratoryCharting labels...")
    eicu_respchart_t0, eicu_respchart_imv, eicu_respchart_labels = scan_eicu_respchart(adult_units)
    eicu_t0 = build_eicu_t0(eicu_patient, eicu_treat_t0, eicu_respchart_t0, eicu_treat_imv, eicu_respcare_imv, eicu_respchart_imv)
    eicu_t0 = add_eicu_outcomes(eicu_t0)
    print("Scanning eICU physiology windows...")
    eicu_t0, eicu_phys_rows = scan_eicu_physiology(eicu_t0)
    eicu_rows = eicu_counts(eicu_patient, eicu_t0, eicu_phys_rows)

    variable_lock = build_variable_lock(mimic_rows, eicu_rows)
    go_rows = build_go_no_go(mimic_rows, eicu_rows)
    support_value_rows = mimic_o2_values + eicu_treat_values + eicu_respchart_labels

    write_csv(
        METADATA_DIR / f"icu_translation_variable_lock_{VERSION}.csv",
        variable_lock,
        ["dataset", "domain", "source_table", "variable_or_item", "role", "lock_status", "notes"],
    )
    write_csv(TABLE_DIR / f"mimiciv_icu_translation_feasibility_counts_{VERSION}.csv", mimic_rows, ["dataset", "metric", "value", "note"])
    write_csv(TABLE_DIR / f"eicu_icu_translation_feasibility_counts_{VERSION}.csv", eicu_rows, ["dataset", "metric", "value", "note"])
    write_csv(
        TABLE_DIR / f"icu_translation_support_value_audit_{VERSION}.csv",
        support_value_rows,
        ["dataset", "source", "value", "n_events", "support_class"],
    )
    write_csv(TABLE_DIR / f"icu_translation_go_no_go_{VERSION}.csv", go_rows, ["criterion", "threshold", "observed", "status", "interpretation"])
    build_log(mimic_rows, eicu_rows, go_rows)
    build_manuscript_note(go_rows)
    print("Wrote ICU translation feasibility V21.0 outputs.")


if __name__ == "__main__":
    main()
