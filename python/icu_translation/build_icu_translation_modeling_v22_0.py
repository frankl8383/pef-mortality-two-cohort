#!/usr/bin/env python3
"""Build V22 ICU translation modeling datasets and transparent baselines.

V22 converts the V21 feasibility design into local-only row-level modeling
sets, then fits simple MIMIC-derived logistic models and evaluates them in a
MIMIC holdout split and eICU external validation. Aggregate outputs only are
written under results/ and metadata/.
"""

from __future__ import annotations

import csv
import gzip
import importlib.util
import math
import os
import re
from collections import defaultdict
from datetime import date
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[2]
V21_PATH = PROJECT_ROOT / "python" / "icu_translation" / "build_icu_translation_feasibility_v21_0.py"
MIMIC_ROOT = Path.home() / "secure_data" / "mimiciv" / "3.1"
EICU_ROOT = Path.home() / "secure_data" / "eicu-crd" / "2.0"
SENSITIVE_DIR = PROJECT_ROOT / "derived_sensitive" / "icu_translation"
TABLE_DIR = PROJECT_ROOT / "results" / "tables"
LOG_DIR = PROJECT_ROOT / "results" / "logs"
METADATA_DIR = PROJECT_ROOT / "metadata"
MANUSCRIPT_DIR = PROJECT_ROOT / "manuscript"
VERSION = "v22_0"
CHUNK_ROWS = int(os.environ.get("ICU_TRANSLATION_CHUNK_ROWS", "750000"))
RANDOM_SEED = 42

MIMIC_ITEM_TO_VAR = {
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
}

COMMON_VARS = ["spo2", "fio2", "respiratory_rate", "heart_rate", "oxygen_flow", "peep"]
FEATURE_SETS = {
    "clinical_basic": ["age", "male", "support_hfnc"],
    "rox_only": ["rox_l2"],
    "dynamic_core": [
        "age",
        "male",
        "support_hfnc",
        "spo2_last",
        "fio2_last",
        "respiratory_rate_last",
        "heart_rate_last",
        "oxygen_flow_max",
        "peep_max",
        "rox_l2",
    ],
}


def load_v21_module():
    spec = importlib.util.spec_from_file_location("icu_v21", V21_PATH)
    if spec is None or spec.loader is None:
        raise RuntimeError("Could not load V21 feasibility module.")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


v21 = load_v21_module()


def ensure_dirs() -> None:
    for directory in [SENSITIVE_DIR, TABLE_DIR, LOG_DIR, METADATA_DIR, MANUSCRIPT_DIR]:
        directory.mkdir(parents=True, exist_ok=True)


def write_csv(path: Path, rows: list[dict[str, object]], fieldnames: list[str] | None = None) -> None:
    if fieldnames is None:
        fieldnames = list(rows[0].keys()) if rows else ["status"]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def normalize_fio2(value: object) -> float:
    try:
        x = float(value)
    except Exception:
        return np.nan
    if not np.isfinite(x) or x <= 0:
        return np.nan
    if x <= 1.0:
        frac = x
    elif x <= 100.0:
        frac = x / 100.0
    else:
        return np.nan
    if frac < 0.20 or frac > 1.0:
        return np.nan
    return frac


def clean_value(var: str, value: object) -> float:
    try:
        x = float(value)
    except Exception:
        return np.nan
    if not np.isfinite(x):
        return np.nan
    if var == "fio2":
        return normalize_fio2(x)
    ranges = {
        "spo2": (50, 100),
        "respiratory_rate": (1, 80),
        "heart_rate": (20, 250),
        "oxygen_flow": (0, 100),
        "peep": (0, 40),
    }
    lo, hi = ranges.get(var, (-np.inf, np.inf))
    if x < lo or x > hi:
        return np.nan
    return x


def update_agg(store: dict[tuple[int, str], dict[str, object]], entity_id: int, var: str, value: float, time_value: object) -> None:
    if not np.isfinite(value):
        return
    key = (int(entity_id), var)
    entry = store.get(key)
    if entry is None:
        store[key] = {
            "count": 1,
            "sum": float(value),
            "min": float(value),
            "max": float(value),
            "last": float(value),
            "last_time": time_value,
        }
        return
    entry["count"] = int(entry["count"]) + 1
    entry["sum"] = float(entry["sum"]) + float(value)
    entry["min"] = min(float(entry["min"]), float(value))
    entry["max"] = max(float(entry["max"]), float(value))
    if pd.isna(entry["last_time"]) or time_value >= entry["last_time"]:
        entry["last"] = float(value)
        entry["last_time"] = time_value


def add_agg_features(df: pd.DataFrame, id_col: str, agg: dict[tuple[int, str], dict[str, object]]) -> pd.DataFrame:
    out = df.copy()
    for var in COMMON_VARS:
        for suffix in ["last", "min", "max", "mean", "count"]:
            out[f"{var}_{suffix}"] = np.nan
    for idx, entity_id in out[id_col].items():
        for var in COMMON_VARS:
            entry = agg.get((int(entity_id), var))
            if entry is None:
                continue
            out.at[idx, f"{var}_last"] = float(entry["last"])
            out.at[idx, f"{var}_min"] = float(entry["min"])
            out.at[idx, f"{var}_max"] = float(entry["max"])
            out.at[idx, f"{var}_mean"] = float(entry["sum"]) / int(entry["count"])
            out.at[idx, f"{var}_count"] = int(entry["count"])
    out["rox_l2"] = (out["spo2_last"] / out["fio2_last"]) / out["respiratory_rate_last"]
    out.loc[~np.isfinite(out["rox_l2"]), "rox_l2"] = np.nan
    return out


def scan_mimic_feature_values(t0: pd.DataFrame) -> tuple[pd.DataFrame, list[dict[str, object]]]:
    candidate = t0[t0["l2_eligible"]].copy()
    stay_ids = set(int(x) for x in candidate["stay_id"])
    if not stay_ids:
        return candidate, []
    t0_index = candidate.set_index("stay_id")["t0"]
    agg: dict[tuple[int, str], dict[str, object]] = {}
    event_counts = defaultdict(int)
    path = MIMIC_ROOT / "icu" / "chartevents.csv.gz"
    usecols = ["stay_id", "charttime", "itemid", "valuenum"]
    scanned = 0
    matched = 0
    in_window = 0
    for chunk in pd.read_csv(path, usecols=usecols, chunksize=CHUNK_ROWS):
        scanned += len(chunk)
        chunk = chunk[chunk["itemid"].isin(MIMIC_ITEM_TO_VAR) & chunk["stay_id"].isin(stay_ids)].copy()
        if chunk.empty:
            continue
        matched += len(chunk)
        chunk["charttime"] = pd.to_datetime(chunk["charttime"], errors="coerce")
        chunk["t0"] = chunk["stay_id"].map(t0_index)
        chunk = chunk[(chunk["charttime"] >= chunk["t0"] - pd.Timedelta(hours=6)) & (chunk["charttime"] <= chunk["t0"] + pd.Timedelta(hours=2))]
        if chunk.empty:
            continue
        in_window += len(chunk)
        chunk["var"] = chunk["itemid"].map(MIMIC_ITEM_TO_VAR)
        for row in chunk.itertuples(index=False):
            var = str(row.var)
            value = clean_value(var, row.valuenum)
            if np.isfinite(value):
                update_agg(agg, int(row.stay_id), var, value, row.charttime)
                event_counts[var] += 1
    out = add_agg_features(candidate, "stay_id", agg)
    rows = [
        {"dataset": "MIMIC-IV", "metric": "feature_scan_rows_total", "value": scanned, "note": "chartevents rows streamed"},
        {"dataset": "MIMIC-IV", "metric": "feature_candidate_rows", "value": matched, "note": "candidate itemids in L2 cohort stays before window filter"},
        {"dataset": "MIMIC-IV", "metric": "feature_rows_in_l2_window", "value": in_window, "note": "candidate item rows in T0-6h to L2"},
    ]
    for var in COMMON_VARS:
        rows.append({"dataset": "MIMIC-IV", "metric": f"clean_feature_events::{var}", "value": event_counts[var], "note": "cleaned numeric events used in aggregation"})
    return out, rows


def build_mimic_modeling_set() -> tuple[pd.DataFrame, list[dict[str, object]]]:
    core = v21.load_mimic_core()
    adult_stays = set(int(x) for x in core.loc[core["adult"], "stay_id"])
    proc_niv, proc_imv = v21.extract_mimic_procedure_events(adult_stays)
    device_t0, _, scan_rows = v21.scan_mimic_o2_device(adult_stays)
    t0 = v21.build_mimic_t0(core, proc_niv, device_t0)
    t0 = v21.add_mimic_outcomes(t0, proc_imv)
    t0, feature_rows = scan_mimic_feature_values(t0)
    out = pd.DataFrame(
        {
            "dataset": "MIMIC-IV",
            "subject_id": t0["subject_id"].astype("int64"),
            "stay_id": t0["stay_id"].astype("int64"),
            "outcome_failure_48h_l2": t0["failure_48h_l2"].astype(int),
            "age": pd.to_numeric(t0["anchor_age"], errors="coerce"),
            "male": t0["gender"].astype(str).str.upper().eq("M").astype(int),
            "support_hfnc": t0["t0_support_class"].astype(str).str.lower().eq("hfnc").astype(int),
            "support_class": t0["t0_support_class"].astype(str),
            "first_careunit": t0["first_careunit"].astype(str),
        }
    )
    feature_cols = [col for col in t0.columns if any(col.startswith(f"{var}_") for var in COMMON_VARS) or col == "rox_l2"]
    out = pd.concat([out.reset_index(drop=True), t0[feature_cols].reset_index(drop=True)], axis=1)
    return out, scan_rows + feature_rows


def eicu_respchart_var(label: object) -> str:
    text = str(label or "").strip().lower()
    if not text:
        return ""
    if "fio2" in text or "o2 percentage" in text:
        return "fio2"
    if "peep" in text:
        return "peep"
    if "lpm o2" in text or "o2 flow" in text:
        return "oxygen_flow"
    if text in {"total rr", "rr (patient)", "vent rate"} or "respiratory rate" in text:
        return "respiratory_rate"
    return ""


def scan_eicu_feature_values(t0: pd.DataFrame) -> tuple[pd.DataFrame, list[dict[str, object]]]:
    candidate = t0[t0["l2_eligible"]].copy()
    unit_ids = set(int(x) for x in candidate["patientunitstayid"])
    if not unit_ids:
        return candidate, []
    t0_index = candidate.set_index("patientunitstayid")["t0_offset"]
    agg: dict[tuple[int, str], dict[str, object]] = {}
    event_counts = defaultdict(int)
    vital_scanned = 0
    vital_matched = 0
    vital_window = 0

    vital_path = EICU_ROOT / "vitalPeriodic.csv.gz"
    vital_cols = ["patientunitstayid", "observationoffset", "sao2", "heartrate", "respiration"]
    for chunk in pd.read_csv(vital_path, usecols=vital_cols, chunksize=CHUNK_ROWS):
        vital_scanned += len(chunk)
        chunk = chunk[chunk["patientunitstayid"].isin(unit_ids)].copy()
        if chunk.empty:
            continue
        vital_matched += len(chunk)
        chunk["t0"] = chunk["patientunitstayid"].map(t0_index)
        chunk["observationoffset"] = pd.to_numeric(chunk["observationoffset"], errors="coerce")
        chunk = chunk[(chunk["observationoffset"] >= chunk["t0"] - 360) & (chunk["observationoffset"] <= chunk["t0"] + 120)]
        if chunk.empty:
            continue
        vital_window += len(chunk)
        for source_col, var in [("sao2", "spo2"), ("heartrate", "heart_rate"), ("respiration", "respiratory_rate")]:
            values = pd.to_numeric(chunk[source_col], errors="coerce")
            for unit, offset, value in zip(chunk["patientunitstayid"], chunk["observationoffset"], values):
                clean = clean_value(var, value)
                if np.isfinite(clean):
                    update_agg(agg, int(unit), var, clean, float(offset))
                    event_counts[var] += 1

    chart_scanned = 0
    chart_matched = 0
    chart_window = 0
    chart_path = EICU_ROOT / "respiratoryCharting.csv.gz"
    chart_cols = ["patientunitstayid", "respchartoffset", "respchartvaluelabel", "respchartvalue"]
    for chunk in pd.read_csv(chart_path, usecols=chart_cols, dtype=str, chunksize=CHUNK_ROWS):
        chart_scanned += len(chunk)
        chunk["patientunitstayid"] = pd.to_numeric(chunk["patientunitstayid"], errors="coerce")
        chunk = chunk[chunk["patientunitstayid"].isin(unit_ids)].copy()
        if chunk.empty:
            continue
        chart_matched += len(chunk)
        chunk["offset"] = pd.to_numeric(chunk["respchartoffset"], errors="coerce")
        chunk["t0"] = chunk["patientunitstayid"].map(t0_index)
        chunk = chunk[(chunk["offset"] >= chunk["t0"] - 360) & (chunk["offset"] <= chunk["t0"] + 120)]
        if chunk.empty:
            continue
        chart_window += len(chunk)
        chunk["var"] = chunk["respchartvaluelabel"].map(eicu_respchart_var)
        chunk = chunk[chunk["var"] != ""]
        values = pd.to_numeric(chunk["respchartvalue"], errors="coerce")
        for unit, offset, var, value in zip(chunk["patientunitstayid"], chunk["offset"], chunk["var"], values):
            clean = clean_value(str(var), value)
            if np.isfinite(clean):
                update_agg(agg, int(unit), str(var), clean, float(offset))
                event_counts[str(var)] += 1

    out = add_agg_features(candidate, "patientunitstayid", agg)
    rows = [
        {"dataset": "eICU-CRD", "metric": "feature_vital_rows_total", "value": vital_scanned, "note": "vitalPeriodic rows streamed"},
        {"dataset": "eICU-CRD", "metric": "feature_vital_candidate_rows", "value": vital_matched, "note": "rows in L2 cohort units before window filter"},
        {"dataset": "eICU-CRD", "metric": "feature_vital_rows_in_l2_window", "value": vital_window, "note": "vital rows in T0-6h to L2"},
        {"dataset": "eICU-CRD", "metric": "feature_respchart_rows_total", "value": chart_scanned, "note": "respiratoryCharting rows streamed"},
        {"dataset": "eICU-CRD", "metric": "feature_respchart_candidate_rows", "value": chart_matched, "note": "rows in L2 cohort units before window filter"},
        {"dataset": "eICU-CRD", "metric": "feature_respchart_rows_in_l2_window", "value": chart_window, "note": "respiratoryCharting rows in T0-6h to L2"},
    ]
    for var in COMMON_VARS:
        rows.append({"dataset": "eICU-CRD", "metric": f"clean_feature_events::{var}", "value": event_counts[var], "note": "cleaned numeric events used in aggregation"})
    return out, rows


def build_eicu_modeling_set() -> tuple[pd.DataFrame, list[dict[str, object]]]:
    patient = v21.load_eicu_patient()
    adult_units = set(int(x) for x in patient.loc[patient["adult"], "patientunitstayid"])
    treatment_t0, treatment_imv, _ = v21.scan_eicu_treatment(adult_units)
    respcare_imv = v21.scan_eicu_respiratory_care(adult_units)
    respchart_t0, respchart_imv, _ = v21.scan_eicu_respchart(adult_units)
    t0 = v21.build_eicu_t0(patient, treatment_t0, respchart_t0, treatment_imv, respcare_imv, respchart_imv)
    t0 = v21.add_eicu_outcomes(t0)
    t0, feature_rows = scan_eicu_feature_values(t0)
    out = pd.DataFrame(
        {
            "dataset": "eICU-CRD",
            "patientunitstayid": t0["patientunitstayid"].astype("int64"),
            "outcome_failure_48h_l2": t0["failure_48h_l2"].astype(int),
            "age": pd.to_numeric(t0["age_num"], errors="coerce"),
            "male": t0["gender"].astype(str).str.lower().str.startswith("male").astype(int),
            "support_hfnc": t0["t0_support_class"].astype(str).str.lower().eq("hfnc").astype(int),
            "support_class": t0["t0_support_class"].astype(str),
            "unittype": t0["unittype"].astype(str),
        }
    )
    feature_cols = [col for col in t0.columns if any(col.startswith(f"{var}_") for var in COMMON_VARS) or col == "rox_l2"]
    out = pd.concat([out.reset_index(drop=True), t0[feature_cols].reset_index(drop=True)], axis=1)
    return out, feature_rows


def stratified_split(y: np.ndarray, train_fraction: float = 0.70, seed: int = RANDOM_SEED) -> np.ndarray:
    rng = np.random.default_rng(seed)
    split = np.array(["test"] * len(y), dtype=object)
    for klass in [0, 1]:
        idx = np.where(y == klass)[0]
        rng.shuffle(idx)
        n_train = int(math.floor(len(idx) * train_fraction))
        split[idx[:n_train]] = "train"
    return split


def sigmoid(x: np.ndarray) -> np.ndarray:
    return 1 / (1 + np.exp(-np.clip(x, -35, 35)))


def fit_logistic(X: np.ndarray, y: np.ndarray, penalty: float = 1.0, max_iter: int = 100) -> np.ndarray:
    beta = np.zeros(X.shape[1])
    penalty_vec = np.ones(X.shape[1]) * penalty
    penalty_vec[0] = 0.0
    for _ in range(max_iter):
        eta = X @ beta
        p = sigmoid(eta)
        w = np.clip(p * (1 - p), 1e-6, None)
        grad = X.T @ (y - p) - penalty_vec * beta
        hess = (X.T * w) @ X + np.diag(penalty_vec)
        try:
            step = np.linalg.solve(hess, grad)
        except np.linalg.LinAlgError:
            step = np.linalg.lstsq(hess, grad, rcond=None)[0]
        beta_new = beta + step
        if np.max(np.abs(step)) < 1e-6:
            beta = beta_new
            break
        beta = beta_new
    return beta


def rank_auc(y: np.ndarray, p: np.ndarray) -> float:
    y = np.asarray(y).astype(int)
    p = np.asarray(p).astype(float)
    n_pos = int(y.sum())
    n_neg = len(y) - n_pos
    if n_pos == 0 or n_neg == 0:
        return np.nan
    ranks = pd.Series(p).rank(method="average").to_numpy()
    sum_pos = ranks[y == 1].sum()
    return float((sum_pos - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg))


def average_precision(y: np.ndarray, p: np.ndarray) -> float:
    y = np.asarray(y).astype(int)
    p = np.asarray(p).astype(float)
    n_pos = int(y.sum())
    if n_pos == 0:
        return np.nan
    order = np.argsort(-p)
    y_sorted = y[order]
    tp = np.cumsum(y_sorted)
    fp = np.cumsum(1 - y_sorted)
    precision = tp / np.maximum(tp + fp, 1)
    recall_step = y_sorted / n_pos
    return float(np.sum(precision * recall_step))


def logit(p: np.ndarray | float) -> np.ndarray | float:
    arr = np.clip(p, 1e-6, 1 - 1e-6)
    return np.log(arr / (1 - arr))


def calibration_slope(y: np.ndarray, p: np.ndarray) -> float:
    lp = np.asarray(logit(p), dtype=float)
    X = np.column_stack([np.ones(len(lp)), lp])
    try:
        beta = fit_logistic(X, y.astype(float), penalty=0.0, max_iter=100)
        return float(beta[1])
    except Exception:
        return np.nan


def metric_row(y: np.ndarray, p: np.ndarray) -> dict[str, object]:
    y = y.astype(int)
    p = np.clip(p.astype(float), 1e-6, 1 - 1e-6)
    obs = float(np.mean(y))
    mean_pred = float(np.mean(p))
    return {
        "n": int(len(y)),
        "events": int(y.sum()),
        "event_rate": obs,
        "mean_predicted_risk": mean_pred,
        "auroc": rank_auc(y, p),
        "auprc": average_precision(y, p),
        "brier": float(np.mean((p - y) ** 2)),
        "calibration_intercept_approx": float(logit(obs) - logit(mean_pred)) if 0 < obs < 1 else np.nan,
        "calibration_slope": calibration_slope(y, p),
    }


def prepare_design(
    df: pd.DataFrame,
    features: list[str],
    train_stats: dict[str, dict[str, float]] | None = None,
    impute_with_missing_indicators: bool = True,
) -> tuple[np.ndarray, list[str], dict[str, dict[str, float]], pd.Series]:
    base = df.copy()
    valid = pd.Series(True, index=base.index)
    for col in ["outcome_failure_48h_l2"] + features:
        if col not in base.columns:
            base[col] = np.nan
    if not impute_with_missing_indicators:
        for col in features:
            valid &= pd.notna(base[col])
    valid &= pd.notna(base["outcome_failure_48h_l2"])
    base = base.loc[valid].copy()

    stats = train_stats or {}
    columns: list[np.ndarray] = [np.ones(len(base))]
    names = ["intercept"]
    for feature in features:
        values = pd.to_numeric(base[feature], errors="coerce").astype(float)
        binary = feature in {"male", "support_hfnc"}
        if train_stats is None:
            median = float(values.median()) if values.notna().any() else 0.0
            mean = float(values.fillna(median).mean())
            sd = float(values.fillna(median).std(ddof=0))
            if not np.isfinite(sd) or sd == 0 or binary:
                sd = 1.0
            stats[feature] = {"median": median, "mean": mean, "sd": sd, "binary": float(binary)}
        med = stats[feature]["median"]
        mean = stats[feature]["mean"]
        sd = stats[feature]["sd"]
        missing = values.isna().astype(float).to_numpy()
        filled = values.fillna(med).to_numpy()
        scaled = filled if bool(stats[feature].get("binary", 0.0)) else (filled - mean) / sd
        columns.append(scaled.astype(float))
        names.append(feature)
        if impute_with_missing_indicators and not bool(stats[feature].get("binary", 0.0)):
            columns.append(missing)
            names.append(f"{feature}_missing")
    X = np.column_stack(columns)
    return X, names, stats, valid


def fit_and_evaluate_models(mimic: pd.DataFrame, eicu: pd.DataFrame) -> tuple[list[dict[str, object]], list[dict[str, object]], pd.DataFrame, pd.DataFrame]:
    mimic = mimic.copy()
    y_all = mimic["outcome_failure_48h_l2"].astype(int).to_numpy()
    mimic["split"] = stratified_split(y_all)
    eicu = eicu.copy()
    eicu["split"] = "external"

    perf_rows: list[dict[str, object]] = []
    coef_rows: list[dict[str, object]] = []
    for model_name, features in FEATURE_SETS.items():
        impute = model_name != "rox_only"
        train_df = mimic[mimic["split"].eq("train")].copy()
        X_train, names, stats, train_valid = prepare_design(train_df, features, train_stats=None, impute_with_missing_indicators=impute)
        y_train = train_df.loc[train_valid, "outcome_failure_48h_l2"].astype(int).to_numpy()
        beta = fit_logistic(X_train, y_train, penalty=1.0)
        for name, value in zip(names, beta):
            coef_rows.append({"model": model_name, "term": name, "coefficient": float(value), "odds_ratio": float(np.exp(np.clip(value, -20, 20)))})

        for eval_name, eval_df in [
            ("mimic_train_apparent", train_df),
            ("mimic_internal_test", mimic[mimic["split"].eq("test")].copy()),
            ("eicu_external", eicu.copy()),
        ]:
            X_eval, _, _, eval_valid = prepare_design(eval_df, features, train_stats=stats, impute_with_missing_indicators=impute)
            y_eval = eval_df.loc[eval_valid, "outcome_failure_48h_l2"].astype(int).to_numpy()
            pred = sigmoid(X_eval @ beta)
            row = metric_row(y_eval, pred)
            row.update({"model": model_name, "evaluation": eval_name, "features": ";".join(features), "imputation": "median_plus_missing_indicators" if impute else "complete_case"})
            perf_rows.append(row)

            # Store predictions locally only.
            pred_col = f"pred_{model_name}"
            if pred_col not in eval_df.columns:
                eval_df[pred_col] = np.nan
            if eval_name.startswith("mimic"):
                mimic.loc[eval_df.loc[eval_valid].index, pred_col] = pred
            else:
                eicu.loc[eval_df.loc[eval_valid].index, pred_col] = pred

    return perf_rows, coef_rows, mimic, eicu


def feature_missingness_rows(df: pd.DataFrame, dataset: str) -> list[dict[str, object]]:
    rows = []
    cols = ["age", "male", "support_hfnc", "rox_l2"]
    for var in COMMON_VARS:
        cols.extend([f"{var}_last", f"{var}_min", f"{var}_max", f"{var}_mean", f"{var}_count"])
    for col in cols:
        if col not in df.columns:
            continue
        missing = int(df[col].isna().sum())
        rows.append(
            {
                "dataset": dataset,
                "feature": col,
                "n": int(len(df)),
                "nonmissing": int(len(df) - missing),
                "missing": missing,
                "missing_pct": missing / len(df) if len(df) else np.nan,
            }
        )
    return rows


def cohort_count_rows(mimic: pd.DataFrame, eicu: pd.DataFrame, scan_rows: list[dict[str, object]]) -> list[dict[str, object]]:
    rows = []
    for dataset, df in [("MIMIC-IV", mimic), ("eICU-CRD", eicu)]:
        rows.append({"dataset": dataset, "metric": "modeling_rows_l2_eligible", "value": int(len(df)), "note": "row-level local modeling set rows"})
        rows.append({"dataset": dataset, "metric": "failure_48h_l2_events", "value": int(df["outcome_failure_48h_l2"].sum()), "note": "primary V22 outcome events"})
        rows.append({"dataset": dataset, "metric": "failure_48h_l2_event_rate", "value": float(df["outcome_failure_48h_l2"].mean()), "note": "primary V22 outcome rate"})
        if dataset == "MIMIC-IV":
            for split, count in df["split"].value_counts().items():
                rows.append({"dataset": dataset, "metric": f"split_rows::{split}", "value": int(count), "note": "stratified 70/30 split"})
                events = int(df.loc[df["split"].eq(split), "outcome_failure_48h_l2"].sum())
                rows.append({"dataset": dataset, "metric": f"split_events::{split}", "value": events, "note": "events within split"})
    rows.extend(scan_rows)
    return rows


def feature_map_rows() -> list[dict[str, object]]:
    return [
        {"feature": "age", "MIMIC-IV source": "hosp/patients.anchor_age", "eICU source": "patient.age", "role": "clinical static predictor", "cleaning": "numeric years; eICU >89 treated as 89"},
        {"feature": "male", "MIMIC-IV source": "hosp/patients.gender", "eICU source": "patient.gender", "role": "clinical static predictor", "cleaning": "binary male indicator"},
        {"feature": "support_hfnc", "MIMIC-IV source": "O2 Delivery Device(s); procedureevents", "eICU source": "treatment; respiratoryCharting", "role": "support type at T0", "cleaning": "HFNC vs NIV-type start; eICU HFNC sparse"},
        {"feature": "spo2_last/min/mean/max", "MIMIC-IV source": "chartevents 220277", "eICU source": "vitalPeriodic.sao2", "role": "oxygenation response", "cleaning": "50-100 percent in T0-6h to L2"},
        {"feature": "fio2_last/min/mean/max", "MIMIC-IV source": "chartevents 223835", "eICU source": "respiratoryCharting FiO2/O2 Percentage", "role": "oxygen intensity", "cleaning": "converted to fraction 0.20-1.00"},
        {"feature": "respiratory_rate_last/min/mean/max", "MIMIC-IV source": "chartevents 220210/224690/224689", "eICU source": "vitalPeriodic.respiration; respiratoryCharting RR labels", "role": "respiratory load", "cleaning": "1-80 breaths/min"},
        {"feature": "heart_rate_last/min/mean/max", "MIMIC-IV source": "chartevents 220045", "eICU source": "vitalPeriodic.heartrate", "role": "physiologic stress", "cleaning": "20-250 bpm"},
        {"feature": "oxygen_flow_last/min/mean/max", "MIMIC-IV source": "chartevents O2 flow items", "eICU source": "respiratoryCharting LPM O2", "role": "support intensity", "cleaning": "0-100 L/min"},
        {"feature": "peep_last/min/mean/max", "MIMIC-IV source": "chartevents 220339/224700", "eICU source": "respiratoryCharting PEEP/PEEP-CPAP", "role": "positive pressure intensity", "cleaning": "0-40 cmH2O"},
        {"feature": "rox_l2", "MIMIC-IV source": "derived from SpO2, FiO2, RR", "eICU source": "derived from SaO2, FiO2, respiration", "role": "ROX-like comparator", "cleaning": "(SpO2 / FiO2 fraction) / respiratory rate using last values"},
    ]


def write_local_modeling_sets(mimic: pd.DataFrame, eicu: pd.DataFrame) -> None:
    mimic.to_csv(SENSITIVE_DIR / f"mimiciv_icu_translation_modeling_{VERSION}.csv.gz", index=False, compression="gzip")
    eicu.to_csv(SENSITIVE_DIR / f"eicu_icu_translation_modeling_{VERSION}.csv.gz", index=False, compression="gzip")


def build_log(perf_rows: list[dict[str, object]], mimic: pd.DataFrame, eicu: pd.DataFrame) -> None:
    perf = pd.DataFrame(perf_rows)
    lines = [
        "# ICU Translation Modeling V22.0",
        "",
        f"- Run date: {date.today().isoformat()}.",
        "- Scope: first transparent MIMIC-derived prediction baselines with eICU external validation.",
        "- Design: train on a stratified 70% MIMIC split; evaluate MIMIC 30% holdout and eICU external validation.",
        "- Primary endpoint: invasive ventilation or death from L2 to T0+48h.",
        "- Boundary: ICU translation layer only; not direct validation of field peak-flow vulnerability.",
        "- Data governance: row-level model files were written only to the local ignored sensitive folder; aggregate outputs are in results/ and metadata/.",
        "",
        "## Modeling Cohorts",
        "",
        f"- MIMIC L2-eligible modeling rows: {len(mimic)}; events: {int(mimic['outcome_failure_48h_l2'].sum())}.",
        f"- eICU L2-eligible modeling rows: {len(eicu)}; events: {int(eicu['outcome_failure_48h_l2'].sum())}.",
        "",
        "## Main Performance Snapshot",
        "",
    ]
    for _, row in perf[perf["evaluation"].isin(["mimic_internal_test", "eicu_external"])].iterrows():
        lines.append(
            f"- {row['model']} / {row['evaluation']}: AUROC {float(row['auroc']):.3f}, "
            f"AUPRC {float(row['auprc']):.3f}, Brier {float(row['brier']):.3f}, "
            f"calibration slope {float(row['calibration_slope']):.3f}."
        )
    lines.extend(
        [
            "",
            "## Interpretation Guardrails",
            "",
            "- These are first-pass transparent baselines, not final clinical models.",
            "- eICU mapping remains broader than exact MIMIC HFNC/NIV mapping; report it as external respiratory-support validation unless manually tightened.",
            "- Model coefficients and performance should be reviewed before adding the ICU layer to the main manuscript.",
            "",
            "PASS: ICU translation modeling V22.0 built.",
        ]
    )
    (LOG_DIR / f"icu_translation_modeling_{VERSION}.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def build_status_note(perf_rows: list[dict[str, object]]) -> None:
    perf = pd.DataFrame(perf_rows)
    lines = [
        "# ICU Translation Model Status Note V22.0",
        "",
        "V22.0 turns the MIMIC/eICU feasibility gate into a first transparent prediction analysis. It should still be read as an ICU translation layer, not as a third-database validation of the community peak-flow score.",
        "",
        "## Performance Summary",
        "",
    ]
    for _, row in perf[perf["evaluation"].isin(["mimic_internal_test", "eicu_external"])].iterrows():
        lines.append(
            f"- {row['model']} in {row['evaluation']}: AUROC {float(row['auroc']):.3f}; "
            f"AUPRC {float(row['auprc']):.3f}; Brier {float(row['brier']):.3f}; "
            f"calibration slope {float(row['calibration_slope']):.3f}."
        )
    lines.extend(
        [
            "",
            "The next scientific decision is whether the dynamic model materially improves over the ROX-like comparator and transports acceptably to eICU after calibration review.",
        ]
    )
    (MANUSCRIPT_DIR / f"icu_translation_model_status_note_{VERSION}.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    ensure_dirs()
    print("Building MIMIC V22 modeling set...")
    mimic, mimic_scan_rows = build_mimic_modeling_set()
    print("Building eICU V22 modeling set...")
    eicu, eicu_scan_rows = build_eicu_modeling_set()
    print("Fitting transparent MIMIC-derived models and evaluating eICU external validation...")
    perf_rows, coef_rows, mimic_pred, eicu_pred = fit_and_evaluate_models(mimic, eicu)
    print("Writing local-only row-level modeling sets...")
    write_local_modeling_sets(mimic_pred, eicu_pred)

    scan_rows = mimic_scan_rows + eicu_scan_rows
    cohort_rows = cohort_count_rows(mimic_pred, eicu_pred, scan_rows)
    missing_rows = feature_missingness_rows(mimic_pred, "MIMIC-IV") + feature_missingness_rows(eicu_pred, "eICU-CRD")

    write_csv(METADATA_DIR / f"icu_translation_feature_map_{VERSION}.csv", feature_map_rows())
    write_csv(TABLE_DIR / f"icu_translation_model_performance_{VERSION}.csv", perf_rows)
    write_csv(TABLE_DIR / f"icu_translation_model_coefficients_{VERSION}.csv", coef_rows)
    write_csv(TABLE_DIR / f"icu_translation_model_cohort_counts_{VERSION}.csv", cohort_rows)
    write_csv(TABLE_DIR / f"icu_translation_feature_missingness_{VERSION}.csv", missing_rows)
    build_log(perf_rows, mimic_pred, eicu_pred)
    build_status_note(perf_rows)
    print("Wrote ICU translation modeling V22.0 outputs.")


if __name__ == "__main__":
    main()
