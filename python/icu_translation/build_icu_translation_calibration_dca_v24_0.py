#!/usr/bin/env python3
"""Build V24 ICU translation calibration, DCA, and mapping diagnostics.

V24 starts from the V22 local-only row-level modeling sets and recomputes the
best V23 transportable model plus the ROX-like comparator. It writes aggregate
diagnostic outputs only. Post-hoc recalibration rows are calibration diagnostics,
not new external validation claims.
"""

from __future__ import annotations

import csv
import importlib.util
import math
from datetime import date
from pathlib import Path

import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[2]
V22_PATH = PROJECT_ROOT / "python" / "icu_translation" / "build_icu_translation_modeling_v22_0.py"
SENSITIVE_DIR = PROJECT_ROOT / "derived_sensitive" / "icu_translation"
TABLE_DIR = PROJECT_ROOT / "results" / "tables"
LOG_DIR = PROJECT_ROOT / "results" / "logs"
MANUSCRIPT_DIR = PROJECT_ROOT / "manuscript"
VERSION = "v24_0"

MODEL_SPECS = {
    "dynamic_last_no_support": {
        "features": [
            "age",
            "male",
            "spo2_last",
            "fio2_last",
            "respiratory_rate_last",
            "heart_rate_last",
            "oxygen_flow_max",
            "peep_max",
            "rox_l2",
        ],
        "impute": True,
        "label": "Best V23 dynamic model without support-type indicator",
    },
    "rox_only": {
        "features": ["rox_l2"],
        "impute": False,
        "label": "ROX-like complete-case comparator",
    },
}

EVALUATIONS = ["mimic_internal_test", "eicu_external"]
THRESHOLDS = [x / 100 for x in range(5, 55, 5)]


def load_v22_module():
    spec = importlib.util.spec_from_file_location("icu_v22", V22_PATH)
    if spec is None or spec.loader is None:
        raise RuntimeError("Could not load V22 modeling module.")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


v22 = load_v22_module()


def ensure_dirs() -> None:
    for directory in [TABLE_DIR, LOG_DIR, MANUSCRIPT_DIR]:
        directory.mkdir(parents=True, exist_ok=True)


def write_csv(path: Path, rows: list[dict[str, object]], fieldnames: list[str] | None = None) -> None:
    if fieldnames is None:
        fieldnames = list(rows[0].keys()) if rows else ["status"]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def load_local_sets() -> tuple[pd.DataFrame, pd.DataFrame]:
    mimic = pd.read_csv(SENSITIVE_DIR / "mimiciv_icu_translation_modeling_v22_0.csv.gz")
    eicu = pd.read_csv(SENSITIVE_DIR / "eicu_icu_translation_modeling_v22_0.csv.gz")
    if "split" not in mimic.columns:
        mimic["split"] = v22.stratified_split(mimic["outcome_failure_48h_l2"].astype(int).to_numpy())
    if "split" not in eicu.columns:
        eicu["split"] = "external"
    return mimic, eicu


def safe_metric_row(y: np.ndarray, p: np.ndarray) -> dict[str, object]:
    y = np.asarray(y).astype(int)
    p = np.clip(np.asarray(p).astype(float), 1e-6, 1 - 1e-6)
    if len(y) == 0:
        return {
            "n": 0,
            "events": 0,
            "event_rate": np.nan,
            "mean_predicted_risk": np.nan,
            "auroc": np.nan,
            "auprc": np.nan,
            "brier": np.nan,
            "calibration_intercept_approx": np.nan,
            "calibration_slope": np.nan,
        }
    if len(np.unique(y)) < 2:
        return {
            "n": int(len(y)),
            "events": int(y.sum()),
            "event_rate": float(y.mean()),
            "mean_predicted_risk": float(p.mean()),
            "auroc": np.nan,
            "auprc": np.nan,
            "brier": float(np.mean((p - y) ** 2)),
            "calibration_intercept_approx": np.nan,
            "calibration_slope": np.nan,
        }
    return v22.metric_row(y, p)


def fit_predictions(mimic: pd.DataFrame, eicu: pd.DataFrame) -> tuple[pd.DataFrame, list[dict[str, object]], list[dict[str, object]]]:
    pred_frames: list[pd.DataFrame] = []
    perf_rows: list[dict[str, object]] = []
    coef_rows: list[dict[str, object]] = []
    keep_base = [
        "dataset",
        "outcome_failure_48h_l2",
        "age",
        "male",
        "support_hfnc",
        "support_class",
        "first_careunit",
        "unittype",
        "rox_l2",
    ]

    for model_name, spec in MODEL_SPECS.items():
        features = list(spec["features"])
        impute = bool(spec["impute"])
        train_df = mimic[mimic["split"].eq("train")].copy()
        X_train, names, stats, train_valid = v22.prepare_design(
            train_df,
            features,
            train_stats=None,
            impute_with_missing_indicators=impute,
        )
        y_train = train_df.loc[train_valid, "outcome_failure_48h_l2"].astype(int).to_numpy()
        beta = v22.fit_logistic(X_train, y_train, penalty=1.0)
        for name, value in zip(names, beta):
            coef_rows.append(
                {
                    "model": model_name,
                    "term": name,
                    "coefficient": float(value),
                    "odds_ratio": float(np.exp(np.clip(value, -20, 20))),
                }
            )

        eval_sources = [
            ("mimic_train_apparent", train_df),
            ("mimic_internal_test", mimic[mimic["split"].eq("test")].copy()),
            ("eicu_external", eicu.copy()),
        ]
        for evaluation, eval_df in eval_sources:
            X_eval, _, _, eval_valid = v22.prepare_design(
                eval_df,
                features,
                train_stats=stats,
                impute_with_missing_indicators=impute,
            )
            y_eval = eval_df.loc[eval_valid, "outcome_failure_48h_l2"].astype(int).to_numpy()
            pred = v22.sigmoid(X_eval @ beta)
            metric = safe_metric_row(y_eval, pred)
            metric.update(
                {
                    "model": model_name,
                    "evaluation": evaluation,
                    "features": ";".join(features),
                    "imputation": "median_plus_missing_indicators" if impute else "complete_case",
                    "model_label": spec["label"],
                }
            )
            perf_rows.append(metric)

            cols = [col for col in keep_base if col in eval_df.columns]
            frame = eval_df.loc[eval_valid, cols].copy()
            frame["source_index"] = eval_df.loc[eval_valid].index.astype(int)
            frame["model"] = model_name
            frame["evaluation"] = evaluation
            frame["prediction"] = pred
            pred_frames.append(frame)

    predictions = pd.concat(pred_frames, ignore_index=True)
    predictions["outcome_failure_48h_l2"] = predictions["outcome_failure_48h_l2"].astype(int)
    return predictions, perf_rows, coef_rows


def calibration_bin_rows(predictions: pd.DataFrame) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for (evaluation, model), sub in predictions[predictions["evaluation"].isin(EVALUATIONS)].groupby(["evaluation", "model"], sort=True):
        sub = sub.sort_values("prediction").reset_index(drop=True)
        if sub.empty:
            continue
        n_bins = min(10, len(sub))
        sub["calibration_bin"] = pd.qcut(np.arange(len(sub)), q=n_bins, labels=False) + 1
        for bin_id, group in sub.groupby("calibration_bin", sort=True):
            y = group["outcome_failure_48h_l2"].to_numpy()
            p = group["prediction"].to_numpy()
            rows.append(
                {
                    "evaluation": evaluation,
                    "model": model,
                    "bin": int(bin_id),
                    "n": int(len(group)),
                    "events": int(y.sum()),
                    "observed_event_rate": float(y.mean()),
                    "mean_predicted_risk": float(p.mean()),
                    "expected_events": float(p.sum()),
                    "min_predicted_risk": float(p.min()),
                    "max_predicted_risk": float(p.max()),
                }
            )
    return rows


def fit_intercept_offset(y: np.ndarray, p: np.ndarray, max_iter: int = 100) -> float:
    y = np.asarray(y).astype(float)
    lp = np.asarray(v22.logit(p), dtype=float)
    alpha = 0.0
    for _ in range(max_iter):
        pred = v22.sigmoid(lp + alpha)
        grad = float(np.sum(y - pred))
        info = float(np.sum(pred * (1 - pred)))
        if info <= 1e-9:
            break
        step = grad / info
        alpha += step
        if abs(step) < 1e-8:
            break
    return float(alpha)


def recalibration_rows(predictions: pd.DataFrame) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    fields = [
        "evaluation",
        "model",
        "recalibration_method",
        "n",
        "events",
        "event_rate",
        "mean_predicted_risk",
        "auroc",
        "auprc",
        "brier",
        "calibration_intercept_approx",
        "calibration_slope",
        "posthoc_intercept_alpha",
        "posthoc_slope_beta",
        "note",
    ]
    for (evaluation, model), sub in predictions[predictions["evaluation"].isin(EVALUATIONS)].groupby(["evaluation", "model"], sort=True):
        y = sub["outcome_failure_48h_l2"].to_numpy(dtype=int)
        p = sub["prediction"].to_numpy(dtype=float)
        original = safe_metric_row(y, p)
        original.update(
            {
                "evaluation": evaluation,
                "model": model,
                "recalibration_method": "original_mimic_trained",
                "posthoc_intercept_alpha": np.nan,
                "posthoc_slope_beta": np.nan,
                "note": "MIMIC-trained model evaluated without target-set recalibration.",
            }
        )
        rows.append({key: original.get(key, np.nan) for key in fields})

        if len(np.unique(y)) < 2:
            continue
        alpha = fit_intercept_offset(y, p)
        p_intercept = v22.sigmoid(np.asarray(v22.logit(p), dtype=float) + alpha)
        intercept_metric = safe_metric_row(y, p_intercept)
        intercept_metric.update(
            {
                "evaluation": evaluation,
                "model": model,
                "recalibration_method": "posthoc_intercept_only",
                "posthoc_intercept_alpha": alpha,
                "posthoc_slope_beta": np.nan,
                "note": "Calibration diagnostic fitted on the evaluation-set outcome labels; not an independent validation estimate.",
            }
        )
        rows.append({key: intercept_metric.get(key, np.nan) for key in fields})

        lp = np.asarray(v22.logit(p), dtype=float)
        X = np.column_stack([np.ones(len(lp)), lp])
        beta = v22.fit_logistic(X, y.astype(float), penalty=0.0, max_iter=100)
        p_slope = v22.sigmoid(X @ beta)
        slope_metric = safe_metric_row(y, p_slope)
        slope_metric.update(
            {
                "evaluation": evaluation,
                "model": model,
                "recalibration_method": "posthoc_intercept_and_slope",
                "posthoc_intercept_alpha": float(beta[0]),
                "posthoc_slope_beta": float(beta[1]),
                "note": "Calibration diagnostic fitted on the evaluation-set outcome labels; not an independent validation estimate.",
            }
        )
        rows.append({key: slope_metric.get(key, np.nan) for key in fields})
    return rows


def net_benefit(y: np.ndarray, p: np.ndarray, threshold: float) -> float:
    y = np.asarray(y).astype(int)
    p = np.asarray(p).astype(float)
    positives = p >= threshold
    tp = int(np.sum(positives & (y == 1)))
    fp = int(np.sum(positives & (y == 0)))
    n = len(y)
    return float(tp / n - fp / n * threshold / (1 - threshold))


def dca_rows(predictions: pd.DataFrame) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for evaluation in EVALUATIONS:
        eval_pred = predictions[predictions["evaluation"].eq(evaluation)].copy()
        dyn = eval_pred[eval_pred["model"].eq("dynamic_last_no_support")]
        rox = eval_pred[eval_pred["model"].eq("rox_only")]
        common_ids = sorted(set(dyn["source_index"]).intersection(set(rox["source_index"])))
        if not common_ids:
            continue
        common = eval_pred[eval_pred["source_index"].isin(common_ids)].copy()
        for threshold in THRESHOLDS:
            for model, group in common.groupby("model", sort=True):
                y = group["outcome_failure_48h_l2"].to_numpy(dtype=int)
                p = group["prediction"].to_numpy(dtype=float)
                rows.append(
                    {
                        "evaluation": evaluation,
                        "analysis_population": "common_rox_complete_case",
                        "strategy": "model",
                        "model": model,
                        "threshold": threshold,
                        "n": int(len(group)),
                        "events": int(y.sum()),
                        "event_rate": float(y.mean()),
                        "net_benefit": net_benefit(y, p, threshold),
                    }
                )
            y_ref = common[common["model"].eq("dynamic_last_no_support")]["outcome_failure_48h_l2"].to_numpy(dtype=int)
            prevalence = float(y_ref.mean())
            rows.append(
                {
                    "evaluation": evaluation,
                    "analysis_population": "common_rox_complete_case",
                    "strategy": "treat_all",
                    "model": "",
                    "threshold": threshold,
                    "n": int(len(y_ref)),
                    "events": int(y_ref.sum()),
                    "event_rate": prevalence,
                    "net_benefit": prevalence - (1 - prevalence) * threshold / (1 - threshold),
                }
            )
            rows.append(
                {
                    "evaluation": evaluation,
                    "analysis_population": "common_rox_complete_case",
                    "strategy": "treat_none",
                    "model": "",
                    "threshold": threshold,
                    "n": int(len(y_ref)),
                    "events": int(y_ref.sum()),
                    "event_rate": prevalence,
                    "net_benefit": 0.0,
                }
            )
    return rows


def subgroup_metric_row(evaluation: str, model: str, group_family: str, subgroup: str, group: pd.DataFrame) -> dict[str, object]:
    y = group["outcome_failure_48h_l2"].to_numpy(dtype=int)
    p = group["prediction"].to_numpy(dtype=float)
    metric = safe_metric_row(y, p)
    if len(group) < 100 or int(y.sum()) < 10 or int(len(y) - y.sum()) < 10:
        flag = "interpret_cautiously_sparse"
    else:
        flag = "adequate_for_descriptive_transport_check"
    metric.update(
        {
            "evaluation": evaluation,
            "model": model,
            "group_family": group_family,
            "subgroup": subgroup,
            "interpretation_flag": flag,
        }
    )
    return metric


def age_group(age: object) -> str:
    try:
        value = float(age)
    except Exception:
        return "age_unknown"
    if not math.isfinite(value):
        return "age_unknown"
    if value < 65:
        return "age_lt65"
    if value < 80:
        return "age_65_79"
    return "age_ge80"


def subgroup_rows(predictions: pd.DataFrame) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for (evaluation, model), sub in predictions[predictions["evaluation"].isin(EVALUATIONS)].groupby(["evaluation", "model"], sort=True):
        sub = sub.copy()
        sub["support_group"] = sub["support_class"].astype(str).str.lower().replace({"nan": "support_unknown"})
        sub["age_group"] = sub["age"].map(age_group)
        sub["sex_group"] = np.where(pd.to_numeric(sub["male"], errors="coerce").eq(1), "male", "female")
        for group_family, column in [
            ("support_class", "support_group"),
            ("age_group", "age_group"),
            ("sex", "sex_group"),
        ]:
            for subgroup, group in sub.groupby(column, sort=True):
                rows.append(subgroup_metric_row(evaluation, model, group_family, str(subgroup), group))
    return rows


def endpoint_mapping_rows(predictions: pd.DataFrame) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    base = predictions[predictions["evaluation"].isin(EVALUATIONS)].copy()
    base["support_group"] = base["support_class"].astype(str).str.lower().replace({"nan": "support_unknown"})
    for (evaluation, model), sub in base.groupby(["evaluation", "model"], sort=True):
        total = len(sub)
        for support, group in sub.groupby("support_group", sort=True):
            metric = safe_metric_row(
                group["outcome_failure_48h_l2"].to_numpy(dtype=int),
                group["prediction"].to_numpy(dtype=float),
            )
            share = len(group) / total if total else np.nan
            if evaluation == "eicu_external" and support == "hfnc" and (len(group) < 500 or share < 0.05):
                flag = "eicu_hfnc_sparse_do_not_claim_exact_hfnc_transport"
            elif evaluation == "eicu_external":
                flag = "external_broader_respiratory_support_mapping"
            else:
                flag = "internal_mimic_support_mapping_check"
            metric.update(
                {
                    "evaluation": evaluation,
                    "model": model,
                    "mapping_dimension": "support_class_at_t0",
                    "mapping_level": support,
                    "support_share_within_model_evaluation": share,
                    "interpretation_flag": flag,
                }
            )
            rows.append(metric)
    return rows


def summarize_key_values(perf_rows: list[dict[str, object]], recal_rows: list[dict[str, object]], dca: list[dict[str, object]], endpoint_rows: list[dict[str, object]]) -> dict[str, object]:
    perf = pd.DataFrame(perf_rows)
    external = perf[(perf["evaluation"].eq("eicu_external")) & (perf["model"].eq("dynamic_last_no_support"))].iloc[0]
    rox_external = perf[(perf["evaluation"].eq("eicu_external")) & (perf["model"].eq("rox_only"))].iloc[0]
    mimic = perf[(perf["evaluation"].eq("mimic_internal_test")) & (perf["model"].eq("dynamic_last_no_support"))].iloc[0]
    recal = pd.DataFrame(recal_rows)
    original_ext = recal[
        recal["evaluation"].eq("eicu_external")
        & recal["model"].eq("dynamic_last_no_support")
        & recal["recalibration_method"].eq("original_mimic_trained")
    ].iloc[0]
    intercept_ext = recal[
        recal["evaluation"].eq("eicu_external")
        & recal["model"].eq("dynamic_last_no_support")
        & recal["recalibration_method"].eq("posthoc_intercept_only")
    ].iloc[0]
    dca_df = pd.DataFrame(dca)
    model_dca = dca_df[
        dca_df["evaluation"].eq("eicu_external")
        & dca_df["strategy"].eq("model")
        & dca_df["model"].isin(["dynamic_last_no_support", "rox_only"])
    ].copy()
    dyn_wins = 0
    checked = 0
    for threshold, group in model_dca.groupby("threshold"):
        if set(group["model"]) == {"dynamic_last_no_support", "rox_only"}:
            checked += 1
            values = dict(zip(group["model"], group["net_benefit"]))
            if values["dynamic_last_no_support"] > values["rox_only"]:
                dyn_wins += 1
    endpoint = pd.DataFrame(endpoint_rows)
    eicu_hfnc = endpoint[
        endpoint["evaluation"].eq("eicu_external")
        & endpoint["model"].eq("dynamic_last_no_support")
        & endpoint["mapping_level"].eq("hfnc")
    ]
    eicu_hfnc_n = int(eicu_hfnc.iloc[0]["n"]) if not eicu_hfnc.empty else 0
    eicu_hfnc_share = float(eicu_hfnc.iloc[0]["support_share_within_model_evaluation"]) if not eicu_hfnc.empty else 0.0
    return {
        "external_auroc": float(external["auroc"]),
        "external_auprc": float(external["auprc"]),
        "external_brier": float(external["brier"]),
        "external_calibration_slope": float(external["calibration_slope"]),
        "external_rox_auroc": float(rox_external["auroc"]),
        "mimic_holdout_auroc": float(mimic["auroc"]),
        "external_original_mean_predicted": float(original_ext["mean_predicted_risk"]),
        "external_original_event_rate": float(original_ext["event_rate"]),
        "external_original_brier": float(original_ext["brier"]),
        "external_intercept_brier": float(intercept_ext["brier"]),
        "external_intercept_alpha": float(intercept_ext["posthoc_intercept_alpha"]),
        "dca_dynamic_wins": dyn_wins,
        "dca_thresholds_checked": checked,
        "eicu_hfnc_n": eicu_hfnc_n,
        "eicu_hfnc_share": eicu_hfnc_share,
    }


def build_log(summary: dict[str, object]) -> None:
    lines = [
        "# ICU Translation Calibration and Decision-Curve Diagnostics V24.0",
        "",
        f"- Run date: {date.today().isoformat()}.",
        "- Scope: calibration bins, post-hoc calibration diagnostics, decision-curve source table, support-mapping sensitivity, and subgroup transportability checks.",
        "- Inputs: previously built local-only ICU modeling sets; outputs here are aggregate tables and text only.",
        "- Boundary: this ICU module tests a dynamic respiratory-support translation layer and does not directly validate the community PEF marker.",
        "",
        "## Key Results",
        "",
        f"- Best V23 dynamic model in MIMIC holdout: AUROC {summary['mimic_holdout_auroc']:.3f}.",
        f"- Best V23 dynamic model in eICU external validation: AUROC {summary['external_auroc']:.3f}, AUPRC {summary['external_auprc']:.3f}, Brier {summary['external_brier']:.3f}, calibration slope {summary['external_calibration_slope']:.3f}.",
        f"- ROX-like eICU comparator AUROC: {summary['external_rox_auroc']:.3f}.",
        f"- eICU event rate versus mean predicted risk for the dynamic model: {summary['external_original_event_rate']:.3f} versus {summary['external_original_mean_predicted']:.3f}.",
        f"- Post-hoc intercept-only recalibration changed dynamic-model eICU Brier from {summary['external_original_brier']:.3f} to {summary['external_intercept_brier']:.3f}; fitted intercept shift {summary['external_intercept_alpha']:.3f}.",
        f"- In common ROX-complete-case decision-curve rows, the dynamic model had higher net benefit than ROX-like comparator at {summary['dca_dynamic_wins']} of {summary['dca_thresholds_checked']} thresholds.",
        f"- eICU HFNC-labeled rows for the dynamic model: {summary['eicu_hfnc_n']} ({summary['eicu_hfnc_share']:.1%} of the eICU evaluation rows), so exact HFNC transport remains sparse.",
        "",
        "## Interpretation",
        "",
        "- The ICU translation layer has a reproducible external signal and acceptable calibration slope, but it remains a supportive translational analysis rather than the main evidence pillar.",
        "- eICU is best described as broader respiratory-support external validation because HFNC-specific labels are sparse.",
        "- Post-hoc recalibration rows are diagnostic only and should not be reported as independent validation performance.",
        "",
        "PASS: ICU translation calibration and DCA diagnostics V24.0 built.",
    ]
    (LOG_DIR / f"icu_translation_calibration_dca_{VERSION}.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def build_status_note(summary: dict[str, object]) -> None:
    lines = [
        "# ICU Translation Calibration Status Note V24.0",
        "",
        "V24.0 strengthens the ICU module by adding calibration-bin, decision-curve, support-mapping, and subgroup transportability diagnostics for the best V23 MIMIC-trained dynamic model and the ROX-like comparator.",
        "",
        f"The dynamic model retained moderate external discrimination in eICU (AUROC {summary['external_auroc']:.3f}; AUPRC {summary['external_auprc']:.3f}) and improved over the ROX-like comparator (AUROC {summary['external_rox_auroc']:.3f}). Calibration slope was close to 1 ({summary['external_calibration_slope']:.3f}), while mean predicted risk differed modestly from the observed event rate ({summary['external_original_mean_predicted']:.3f} versus {summary['external_original_event_rate']:.3f}).",
        "",
        "This is now a credible supportive ICU translation layer, but not a direct third-database validation of the community peak-flow marker and not a stand-alone clinical decision tool. The main manuscript should continue to place CHARLS and NHANES as the primary evidence, with MIMIC/eICU either in a translational supplement or a carefully bounded secondary section.",
    ]
    (MANUSCRIPT_DIR / f"icu_translation_calibration_status_note_{VERSION}.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    ensure_dirs()
    print("Loading V22 local-only ICU modeling sets...")
    mimic, eicu = load_local_sets()
    print("Refitting V24 dynamic and ROX-like models with the V22 MIMIC split...")
    predictions, perf_rows, coef_rows = fit_predictions(mimic, eicu)
    print("Building aggregate calibration, recalibration, DCA, and subgroup tables...")
    cal_bins = calibration_bin_rows(predictions)
    recal = recalibration_rows(predictions)
    dca = dca_rows(predictions)
    subgroup = subgroup_rows(predictions)
    endpoint = endpoint_mapping_rows(predictions)
    summary = summarize_key_values(perf_rows, recal, dca, endpoint)

    write_csv(TABLE_DIR / f"icu_translation_model_performance_{VERSION}.csv", perf_rows)
    write_csv(TABLE_DIR / f"icu_translation_model_coefficients_{VERSION}.csv", coef_rows)
    write_csv(TABLE_DIR / f"icu_translation_calibration_bins_{VERSION}.csv", cal_bins)
    write_csv(TABLE_DIR / f"icu_translation_recalibration_{VERSION}.csv", recal)
    write_csv(TABLE_DIR / f"icu_translation_dca_{VERSION}.csv", dca)
    write_csv(TABLE_DIR / f"icu_translation_subgroup_transport_{VERSION}.csv", subgroup)
    write_csv(TABLE_DIR / f"icu_translation_endpoint_mapping_sensitivity_{VERSION}.csv", endpoint)
    build_log(summary)
    build_status_note(summary)
    print("Wrote ICU translation calibration and DCA diagnostics V24.0 outputs.")


if __name__ == "__main__":
    main()
