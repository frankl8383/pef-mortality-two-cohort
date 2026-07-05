#!/usr/bin/env python3
"""Compare V23 ICU translation model variants using V22 local model sets."""

from __future__ import annotations

import csv
import importlib.util
from datetime import date
from pathlib import Path

import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[2]
V22_PATH = PROJECT_ROOT / "python" / "icu_translation" / "build_icu_translation_modeling_v22_0.py"
SENSITIVE_DIR = PROJECT_ROOT / "derived_sensitive" / "icu_translation"
TABLE_DIR = PROJECT_ROOT / "results" / "tables"
LOG_DIR = PROJECT_ROOT / "results" / "logs"
MANUSCRIPT_DIR = PROJECT_ROOT / "manuscript"
VERSION = "v23_0"


def load_v22_module():
    spec = importlib.util.spec_from_file_location("icu_v22", V22_PATH)
    if spec is None or spec.loader is None:
        raise RuntimeError("Could not load V22 modeling module.")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


v22 = load_v22_module()

MODEL_VARIANTS = {
    "clinical_basic": ["age", "male", "support_hfnc"],
    "rox_only": ["rox_l2"],
    "dynamic_last": [
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
    "dynamic_last_no_support": [
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
    "dynamic_worst": [
        "age",
        "male",
        "support_hfnc",
        "spo2_min",
        "fio2_max",
        "respiratory_rate_max",
        "heart_rate_max",
        "oxygen_flow_max",
        "peep_max",
        "rox_l2",
    ],
    "dynamic_worst_no_support": [
        "age",
        "male",
        "spo2_min",
        "fio2_max",
        "respiratory_rate_max",
        "heart_rate_max",
        "oxygen_flow_max",
        "peep_max",
        "rox_l2",
    ],
    "dynamic_vitals_worst": [
        "age",
        "male",
        "spo2_min",
        "fio2_max",
        "respiratory_rate_max",
        "heart_rate_max",
        "rox_l2",
    ],
}


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
    return mimic, eicu


def comparison_rows(perf_rows: list[dict[str, object]]) -> list[dict[str, object]]:
    perf = pd.DataFrame(perf_rows)
    rows: list[dict[str, object]] = []
    for evaluation in ["mimic_internal_test", "eicu_external"]:
        sub = perf[perf["evaluation"].eq(evaluation)].copy()
        rox = sub[sub["model"].eq("rox_only")]
        if rox.empty:
            continue
        rox_row = rox.iloc[0]
        for _, row in sub.iterrows():
            rows.append(
                {
                    "evaluation": evaluation,
                    "model": row["model"],
                    "n": int(row["n"]),
                    "events": int(row["events"]),
                    "auroc": float(row["auroc"]),
                    "auroc_delta_vs_rox": float(row["auroc"]) - float(rox_row["auroc"]),
                    "auprc": float(row["auprc"]),
                    "auprc_delta_vs_rox": float(row["auprc"]) - float(rox_row["auprc"]),
                    "brier": float(row["brier"]),
                    "brier_delta_vs_rox": float(row["brier"]) - float(rox_row["brier"]),
                    "calibration_slope": float(row["calibration_slope"]),
                }
            )
    return rows


def build_log(perf_rows: list[dict[str, object]], comp_rows: list[dict[str, object]]) -> None:
    perf = pd.DataFrame(perf_rows)
    external = perf[perf["evaluation"].eq("eicu_external")].copy().sort_values("auroc", ascending=False)
    internal = perf[perf["evaluation"].eq("mimic_internal_test")].copy().sort_values("auroc", ascending=False)
    lines = [
        "# ICU Translation Model Variants V23.0",
        "",
        f"- Run date: {date.today().isoformat()}.",
        "- Scope: fast model-variant comparison using V22 local row-level modeling sets.",
        "- No raw MIMIC/eICU tables were rescanned.",
        "- Boundary: model-selection support only; not a final clinical model.",
        "",
        "## Best External AUROC",
        "",
    ]
    for _, row in external.head(5).iterrows():
        lines.append(
            f"- {row['model']}: eICU AUROC {float(row['auroc']):.3f}, "
            f"AUPRC {float(row['auprc']):.3f}, Brier {float(row['brier']):.3f}, "
            f"calibration slope {float(row['calibration_slope']):.3f}."
        )
    lines.extend(["", "## Best Internal AUROC", ""])
    for _, row in internal.head(5).iterrows():
        lines.append(
            f"- {row['model']}: MIMIC holdout AUROC {float(row['auroc']):.3f}, "
            f"AUPRC {float(row['auprc']):.3f}, Brier {float(row['brier']):.3f}, "
            f"calibration slope {float(row['calibration_slope']):.3f}."
        )
    best = external.iloc[0]
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            f"- Best current external variant: {best['model']} with eICU AUROC {float(best['auroc']):.3f}.",
            "- If no variant clearly exceeds the ROX-like comparator externally, the ICU module should be described as promising but still under refinement.",
            "",
            "PASS: ICU translation model variants V23.0 built.",
        ]
    )
    (LOG_DIR / f"icu_translation_model_variants_{VERSION}.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def build_status_note(perf_rows: list[dict[str, object]]) -> None:
    perf = pd.DataFrame(perf_rows)
    external = perf[perf["evaluation"].eq("eicu_external")].copy().sort_values("auroc", ascending=False)
    best = external.iloc[0]
    lines = [
        "# ICU Translation Model Variant Status Note V23.0",
        "",
        f"The strongest current external-validation variant is `{best['model']}` with eICU AUROC {float(best['auroc']):.3f}, AUPRC {float(best['auprc']):.3f}, Brier {float(best['brier']):.3f}, and calibration slope {float(best['calibration_slope']):.3f}.",
        "",
        "This improves the ICU translation layer from feasibility-only to a real, externally tested signal, but it remains below the cleaner high-impact threshold that would make the ICU module a major centerpiece without further calibration and endpoint-mapping review.",
    ]
    (MANUSCRIPT_DIR / f"icu_translation_model_variant_status_note_{VERSION}.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    MANUSCRIPT_DIR.mkdir(parents=True, exist_ok=True)
    mimic, eicu = load_local_sets()
    original_sets = v22.FEATURE_SETS
    try:
        v22.FEATURE_SETS = MODEL_VARIANTS
        perf_rows, coef_rows, _, _ = v22.fit_and_evaluate_models(mimic, eicu)
    finally:
        v22.FEATURE_SETS = original_sets
    comp_rows = comparison_rows(perf_rows)
    write_csv(TABLE_DIR / f"icu_translation_model_variant_performance_{VERSION}.csv", perf_rows)
    write_csv(TABLE_DIR / f"icu_translation_model_variant_coefficients_{VERSION}.csv", coef_rows)
    write_csv(TABLE_DIR / f"icu_translation_model_variant_comparison_{VERSION}.csv", comp_rows)
    build_log(perf_rows, comp_rows)
    build_status_note(perf_rows)
    print("Wrote ICU translation model variants V23.0 outputs.")


if __name__ == "__main__":
    main()
