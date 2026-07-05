#!/usr/bin/env python3
"""Build a safe handoff packet for MIMIC-IV V19/V19.1 feasibility work."""

from __future__ import annotations

import csv
import re
import zipfile
from datetime import date
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[2]
TABLE_DIR = PROJECT_ROOT / "results" / "tables"
LOG_DIR = PROJECT_ROOT / "results" / "logs"
MANUSCRIPT_DIR = PROJECT_ROOT / "manuscript"
HANDOFF_DIR = PROJECT_ROOT / "handoff"


FILES = [
    "results/logs/mimiciv_download_completion_monitor_v19_0.md",
    "results/tables/mimiciv_download_completion_status_v19_0.csv",
    "results/tables/mimiciv_download_completion_monitor_v19_0.csv",
    "results/logs/mimiciv_integrity_schema_feasibility_v0_1.md",
    "results/logs/mimiciv_integrity_schema_feasibility_v0_1_validation.md",
    "results/logs/mimiciv_schema_cohort_endpoint_lock_v19_0.md",
    "results/logs/mimiciv_schema_cohort_endpoint_lock_v19_0_validation.md",
    "results/tables/mimiciv_key_variable_lock_v19_0.csv",
    "results/tables/mimiciv_cohort_feasibility_counts_v19_0.csv",
    "results/tables/mimiciv_respiratory_item_candidates_v19_0.csv",
    "results/tables/mimiciv_icd_respiratory_code_candidates_v19_0.csv",
    "results/tables/mimiciv_endpoint_feasibility_v19_0.csv",
    "manuscript/mimiciv_integration_decision_memo_v19_0.md",
    "results/tables/mimiciv_manuscript_integration_decision_v19_0.csv",
    "results/logs/mimiciv_integration_decision_memo_v19_0_validation.md",
    "manuscript/mimiciv_minimal_model_design_v19_1.md",
    "results/tables/mimiciv_minimal_model_design_v19_1.csv",
    "results/logs/mimiciv_minimal_model_design_v19_1_validation.md",
    "manuscript/target_journal_strategy_mimic_update_v19_1.md",
    "results/tables/target_journal_strategy_mimic_update_v19_1.csv",
    "results/logs/target_journal_strategy_mimic_update_v19_1_validation.md",
]


def write_csv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def sensitive_hits(path: Path) -> list[str]:
    text = path.read_text(errors="ignore")
    markers = [r"${SECURE_DATA_ROOT}", r"password", r"token", r"lzhfrankl", r"PHYSIONET_USER"]
    return [marker for marker in markers if re.search(marker, text, flags=re.IGNORECASE)]


def main() -> None:
    HANDOFF_DIR.mkdir(parents=True, exist_ok=True)
    file_rows = []
    for rel in FILES:
        path = PROJECT_ROOT / rel
        if not path.exists():
            raise SystemExit(f"Missing handoff file: {rel}")
        hits = sensitive_hits(path) if path.suffix.lower() in {".md", ".csv", ".txt"} else []
        file_rows.append(
            {
                "file": rel,
                "exists": "yes",
                "size_bytes": path.stat().st_size,
                "sensitive_scan": "pass" if not hits else "fail:" + ";".join(hits),
            }
        )
    write_csv(TABLE_DIR / "mimiciv_v19_1_handoff_file_index.csv", file_rows, ["file", "exists", "size_bytes", "sensitive_scan"])

    handoff = [
        "# Morning Handoff V19.1",
        "",
        "## 一句话结论",
        "",
        "MIMIC-IV 下载和完整性门已经通过，schema/cohort/respiratory endpoint feasibility 也通过；MIMIC 有 ICU charted Peak Exp Flow Rate 候选项，但尚不能等同于 CHARLS/NHANES 标准化 PEF，因此当前不应写成三库验证，最佳定位是可选 ICU translation layer。",
        "",
        "## 已完成",
        "",
        "- V19.0 download completion monitor: 33/33 files pass, 0 missing, 0 gzip failures, 0 SHA failures.",
        "- V0.1 MIMIC integrity/schema feasibility rerun: PASS in completed-download state.",
        "- V19.0 schema/cohort/endpoint lock: adult ICU denominator, linkage, mortality, LOS and respiratory support endpoint metadata are feasible.",
        "- V19.0 integration decision memo: keep current main evidence as CHARLS + NHANES; do not add MIMIC to abstract/Results/conclusion yet.",
        "- V19.1 minimal model design: design-only, no model executed; run only after endpoint/time-window lock.",
        "- V19.1 target-journal addendum: Respiratory Research remains the default balanced first route; MIMIC feasibility alone does not change the route.",
        "",
        "## 当前最应该阅读的 5 个文件",
        "",
        "1. `manuscript/mimiciv_integration_decision_memo_v19_0.md`",
        "2. `results/logs/mimiciv_schema_cohort_endpoint_lock_v19_0.md`",
        "3. `manuscript/mimiciv_minimal_model_design_v19_1.md`",
        "4. `manuscript/target_journal_strategy_mimic_update_v19_1.md`",
        "5. `results/tables/mimiciv_endpoint_feasibility_v19_0.csv`",
        "",
        "## 下一步",
        "",
        "- 若继续 MIMIC：先做 final endpoint lock table，明确 item IDs、ICD code sets、baseline window、outcome window 和 baseline ventilated exclusion。",
        "- 若准备投稿：暂时仍按 CHARLS + NHANES 主稿推进，MIMIC 只写为 pending/optional translation layer，除非后续模型通过并科学连贯。",
        "- 作者仍需填写 V18 author response worksheet 中的 12 个 P0 作者信息/声明项。",
        "",
        "PASS: Morning handoff V19.1 build complete.",
    ]
    (LOG_DIR / "morning_handoff_v19_1.md").write_text("\n".join(handoff) + "\n")

    zip_path = HANDOFF_DIR / "respiratory_vulnerability_mimic_v19_1.zip"
    with zipfile.ZipFile(zip_path, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        for rel in FILES:
            zf.write(PROJECT_ROOT / rel, arcname=rel)
        zf.write(LOG_DIR / "morning_handoff_v19_1.md", arcname="results/logs/morning_handoff_v19_1.md")
        zf.write(TABLE_DIR / "mimiciv_v19_1_handoff_file_index.csv", arcname="results/tables/mimiciv_v19_1_handoff_file_index.csv")

    failures = [row for row in file_rows if not str(row["sensitive_scan"]).startswith("pass")]
    if failures:
        raise SystemExit("Sensitive scan failed for handoff files.")

    validation = [
        "# MIMIC-IV V19.1 Handoff Validation",
        "",
        f"- Run date: {date.today().isoformat()}.",
        f"- Files indexed: {len(file_rows)}.",
        f"- ZIP: `handoff/{zip_path.name}`.",
        "- Sensitive-marker scan passed for indexed text files.",
        "- No raw MIMIC data files are included.",
        "",
        "PASS: MIMIC-IV V19.1 handoff validation complete.",
    ]
    (LOG_DIR / "mimiciv_v19_1_handoff_validation.md").write_text("\n".join(validation) + "\n")
    print(f"Wrote handoff/{zip_path.name}")


if __name__ == "__main__":
    main()
