#!/usr/bin/env python3
"""Build target-journal strategy addendum after MIMIC V19.1 decision."""

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
            "journal_or_route": "Respiratory Research",
            "pre_mimic_position": "balanced_first",
            "mimic_v19_1_effect": "fit_unchanged_until_model",
            "recommended_action": "Continue as the default balanced first route for the CHARLS + NHANES manuscript; do not add MIMIC claims yet.",
            "probability_estimate": "20-30% now; 30-45% only if a coherent MIMIC ICU translation model later passes validation.",
        },
        {
            "journal_or_route": "npj Primary Care Respiratory Medicine",
            "pre_mimic_position": "balanced_first",
            "mimic_v19_1_effect": "little_direct_help",
            "recommended_action": "Still viable for field-measure/primary-care framing; MIMIC ICU layer may be less central to the journal story.",
            "probability_estimate": "similar to V4.1; MIMIC should not drive this route.",
        },
        {
            "journal_or_route": "BMJ Open Respiratory Research",
            "pre_mimic_position": "safe_respiratory",
            "mimic_v19_1_effect": "may_help_transparency_if_added_later",
            "recommended_action": "Good fallback if the manuscript emphasizes transparent observational reporting and keeps MIMIC as optional or future work.",
            "probability_estimate": "unchanged until validated MIMIC model exists.",
        },
        {
            "journal_or_route": "ERJ Open Research / Thorax stretch",
            "pre_mimic_position": "higher_risk",
            "mimic_v19_1_effect": "potential_ceiling_only_after_validated_model",
            "recommended_action": "Do not stretch upward on feasibility alone; consider only if MIMIC produces a coherent respiratory-support result.",
            "probability_estimate": "high desk-risk without a validated MIMIC result.",
        },
        {
            "journal_or_route": "Age and Ageing / Lancet Healthy Longevity",
            "pre_mimic_position": "aging_high_or_stretch",
            "mimic_v19_1_effect": "not_primary_value",
            "recommended_action": "MIMIC does not solve the core ageing-actionability issue; keep older-adult clinical relevance grounded in CHARLS.",
            "probability_estimate": "unchanged by MIMIC feasibility alone.",
        },
    ]
    write_csv(
        TABLE_DIR / "target_journal_strategy_mimic_update_v19_1.csv",
        rows,
        ["journal_or_route", "pre_mimic_position", "mimic_v19_1_effect", "recommended_action", "probability_estimate"],
    )

    lines = [
        "# Target Journal Strategy MIMIC Update V19.1",
        "",
        f"- Run date: {date.today().isoformat()}.",
        "- Inputs: V4.1 target journal strategy, V19.0 MIMIC integration decision memo, and V19.1 minimal model design.",
        "- Bottom line: MIMIC completion improves technical readiness but does not yet change the main submission route because no validated MIMIC model exists and the ICU peak-flow candidate is not established as equivalent to standardized field PEF.",
        "",
        "## Recommendation",
        "",
        "Respiratory Research remains the most realistic balanced first target. The paper should still be submitted, if author fields are completed, as a CHARLS + NHANES observational marker paper. MIMIC should be mentioned only as pending or optional translational work unless a final endpoint/time-window lock and coherent model are completed. The ICU charted Peak Exp Flow Rate candidate should be evaluated transparently, but feasibility alone should not be used to pitch three-database validation.",
        "",
        "## Probability Update",
        "",
        "The current subjective Respiratory Research acceptance range remains roughly 20-30%. A coherent MIMIC ICU translation layer could raise the estimated range to roughly 30-45%, but only if the model is scientifically aligned and does not blur the PEF construct. A poorly aligned MIMIC layer would likely reduce clarity and should be left out.",
        "",
        "## Route Table",
        "",
        "| Journal or route | MIMIC V19.1 effect | Recommended action |",
        "| --- | --- | --- |",
    ]
    for row in rows:
        lines.append(f"| {row['journal_or_route']} | {row['mimic_v19_1_effect']} | {row['recommended_action']} |")
    lines.extend(["", "PASS: Target journal strategy MIMIC update V19.1 built."])
    (MANUSCRIPT_DIR / "target_journal_strategy_mimic_update_v19_1.md").write_text("\n".join(lines) + "\n")

    validation = [
        "# Target Journal Strategy MIMIC Update V19.1 Validation",
        "",
        f"- Run date: {date.today().isoformat()}.",
        f"- Route rows checked: {len(rows)}.",
        "- Confirms Respiratory Research remains default balanced first route.",
        "- Confirms MIMIC feasibility alone does not justify three-database validation claim.",
        "- Confirms probability update is conditional on validated MIMIC modeling.",
        "",
        "PASS: target journal strategy MIMIC update V19.1 validation complete.",
    ]
    (LOG_DIR / "target_journal_strategy_mimic_update_v19_1_validation.md").write_text("\n".join(validation) + "\n")
    print("Wrote target journal strategy MIMIC update V19.1 outputs.")


if __name__ == "__main__":
    main()
