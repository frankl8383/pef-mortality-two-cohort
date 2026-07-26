#!/usr/bin/env python3
"""Reproduce manuscript Tables 1-3 and Figures 1-4 from aggregate CSVs.

This script is intentionally limited to the disclosure-safe display layer. It
does not open participant-level data, run imputation, or refit cohort models.
Before rendering, it verifies the seven input files and their main numerical
and interpretation contracts.
"""

from __future__ import annotations

import hashlib
import shutil
from pathlib import Path

import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
from matplotlib.ticker import FuncFormatter, LogLocator, NullFormatter
from PIL import Image


ROOT = Path(__file__).resolve().parent
DATA_DIR = ROOT / "data"
OUT_DIR = ROOT / "output"
FIG_DIR = OUT_DIR / "figures"
TABLE_DIR = OUT_DIR / "tables"

INPUTS: dict[str, tuple[str, str, int, tuple[str, ...]]] = {
    "figure1": (
        "figure1_source_data.csv",
        "f38df6d2b7a612d42189ce11a80277e1fd312dc42f596fbf8cf75bba75b46a54",
        16,
        (
            "cohort",
            "step_order",
            "step_id",
            "step_label",
            "n",
            "deaths",
            "n_display",
            "deaths_display",
            "record_type",
        ),
    ),
    "figure2": (
        "figure2_source_data.csv",
        "29d9ba0f3fd526f05c4588067d474a16f7492d4553e56b322a73cbab0c2a037b",
        6,
        (
            "cohort",
            "model_tier",
            "label",
            "n",
            "deaths",
            "hr",
            "ci_low",
            "ci_high",
            "primary",
            "analysis_id",
            "scale",
        ),
    ),
    "figure3": (
        "figure3_source_data.csv",
        "a7c0b954f01f6492f9ee5db74e1485f7a82afde47e6d287b45d3c9da0f0a2501",
        8,
        (
            "panel",
            "analysis_id",
            "term",
            "label",
            "scale",
            "n",
            "deaths",
            "beta",
            "beta_ci_low",
            "beta_ci_high",
            "effect_measure",
            "effect",
            "effect_ci_low",
            "effect_ci_high",
            "p_value",
            "holm_p_value",
        ),
    ),
    "figure4": (
        "figure4_source_data.csv",
        "c33a131b80119b569f882ae487e4c6b0929041ae485629d1117cc61521b58ddc",
        7,
        (
            "panel",
            "cohort",
            "analysis_id",
            "term",
            "variant",
            "label",
            "scale",
            "model_n",
            "model_deaths",
            "weighting_target_n",
            "weighting_target_deaths",
            "effect",
            "effect_ci_low",
            "effect_ci_high",
            "p_value",
        ),
    ),
    "table1": (
        "table1_source_data.csv",
        "10c373a06c374fca42d047a49f49510e6ca553f0ddabc21d5aa201aa820fe389",
        43,
        (
            "cohort",
            "group_id",
            "group_label",
            "group_n_unweighted",
            "group_events_unweighted",
            "variable_order",
            "variable_id",
            "variable_label",
            "level_order",
            "level_id",
            "level_label",
            "variable_type",
            "display_role",
            "summary_rule",
            "estimate",
            "standard_error",
            "standard_deviation",
            "ci_low",
            "ci_high",
            "weighted_q25",
            "weighted_median",
            "weighted_q75",
            "unit",
            "n_observed_unweighted",
            "n_observed_display",
            "missing_n_unweighted",
            "missing_n_display",
            "display_value",
            "descriptive_source",
            "association_model_handling",
            "denominator_rule",
            "design_rule",
            "disclosure_status",
            "reliability_status",
        ),
    ),
    "table2": (
        "table2_source_data.csv",
        "710ca9ad7c986b6f1fded4dcc34d5e042a301d5b9b7b8cc068ea5d25f7dfd0b4",
        12,
        (
            "cohort",
            "scale",
            "model_tier",
            "n",
            "deaths",
            "hr",
            "ci_low",
            "ci_high",
            "analysis_id",
        ),
    ),
    "table3": (
        "table3_source_data.csv",
        "61e8880616232d9fb630f6705c888449744417acbb7a2e54ef35b9adec78a588",
        8,
        (
            "section",
            "analysis_id",
            "term",
            "label",
            "scale",
            "n",
            "deaths",
            "beta",
            "beta_ci_low",
            "beta_ci_high",
            "effect_measure",
            "effect",
            "effect_ci_low",
            "effect_ci_high",
            "p_value",
            "holm_p_value",
        ),
    ),
}

CHARLS = "#2F6B9A"
NHANES = "#C56A2D"
PEF = "#2F6B9A"
GLI = "#6C7784"
PRIMARY = "#243746"
GRID = "#D7DCE2"
TEXT = "#111111"
LIGHT_BLUE = "#EEF4F8"
LIGHT_ORANGE = "#FAF1E9"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def assert_close(observed: object, expected: float, label: str) -> None:
    if not np.isclose(float(observed), expected, rtol=0, atol=1e-12):
        raise RuntimeError(
            f"{label} mismatch: expected {expected:.15g}, "
            f"found {float(observed):.15g}"
        )


def one_row(
    frame: pd.DataFrame,
    *,
    analysis_id: str,
    term: str | None = None,
    variant: str | None = None,
) -> pd.Series:
    selected = frame.loc[frame["analysis_id"] == analysis_id]
    if term is not None:
        selected = selected.loc[selected["term"] == term]
    if variant is not None:
        selected = selected.loc[selected["variant"] == variant]
    if len(selected) != 1:
        raise RuntimeError(
            f"Expected one row for analysis_id={analysis_id}, term={term}, "
            f"variant={variant}; found {len(selected)}"
        )
    return selected.iloc[0]


def verify_input_files() -> None:
    for name, (filename, expected_hash, _, _) in INPUTS.items():
        path = DATA_DIR / filename
        if not path.is_file():
            raise FileNotFoundError(f"Missing aggregate input: data/{filename}")
        observed_hash = sha256(path)
        if observed_hash != expected_hash:
            raise RuntimeError(
                f"Input hash mismatch for {name}: expected {expected_hash}, "
                f"found {observed_hash}"
            )


def load_frames() -> dict[str, pd.DataFrame]:
    frames: dict[str, pd.DataFrame] = {}
    for name, (filename, _, _, _) in INPUTS.items():
        kwargs: dict[str, object] = {}
        if name == "table1":
            kwargs = {"dtype": str, "keep_default_na": False}
        frames[name] = pd.read_csv(DATA_DIR / filename, **kwargs)
    return frames


def validate_schema_and_privacy(frames: dict[str, pd.DataFrame]) -> None:
    prohibited_columns = {
        "participant_id",
        "person_id",
        "respondent_id",
        "household_id",
        "medical_record_number",
        "seqn",
    }
    for name, frame in frames.items():
        expected_rows = INPUTS[name][2]
        expected_columns = INPUTS[name][3]
        if len(frame) != expected_rows:
            raise RuntimeError(
                f"{name} must contain {expected_rows} aggregate rows; "
                f"found {len(frame)}"
            )
        if tuple(frame.columns) != expected_columns:
            raise RuntimeError(
                f"{name} schema mismatch.\nExpected: {expected_columns}\n"
                f"Found: {tuple(frame.columns)}"
            )
        column_hits = prohibited_columns.intersection(
            {str(column).casefold() for column in frame.columns}
        )
        if column_hits:
            raise RuntimeError(
                f"{name} contains participant identifier columns: "
                f"{sorted(column_hits)}"
            )
        for column in frame.columns:
            values = frame[column].dropna().astype(str)
            if values.str.contains(
                r"/Users/|/home/|[A-Za-z]:\\\\", regex=True
            ).any():
                raise RuntimeError(f"{name}.{column} contains a local path")


def validate_numeric_contracts(frames: dict[str, pd.DataFrame]) -> None:
    flow = frames["figure1"]
    if set(flow["record_type"]) != {
        "sequential_flow",
        "nested_analysis_sample",
    }:
        raise RuntimeError("Figure 1 record types are incomplete or mislabeled")
    sequential = flow.loc[flow["record_type"] == "sequential_flow"].copy()
    for cohort in ("CHARLS", "NHANES"):
        part = sequential.loc[sequential["cohort"] == cohort].sort_values(
            "step_order"
        )
        if not part["n"].astype(float).is_monotonic_decreasing:
            raise RuntimeError(f"{cohort} Figure 1 sample counts are not monotone")
        if not part["deaths"].astype(float).is_monotonic_decreasing:
            raise RuntimeError(f"{cohort} Figure 1 death counts are not monotone")
    final_flow = {
        ("CHARLS", "subsequent_vital_status_known"): (12555, 1735),
        ("NHANES", "pef_A_only"): (6719, 925),
    }
    for (cohort, step_id), (expected_n, expected_deaths) in final_flow.items():
        row = flow.loc[
            (flow["cohort"] == cohort) & (flow["step_id"] == step_id)
        ]
        if len(row) != 1:
            raise RuntimeError(f"Missing Figure 1 endpoint: {cohort}/{step_id}")
        assert_close(row.iloc[0]["n"], expected_n, f"{cohort} final n")
        assert_close(
            row.iloc[0]["deaths"], expected_deaths, f"{cohort} final deaths"
        )
    nested = {
        "complete_three_trials": ("12,512", "1,722"),
        "gli_common": ("6,540", "885"),
        "gli_reference_2022": ("5,304-5,309", "582-584"),
    }
    for step_id, expected in nested.items():
        row = flow.loc[flow["step_id"] == step_id]
        if len(row) != 1:
            raise RuntimeError(f"Missing Figure 1 nested sample: {step_id}")
        observed = (str(row.iloc[0]["n_display"]), str(row.iloc[0]["deaths_display"]))
        if observed != expected:
            raise RuntimeError(
                f"Nested sample mismatch for {step_id}: expected {expected}, "
                f"found {observed}"
            )

    figure2 = frames["figure2"]
    bool_values = (
        figure2["primary"]
        .astype(str)
        .str.strip()
        .str.casefold()
        .map({"true": True, "false": False})
    )
    if bool_values.isna().any():
        raise RuntimeError("Figure 2 contains an invalid primary indicator")
    figure2["primary"] = bool_values.astype(bool)
    for cohort in ("CHARLS", "NHANES"):
        part = figure2.loc[figure2["cohort"] == cohort]
        if set(part["model_tier"]) != {"A0", "A1", "A2"}:
            raise RuntimeError(f"Figure 2 adjustment tiers incomplete for {cohort}")
        if int(part["primary"].sum()) != 1:
            raise RuntimeError(f"Figure 2 must mark one A1 row for {cohort}")
    figure2_anchors = {
        ("CHARLS", "A1"): 1.4254815272918,
        ("NHANES", "A1"): 1.452707152858889,
    }
    for (cohort, tier), expected in figure2_anchors.items():
        row = figure2.loc[
            (figure2["cohort"] == cohort)
            & (figure2["model_tier"] == tier)
        ].iloc[0]
        assert_close(row["hr"], expected, f"Figure 2 {cohort} {tier} HR")

    for name in ("figure2", "figure3", "figure4", "table2", "table3"):
        frame = frames[name]
        if name in {"figure2", "table2"}:
            effect, low, high = "hr", "ci_low", "ci_high"
        else:
            effect, low, high = "effect", "effect_ci_low", "effect_ci_high"
        numeric = frame[[effect, low, high]].apply(
            pd.to_numeric, errors="raise"
        )
        if not np.isfinite(numeric.to_numpy()).all():
            raise RuntimeError(f"{name} contains non-finite effect estimates")
        if not (
            (numeric[low] <= numeric[effect])
            & (numeric[effect] <= numeric[high])
        ).all():
            raise RuntimeError(f"{name} contains an unordered confidence interval")

    figure3 = frames["figure3"]
    expected_figure3_ids = {
        "v45_g1",
        "v45_g2",
        "v45_g3",
        "v45_g4",
        "v45_g6a",
        "v45_g6b",
        "v45_gli_ref_2022",
    }
    if set(figure3["analysis_id"]) != expected_figure3_ids:
        raise RuntimeError("Figure 3 analysis set is incomplete")
    if (figure3["effect_measure"] != "HR").any():
        raise RuntimeError("Figure 3 must contain exposure HRs only")
    if set(figure3["panel"]) != {"A", "B"}:
        raise RuntimeError("Figure 3 panel labels are invalid")
    assert_close(
        one_row(figure3, analysis_id="v45_g6a")["effect"],
        1.28954518174239,
        "G6a conditional PEF HR",
    )
    ref_row = one_row(figure3, analysis_id="v45_gli_ref_2022")
    assert_close(ref_row["effect"], 1.25568516880495, "GLI reference-range PEF HR")
    if str(ref_row["n"]) != "5304-5309" or str(ref_row["deaths"]) != "582-584":
        raise RuntimeError("GLI reference-range denominator must remain a range")

    figure4 = frames["figure4"]
    expected_ipw = {
        "primary_original",
        "primary_uncapped",
        "primary_cap_p99",
        "primary_cap_p975",
    }
    observed_ipw = set(
        figure4.loc[figure4["cohort"] == "NHANES", "variant"]
    )
    if observed_ipw != expected_ipw:
        raise RuntimeError("Figure 4 NHANES observation-IPW variants are incomplete")
    charls_trial_ids = {
        "v45_c_trial_max_a1": 1.31945523717765,
        "v45_c_trial_mean_a1": 1.35827893527892,
        "v45_c_trial_median_a1": 1.34753895505856,
    }
    for analysis_id, expected in charls_trial_ids.items():
        row = one_row(figure4, analysis_id=analysis_id)
        assert_close(row["effect"], expected, f"{analysis_id} HR")
        if int(row["model_n"]) != 12512 or int(row["model_deaths"]) != 1722:
            raise RuntimeError(f"{analysis_id} denominator mismatch")
    original_ipw = one_row(
        figure4, analysis_id="v45_ipw_a1", variant="primary_original"
    )
    assert_close(
        original_ipw["effect"], 1.45529466911348, "Observation-IPW anchor HR"
    )
    if (
        int(original_ipw["weighting_target_n"]) != 7926
        or int(original_ipw["weighting_target_deaths"]) != 1228
    ):
        raise RuntimeError("Observation-IPW target denominator mismatch")

    table1 = frames["table1"]
    if set(table1["cohort"]) != {"CHARLS", "NHANES"}:
        raise RuntimeError("Table 1 must contain both cohorts")
    table1_totals = {
        "CHARLS": ("12555", "1735"),
        "NHANES": ("6719", "925"),
    }
    for cohort, expected in table1_totals.items():
        part = table1.loc[table1["cohort"] == cohort]
        observed_n = set(part["group_n_unweighted"])
        observed_deaths = set(part["group_events_unweighted"])
        if observed_n != {expected[0]} or observed_deaths != {expected[1]}:
            raise RuntimeError(f"Table 1 totals inconsistent for {cohort}")
    if not set(table1["reliability_status"]).issubset({"OK", "SUPPRESSED"}):
        raise RuntimeError("Table 1 contains an invalid reliability status")
    if int((table1["reliability_status"] == "SUPPRESSED").sum()) != 1:
        raise RuntimeError("Table 1 registered cell-suppression count changed")
    suppressed = table1["disclosure_status"].str.contains(
        "SUPPRESSION", regex=False
    )
    if not suppressed.any():
        raise RuntimeError("Table 1 disclosure-suppression contract is missing")
    if not (table1.loc[suppressed, "missing_n_display"] == "Suppressed").all():
        raise RuntimeError("Table 1 suppressed missingness cells were disclosed")

    table2 = frames["table2"]
    combinations = table2[
        ["cohort", "scale", "model_tier"]
    ].drop_duplicates()
    if len(combinations) != 12:
        raise RuntimeError("Table 2 cohort-scale-tier combinations are not unique")
    sd_rows = table2.loc[table2["scale"] == "1 sex-specific SD lower"]
    merged = figure2.merge(
        sd_rows,
        on=["cohort", "model_tier"],
        suffixes=("_figure", "_table"),
        validate="one_to_one",
    )
    if len(merged) != 6 or not np.allclose(
        merged["hr_figure"], merged["hr_table"], rtol=0, atol=1e-12
    ):
        raise RuntimeError("Figure 2 and Table 2 standardized HRs disagree")

    table3 = frames["table3"]
    expected_table3_ids = expected_figure3_ids | {
        "v45_paired_g6a_minus_g1"
    }
    if set(table3["analysis_id"]) != expected_table3_ids:
        raise RuntimeError("Table 3 analysis set is incomplete")
    paired = one_row(table3, analysis_id="v45_paired_g6a_minus_g1")
    if (
        paired["effect_measure"] != "Ratio of HRs"
        or str(paired["scale"]) != "Ratio of HRs"
    ):
        raise RuntimeError("Paired coefficient contrast is mislabeled")
    assert_close(paired["effect"], 0.883315828886649, "Paired ratio of HRs")
    if not (
        float(paired["effect_ci_low"])
        < 1.0
        < float(paired["effect_ci_high"])
    ):
        raise RuntimeError("Paired ratio-of-HRs CI must include no change")
    if float(paired["p_value"]) <= 0.05:
        raise RuntimeError("Paired coefficient contrast must not be labeled significant")
    for analysis_id in expected_figure3_ids:
        table_row = one_row(table3, analysis_id=analysis_id)
        figure_rows = figure3.loc[figure3["analysis_id"] == analysis_id]
        if not np.allclose(
            figure_rows["effect"].astype(float),
            float(table_row["effect"]),
            rtol=0,
            atol=1e-12,
        ):
            raise RuntimeError(
                f"Figure 3 and Table 3 disagree for {analysis_id}"
            )


def prepare_output() -> None:
    if OUT_DIR.exists():
        shutil.rmtree(OUT_DIR)
    FIG_DIR.mkdir(parents=True)
    TABLE_DIR.mkdir(parents=True)


def set_plot_style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": [
                "Arial",
                "Helvetica",
                "DejaVu Sans",
                "sans-serif",
            ],
            "font.size": 7.4,
            "axes.labelsize": 7.4,
            "axes.titlesize": 8.6,
            "xtick.labelsize": 6.9,
            "ytick.labelsize": 6.9,
            "axes.spines.right": False,
            "axes.spines.top": False,
            "axes.linewidth": 0.7,
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "text.color": TEXT,
            "axes.labelcolor": TEXT,
            "axes.edgecolor": TEXT,
            "xtick.color": TEXT,
            "ytick.color": TEXT,
        }
    )


def fmt_hr(effect: float, low: float, high: float) -> str:
    return f"{effect:.3f} ({low:.3f}-{high:.3f})"


def fmt_p(value: object) -> str:
    if value is None or pd.isna(value) or str(value).strip() == "":
        return ""
    numeric = float(value)
    return "<0.001" if numeric < 0.001 else f"{numeric:.3f}"


def count_display(value: object) -> str:
    text = str(value).strip()
    if not text or text.lower() == "nan":
        return ""
    try:
        return "--".join(f"{int(float(part)):,}" for part in text.split("-"))
    except ValueError:
        return text


def tex_escape(value: object) -> str:
    text = "" if value is None else str(value)
    for old, new in (
        ("\\", r"\textbackslash{}"),
        ("&", r"\&"),
        ("%", r"\%"),
        ("_", r"\_"),
        ("#", r"\#"),
    ):
        text = text.replace(old, new)
    return text


def save_figure(fig: plt.Figure, stem: str) -> list[Path]:
    outputs: list[Path] = []
    for suffix, kwargs in (
        (".svg", {}),
        (".pdf", {}),
        (".tiff", {"dpi": 600}),
        (".png", {"dpi": 300}),
    ):
        path = FIG_DIR / f"{stem}{suffix}"
        fig.savefig(path, facecolor="white", bbox_inches="tight", **kwargs)
        if suffix == ".tiff":
            with Image.open(path) as image:
                image.convert("RGB").save(
                    path,
                    dpi=(600, 600),
                    compression="tiff_lzw",
                )
        outputs.append(path)
    plt.close(fig)
    return outputs


def flow_box(
    ax: plt.Axes,
    x: float,
    y: float,
    width: float,
    height: float,
    label: str,
    count: str,
    color: str,
    *,
    primary: bool,
) -> None:
    face = color if primary else (LIGHT_BLUE if color == CHARLS else LIGHT_ORANGE)
    text_color = "white" if primary else TEXT
    patch = FancyBboxPatch(
        (x, y),
        width,
        height,
        boxstyle="round,pad=0.006,rounding_size=0.01",
        edgecolor=color,
        facecolor=face,
        linewidth=0.9 if primary else 0.7,
    )
    ax.add_patch(patch)
    ax.text(
        x + 0.025 * width,
        y + 0.60 * height,
        label,
        ha="left",
        va="center",
        fontsize=6.3,
        color=text_color,
        linespacing=1.1,
    )
    ax.text(
        x + 0.975 * width,
        y + 0.22 * height,
        count,
        ha="right",
        va="center",
        fontsize=6.5,
        fontweight="bold",
        color=text_color,
    )


def flow_arrow(
    ax: plt.Axes,
    x: float,
    y_top: float,
    y_bottom: float,
    color: str,
    *,
    dashed: bool = False,
) -> None:
    ax.add_patch(
        FancyArrowPatch(
            (x, y_top),
            (x, y_bottom),
            arrowstyle="-|>",
            mutation_scale=7,
            linewidth=0.7,
            color=color,
            linestyle="--" if dashed else "-",
        )
    )


def build_figure1(source: pd.DataFrame) -> list[Path]:
    sequential = source.loc[source["record_type"] == "sequential_flow"].copy()
    nested = source.loc[
        source["record_type"] == "nested_analysis_sample"
    ].copy()
    fig, axes = plt.subplots(1, 2, figsize=(6.69, 6.35))
    for panel, (ax, cohort, color) in enumerate(
        ((axes[0], "CHARLS", CHARLS), (axes[1], "NHANES", NHANES))
    ):
        frame = sequential.loc[sequential["cohort"] == cohort].sort_values(
            "step_order"
        )
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.axis("off")
        ax.text(
            0,
            1.01,
            chr(ord("a") + panel),
            fontsize=10,
            fontweight="bold",
            va="bottom",
        )
        ax.text(
            0.08,
            1.01,
            cohort,
            fontsize=9,
            color=color,
            fontweight="bold",
            va="bottom",
        )
        height = 0.077 if cohort == "NHANES" else 0.084
        ys = np.linspace(0.87, 0.35, len(frame))
        for index, (record, y) in enumerate(zip(frame.itertuples(), ys)):
            flow_box(
                ax,
                0.10,
                float(y),
                0.80,
                height,
                str(record.step_label),
                f"n={int(record.n):,}",
                color,
                primary=index == len(frame) - 1,
            )
            if index < len(frame) - 1:
                flow_arrow(
                    ax,
                    0.50,
                    float(y),
                    float(ys[index + 1] + height),
                    color,
                )
        final = frame.iloc[-1]
        ax.text(
            0.50,
            float(ys[-1] - 0.020),
            f"Deaths={int(final['deaths']):,}",
            ha="center",
            va="top",
            fontsize=6.4,
            color=color,
            fontweight="bold",
        )
        branches = nested.loc[nested["cohort"] == cohort].reset_index(drop=True)
        positions = (
            [(0.12, 0.12, 0.76), (0.12, 0.01, 0.76)]
            if len(branches) == 2
            else [(0.20, 0.065, 0.60)]
        )
        previous_bottom = float(ys[-1])
        for branch, (x, branch_y, branch_width) in zip(
            branches.itertuples(), positions
        ):
            flow_arrow(
                ax,
                0.50,
                previous_bottom,
                branch_y + 0.085,
                color,
                dashed=True,
            )
            label = str(branch.step_label)
            if len(label) > 35:
                label = label.replace(" within ", "\nwithin ")
            flow_box(
                ax,
                x,
                branch_y,
                branch_width,
                0.085,
                label,
                f"n={branch.n_display}; deaths={branch.deaths_display}",
                color,
                primary=False,
            )
            previous_bottom = branch_y
    fig.text(
        0.5,
        0.005,
        "Nested samples are shown beneath each original cohort; "
        "reference-range counts vary across imputations.",
        ha="center",
        fontsize=6.5,
        color="#4B5563",
    )
    fig.subplots_adjust(wspace=0.16, bottom=0.06, top=0.95)
    return save_figure(fig, "figure1_cohort_flow")


def forest_axis(
    ax: plt.Axes,
    labels: list[str],
    effects: np.ndarray,
    lows: np.ndarray,
    highs: np.ndarray,
    colors: list[str],
    title: str,
    xlabel: str,
    xlim: tuple[float, float],
) -> None:
    positions = np.arange(len(labels))[::-1]
    for y, effect, low, high, color in zip(
        positions, effects, lows, highs, colors
    ):
        ax.errorbar(
            effect,
            y,
            xerr=[[effect - low], [high - effect]],
            fmt="o",
            color=color,
            ecolor=color,
            markersize=4.1,
            elinewidth=1.0,
            capsize=2.2,
            markeredgecolor="white",
            markeredgewidth=0.35,
            zorder=3,
        )
        ax.text(
            1.025,
            y,
            fmt_hr(effect, low, high),
            ha="left",
            va="center",
            fontsize=6.5,
            clip_on=False,
            transform=ax.get_yaxis_transform(),
        )
    ax.axvline(1.0, color="#6B7280", linewidth=0.8, linestyle="--")
    ax.set_xscale("log")
    ax.set_xlim(*xlim)
    ax.set_ylim(-0.7, len(labels) - 0.3)
    ax.set_yticks(positions)
    ax.set_yticklabels(labels)
    ax.set_title(title, loc="left", fontweight="bold", pad=8)
    ax.set_xlabel(xlabel)
    ax.xaxis.set_major_locator(LogLocator(base=10, subs=(1.0, 1.25, 1.5, 2.0)))
    ax.xaxis.set_major_formatter(FuncFormatter(lambda value, _: f"{value:g}"))
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.grid(axis="x", color=GRID, linewidth=0.55)
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="y", length=0)


def build_figure2(source: pd.DataFrame) -> list[Path]:
    fig, axes = plt.subplots(2, 1, figsize=(6.69, 5.45), sharex=False)
    for panel, (ax, cohort, color) in enumerate(
        ((axes[0], "CHARLS", CHARLS), (axes[1], "NHANES", NHANES))
    ):
        subset = source.loc[source["cohort"] == cohort].copy()
        forest_axis(
            ax,
            subset["label"].tolist(),
            subset["hr"].to_numpy(float),
            subset["ci_low"].to_numpy(float),
            subset["ci_high"].to_numpy(float),
            [color] * len(subset),
            cohort,
            "Hazard ratio (95% CI)",
            (0.95, 2.00),
        )
        ax.text(
            -0.10,
            1.05,
            chr(ord("a") + panel),
            transform=ax.transAxes,
            fontsize=10,
            fontweight="bold",
        )
    fig.text(
        0.5,
        0.01,
        "Per one original cohort-specific, sex-specific SD lower PEF; "
        "cohorts were modeled separately.",
        ha="center",
        fontsize=6.5,
        color="#4B5563",
    )
    fig.subplots_adjust(
        left=0.36,
        right=0.72,
        bottom=0.14,
        top=0.93,
        hspace=0.58,
    )
    return save_figure(fig, "figure2_adjustment_ladder")


def build_figure3(source: pd.DataFrame) -> list[Path]:
    fig, axes = plt.subplots(2, 1, figsize=(6.69, 7.0))
    panels = (
        (
            "A",
            "Separate A1 models in the GLI common sample",
            [PEF, GLI, GLI, GLI],
            "Different PEF and GLI standardization systems",
        ),
        (
            "B",
            "PEF before and after spirometric adjustment",
            [PEF, PEF, PEF, PRIMARY],
            "Per one original sex-specific PEF SD",
        ),
    )
    for index, (panel, title, colors, xlabel) in enumerate(panels):
        subset = source.loc[source["panel"] == panel]
        forest_axis(
            axes[index],
            subset["label"].tolist(),
            subset["effect"].to_numpy(float),
            subset["effect_ci_low"].to_numpy(float),
            subset["effect_ci_high"].to_numpy(float),
            colors,
            title,
            xlabel,
            (0.92, 1.82),
        )
        axes[index].text(
            -0.10,
            1.05,
            panel.lower(),
            transform=axes[index].transAxes,
            fontsize=10,
            fontweight="bold",
        )
        axes[index].text(
            0.0,
            -0.24,
            (
                "n=6,540; deaths=885"
                if panel == "A"
                else "Reference-range n and deaths vary across imputations"
            ),
            transform=axes[index].transAxes,
            fontsize=6.2,
            color="#4B5563",
        )
    fig.subplots_adjust(
        left=0.42,
        right=0.72,
        bottom=0.10,
        top=0.94,
        hspace=0.60,
    )
    return save_figure(fig, "figure3_gli_comparison")


def build_figure4(source: pd.DataFrame) -> list[Path]:
    fig, axes = plt.subplots(2, 1, figsize=(6.69, 6.4))
    panel_contracts = (
        (
            "A",
            "NHANES strict-quality selection sensitivity",
            NHANES,
            "Per one original sex-specific PEF SD",
            "Model n=6,719; deaths=925; weighting target n=7,926",
        ),
        (
            "B",
            "CHARLS three-trial summary sensitivity",
            CHARLS,
            "Per 100 L/min lower PEF",
            "Complete-three sample n=12,512; deaths=1,722",
        ),
    )
    for index, (panel, title, color, xlabel, note) in enumerate(panel_contracts):
        subset = source.loc[source["panel"] == panel]
        forest_axis(
            axes[index],
            subset["label"].tolist(),
            subset["effect"].to_numpy(float),
            subset["effect_ci_low"].to_numpy(float),
            subset["effect_ci_high"].to_numpy(float),
            [color] * len(subset),
            title,
            xlabel,
            (0.95, 1.75),
        )
        axes[index].text(
            -0.10,
            1.05,
            panel.lower(),
            transform=axes[index].transAxes,
            fontsize=10,
            fontweight="bold",
        )
        axes[index].text(
            0.0,
            -0.24,
            note,
            transform=axes[index].transAxes,
            fontsize=6.2,
            color="#4B5563",
        )
    fig.subplots_adjust(
        left=0.40,
        right=0.72,
        bottom=0.11,
        top=0.94,
        hspace=0.60,
    )
    return save_figure(fig, "figure4_selection_measurement")


def table1_rows(table1: pd.DataFrame) -> list[str]:
    lines: list[str] = []
    for cohort, panel in (("CHARLS", "A"), ("NHANES", "B")):
        part = table1.loc[table1["cohort"] == cohort].copy()
        n = int(part["group_n_unweighted"].iloc[0])
        deaths = int(part["group_events_unweighted"].iloc[0])
        lines.append(
            rf"\multicolumn{{3}}{{l}}{{\textbf{{Panel {panel}. {cohort} "
            rf"(N={n:,}; deaths={deaths:,})}}}} \\"
        )
        previous_variable: str | None = None
        for row in part.itertuples(index=False):
            is_multilevel = str(row.variable_type) == "categorical"
            new_variable = str(row.variable_id) != previous_variable
            if is_multilevel:
                if new_variable:
                    characteristic = (
                        f"{row.variable_label} --- {row.level_label}"
                    )
                else:
                    characteristic = rf"\quad {tex_escape(row.level_label)}"
            else:
                characteristic = str(row.variable_label)
            characteristic = characteristic.replace(
                "Body-mass index", "Body mass index"
            )
            if not (is_multilevel and not new_variable):
                characteristic = tex_escape(characteristic)
            characteristic = characteristic.replace("kg/m2", r"kg/m$^2$")
            overall = tex_escape(row.display_value)
            missing = (
                tex_escape(row.missing_n_display)
                if new_variable or not is_multilevel
                else ""
            )
            lines.append(f"{characteristic} & {overall} & {missing} \\\\")
            previous_variable = str(row.variable_id)
        lines.append(r"\addlinespace")
    return lines


def write_table1(table1: pd.DataFrame) -> Path:
    lines = [
        r"\begin{longtable}{p{0.54\textwidth}p{0.25\textwidth}p{0.12\textwidth}}",
        r"\caption{Baseline characteristics of the CHARLS and NHANES primary mortality populations}\label{tab:baseline}\\",
        r"\toprule",
        r"Characteristic & Overall & Missing, n \\",
        r"\midrule",
        r"\endfirsthead",
        r"\multicolumn{3}{l}{\small\itshape Table \thetable\ continued} \\",
        r"\toprule",
        r"Characteristic & Overall & Missing, n \\",
        r"\midrule",
        r"\endhead",
        r"\midrule",
        r"\multicolumn{3}{r}{\small Continued on next page} \\",
        r"\endfoot",
        r"\bottomrule",
        r"\endlastfoot",
        *table1_rows(table1),
        r"\multicolumn{3}{p{0.95\textwidth}}{\footnotesize Values are design-weighted mean (SD), median [Q1, Q3], or unweighted n (weighted percentage). Percentages use nonmissing observations; summaries precede imputation. CHARLS used biomarker weights with community primary sampling units; NHANES used the Mobile Examination Center design with the strict PEF-quality sample as a domain. ``Suppressed'' denotes disclosure control. Cohort-specific covariates were not fully harmonized.} \\",
        r"\end{longtable}",
    ]
    path = TABLE_DIR / "table1_baseline.tex"
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def write_table2(table2: pd.DataFrame) -> Path:
    roles = {
        "A0": "Demographic/body-size",
        "A1": "Primary model",
        "A2": "Extended health/function",
    }
    lines = [
        r"\begin{longtable}{llllrrl}",
        r"\caption{Associations of lower baseline peak expiratory flow with all-cause mortality across adjustment tiers}\label{tab:primary}\\",
        r"\toprule",
        r"Cohort & Exposure scale & Tier & Role & n & Deaths & HR (95\% CI) \\",
        r"\midrule",
        r"\endfirsthead",
        r"\toprule",
        r"Cohort & Exposure scale & Tier & Role & n & Deaths & HR (95\% CI) \\",
        r"\midrule",
        r"\endhead",
        r"\bottomrule",
        r"\endfoot",
    ]
    for row in table2.itertuples(index=False):
        role = roles[str(row.model_tier)]
        result = fmt_hr(float(row.hr), float(row.ci_low), float(row.ci_high))
        if row.model_tier == "A1":
            tier_text = r"\textbf{A1}"
            role_text = rf"\textbf{{{role}}}"
            result_text = rf"\textbf{{{result}}}"
        else:
            tier_text = str(row.model_tier)
            role_text = role
            result_text = result
        lines.append(
            f"{row.cohort} & {tex_escape(row.scale)} & {tier_text} & "
            f"{role_text} & {int(row.n):,} & {int(row.deaths):,} & "
            f"{result_text} \\\\"
        )
    lines.extend(
        [
            r"\addlinespace",
            r"\multicolumn{7}{p{0.96\textwidth}}{\footnotesize HR, hazard ratio; CI, confidence interval. Standardized results are per one original cohort-specific, sex-specific, design-weighted SD lower raw PEF. CHARLS estimates are interval-hazard ratios from survey-weighted grouped-time complementary-log-log models. NHANES estimates are from survey-weighted Cox models. Unadjusted estimates were not used because age, sex, and body size are structural determinants of PEF; A0 was the narrowest specified model. A1 was the original primary model. A2 additionally included chronic disease and physical-function measures.} \\",
            r"\end{longtable}",
        ]
    )
    path = TABLE_DIR / "table2_primary_associations.tex"
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def write_table3(table3: pd.DataFrame) -> Path:
    focal_labels = {
        "E0": "Lower PEF",
        "L_FEV1": "Lower FEV1 z score",
        "L_FVC": "Lower FVC z score",
        "L_RATIO": "Lower FEV1/FVC z score",
        "G6a_minus_G1": "Conditional/separate PEF",
    }
    model_labels = {
        "v45_g1": "Separate PEF model",
        "v45_g2": "Separate FEV1 model",
        "v45_g3": "Separate FVC model",
        "v45_g4": "Separate FEV1/FVC model",
        "v45_g6a": "PEF + FEV1 + FEV1/FVC z scores",
        "v45_g6b": "PEF + FVC + FEV1/FVC z scores",
        "v45_gli_ref_2022": "PEF: all three GLI indices in range",
        "v45_paired_g6a_minus_g1": "Conditional versus separate PEF model",
    }
    lines = [
        r"\begin{longtable}{p{0.25\textwidth}p{0.25\textwidth}p{0.12\textwidth}rrp{0.15\textwidth}ll}",
        r"\caption{NHANES same-sample comparisons with GLI-standardized spirometric indices}\label{tab:gli}\\",
        r"\toprule",
        r"Model or contrast & Exposure or contrast & Scale & n & Deaths & Effect (95\% CI) & Nominal P & Holm P \\",
        r"\midrule",
        r"\endfirsthead",
        r"\toprule",
        r"Model or contrast & Exposure or contrast & Scale & n & Deaths & Effect (95\% CI) & Nominal P & Holm P \\",
        r"\midrule",
        r"\endhead",
        r"\bottomrule",
        r"\endfoot",
    ]
    for section, subset in table3.groupby("section", sort=False):
        lines.append(
            rf"\multicolumn{{8}}{{l}}{{\textbf{{{tex_escape(section)}}}}} \\"
        )
        for row in subset.itertuples(index=False):
            effect_label = (
                "HR" if row.effect_measure == "HR" else "Ratio of HRs"
            )
            estimate = (
                f"{effect_label} "
                f"{fmt_hr(float(row.effect), float(row.effect_ci_low), float(row.effect_ci_high))}"
            )
            lines.append(
                f"{tex_escape(model_labels[str(row.analysis_id)])} & "
                f"{tex_escape(focal_labels[str(row.term)])} & "
                f"{tex_escape(row.scale)} & {count_display(row.n)} & "
                f"{count_display(row.deaths)} & "
                f"{estimate} & {fmt_p(row.p_value)} & "
                f"{fmt_p(row.holm_p_value)} \\\\"
            )
        lines.append(r"\addlinespace")
    lines.extend(
        [
            r"\multicolumn{8}{p{0.98\textwidth}}{\footnotesize GLI, Global Lung Function Initiative; HR, hazard ratio; CI, confidence interval; PEF, peak expiratory flow. All models used NHANES survey weights and the A1 covariates. PEF effects are per one original sex-specific SD lower PEF; GLI effects are per one-unit lower GLI z score. Their point estimates therefore should not be interpreted as a ranking of measurements. The reference-range row required FEV1, FVC, and FEV1/FVC to be at or above the GLI Global 2022 lower limit in each completed dataset; n and deaths are ranges across imputations. The paired contrast is the conditional-model PEF HR divided by the separate-model PEF HR, not a new exposure HR. Holm P values apply only to the two focal PEF tests; all other P values are descriptive and unadjusted for multiplicity.} \\",
            r"\end{longtable}",
        ]
    )
    path = TABLE_DIR / "table3_nhanes_gli.tex"
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def validate_outputs(outputs: list[Path]) -> None:
    expected_names = {
        "table1_baseline.tex",
        "table2_primary_associations.tex",
        "table3_nhanes_gli.tex",
        *{
            f"{stem}{suffix}"
            for stem in (
                "figure1_cohort_flow",
                "figure2_adjustment_ladder",
                "figure3_gli_comparison",
                "figure4_selection_measurement",
            )
            for suffix in (".tiff", ".png", ".pdf", ".svg")
        },
    }
    observed_names = {path.name for path in outputs}
    if observed_names != expected_names:
        raise RuntimeError(
            f"Output set mismatch: expected {sorted(expected_names)}, "
            f"found {sorted(observed_names)}"
        )
    for path in outputs:
        if not path.is_file() or path.stat().st_size == 0:
            raise RuntimeError(f"Missing or empty output: {path}")
    for path in FIG_DIR.glob("*.tiff"):
        with Image.open(path) as image:
            if image.mode != "RGB":
                raise RuntimeError(f"{path.name} must be RGB; found {image.mode}")
            if min(image.size) < 2400:
                raise RuntimeError(
                    f"{path.name} is unexpectedly small: {image.size}"
                )
            compression = str(image.info.get("compression", "")).casefold()
            if compression not in {"tiff_lzw", "lzw"}:
                raise RuntimeError(
                    f"{path.name} must use LZW compression; found "
                    f"{compression or 'unknown'}"
                )
            dpi = image.info.get("dpi")
            if (
                not dpi
                or len(dpi) != 2
                or any(abs(float(value) - 600.0) > 1.0 for value in dpi)
            ):
                raise RuntimeError(
                    f"{path.name} must have 600-dpi metadata; found {dpi}"
                )
    table_text = {
        path.name: path.read_text(encoding="utf-8")
        for path in TABLE_DIR.glob("*.tex")
    }
    required_fragments = {
        "table1_baseline.tex": ("N=12,555", "N=6,719", "Suppressed"),
        "table2_primary_associations.tex": (
            "1.425 (1.313-1.548)",
            "1.453 (1.312-1.608)",
        ),
        "table3_nhanes_gli.tex": (
            "1.290 (1.066-1.560)",
            "0.883 (0.765-1.021)",
            "Ratio of HRs",
        ),
    }
    for filename, fragments in required_fragments.items():
        for fragment in fragments:
            if fragment not in table_text[filename]:
                raise RuntimeError(
                    f"{filename} is missing required content: {fragment}"
                )


def write_manifest(outputs: list[Path]) -> Path:
    records: list[dict[str, str]] = []
    for name, (filename, expected_hash, _, _) in INPUTS.items():
        path = DATA_DIR / filename
        records.append(
            {
                "asset_type": "input",
                "asset_id": name,
                "relative_path": str(path.relative_to(ROOT)),
                "sha256": sha256(path),
                "expected_sha256": expected_hash,
                "status": "VERIFIED",
            }
        )
    for path in sorted(outputs):
        records.append(
            {
                "asset_type": "output",
                "asset_id": path.stem,
                "relative_path": str(path.relative_to(ROOT)),
                "sha256": sha256(path),
                "expected_sha256": "",
                "status": "GENERATED_AND_CHECKED",
            }
        )
    manifest = OUT_DIR / "manifest.csv"
    pd.DataFrame(records).to_csv(manifest, index=False)
    return manifest


def main() -> None:
    verify_input_files()
    frames = load_frames()
    validate_schema_and_privacy(frames)
    validate_numeric_contracts(frames)
    prepare_output()
    set_plot_style()

    outputs: list[Path] = []
    outputs.extend(build_figure1(frames["figure1"]))
    outputs.extend(build_figure2(frames["figure2"]))
    outputs.extend(build_figure3(frames["figure3"]))
    outputs.extend(build_figure4(frames["figure4"]))
    outputs.append(write_table1(frames["table1"]))
    outputs.append(write_table2(frames["table2"]))
    outputs.append(write_table3(frames["table3"]))
    validate_outputs(outputs)
    manifest = write_manifest(outputs)
    print(
        "REPRODUCTION PASS: verified 7 aggregate inputs, generated and "
        f"checked {len(outputs)} display files, and wrote "
        f"{manifest.relative_to(ROOT)}"
    )


if __name__ == "__main__":
    main()
