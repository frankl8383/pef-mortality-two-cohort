#!/usr/bin/env python3
"""Build manuscript displays from disclosure-safe aggregate-only sources.

This script never opens participant-level data. It verifies the SHA-256 of
every aggregate input before producing source-data CSVs, LaTeX tables, figures,
and a build manifest.
"""

from __future__ import annotations

import csv
import hashlib
import math
import textwrap
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from matplotlib.ticker import FuncFormatter, LogLocator, NullFormatter
from PIL import Image


PROJECT = Path(__file__).resolve().parents[1]
OUT = PROJECT / "rebuild"
FIG_DIR = OUT / "figures"
TABLE_DIR = OUT / "tables"
SOURCE_DIR = OUT / "source_data"

SOURCES = {
    "charls_registry": (
        PROJECT
        / "results/v44_1_1/stage3_charls_primary_association_registry_v44_1_1.csv",
        "ac05ec1528c413a6355684e8faca408d12bce393f1d5fb4d3748bf86495d00bd",
    ),
    "charls_anthro": (
        PROJECT
        / "results/v44_1_1/stage3_charls_anthropometric_sensitivity_registry_v44_1_1.csv",
        "f123e879a27d832c8b28a5f6407766bbe1927bf0234ed6f2481606f0fd0b5d41",
    ),
    "nhanes_registry": (
        PROJECT / "results/v43/stage3_nhanes_primary_association_registry_v43.csv",
        "e6257e8483e416ad43b39624f431bbdbda141186eb4e87e00bcb4f68374e1b76",
    ),
    "charls_flow": (
        PROJECT / "results/v43/charls_sequential_flow_v43.csv",
        "61f8013b4139f10e1a40181bd45eb3192a4763dfff81a174c61a81aa843439e9",
    ),
    "nhanes_flow": (
        PROJECT / "results/v43/nhanes_cohort_flow_v43.csv",
        "8bd7b6d7531d9ed74efab3068a9efaaf3840bc2adfb868b86b843957b10b3843",
    ),
    "table1": (
        PROJECT / "results/v44_2/table1_baseline_characteristics_source_v44_2.csv",
        "2172d7cdaf5dfc57191a89a250b4892959ec81f9d81af860640210c157b3f4c4",
    ),
}

CHARLS = "#2F6B9A"
NHANES = "#D17B32"
PRIMARY_DARK = "#243746"
SUPPORT = "#8C98A4"
GRID = "#D8DDE3"
TEXT = "#1F2933"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def verify_sources() -> None:
    for name, (path, expected) in SOURCES.items():
        if not path.exists():
            raise FileNotFoundError(f"Missing frozen aggregate source: {name}")
        observed = sha256(path)
        if observed != expected:
            raise RuntimeError(
                f"Frozen aggregate hash mismatch for {name}: "
                f"expected {expected}, observed {observed}"
            )


def prepare_dirs() -> None:
    for path in (FIG_DIR, TABLE_DIR, SOURCE_DIR):
        path.mkdir(parents=True, exist_ok=True)


def row_by_id(frame: pd.DataFrame, analysis_id: str) -> pd.Series:
    rows = frame.loc[frame["analysis_id"] == analysis_id]
    if len(rows) != 1:
        raise RuntimeError(
            f"Expected exactly one aggregate row for {analysis_id}; found {len(rows)}"
        )
    return rows.iloc[0]


def hr_fields(row: pd.Series, nhanes: bool = False) -> tuple[float, float, float]:
    if nhanes:
        return (
            float(row["exp_estimate"]),
            float(row["exp_ci_low"]),
            float(row["exp_ci_high"]),
        )
    return (
        float(row["hazard_ratio"]),
        float(row["hr_ci_low"]),
        float(row["hr_ci_high"]),
    )


def fmt_int(value: float | int | str) -> str:
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return "--"
    return f"{int(float(value)):,}"


def fmt_hr(hr: float, low: float, high: float) -> str:
    return f"{hr:.3f} ({low:.3f}--{high:.3f})"


def tex_escape(value: object) -> str:
    text = "" if value is None else str(value)
    replacements = [
        ("\\", r"\textbackslash{}"),
        ("&", r"\&"),
        ("%", r"\%"),
        ("$", r"\$"),
        ("#", r"\#"),
        ("_", r"\_"),
        ("{", r"\{"),
        ("}", r"\}"),
        ("~", r"\textasciitilde{}"),
        ("^", r"\textasciicircum{}"),
    ]
    for old, new in replacements:
        text = text.replace(old, new)
    text = text.replace("≥", r"$\geq$")
    text = text.replace("–", "--")
    text = text.replace("²", r"$^2$")
    return text


def write_source_csv(frame: pd.DataFrame, filename: str) -> Path:
    path = SOURCE_DIR / filename
    frame.to_csv(path, index=False, quoting=csv.QUOTE_MINIMAL)
    return path


def save_figure(fig: plt.Figure, stem: str) -> list[Path]:
    outputs = []
    for suffix, kwargs in (
        (".svg", {}),
        (".pdf", {}),
        (".tiff", {"dpi": 600}),
        (".png", {"dpi": 300}),
    ):
        path = FIG_DIR / f"{stem}{suffix}"
        fig.savefig(path, facecolor="white", **kwargs)
        if suffix == ".tiff":
            with Image.open(path) as raster:
                raster.convert("RGB").save(
                    path,
                    dpi=(600, 600),
                    compression="tiff_lzw",
                )
        outputs.append(path)
    plt.close(fig)
    return outputs


def set_plot_style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "font.size": 7.2,
            "axes.labelsize": 7.2,
            "axes.titlesize": 8.4,
            "xtick.labelsize": 6.8,
            "ytick.labelsize": 6.8,
            "axes.spines.right": False,
            "axes.spines.top": False,
            "axes.linewidth": 0.7,
            "legend.frameon": False,
            "text.color": TEXT,
            "axes.labelcolor": TEXT,
            "axes.edgecolor": TEXT,
            "xtick.color": TEXT,
            "ytick.color": TEXT,
        }
    )


def draw_flow_box(
    ax: plt.Axes,
    x: float,
    y: float,
    width: float,
    height: float,
    label: str,
    n: int,
    color: str,
    excluded: int | None = None,
    primary: bool = False,
) -> None:
    face = "#EFF5FA" if color == CHARLS else "#FBF2E9"
    if primary:
        face = color
        txt_color = "white"
        edge = color
    else:
        txt_color = TEXT
        edge = color
    box = FancyBboxPatch(
        (x, y),
        width,
        height,
        boxstyle="round,pad=0.008,rounding_size=0.012",
        linewidth=1.0 if primary else 0.75,
        edgecolor=edge,
        facecolor=face,
    )
    ax.add_patch(box)
    wrapped = "\n".join(textwrap.wrap(label, width=36))
    ax.text(
        x + width * 0.04,
        y + height * 0.58,
        wrapped,
        ha="left",
        va="center",
        fontsize=6.3,
        color=txt_color,
        linespacing=1.05,
    )
    ax.text(
        x + width * 0.96,
        y + height * 0.21,
        f"n={n:,}",
        ha="right",
        va="center",
        fontsize=6.5,
        fontweight="bold",
        color=txt_color,
    )
    if excluded is not None and excluded > 0:
        ax.text(
            x - 0.018,
            y + height * 0.5,
            f"−{excluded:,}",
            ha="right",
            va="center",
            fontsize=5.8,
            color="#6B7280",
        )


def draw_vertical_arrow(
    ax: plt.Axes,
    x: float,
    y_top: float,
    y_bottom: float,
    color: str,
    dashed: bool = False,
) -> None:
    arrow = FancyArrowPatch(
        (x, y_top),
        (x, y_bottom),
        arrowstyle="-|>",
        mutation_scale=7,
        linewidth=0.75,
        color=color,
        linestyle="--" if dashed else "-",
    )
    ax.add_patch(arrow)


def build_figure1(
    charls_flow: pd.DataFrame,
    nhanes_flow: pd.DataFrame,
    charls_registry: pd.DataFrame,
    nhanes_registry: pd.DataFrame,
) -> tuple[pd.DataFrame, list[Path]]:
    c = charls_flow.loc[
        charls_flow["flow_branch"] == "primary_P_C_upper_bound"
    ].sort_values("step_order").copy()
    c.loc[c.index[0], "step_label"] = "Harmonized CHARLS core records"
    n_core = nhanes_flow.loc[nhanes_flow["flow_branch"] == "core"].copy()
    n_quality = nhanes_flow.loc[nhanes_flow["flow_branch"] == "pef_quality"].copy()
    n_selected = pd.concat(
        [
            n_core.loc[n_core["step_id"].isin(["all_demo", "age_45_79", "age_45_79_mec"])],
            n_core.loc[n_core["step_id"] == "target_linked"],
            n_core.loc[n_core["step_id"] == "safe_target"],
            n_quality.loc[n_quality["step_id"] == "pef_nonmissing"],
            n_quality.loc[n_quality["step_id"] == "pef_A_only"],
        ],
        ignore_index=True,
    )
    n_selected["step_order"] = range(1, len(n_selected) + 1)
    n_selected["excluded_from_previous"] = (
        n_selected["n"].shift(1) - n_selected["n"]
    )
    n_selected.loc[n_selected.index[0], "excluded_from_previous"] = np.nan
    n_selected.loc[
        n_selected["step_id"] == "age_45_79_mec", "step_label"
    ] = "MEC examined and mortality-file matched"
    n_selected.loc[
        n_selected["step_id"] == "all_demo", "step_label"
    ] = "All NHANES 2007–2012 demographic records"
    n_selected.loc[
        n_selected["step_id"] == "target_linked", "step_label"
    ] = "Mortality-linkage eligible (ELIGSTAT=1)"
    n_selected.loc[
        n_selected["step_id"] == "safe_target", "step_label"
    ] = "Spirometry status available and no official safety exclusion"
    n_selected.loc[
        n_selected["step_id"] == "pef_A_only", "step_label"
    ] = "Primary sample: grade-A spirometry"

    c_cc = row_by_id(charls_registry, "v44_s3_charls_pc_e0_a1_complete_case")
    c_land = row_by_id(
        charls_registry, "v44_s3_charls_pc_e0_a1_w2_landmark_conditional_survivor"
    )
    n_cc = row_by_id(nhanes_registry, "v43_s3_nhanes_pn_e0_a1_complete_case")
    n_land = row_by_id(
        nhanes_registry, "v43_s3_nhanes_pn_e0_two_year_landmark_a1"
    )
    branch_rows = pd.DataFrame(
        [
            {
                "cohort": "CHARLS",
                "branch": "Primary multiple-imputation analysis",
                "n": 12555,
                "deaths": 1735,
                "branch_type": "primary",
            },
            {
                "cohort": "CHARLS",
                "branch": "Complete-case sensitivity",
                "n": int(c_cc["n"]),
                "deaths": int(c_cc["events"]),
                "branch_type": "sensitivity",
            },
            {
                "cohort": "CHARLS",
                "branch": "Wave 2 conditional-survivor sensitivity",
                "n": int(c_land["n"]),
                "deaths": int(c_land["events"]),
                "branch_type": "sensitivity",
            },
            {
                "cohort": "NHANES",
                "branch": "Primary multiple-imputation analysis",
                "n": 6719,
                "deaths": 925,
                "branch_type": "primary",
            },
            {
                "cohort": "NHANES",
                "branch": "Complete-case sensitivity",
                "n": int(n_cc["n"]),
                "deaths": int(n_cc["events"]),
                "branch_type": "sensitivity",
            },
            {
                "cohort": "NHANES",
                "branch": "Two-year conditional-survivor sensitivity",
                "n": int(n_land["n"]),
                "deaths": int(n_land["events"]),
                "branch_type": "sensitivity",
            },
            {
                "cohort": "NHANES",
                "branch": "Quality A plus at least 3 acceptable curves",
                "n": 6634,
                "deaths": 905,
                "branch_type": "quality",
            },
            {
                "cohort": "NHANES",
                "branch": "Quality A or C",
                "n": 6971,
                "deaths": 967,
                "branch_type": "quality",
            },
            {
                "cohort": "NHANES",
                "branch": "Quality A/B/C",
                "n": 7083,
                "deaths": 989,
                "branch_type": "quality",
            },
        ]
    )
    flow_source = pd.concat(
        [
            c.assign(cohort="CHARLS")[
                [
                    "cohort",
                    "step_order",
                    "step_id",
                    "step_label",
                    "n",
                    "excluded_from_previous",
                    "deaths",
                ]
            ],
            n_selected.assign(cohort="NHANES")[
                [
                    "cohort",
                    "step_order",
                    "step_id",
                    "step_label",
                    "n",
                    "excluded_from_previous",
                    "deaths",
                ]
            ],
        ],
        ignore_index=True,
    )
    flow_source["record_type"] = "sequential_flow"
    branch_export = branch_rows.copy()
    branch_export["record_type"] = "analysis_branch"
    source = pd.concat([flow_source, branch_export], ignore_index=True, sort=False)
    write_source_csv(source, "figure1_source_data.csv")

    fig, axes = plt.subplots(1, 2, figsize=(6.69, 6.35))
    panels = [
        (axes[0], c, "a", "CHARLS", CHARLS),
        (axes[1], n_selected, "b", "NHANES", NHANES),
    ]
    for ax, frame, panel, title, color in panels:
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.axis("off")
        ax.text(
            0.0,
            1.01,
            panel,
            transform=ax.transAxes,
            fontsize=10,
            fontweight="bold",
            va="bottom",
        )
        ax.text(
            0.08,
            1.01,
            title,
            transform=ax.transAxes,
            fontsize=9,
            fontweight="bold",
            color=color,
            va="bottom",
        )
        x, width, height = 0.12, 0.78, 0.082
        ys = np.linspace(0.88, 0.30, len(frame))
        rows = list(frame.itertuples(index=False))
        for idx, (row, y) in enumerate(zip(rows, ys)):
            draw_flow_box(
                ax,
                x,
                y,
                width,
                height,
                str(row.step_label),
                int(row.n),
                color,
                excluded=(
                    None
                    if pd.isna(row.excluded_from_previous)
                    else int(row.excluded_from_previous)
                ),
                primary=idx == len(rows) - 1,
            )
            if idx < len(rows) - 1:
                next_y = ys[idx + 1] + height
                draw_vertical_arrow(
                    ax,
                    x + width / 2,
                    y,
                    next_y,
                    color,
                )
        terminal = rows[-1]
        ax.text(
            x + width / 2,
            ys[-1] - 0.026,
            f"Primary mortality population: n={int(terminal.n):,}; "
            f"deaths={int(terminal.deaths):,}",
            ha="center",
            va="top",
            fontsize=6.5,
            fontweight="bold",
            color=color,
        )
        cohort_branches = branch_rows.loc[
            (branch_rows["cohort"] == title)
            & (branch_rows["branch_type"] != "quality")
        ]
        quality_branches = branch_rows.loc[
            (branch_rows["cohort"] == title)
            & (branch_rows["branch_type"] == "quality")
        ]
        branch_y = 0.145
        n_branches = len(cohort_branches)
        branch_width = 0.27
        gap = (0.96 - n_branches * branch_width) / max(n_branches - 1, 1)
        branch_xs = [0.02 + i * (branch_width + gap) for i in range(n_branches)]
        draw_vertical_arrow(
            ax,
            x + width / 2,
            ys[-1],
            branch_y + 0.095,
            color,
            dashed=True,
        )
        for bx, branch in zip(branch_xs, cohort_branches.itertuples(index=False)):
            face = "#FFFFFF"
            edge = color if branch.branch_type == "primary" else SUPPORT
            patch = FancyBboxPatch(
                (bx, branch_y),
                branch_width,
                0.092,
                boxstyle="round,pad=0.006,rounding_size=0.01",
                linewidth=0.75,
                edgecolor=edge,
                facecolor=face,
                linestyle="-" if branch.branch_type == "primary" else "--",
            )
            ax.add_patch(patch)
            label = "\n".join(textwrap.wrap(branch.branch, width=21))
            ax.text(
                bx + branch_width / 2,
                branch_y + 0.050,
                label,
                ha="center",
                va="center",
                fontsize=5.4,
                linespacing=1.0,
            )
            ax.text(
                bx + branch_width / 2,
                branch_y + 0.012,
                f"n={branch.n:,}; deaths={branch.deaths:,}",
                ha="center",
                va="center",
                fontsize=5.3,
                color=edge,
                fontweight="bold" if branch.branch_type == "primary" else "normal",
            )
        if title == "NHANES" and len(quality_branches):
            quality_y = 0.012
            quality_box = FancyBboxPatch(
                (0.05, quality_y),
                0.90,
                0.082,
                boxstyle="round,pad=0.006,rounding_size=0.01",
                linewidth=0.75,
                edgecolor=SUPPORT,
                facecolor="white",
                linestyle="--",
            )
            ax.add_patch(quality_box)
            quality_text = "  |  ".join(
                [
                    f"{row.branch.replace('Quality ', '')}: n={row.n:,}"
                    for row in quality_branches.itertuples(index=False)
                ]
            )
            ax.text(
                0.50,
                quality_y + 0.052,
                "Alternative spirometry-quality samples",
                ha="center",
                va="center",
                fontsize=5.7,
                fontweight="bold",
                color=NHANES,
            )
            ax.text(
                0.50,
                quality_y + 0.026,
                quality_text,
                ha="center",
                va="center",
                fontsize=5.0,
                color="#6B7280",
            )
    fig.subplots_adjust(left=0.02, right=0.98, bottom=0.02, top=0.96, wspace=0.12)
    return source, save_figure(fig, "figure1_cohort_flow")


def forest_panel(
    ax: plt.Axes,
    frame: pd.DataFrame,
    cohort: str,
    color: str,
    label_col: str,
    highlight_col: str,
    xlim: tuple[float, float],
) -> None:
    frame = frame.reset_index(drop=True)
    y = np.arange(len(frame))[::-1]
    ax.axvline(1.0, color="#6B7280", linewidth=0.8, linestyle="--", zorder=0)
    for idx, row in frame.iterrows():
        yi = y[idx]
        primary = bool(row[highlight_col])
        point_color = color if primary else SUPPORT
        marker = "s" if primary else "o"
        ax.errorbar(
            row["hr"],
            yi,
            xerr=np.array(
                [[row["hr"] - row["ci_low"]], [row["ci_high"] - row["hr"]]]
            ),
            fmt=marker,
            markersize=5.2 if primary else 4.0,
            markerfacecolor=point_color,
            markeredgecolor=point_color,
            ecolor=point_color,
            elinewidth=1.2 if primary else 0.9,
            capsize=2.2,
            capthick=0.9,
            zorder=3,
        )
        weight = "bold" if primary else "normal"
        ax.text(
            -0.02,
            yi,
            str(row[label_col]),
            transform=ax.get_yaxis_transform(),
            ha="right",
            va="center",
            fontsize=6.6,
            fontweight=weight,
            clip_on=False,
        )
        ax.text(
            1.02,
            yi,
            f"{int(row['n']):,}/{int(row['deaths']):,}   "
            f"{row['hr']:.3f} ({row['ci_low']:.3f}–{row['ci_high']:.3f})",
            transform=ax.get_yaxis_transform(),
            ha="left",
            va="center",
            fontsize=6.2,
            fontweight=weight,
            clip_on=False,
        )
    ax.set_xscale("log")
    ax.set_xlim(*xlim)
    ax.set_ylim(-0.7, len(frame) - 0.3)
    ax.set_yticks([])
    ax.grid(axis="x", color=GRID, linewidth=0.55)
    ax.xaxis.set_major_locator(LogLocator(base=10, subs=(1.0, 1.2, 1.4, 1.6, 1.8, 2.0)))
    ax.xaxis.set_major_formatter(
        FuncFormatter(lambda value, _: f"{value:.1f}" if 0.9 <= value <= 2.1 else "")
    )
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.set_title(cohort, loc="left", color=color, fontweight="bold", pad=6)
    ax.spines["left"].set_visible(False)
    ax.spines["bottom"].set_color("#9AA3AC")
    ax.text(
        1.02,
        len(frame) - 0.03,
        "n/deaths   HR (95% CI)",
        transform=ax.get_yaxis_transform(),
        ha="left",
        va="bottom",
        fontsize=5.9,
        color="#6B7280",
        clip_on=False,
    )


def build_figure2(
    charls_registry: pd.DataFrame, nhanes_registry: pd.DataFrame
) -> tuple[pd.DataFrame, list[Path]]:
    records = []
    for cohort, registry, prefix, nhanes, color in (
        ("CHARLS", charls_registry, "v44_s3_charls_pc_e0_", False, CHARLS),
        ("NHANES", nhanes_registry, "v43_s3_nhanes_pn_e0_", True, NHANES),
    ):
        for tier, label in (
            ("a0", "A0: demographic/body-size"),
            ("a1", "A1: primary model"),
            ("a2", "A2: extended health/function"),
        ):
            row = row_by_id(registry, f"{prefix}{tier}")
            hr, low, high = hr_fields(row, nhanes=nhanes)
            records.append(
                {
                    "cohort": cohort,
                    "model_tier": tier.upper(),
                    "label": label,
                    "n": int(row["n"]),
                    "deaths": int(row["events"]),
                    "hr": hr,
                    "ci_low": low,
                    "ci_high": high,
                    "primary": tier == "a1",
                    "analysis_id": row["analysis_id"],
                    "cohort_color": color,
                }
            )
    source = pd.DataFrame(records)
    write_source_csv(source.drop(columns="cohort_color"), "figure2_source_data.csv")

    fig, axes = plt.subplots(2, 1, figsize=(6.69, 4.45), sharex=True)
    for idx, (ax, cohort, color) in enumerate(
        [(axes[0], "CHARLS", CHARLS), (axes[1], "NHANES", NHANES)]
    ):
        forest_panel(
            ax,
            source.loc[source["cohort"] == cohort],
            cohort,
            color,
            "label",
            "primary",
            (0.92, 2.02),
        )
        ax.text(
            -0.42,
            1.08,
            chr(ord("a") + idx),
            transform=ax.transAxes,
            fontsize=10,
            fontweight="bold",
            va="top",
            clip_on=False,
        )
    axes[1].set_xlabel("Hazard ratio for all-cause mortality (log scale)")
    fig.subplots_adjust(left=0.33, right=0.67, top=0.94, bottom=0.13, hspace=0.40)
    return source, save_figure(fig, "figure2_adjustment_ladder")


def build_figure3(
    charls_registry: pd.DataFrame,
    charls_anthro: pd.DataFrame,
    nhanes_registry: pd.DataFrame,
) -> tuple[pd.DataFrame, list[Path]]:
    specs = [
        (
            "CHARLS",
            charls_registry,
            "v44_s3_charls_pc_e0_a1",
            "Primary: WHO-block QC; PEF boundaries retained",
            False,
            True,
        ),
        (
            "CHARLS",
            charls_registry,
            "v44_s3_charls_pc_e0_a1_boundary_excluded",
            "PEF boundary values 30/890 L/min excluded",
            False,
            False,
        ),
        (
            "CHARLS",
            charls_registry,
            "v44_s3_charls_pc_e0_a1_strict_measurement",
            "Standing, full effort, and at least 2 valid trials",
            False,
            False,
        ),
        (
            "CHARLS",
            charls_anthro,
            "v44_s3_charls_pc_e0_a1_who_component_no_later_aux",
            "WHO component-only anthropometric invalidation",
            False,
            False,
        ),
        (
            "CHARLS",
            charls_anthro,
            "v44_s3_charls_pc_e0_a1_pcornet_block_no_later_aux",
            "Stricter PCORnet-style block bounds",
            False,
            False,
        ),
        (
            "NHANES",
            nhanes_registry,
            "v43_s3_nhanes_pn_e0_a1",
            "Grade A (primary sample)",
            True,
            True,
        ),
        (
            "NHANES",
            nhanes_registry,
            "v43_s3_nhanes_pn_e0_pef_A_acc3_a1",
            "Grade A with at least 3 acceptable curves",
            True,
            False,
        ),
        (
            "NHANES",
            nhanes_registry,
            "v43_s3_nhanes_pn_e0_pef_A_plus_C_a1",
            "Grades A or C",
            True,
            False,
        ),
        (
            "NHANES",
            nhanes_registry,
            "v43_s3_nhanes_pn_e0_pef_ABC_a1",
            "Grades A, B, or C",
            True,
            False,
        ),
        (
            "NHANES",
            nhanes_registry,
            "v43_s3_nhanes_pn_e0_gt900_influence_exclusion_a1",
            "Values >900 L/min excluded, influence diagnostic",
            True,
            False,
        ),
    ]
    records = []
    for cohort, registry, aid, label, nhanes, primary in specs:
        row = row_by_id(registry, aid)
        hr, low, high = hr_fields(row, nhanes=nhanes)
        records.append(
            {
                "cohort": cohort,
                "analysis": label,
                "n": int(row["n"]),
                "deaths": int(row["events"]),
                "hr": hr,
                "ci_low": low,
                "ci_high": high,
                "primary": primary,
                "analysis_id": aid,
                "scale": "1 sex-specific within-cohort weighted SD lower raw PEF",
            }
        )
    source = pd.DataFrame(records)
    write_source_csv(source, "figure3_source_data.csv")
    fig, axes = plt.subplots(2, 1, figsize=(6.69, 6.20), sharex=True)
    for idx, (ax, cohort, color) in enumerate(
        [(axes[0], "CHARLS", CHARLS), (axes[1], "NHANES", NHANES)]
    ):
        forest_panel(
            ax,
            source.loc[source["cohort"] == cohort],
            cohort,
            color,
            "analysis",
            "primary",
            (0.92, 1.82),
        )
        ax.text(
            -0.60,
            1.06,
            chr(ord("a") + idx),
            transform=ax.transAxes,
            fontsize=10,
            fontweight="bold",
            va="top",
            clip_on=False,
        )
    axes[1].set_xlabel("Hazard ratio per sex-specific SD lower PEF (log scale)")
    fig.subplots_adjust(left=0.43, right=0.66, top=0.96, bottom=0.10, hspace=0.35)
    return source, save_figure(fig, "figure3_measurement_robustness")


def table1_characteristic_rows(table1: pd.DataFrame) -> list[str]:
    lines: list[str] = []
    for cohort, panel in (("CHARLS", "A"), ("NHANES", "B")):
        part = table1.loc[table1["cohort"] == cohort].copy()
        if cohort == "CHARLS":
            n, deaths = 12555, 1735
        else:
            n, deaths = 6719, 925
        lines.append(
            rf"\multicolumn{{3}}{{l}}{{\textbf{{Panel {panel}. {cohort} "
            rf"(N={n:,}; deaths={deaths:,})}}}} \\"
        )
        previous_variable = None
        for row in part.itertuples(index=False):
            is_multilevel = str(row.variable_type) == "categorical"
            new_variable = row.variable_id != previous_variable
            if is_multilevel:
                if new_variable:
                    characteristic = f"{row.variable_label} --- {row.level_label}"
                else:
                    characteristic = tex_escape(str(row.level_label))
                    characteristic = rf"\quad {characteristic}"
            else:
                characteristic = str(row.variable_label)
            characteristic = characteristic.replace("Body-mass index", "Body mass index")
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
            previous_variable = row.variable_id
        lines.append(r"\addlinespace")
    return lines


def write_table1(table1: pd.DataFrame) -> Path:
    if len(table1) != 43:
        raise RuntimeError(f"Table 1 aggregate source must contain 43 rows; found {len(table1)}")
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
    ]
    lines.extend(table1_characteristic_rows(table1))
    lines.extend(
        [
            r"\multicolumn{3}{p{0.95\textwidth}}{\footnotesize Values are design-weighted mean (weighted SD), weighted median [Q1, Q3], or unweighted n (design-weighted percentage), as appropriate. Percentages use nonmissing observations for the variable. Descriptions are observed, pre-imputation summaries. CHARLS estimates use the biomarker weight and community primary sampling unit; NHANES estimates use the Mobile Examination Center design with the grade-A spirometry sample specified as a domain. Cells marked ``Suppressed'' follow the study disclosure rules described in the Supplementary Methods. Cohort-specific disease, socioeconomic, depressive-symptom, and function variables were not fully harmonized.} \\",
            r"\end{longtable}",
        ]
    )
    path = TABLE_DIR / "table1_baseline.tex"
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    table1.to_csv(SOURCE_DIR / "table1_source_data.csv", index=False)
    return path


def table2_records(
    charls_registry: pd.DataFrame, nhanes_registry: pd.DataFrame
) -> pd.DataFrame:
    records = []
    for cohort, registry, prefix, nhanes in (
        ("CHARLS", charls_registry, "v44_s3_charls_pc_", False),
        ("NHANES", nhanes_registry, "v43_s3_nhanes_pn_", True),
    ):
        for exposure, scale in (
            ("e0", "1 sex-specific SD lower"),
            ("e0b", "100 L/min lower"),
        ):
            for tier in ("a0", "a1", "a2"):
                row = row_by_id(registry, f"{prefix}{exposure}_{tier}")
                hr, low, high = hr_fields(row, nhanes=nhanes)
                records.append(
                    {
                        "cohort": cohort,
                        "scale": scale,
                        "model_tier": tier.upper(),
                        "n": int(row["n"]),
                        "deaths": int(row["events"]),
                        "hr": hr,
                        "ci_low": low,
                        "ci_high": high,
                        "analysis_id": row["analysis_id"],
                    }
                )
    return pd.DataFrame(records)


def write_table2(table2: pd.DataFrame) -> Path:
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
    roles = {
        "A0": "Demographic/body-size",
        "A1": "Primary model",
        "A2": "Extended health/function",
    }
    for row in table2.itertuples(index=False):
        role = roles[row.model_tier]
        if row.model_tier == "A1":
            tier = r"\textbf{A1}"
            role = rf"\textbf{{{role}}}"
            result = rf"\textbf{{{fmt_hr(row.hr, row.ci_low, row.ci_high)}}}"
        else:
            tier = row.model_tier
            result = fmt_hr(row.hr, row.ci_low, row.ci_high)
        lines.append(
            f"{row.cohort} & {tex_escape(row.scale)} & {tier} & {role} & "
            f"{row.n:,} & {row.deaths:,} & {result} \\\\"
        )
    lines.extend(
        [
            r"\addlinespace",
            r"\multicolumn{7}{p{0.96\textwidth}}{\footnotesize HR, hazard ratio; CI, confidence interval. Standardized results are per one cohort-specific, sex-specific, design-weighted SD lower raw PEF. CHARLS estimates are interval-hazard ratios from survey-weighted grouped-time complementary-log-log models. NHANES estimates are from survey-weighted Cox models. A1 was the primary model. A2 additionally included chronic disease and physical-function measures.} \\",
            r"\end{longtable}",
        ]
    )
    path = TABLE_DIR / "table2_primary_associations.tex"
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    table2.to_csv(SOURCE_DIR / "table2_source_data.csv", index=False)
    return path


def build_table3_records(
    charls_registry: pd.DataFrame,
    charls_anthro: pd.DataFrame,
    nhanes_registry: pd.DataFrame,
) -> pd.DataFrame:
    specifications = [
        ("CHARLS", charls_registry, "v44_s3_charls_pc_e0_a1", "Primary model", "Multiple imputation", False),
        ("CHARLS", charls_registry, "v44_s3_charls_pc_e0_a1_complete_case", "Complete case", "Missing data", False),
        ("CHARLS", charls_registry, "v44_s3_charls_pc_e0_a1_w2_landmark_conditional_survivor", "Wave 2 conditional-survivor", "Early follow-up", False),
        ("CHARLS", charls_registry, "v44_s3_charls_pc_e0_a1_boundary_excluded", "PEF boundary values excluded", "PEF measurement", False),
        ("CHARLS", charls_registry, "v44_s3_charls_pc_e0_a1_strict_measurement", "Strict PEF measurement", "PEF measurement", False),
        ("CHARLS", charls_anthro, "v44_s3_charls_pc_e0_a1_who_component_no_later_aux", "WHO component-only routing", "Anthropometry", False),
        ("CHARLS", charls_anthro, "v44_s3_charls_pc_e0_a1_pcornet_block_no_later_aux", "PCORnet-style block bounds", "Anthropometry", False),
        ("NHANES", nhanes_registry, "v43_s3_nhanes_pn_e0_a1", "Primary model", "Multiple imputation", True),
        ("NHANES", nhanes_registry, "v43_s3_nhanes_pn_e0_a1_complete_case", "Complete case", "Missing data", True),
        ("NHANES", nhanes_registry, "v43_s3_nhanes_pn_e0_two_year_landmark_a1", "Two-year conditional-survivor", "Early follow-up", True),
        ("NHANES", nhanes_registry, "v43_s3_nhanes_pn_e0_pef_A_acc3_a1", "Grade A with at least 3 acceptable curves", "Spirometry quality", True),
        ("NHANES", nhanes_registry, "v43_s3_nhanes_pn_e0_pef_A_plus_C_a1", "Grades A or C", "Spirometry quality", True),
        ("NHANES", nhanes_registry, "v43_s3_nhanes_pn_e0_pef_ABC_a1", "Grades A, B, or C", "Spirometry quality", True),
        ("NHANES", nhanes_registry, "v43_s3_nhanes_pn_e0_gt900_influence_exclusion_a1", "Values >900 L/min excluded", "High PEF values", True),
        ("NHANES", nhanes_registry, "v43_s3_nhanes_pn_e0_leave_cycle_E_a1", "Cycle E excluded", "Survey cycle", True),
        ("NHANES", nhanes_registry, "v43_s3_nhanes_pn_e0_leave_cycle_F_a1", "Cycle F excluded", "Survey cycle", True),
        ("NHANES", nhanes_registry, "v43_s3_nhanes_pn_e0_leave_cycle_G_a1", "Cycle G excluded", "Survey cycle", True),
    ]
    records = []
    for cohort, registry, aid, analysis, assumption, nhanes in specifications:
        row = row_by_id(registry, aid)
        hr, low, high = hr_fields(row, nhanes=nhanes)
        records.append(
            {
                "cohort": cohort,
                "assumption": assumption,
                "analysis": analysis,
                "n": int(row["n"]),
                "deaths": int(row["events"]),
                "result": fmt_hr(hr, low, high),
                "hr": hr,
                "ci_low": low,
                "ci_high": high,
                "analysis_id": aid,
                "status": "ESTIMATED",
            }
        )
    records.extend(
        [
            {
                "cohort": "CHARLS",
                "assumption": "Anthropometric auxiliaries",
                "analysis": "Later-wave auxiliary variant",
                "n": np.nan,
                "deaths": np.nan,
                "result": "Not stably estimable; no mortality model fitted",
                "hr": np.nan,
                "ci_low": np.nan,
                "ci_high": np.nan,
                "analysis_id": "who_block_later_aux",
                "status": "NOT_STABLY_ESTIMABLE",
            },
            {
                "cohort": "NHANES",
                "assumption": "Older-adult routing boundary",
                "analysis": "Age 60–79, A1",
                "n": 3346,
                "deaths": 701,
                "result": "Not estimable; no coefficient reported",
                "hr": np.nan,
                "ci_low": np.nan,
                "ci_high": np.nan,
                "analysis_id": "v43_s3_nhanes_pn_e0_age60_79_a1",
                "status": "NOT_ESTIMABLE",
            },
            {
                "cohort": "NHANES",
                "assumption": "Older-adult routing boundary",
                "analysis": "Age 60–79, A2",
                "n": 3346,
                "deaths": 701,
                "result": "Not attempted after the paired A1 model was not estimable",
                "hr": np.nan,
                "ci_low": np.nan,
                "ci_high": np.nan,
                "analysis_id": "v43_s3_nhanes_pn_e0_age60_79_a2",
                "status": "NOT_ATTEMPTED_PANEL_DEPENDENCY",
            },
        ]
    )
    return pd.DataFrame(records)


def write_table3(table3: pd.DataFrame) -> Path:
    lines = [
        r"\begin{longtable}{llp{0.32\textwidth}rrp{0.25\textwidth}}",
        r"\caption{Sensitivity analyses of measurement and analytic choices}\label{tab:sensitivity}\\",
        r"\toprule",
        r"Cohort & Assumption & Analysis & n & Deaths & Result \\",
        r"\midrule",
        r"\endfirsthead",
        r"\toprule",
        r"Cohort & Assumption & Analysis & n & Deaths & Result \\",
        r"\midrule",
        r"\endhead",
        r"\bottomrule",
        r"\endfoot",
    ]
    display_rows = table3.loc[table3["status"] == "ESTIMATED"]
    for row in display_rows.itertuples(index=False):
        lines.append(
            f"{row.cohort} & {tex_escape(row.assumption)} & "
            f"{tex_escape(row.analysis)} & {fmt_int(row.n)} & "
            f"{fmt_int(row.deaths)} & {tex_escape(row.result)} \\\\"
        )
    lines.extend(
        [
            r"\addlinespace",
            r"\multicolumn{6}{p{0.96\textwidth}}{\footnotesize Results are A1 hazard ratios (95\% confidence intervals) per sex-specific within-cohort weighted SD lower raw PEF. Complete-case analyses use the common A2-complete sample. The NHANES A/B/C definition is the broadest spirometry-quality sensitivity. WHO and PCORnet-style ranges are operational anthropometric data ranges.} \\",
            r"\end{longtable}",
        ]
    )
    path = TABLE_DIR / "table3_sensitivity.tex"
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    table3.to_csv(SOURCE_DIR / "table3_source_data.csv", index=False)
    return path


def write_manifest(outputs: list[Path]) -> Path:
    rows = []
    for name, (path, expected) in SOURCES.items():
        rows.append(
            {
                "asset_type": "input",
                "asset_id": name,
                "relative_path": str(path.relative_to(PROJECT)),
                "sha256": sha256(path),
                "expected_sha256": expected,
                "status": "PASS",
            }
        )
    for path in sorted(set(outputs)):
        rows.append(
            {
                "asset_type": "output",
                "asset_id": path.stem,
                "relative_path": str(path.relative_to(PROJECT)),
                "sha256": sha256(path),
                "expected_sha256": "",
                "status": "GENERATED_FROM_VERIFIED_AGGREGATES",
            }
        )
    manifest = pd.DataFrame(rows)
    path = OUT / "display_build_manifest.csv"
    manifest.to_csv(path, index=False)
    return path


def main() -> None:
    verify_sources()
    prepare_dirs()
    set_plot_style()
    charls_registry = pd.read_csv(SOURCES["charls_registry"][0])
    charls_anthro = pd.read_csv(SOURCES["charls_anthro"][0])
    nhanes_registry = pd.read_csv(SOURCES["nhanes_registry"][0])
    charls_flow = pd.read_csv(SOURCES["charls_flow"][0])
    nhanes_flow = pd.read_csv(SOURCES["nhanes_flow"][0])
    table1 = pd.read_csv(SOURCES["table1"][0], dtype=str, keep_default_na=False)

    outputs: list[Path] = []
    _, fig1 = build_figure1(
        charls_flow, nhanes_flow, charls_registry, nhanes_registry
    )
    outputs.extend(fig1)
    _, fig2 = build_figure2(charls_registry, nhanes_registry)
    outputs.extend(fig2)
    _, fig3 = build_figure3(charls_registry, charls_anthro, nhanes_registry)
    outputs.extend(fig3)
    outputs.append(write_table1(table1))
    table2 = table2_records(charls_registry, nhanes_registry)
    outputs.append(write_table2(table2))
    table3 = build_table3_records(charls_registry, charls_anthro, nhanes_registry)
    outputs.append(write_table3(table3))
    outputs.extend(sorted(SOURCE_DIR.glob("*.csv")))
    manifest = write_manifest(outputs)
    print(
        "DISPLAY BUILD PASS; verified 6 frozen aggregate inputs; "
        f"generated {len(outputs)} display/source artifacts and {manifest.relative_to(PROJECT)}"
    )


if __name__ == "__main__":
    main()
