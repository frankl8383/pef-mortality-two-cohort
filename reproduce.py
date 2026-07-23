#!/usr/bin/env python3
"""Reproduce manuscript Tables 1-3 and Figures 1-3 from aggregate data.

The six inputs are disclosure-safe, manuscript-facing source-data files. This
script does not open participant-level data, run multiple imputation, or refit
the cohort models. It verifies every input before producing the display files
and an auditable SHA-256 manifest.
"""

from __future__ import annotations

import hashlib
import math
import textwrap
from pathlib import Path

import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from matplotlib.ticker import FuncFormatter, LogLocator, NullFormatter
from PIL import Image


ROOT = Path(__file__).resolve().parent
DATA_DIR = ROOT / "data"
OUT = ROOT / "output"
FIG_DIR = OUT / "figures"
TABLE_DIR = OUT / "tables"

INPUTS = {
    "figure1": (
        DATA_DIR / "figure1_source_data.csv",
        "ce799bc0351439d300647aa135fd662a412baa2b839f2c617144afa9fc5ca64e",
    ),
    "figure2": (
        DATA_DIR / "figure2_source_data.csv",
        "aa9f9f4c30b40f62be682f21caf0098b87b0d788f6b0a2af5bf6b48f1a9a2671",
    ),
    "figure3": (
        DATA_DIR / "figure3_source_data.csv",
        "370bc7275f640fbd3fc434a5248b9334d4970ae409b7ae3af6f89d79406d7b64",
    ),
    "table1": (
        DATA_DIR / "table1_source_data.csv",
        "2172d7cdaf5dfc57191a89a250b4892959ec81f9d81af860640210c157b3f4c4",
    ),
    "table2": (
        DATA_DIR / "table2_source_data.csv",
        "710ca9ad7c986b6f1fded4dcc34d5e042a301d5b9b7b8cc068ea5d25f7dfd0b4",
    ),
    "table3": (
        DATA_DIR / "table3_source_data.csv",
        "6d4b94b67b3b0202b1dcbd57f21cee6c8a2dc721ffe1b22bbeb4af98905f9d92",
    ),
}

EXPECTED_TABLES = {
    "table1_baseline.tex":
        "dee6968c47224f6ff07b4f7c600959cadafa70932771a6ecdcc0c7680eea72a5",
    "table2_primary_associations.tex":
        "f62925ad434ccaff1e8406b0f09c7bf1dc69f9776ed1e897998cf523ccdbbf8a",
    "table3_sensitivity.tex":
        "c3fcf286c59be719a2a993c940935cacbb9799d187bd1c86bd08c54be188b82f",
}

EXPECTED_TIFFS = {
    "figure1_cohort_flow.tiff": (4014, 3810),
    "figure2_adjustment_ladder.tiff": (4014, 2670),
    "figure3_measurement_robustness.tiff": (4014, 3720),
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


def verify_inputs() -> None:
    for name, (path, expected) in INPUTS.items():
        if not path.exists():
            raise FileNotFoundError(f"Missing aggregate input: {path.name}")
        observed = sha256(path)
        if observed != expected:
            raise RuntimeError(
                f"Aggregate input hash mismatch for {name}: "
                f"expected {expected}, observed {observed}"
            )


def prepare_dirs() -> None:
    for path in (FIG_DIR, TABLE_DIR):
        path.mkdir(parents=True, exist_ok=True)


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


def build_figure1(source: pd.DataFrame) -> list[Path]:
    sequential = source.loc[source["record_type"] == "sequential_flow"].copy()
    c = sequential.loc[sequential["cohort"] == "CHARLS"].sort_values("step_order")
    n_selected = sequential.loc[
        sequential["cohort"] == "NHANES"
    ].sort_values("step_order")
    branch_rows = source.loc[source["record_type"] == "analysis_branch"].copy()

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
    return save_figure(fig, "figure1_cohort_flow")


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


def build_figure2(source: pd.DataFrame) -> list[Path]:
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
    return save_figure(fig, "figure2_adjustment_ladder")


def build_figure3(source: pd.DataFrame) -> list[Path]:
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
    return save_figure(fig, "figure3_measurement_robustness")


def table1_characteristic_rows(table1: pd.DataFrame) -> list[str]:
    lines: list[str] = []
    for cohort, panel in (("CHARLS", "A"), ("NHANES", "B")):
        part = table1.loc[table1["cohort"] == cohort].copy()
        cohort_n = pd.to_numeric(part["group_n_unweighted"], errors="raise").unique()
        cohort_deaths = pd.to_numeric(
            part["group_events_unweighted"], errors="raise"
        ).unique()
        if len(cohort_n) != 1 or len(cohort_deaths) != 1:
            raise RuntimeError(f"Table 1 {cohort} panel totals are inconsistent")
        n, deaths = int(cohort_n[0]), int(cohort_deaths[0])
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
    return path


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
    return path


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
    return path


def validate_frames(frames: dict[str, pd.DataFrame]) -> None:
    expected_rows = {
        "figure1": 23,
        "figure2": 6,
        "figure3": 10,
        "table1": 43,
        "table2": 12,
        "table3": 20,
    }
    for name, expected in expected_rows.items():
        observed = len(frames[name])
        if observed != expected:
            raise RuntimeError(
                f"{name} must contain {expected} aggregate rows; found {observed}"
            )

    for name in ("figure2", "figure3"):
        primary = (
            frames[name]["primary"]
            .astype(str)
            .str.strip()
            .str.lower()
            .map({"true": True, "false": False})
        )
        if primary.isna().any():
            raise RuntimeError(f"{name} contains an invalid primary indicator")
        frames[name]["primary"] = primary.astype(bool)

    record_types = set(frames["figure1"]["record_type"].dropna())
    if record_types != {"sequential_flow", "analysis_branch"}:
        raise RuntimeError("Figure 1 record types are incomplete")

    expected_tiers = {"A0", "A1", "A2"}
    for cohort in ("CHARLS", "NHANES"):
        observed_tiers = set(
            frames["figure2"].loc[
                frames["figure2"]["cohort"] == cohort, "model_tier"
            ]
        )
        if observed_tiers != expected_tiers:
            raise RuntimeError(f"Figure 2 tiers are incomplete for {cohort}")

    primary_counts = (
        frames["figure3"].groupby("cohort")["primary"].sum().astype(int).to_dict()
    )
    if primary_counts != {"CHARLS": 1, "NHANES": 1}:
        raise RuntimeError("Figure 3 must contain one primary row per cohort")

    table1 = frames["table1"]
    empty_education = (
        (table1["cohort"] == "CHARLS")
        & (table1["variable_id"] == "education")
        & (table1["level_id"] == "none")
    )
    if empty_education.any():
        raise RuntimeError("Table 1 contains an unobserved education category")

    table2 = frames["table2"]
    combinations = table2[
        ["cohort", "scale", "model_tier"]
    ].drop_duplicates()
    if len(combinations) != 12:
        raise RuntimeError("Table 2 cohort-scale-tier combinations are not unique")

    table3 = frames["table3"]
    if int((table3["status"] == "ESTIMATED").sum()) != 17:
        raise RuntimeError("Table 3 must contain 17 estimated rows")
    if int((table3["status"] != "ESTIMATED").sum()) != 3:
        raise RuntimeError("Table 3 must contain three non-estimated rows")

    prohibited = {
        "participant_id",
        "person_id",
        "respondent_id",
        "household_id",
        "medical_record_number",
        "seqn",
    }
    for name, frame in frames.items():
        hits = prohibited.intersection(
            {str(column).lower() for column in frame.columns}
        )
        if hits:
            raise RuntimeError(
                f"{name} contains participant identifiers: {sorted(hits)}"
            )


def validate_outputs(outputs: list[Path]) -> None:
    missing = [str(path) for path in outputs if not path.is_file()]
    if missing:
        raise RuntimeError(f"Expected outputs were not created: {missing}")

    for filename, expected in EXPECTED_TABLES.items():
        path = TABLE_DIR / filename
        observed = sha256(path)
        if observed != expected:
            raise RuntimeError(
                f"Table hash mismatch for {filename}: "
                f"expected {expected}, observed {observed}"
            )

    for filename, expected_size in EXPECTED_TIFFS.items():
        path = FIG_DIR / filename
        with Image.open(path) as image:
            if image.mode != "RGB":
                raise RuntimeError(f"{filename} must be RGB; found {image.mode}")
            if image.size != expected_size:
                raise RuntimeError(
                    f"{filename} must be {expected_size}; found {image.size}"
                )
            compression = str(image.info.get("compression", "")).lower()
            if compression not in {"tiff_lzw", "lzw"}:
                raise RuntimeError(
                    f"{filename} must use LZW compression; "
                    f"found {compression or 'unknown'}"
                )
            dpi = image.info.get("dpi")
            if (
                not dpi
                or len(dpi) != 2
                or any(abs(float(value) - 600.0) > 1.0 for value in dpi)
            ):
                raise RuntimeError(
                    f"{filename} must have 600-dpi metadata; found {dpi}"
                )


def write_manifest(outputs: list[Path]) -> Path:
    rows = []
    for name, (path, expected) in INPUTS.items():
        rows.append(
            {
                "asset_type": "input",
                "asset_id": name,
                "relative_path": str(path.relative_to(ROOT)),
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
                "relative_path": str(path.relative_to(ROOT)),
                "sha256": sha256(path),
                "expected_sha256": "",
                "status": "GENERATED",
            }
        )
    manifest = pd.DataFrame(rows)
    path = OUT / "manifest.csv"
    manifest.to_csv(path, index=False)
    return path


def main() -> None:
    verify_inputs()
    prepare_dirs()
    set_plot_style()
    frames = {
        "figure1": pd.read_csv(INPUTS["figure1"][0]),
        "figure2": pd.read_csv(INPUTS["figure2"][0]),
        "figure3": pd.read_csv(INPUTS["figure3"][0]),
        "table1": pd.read_csv(
            INPUTS["table1"][0], dtype=str, keep_default_na=False
        ),
        "table2": pd.read_csv(INPUTS["table2"][0]),
        "table3": pd.read_csv(INPUTS["table3"][0]),
    }
    validate_frames(frames)

    outputs: list[Path] = []
    outputs.extend(build_figure1(frames["figure1"]))
    outputs.extend(build_figure2(frames["figure2"]))
    outputs.extend(build_figure3(frames["figure3"]))
    outputs.append(write_table1(frames["table1"]))
    outputs.append(write_table2(frames["table2"]))
    outputs.append(write_table3(frames["table3"]))
    validate_outputs(outputs)
    manifest = write_manifest(outputs)
    print(
        "REPRODUCTION PASS; verified 6 aggregate inputs and generated "
        f"{len(outputs)} table/figure files plus {manifest.relative_to(ROOT)}"
    )


if __name__ == "__main__":
    main()
