"""
Publication figure generator — extracted from analysis session lineage.
Figures were rendered with matplotlib in a Python 3.11 environment.
Inputs are the aggregate summary CSVs in ../summary_data/ (no row-level data).
Local absolute paths in the original session have been left as-is inside string
literals ONLY where they reference artifact inputs; replace with your local copy
of the corresponding summary CSV. No controlled/row-level data is required.
"""

# skill:figure-style kernel.py (auto-injected on skill load)
META_GREY = "#888888"


def apply_figure_style(*, frame="open", font=None, sizes=(8, 7, 6), grid=False):
    import matplotlib as mpl
    if frame not in ("open", "boxed", "none"):
        raise ValueError(f"frame must be 'open'|'boxed'|'none', got {frame!r}")
    try:
        import os, sys, glob, matplotlib.font_manager as fm
        fdir = os.path.join(os.environ.get("CONDA_PREFIX") or sys.prefix, "fonts")
        if os.path.isdir(fdir):
            known = {f.fname for f in fm.fontManager.ttflist}
            for f in glob.glob(os.path.join(fdir, "*.ttf")):
                if f not in known:
                    fm.fontManager.addfont(f)
    except Exception:
        pass
    base, secondary, tick = sizes
    boxed = (frame == "boxed")
    rc = {
        "font.family": "sans-serif",
        "font.size": base,
        "axes.labelsize": base,
        "axes.titlesize": base,
        "legend.fontsize": secondary,
        "xtick.labelsize": tick,
        "ytick.labelsize": tick,
        "axes.linewidth": 0.6,
        "xtick.direction": "out", "ytick.direction": "out",
        "xtick.major.size": 3, "ytick.major.size": 3,
        "xtick.major.width": 0.6, "ytick.major.width": 0.6,
        "axes.spines.top": boxed, "axes.spines.right": boxed,
        "axes.spines.left": frame != "none", "axes.spines.bottom": frame != "none",
        "axes.grid": bool(grid),
        "legend.frameon": False,
        "figure.dpi": 200,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
        "axes.titleweight": "normal",
        "axes.titlelocation": "left",
        "axes.labelweight": "normal",
        "lines.linewidth": 1.2,
        "patch.linewidth": 0.6,
        "pdf.fonttype": 42, "ps.fonttype": 42,
    }
    if font:
        rc["font.sans-serif"] = [font, "DejaVu Sans"]
    mpl.rcParams.update(rc)


import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from survey import svydesign
import os

set.seed = None  # not needed in Python path

# Load the data
import pandas as pd

d0_path = "${HOME}/.claude-science/orgs/856a475a-d0be-450c-b905-47bbd7e2d546/artifacts/proj_05e754e6a4d7/602fecae-8210-4210-baea-dbdb0ab98075/v711cf3e6_incremental_value_M2base_models_v1_0.csv"

# The figure is produced by the R code block that uses ggplot2
# But since we need Python, we reproduce from the CSV data

# The R code reads the RDS and computes models; the figure uses models_tbl
# We need to load the models CSV and reproduce the ggplot in matplotlib

models_df = pd.read_csv(d0_path)

# The R code builds specs with these ids and labels
specs_labels = {
    "M0": "Base covariates",
    "M1": "Base + raw PEF (per SD)",
    "M2": "Base + GLI FEV1 z",
    "M3": "Base + residualized PEF marker",
    "M4": "Base + GLI FEV1 z + PEF marker",
    "M5": "Base + full GLI panel",
    "M6": "Base + full GLI panel + PEF marker",
}

# The models CSV should have id, C, C_lo, C_hi columns
# merge label if not present
if 'label' not in models_df.columns:
    models_df['label'] = models_df['id'].map(specs_labels)

apply_figure_style()

plot_df = models_df.copy()
plot_df['id'] = pd.Categorical(plot_df['id'], categories=list(reversed(list(specs_labels.keys()))), ordered=True)
plot_df = plot_df.sort_values('id')

fig, ax = plt.subplots(figsize=(7.6, 4.2))

y_positions = {mid: i for i, mid in enumerate(reversed(list(specs_labels.keys())))}

for _, row in plot_df.iterrows():
    yi = y_positions[row['id']]
    ax.errorbar(row['C'], yi,
                xerr=[[row['C'] - row['C_lo']], [row['C_hi'] - row['C']]],
                fmt='o', color='#1f4e79', ms=4, elinewidth=0.8, capsize=0)
    ax.text(row['C'], yi + 0.30, row['label'],
            ha='left', va='bottom', fontsize=2.9 * 3, color='grey')

order = list(reversed(list(specs_labels.keys())))
ax.set_yticks(list(range(len(order))))
ax.set_yticklabels(order)
ax.set_xlabel("Weighted Harrell's C (95% CI)")

n_used = int(models_df['n'].iloc[0]) if 'n' in models_df.columns else 7035
n_events = int(models_df['events'].iloc[0]) if 'events' in models_df.columns else 974

ax.set_title("Incremental discrimination for all-cause mortality")
ax.set_title(f"NHANES 2007-2012 linked mortality; common sample N={n_used}, deaths={n_events}",
             loc='left', fontsize=8, color='grey', pad=2)

ax.set_ylim(-0.7, len(order) - 0.15)
ax.set_xlim(ax.get_xlim())

fig.tight_layout()
fig.savefig("incremental_value_cindex_v1_0.png", dpi=600, bbox_inches='tight')
print("saved incremental_value_cindex_v1_0.png")