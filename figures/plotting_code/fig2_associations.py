"""
Publication figure generator — extracted from analysis session lineage.
Figures were rendered with matplotlib in a Python 3.11 environment.
Inputs are the aggregate summary CSVs in ../summary_data/ (no row-level data).
Local absolute paths in the original session have been left as-is inside string
literals ONLY where they reference artifact inputs; replace with your local copy
of the corresponding summary CSV. No controlled/row-level data is required.
"""

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import numpy as np

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


apply_figure_style(sizes=(8, 7, 6))

rows = [
    ("CHARLS", "incident chronic lung disease", "HR", 1.20, 1.07, 1.35, "P=0.002", "A"),
    ("CHARLS", "incident CLD, frailty-adjusted", "HR", 1.16, 1.04, 1.31, "P=0.010", "A"),
    ("CHARLS", "all-cause mortality", "HR", 1.33, 1.24, 1.43, "P<0.001", "A"),
    ("CHARLS", "all-cause mortality, frailty-adjusted", "HR", 1.24, 1.16, 1.33, "P<0.001", "A"),
    ("NHANES", "linked all-cause mortality", "HR", 1.38, 1.27, 1.51, "P<0.001", "A"),
    ("NHANES", "mortality, beyond full GLI panel", "HR", 1.33, 1.15, 1.54, "P<0.001", "A"),
    ("NHANES", "GLI airflow obstruction", "OR", 4.57, 3.95, 5.29, "P<0.001", "B"),
    ("NHANES", "PRISm-z", "OR", 3.05, 2.59, 3.60, "P<0.001", "B"),
]
A = [r for r in rows if r[7] == "A"]
B = [r for r in rows if r[7] == "B"]
ordered = A + B
n = len(ordered)
ys = []
y = 0.0
for i in range(n):
    if i == len(A):
        y -= 0.8
    ys.append(y)
    y -= 1.0
ys = np.array(ys)

fig = plt.figure(figsize=(9.6, 4.9))
gs = GridSpec(1, 3, width_ratios=[3.4, 4.0, 2.5], wspace=0.0)
axL = fig.add_subplot(gs[0])
axM = fig.add_subplot(gs[1])
axR = fig.add_subplot(gs[2])
for ax in (axL, axR, axM):
    ax.set_ylim(ys.min() - 1.1, ys.max() + 1.5)
axL.axis("off")
axR.axis("off")
axL.set_xlim(0, 1)
axR.set_xlim(0, 1)

colA = "#1f4e79"
colB = "#9c5a2b"
yb_top = ys[len(A):].max() + 0.5
yb_bot = ys[len(A):].min() - 0.5
axM.axhspan(yb_bot, yb_top, color="#f3e9e0", alpha=0.8, zorder=0)
axM.axvline(1.0, color="grey", lw=1.0, ls="--", zorder=1)

for (coh, desc, meas, est, lo, hi, ptxt, grp), yy in zip(ordered, ys):
    c = colA if grp == "A" else colB
    mk = "o" if meas == "HR" else "s"
    axM.plot([lo, hi], [yy, yy], color=c, lw=1.7, zorder=2)
    axM.plot([est], [yy], mk, color=c, ms=7.5, zorder=3)
    axL.text(0.02, yy, coh, ha="left", va="center", fontsize=8.4, fontweight="bold", color=c)
    axL.text(0.26, yy, desc, ha="left", va="center", fontsize=8.1, color="#222")
    axR.text(0.98, yy, f"{meas} {est:.2f} ({lo:.2f}\u2013{hi:.2f}); {ptxt}", ha="right", va="center", fontsize=7.9, color="#222")

axM.set_xscale("log")
axM.set_xlim(0.9, 6.6)
axM.set_xticks([1, 1.5, 2, 3, 4, 5, 6])
axM.set_xticklabels(["1.0", "1.5", "2.0", "3.0", "4.0", "5.0", "6.0"], fontsize=8)
axM.set_yticks([])
for s in ["left", "right", "top"]:
    axM.spines[s].set_visible(False)
axM.set_xlabel("Effect per 1-SD higher marker (log scale)\nHR = circles, OR = squares", fontsize=8.4)

axL.text(0.02, ys[:len(A)].max() + 0.95, "Independent evidence", ha="left", va="bottom", fontsize=8.9, fontweight="bold", color=colA)
axL.text(0.02, ys[len(A):].max() + 0.60, "Construct convergence", ha="left", va="bottom", fontsize=8.6, fontweight="bold", color=colB)

fig.suptitle("Figure 2. Cross-dataset associations for lower-than-expected peak expiratory flow",
             x=0.02, ha="left", fontsize=10.5, fontweight="bold")
fig.subplots_adjust(left=0.02, right=0.99, top=0.90, bottom=0.14)
fig.savefig("figure2_associations_v39.png", dpi=600, bbox_inches="tight")
fig.savefig("figure2_associations_v39.pdf", bbox_inches="tight")
print("saved clean figure2")