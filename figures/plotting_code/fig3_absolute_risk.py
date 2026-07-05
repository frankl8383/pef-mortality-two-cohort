"""
Publication figure generator — extracted from analysis session lineage.
Figures were rendered with matplotlib in a Python 3.11 environment.
Inputs are the aggregate summary CSVs in ../summary_data/ (no row-level data).
Local absolute paths in the original session have been left as-is inside string
literals ONLY where they reference artifact inputs; replace with your local copy
of the corresponding summary CSV. No controlled/row-level data is required.
"""

import matplotlib.pyplot as plt
import numpy as np

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


apply_figure_style()

# verified from nhanes_v38_1_absolute_mortality_risk.csv
q=["Q1","Q2","Q3","Q4"]
y5=[1.8,2.5,3.0,7.3]; y5lo=[1.1,1.6,2.0,6.0]; y5hi=[2.5,3.4,3.9,8.6]
y10=[12.2,12.3,14.2,27.2]; y10lo=[8.9,9.6,11.1,23.6]; y10hi=[15.6,14.9,17.2,30.7]
rd5=[0.0,0.7,1.2,5.5]; rd10=[0.0,0.1,1.9,14.9]

fig,(ax1,ax2)=plt.subplots(1,2,figsize=(9.6,4.4))
x=np.arange(4); col="#1f4e79"; colhi="#9c5a2b"
bars_col=[col,col,col,colhi]

for ax,(yy,lo,hi,rd,lab,ymax) in zip((ax1,ax2),
        [(y5,y5lo,y5hi,rd5,"5-year",10.0),(y10,y10lo,y10hi,rd10,"10-year",34.0)]):
    err=[np.array(yy)-np.array(lo),np.array(hi)-np.array(yy)]
    ax.bar(x,yy,yerr=err,capsize=4,color=bars_col,edgecolor="black",linewidth=0.6,width=0.62,
           error_kw=dict(lw=1.0,ecolor="#444"))
    for xi,yi,hii,rdi in zip(x,yy,hi,rd):
        ax.text(xi,hii+ymax*0.02,f"{yi:.1f}%",ha="center",va="bottom",fontsize=8.2,fontweight="bold")
        if rdi>0:
            ax.text(xi,yi/2,f"+{rdi:.1f} pp",ha="center",va="center",fontsize=7.0,color="white",fontweight="bold")
    ax.set_xticks(x); ax.set_xticklabels(q,fontsize=9)
    ax.set_ylim(0,ymax); ax.set_title(f"{lab} all-cause mortality",fontsize=9.5,fontweight="bold")
    ax.set_xlabel("Lower-than-expected PEF marker quartile\n(Q1 = highest expected reserve \u2192 Q4 = lowest)",fontsize=8.2)
    ax.spines[["top","right"]].set_visible(False)
    ax.grid(axis="y",alpha=0.25,lw=0.6)
ax1.set_ylabel("Weighted mortality risk (%)",fontsize=9)
fig.suptitle("Figure 3. Absolute mortality risk across lower-than-expected PEF marker quartiles (NHANES 2007\u20132012)",
             x=0.02,ha="left",fontsize=10.2,fontweight="bold")
fig.text(0.02,0.005,"Survey-weighted fixed-time risks with 95% CI; risk difference vs Q1 annotated. N=7,035 (5-year and 10-year fixed-time samples per Supplementary Table S1).",
         fontsize=6.6,color="#555")
fig.tight_layout(rect=[0,0.02,1,0.94])
fig.savefig("figure3_absolute_risk_v39.png",dpi=600,bbox_inches="tight")
fig.savefig("figure3_absolute_risk_v39.pdf",bbox_inches="tight")
print("saved figure3_absolute_risk_v39")