"""
Publication figure generator — extracted from analysis session lineage.
Figures were rendered with matplotlib in a Python 3.11 environment.
Inputs are the aggregate summary CSVs in ../summary_data/ (no row-level data).
Local absolute paths in the original session have been left as-is inside string
literals ONLY where they reference artifact inputs; replace with your local copy
of the corresponding summary CSV. No controlled/row-level data is required.
"""

import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

fig=plt.figure(figsize=(13.2,6.7))
gsL=fig.add_axes([0.0,0.0,0.42,0.93]); gsR=fig.add_axes([0.44,0.0,0.56,0.93])
for ax in (gsL,gsR): ax.set_xlim(0,10); ax.set_ylim(0,10); ax.axis("off")
def box(ax,x,y,w,h,lab,face="#eaf1f8",edge="#1f4e79",fs=7.2):
    ax.add_patch(FancyBboxPatch((x-w/2,y-h/2),w,h,boxstyle="round,pad=0.02,rounding_size=0.08",lw=1.1,edgecolor=edge,facecolor=face))
    ax.text(x,y,lab,ha="center",va="center",fontsize=fs,linespacing=1.3)
def arr(ax,x1,y1,x2,y2,c="#333",lw=1.0,ms=11):
    ax.add_patch(FancyArrowPatch((x1,y1),(x2,y2),arrowstyle="-|>",mutation_scale=ms,lw=lw,color=c))
gsL.text(0.2,9.85,"A  CHARLS (China)",fontweight="bold",fontsize=9.5,ha="left",va="top")
box(gsL,5,8.7,6.4,1.0,"Baseline respiratory-score eligible cohort\n(CHARLS Wave 1)\nn = 12,633")
arr(gsL,5,8.2,5,7.35)
gsL.text(6.35,7.78,"prospective follow-up\nto Wave 5 / 2020",ha="left",va="center",fontsize=6.2,style="italic",color="#666",linespacing=1.2)
box(gsL,2.7,6.0,4.5,1.5,"Incident chronic lung disease\nanalytic cohort\n\nn = 11,066\n1,647 incident events",fs=6.9)
box(gsL,7.3,6.0,4.5,1.5,"All-cause death\nanalytic cohort\n\nn = 12,569\n1,681 deaths",fs=6.9)
arr(gsL,5,7.3,2.9,6.8); arr(gsL,5,7.3,7.1,6.8)
box(gsL,2.7,3.7,4.7,1.35,"Excluded from incident-CLD analysis:\nprevalent baseline chronic lung\ndisease, or missing Wave 5 status\n/ model covariates",face="#f7f7f7",edge="#999",fs=6.3)
arr(gsL,2.7,5.25,2.7,4.4,c="#999",lw=0.85,ms=9)
box(gsL,5,1.45,8.4,1.15,"Both endpoints \u2014 survey-weighted Cox covariates:\nresidualized PEF marker + age (decade) + sex + ever-smoking + BMI\n(frailty-proxy count added in sensitivity model)",face="#fff7e6",edge="#d9a441",fs=6.7)
arr(gsL,2.7,3.02,3.5,2.05,c="#d9a441",lw=0.8,ms=8); arr(gsL,7.3,5.25,6.5,2.05,c="#d9a441",lw=0.8,ms=8)
gsR.text(0.2,9.85,"B  NHANES (United States)",fontweight="bold",fontsize=9.5,ha="left",va="top")
flow=[("Merged NHANES 2007\u20132012 person rows","30,442"),("Adults aged \u2265 45 years","10,270"),
      ("Adults 45+ with positive MEC weight","9,874"),("Adults 45+ with PEF available","7,108"),
      ("Primary PEF residual model\n(residualized on age, sex, height)","7,062"),
      ("All-cause mortality analytic sample\n(complete covariates; incremental analysis)","7,035")]
excl={0:["Excluded: age < 45 y (\u221220,172)"],1:["Excluded: zero/negative MEC","examination weight (\u2212396)"],
      2:["Excluded: PEF not available (\u22122,766)"],3:["Excluded: PEF quality / model","criteria not met (\u221246)"],
      4:["Excluded: missing model covariates (\u221227)"]}
n=len(flow); top=9.0; bot=1.95; gap=(top-bot)/(n-1); ys=[top-i*gap for i in range(n)]
xm=3.0; bw=5.0; bh=gap*0.56
for i,(lab,ns) in enumerate(flow): box(gsR,xm,ys[i],bw,bh,f"{lab}\nn = {ns}",fs=6.5)
for i in range(n-1): arr(gsR,xm,ys[i]-bh/2,xm,ys[i+1]+bh/2,ms=10,lw=0.95)
for idx,lines in excl.items():
    ymid=(ys[idx]+ys[idx+1])/2; xe=7.15; ew=4.7; eh=max(0.5,0.3*len(lines)+0.3)
    gsR.add_patch(FancyBboxPatch((xe-ew/2,ymid-eh/2),ew,eh,boxstyle="round,pad=0.02,rounding_size=0.05",lw=0.85,edgecolor="#999",facecolor="#f7f7f7"))
    gsR.text(xe,ymid,"\n".join(lines),ha="center",va="center",fontsize=6.1,linespacing=1.2)
    arr(gsR,xm+bw/2,ymid,xe-ew/2,ymid,c="#999",lw=0.75,ms=8)
gsR.add_patch(FancyBboxPatch((0.5,0.2),9.0,1.1,boxstyle="round,pad=0.02,rounding_size=0.05",lw=1.0,edgecolor="#d9a441",facecolor="#fff7e6"))
gsR.text(5,0.75,"Mortality models (survey-weighted Cox): M1 age + sex + race/ethnicity + BMI + smoking;\nM2 (primary) + education + income-poverty ratio.  Objective phenotypes (obstruction, PRISm-z): age + sex + race/ethnicity + BMI + smoking",ha="center",va="center",fontsize=6.0,linespacing=1.35)
fig.suptitle("Figure 1. Participant inclusion and exclusion (STROBE)",fontweight="bold",fontsize=11,x=0.01,ha="left",y=0.995)
fig.savefig("figure1_strobe_v39.png",dpi=600,bbox_inches="tight")
fig.savefig("figure1_strobe_v39.pdf",bbox_inches="tight")
print("saved fixed figure1")