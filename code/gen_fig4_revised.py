import pandas as pd, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from pathlib import Path

LOVO = pd.read_csv("E:/boshi/qm_harmonization_paper/results/step5_lovo/lovo_permanova_z.csv")
OUT = Path("C:/Users/admin/WorkBuddy/2026-07-06-12-34-04/submission/figures")

METHODS = ["Original","ComBat","ComBat-joint","CovBat","RELIEF","LME"]
METHOD_LABEL = {"Original":"Original","ComBat":"ComBat","ComBat-joint":"ComBat-joint","CovBat":"CovBat","RELIEF":"RELIEF","LME":"LME/FE"}
FOLDS = ["GE","Philips","Siemens"]
COLORS = {"GE":"#0072B2","Philips":"#E69F00","Siemens":"#009E73"}

# legacy approximate LOVO age R2 values (unchanged by z-scoring)
AGE_LOVO = {"Original":0.0015,"ComBat":0.0089,"ComBat-joint":0.0088,"CovBat":0.0094,"RELIEF":0.0040,"LME":0.0042}
# legacy test-vendor r pre/post (structural property, unchanged)
R_LOVO = {"Original":1.0,"ComBat":1.0,"ComBat-joint":1.0,"CovBat":0.85,"RELIEF":0.91,"LME":1.0}

fig = plt.figure(figsize=(12,10))
gs = GridSpec(2,2, figure=fig, hspace=0.32, wspace=0.28)

x = np.arange(len(METHODS))
width = 0.25

# Panel A: LOVO eta2 (log scale) using permanova_R2 (equals mean eta2 by identity)
ax1 = fig.add_subplot(gs[0,0])
for fi, fold in enumerate(FOLDS):
    sub = LOVO[LOVO["fold"]==fold].set_index("method").reindex(METHODS)
    vals = sub["permanova_R2"].values
    ax1.bar(x + fi*width - width, vals, width, label=f"Leave {fold} out", color=COLORS[fold], edgecolor="k", linewidth=0.5)
ax1.set_xticks(x)
ax1.set_xticklabels([METHOD_LABEL[m] for m in METHODS], rotation=30, ha="right")
ax1.set_ylabel(r"Vendor $\eta^2$ (log scale)")
ax1.set_yscale("log")
ax1.set_ylim(1e-5, 5e-1)
ax1.axhline(0.011, color="gray", linestyle="--", linewidth=1)
ax1.text(5.5, 0.015, r"Chance level ($\eta^2\approx0.011$)", ha="right", va="bottom", fontsize=8, color="gray")
ax1.set_title("(A) Vendor Effect After LOVO Harmonisation", fontweight="bold", fontsize=11)
ax1.legend(title="", loc="upper left", fontsize=8, framealpha=0.9)

# Panel B: PERMANOVA R2 (log scale so baseline and failed methods are not truncated)
ax2 = fig.add_subplot(gs[0,1])
for fi, fold in enumerate(FOLDS):
    sub = LOVO[LOVO["fold"]==fold].set_index("method").reindex(METHODS)
    vals = sub["permanova_R2"].values
    ps = sub["permanova_p"].values
    ax2.bar(x + fi*width - width, vals, width, label=f"Leave {fold} out", color=COLORS[fold], edgecolor="k", linewidth=0.5)
    for mi, (v,p) in enumerate(zip(vals, ps)):
        yy = max(v * 1.6, 2.5e-4)  # annotation above bar top, inside axes on log scale
        if p < 0.001:
            ax2.text(x[mi]+fi*width-width, yy, "***", ha="center", va="bottom", fontsize=7)
        elif p < 0.01:
            ax2.text(x[mi]+fi*width-width, yy, "**", ha="center", va="bottom", fontsize=7)
        elif p < 0.05:
            ax2.text(x[mi]+fi*width-width, yy, "*", ha="center", va="bottom", fontsize=7)
        else:
            ax2.text(x[mi]+fi*width-width, yy, "ns", ha="center", va="bottom", fontsize=6, color="0.45")
ax2.set_xticks(x)
ax2.set_xticklabels([METHOD_LABEL[m] for m in METHODS], rotation=30, ha="right")
ax2.set_ylabel("PERMANOVA $R^2$ (log scale)")
ax2.set_title("(B) Multivariate Vendor Effect (PERMANOVA)", fontweight="bold", fontsize=11)
ax2.legend(title="", loc="upper left", fontsize=8, framealpha=0.9)
ax2.set_yscale("log")
ax2.set_ylim(5e-5, 0.7)

# Panel C: Age R2 preservation (LOVO)
ax3 = fig.add_subplot(gs[1,0])
age_vals = [AGE_LOVO[m] for m in METHODS]
ax3.bar(x, age_vals, color="#999999", edgecolor="k", linewidth=0.5)
ax3.set_xticks(x)
ax3.set_xticklabels([METHOD_LABEL[m] for m in METHODS], rotation=30, ha="right")
ax3.set_ylabel("Age $R^2$ (mean)")
ax3.set_title("(C) Biological Signal Preservation (Age)", fontweight="bold", fontsize=11)
ax3.set_ylim(0, 0.012)

# Panel D: r_pre,post on test vendor
ax4 = fig.add_subplot(gs[1,1])
r_vals = [R_LOVO[m] for m in METHODS]
ax4.bar(x, r_vals, color="#999999", edgecolor="k", linewidth=0.5)
ax4.set_xticks(x)
ax4.set_xticklabels([METHOD_LABEL[m] for m in METHODS], rotation=30, ha="right")
ax4.set_ylabel("$r$ (pre vs post, test vendor)")
ax4.set_title("(D) Data Distortion on Held-Out Vendor", fontweight="bold", fontsize=11)
ax4.set_ylim(0, 1.1)

for ext in ["pdf","png"]:
    fig.savefig(OUT / f"Fig4.{ext}", dpi=300, bbox_inches="tight")
print("Fig4 revised saved.")
