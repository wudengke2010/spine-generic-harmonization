#!/usr/bin/env python3
"""Regenerate HBM main-text figures Fig2/Fig3 from complete-case results (v1.3.0).

Replaces ad-hoc inline code used in an earlier session so the exact
submission figures are reproducible from versioned scripts + CSVs.

Fig2 (submission) — Harmonisation performance across endpoints:
  A. Univariate vendor eta2 (log10), grouped bars HC + ALL, baseline lines
  B. PERMANOVA R2 (log10), grouped bars HC + ALL, baseline lines
  C. Age semi-partial R2 (grouped bars, Original baseline shown)
  D. Sex semi-partial R2 (grouped bars, Original baseline shown)

Fig3 (submission) — Pareto frontier:
  A. HC cohort: -log10(PERMANOVA R2) vs r_pre_post
  B. ALL cohort: same

Data (complete-case unified sample, HC N=188 / ALL N=246):
  E:/boshi/qm_harmonization_paper/results/step4_eval/cc_summary.csv
  E:/boshi/qm_harmonization_paper/results/step4_eval/permanova_results.csv
"""
from __future__ import annotations
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 9, "axes.titlesize": 10, "axes.labelsize": 9,
    "xtick.labelsize": 8, "ytick.labelsize": 8, "legend.fontsize": 8,
    "axes.linewidth": 0.8, "pdf.fonttype": 42, "ps.fonttype": 42,
    "figure.dpi": 300, "savefig.dpi": 300,
    "axes.spines.top": False, "axes.spines.right": False,
})

EVAL_DIR = Path("E:/boshi/qm_harmonization_paper/results/step4_eval")
OUT_DIRS = [Path("C:/Users/admin/WorkBuddy/2026-07-06-12-34-04/submission/figures"),
            Path("C:/Users/admin/WorkBuddy/2026-07-06-12-34-04/repo_push_tmp/figures")]

# Wong 2011 color-blind-safe palette
WONG = {"blue": "#0072B2", "sky": "#56B4E9", "green": "#009E73",
        "orange": "#E69F00", "vermilion": "#D55E00", "pink": "#CC79A7",
        "grey": "#999999"}

METHODS = ["LME", "ComBat", "ComBat-joint", "RELIEF", "CovBat"]
M_LABEL = {"LME": "LME/FE", "ComBat": "ComBat", "ComBat-joint": "ComBat-joint",
           "RELIEF": "RELIEF", "CovBat": "CovBat"}
M_COLOR = {"LME": WONG["pink"], "ComBat": WONG["blue"],
           "ComBat-joint": WONG["sky"], "RELIEF": WONG["orange"],
           "CovBat": WONG["green"]}
COHORT_COLOR = {"HC": WONG["vermilion"], "ALL": WONG["grey"]}


def load():
    cc = pd.read_csv(EVAL_DIR / "cc_summary.csv")
    perm = pd.read_csv(EVAL_DIR / "permanova_results.csv")
    return cc, perm


def grouped_bar_panel(ax, cc, metric_col, baseline_source, title, zoomed):
    """Grouped bars: methods x (HC, ALL); baseline horizontal lines."""
    x = np.arange(len(METHODS))
    w = 0.38
    for i, cohort in enumerate(["HC", "ALL"]):
        vals = []
        for m in METHODS:
            row = cc[(cc["cohort"] == cohort) & (cc["method"] == m)]
            vals.append(row[metric_col].iloc[0])
        bars = ax.bar(x + (i - 0.5) * w, vals, width=w,
                      color=COHORT_COLOR[cohort], alpha=0.85,
                      label=f"{cohort} (N={188 if cohort=='HC' else 246})")
        for b, v in zip(bars, vals):
            ax.text(b.get_x() + b.get_width() / 2,
                    b.get_height() * 1.15 if zoomed else b.get_height(),
                    f"{v:.1e}".replace("e-0", "e-") if zoomed else f"{v:.4f}",
                    ha="center", va="bottom", fontsize=6, rotation=0)
    # baseline lines
    if isinstance(baseline_source, dict):
        for cohort in ["HC", "ALL"]:
            bval = baseline_source[cohort]
            ax.axhline(bval, color=COHORT_COLOR[cohort], ls="--", lw=0.8, alpha=0.5)
    ax.set_xticks(x)
    ax.set_xticklabels([M_LABEL[m] for m in METHODS], rotation=15)
    ax.set_title(title, fontsize=9)
    return ax


def gen_fig2(cc, perm):
    fig, axes = plt.subplots(2, 2, figsize=(12, 8.5))

    # Panel A: univariate eta2 (log scale, zoomed)
    ax = axes[0, 0]
    grouped_bar_panel(ax, cc, "eta2_mean",
                      {"HC": cc[(cc.cohort == "HC") & (cc.method == "Original")]["eta2_mean"].iloc[0],
                       "ALL": cc[(cc.cohort == "ALL") & (cc.method == "Original")]["eta2_mean"].iloc[0]},
                      "Univariate vendor effect (zoomed)", zoomed=True)
    ax.set_yscale("log")
    ax.set_ylim(5e-5, 3e-3)
    ax.set_ylabel(r"log$_{10}$($\eta^2$)")
    # baseline annotations
    ax.text(0.99, cc[(cc.cohort == "HC") & (cc.method == "Original")]["eta2_mean"].iloc[0],
            r"  Baseline HC $\eta^2$=0.315", transform=ax.get_yaxis_transform(),
            ha="right", va="bottom", fontsize=6.5, color=COHORT_COLOR["HC"])
    ax.text(0.99, cc[(cc.cohort == "ALL") & (cc.method == "Original")]["eta2_mean"].iloc[0],
            r"  Baseline ALL $\eta^2$=0.336", transform=ax.get_yaxis_transform(),
            ha="right", va="bottom", fontsize=6.5, color=COHORT_COLOR["ALL"])
    ax.text(0.02, 0.98, "A", transform=ax.transAxes, fontsize=12,
            fontweight="bold", va="top")

    # Panel B: PERMANOVA R2 (log scale, zoomed)
    ax = axes[0, 1]
    perm_p = perm.set_index(["cohort", "method"])["R2"]
    # For complete-case fits the PERMANOVA R2 equals mean univariate eta2
    # (z-score Euclidean identity); pull directly from permanova_results.csv.
    cc_perm = perm.copy()
    grouped_bar_panel(ax, cc_perm.assign(
        eta2_mean=[perm_p.get((r.cohort, r.method), np.nan) for r in cc_perm.itertuples()]
    ) if False else cc_perm, "R2", None, "Multivariate vendor effect (zoomed)", zoomed=True)
    ax.set_yscale("log")
    ax.set_ylim(5e-5, 3e-3)
    ax.set_ylabel(r"log$_{10}$(PERMANOVA $R^2$)")
    for cohort, lab in [("HC", "Baseline HC $R^2$=0.315"), ("ALL", "Baseline ALL $R^2$=0.336")]:
        bval = perm[(perm.cohort == cohort) & (perm.method == "Original")]["R2"].iloc[0]
        ax.axhline(bval, color=COHORT_COLOR[cohort], ls="--", lw=0.8, alpha=0.5)
        ax.text(0.99, bval, f"  {lab}", transform=ax.get_yaxis_transform(),
                ha="right", va="bottom", fontsize=6.5, color=COHORT_COLOR[cohort])
    ax.text(0.02, 0.98, "B", transform=ax.transAxes, fontsize=12,
            fontweight="bold", va="top")

    # Panel C: age R2
    ax = axes[1, 0]
    grouped_bar_panel(ax, cc, "age_R2_mean", None,
                      "Age semi-partial $R^2$", zoomed=False)
    ax.set_ylabel(r"Age $R^2$ (mean)")
    for cohort in ["HC", "ALL"]:
        bval = cc[(cc.cohort == cohort) & (cc.method == "Original")]["age_R2_mean"].iloc[0]
        ax.axhline(bval, color=COHORT_COLOR[cohort], ls="--", lw=0.8, alpha=0.5)
    ax.text(0.99, 0.97, "best", transform=ax.transAxes, ha="right",
            va="top", fontsize=7, style="italic", color="#555")
    ax.text(0.02, 0.98, "C", transform=ax.transAxes, fontsize=12,
            fontweight="bold", va="top")

    # Panel D: sex R2
    ax = axes[1, 1]
    grouped_bar_panel(ax, cc, "sex_R2_mean", None,
                      "Sex semi-partial $R^2$", zoomed=False)
    ax.set_ylabel(r"Sex $R^2$ (mean)")
    for cohort in ["HC", "ALL"]:
        bval = cc[(cc.cohort == cohort) & (cc.method == "Original")]["sex_R2_mean"].iloc[0]
        ax.axhline(bval, color=COHORT_COLOR[cohort], ls="--", lw=0.8, alpha=0.5)
    ax.text(0.99, 0.97, "best", transform=ax.transAxes, ha="right",
            va="top", fontsize=7, style="italic", color="#555")
    ax.text(0.02, 0.98, "D", transform=ax.transAxes, fontsize=12,
            fontweight="bold", va="top")

    axes[0, 0].legend(loc="upper left", bbox_to_anchor=(0.0, -0.12),
                      ncol=2, frameon=False)
    fig.suptitle("Harmonisation performance across evaluation endpoints",
                 fontsize=12, fontweight="bold", y=0.995)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    for d in OUT_DIRS:
        fig.savefig(d / "Fig2.pdf", dpi=300, bbox_inches="tight")
        fig.savefig(d / "Fig2.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("Fig2 (submission) saved.")


def gen_fig3(cc, perm):
    fig, axes = plt.subplots(1, 2, figsize=(12, 5.2))
    for ax, cohort, n in [(axes[0], "HC", 188), (axes[1], "ALL", 246)]:
        xs, ys, names = [], [], []
        for m in METHODS:
            r2 = perm[(perm.cohort == cohort) & (perm.method == m)]["R2"].iloc[0]
            r = cc[(cc.cohort == cohort) & (cc.method == m)]["r_pre_post_mean"].iloc[0]
            xs.append(-np.log10(r2))
            ys.append(r)
            names.append(M_LABEL[m])
        for x, y, m in zip(xs, ys, METHODS):
            ax.scatter(x, y, s=110, color=M_COLOR[m], zorder=3, edgecolor="white", linewidth=1)
        # Pareto frontier (upper-left envelope)
        pts = sorted(zip(xs, ys), key=lambda t: t[1], reverse=True)
        front = [pts[0]]
        for p in pts[1:]:
            if p[0] < front[-1][0]:
                front.append(p)
        ax.plot([p[0] for p in front], [p[1] for p in front],
                color="#555", ls="--", lw=0.9, zorder=2)
        for x, y, m in zip(xs, ys, METHODS):
            ax.annotate(M_LABEL[m], (x, y), textcoords="offset points",
                        xytext=(7, -3), fontsize=8)
        ax.set_xlabel(r"$-$log$_{10}$(PERMANOVA $R^2$)")
        ax.set_ylabel(r"$r_{\mathrm{pre,post}}$ (subject preservation)")
        ax.set_title(f"{'A' if cohort=='HC' else 'B'}  Pareto frontier — {cohort} cohort (N={n})",
                     fontsize=10, loc="left")
    fig.suptitle("Central trade-off: multivariate decontamination vs. subject-level preservation",
                 fontsize=12, fontweight="bold", y=1.0)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    for d in OUT_DIRS:
        fig.savefig(d / "Fig3.pdf", dpi=300, bbox_inches="tight")
        fig.savefig(d / "Fig3.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("Fig3 (submission) saved.")


if __name__ == "__main__":
    cc, perm = load()
    gen_fig2(cc, perm)
    gen_fig3(cc, perm)
    print("DONE.")
