#!/usr/bin/env python3
"""
Regenerate ALL 7 figures from PRE-COMPUTED analysis results.
Uses biomarkers_master.csv (267 subjects, 44 sites, 8 biomarkers)
and all pre-computed harmonized outputs + eval metrics.

Data sources:
  E:/boshi/qm_harmonization_paper/results/biomarkers_master.csv
  E:/boshi/qm_harmonization_paper/results/step4_eval/*.csv
  E:/qm_harmonization_paper/derivatives/all-files/results_design_C.csv
"""
from __future__ import annotations
import os, sys, warnings
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle, Patch
from matplotlib.lines import Line2D
from scipy import stats
warnings.filterwarnings("ignore")

# ============================================================
# Configuration
# ============================================================
MASTER_CSV = Path("E:/boshi/qm_harmonization_paper/results/biomarkers_master.csv")
EVAL_DIR   = Path("E:/boshi/qm_harmonization_paper/results/step4_eval")
DESIGN_C   = Path("E:/qm_harmonization_paper/derivatives/all-files/results_design_C.csv")
OUT_DIR    = Path("C:/Users/admin/WorkBuddy/2026-07-06-12-34-04/figures")
OUT_DIR.mkdir(parents=True, exist_ok=True)

BIOMARKERS = ["T2w_CSA", "GM_CSA_mm2", "MTR", "MTsat", "FA", "MD", "AD", "RD"]
BM_LABEL   = {"T2w_CSA": "T2w CSA", "GM_CSA_mm2": "T2* GM CSA",
              "MTR": "MTR", "MTsat": "MT$_{sat}$",
              "FA": "FA", "MD": "MD", "AD": "AD", "RD": "RD"}
METHOD_ORDER = ["Original", "LME", "ComBat", "ComBat-joint", "RELIEF", "CovBat"]
METHOD_SHORT = {"Original": "Orig", "LME": "LME", "ComBat": "CBT",
                "ComBat-joint": "CBT-J", "RELIEF": "RLF", "CovBat": "CVB"}

# Wong 2011 palette
WONG = {"blue":"#0072B2","sky":"#56B4E9","green":"#009E73","orange":"#E69F00",
        "vermilion":"#D55E00","pink":"#CC79A7","yellow":"#F0E442","grey":"#999999"}
VENDOR_COLOR = {"GE": WONG["blue"], "Philips": WONG["green"], "Siemens": WONG["vermilion"]}
METHOD_COLOR = {"LME": WONG["pink"], "ComBat": WONG["blue"],
                "ComBat-joint": WONG["sky"], "RELIEF": WONG["orange"],
                "CovBat": WONG["green"], "Original": WONG["grey"]}

# SCI-quality rcParams
mpl.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 9, "axes.titlesize": 10, "axes.labelsize": 9,
    "xtick.labelsize": 8, "ytick.labelsize": 8, "legend.fontsize": 7.5,
    "axes.linewidth": 0.8, "xtick.major.width": 0.8, "ytick.major.width": 0.8,
    "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
    "figure.dpi": 300, "savefig.dpi": 300,
})

# ============================================================
# Data loading
# ============================================================
def load_all():
    """Load master CSV + pre-computed eval results."""
    print("[1/4] Loading master CSV...")
    master = pd.read_csv(MASTER_CSV)
    master.rename(columns={"participant_id": "subject", "Manufacturer": "vendor",
                           "Age": "age", "Sex": "sex", "Pathology": "pathology"}, inplace=True)
    print(f"  Master: {len(master)} subjects, {master['site'].nunique()} sites, "
          f"{master['vendor'].nunique()} vendors")

    print("[2/4] Loading comparison_long (per-biomarker metrics)...")
    comp = pd.read_csv(EVAL_DIR / "comparison_long.csv")
    # Aggregate mean across biomarkers per (cohort, method)
    agg = comp.groupby(["cohort", "method"]).agg(
        vendor_eta2_mean=("vendor_eta2", "mean"),
        r_pre_post_mean=("r_pre_post", "mean"),
        age_R2_mean=("age_R2", "mean"),
        sex_R2_mean=("sex_R2", "mean"),
        vendor_eta2_std=("vendor_eta2", "std"),
        r_pre_post_std=("r_pre_post", "std"),
    ).reset_index()
    print(f"  Aggregated: {len(agg)} rows ({agg['cohort'].nunique()} cohorts × "
          f"{agg['method'].nunique()} methods)")

    print("[3/4] Loading permanova_results...")
    perm = pd.read_csv(EVAL_DIR / "permanova_results.csv")
    print(f"  PERMANOVA: {len(perm)} rows")

    print("[4/4] Loading UMAP embeddings...")
    umap = pd.read_csv(EVAL_DIR / "umap_embeddings.csv")
    print(f"  UMAP: {umap.shape[0]} rows, {umap['cohort'].nunique()} cohorts × "
          f"{umap['method'].nunique()} methods")

    return master, agg, perm, umap


def eta2(y, group):
    """Compute eta-squared (ANOVA effect size)."""
    grand = y.mean()
    ss_total = ((y - grand)**2).sum()
    if ss_total == 0:
        return 0.0
    ss_between = 0.0
    for g in np.unique(group):
        yg = y[group == g]
        ss_between += len(yg) * (yg.mean() - grand)**2
    return float(ss_between / ss_total)


# ============================================================
# Fig1a: Dataset & Biomarkers schematic
# ============================================================
def gen_fig1a(master):
    print("\n[Fig1a] Dataset & Biomarkers schematic")
    fig, axes = plt.subplots(1, 2, figsize=(8.0, 4.0),
                             gridspec_kw={"width_ratios": [1.2, 1]})
    axL, axR = axes

    pathology_counts = master["pathology"].value_counts()
    hc_n = pathology_counts.get("HC", 0)
    all_n = len(master)
    n_sites = master["site"].nunique()

    # -- Left panel: Cohort breakdown --
    axL.set_xlim(0, 10)
    axL.set_ylim(0, 10)
    axL.axis("off")
    axL.text(5, 9.8, "spine-generic multi-subject dataset", ha="center", fontsize=11, fontweight="bold")

    # Main box
    main = FancyBboxPatch((1.0, 1.0), 8.0, 8.2, boxstyle="round,pad=0.15",
                          facecolor="#E8F4FD", edgecolor="#0072B2", linewidth=1.5)
    axL.add_patch(main)
    axL.text(5, 9.2, f"N = {all_n} subjects, {n_sites} sites, 3 vendors",
             ha="center", fontsize=9, fontweight="bold", color="#0072B2")

    # HC sub-box
    hc_complete = master[master["pathology"] == "HC"][BIOMARKERS].notna().all(axis=1).sum()
    hc_vendor_counts = master[master["pathology"] == "HC"]["vendor"].value_counts()
    v_str = ", ".join([f"{v}={c}" for v, c in sorted(hc_vendor_counts.items())])

    hc_box = FancyBboxPatch((1.5, 5.8), 7.0, 2.8, boxstyle="round,pad=0.1",
                            facecolor="#D4EFDF", edgecolor="#009E73", linewidth=1.2)
    axL.add_patch(hc_box)
    axL.text(5, 8.0, f"HC cohort: N = {hc_n}", ha="center", fontsize=9,
             fontweight="bold", color="#009E73")
    axL.text(5, 7.2, f"Complete 8-biomarker: N = {hc_complete}",
             ha="center", fontsize=8, color="#333")
    axL.text(5, 6.5, f"Vendors: {v_str}", ha="center", fontsize=7, color="#555")

    # Full cohort sub-box
    all_complete = master[BIOMARKERS].notna().all(axis=1).sum()
    path_str = f"HC + MildCompression({pathology_counts.get('MildCompression', 0)}) + DCM({pathology_counts.get('DCM', 0)})"

    all_box = FancyBboxPatch((1.5, 2.5), 7.0, 2.8, boxstyle="round,pad=0.1",
                             facecolor="#FDEBD0", edgecolor="#E69F00", linewidth=1.2)
    axL.add_patch(all_box)
    axL.text(5, 4.7, f"Full cohort: N = {all_n}", ha="center", fontsize=9,
             fontweight="bold", color="#E69F00")
    axL.text(5, 3.9, f"Complete 8-biomarker: N = {all_complete}",
             ha="center", fontsize=8, color="#333")
    axL.text(5, 3.2, path_str, ha="center", fontsize=7, color="#555")

    # -- Right panel: Biomarker list --
    axR.set_xlim(0, 10)
    axR.set_ylim(0, 10)
    axR.axis("off")
    axR.text(5, 9.8, "8 qMRI Biomarkers (C2-C3)", ha="center", fontsize=11, fontweight="bold")

    bio_info = [
        ("T2w CSA", "SC cross-sectional area (T2w)", "CSA", WONG["blue"]),
        ("T2* GM CSA", "Gray matter CSA (T2*w)", "CSA", WONG["vermilion"]),
        ("MTR", "Magnetization transfer ratio", "MT", WONG["green"]),
        ("MTsat", "MT saturation", "MT", WONG["orange"]),
        ("FA", "Fractional anisotropy", "DTI", WONG["pink"]),
        ("MD", "Mean diffusivity", "DTI", WONG["sky"]),
        ("AD", "Axial diffusivity", "DTI", WONG["grey"]),
        ("RD", "Radial diffusivity", "DTI", WONG["yellow"]),
    ]

    for i, (name, desc, modality, color) in enumerate(bio_info):
        y_base = 8.5 - i * 1.0
        axR.add_patch(Rectangle((0.5, y_base - 0.25), 0.7, 0.6,
                                facecolor=color, edgecolor="black", linewidth=0.5))
        axR.text(1.5, y_base, name, fontsize=9, fontweight="bold", color=color)
        axR.text(1.5, y_base - 0.4, f"{desc} ({modality})", fontsize=7, color="#555")

    fig.tight_layout(pad=0.3, rect=[0, 0, 1, 0.97])
    for fmt in ["png", "pdf"]:
        fig.savefig(OUT_DIR / f"Fig1a_dataset_biomarkers.{fmt}", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("  Fig1a saved.")


# ============================================================
# Fig1b: Methods & Evaluation schematic
# ============================================================
def gen_fig1b():
    """Methods & Evaluation framework — clear 2-panel flow diagram."""
    print("[Fig1b] Methods & Evaluation schematic")
    fig, ax = plt.subplots(figsize=(14.0, 8.0))
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 100)
    ax.axis("off")

    # === Panel A: Harmonization pipeline =====================================
    ax.text(1, 96, "A", fontsize=14, fontweight="bold", va="top")
    ax.text(5, 96, "Harmonization pipeline: input data, five harmonization methods, and outputs",
            fontsize=11, fontweight="bold", va="top")

    # Input Data
    ix, iy, iw, ih = 1, 56, 11, 32
    ax.add_patch(FancyBboxPatch((ix, iy), iw, ih, boxstyle="round,pad=0.4",
                                facecolor="#F0F0F0", edgecolor="#888", linewidth=1.5))
    ax.text(ix + iw/2, iy + ih - 2.5, "Input Data",
            ha="center", va="top", fontsize=9, fontweight="bold", color="#333")
    ax.text(ix + iw/2, iy + ih/2, "8 biomarkers\nC2-C3 aggregated\nN = 267 subjects\n44 sites / 3 vendors",
            ha="center", va="center", fontsize=7.5, color="#444", linespacing=1.4)

    # Pipeline arrows between Input and LME
    ax.annotate("", xy=(13, 72), xytext=(12, 72),
                arrowprops=dict(arrowstyle="->", color="#555", lw=2.0))

    # Method boxes — all five methods in parallel
    methods = [
        ("LME", "Linear Mixed Effects\n(site random intercept)", WONG["pink"]),
        ("ComBat", "Empirical Bayes ComBat\n(per metric, no EB)", WONG["blue"]),
        ("ComBat-\njoint", "Joint EB ComBat\n(all 8 metrics)", WONG["sky"]),
        ("RELIEF", "ICA + BH-FDR\ncomponent removal", WONG["orange"]),
        ("CovBat", "PCA + per-PC ComBat\n(covariance adjustment)", WONG["green"]),
    ]
    m_w, m_h, m_y = 12.5, 32, 56
    m_gap = 1.5
    method_positions = []
    for i, (name, desc, color) in enumerate(methods):
        x0 = 13 + i * (m_w + m_gap)
        method_positions.append((x0, m_w))
        # Semi-transparent background
        ax.add_patch(FancyBboxPatch((x0, m_y), m_w, m_h, boxstyle="round,pad=0.4",
                                    facecolor=color, edgecolor="white", alpha=0.12))
        # Colored border
        ax.add_patch(FancyBboxPatch((x0, m_y), m_w, m_h, boxstyle="round,pad=0.4",
                                    facecolor="none", edgecolor=color, linewidth=2.0))
        # Method name
        ax.text(x0 + m_w/2, m_y + m_h - 2.5, name,
                ha="center", va="top", fontsize=9, fontweight="bold", color=color,
                linespacing=1.1)
        # Method description
        ax.text(x0 + m_w/2, m_y + m_h/2 - 1.0, desc,
                ha="center", va="center", fontsize=7.5, color="#333",
                linespacing=1.4)

        # Arrow to next method
        if i < len(methods) - 1:
            ax.annotate("", xy=(x0 + m_w + m_gap - 0.3, 72),
                        xytext=(x0 + m_w + 0.3, 72),
                        arrowprops=dict(arrowstyle="->", color="#AAA", lw=1.5))

    # Last arrow to Output
    last_x, last_w = method_positions[-1]
    ax.annotate("", xy=(83, 72), xytext=(last_x + last_w + 0.5, 72),
                arrowprops=dict(arrowstyle="->", color="#555", lw=2.0))

    # Output Data
    ox, oy, ow, oh = 83, 56, 11, 32
    ax.add_patch(FancyBboxPatch((ox, oy), ow, oh, boxstyle="round,pad=0.4",
                                facecolor="#E8F5E8", edgecolor="#009E73", linewidth=1.5))
    ax.text(ox + ow/2, oy + oh - 2.5, "Output Data",
            ha="center", va="top", fontsize=9, fontweight="bold", color="#007050")
    ax.text(ox + ow/2, oy + oh/2, "5 harmonized\n8-biomarker\ndatasets",
            ha="center", va="center", fontsize=7.5, color="#444", linespacing=1.4)

    # === Connection: pipeline -> evaluation =====================================
    ax.plot([5, 95], [53, 53], color="#888", linewidth=0.8, linestyle="--")

    ax.text(50, 48, "All five harmonized outputs are evaluated by four orthogonal endpoints",
            ha="center", va="center", fontsize=8.5, color="#333", fontweight="bold")
    ax.annotate("", xy=(50, 44), xytext=(50, 46.5),
                arrowprops=dict(arrowstyle="->", color="#555", lw=2.0))

    # === Panel B: Evaluation endpoints =========================================
    ax.text(1, 42, "B", fontsize=14, fontweight="bold", va="top")
    ax.text(5, 42, "Orthogonal evaluation endpoints",
            fontsize=11, fontweight="bold", va="top")

    evals = [
        ("Vendor \u03b7\u00b2", "Univariate vendor\neffect removal", WONG["vermilion"]),
        ("PERMANOVA R\u00b2", "Multivariate vendor\neffect removal", WONG["blue"]),
        ("r_pre,post", "Per-subject rank\npreservation", WONG["green"]),
        ("Age + Sex R\u00b2", "Biological signal\nrecovery", WONG["orange"]),
    ]
    ev_w, ev_h, ev_y = 18, 28, 5
    ev_gap = 4
    for i, (name, desc, color) in enumerate(evals):
        x0 = 5 + i * (ev_w + ev_gap)
        ax.add_patch(FancyBboxPatch((x0, ev_y), ev_w, ev_h, boxstyle="round,pad=0.4",
                                    facecolor=color, edgecolor=color, alpha=0.08, linewidth=1.5))
        ax.annotate("", xy=(x0 + ev_w/2, ev_y + ev_h),
                    xytext=(x0 + ev_w/2, 44),
                    arrowprops=dict(arrowstyle="->", color=color, lw=1.2))
        ax.text(x0 + ev_w/2, ev_y + ev_h - 2.5, name,
                ha="center", va="top", fontsize=9, fontweight="bold", color=color)
        ax.text(x0 + ev_w/2, ev_y + ev_h/2 - 0.5, desc,
                ha="center", va="center", fontsize=8, color="#333", linespacing=1.45)

    for fmt in ["png", "pdf"]:
        fig.savefig(OUT_DIR / f"Fig1b_methods_evaluation.{fmt}", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("  Fig1b saved.")



# ============================================================
# Fig2: Baseline vendor effects
# ============================================================
def gen_fig2(master, agg, perm, umap):
    print("[Fig2] Baseline vendor effects")
    hc = master[master["pathology"] == "HC"].copy()

    fig = plt.figure(figsize=(8.0, 7.5))
    gs = fig.add_gridspec(2, 2, hspace=0.75, wspace=0.45)

    # -- A: Representative biomarker boxplots --
    axA = fig.add_subplot(gs[0, 0])
    axA.text(-0.15, 1.18, "A", transform=axA.transAxes, fontsize=12, fontweight="bold", va="top")
    axA.text(0.02, 1.08, "Per-vendor distribution of 4 biomarkers (HC)",
             transform=axA.transAxes, fontsize=9, fontweight="bold", va="top")
    axA.axis("off")

    vendors_list = ["Siemens", "GE", "Philips"]
    vendor_short = ["Si", "GE", "Ph"]
    sub_bios = ["T2w_CSA", "GM_CSA_mm2", "FA", "MTR"]
    sub_names = ["T2w CSA", "T2* GM CSA", "FA", "MTR"]
    colors_v = [VENDOR_COLOR[v] for v in vendors_list]

    for idx, (bio, bname) in enumerate(zip(sub_bios, sub_names)):
        bx_pos = [0.11 + idx * 0.22, 0.12, 0.19, 0.72]  # [left, bottom, width, height]
        ax_sub = fig.add_axes(bx_pos)
        data_list = [hc[hc["vendor"] == v][bio].dropna().values for v in vendors_list]

        bp = ax_sub.boxplot(data_list, positions=[0, 1, 2], widths=0.5,
                           patch_artist=True, showfliers=False,
                           medianprops={"color": "black", "linewidth": 0.8},
                           whiskerprops={"linewidth": 0.6},
                           capprops={"linewidth": 0.6})
        for patch, c in zip(bp["boxes"], colors_v):
            patch.set_facecolor(c)
            patch.set_alpha(0.7)

        ax_sub.set_title(bname, fontsize=7.5, fontweight="bold")
        ax_sub.set_xticks([0, 1, 2])
        ax_sub.set_xticklabels(vendor_short, fontsize=6)
        ax_sub.tick_params(labelsize=6)

        # ANOVA p-value
        try:
            fv, pv = stats.f_oneway(*data_list)
            label = f"p={pv:.3f}" if pv >= 0.001 else f"p={pv:.2e}"
            ax_sub.text(0.98, 0.05, label, transform=ax_sub.transAxes,
                       ha="right", fontsize=6, color="#555")
        except:
            pass

    # -- B: UMAP --
    axB = fig.add_subplot(gs[0, 1])
    axB.text(-0.15, 1.18, "B", transform=axB.transAxes, fontsize=12, fontweight="bold", va="top")
    axB.text(0.02, 1.08, "UMAP of standardized 8-biomarker vectors (HC, Original)",
             transform=axB.transAxes, fontsize=9, fontweight="bold", va="top")

    umap_orig = umap[(umap["method"] == "Original") & (umap["cohort"] == "HC")]
    for v in vendors_list:
        mask = umap_orig["vendor"] == v
        axB.scatter(umap_orig.loc[mask, "x"], umap_orig.loc[mask, "y"],
                   s=18, c=VENDOR_COLOR[v], label=v, alpha=0.85,
                   edgecolors="white", linewidths=0.3)
    axB.set_xlabel("UMAP1", fontsize=8)
    axB.set_ylabel("UMAP2", fontsize=8)
    axB.legend(fontsize=7, loc="lower left", markerscale=0.8)
    axB.tick_params(labelsize=7)

    # -- C: Vendor η² text --
    axC = fig.add_subplot(gs[1, 0])
    axC.axis("off")
    axC.text(-0.15, 1.18, "C", transform=axC.transAxes, fontsize=12, fontweight="bold", va="top")
    axC.text(-0.02, 1.05, "Key Observations", transform=axC.transAxes,
             fontsize=9, fontweight="bold", va="top")

    obs_lines = []
    for bio in BIOMARKERS:
        vals = hc[bio].dropna()
        groups = hc.loc[vals.index, "vendor"].astype(str).values
        e2 = eta2(vals.values, groups)
        obs_lines.append(f"  {BM_LABEL[bio]:14s} η² = {e2:.3f}")

    # Get PERMANOVA R² for Original
    perm_orig = perm[(perm["method"] == "Original") & (perm["cohort"] == "HC")]
    pr2 = perm_orig["R2"].values[0] if len(perm_orig) > 0 else 0
    p_perm = perm_orig["p"].values[0] if len(perm_orig) > 0 else 1

    axC.text(0.02, 0.88, "\n".join(obs_lines), transform=axC.transAxes,
             fontsize=7.5, fontfamily="monospace", va="top", color="#333", linespacing=1.4)
    axC.text(0.02, 0.18, f"Multivariate PERMANOVA R² = {pr2:.3f} (p = {p_perm:.4f})",
             transform=axC.transAxes, fontsize=8, fontweight="bold", color="#D55E00")
    axC.text(0.02, 0.08, "T2* GM CSA and MTsat show largest vendor effects (η² > 0.15)",
             transform=axC.transAxes, fontsize=7.5, color="#555")

    # -- D: Baseline η² bar chart --
    axD = fig.add_subplot(gs[1, 1])
    axD.text(-0.15, 1.18, "D", transform=axD.transAxes, fontsize=12, fontweight="bold", va="top")
    axD.text(-0.02, 1.08, "Baseline vendor η², with PERMANOVA R² reference",
             transform=axD.transAxes, fontsize=9, fontweight="bold", va="top")

    eta_vals = []
    for bio in BIOMARKERS:
        vals = hc[bio].dropna()
        groups = hc.loc[vals.index, "vendor"].values
        eta_vals.append(eta2(vals.values, groups))

    x_pos = np.arange(len(BIOMARKERS))
    bars = axD.bar(x_pos, eta_vals, color=[VENDOR_COLOR["Siemens"]]*len(BIOMARKERS),
                  edgecolor="black", linewidth=0.5, alpha=0.7)
    axD.axhline(y=pr2, color="#D55E00", linestyle="--", linewidth=1.2,
               label=f"PERMANOVA R² = {pr2:.3f}")
    axD.set_xticks(x_pos)
    axD.set_xticklabels([BM_LABEL[b] for b in BIOMARKERS], rotation=30, ha="right", fontsize=7)
    axD.set_ylabel("Vendor η²", fontsize=8)
    axD.legend(fontsize=7, loc="upper right")
    axD.tick_params(labelsize=7)

    for xi, vi in zip(x_pos, eta_vals):
        axD.text(xi, vi + 0.01, f"{vi:.3f}", ha="center", fontsize=6.5, color="#333")

    axD.set_title("Baseline vendor effects — unharmonized qMRI biomarkers",
                  fontsize=9, fontweight="bold")

    for fmt in ["png", "pdf"]:
        fig.savefig(OUT_DIR / f"Fig2_baseline_effects.{fmt}", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("  Fig2 saved.")


# ============================================================
# Fig3: Harmonization performance (grouped bar)
# ============================================================
def gen_fig3(agg, perm):
    print("[Fig3] Harmonization performance")

    metrics_plot = ["vendor_eta2_mean", "permanova_R2", "r_pre_post_mean", "age_sex_R2"]
    metric_labels = ["Vendor η² (mean across biomarkers)", "PERMANOVA R²",
                     "r_pre,post (mean across biomarkers)", "Δ Age + Sex R²"]
    cohorts = ["HC", "ALL"]

    methods_plot = [m for m in METHOD_ORDER[1:] if m in agg["method"].values]
    colors_c = {"HC": WONG["blue"], "ALL": WONG["orange"]}

    fig, axes = plt.subplots(2, 2, figsize=(8.0, 6.0))
    axes = axes.flatten()

    for idx, (metric, mlabel) in enumerate(zip(metrics_plot, metric_labels)):
        ax = axes[idx]
        ax.text(-0.15, 1.05, chr(65 + idx), transform=ax.transAxes,
                fontsize=11, fontweight="bold", va="top")

        x = np.arange(len(methods_plot))
        width = 0.32

        for c_idx, cohort in enumerate(cohorts):
            vals = []
            if metric == "permanova_R2":
                # From permanova_results
                for m in methods_plot:
                    row = perm[(perm["method"] == m) & (perm["cohort"] == cohort)]
                    vals.append(row["R2"].values[0] if len(row) > 0 else 0)
            elif metric == "age_sex_R2":
                # Sum of age_R2 + sex_R2 from agg
                sub = agg[agg["cohort"] == cohort]
                for m in methods_plot:
                    row = sub[sub["method"] == m]
                    if len(row) > 0:
                        vals.append(row["age_R2_mean"].values[0] + row["sex_R2_mean"].values[0])
                    else:
                        vals.append(0)
                # Also get Original
                orig = agg[(agg["cohort"] == cohort) & (agg["method"] == "Original")]
                bl = orig["age_R2_mean"].values[0] + orig["sex_R2_mean"].values[0] if len(orig) > 0 else 0
            else:
                sub = agg[agg["cohort"] == cohort]
                for m in methods_plot:
                    row = sub[sub["method"] == m]
                    vals.append(row[metric].values[0] if len(row) > 0 else 0)

            bars = ax.bar(x + c_idx * width - width/2, vals, width,
                         label=cohort, color=colors_c[cohort],
                         edgecolor="black", linewidth=0.3, alpha=0.85)

            # Annotate
            for xi, vi in zip(x, vals):
                if abs(vi) > 0.0005:
                    rot = 90 if metric in ["permanova_R2", "vendor_eta2_mean"] else 0
                    ax.text(xi + c_idx * width - width/2, vi + 0.002,
                           f"{vi:.3f}", ha="center", fontsize=6, color="#333", rotation=rot)

        # Baseline line for eta2 and permanova
        if metric in ["vendor_eta2_mean", "permanova_R2"]:
            if metric == "permanova_R2":
                bl_row = perm[(perm["method"] == "Original") & (perm["cohort"] == "HC")]
                bl = bl_row["R2"].values[0] if len(bl_row) > 0 else 0
            else:
                bl_row = agg[(agg["method"] == "Original") & (agg["cohort"] == "HC")]
                bl = bl_row[metric].values[0] if len(bl_row) > 0 else 0
            ax.axhline(y=bl, color="#999", linestyle="--", linewidth=0.8)
            ax.text(len(methods_plot)-0.3, bl, "baseline", fontsize=6.5, color="#999", va="bottom")

        ax.set_xticks(x)
        ax.set_xticklabels([METHOD_SHORT.get(m, m) for m in methods_plot],
                          fontsize=7, rotation=15, ha="right")
        ax.set_title(mlabel, fontsize=8.5, fontweight="bold")
        ax.tick_params(labelsize=7)
        if idx == 0:
            ax.legend(fontsize=7, loc="upper right", ncol=2, framealpha=0.8)

    fig.suptitle("Harmonization performance across orthogonal evaluation endpoints",
                 fontsize=11, fontweight="bold", y=0.99)
    fig.tight_layout(pad=0.4, rect=[0, 0.04, 1, 0.95])

    for fmt in ["png", "pdf"]:
        fig.savefig(OUT_DIR / f"Fig3_performance.{fmt}", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("  Fig3 saved.")


# ============================================================
# Fig4: Decontamination vs preservation trade-off
# ============================================================
def gen_fig4(agg, perm):
    print("[Fig4] Trade-off scatter")

    fig = plt.figure(figsize=(8.5, 5.8))
    gs = fig.add_gridspec(2, 2, height_ratios=[12, 1.8],
                          hspace=0.55, wspace=0.40)

    cohorts = ["HC", "ALL"]
    methods_plot = [m for m in METHOD_ORDER[1:] if m in agg["method"].values]

    for idx, cohort in enumerate(cohorts):
        ax = fig.add_subplot(gs[0, idx])
        ax.text(-0.15, 1.05, chr(65 + idx), transform=ax.transAxes,
                fontsize=11, fontweight="bold", va="top")
        ax.set_title(f"{cohort} cohort", fontsize=9, fontweight="bold")

        sub_agg = agg[agg["cohort"] == cohort]
        sub_perm = perm[perm["cohort"] == cohort]

        for m in methods_plot:
            # PERMANOVA R² from permanova_results
            pr = sub_perm[sub_perm["method"] == m]
            r2 = pr["R2"].values[0] if len(pr) > 0 else 0

            # r_pre_post from agg
            ar = sub_agg[sub_agg["method"] == m]
            rpp = ar["r_pre_post_mean"].values[0] if len(ar) > 0 else 0

            x_log = -np.log10(max(r2, 1e-5))
            ax.scatter(x_log, rpp, s=120, c=METHOD_COLOR.get(m, "#999"),
                      edgecolors="black", linewidths=0.6, zorder=3,
                      label=METHOD_SHORT.get(m, m))

            # Label with offset
            dx = 0.08
            ax.annotate(METHOD_SHORT.get(m, m), xy=(x_log, rpp),
                       xytext=(x_log + dx, rpp + 0.006),
                       fontsize=7.5, fontweight="bold", color=METHOD_COLOR.get(m, "#333"))

        ax.set_xlabel("–log₁₀(R²_PERMANOVA)", fontsize=8)
        ax.set_ylabel("r_pre,post", fontsize=8)
        ax.tick_params(labelsize=7)

        # Add Original reference
        orig_perm = sub_perm[sub_perm["method"] == "Original"]
        orig_agg = sub_agg[sub_agg["method"] == "Original"]
        if len(orig_perm) > 0 and len(orig_agg) > 0:
            r2_o = orig_perm["R2"].values[0]
            rpp_o = orig_agg["r_pre_post_mean"].values[0]
            ax.axhline(y=rpp_o, color="#999", linestyle=":", linewidth=0.8, alpha=0.6)
            ax.text(0.98, rpp_o, " r=1.0", transform=ax.get_yaxis_transform(),
                   fontsize=6.5, color="#999", va="center", ha="left")

    # Legend row
    ax_leg = fig.add_subplot(gs[1, :])
    ax_leg.axis("off")
    handles = []
    for m in methods_plot:
        handles.append(Line2D([0], [0], marker="o", color="w",
                              markerfacecolor=METHOD_COLOR[m], markersize=7,
                              label=m, markeredgecolor="black", markeredgewidth=0.4))
    ax_leg.legend(handles=handles, loc="center", ncol=len(methods_plot),
                  fontsize=8, frameon=False)

    for fmt in ["png", "pdf"]:
        fig.savefig(OUT_DIR / f"Fig4_tradeoff.{fmt}", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("  Fig4 saved.")


# ============================================================
# Fig5: r_pre,post - per-biomarker bar chart
# ============================================================
def gen_fig5(agg, comp):
    """Show r_pre,post per biomarker per method (bar chart)."""
    print("[Fig5] r_pre,post per biomarker")

    fig, axes = plt.subplots(2, 2, figsize=(10.0, 7.5))
    axes = axes.flatten()
    methods_plot = [m for m in METHOD_ORDER[1:] if m in comp["method"].values]

    panels = [
        ("HC", "4 representative biomarkers", ["T2w_CSA", "GM_CSA_mm2", "FA", "MTR"]),
        ("HC", "4 DTI biomarkers", ["FA", "MD", "AD", "RD"]),
        ("ALL", "4 representative biomarkers", ["T2w_CSA", "GM_CSA_mm2", "FA", "MTR"]),
        ("ALL", "4 DTI biomarkers", ["FA", "MD", "AD", "RD"]),
    ]

    for idx, (cohort, title, bios) in enumerate(panels):
        ax = axes[idx]
        ax.text(-0.15, 1.06, chr(65 + idx), transform=ax.transAxes,
                fontsize=11, fontweight="bold", va="top")

        n_bios = len(bios)
        n_methods = len(methods_plot)
        x = np.arange(n_bios)
        width = 0.7 / n_methods

        for j, m in enumerate(methods_plot):
            vals = []
            for bio in bios:
                row = comp[(comp["cohort"] == cohort) & (comp["method"] == m) & (comp["biomarker"] == bio)]
                vals.append(row["r_pre_post"].values[0] if len(row) > 0 else 0)

            bars = ax.bar(x + (j - n_methods/2 + 0.5) * width, vals, width,
                         color=METHOD_COLOR.get(m, "#999"), alpha=0.85,
                         edgecolor="white", linewidth=0.3)

        ax.set_xticks(x)
        ax.set_xticklabels([BM_LABEL.get(b, b) for b in bios], fontsize=7.5)
        ax.set_ylabel("r_pre,post", fontsize=8)
        ax.set_title(f"{cohort} — {title}", fontsize=8.5, fontweight="bold")
        ax.tick_params(labelsize=7)
        ax.set_ylim(0.4, 1.05)

        # Add legend only to top-left
        if idx == 0:
            leg_handles = [Patch(facecolor=METHOD_COLOR[m], label=METHOD_SHORT.get(m, m),
                                edgecolor="white", linewidth=0.3) for m in methods_plot]
            ax.legend(handles=leg_handles, fontsize=6.5, loc="lower left",
                     ncol=3, framealpha=0.8)

    fig.suptitle("Per-subject biomarker preservation (r_pre,post) after harmonization",
                 fontsize=11, fontweight="bold", y=0.99)
    fig.tight_layout(pad=0.4, rect=[0, 0.02, 1, 0.96])

    for fmt in ["png", "pdf"]:
        fig.savefig(OUT_DIR / f"Fig5_rpre_post_marginals.{fmt}", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("  Fig5 saved.")


# ============================================================
# FigS: Design C robustness
# ============================================================
def gen_figs():
    print("[FigS] Design C robustness")

    if not DESIGN_C.exists():
        print(f"  WARNING: {DESIGN_C} not found, creating placeholder")
        fig, axes = plt.subplots(1, 3, figsize=(9, 3.5))
        for ax in axes:
            ax.text(0.5, 0.5, "Design C data not available", ha="center", va="center", fontsize=9)
        fig.tight_layout()
        for fmt in ["png", "pdf"]:
            fig.savefig(OUT_DIR / f"FigS_design_C_robustness.{fmt}", dpi=300, bbox_inches="tight")
        plt.close(fig)
        return

    df = pd.read_csv(DESIGN_C)
    print(f"  Design C: {len(df)} rows, columns: {list(df.columns)}")

    # Identify metric columns and method column
    if "method" in df.columns:
        df["method_clean"] = df["method"]
    elif any(df.columns.str.contains("method", case=False)):
        mcol = [c for c in df.columns if "method" in c.lower()][0]
        df["method_clean"] = df[mcol]
    else:
        print("  Cannot identify method column, skipping")
        return

    methods_here = [m for m in METHOD_ORDER if m in df["method_clean"].values]

    # Find metric columns (numeric, likely to be metrics)
    numeric_cols = df.select_dtypes(include=[np.number]).columns.tolist()
    metric_cols = [c for c in numeric_cols if c not in ["seed", "n_perm", "N"]]

    fig, axes = plt.subplots(1, 3, figsize=(9.0, 3.5))
    axes = axes.flatten()
    panel_labels = ["A", "B", "C"]

    if len(metric_cols) >= 3:
        metric_cols = metric_cols[:3]
    metrics = metric_cols[:3] if metric_cols else ["vendor_eta2", "permanova_r2", "r_pre_post"]
    metric_labels_map = {"vendor_eta2": "Vendor η²", "permanova_r2": "PERMANOVA R²",
                         "r_pre_post": "r_pre,post", "age_R2": "Age R²", "sex_R2": "Sex R²"}

    for idx, metric in enumerate(metrics):
        ax = axes[idx]
        ax.text(-0.12, 1.05, panel_labels[idx], transform=ax.transAxes,
                fontsize=11, fontweight="bold", va="top")

        if metric in df.columns:
            data_by_method = [df[df["method_clean"] == m][metric].dropna().values
                            for m in methods_here if m in df["method_clean"].values]

            if data_by_method:
                labels = [m for m in methods_here if m in df["method_clean"].values]
                bp = ax.boxplot(data_by_method, tick_labels=[METHOD_SHORT.get(m, m) for m in labels],
                               patch_artist=True, showfliers=True, widths=0.55,
                               flierprops={"markersize": 2, "alpha": 0.4})
                for patch, m in zip(bp["boxes"], labels):
                    patch.set_facecolor(METHOD_COLOR.get(m, "#999"))
                    patch.set_alpha(0.6)

                # Also strip plots
                for i, vals in enumerate(data_by_method):
                    jitter = np.random.normal(0, 0.06, len(vals))
                    ax.scatter(np.full(len(vals), i+1) + jitter, vals,
                              s=3, alpha=0.3, color="black", zorder=3)

        ax.set_title(metric_labels_map.get(metric, metric), fontsize=9, fontweight="bold")
        ax.tick_params(labelsize=7)

    fig.suptitle("Design C: Robustness to balanced sub-sampling (20 random seeds)",
                 fontsize=10, fontweight="bold", y=1.02)
    fig.tight_layout(pad=0.4, rect=[0, 0, 1, 0.94])

    for fmt in ["png", "pdf"]:
        fig.savefig(OUT_DIR / f"FigS_design_C_robustness.{fmt}", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("  FigS saved.")


# ============================================================
# Main
# ============================================================
if __name__ == "__main__":
    print("=" * 60)
    print("Regenerating ALL 7 Figures (v2 — from pre-computed data)")
    print("=" * 60)

    master, agg, perm, umap = load_all()

    gen_fig1a(master)
    gen_fig1b()
    gen_fig2(master, agg, perm, umap)
    gen_fig3(agg, perm)
    gen_fig4(agg, perm)
    gen_fig5(agg, pd.read_csv(EVAL_DIR / "comparison_long.csv"))
    gen_figs()

    print("\n" + "=" * 60)
    print(f"DONE! All 7 figures saved to: {OUT_DIR}")
    for f in sorted(OUT_DIR.glob("*.png")):
        size_kb = f.stat().st_size / 1024
        print(f"  {f.name:45s} {size_kb:6.0f} KB")
