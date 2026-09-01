#!/usr/bin/env python3
"""
High-quality SCI figure generation for spine-generic harmonization paper.

Fig1a/Fig1b: SVG schematics (draw.io style) via generate_svg_figures.py
Fig2-Fig5, FigS: matplotlib statistical charts with REAL data

Data sources:
  E:/boshi/qm_harmonization_paper/results/biomarkers_master.csv
  E:/boshi/qm_harmonization_paper/results/step4_eval/comparison_long.csv
  E:/boshi/qm_harmonization_paper/results/step4_eval/permanova_results.csv
  E:/boshi/qm_harmonization_paper/results/step4_eval/umap_embeddings.csv
  E:/qm_harmonization_paper/derivatives/all-files/results_design_C.csv
"""
from __future__ import annotations
import warnings
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
from matplotlib.patches import FancyBboxPatch, Patch, Rectangle
from matplotlib.lines import Line2D
from matplotlib.colors import LinearSegmentedColormap
from scipy import stats

warnings.filterwarnings("ignore")

# ============================================================
# Paths
# ============================================================
MASTER_CSV = Path("E:/boshi/qm_harmonization_paper/results/biomarkers_master.csv")
EVAL_DIR   = Path("E:/boshi/qm_harmonization_paper/results/step4_eval")
DESIGN_C   = Path("E:/qm_harmonization_paper/derivatives/all-files/results_design_C.csv")
OUT_DIR    = Path("C:/Users/admin/WorkBuddy/2026-07-06-12-34-04/figures")
OUT_DIR.mkdir(parents=True, exist_ok=True)

# ============================================================
# Constants
# ============================================================
BIOMARKERS = ["T2w_CSA", "GM_CSA_mm2", "MTR", "MTsat", "FA", "MD", "AD", "RD"]
BM_LABEL = {
    "T2w_CSA": "T2w CSA", "GM_CSA_mm2": "T2* GM CSA",
    "MTR": "MTR", "MTsat": "MT$_{sat}$",
    "FA": "FA", "MD": "MD", "AD": "AD", "RD": "RD",
}
BM_FAMILY = {
    "T2w_CSA": "T2w", "GM_CSA_mm2": "T2*",
    "MTR": "MT", "MTsat": "MT",
    "FA": "DTI", "MD": "DTI", "AD": "DTI", "RD": "DTI",
}
METHOD_ORDER = ["Original", "LME", "ComBat", "ComBat-joint", "RELIEF", "CovBat"]
METHOD_LABEL = {
    "Original": "Original", "LME": "LME/FE", "ComBat": "ComBat",
    "ComBat-joint": "ComBat-joint", "RELIEF": "RELIEF", "CovBat": "CovBat",
}

# Wong 2011 color-blind-safe palette
WONG = {
    "blue": "#0072B2", "sky": "#56B4E9", "green": "#009E73",
    "orange": "#E69F00", "vermilion": "#D55E00", "pink": "#CC79A7",
    "yellow": "#F0E442", "grey": "#999999",
}
VENDOR_COLOR = {"GE": WONG["blue"], "Philips": WONG["green"], "Siemens": WONG["vermilion"]}
METHOD_COLOR = {
    "LME": WONG["pink"], "ComBat": WONG["blue"],
    "ComBat-joint": WONG["sky"], "RELIEF": WONG["orange"],
    "CovBat": WONG["green"], "Original": WONG["grey"],
}
FAMILY_COLOR = {
    "T2w": WONG["blue"], "T2*": WONG["vermilion"],
    "MT": WONG["green"], "DTI": WONG["orange"],
}

# SCI-quality global style
mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 9,
    "axes.titlesize": 10,
    "axes.labelsize": 9,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "legend.fontsize": 8,
    "axes.linewidth": 0.8,
    "xtick.major.width": 0.8,
    "ytick.major.width": 0.8,
    "xtick.major.size": 3,
    "ytick.major.size": 3,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "svg.fonttype": "none",
    "figure.dpi": 300,
    "savefig.dpi": 300,
    "axes.spines.top": False,
    "axes.spines.right": False,
})


# ============================================================
# Data loading
# ============================================================
def load_all():
    master = pd.read_csv(MASTER_CSV)
    master.rename(columns={
        "participant_id": "subject", "Manufacturer": "vendor",
        "Age": "age", "Sex": "sex", "Pathology": "pathology",
    }, inplace=True)
    master["cohort"] = master["pathology"].apply(lambda x: "HC" if x == "HC" else "Path")
    comp = pd.read_csv(EVAL_DIR / "comparison_long.csv")
    perm = pd.read_csv(EVAL_DIR / "permanova_results.csv")
    umap = pd.read_csv(EVAL_DIR / "umap_embeddings.csv")
    return master, comp, perm, umap


# ============================================================
# Fig2: Baseline vendor effects
# ============================================================
def gen_fig2(master, comp, perm, umap):
    """Baseline vendor effects before harmonization.

    Layout (GridSpec nested, NO manual add_axes):
      Row 0: [Panel A: 4 violins | Panel B: 2 UMAPs]
      Row 1: [Panel C: full-width eta2 bar chart]
    """
    print("[Fig2] Baseline vendor effects")
    vendors = ["GE", "Philips", "Siemens"]
    v_colors = [VENDOR_COLOR[v] for v in vendors]
    # Complete-case HC analysis sample (N=188): all 8 biomarkers + age + sex + vendor
    # present — identical filter to step3_refit_completecase.py / step4_eval_completecase.py,
    # so violins / ANOVA p / eta2 annotations share the same subjects as Table 2.
    CC_BIOMARKERS = ["T2w_CSA", "GM_CSA_mm2", "MTR", "MTsat",
                     "FA", "MD", "AD", "RD"]
    hc = master[master["pathology"] == "HC"].dropna(
        subset=CC_BIOMARKERS + ["age", "sex", "vendor"])
    hc_orig = comp[(comp["cohort"] == "HC") & (comp["method"] == "Original")].set_index("biomarker")

    fig = plt.figure(figsize=(13, 9))

    # Master grid: 2 rows, 1 col
    outer = GridSpec(2, 1, figure=fig, height_ratios=[1.0, 0.85], hspace=0.38)

    # Row 0: split into Panel A (left, 4 violins) and Panel B (right, 2 UMAPs)
    row0 = GridSpecFromSubplotSpec(1, 2, subplot_spec=outer[0, 0],
                                   width_ratios=[4, 3], wspace=0.25)

    # Panel A: 4 violin subplots
    gs_a = GridSpecFromSubplotSpec(1, 4, subplot_spec=row0[0, 0], wspace=0.55)
    rep_bios = ["T2w_CSA", "FA", "MTsat", "AD"]
    rep_names = ["T2w CSA", "FA", "MT$_{sat}$", "AD"]

    for idx, (bio, bname) in enumerate(zip(rep_bios, rep_names)):
        ax = fig.add_subplot(gs_a[0, idx])
        data_list = [hc[hc["vendor"] == v][bio].dropna().values for v in vendors]

        parts = ax.violinplot(data_list, positions=[0, 1, 2], widths=0.7, showextrema=False)
        for pc, c in zip(parts["bodies"], v_colors):
            pc.set_facecolor(c)
            pc.set_alpha(0.55)
            pc.set_edgecolor(c)
            pc.set_linewidth(0.8)

        bp = ax.boxplot(data_list, positions=[0, 1, 2], widths=0.25,
                        patch_artist=True, showfliers=False,
                        medianprops={"color": "black", "linewidth": 0.8},
                        whiskerprops={"linewidth": 0.6},
                        capprops={"linewidth": 0.6},
                        boxprops={"linewidth": 0.5})
        for patch, c in zip(bp["boxes"], v_colors):
            patch.set_facecolor(c)
            patch.set_alpha(0.8)

        ax.set_xticks([0, 1, 2])
        ax.set_xticklabels(["GE", "Phil", "Siem"], fontsize=7)
        ax.tick_params(labelsize=7)
        eta2_val = hc_orig.loc[bio, "vendor_eta2"] if bio in hc_orig.index else 0
        ax.set_title(f"{bname}\n($\\eta^2$={eta2_val:.2f})", fontsize=8, fontweight="bold", pad=6)
        try:
            _, pv = stats.f_oneway(*data_list)
            p_str = f"p={pv:.1e}" if pv < 0.001 else f"p={pv:.3f}"
            ax.text(0.97, 0.05, p_str, transform=ax.transAxes,
                    ha="right", va="bottom", fontsize=7, color="#444",
                    bbox=dict(boxstyle="round,pad=0.12", facecolor="white", alpha=0.8, edgecolor="none"))
        except Exception:
            pass
        if idx == 0:
            ax.set_ylabel("Value (a.u.)", fontsize=8)
        else:
            # Hide redundant y-axis tick labels on inner subplots
            ax.set_yticklabels([])

    # Panel A title + letter
    fig.text(0.01, 0.96, "A", fontsize=13, fontweight="bold", va="top")
    fig.text(0.03, 0.96, "Per-vendor distributions of representative biomarkers (HC)",
             fontsize=10, fontweight="bold", va="top")

    # Panel B: 2 UMAP subplots
    gs_b = GridSpecFromSubplotSpec(1, 2, subplot_spec=row0[0, 1], wspace=0.40)
    for ci, (cohort, sub_title) in enumerate([("HC", "HC"), ("ALL", "ALL")]):
        ax = fig.add_subplot(gs_b[0, ci])
        umap_coh = umap[(umap["method"] == "Original") & (umap["cohort"] == cohort)]
        for v in vendors:
            mask = umap_coh["vendor"] == v
            ax.scatter(umap_coh.loc[mask, "x"], umap_coh.loc[mask, "y"],
                       s=12, c=VENDOR_COLOR[v], label=v if ci == 0 else None,
                       alpha=0.75, edgecolors="white", linewidths=0.2)
        n_val = len(umap_coh)
        ax.set_title(f"{sub_title} (N={n_val})", fontsize=8.5, fontweight="bold", pad=6)
        ax.set_xlabel("UMAP-1", fontsize=8)
        if ci == 0:
            ax.set_ylabel("UMAP-2", fontsize=8)
        ax.tick_params(labelsize=7)
        if ci == 0:
            ax.legend(fontsize=7, loc="lower left", markerscale=0.6, framealpha=0.85)

    # Panel B title + letter
    fig.text(0.58, 0.96, "B", fontsize=13, fontweight="bold", va="top")
    fig.text(0.60, 0.96, "UMAP: vendor clustering (Original)",
             fontsize=10, fontweight="bold", va="top")

    # Panel C: full-width bar chart
    axC = fig.add_subplot(outer[1, 0])
    axC.set_title("Baseline vendor $\\eta^2$ per biomarker: HC vs ALL",
                  fontsize=10, fontweight="bold", pad=10)

    hc_eta = [comp[(comp["cohort"] == "HC") & (comp["method"] == "Original") &
                   (comp["biomarker"] == b)]["vendor_eta2"].values[0] for b in BIOMARKERS]
    all_eta = [comp[(comp["cohort"] == "ALL") & (comp["method"] == "Original") &
                    (comp["biomarker"] == b)]["vendor_eta2"].values[0] for b in BIOMARKERS]

    x = np.arange(len(BIOMARKERS))
    w = 0.35

    bars_hc = axC.bar(x - w/2, hc_eta, w, label="HC (N=188)", color=WONG["blue"],
                      edgecolor="black", linewidth=0.4, alpha=0.85)
    bars_all = axC.bar(x + w/2, all_eta, w, label="ALL (N=246)", color=WONG["orange"],
                       edgecolor="black", linewidth=0.4, alpha=0.85)

    perm_hc = perm[(perm["cohort"] == "HC") & (perm["method"] == "Original")]["R2"].values[0]
    perm_all = perm[(perm["cohort"] == "ALL") & (perm["method"] == "Original")]["R2"].values[0]

    axC.axhline(y=perm_hc, color=WONG["blue"], linestyle="--", linewidth=1.0, alpha=0.7)
    axC.axhline(y=perm_all, color=WONG["orange"], linestyle="--", linewidth=1.0, alpha=0.7)

    # Use legend for PERMANOVA reference lines instead of right-side text to avoid bar overlap
    from matplotlib.lines import Line2D
    perm_handles = [
        Line2D([0], [0], color=WONG["blue"], linestyle="--", linewidth=1.5, label=f"HC PERMANOVA $R^2$={perm_hc:.3f}"),
        Line2D([0], [0], color=WONG["orange"], linestyle="--", linewidth=1.5, label=f"ALL PERMANOVA $R^2$={perm_all:.3f}"),
    ]
    y_max = max(max(hc_eta), max(all_eta)) * 1.22
    axC.set_ylim(0, y_max * 1.12)

    # Value labels above bars
    for xi, vi in zip(x - w/2, hc_eta):
        axC.text(xi, vi + y_max * 0.015, f"{vi:.2f}", ha="center", fontsize=6.5, color="#333")
    for xi, vi in zip(x + w/2, all_eta):
        axC.text(xi, vi + y_max * 0.015, f"{vi:.2f}", ha="center", fontsize=6.5, color="#333")

    axC.set_xticks(x)
    axC.set_xticklabels([BM_LABEL[b] for b in BIOMARKERS], fontsize=9)
    axC.set_ylabel("Vendor $\\eta^2$ (one-way ANOVA)", fontsize=9.5)
    # Two-column legend: cohort bars + PERMANOVA lines
    handles_bars, labels_bars = axC.get_legend_handles_labels()
    axC.legend(handles=handles_bars + perm_handles, fontsize=7.5, loc="upper left",
               framealpha=0.9, ncol=1)
    for tick, bio in zip(axC.get_xticklabels(), BIOMARKERS):
        tick.set_color(FAMILY_COLOR[BM_FAMILY[bio]])

    # Panel C letter
    fig.text(0.01, 0.46, "C", fontsize=13, fontweight="bold", va="top",
             transform=fig.transFigure)

    for fmt in ["png", "pdf"]:
        fig.savefig(OUT_DIR / f"Fig2_baseline_effects.{fmt}", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("  Fig2 saved.")


# ============================================================
# Fig3: Harmonization performance
# ============================================================
def gen_fig3(comp, perm):
    """4-panel performance comparison.

    Panels A & B: ZOOMED log-scale bars (y-axis [-4.2, -2.0]) with value labels,
                   baseline shown as top annotation instead of dashed line.
    Panels C & D: Linear bars with value labels and star markers.
    """
    print("[Fig3] Harmonization performance")
    fig, axes = plt.subplots(2, 2, figsize=(13, 10))
    axes = axes.flatten()

    cohorts = ["HC", "ALL"]
    cohort_colors = {"HC": WONG["blue"], "ALL": WONG["orange"]}
    methods_plot = ["LME", "ComBat", "ComBat-joint", "RELIEF", "CovBat"]

    agg = comp.groupby(["cohort", "method"]).agg(
        eta2_mean=("vendor_eta2", "mean"),
        age_r2_mean=("age_R2", "mean"),
        sex_r2_mean=("sex_R2", "mean"),
    ).reset_index()

    x = np.arange(len(methods_plot))
    w = 0.32

    # ---- Panels A & B: zoomed log bars ----
    def add_log_bars(ax, values_getter, baseline_getter, ylabel, title, panel_letter):
        """Zoomed log-scale bar chart with value labels and baseline annotation."""
        ax.text(-0.08, 1.06, panel_letter, transform=ax.transAxes,
                fontsize=12, fontweight="bold", va="top")
        ax.set_title(title, fontsize=9.5, fontweight="bold", pad=8)

        # Collect all post-harmonization values
        all_log_vals = []
        for ci, coh in enumerate(cohorts):
            raw_vals = [values_getter(coh, m) for m in methods_plot]
            log_vals = [np.log10(max(v, 1e-7)) for v in raw_vals]
            all_log_vals.extend(log_vals)
            ax.bar(x + (ci - 0.5) * w, log_vals, w, label=coh, color=cohort_colors[coh],
                   edgecolor="black", linewidth=0.4, alpha=0.85)
            # Value labels: actual scientific notation
            for xi, vi_raw in zip(x + (ci - 0.5) * w, raw_vals):
                ax.text(xi, np.log10(max(vi_raw, 1e-7)) + 0.03, f"{vi_raw:.1e}",
                        ha="center", fontsize=6.5, color="#333", rotation=0)

        # Zoomed y-axis to amplify small differences
        ax.set_ylim(-4.3, -1.8)

        # Baseline annotation at top (arrows pointing down to the axis edge)
        for ci, coh in enumerate(cohorts):
            bl = baseline_getter(coh)
            bl_log = np.log10(max(bl, 1e-7))
            color = cohort_colors[coh]
            # Draw a small triangle + text at the top edge to indicate baseline
            bx = 4.4 if ci == 0 else 4.4
            ax.annotate(f"Baseline {coh}\n$\\eta^2$={bl:.3f}" if "\\eta" in ylabel
                        else f"Baseline {coh}\n$R^2$={bl:.3f}",
                        xy=(bx, -1.85), xytext=(bx, -1.82),
                        fontsize=7, color=color, ha="right", va="top", fontweight="bold",
                        bbox=dict(boxstyle="round,pad=0.2", facecolor="white",
                                  edgecolor=color, alpha=0.9, linewidth=0.8))

        ax.set_xticks(x)
        ax.set_xticklabels([METHOD_LABEL[m] for m in methods_plot], fontsize=8)
        ax.set_ylabel(ylabel, fontsize=9)
        ax.tick_params(labelsize=7.5)
        # Add horizontal grid lines for readability
        ax.yaxis.grid(True, linestyle=":", linewidth=0.4, alpha=0.6)
        ax.set_axisbelow(True)

    # ---- Panels C & D: linear bars ----
    def add_linear_bars(ax, values_getter, baseline_getter, ylabel, title, panel_letter,
                        fmt_str="{:.4f}"):
        ax.text(-0.08, 1.06, panel_letter, transform=ax.transAxes,
                fontsize=12, fontweight="bold", va="top")
        ax.set_title(title, fontsize=9.5, fontweight="bold", pad=8)

        for coh in cohorts:
            bl = baseline_getter(coh)
            ax.axhline(y=bl, color=cohort_colors[coh], linestyle="--",
                       linewidth=1.0, alpha=0.5)

        for ci, coh in enumerate(cohorts):
            vals = [values_getter(coh, m) for m in methods_plot]
            ax.bar(x + (ci - 0.5) * w, vals, w, label=coh, color=cohort_colors[coh],
                   edgecolor="black", linewidth=0.3, alpha=0.85)
            y_range = max(vals) - min(vals) if max(vals) > min(vals) else 0.01
            for xi, vi in zip(x + (ci - 0.5) * w, vals):
                ax.text(xi, vi + y_range * 0.04, fmt_str.format(vi),
                        ha="center", fontsize=6.5, color="#333")

        ax.set_xticks(x)
        ax.set_xticklabels([METHOD_LABEL[m] for m in methods_plot], fontsize=8)
        ax.set_ylabel(ylabel, fontsize=9)
        ax.tick_params(labelsize=7.5)
        ax.yaxis.grid(True, linestyle=":", linewidth=0.4, alpha=0.6)
        ax.set_axisbelow(True)

    # Panel A: zoomed log eta2
    add_log_bars(axes[0],
                 lambda c, m: agg[(agg["cohort"] == c) & (agg["method"] == m)]["eta2_mean"].values[0],
                 lambda c: agg[(agg["cohort"] == c) & (agg["method"] == "Original")]["eta2_mean"].values[0],
                 "log$_{10}$($\\eta^2$)", "Univariate vendor effect (zoomed)", "A")

    # Panel B: zoomed log PERMANOVA R2
    add_log_bars(axes[1],
                 lambda c, m: perm[(perm["cohort"] == c) & (perm["method"] == m)]["R2"].values[0],
                 lambda c: perm[(perm["cohort"] == c) & (perm["method"] == "Original")]["R2"].values[0],
                 "log$_{10}$(PERMANOVA $R^2$)", "Multivariate vendor effect (zoomed)", "B")

    # Panel C: Age R2 (linear)
    add_linear_bars(axes[2],
                    lambda c, m: agg[(agg["cohort"] == c) & (agg["method"] == m)]["age_r2_mean"].values[0],
                    lambda c: agg[(agg["cohort"] == c) & (agg["method"] == "Original")]["age_r2_mean"].values[0],
                    "Age $R^2$ (mean)", "Age semi-partial $R^2$", "C")

    # Star for best method in Panel C (CovBat HC)
    best_idx_c = methods_plot.index("CovBat")
    best_y_c = agg[(agg["cohort"] == "HC") & (agg["method"] == "CovBat")]["age_r2_mean"].values[0]
    star_y_c = best_y_c + 0.0045
    axes[2].scatter(best_idx_c - w/2, star_y_c, s=130, marker="*",
                    color=WONG["vermilion"], zorder=5, edgecolors="black", linewidths=0.4,
                    clip_on=False)
    axes[2].annotate("best", xy=(best_idx_c - w/2, star_y_c),
                     xytext=(best_idx_c - w/2, star_y_c + 0.0035),
                     fontsize=7, color=WONG["vermilion"], ha="center", fontweight="bold")

    # Panel D: Sex R2 (linear)
    add_linear_bars(axes[3],
                    lambda c, m: agg[(agg["cohort"] == c) & (agg["method"] == m)]["sex_r2_mean"].values[0],
                    lambda c: agg[(agg["cohort"] == c) & (agg["method"] == "Original")]["sex_r2_mean"].values[0],
                    "Sex $R^2$ (mean)", "Sex semi-partial $R^2$", "D")

    # Star for best method in Panel D (RELIEF HC)
    best_idx_d = methods_plot.index("RELIEF")
    best_y_d = agg[(agg["cohort"] == "HC") & (agg["method"] == "RELIEF")]["sex_r2_mean"].values[0]
    star_y_d = best_y_d + 0.0025
    axes[3].scatter(best_idx_d - w/2, star_y_d, s=130, marker="*",
                    color=WONG["vermilion"], zorder=5, edgecolors="black", linewidths=0.4,
                    clip_on=False)
    axes[3].annotate("best", xy=(best_idx_d - w/2, star_y_d),
                     xytext=(best_idx_d - w/2, star_y_d + 0.0025),
                     fontsize=7, color=WONG["vermilion"], ha="center", fontweight="bold")

    # Single figure-level legend
    handles = [Patch(facecolor=cohort_colors[c], edgecolor="black",
                     label=c + " (N=188)" if c == "HC" else c + " (N=246)") for c in cohorts]
    handles += [Line2D([0], [0], color=cohort_colors[c], linestyle="--",
                       linewidth=1.5, label=f"Baseline {c}") for c in cohorts]
    handles += [Line2D([0], [0], marker="*", color="w", markerfacecolor=WONG["vermilion"],
                       markersize=10, label="Best method", markeredgecolor="black", markeredgewidth=0.4)]
    fig.legend(handles=handles, loc="upper center", ncol=5, fontsize=8.5,
               frameon=False, bbox_to_anchor=(0.5, 1.01))

    fig.suptitle("Harmonization performance across evaluation endpoints",
                 fontsize=12, fontweight="bold", y=1.05)
    fig.subplots_adjust(hspace=0.45, wspace=0.30, top=0.90, bottom=0.06, left=0.08, right=0.97)
    for fmt in ["png", "pdf"]:
        fig.savefig(OUT_DIR / f"Fig3_performance.{fmt}", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("  Fig3 saved.")


# ============================================================
# Fig4: Trade-off scatter
# ============================================================
def gen_fig4(comp, perm):
    """Pareto frontier: x = -log10(PERMANOVA R2), y = r_pre,post."""
    print("[Fig4] Trade-off scatter (Pareto frontier)")
    fig, axes = plt.subplots(1, 2, figsize=(12, 5.5))

    cohorts = ["HC", "ALL"]
    methods_plot = ["LME", "ComBat", "ComBat-joint", "RELIEF", "CovBat"]

    agg = comp.groupby(["cohort", "method"]).agg(
        r_mean=("r_pre_post", "mean"),
        age_r2_mean=("age_R2", "mean"),
    ).reset_index()

    for idx, cohort in enumerate(cohorts):
        ax = axes[idx]
        ax.text(-0.08, 1.06, chr(65 + idx), transform=ax.transAxes,
                fontsize=12, fontweight="bold", va="top")
        ax.set_title(f"{cohort} cohort (N={'188' if cohort == 'HC' else '246'})",
                     fontsize=10, fontweight="bold", pad=8)

        bl_age = agg[(agg["cohort"] == cohort) & (agg["method"] == "Original")]["age_r2_mean"].values[0]

        xs, ys, sizes, colors_list = [], [], [], []
        for m in methods_plot:
            r2 = perm[(perm["cohort"] == cohort) & (perm["method"] == m)]["R2"].values[0]
            rpp = agg[(agg["cohort"] == cohort) & (agg["method"] == m)]["r_mean"].values[0]
            age_r2 = agg[(agg["cohort"] == cohort) & (agg["method"] == m)]["age_r2_mean"].values[0]
            age_gain = max(age_r2 - bl_age, 0.0001)
            xs.append(-np.log10(max(r2, 1e-6)))
            ys.append(rpp)
            sizes.append(120 + age_gain * 6000)
            colors_list.append(METHOD_COLOR[m])

        ax.scatter(xs, ys, s=sizes, c=colors_list, edgecolors="black",
                   linewidths=0.6, zorder=3, alpha=0.85)

        # Pareto frontier line
        order = sorted(range(len(xs)), key=lambda i: xs[i])
        ax.plot([xs[i] for i in order], [ys[i] for i in order],
                linestyle="--", color="#aaa", linewidth=1.0, zorder=2)

        # Labels: position to avoid overlap, use leader lines
        label_offsets = {
            "LME":          (-0.30,  0.008, "right"),
            "ComBat":       ( 0.12,  0.010, "left"),
            "ComBat-joint": ( 0.15, -0.015, "left"),
            "RELIEF":       ( 0.15,  0.012, "left"),
            "CovBat":       (-0.15, -0.015, "right"),
        }
        for xi, yi, m in zip(xs, ys, methods_plot):
            dx, dy, ha = label_offsets[m]
            ax.annotate(METHOD_LABEL[m], xy=(xi, yi), xytext=(xi + dx, yi + dy),
                        fontsize=8.5, fontweight="bold", color=METHOD_COLOR[m],
                        ha=ha,
                        arrowprops=dict(arrowstyle="-", color=METHOD_COLOR[m], lw=0.5,
                                        connectionstyle="arc3,rad=0.1"))

        ax.set_xlabel("$-$log$_{10}$(PERMANOVA $R^2$)", fontsize=9)
        ax.set_ylabel("r$_{pre,post}$ (subject preservation)", fontsize=9)
        ax.tick_params(labelsize=7.5)
        y_pad = (max(ys) - min(ys)) * 0.25
        ax.set_ylim(min(ys) - y_pad, max(ys) + y_pad)
        x_pad = (max(xs) - min(xs)) * 0.15
        ax.set_xlim(min(xs) - x_pad, max(xs) + x_pad * 1.3)
        ax.text(0.97, 0.03, "Pareto frontier", transform=ax.transAxes,
                fontsize=8, color="#999", ha="right", style="italic")

    handles = [Line2D([0], [0], marker="o", color="w", markerfacecolor=METHOD_COLOR[m],
                      markersize=9, label=METHOD_LABEL[m], markeredgecolor="black",
                      markeredgewidth=0.4) for m in methods_plot]
    fig.legend(handles=handles, loc="lower center", ncol=5, fontsize=9,
               frameon=False, bbox_to_anchor=(0.5, -0.04))
    fig.suptitle("Central trade-off: multivariate decontamination vs. subject-level preservation",
                 fontsize=11, fontweight="bold", y=1.02)
    fig.subplots_adjust(wspace=0.30, top=0.88, bottom=0.14, left=0.08, right=0.97)
    for fmt in ["png", "pdf"]:
        fig.savefig(OUT_DIR / f"Fig4_tradeoff.{fmt}", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("  Fig4 saved.")


# ============================================================
# Fig5: r_pre,post heatmap with marginals
# ============================================================
def gen_fig5(comp):
    """Heatmap of r_pre,post (5 methods x 8 biomarkers) for HC and ALL."""
    print("[Fig5] r_pre,post heatmap with marginals")
    methods_plot = ["LME", "ComBat", "ComBat-joint", "RELIEF", "CovBat"]

    hc_data = comp[(comp["cohort"] == "HC") & (comp["method"].isin(methods_plot))].pivot_table(
        index="method", columns="biomarker", values="r_pre_post")
    hc_data = hc_data.reindex(index=methods_plot, columns=BIOMARKERS)

    all_data = comp[(comp["cohort"] == "ALL") & (comp["method"].isin(methods_plot))].pivot_table(
        index="method", columns="biomarker", values="r_pre_post")
    all_data = all_data.reindex(index=methods_plot, columns=BIOMARKERS)

    eta2_hc = comp[(comp["cohort"] == "HC") & (comp["method"] == "Original")].set_index("biomarker")["vendor_eta2"]
    eta2_hc = eta2_hc.reindex(BIOMARKERS)

    cmap_r = LinearSegmentedColormap.from_list("r_cmap", ["#FFFFFF", "#A8D5A8", "#009E73"])

    fig = plt.figure(figsize=(16, 8.5))
    # 3 rows x 6 cols: [heatmap | gap | right_marg] x 2 cohorts
    gs = GridSpec(3, 6, figure=fig,
                  width_ratios=[5, 0.35, 3.2, 5, 0.35, 3.2],
                  height_ratios=[0.6, 3.2, 1.0],
                  hspace=0.28, wspace=0.12)

    for ci, (cohort, data) in enumerate([("HC", hc_data), ("ALL", all_data)]):
        col_offset = ci * 3

        # Top marginal: per-biomarker mean r
        ax_top = fig.add_subplot(gs[0, col_offset])
        bm_means = data.mean(axis=0)
        bars = ax_top.bar(np.arange(len(BIOMARKERS)), bm_means,
                          edgecolor="black", linewidth=0.3, alpha=0.75, width=0.7)
        for bar, bio in zip(bars, BIOMARKERS):
            bar.set_color(FAMILY_COLOR[BM_FAMILY[bio]])
        ax_top.set_ylim(0.7, 1.05)
        ax_top.set_xticks([])
        ax_top.set_yticks([0.8, 1.0])
        ax_top.tick_params(labelsize=7)
        ax_top.set_ylabel("Mean r", fontsize=8)
        ax_top.set_title(f"{cohort} cohort (N={'188' if cohort == 'HC' else '246'})",
                         fontsize=10, fontweight="bold", pad=8)
        ax_top.spines["bottom"].set_visible(False)

        # Heatmap
        ax_hm = fig.add_subplot(gs[1, col_offset])
        im = ax_hm.imshow(data.values, aspect="auto", cmap=cmap_r,
                          vmin=0.4, vmax=1.0, interpolation="nearest")
        for i in range(data.shape[0]):
            for j in range(data.shape[1]):
                val = data.values[i, j]
                color = "white" if val < 0.72 else "black"
                ax_hm.text(j, i, f"{val:.2f}", ha="center", va="center",
                           fontsize=7.5, color=color, fontweight="bold")
        ax_hm.set_xticks(np.arange(len(BIOMARKERS)))
        ax_hm.set_xticklabels([BM_LABEL[b] for b in BIOMARKERS], fontsize=8, rotation=0)
        ax_hm.set_yticks(np.arange(len(methods_plot)))
        ax_hm.set_yticklabels([METHOD_LABEL[m] for m in methods_plot], fontsize=8.5)
        for tick, bio in zip(ax_hm.get_xticklabels(), BIOMARKERS):
            tick.set_color(FAMILY_COLOR[BM_FAMILY[bio]])
        ax_hm.tick_params(axis="x", pad=6)

        # Right marginal: per-method mean r with rank (horizontal bars)
        ax_right = fig.add_subplot(gs[1, col_offset + 2])
        method_means = data.mean(axis=1)
        ranks = method_means.rank(ascending=False, method="min").astype(int)
        y_pos = np.arange(len(methods_plot))
        ax_right.barh(y_pos, method_means - 0.7, left=0.7,
                      color=[METHOD_COLOR[m] for m in methods_plot],
                      edgecolor="black", linewidth=0.3, alpha=0.85, height=0.55)
        # Rank + value labels to the right of bars
        for yi, (val, rank, m) in enumerate(zip(method_means, ranks, methods_plot)):
            ax_right.text(1.02, yi, f"#{rank}  r={val:.3f}",
                          va="center", ha="left", fontsize=8, color="#333")
        ax_right.set_xlim(0.65, 1.38)
        ax_right.set_yticks([])
        ax_right.set_xlabel("Mean r", fontsize=8)
        ax_right.tick_params(labelsize=7)
        ax_right.invert_yaxis()
        ax_right.set_title("Rank", fontsize=9, fontweight="bold", pad=8)
        ax_right.spines["left"].set_visible(False)

        # Bottom strip: baseline vendor eta2
        ax_bot = fig.add_subplot(gs[2, col_offset])
        cmap_eta = plt.cm.magma
        for j, (bio, e2) in enumerate(zip(BIOMARKERS, eta2_hc.values)):
            ax_bot.add_patch(plt.Rectangle((j - 0.4, 0), 0.8, 1,
                                           facecolor=cmap_eta(e2 / eta2_hc.max()),
                                           edgecolor="black", linewidth=0.3))
            text_color = "white" if e2 / eta2_hc.max() > 0.35 else "black"
            ax_bot.text(j, 0.5, f"{e2:.2f}", ha="center", va="center",
                        fontsize=8, color=text_color, fontweight="bold")
        ax_bot.set_xlim(-0.5, len(BIOMARKERS) - 0.5)
        ax_bot.set_ylim(0, 1)
        ax_bot.set_xticks([])
        ax_bot.set_yticks([])
        ax_bot.set_ylabel("Baseline\nvendor $\\eta^2$", fontsize=8)
        ax_bot.spines["top"].set_visible(False)
        ax_bot.spines["bottom"].set_visible(False)

    # Colorbar — placed far right to avoid both right marginals
    ax_cb = fig.add_axes([0.965, 0.30, 0.012, 0.35])
    cb = fig.colorbar(im, cax=ax_cb, orientation="vertical")
    cb.set_label("r$_{pre,post}$", fontsize=9)
    cb.ax.tick_params(labelsize=7.5)

    fig.suptitle("Subject-level preservation (r$_{pre,post}$) by biomarker and method",
                 fontsize=12, fontweight="bold", y=0.98)
    for fmt in ["png", "pdf"]:
        fig.savefig(OUT_DIR / f"Fig5_rpre_post_marginals.{fmt}", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("  Fig5 saved.")


# ============================================================
# FigS: Design C robustness
# ============================================================
def gen_figs():
    """Design C balanced subsample robustness analysis."""
    print("[FigS] Design C robustness")
    if not DESIGN_C.exists():
        print(f"  WARNING: {DESIGN_C} not found")
        return

    df = pd.read_csv(DESIGN_C)
    methods_here = [m for m in ["Original", "LME", "ComBat", "ComBat-joint", "RELIEF", "CovBat"]
                    if m in df["method"].values]
    methods_no_orig = [m for m in methods_here if m != "Original"]

    fig, axes = plt.subplots(1, 2, figsize=(11, 5))

    # Panel A: r_pre_post
    ax = axes[0]
    ax.text(-0.08, 1.06, "A", transform=ax.transAxes, fontsize=12, fontweight="bold", va="top")
    ax.set_title("Subject preservation (r$_{pre,post}$)\nunder balanced subsampling",
                 fontsize=9.5, fontweight="bold", pad=8)
    data_r = [df[df["method"] == m]["r_pre_post"].dropna().values for m in methods_no_orig]
    bp = ax.boxplot(data_r, tick_labels=[METHOD_LABEL[m] for m in methods_no_orig],
                    patch_artist=True, showfliers=True, widths=0.5,
                    flierprops={"markersize": 3, "alpha": 0.3, "marker": "o"},
                    medianprops={"color": "black", "linewidth": 1.0},
                    whiskerprops={"linewidth": 0.6},
                    capprops={"linewidth": 0.6})
    for patch, m in zip(bp["boxes"], methods_no_orig):
        patch.set_facecolor(METHOD_COLOR[m])
        patch.set_alpha(0.55)
    for i, vals in enumerate(data_r):
        jitter = np.random.normal(0, 0.04, len(vals))
        ax.scatter(np.full(len(vals), i + 1) + jitter, vals, s=8, alpha=0.25, color="black", zorder=3)

    # Mean labels: placed above each box with alternating vertical offsets
    label_offsets_y = [0.012, -0.015, 0.018, -0.010, 0.008]
    for i, vals in enumerate(data_r):
        mean_v = np.mean(vals)
        off_idx = i % len(label_offsets_y)
        label_y = mean_v + label_offsets_y[off_idx]
        label_y = max(0.80, min(label_y, 1.04))
        ax.text(i + 1, label_y, f"$\\mu$={mean_v:.3f}", fontsize=7, color="#333",
                va="center", ha="center",
                bbox=dict(boxstyle="round,pad=0.15", facecolor="white", alpha=0.8, edgecolor="none"))
    ax.set_ylabel("r$_{pre,post}$", fontsize=9.5)
    ax.set_ylim(0.76, 1.06)
    ax.tick_params(labelsize=8, axis="x", rotation=15)

    # Panel B: eta2
    ax = axes[1]
    ax.text(-0.08, 1.06, "B", transform=ax.transAxes, fontsize=12, fontweight="bold", va="top")
    ax.set_title("Vendor effect ($\\eta^2$) under\nbalanced subsampling",
                 fontsize=9.5, fontweight="bold", pad=8)
    data_e = [df[df["method"] == m]["eta2"].dropna().values for m in methods_no_orig]
    bp2 = ax.boxplot(data_e, tick_labels=[METHOD_LABEL[m] for m in methods_no_orig],
                     patch_artist=True, showfliers=True, widths=0.5,
                     flierprops={"markersize": 3, "alpha": 0.3, "marker": "o"},
                     medianprops={"color": "black", "linewidth": 1.0},
                     whiskerprops={"linewidth": 0.6},
                     capprops={"linewidth": 0.6})
    for patch, m in zip(bp2["boxes"], methods_no_orig):
        patch.set_facecolor(METHOD_COLOR[m])
        patch.set_alpha(0.55)
    for i, vals in enumerate(data_e):
        jitter = np.random.normal(0, 0.04, len(vals))
        ax.scatter(np.full(len(vals), i + 1) + jitter, vals, s=8, alpha=0.25, color="black", zorder=3)
    ax.set_ylabel("Vendor $\\eta^2$", fontsize=9.5)
    ax.set_yscale("log")
    ax.tick_params(labelsize=8, axis="x", rotation=15)
    orig_eta = df[df["method"] == "Original"]["eta2"].mean()
    ax.axhline(y=orig_eta, color=WONG["grey"], linestyle="--", linewidth=1.0)
    ax.text(0.02, orig_eta, f" baseline $\\eta^2$={orig_eta:.3f}",
            fontsize=7.5, color="#555",
            transform=ax.get_yaxis_transform(), va="center", ha="left")

    fig.suptitle("Design C: Robustness to balanced vendor composition (20 random seeds, n$\\approx$98/seed)",
                 fontsize=11, fontweight="bold", y=1.03)
    fig.subplots_adjust(wspace=0.30, top=0.82, bottom=0.12, left=0.08, right=0.97)
    for fmt in ["png", "pdf"]:
        fig.savefig(OUT_DIR / f"FigS_design_C_robustness.{fmt}", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("  FigS saved.")


# ============================================================
# Main
# ============================================================
if __name__ == "__main__":
    print("=" * 60)
    print("Generating ALL statistical figures (Fig2-Fig5, FigS)")
    print("with REAL data from pre-computed analysis results")
    print("=" * 60)

    master, comp, perm, umap = load_all()
    gen_fig2(master, comp, perm, umap)
    gen_fig3(comp, perm)
    gen_fig4(comp, perm)
    gen_fig5(comp)
    gen_figs()

    print("\n" + "=" * 60)
    print(f"DONE! All figures saved to: {OUT_DIR}")
    for f in sorted(OUT_DIR.glob("Fig*.png")):
        size_kb = f.stat().st_size / 1024
        print(f"  {f.name:45s} {size_kb:6.0f} KB")
