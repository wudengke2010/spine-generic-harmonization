"""
Supplementary Fig S4 — Cross-metric covariance structure preservation
======================================================================
Following the covariance-harmonisation emphasis of Chen et al. (HBM 2022):

Panel (a): 8x8 Pearson correlation matrices of the biomarkers (HC cohort,
           complete cases N=188) before and after each method
           (Original, ComBat, ComBat-joint, CovBat, RELIEF, LME/FE).
Panel (b): 6x6 Frobenius distance matrices between the methods' correlation
           matrices, for HC (N=188) and ALL (N=246) cohorts.

Frobenius distance on the upper-triangle-excluding-diagonal correlation
entries: d(A,B) = ||C_A - C_B||_F, lower = more similar covariance structure.

Outputs:
  results/step4_eval/corr_matrices_HC.csv / corr_matrices_ALL.csv (long)
  results/step4_eval/frobenius_distance_HC.csv / frobenius_distance_ALL.csv
  submission/supplementary/FigS4_covariance.pdf / .png
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
    "font.size": 10, "axes.titlesize": 11, "axes.labelsize": 10,
    "xtick.labelsize": 9, "ytick.labelsize": 9,
    "axes.linewidth": 0.8, "xtick.major.width": 0.8, "ytick.major.width": 0.8,
    "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
})

RES = Path("E:/boshi/qm_harmonization_paper/results")
OUT_CSV = RES / "step4_eval"
OUT_FIG = Path("C:/Users/admin/WorkBuddy/2026-07-06-12-34-04/submission/supplementary")

BIOMARKERS = ["T2w_CSA", "GM_CSA_mm2", "MTR", "MTsat",
              "FA", "MD", "AD", "RD"]
BM_LABEL = {
    "T2w_CSA": "Cord CSA", "GM_CSA_mm2": "GM CSA",
    "MTR": "MTR", "MTsat": "MTsat",
    "FA": "FA", "MD": "MD", "AD": "AD", "RD": "RD",
}
# short horizontal labels for panel (b) frobenius heatmaps to avoid vertical crowding
METHOD_SHORT = {
    "Original": "Orig",
    "ComBat": "ComBat",
    "ComBat-joint": "ComBat-j",
    "CovBat": "CovBat",
    "RELIEF": "RELIEF",
    "LME/FE": "LME/FE",
}
METHODS = [  # (display, csv, suffix)
    ("Original",    RES / "step3_combat/biomarkers_combat_HC.csv", ""),
    ("ComBat",      RES / "step3_combat/biomarkers_combat_HC.csv", "_combat"),
    ("ComBat-joint", RES / "step3_combat/biomarkers_combat_joint_HC.csv", "_combatJ"),
    ("CovBat",      RES / "step3_covbat/biomarkers_covbat_HC.csv", "_covbat"),
    ("RELIEF",      RES / "step3_relief/biomarkers_relief_HC.csv", "_relief"),
    ("LME/FE",      RES / "step3_lme/biomarkers_lme_HC.csv", "_lme"),
]
# accent colour per method (muted, print-safe) used for titles and highlights
METHOD_COLORS = {
    "Original":     "#404040",
    "ComBat":       "#B2182B",
    "ComBat-joint": "#D6604D",
    "CovBat":       "#2166AC",
    "RELIEF":       "#4393C3",
    "LME/FE":       "#1B7837",
}
# key method pairs to highlight in the Frobenius panel
HIGHLIGHT_PAIRS = [("CovBat", "RELIEF"), ("LME/FE", "Original")]


def corr_matrix(cohort_label: str) -> tuple[dict[str, np.ndarray], int]:
    """Compute 8x8 correlation matrices per method for one cohort."""
    mats, n_used = {}, None
    for name, csv, suf in METHODS:
        if cohort_label == "ALL":
            csv = Path(str(csv).replace("_HC.csv", "_ALL.csv"))
        df = pd.read_csv(csv)
        cols = [m + suf for m in BIOMARKERS]
        sub = df[cols].dropna()
        C = np.corrcoef(sub.values, rowvar=False)
        mats[name] = C
        if name == "Original":
            n_used = len(sub)
    return mats, n_used


def long_corr(mats: dict[str, np.ndarray]) -> pd.DataFrame:
    rows = []
    for method, C in mats.items():
        for i, a in enumerate(BIOMARKERS):
            for j, b in enumerate(BIOMARKERS):
                rows.append({"method": method, "biomarker1": a,
                             "biomarker2": b, "r": C[i, j]})
    return pd.DataFrame(rows)


def frobenius(mats: dict[str, np.ndarray]) -> pd.DataFrame:
    names = list(mats)
    D = pd.DataFrame(0.0, index=names, columns=names)
    for a in names:
        for b in names:
            D.loc[a, b] = np.linalg.norm(mats[a] - mats[b], ord="fro")
    return D


def draw_correlation_grid(ax_arr, mats: dict[str, np.ndarray]) -> None:
    names = list(mats)
    n = len(BIOMARKERS)
    for k, name in enumerate(names):
        ax = ax_arr[k]
        C = mats[name]
        im = ax.imshow(C, vmin=-1, vmax=1, cmap="RdBu_r",
                       interpolation="nearest")
        ax.set_title(name, fontweight="bold", fontsize=11, pad=8,
                     color=METHOD_COLORS[name])
        ax.set_xticks(range(n))
        ax.set_xticklabels([BM_LABEL[b] for b in BIOMARKERS], rotation=90,
                           ha="center", va="top", fontsize=8)
        ax.set_yticks(range(n))
        ax.set_yticklabels([BM_LABEL[b] for b in BIOMARKERS], fontsize=8)
        ax.tick_params(length=0)
        # white gridlines between cells
        ax.set_xticks(np.arange(-0.5, n), minor=True)
        ax.set_yticks(np.arange(-0.5, n), minor=True)
        ax.grid(which="minor", color="white", linewidth=0.7)
        ax.tick_params(which="minor", length=0)
        # Omit cell numbers: the colour map conveys the pattern; exact
        # values are available in the supplementary CSV files.
        # keep tick labels only on outer edges to reduce clutter
        if k % 3 != 0:
            ax.set_yticklabels([])
        if k < 3:
            ax.set_xticklabels([])
    return im


def draw_frobenius(ax, D: pd.DataFrame, cohort: str, vmax: float):
    names = list(D.index)
    M = D.values
    im = ax.imshow(M, vmin=0, vmax=vmax, cmap="viridis",
                   interpolation="nearest")
    n = len(names)
    ax.set_xticks(range(n))
    ax.set_xticklabels([METHOD_SHORT[n] for n in names], rotation=0,
                       ha="center", fontsize=8)
    ax.set_yticks(range(n))
    ax.set_yticklabels([METHOD_SHORT[n] for n in names], fontsize=8)
    for tick, name in zip(ax.get_yticklabels(), names):
        tick.set_color(METHOD_COLORS[name])
    for tick, name in zip(ax.get_xticklabels(), names):
        tick.set_color(METHOD_COLORS[name])
    ax.set_title(f"{cohort}", fontweight="bold", loc="left", fontsize=11,
                 pad=8)
    ax.tick_params(length=0)
    # Move ticks to top/right so labels do not overlap first row/column cells
    ax.xaxis.tick_top()
    ax.yaxis.tick_right()
    for i in range(n):
        for j in range(n):
            v = M[i, j]
            color = "white" if v > vmax * 0.55 else "black"
            txt = "\u2013" if i == j else f"{v:.2f}"
            ax.text(j, i, txt, ha="center", va="center",
                    fontsize=9.5, color=color)
    ax.set_xticks(np.arange(-0.5, n), minor=True)
    ax.set_yticks(np.arange(-0.5, n), minor=True)
    ax.grid(which="minor", color="white", linewidth=0.8)
    ax.tick_params(which="minor", length=0)
    # highlight key method pairs with a crisp outline
    for a, b in HIGHLIGHT_PAIRS:
        if a in names and b in names:
            ia, ib = names.index(a), names.index(b)
            for (r_, c_) in [(ia, ib), (ib, ia)]:
                ax.add_patch(plt.Rectangle((c_ - 0.5, r_ - 0.5), 1, 1,
                             fill=False, edgecolor="white", linewidth=2.2,
                             zorder=3))
    return im


def main() -> None:
    OUT_CSV.mkdir(parents=True, exist_ok=True)
    OUT_FIG.mkdir(parents=True, exist_ok=True)

    mats_HC, n_hc = corr_matrix("HC")
    mats_ALL, n_all = corr_matrix("ALL")
    D_HC = frobenius(mats_HC)
    D_ALL = frobenius(mats_ALL)

    long_corr(mats_HC).to_csv(OUT_CSV / "corr_matrices_HC.csv", index=False)
    long_corr(mats_ALL).to_csv(OUT_CSV / "corr_matrices_ALL.csv", index=False)
    D_HC.to_csv(OUT_CSV / "frobenius_distance_HC.csv")
    D_ALL.to_csv(OUT_CSV / "frobenius_distance_ALL.csv")

    print(f"HC  N={n_hc}   ALL N={n_all}")
    print("\nFrobenius distance from Original (HC):")
    print(D_HC["Original"].round(3).to_string())
    print("\nFrobenius distance from Original (ALL):")
    print(D_ALL["Original"].round(3).to_string())

    # ---- figure ---------------------------------------------------------
    # 6-column mosaic: panel (a) 2 rows x 3 cols (each matrix spans 2 cols);
    # panel (b) each heatmap spans 3 cols for relaxed labels.
    mosaic = [
        ["hdr_a", "hdr_a", "hdr_a", "hdr_a", "hdr_a", "hdr_a"],
        ["a0", "a0", "a1", "a1", "a2", "a2"],
        ["a3", "a3", "a4", "a4", "a5", "a5"],
        ["cb_a", "cb_a", "cb_a", "cb_a", "cb_a", "cb_a"],
        ["hdr_b", "hdr_b", "hdr_b", "hdr_b", "hdr_b", "hdr_b"],
        ["b_hc", "b_hc", "b_hc", "b_all", "b_all", "b_all"],
        ["cb_b", "cb_b", "cb_b", "cb_b", "cb_b", "cb_b"],
    ]
    axes = plt.figure(figsize=(14.0, 13.5), constrained_layout=True
                      ).subplot_mosaic(
        mosaic,
        height_ratios=[0.32, 1.4, 1.4, 0.20, 0.32, 1.5, 0.20],
        empty_sentinel=None,
    )
    fig = axes["a0"].figure
    fig.get_layout_engine().set(w_pad=0.14, h_pad=0.06,
                                wspace=0.05, hspace=0.05)

    # panel headers as text-only mosaic cells (auto-centred by layout engine)
    for key, label in [
        ("hdr_a", "(a) Biomarker correlation matrices "
                  f"(HC cohort, N={n_hc})"),
        ("hdr_b", "(b) Frobenius distance between correlation matrices"),
    ]:
        axh = axes[key]
        axh.axis("off")
        axh.text(0.5, 0.5, label, transform=axh.transAxes,
                 ha="center", va="center", fontsize=12, fontweight="bold")

    # (a) correlation grids
    ax_arr = [axes[k] for k in ["a0", "a1", "a2", "a3", "a4", "a5"]]
    im_a = draw_correlation_grid(ax_arr, mats_HC)

    cbar_a = fig.colorbar(im_a, cax=axes["cb_a"], orientation="horizontal",
                          aspect=38, shrink=0.82)
    cbar_a.outline.set_linewidth(0.8)
    cbar_a.outline.set_edgecolor("0.3")
    axes["cb_a"].set_xlabel("Pearson r", fontsize=9)
    axes["cb_a"].tick_params(labelsize=8)

    # (b) Frobenius heatmaps
    vmax = float(max(D_HC.values.max(), D_ALL.values.max()))
    ax_hc, ax_all = axes["b_hc"], axes["b_all"]
    im_b = draw_frobenius(ax_hc, D_HC, f"HC (N={n_hc})", vmax)
    draw_frobenius(ax_all, D_ALL, f"ALL (N={n_all})", vmax)

    cbar_b = fig.colorbar(im_b, cax=axes["cb_b"], orientation="horizontal",
                          aspect=38, shrink=0.82)
    cbar_b.outline.set_linewidth(0.8)
    cbar_b.outline.set_edgecolor("0.3")
    axes["cb_b"].set_xlabel(
        "||C$_{method A}$ − C$_{method B}$||$_F$", fontsize=9)
    axes["cb_b"].tick_params(labelsize=8)

    fig.savefig(OUT_FIG / "FigS4_covariance.pdf", bbox_inches="tight")
    fig.savefig(OUT_FIG / "FigS4_covariance.png", dpi=300, bbox_inches="tight")
    fig.savefig(OUT_FIG / "FigS4_covariance.svg", bbox_inches="tight")
    plt.close(fig)
    print(f"\n[done] -> {OUT_FIG/'FigS4_covariance.pdf'}")


if __name__ == "__main__":
    main()
