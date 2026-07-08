"""
Step 5.8 — Figure 1a + 1b: Study Design Schematic (split into two landscape panels)
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Rectangle

mpl.rcParams.update({
    "font.family":      "sans-serif",
    "font.sans-serif":  ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size":        9,
    "pdf.fonttype":     42,
    "ps.fonttype":      42,
    "svg.fonttype":     "none",
})

WONG = {
    "blue":      "#0072B2",
    "sky":       "#56B4E9",
    "green":     "#009E73",
    "orange":    "#E69F00",
    "vermilion": "#D55E00",
    "pink":      "#CC79A7",
    "yellow":    "#F0E442",
    "grey":      "#999999",
}
VENDOR_COLOR = {"GE":      WONG["blue"],
                "Philips": WONG["green"],
                "Siemens": WONG["vermilion"]}
METHOD_COLOR = {"LME":          WONG["pink"],
                "ComBat":       WONG["blue"],
                "ComBat-joint": WONG["sky"],
                "RELIEF":       WONG["orange"],
                "CovBat":       WONG["green"]}
ACQ_COLOR    = {"T2w": "#5b5b5b", "T2*": "#5b5b5b",
                "MT":  "#5b5b5b", "DTI": "#5b5b5b"}


def panel_frame(ax, x, y, w, h, letter, title, *, subtitle=None):
    bg = FancyBboxPatch((x, y), w, h,
                        boxstyle="round,pad=0.4,rounding_size=1.0",
                        linewidth=0.8, edgecolor="#bbbbbb",
                        facecolor="#f7f7f7", zorder=1)
    ax.add_patch(bg)
    ax.text(x + 1.4, y + h - 1.6, letter,
            ha="left", va="top",
            fontsize=13, fontweight="bold", color="#222", zorder=5)
    ax.text(x + 5.0, y + h - 2.1, title,
            ha="left", va="top",
            fontsize=10.5, fontweight="bold", color="#222", zorder=5)
    if subtitle:
        ax.text(x + 5.0, y + h - 4.0, subtitle,
                ha="left", va="top",
                fontsize=8.5, style="italic", color="#555", zorder=5)


def chip(ax, x, y, w, h, label, *, face="white", edge="#444",
         text_color="black", fontsize=8.5, fontweight="normal",
         radius=0.5, zorder=3):
    box = FancyBboxPatch((x, y), w, h,
                         boxstyle=f"round,pad=0.12,rounding_size={radius}",
                         linewidth=0.7, edgecolor=edge,
                         facecolor=face, zorder=zorder)
    ax.add_patch(box)
    ax.text(x + w / 2, y + h / 2, label,
            ha="center", va="center",
            color=text_color, fontsize=fontsize,
            fontweight=fontweight, zorder=zorder + 1,
            linespacing=1.25)


def vertical_arrow(ax, x_center, y_top, y_bot, *,
                   color="#222", lw=1.6, mutation=18):
    a = FancyArrowPatch((x_center, y_top), (x_center, y_bot),
                        arrowstyle="-|>", mutation_scale=mutation,
                        color=color, linewidth=lw, zorder=2,
                        shrinkA=2, shrinkB=2)
    ax.add_patch(a)


def draw_panel_A(ax, x0, y0, w, h):
    panel_frame(ax, x0, y0, w, h, "A", "Dataset & cohorts",
                subtitle="spine-generic-multi-subject  ·  open BIDS")

    y_top_up = y0 + h - 6.0
    y_mid    = y0 + h * 0.46
    up_h     = y_top_up - y_mid

    left_x = x0 + 3.0
    ax.text(left_x, y_top_up - 1.5, "267 subjects",
            ha="left", va="top",
            fontsize=14, fontweight="bold", color="#222", zorder=5)
    ax.text(left_x, y_top_up - 5.6, "42 sites  ·  13 scanner models",
            ha="left", va="top", fontsize=9, color="#444", zorder=5)

    chip_w = (w - 6.0) * 0.26
    chip_h = up_h - 1.2
    gap    = 1.2
    chip_y = y_mid + 0.6
    chip_x2 = x0 + w - 3.0 - chip_w
    chip_x1 = chip_x2 - gap - chip_w
    chip(ax, chip_x1, chip_y, chip_w, chip_h,
         "HC cohort\nN = 188",
         face=WONG["yellow"], edge="#444",
         fontsize=9, fontweight="bold")
    chip(ax, chip_x2, chip_y, chip_w, chip_h,
         "ALL cohort\nN = 246\n(HC + MildComp + DCM)",
         face="#ffffff", edge="#444", fontsize=8.0)

    bar_x0 = x0 + 3.0
    bar_x1 = x0 + w - 3.0
    bar_h  = 3.6
    bar_y  = y0 + 2.6
    ax.text((bar_x0 + bar_x1) / 2, bar_y + bar_h + 1.6,
            "Vendor composition (n)",
            ha="center", va="bottom",
            fontsize=9, fontweight="bold", color="#333", zorder=5)

    counts = [("GE", 37), ("Philips", 50), ("Siemens", 180)]
    total  = sum(c for _, c in counts)
    xc = bar_x0
    for v, n in counts:
        seg_w = (bar_x1 - bar_x0) * n / total
        ax.add_patch(Rectangle((xc, bar_y), seg_w, bar_h,
                               facecolor=VENDOR_COLOR[v],
                               edgecolor="black", linewidth=0.6, zorder=3))
        ax.text(xc + seg_w / 2, bar_y + bar_h / 2,
                f"{v}   n = {n}",
                ha="center", va="center",
                fontsize=8.4, color="white", fontweight="bold",
                zorder=4)
        xc += seg_w


def draw_panel_B(ax, x0, y0, w, h):
    panel_frame(ax, x0, y0, w, h, "B", "Biomarker extraction at C2-C3",
                subtitle="Spinal Cord Toolbox v6.4  ·  8 per-subject scalars")
    groups = [
        ("T2w", ["Cord CSA"]),
        ("T2*", ["GM CSA"]),
        ("MT",  ["MTR", "MTsat"]),
        ("DTI", ["FA", "MD", "AD", "RD"]),
    ]
    pad_x = 3.0
    region_w = w - 2 * pad_x
    col_w = (region_w - 3 * 1.2) / 4
    top_y = y0 + h - 6.5
    bot_y = y0 + 2.0
    for i, (acq, biomarkers) in enumerate(groups):
        cx = x0 + pad_x + i * (col_w + 1.2)
        chip(ax, cx, top_y - 3.0, col_w, 2.8,
             acq, face=ACQ_COLOR[acq], text_color="white",
             fontsize=9.5, fontweight="bold")
        bm_top = top_y - 4.0
        bm_bot = bot_y
        n_bm = len(biomarkers)
        gap = 0.5
        bm_h = (bm_top - bm_bot - (n_bm - 1) * gap) / n_bm
        bm_h = max(min(bm_h, 2.4), 1.4)
        stack_h = n_bm * bm_h + (n_bm - 1) * gap
        first_y = bm_top - (bm_top - bm_bot - stack_h) / 2 - bm_h
        for j, bm in enumerate(biomarkers):
            by = first_y - j * (bm_h + gap)
            chip(ax, cx, by, col_w, bm_h,
                 bm, face="#ffffff", edge="#444", fontsize=8.5)


def draw_panel_C(ax, x0, y0, w, h):
    panel_frame(ax, x0, y0, w, h, "C", "Harmonization methods",
                subtitle="ordered by modeled vendor components")
    methods = [
        ("LME",          "mean only",        "(random intercept)"),
        ("ComBat",       "+ scale",          "(per-metric, EB)"),
        ("ComBat-joint", "+ joint EB",       "(across metrics)"),
        ("CovBat",       "+ covariance",     "(PCA + ComBat)"),
        ("RELIEF",       "+ latent factors", "(ICA + BH-FDR)"),
    ]
    pad_x = 3.0
    region_w = w - 2 * pad_x
    n = len(methods)
    gap = 1.0
    box_w = (region_w - (n - 1) * gap) / n
    header_h = 3.0
    desc_h   = 5.2
    block_h  = header_h + desc_h + 0.3
    block_y  = y0 + (h - block_h) / 2 - 0.6
    for i, (name, line1, line2) in enumerate(methods):
        bx = x0 + pad_x + i * (box_w + gap)
        col = METHOD_COLOR[name]
        chip(ax, bx, block_y + desc_h + 0.3, box_w, header_h,
             name, face=col, text_color="white",
             fontsize=9.5, fontweight="bold")
        chip(ax, bx, block_y, box_w, desc_h,
             f"{line1}\n{line2}",
             face="white", edge="#444", fontsize=8.0)
    arr_y = block_y - 2.0
    a = FancyArrowPatch((x0 + pad_x + 0.5,            arr_y),
                        (x0 + pad_x + region_w - 0.5, arr_y),
                        arrowstyle="-|>", mutation_scale=16,
                        color="#666", linewidth=1.3, zorder=2)
    ax.add_patch(a)
    ax.text(x0 + pad_x + region_w / 2, arr_y - 1.4,
            "Increasing model complexity",
            ha="center", va="top",
            fontsize=8.5, color="#444", style="italic", zorder=4)


def draw_panel_D(ax, x0, y0, w, h):
    panel_frame(ax, x0, y0, w, h, "D", "Evaluation framework",
                subtitle="three orthogonal endpoints")
    rows = [
        ("Vendor removal",       r"$\eta^{2}$  ·  ICC  ·  PERMANOVA $R^{2}$",
         WONG["vermilion"]),
        ("Biology recovery",     r"age $R^{2}$  ·  sex $R^{2}$",
         WONG["green"]),
        ("Subject preservation", r"within-vendor $r_{\mathrm{pre,post}}$",
         WONG["pink"]),
    ]
    pad_x = 3.0
    region_w = w - 2 * pad_x
    top_y = y0 + h - 6.0
    bot_y = y0 + 1.8
    region_h = top_y - bot_y
    gap = 0.6
    n = len(rows)
    row_h = (region_h - (n - 1) * gap) / n
    for i, (name, metrics, col) in enumerate(rows):
        ry = top_y - (i + 1) * row_h - i * gap
        ax.add_patch(Rectangle((x0 + pad_x, ry), 1.7, row_h,
                               facecolor=col, edgecolor="none", zorder=3))
        chip(ax, x0 + pad_x + 2.2, ry, region_w - 2.2, row_h,
             f"{name}     {metrics}",
             face="white", edge="#444", fontsize=9.0)


def draw_fig1a(out_dir: Path):
    """Fig 1a: Panels A + B (Dataset & Biomarkers) — landscape"""
    fig = plt.figure(figsize=(7.0, 4.0))
    ax = fig.add_subplot(111)
    ax.set_xlim(0, 100); ax.set_ylim(60, 126)
    ax.set_aspect("equal")
    ax.axis("off")

    fig.text(0.5, 0.975,
             "Study design: data and biomarker extraction",
             ha="center", va="top",
             fontsize=11, fontweight="bold")

    margin = 4
    pw = 100 - 2 * margin
    h_A, h_B = 28, 22
    arrow_gap = 3.5

    y_top = 122
    y_A = y_top - h_A
    y_B = y_A - arrow_gap - h_B

    draw_panel_A(ax, margin, y_A, pw, h_A)
    draw_panel_B(ax, margin, y_B, pw, h_B)

    xc = 50
    vertical_arrow(ax, xc, y_A - 0.3, y_B + h_B + 0.3)

    fig.tight_layout(pad=0.5)
    fig.savefig(out_dir / "Fig1a_dataset_biomarkers.pdf", dpi=300)
    fig.savefig(out_dir / "Fig1a_dataset_biomarkers.png", dpi=300)
    plt.close(fig)
    print(f"[saved] {out_dir / 'Fig1a_dataset_biomarkers.pdf'}")
    print(f"[saved] {out_dir / 'Fig1a_dataset_biomarkers.png'}")


def draw_fig1b(out_dir: Path):
    """Fig 1b: Panels C + D (Harmonization & Evaluation) — landscape"""
    fig = plt.figure(figsize=(7.0, 3.5))
    ax = fig.add_subplot(111)
    ax.set_xlim(0, 100); ax.set_ylim(10, 65)
    ax.set_aspect("equal")
    ax.axis("off")

    fig.text(0.5, 0.975,
             "Study design: harmonization methods and evaluation",
             ha="center", va="top",
             fontsize=11, fontweight="bold")

    margin = 4
    pw = 100 - 2 * margin
    h_C, h_D = 24, 20
    arrow_gap = 3.5

    y_top = 61
    y_C = y_top - h_C
    y_D = y_C - arrow_gap - h_D

    draw_panel_C(ax, margin, y_C, pw, h_C)
    draw_panel_D(ax, margin, y_D, pw, h_D)

    xc = 50
    vertical_arrow(ax, xc, y_C - 0.3, y_D + h_D + 0.3)

    fig.tight_layout(pad=0.5)
    fig.savefig(out_dir / "Fig1b_methods_evaluation.pdf", dpi=300)
    fig.savefig(out_dir / "Fig1b_methods_evaluation.png", dpi=300)
    plt.close(fig)
    print(f"[saved] {out_dir / 'Fig1b_methods_evaluation.pdf'}")
    print(f"[saved] {out_dir / 'Fig1b_methods_evaluation.png'}")


def main(out_dir: Path):
    out_dir.mkdir(parents=True, exist_ok=True)
    draw_fig1a(out_dir)
    draw_fig1b(out_dir)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()
    main(args.out)
