# -*- coding: utf-8 -*-
"""Graphical Abstract for the NMR in Biomedicine submission (revision v1.7.x).
Left: multi-site spinal-cord QMRI -> vendor effects -> harmonisation schematic.
Right: five-method trade-off (vendor-effect removal vs subject-level preservation,
       HC cohort ordering), highlighting ComBat-joint as the recommended balance.
Outputs a clean 16:9 PNG (and PDF) suitable for the Table of Contents.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Ellipse
import os

OUT = "C:/Users/admin/WorkBuddy/2026-07-06-12-34-04/submission/figures/"
os.makedirs(OUT, exist_ok=True)

# ---- palette (jargon-free, readable at thumbnail) ----
TEAL   = "#0E6E6B"   # primary pipeline / methods
AMBER  = "#E08A1E"   # recommended-balance accent
RED    = "#C0392B"   # vendor effects
GRAY   = "#9AA0A6"
DGRAY  = "#3C4043"

fig = plt.figure(figsize=(16, 9), dpi=100)
fig.patch.set_facecolor("white")

# ============================================================
# LEFT  ~ schematic pipeline (0.00 - 0.60)
# ============================================================
ax0 = fig.add_axes([0.02, 0.06, 0.60, 0.88])
ax0.set_xlim(0, 10); ax0.set_ylim(0, 10); ax0.axis("off")

# --- header ---
ax0.text(0.4, 9.55, "Multi-site spinal cord quantitative MRI", fontsize=15,
         fontweight="bold", color=DGRAY, va="center")
ax0.text(0.4, 9.0, "43 sites  ·  3 vendors  ·  267 subjects", fontsize=12,
         color=GRAY, va="center", style="italic")

# --- top row: MR thumbnails per vendor ---
def mri_thumb(x, y, label, color):
    # spinal cord cross-section glyph (grey matter butterfly)
    ax0.add_patch(Ellipse((x, y), 1.3, 1.9, fc="#F1F3F4", ec=color, lw=1.6))
    ax0.add_patch(Ellipse((x-0.28, y+0.05), 0.32, 0.5, fc="#C6CACF", ec="none"))
    ax0.add_patch(Ellipse((x+0.28, y-0.05), 0.32, 0.5, fc="#C6CACF", ec="none"))
    ax0.text(x, y-1.15, label, fontsize=10.5, color=color, ha="center", fontweight="bold")

mri_thumb(1.6, 7.3, "Siemens", TEAL)
mri_thumb(4.0, 7.3, "GE",      TEAL)
mri_thumb(6.4, 7.3, "Philips", TEAL)

# --- vendor-effect arrow ---
ax0.annotate("", xy=(7.9, 6.6), xytext=(7.9, 7.9),
             arrowprops=dict(arrowstyle="-|>", color=RED, lw=2.4))
ax0.text(8.35, 7.25, "Vendor\nscanner\neffects", fontsize=11, color=RED,
         ha="left", va="center", fontweight="bold")

# --- harmonisation methods row ---
ax0.text(0.4, 5.4, "Harmonisation methods", fontsize=13.5, fontweight="bold", color=DGRAY)
methods = ["ComBat", "ComBat-joint", "CovBat", "RELIEF", "LME/FE"]
for i, m in enumerate(methods):
    x = 1.0 + i * 1.55
    w = 1.4
    fc = AMBER if m == "ComBat-joint" else "#FFFFFF"
    ec = AMBER if m == "ComBat-joint" else GRAY
    ax0.add_patch(FancyBboxPatch((x, 4.35), w-0.15, 0.75,
                  boxstyle="round,pad=0.02,rounding_size=0.12",
                  fc=fc, ec=ec, lw=2.0))
    ax0.text(x + (w-0.15)/2, 4.72, m, fontsize=10.5, color=DGRAY,
             ha="center", va="center", fontweight="bold")

# --- harmonised output ---
ax0.annotate("", xy=(9.0, 3.55), xytext=(2.0, 3.55),
             arrowprops=dict(arrowstyle="-|>", color=TEAL, lw=2.4))
ax0.text(0.4, 2.7, "Harmonised biomarkers", fontsize=13.5, fontweight="bold", color=DGRAY)
ax0.text(0.4, 2.05, "Vendor effects reduced, but residual vendor signatures remain.",
         fontsize=11.5, color=GRAY)

# --- bottom: preservation note ---
ax0.text(0.4, 1.05, "Subject-ranking preservation and age/sex associations were evaluated",
         fontsize=11.5, color=AMBER, fontweight="bold")

# ============================================================
# RIGHT ~ trade-off scatter (0.62 - 0.98)
# Positions follow the HC-cohort ordering: x = multivariate vendor
# removal (-log10 PERMANOVA R2), y = subject-level preservation (r).
# ============================================================
ax1 = fig.add_axes([0.63, 0.14, 0.35, 0.74])
ax1.set_facecolor("white")

# stylized schematic (qualitative ordering, not exact values)
# (x, y, color, marker, label_x, label_y, ha)
pts = {
    "Original":    (0.06, 1.00, GRAY,  "o", 0.06, 1.075, "center"),
    "LME/FE":      (0.16, 0.86, TEAL,  "o", 0.16, 0.94,  "center"),
    "ComBat-joint":(0.46, 0.80, AMBER, "D", 0.46, 0.87,  "center"),
    "RELIEF":      (0.54, 0.69, TEAL,  "o", 0.54, 0.585, "center"),
    "ComBat":      (0.60, 0.66, TEAL,  "o", 0.68, 0.66,  "left"),
    "CovBat":      (0.88, 0.55, TEAL,  "o", 0.88, 0.46,  "center"),
}
for name, (x, y, c, m, lx, ly, ha) in pts.items():
    ax1.scatter([x], [y], s=210, c=c, marker=m, zorder=3,
                edgecolor="white", linewidth=1.4)
    ax1.text(lx, ly, name, fontsize=10.5, ha=ha, va="center",
             color=DGRAY, fontweight="bold" if name == "ComBat-joint" else "normal")

# non-dominated envelope through all five methods (HC cohort)
ax1.plot([0.06, 0.16, 0.46, 0.54, 0.60, 0.88],
         [1.00, 0.86, 0.80, 0.69, 0.66, 0.55],
         "--", color=GRAY, lw=1.6, zorder=2)

# axes labels
ax1.set_xlabel("Vendor-effect removal  →", fontsize=12.5, color=DGRAY, fontweight="bold")
ax1.set_ylabel("Subject-level\npreservation  →", fontsize=12.5, color=DGRAY, fontweight="bold")
ax1.set_xlim(-0.02, 1.15); ax1.set_ylim(0.30, 1.22)
ax1.set_xticks([]); ax1.set_yticks([])
for s in ax1.spines.values():
    s.set_color(GRAY); s.set_linewidth(1.0)

ax1.set_title("Decontamination–preservation trade-off", fontsize=13,
              color=DGRAY, fontweight="bold", pad=10)

# ============================================================
# bottom banner
# ============================================================
fig.text(0.5, 0.035,
         "ComBat-joint offered a favourable balance in this dataset; LME/FE failed to centre "
         "the held-out vendor under the specified transductive rule",
         fontsize=12.5, color=DGRAY, ha="center", fontweight="bold")

for ext in ("png", "pdf"):
    fig.savefig(OUT + f"GraphicalAbstract.{ext}", dpi=300, bbox_inches="tight",
                facecolor="white")
plt.close(fig)
print("saved GraphicalAbstract.png / .pdf")
