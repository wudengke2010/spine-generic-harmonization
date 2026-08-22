#!/usr/bin/env python3
"""Patch the existing paper PDF by replacing the old Figure 1 images on page 8
with the newly generated Fig1a/Fig1b PNGs, preserving aspect ratio and fitting
them into the original figure area.
"""

from pathlib import Path
import fitz

WORKSPACE = Path("C:/Users/admin/WorkBuddy/2026-07-06-12-34-04")
SRC_PDF = WORKSPACE / "spine_generic_harmonization_paper.pdf"
OUT_PDF = WORKSPACE / "spine_generic_harmonization_paper_patched.pdf"
FIG1A = WORKSPACE / "figures" / "Fig1a_dataset_biomarkers.png"
FIG1B = WORKSPACE / "figures" / "Fig1b_methods_evaluation.png"

PAGE_INDEX = 7  # page 8 (0-indexed)

# Original figure area on page 8 (from PyMuPDF inspection)
FIG_TOP = 70.87
FIG_BOTTOM = 555.54
FIG_LEFT = 70.87
FIG_RIGHT = 524.41
FIG_WIDTH = FIG_RIGHT - FIG_LEFT  # 453.54
FIG_HEIGHT = FIG_BOTTOM - FIG_TOP  # 484.67

# New figure pixel dimensions and aspect ratios
fig1a_w, fig1a_h = 1000, 680
fig1b_w, fig1b_h = 1100, 560
ar1a = fig1a_w / fig1a_h
ar1b = fig1b_w / fig1b_h

# Determine display width that fits both figures within original height, with a 10 pt gap
gap = 10
# W/ar1a + W/ar1b + gap = FIG_HEIGHT  =>  W*(1/ar1a + 1/ar1b) = FIG_HEIGHT - gap
W = (FIG_HEIGHT - gap) / (1/ar1a + 1/ar1b)
# Center the figure column within the original text width
center_x = (FIG_LEFT + FIG_RIGHT) / 2
x0 = center_x - W / 2
x1 = center_x + W / 2

h1a = W / ar1a
h1b = W / ar1b

y1a_top = FIG_TOP
y1a_bottom = y1a_top + h1a
y1b_top = y1a_bottom + gap
y1b_bottom = y1b_top + h1b

print(f"Display width: {W:.1f} pt")
print(f"Fig1a rect: ({x0:.1f}, {y1a_top:.1f}, {x1:.1f}, {y1a_bottom:.1f})")
print(f"Fig1b rect: ({x0:.1f}, {y1b_top:.1f}, {x1:.1f}, {y1b_bottom:.1f})")

# Open PDF
doc = fitz.open(SRC_PDF)
page = doc[PAGE_INDEX]

# Draw white rectangle over the original figure area to erase old images
page.draw_rect(fitz.Rect(FIG_LEFT, FIG_TOP, FIG_RIGHT, FIG_BOTTOM),
               color=(1, 1, 1), fill=(1, 1, 1), width=0)

# Insert new images
page.insert_image(fitz.Rect(x0, y1a_top, x1, y1a_bottom), filename=str(FIG1A))
page.insert_image(fitz.Rect(x0, y1b_top, x1, y1b_bottom), filename=str(FIG1B))

# Save
OUT_PDF.unlink(missing_ok=True)
doc.save(OUT_PDF, garbage=4, deflate=True)
doc.close()
print(f"Saved patched PDF: {OUT_PDF}")
