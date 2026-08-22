#!/usr/bin/env python3
"""Patch the paper PDF by replacing Fig2–FigS with newly generated matplotlib figures.

Approach: For each figure page, draw a white rectangle over the old image area,
then insert the new PNG centered within the original rect, preserving aspect ratio.

Fig1a/Fig1b were already patched in a previous session — we skip page 8.
"""

from pathlib import Path
import fitz  # PyMuPDF

WORKSPACE = Path("C:/Users/admin/WorkBuddy/2026-07-06-12-34-04")
SRC_PDF = WORKSPACE / "spine_generic_harmonization_paper.pdf"
OUT_PDF = WORKSPACE / "spine_generic_harmonization_paper_final.pdf"
FIG_DIR = WORKSPACE / "figures"

# (page_index, figure_file, original_rect)
# original_rect = (x0, y0, x1, y1) from PDF inspection
FIGURES = [
    # Page 12 (idx 11) — Fig2
    (11, "Fig2_baseline_effects.png", (70.9, 70.9, 524.4, 469.2)),
    # Page 15 (idx 14) — Fig3
    (14, "Fig3_performance.png", (70.9, 176.8, 524.4, 502.5)),
    # Page 17 (idx 16) — Fig4
    (16, "Fig4_tradeoff.png", (82.2, 237.1, 513.1, 515.2)),
    # Page 18 (idx 17) — Fig5
    (17, "Fig5_rpre_post_marginals.png", (70.9, 115.7, 524.4, 444.8)),
    # Page 20 (idx 19) — FigS (Figure 6 / Supplementary)
    (19, "FigS_design_C_robustness.png", (70.9, 330.4, 524.4, 515.1)),
]


def fit_image_in_rect(img_w, img_h, rect_w, rect_h):
    """Compute display dimensions that fit image within rect, preserving aspect ratio."""
    ar = img_w / img_h
    rect_ar = rect_w / rect_h
    if ar > rect_ar:
        # Image is wider than rect — fit to width
        disp_w = rect_w
        disp_h = rect_w / ar
    else:
        # Image is taller than rect — fit to height
        disp_h = rect_h
        disp_w = rect_h * ar
    return disp_w, disp_h


def main():
    doc = fitz.open(SRC_PDF)
    print(f"Opened: {SRC_PDF.name} ({len(doc)} pages)")

    for page_idx, fig_name, (rx0, ry0, rx1, ry1) in FIGURES:
        fig_path = FIG_DIR / fig_name
        if not fig_path.exists():
            print(f"  WARNING: {fig_name} not found, skipping")
            continue

        page = doc[page_idx]
        rect_w = rx1 - rx0
        rect_h = ry1 - ry0

        # Get new image dimensions
        info = fitz.Pixmap(str(fig_path))
        img_w, img_h = info.width, info.height
        info = None  # free

        # Compute fitted display size
        disp_w, disp_h = fit_image_in_rect(img_w, img_h, rect_w, rect_h)

        # Center within original rect
        cx = (rx0 + rx1) / 2
        cy = (ry0 + ry1) / 2
        dx0 = cx - disp_w / 2
        dy0 = cy - disp_h / 2
        dx1 = cx + disp_w / 2
        dy1 = cy + disp_h / 2

        # Step 1: Draw white rectangle over old image area
        page.draw_rect(fitz.Rect(rx0, ry0, rx1, ry1),
                       color=(1, 1, 1), fill=(1, 1, 1), width=0)

        # Step 2: Insert new image
        page.insert_image(fitz.Rect(dx0, dy0, dx1, dy1), filename=str(fig_path))

        print(f"  Page {page_idx + 1}: {fig_name} "
              f"({img_w}x{img_h} -> {disp_w:.1f}x{disp_h:.1f} pt, "
              f"rect=({dx0:.1f}, {dy0:.1f}, {dx1:.1f}, {dy1:.1f}))")

    # Save
    OUT_PDF.unlink(missing_ok=True)
    doc.save(OUT_PDF, garbage=4, deflate=True)
    doc.close()
    print(f"\nSaved: {OUT_PDF}")
    print(f"Size: {OUT_PDF.stat().st_size / 1024:.0f} KB")


if __name__ == "__main__":
    main()
