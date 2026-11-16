#!/usr/bin/env python3
"""Generate Fig1a and Fig1b as clean, hand-drawn SVG schematic diagrams.

These replace the matplotlib-drawn schematics that had layout problems.
SVG is converted to PDF and PNG via cairosvg.
"""

from pathlib import Path
from svglib.svglib import svg2rlg
from reportlab.graphics import renderPDF

# PyMuPDF for PNG preview (optional; PDF is the primary deliverable)
try:
    import fitz
    _HAS_FITZ = True
except Exception:
    _HAS_FITZ = False

OUT_DIR = Path("C:/Users/admin/WorkBuddy/2026-07-06-12-34-04/figures")
OUT_DIR.mkdir(parents=True, exist_ok=True)

# ── Color palette ──
C_GE       = "#2958A3"
C_PHILIPS  = "#6AAD4A"
C_SIEMENS  = "#D55E00"
C_PINK     = "#CC79A7"
C_BLUE     = "#0072B2"
C_SKY      = "#56B4E9"
C_ORANGE   = "#E69F00"
C_GREEN    = "#009E73"
C_RED      = "#D55E00"
C_GREY_BG  = "#F5F5F5"
C_GREY_TXT = "#555555"
C_DARK     = "#222222"
C_WHITE    = "#FFFFFF"
C_BORDER   = "#CCCCCC"

METHOD_COLORS = {
    "LME/FE":       C_PINK,
    "ComBat":       C_BLUE,
    "ComBat-joint": C_SKY,
    "RELIEF":       C_ORANGE,
    "CovBat":       C_GREEN,
}


def svg_wrap(body, width, height):
    return f'''<?xml version="1.0" encoding="UTF-8"?>
<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}"
     viewBox="0 0 {width} {height}">
  <rect width="{width}" height="{height}" fill="white"/>
{body}
</svg>'''

def rect(x, y, w, h, rx=6, fill="white", stroke=C_BORDER, sw=1.5):
    return f'  <rect x="{x}" y="{y}" width="{w}" height="{h}" rx="{rx}" fill="{fill}" stroke="{stroke}" stroke-width="{sw}"/>'

def text(x, y, content, size=12, color=C_DARK, anchor="middle", bold=False, font="Arial, sans-serif"):
    fw = 'font-weight="bold" ' if bold else ''
    return f'  <text x="{x}" y="{y}" {fw}font-family="{font}" font-size="{size}" fill="{color}" text-anchor="{anchor}">{content}</text>'

def multiline(x, y, lines, size=11, color=C_GREY_TXT, anchor="middle", lh=16):
    return "\n".join([text(x, y + i*lh, line, size=size, color=color, anchor=anchor) for i, line in enumerate(lines)])

def arrow(x1, y1, x2, y2, color=C_DARK, sw=1.8):
    """Simple arrow with triangle head."""
    import math
    dx = x2 - x1
    dy = y2 - y1
    length = math.hypot(dx, dy)
    if length < 1: return ""
    ux, uy = dx/length, dy/length
    hx = x2 - ux*10
    hy = y2 - uy*10
    return f'''  <line x1="{x1}" y1="{y1}" x2="{hx}" y2="{hy}" stroke="{color}" stroke-width="{sw}" stroke-linecap="round"/>
  <polygon points="{x2},{y2} {hx-uy*5},{hy+ux*5} {hx+uy*5},{hy-ux*5}" fill="{color}"/>'''


def label_letter(x, y, letter, color=C_DARK, size=20, bg_color=C_GREY_BG, r=18):
    """Draw a circular letter label (A, B, C, D)."""
    out = []
    out.append(f'  <circle cx="{x}" cy="{y}" r="{r}" fill="{bg_color}" stroke="{color}" stroke-width="2"/>')
    out.append(text(x, y + 7, letter, size=size, color=color, bold=True))
    return "\n".join(out)


def stacked_bar(x, y, w, h, segments, color_map, stroke_width=1.0, labels=True):
    """Draw a horizontal stacked bar chart.
    segments: list of (label, value, color)
    """
    total = sum(v for _, v, _ in segments)
    out = []
    cx = x
    for label, value, color in segments:
        seg_w = w * value / total
        out.append(f'  <rect x="{cx}" y="{y}" width="{seg_w}" height="{h}" fill="{color}" stroke="white" stroke-width="1"/>')
        # Label with count and percentage
        pct = value / total * 100
        if seg_w > 30 and labels:
            out.append(text(cx + seg_w/2, y + h/2 - 4, f"{value}", size=10, color="white", bold=True))
            out.append(text(cx + seg_w/2, y + h/2 + 8, f"({pct:.0f}%)", size=9, color="white"))
        cx += seg_w
    # Border
    out.append(f'  <rect x="{x}" y="{y}" width="{w}" height="{h}" rx="4" fill="none" stroke="{C_DARK}" stroke-width="{stroke_width}"/>')
    return "\n".join(out)


def save_and_convert(svg, name, W, H):
    """Save SVG and convert to PDF (required) and PNG preview (optional, via PyMuPDF)."""
    svg_path = OUT_DIR / f"{name}.svg"
    svg_path.write_text(svg, encoding="utf-8")
    drawing = svg2rlg(str(svg_path))
    if drawing is None:
        raise RuntimeError(f"svg2rlg failed to parse {svg_path}")
    # PDF
    pdf_path = OUT_DIR / f"{name}.pdf"
    renderPDF.drawToFile(drawing, str(pdf_path))
    # PNG preview (2x size)
    png_path = OUT_DIR / f"{name}.png"
    if _HAS_FITZ:
        doc = fitz.open(str(pdf_path))
        page = doc[0]
        mat = fitz.Matrix(2, 2)
        pix = page.get_pixmap(matrix=mat, alpha=False)
        pix.save(str(png_path))
        doc.close()
        print(f"  {name} saved (SVG + PDF + PNG)")
    else:
        print(f"  {name} saved (SVG + PDF); PNG preview skipped (PyMuPDF not available)")


# ═══════════════════════════════════════════════════════════════
#  FIG 1a: Dataset & Biomarkers Overview
# ═══════════════════════════════════════════════════════════════
def gen_fig1a():
    W, H = 1000, 680
    parts = []

    # Title
    parts.append(text(W/2, 32, "Study Cohort and Quantitative MRI Biomarkers", size=18, color=C_DARK, bold=True))

    # ── Panel A: Cohort overview ──
    py = 60
    ph = 230
    panel_left = 30
    panel_right = W - 30
    panel_w = panel_right - panel_left
    # Background panel
    parts.append(rect(panel_left, py, panel_w, ph, rx=10, fill=C_GREY_BG, stroke=C_BORDER))
    parts.append(label_letter(60, py + 32, "a"))
    parts.append(text(92, py + 34, "Study Cohort", size=15, color=C_DARK, bold=True, anchor="start"))

    # Main dataset box (left)
    bx, by = 60, py + 55
    bw, bh = 280, 150
    parts.append(rect(bx, by, bw, bh, rx=8, fill=C_WHITE, stroke=C_BORDER, sw=1.5))
    parts.append(text(bx + bw/2, by + 28, "spine-generic dataset", size=14, color=C_DARK, bold=True))
    parts.append(text(bx + bw/2, by + 50, "release r20231212", size=11, color=C_GREY_TXT))
    parts.append(multiline(bx + bw/2, by + 72,
                           ["267 subjects", "43 sites", "13 scanner models", "3 vendors", "(full dataset)"],
                           size=11, color=C_DARK, lh=18))

    # Key stats (middle)
    stats = [
        ("HC", "188", "complete cases"),
        ("ALL", "246", "complete cases"),
    ]
    sx = 380
    for i, (name, num, note) in enumerate(stats):
        x0 = sx + i * 130
        parts.append(f'  <rect x="{x0}" y="{by+10}" width="115" height="130" rx="8" fill="white" stroke="{C_BORDER}" stroke-width="1.2"/>')
        parts.append(text(x0+57, by+36, name, size=11, color=C_GREY_TXT, bold=True))
        parts.append(text(x0+57, by+68, num, size=26, color=C_DARK, bold=True))
        parts.append(text(x0+57, by+95, note, size=10, color=C_GREY_TXT))

    # Vendor stacked bar (right)
    bar_x = 660
    bar_y = by + 30
    bar_w = 290
    bar_h = 28
    parts.append(text(bar_x + bar_w/2, by + 18, "Vendor composition", size=12, color=C_DARK, bold=True))
    segments = [
        ("GE", 37, C_GE),
        ("Philips", 50, C_PHILIPS),
        ("Siemens", 180, C_SIEMENS),
    ]
    parts.append(stacked_bar(bar_x, bar_y, bar_w, bar_h, segments, None))
    # Legend below bar
    leg_y = bar_y + 45
    parts.append(f'  <rect x="{bar_x}" y="{leg_y}" width="14" height="14" rx="3" fill="{C_GE}"/>')
    parts.append(text(bar_x + 20, leg_y + 11, "GE", size=10, color=C_DARK, anchor="start"))
    parts.append(f'  <rect x="{bar_x + 70}" y="{leg_y}" width="14" height="14" rx="3" fill="{C_PHILIPS}"/>')
    parts.append(text(bar_x + 90, leg_y + 11, "Philips", size=10, color=C_DARK, anchor="start"))
    parts.append(f'  <rect x="{bar_x + 170}" y="{leg_y}" width="14" height="14" rx="3" fill="{C_SIEMENS}"/>')
    parts.append(text(bar_x + 190, leg_y + 11, "Siemens", size=10, color=C_DARK, anchor="start"))

    # Cohort definition note
    parts.append(text(bar_x + bar_w/2, by + 125, "Two pre-specified cohorts: HC (healthy controls) and ALL (full cohort)", size=10, color=C_GREY_TXT))

    # ── Panel B: 8 Biomarker cards ──
    by2 = py + ph + 50
    parts.append(label_letter(60, by2 - 12, "b"))
    parts.append(text(92, by2 - 10, "Eight qMRI Biomarkers (seven sampled at C2–C3)", size=15, color=C_DARK, bold=True, anchor="start"))

    biomarkers = [
        ("T2w CSA", "SC cross-sectional area", "T2w", "Morphometry", "#1B9E77"),
        ("T2* GM CSA", "Gray matter CSA", "T2*w", "GM morphometry", "#D95F02"),
        ("MTR", "Magnetization transfer ratio", "MT", "Myelin", "#7570B3"),
        ("MTsat", "MT saturation index", "MT", "Macromolecules", "#E7298A"),
        ("FA", "Fractional anisotropy", "DTI", "White-matter coherence", "#66A61E"),
        ("MD", "Mean diffusivity", "DTI", "Tissue microstructure", "#E6AB02"),
        ("AD", "Axial diffusivity", "DTI", "Axonal integrity", "#A6761D"),
        ("RD", "Radial diffusivity", "DTI", "Myelin condition", "#666666"),
    ]

    card_w, card_h = 210, 110
    gap_x, gap_y = 22, 18
    cols = 4
    start_x = (W - (cols*card_w + (cols-1)*gap_x)) / 2

    for i, (abbr, full, family, category, color) in enumerate(biomarkers):
        row, col = i // cols, i % cols
        cx = start_x + col*(card_w + gap_x)
        cy = by2 + 18 + row*(card_h + gap_y)

        # Card background
        parts.append(f'  <rect x="{cx}" y="{cy}" width="{card_w}" height="{card_h}" rx="8" fill="white" stroke="{color}" stroke-width="1.8"/>')
        # Color bar at top
        parts.append(f'  <rect x="{cx+8}" y="{cy+8}" width="{card_w-16}" height="5" rx="3" fill="{color}" opacity="0.7"/>')
        # Abbreviation
        parts.append(text(cx+card_w/2, cy+33, abbr, size=15, color=color, bold=True))
        # Full name
        parts.append(text(cx+card_w/2, cy+52, full, size=11, color=C_GREY_TXT))
        # Family tag
        parts.append(f'  <rect x="{cx+10}" y="{cy+card_h-34}" width="{card_w-20}" height="20" rx="10" fill="{color}" opacity="0.12"/>')
        parts.append(text(cx+card_w/2, cy+card_h-20, f"{family} · {category}", size=10, color=color))

    body = "\n".join(parts)
    svg = svg_wrap(body, W, H)
    save_and_convert(svg, "Fig1a_dataset_biomarkers", W, H)


# ═══════════════════════════════════════════════════════════════
#  FIG 1b: Harmonization & Evaluation Framework
# ═══════════════════════════════════════════════════════════════
def gen_fig1b():
    W, H = 1100, 560
    parts = []

    # Title
    parts.append(text(W/2, 32, "Harmonization and Evaluation Framework", size=18, color=C_DARK, bold=True))

    # ── Panel C: Harmonization Pipeline ──
    py = 65
    ph = 110

    parts.append(label_letter(40, py - 15, "C"))
    parts.append(text(72, py - 13, "Harmonization Pipeline", size=15, color=C_DARK, bold=True, anchor="start"))

    # Input box
    ix, iw = 40, 110
    parts.append(rect(ix, py, iw, ph, rx=8, fill=C_GREY_BG, stroke=C_BORDER, sw=1.5))
    parts.append(text(ix+iw/2, py+25, "Input", size=13, color=C_DARK, bold=True))
    parts.append(multiline(ix+iw/2, py+45,
                           ["8 biomarkers", "267 subjects", "44 sites"],
                           size=10, color=C_GREY_TXT, lh=15))

    # Arrow Input → first method
    parts.append(arrow(ix+iw, py+ph/2, 185, py+ph/2, C_DARK, 2.0))

    # Method boxes
    methods = [
        ("LME/FE", "Linear Mixed /", "Fixed Effects"),
        ("ComBat", "Empirical Bayes", "per metric"),
        ("ComBat-joint", "Joint EB", "all 8 metrics"),
        ("RELIEF", "ICA + BH-FDR", "components"),
        ("CovBat", "PCA + ComBat", "covariance"),
    ]
    mw, mg = 130, 16
    mh = ph
    mx_start = 190

    for i, (name, line1, line2) in enumerate(methods):
        mx = mx_start + i*(mw + mg)
        color = METHOD_COLORS[name]
        # Box
        parts.append(rect(mx, py, mw, mh, rx=8, fill="white", stroke=color, sw=2.0))
        # Top accent
        parts.append(f'  <rect x="{mx+2}" y="{py+2}" width="{mw-4}" height="5" rx="3" fill="{color}" opacity="0.6"/>')
        # Name
        parts.append(text(mx+mw/2, py+28, name, size=13, color=color, bold=True))
        # Description
        parts.append(text(mx+mw/2, py+50, line1, size=10, color=C_GREY_TXT))
        parts.append(text(mx+mw/2, py+64, line2, size=10, color=C_GREY_TXT))

        # Arrow to next
        if i < len(methods) - 1:
            parts.append(arrow(mx+mw, py+ph/2, mx+mw+mg, py+ph/2, "#AAAAAA", 1.5))

    # Arrow last method → Output
    last_x = mx_start + (len(methods)-1)*(mw+mg) + mw
    parts.append(arrow(last_x, py+ph/2, last_x+25, py+ph/2, C_DARK, 2.0))

    # Output box
    ox = last_x + 25
    ow = 100
    parts.append(rect(ox, py, ow, ph, rx=8, fill="#E8F5E9", stroke=C_GREEN, sw=2.0))
    parts.append(text(ox+ow/2, py+25, "Output", size=13, color=C_GREEN, bold=True))
    parts.append(multiline(ox+ow/2, py+45, ["5 harmonized", "datasets"], size=10, color=C_GREY_TXT, lh=15))

    # ── Vertical link ──
    mid_y = py + ph + 45
    parts.append(arrow(W/2, py+ph+5, W/2, mid_y, C_DARK, 2.5))
    # Label on arrow
    parts.append(f'  <rect x="{W/2-140}" y="{py+ph+12}" width="280" height="24" rx="12" fill="white" stroke="{C_BORDER}" stroke-width="1"/>')
    parts.append(text(W/2, py+ph+28, "All outputs are evaluated by four orthogonal endpoints", size=11, color=C_DARK, bold=True))

    # ── Panel D: Evaluation Endpoints ──
    parts.append(label_letter(40, mid_y + 22, "D"))
    parts.append(text(72, mid_y + 24, "Orthogonal Evaluation Endpoints", size=15, color=C_DARK, bold=True, anchor="start"))

    # ── Evaluation row ──
    ey = mid_y + 50
    eh = 130
    evals = [
        ("Vendor η²", "Univariate vendor effect", "Lower = better", C_RED,
         "Residual per-biomarker vendor bias"),
        ("PERMANOVA R²", "Multivariate vendor effect", "Lower = better", C_BLUE,
         "Joint 8-biomarker vendor structure"),
        ("r_pre,post", "Per-subject rank", "Higher = better", C_GREEN,
         "Pre- vs post-harmonization correlation"),
        ("Age + Sex R²", "Biological signal", "Higher = better", C_ORANGE,
         "Age/sex variance recovery"),
    ]
    ew, eg = 240, 20
    ex_start = (W - (4*ew + 3*eg)) / 2

    for i, (name, line1, compare, color, long_desc) in enumerate(evals):
        ex = ex_start + i*(ew + eg)

        # Arrow from central vertical line to each box
        parts.append(f'  <line x1="{W/2}" y1="{mid_y+5}" x2="{ex+ew/2}" y2="{ey-10}" stroke="{color}" stroke-width="1.2" opacity="0.5" stroke-dasharray="4,4"/>')
        parts.append(arrow(ex+ew/2-0.1, ey-10, ex+ew/2, ey, color, 1.5))

        # Box
        parts.append(rect(ex, ey, ew, eh, rx=10, fill="white", stroke=color, sw=2.0))
        # Top accent
        parts.append(f'  <rect x="{ex+2}" y="{ey+2}" width="{ew-4}" height="5" rx="3" fill="{color}" opacity="0.6"/>')

        # Name
        parts.append(text(ex+ew/2, ey+28, name, size=16, color=color, bold=True))
        # Short description
        parts.append(text(ex+ew/2, ey+50, line1, size=11, color=C_GREY_TXT))
        # Compare
        parts.append(f'  <rect x="{ex+ew/2-55}" y="{ey+60}" width="110" height="18" rx="9" fill="{color}" opacity="0.12"/>')
        parts.append(text(ex+ew/2, ey+72, compare, size=10, color=color))
        # Divider
        parts.append(f'  <line x1="{ex+30}" y1="{ey+85}" x2="{ex+ew-30}" y2="{ey+85}" stroke="{C_BORDER}" stroke-width="0.8"/>')
        # Long description
        parts.append(text(ex+ew/2, ey+103, long_desc, size=10, color="#777777"))

    body = "\n".join(parts)
    svg = svg_wrap(body, W, H)
    save_and_convert(svg, "Fig1b_methods_evaluation", W, H)


if __name__ == "__main__":
    print("Generating Fig1a...")
    gen_fig1a()
    print("Generating Fig1b...")
    gen_fig1b()
    print("Done.")
