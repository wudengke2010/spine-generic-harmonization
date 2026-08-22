#!/usr/bin/env python3
"""Visual quality diagnostic for SCI figures."""
import os
from pathlib import Path
import numpy as np
from PIL import Image, ImageStat

FIG_DIR = Path("C:/Users/admin/WorkBuddy/2026-07-06-12-34-04/figures")
OUT_FILE = Path("C:/Users/admin/WorkBuddy/2026-07-06-12-34-04/figure_quality_diagnosis.md")

FIGURES = {
    "Fig1a_dataset_biomarkers": "Study design: data and biomarkers",
    "Fig1b_methods_evaluation": "Study design: methods and evaluation",
    "Fig2_baseline_effects": "Baseline vendor effects",
    "Fig3_performance": "Harmonization performance",
    "Fig4_tradeoff": "Trade-off scatter",
    "Fig5_rpre_post_marginals": "r_pre,post heatmap",
    "FigS_design_C_robustness": "Design C robustness",
}

COLUMNS = [
    "Figure", "Width(px)", "Height(px)", "DPI_est",
    "Width(in)", "Height(in)", "AspectRatio", "FileSize(MB)",
    "WhiteEdge_Top", "WhiteEdge_Bottom", "WhiteEdge_Left", "WhiteEdge_Right",
    "TextDensity_Center", "MeanBrightness", "Issues"
]


def estimate_dpi(img):
    # Try to read DPI from metadata
    try:
        dpi = img.info.get("dpi")
        if dpi and dpi[0] and dpi[0] > 0:
            return dpi[0]
    except Exception:
        pass
    # Fallback: assume saved at 300 dpi from matplotlib
    return 300


def edge_margin(img_arr, threshold=250, scan_lines=30):
    """Return margin (in pixels) where image is mostly white near each edge."""
    # For grayscale luminance
    if len(img_arr.shape) == 3:
        gray = np.mean(img_arr, axis=2)
    else:
        gray = img_arr
    h, w = gray.shape

    def white_fraction(col_or_row):
        return np.sum(col_or_row > threshold) / len(col_or_row)

    top = 0
    for i in range(min(scan_lines, h)):
        if white_fraction(gray[i, :]) > 0.95:
            top += 1
        else:
            break
    bottom = 0
    for i in range(min(scan_lines, h)):
        if white_fraction(gray[h-1-i, :]) > 0.95:
            bottom += 1
        else:
            break
    left = 0
    for i in range(min(scan_lines, w)):
        if white_fraction(gray[:, i]) > 0.95:
            left += 1
        else:
            break
    right = 0
    for i in range(min(scan_lines, w)):
        if white_fraction(gray[:, w-1-i]) > 0.95:
            right += 1
        else:
            break
    return top, bottom, left, right


def text_density_center(img_arr, frac=0.7):
    """Estimate text density in central region by edge/ink pixel fraction."""
    if len(img_arr.shape) == 3:
        gray = np.mean(img_arr, axis=2)
    else:
        gray = img_arr
    h, w = gray.shape
    y0, y1 = int(h*(1-frac)/2), int(h*(1+frac)/2)
    x0, x1 = int(w*(1-frac)/2), int(w*(1+frac)/2)
    center = gray[y0:y1, x0:x1]
    # Fraction of pixels that are not near-white
    density = np.mean((center < 200) & (center > 30))
    return density


def diagnose(path, label):
    img = Image.open(path).convert("RGB")
    arr = np.array(img)
    w, h = img.size
    dpi = estimate_dpi(img)
    size_mb = os.path.getsize(path) / (1024*1024)
    top, bottom, left, right = edge_margin(arr)
    density = text_density_center(arr)
    brightness = np.mean(arr)

    width_in = w / dpi
    height_in = h / dpi
    aspect = w / h

    issues = []
    if w < 1600 or h < 1200:
        issues.append("low resolution")
    if top < 5 or bottom < 5 or left < 5 or right < 5:
        issues.append("content near edge/clipping risk")
    if density > 0.35:
        issues.append("high text density, possible overlap")
    elif density < 0.05:
        issues.append("very low content density")

    issues_str = "; ".join(issues) if issues else "OK"
    return [
        label, w, h, dpi,
        f"{width_in:.2f}", f"{height_in:.2f}", f"{aspect:.2f}", f"{size_mb:.2f}",
        top, bottom, left, right,
        f"{density:.3f}", f"{brightness:.1f}", issues_str
    ]


def main():
    rows = []
    for name, label in FIGURES.items():
        png = FIG_DIR / f"{name}.png"
        if not png.exists():
            rows.append([label] + ["-"] * (len(COLUMNS)-2) + ["file missing"])
            continue
        rows.append(diagnose(png, label))

    # Markdown table
    lines = ["# Figure Quality Diagnosis Report", ""]
    lines.append("| " + " | ".join(COLUMNS) + " |")
    lines.append("| " + " | ".join(["---" for _ in COLUMNS]) + " |")
    for row in rows:
        lines.append("| " + " | ".join(str(x) for x in row) + " |")
    lines.append("")
    lines.append("""
## 说明
- **DPI_est**: 估计DPI（matplotlib默认300）。
- **WhiteEdge_X**: 距离边缘多少像素内仍是白色/空白，值越小说明内容越贴近边缘，存在被裁切风险。
- **TextDensity_Center**: 中央区域内非白色像素占比，越高说明内容越密集，>0.35 表示可能文字/图形重叠。
- **Issues**: 自动检测到的潜在问题。
""")

    OUT_FILE.write_text("\n".join(lines), encoding="utf-8")
    print(f"Report written to: {OUT_FILE}")

    # Also print a concise summary
    print("\n" + "\n".join(lines))


if __name__ == "__main__":
    main()
