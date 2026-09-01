"""
Generate Supplementary Materials PDF (Supplementary Table S1).
- Page 1: Title page with manuscript title and author list
- Pages 2+: Supplementary Table S1 (96 rows, 11 columns)
"""
import csv
import os
from fpdf import FPDF

CSV_PATH = r'C:\Users\admin\WorkBuddy\2026-07-06-12-34-04\supplementary_table_S1.csv'
PDF_PATH = r'C:\Users\admin\WorkBuddy\2026-07-06-12-34-04\supplementary_materials.pdf'

# Read CSV
with open(CSV_PATH, 'r', encoding='utf-8') as f:
    reader = csv.reader(f)
    rows = list(reader)

header = rows[0]
data = rows[1:]

# Column display names (shorter for table fit)
col_names = ['Cohort', 'Method', 'Biomarker', 'N',
             'Vendor F', 'Vendor p', 'Vendor eta2', 'Vendor ICC',
             'Age R2', 'Sex R2', 'r(pre,post)']

# Format numeric values
def fmt(val):
    try:
        f = float(val)
        if abs(f) < 0.001 and f != 0:
            return f'{f:.2e}'
        elif abs(f) >= 100:
            return f'{f:.1f}'
        elif abs(f) >= 10:
            return f'{f:.2f}'
        elif abs(f) >= 1:
            return f'{f:.3f}'
        else:
            return f'{f:.4f}'
    except:
        return val

# Create PDF
pdf = FPDF(orientation='L', unit='mm', format='A4')
pdf.set_auto_page_break(auto=True, margin=15)
pdf.set_margins(10, 10, 10)
page_w = 297 - 20  # A4 landscape width minus margins

# ===== PAGE 1: Title Page =====
pdf.add_page()
pdf.ln(40)

pdf.set_font('Helvetica', 'B', 16)
title = "Supplementary Materials"
pdf.cell(page_w, 8, title, align='C', new_x="LMARGIN", new_y="NEXT")
pdf.ln(10)

pdf.set_font('Helvetica', '', 12)
subtitle = ("Scanner Vendor Harmonisation of Spinal Cord Quantitative MRI "
            "Biomarkers: A Five-Method Benchmark")
pdf.multi_cell(page_w, 6, subtitle, align='C')
pdf.ln(8)

pdf.set_font('Helvetica', '', 11)
pdf.cell(page_w, 6, "Yilin Zhu1,* & Dengke Wu2", align='C', new_x="LMARGIN", new_y="NEXT")
pdf.ln(3)
pdf.set_font('Helvetica', '', 9)
pdf.multi_cell(page_w, 5, "1 Changsha Hospital of Traditional Chinese Medicine (Changsha No. 8 Hospital), Changsha, Hunan, China", align='C')
pdf.multi_cell(page_w, 5, "2 Department of Emergency Medicine, and Emergency Medicine and Difficult Diseases Institute, The Second Xiangya Hospital of Central South University, Changsha 410011, Hunan, China", align='C')
pdf.ln(3)
pdf.set_font('Helvetica', 'I', 9)
pdf.cell(page_w, 5, "* Corresponding author: Dengke Wu (wudk2010@csu.edu.cn)", align='C', new_x="LMARGIN", new_y="NEXT")
pdf.ln(15)

pdf.set_font('Helvetica', 'B', 12)
pdf.cell(page_w, 7, "Supplementary Table S1", align='C', new_x="LMARGIN", new_y="NEXT")
pdf.ln(3)
pdf.set_font('Helvetica', '', 10)
desc = ("Per-biomarker evaluation table: all metrics for 8 biomarkers x 5 methods x 2 cohorts "
        "(96 rows, 11 columns). Metrics include vendor F-statistic, vendor p-value, "
        "vendor eta-squared, vendor ICC, age semi-partial R-squared, "
        "sex semi-partial R-squared, and pre-post Pearson r.")
pdf.multi_cell(page_w, 5, desc, align='C')
pdf.ln(5)

pdf.set_font('Helvetica', 'I', 9)
abbrev = ("Abbreviations: ComBat = Combatting batch effects; CovBat = Covariance Batch harmonization; "
          "LME/FE = linear mixed-effects / fixed-effects hybrid model; RELIEF = Removal of Latent "
          "Inter-scanner Effects via independent components; ICC = intraclass correlation coefficient; "
          "eta2 = eta-squared; R2 = semi-partial R-squared.")
pdf.multi_cell(page_w, 5, abbrev, align='C')

# ===== PAGE 2+: Data Table =====
pdf.add_page()

# Column widths (total ~277mm for A4 landscape with 10mm margins)
col_widths = [16, 22, 28, 10, 22, 22, 24, 24, 22, 22, 28]
total_w = sum(col_widths)

# Header row
pdf.set_font('Helvetica', 'B', 7.5)
pdf.set_fill_color(200, 210, 220)
pdf.set_text_color(0)
for i, name in enumerate(col_names):
    pdf.cell(col_widths[i], 6, name, border=1, align='C', fill=True)
pdf.ln()

# Data rows
pdf.set_font('Helvetica', '', 7)
for row_idx, row in enumerate(data):
    # Alternate row colors
    if row_idx % 2 == 0:
        pdf.set_fill_color(248, 248, 252)
    else:
        pdf.set_fill_color(255, 255, 255)

    for i, val in enumerate(row):
        formatted = fmt(val) if i >= 3 else val
        align = 'C' if i in [0, 1, 3] else 'R'
        pdf.cell(col_widths[i], 5, formatted, border=1, align=align, fill=True)
    pdf.ln()

# Footer note
pdf.ln(5)
pdf.set_font('Helvetica', 'I', 8)
pdf.set_text_color(100)
pdf.multi_cell(0, 4,
    "Note: 'Original' rows show baseline metrics before harmonization. "
    "Vendor F and p-value from one-way ANOVA (vendor as factor). "
    "eta2 = effect size (vendor variance / total variance). "
    "ICC = intraclass correlation coefficient (vendor as random effect). "
    "Age and Sex R2 = semi-partial R-squared from multiple regression "
    "(age and sex as predictors, vendor as covariate). "
    "r(pre,post) = Pearson correlation between pre- and post-harmonization "
    "values within each vendor, averaged across vendors.")

# Save
pdf.output(PDF_PATH)
print(f"PDF generated: {PDF_PATH}")
print(f"Pages: {len(pdf.pages)}")
print(f"Data rows: {len(data)}")
print(f"Columns: {len(header)}")
