# Spine-Generic QMRI Harmonization

Companion code and data for:

> Zhu Y, Wu D. **Scanner Vendor Harmonisation of Spinal Cord Quantitative MRI Biomarkers: A Five-Method Benchmark.** (under review).

## Overview

This repository contains all processing scripts, harmonization implementations, evaluation code, figure-generation pipelines, and biomarker tables supporting the above manuscript.

## Dataset

The analysis uses the **spine-generic multi-subject dataset** (release r20231212, archived at Zenodo DOI: 10.5281/zenodo.4299140), publicly available at:
- https://github.com/spine-generic/data-multi-subject

- 267 subjects, 44 sites, 3 vendors (GE, Philips, Siemens)
- 8 qMRI biomarkers at C2–C3: T2w CSA, T2\* GM CSA, MTR, MTsat, FA, MD, AD, RD

## Repository Structure

```
.
├── code/                           # Analysis and figure-generation scripts
│   ├── lovo_validation.py          # Leave-one-vendor-out (LOVO) external validation
│   ├── generate_lovo_figure.py     # LOVO validation figure generation
│   ├── generate_all_figures_final.py  # Final manuscript figure pipeline (Figs 1–5, FigS)
│   ├── fix_lovo_rprepost.py        # LOVO R-pre/R-post statistic computation
│   ├── analyze_spine_generic.py    # [DEPRECATED] Early exploratory T2*/T2w analysis
│   ├── resize_figures.py           # Figure DPI/format conversion
│   ├── step5_8_fig1a_fig1b_split.py  # Figure 1 generation
│   ├── verify_vendor_mapping.py    # Vendor mapping cross-check utility
│   ├── verify_master_csv.py        # Master CSV integrity check
│   ├── check_paper_stats.py        # Paper statistics validator
│   ├── check_sample_size.py        # Cohort sample-size audit
│   ├── check_contamination.py      # Harmonization residual-contamination check
│   ├── check_figure_quality.py     # Figure quality diagnostics
│   ├── patch_figure1_pdf.py        # Figure 1 vector patching
│   ├── patch_all_figures_pdf.py    # Batch vector patching of figures
│   ├── regenerate_all_figures.py   # Full figure regeneration (v1)
│   ├── regenerate_figures_v2.py    # Full figure regeneration (v2)
│   ├── sync_html_md.py             # LaTeX → HTML/MD synchronization
│   ├── check_abstract.py           # Abstract word-count validator
│   ├── fix_paper.py                # Automated paper fixes
│   └── make_supplementary_pdf.py   # Supplementary PDF generator
├── data/                           # Raw and harmonized biomarker tables
│   ├── supplementary_table_S1.csv  # Per-biomarker evaluation table (96 rows)
│   ├── cv_gm_by_site.csv           # GM CSA CV by site
│   ├── cv_sc_by_site.csv           # SC CSA CV by site
│   ├── descriptive_stats_by_site.csv
│   ├── descriptive_stats_by_vendor.csv
│   ├── descriptive_stats_t2star.csv
│   ├── descriptive_stats_t2w.csv
│   ├── merged_csa_data.csv
│   ├── pairwise_vendor_comparisons.csv
│   ├── regression_results.txt
│   ├── t2star_csa_data.csv
│   ├── t2w_by_vertebral_level.csv
│   ├── t2w_csa_data.csv
│   ├── voxel_resolution_by_vendor.csv
│   ├── lovo_vs_internal_comparison.csv  # LOVO vs internal-validation results
│   ├── results_design_C.csv        # Design C balanced-subsample bootstrap results
│   └── summary_design_C.csv        # Design C summary statistics
├── figures/                        # All manuscript figures (PNG + PDF)
│   ├── Fig1_study_design_v3.{png,pdf}
│   ├── Fig2_baseline_effects.{png,pdf}
│   ├── Fig3_performance.{png,pdf}
│   ├── Fig4_tradeoff.{png,pdf}
│   ├── Fig5_rpre_post_marginals.{png,pdf}
│   ├── fig_lovo_validation.{png,pdf}
│   └── FigS_design_C_robustness.{png,pdf}
├── paper/                          # Manuscript source files
│   ├── spine_generic_harmonization_paper.tex   # Final manuscript (LaTeX source)
│   └── spine_generic_harmonization_paper.pdf
├── LICENSE                         # MIT License
└── README.md                       # This file
```

## Requirements

- Python 3.13+
- Packages: numpy, scipy, scikit-learn, pandas, matplotlib, statsmodels, lifelines
- Spinal Cord Toolbox (SCT) v7.4 — for raw biomarker extraction
- LaTeX (pdflatex) — for manuscript compilation

## Harmonization Methods

| Method | Description |
|--------|-------------|
| ComBat | Per-metric empirical Bayes batch correction |
| ComBat-joint | Joint estimation across all biomarkers |
| CovBat | Covariance batch harmonization |
| RELIEF | Removal of latent inter-scanner effects via ICA |
| LME/FE | Linear mixed-effects / fixed-effects hybrid model |

## Citation

If you use this code or data, please cite:

```bibtex
@misc{zhu2026harmonization,
  author = {Zhu, Yilin and Wu, Dengke},
  title  = {Spine-Generic QMRI Harmonization},
  year   = {2026},
  doi    = {10.5281/zenodo.21264267},
  url    = {https://github.com/wudengke2010/spine-generic-harmonization}
}
```

## License

- **Code**: MIT License (see [LICENSE](LICENSE))
- **Data**: CC BY 4.0 (derived from spine-generic, which is openly available)
- **Figures**: CC BY 4.0

## Contact

- **Yilin Zhu** — Changsha Hospital of Traditional Chinese Medicine (Changsha No. 8 Hospital), Changsha City, Hunan Province, China
- **Dengke Wu** (corresponding) — Department of Emergency Medicine, and Emergency Medicine and Difficult Diseases Institute, The Second Xiangya Hospital of Central South University, Changsha 410011, Hunan, China — wudk2010@csu.edu.cn
