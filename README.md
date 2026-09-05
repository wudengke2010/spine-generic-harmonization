# Spine-Generic QMRI Harmonization

Companion code and data for:

> Zhu Y, Wu D. **Scanner Vendor Harmonisation of Spinal Cord Quantitative MRI Biomarkers: A Five-Method Benchmark.** (under review).

## Overview

This repository contains all processing scripts, harmonization implementations, evaluation code, figure-generation pipelines, and biomarker tables supporting the above manuscript.

## Dataset

The analysis uses the **spine-generic multi-subject dataset** (release r20231212, archived at Zenodo DOI: 10.5281/zenodo.4299140), publicly available at:
- https://github.com/spine-generic/data-multi-subject

- 267 subjects, 43 sites (institution level, participants.tsv), 3 vendors (GE, Philips, Siemens)
- 8 qMRI biomarkers at C2–C3: T2w CSA, T2\* GM CSA, MTR, MTsat, FA, MD, AD, RD

**Per-subject measurements** (`data/per_subject/`) are the direct SCT v7.4 output computed from the spine-generic derivatives: `csa_t2star/` holds one summary CSV per subject (SC/GM/WM CSA), and `csa_t2w/` holds the slice-level `sct_process_segmentation` output (per-slice area, angles, diameters). These are the raw inputs consumed by the analysis scripts in `code/`.

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
│   ├── ml_vendor_auc.py            # Machine-learnability of the vendor signature (RF, repeated stratified 5-fold CV)
│   ├── ml_fold_internal.py         # Fold-internal (leak-free) harmonisation ML analysis (v1.4.0)
│   ├── lovo_permanova_z.py         # LOVO PERMANOVA with z-score standardisation (v1.4.0)
│   ├── step6_design_C_via_paper_scripts.py  # Design C balanced-subsample analysis (Siemens downsampled to 40, 20 seeds)
│   ├── gen_fig4_revised.py         # Figure 4 (LOVO z-standardised, log scale) regeneration
│   ├── bootstrap_ci_table2.py      # Stratified bootstrap 95% CIs for Table 2 effect sizes
│   ├── covbat_pc_report.py         # CovBat PC retention analysis (K, cumulative variance)
│   ├── cov_frobenius_fig.py        # FigS4: cross-metric correlation matrices + Frobenius distance heatmaps
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
│   ├── summary_design_C.csv        # Design C summary statistics
│   ├── ml_auc_summary.csv          # ML vendor-signature AUC, ALL cohort (N=246)
│   ├── ml_auc_summary_HC.csv       # ML vendor-signature AUC, HC validation cohort (N=188)
│   ├── ml_leakfree_summary.csv     # Fold-internal (leak-free) ML vendor AUC summary (v1.4.0, primary ML endpoint)
│   ├── ml_leakfree_long.csv        # Per-repetition fold-internal ML results
│   ├── ml_leakfree_paired.csv      # Paired bootstrap tests between methods
│   ├── ml_leakfree_vendorclass.csv # Per-vendor classification AUC
│   ├── ml_leakfree_hmci.csv        # Pathology-label HM AUC with CIs
│   ├── cc_fit_summary.csv          # Complete-case fit summary per method x cohort x biomarker (LME convergence, eta2)
│   ├── lovo_permanova_z.csv        # LOVO PERMANOVA with z-score standardisation (v1.4.0)
│   ├── lovo_permanova_z_identity_check.csv  # Identity check for the z-standardised LOVO analysis
│   ├── age_R2_significance.csv     # Age-association R^2 and p-values per method/cohort
│   ├── comparison_long_cc.csv      # Complete-case long-format comparison table
│   ├── ml_vendor_detection.csv     # Per-run RF vendor-detection AUC traces
│   ├── bootstrap_ci_table2.csv     # Bootstrap 95% CIs for Table 2 effect sizes
│   ├── covbat_pc_retention.csv     # CovBat retained PC count (K) per cohort
│   ├── covbat_pc_variance.csv      # Per-PC explained variance of CovBat residuals
│   ├── corr_matrices_HC.csv        # 8x8 cross-metric correlations per method (HC)
│   ├── corr_matrices_ALL.csv       # 8x8 cross-metric correlations per method (ALL)
│   ├── frobenius_distance_HC.csv   # Pairwise Frobenius distances between methods (HC)
│   ├── frobenius_distance_ALL.csv  # Pairwise Frobenius distances between methods (ALL)
│   └── per_subject/                # Per-subject SCT output (raw measurements)
│       ├── csa_t2star/             # 267 files: sub-*_csa_t2star.csv (subject-level T2* CSA)
│       └── csa_t2w/                # 263 files: sub-*_csa.csv (slice-level T2w SCT output)
├── figures/                        # All manuscript figures (PNG + PDF)
│   ├── Fig1_study_design_v3.{png,pdf}
│   ├── Fig2_baseline_effects.{png,pdf}
│   ├── Fig3_performance.{png,pdf}
│   ├── Fig4_tradeoff.{png,pdf}
│   ├── Fig5_rpre_post_marginals.{png,pdf}
│   ├── fig_lovo_validation.{png,pdf}
│   ├── FigS_design_C_robustness.{png,pdf}
│   └── FigS4_covariance.{pdf,png,svg}  # Cross-metric covariance + Frobenius distance (Supplementary Fig. S4)
├── paper/                          # Manuscript source files (HBM submission version)
│   ├── main_text.tex               # Main manuscript (LaTeX source)
│   ├── main_text.pdf               # Compiled manuscript (19 pages)
│   ├── figures/                    # Main-text figures (Fig1–Fig4, vector PDF)
│   └── supplementary/              # ESM.tex/pdf + supplementary figures FigS1–S4
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

## Version History

- **v1.4.2** — Third-revision (minor) submission package. All analyses unified on the complete-case cohorts (N=188 HC / N=246 ALL). Fold-internal (leak-free) machine-learnability promoted to primary ML endpoint (`ml_fold_internal.py`, `ml_leakfree_*.csv`); LOVO PERMANOVA z-score standardisation (`lovo_permanova_z.py`, `lovo_permanova_z.csv`); Design C protocol description corrected (Siemens restricted to 40 per seed, all GE 25 + Philips 33 retained, n≈98 per seed, 20 random seeds) and FigS2 regenerated from the authoritative results (LME/FE mean r = 0.843); site counts harmonised to institution level (43 dataset / 40 analysis); numerical corrections (RELIEF LOVO ratio 428; non-diffusion age p ≥ 0.18); ESM Table S8 (per-vendor AUC) rendering fixed; paper/ restructured to the submission manuscript (main_text + supplementary).
- **v1.4.0 / v1.4.1** — Second-round (major) revision: complete-case unified refit of all methods, leakage-free ML analysis, LOVO z-standardised PERMANOVA, ML promoted to primary endpoint, Figures 1/2/4 redrawn, Table 1 site counts corrected, ESM Table S8 added.
- **v1.3.0** — Complete-case unified analysis (N=188 HC / N=246 ALL): refit ComBat/LME, paired bootstrap (8/10 HC, 9/10 ALL), LOVO per-fold PERMANOVA, numeric audit fixes (ESM Table S6 Philips LME 0.93→0.92).
- **v1.2.0** — Added ML vendor-signature analysis (`ml_vendor_auc.py`), bootstrap 95% CIs for Table 2 (`bootstrap_ci_table2.py`), CovBat PC retention report (`covbat_pc_report.py`), and covariance/Frobenius supplementary figure (`cov_frobenius_fig.py`, FigS4). Manuscript updated to the submission version (machine-learnability endpoint, uncertainty quantification, CovBat PC retention).
- **v1.1.0** — LOVO validation code, Design C data, final manuscript, per-subject SCT measurements.

## Citation

If you use this code or data, please cite:

```bibtex
@misc{zhu2026harmonization,
  author = {Zhu, Yilin and Wu, Dengke},
  title  = {Spine-Generic QMRI Harmonization},
  year   = {2026},
  doi    = {10.5281/zenodo.21264266},
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
