# Spine-Generic QMRI Harmonization

Companion code and data for:

> Zhu Y, Wu D. **Scanner Vendor Harmonisation of Spinal Cord Quantitative MRI Biomarkers: A Five-Method Benchmark.** (under review).

## Overview

This repository contains all processing scripts, harmonization implementations, evaluation code, figure-generation pipelines, and biomarker tables supporting the above manuscript.

## Dataset

The analysis uses the **spine-generic multi-subject dataset** (release r20231212, archived at Zenodo DOI: 10.5281/zenodo.4299140), publicly available at:
- https://github.com/spine-generic/data-multi-subject

- 267 subjects, 43 sites (institution level, participants.tsv), 3 vendors (GE, Philips, Siemens)
- 8 qMRI biomarkers: T2w CSA, T2\* GM CSA, MTR, MTsat, FA, MD, AD, RD (seven sampled at C2–C3; GM-CSA averaged over the full cord — see manuscript Section 2.1)

**Per-subject measurements** (`data/per_subject/`) are the direct SCT v7.4 output computed from the spine-generic derivatives: `csa_t2star/` holds one summary CSV per subject (SC/GM/WM CSA), and `csa_t2w/` holds the slice-level `sct_process_segmentation` output (per-slice area, angles, diameters). These are the raw inputs consumed by the analysis scripts in `code/`.

## Repository Structure

```
.
├── code/                           # Analysis and figure-generation scripts
│   ├── lovo_validation.py          # Leave-one-vendor-out (LOVO) transductive validation
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
│   ├── lovo_summary.csv           # LOVO validation results by method and held-out vendor (v1.7.3 rerun)
│   ├── lovo_testvendor.csv         # Held-out-vendor shifts and correlations per fold/method/biomarker
│   ├── lovo_per_vendor.csv         # Per-vendor LOVO detail (fold-wise r ranges)
│   ├── results_design_C.csv        # Design C balanced-subsample bootstrap results
│   ├── summary_design_C.csv        # Design C summary statistics
│   ├── ml_auc_summary.csv          # ML vendor-signature AUC, ALL cohort (N=246)
│   ├── ml_auc_summary_HC.csv       # ML vendor-signature AUC, HC validation cohort (N=188)
│   ├── ml_leakfree_summary.csv     # Fold-internal (leak-free) ML vendor AUC summary (v1.4.0, primary ML endpoint)
│   ├── ml_leakfree_long.csv        # Per-repetition fold-internal ML results
│   ├── ml_leakfree_paired.csv      # Paired bootstrap tests between methods
│   ├── ml_leakfree_vendorclass.csv # Per-vendor classification AUC
│   ├── ml_leakfree_hmci.csv        # Compression-label AUC with CIs (exploratory; age-confounded)
│   ├── r1_age_confound.csv         # Age-confounding sensitivity of compression-label AUC (S9)
│   ├── r1_age_confound_infold.csv  # AUC after in-fold age/sex residualisation (S9)
│   ├── r1_relief_lovo_match.csv    # RELIEF LOVO under published vs distribution-matched protocols (S5)
│   ├── r1_site_harmonization.csv   # ComBat with site/scanner-model batch variables (S10)
│   ├── comparison_long.csv         # Per-biomarker metrics, all methods/cohorts, unified complete-case
│   ├── cc_summary.csv              # Method-level summary of comparison_long.csv
│   ├── demographics_completecase.csv # Complete-case demographics (Table 1 / S1)
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
│   ├── FigS4_covariance.{pdf,png,svg}  # Cross-metric correlation matrices + Frobenius distance (Supplementary Fig. S4)
│   └── GraphicalAbstract.pdf       # Graphical abstract (Table of Contents)
├── paper/                          # Manuscript source files (NMR in Biomedicine submission version)
│   ├── main_text.tex               # Main manuscript (LaTeX source)
│   ├── main_text.pdf               # Compiled manuscript (28 pages)
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
| ComBat | Per-metric location-and-scale correction without empirical-Bayes pooling (eb=False) |
| ComBat-joint | Joint estimation with a single EB prior shared across biomarkers (eb=True); covariance not explicitly modelled |
| CovBat | Covariance batch harmonization |
| RELIEF | Removal of latent inter-scanner effects via ICA |
| LME/FE | Linear mixed-effects / fixed-effects hybrid model |

## Version History

- **v1.7.9** — Data-consistency pass after a full numerical audit (285 checks against the authoritative results): ESM S5 ComBat description corrected to match the main text and the code (per-metric ComBat uses `eb=False` in both the internal and LOVO analyses; the stale claim that the internal analysis used `eb=True`, and that the LOVO ratio of 0.5 was an estimator effect, was removed); the Cover Letter statement on AI assistance aligned with the expanded disclosure in the Acknowledgements; the stale "54 derived output tables" count (correct at v1.7.3) removed from the Declarations, Cover Letter and Title Page; and the obsolete `data/lovo_vs_internal_comparison.csv` (v1.7.2 LOVO values that contradicted Table 3) deleted from the archived repository. No data, results, or conclusions changed.
- **v1.7.8** — Figure 4 layout fix: panel (a) and (b) legends moved to the upper-centre position so they no longer overlap the tall Original bars; the chance-level annotation in panel (a) relocated to the free mid-panel region. Cosmetic only — no data, numerical results, or conclusions changed.
- **v1.7.7** — Pre-resubmission consistency pass: ComBat EB description corrected to match the code (per-metric ComBat uses eb=False in both the internal and LOVO analyses; the earlier claim that the internal analysis used eb=True was wrong); ComBat-joint re-described (single EB prior shared across biomarkers; cross-metric covariance not explicitly modelled); "covariance structure" terminology aligned to "correlation structure" where correlation matrices are meant; graphical abstract redrawn to match the revised conclusions; Declarations consent-for-publication wording corrected; repository data cleaned (stale lovo_vs_internal_comparison.csv removed; manuscript CSVs documented; correspondence synced to the resubmission versions).
- **v1.7.6** — Word-by-word final check: LOVO fold-wise r range corrected to 0.84–0.95; LOVO PERMANOVA floor stated as p>0.99; missingness double-count clarified (21 unique subjects); introduction citations redistributed; ESM PERMANOVA citation fixed; Fig 1 British spellings.
- **v1.7.5** — Table 3 LOVO columns resynced to the protocol-symmetric v1.7.3 rerun (ComBat 1.84e-04, ratio 0.5; RELIEF 8.86e-04, ratio 2.0; LME/FE 1.99e-01, ratio 250); data/lovo_*.csv refreshed; figure panel titles aligned with captions.
- **v1.7.4** — Terminology pass: compression labels replace pathology labels throughout; "improved" softened to "numerically higher"; LOVO described as transductive validation; conclusions hedged.
- **v1.7.3** — LOVO protocol symmetry: RELIEF held-out-vendor pathway switched to test-vendor distribution matching; LOVO results rerun; tracked-changes PDF generated.
- **v1.5.1** — Fifth-revision (minor) submission package. Figure 3 (Pareto frontier) regenerated with the leader-line method labels removed (single figure-level legend; the previous redraw had introduced label overflow and a duplicated "CovBat" overprinting the x-axis title) and an in-panel "compressed x-axis range" note in panel (b) plus a matching caption sentence. Measured baseline η² of both Design C subsample designs reported with seed-wise ranges and the explicit √(1−η²) pooled-r ceiling (reduced-imbalance 0.374 → ≈0.79; true-equal-n 0.389 → ≈0.78), plus a per-seed verification that LME/FE pooled r never exceeded the per-seed ceiling in any of the 20 seeds of either design (ESM S6). Long-sentence restructuring in Methods (LOVO), Results (covariance, Design C, LOVO) and a three-way split of the Limitations paragraph.
- **v1.5.0** — Fourth-revision (minor) submission package. Design C subsample analysis rebuilt on the paper-canonical pooled-r metric (`step6_design_C_via_paper_scripts.py` v2; the previously reported within-vendor-based LME/FE r = 0.843 was not reproducible from the canonical code and is superseded): pooled-r hierarchy retained under both the reduced-imbalance design (n≈98/seed; LME/FE 0.779, ComBat-joint 0.773, RELIEF 0.760, CovBat 0.749, ComBat 0.745) and a new true-equal-n sensitivity design (`step6_design_C_equaln.py`, 25 subjects/vendor, n=75; LME/FE 0.772, ComBat-family 0.732–0.768); terminology "balanced conditions" → "reduced-imbalance subsample"; Frobenius distances normalised by ‖R_orig‖_F (HC 38.1–41.3%, ALL 40.2–43.0%); LOVO r = 0.85–0.97 disambiguated as per-fold range with fold means 0.90/0.95; Fig 1A site/scanner-model scope annotated (43 sites · 13 models, full dataset); Fig 4B baseline asterisks removed; Fig 3 caption panel letters lower-cased; Word manuscript regenerated (caption numbering prefixes, numbered section headings, recomputed app.xml statistics, pandoc-residue styles removed).
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
