# -*- coding: utf-8 -*-
"""
Regenerate step5_lovo/lovo_summary.csv (and lovo_results_long.csv) on the
AUTHORITATIVE complete-case sample (n=188), derived directly from the saved
lovo_combined_{fold}_{method}.csv files.

Background
----------
The shipped lovo_summary.csv was aggregated from an older per-biomarker run that
used pairwise deletion and stale harmonised columns. It disagrees with
lovo_permanova_z.csv for ComBat (1.84e-4 vs 2.46e-3, a factor of 13) and LME
(0.2169 vs 0.2008), and therefore with main-text Table 2.

Because z-scored Euclidean PERMANOVA R^2 is mathematically identical to the mean
per-metric eta^2 computed on the same complete-case sample, regenerating from the
combined files must reproduce lovo_permanova_z.csv exactly. That identity is used
below as the pass/fail check.
"""
import io
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

OUT_DIR = Path("E:/boshi/qm_harmonization_paper/results/step5_lovo")
BIOMARKERS = ["T2w_CSA", "GM_CSA_mm2", "MTR", "MTsat", "FA", "MD", "AD", "RD"]
BATCH_COL = "Manufacturer"
VENDORS = ["GE", "Philips", "Siemens"]
METHODS = ["Original", "ComBat", "ComBat-joint", "CovBat", "RELIEF", "LME"]
SUFFIX = {"Original": "", "ComBat": "_combat", "ComBat-joint": "_combatJ",
          "CovBat": "_covbat", "RELIEF": "_relief", "LME": "_lme"}


# --- helpers copied verbatim from code/lovo_validation.py -------------------
def vendor_stats(y, batch):
    sub = pd.concat([y.rename("y"), batch.rename("b")], axis=1).dropna()
    if sub["b"].nunique() < 2 or len(sub) < 4:
        return dict(F=np.nan, p=np.nan, eta2=np.nan, n=len(sub))
    groups = [g["y"].values for _, g in sub.groupby("b") if len(g) > 1]
    if len(groups) < 2:
        return dict(F=np.nan, p=np.nan, eta2=np.nan, n=len(sub))
    F, p = stats.f_oneway(*groups)
    grand = sub["y"].mean()
    SSB = sum(len(g) * (g.mean() - grand) ** 2 for g in groups)
    SST = ((sub["y"] - grand) ** 2).sum()
    return dict(F=F, p=p, eta2=SSB / SST if SST > 0 else np.nan, n=len(sub))


def bio_partial_R2(y, age, sex):
    sub = pd.concat([y.rename("y"), age.rename("Age"), sex.rename("Sex")],
                    axis=1).dropna()
    if len(sub) < 10:
        return dict(R2_age=np.nan, R2_sex=np.nan)
    sex_d = pd.get_dummies(sub["Sex"], drop_first=True).astype(float)
    age_v = sub["Age"].astype(float).values
    y_v = sub["y"].astype(float).values
    X_full = sm.add_constant(np.column_stack([age_v, sex_d.values]))
    X_ageonly = sm.add_constant(age_v)
    X_sexonly = sm.add_constant(sex_d.values)
    r2_full = sm.OLS(y_v, X_full).fit().rsquared
    return dict(
        R2_age=max(r2_full - sm.OLS(y_v, X_sexonly).fit().rsquared, 0.0),
        R2_sex=max(r2_full - sm.OLS(y_v, X_ageonly).fit().rsquared, 0.0),
    )


def main():
    long_rows = []
    for fold in VENDORS:
        for method in METHODS:
            suf = SUFFIX[method]
            path = OUT_DIR / f"lovo_combined_{fold}_{method}.csv"
            df = pd.read_csv(path)
            harm = [b + suf for b in BIOMARKERS]
            # >>> authoritative sample: complete cases over the 8 harmonised cols
            keep = df[harm + [BATCH_COL, "Age", "Sex"]].dropna().index
            sub = df.loc[keep]
            if len(sub) != 188:
                raise RuntimeError(f"{fold}/{method}: n={len(sub)} (expected 188)")
            for b in BIOMARKERS:
                col = b + suf
                vs = vendor_stats(sub[col], sub[BATCH_COL])
                te_mask = sub[BATCH_COL] == fold
                orig_te = sub.loc[te_mask, b]
                harm_te = sub.loc[te_mask, col]
                n_min = min(len(orig_te), len(harm_te))
                r_pp = (stats.pearsonr(orig_te.values[:n_min],
                                       harm_te.values[:n_min])[0]
                        if n_min >= 3 else np.nan)
                bio = bio_partial_R2(sub[col], sub["Age"], sub["Sex"])
                long_rows.append(dict(fold=fold, method=method, biomarker=b,
                                      n=vs["n"], vendor_eta2=vs["eta2"],
                                      vendor_F=vs["F"], vendor_p=vs["p"],
                                      r_pre_post_test=r_pp,
                                      age_R2=bio["R2_age"], sex_R2=bio["R2_sex"]))

    long_df = pd.DataFrame(long_rows)
    summary = (long_df.groupby(["fold", "method"])
               .agg(vendor_eta2_mean=("vendor_eta2", "mean"),
                    vendor_eta2_std=("vendor_eta2", "std"),
                    r_pre_post_test_mean=("r_pre_post_test", "mean"),
                    age_R2_mean=("age_R2", "mean"),
                    sex_R2_mean=("sex_R2", "mean"),
                    n_mean=("n", "mean")).reset_index())

    # ---------- verification: identity vs lovo_permanova_z.csv ----------
    z = pd.read_csv(OUT_DIR / "lovo_permanova_z.csv")
    print("=" * 74)
    print("VERIFICATION  mean per-metric eta2  vs  z-scored PERMANOVA R2")
    print("=" * 74)
    ok = True
    for m in METHODS:
        s = summary[summary.method == m].set_index("fold").vendor_eta2_mean
        r2 = z[z.method == m].set_index("fold").permanova_R2
        diff = (s - r2).abs().max()
        flag = "OK " if diff < 1e-6 else "FAIL"
        if diff >= 1e-6:
            ok = False
        print(f"  {flag} {m:13s} max|eta2 - R2| = {diff:.2e}")
    print()
    mean_across = summary.groupby("method").vendor_eta2_mean.mean()
    table2 = {"ComBat": 2.46e-3, "ComBat-joint": 4.42e-4, "CovBat": 1.01e-3,
              "RELIEF": 1.87e-1, "LME": 2.01e-1, "Original": 3.15e-1}
    print("VERIFICATION  vs main-text Table 2 (LOVO column)")
    for m in ["ComBat", "ComBat-joint", "CovBat", "RELIEF", "LME", "Original"]:
        mine = mean_across[m]
        t2 = table2[m]
        print("  %-13s regenerated %.4e   Table2 %.4e   reldiff %.1f%%"
              % (m, mine, t2, 100 * abs(mine - t2) / t2))
    print()
    print("old lovo_summary means (pairwise-deletion, stale):")
    old = pd.read_csv(OUT_DIR / "lovo_summary.csv").groupby("method").vendor_eta2_mean.mean()
    print(old.to_string())

    if not ok:
        print("\n*** IDENTITY FAILED -> files NOT written ***")
        return

    # ---------- write ----------
    for f in ["lovo_results_long.csv", "lovo_summary.csv"]:
        bak = OUT_DIR / f"{f}.bak_stale"
        if not bak.exists():
            (OUT_DIR / f).replace(bak)
    long_df.to_csv(OUT_DIR / "lovo_results_long.csv", index=False)
    summary.to_csv(OUT_DIR / "lovo_summary.csv", index=False)
    print("\n[saved] lovo_summary.csv (%d rows), lovo_results_long.csv (%d rows)"
          % (len(summary), len(long_df)))
    print("[backup] originals kept as <name>.bak_stale")

if __name__ == "__main__":
    main()
