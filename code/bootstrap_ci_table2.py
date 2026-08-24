"""
Bootstrap 95% CI for Table 2 (HC cohort).
==========================================
Faithfully replicates the original pipeline's N strategy:
  - vendor eta^2, r_pre_post, Age R^2: AVAILABLE-CASE per biomarker
    (N=203/194/190/198 etc., matching step4_eval_comparison.py)
  - PERMANOVA R^2: COMPLETE-CASE (N=188), z-scored ddof=1
    (matching step4_eval_multivariate.py)

Bootstrap: stratified nonparametric (B=1000), resampling subjects
WITHIN each scanner-manufacturer stratum (preserves vendor group sizes).
CI method: basic (reverse-percentile), which corrects for median bias
in near-zero effect sizes:
  CI = [2*theta_hat - q_{97.5}, 2*theta_hat - q_{2.5}]

Outputs:
  - E:/boshi/qm_harmonization_paper/results/step6_mlauc/bootstrap_ci_table2.csv
"""
import warnings
warnings.filterwarnings('ignore')

import numpy as np
import pandas as pd
import statsmodels.api as sm

# -- Config ---------------------------------------------------------
BASE = "E:/boshi/qm_harmonization_paper/results"
BIOM = ["T2w_CSA", "GM_CSA_mm2", "MTR", "MTsat", "FA", "MD", "AD", "RD"]
BATCH = "Manufacturer"
N_BOOT = 1000
SEED = 42

METHODS = [
    ("Original",      "",          None,             None),
    ("ComBat",        "_combat",   "step3_combat",   "biomarkers_combat_HC.csv"),
    ("ComBat-joint",  "_combatJ",  "step3_combat",   "biomarkers_combat_joint_HC.csv"),
    ("CovBat",        "_covbat",   "step3_covbat",   "biomarkers_covbat_HC.csv"),
    ("RELIEF",        "_relief",   "step3_relief",   "biomarkers_relief_HC.csv"),
    ("LME",           "_lme",      "step3_lme",      "biomarkers_lme_HC.csv"),
]

# -- Helper functions ------------------------------------------------

def anova_eta2(y, groups):
    """One-way ANOVA eta^2 = SSB/SST. Handles NaN by available-case."""
    mask = np.isfinite(y)
    y_f, g_f = y[mask], groups[mask]
    if len(y_f) < 4:
        return np.nan
    grand = y_f.mean()
    SST = ((y_f - grand) ** 2).sum()
    if SST <= 0:
        return np.nan
    SSB = 0.0
    for g in np.unique(g_f):
        idx = g_f == g
        n_g = idx.sum()
        if n_g > 0:
            SSB += n_g * (y_f[idx].mean() - grand) ** 2
    return SSB / SST


def permanova_r2_complete(X, groups):
    """PERMANOVA R^2 on complete cases, z-scored ddof=1 (matches original)."""
    # Complete cases only
    mask = np.all(np.isfinite(X), axis=1)
    Xc, gc = X[mask], groups[mask]
    N = len(Xc)
    if N < 10:
        return np.nan
    # z-score ddof=1 (matching step4_eval_multivariate.py L112)
    Z = (Xc - Xc.mean(axis=0)) / (Xc.std(axis=0, ddof=1) + 1e-12)
    # squared Euclidean distance
    diffs = Z[:, None, :] - Z[None, :, :]
    D2 = (diffs ** 2).sum(-1)
    SST = D2.sum() / (2 * N)
    if SST <= 0:
        return np.nan
    SSW = 0.0
    for u in np.unique(gc):
        idx = np.where(gc == u)[0]
        n_u = len(idx)
        if n_u >= 2:
            SSW += D2[np.ix_(idx, idx)].sum() / (2 * n_u)
    return (SST - SSW) / SST


def age_semi_r2(y, age, sex):
    """Semi-partial R^2 for Age: R^2(Age+Sex) - R^2(Sex). Available-case."""
    mask = np.isfinite(y) & np.isfinite(age)
    y_f = y[mask]; age_f = age[mask]; sex_f = sex[mask]
    if mask.sum() < 10:
        return np.nan
    sex_d = pd.get_dummies(pd.Series(sex_f), drop_first=True).astype(float).values
    age_v = age_f.astype(float)
    y_v = y_f.astype(float)
    if sex_d.shape[1] == 0:
        X_full = sm.add_constant(age_v)
        return max(sm.OLS(y_v, X_full).fit().rsquared, 0.0)
    X_full = sm.add_constant(np.column_stack([age_v, sex_d]))
    X_sex = sm.add_constant(sex_d)
    R2_full = sm.OLS(y_v, X_full).fit().rsquared
    R2_sex = sm.OLS(y_v, X_sex).fit().rsquared
    return max(R2_full - R2_sex, 0.0)


def compute_all_metrics(harm, orig, vendor, age, sex, name):
    """Compute all 4 Table-2 metrics on one (possibly resampled) dataset.

    harm: (N, 8) harmonized values (may contain NaN)
    orig: (N, 8) original values (may contain NaN)
    vendor/age/sex: (N,) covariates
    """
    eta2_vals, r_vals, age_r2_vals = [], [], []
    for j in range(harm.shape[1]):
        y = harm[:, j]
        eta2_vals.append(anova_eta2(y, vendor))
        if name == "Original":
            r_vals.append(1.0)
        else:
            o = orig[:, j]
            mask = np.isfinite(y) & np.isfinite(o)
            if mask.sum() >= 3 and y[mask].std() > 0 and o[mask].std() > 0:
                r_vals.append(np.corrcoef(o[mask], y[mask])[0, 1])
            else:
                r_vals.append(np.nan)
        age_r2_vals.append(age_semi_r2(y, age, sex))

    return {
        "eta2": np.nanmean(eta2_vals),
        "perma_r2": permanova_r2_complete(harm, vendor),
        "r_prepost": np.nanmean(r_vals),
        "age_r2": np.nanmean(age_r2_vals),
    }


def basic_ci(point, boot_vals, alpha=0.05):
    """Basic (reverse-percentile) bootstrap CI."""
    lo_q = np.percentile(boot_vals, 100 * (1 - alpha / 2))
    hi_q = np.percentile(boot_vals, 100 * (alpha / 2))
    return 2 * point - lo_q, 2 * point - hi_q


# -- Load data (all HC subjects, N=203) -------------------------------
print("Loading data (all HC subjects, available-case)...")
master = pd.read_csv(f"{BASE}/biomarkers_master.csv")
hc = master[master["Pathology"].fillna("HC").str.upper() == "HC"].reset_index(drop=True)
N = len(hc)
print(f"HC subjects: N={N}")

method_data = {}
for name, suf, subdir, fname in METHODS:
    if name == "Original":
        harm = hc[BIOM].to_numpy(dtype=float)
    else:
        h = pd.read_csv(f"{BASE}/{subdir}/{fname}")
        merged = hc[["participant_id"]].merge(
            h[["participant_id"] + [b + suf for b in BIOM]],
            on="participant_id", how="left"
        )
        harm_cols = [b + suf for b in BIOM]
        harm = merged[harm_cols].to_numpy(dtype=float)

    method_data[name] = {
        "harm": harm,
        "orig": hc[BIOM].to_numpy(dtype=float),
        "vendor": hc[BATCH].to_numpy(),
        "age": hc["Age"].to_numpy(dtype=float),
        "sex": hc["Sex"].to_numpy(),
    }
    nan_count = np.sum(~np.isfinite(harm))
    print(f"  {name:12s}: shape {harm.shape}, NaN cells: {nan_count}")

# -- Point estimates (should match Table 2) ---------------------------
print("\nPoint estimates (available-case for eta2/r/ageR2, complete-case for PERMANOVA):")
point_est = {}
for name in [m[0] for m in METHODS]:
    d = method_data[name]
    point_est[name] = compute_all_metrics(
        d["harm"], d["orig"], d["vendor"], d["age"], d["sex"], name
    )
    pe = point_est[name]
    print(f"  {name:12s}: eta2={pe['eta2']:.4e}  perma={pe['perma_r2']:.4e}  "
          f"r={pe['r_prepost']:.4f}  ageR2={pe['age_r2']:.4e}")

# -- Stratified bootstrap ---------------------------------------------
print(f"\nRunning {N_BOOT} stratified bootstrap iterations (within-vendor)...")

vendor_all = method_data["Original"]["vendor"]
strata = {}
for g in np.unique(vendor_all):
    strata[g] = np.where(vendor_all == g)[0]
print(f"  Strata sizes: { {g: len(idx) for g, idx in strata.items()} }")

rng = np.random.default_rng(SEED)

boot_results = {name: {"eta2": [], "perma_r2": [], "r_prepost": [], "age_r2": []}
                for name in [m[0] for m in METHODS]}

for b in range(N_BOOT):
    # Stratified resample: same vendor group sizes, new subjects within each
    idx = np.concatenate([
        rng.choice(strata[g], size=len(strata[g]), replace=True)
        for g in sorted(strata.keys())
    ])

    for name in [m[0] for m in METHODS]:
        d = method_data[name]
        m_boot = compute_all_metrics(
            d["harm"][idx], d["orig"][idx], d["vendor"][idx],
            d["age"][idx], d["sex"][idx], name
        )
        for metric in ["eta2", "perma_r2", "r_prepost", "age_r2"]:
            boot_results[name][metric].append(m_boot[metric])

    if (b + 1) % 100 == 0:
        print(f"  ...{b+1}/{N_BOOT}")

# -- Summarize with basic CI -----------------------------------------
print("\n\nBootstrap 95% CI results (basic / reverse-percentile):")
print("=" * 130)

rows = []
for name in [m[0] for m in METHODS]:
    pe = point_est[name]
    line = f"{name:12s} "
    for metric in ["eta2", "perma_r2", "r_prepost", "age_r2"]:
        arr = np.array(boot_results[name][metric])
        arr = arr[np.isfinite(arr)]
        pt = pe[metric]

        # Basic CI
        ci_lo_b, ci_hi_b = basic_ci(pt, arr)
        # Percentile CI (for comparison)
        ci_lo_p, ci_hi_p = np.percentile(arr, [2.5, 97.5])

        # Clip to valid ranges
        if metric in ("eta2", "perma_r2", "age_r2"):
            ci_lo_b = max(ci_lo_b, 0.0)
            ci_lo_p = max(ci_lo_p, 0.0)
            ci_hi_b = min(ci_hi_b, 1.0)
            ci_hi_p = min(ci_hi_p, 1.0)
        if metric == "r_prepost":
            ci_lo_b = max(ci_lo_b, -1.0)
            ci_hi_b = min(ci_hi_b, 1.0)

        covered_b = ci_lo_b <= pt <= ci_hi_b

        rows.append({
            "method": name,
            "metric": metric,
            "point": pt,
            "ci_lo_basic": ci_lo_b,
            "ci_hi_basic": ci_hi_b,
            "ci_lo_percentile": ci_lo_p,
            "ci_hi_percentile": ci_hi_p,
            "covered_by_basic": covered_b,
            "n_boot": len(arr),
        })

        fmt = ".4f" if metric == "r_prepost" else ".3e"
        line += f"  {metric}={pt:{fmt}} [{ci_lo_b:{fmt}}, {ci_hi_b:{fmt}}]"
        if not covered_b:
            line += " (!)"

    print(line)

df_out = pd.DataFrame(rows)
out_path = f"{BASE}/step6_mlauc/bootstrap_ci_table2.csv"
df_out.to_csv(out_path, index=False)
print(f"\nSaved: {out_path}")
print(f"Total rows: {len(df_out)}")

n_basic = df_out["covered_by_basic"].sum()
print(f"Coverage: basic CI {n_basic}/{len(df_out)}")
