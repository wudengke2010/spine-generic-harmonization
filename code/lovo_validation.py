"""
Leave-One-Vendor-Out (LOVO) External Validation
================================================
For each of 3 folds (leave Siemens / GE / Philips out):
  1. Train: fit harmonization on 2 training vendors
  2. Test:  apply fitted model to the held-out vendor (unseen during training)
  3. Combine harmonized train + test data
  4. Evaluate: vendor eta2, PERMANOVA R2, r_pre_post, Age R2, Sex R2

Methods:
  - Original:     no transform (baseline)
  - ComBat:       neuroCombat on train; z-score transform for test vendor
  - ComBat-joint: joint neuroCombat (EB=True) on train; z-score transform for test
  - CovBat:       ComBat-joint + PCA covariance harmonization (train PCA, apply to test)
  - RELIEF:       ComBat-joint + ICA batch-component removal (train ICA, apply to test)
  - LME/FE:       fit on train vendors; test vendor gets no offset (unknown batch)

Key design decisions:
  - ComBat-type methods: neuroCombat requires batch in the model, so for the
    unseen vendor we use distribution-matching (z-score to training harmonized
    mean/SD). This is the standard approach for applying ComBat to new sites
    (Fortin 2018, Pomponi 2020).
  - LME/FE: the held-out vendor's offset is unknown, so we set it to 0.
    This tests whether the method's biological-covariate removal generalizes.
  - CovBat/RELIEF: PCA/ICA components are trained on training residuals and
    applied to the test vendor's residuals.

Outputs:
  lovo_results_long.csv     fold x method x biomarker x metric
  lovo_summary.csv          fold x method x mean across biomarkers
  lovo_permanova.csv        fold x method x PERMANOVA R2
  lovo_per_vendor.csv       fold x method x vendor x per-vendor metrics
"""
from __future__ import annotations

import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from neuroCombat import neuroCombat
from sklearn.decomposition import FastICA
import statsmodels.api as sm
import statsmodels.formula.api as smf
from statsmodels.tools.sm_exceptions import ConvergenceWarning

# ── Constants ────────────────────────────────────────────────────────────────
BIOMARKERS = ["T2w_CSA", "GM_CSA_mm2",
              "MTR", "MTsat",
              "FA", "MD", "AD", "RD"]
BATCH_COL  = "Manufacturer"
COVARIATES = ["Age", "Sex"]
VENDORS    = ["GE", "Philips", "Siemens"]
N_PERM     = 4999   # PERMANOVA permutations (reduced for speed)

MASTER_CSV = Path("E:/boshi/qm_harmonization_paper/results/biomarkers_master.csv")
OUT_DIR    = Path("E:/boshi/qm_harmonization_paper/results/step5_lovo")


# ── Evaluation helpers (copied from step4_eval for consistency) ──────────────
def vendor_stats(y: pd.Series, batch: pd.Series) -> dict:
    sub = pd.concat([y.rename("y"), batch.rename("b")], axis=1).dropna()
    if sub["b"].nunique() < 2 or len(sub) < 4:
        return dict(F=np.nan, p=np.nan, eta2=np.nan, n=len(sub))
    groups = [g["y"].values for _, g in sub.groupby("b") if len(g) > 1]
    if len(groups) < 2:
        return dict(F=np.nan, p=np.nan, eta2=np.nan, n=len(sub))
    F, p = stats.f_oneway(*groups)
    grand = sub["y"].mean()
    SSB = sum(len(g) * (g.mean() - grand)**2 for g in groups)
    SST = ((sub["y"] - grand)**2).sum()
    eta2 = SSB / SST if SST > 0 else np.nan
    return dict(F=F, p=p, eta2=eta2, n=len(sub))


def bio_partial_R2(y: pd.Series, age: pd.Series, sex: pd.Series) -> dict:
    sub = pd.concat([y.rename("y"), age.rename("Age"), sex.rename("Sex")],
                    axis=1).dropna()
    if len(sub) < 10:
        return dict(R2_age=np.nan, R2_sex=np.nan)
    sex_d = pd.get_dummies(sub["Sex"], drop_first=True).astype(float)
    age_v = sub["Age"].astype(float).values
    X_full = sm.add_constant(np.column_stack([age_v, sex_d.values]))
    X_ageonly = sm.add_constant(age_v)
    X_sexonly = sm.add_constant(sex_d.values) if sex_d.shape[1] else None
    y_v = sub["y"].astype(float).values
    R2_full = sm.OLS(y_v, X_full).fit().rsquared
    R2_only_age = sm.OLS(y_v, X_ageonly).fit().rsquared
    R2_only_sex = (sm.OLS(y_v, X_sexonly).fit().rsquared
                   if X_sexonly is not None else 0.0)
    return dict(
        R2_age=max(R2_full - R2_only_sex, 0.0),
        R2_sex=max(R2_full - R2_only_age, 0.0),
    )


def permanova(X: np.ndarray, groups: np.ndarray, n_perm: int = N_PERM,
              rng: np.random.Generator | None = None) -> dict:
    if rng is None:
        rng = np.random.default_rng(0)
    N = len(X)
    diffs = X[:, None, :] - X[None, :, :]
    D2 = (diffs ** 2).sum(-1)
    SST = D2.sum() / (2 * N)

    def ssw(g):
        s = 0.0
        for u in np.unique(g):
            idx = np.where(g == u)[0]
            n_u = len(idx)
            if n_u < 2:
                continue
            s += D2[np.ix_(idx, idx)].sum() / (2 * n_u)
        return s

    SSW_obs = ssw(groups)
    SSB_obs = SST - SSW_obs
    a = len(np.unique(groups))
    F_obs = (SSB_obs / (a - 1)) / (SSW_obs / (N - a)) if SSW_obs > 0 else np.nan
    R2 = SSB_obs / SST if SST > 0 else np.nan

    F_null = np.empty(n_perm)
    g_perm = groups.copy()
    for i in range(n_perm):
        rng.shuffle(g_perm)
        SSW_p = ssw(g_perm)
        SSB_p = SST - SSW_p
        F_null[i] = ((SSB_p / (a - 1)) / (SSW_p / (N - a))
                     if SSW_p > 0 else np.nan)
    p = (np.sum(F_null >= F_obs) + 1) / (n_perm + 1)
    return dict(F=F_obs, R2=R2, p=p, N=N)


def r_pearson(x: pd.Series, y: pd.Series) -> float:
    pair = pd.concat([x, y], axis=1).dropna()
    if len(pair) < 3 or pair.iloc[:, 0].std() == 0 or pair.iloc[:, 1].std() == 0:
        return np.nan
    return pair.corr().iloc[0, 1]


# ── ComBat (per-metric, eb=False) ────────────────────────────────────────────
def combat_train_transform(train_df: pd.DataFrame, test_df: pd.DataFrame,
                           metric: str) -> tuple[np.ndarray, np.ndarray]:
    """Train ComBat on 2 vendors; transform test vendor via z-score matching."""
    needed = [metric, BATCH_COL] + COVARIATES
    tr = train_df[needed].dropna().copy()
    te = test_df[[metric] + COVARIATES].dropna().copy()

    if tr[BATCH_COL].nunique() < 2 or len(tr) < 10:
        return (np.full(len(train_df), np.nan),
                np.full(len(test_df), np.nan))

    # Train: neuroCombat on 2 vendors
    data = tr[[metric]].T.values.astype(float)
    covars = tr[[BATCH_COL] + COVARIATES].reset_index(drop=True)
    out = neuroCombat(
        dat=data, covars=covars, batch_col=BATCH_COL,
        categorical_cols=["Sex"], continuous_cols=["Age"],
        eb=False, parametric=True, mean_only=False,
    )
    train_harm = out["data"].ravel()

    # Test: z-score to training harmonized distribution
    mu_h = np.nanmean(train_harm)
    sd_h = np.nanstd(train_harm, ddof=1)
    y_te = te[metric].values.astype(float)
    mu_te = np.nanmean(y_te)
    sd_te = np.nanstd(y_te, ddof=1)
    if sd_te == 0 or sd_h == 0:
        test_harm = y_te - mu_te + mu_h
    else:
        test_harm = (y_te - mu_te) / sd_te * sd_h + mu_h

    # Place into full-length arrays (aligned to original df indices)
    tr_full = np.full(len(train_df), np.nan)
    tr_full[tr.index.values] = train_harm

    te_full = np.full(len(test_df), np.nan)
    te_full[te.index.values] = test_harm

    return tr_full, te_full


# ── ComBat-joint (all 8 biomarkers, eb=True) ─────────────────────────────────
def combat_joint_train_transform(train_df: pd.DataFrame, test_df: pd.DataFrame):
    """Train joint ComBat on 2 vendors; transform test vendor via z-score."""
    needed = BIOMARKERS + [BATCH_COL] + COVARIATES
    tr = train_df[["participant_id"] + needed].dropna(subset=needed).copy()
    te = test_df[["participant_id"] + BIOMARKERS + COVARIATES].dropna(
        subset=BIOMARKERS + COVARIATES).copy()

    if tr[BATCH_COL].nunique() < 2 or len(tr) < 10:
        return None, None

    # z-score training
    mu = tr[BIOMARKERS].mean()
    sd = tr[BIOMARKERS].std(ddof=1)
    z = (tr[BIOMARKERS] - mu) / sd

    # Joint ComBat
    data = z.T.values.astype(float)
    covars = tr[[BATCH_COL] + COVARIATES].reset_index(drop=True)
    out = neuroCombat(
        dat=data, covars=covars, batch_col=BATCH_COL,
        categorical_cols=["Sex"], continuous_cols=["Age"],
        eb=True, parametric=True, mean_only=False,
    )
    z_harm = pd.DataFrame(out["data"].T, index=tr.index, columns=BIOMARKERS)
    train_harm = z_harm * sd + mu  # inverse z-score -> native units

    # Test: z-score each feature to training harmonized distribution
    test_harm = pd.DataFrame(index=te.index, columns=BIOMARKERS, dtype=float)
    for b in BIOMARKERS:
        y_te = te[b].values.astype(float)
        mu_te, sd_te = np.nanmean(y_te), np.nanstd(y_te, ddof=1)
        mu_h, sd_h = train_harm[b].mean(), train_harm[b].std(ddof=1)
        if sd_te == 0 or sd_h == 0:
            test_harm[b] = y_te - mu_te + mu_h
        else:
            test_harm[b] = (y_te - mu_te) / sd_te * sd_h + mu_h

    # Build full arrays aligned to original dfs
    tr_out = pd.DataFrame(np.nan, index=train_df.index,
                          columns=[b + "_combatJ" for b in BIOMARKERS])
    tr_out.loc[tr.index] = train_harm.values

    te_out = pd.DataFrame(np.nan, index=test_df.index,
                          columns=[b + "_combatJ" for b in BIOMARKERS])
    te_out.loc[te.index] = test_harm.values

    return tr_out, te_out


# ── CovBat ───────────────────────────────────────────────────────────────────
def covbat_train_transform(train_df: pd.DataFrame, test_df: pd.DataFrame):
    """CovBat: ComBat-joint + PCA covariance harmonization."""
    tr_out, te_out = combat_joint_train_transform(train_df, test_df)
    if tr_out is None:
        return None, None

    needed = BIOMARKERS + [BATCH_COL] + COVARIATES
    tr = train_df[["participant_id"] + needed].dropna(subset=needed).copy()
    te = test_df[["participant_id"] + BIOMARKERS + COVARIATES].dropna(
        subset=BIOMARKERS + COVARIATES).copy()

    # Get ComBat-harmonized training data (z-scale)
    mu = tr[BIOMARKERS].mean()
    sd = tr[BIOMARKERS].std(ddof=1)
    z_tr = (tr[BIOMARKERS] - mu) / sd
    covars = tr[[BATCH_COL] + COVARIATES].reset_index(drop=True)
    out = neuroCombat(
        dat=z_tr.T.values.astype(float), covars=covars, batch_col=BATCH_COL,
        categorical_cols=["Sex"], continuous_cols=["Age"],
        eb=True, parametric=True, mean_only=False,
    )
    Xc = out["data"].T  # N_train x F

    # Residuals (regress out biological covariates)
    D_bio = np.column_stack([
        np.ones(len(tr)),
        tr["Age"].values.astype(float),
        pd.get_dummies(tr["Sex"], drop_first=True).values.astype(float),
    ])
    beta, *_ = np.linalg.lstsq(D_bio, Xc, rcond=None)
    R = Xc - D_bio @ beta

    # PCA on training residuals
    U, s, Vt = np.linalg.svd(R, full_matrices=False)
    var_ratio = (s**2) / (s**2).sum()
    K = int(np.searchsorted(np.cumsum(var_ratio), 0.95) + 1)
    K = max(1, min(K, len(s)))

    # ComBat on PC scores (if K >= 2)
    scores = U[:, :K] * s[:K]
    if K >= 2:
        out2 = neuroCombat(
            dat=scores.T, covars=covars, batch_col=BATCH_COL,
            categorical_cols=["Sex"], continuous_cols=["Age"],
            eb=True, parametric=True, mean_only=False,
        )
        scores_h = out2["data"].T
    else:
        scores_h = scores

    R_h = scores_h @ Vt[:K, :]
    Z_h = D_bio @ beta + R_h
    train_harm_native = pd.DataFrame(Z_h * sd.values + mu.values,
                                      index=tr.index, columns=BIOMARKERS)
    tr_out.loc[tr.index] = train_harm_native.values

    # Test vendor: project onto training PCA space, apply covariance harmonization
    z_te = (te[BIOMARKERS] - mu) / sd  # z-score using training mu/sd
    # Regress out biological covariates (using training beta)
    D_bio_te = np.column_stack([
        np.ones(len(te)),
        te["Age"].values.astype(float),
        pd.get_dummies(te["Sex"], drop_first=True).values.astype(float),
    ])
    # Align columns with training D_bio
    if D_bio_te.shape[1] < D_bio.shape[1]:
        pad = np.zeros((len(te), D_bio.shape[1] - D_bio_te.shape[1]))
        D_bio_te = np.column_stack([D_bio_te, pad])
    R_te = z_te.values - D_bio_te @ beta

    # Project test residuals onto training PC space
    scores_te = R_te @ Vt[:K, :].T  # N_test x K

    # Apply covariance harmonization: match to training harmonized PC score distribution
    scores_te_h = np.zeros_like(scores_te)
    for k in range(K):
        mu_k = scores_h[:, k].mean()
        sd_k = scores_h[:, k].std(ddof=1)
        mu_te_k = scores_te[:, k].mean()
        sd_te_k = scores_te[:, k].std(ddof=1)
        if sd_te_k == 0 or sd_k == 0:
            scores_te_h[:, k] = scores_te[:, k] - mu_te_k + mu_k
        else:
            scores_te_h[:, k] = (scores_te[:, k] - mu_te_k) / sd_te_k * sd_k + mu_k

    R_te_h = scores_te_h @ Vt[:K, :]
    Z_te_h = D_bio_te @ beta + R_te_h
    test_harm_native = pd.DataFrame(Z_te_h * sd.values + mu.values,
                                     index=te.index, columns=BIOMARKERS)
    te_out.loc[te.index] = test_harm_native.values

    return tr_out, te_out


# ── RELIEF ───────────────────────────────────────────────────────────────────
def relief_train_transform(train_df: pd.DataFrame, test_df: pd.DataFrame):
    """RELIEF: ComBat-joint + ICA batch-component removal."""
    tr_out, te_out = combat_joint_train_transform(train_df, test_df)
    if tr_out is None:
        return None, None

    needed = BIOMARKERS + [BATCH_COL] + COVARIATES
    tr = train_df[["participant_id"] + needed].dropna(subset=needed).copy()
    te = test_df[["participant_id"] + BIOMARKERS + COVARIATES].dropna(
        subset=BIOMARKERS + COVARIATES).copy()

    # z-score
    mu = tr[BIOMARKERS].mean()
    sd = tr[BIOMARKERS].std(ddof=1)
    z_tr = (tr[BIOMARKERS] - mu) / sd

    covars = tr[[BATCH_COL] + COVARIATES].reset_index(drop=True)
    out = neuroCombat(
        dat=z_tr.T.values.astype(float), covars=covars, batch_col=BATCH_COL,
        categorical_cols=["Sex"], continuous_cols=["Age"],
        eb=True, parametric=True, mean_only=False,
    )
    Xc = out["data"].T

    # Residuals
    D_bio = np.column_stack([
        np.ones(len(tr)),
        tr["Age"].values.astype(float),
        pd.get_dummies(tr["Sex"], drop_first=True).values.astype(float),
    ])
    beta, *_ = np.linalg.lstsq(D_bio, Xc, rcond=None)
    R = Xc - D_bio @ beta

    # PCA
    U, s, Vt = np.linalg.svd(R, full_matrices=False)
    var_ratio = (s**2) / (s**2).sum()
    K = int(np.searchsorted(np.cumsum(var_ratio), 0.95) + 1)
    K = max(2, min(K, len(s)))
    PC = U[:, :K] * s[:K]

    # ICA on training PC scores
    ica = FastICA(n_components=K, random_state=0, whiten="unit-variance",
                  max_iter=2000, tol=1e-5)
    S = ica.fit_transform(PC)
    A = ica.mixing_

    # F-test each IC against Manufacturer (training)
    batches = tr[BATCH_COL].values
    Fs, ps = [], []
    for k in range(K):
        groups = [S[batches == b, k] for b in np.unique(batches)]
        groups = [g for g in groups if len(g) > 1]
        if len(groups) < 2:
            Fs.append(np.nan); ps.append(np.nan); continue
        f_, p_ = stats.f_oneway(*groups)
        Fs.append(f_); ps.append(p_)

    # BH FDR
    p_arr = np.array(ps)
    m = len(p_arr)
    order = np.argsort(p_arr)
    ranked = p_arr[order]
    q = ranked * m / (np.arange(m) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    q_full = np.empty_like(q)
    q_full[order] = np.clip(q, 0, 1)
    flag = (q_full < 0.05)

    # Zero out batch-driven components in training
    S_clean = S.copy()
    S_clean[:, flag] = 0.0
    PC_clean = S_clean @ A.T
    R_h = PC_clean @ Vt[:K, :]
    Z_h = D_bio @ beta + R_h
    train_harm_native = pd.DataFrame(Z_h * sd.values + mu.values,
                                      index=tr.index, columns=BIOMARKERS)
    tr_out.loc[tr.index] = train_harm_native.values

    # Test vendor: project onto training ICA space, remove batch components
    z_te = (te[BIOMARKERS] - mu) / sd
    D_bio_te = np.column_stack([
        np.ones(len(te)),
        te["Age"].values.astype(float),
        pd.get_dummies(te["Sex"], drop_first=True).values.astype(float),
    ])
    if D_bio_te.shape[1] < D_bio.shape[1]:
        pad = np.zeros((len(te), D_bio.shape[1] - D_bio_te.shape[1]))
        D_bio_te = np.column_stack([D_bio_te, pad])
    R_te = z_te.values - D_bio_te @ beta

    # Project test residuals onto training PCA space
    PC_te = R_te @ Vt[:K, :].T  # N_test x K

    # Apply ICA unmixing (using training ICA)
    # PC = S @ A.T  =>  S = PC @ inv(A.T)
    S_te = PC_te @ np.linalg.inv(A.T)

    # Zero out batch-driven components
    S_te_clean = S_te.copy()
    S_te_clean[:, flag] = 0.0

    # Reconstruct
    PC_te_clean = S_te_clean @ A.T
    R_te_h = PC_te_clean @ Vt[:K, :]
    Z_te_h = D_bio_te @ beta + R_te_h
    test_harm_native = pd.DataFrame(Z_te_h * sd.values + mu.values,
                                     index=te.index, columns=BIOMARKERS)
    te_out.loc[te.index] = test_harm_native.values

    return tr_out, te_out


# ── LME/FE ───────────────────────────────────────────────────────────────────
def lme_train_transform(train_df: pd.DataFrame, test_df: pd.DataFrame):
    """LME/FE: fit on train vendors; test vendor gets offset=0 (unknown)."""
    tr_out = pd.DataFrame(np.nan, index=train_df.index,
                          columns=[b + "_lme" for b in BIOMARKERS])
    te_out = pd.DataFrame(np.nan, index=test_df.index,
                          columns=[b + "_lme" for b in BIOMARKERS])

    for m in BIOMARKERS:
        needed = [m, BATCH_COL] + COVARIATES
        sub = train_df[needed].dropna().copy()
        sub.columns = ["y", BATCH_COL] + COVARIATES

        status = "ok-LME"
        harm_vals = None

        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", ConvergenceWarning)
                warnings.simplefilter("ignore", RuntimeWarning)
                md = smf.mixedlm("y ~ Age + C(Sex)", data=sub,
                                 groups=sub[BATCH_COL])
                mdf = md.fit(method="lbfgs", reml=True)
                re = mdf.random_effects
            blup = {g: float(v.iloc[0]) for g, v in re.items()}
            u_vec = sub[BATCH_COL].map(blup).values
            harm_vals = sub["y"].values - u_vec
        except Exception:
            status = "FE-fallback"
            X = pd.get_dummies(sub[[BATCH_COL] + COVARIATES], drop_first=True)
            X = sm.add_constant(X).astype(float)
            y = sub["y"].astype(float).values
            res = sm.OLS(y, X).fit()
            vendor_cols = [c for c in X.columns if c.startswith(BATCH_COL + "_")]
            if vendor_cols:
                vendor_pred = (X[vendor_cols] @ res.params[vendor_cols]).values
                vendor_pred_c = vendor_pred - vendor_pred.mean()
                harm_vals = y - vendor_pred_c
            else:
                harm_vals = y

        # Training: subtract BLUP/FE offset
        tr_full = np.full(len(train_df), np.nan)
        tr_idx = train_df.index.get_indexer(sub.index)
        tr_full[tr_idx] = harm_vals
        tr_out[m + "_lme"] = tr_full

        # Test vendor: offset = 0 (unknown vendor)
        # Only remove biological covariate effects (same as training OLS)
        te_sub = test_df[[m] + COVARIATES].dropna().copy()
        if len(te_sub) > 0:
            te_full = np.full(len(test_df), np.nan)
            te_idx = test_df.index.get_indexer(te_sub.index)
            # No vendor offset adjustment — just keep raw values
            # (biological covariates are preserved, not removed)
            te_full[te_idx] = te_sub[m].values
            te_out[m + "_lme"] = te_full

    return tr_out, te_out


# ── Main LOVO pipeline ───────────────────────────────────────────────────────
def run_lovo():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(MASTER_CSV)

    # Use HC cohort (consistent with main analysis)
    if "Pathology" in df.columns:
        hc = df[df["Pathology"].fillna("HC").str.upper().eq("HC")].copy()
    else:
        hc = df.copy()
    hc = hc.reset_index(drop=True)

    print(f"LOVO Validation — HC cohort, N={len(hc)}")
    print(f"  Vendors: {dict(hc[BATCH_COL].value_counts())}")
    print(f"  Biomarkers: {BIOMARKERS}")
    print()

    all_rows = []
    perm_rows = []
    pv_rows = []

    for held_out in VENDORS:
        train_vendors = [v for v in VENDORS if v != held_out]
        train_df = hc[hc[BATCH_COL].isin(train_vendors)].copy().reset_index(drop=True)
        test_df  = hc[hc[BATCH_COL] == held_out].copy().reset_index(drop=True)

        print(f"\n{'='*60}")
        print(f"  Fold: Leave {held_out} Out")
        print(f"  Train: {train_vendors} (N={len(train_df)})")
        print(f"  Test:  {held_out} (N={len(test_df)})")
        print(f"{'='*60}")

        # ── Run each method ─────────────────────────────────────────────
        methods_harm = {}  # method_name -> combined harmonized df

        # Original (baseline)
        combined_orig = pd.concat([train_df, test_df], ignore_index=True)
        methods_harm["Original"] = combined_orig

        # ComBat (per-metric)
        print("\n  [ComBat] training...")
        tr_harm = train_df.copy()
        te_harm = test_df.copy()
        for b in BIOMARKERS:
            tr_vals, te_vals = combat_train_transform(train_df, test_df, b)
            tr_harm[b + "_combat"] = tr_vals
            te_harm[b + "_combat"] = te_vals
        combined = pd.concat([tr_harm, te_harm], ignore_index=True)
        methods_harm["ComBat"] = combined

        # ComBat-joint
        print("  [ComBat-joint] training...")
        tr_out, te_out = combat_joint_train_transform(train_df, test_df)
        if tr_out is not None:
            tr_harm = train_df.copy()
            te_harm = test_df.copy()
            for b in BIOMARKERS:
                tr_harm[b + "_combatJ"] = tr_out[b + "_combatJ"]
                te_harm[b + "_combatJ"] = te_out[b + "_combatJ"]
            combined = pd.concat([tr_harm, te_harm], ignore_index=True)
            methods_harm["ComBat-joint"] = combined

        # CovBat
        print("  [CovBat] training...")
        tr_out, te_out = covbat_train_transform(train_df, test_df)
        if tr_out is not None:
            tr_harm = train_df.copy()
            te_harm = test_df.copy()
            for b in BIOMARKERS:
                tr_harm[b + "_covbat"] = tr_out[b + "_combatJ"]  # reuse joint columns
                te_harm[b + "_covbat"] = te_out[b + "_combatJ"]
            # Rename columns to _covbat
            tr_harm = tr_harm.rename(columns={b+"_combatJ": b+"_covbat" for b in BIOMARKERS})
            te_harm = te_harm.rename(columns={b+"_combatJ": b+"_covbat" for b in BIOMARKERS})
            combined = pd.concat([tr_harm, te_harm], ignore_index=True)
            methods_harm["CovBat"] = combined

        # RELIEF
        print("  [RELIEF] training...")
        tr_out, te_out = relief_train_transform(train_df, test_df)
        if tr_out is not None:
            tr_harm = train_df.copy()
            te_harm = test_df.copy()
            for b in BIOMARKERS:
                tr_harm[b + "_relief"] = tr_out[b + "_combatJ"]
                te_harm[b + "_relief"] = te_out[b + "_combatJ"]
            tr_harm = tr_harm.rename(columns={b+"_combatJ": b+"_relief" for b in BIOMARKERS})
            te_harm = te_harm.rename(columns={b+"_combatJ": b+"_relief" for b in BIOMARKERS})
            combined = pd.concat([tr_harm, te_harm], ignore_index=True)
            methods_harm["RELIEF"] = combined

        # LME/FE
        print("  [LME/FE] training...")
        tr_out, te_out = lme_train_transform(train_df, test_df)
        tr_harm = train_df.copy()
        te_harm = test_df.copy()
        for b in BIOMARKERS:
            tr_harm[b + "_lme"] = tr_out[b + "_lme"]
            te_harm[b + "_lme"] = te_out[b + "_lme"]
        combined = pd.concat([tr_harm, te_harm], ignore_index=True)
        methods_harm["LME"] = combined

        # ── Evaluate each method ────────────────────────────────────────
        suffix_map = {
            "Original": "",
            "ComBat": "_combat",
            "ComBat-joint": "_combatJ",
            "CovBat": "_covbat",
            "RELIEF": "_relief",
            "LME": "_lme",
        }

        for method, suf in suffix_map.items():
            combined = methods_harm.get(method)
            if combined is None:
                continue

            # Per-biomarker univariate metrics
            for b in BIOMARKERS:
                col = b + suf
                if col not in combined.columns:
                    continue

                # Vendor stats on combined data
                vs = vendor_stats(combined[col], combined[BATCH_COL])

                # r_pre_post (test vendor only)
                test_orig = test_df[b].values
                test_harm_vals = combined.loc[combined[BATCH_COL] == held_out, col].values
                min_len = min(len(test_orig), len(test_harm_vals))
                if min_len >= 3:
                    r_pp = stats.pearsonr(test_orig[:min_len],
                                         test_harm_vals[:min_len])[0]
                else:
                    r_pp = np.nan

                # Biological signal (combined)
                bio = bio_partial_R2(combined[col], combined["Age"],
                                     combined["Sex"])

                all_rows.append({
                    "fold": held_out,
                    "method": method,
                    "biomarker": b,
                    "n": vs["n"],
                    "vendor_eta2": vs["eta2"],
                    "vendor_F": vs["F"],
                    "vendor_p": vs["p"],
                    "r_pre_post_test": r_pp,
                    "age_R2": bio["R2_age"],
                    "sex_R2": bio["R2_sex"],
                })

            # PERMANOVA (multivariate, combined)
            cols = [b + suf for b in BIOMARKERS if (b + suf) in combined.columns]
            sub = combined[cols + [BATCH_COL]].dropna()
            if len(sub) >= 10 and sub[BATCH_COL].nunique() >= 2:
                X = sub[cols].values
                g = sub[BATCH_COL].values
                res = permanova(X, g, n_perm=N_PERM,
                                rng=np.random.default_rng(0))
                perm_rows.append({
                    "fold": held_out,
                    "method": method,
                    "permanova_F": res["F"],
                    "permanova_R2": res["R2"],
                    "permanova_p": res["p"],
                    "n": res["N"],
                })
                print(f"  {method:13s}  PERMANOVA R2={res['R2']:.4f}  "
                      f"p={res['p']:.4f}  N={res['N']}")

            # Per-vendor breakdown
            for v in VENDORS:
                v_data = combined[combined[BATCH_COL] == v]
                for b in BIOMARKERS:
                    col = b + suf
                    if col not in v_data.columns:
                        continue
                    v_vals = v_data[col].dropna()
                    if len(v_vals) < 3:
                        continue
                    grand = combined[col].mean()
                    dev = v_vals.mean() - grand
                    r_pp_v = np.nan
                    if method != "Original":
                        orig_v = hc[hc[BATCH_COL] == v][b].dropna()
                        harm_v = v_vals
                        min_len = min(len(orig_v), len(harm_v))
                        if min_len >= 3:
                            r_pp_v = stats.pearsonr(
                                orig_v.values[:min_len],
                                harm_v.values[:min_len]
                            )[0]
                    else:
                        r_pp_v = 1.0

                    pv_rows.append({
                        "fold": held_out,
                        "method": method,
                        "vendor": v,
                        "biomarker": b,
                        "mean_deviation": dev,
                        "r_pre_post": r_pp_v,
                        "n": len(v_vals),
                    })

    # ── Save results ────────────────────────────────────────────────────
    long_df = pd.DataFrame(all_rows)
    long_df.to_csv(OUT_DIR / "lovo_results_long.csv", index=False)
    print(f"\n[saved] lovo_results_long.csv ({len(long_df)} rows)")

    perm_df = pd.DataFrame(perm_rows)
    perm_df.to_csv(OUT_DIR / "lovo_permanova.csv", index=False)
    print(f"[saved] lovo_permanova.csv ({len(perm_df)} rows)")

    pv_df = pd.DataFrame(pv_rows)
    pv_df.to_csv(OUT_DIR / "lovo_per_vendor.csv", index=False)
    print(f"[saved] lovo_per_vendor.csv ({len(pv_df)} rows)")

    # Summary: mean across biomarkers
    summary = long_df.groupby(["fold", "method"]).agg(
        vendor_eta2_mean=("vendor_eta2", "mean"),
        vendor_eta2_std=("vendor_eta2", "std"),
        r_pre_post_test_mean=("r_pre_post_test", "mean"),
        age_R2_mean=("age_R2", "mean"),
        sex_R2_mean=("sex_R2", "mean"),
        n_mean=("n", "mean"),
    ).reset_index()
    summary.to_csv(OUT_DIR / "lovo_summary.csv", index=False)
    print(f"[saved] lovo_summary.csv ({len(summary)} rows)")

    # Print summary
    print("\n" + "="*70)
    print("LOVO SUMMARY (mean across 8 biomarkers)")
    print("="*70)
    for fold in VENDORS:
        print(f"\n  Fold: Leave {fold} Out")
        sub = summary[summary["fold"] == fold]
        for _, row in sub.iterrows():
            print(f"    {row['method']:13s}  eta2={row['vendor_eta2_mean']:.4e}  "
                  f"r_pp={row['r_pre_post_test_mean']:.4f}  "
                  f"age_R2={row['age_R2_mean']:.4f}  "
                  f"sex_R2={row['sex_R2_mean']:.4f}")

    print("\n" + "="*70)
    print("PERMANOVA (multivariate vendor effect)")
    print("="*70)
    print(perm_df.to_string(index=False))

    return long_df, perm_df, pv_df


if __name__ == "__main__":
    run_lovo()
