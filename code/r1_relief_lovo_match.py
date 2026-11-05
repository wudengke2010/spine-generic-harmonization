"""R1-4: RELIEF LOVO fairness check.

The published LOVO implementation applies test-vendor distribution matching for
ComBat (per-biomarker) and CovBat (per-PC), but NOT for RELIEF, whose test vendor
is merely projected onto the training PCA/ICA basis and reconstructed with the
training mu/sd.  That asymmetry could by itself explain why RELIEF "fails" LOVO.

This script re-runs the RELIEF fold evaluating two test-vendor protocols:
  (A) as published  : projection + reconstruction, no matching (reproduces the paper)
  (B) matched       : per-IC z-score matching to the harmonised training scores,
                      exactly as CovBat does in PC space
and reports the pooled-cohort vendor eta^2 (mean over the eight biomarkers) and
the test-vendor pre/post correlation for each.

Output: results/step5_lovo/r1_relief_lovo_match.csv
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.decomposition import FastICA
from neuroCombat import neuroCombat

sys.path.insert(0, str(Path(__file__).resolve().parent))
import lovo_validation as lv  # noqa: E402

BM = lv.BIOMARKERS
BATCH = lv.BATCH_COL
COV = lv.COVARIATES
OUT = Path("E:/boshi/qm_harmonization_paper/results/step5_lovo/r1_relief_lovo_match.csv")


def eta2_mean(X, labels):
    """Mean one-way eta^2 across columns of X (n x p) given vendor labels."""
    vals = []
    X = np.asarray(X, dtype=float)
    labels = np.asarray(labels)
    u = np.unique(labels)
    for j in range(X.shape[1]):
        col = X[:, j]
        groups = [col[labels == v] for v in u]
        groups = [g for g in groups if len(g) > 1]
        allv = col
        ssb = sum(len(g) * (g.mean() - allv.mean()) ** 2 for g in groups)
        sst = ((allv - allv.mean()) ** 2).sum()
        vals.append(ssb / sst)
    return float(np.mean(vals))


def rprepost_mean(pre, post):
    rs = []
    for i in range(pre.shape[1]):
        x, y = pre[:, i], post[:, i]
        if x.std() > 0 and y.std() > 0:
            rs.append(np.corrcoef(x, y)[0, 1])
    return float(np.mean(rs))


def relief_transform(train_df, test_df, match):
    """Copy of lovo_validation.relief_train_transform, with optional per-IC matching."""
    tr = train_df[["participant_id"] + BM + [BATCH] + COV].dropna(subset=BM + [BATCH] + COV).copy()
    te = test_df[["participant_id"] + BM + COV].dropna(subset=BM + COV).copy()

    mu = tr[BM].mean()
    sd = tr[BM].std(ddof=1)
    z_tr = (tr[BM] - mu) / sd

    covars = tr[[BATCH] + COV].reset_index(drop=True)
    out = neuroCombat(dat=z_tr.T.values.astype(float), covars=covars, batch_col=BATCH,
                      categorical_cols=["Sex"], continuous_cols=["Age"], eb=True,
                      parametric=True, mean_only=False)
    Xc = out["data"].T

    D_bio = np.column_stack([np.ones(len(tr)), tr["Age"].values.astype(float),
                             pd.get_dummies(tr["Sex"], drop_first=True).values.astype(float)])
    beta, *_ = np.linalg.lstsq(D_bio, Xc, rcond=None)
    R = Xc - D_bio @ beta

    U, s, Vt = np.linalg.svd(R, full_matrices=False)
    var_ratio = (s ** 2) / (s ** 2).sum()
    K = int(np.searchsorted(np.cumsum(var_ratio), 0.95) + 1)
    K = max(2, min(K, len(s)))
    PC = U[:, :K] * s[:K]

    ica = FastICA(n_components=K, random_state=0, whiten="unit-variance", max_iter=2000, tol=1e-5)
    S = ica.fit_transform(PC)
    A = ica.mixing_

    batches = tr[BATCH].values
    Fs, ps = [], []
    for k in range(K):
        groups = [S[batches == b, k] for b in np.unique(batches)]
        groups = [g for g in groups if len(g) > 1]
        f_, p_ = stats.f_oneway(*groups)
        Fs.append(f_); ps.append(p_)
    p_arr = np.array(ps); m = len(p_arr)
    order = np.argsort(p_arr); ranked = p_arr[order]
    q = ranked * m / (np.arange(m) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    q_full = np.empty_like(q); q_full[order] = np.clip(q, 0, 1)
    flag = q_full < 0.05

    S_clean = S.copy(); S_clean[:, flag] = 0.0
    PC_clean = S_clean @ A.T
    R_h = PC_clean @ Vt[:K, :]
    Z_h = D_bio @ beta + R_h
    train_harm = pd.DataFrame(Z_h * sd.values + mu.values, index=tr.index, columns=BM)

    # ---- test vendor ----
    z_te = (te[BM] - mu) / sd
    D_te = np.column_stack([np.ones(len(te)), te["Age"].values.astype(float),
                            pd.get_dummies(te["Sex"], drop_first=True).values.astype(float)])
    if D_te.shape[1] < D_bio.shape[1]:
        D_te = np.column_stack([D_te, np.zeros((len(te), D_bio.shape[1] - D_te.shape[1]))])
    R_te = z_te.values - D_te @ beta
    PC_te = R_te @ Vt[:K, :].T
    S_te = PC_te @ np.linalg.inv(A.T)
    S_te_clean = S_te.copy(); S_te_clean[:, flag] = 0.0

    if match:
        # mirror CovBat: match each IC of the test vendor to the harmonised training ICs
        S_te_m = np.zeros_like(S_te_clean)
        for k in range(K):
            mu_k, sd_k = S_clean[:, k].mean(), S_clean[:, k].std(ddof=1)
            mu_t, sd_t = S_te_clean[:, k].mean(), S_te_clean[:, k].std(ddof=1)
            if sd_t == 0 or sd_k == 0:
                S_te_m[:, k] = S_te_clean[:, k] - mu_t + mu_k
            else:
                S_te_m[:, k] = (S_te_clean[:, k] - mu_t) / sd_t * sd_k + mu_k
        S_te_clean = S_te_m

    PC_te_c = S_te_clean @ A.T
    R_te_h = PC_te_c @ Vt[:K, :]
    Z_te_h = D_te @ beta + R_te_h
    test_harm = pd.DataFrame(Z_te_h * sd.values + mu.values, index=te.index, columns=BM)

    return train_harm, test_harm, K, int(flag.sum())


def main():
    df = pd.read_csv(lv.MASTER_CSV)
    hc = df[df["Pathology"].fillna("HC").str.upper().eq("HC")].copy()
    hc = hc.dropna(subset=BM + COV + [BATCH]).reset_index(drop=True)
    print(f"HC complete-case N={len(hc)}  vendors={dict(hc[BATCH].value_counts())}")

    rows = []
    for held in lv.VENDORS:
        train_df = hc[hc[BATCH] != held].reset_index(drop=True)
        test_df = hc[hc[BATCH] == held].reset_index(drop=True)

        # baseline (unharmonised) for reference
        base = pd.concat([train_df[BM], test_df[BM]]).reset_index(drop=True)
        labels = pd.concat([train_df[BATCH], test_df[BATCH]]).values

        for match in (False, True):
            tr_h, te_h, K, nflag = relief_transform(train_df, test_df, match=match)
            comb = pd.concat([tr_h, te_h]).reset_index(drop=True)
            pre = base.values
            post = comb.values
            rows.append({
                "held_out": held, "protocol": "matched" if match else "as_published",
                "K": K, "n_flagged": nflag,
                "eta2": eta2_mean(post, labels),
                "r_prepost_test_vendor": rprepost_mean(
                    test_df[BM].values, te_h.loc[test_df.index].values),
                "n_train": len(train_df), "n_test": len(test_df),
            })
            print(f"  {held:8s} {'matched' if match else 'as_published':13s}"
                  f" eta2={rows[-1]['eta2']:.4f}  r_test={rows[-1]['r_prepost_test_vendor']:.3f}")

    out = pd.DataFrame(rows)
    out.to_csv(OUT, index=False)
    print("\n== summary ==")
    print(out.groupby("protocol")[["eta2", "r_prepost_test_vendor"]].agg(["mean", "max"]).round(4))
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()
