"""
Re-compute LOVO PERMANOVA with per-dimension z-scoring (v1.4.0, review-2 fix).

Rationale (reviewer must-fix #2):
  The internal PERMANOVA (step4_eval_multivariate.py) z-scores each biomarker
  within the cohort before computing the Euclidean PERMANOVA, which yields
  R^2 = mean(per-dimension eta^2) (mathematical identity, see main_text P33).
  The original LOVO PERMANOVA (lovo_validation.py line ~701) fed RAW
  (un-standardised) values into the same statistic, so LOVO R^2 was on a
  different scale than the internal R^2 (Original: 0.116 vs 0.315).

This script re-computes the LOVO PERMANOVA on the SAME pooled fold samples
(results/step5_lovo/lovo_combined_{vendor}_{method}.csv), applying the SAME
convention as the internal evaluation:
    Z = (X - mean) / std(ddof=1)  per dimension, over the pooled fold sample
    PERMANOVA: Euclidean, 9,999 permutations, seed 0, grouping = Manufacturer

Outputs:
  results/step5_lovo/lovo_permanova_z.csv   (fold x method -> F, R2, p, n)
  Also prints verification that R2 == mean(per-dimension eta^2).
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path("E:/boshi/qm_harmonization_paper")
LOVO_DIR = ROOT / "results" / "step5_lovo"

BIOMARKERS = ["T2w_CSA", "GM_CSA_mm2", "MTR", "MTsat", "FA", "MD", "AD", "RD"]
BATCH_COL = "Manufacturer"

METHODS = {
    "Original":     "",
    "ComBat":       "_combat",
    "ComBat-joint": "_combatJ",
    "CovBat":       "_covbat",
    "RELIEF":       "_relief",
    "LME":          "_lme",
}

FOLDS = ["GE", "Philips", "Siemens"]


def permanova(X: np.ndarray, groups: np.ndarray, n_perm: int = 9999,
              rng: np.random.Generator | None = None) -> dict:
    """Custom Euclidean PERMANOVA (identical to step4_eval_multivariate.py)."""
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
        F_null[i] = (SST - ssw(g_perm)) / (a - 1) / (ssw(g_perm) / (N - a))
    p = (np.sum(F_null >= F_obs) + 1) / (n_perm + 1)
    return {"F": F_obs, "R2": R2, "p": p, "n": N}


def eta2_1d(x: np.ndarray, g: np.ndarray) -> float:
    """One-way ANOVA eta^2 for a single dimension (scale-invariant)."""
    x = np.asarray(x, float)
    grand = x.mean()
    sst = ((x - grand) ** 2).sum()
    ssb = 0.0
    for u in np.unique(g):
        idx = g == u
        ssb += idx.sum() * (x[idx].mean() - grand) ** 2
    return ssb / sst if sst > 0 else np.nan


def main() -> None:
    rows = []
    ident_rows = []
    for fold in FOLDS:
        for method, suf in METHODS.items():
            f = LOVO_DIR / f"lovo_combined_{fold}_{method}.csv"
            df = pd.read_csv(f)
            cols = [b + suf for b in BIOMARKERS]
            df = df.dropna(subset=cols + [BATCH_COL]).reset_index(drop=True)
            X = df[cols].values.astype(float)
            g = df[BATCH_COL].values
            # per-dimension z-scoring over the pooled fold sample (internal convention)
            Z = (X - X.mean(axis=0)) / X.std(axis=0, ddof=1)
            res = permanova(Z, g, n_perm=9999, rng=np.random.default_rng(0))
            # identity check: R2 == mean per-dimension eta^2
            etas = [eta2_1d(Z[:, j], g) for j in range(Z.shape[1])]
            ident_rows.append({
                "fold": fold, "method": method,
                "permanova_R2": res["R2"],
                "mean_eta2_dim": float(np.mean(etas)),
                "abs_diff": abs(res["R2"] - np.mean(etas)),
            })
            rows.append({"fold": fold, "method": method,
                         "permanova_F": res["F"],
                         "permanova_R2": res["R2"],
                         "permanova_p": res["p"], "n": res["n"]})
            print(f"{fold:8s} {method:12s} F={res['F']:8.3f}  R2={res['R2']:.6g}  "
                  f"p={res['p']:.4f}  mean_eta2={np.mean(etas):.6g}")

    out = pd.DataFrame(rows)
    out.to_csv(LOVO_DIR / "lovo_permanova_z.csv", index=False)
    ident = pd.DataFrame(ident_rows)
    ident.to_csv(LOVO_DIR / "lovo_permanova_z_identity_check.csv", index=False)

    print("\n=== identity check: max |R2 - mean(eta2_dim)| ===")
    print(ident["abs_diff"].max())

    print("\n=== 3-fold mean R2 per method (for Fig 4B) ===")
    print(out.groupby("method")["permanova_R2"].mean().to_string())
    print("\n=== 3-fold min/max p per method ===")
    print(out.groupby("method")["permanova_p"].agg(["min", "max"]).to_string())


if __name__ == "__main__":
    main()
