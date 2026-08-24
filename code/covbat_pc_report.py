"""
CovBat PC retention report (Chen et al. 2022 recommendation)
============================================================
Faithfully replicates the PCA stage of code/step3_covbat.py up to (and
including) the K-selection step, then writes a per-PC variance table for
the HC and ALL cohorts.

Steps replicated (identical to step3_covbat.py):
  1. complete-case subset (biomarkers + Manufacturer + Age + Sex)
  2. per-feature z-score (ddof=1)
  3. ComBat (eb=True, parametric, mean_only=False) with Age+Sex covariates
  4. residuals R = Xc - D_bio @ beta_hat (intercept + Age + Sex only)
  5. SVD of R -> variance ratio, cumulative variance, K at 95%

Outputs:
  results/step3_covbat/covbat_pc_retention.csv   (summary per cohort)
  results/step3_covbat/covbat_pc_variance.csv    (per-PC variance, all PCs)
"""
from pathlib import Path

import numpy as np
import pandas as pd
from neuroCombat import neuroCombat

BIOMARKERS = ["T2w_CSA", "GM_CSA_mm2",
              "MTR", "MTsat",
              "FA", "MD", "AD", "RD"]
BATCH_COL  = "Manufacturer"
COVARIATES = ["Age", "Sex"]
PC_VAR_KEEP = 0.95

MASTER = Path("E:/boshi/qm_harmonization_paper/results/biomarkers_master.csv")
OUT_DIR = Path("E:/boshi/qm_harmonization_paper/results/step3_covbat")


def pc_retention(df: pd.DataFrame, label: str) -> tuple[dict, pd.DataFrame]:
    needed = BIOMARKERS + [BATCH_COL] + COVARIATES
    sub = df[needed].dropna().copy()
    n = len(sub)

    # z-score per feature
    mu = sub[BIOMARKERS].mean()
    sd = sub[BIOMARKERS].std(ddof=1)
    Z = ((sub[BIOMARKERS] - mu) / sd).values

    covars = sub[[BATCH_COL] + COVARIATES].reset_index(drop=True)

    # ComBat on features
    out = neuroCombat(
        dat=Z.T, covars=covars, batch_col=BATCH_COL,
        categorical_cols=["Sex"], continuous_cols=["Age"],
        eb=True, parametric=True, mean_only=False,
    )
    Xc = out["data"].T

    # residuals on biological design only
    D_bio = np.column_stack([
        np.ones(n),
        sub["Age"].values.astype(float),
        pd.get_dummies(sub["Sex"], drop_first=True).values.astype(float),
    ])
    beta, *_ = np.linalg.lstsq(D_bio, Xc, rcond=None)
    R = Xc - D_bio @ beta

    # PCA via SVD
    U, s, Vt = np.linalg.svd(R, full_matrices=False)
    var_ratio = (s ** 2) / (s ** 2).sum()
    cum = np.cumsum(var_ratio)
    K = int(np.searchsorted(cum, PC_VAR_KEEP) + 1)
    K = max(1, min(K, len(s)))

    summary = {
        "cohort": label,
        "n_complete": n,
        "n_pc_total": len(s),
        "k_kept": K,
        "pc_var_threshold": PC_VAR_KEEP,
        "cum_var_kept": float(cum[K - 1]),
        "var_pc1": float(var_ratio[0]),
        "var_pc2": float(var_ratio[1]) if len(s) > 1 else np.nan,
        "var_pc3": float(var_ratio[2]) if len(s) > 2 else np.nan,
        "var_pc4": float(var_ratio[3]) if len(s) > 3 else np.nan,
        "var_pc5": float(var_ratio[4]) if len(s) > 4 else np.nan,
    }

    per_pc = pd.DataFrame({
        "cohort": label,
        "pc": np.arange(1, len(s) + 1),
        "var_ratio": var_ratio,
        "cum_var": cum,
        "kept": np.arange(1, len(s) + 1) <= K,
    })
    return summary, per_pc


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(MASTER)
    hc = df[df["Pathology"].fillna("HC").str.upper().eq("HC")].copy()

    summaries, per_pcs = [], []
    for cohort_df, label in [(hc, "HC"), (df, "ALL")]:
        s, p = pc_retention(cohort_df, label)
        summaries.append(s)
        per_pcs.append(p)
        print(f"[{label}] N={s['n_complete']}  K={s['k_kept']}/{s['n_pc_total']}"
              f"  cum_var={s['cum_var_kept']:.4f}"
              f"  PC1={s['var_pc1']:.3f} PC2={s['var_pc2']:.3f}"
              f" PC3={s['var_pc3']:.3f}")

    pd.DataFrame(summaries).to_csv(OUT_DIR / "covbat_pc_retention.csv", index=False)
    pd.concat(per_pcs, ignore_index=True).to_csv(
        OUT_DIR / "covbat_pc_variance.csv", index=False)
    print(f"\n[done] -> {OUT_DIR/'covbat_pc_retention.csv'}")
    print(pd.DataFrame(summaries).to_string(index=False))


if __name__ == "__main__":
    main()
