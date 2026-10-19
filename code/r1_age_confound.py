# -*- coding: utf-8 -*-
"""
Reviewer 1 (NMRB 7830229) — age-confound sensitivity for the pathology
classification endpoint.

R1 points out that the 'pathology' label is a radiological assessment (mostly
mild cervical compression, Valošek et al. 2024, doi:10.1162/imag_a_00075),
that mild compression is age-related, and that the pathology group is
essentially mild compression (only 2 DCM). We therefore quantify how much of
the pathology classification signal is explained by age:

  A1. Univariate age difference between groups (MWU).
  A2. AUC of age alone for pathology (logistic-style: use age directly).
  A3. RF pathology AUC (same protocol as ml_fold_internal.py: 8 biomarkers,
      5-fold x 20 repeats, harmonisation refit per fold) computed on
      biomarkers residualised on age+sex -- for Original and ComBat.
  A4. For reference, pathology AUC WITHOUT residualisation is taken from
      ml_leakfree_summary.csv (already computed).

Outputs -> results/step6_mlauc/r1_age_confound.csv  (and stdout)
"""
from __future__ import annotations
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score

BASE = "E:/boshi/qm_harmonization_paper/results"
BIOM = ["T2w_CSA", "GM_CSA_mm2", "MTR", "MTsat", "FA", "MD", "AD", "RD"]
N_FOLDS, N_REPEATS, SEED = 5, 20, 42

# --- canonical ComBat (identical math to ml_fold_internal.py) ----------------
from numpy.linalg import inv

def _design(batch, age, sex, levels):
    n = len(batch)
    B = np.zeros((n, len(levels)))
    for i, b in enumerate(batch):
        B[i, levels.index(b)] = 1.0
    s_levels = sorted(np.unique(sex))
    S = np.column_stack([(sex == s_levels[j]).astype(float)
                         for j in range(1, len(s_levels))]) if len(s_levels) > 1 else np.zeros((n, 0))
    return np.column_stack([B, S, np.asarray(age, float).reshape(-1, 1)])

def combat_fit(X, batch, age, sex, eb=False):
    levels = sorted(np.unique(batch))
    D = _design(batch, age, sex, levels)
    n_batch = len(levels)
    B_hat = inv(D.T @ D) @ D.T @ X
    ns = np.array([(np.asarray(batch) == b).sum() for b in levels], float)
    grand = (ns / ns.sum()) @ B_hat[:n_batch]
    resid = X - D @ B_hat
    var_pooled = (resid ** 2).mean(axis=0)
    D_nb = D.copy(); D_nb[:, :n_batch] = 0.0
    mod_mean = D_nb @ B_hat
    s_data = (X - grand - mod_mean) / np.sqrt(var_pooled)
    gamma_hat = np.zeros((n_batch, X.shape[1])); delta_hat = np.zeros_like(gamma_hat)
    for i, b in enumerate(levels):
        m = np.asarray(batch) == b
        gamma_hat[i] = s_data[m].mean(axis=0)
        delta_hat[i] = s_data[m].var(axis=0, ddof=1)
    return {"levels": levels, "B_hat": B_hat, "n_batch": n_batch, "grand": grand,
            "var_pooled": var_pooled, "gamma_star": gamma_hat, "delta_star": delta_hat}

def combat_apply(p, X, batch, age, sex):
    levels = p["levels"]
    D = _design(batch, age, sex, levels)
    D_nb = D.copy(); D_nb[:, :p["n_batch"]] = 0.0
    mod_mean = D_nb @ p["B_hat"]
    s = (X - p["grand"] - mod_mean) / np.sqrt(p["var_pooled"])
    out = np.empty_like(s)
    bidx = {b: i for i, b in enumerate(levels)}
    for j, b in enumerate(batch):
        out[j] = (s[j] - p["gamma_star"][bidx[b]]) / np.sqrt(p["delta_star"][bidx[b]])
    return out * np.sqrt(p["var_pooled"]) + p["grand"] + mod_mean

def residualise(X, age, sex):
    """OLS-residualise each biomarker on [1, age, sex]."""
    n = X.shape[0]
    s_levels = sorted(np.unique(sex))
    S = np.column_stack([(np.asarray(sex) == s_levels[j]).astype(float)
                         for j in range(1, len(s_levels))])
    D = np.column_stack([np.ones(n), np.asarray(age, float), S])
    beta, *_ = np.linalg.lstsq(D, X, rcond=None)
    return X - D @ beta

def rf_auc(X, y, meta=None, method=None):
    """Same CV protocol as ml_fold_internal.run_task."""
    aucs = []
    y = np.asarray(y)
    for r in range(N_REPEATS):
        skf = StratifiedKFold(n_splits=N_FOLDS, shuffle=True, random_state=SEED + r)
        p = np.zeros(len(y))
        for tr, te in skf.split(X, y):
            Xtr, Xte = X[tr], X[te]
            if method == "ComBat":
                mtr = {k: np.asarray(v)[tr] for k, v in meta.items()}
                mte = {k: np.asarray(v)[te] for k, v in meta.items()}
                prm = combat_fit(Xtr, mtr["vendor"], mtr["age"], mtr["sex"])
                Htr = combat_apply(prm, Xtr, mtr["vendor"], mtr["age"], mtr["sex"])
                Hte = combat_apply(prm, Xte, mte["vendor"], mte["age"], mte["sex"])
            else:
                Htr, Hte = Xtr, Xte
            rf = RandomForestClassifier(n_estimators=500, min_samples_leaf=5,
                                        class_weight="balanced", random_state=SEED + r)
            rf.fit(Htr, y[tr])
            p[te] = rf.predict_proba(Hte)[:, 1]
        aucs.append(roc_auc_score(y, p))
    return float(np.mean(aucs)), float(np.std(aucs))

def main():
    master = pd.read_csv(f"{BASE}/biomarkers_master.csv")
    cc = master.dropna(subset=BIOM + ["Age", "Sex", "Manufacturer"]).reset_index(drop=True)
    y = (cc["Pathology"] != "HC").astype(int).values
    age = cc["Age"].values.astype(float)
    sex = cc["Sex"].values
    vendor = cc["Manufacturer"].values
    X = cc[BIOM].to_numpy()
    print(f"ALL complete-case N={len(cc)}, positives={y.sum()}")

    rows = []
    # A1: age difference
    u = stats.mannwhitneyu(age[y == 1], age[y == 0])
    rows.append({"analysis": "A1 age MWU (pathology vs HC)",
                 "value": f"median {np.median(age[y==1]):.1f} vs {np.median(age[y==0]):.1f}; p={u.pvalue:.2e}",
                 "auc_mean": np.nan, "auc_sd": np.nan})
    # A2: age alone AUC
    auc_age = roc_auc_score(y, age)
    rows.append({"analysis": "A2 age alone", "value": "", "auc_mean": auc_age, "auc_sd": np.nan})
    # sex alone for completeness
    auc_sex = roc_auc_score(y, (sex == "M").astype(float))
    rows.append({"analysis": "A2b sex alone", "value": "", "auc_mean": max(auc_sex, 1 - auc_sex), "auc_sd": np.nan})
    # age+sex logistic-style combined (linear score)
    D = np.column_stack([np.ones(len(y)), age, (sex == "M").astype(float)])
    beta, *_ = np.linalg.lstsq(np.column_stack([D]), y.astype(float), rcond=None)
    rows.append({"analysis": "A2c age+sex linear score", "value": "", "auc_mean": roc_auc_score(y, D @ beta), "auc_sd": np.nan})

    # A3: residualised biomarkers
    Xr = residualise(X, age, sex)
    meta = {"vendor": vendor, "age": age, "sex": sex}
    m, s = rf_auc(Xr, y)
    rows.append({"analysis": "A3 Original, biomarkers residualised on age+sex", "value": "", "auc_mean": m, "auc_sd": s})
    m, s = rf_auc(Xr, y, meta, method="ComBat")
    rows.append({"analysis": "A3 ComBat, biomarkers residualised on age+sex", "value": "", "auc_mean": m, "auc_sd": s})

    # A4 reference (from existing ml_leakfree_summary)
    try:
        summ = pd.read_csv(f"{BASE}/step6_mlauc/ml_leakfree_summary.csv")
        ref = summ[(summ.cohort == "ALL") & (summ.task == "pathology")]
        for _, rr in ref.iterrows():
            rows.append({"analysis": f"A4 reference (no residualisation): {rr['method']}",
                         "value": "", "auc_mean": rr["auc_mean"], "auc_sd": rr["auc_sd"]})
    except Exception as e:
        print("ref read failed:", e)

    out = pd.DataFrame(rows)
    out.to_csv(f"{BASE}/step6_mlauc/r1_age_confound.csv", index=False)
    print(out.to_string(index=False))

if __name__ == "__main__":
    main()
