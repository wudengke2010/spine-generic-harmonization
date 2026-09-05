# -*- coding: utf-8 -*-
"""
Step 6b — Leakage-free ML evaluation: harmonisation refit WITHIN each CV fold
==============================================================================
Addresses reviewer point 5 (round 2): the original ml_vendor_auc.py applied
cross-validation to FULL-SAMPLE harmonised features, so the harmonisation
models saw the test subjects during fitting. Here every fold refits the
harmonisation on training subjects only and applies the fitted transform to
the held-out fold (vendor label of a new subject is assumed known at
deployment, standard for ComBat apply).

Methods and their batch variable replicate the v1.3.0 pipeline exactly:
  ComBat        per-metric, batch=Manufacturer, covars Age+Sex, eb=False
  ComBat-joint  per-feature z-score (train mu/sd) -> ComBat eb=True jointly
  CovBat        ComBat-joint + PCA(95%) on residuals + per-PC ComBat eb=True
  RELIEF        ComBat-joint + PCA(95%) + FastICA + BH-FDR(q<0.05) zeroing
  LME/FE        per-metric mixedlm(y~Age+C(Sex), groups=vendor), FE fallback

RF protocol identical to ml_vendor_auc.py: 500 trees, min_samples_leaf=5,
balanced class weights, RepeatedStratifiedKFold(5 folds x 20 repeats, seed 42).

Extra outputs requested by the reviewer:
  - Hanley-McNeil SE and 95% CI for point AUCs (pathology task, per-class OvR)
  - paired comparison Original vs each method (per-repeat paired t + HM z-test)
  - per-class one-vs-rest vendor AUCs
  - LME/FE fallback rate per fold

Outputs -> results/step6_mlauc/:
  ml_leakfree_long.csv        method x task x repeat AUCs
  ml_leakfree_summary.csv     method x task mean +/- SD (+ point AUC, HM CI)
  ml_leakfree_vendorclass.csv per-class OvR AUC (ALL, point estimates)
  ml_leakfree_paired.csv      pathology paired tests Original vs method
  ml_leakfree_hmci.csv        HM SE/CI table
"""
from __future__ import annotations

import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
from numpy.linalg import inv
from scipy import stats
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold
from sklearn.decomposition import FastICA
import statsmodels.formula.api as smf

BASE = "E:/boshi/qm_harmonization_paper/results"
OUT = f"{BASE}/step6_mlauc"
BIOM = ["T2w_CSA", "GM_CSA_mm2", "MTR", "MTsat", "FA", "MD", "AD", "RD"]
N_FOLDS, N_REPEATS, SEED = 5, 20, 42
PC_VAR_KEEP, FDR_ALPHA = 0.95, 0.05

# ----------------------------------------------------------------------------
# ComBat core (exact replication of neuroCombat math, fit/apply split)
# ----------------------------------------------------------------------------
def _design(batch, age, sex, levels):
    """batch one-hot (full) + sex dummy (drop first, sorted levels) + age."""
    n = len(batch)
    B = np.zeros((n, len(levels)))
    for i, b in enumerate(batch):
        B[i, levels.index(b)] = 1.0
    s_levels = sorted(np.unique(sex))
    S = np.column_stack([(sex == s_levels[j]).astype(float)
                         for j in range(1, len(s_levels))]) if len(s_levels) > 1 else np.zeros((n, 0))
    return np.column_stack([B, S, np.asarray(age, float).reshape(-1, 1)])


def _aprior(d):
    m, s2 = np.mean(d), np.var(d, ddof=1)
    return (2 * s2 + m ** 2) / s2


def _bprior(d):
    m, s2 = np.mean(d), np.var(d, ddof=1)
    return (m * s2 + m ** 3) / s2


def _it_sol(sdat, g_hat, d_hat, g_bar, t2, a, b, conv=1e-4):
    """Parametric EB iteration (replicates neuroCombat it_sol)."""
    n = sdat.shape[1]
    g_old, d_old = g_hat.copy(), d_hat.copy()
    change = 1.0
    while change > conv:
        g_new = (t2 * n * g_hat + d_old * g_bar) / (t2 * n + d_old)
        sum2 = ((sdat - g_new.reshape(-1, 1)) ** 2).sum(axis=1)
        d_new = (0.5 * sum2 + b) / (n / 2.0 + a - 1.0)
        change = max(np.max(np.abs(g_new - g_old) / np.abs(g_old)),
                     np.max(np.abs(d_new - d_old) / np.abs(d_old)))
        g_old, d_old = g_new, d_new
    return g_new, d_new


def combat_fit(X, batch, age, sex, eb=True):
    """X: N x F. Returns parameter dict. batch=vendor (3 levels)."""
    levels = sorted(np.unique(batch))
    D = _design(batch, age, sex, levels)
    n_batch = len(levels)
    # B_hat (p x F)
    B_hat = inv(D.T @ D) @ D.T @ X
    # batch sizes
    ns = np.array([(np.asarray(batch) == b).sum() for b in levels], float)
    grand = (ns / ns.sum()) @ B_hat[:n_batch]              # F
    # pooled variance (divide by N, as neuroCombat)
    resid = X - D @ B_hat
    var_pooled = (resid ** 2).mean(axis=0)                 # F
    var_pooled[var_pooled == 0] = np.median(var_pooled[var_pooled != 0]) if (var_pooled != 0).any() else 1.0
    # mod_mean (covariate part only)
    D_nb = D.copy()
    D_nb[:, :n_batch] = 0.0
    mod_mean = D_nb @ B_hat                                # N x F
    s_data = (X - grand - mod_mean) / np.sqrt(var_pooled)  # N x F
    # gamma/delta per batch
    gamma_hat = np.zeros((n_batch, X.shape[1]))
    delta_hat = np.zeros((n_batch, X.shape[1]))
    for i, b in enumerate(levels):
        m = np.asarray(batch) == b
        gamma_hat[i] = s_data[m].mean(axis=0)
        d = s_data[m].var(axis=0, ddof=1)
        d[d == 0] = 1.0
        delta_hat[i] = d
    if eb:
        gamma_star = np.zeros_like(gamma_hat)
        delta_star = np.zeros_like(delta_hat)
        # EB priors borrow strength ACROSS FEATURES (neuroCombat: gamma_bar and
        # t2 are per-batch scalars computed over the feature axis)
        gamma_bar = gamma_hat.mean(axis=1)                 # (n_batch,)
        t2 = gamma_hat.var(axis=1, ddof=1)                 # (n_batch,)
        for i, b in enumerate(levels):
            m = np.asarray(batch) == b
            # neuroCombat: aprior/bprior on the per-feature delta_hat vector of
            # batch i -> a SCALAR prior shared by all features of the batch
            a_p = _aprior(delta_hat[i])
            b_p = _bprior(delta_hat[i])
            g, d = _it_sol(s_data[m].T, gamma_hat[i], delta_hat[i],
                           gamma_bar[i], t2[i], a_p, b_p)
            gamma_star[i], delta_star[i] = g, d
    else:
        gamma_star, delta_star = gamma_hat, delta_hat
    return {"levels": levels, "B_hat": B_hat, "n_batch": n_batch,
            "grand": grand, "var_pooled": var_pooled,
            "gamma_star": gamma_star, "delta_star": delta_star}


def combat_apply(params, X, batch, age, sex):
    levels = params["levels"]
    D = _design(batch, age, sex, levels)
    n_batch = params["n_batch"]
    D_nb = D.copy()
    D_nb[:, :n_batch] = 0.0
    mod_mean = D_nb @ params["B_hat"]
    s = (X - params["grand"] - mod_mean) / np.sqrt(params["var_pooled"])
    out = np.empty_like(s)
    bidx = {b: i for i, b in enumerate(levels)}
    for j, b in enumerate(batch):
        out[j] = (s[j] - params["gamma_star"][bidx[b]]) / np.sqrt(
            params["delta_star"][bidx[b]])
    return out * np.sqrt(params["var_pooled"]) + params["grand"] + mod_mean


# ----------------------------------------------------------------------------
# Five harmonisers: fit on train, apply to test
# ----------------------------------------------------------------------------
def harm_combat(Xtr, mtr, Xte, mte):
    """Per-metric ComBat, eb=False (identical to joint fit with eb=False)."""
    p = combat_fit(Xtr, mtr["vendor"], mtr["age"], mtr["sex"], eb=False)
    return combat_apply(p, Xtr, mtr["vendor"], mtr["age"], mtr["sex"]), \
           combat_apply(p, Xte, mte["vendor"], mte["age"], mte["sex"])


def harm_combat_joint(Xtr, mtr, Xte, mte):
    mu, sd = Xtr.mean(axis=0), Xtr.std(axis=0, ddof=1)
    Ztr, Zte = (Xtr - mu) / sd, (Xte - mu) / sd
    p = combat_fit(Ztr, mtr["vendor"], mtr["age"], mtr["sex"], eb=True)
    Htr = combat_apply(p, Ztr, mtr["vendor"], mtr["age"], mtr["sex"])
    Hte = combat_apply(p, Zte, mte["vendor"], mte["age"], mte["sex"])
    return Htr * sd + mu, Hte * sd + mu


def _bio_design(age, sex):
    s_levels = sorted(np.unique(sex))
    S = np.column_stack([(np.asarray(sex) == s_levels[j]).astype(float)
                         for j in range(1, len(s_levels))]) if len(s_levels) > 1 else np.zeros((len(sex), 0))
    return np.column_stack([np.ones(len(age)), np.asarray(age, float).reshape(-1, 1), S])


def harm_covbat(Xtr, mtr, Xte, mte):
    mu, sd = Xtr.mean(axis=0), Xtr.std(axis=0, ddof=1)
    Ztr, Zte = (Xtr - mu) / sd, (Xte - mu) / sd
    # step 1: ComBat eb=True
    p1 = combat_fit(Ztr, mtr["vendor"], mtr["age"], mtr["sex"], eb=True)
    Ctr = combat_apply(p1, Ztr, mtr["vendor"], mtr["age"], mtr["sex"])
    Cte = combat_apply(p1, Zte, mte["vendor"], mte["age"], mte["sex"])
    # step 2: bio residuals (beta from train)
    Db_tr = _bio_design(mtr["age"], mtr["sex"])
    Db_te = _bio_design(mte["age"], mte["sex"])
    beta, *_ = np.linalg.lstsq(Db_tr, Ctr, rcond=None)
    Rtr, Rte = Ctr - Db_tr @ beta, Cte - Db_te @ beta
    # step 3: PCA on train residuals
    U, s, Vt = np.linalg.svd(Rtr, full_matrices=False)
    var_ratio = (s ** 2) / (s ** 2).sum()
    K = int(np.searchsorted(np.cumsum(var_ratio), PC_VAR_KEEP) + 1)
    K = max(1, min(K, len(s)))
    Str = U[:, :K] * s[:K]                     # N_tr x K
    Ste = Rte @ Vt[:K, :].T                    # N_te x K
    # step 4: ComBat on PC scores
    p2 = combat_fit(Str, mtr["vendor"], mtr["age"], mtr["sex"], eb=True)
    SHtr = combat_apply(p2, Str, mtr["vendor"], mtr["age"], mtr["sex"])
    SHte = combat_apply(p2, Ste, mte["vendor"], mte["age"], mte["sex"])
    # step 5: reconstruct
    Rh_tr, Rh_te = SHtr @ Vt[:K, :], SHte @ Vt[:K, :]
    return (Db_tr @ beta + Rh_tr) * sd + mu, (Db_te @ beta + Rh_te) * sd + mu


def _bh_fdr(p):
    p = np.asarray(p, float)
    m = len(p)
    order = np.argsort(p)
    ranked = p[order]
    q = ranked * m / (np.arange(m) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty_like(q)
    out[order] = np.clip(q, 0, 1)
    return out


def harm_relief(Xtr, mtr, Xte, mte):
    mu, sd = Xtr.mean(axis=0), Xtr.std(axis=0, ddof=1)
    Ztr, Zte = (Xtr - mu) / sd, (Xte - mu) / sd
    p1 = combat_fit(Ztr, mtr["vendor"], mtr["age"], mtr["sex"], eb=True)
    Ctr = combat_apply(p1, Ztr, mtr["vendor"], mtr["age"], mtr["sex"])
    Cte = combat_apply(p1, Zte, mte["vendor"], mte["age"], mte["sex"])
    Db_tr = _bio_design(mtr["age"], mtr["sex"])
    Db_te = _bio_design(mte["age"], mte["sex"])
    beta, *_ = np.linalg.lstsq(Db_tr, Ctr, rcond=None)
    Rtr, Rte = Ctr - Db_tr @ beta, Cte - Db_te @ beta
    U, s, Vt = np.linalg.svd(Rtr, full_matrices=False)
    var_ratio = (s ** 2) / (s ** 2).sum()
    K = int(np.searchsorted(np.cumsum(var_ratio), PC_VAR_KEEP) + 1)
    K = max(2, min(K, len(s)))
    PCtr = U[:, :K] * s[:K]
    PCte = Rte @ Vt[:K, :].T
    ica = FastICA(n_components=K, random_state=0, whiten="unit-variance",
                  max_iter=2000, tol=1e-5)
    Str = ica.fit_transform(PCtr)
    A = ica.mixing_
    Ste = ica.transform(PCte)
    # F-test each IC vs vendor on TRAIN
    Fs, ps = [], []
    for k in range(K):
        groups = [Str[np.asarray(mtr["vendor"]) == b, k]
                  for b in np.unique(mtr["vendor"])]
        groups = [g for g in groups if len(g) > 1]
        if len(groups) < 2:
            Fs.append(np.nan); ps.append(np.nan); continue
        f_, p_ = stats.f_oneway(*groups)
        Fs.append(f_); ps.append(p_)
    qs = _bh_fdr(np.nan_to_num(np.array(ps), nan=1.0))
    flag = qs < FDR_ALPHA
    Str_c, Ste_c = Str.copy(), Ste.copy()
    Str_c[:, flag] = 0.0
    Ste_c[:, flag] = 0.0
    Rh_tr = (Str_c @ A.T) @ Vt[:K, :]
    Rh_te = (Ste_c @ A.T) @ Vt[:K, :]
    return (Db_tr @ beta + Rh_tr) * sd + mu, (Db_te @ beta + Rh_te) * sd + mu


def harm_lme(Xtr, mtr, Xte, mte):
    """Per-metric mixedlm(y~Age+C(Sex), groups=vendor) with FE fallback."""
    n_f = Xtr.shape[1]
    levels = sorted(np.unique(mtr["vendor"]))
    Htr = np.empty_like(Xtr)
    Hte = np.empty_like(Xte)
    n_fallback = 0
    for f in range(n_f):
        d = pd.DataFrame({"y": Xtr[:, f], "Age": mtr["age"], "Sex": mtr["sex"],
                          "b": mtr["vendor"]})
        ok = False
        try:
            md = smf.mixedlm("y ~ Age + C(Sex)", data=d, groups=d["b"])
            mdf = md.fit(method="lbfgs", reml=True)
            blup = {g: float(v.iloc[0]) for g, v in mdf.random_effects.items()}
            if set(levels).issubset(blup.keys()) and all(np.isfinite(v) for v in blup.values()):
                Htr[:, f] = Xtr[:, f] - d["b"].map(blup).values
                Hte[:, f] = Xte[:, f] - pd.Series(mte["vendor"]).map(blup).values
                ok = True
        except Exception:
            ok = False
        if not ok:
            n_fallback += 1
            Xd = pd.get_dummies(d[["b", "Age", "Sex"]], drop_first=True).astype(float)
            Xd = pd.concat([pd.Series(1.0, index=Xd.index, name="const"), Xd], axis=1)
            import statsmodels.api as sm
            res = sm.OLS(d["y"].values, Xd).fit()
            offs = {v: float(res.params.get(f"b_{v}", 0.0)) for v in levels}
            mean_off = np.mean(list(offs.values()))
            offs = {v: o - mean_off for v, o in offs.items()}
            Htr[:, f] = Xtr[:, f] - d["b"].map(offs).values
            Hte[:, f] = Xte[:, f] - pd.Series(mte["vendor"]).map(offs).values
    return Htr, Hte, n_fallback


HARMONIZERS = {
    "ComBat": harm_combat,
    "ComBat-joint": harm_combat_joint,
    "CovBat": harm_covbat,
    "RELIEF": harm_relief,
    "LME": None,  # special-cased (returns fallback count)
}


# ----------------------------------------------------------------------------
# Hanley-McNeil SE / CI
# ----------------------------------------------------------------------------
def hm_se(auc, n_pos, n_neg):
    q1 = auc / (2 - auc)
    q2 = 2 * auc ** 2 / (1 + auc)
    v = (auc * (1 - auc) + (n_pos - 1) * (q1 - auc ** 2)
         + (n_neg - 1) * (q2 - auc ** 2)) / (n_pos * n_neg)
    return np.sqrt(v)


# ----------------------------------------------------------------------------
# Main evaluation
# ----------------------------------------------------------------------------
def run_task(X, meta, y, method, task, seed):
    """Returns per-repeat AUC list, pooled predictions, and fallback count."""
    classes = np.unique(y)
    multiclass = len(classes) > 2
    y_arr = np.asarray(y)
    rng_seed = seed
    aucs, fallbacks = [], 0
    pooled = np.zeros((len(y_arr), len(classes))) if multiclass else np.zeros(len(y_arr))
    for r in range(N_REPEATS):
        skf = StratifiedKFold(n_splits=N_FOLDS, shuffle=True,
                              random_state=seed + r)
        if multiclass:
            p = np.zeros((len(y_arr), len(classes)))
        else:
            p = np.zeros(len(y_arr))
        for tr, te in skf.split(X, y_arr):
            Xtr, Xte = X[tr], X[te]
            mtr = {k: np.asarray(v)[tr] for k, v in meta.items()}
            mte = {k: np.asarray(v)[te] for k, v in meta.items()}
            if method == "Original":
                Htr, Hte = Xtr, Xte
            elif method == "LME":
                Htr, Hte, fb = harm_lme(Xtr, mtr, Xte, mte)
                fallbacks += fb
            else:
                Htr, Hte = HARMONIZERS[method](Xtr, mtr, Xte, mte)
            rf = RandomForestClassifier(
                n_estimators=500, min_samples_leaf=5, class_weight="balanced",
                random_state=seed + r)
            rf.fit(Htr, y_arr[tr])
            if multiclass:
                p[te] = rf.predict_proba(Hte)
            else:
                p[te] = rf.predict_proba(Hte)[:, 1]
        pooled += p
        if multiclass:
            aucs.append(roc_auc_score(y_arr, p, multi_class="ovr",
                                      average="macro", labels=classes))
        else:
            aucs.append(roc_auc_score(y_arr, p))
    pooled /= N_REPEATS
    return aucs, pooled, fallbacks


from sklearn.metrics import roc_auc_score  # noqa: E402


def main():
    import os
    os.makedirs(OUT, exist_ok=True)
    master = pd.read_csv(f"{BASE}/biomarkers_master.csv")
    cc = master.dropna(subset=BIOM + ["Age", "Sex", "Manufacturer"]).reset_index(drop=True)
    hc = cc[cc["Pathology"].fillna("HC").str.upper().eq("HC")].reset_index(drop=True)
    print(f"ALL complete-case N={len(cc)}  HC N={len(hc)}")

    X_all = cc[BIOM].to_numpy()
    meta_all = {"vendor": cc["Manufacturer"].values, "age": cc["Age"].values,
                "sex": cc["Sex"].values}
    X_hc = hc[BIOM].to_numpy()
    meta_hc = {"vendor": hc["Manufacturer"].values, "age": hc["Age"].values,
               "sex": hc["Sex"].values}

    tasks = [
        ("ALL", "vendor", X_all, meta_all, cc["Manufacturer"].values),
        ("ALL", "pathology", X_all, meta_all, (cc["Pathology"] != "HC").astype(int).values),
        ("ALL", "sex", X_all, meta_all, (cc["Sex"] == "M").astype(int).values),
        ("HC", "vendor", X_hc, meta_hc, hc["Manufacturer"].values),
    ]

    methods = ["Original", "ComBat", "ComBat-joint", "CovBat", "RELIEF", "LME"]
    rows_long, summary_rows, pooled_store = [], [], {}

    for cohort, task, X, meta, y in tasks:
        for method in methods:
            aucs, pooled, fb = run_task(X, meta, y, method, task, SEED)
            pooled_store[(cohort, task, method)] = pooled
            m, s = float(np.mean(aucs)), float(np.std(aucs))
            # point AUC from pooled (mean over repeats) predictions
            classes = np.unique(y)
            if len(classes) > 2:
                point = roc_auc_score(y, pooled, multi_class="ovr",
                                      average="macro", labels=classes)
            else:
                point = roc_auc_score(y, pooled)
            summary_rows.append({"cohort": cohort, "method": method, "task": task,
                                 "auc_mean": round(m, 4), "auc_sd": round(s, 4),
                                 "auc_pooled": round(float(point), 4),
                                 "n": len(y), "lme_fallbacks": fb})
            for r, v in enumerate(aucs):
                rows_long.append({"cohort": cohort, "method": method,
                                  "task": task, "repeat": r, "auc": v})
            print(f"{cohort:3s} {task:9s} {method:12s} "
                  f"AUC={m:.3f}±{s:.3f}  pooled={point:.3f}"
                  + (f"  fb={fb}" if method == "LME" else ""))

    # permutation reference (ALL vendor, Original features, shuffled labels)
    rng = np.random.default_rng(SEED)
    perm = []
    yv = cc["Manufacturer"].values
    classes = np.unique(yv)
    for r in range(N_REPEATS):
        y_perm = rng.permutation(yv)
        skf = StratifiedKFold(n_splits=N_FOLDS, shuffle=True, random_state=SEED + r)
        p = np.zeros((len(yv), len(classes)))
        for tr, te in skf.split(X_all, y_perm):
            rf = RandomForestClassifier(n_estimators=500, min_samples_leaf=5,
                                        class_weight="balanced",
                                        random_state=SEED + r)
            rf.fit(X_all[tr], y_perm[tr])
            p[te] = rf.predict_proba(X_all[te])
        perm.append(roc_auc_score(y_perm, p, multi_class="ovr",
                                  average="macro", labels=classes))
    summary_rows.append({"cohort": "ALL", "method": "Permutation", "task": "vendor",
                         "auc_mean": round(float(np.mean(perm)), 4),
                         "auc_sd": round(float(np.std(perm)), 4),
                         "auc_pooled": np.nan, "n": len(yv), "lme_fallbacks": 0})
    for r, v in enumerate(perm):
        rows_long.append({"cohort": "ALL", "method": "Permutation", "task": "vendor",
                          "repeat": r, "auc": v})
    print(f"Permutation vendor AUC = {np.mean(perm):.3f}±{np.std(perm):.3f}")

    pd.DataFrame(rows_long).to_csv(f"{OUT}/ml_leakfree_long.csv", index=False)
    pd.DataFrame(summary_rows).to_csv(f"{OUT}/ml_leakfree_summary.csv", index=False)

    # ---- per-class OvR AUC (ALL vendor, pooled predictions) ----
    vc_rows = []
    classes = np.unique(yv)
    for method in methods:
        pooled = pooled_store[("ALL", "vendor", method)]
        for i, c in enumerate(classes):
            y_bin = (yv == c).astype(int)
            auc = roc_auc_score(y_bin, pooled[:, i])
            n_pos, n_neg = y_bin.sum(), (1 - y_bin).sum()
            se = hm_se(auc, int(n_pos), int(n_neg))
            vc_rows.append({"method": method, "vendor": c, "auc": round(auc, 4),
                            "hm_se": round(se, 4),
                            "ci95_lo": round(auc - 1.96 * se, 4),
                            "ci95_hi": round(auc + 1.96 * se, 4),
                            "n_pos": int(n_pos), "n_neg": int(n_neg)})
    pd.DataFrame(vc_rows).to_csv(f"{OUT}/ml_leakfree_vendorclass.csv", index=False)

    # ---- pathology: HM CI + paired tests Original vs each method ----
    y_path = (cc["Pathology"] != "HC").astype(int).values
    n_pos, n_neg = int(y_path.sum()), int((1 - y_path).sum())
    hm_rows, paired_rows = [], []
    long_df = pd.DataFrame(rows_long)
    for method in methods:
        pooled = pooled_store[("ALL", "pathology", method)]
        point = roc_auc_score(y_path, pooled)
        se = hm_se(point, n_pos, n_neg)
        hm_rows.append({"method": method, "auc": round(point, 4),
                        "hm_se": round(se, 4),
                        "ci95_lo": round(point - 1.96 * se, 4),
                        "ci95_hi": round(point + 1.96 * se, 4),
                        "n_pos": n_pos, "n_neg": n_neg})
        if method != "Original":
            a0 = long_df[(long_df.cohort == "ALL") & (long_df.task == "pathology")
                         & (long_df.method == "Original")]["auc"].values
            a1 = long_df[(long_df.cohort == "ALL") & (long_df.task == "pathology")
                         & (long_df.method == method)]["auc"].values
            t, p_t = stats.ttest_rel(a1, a0)
            try:
                w, p_w = stats.wilcoxon(a1, a0)
            except Exception:
                p_w = np.nan
            paired_rows.append({"comparison": f"{method} vs Original",
                                "mean_diff": round(float(np.mean(a1 - a0)), 4),
                                "t": round(float(t), 3), "p_ttest": round(float(p_t), 5),
                                "p_wilcoxon": round(float(p_w), 5) if p_w == p_w else np.nan})
    # HM z-test for point AUC difference vs Original
    orig = hm_rows[0]
    for row in hm_rows[1:]:
        se_d = np.sqrt(row["hm_se"] ** 2 + orig["hm_se"] ** 2)
        z = (row["auc"] - orig["auc"]) / se_d
        paired_rows.append({"comparison": f"z-test {row['method']} vs Original",
                            "mean_diff": np.nan, "t": round(float(z), 3),
                            "p_ttest": round(float(2 * (1 - stats.norm.cdf(abs(z)))), 5),
                            "p_wilcoxon": np.nan})
    pd.DataFrame(hm_rows).to_csv(f"{OUT}/ml_leakfree_hmci.csv", index=False)
    pd.DataFrame(paired_rows).to_csv(f"{OUT}/ml_leakfree_paired.csv", index=False)

    print("\n=== pathology HM CI ===")
    print(pd.DataFrame(hm_rows).to_string(index=False))
    print("\n=== paired tests ===")
    print(pd.DataFrame(paired_rows).to_string(index=False))
    print("\nSaved to", OUT)


if __name__ == "__main__":
    main()
