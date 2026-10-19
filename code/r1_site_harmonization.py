# -*- coding: utf-8 -*-
"""
Reviewer 1 (NMRB 7830229) — "As there are only 3 vendors, the authors should
consider harmonizing within a vendor with scanner models or sites."

We refit per-metric ComBat and joint ComBat on the HC complete-case sample
(N=188) using, as the batch variable:
  S1. site (40 levels; sites are perfectly nested within vendors)
  S2. scanner model (12 levels; models are nested within vendors but some
      models are shared across sites)
and compare with the canonical vendor-batch fits (3 levels).

Endpoints (identical definitions to the main analysis):
  vendor eta^2 (ANOVA, mean over 8 biomarkers)   [target: vendor removal]
  site eta^2 (ANOVA over 40 sites, mean)         [target: site removal]
  PERMANOVA R^2 (z-scored Euclidean; equals mean per-metric eta^2)
  age / sex semi-partial R^2 (Type II),          [biological recovery]
  r_pre/post (subject-level preservation),       [preservation]
  plus per-vendor residual means to check within-vendor structure.

Note on EB: with 40 site batches (n=2-8 per site), empirical-Bayes shrinkage
across batches is essential; we therefore report eb=True for the site and
model fits. For the vendor-batch reference we use the canonical settings of
the main analysis (ComBat per-metric eb=False; ComBat-joint eb=True).

Outputs -> results/step4_eval/r1_site_harmonization.csv
"""
from __future__ import annotations
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
from numpy.linalg import inv
from scipy import stats

BASE = "E:/boshi/qm_harmonization_paper/results"
BIOM = ["T2w_CSA", "GM_CSA_mm2", "MTR", "MTsat", "FA", "MD", "AD", "RD"]

# ---------------------------------------------------------------- ComBat ----
def _design(batch, age, sex, levels):
    n = len(batch)
    B = np.zeros((n, len(levels)))
    for i, b in enumerate(batch):
        B[i, levels.index(b)] = 1.0
    s_levels = sorted(np.unique(sex))
    S = np.column_stack([(sex == s_levels[j]).astype(float)
                         for j in range(1, len(s_levels))]) if len(s_levels) > 1 else np.zeros((n, 0))
    return np.column_stack([B, S, np.asarray(age, float).reshape(-1, 1)])

def _aprior(d):
    m, s2 = np.mean(d), np.var(d, ddof=1); return (2 * s2 + m ** 2) / s2

def _bprior(d):
    m, s2 = np.mean(d), np.var(d, ddof=1); return (m * s2 + m ** 3) / s2

def _it_sol(sdat, g_hat, d_hat, g_bar, t2, a, b, conv=1e-4):
    g_old, d_old = g_hat.copy(), d_hat.copy()
    change = 1.0
    while change > conv:
        g_new = (t2 * len(g_hat) * 0 + t2 * sdat.shape[1] * g_hat + d_old * g_bar) / (t2 * sdat.shape[1] + d_old)
        sum2 = ((sdat - g_new.reshape(-1, 1)) ** 2).sum(axis=1)
        d_new = (0.5 * sum2 + b) / (sdat.shape[1] / 2.0 + a - 1.0)
        change = max(np.max(np.abs(g_new - g_old) / np.abs(g_old)),
                     np.max(np.abs(d_new - d_old) / np.abs(d_old)))
        g_old, d_old = g_new, d_new
    return g_new, d_new

def combat_fit(X, batch, age, sex, eb=True):
    levels = sorted(np.unique(batch))
    D = _design(batch, age, sex, levels)
    n_batch = len(levels)
    B_hat = inv(D.T @ D) @ D.T @ X
    ns = np.array([(np.asarray(batch) == b).sum() for b in levels], float)
    grand = (ns / ns.sum()) @ B_hat[:n_batch]
    resid = X - D @ B_hat
    var_pooled = (resid ** 2).mean(axis=0)
    var_pooled[var_pooled == 0] = 1.0
    D_nb = D.copy(); D_nb[:, :n_batch] = 0.0
    mod_mean = D_nb @ B_hat
    s_data = (X - grand - mod_mean) / np.sqrt(var_pooled)
    gamma_hat = np.zeros((n_batch, X.shape[1])); delta_hat = np.zeros_like(gamma_hat)
    for i, b in enumerate(levels):
        m = np.asarray(batch) == b
        gamma_hat[i] = s_data[m].mean(axis=0)
        d = s_data[m].var(axis=0, ddof=1)
        d[d == 0] = 1.0
        delta_hat[i] = d
    if eb:
        gamma_star = np.zeros_like(gamma_hat); delta_star = np.zeros_like(delta_hat)
        gamma_bar = gamma_hat.mean(axis=1); t2 = gamma_hat.var(axis=1, ddof=1)
        for i, b in enumerate(levels):
            m = np.asarray(batch) == b
            a_p, b_p = _aprior(delta_hat[i]), _bprior(delta_hat[i])
            g, d = _it_sol(s_data[m].T, gamma_hat[i], delta_hat[i],
                           gamma_bar[i], t2[i], a_p, b_p)
            gamma_star[i], delta_star[i] = g, d
    else:
        gamma_star, delta_star = gamma_hat, delta_hat
    return {"levels": levels, "B_hat": B_hat, "n_batch": n_batch, "grand": grand,
            "var_pooled": var_pooled, "gamma_star": gamma_star, "delta_star": delta_star}

def combat_apply(p, X, batch, age, sex):
    D = _design(batch, age, sex, p["levels"])
    D_nb = D.copy(); D_nb[:, :p["n_batch"]] = 0.0
    mod_mean = D_nb @ p["B_hat"]
    s = (X - p["grand"] - mod_mean) / np.sqrt(p["var_pooled"])
    out = np.empty_like(s)
    bidx = {b: i for i, b in enumerate(p["levels"])}
    for j, b in enumerate(batch):
        out[j] = (s[j] - p["gamma_star"][bidx[b]]) / np.sqrt(p["delta_star"][bidx[b]])
    return out * np.sqrt(p["var_pooled"]) + p["grand"] + mod_mean

def combat_joint(X, batch, age, sex, eb=True):
    mu, sd = X.mean(axis=0), X.std(axis=0, ddof=1)
    Z = (X - mu) / sd
    p = combat_fit(Z, batch, age, sex, eb=eb)
    return combat_apply(p, Z, batch, age, sex) * sd + mu

# ------------------------------------------------------------- endpoints ----
def eta2(X, groups):
    vals = []
    groups = np.asarray(groups)
    for f in range(X.shape[1]):
        y = X[:, f]
        ss_t = ((y - y.mean()) ** 2).sum()
        ss_b = sum(((y[groups == g].mean() - y.mean()) ** 2) * (groups == g).sum()
                   for g in np.unique(groups))
        vals.append(ss_b / ss_t)
    return np.array(vals)

def semi_partial_R2(X, age, sex):
    """Type-II semi-partial R^2 for age and sex, mean over biomarkers."""
    n = X.shape[0]
    s_levels = sorted(np.unique(sex))
    S = np.column_stack([(np.asarray(sex) == s_levels[j]).astype(float)
                         for j in range(1, len(s_levels))])
    age = np.asarray(age, float).reshape(-1, 1)
    Dfull = np.column_stack([np.ones(n), age, S])
    Dage = np.column_stack([np.ones(n), S])
    Dsex = np.column_stack([np.ones(n), age])
    r2_age, r2_sex = [], []
    for f in range(X.shape[1]):
        y = X[:, f]
        sst = ((y - y.mean()) ** 2).sum()
        rss_full = ((y - Dfull @ np.linalg.lstsq(Dfull, y, rcond=None)[0]) ** 2).sum()
        rss_age = ((y - Dage @ np.linalg.lstsq(Dage, y, rcond=None)[0]) ** 2).sum()
        rss_sex = ((y - Dsex @ np.linalg.lstsq(Dsex, y, rcond=None)[0]) ** 2).sum()
        r2_age.append((rss_age - rss_full) / sst)
        r2_sex.append((rss_sex - rss_full) / sst)
    return float(np.mean(r2_age)), float(np.mean(r2_sex))

def r_prepost(Xpre, Xpost):
    """Main-analysis definition (step4_eval_completecase.py): for each
    biomarker, Pearson r across subjects between pre- and post-harmonisation
    values (subject-ranking preservation), averaged over biomarkers."""
    rs = []
    for f in range(Xpre.shape[1]):
        x, y = Xpre[:, f], Xpost[:, f]
        if x.std() > 0 and y.std() > 0:
            rs.append(np.corrcoef(x, y)[0, 1])
    return float(np.mean(rs))

def main():
    master = pd.read_csv(f"{BASE}/biomarkers_master.csv")
    hc = master.dropna(subset=BIOM + ["Age", "Sex", "Manufacturer"])
    hc = hc[hc["Pathology"].fillna("HC").str.upper().eq("HC")].reset_index(drop=True)
    print(f"HC complete-case N={len(hc)}, sites={hc['Site'].nunique()}, "
          f"models={hc['ScannerModel'].nunique()}, vendors={hc['Manufacturer'].nunique()}")
    X = hc[BIOM].to_numpy()
    vendor = hc["Manufacturer"].values
    site = hc["Site"].values
    model = hc["ScannerModel"].values
    age = hc["Age"].values.astype(float)
    sex = hc["Sex"].values

    fits = {
        "Original": X,
        "ComBat (batch=vendor, eb=F)": combat_apply(
            combat_fit(X, vendor, age, sex, eb=False), X, vendor, age, sex),
        "ComBat (batch=site, eb=T)": combat_apply(
            combat_fit(X, site, age, sex, eb=True), X, site, age, sex),
        "ComBat (batch=model, eb=T)": combat_apply(
            combat_fit(X, model, age, sex, eb=True), X, model, age, sex),
        "ComBat-joint (batch=vendor, eb=T)": combat_joint(X, vendor, age, sex, eb=True),
        "ComBat-joint (batch=site, eb=T)": combat_joint(X, site, age, sex, eb=True),
        "ComBat-joint (batch=model, eb=T)": combat_joint(X, model, age, sex, eb=True),
    }

    rows = []
    for name, H in fits.items():
        ev = eta2(H, vendor); es = eta2(H, site)
        r2a, r2s = semi_partial_R2(H, age, sex)
        rows.append({
            "fit": name,
            "vendor_eta2_mean": round(float(ev.mean()), 6),
            "site_eta2_mean": round(float(es.mean()), 6),
            "permanova_R2_vendor": round(float(ev.mean()), 6),  # identity (z-scored Euclidean)
            "age_R2": round(r2a, 5),
            "sex_R2": round(r2s, 5),
            "r_prepost": round(r_prepost(X, H), 4),
        })
        print(f"{name:38s} vendor_eta2={ev.mean():.6f}  site_eta2={es.mean():.4f}  "
              f"ageR2={r2a:.4f}  sexR2={r2s:.4f}  r={rows[-1]['r_prepost']:.4f}")
    out = pd.DataFrame(rows)
    out.to_csv(f"{BASE}/step4_eval/r1_site_harmonization.csv", index=False)
    print("\nSaved:", f"{BASE}/step4_eval/r1_site_harmonization.csv")

if __name__ == "__main__":
    main()
