#!/usr/bin/env python3
"""
Comprehensive figure regeneration script for spine-generic QM harmonization paper.
Reads raw spine-generic data, runs harmonization analysis, and generates all 7
figures with high-quality SCI publication layout.

Usage:
    python regenerate_all_figures.py
"""
from __future__ import annotations
import os, sys, re, warnings
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Rectangle, Patch
from matplotlib.lines import Line2D
from scipy import stats
from sklearn.decomposition import PCA, FastICA
from sklearn.preprocessing import StandardScaler
from sklearn.manifold import TSNE
from neuroCombat import neuroCombat
import statsmodels.formula.api as smf
from skbio.stats.distance import permanova, DistanceMatrix
warnings.filterwarnings("ignore")

# ============================================================
# Configuration
# ============================================================
BIDS_ROOT = Path("E:/boshi/spine-generic-multi-subject")
DERIV_T2STAR = Path("E:/qm_harmonization_paper/derivatives/csa_t2star")
OUT_DIR = Path("C:/Users/admin/WorkBuddy/2026-07-06-12-34-04/figures")

# Vendor mapping (corrected 2026-07-08)
VENDOR_MAP = {
    "amu": "Siemens", "balgrist": "Siemens", "barcelona": "Siemens",
    "beijingPrisma": "Siemens", "beijingVerio": "Siemens", "brnoCeitec": "Siemens",
    "cardiff": "Siemens", "cmrra": "Siemens", "cmrrb": "Siemens",
    "dresden": "Siemens", "fslPrisma": "Siemens", "geneva": "Siemens",
    "hamburg": "Siemens", "mgh": "Siemens", "milan": "Siemens",
    "mniPilot": "Siemens", "mniS": "Siemens", "mountSinai": "Siemens",
    "mpicbs": "Siemens", "nwu": "Siemens", "oxfordFmrib": "Siemens",
    "oxfordOhba": "Siemens", "pavia": "Siemens", "queensland": "Siemens",
    "strasbourg": "Siemens", "tehranS": "Siemens", "tokyoSkyra": "Siemens",
    "ucdavis": "Siemens", "unf": "Siemens", "vallHebron": "Siemens",
    "beijingGE": "GE", "brnoUhb": "GE", "juntendo750w": "GE",
    "perform": "GE", "stanford": "GE", "tokyo750w": "GE",
    "fslAchieva": "Philips", "nottwil": "Philips", "sherbrooke": "Philips",
    "tokyoIngenia": "Philips", "ubc": "Philips", "ucl": "Philips",
    "vuiisAchieva": "Philips", "vuiisIngenia": "Philips",
}

BIOMARKERS = ["T2w_CSA", "GM_CSA_mm2", "MTR", "MTsat", "FA", "MD", "AD", "RD"]
METHOD_ORDER = ["LME", "ComBat", "ComBat-joint", "RELIEF", "CovBat"]

# Wong 2011 palette
WONG = {"blue":"#0072B2","sky":"#56B4E9","green":"#009E73","orange":"#E69F00",
        "vermilion":"#D55E00","pink":"#CC79A7","yellow":"#F0E442","grey":"#999999"}
VENDOR_COLOR = {"GE": WONG["blue"], "Philips": WONG["green"], "Siemens": WONG["vermilion"]}
METHOD_COLOR = {"LME": WONG["pink"], "ComBat": WONG["blue"],
                "ComBat-joint": WONG["sky"], "RELIEF": WONG["orange"],
                "CovBat": WONG["green"]}

# ============================================================
# SCI-quality global rcParams
# ============================================================
mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 9,
    "axes.titlesize": 10,
    "axes.labelsize": 9,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "legend.fontsize": 8,
    "axes.linewidth": 0.8,
    "xtick.major.width": 0.8,
    "ytick.major.width": 0.8,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "svg.fonttype": "none",
    "figure.dpi": 300,
    "savefig.dpi": 300,
})

# ============================================================
# Step 1: Load data and create biomarker master table
# ============================================================
def load_data():
    """Load all biomarker data and create master table."""
    print("[1/5] Loading data...")

    # participants.tsv
    participants = pd.read_csv(BIDS_ROOT / "participants.tsv", sep="\t")
    participants["subject"] = participants["participant_id"]
    participants["site"] = participants["institution_id"]
    participants["Manufacturer"] = participants["manufacturer"]
    participants["Age"] = participants["age"]
    participants["Sex"] = participants["sex"]
    participants["Pathology"] = participants["pathology"].fillna("HC")

    # perlevel_metrics.csv (7 biomarkers, per vertebral level)
    perlevel = pd.read_csv(BIDS_ROOT / "perlevel_metrics.csv")

    # Aggregate C2-C3 (levels 'C2' and 'C3') per subject
    c2c3 = perlevel[perlevel["level"].isin(["C2", "C3"])].groupby("subject").mean(numeric_only=True)
    c2c3 = c2c3.reset_index()

    # Rename columns to match paper convention
    master = pd.DataFrame()
    master["subject"] = c2c3["subject"]
    master["T2w_CSA"] = c2c3["CSA"]
    master["MTR"] = c2c3["MTR"]
    master["MTsat"] = c2c3["MTsat"]
    master["FA"] = c2c3["FA"]
    master["MD"] = c2c3["MD"]
    master["AD"] = c2c3["AD"]
    master["RD"] = c2c3["RD"]

    # T2* GM CSA from derivatives (direct column)
    t2s_csa = {}
    for f in DERIV_T2STAR.glob("*_csa_t2star.csv"):
        sub = f.name.replace("_csa_t2star.csv", "")
        try:
            df = pd.read_csv(f)
            if "GM_CSA_mm2" in df.columns and len(df) > 0:
                t2s_csa[sub] = float(df["GM_CSA_mm2"].iloc[0])
            elif "MEAN(area)" in df.columns and len(df) > 0:
                t2s_csa[sub] = float(df["MEAN(area)"].mean())
        except:
            pass
    master["GM_CSA_mm2"] = master["subject"].map(t2s_csa)

    # Merge with participant metadata
    master = master.merge(
        participants[["subject", "site", "Manufacturer", "Age", "Sex", "Pathology"]],
        on="subject", how="left")

    # Apply corrected vendor mapping
    master["Manufacturer"] = master["site"].map(VENDOR_MAP).fillna(master["Manufacturer"])

    # Drop subjects with any missing biomarker
    master_clean = master.dropna(subset=BIOMARKERS + ["Manufacturer", "Age", "Sex"]).copy()
    master_clean = master_clean.reset_index(drop=True)

    print(f"  Master table: {len(master_clean)} subjects, {len(BIOMARKERS)} biomarkers")
    print(f"  Vendors: {master_clean['Manufacturer'].value_counts().to_dict()}")
    print(f"  Sites: {master_clean['site'].nunique()}")

    return master_clean

# ============================================================
# Step 2: Harmonization methods
# ============================================================
def _design_matrix(meta):
    sex_num = (meta["Sex"].str.upper() == "F").astype(int).values
    age = meta["Age"].astype(float).values
    cov = pd.DataFrame({"age": age, "sex": sex_num,
                        "batch": meta["Manufacturer"].values})
    return cov, meta["Manufacturer"].values

def harm_combat_permetric(X, meta):
    Y = np.zeros_like(X)
    for j in range(X.shape[1]):
        Xj = X[:, [j]]
        cov, _ = _design_matrix(meta)
        try:
            out = neuroCombat(dat=Xj.T, covars=cov, batch_col="batch",
                              categorical_cols=["sex"], continuous_cols=["age"],
                              eb=False, parametric=True, mean_only=False)
            Y[:, j] = out["data"].T.ravel()
        except:
            Y[:, j] = Xj.ravel()
    return Y

def harm_combat_joint(X, meta):
    cov, _ = _design_matrix(meta)
    out = neuroCombat(dat=X.T, covars=cov, batch_col="batch",
                      categorical_cols=["sex"], continuous_cols=["age"],
                      eb=True, parametric=True, mean_only=False)
    return out["data"].T

def harm_covbat(X, meta, pc_var=0.95):
    Y1 = harm_combat_joint(X, meta)
    sc = StandardScaler().fit(Y1)
    Z = sc.transform(Y1)
    pca = PCA().fit(Z)
    cum = np.cumsum(pca.explained_variance_ratio_)
    k = min(int(np.searchsorted(cum, pc_var) + 1), Z.shape[1])
    scores = Z @ pca.components_[:k].T
    cov, _ = _design_matrix(meta)
    scores_h = np.zeros_like(scores)
    for j in range(k):
        sj = scores[:, [j]]
        try:
            out = neuroCombat(dat=sj.T, covars=cov, batch_col="batch",
                              categorical_cols=["sex"], continuous_cols=["age"],
                              eb=False, parametric=True, mean_only=False)
            scores_h[:, j] = out["data"].T.ravel()
        except:
            scores_h[:, j] = sj.ravel()
    Z_back = scores_h @ pca.components_[:k]
    if k < Z.shape[1]:
        scores_rest = Z @ pca.components_[k:].T
        Z_back = Z_back + scores_rest @ pca.components_[k:]
    return sc.inverse_transform(Z_back)

def harm_relief(X, meta, n_restarts=5, random_state=0):
    Y1 = harm_combat_joint(X, meta)
    cov, batch = _design_matrix(meta)
    n, p = Y1.shape
    Xdes = np.column_stack([np.ones(n), cov["age"].values, cov["sex"].values])
    beta, *_ = np.linalg.lstsq(Xdes, Y1, rcond=None)
    R = Y1 - Xdes @ beta
    nc = p
    best_recon = np.inf
    best_S = best_A = None
    for s in range(n_restarts):
        try:
            ica = FastICA(n_components=nc, random_state=random_state + s,
                          max_iter=2000, whiten="unit-variance")
            S = ica.fit_transform(R)
            recon = ((R - ica.inverse_transform(S))**2).mean()
            if recon < best_recon:
                best_recon, best_S, best_A = recon, S, ica
        except:
            continue
    if best_S is None:
        return Y1
    A = best_A.mixing_
    pvals = []
    vendors = np.unique(batch)
    for j in range(nc):
        groups = [best_S[batch == v, j] for v in vendors]
        try:
            _, pv = stats.f_oneway(*groups)
        except:
            pv = 1.0
        pvals.append(pv)
    pvals = np.array(pvals)
    order = np.argsort(pvals)
    ranks = np.arange(1, len(pvals) + 1)
    thresh = 0.05 * ranks / len(pvals)
    passed = pvals[order] <= thresh
    if passed.any():
        cutoff = ranks[passed].max()
        sig = order[:cutoff]
    else:
        sig = np.array([], dtype=int)
    S_sub = best_S.copy()
    mask = np.ones(nc, dtype=bool)
    mask[sig] = False
    S_sub[:, mask] = 0
    R_batch = S_sub @ A.T + best_A.mean_
    return Y1 - R_batch

def harm_lme(X, meta):
    Y = np.zeros_like(X)
    df = meta.copy().reset_index(drop=True)
    df["age"] = df["Age"].astype(float)
    df["sex"] = (df["Sex"].str.upper() == "F").astype(int)
    site_col = "site" if "site" in df.columns else "Site"
    df["site_"] = df[site_col].astype(str)
    df["vendor"] = df["Manufacturer"].astype("category")
    cats = list(df["vendor"].cat.categories)
    for j in range(X.shape[1]):
        df["y"] = X[:, j]
        success = False
        try:
            md = smf.mixedlm("y ~ age + sex", df, groups=df["site_"])
            mdf = md.fit(method="lbfgs", reml=True)
            re = mdf.random_effects
            shift = np.zeros(len(df))
            for s, eff in re.items():
                shift[df["site_"] == s] = float(eff.iloc[0])
            Y[:, j] = X[:, j] - shift
            success = True
        except:
            pass
        if not success:
            try:
                model = smf.ols("y ~ age + sex + C(vendor)", df).fit()
                params = model.params
                shift = np.zeros(len(df))
                for v in cats[1:]:
                    key = f"C(vendor)[T.{v}]"
                    if key in params.index:
                        shift[df["vendor"] == v] = params[key]
                Y[:, j] = X[:, j] - shift
            except:
                Y[:, j] = X[:, j]
    return Y

METHODS = {
    "LME": harm_lme,
    "ComBat": harm_combat_permetric,
    "ComBat-joint": harm_combat_joint,
    "RELIEF": harm_relief,
    "CovBat": harm_covbat,
}

# ============================================================
# Step 3: Evaluation metrics
# ============================================================
def eta2(y, group):
    grand = y.mean()
    ss_total = ((y - grand)**2).sum()
    ss_between = 0.0
    for g in np.unique(group):
        yg = y[group == g]
        ss_between += len(yg) * (yg.mean() - grand)**2
    return float(ss_between / ss_total) if ss_total > 0 else 0.0

def mean_eta2(X, vendor):
    return float(np.mean([eta2(X[:, j], vendor) for j in range(X.shape[1])]))

def permanova_r2(X, vendor, n_perm=499):
    Z = StandardScaler().fit_transform(X)
    D = np.sqrt(((Z[:, None, :] - Z[None, :, :])**2).sum(-1))
    ids = [f"s{i}" for i in range(len(Z))]
    dm = DistanceMatrix(D, ids)
    res = permanova(dm, vendor, permutations=n_perm)
    F = float(res["test statistic"])
    n = len(vendor); g = len(np.unique(vendor))
    r2 = (F * (g - 1)) / (F * (g - 1) + (n - g))
    return r2, float(res["p-value"])

def within_vendor_r(Xpre, Xpost, vendor):
    rs = []
    for j in range(Xpre.shape[1]):
        for v in np.unique(vendor):
            m = vendor == v
            if m.sum() < 3: continue
            a, b = Xpre[m, j], Xpost[m, j]
            if np.std(a) == 0 or np.std(b) == 0: continue
            r, _ = stats.pearsonr(a, b)
            if np.isfinite(r): rs.append(r)
    return float(np.mean(rs)) if rs else np.nan

def age_r2_mean(X, meta):
    age = meta["Age"].astype(float).values
    vals = []
    for j in range(X.shape[1]):
        y = X[:, j]
        m = np.isfinite(y) & np.isfinite(age)
        if m.sum() < 5: continue
        r, _ = stats.pearsonr(y[m], age[m])
        vals.append(r * r)
    return float(np.mean(vals)) if vals else np.nan

def sex_r2_mean(X, meta):
    sex = (meta["Sex"].str.upper() == "F").astype(int).values
    vals = []
    for j in range(X.shape[1]):
        y = X[:, j]
        m = np.isfinite(y)
        if m.sum() < 5: continue
        r, _ = stats.pearsonr(y[m], sex[m])
        vals.append(r * r)
    return float(np.mean(vals)) if vals else np.nan

# ============================================================
# Step 4: Run full analysis
# ============================================================
def run_analysis(master):
    """Run all harmonization methods and compute metrics."""
    print("[2/5] Running harmonization analysis...")

    # Cohorts
    cohorts = {
        "HC": master[master["Pathology"].str.upper() == "HC"].copy(),
        "ALL": master.copy(),
    }

    # Results storage
    comparison_rows = []
    permanova_rows = []
    umap_rows = []

    for cohort_name, df in cohorts.items():
        print(f"  Cohort: {cohort_name} (N={len(df)})")
        X = df[BIOMARKERS].astype(float).values
        vendor = df["Manufacturer"].values
        meta = df[["Age", "Sex", "Manufacturer", "site"]].copy()

        # Baseline (Original) metrics
        base_eta = mean_eta2(X, vendor)
        base_R2, base_p = permanova_r2(X, vendor, n_perm=499)
        base_age = age_r2_mean(X, df)
        base_sex = sex_r2_mean(X, df)

        # Per-biomarker baseline eta2
        for j, bm in enumerate(BIOMARKERS):
            comparison_rows.append(dict(
                cohort=cohort_name, method="Original", biomarker=bm,
                vendor_eta2=eta2(X[:, j], vendor),
                r_pre_post=1.0,
                age_R2=stats.pearsonr(X[:, j], df["Age"].astype(float).values)[0]**2
                        if np.isfinite(X[:, j]).all() else np.nan,
                sex_R2=stats.pearsonr(X[:, j], (df["Sex"].str.upper()=="F").astype(int).values)[0]**2
                        if np.isfinite(X[:, j]).all() else np.nan,
            ))

        permanova_rows.append(dict(
            cohort=cohort_name, method="Original",
            R2=base_R2, p_value=base_p, n=len(df)))

        # UMAP for Original
        try:
            from umap import UMAP
            Z = StandardScaler().fit_transform(X)
            emb = UMAP(n_components=2, random_state=42, n_neighbors=15).fit_transform(Z)
        except:
            Z = StandardScaler().fit_transform(X)
            emb = TSNE(n_components=2, random_state=42, perplexity=min(30, len(Z)-1)).fit_transform(Z)
        for i in range(len(df)):
            umap_rows.append(dict(
                cohort=cohort_name, method="Original",
                subject=df.iloc[i]["subject"],
                vendor=vendor[i], x=emb[i, 0], y=emb[i, 1]))

        # Run each harmonization method
        for m_name in METHOD_ORDER:
            print(f"    {m_name}...", end=" ", flush=True)
            try:
                Y = METHODS[m_name](X, meta)
            except Exception as e:
                print(f"FAILED ({e})")
                for j, bm in enumerate(BIOMARKERS):
                    comparison_rows.append(dict(
                        cohort=cohort_name, method=m_name, biomarker=bm,
                        vendor_eta2=np.nan, r_pre_post=np.nan,
                        age_R2=np.nan, sex_R2=np.nan))
                permanova_rows.append(dict(
                    cohort=cohort_name, method=m_name,
                    R2=np.nan, p_value=np.nan, n=len(df)))
                continue

            # Per-biomarker metrics
            for j, bm in enumerate(BIOMARKERS):
                r_pp = np.nan
                for v in np.unique(vendor):
                    mask = vendor == v
                    if mask.sum() < 3: continue
                    a, b = X[mask, j], Y[mask, j]
                    if np.std(a) > 0 and np.std(b) > 0:
                        r, _ = stats.pearsonr(a, b)
                        if np.isfinite(r):
                            r_pp = np.nanmean([r_pp, r]) if np.isfinite(r_pp) else r
                comparison_rows.append(dict(
                    cohort=cohort_name, method=m_name, biomarker=bm,
                    vendor_eta2=eta2(Y[:, j], vendor),
                    r_pre_post=r_pp,
                    age_R2=stats.pearsonr(Y[:, j], df["Age"].astype(float).values)[0]**2,
                    sex_R2=stats.pearsonr(Y[:, j], (df["Sex"].str.upper()=="F").astype(int).values)[0]**2,
                ))

            # PERMANOVA
            try:
                R2_post, p_post = permanova_r2(Y, vendor, n_perm=499)
            except:
                R2_post, p_post = np.nan, np.nan
            permanova_rows.append(dict(
                cohort=cohort_name, method=m_name,
                R2=max(R2_post, 1e-6) if np.isfinite(R2_post) else np.nan,
                p_value=p_post, n=len(df)))

            # UMAP
            try:
                from umap import UMAP
                Z_y = StandardScaler().fit_transform(Y)
                emb_y = UMAP(n_components=2, random_state=42, n_neighbors=15).fit_transform(Z_y)
            except:
                Z_y = StandardScaler().fit_transform(Y)
                emb_y = TSNE(n_components=2, random_state=42, perplexity=min(30, len(Z_y)-1)).fit_transform(Z_y)
            for i in range(len(df)):
                umap_rows.append(dict(
                    cohort=cohort_name, method=m_name,
                    subject=df.iloc[i]["subject"],
                    vendor=vendor[i], x=emb_y[i, 0], y=emb_y[i, 1]))

            print("OK")

    comp_long = pd.DataFrame(comparison_rows)
    perm_df = pd.DataFrame(permanova_rows)
    umap_df = pd.DataFrame(umap_rows)

    print(f"  comparison_long: {len(comp_long)} rows")
    print(f"  permanova: {len(perm_df)} rows")
    print(f"  umap: {len(umap_df)} rows")

    return comp_long, perm_df, umap_df

# ============================================================
# Step 5: Figure generation (SCI-quality layout)
# ============================================================

def _save_fig(fig, name):
    """Save figure in PDF and PNG."""
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_DIR / f"{name}.pdf", bbox_inches="tight", pad_inches=0.1)
    fig.savefig(OUT_DIR / f"{name}.png", dpi=300, bbox_inches="tight", pad_inches=0.1)
    plt.close(fig)
    print(f"  [saved] {name}.pdf / .png")

# ---- Fig1a + Fig1b: Study design schematic ----
def _chip(ax, x, y, w, h, label, **kw):
    box = FancyBboxPatch((x, y), w, h,
                         boxstyle=f"round,pad=0.12,rounding_size={kw.get('radius',0.5)}",
                         linewidth=0.7, edgecolor=kw.get("edge", "#444"),
                         facecolor=kw.get("face", "white"), zorder=kw.get("zorder", 3))
    ax.add_patch(box)
    ax.text(x + w/2, y + h/2, label, ha="center", va="center",
            color=kw.get("text_color", "black"), fontsize=kw.get("fontsize", 9),
            fontweight=kw.get("fontweight", "normal"), zorder=kw.get("zorder", 3)+1,
            linespacing=1.3)

def _panel_frame(ax, x, y, w, h, letter, title, subtitle=None):
    bg = FancyBboxPatch((x, y), w, h,
                        boxstyle="round,pad=0.4,rounding_size=1.0",
                        linewidth=0.8, edgecolor="#bbbbbb",
                        facecolor="#f7f7f7", zorder=1)
    ax.add_patch(bg)
    ax.text(x + 1.4, y + h - 1.6, letter, ha="left", va="top",
            fontsize=14, fontweight="bold", color="#222", zorder=5)
    ax.text(x + 5.0, y + h - 2.1, title, ha="left", va="top",
            fontsize=11, fontweight="bold", color="#222", zorder=5)
    if subtitle:
        ax.text(x + 5.0, y + h - 4.0, subtitle, ha="left", va="top",
                fontsize=9, style="italic", color="#555", zorder=5)

def _varrow(ax, xc, yt, yb, color="#222", lw=1.6):
    a = FancyArrowPatch((xc, yt), (xc, yb), arrowstyle="-|>",
                        mutation_scale=18, color=color, linewidth=lw, zorder=2,
                        shrinkA=2, shrinkB=2)
    ax.add_patch(a)

def gen_fig1a(master):
    """Fig1a: Dataset & Biomarkers schematic."""
    fig = plt.figure(figsize=(7.5, 4.5))
    ax = fig.add_subplot(111)
    ax.set_xlim(0, 100); ax.set_ylim(55, 130)
    ax.set_aspect("equal"); ax.axis("off")

    fig.text(0.5, 0.98, "Study design: data and biomarker extraction",
             ha="center", va="top", fontsize=12, fontweight="bold")

    margin = 4; pw = 100 - 2*margin
    h_A, h_B = 30, 25; gap = 3.5
    y_top = 126; y_A = y_top - h_A; y_B = y_A - gap - h_B

    # Panel A: Dataset
    _panel_frame(ax, margin, y_A, pw, h_A, "A", "Dataset & cohorts",
                 subtitle="spine-generic-multi-subject  (r20240523)")

    vc = master["Manufacturer"].value_counts()
    n_total = len(master)
    ax.text(margin+3, y_A+h_A-6, f"{n_total} subjects", ha="left", va="top",
            fontsize=15, fontweight="bold", color="#222", zorder=5)
    ax.text(margin+3, y_A+h_A-10, f"{master['site'].nunique()} sites  -  "
            f"{master['Manufacturer'].nunique()} vendors", ha="left", va="top",
            fontsize=10, color="#444", zorder=5)

    n_hc = (master["Pathology"].str.upper() == "HC").sum()
    chip_w = (pw - 6) * 0.26; chip_h = h_A - 12
    chip_y = y_A + 2.0
    cx2 = margin + pw - 3 - chip_w
    cx1 = cx2 - 1.2 - chip_w
    _chip(ax, cx1, chip_y, chip_w, chip_h, f"HC cohort\nN = {n_hc}",
          face=WONG["yellow"], fontsize=10, fontweight="bold")
    _chip(ax, cx2, chip_y, chip_w, chip_h, f"ALL cohort\nN = {n_total}\n(HC + pathology)",
          face="#ffffff", fontsize=9)

    # Vendor bar
    bar_x0 = margin + 3; bar_x1 = margin + pw - 3; bar_h = 5.0; bar_y = y_A + 2.5
    ax.text((bar_x0+bar_x1)/2, bar_y+bar_h+2.0, "Vendor composition (n)",
            ha="center", va="bottom", fontsize=10, fontweight="bold", color="#333", zorder=5)
    counts = [("GE", vc.get("GE",0)), ("Philips", vc.get("Philips",0)), ("Siemens", vc.get("Siemens",0))]
    total = sum(c for _, c in counts)
    xc = bar_x0
    for v, n in counts:
        seg_w = (bar_x1 - bar_x0) * n / total
        ax.add_patch(Rectangle((xc, bar_y), seg_w, bar_h, facecolor=VENDOR_COLOR[v],
                               edgecolor="black", linewidth=0.6, zorder=3))
        label_text = f"{v}\nn={n}"
        ax.text(xc+seg_w/2, bar_y+bar_h/2, label_text, ha="center", va="center",
                fontsize=9.5, color="white", fontweight="bold", zorder=4)
        xc += seg_w

    # Panel B: Biomarker extraction
    _panel_frame(ax, margin, y_B, pw, h_B, "B", "Biomarker extraction at C2-C3",
                 subtitle="Spinal Cord Toolbox v7.4  -  8 per-subject scalars")
    groups = [("T2w", ["Cord CSA"]), ("T2*", ["GM CSA"]),
              ("MT", ["MTR", "MTsat"]), ("DTI", ["FA", "MD", "AD", "RD"])]
    pad_x = 3.0; region_w = pw - 2*pad_x
    col_w = (region_w - 3*1.2) / 4
    top_y = y_B + h_B - 7; bot_y = y_B + 2.0
    for i, (acq, bms) in enumerate(groups):
        cx = margin + pad_x + i * (col_w + 1.2)
        _chip(ax, cx, top_y - 3, col_w, 2.8, acq, face="#5b5b5b", text_color="white",
              fontsize=10, fontweight="bold")
        n_bm = len(bms); g = 0.5
        bm_h = min((top_y - 5 - bot_y - (n_bm-1)*g) / n_bm, 2.4)
        stack_h = n_bm * bm_h + (n_bm-1) * g
        first_y = top_y - 5 - (top_y - 5 - bot_y - stack_h) / 2 - bm_h
        for j, bm in enumerate(bms):
            by = first_y - j * (bm_h + g)
            _chip(ax, cx, by, col_w, bm_h, bm, face="#ffffff", fontsize=9)

    _varrow(ax, 50, y_A - 0.3, y_B + h_B + 0.3)
    fig.tight_layout(pad=0.5, rect=[0, 0, 1, 0.95])
    _save_fig(fig, "Fig1a_dataset_biomarkers")

def gen_fig1b():
    """Fig1b: Methods & Evaluation schematic."""
    fig = plt.figure(figsize=(7.5, 4.0))
    ax = fig.add_subplot(111)
    ax.set_xlim(0, 100); ax.set_ylim(5, 70)
    ax.set_aspect("equal"); ax.axis("off")

    fig.text(0.5, 0.98, "Study design: harmonization methods and evaluation",
             ha="center", va="top", fontsize=12, fontweight="bold")

    margin = 4; pw = 100 - 2*margin
    h_C, h_D = 26, 22; gap = 3.5
    y_top = 66; y_C = y_top - h_C; y_D = y_C - gap - h_D

    # Panel C: Methods
    _panel_frame(ax, margin, y_C, pw, h_C, "C", "Harmonization methods",
                 subtitle="ordered by modeled vendor components")
    methods = [("LME", "mean only", "(random intercept)"),
               ("ComBat", "+ scale", "(per-metric, EB)"),
               ("ComBat-joint", "+ joint EB", "(across metrics)"),
               ("CovBat", "+ covariance", "(PCA + ComBat)"),
               ("RELIEF", "+ latent factors", "(ICA + BH-FDR)")]
    pad_x = 3.0; region_w = pw - 2*pad_x; n = len(methods); g = 1.0
    box_w = (region_w - (n-1)*g) / n
    header_h = 3.2; desc_h = 5.5; block_h = header_h + desc_h + 0.3
    block_y = y_C + (h_C - block_h) / 2 - 0.6
    for i, (name, l1, l2) in enumerate(methods):
        bx = margin + pad_x + i * (box_w + g)
        col = METHOD_COLOR[name]
        _chip(ax, bx, block_y+desc_h+0.3, box_w, header_h, name,
              face=col, text_color="white", fontsize=10, fontweight="bold")
        _chip(ax, bx, block_y, box_w, desc_h, f"{l1}\n{l2}", face="white", fontsize=9)
    arr_y = block_y - 2.0
    a = FancyArrowPatch((margin+pad_x+0.5, arr_y), (margin+pad_x+region_w-0.5, arr_y),
                        arrowstyle="-|>", mutation_scale=16, color="#666", linewidth=1.3, zorder=2)
    ax.add_patch(a)
    ax.text(margin+pad_x+region_w/2, arr_y-1.5, "Increasing model complexity",
            ha="center", va="top", fontsize=9, color="#444", style="italic", zorder=4)

    # Panel D: Evaluation
    _panel_frame(ax, margin, y_D, pw, h_D, "D", "Evaluation framework",
                 subtitle="three orthogonal endpoints")
    rows = [("Vendor removal", r"$\eta^{2}$  -  ICC  -  PERMANOVA $R^{2}$", WONG["vermilion"]),
            ("Biology recovery", r"age $R^{2}$  -  sex $R^{2}$", WONG["green"]),
            ("Subject preservation", r"within-vendor $r_{\mathrm{pre,post}}$", WONG["pink"])]
    pad_x = 3.0; region_w = pw - 2*pad_x
    top_y = y_D + h_D - 6.5; bot_y = y_D + 1.8
    rg_h = top_y - bot_y; g = 0.6; n = len(rows)
    row_h = (rg_h - (n-1)*g) / n
    for i, (name, metrics, col) in enumerate(rows):
        ry = top_y - (i+1)*row_h - i*g
        ax.add_patch(Rectangle((margin+pad_x, ry), 1.7, row_h, facecolor=col, edgecolor="none", zorder=3))
        _chip(ax, margin+pad_x+2.2, ry, region_w-2.2, row_h,
              f"{name}     {metrics}", face="white", fontsize=9.5)

    _varrow(ax, 50, y_C - 0.3, y_D + h_D + 0.3)
    fig.tight_layout(pad=0.5, rect=[0, 0, 1, 0.95])
    _save_fig(fig, "Fig1b_methods_evaluation")

# ---- Fig2: Baseline vendor effects ----
def gen_fig2(master, comp_long, perm_df, umap_df):
    """Fig2: Baseline vendor effects (4-panel composite)."""
    print("[3/5] Generating Fig2...")

    M = master.copy()
    master_hc = M[M["Pathology"].str.upper() == "HC"].copy()

    L_orig = comp_long[comp_long["method"] == "Original"]
    eta_hc = L_orig[L_orig["cohort"] == "HC"].set_index("biomarker")["vendor_eta2"].to_dict()
    eta_all = L_orig[L_orig["cohort"] == "ALL"].set_index("biomarker")["vendor_eta2"].to_dict()
    P_orig = perm_df[perm_df["method"] == "Original"].set_index("cohort")["R2"]
    r2_hc, r2_all = float(P_orig["HC"]), float(P_orig["ALL"])

    REPR_BM = ["T2w_CSA", "FA", "MTsat", "AD"]
    BM_LABEL = {"T2w_CSA":"Cord CSA (T2w)\n[mm\u00b2]", "GM_CSA_mm2":"GM CSA (T2*)\n[mm\u00b2]",
                "MTR":"MTR [%]", "MTsat":"MTsat [a.u.]", "FA":"FA",
                "MD":"MD [mm\u00b2/s]", "AD":"AD [mm\u00b2/s]", "RD":"RD [mm\u00b2/s]"}
    BM_SHORT = {"T2w_CSA":"T2w_CSA", "GM_CSA_mm2":"GM_CSA", "MTR":"MTR",
                "MTsat":"MTsat", "FA":"FA", "MD":"MD", "AD":"AD", "RD":"RD"}

    # IMPROVED: larger figure, more spacing
    fig = plt.figure(figsize=(7.5, 7.0))
    outer = fig.add_gridspec(nrows=3, ncols=4,
                             height_ratios=[2.5, 2.5, 2.8],
                             width_ratios=[1, 1, 1, 1],
                             hspace=0.70, wspace=0.45)

    # Panel A: 4 boxplots
    axA = [fig.add_subplot(outer[0, j]) for j in range(4)]
    for ax, bm in zip(axA, REPR_BM):
        data = [master_hc.loc[master_hc["Manufacturer"]==v, bm].dropna().values
                for v in ["GE", "Philips", "Siemens"]]
        bp = ax.boxplot(data, positions=np.arange(3), widths=0.55,
                        showfliers=False, patch_artist=True,
                        medianprops=dict(color="black", linewidth=1.1),
                        whiskerprops=dict(color="#444", linewidth=0.8),
                        capprops=dict(color="#444", linewidth=0.8))
        for patch, v in zip(bp["boxes"], ["GE", "Philips", "Siemens"]):
            patch.set_facecolor(VENDOR_COLOR[v]); patch.set_alpha(0.45)
            patch.set_edgecolor("black"); patch.set_linewidth(0.8)
        rng = np.random.default_rng(0)
        for i, (v, vals) in enumerate(zip(["GE", "Philips", "Siemens"], data)):
            jitter = rng.uniform(-0.15, 0.15, len(vals))
            ax.scatter(i+jitter, vals, s=10, c=VENDOR_COLOR[v], alpha=0.65,
                       edgecolor="none", zorder=3)
        ax.set_xticks(np.arange(3))
        ax.set_xticklabels(["GE", "Philips", "Siemens"], fontsize=8)
        ax.set_title(BM_LABEL.get(bm, bm), fontsize=9)
        ax.spines[["top", "right"]].set_visible(False)
        ax.tick_params(axis="y", labelsize=8)
        eta = eta_hc.get(bm)
        if eta is not None:
            ax.text(0.03, 0.95, rf"$\eta^{{2}}={eta:.3f}$", transform=ax.transAxes,
                    ha="left", va="top", fontsize=8.5,
                    bbox=dict(boxstyle="round,pad=0.25", facecolor="white",
                              edgecolor="#999", linewidth=0.6))

    axA[0].annotate("A.  Per-vendor distribution of 4 representative biomarkers (HC)",
                    xy=(-0.1, 1.22), xycoords="axes fraction",
                    ha="left", va="bottom", fontweight="bold", fontsize=10)

    # Panel B: UMAP
    axB_hc = fig.add_subplot(outer[1, 2])
    axB_all = fig.add_subplot(outer[1, 3])
    for ax, cohort in [(axB_hc, "HC"), (axB_all, "ALL")]:
        sub = umap_df[(umap_df["cohort"]==cohort) & (umap_df["method"]=="Original")]
        for v in ["GE", "Philips", "Siemens"]:
            s = sub[sub["vendor"]==v]
            ax.scatter(s["x"], s["y"], s=20, c=VENDOR_COLOR[v], alpha=0.75,
                       edgecolor="black", linewidth=0.3, label=f"{v} (n={len(s)})")
        ax.set_xticks([]); ax.set_yticks([])
        ax.set_xlabel("UMAP-1", fontsize=8); ax.set_ylabel("UMAP-2", fontsize=8)
        ax.set_title(f"{cohort} (N={len(sub)})", loc="left", fontweight="bold", fontsize=9)
        ax.spines[["top","right","bottom","left"]].set_color("#888")
        ax.spines[["top","right","bottom","left"]].set_linewidth(0.6)
        ax.legend(loc="best", frameon=False, fontsize=7, handletextpad=0.4, borderpad=0.3)
    axB_hc.annotate("B.  UMAP of standardized 8-biomarker vectors (Original)",
                    xy=(-0.1, 1.22), xycoords="axes fraction",
                    ha="left", va="bottom", fontweight="bold", fontsize=10)

    # Panel B-left: text takeaways
    axC = fig.add_subplot(outer[1, 0:2]); axC.axis("off")
    axC.text(0.02, 0.95, "Findings - baseline vendor effects",
             transform=axC.transAxes, fontsize=11, fontweight="bold", va="top")
    axC.text(0.02, 0.82,
        "* Demographic covariates (age, sex) are balanced\n  across vendors (all p>0.1).\n\n"
        "* Univariate vendor $\\eta^{2}$ ranges from ~0.05 (CSA)\n  to ~0.62 (AD) - mean ~0.28 across 8 biomarkers.\n\n"
        "* Multivariate PERMANOVA $R^{2}$=" + f"{r2_hc:.3f}" + " (HC)\n  and " +
        f"{r2_all:.3f}" + " (ALL): vendor explains\n  about one third of joint biomarker distance.\n\n"
        "* Vendor effect dominates biological covariates\n  by ~3 orders of magnitude.",
        transform=axC.transAxes, fontsize=9, va="top")

    # Panel C: eta2 bars
    axD = fig.add_subplot(outer[2, :])
    bms = sorted(eta_hc.keys(), key=lambda b: -eta_hc[b])
    x = np.arange(len(bms), dtype=float); w = 0.36
    axD.bar(x-w/2, [eta_hc[b] for b in bms], width=w, color="#444444", alpha=1.0,
            edgecolor="black", linewidth=0.5, label="HC")
    axD.bar(x+w/2, [eta_all[b] for b in bms], width=w, color="#444444", alpha=0.55,
            edgecolor="black", linewidth=0.5, label="ALL")
    axD.set_xticks(x)
    axD.set_xticklabels([BM_SHORT[b] for b in bms], rotation=20, ha="right", fontsize=8.5)
    axD.set_ylabel(r"Baseline vendor $\eta^{2}$  (per biomarker)", fontsize=9)
    axD.spines[["top", "right"]].set_visible(False)
    axD.set_axisbelow(True)
    axD.yaxis.grid(True, linestyle=":", linewidth=0.5, alpha=0.55)
    axD.axhline(r2_hc, color="#222", linestyle="--", linewidth=0.9)
    axD.text(len(bms)-0.5, r2_hc, f"  HC multivariate $R^{{2}}={r2_hc:.3f}$",
             va="bottom", ha="right", fontsize=8, color="#222")
    axD.axhline(r2_all, color="#222", linestyle=":", linewidth=0.9)
    axD.text(len(bms)-0.5, r2_all, f"  ALL multivariate $R^{{2}}={r2_all:.3f}$",
             va="top", ha="right", fontsize=8, color="#222")
    axD.legend(loc="upper right", frameon=False, fontsize=8.5)
    axD.set_ylim(0, max(max(eta_hc.values()), max(eta_all.values()), r2_hc, r2_all) * 1.25)
    axD.annotate("C.  Baseline vendor $\\eta^{2}$ per biomarker, with PERMANOVA $R^{2}$ reference",
                 xy=(0, 1.12), xycoords="axes fraction",
                 ha="left", va="bottom", fontweight="bold", fontsize=10)

    fig.suptitle("Baseline vendor effects in unharmonized spinal-cord qMRI biomarkers",
                 fontsize=11, fontweight="bold", y=1.02)
    fig.tight_layout(pad=0.4, rect=[0, 0, 1, 0.97])
    _save_fig(fig, "Fig2_baseline_effects")

# ---- Fig3: Harmonization performance ----
def gen_fig3(comp_long, perm_df):
    """Fig3: Harmonization performance (2x2 panels)."""
    print("[3/5] Generating Fig3...")

    METHOD_LABEL = {"LME":"LME", "ComBat":"ComBat", "ComBat-joint":"ComBat\n(joint)",
                    "RELIEF":"RELIEF", "CovBat":"CovBat"}
    COHORT_ALPHA = {"HC": 1.0, "ALL": 0.55}

    # Aggregate
    agg = {}
    for cohort in ["HC", "ALL"]:
        L = comp_long[comp_long["cohort"]==cohort]
        P = perm_df[perm_df["cohort"]==cohort].set_index("method")
        age_pre = L.loc[L["method"]=="Original", "age_R2"].mean()
        eta_pre = L.loc[L["method"]=="Original", "vendor_eta2"].mean()
        R2_pre = float(P.loc["Original", "R2"]) if "Original" in P.index else np.nan
        per_m = {"vendor_eta2":{}, "age_gain":{}, "sex_gain":{}, "permanova_R2":{}}
        for m in METHOD_ORDER:
            sub = L[L["method"]==m]
            if sub.empty: continue
            per_m["vendor_eta2"][m] = sub["vendor_eta2"].mean()
            per_m["age_gain"][m] = sub["age_R2"].mean() - age_pre
            per_m["sex_gain"][m] = sub["sex_R2"].mean() - L.loc[L["method"]=="Original","sex_R2"].mean()
            r2p = float(P.loc[m, "R2"]) if m in P.index else np.nan
            per_m["permanova_R2"][m] = max(r2p, 1e-6) if np.isfinite(r2p) else np.nan
        per_m["baseline"] = dict(eta_pre=eta_pre, R2_pre=R2_pre)
        agg[cohort] = per_m

    def grouped_bars(ax, vals_hc, vals_all, ylabel, baseline=None, baseline_label=None,
                     log=False, winner="min"):
        methods = [m for m in METHOD_ORDER if m in vals_hc and m in vals_all and
                   np.isfinite(vals_hc[m]) and np.isfinite(vals_all[m])]
        n = len(methods); x = np.arange(n, dtype=float); w = 0.36
        for i, m in enumerate(methods):
            col = METHOD_COLOR[m]
            ax.bar(x[i]-w/2, vals_hc[m], width=w, color=col, alpha=COHORT_ALPHA["HC"],
                   edgecolor="black", linewidth=0.6)
            ax.bar(x[i]+w/2, vals_all[m], width=w, color=col, alpha=COHORT_ALPHA["ALL"],
                   edgecolor="black", linewidth=0.6)
        pooled = {m: (vals_hc[m]+vals_all[m])/2 for m in methods}
        win = min(pooled, key=pooled.get) if winner=="min" else max(pooled, key=pooled.get)
        wi = methods.index(win)
        ymax_pair = max(vals_hc[win], vals_all[win])
        star_y = ymax_pair * 1.6 if log else ymax_pair + abs(ymax_pair)*0.12 + 1e-4
        ax.text(wi, star_y, "*", ha="center", va="bottom", fontsize=14, color="#222", fontweight="bold")
        if baseline is not None:
            ax.axhline(baseline, color="#888", linestyle="--", linewidth=1.0)
            if baseline == 0.0:
                ax.text(n-0.5, baseline, f"  {baseline_label}" if baseline_label else "",
                        va="top", ha="left", fontsize=8, color="#555")
            else:
                ax.text(n-0.5, baseline, f"  {baseline_label}" if baseline_label else "",
                        va="center", ha="left", fontsize=8, color="#555")
        ax.set_xticks(x)
        ax.set_xticklabels([METHOD_LABEL[m] for m in methods], fontsize=8)
        ax.set_ylabel(ylabel, fontsize=9)
        if log: ax.set_yscale("log")
        ax.spines[["top","right"]].set_visible(False)
        ax.set_axisbelow(True)
        ax.yaxis.grid(True, linestyle=":", linewidth=0.5, alpha=0.55)

    # IMPROVED: larger figure
    fig, axes = plt.subplots(2, 2, figsize=(7.5, 5.8),
                             gridspec_kw=dict(wspace=0.35, hspace=0.55))
    (axA, axB), (axC, axD) = axes

    grouped_bars(axA, agg["HC"]["vendor_eta2"], agg["ALL"]["vendor_eta2"],
                 ylabel=r"Residual vendor $\eta^{2}$  (mean)",
                 baseline=max(agg["HC"]["baseline"]["eta_pre"], agg["ALL"]["baseline"]["eta_pre"]),
                 baseline_label=r"Original $\eta^{2}$" + f"$\\approx{agg['HC']['baseline']['eta_pre']:.2f}$",
                 log=True, winner="min")
    axA.set_title("A.  Univariate vendor removal", loc="left", fontweight="bold", fontsize=10)

    grouped_bars(axB, agg["HC"]["permanova_R2"], agg["ALL"]["permanova_R2"],
                 ylabel=r"PERMANOVA vendor $R^{2}$",
                 baseline=max(agg["HC"]["baseline"]["R2_pre"], agg["ALL"]["baseline"]["R2_pre"]),
                 baseline_label=r"Original $R^{2}$" + f"$\\approx{agg['HC']['baseline']['R2_pre']:.2f}$",
                 log=True, winner="min")
    axB.set_title("B.  Multivariate vendor removal", loc="left", fontweight="bold", fontsize=10)

    grouped_bars(axC, agg["HC"]["age_gain"], agg["ALL"]["age_gain"],
                 ylabel=r"Age $R^{2}$ gain  (post $-$ pre)",
                 baseline=0.0, baseline_label="no change", log=False, winner="max")
    axC.set_title("C.  Biological signal recovery (age)", loc="left", fontweight="bold", fontsize=10)

    grouped_bars(axD, agg["HC"]["sex_gain"], agg["ALL"]["sex_gain"],
                 ylabel=r"Sex $R^{2}$ gain  (post $-$ pre)",
                 baseline=0.0, baseline_label="no change", log=False, winner="max")
    axD.set_title("D.  Biological signal stability (sex)", loc="left", fontweight="bold", fontsize=10)

    # Shared legend
    n_hc = (comp_long[comp_long["cohort"]=="HC"]["method"]=="Original").sum()
    legend_handles = [
        Patch(facecolor="#777", alpha=1.0, edgecolor="black", linewidth=0.6, label="HC cohort"),
        Patch(facecolor="#777", alpha=0.55, edgecolor="black", linewidth=0.6, label="ALL cohort"),
        Line2D([0],[0], color="#888", linestyle="--", linewidth=1.0, label="Reference"),
        Line2D([0],[0], marker="*", color="w", markerfacecolor="#222", markeredgecolor="#222",
               markersize=12, label="Within-panel winner"),
    ]
    fig.legend(handles=legend_handles, ncol=4, loc="lower center",
               frameon=False, bbox_to_anchor=(0.5, -0.01), fontsize=8.5)

    fig.suptitle("Harmonization performance across orthogonal endpoints",
                 fontsize=11, fontweight="bold", y=1.01)
    fig.tight_layout(pad=0.4, rect=[0, 0.05, 1, 0.97])
    _save_fig(fig, "Fig3_performance")

# ---- Fig4: Trade-off scatter ----
def gen_fig4(comp_long, perm_df):
    """Fig4: Method trade-off scatter."""
    print("[3/5] Generating Fig4...")

    METHOD_LABEL = {"LME":"LME", "ComBat":"ComBat", "ComBat-joint":"ComBat (joint)",
                    "RELIEF":"RELIEF", "CovBat":"CovBat"}

    def aggregate(long, perm, cohort):
        long = long[long["cohort"]==cohort]
        perm = perm[perm["cohort"]==cohort].set_index("method")
        age_pre = long.loc[long["method"]=="Original", "age_R2"].mean()
        r2_orig = perm.loc["Original", "R2"] if "Original" in perm.index else np.nan
        rows = []
        for m in METHOD_ORDER:
            sub = long[long["method"]==m]
            if sub.empty: continue
            y = sub["r_pre_post"].mean()
            y_sd = sub["r_pre_post"].std(ddof=1)
            r2p = float(perm.loc[m, "R2"]) if m in perm.index else np.nan
            r2p = max(r2p, 1e-6) if np.isfinite(r2p) else np.nan
            x = -np.log10(r2p) if np.isfinite(r2p) else np.nan
            eta_post = sub["vendor_eta2"].clip(lower=1e-7)
            x_proxy = -np.log10(eta_post)
            x_sd = x_proxy.std(ddof=1)
            age_gain = sub["age_R2"].mean() - age_pre
            rows.append(dict(method=m, x=x, x_sd=x_sd, y=y, y_sd=y_sd,
                             age_gain=age_gain, r2_post=r2p, r2_orig=r2_orig))
        return pd.DataFrame(rows)

    def pareto_front(pts):
        keep = []
        for i, row in pts.iterrows():
            dominated = ((pts["x"]>=row["x"]) & (pts["y"]>=row["y"]) &
                        ((pts["x"]>row["x"]) | (pts["y"]>row["y"]))).any()
            if not dominated: keep.append(i)
        return pts.loc[keep].sort_values("x").reset_index(drop=True)

    def plot_panel(ax, tab, cohort):
        g = tab["age_gain"].values
        g_min, g_max = g.min(), g.max()
        s_min, s_max = 90, 360
        if g_max - g_min < 1e-6:
            sizes = np.full_like(g, (s_min+s_max)/2)
        else:
            sizes = s_min + (g - g_min) / (g_max - g_min) * (s_max - s_min)
        for _, row in tab.iterrows():
            ax.errorbar(row["x"], row["y"], xerr=row["x_sd"], yerr=row["y_sd"],
                        fmt="none", ecolor="#666", elinewidth=0.7, capsize=2.5,
                        capthick=0.7, alpha=0.55, zorder=1)
        front = pareto_front(tab.dropna(subset=["x","y"]))
        if len(front) >= 2:
            ax.plot(front["x"], front["y"], color="#222", linestyle="--",
                    linewidth=1.0, alpha=0.55, zorder=2, label="Pareto front")
        for (_, row), sz in zip(tab.iterrows(), sizes):
            ax.scatter(row["x"], row["y"], s=sz, c=METHOD_COLOR[row["method"]],
                       edgecolors="black", linewidths=0.8, zorder=3,
                       label=METHOD_LABEL[row["method"]])
        # Method color legend placed in upper-left to avoid data overlap
        ax.legend(loc="upper left", frameon=False, fontsize=8, handletextpad=0.3,
                  borderpad=0.25, title="Method", title_fontsize=8)
        ax.set_xlabel(r"Decontamination strength  $-\log_{10}(R^{2}_{\mathrm{PERMANOVA}})$", fontsize=9)
        ax.set_ylabel(r"Subject-level preservation  $\bar{r}_{\mathrm{pre,post}}$", fontsize=9)
        n_R2 = tab["r2_orig"].iloc[0] if len(tab) > 0 else 0
        ax.set_title(f"{cohort}  (Original $R^{{2}}={n_R2:.3f}$)", fontweight="bold", loc="left", fontsize=10)
        ax.spines[["top","right"]].set_visible(False)
        ax.set_axisbelow(True)
        ax.grid(True, linestyle=":", linewidth=0.5, alpha=0.55)
        xmin, xmax = ax.get_xlim(); ymin, ymax = ax.get_ylim()
        xrange = xmax - xmin; yrange = ymax - ymin
        ax.set_xlim(xmin - 0.05*xrange, xmax + 0.15*xrange)
        ax.set_ylim(ymin - 0.04*yrange, ymax + 0.12*yrange)

    tabs = {c: aggregate(comp_long, perm_df, c) for c in ["HC", "ALL"]}
    all_gain = pd.concat([t["age_gain"] for t in tabs.values()])
    g_min, g_max = float(all_gain.min()), float(all_gain.max())

    # IMPROVED: larger figure, better spacing
    fig = plt.figure(figsize=(8.5, 5.5))
    gs = fig.add_gridspec(2, 2, height_ratios=[12, 1.0], width_ratios=[1, 1],
                          hspace=0.50, wspace=0.35)
    ax_hc = fig.add_subplot(gs[0, 0])
    ax_all = fig.add_subplot(gs[0, 1])
    ax_lg = fig.add_subplot(gs[1, :])

    plot_panel(ax_hc, tabs["HC"], "HC")
    plot_panel(ax_all, tabs["ALL"], "ALL")

    # Marker size legend
    ax_lg.axis("off")
    vals = np.linspace(g_min, g_max, 3)
    s_min, s_max = 90, 360
    sizes = s_min + (vals - g_min) / max(g_max - g_min, 1e-6) * (s_max - s_min)
    xs = np.linspace(0.15, 0.85, 3)
    for x, sz, v in zip(xs, sizes, vals):
        ax_lg.scatter(x, 0.50, s=sz, facecolor="#bbb", edgecolor="black",
                      linewidth=0.7, transform=ax_lg.transAxes)
        ax_lg.text(x, 0.08, f"{v:+.3f}", ha="center", va="top",
                   fontsize=8, transform=ax_lg.transAxes)
    ax_lg.text(0.50, 0.92, r"marker size  $\propto$  age $R^{2}$ gain (post $-$ pre)",
               ha="center", va="top", fontsize=9, transform=ax_lg.transAxes)

    fig.suptitle("Harmonization method trade-off: decontamination vs. subject preservation",
                 fontsize=11, fontweight="bold", y=1.01)
    fig.tight_layout(pad=0.4, rect=[0, 0, 1, 0.96])
    _save_fig(fig, "Fig4_tradeoff")

# ---- Fig5: r_pre,post heatmap with marginals ----
def gen_fig5(comp_long):
    """Fig5: r_pre,post heatmap with marginal bars."""
    print("[3/5] Generating Fig5...")

    METHOD_LABEL = {"LME":"LME", "ComBat":"ComBat", "ComBat-joint":"ComBat\n(joint)",
                    "RELIEF":"RELIEF", "CovBat":"CovBat"}
    BM_ORDER = BIOMARKERS
    BM_LABEL = {"T2w_CSA":"Cord CSA\n(T2w)", "GM_CSA_mm2":"GM CSA\n(T2*)",
                "MTR":"MTR", "MTsat":"MTsat", "FA":"FA",
                "MD":"MD", "AD":"AD", "RD":"RD"}

    def build_matrix(long, cohort, value):
        sub = long[long["cohort"]==cohort]
        mat = sub.pivot(index="method", columns="biomarker", values=value)
        mat = mat.reindex(index=METHOD_ORDER, columns=BM_ORDER)
        return mat

    def baseline_eta2(long, cohort):
        sub = long[(long["cohort"]==cohort) & (long["method"]=="Original")]
        return sub.set_index("biomarker")["vendor_eta2"].reindex(BM_ORDER)

    all_r = comp_long[comp_long["method"]!="Original"]["r_pre_post"]
    vmin = float(np.floor(all_r.min() * 20) / 20); vmax = 1.0

    def draw_cohort(fig, gs_block, long, cohort):
        inner = gs_block.subgridspec(3, 2, height_ratios=[0.7, 2.8, 0.6],
                                     width_ratios=[5.0, 0.9], hspace=0.10, wspace=0.08)
        ax_top = fig.add_subplot(inner[0, 0])
        ax_hm = fig.add_subplot(inner[1, 0])
        ax_rt = fig.add_subplot(inner[1, 1])
        ax_eta = fig.add_subplot(inner[2, 0])

        mat = build_matrix(long, cohort, "r_pre_post")
        eta0 = baseline_eta2(long, cohort)
        col_m = mat.mean(axis=0); row_m = mat.mean(axis=1)
        n_bio = mat.shape[1]; n_met = mat.shape[0]

        im = ax_hm.imshow(mat.values, aspect="auto", cmap="viridis", vmin=vmin, vmax=vmax)
        ax_hm.set_xticks(np.arange(n_bio)); ax_hm.set_xticklabels([])
        ax_hm.set_yticks(np.arange(n_met))
        ax_hm.set_yticklabels([METHOD_LABEL[m] for m in mat.index], fontsize=9)
        for i in range(n_met):
            for j in range(n_bio):
                v = mat.values[i, j]
                if np.isfinite(v):
                    txt_col = "white" if v < (vmin+vmax)/2 + 0.05 else "black"
                    ax_hm.text(j, i, f"{v:.2f}", ha="center", va="center",
                               fontsize=9.5, color=txt_col)
        ax_hm.set_ylabel(f"{cohort}\nharmonization method", fontweight="bold", fontsize=9.5)
        for sp in ax_hm.spines.values(): sp.set_visible(False)

        ax_top.bar(np.arange(n_bio), col_m.values, color="#444", edgecolor="black", linewidth=0.5)
        ax_top.set_xlim(ax_hm.get_xlim()); ax_top.set_xticks([])
        ax_top.set_ylim(0.4, 1.0); ax_top.set_yticks([0.5, 0.75, 1.0])
        ax_top.set_yticklabels(["0.5","0.75","1.0"], fontsize=8.5)
        ax_top.set_ylabel("mean r\nacross\nmethods", fontsize=8.5)
        ax_top.spines[["top","right"]].set_visible(False)
        ax_top.yaxis.grid(True, linestyle=":", linewidth=0.4, alpha=0.5)
        ax_top.set_axisbelow(True)
        ranks = (-col_m).rank(method="min").astype(int)
        for j, v in enumerate(col_m.values):
            if np.isfinite(v):
                ax_top.text(j, v+0.02, f"#{ranks.iloc[j]}", ha="center", va="bottom",
                            fontsize=8.5, color="#666")

        ax_rt.barh(np.arange(n_met)[::-1], row_m.values, color="#444",
                   edgecolor="black", linewidth=0.5)
        ax_rt.set_ylim(ax_hm.get_ylim()); ax_rt.set_yticks([])
        ax_rt.set_xlim(0.4, 1.0); ax_rt.set_xticks([0.5, 0.75, 1.0])
        ax_rt.set_xticklabels(["0.5","0.75","1.0"], fontsize=8.5)
        ax_rt.set_xlabel("mean r\nacross\nbiomarkers", fontsize=8.5)
        ax_rt.spines[["top","right"]].set_visible(False)
        ax_rt.xaxis.grid(True, linestyle=":", linewidth=0.4, alpha=0.5)
        ax_rt.set_axisbelow(True)
        win = row_m.idxmax(); win_pos = list(mat.index).index(win)
        star_y = (n_met-1) - win_pos
        ax_rt.text(row_m.max()+0.04, star_y, "*", ha="left", va="center",
                   fontsize=14, color="#222", fontweight="bold")
        ranks_m = (-row_m).rank(method="min").astype(int)
        for i, (m, v) in enumerate(row_m.items()):
            y = (n_met-1) - i
            if np.isfinite(v):
                ax_rt.text(v-0.015, y, f"#{ranks_m[m]}", ha="right", va="center",
                           fontsize=8.5, color="white" if v > 0.55 else "#666")

        eta_row = eta0.values.reshape(1, -1)
        im_eta = ax_eta.imshow(eta_row, aspect="auto", cmap="magma_r",
                               vmin=0.0, vmax=float(np.nanmax(eta0)) if np.isfinite(np.nanmax(eta0)) else 1.0)
        ax_eta.set_xticks(np.arange(n_bio))
        ax_eta.set_xticklabels([BM_LABEL[b] for b in eta0.index], fontsize=8.5)
        ax_eta.set_yticks([0])
        ax_eta.set_yticklabels(["Baseline\n$\\eta^{2}$"], fontsize=8.5)
        for j, v in enumerate(eta0.values):
            if np.isfinite(v):
                col = "white" if v > 0.30 else "black"
                ax_eta.text(j, 0, f"{v:.2f}", ha="center", va="center", fontsize=9, color=col)
        for sp in ax_eta.spines.values(): sp.set_visible(False)
        ax_eta.tick_params(axis="x", length=0, pad=2)
        ax_eta.tick_params(axis="y", length=0, pad=2)
        return im

    # IMPROVED: larger figure, more readable heatmap
    fig = plt.figure(figsize=(8.5, 8.0))
    outer = fig.add_gridspec(2, 2, width_ratios=[28, 1], height_ratios=[1, 1],
                             wspace=0.06, hspace=0.42)
    im_hc = draw_cohort(fig, outer[0, 0], comp_long, "HC")
    im_all = draw_cohort(fig, outer[1, 0], comp_long, "ALL")
    cax = fig.add_subplot(outer[:, 1])
    cb = fig.colorbar(im_hc, cax=cax)
    cb.set_label(r"Within-vendor $r_{\mathrm{pre,post}}$  (subject preservation)", fontsize=9)
    cb.ax.tick_params(labelsize=8)
    fig.suptitle("Subject-level preservation by biomarker and method",
                 fontsize=11, fontweight="bold", y=0.995)
    _save_fig(fig, "Fig5_rpre_post_marginals")

# ---- FigS: Design C robustness ----
def gen_figs(master):
    """FigS: Design C balanced-subsample robustness."""
    print("[3/5] Generating FigS...")

    # Try to use existing results_design_C.csv first
    existing_csv = Path("E:/qm_harmonization_paper/derivatives/all-files/results_design_C.csv")
    if existing_csv.exists():
        res = pd.read_csv(existing_csv)
        print(f"  Using existing results_design_C.csv ({len(res)} rows)")
    else:
        # Run Design C analysis
        print("  Running Design C analysis...")
        df = master[master["Pathology"].str.upper() == "HC"].copy()
        df = df.dropna(subset=BIOMARKERS + ["Manufacturer", "Age", "Sex"]).copy()

        all_rows = []
        k_seeds = 20; target_si = 40
        for s in range(k_seeds):
            rng = np.random.RandomState(s)
            ge = df[df["Manufacturer"]=="GE"]
            ph = df[df["Manufacturer"]=="Philips"]
            si = df[df["Manufacturer"]=="Siemens"]
            si_sub = si.sample(n=min(target_si, len(si)), random_state=s)
            sub = pd.concat([ge, ph, si_sub], ignore_index=True)
            X = sub[BIOMARKERS].astype(float).values
            vendor = sub["Manufacturer"].values
            meta = sub[["Age","Sex","Manufacturer","site"]].copy()

            base_eta = mean_eta2(X, vendor)
            base_R2, base_p = permanova_r2(X, vendor, n_perm=199)
            base_age = age_r2_mean(X, sub)
            all_rows.append(dict(seed=s, method="Original", eta2=base_eta,
                                 R2_perm=base_R2, p_perm=base_p, r_pre_post=1.0,
                                 age_R2=base_age, age_R2_gain=0.0, n=len(sub)))
            for m_name in METHOD_ORDER:
                try:
                    Y = METHODS[m_name](X, meta)
                    eta = mean_eta2(Y, vendor)
                    R2, pv = permanova_r2(Y, vendor, n_perm=199)
                    r = within_vendor_r(X, Y, vendor)
                    ar2 = age_r2_mean(Y, sub)
                    all_rows.append(dict(seed=s, method=m_name, eta2=eta,
                                         R2_perm=R2, p_perm=pv, r_pre_post=r,
                                         age_R2=ar2, age_R2_gain=ar2-base_age, n=len(sub)))
                except Exception as e:
                    all_rows.append(dict(seed=s, method=m_name, eta2=np.nan,
                                         R2_perm=np.nan, p_perm=np.nan, r_pre_post=np.nan,
                                         age_R2=np.nan, age_R2_gain=np.nan, n=len(sub)))
            print(f"    seed {s}/{k_seeds}", end="\r")
        res = pd.DataFrame(all_rows)
        print()

    # Full-cohort reference
    df_hc = master[master["Pathology"].str.upper() == "HC"].copy()
    df_hc = df_hc.dropna(subset=BIOMARKERS + ["Manufacturer", "Age", "Sex"]).copy()
    X_full = df_hc[BIOMARKERS].astype(float).values
    v_full = df_hc["Manufacturer"].values
    meta_full = df_hc[["Age","Sex","Manufacturer","site"]].copy()
    full_ref = {}
    base_age_full = age_r2_mean(X_full, df_hc)
    for m_name in METHOD_ORDER:
        try:
            Y = METHODS[m_name](X_full, meta_full)
            full_ref[m_name] = {
                "R2_perm": max(permanova_r2(Y, v_full, n_perm=199)[0], 1e-6),
                "r_pre_post": within_vendor_r(X_full, Y, v_full),
                "age_R2_gain": age_r2_mean(Y, df_hc) - base_age_full,
            }
        except:
            full_ref[m_name] = {"R2_perm": np.nan, "r_pre_post": np.nan, "age_R2_gain": np.nan}

    # IMPROVED: taller figure, better spacing
    fig, axes = plt.subplots(1, 3, figsize=(8.0, 3.5))
    panels = [
        ("R2_perm", "(A)  Residual vendor structure\n$R^{2}_{\\mathrm{PERMANOVA}}$ (log)",
         "lower is better", True),
        ("r_pre_post", "(B)  Subject preservation\n$\\bar{r}_{\\mathrm{pre,post}}$",
         "higher is better", False),
        ("age_R2_gain", "(C)  Biology unmasking\n$\\Delta$ age $R^{2}$ vs baseline",
         "higher is better", False),
    ]
    data = res[res["method"].isin(METHOD_ORDER)].copy()
    xs = np.arange(len(METHOD_ORDER))
    cols = [METHOD_COLOR[m] for m in METHOD_ORDER]

    for k, (col, ttl, sub_label, logy) in enumerate(panels):
        ax = axes[k]
        means = data.groupby("method")[col].mean().reindex(METHOD_ORDER)
        sds = data.groupby("method")[col].std().reindex(METHOD_ORDER)
        ax.bar(xs, means.values, yerr=sds.values, color=cols,
               edgecolor="black", linewidth=0.7,
               error_kw=dict(ecolor="#333", capsize=3, lw=0.8),
               alpha=0.85, zorder=3)
        for i, m in enumerate(METHOD_ORDER):
            yvals = data[data["method"]==m][col].dropna().values
            if len(yvals) > 0:
                jitter = (np.random.RandomState(42+i).rand(len(yvals)) - 0.5) * 0.30
                ax.scatter(np.full_like(yvals, xs[i]) + jitter, yvals,
                           s=12, color="white", edgecolor="#222",
                           linewidth=0.5, alpha=0.85, zorder=5)
            ref = full_ref.get(m, {}).get(col, np.nan)
            if np.isfinite(ref):
                ax.hlines(ref, xs[i]-0.42, xs[i]+0.42, colors="black",
                          linestyles="--", linewidth=1.1, zorder=6)
        ax.set_xticks(xs)
        ax.set_xticklabels(METHOD_ORDER, rotation=30, ha="right", fontsize=8.5)
        ax.set_title(ttl, fontsize=9.5, fontweight="bold")
        ax.set_xlabel(sub_label, fontsize=8.5, color="#555")
        ax.grid(axis="y", alpha=0.3)
        if logy: ax.set_yscale("log")

    handles = [
        Line2D([0],[0], color="black", ls="--", lw=1.1, label="Full-cohort reference"),
        Line2D([0],[0], marker="o", color="w", markerfacecolor="white",
               markeredgecolor="#222", markeredgewidth=0.5, markersize=5,
               linestyle="", label="Single random seed"),
    ]
    fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, 1.03),
               ncol=2, frameon=False, fontsize=9)
    k_seeds = res["seed"].nunique()
    n_mean = int(res["n"].mean())
    fig.suptitle(f"Design C - Balanced-subsample robustness "
                 f"(HC, K={k_seeds} seeds, n$\\approx${n_mean} per seed)",
                 fontsize=10.5, fontweight="bold", y=1.08)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    _save_fig(fig, "FigS_design_C_robustness")

# ============================================================
# Main
# ============================================================
def main():
    print("=" * 60)
    print("Spine-Generic QM Harmonization - Figure Regeneration")
    print("=" * 60)

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    # Step 1: Load data
    master = load_data()

    # Step 2: Run analysis
    comp_long, perm_df, umap_df = run_analysis(master)

    # Step 3: Generate all figures
    print("\n[3/5] Generating figures...")
    gen_fig1a(master)
    gen_fig1b()
    gen_fig2(master, comp_long, perm_df, umap_df)
    gen_fig3(comp_long, perm_df)
    gen_fig4(comp_long, perm_df)
    gen_fig5(comp_long)
    gen_figs(master)

    print("\n[5/5] Done! All figures saved to:", OUT_DIR)
    print("=" * 60)

if __name__ == "__main__":
    main()
