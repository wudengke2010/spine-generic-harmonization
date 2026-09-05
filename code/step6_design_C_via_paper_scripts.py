"""
Step 6.C v2 — Design C: Balanced-subsample robustness check
            (wrapper around the paper's own step3_*.py scripts)
============================================================
This wrapper guarantees that the harmonization implementations
exactly match the paper (no re-implementation drift).  For each
random seed, it:
    1. Reads biomarkers_master.csv and restricts to the chosen cohort
    2. Subsamples Siemens down to `--siemens` subjects, keeping all
       GE and Philips → balanced subset with n ≈ 25 + 33 + 40 ≈ 98 (HC)
    3. Writes the subset to a tmp CSV
    4. Calls each of the 5 paper step3_*.py scripts via subprocess
       with that tmp CSV as --master  (each script re-runs the
       paper-canonical implementation on the subsample)
    5. Reads each harmonized CSV, extracts the per-method columns,
       and computes:
          - mean vendor η² across 8 biomarkers
          - PERMANOVA R² (vendor)
          - within-vendor r_{pre,post}
          - mean age semi-partial R² gain vs baseline
Aggregate over K seeds → mean ± SD.

Output (--out dir):
  results_design_C.csv
  summary_design_C.csv
  FigS_design_C_robustness.png/pdf

Local run:
  python3 step6_design_C_via_paper_scripts.py \
      --master  ~/qm_harmonization_paper/results/biomarkers_master.csv \
      --code    ~/qm_harmonization_paper/code \
      --out     ~/qm_harmonization_paper/paper \
      --cohort  HC --k 20 --siemens 40
"""

from __future__ import annotations
import argparse
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy import stats
from sklearn.preprocessing import StandardScaler
from skbio.stats.distance import permanova, DistanceMatrix

# ---------- SCI rcParams ----------
mpl.rcParams.update({
    "font.family":     "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size":       9,
    "pdf.fonttype":    42,
    "ps.fonttype":     42,
})
WONG = {"pink": "#CC79A7", "blue": "#0072B2", "sky": "#56B4E9",
        "orange": "#E69F00", "green": "#009E73"}
METHOD_COLOR = {"LME": WONG["pink"], "ComBat": WONG["blue"],
                "ComBat-joint": WONG["sky"], "RELIEF": WONG["orange"],
                "CovBat": WONG["green"]}
METHOD_ORDER = ["LME", "ComBat", "ComBat-joint", "RELIEF", "CovBat"]

BIOMARKERS = ["T2w_CSA", "GM_CSA_mm2", "MTR", "MTsat",
              "FA", "MD", "AD", "RD"]

# ----------------------------------------------------------------------------
# Paper-canonical method config: maps method name → step3 script + output CSV
# + harmonized column suffix.  Confirmed against the user's local layout.
# ----------------------------------------------------------------------------
METHOD_SPECS = {
    "ComBat":       dict(script="step3_combat.py",
                         out_csv_HC ="biomarkers_combat_HC.csv",
                         out_csv_ALL="biomarkers_combat_ALL.csv",
                         suffix="combat"),
    "ComBat-joint": dict(script="step3_combat_joint.py",
                         out_csv_HC ="biomarkers_combat_joint_HC.csv",
                         out_csv_ALL="biomarkers_combat_joint_ALL.csv",
                         suffix="combatJ"),
    "CovBat":       dict(script="step3_covbat.py",
                         out_csv_HC ="biomarkers_covbat_HC.csv",
                         out_csv_ALL="biomarkers_covbat_ALL.csv",
                         suffix="covbat"),
    "RELIEF":       dict(script="step3_relief.py",
                         out_csv_HC ="biomarkers_relief_HC.csv",
                         out_csv_ALL="biomarkers_relief_ALL.csv",
                         suffix="relief"),
    "LME":          dict(script="step3_lme.py",
                         out_csv_HC ="biomarkers_lme_HC.csv",
                         out_csv_ALL="biomarkers_lme_ALL.csv",
                         suffix="lme"),
}


# ---------------------------------------------------------------------------
# Metric helpers (paper-equivalent)
# ---------------------------------------------------------------------------
def eta2(y: np.ndarray, group: np.ndarray) -> float:
    grand = y.mean()
    ss_total = ((y - grand) ** 2).sum()
    ss_between = sum(len(y[group == g]) * (y[group == g].mean() - grand) ** 2
                     for g in np.unique(group))
    return float(ss_between / ss_total) if ss_total > 0 else 0.0


def mean_eta2(X: np.ndarray, vendor: np.ndarray) -> float:
    return float(np.mean([eta2(X[:, j], vendor) for j in range(X.shape[1])]))


def permanova_r2(X: np.ndarray, vendor: np.ndarray,
                 n_perm: int = 999) -> Tuple[float, float]:
    Z = StandardScaler().fit_transform(X)
    D = np.sqrt(((Z[:, None, :] - Z[None, :, :]) ** 2).sum(-1))
    ids = [f"s{i}" for i in range(len(Z))]
    dm = DistanceMatrix(D, ids)
    res = permanova(dm, vendor, permutations=n_perm)
    F = float(res["test statistic"])
    n = len(vendor); g = len(np.unique(vendor))
    r2 = (F * (g - 1)) / (F * (g - 1) + (n - g))
    return r2, float(res["p-value"])


def within_vendor_r(Xpre: np.ndarray, Xpost: np.ndarray,
                    vendor: np.ndarray) -> float:
    rs = []
    for j in range(Xpre.shape[1]):
        for v in np.unique(vendor):
            m = vendor == v
            if m.sum() < 3: continue
            a, b = Xpre[m, j], Xpost[m, j]
            if np.std(a) == 0 or np.std(b) == 0: continue
            if not (np.isfinite(a).all() and np.isfinite(b).all()): continue
            r, _ = stats.pearsonr(a, b)
            if np.isfinite(r): rs.append(r)
    return float(np.mean(rs)) if rs else np.nan


def age_r2_mean(X: np.ndarray, meta: pd.DataFrame) -> float:
    age = meta["Age"].astype(float).values
    vals = []
    for j in range(X.shape[1]):
        y = X[:, j]
        m = np.isfinite(y) & np.isfinite(age)
        if m.sum() < 5: continue
        r, _ = stats.pearsonr(y[m], age[m])
        vals.append(r * r)
    return float(np.mean(vals)) if vals else np.nan


# ---------------------------------------------------------------------------
# Subprocess: call a step3 script and read its harmonized CSV
# ---------------------------------------------------------------------------
def run_paper_script(script_path: Path, master_csv: Path,
                     out_dir: Path) -> bool:
    """Run a step3_*.py with --master and --out.  Return True on success."""
    out_dir.mkdir(parents=True, exist_ok=True)
    try:
        result = subprocess.run(
            [sys.executable, str(script_path),
             "--master", str(master_csv),
             "--out",    str(out_dir)],
            capture_output=True, text=True, timeout=600,
        )
        if result.returncode != 0:
            print(f"  [WARN] {script_path.name} exit={result.returncode}: "
                  f"{result.stderr.strip()[:200]}")
            return False
        return True
    except subprocess.TimeoutExpired:
        print(f"  [WARN] {script_path.name} timed out (10 min)")
        return False
    except Exception as e:
        print(f"  [WARN] {script_path.name} raised {type(e).__name__}: {e}")
        return False


def read_harmonized(csv_path: Path, suffix: str) -> np.ndarray | None:
    """Read harmonized CSV; extract the 8 BIOMARKERS by <bm>_<suffix> col."""
    if not csv_path.exists():
        return None
    df = pd.read_csv(csv_path)
    cols = [f"{bm}_{suffix}" for bm in BIOMARKERS]
    missing = [c for c in cols if c not in df.columns]
    if missing:
        print(f"  [WARN] missing cols in {csv_path.name}: {missing}")
        return None
    return df[cols].astype(float).values


# ---------------------------------------------------------------------------
def one_seed(df_subset: pd.DataFrame, seed: int,
             code_dir: Path, cohort: str,
             tmp_root: Path) -> List[Dict]:
    """Run the full 5-method evaluation on one subsampled subset."""
    # write tmp master CSV
    tmp_dir = tmp_root / f"seed{seed:03d}"
    tmp_dir.mkdir(parents=True, exist_ok=True)
    tmp_master = tmp_dir / "biomarkers_master.csv"
    df_subset.to_csv(tmp_master, index=False)

    X = df_subset[BIOMARKERS].astype(float).values
    vendor = df_subset["Manufacturer"].values

    rows = []
    base_eta = mean_eta2(X, vendor)
    base_R2, base_p = permanova_r2(X, vendor, n_perm=999)
    base_age = age_r2_mean(X, df_subset)
    rows.append(dict(seed=seed, method="Original",
                     eta2=base_eta, R2_perm=base_R2, p_perm=base_p,
                     r_pre_post=1.0, age_R2=base_age,
                     age_R2_gain=0.0, n=len(df_subset)))

    for method, spec in METHOD_SPECS.items():
        script_path = code_dir / spec["script"]
        method_out  = tmp_dir / method.replace("-", "_")
        ok = run_paper_script(script_path, tmp_master, method_out)
        if not ok:
            rows.append(dict(seed=seed, method=method,
                             eta2=np.nan, R2_perm=np.nan, p_perm=np.nan,
                             r_pre_post=np.nan, age_R2=np.nan,
                             age_R2_gain=np.nan, n=len(df_subset),
                             error="script_failed"))
            continue
        out_key = "out_csv_HC" if cohort == "HC" else "out_csv_ALL"
        Y = read_harmonized(method_out / spec[out_key], spec["suffix"])
        if Y is None or Y.shape != X.shape:
            rows.append(dict(seed=seed, method=method,
                             eta2=np.nan, R2_perm=np.nan, p_perm=np.nan,
                             r_pre_post=np.nan, age_R2=np.nan,
                             age_R2_gain=np.nan, n=len(df_subset),
                             error="csv_read_failed"))
            continue
        # mask rows where the harmonized output has NaN in any biomarker
        # (paper's step3 sometimes drops subjects per-metric; we use the
        # intersection so all metrics share the same n)
        ok_rows = np.isfinite(Y).all(axis=1)
        Xc = X[ok_rows]; Yc = Y[ok_rows]; vc = vendor[ok_rows]
        eta = mean_eta2(Yc, vc)
        R2, pv = permanova_r2(Yc, vc, n_perm=999)
        r = within_vendor_r(Xc, Yc, vc)
        ar2 = age_r2_mean(Yc, df_subset.iloc[ok_rows])
        rows.append(dict(seed=seed, method=method,
                         eta2=eta, R2_perm=R2, p_perm=pv,
                         r_pre_post=r, age_R2=ar2,
                         age_R2_gain=ar2 - base_age,
                         n=int(ok_rows.sum())))

    # cleanup tmp dir to save disk
    shutil.rmtree(tmp_dir, ignore_errors=True)
    return rows


# ---------------------------------------------------------------------------
# Plotting (3-panel)
# ---------------------------------------------------------------------------
def make_figure(res: pd.DataFrame, full_ref: Dict[str, Dict[str, float]],
                out_dir: Path, cohort: str) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(11.5, 3.8))
    panels = [
        ("R2_perm",     "(A)  Residual vendor structure\n"
                        "$R^{2}_{\\mathrm{PERMANOVA}}$ (log)",
         "lower is better", True),
        ("r_pre_post",  "(B)  Subject preservation\n"
                        "$\\bar{r}_{\\mathrm{pre,post}}$",
         "higher is better", False),
        ("age_R2_gain", "(C)  Biology unmasking\n"
                        "$\\Delta$ age $R^{2}$ vs baseline",
         "higher is better", False),
    ]
    data = res[res["method"].isin(METHOD_ORDER)].copy()
    xs = np.arange(len(METHOD_ORDER))
    cols = [METHOD_COLOR[m] for m in METHOD_ORDER]

    for k, (col, ttl, sub, logy) in enumerate(panels):
        ax = axes[k]
        means = data.groupby("method")[col].mean().reindex(METHOD_ORDER)
        sds   = data.groupby("method")[col].std().reindex(METHOD_ORDER)
        ax.bar(xs, means.values, yerr=sds.values, color=cols,
               edgecolor="black", linewidth=0.7,
               error_kw=dict(ecolor="#333", capsize=2.5, lw=0.8),
               alpha=0.85, zorder=3)
        for i, m in enumerate(METHOD_ORDER):
            yvals = data[data["method"] == m][col].dropna().values
            jitter = (np.random.RandomState(42 + i).rand(len(yvals)) - 0.5) * 0.30
            ax.scatter(np.full_like(yvals, xs[i]) + jitter, yvals,
                       s=10, color="white", edgecolor="#222",
                       linewidth=0.5, alpha=0.85, zorder=5)
            ref = full_ref.get(m, {}).get(col, np.nan)
            if np.isfinite(ref):
                ax.hlines(ref, xs[i] - 0.42, xs[i] + 0.42,
                          colors="black", linestyles="--", linewidth=1.1,
                          zorder=6)
        ax.set_xticks(xs)
        ax.set_xticklabels(METHOD_ORDER, rotation=30, ha="right", fontsize=8.5)
        ax.set_title(ttl, fontsize=9.5, fontweight="bold")
        ax.set_xlabel(sub, fontsize=8, color="#555")
        ax.grid(axis="y", alpha=0.3)
        if logy:
            ax.set_yscale("log")

    handles = [
        plt.Line2D([0], [0], color="black", ls="--", lw=1.1,
                   label="Full-cohort reference"),
        plt.Line2D([0], [0], marker="o", color="w",
                   markerfacecolor="white", markeredgecolor="#222",
                   markeredgewidth=0.5, markersize=4,
                   linestyle="", label="Single random seed"),
    ]
    fig.legend(handles=handles, loc="upper center",
               bbox_to_anchor=(0.5, 1.02), ncol=2, frameon=False, fontsize=8.5)
    fig.suptitle(
        f"Design C — Balanced-subsample robustness "
        f"({cohort} cohort, K = {res['seed'].nunique()} seeds; "
        f"$n \\approx$ {int(res['n'].median())} per seed)",
        fontsize=10, fontweight="bold", y=1.07)
    fig.tight_layout()
    out_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_dir / "FigS_design_C_robustness.png", dpi=300,
                bbox_inches="tight")
    fig.savefig(out_dir / "FigS_design_C_robustness.pdf",
                bbox_inches="tight")
    plt.close(fig)
    print(f"[saved] {out_dir / 'FigS_design_C_robustness.png'}")


# ---------------------------------------------------------------------------
def main(master_csv: Path, code_dir: Path, out_dir: Path,
         cohort: str, k_seeds: int, target_siemens: int) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)

    # sanity-check the 5 scripts exist
    missing = [m for m, s in METHOD_SPECS.items()
               if not (code_dir / s["script"]).exists()]
    if missing:
        sys.exit(f"[ERROR] missing scripts in {code_dir}: {missing}")

    df_all = pd.read_csv(master_csv)
    # column-name compatibility: paper's scripts use Manufacturer / Age / Sex
    needed = BIOMARKERS + ["Manufacturer", "Age", "Sex", "Pathology"]
    df = df_all.dropna(subset=BIOMARKERS + ["Manufacturer", "Age", "Sex"]).copy()
    if cohort == "HC":
        df = df[df["Pathology"].astype(str).str.upper() == "HC"].copy()

    counts = df["Manufacturer"].value_counts().to_dict()
    print(f"[cohort={cohort}] vendor counts: {counts}")
    print(f"[design C] subsampling Siemens → {target_siemens} × K = {k_seeds} seeds")
    print(f"[design C] calling 5 paper step3_*.py scripts per seed from {code_dir}")

    tmp_root = Path(tempfile.mkdtemp(prefix="design_C_"))
    print(f"[design C] tmp dir: {tmp_root}")

    all_rows = []
    for s in range(k_seeds):
        ge = df[df["Manufacturer"] == "GE"]
        ph = df[df["Manufacturer"] == "Philips"]
        si = df[df["Manufacturer"] == "Siemens"].sample(
                 n=min(target_siemens, (df["Manufacturer"] == "Siemens").sum()),
                 random_state=s)
        sub = pd.concat([ge, ph, si], ignore_index=True)
        print(f"\n[seed {s:02d}] subset n={len(sub)} "
              f"(GE={len(ge)}, Philips={len(ph)}, Siemens={len(si)})")
        rows = one_seed(sub, seed=s, code_dir=code_dir,
                        cohort=cohort, tmp_root=tmp_root)
        all_rows.extend(rows)
        r_avg = np.nanmean([r["r_pre_post"] for r in rows
                            if r["method"] != "Original"])
        print(f"  -> mean r_pre_post across methods = {r_avg:.3f}")

    res = pd.DataFrame(all_rows)
    res.to_csv(out_dir / "results_design_C.csv", index=False)
    print(f"\n[saved] {out_dir / 'results_design_C.csv'}")

    # full-cohort reference (no subsampling)
    print("\n[design C] computing full-cohort reference (no subsampling) …")
    n_si = (df["Manufacturer"] == "Siemens").sum()
    full_rows = one_seed(df, seed=99999, code_dir=code_dir,
                         cohort=cohort, tmp_root=tmp_root)
    full_ref = {r["method"]: {k: r[k] for k in
                              ("eta2", "R2_perm", "r_pre_post",
                               "age_R2", "age_R2_gain")}
                for r in full_rows}

    make_figure(res, full_ref, out_dir, cohort)

    print("\n=== Design C summary (mean ± SD across seeds) ===")
    summary = res[res["method"].isin(METHOD_ORDER)].groupby("method").agg(
        eta2_mean=("eta2", "mean"), eta2_sd=("eta2", "std"),
        R2_mean=("R2_perm", "mean"), R2_sd=("R2_perm", "std"),
        r_mean=("r_pre_post", "mean"), r_sd=("r_pre_post", "std"),
        age_gain_mean=("age_R2_gain", "mean"),
        age_gain_sd=("age_R2_gain", "std"),
    ).reindex(METHOD_ORDER)
    print(summary.to_string())
    summary.to_csv(out_dir / "summary_design_C.csv")
    print(f"[saved] {out_dir / 'summary_design_C.csv'}")

    # cleanup
    shutil.rmtree(tmp_root, ignore_errors=True)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--master",  required=True, type=Path,
                    help="path to biomarkers_master.csv")
    ap.add_argument("--code",    required=True, type=Path,
                    help="directory containing step3_*.py scripts")
    ap.add_argument("--out",     required=True, type=Path,
                    help="output directory for results & figure")
    ap.add_argument("--cohort",  default="HC", choices=["HC", "ALL"])
    ap.add_argument("--k",       type=int, default=20,
                    help="number of random seeds")
    ap.add_argument("--siemens", type=int, default=40,
                    help="target Siemens subsample size")
    args = ap.parse_args()
    main(args.master, args.code, args.out,
         args.cohort, args.k, args.siemens)
