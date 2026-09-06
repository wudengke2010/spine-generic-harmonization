"""
Step 6.C v3 — Design C: TRUE-EQUAL-n sensitivity analysis
            (wrapper around the paper's own step3_*.py scripts)
============================================================
Per reviewer request ("add a true equal-n sensitivity analysis"):
Instead of Siemens-subsetting to 40 (GE 25 + Philips 33 + Siemens 40 = 98,
residual imbalance ratio max/min = 1.60), we subsample EVERY vendor to the
same n.  In the HC complete-case population (GE 25 / Philips 33 / Siemens 130)
this gives 25 + 25 + 25 = 75 subjects per seed → balance ratio 1.00.

For each random seed:
  1. reads biomarkers_master.csv, restricts to the HC cohort (complete-case
     by the same dropna the wrapper always applies)
  2. subsamples GE/Philips/Siemens to `--ge`, `--philips`, `--siemens`
     (default 25 / 25 / 25) with reproducibility = seed
  3. writes the subset to a tmp CSV and calls each of the 5 paper step3_*.py
     scripts (identical canonical implementations as the main analysis)
  4. recomputes vendor eta², PERMANOVA R², pooled r_{pre,post}
     (paper-canonical metric, identical to Table 2), and age R²
Aggregate over K seeds → mean ± SD, plus a figure.

Output (--out dir):  results_design_C_equaln.csv / summary_design_C_equaln.csv
                     FigS_design_C_equaln.png/pdf

Local run:
  python3 step6_design_C_equaln.py \
      --master  results/biomarkers_master.csv \
      --code    code \
      --out     paper/equaln \
      --cohort  HC --k 20 --ge 25 --philips 25 --siemens 25
"""
from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

# Reuse the paper-canonical helpers / methods from the Siemens-subset wrapper.
sys.path.insert(0, str(Path(__file__).resolve().parent))
from step6_design_C_via_paper_scripts import (  # noqa: E402
    one_seed,
    METHOD_SPECS,
    METHOD_ORDER,
    METHOD_COLOR,
)

mpl.rcParams.update({
    "font.family":     "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size":       9,
    "pdf.fonttype":    42,
    "ps.fonttype":     42,
})
BIOMARKERS = ["T2w_CSA", "GM_CSA_mm2", "MTR", "MTsat",
              "FA", "MD", "AD", "RD"]


def make_figure(res: pd.DataFrame, full_ref: dict, out_dir: Path,
                cohort: str, n_targets: dict) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(11.5, 3.8))
    panels = [
        ("R2_perm", "(a) Residual vendor structure\n$R^{2}_{\\mathrm{PERMANOVA}}$ (log)",
         "lower is better", True),
        ("r_pre_post", "(b) Subject preservation (pooled)\n$\\bar{r}_{\\mathrm{pre,post}}$",
         "higher is better", False),
        ("age_R2_gain", "(c) Biology unmasking\n$\\Delta$ age $R^{2}$ vs baseline",
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
    balance = "balanced" if n_targets["ge"] == n_targets["philips"] == n_targets["siemens"] \
        else "reduced-imbalance"
    fig.suptitle(
        f"Design C — True equal-n sensitivity ({balance}; {cohort} cohort, "
        f"K = {res['seed'].nunique()} seeds; "
        f"$n \\approx$ {int(res['n'].median())} per seed)",
        fontsize=10, fontweight="bold", y=1.07)
    fig.tight_layout()
    out_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_dir / "FigS_design_C_equaln.png", dpi=300, bbox_inches="tight")
    fig.savefig(out_dir / "FigS_design_C_equaln.pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"[saved] {out_dir / 'FigS_design_C_equaln.png'}")


def main(master_csv: Path, code_dir: Path, out_dir: Path,
         cohort: str, k_seeds: int,
         ge_target: int, ph_target: int, si_target: int) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    missing = [m for m, s in METHOD_SPECS.items()
               if not (code_dir / s["script"]).exists()]
    if missing:
        sys.exit(f"[ERROR] missing scripts in {code_dir}: {missing}")

    df_all = pd.read_csv(master_csv)
    needed = BIOMARKERS + ["Manufacturer", "Age", "Sex", "Pathology"]
    df = df_all.dropna(subset=BIOMARKERS + ["Manufacturer", "Age", "Sex"]).copy()
    if cohort == "HC":
        df = df[df["Pathology"].astype(str).str.upper() == "HC"].copy()

    counts = df["Manufacturer"].value_counts().to_dict()
    print(f"[cohort={cohort}] complete-case vendor counts: {counts}")
    n_ge = min(ge_target, counts.get("GE", 0))
    n_ph = min(ph_target, counts.get("Philips", 0))
    n_si = min(si_target, counts.get("Siemens", 0))
    print(f"[design C equal-n] subsampling GE={n_ge}, Philips={n_ph}, "
          f"Siemens={n_si} (n≈{n_ge+n_ph+n_si}) × K = {k_seeds} seeds")

    tmp_root = Path(tempfile.mkdtemp(prefix="design_C_equaln_"))
    all_rows = []
    for s in range(k_seeds):
        ge = df[df["Manufacturer"] == "GE"]
        ph = df[df["Manufacturer"] == "Philips"]
        si = df[df["Manufacturer"] == "Siemens"]
        ge_s = ge.sample(n=n_ge, random_state=s) if len(ge) > n_ge else ge
        ph_s = ph.sample(n=n_ph, random_state=s) if len(ph) > n_ph else ph
        si_s = si.sample(n=n_si, random_state=s) if len(si) > n_si else si
        sub = pd.concat([ge_s, ph_s, si_s], ignore_index=True)
        print(f"\n[seed {s:02d}] subset n={len(sub)} "
              f"(GE={len(ge_s)}, Philips={len(ph_s)}, Siemens={len(si_s)})")
        rows = one_seed(sub, seed=s, code_dir=code_dir,
                        cohort=cohort, tmp_root=tmp_root)
        all_rows.extend(rows)
        r_avg = np.nanmean([r["r_pre_post"] for r in rows
                            if r["method"] != "Original"])
        print(f"  -> mean r_pre_post across methods = {r_avg:.3f}")

    res = pd.DataFrame(all_rows)
    res.to_csv(out_dir / "results_design_C_equaln.csv", index=False)
    print(f"\n[saved] {out_dir / 'results_design_C_equaln.csv'}")

    print("\n[design C equal-n] computing full-cohort reference …")
    full_rows = one_seed(df, seed=99999, code_dir=code_dir,
                         cohort=cohort, tmp_root=tmp_root)
    full_ref = {r["method"]: {k: r[k] for k in
                              ("eta2", "R2_perm", "r_pre_post", "age_R2", "age_R2_gain")}
                for r in full_rows}

    targets = {"ge": n_ge, "philips": n_ph, "siemens": n_si}
    make_figure(res, full_ref, out_dir, cohort, targets)

    print("\n=== Design C equal-n summary (mean ± SD across seeds) ===")
    summary = res[res["method"].isin(METHOD_ORDER)].groupby("method").agg(
        eta2_mean=("eta2", "mean"), eta2_sd=("eta2", "std"),
        R2_mean=("R2_perm", "mean"), R2_sd=("R2_perm", "std"),
        r_mean=("r_pre_post", "mean"), r_sd=("r_pre_post", "std"),
        age_gain_mean=("age_R2_gain", "mean"),
        age_gain_sd=("age_R2_gain", "std"),
    ).reindex(METHOD_ORDER)
    print(summary.to_string())
    summary.to_csv(out_dir / "summary_design_C_equaln.csv")
    print(f"[saved] {out_dir / 'summary_design_C_equaln.csv'}")

    shutil.rmtree(tmp_root, ignore_errors=True)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--master",  required=True, type=Path)
    ap.add_argument("--code",    required=True, type=Path)
    ap.add_argument("--out",     required=True, type=Path)
    ap.add_argument("--cohort",  default="HC", choices=["HC", "ALL"])
    ap.add_argument("--k",       type=int, default=20)
    ap.add_argument("--ge",      type=int, default=25)
    ap.add_argument("--philips", type=int, default=25)
    ap.add_argument("--siemens", type=int, default=25)
    args = ap.parse_args()
    main(args.master, args.code, args.out, args.cohort, args.k,
         args.ge, args.philips, args.siemens)
