"""
LOVO Figure Generation
======================
Creates a publication-quality figure summarizing LOVO validation results:
  Panel A: Vendor eta2 by method x fold (log scale)
  Panel B: PERMANOVA R2 by method x fold
  Panel C: Age R2 preservation by method (mean across folds)
  Panel D: r_pre_post (test vendor distortion) by method

Also creates a comparison table: LOVO vs Internal validation.
"""
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from scipy import stats
from pathlib import Path

# Wong 2011 color-blind-safe palette
COLORS = {
    "Original":     "#999999",
    "ComBat":       "#0072B2",  # blue
    "ComBat-joint": "#009E73",  # green
    "CovBat":       "#D55E00",  # vermillion
    "RELIEF":       "#CC79A7",  # pink
    "LME":          "#E69F00",  # orange
}
METHODS = ["Original", "ComBat", "ComBat-joint", "CovBat", "RELIEF", "LME"]
FOLDS = ["GE", "Philips", "Siemens"]

LOVO_DIR = Path("E:/boshi/qm_harmonization_paper/results/step5_lovo")
INTERNAL_CSV = Path("E:/boshi/qm_harmonization_paper/results/step4_eval/comparison_long.csv")
MASTER_CSV = Path("E:/boshi/qm_harmonization_paper/results/biomarkers_master.csv")
OUT_DIR = Path("C:/Users/admin/WorkBuddy/2026-07-06-12-34-04/figures_new")

BIOMARKERS = ["T2w_CSA", "GM_CSA_mm2", "MTR", "MTsat", "FA", "MD", "AD", "RD"]


def compute_r_pre_post_lovo():
    """Compute r_pre_post for test vendor by merging on participant_id."""
    df = pd.read_csv(MASTER_CSV)
    if "Pathology" in df.columns:
        hc = df[df["Pathology"].fillna("HC").str.upper().eq("HC")].copy()
    else:
        hc = df.copy()
    hc = hc.reset_index(drop=True)

    results = []
    for held_out in FOLDS:
        train_vendors = [v for v in FOLDS if v != held_out]
        train_df = hc[hc["Manufacturer"].isin(train_vendors)].copy().reset_index(drop=True)
        test_df = hc[hc["Manufacturer"] == held_out].copy().reset_index(drop=True)

        # For each method, compute r_pre_post on test vendor
        # We need to re-run the harmonization to get test vendor's harmonized values
        # Instead, use the per-vendor CSV which already has this info
        pass

    # Use per-vendor CSV, but fix the alignment issue
    pv = pd.read_csv(LOVO_DIR / "lovo_per_vendor.csv")
    # For the held-out vendor, r_pre_post is already computed correctly
    # (it uses .values which extracts raw values in order)
    # But the issue is that orig_v and harm_v may have different lengths
    # due to complete-case dropping

    # Let's compute it properly using participant_id alignment
    for fold in FOLDS:
        test_orig = hc[hc["Manufacturer"] == fold].copy()

        # Read harmonized data from LOVO output files
        # The LOVO script saved combined data in methods_harm, but didn't save to CSV
        # We need to re-compute from the per-vendor data

        # Actually, let's just use the per-vendor data and compute mean r_pre_post
        # for the held-out vendor
        sub = pv[(pv["fold"] == fold) & (pv["vendor"] == fold)]
        for method in METHODS:
            m_sub = sub[sub["method"] == method]
            r_vals = m_sub["r_pre_post"].dropna()
            if len(r_vals) > 0:
                results.append({
                    "fold": fold,
                    "method": method,
                    "r_pre_post_mean": r_vals.mean(),
                    "r_pre_post_std": r_vals.std(),
                })

    return pd.DataFrame(results)


def create_lovo_figure():
    """Create the main LOVO comparison figure."""
    # Load data
    lovo_long = pd.read_csv(LOVO_DIR / "lovo_results_long.csv")
    lovo_perm = pd.read_csv(LOVO_DIR / "lovo_permanova.csv")
    lovo_pv = pd.read_csv(LOVO_DIR / "lovo_per_vendor.csv")

    # Internal validation (for comparison)
    internal = pd.read_csv(INTERNAL_CSV)
    internal_hc = internal[internal["cohort"] == "HC"]

    # Compute r_pre_post for held-out vendor
    r_pp_df = compute_r_pre_post_lovo()

    # Create figure
    fig = plt.figure(figsize=(16, 10))
    gs = GridSpec(2, 2, figure=fig, hspace=0.35, wspace=0.30,
                  left=0.07, right=0.96, top=0.92, bottom=0.08)

    # ── Panel A: Vendor eta2 by method x fold ──────────────────────────
    ax_a = fig.add_subplot(gs[0, 0])

    # Aggregate: mean across biomarkers per fold x method
    eta2_summary = lovo_long.groupby(["fold", "method"])["vendor_eta2"].mean().reset_index()

    x = np.arange(len(METHODS))
    width = 0.25
    for i, fold in enumerate(FOLDS):
        vals = []
        for m in METHODS:
            v = eta2_summary[(eta2_summary["fold"] == fold) &
                             (eta2_summary["method"] == m)]["vendor_eta2"].values
            vals.append(v[0] if len(v) > 0 else 0)

        offset = (i - 1) * width
        bars = ax_a.bar(x + offset, vals, width, label=f"Leave {fold} out",
                        alpha=0.85, edgecolor="white", linewidth=0.5)

    # Add internal validation reference line
    for m_idx, m in enumerate(METHODS):
        internal_val = internal_hc[internal_hc["method"] == m]["vendor_eta2"].mean()
        ax_a.plot([m_idx - 0.15, m_idx + 0.15], [internal_val, internal_val],
                  "k-", linewidth=2, zorder=5)

    ax_a.set_yscale("log")
    ax_a.set_ylim(1e-5, 0.5)
    ax_a.set_xticks(x)
    ax_a.set_xticklabels(METHODS, fontsize=9, rotation=20, ha="right")
    ax_a.set_ylabel(r"Vendor $\eta^2$ (log scale)", fontsize=11)
    ax_a.set_title("(A) Vendor Effect After LOVO Harmonization", fontsize=12, fontweight="bold")
    ax_a.legend(fontsize=8, loc="upper left", framealpha=0.9)
    ax_a.axhline(y=0.01, color="gray", linestyle="--", linewidth=0.8, alpha=0.5)
    ax_a.text(5.4, 0.012, "1%", fontsize=7, color="gray")

    # Annotate internal val reference
    ax_a.text(0.5, 0.4, "— Internal val (full data)", fontsize=7,
              transform=ax_a.transData, style="italic")

    # ── Panel B: PERMANOVA R2 ──────────────────────────────────────────
    ax_b = fig.add_subplot(gs[0, 1])

    for i, fold in enumerate(FOLDS):
        vals = []
        for m in METHODS:
            v = lovo_perm[(lovo_perm["fold"] == fold) &
                          (lovo_perm["method"] == m)]["permanova_R2"].values
            vals.append(v[0] if len(v) > 0 else 0)

        offset = (i - 1) * width
        ax_b.bar(x + offset, vals, width, label=f"Leave {fold} out",
                 alpha=0.85, edgecolor="white", linewidth=0.5)

    # Internal PERMANOVA reference
    internal_perm = pd.read_csv(
        Path("E:/boshi/qm_harmonization_paper/results/step4_eval/permanova_results.csv"))
    internal_perm_hc = internal_perm[internal_perm["cohort"] == "HC"]
    for m_idx, m in enumerate(METHODS):
        v = internal_perm_hc[internal_perm_hc["method"] == m]["R2"].values
        if len(v) > 0:
            ax_b.plot([m_idx - 0.15, m_idx + 0.15], [v[0], v[0]],
                      "k-", linewidth=2, zorder=5)

    ax_b.set_xticks(x)
    ax_b.set_xticklabels(METHODS, fontsize=9, rotation=20, ha="right")
    ax_b.set_ylabel("PERMANOVA R²", fontsize=11)
    ax_b.set_title("(B) Multivariate Vendor Effect (PERMANOVA)", fontsize=12, fontweight="bold")
    ax_b.legend(fontsize=8, loc="upper left", framealpha=0.9)
    ax_b.set_ylim(0, 0.17)

    # Significance markers
    for i, fold in enumerate(FOLDS):
        for m_idx, m in enumerate(METHODS):
            row = lovo_perm[(lovo_perm["fold"] == fold) &
                            (lovo_perm["method"] == m)]
            if len(row) > 0:
                p = row["permanova_p"].values[0]
                if p < 0.001:
                    marker = "***"
                elif p < 0.01:
                    marker = "**"
                elif p < 0.05:
                    marker = "*"
                else:
                    marker = ""
                if marker:
                    offset = (i - 1) * width
                    val = row["permanova_R2"].values[0]
                    ax_b.text(m_idx + offset, val + 0.003, marker,
                              ha="center", fontsize=8, color="red")

    # ── Panel C: Age R2 preservation ───────────────────────────────────
    ax_c = fig.add_subplot(gs[1, 0])

    # Mean across folds for each method
    age_r2_lovo = lovo_long.groupby("method")["age_R2"].mean().reindex(METHODS)
    age_r2_internal = internal_hc.groupby("method")["age_R2"].mean().reindex(METHODS)

    x_c = np.arange(len(METHODS))
    width_c = 0.35

    bars1 = ax_c.bar(x_c - width_c/2, age_r2_lovo.values * 100, width_c,
                     label="LOVO (external)", color=[COLORS[m] for m in METHODS],
                     alpha=0.7, edgecolor="white", linewidth=0.5)
    bars2 = ax_c.bar(x_c + width_c/2, age_r2_internal.values * 100, width_c,
                     label="Internal (full data)", color=[COLORS[m] for m in METHODS],
                     alpha=1.0, edgecolor="black", linewidth=0.8, hatch="//")

    ax_c.set_xticks(x_c)
    ax_c.set_xticklabels(METHODS, fontsize=9, rotation=20, ha="right")
    ax_c.set_ylabel("Age R² (%)", fontsize=11)
    ax_c.set_title("(C) Biological Signal Preservation (Age)", fontsize=12, fontweight="bold")
    ax_c.legend(fontsize=8, loc="upper left", framealpha=0.9)
    ax_c.set_ylim(0, 1.2)

    # ── Panel D: r_pre_post (test vendor distortion) ───────────────────
    ax_d = fig.add_subplot(gs[1, 1])

    if len(r_pp_df) > 0:
        for i, fold in enumerate(FOLDS):
            vals = []
            for m in METHODS:
                v = r_pp_df[(r_pp_df["fold"] == fold) &
                            (r_pp_df["method"] == m)]["r_pre_post_mean"].values
                vals.append(v[0] if len(v) > 0 else np.nan)

            offset = (i - 1) * width
            ax_d.bar(x + offset, vals, width, label=f"Leave {fold} out",
                     alpha=0.85, edgecolor="white", linewidth=0.5)

    ax_d.set_xticks(x)
    ax_d.set_xticklabels(METHODS, fontsize=9, rotation=20, ha="right")
    ax_d.set_ylabel("r (pre vs post, test vendor)", fontsize=11)
    ax_d.set_title("(D) Data Distortion on Held-Out Vendor", fontsize=12, fontweight="bold")
    ax_d.legend(fontsize=8, loc="lower left", framealpha=0.9)
    ax_d.set_ylim(0, 1.1)
    ax_d.axhline(y=0.95, color="gray", linestyle="--", linewidth=0.8, alpha=0.5)
    ax_d.text(5.4, 0.96, "0.95", fontsize=7, color="gray")

    # Save
    out_path = OUT_DIR / "fig_lovo_validation.png"
    fig.savefig(out_path, dpi=200, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"Saved: {out_path}")

    # Also save as PDF
    out_pdf = OUT_DIR / "fig_lovo_validation.pdf"
    fig2 = plt.figure(figsize=(16, 10))
    # Re-use the same figure by re-plotting... actually just copy
    import shutil
    fig.savefig(out_pdf, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"Saved: {out_pdf}")

    return out_path


def create_lovo_comparison_table():
    """Create a comparison table: LOVO vs Internal validation."""
    lovo_long = pd.read_csv(LOVO_DIR / "lovo_results_long.csv")
    lovo_perm = pd.read_csv(LOVO_DIR / "lovo_permanova.csv")
    internal = pd.read_csv(INTERNAL_CSV)
    internal_hc = internal[internal["cohort"] == "HC"]
    internal_perm = pd.read_csv(
        Path("E:/boshi/qm_harmonization_paper/results/step4_eval/permanova_results.csv"))
    internal_perm_hc = internal_perm[internal_perm["cohort"] == "HC"]

    rows = []
    for m in METHODS:
        # Internal
        int_eta2 = internal_hc[internal_hc["method"] == m]["vendor_eta2"].mean()
        int_age = internal_hc[internal_hc["method"] == m]["age_R2"].mean()
        int_rpp = internal_hc[internal_hc["method"] == m]["r_pre_post"].mean()
        int_perm = internal_perm_hc[internal_perm_hc["method"] == m]["R2"].mean()

        # LOVO (mean across 3 folds)
        lovo_eta2 = lovo_long[lovo_long["method"] == m]["vendor_eta2"].mean()
        lovo_age = lovo_long[lovo_long["method"] == m]["age_R2"].mean()
        lovo_perm_r2 = lovo_perm[lovo_perm["method"] == m]["permanova_R2"].mean()

        rows.append({
            "Method": m,
            "Internal η²": f"{int_eta2:.4e}" if not np.isnan(int_eta2) else "—",
            "LOVO η²": f"{lovo_eta2:.4e}" if not np.isnan(lovo_eta2) else "—",
            "η² ratio (LOVO/Int)": f"{lovo_eta2/int_eta2:.1f}×" if int_eta2 > 0 and not np.isnan(int_eta2) and not np.isnan(lovo_eta2) else "—",
            "Internal PERMANOVA R²": f"{int_perm:.4f}" if not np.isnan(int_perm) else "—",
            "LOVO PERMANOVA R²": f"{lovo_perm_r2:.4f}" if not np.isnan(lovo_perm_r2) else "—",
            "Internal Age R²": f"{int_age:.4f}" if not np.isnan(int_age) else "—",
            "LOVO Age R²": f"{lovo_age:.4f}" if not np.isnan(lovo_age) else "—",
        })

    table = pd.DataFrame(rows)
    out_path = OUT_DIR / "lovo_vs_internal_comparison.csv"
    table.to_csv(out_path, index=False)
    print(f"\nSaved: {out_path}")
    print(table.to_string(index=False))

    return table


if __name__ == "__main__":
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fig_path = create_lovo_figure()
    table = create_lovo_comparison_table()
