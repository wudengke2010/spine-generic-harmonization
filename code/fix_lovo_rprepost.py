"""
Fix LOVO r_pre_post alignment and save combined harmonized data.

This script reads the original master data, re-applies the LOVO
harmonization (training on 2 vendors, testing on 1), and:
  1. Saves combined harmonized DataFrames to CSV
  2. Computes r_pre_post correctly aligned by participant_id
  3. Updates lovo_results_long.csv and lovo_summary.csv
"""
import pandas as pd
import numpy as np
from scipy import stats
from pathlib import Path

import lovo_validation as lovo

MASTER_CSV = Path("E:/boshi/qm_harmonization_paper/results/biomarkers_master.csv")
OUT_DIR = Path("E:/boshi/qm_harmonization_paper/results/step5_lovo")
BIOMARKERS = lovo.BIOMARKERS
VENDORS = lovo.VENDORS


def main():
    df = pd.read_csv(MASTER_CSV)
    if "Pathology" in df.columns:
        hc = df[df["Pathology"].fillna("HC").str.upper().eq("HC")].copy()
    else:
        hc = df.copy()
    hc = hc.reset_index(drop=True)

    all_rows = []
    perm_rows = list(pd.read_csv(OUT_DIR / "lovo_permanova.csv").to_dict("records"))
    pv_rows = list(pd.read_csv(OUT_DIR / "lovo_per_vendor.csv").to_dict("records"))

    # We'll keep the existing perm_rows and pv_rows, but recompute r_pre_post
    # in all_rows and update pv_rows

    for held_out in VENDORS:
        train_vendors = [v for v in VENDORS if v != held_out]
        train_df = hc[hc["Manufacturer"].isin(train_vendors)].copy().reset_index(drop=True)
        test_df = hc[hc["Manufacturer"] == held_out].copy().reset_index(drop=True)

        print(f"\nFold: Leave {held_out} Out")

        # Run harmonizations
        methods = {}
        methods["Original"] = pd.concat([train_df, test_df], ignore_index=True)

        tr_harm = train_df.copy()
        te_harm = test_df.copy()
        for b in BIOMARKERS:
            tr_vals, te_vals = lovo.combat_train_transform(train_df, test_df, b)
            tr_harm[b + "_combat"] = tr_vals
            te_harm[b + "_combat"] = te_vals
        methods["ComBat"] = pd.concat([tr_harm, te_harm], ignore_index=True)

        tr_out, te_out = lovo.combat_joint_train_transform(train_df, test_df)
        tr_harm = train_df.copy()
        te_harm = test_df.copy()
        for b in BIOMARKERS:
            tr_harm[b + "_combatJ"] = tr_out[b + "_combatJ"]
            te_harm[b + "_combatJ"] = te_out[b + "_combatJ"]
        methods["ComBat-joint"] = pd.concat([tr_harm, te_harm], ignore_index=True)

        tr_out, te_out = lovo.covbat_train_transform(train_df, test_df)
        tr_harm = train_df.copy()
        te_harm = test_df.copy()
        for b in BIOMARKERS:
            tr_harm[b + "_covbat"] = tr_out[b + "_combatJ"]
            te_harm[b + "_covbat"] = te_out[b + "_combatJ"]
        tr_harm = tr_harm.rename(columns={b+"_combatJ": b+"_covbat" for b in BIOMARKERS})
        te_harm = te_harm.rename(columns={b+"_combatJ": b+"_covbat" for b in BIOMARKERS})
        methods["CovBat"] = pd.concat([tr_harm, te_harm], ignore_index=True)

        tr_out, te_out = lovo.relief_train_transform(train_df, test_df)
        tr_harm = train_df.copy()
        te_harm = test_df.copy()
        for b in BIOMARKERS:
            tr_harm[b + "_relief"] = tr_out[b + "_combatJ"]
            te_harm[b + "_relief"] = te_out[b + "_combatJ"]
        tr_harm = tr_harm.rename(columns={b+"_combatJ": b+"_relief" for b in BIOMARKERS})
        te_harm = te_harm.rename(columns={b+"_combatJ": b+"_relief" for b in BIOMARKERS})
        methods["RELIEF"] = pd.concat([tr_harm, te_harm], ignore_index=True)

        tr_out, te_out = lovo.lme_train_transform(train_df, test_df)
        tr_harm = train_df.copy()
        te_harm = test_df.copy()
        for b in BIOMARKERS:
            tr_harm[b + "_lme"] = tr_out[b + "_lme"]
            te_harm[b + "_lme"] = te_out[b + "_lme"]
        methods["LME"] = pd.concat([tr_harm, te_harm], ignore_index=True)

        # Save combined harmonized data
        for method_name, combined in methods.items():
            out_csv = OUT_DIR / f"lovo_combined_{held_out}_{method_name}.csv"
            combined.to_csv(out_csv, index=False)

        suffix_map = {
            "Original": "", "ComBat": "_combat", "ComBat-joint": "_combatJ",
            "CovBat": "_covbat", "RELIEF": "_relief", "LME": "_lme",
        }

        for method, suf in suffix_map.items():
            combined = methods[method]

            for b in BIOMARKERS:
                col = b + suf
                if col not in combined.columns:
                    continue

                vs = lovo.vendor_stats(combined[col], combined["Manufacturer"])
                bio = lovo.bio_partial_R2(combined[col], combined["Age"], combined["Sex"])

                # r_pre_post: compute correctly aligned by participant_id
                # For test vendor
                test_orig = test_df[["participant_id", b]].copy()
                test_orig = test_orig.rename(columns={b: "original"})
                test_harm = combined.loc[combined["Manufacturer"] == held_out,
                                         ["participant_id", col]].copy()
                test_harm = test_harm.rename(columns={col: "harmonized"})
                merged = test_orig.merge(test_harm, on="participant_id", how="inner")
                merged = merged.dropna()
                if len(merged) >= 3 and merged["original"].std() > 0 and merged["harmonized"].std() > 0:
                    r_pp = stats.pearsonr(merged["original"], merged["harmonized"])[0]
                else:
                    r_pp = np.nan

                all_rows.append({
                    "fold": held_out,
                    "method": method,
                    "biomarker": b,
                    "n": vs["n"],
                    "vendor_eta2": vs["eta2"],
                    "vendor_F": vs["F"],
                    "vendor_p": vs["p"],
                    "r_pre_post_test": r_pp,
                    "age_R2": bio["R2_age"],
                    "sex_R2": bio["R2_sex"],
                })

                # Update per-vendor r_pre_post for all vendors
                for v in VENDORS:
                    v_orig = hc[hc["Manufacturer"] == v][["participant_id", b]].copy()
                    v_orig = v_orig.rename(columns={b: "original"})
                    v_harm = combined.loc[combined["Manufacturer"] == v,
                                          ["participant_id", col]].copy()
                    v_harm = v_harm.rename(columns={col: "harmonized"})
                    v_merged = v_orig.merge(v_harm, on="participant_id", how="inner")
                    v_merged = v_merged.dropna()
                    if len(v_merged) >= 3 and v_merged["original"].std() > 0 and v_merged["harmonized"].std() > 0:
                        r_pp_v = stats.pearsonr(v_merged["original"], v_merged["harmonized"])[0]
                    else:
                        r_pp_v = np.nan

                    # Update pv_rows
                    for row in pv_rows:
                        if (row["fold"] == held_out and row["method"] == method
                                and row["vendor"] == v and row["biomarker"] == b):
                            row["r_pre_post"] = r_pp_v

    # Save updated files
    long_df = pd.DataFrame(all_rows)
    long_df.to_csv(OUT_DIR / "lovo_results_long.csv", index=False)
    print(f"\nUpdated: {OUT_DIR / 'lovo_results_long.csv'} ({len(long_df)} rows)")

    pv_df = pd.DataFrame(pv_rows)
    pv_df.to_csv(OUT_DIR / "lovo_per_vendor.csv", index=False)
    print(f"Updated: {OUT_DIR / 'lovo_per_vendor.csv'} ({len(pv_df)} rows)")

    summary = long_df.groupby(["fold", "method"]).agg(
        vendor_eta2_mean=("vendor_eta2", "mean"),
        vendor_eta2_std=("vendor_eta2", "std"),
        r_pre_post_test_mean=("r_pre_post_test", "mean"),
        r_pre_post_test_std=("r_pre_post_test", "std"),
        age_R2_mean=("age_R2", "mean"),
        sex_R2_mean=("sex_R2", "mean"),
        n_mean=("n", "mean"),
    ).reset_index()
    summary.to_csv(OUT_DIR / "lovo_summary.csv", index=False)
    print(f"Updated: {OUT_DIR / 'lovo_summary.csv'} ({len(summary)} rows)")

    print("\nCorrected LOVO r_pre_post (held-out vendor):")
    for fold in VENDORS:
        print(f"\n  Fold: Leave {fold} Out")
        sub = summary[summary["fold"] == fold]
        for m in ["ComBat", "ComBat-joint", "CovBat", "RELIEF", "LME"]:
            row = sub[sub["method"] == m]
            if len(row) > 0:
                print(f"    {m:13s}  r = {row['r_pre_post_test_mean'].values[0]:.4f}  "
                      f"eta2 = {row['vendor_eta2_mean'].values[0]:.4e}")


if __name__ == "__main__":
    main()
