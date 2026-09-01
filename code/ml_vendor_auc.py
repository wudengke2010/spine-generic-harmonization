# -*- coding: utf-8 -*-
"""
Step 6: ML vendor-detection AUC analysis (Chen et al. 2022, HBM 43:1179 paradigm).

Random-forest classification of scanner manufacturer (Siemens/GE/Philips) from the
8 spinal cord qMRI biomarkers, before and after each harmonisation method.
Vendor-detectability AUC -> 0.5 indicates removal of vendor signature.
Biology-preservation: pathology (HC vs disease) and sex classification AUC.

Design:
  - Complete cases (all 8 biomarkers non-missing), N=246
  - Repeated stratified 5-fold CV (n_repeats=20), RF with 500 trees
  - Macro one-vs-rest AUC (vendor, 3-class), binary AUC (pathology, sex)
  - Permutation reference (shuffled vendor labels) for chance-level calibration

Outputs: E:/boshi/qm_harmonization_paper/results/step6_mlauc/
  ml_vendor_detection.csv   long-format results
  ml_auc_summary.csv        method x task summary (mean +/- SD)
"""
import warnings
warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import RepeatedStratifiedKFold, cross_val_predict
from sklearn.metrics import roc_auc_score

BASE = "E:/boshi/qm_harmonization_paper/results"
BIOM = ["T2w_CSA", "GM_CSA_mm2", "MTR", "MTsat", "FA", "MD", "AD", "RD"]
SUFFIX = {
    "Original": None,
    "ComBat": "combat",
    "ComBat-joint": "combatJ",
    "CovBat": "covbat",
    "RELIEF": "relief",
    "LME": "lme",
}
FILE = {
    "Original": f"{BASE}/biomarkers_master.csv",
    # v1.3.0: ComBat and LME/FE use the complete-case refits (step3_cc), so the
    # ML analysis evaluates the same fitted models as the univariate panel.
    "ComBat": f"{BASE}/step3_cc/biomarkers_combat_ALL.csv",
    "ComBat-joint": f"{BASE}/step3_combat/biomarkers_combat_joint_ALL.csv",
    "CovBat": f"{BASE}/step3_covbat/biomarkers_covbat_ALL.csv",
    "RELIEF": f"{BASE}/step3_relief/biomarkers_relief_ALL.csv",
    "LME": f"{BASE}/step3_cc/biomarkers_lme_ALL.csv",
}
FILE_HC = {
    "Original": f"{BASE}/biomarkers_master.csv",
    "ComBat": f"{BASE}/step3_cc/biomarkers_combat_HC.csv",
    "ComBat-joint": f"{BASE}/step3_combat/biomarkers_combat_joint_HC.csv",
    "CovBat": f"{BASE}/step3_covbat/biomarkers_covbat_HC.csv",
    "RELIEF": f"{BASE}/step3_relief/biomarkers_relief_HC.csv",
    "LME": f"{BASE}/step3_cc/biomarkers_lme_HC.csv",
}
N_REPEATS = 20
N_FOLDS = 5
SEED = 42


def load_method(name, cohort="ALL"):
    fdict = FILE if cohort == "ALL" else FILE_HC
    df = pd.read_csv(fdict[name])
    if SUFFIX[name] is None:
        cols = BIOM
    else:
        cols = [f"{b}_{SUFFIX[name]}" for b in BIOM]
    df = df[["participant_id", "Manufacturer", "Sex", "Pathology"] + cols].copy()
    df.columns = ["participant_id", "Manufacturer", "Sex", "Pathology"] + BIOM
    return df


def rf_auc_matrix(X, y_vendor, y_path, y_sex, seed):
    """Run repeated stratified CV; return dict of AUCs (macro OVR for vendor)."""
    classes = np.sort(np.unique(y_vendor))
    rf = RandomForestClassifier(
        n_estimators=500, min_samples_leaf=5, class_weight="balanced",
        random_state=seed,
    )
    cv = RepeatedStratifiedKFold(n_splits=N_FOLDS, n_repeats=N_REPEATS, random_state=seed)

    # vendor (multiclass): collect out-of-fold predicted probabilities over all repeats
    probs = np.zeros((len(y_vendor), len(classes), N_REPEATS))
    for r in range(N_REPEATS):
        skf = RepeatedStratifiedKFold(n_splits=N_FOLDS, n_repeats=1, random_state=seed + r)
        parts = list(RepeatedStratifiedKFold(
            n_splits=N_FOLDS, n_repeats=1, random_state=seed + r).split(X, y_vendor))
        # single pass over folds (one repeat)
        p = np.zeros((len(y_vendor), len(classes)))
        for tr, te in parts:
            rf.set_params(random_state=seed + r)
            rf.fit(X[tr], y_vendor[tr])
            p[te] = rf.predict_proba(X[te])
        probs[:, :, r] = p
    auc_vendor = [roc_auc_score(y_vendor, probs[:, :, r], multi_class="ovr",
                                average="macro", labels=classes) for r in range(N_REPEATS)]

    # binary tasks
    def binary_auc(y):
        out = []
        for r in range(N_REPEATS):
            parts = list(RepeatedStratifiedKFold(
                n_splits=N_FOLDS, n_repeats=1, random_state=seed + r).split(X, y))
            p = np.zeros(len(y))
            for tr, te in parts:
                rf.set_params(random_state=seed + r)
                rf.fit(X[tr], y[tr])
                p[te] = rf.predict_proba(X[te])[:, 1]
            out.append(roc_auc_score(y, p))
        return out

    y_path_bin = (y_path != "HC").astype(int)
    y_sex_bin = (y_sex == "M").astype(int)
    return {
        "vendor": auc_vendor,
        "pathology": binary_auc(y_path_bin),
        "sex": binary_auc(y_sex_bin),
    }


def main():
    import os
    outdir = f"{BASE}/step6_mlauc"
    os.makedirs(outdir, exist_ok=True)

    # complete-case mask from master (same subjects across all method files)
    master = pd.read_csv(FILE["Original"])
    mask = master[BIOM].notna().all(axis=1)
    keep_ids = set(master.loc[mask, "participant_id"])

    rows_long = []
    summary = {}

    for name in SUFFIX:
        df = load_method(name)
        df = df[df["participant_id"].isin(keep_ids)].reset_index(drop=True)
        assert len(df) == len(keep_ids), f"{name}: subject mismatch"
        X = df[BIOM].to_numpy()
        y_vendor = df["Manufacturer"].to_numpy()
        y_path = df["Pathology"].to_numpy()
        y_sex = df["Sex"].to_numpy()
        res = rf_auc_matrix(X, y_vendor, y_path, y_sex, SEED)
        summary[name] = {k: (float(np.mean(v)), float(np.std(v))) for k, v in res.items()}
        for task, vals in res.items():
            for r, v in enumerate(vals):
                rows_long.append({"method": name, "task": task, "repeat": r, "auc": v})
        print(f"{name:12s} vendor={summary[name]['vendor'][0]:.3f}+-{summary[name]['vendor'][1]:.3f}  "
              f"pathology={summary[name]['pathology'][0]:.3f}+-{summary[name]['pathology'][1]:.3f}  "
              f"sex={summary[name]['sex'][0]:.3f}+-{summary[name]['sex'][1]:.3f}")

    # permutation reference: shuffled vendor labels on Original features
    df = load_method("Original")
    df = df[df["participant_id"].isin(keep_ids)].reset_index(drop=True)
    X = df[BIOM].to_numpy()
    rng = np.random.default_rng(SEED)
    perm_aucs = []
    for r in range(N_REPEATS):
        y_perm = rng.permutation(df["Manufacturer"].to_numpy())
        classes = np.sort(np.unique(y_perm))
        parts = list(RepeatedStratifiedKFold(
            n_splits=N_FOLDS, n_repeats=1, random_state=SEED + r).split(X, y_perm))
        rf = RandomForestClassifier(n_estimators=500, min_samples_leaf=5,
                                    class_weight="balanced", random_state=SEED + r)
        p = np.zeros((len(y_perm), len(classes)))
        for tr, te in parts:
            rf.fit(X[tr], y_perm[tr])
            p[te] = rf.predict_proba(X[te])
        perm_aucs.append(roc_auc_score(y_perm, p, multi_class="ovr", average="macro",
                                       labels=classes))
    print(f"{'Permutation':12s} vendor={np.mean(perm_aucs):.3f}+-{np.std(perm_aucs):.3f}")
    for r, v in enumerate(perm_aucs):
        rows_long.append({"method": "Permutation", "task": "vendor", "repeat": r, "auc": v})

    pd.DataFrame(rows_long).to_csv(f"{outdir}/ml_vendor_detection.csv", index=False)

    sm = []
    for name, tasks in summary.items():
        for task, (m, s) in tasks.items():
            sm.append({"method": name, "task": task, "auc_mean": round(m, 4),
                       "auc_sd": round(s, 4), "n": len(keep_ids)})
    sm.append({"method": "Permutation", "task": "vendor", "auc_mean": round(float(np.mean(perm_aucs)), 4),
               "auc_sd": round(float(np.std(perm_aucs)), 4), "n": len(keep_ids)})
    sm = pd.DataFrame(sm)
    sm.to_csv(f"{outdir}/ml_auc_summary.csv", index=False)
    print("\nSaved:", outdir)
    print(sm.pivot(index="method", columns="task", values="auc_mean").round(3))

    # ---- HC cohort: vendor task only (complete-case N=188) ----
    master_hc = pd.read_csv(FILE_HC["Original"])
    hc_mask = (master_hc["Pathology"].astype(str).str.upper() == "HC") & master_hc[BIOM].notna().all(axis=1)
    hc_ids = set(master_hc.loc[hc_mask, "participant_id"])
    hc_rows, hc_summary = [], {}
    for name in SUFFIX:
        df = load_method(name, cohort="HC")
        df = df[df["participant_id"].isin(hc_ids)].reset_index(drop=True)
        assert len(df) == len(hc_ids), f"{name} HC: subject mismatch"
        X = df[BIOM].to_numpy()
        y_vendor = df["Manufacturer"].to_numpy()
        classes = np.sort(np.unique(y_vendor))
        aucs = []
        for r in range(N_REPEATS):
            parts = list(RepeatedStratifiedKFold(
                n_splits=N_FOLDS, n_repeats=1, random_state=SEED + r).split(X, y_vendor))
            rf = RandomForestClassifier(n_estimators=500, min_samples_leaf=5,
                                        class_weight="balanced", random_state=SEED + r)
            p = np.zeros((len(y_vendor), len(classes)))
            for tr, te in parts:
                rf.fit(X[tr], y_vendor[tr])
                p[te] = rf.predict_proba(X[te])
            aucs.append(roc_auc_score(y_vendor, p, multi_class="ovr",
                                      average="macro", labels=classes))
        hc_summary[name] = (float(np.mean(aucs)), float(np.std(aucs)))
        print(f"HC {name:12s} vendor={hc_summary[name][0]:.3f}+-{hc_summary[name][1]:.3f}")
        for r, v in enumerate(aucs):
            hc_rows.append({"method": name, "task": "vendor_HC", "repeat": r, "auc": v})
    pd.DataFrame(hc_rows).to_csv(f"{outdir}/ml_vendor_detection_HC.csv", index=False)
    hc_sm = pd.DataFrame([
        {"method": k, "task": "vendor_HC", "auc_mean": round(v[0], 4), "auc_sd": round(v[1], 4)}
        for k, v in hc_summary.items()])
    hc_sm.to_csv(f"{outdir}/ml_auc_summary_HC.csv", index=False)
    print("Saved HC summary:", f"{outdir}/ml_auc_summary_HC.csv")


if __name__ == "__main__":
    main()
