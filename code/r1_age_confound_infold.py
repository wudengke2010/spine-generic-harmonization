"""R1-11: is the age/sex residualisation inside the CV folds?

The published sensitivity analysis (r1_age_confound.py) residualises every
biomarker on [1, age, sex] ONCE on the full sample and then runs the same
repeated stratified 5-fold CV.  Residualising on age and sex does not use the
outcome label, so this is not label leakage, but the coefficients are estimated
with test-fold subjects included.

This script repeats A1 (Original) and A2 (ComBat) with the residualisation
refitted inside each training fold, and compares the AUCs with the published
full-sample-residualisation values.

Output: results/step6_mlauc/r1_age_confound_infold.csv
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold

sys.path.insert(0, str(Path(__file__).resolve().parent))
import r1_age_confound as base  # noqa: E402

BIOM = base.BIOM
N_REPEATS, N_FOLDS, SEED = base.N_REPEATS, base.N_FOLDS, base.SEED
OUT = Path("E:/boshi/qm_harmonization_paper/results/step6_mlauc/r1_age_confound_infold.csv")


def design(age, sex):
    lv = sorted(np.unique(sex))
    S = np.column_stack([(np.asarray(sex) == lv[j]).astype(float) for j in range(1, len(lv))])
    return np.column_stack([np.ones(len(age)), np.asarray(age, float), S])


def resid_fit(D, X):
    beta, *_ = np.linalg.lstsq(D, X, rcond=None)
    return beta


def rf_auc_infold(X, y, age, sex, meta=None, method=None):
    aucs = []
    y = np.asarray(y)
    for r in range(N_REPEATS):
        skf = StratifiedKFold(n_splits=N_FOLDS, shuffle=True, random_state=SEED + r)
        p = np.zeros(len(y))
        for tr, te in skf.split(X, y):
            # residualisation coefficients from the TRAINING fold only
            beta = resid_fit(design(age[tr], sex[tr]), X[tr])
            Rtr = X[tr] - design(age[tr], sex[tr]) @ beta
            Rte = X[te] - design(age[te], sex[te]) @ beta

            if method == "ComBat":
                mtr = {k: np.asarray(v)[tr] for k, v in meta.items()}
                mte = {k: np.asarray(v)[te] for k, v in meta.items()}
                prm = base.combat_fit(Rtr, mtr["vendor"], mtr["age"], mtr["sex"])
                Htr = base.combat_apply(prm, Rtr, mtr["vendor"], mtr["age"], mtr["sex"])
                Hte = base.combat_apply(prm, Rte, mte["vendor"], mte["age"], mte["sex"])
            else:
                Htr, Hte = Rtr, Rte

            rf = RandomForestClassifier(n_estimators=500, min_samples_leaf=5,
                                        class_weight="balanced", random_state=SEED + r)
            rf.fit(Htr, y[tr])
            p[te] = rf.predict_proba(Hte)[:, 1]
        aucs.append(roc_auc_score(y, p))
    return float(np.mean(aucs)), float(np.std(aucs))


def main():
    master = pd.read_csv(f"{base.BASE}/biomarkers_master.csv")
    cc = master.dropna(subset=BIOM + ["Age", "Sex", "Manufacturer"]).reset_index(drop=True)
    y = (cc["Pathology"] != "HC").astype(int).values
    age = cc["Age"].values.astype(float)
    sex = cc["Sex"].values
    X = cc[BIOM].values.astype(float)
    meta = {"vendor": cc["Manufacturer"].values, "age": age, "sex": sex}
    print(f"N={len(cc)}  positives={y.sum()}")

    rows = []
    m, s = rf_auc_infold(X, y, age, sex)
    rows.append({"analysis": "Original, in-fold residualisation", "auc_mean": round(m, 4), "auc_sd": round(s, 4)})
    print(f"Original in-fold : {m:.4f} +- {s:.4f}")

    m2, s2 = rf_auc_infold(X, y, age, sex, meta=meta, method="ComBat")
    rows.append({"analysis": "ComBat, in-fold residualisation", "auc_mean": round(m2, 4), "auc_sd": round(s2, 4)})
    print(f"ComBat   in-fold : {m2:.4f} +- {s2:.4f}")

    # published values for side-by-side comparison
    pub = pd.read_csv(Path(base.BASE) / "step6_mlauc/r1_age_confound.csv")
    print("\n== published ==")
    print(pub.to_string(index=False))

    out = pd.DataFrame(rows)
    out["delta_vs_original"] = out["auc_mean"] - out["auc_mean"].iloc[0]
    out.to_csv(OUT, index=False)
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()
