# -*- coding: utf-8 -*-
import io, sys
import numpy as np, pandas as pd
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
R = r"E:/boshi/qm_harmonization_paper/results"
BM = ['T2w_CSA','GM_CSA_mm2','MTR','MTsat','FA','MD','AD','RD']

print("########## 1. ||R_orig||_F ##########")
for coh in ['HC','ALL']:
    d = pd.read_csv(f"{R}/step4_eval/corr_matrices_{coh}.csv")
    for m in ['Original','ComBat','ComBat-joint','CovBat','RELIEF','LME/FE']:
        sub = d[d.method==m]
        P = sub.pivot(index='biomarker1', columns='biomarker2', values='r').loc[BM, BM]
        print(f"  [{coh}] {m:14s} ||R||_F = {np.linalg.norm(P.values):.4f}")

print("\n########## 2. Frobenius % relative ##########")
for coh in ['HC','ALL']:
    d = pd.read_csv(f"{R}/step4_eval/corr_matrices_{coh}.csv")
    mats = {m: d[d.method==m].pivot(index='biomarker1',columns='biomarker2',values='r').loc[BM,BM].values for m in ['Original','ComBat','ComBat-joint','CovBat','RELIEF','LME/FE']}
    n0 = np.linalg.norm(mats['Original'])
    print(f"  [{coh}] ||R_orig||_F = {n0:.4f}")
    for m in ['ComBat','ComBat-joint','CovBat','RELIEF','LME/FE']:
        dist = np.linalg.norm(mats[m]-mats['Original'])
        print(f"      {m:14s} dist={dist:.4f}  pct_of_orig={100*dist/n0:.1f}%")
    print("      --- pairwise between methods ---")
    ms = ['ComBat','ComBat-joint','CovBat','RELIEF','LME/FE']
    for i in range(len(ms)):
        for j in range(i+1,len(ms)):
            print(f"      {ms[i]:14s} vs {ms[j]:14s} = {np.linalg.norm(mats[ms[i]]-mats[ms[j]]):.4f}")

print("\n########## 3. bootstrap paired r ##########")
try:
    b = pd.read_csv(f"{R}/step4_eval/bootstrap_paired_r.csv"); print(b.to_string())
except Exception as e: print("ERR", e)
try:
    b2 = pd.read_csv(f"{R}/step4_eval/bootstrap_r_ci_cc.csv"); print(b2.to_string())
except Exception as e: print("ERR", e)
try:
    b3 = pd.read_csv(f"{R}/step4_eval/bootstrap_eta2_baseline_cc.csv"); print(b3.to_string())
except Exception as e: print("ERR", e)

print("\n########## 4. demographics / per-site ##########")
try:
    dm = pd.read_csv(f"{R}/step4_eval/demographics_completecase.csv")
    print(dm.columns.tolist()); print(dm.head(5).to_string())
except Exception as e: print("ERR", e)
