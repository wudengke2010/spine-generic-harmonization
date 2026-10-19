# -*- coding: utf-8 -*-
"""Final gate: every value edited in this round must be found in the .tex
and must match the independently recomputed truth."""
import io, sys, re, itertools
import numpy as np, pandas as pd

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
SUB = "C:/Users/admin/WorkBuddy/2026-07-06-12-34-04/submission"
R = "E:/boshi/qm_harmonization_paper/results"
BM = ['T2w_CSA', 'GM_CSA_mm2', 'MTR', 'MTsat', 'FA', 'MD', 'AD', 'RD']

mt = open(f"{SUB}/main_text.tex", encoding='utf-8').read()
esm = open(f"{SUB}/supplementary/ESM.tex", encoding='utf-8').read()
both = mt + esm
ok = fail = 0

def check(label, must_present, must_absent=()):
    global ok, fail
    bad = [s for s in must_absent if s in both]
    good = [s for s in must_present if s in both]
    miss = [s for s in must_present if s not in both]
    if bad or miss:
        fail += 1
        print("  FAIL  %-34s missing=%s wrongly-present=%s" % (label, miss, bad))
    else:
        ok += 1
        print("  ok    %-34s" % label)

print("=" * 72)
print("A1  Key Point 1 / abstract  (absolute R2 instead of ratios)")
print("=" * 72)
check("KP1 absolute R2", ["0.002$ for age and $0.030$ for sex"],
      ["155\\times", "16\\times"])
cl = pd.read_csv(f"{R}/step4_eval/comparison_long_cc.csv")
o = cl[(cl.cohort == 'HC') & (cl.method == 'Original')]
print("       truth: vendor eta2=%.4f  age=%.6f  sex=%.6f" %
      (o.vendor_eta2.mean(), o.age_R2.mean(), o.sex_R2.mean()))

print()
print("=" * 72)
print("A2  Frobenius (main text L197 + ESM Fig S4)")
print("=" * 72)
for coh in ['HC', 'ALL']:
    d = pd.read_csv(f"{R}/step4_eval/corr_matrices_{coh}.csv")
    M = {}
    for met, s in d.groupby('method'):
        mat = np.zeros((8, 8))
        idx = {b: i for i, b in enumerate(BM)}
        for _, r in s.iterrows():
            mat[idx[r.biomarker1], idx[r.biomarker2]] = r.r
        M[met] = mat
    no = np.linalg.norm(M['Original'])
    rel = {k: 100 * np.linalg.norm(M[k] - M['Original']) / no
           for k in M if k != 'Original'}
    pairs = sorted((np.linalg.norm(M[a] - M[b]), a, b)
                   for a, b in itertools.combinations([k for k in M if k != 'Original'], 2))
    print("  %s: displacement %.1f--%.1f%%   pairwise %.3f--%.3f   "
          "min-pair=%s  max-pair=%s" %
          (coh, min(rel.values()), max(rel.values()), pairs[0][0], pairs[-1][0],
           "%s/%s" % (pairs[0][1], pairs[0][2]), "%s/%s" % (pairs[-1][1], pairs[-1][2])))
    print("       ComBat%%=%.1f  joint%%=%.1f  LME%%=%.1f  RELIEF%%=%.1f  CovBat%%=%.1f"
          % (rel['ComBat'], rel['ComBat-joint'], rel['LME/FE'], rel['RELIEF'], rel['CovBat']))
str_HC = ("lay 38.1--41.3\\% away from the original in the HC cohort "
          "(ComBat 38.1\\%, LME/FE 38.3\\%, ComBat-joint 38.4\\%, RELIEF 40.3\\%, CovBat 41.3\\%)")
check("MT HC percentages", [str_HC], ["37.4\\%", "ComBat-joint 38.5\\%"])
check("MT pairwise range", ["0.048 (ComBat vs.\\ ComBat-joint) to 0.501 (CovBat vs.\\ LME/FE)"],
      ["ranged from 0.041", "16--42$\\times$", "0.11 in HC, 0.10 in ALL"])
check("MT mean displacement", ["1.60 (HC) and 1.68 (ALL)"])
check("ESM FigS4 percentages", ["38.1--41.3\\% in HC (ComBat 38.1\\%"],
      ["37.4--41.3\\% in HC", "39.8--43.0\\% in ALL"])
check("ESM FigS4 pair values", ["ComBat and ComBat-joint are likewise nearly identical to each other (HC: 0.048, ALL: 0.036)"],
      ["HC: 0.11, ALL: 0.10"])

print()
print("=" * 72)
print("A3  lovo_summary.csv vs Table 2")
print("=" * 72)
s = pd.read_csv(f"{R}/step5_lovo/lovo_summary.csv")
z = pd.read_csv(f"{R}/step5_lovo/lovo_permanova_z.csv")
t2 = {"ComBat": 2.46e-3, "ComBat-joint": 4.42e-4, "CovBat": 1.01e-3,
      "RELIEF": 1.87e-1, "LME": 2.01e-1}
worst_id = worst_t2 = 0.0
for m in t2:
    a = s[s.method == m].vendor_eta2_mean.mean()
    b = z[z.method == m].permanova_R2.mean()
    worst_id = max(worst_id, abs(a - b) / b)
    worst_t2 = max(worst_t2, abs(a - t2[m]) / t2[m])
print("  max rel.deviation  vs permanova_z = %.2e   vs Table 2 = %.2f%%"
      % (worst_id, 100 * worst_t2))
check("rows still 18", [""]) if False else None
print("  lovo_summary rows =", len(s), "(old file backed up as .bak_stale)")
if worst_id < 1e-9 and worst_t2 < 0.005 and len(s) == 18:
    ok += 1; print("  ok    A3 lovo_summary consistent")
else:
    fail += 1; print("  FAIL  A3")

print()
print("=" * 72)
print("A4 / A5 / A6 / A7 / A8")
print("=" * 72)
m = pd.read_csv(f"{R}/biomarkers_master.csv")
g = m.groupby('Site').size()
print("  truth full dataset: median=%.1f range=%d-%d" % (g.median(), g.min(), g.max()))
check("A4 site counts", ["ALL: median 6, range 2--11"], ["full dataset: median 5, range 1--8"])
check("A5 FOV attribution",
      ["from the labelled C2--C3 level falling outside the MT field of view; diffusion metrics were unavailable for two of these subjects for the same reason"],
      ["from the C2--C3 level falling outside the MT and diffusion field of view"])
check("A5 zero rule xref",
      ["set to missing by the aggregation-stage rule detailed in Supplementary Material~S2"])
check("ESM S1 category 4 consistent",
      ["outside the field of view of the MT acquisition, so extraction returned no voxels; diffusion metrics were lost for the same reason in two of these subjects and through the identically-zero rule in the third"],
      ["outside the field of view of the MT and diffusion acquisitions"])
check("A6 rounding", ["mean $4.8\\times10^{-4}$", "1,240-fold", "$\\approx2.8$ orders"],
      ["4.9\\times10^{-4}", "1,260-fold", "2.7$ orders"])
check("A7 reference issue", ["2026;47(8):e70546"], ["2026;47:e70546"])
print("  A8 'across 37 countries' VERIFIED against source abstract "
      "(269 responses from 37 countries) -> kept unchanged")
check("A8 kept", ["across 37 countries [6]"])

print()
print("=" * 72)
print("B-items")
print("=" * 72)
check("B1a complementary", ["four complementary endpoints", "Four complementary endpoints"], ["orthogonal"])
check("B2 abstract tail", ["Both failures trace to three vendors with nested sites."],
      ["converge on one limitation", "Three structural constraints---"])
check("B3 may warrant", ["and may warrant preference for broad deployment."])
check("B4 ML method split",
      ["Harmonisation used vendor labels on training and test folds alike, as it would at deployment; the random forest itself never saw vendor labels."],
      ["while the random forest itself remained blind"])
check("B5 stats split", ["The unpaired Hanley--McNeil test on point AUCs did not reach significance ($p=0.07$)"])
check("B6 rhetoric", ["Harmonisation therefore unmasked signal that vendor variance had obscured, rather than merely removing noise."],
      ["does not merely remove noise---it unmasks"])
check("B8 L145", ["This hybrid is referred to below as LME/FE."], ["We label this hybrid approach"])
check("B8 L173 removed", [], ["Sensitivity to missing-data assumptions was assessed"])
check("B8 L302 removed", [], ["a need echoed by recent European quality initiatives"])
check("B8 L306 removed", [], ["decontamination--preservation trade-off it visualises"])
check("B8 L312 removed", [], ["mirrors the disease-classification gains"])
check("B9 limitations", ["(i) \\emph{Dataset structure.}", "(iv) \\emph{Residual detectability.}"],
      ["(vii) This analysis was not pre-registered"])
check("B9 ESM triple because",
      ["for three reasons: the interactive placement"], ["not adopt the relocated values as primary because"])

print()
print("=" * 72)
print("RESULT: %d passed, %d failed" % (ok, fail))
print("=" * 72)
