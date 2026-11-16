# -*- coding: utf-8 -*-
"""v1.7.4: wording corrections from final pre-submission check (8 items).

Each replacement is asserted; a miss aborts the run so nothing is silently skipped.
"""
import io, re, sys
from pathlib import Path

MAIN = Path("submission/main_text.tex")
ESM = Path("submission/supplementary/ESM.tex")
ABS = Path("submission/submission_package/Abstract_Summary.txt")

def apply(path, pairs):
    s = path.read_text(encoding="utf-8")
    for old, new, n in pairs:
        c = s.count(old)
        if c != n:
            print(f"FAIL {path.name}: expected {n} got {c} :: {old[:70]!r}")
            sys.exit(1)
        s = s.replace(old, new)
    path.write_text(s, encoding="utf-8")
    left = len(re.findall("pathology", s, re.I))
    print(f"OK  {path.name}: {len(pairs)} replacements, residual 'pathology' = {left}")

main_pairs = [
    # --- item 1: stop calling the cohort / labels "pathology" ---
    ("and pathology-inclusive ($N=246$) complete-case cohorts",
     "and full-cohort ($N=246$) complete-case samples", 1),
    ("Pathology distribution was balanced across vendors",
     "The distribution of these labels was balanced across vendors", 1),
    ("The pathology labels shipped in the dataset's",
     "The compression labels shipped in the dataset's", 1),
    ("the pathology group is essentially radiological mild compression",
     "the labelled group is essentially radiological mild compression", 1),
    ("applied to the radiological pathology labels",
     "applied to the radiological compression labels", 1),
    ("Pathology AUC (exploratory)", "Compression-label AUC (exploratory)", 1),
    ("the ALL-cohort pathology/sex classification indices",
     "the ALL-cohort compression-label and sex classification indices", 1),
    ("Classification of the exploratory radiological pathology labels",
     "Classification of the exploratory radiological compression labels", 1),
    ("In the pathology-inclusive cohort no such frontier exists",
     "In the full cohort no such frontier exists", 1),
    ("in the pathology-inclusive cohort, ComBat-joint weakly dominated",
     "in the full cohort, ComBat-joint weakly dominated", 1),
    ("pathology-inclusive cohort ($N=246$)", "full cohort ($N=246$)", 1),
    ("pathology-inclusive cohort; HC 0.67", "full cohort; HC 0.67", 1),
    ("in the pathology-inclusive cohort the decontamination axis",
     "in the full cohort the decontamination axis", 1),
    ("(v) \\emph{Exploratory pathology endpoint.} The pathology labels",
     "(v) \\emph{Exploratory compression-label endpoint.} The compression labels", 1),

    # --- item 2: do not call a non-significant change an "improvement" ---
    ("improved from AUC $0.520\\pm0.024$ (unharmonised)",
     "was numerically higher after harmonisation, from AUC $0.520\\pm0.024$ (unharmonised)", 1),
    ("classification of radiological cord-compression labels improved",
     "classification of radiological cord-compression labels was numerically higher but not statistically significant", 1),
    ("the improvement persisted (Original $0.563\\pm0.022$",
     "the numerical increase persisted (Original $0.563\\pm0.022$", 1),
    ("Although the post-harmonisation improvement persisted after age--sex residualisation",
     "Although the numerical increase persisted after age--sex residualisation", 1),

    # --- item 3: residualisation hedge ---
    ("so it is not attributable to the between-group age and sex difference, although residual confounding by unmeasured age-related change cannot be excluded.",
     "so the increase is not fully explained by the modelled linear age and sex effects, although residual confounding by non-linear age effects, site or vendor composition, and unmeasured age-related change cannot be excluded.", 1),

    # --- item 4: figure captions should not claim "biological signal" ---
    ("Age-related biological signal (semi-partial $R^2$) before and after harmonisation",
     "Age-associated variance (semi-partial $R^2$) before and after harmonisation", 1),
    ("Sex-related biological signal (semi-partial $R^2$) before and after harmonisation",
     "Sex-associated variance (semi-partial $R^2$) before and after harmonisation", 1),
    ("Age-related biological signal (semi-partial $R^2$) recovered on the held-out vendor",
     "Age-associated variance (semi-partial $R^2$) recovered on the held-out vendor", 1),

    # --- item 5: LOVO is transductive, not external validation ---
    ("leave-one-vendor-out external validation to test generalisability",
     "leave-one-vendor-out (transductive) validation to test generalisability", 1),
    ("LOVO columns show leave-one-vendor-out external validation",
     "LOVO columns show held-out-vendor (transductive) validation", 1),
    ("\\caption{Leave-one-vendor-out (LOVO) external validation.",
     "\\caption{Leave-one-vendor-out (LOVO) transductive validation.", 1),
    ("Fourth, LOVO external validation isolates a single genuine failure.",
     "Fourth, the LOVO transductive evaluation isolates a single genuine failure.", 1),
    ("(iii) \\emph{Scope of external validation.}",
     "(iii) \\emph{Scope of the held-out-vendor validation.}", 1),
    ("and our LOVO results independently support ComBat-family harmonisation",
     "and our LOVO results provide preliminary within-dataset support for ComBat-family harmonisation", 1),

    # --- item 6: no formal "per unit" metric exists ---
    ("they cost more subject-level distortion per unit of vendor removal.",
     "in this dataset the finer-grained definitions gave lower subject-ranking preservation; given the sparse per-site sample sizes and the nested design, these comparisons should be read as exploratory.", 1),

    # --- item 7: conclusion wording ---
    ("ComBat-joint provides the best overall balance in our sample",
     "ComBat-joint offered the most favourable descriptive trade-off in our sample", 1),
    ("ComBat-joint may therefore be appropriate for broad multi-centre deployment",
     "ComBat-joint may therefore be considered for multicentre datasets of similar structure", 1),

    # --- item 8: GM-CSA does not share the C2-C3 window ---
    ("so that all metrics cover identical anatomy",
     "so that these seven biomarkers cover the same nominal vertebral window", 1),
]

esm_pairs = [
    ("64 had pathology (62 mild cervical compression; 2 degenerative cervical myelopathy)",
     "64 carried radiological compression labels (62 mild cervical compression; 2 degenerative cervical myelopathy)", 1),
    ("Pathology (HC/MC/DCM)", "Compression label (HC/MC/DCM)", 1),
    ("(vendor, pathology, and sex AUC by method and cohort)",
     "(vendor, compression-label, and sex AUC by method and cohort)", 1),
    ("Paired comparison of pathology AUC", "Paired comparison of compression-label AUC", 1),
    ("Age-confounding sensitivity analysis of the pathology-label classification",
     "Age-confounding sensitivity analysis of the compression-label classification", 1),
    ("Pathology classification after in-fold resid",
     "Compression-label classification after in-fold resid", 1),
    ("S9. Age Confounding of the Radiological Pathology Labels",
     "S9. Age Confounding of the Radiological Compression Labels", 1),
    ("The pathology labels are radiological assessments",
     "The compression labels are radiological assessments", 1),
    ("the exploratory pathology classification (main text, Section",
     "the exploratory compression-label classification (main text, Section", 1),
    ("\\textbf{Table S9.} Pathology-label classification before and af",
     "\\textbf{Table S9.} Compression-label classification before and af", 1),
]

apply(MAIN, main_pairs)
apply(ESM, esm_pairs)

s = ABS.read_text(encoding="utf-8")
n = len(re.findall("pathology", s, re.I))
print(f"Abstract_Summary.txt: {n} occurrence(s) of pathology")
for m in re.finditer("pathology", s, re.I):
    print("   ...", s[max(0, m.start()-60):m.end()+40].replace("\n", " "))
