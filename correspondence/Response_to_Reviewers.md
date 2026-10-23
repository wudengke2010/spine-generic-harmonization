# Response to Reviewers

**Manuscript:** 7830229 — "Scanner Vendor Harmonisation of Spinal Cord Quantitative MRI Biomarkers: A Five-Method Benchmark"
**Journal:** NMR in Biomedicine
**Authors:** Yilin Zhu and Dengke Wu
**Date:** 23 October 2026

Dear Editor and Reviewer,

We thank the editorial team and Reviewer 1 for the careful and constructive assessment. The comments identified real weaknesses — an unclear endpoint framework, an over-interpreted exploratory pathology analysis, and insufficient methodological detail for reproducibility — and addressing them has substantially improved the manuscript. All changes are marked in the revised files; section, table and page numbers below refer to the revised manuscript. Two new analyses were performed for this revision: an age-confounding sensitivity analysis of the pathology-label classification (new Supplementary Material S9) and site-/scanner-model-level harmonisation (new Supplementary Material S10 and main-text Section 3.8). No previously reported number was altered by these additions.

---

## Reviewer 1

### General comments

> **G1.** *The paper is very hard to follow. It is difficult to know what each metric means, what direction is "better", and how methods rank against each other. I suggest a table listing each endpoint, what it measures, which direction is desirable, and its main caveat. I also suggest fewer, more consistently named metrics, and an abstract and key points with fewer numbers and plainer language. The conclusion drawn in the abstract and conclusion are not easily findable in the main text. The authors should consider rewriting to clarify intermediate results and conclusion.*

We agree and have restructured the paper around this comment:

1. **New Table 2 ("Evaluation endpoints at a glance")** lists every endpoint, what it measures, the desirable direction, and its main caveat (Methods, Section 2.4).
2. **Metric naming has been unified.** We renamed "biological-signal recovery" to the more precise **biological-covariate recovery** throughout, and corrected the definition of subject-level preservation: the reported $r_{pre,post}$ is the per-biomarker correlation across subjects between pre- and post-harmonisation values (i.e., preservation of subject ranking), averaged over the eight biomarkers. The previous wording ("per-subject correlation of 8-dimensional biomarker vectors") did not match the implemented quantity; the definition is now stated identically in Methods 2.4 and Supplementary S4.4, and is consistent with the $\sqrt{1-\eta^2}$ ceiling argument in S6. No numerical value changed.
3. **Abstract rewritten**: 295 → 250 words, most numbers removed, plain statements of what was found.
4. **Key Points rewritten** in plainer language with minimal numbers.
5. **Conclusions made traceable**: each Results subsection now ends with a one-sentence mini-conclusion ("In summary, ..."), and the Abstract/Conclusion statements map one-to-one onto Results Sections 3.1–3.8.

> **G2.** *The pathology classification is a general flaw. It is the main evidence for "biological-signal preservation", but it relies on one dataset, internal cross-validation, and the "pathology" is only a radiological assessment from MRI image, "mild compression" is also normal with aging, thus should not be considered as a pathology, no diagnostic were done. It should be validated on an independent external dataset with well-characterised disease. The authors should also give references showing that these biomarkers can discriminate the condition in question.*

We accept this criticism and have made four changes:

1. **The endpoint is demoted to exploratory and no longer carries the "biological-signal preservation" claim.** The endpoint is now presented as "exploratory classification of radiological compression labels" (Methods 2.4, Table 2, Results 3.3, Discussion sixth finding, new Limitation (v)). The Abstract now describes it as "an exploratory endpoint that is substantially confounded by age", and the biological-signal claims rest on the age/sex covariate-recovery analyses instead.
2. **Label provenance is stated explicitly** (Section 2.1): the labels shipped in `participants.tsv` are radiological assessments by two experienced radiologists screening for spinal cord compression (mostly mild), without clinical diagnosis, exactly as described in Valošek et al. 2024 (now cited at this point, ref [3], DOI 10.1162/imag_a_00075). With only two DCM cases the group is essentially mild compression; we say so.
3. **New age-confounding sensitivity analysis (Supplementary S9).** The labelled group is older (median 34.0 vs 28.0 years, $p=2\times10^{-7}$), and age alone discriminates the labels (AUC 0.725) — above every biomarker-based classifier. Critically, after residualising all biomarkers on age and sex, the post-harmonisation improvement *persists* (Original AUC 0.565 → ComBat 0.655, same leakage-free protocol), so the improvement is not an artefact of the age difference. We report this transparently and draw no disease-related conclusion.
4. **References added** showing these biomarkers can discriminate radiological compression / early degeneration: Martin et al. 2018 (BMJ Open 8:e019809; multiparametric qMRI — CSA, FA, MTR, T2\* WM/GM — detects subclinical tissue injury in asymptomatic cervical cord compression) and Grabher et al. 2016 (Sci Rep 6:24636; tract-specific degeneration rostral to the stenotic level in cervical spondylotic myelopathy) — the latter also speaks to why sampling at C2–C3 rather than at the compressed level can still carry signal.

External validation on a clinically characterised disease cohort: we agree this is the correct standard, but no open multi-vendor spinal cord qMRI dataset with clinically diagnosed disease exists to our knowledge, and acquiring one is beyond the scope of this methodological benchmark. Instead of overstating, we have removed all disease-related interpretation, added a dedicated Limitation (v), and state explicitly that external validation on clinically characterised cohorts is required before any clinical interpretation.

> **G3.** *As there are only 3 vendors, the authors should consider harmonizing within a vendor with scanner models or sites.*

Thank you for this concrete suggestion — we performed the analysis (new Section 3.8 and Supplementary S10). Per-metric ComBat with empirical Bayes was refitted on the HC complete-case sample (N=188) using **site** (40 levels) or **scanner model** (12 levels) as the batch variable. Findings: site-level harmonisation removes the vendor effect to the chance level (η²=0.012 vs chance expectation ≈0.011; site η² 0.54→0.045) and recovers slightly more age-related variance, but at the cost of clearly lower subject-level preservation (r=0.703 vs 0.789 for vendor-batch ComBat); scanner-model harmonisation is intermediate (vendor η²=0.0022, r=0.763). We conclude that finer-grained batch definitions are viable when the vendor count is small, but they neither resolve the vendor–site nesting nor offer a transfer route to unseen vendors/sites, and they cost more subject-level distortion per unit of vendor removal. This is now discussed in Section 3.8 and in the Discussion's practical-application paragraph.

### Methods

> **2.1 — Pathology labels were assigned from the images only (radiological assessment), not clinically. This should be stated explicitly.**

Done. Section 2.1 now states: the labels are radiological assessments only, made by two radiologists screening for spinal cord compression without clinical diagnosis, citing Valošek et al. 2024 (ref [3]) where this assessment is documented.

> **2.1 — Mild cervical compression is not necessarily a pathology. It can reflect normal ageing, so we should not expect all qMRI metrics to reflect it. With only 2 DCM cases, the "pathology" group is essentially mild compression. Please cite https://doi.org/10.1162/imag_a_00075, where this assessment was performed. There was no diagnosis, rather radiological assessment from T2w scans. In addition, metrics were performed at C2-C3, not at the levels of compression, please provide references showing that these biomarkers can discriminate the condition in question.**

All points addressed in Section 2.1 and the new Limitation (v): mild compression is common with ageing and not necessarily pathological; the group is essentially mild compression (2 DCM cases); Valošek et al. 2024 (DOI 10.1162/imag_a_00075) is cited as the source of the radiological assessment; the C2–C3 extraction window is acknowledged as not targeting the compressed level. References supporting biomarker discrimination of this condition were added (Martin 2018; Grabher 2016), and the age confound is quantified in the new S9. We no longer claim the pathology endpoint demonstrates "biological-signal preservation".

> **2.1 — The dataset has 267 participants, but fewer were analysed (the original paper used 260). The reason is only given in Section 2.8 and should be stated where the data are first presented.**

Done. Section 2.1 now states, where the data are first presented: release r20231212 contains 267 subjects at 43 sites, adding seven subjects and one site to the 260 participants at 42 sites of the original data descriptor; 21 subjects (7.9%) were missing at least one biomarker — almost always because a contrast was not acquired at the site — leaving 246 complete cases, of whom 188 are healthy controls. Section 2.8 retains the detailed missingness mechanism.

> **2.2 — SCT v7.4 does not correspond to a released version. If the master branch was used, please give the exact git commit.**

Correct, and thank you for catching this. The installed build reports `7.4.dev0`, i.e., a development snapshot of the master branch after the 7.3 release. We verified the exact snapshot against the GitHub repository (file-content hash match): **master commit `93923f5e363d9400fbadab07459f125a08b21ad8`, dated 2 June 2026, sixteen commits after the 7.3 release (8 May 2026)**. This is now reported in Section 2.2 and Supplementary S2.

> **2.2 — Please provide a github link in this section too.**

Done: the analysis code repository URL (https://github.com/wudengke2010/spine-generic-harmonization) is now given in Section 2.2 as well as in Supplementary S8.

> **2.2 — SCT has changed since the original spine-generic pipeline. Were any changes made to the pipeline?**

Yes, and we now disclose them explicitly. The revised Section 2.2 lists four deliberate differences, and the new Supplementary Table S2a compares the two pipelines step by step: (i) the seven level-restricted biomarkers use a common C2–C3 window (the reference pipeline samples MTR/MTsat and DTI over C2–C5 and GM-CSA over C3–C4); (ii) segmentation uses the current contrast-agnostic nnUNet models (`sct_deepseg spinalcord`; `sct_deepseg graymatter`, model release r20250420) instead of the legacy `sct_deepseg_sc`/`sct_deepseg_gm`; (iii) vertebral labelling is fully automatic on T2w without the shipped manual disc labels (the consequences are audited in Supplementary S2.2); (iv) MT and diffusion metrics are averaged within the whole-cord mask rather than the PAM50 white-matter atlas, and we additionally extract AD, which the reference pipeline does not. Because all methods and cohorts are processed identically, these choices affect absolute biomarker values but not the benchmark's method comparisons.

> **2.2 — The authors state that the metrics were extracted at C2-C3 vertebral levels, although the spine-generic pipeline compute CSA-GM at C3 and C4, MTR and DTI metrics between C2 and C5.**

Acknowledged and now explicitly documented: our common C2–C3 window differs from the reference pipeline (GM-CSA at C3–C4; MTR/MTsat and DTI at C2–C5), and GM-CSA is averaged over every segmented cervical slice because the T2\* acquisition does not span the vertebral levels resolved on T2w. See Section 2.2 and Table S2a. The rationale for a common window — identical anatomy across the seven level-restricted biomarkers — is stated.

> **2.2 — Were the derivatives and labels shipped with the dataset used, or were new ones generated?**

New ones were generated. All segmentations, vertebral labels and metric maps were produced de novo with the pipeline described above; the manual disc labels and segmentations shipped in the dataset's `derivatives/labels` directory were not used, so that every subject was processed by one uniform, fully automatic pipeline. This is now stated in Section 2.2 and Supplementary S2, and the cost of fully automatic labelling (15 subjects with displaced labels) is quantified in the S2.2 audit with a sensitivity analysis.

> **2.2 — This section needs to be substantially extended.**

The section has been substantially extended (roughly doubled): SCT snapshot and commit, code repository link, de novo derivatives statement, the four pipeline differences, the extraction-window rationale, and the pointer to the new step-by-step Table S2a.

> **2.9 — This belongs in the Acknowledgements and does not need to be in the Methods. It is currently duplicated in the Declarations.**

Agreed. The AI-assisted-tools statement has been removed from Methods (old Section 2.9) and now appears once, in the new Acknowledgements section. The Declarations file has been updated to point to the Acknowledgements instead of repeating the statement.

### Results

> **Was the pathology classification performed on all biomarkers combined? Please state this explicitly.**

Yes — the classifier used all eight biomarkers jointly as input features (the same feature set as the vendor and sex tasks). This is now stated explicitly in Methods 2.4 and Results 3.3 ("all eight biomarkers combined as input features").

### Discussion

> **The authors could discuss how these methods could be applied in practice, also to external dataset, if new data is added, etc.**

A dedicated practical-application paragraph has been added to the Discussion. It covers: the fit–store–apply deployment pattern (fitted vendor parameters applied to incoming subjects with known vendor labels; affine maps for the ComBat family, projection onto the reference PCA/ICA basis for CovBat/RELIEF); ComBat-Predict for new sites without re-accessing reference data, with the caveat that our LOVO results validate this route only for the ComBat family; site-/model-level harmonisation as a finer-grained alternative where several sites per vendor exist (Section 3.8), including its preservation cost and lack of transfer to unseen sites; and the recommendation to re-estimate and report the residual vendor signature (leakage-free vendor AUC) whenever new data accrue, rather than assuming vendor effects have been eliminated.

---

We hope the revised manuscript is now clearer, fully reproducible in its methodological details, and appropriately cautious in its claims, and we thank the reviewer again for the detailed and fair criticism.

Yours faithfully,
Dengke Wu, MD (on behalf of both authors)
Department of Emergency Medicine, The Second Xiangya Hospital of Central South University
wudk2010@csu.edu.cn
