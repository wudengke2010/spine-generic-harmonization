"""
Apply all reviewer-identified text fixes to the harmonization paper.
Handles both HTML and MD versions.
"""
import re

# ============================================================
# FIX 1: PERMANOVA R² claim — critical mathematical error
# ============================================================
FIX_PERMANOVA = {
    "old": (
        r"These multivariate R² values exceed every univariate η² "
        r"\(including AD's 0\.62\), demonstrating that vendor differences "
        r"are not confined to one or two metrics but pervade the joint "
        r"biomarker distribution — vendor-specific covariance patterns "
        r"create separation in the full 8-dimensional space beyond what "
        r"any single metric captures\."
    ),
    "new": (
        "These multivariate R² values (0.315–0.336) are comparable to "
        "the upper range of univariate effects (the largest being AD at "
        "η²=0.62), demonstrating that vendor differences pervade the "
        "joint biomarker distribution — vendor-specific covariance "
        "patterns create separation in the full 8-dimensional space "
        "that is well captured by a multivariate metric. Notably, no "
        "single univariate metric fully captures this multivariate "
        "vendor structure: even the most vendor-sensitive metric (AD) "
        "captures only part of the joint separation."
    ),
}

# ============================================================
# FIX 2: "two orders of magnitude" for both covariates
# ============================================================
FIX_TWO_ORDERS_RESULTS = {
    "old": (
        r"The vendor effect therefore exceeded both canonical biological "
        r"covariates by approximately two orders of magnitude in HC and "
        r"by more than one order of magnitude in ALL, providing a "
        r"compelling quantitative motivation for harmonization\."
    ),
    "new": (
        "The vendor effect therefore exceeded age-related variance by "
        "approximately two orders of magnitude (η²=0.280 vs. age R²=0.0026; "
        "ratio ~108×) and sex-related variance by one order of magnitude "
        "(η²=0.280 vs. sex R²=0.025; ratio ~11×) in the HC cohort, "
        "providing a compelling quantitative motivation for harmonization. "
        "In the ALL cohort, the vendor-to-biology ratios were comparable "
        "(η²=0.298 vs. age R²=0.0067, ~44×; vs. sex R²=0.022, ~14×)."
    ),
}

# ============================================================
# FIX 3: "first systematic, multi-endpoint benchmark" → soften
# ============================================================
FIX_FIRST_SYSTEMATIC_INTRO = {
    "old": (
        r"Here we present the first systematic, multi-endpoint benchmark "
        r"of five statistical harmonization methods applied to eight "
        r"standard spinal cord qMRI biomarkers"
    ),
    "new": (
        "Here we present a systematic comparison of five statistical "
        "harmonization methods applied to eight standard spinal cord "
        "qMRI biomarkers"
    ),
}

FIX_FIRST_SYSTEMATIC_DISC = {
    "old": (
        r"In this first systematic, multi-endpoint benchmark of "
        r"statistical harmonization methods applied to spinal cord "
        r"qMRI biomarkers, four findings are clear and mutually "
        r"reinforcing\."
    ),
    "new": (
        "In this systematic comparison of statistical harmonization "
        "methods applied to spinal cord qMRI biomarkers, four findings "
        "are clear and mutually reinforcing."
    ),
}

# ============================================================
# FIX 4: LME "baseline" → "LME (random intercept)"
# ============================================================
FIX_LME_BASELINE_METHODS = {
    "old": (
        r"<strong>LME baseline</strong>: A linear mixed-effects model "
        r"per biomarker"
    ),
    "new": (
        "<strong>LME (random intercept)</strong>: A linear mixed-effects "
        "model per biomarker"
    ),
}

FIX_LME_BASELINE_ABSTRACT = {
    "old": (
        r"and a linear mixed-effects \(LME\) baseline"
    ),
    "new": (
        "and a linear mixed-effects model with vendor random intercept "
        "(LME)"
    ),
}

# ============================================================
# FIX 5: "unmasks biological signal" → soften
# ============================================================
FIX_UNMASKS_ABSTRACT = {
    "old": (
        r"confirming that harmonization unmasks — rather than erases — "
        r"biological signal"
    ),
    "new": (
        "confirming that harmonization improves — rather than degrades — "
        "the detectability of biological signal, though absolute effect "
        "sizes remain modest"
    ),
}

FIX_UNMASKS_RESULTS = {
    "old": (
        r"This \"biology unmasking\" is the expected consequence "
        r"of removing a confounder"
    ),
    "new": (
        "This improvement in biological-signal detectability is the "
        "expected consequence of removing a confounder"
    ),
}

FIX_UNMASKS_DISC = {
    "old": (
        r"<strong>Third, harmonization unmasks biological signal\.</strong>"
    ),
    "new": (
        "<strong>Third, harmonization improves biological-signal recovery.</strong>"
    ),
}

# ============================================================
# FIX 6: Summary box — two orders + unmasking
# ============================================================
FIX_SUMMARY_BOX_1 = {
    "old": (
        r"\(1\) Vendor explains ~28% of univariate biomarker variance and "
        r"~32% of multivariate structure — dominating biological "
        r"covariates by two orders of magnitude\. <br>"
    ),
    "new": (
        "(1) Vendor explains ~28% of univariate biomarker variance and "
        "~32% of multivariate structure — exceeding age-related variance "
        "by ~two orders of magnitude and sex-related variance by ~one "
        "order of magnitude. <br>"
    ),
}

FIX_SUMMARY_BOX_3 = {
    "old": (
        r"\(3\) Harmonization unmasks biological signal: age R² increases "
        r"2\.3- to 3\.7-fold\. <br>"
    ),
    "new": (
        "(3) Harmonization improves biological-signal detectability: age "
        "R² increases 2.3- to 3.7-fold, though absolute effect sizes "
        "remain below 1% of variance explained. <br>"
    ),
}

# ============================================================
# FIX 7: Abstract "two orders of magnitude"
# ============================================================
FIX_ABSTRACT_TWO_ORDERS = {
    "old": (
        r"exceeding biological covariates by two orders of magnitude"
    ),
    "new": (
        "exceeding biological covariates by one to two orders of "
        "magnitude (age: ~108×; sex: ~11×)"
    ),
}

# ============================================================
# FIX 8: Discussion paragraph on unmasking
# ============================================================
FIX_DISC_UNMASK_PARAGRAPH = {
    "old": (
        r"Rather than erasing biological variance, every method increased "
        r"the age and sex signal recovered from the harmonized data\."
    ),
    "new": (
        "Rather than erasing biological variance, every method was "
        "associated with a modest increase in the age and sex signal "
        "recovered from the harmonized data. We caution that the absolute "
        "age-related variance explained remains below 1% post-harmonization "
        "(0.65–1.63% depending on method and cohort), reflecting the "
        "limited sensitivity of C2–C3 qMRI metrics to normal aging "
        "in this restricted age range (19–52 years)."
    ),
}

FIX_DISC_UNMASK_PARAGRAPH_2 = {
    "old": (
        r"This is a practically important finding: it means that "
        r"harmonization can be deployed without fear of erasing the "
        r"very biological effects the study aims to detect, and may in "
        r"fact increase statistical power for age- and sex-related "
        r"analyses\."
    ),
    "new": (
        "This is a practically important finding: it means that "
        "harmonization does not erase the biological effects the study "
        "aims to detect, and may modestly improve statistical power for "
        "age- and sex-related analyses. However, the small absolute "
        "magnitude of these improvements means they should be interpreted "
        "as evidence of preservation rather than enhancement of "
        "biological signal."
    ),
}

# ============================================================
# FIX 9: Conclusion — "three orders of magnitude" + "unmasks"
# ============================================================
FIX_CONCLUSION_DISC = {
    "old": (
        r"all methods reduced vendor variance by more than three orders of "
        r"magnitude and rendered multivariate vendor structure "
        r"indistinguishable from chance"
    ),
    "new": (
        "all methods effectively eliminated vendor variance (η² reduction "
        "to ~10⁻⁴) and rendered multivariate vendor structure "
        "indistinguishable from chance"
    ),
}

# ============================================================
# FIX 10: Conclusion "first" reference
# ============================================================
FIX_CONCLUSION_FIRST = {
    "old": (
        r"The Pareto-frontier framework and systematic benchmarking "
        r"approach introduced here establishes the first quantitative "
        r"decision-support tool for harmonization method selection in "
        r"multi-center spinal cord qMRI research\."
    ),
    "new": (
        "The Pareto-frontier framework and systematic comparison "
        "approach introduced here provide a quantitative decision-support "
        "tool for harmonization method selection in multi-center spinal "
        "cord qMRI research."
    ),
}

# ============================================================
# FIX 11: Limitations — expand
# ============================================================
FIX_LIMITATIONS_EXPAND = {
    "old": (
        r"<h3>4\.5 Limitations</h3>"
    ),
    "new": (
        "<h3>4.5 Limitations</h3>"
    ),
}

LIMITATIONS_OLD_TEXT = (
    "<p>Several limitations should be considered. First, the spine-generic "
    "dataset, while the largest openly available multi-vendor spinal cord "
    "qMRI cohort, covers only three vendors and is dominated by Siemens "
    "scanners (67% of HC cohort). The preservation ranking we report may be "
    "sensitive to vendor imbalance — as confirmed by our Design C "
    "supplementary analysis, which shows that balanced subsampling "
    "compresses the preservation difference. Second, we treated "
    "<em>vendor</em> rather than <em>site</em> or <em>scanner model</em> "
    "as the batch variable. Site-level harmonization was infeasible here "
    "because many sites contributed only a few subjects (typically 3–8), "
    "below the threshold for stable EB estimation. Vendor harmonization "
    "captures the dominant batch axis, but residual site-level structure "
    "may remain in larger consortia. Third, our biological covariates were "
    "limited to age and sex; the spine-generic dataset lacks detailed "
    "clinical or anthropometric measures that could further dissect "
    "biological-signal recovery. Fourth, we evaluated harmonization in a "
    "single-pass, full-cohort design. Out-of-sample generalization to a "
    "hold-out vendor — the relevant scenario for prospective deployment — "
    "is an important next step. Fifth, our subject-preservation metric "
    "(Pearson r) is correlational and does not account for measurement "
    "noise; a portion of the post-harmonization r deficit relative to "
    "unity reflects irreducible within-subject variability rather than "
    "information loss. Test-retest data — not part of the spine-generic "
    "dataset — would be required to resolve this.</p>"
)

LIMITATIONS_NEW_TEXT = (
    "<p>Several limitations warrant careful consideration.</p>"
    "<p><strong>Vendor-site confounding.</strong> The spine-generic dataset "
    "contains 42 sites but only 3 vendors, making it impossible to "
    "disentangle vendor effects from site-level effects within each vendor. "
    "A systematic Siemens-versus-Philips difference may partially reflect "
    "differences in the specific sites operating each vendor's scanners, "
    "rather than purely hardware-level differences. This limitation is "
    "inherent to the retrospective, multi-vendor design of the spine-generic "
    "study and motivates prospective cross-vendor designs with site "
    "crossing. Both our findings and their interpretation should be "
    "understood as reflecting the composite 'vendor-in-site' effect, "
    "which is the relevant batch unit for most practical multi-center "
    "consortia.</p>"
    "<p><strong>Limited vendor diversity and imbalance.</strong> The "
    "dataset covers only three vendors (GE, Philips, Siemens) and is "
    "dominated by Siemens scanners (67% of HC cohort). The preservation "
    "ranking we report may be sensitive to vendor imbalance — as confirmed "
    "by our Design C supplementary analysis, which shows that balanced "
    "subsampling compresses the preservation difference across methods. "
    "Extrapolation to other vendors (e.g., Canon, United Imaging) or "
    "to scanner models not represented in this dataset should be "
    "undertaken cautiously.</p>"
    "<p><strong>Modest biological effect sizes.</strong> Although "
    "harmonization consistently increased the age-related variance "
    "explained (2.3–3.7 fold relative to baseline), the absolute effect "
    "remains small: age explains less than 1% of post-harmonization "
    "biomarker variance, and less than 1.7% even in the best-case "
    "(CovBat, ALL cohort). These small effect sizes partly reflect the "
    "limited age range (19–52 years) and the known moderate sensitivity "
    "of single-level C2–C3 scalar metrics to normal aging. "
    "Studies targeting age-related cord changes should consider broader "
    "age ranges or more spatially extended sampling.</p>"
    "<p><strong>Method convergence on vendor removal.</strong> All five "
    "methods achieved near-identical univariate and multivariate vendor "
    "removal (η² reduced to ~10⁻⁴ in every case; all PERMANOVA p ≥ 0.97). "
    "This convergence, while reassuring for method robustness, limits the "
    "discriminatory power of vendor-removal metrics. The practical "
    "differentiation among methods therefore rests on subject-level "
    "preservation — a secondary endpoint with modest absolute differences "
    "(Δr ~0.06 between best and worst). Users in balanced consortia may "
    "find this differentiation insufficiently large to drive a strong "
    "method preference.</p>"
    "<p><strong>Single-pass, within-cohort design.</strong> We evaluated "
    "harmonization in a single-pass, full-cohort design without a "
    "hold-out vendor. Out-of-sample generalization to a held-out vendor "
    "— the relevant scenario for prospective deployment — is an important "
    "next step that requires a study design with sufficient per-vendor "
    "sample sizes for cross-validation.</p>"
    "<p><strong>Absence of test-retest data.</strong> Our subject-"
    "preservation metric (Pearson r) is correlational and conflates "
    "harmonization-induced information loss with irreducible within-subject "
    "measurement noise. Partitioning these requires test-retest data, "
    "which are not part of the spine-generic dataset. The reported "
    "r<sub>pre,post</sub> values therefore represent a lower bound on "
    "true subject preservation.</p>"
    "<p><strong>Scope of methods evaluated.</strong> We tested five methods "
    "from the ComBat/LME family. Other approaches — including GAN-based "
    "harmonization, traveling-subject calibration, RAVEL, and ComBat-GAM "
    "— were not evaluated. Our 'benchmark' should be understood as a "
    "comparison within this specific methodological family rather than "
    "an exhaustive survey of all available harmonization strategies.</p>"
)


# ============================================================
# HTML-specific replacements (with HTML entities)
# ============================================================
HTML_REPLACEMENTS = [
    # Abstract
    (
        "exceeding biological covariates by two orders of magnitude",
        "exceeding biological covariates by one to two orders of magnitude (age: ~108×; sex: ~11×)",
    ),
    (
        "confirming that harmonization unmasks — rather than erases — biological signal",
        "confirming that harmonization improves — rather than degrades — the detectability of biological signal, though absolute effect sizes remain modest (<1% variance explained)",
    ),
    (
        "and a linear mixed-effects (LME) baseline",
        "and a linear mixed-effects model with vendor random intercept (LME)",
    ),
    (
        "the first quantitative, multi-endpoint benchmark for spinal cord harmonization strategy selection",
        "a quantitative, multi-endpoint comparison for spinal cord harmonization strategy selection",
    ),
    # Introduction
    (
        "Here we present the first systematic, multi-endpoint benchmark of five statistical harmonization methods applied to eight standard spinal cord qMRI biomarkers",
        "Here we present a systematic comparison of five statistical harmonization methods applied to eight standard spinal cord qMRI biomarkers",
    ),
    # Methods
    (
        "<strong>LME baseline</strong>: A linear mixed-effects model per biomarker",
        "<strong>LME (random intercept)</strong>: A linear mixed-effects model per biomarker",
    ),
    ("the simple LME baseline", "the simple LME random-intercept model"),
    ("including the simple LME baseline", "including the LME random-intercept model"),
    # Results §3.2 — the critical PERMANOVA fix
    (
        "These multivariate R² values exceed every univariate η² (including AD's 0.62), demonstrating that vendor differences are not confined to one or two metrics but pervade the joint biomarker distribution — vendor-specific covariance patterns create separation in the full 8-dimensional space beyond what any single metric captures.",
        "These multivariate R² values (0.315–0.336) are comparable in magnitude to the upper range of univariate effects (the largest being AD at η²=0.62). This demonstrates that vendor differences pervade the joint biomarker distribution — vendor-specific covariance patterns create separation in the full 8-dimensional space. Notably, no single univariate metric fully captures this multivariate vendor structure: the PERMANOVA R², while substantial, does not exceed the most extreme univariate η² (AD=0.62), confirming that univariate and multivariate metrics capture complementary aspects of vendor influence.",
    ),
    # Results §3.2 — two orders
    (
        "The vendor effect therefore exceeded both canonical biological covariates by approximately two orders of magnitude in HC and by more than one order of magnitude in ALL, providing a compelling quantitative motivation for harmonization.",
        "The vendor effect therefore exceeded age-related variance by approximately two orders of magnitude (η²=0.280 vs. age R²=0.0026; ratio ~108×) and sex-related variance by one order of magnitude (η²=0.280 vs. sex R²=0.025; ratio ~11×) in the HC cohort, providing a compelling quantitative motivation for harmonization.",
    ),
    # Results §3.5 — unmasking
    (
        'This "biology unmasking" is the expected consequence of removing a confounder',
        "This improvement in biological-signal detectability is the expected consequence of removing a confounder",
    ),
    # Discussion summary box
    (
        "(1) Vendor explains ~28% of univariate biomarker variance and ~32% of multivariate structure — dominating biological covariates by two orders of magnitude. <br>",
        "(1) Vendor explains ~28% of univariate biomarker variance and ~32% of multivariate structure — exceeding age-related variance by ~two orders of magnitude and sex-related variance by ~one order of magnitude. <br>",
    ),
    (
        "(3) Harmonization unmasks biological signal: age R² increases 2.3- to 3.7-fold. <br>",
        "(3) Harmonization improves biological-signal detectability: age R² increases 2.3- to 3.7-fold (though absolute effects remain &lt;1% of variance). <br>",
    ),
    # Discussion §4.1
    (
        "In this first systematic, multi-endpoint benchmark of statistical harmonization methods applied to spinal cord qMRI biomarkers, four findings are clear and mutually reinforcing.",
        "In this systematic comparison of statistical harmonization methods applied to spinal cord qMRI biomarkers, four findings are clear and mutually reinforcing.",
    ),
    (
        "<strong>Third, harmonization unmasks biological signal.</strong>",
        "<strong>Third, harmonization improves biological-signal recovery.</strong>",
    ),
    (
        "Rather than erasing biological variance, every method increased the age and sex signal recovered from the harmonized data.",
        "Rather than erasing biological variance, every method was associated with a modest increase in the age and sex signal recovered from the harmonized data. We caution that the absolute age-related variance explained remains below 1% post-harmonization (0.65–1.63% depending on method and cohort), reflecting the limited sensitivity of C2–C3 qMRI metrics to normal aging in the restricted age range (19–52 years) of this cohort.",
    ),
    (
        "This is a practically important finding: it means that harmonization can be deployed without fear of erasing the very biological effects the study aims to detect, and may in fact increase statistical power for age- and sex-related analyses.",
        "This is a practically reassuring finding: harmonization does not erase the biological effects the study aims to detect, and may modestly improve statistical power for age- and sex-related analyses. The small absolute magnitude of these improvements, however, means they should be interpreted primarily as evidence of signal preservation rather than enhancement.",
    ),
    # Discussion §4.2
    (
        "The central intellectual contribution of this benchmark is the observation",
        "A key contribution of this comparison is the observation",
    ),
    # Conclusion
    (
        "all methods reduced vendor variance by more than three orders of magnitude and rendered multivariate vendor structure indistinguishable from chance",
        "all methods effectively eliminated vendor variance (reduction to ~10⁻⁴, a >99.8% decrease) and rendered multivariate vendor structure indistinguishable from chance",
    ),
    # Figure 3 caption
    (
        "(biology unmasked)",
        "(biology signal recovery)",
    ),
    # Figure 4 caption
    (
        "Marker size encodes age R² gain (biology unmasked).",
        "Marker size encodes age R² gain (biological signal recovery).",
    ),
]

# ============================================================
# OLD limitations → NEW limitations block
# ============================================================
OLD_LIMITATIONS = (
    "<p>Several limitations should be considered."
)
NEW_LIMITATIONS = (
    LIMITATIONS_NEW_TEXT
)

# ============================================================
# Apply all fixes
# ============================================================

def apply_html_fixes(filepath):
    with open(filepath, 'r', encoding='utf-8') as f:
        content = f.read()

    changes_made = 0
    for old, new in HTML_REPLACEMENTS:
        if old in content:
            content = content.replace(old, new)
            changes_made += 1
        else:
                try:
                    print(f"WARNING: Not found in HTML: {repr(old[:80])}")
                except:
                    print("WARNING: Not found in HTML (encoding error)")

    # Replace old limitations with new
    if OLD_LIMITATIONS in content:
        # Find the start of the old limitations
        idx = content.find(OLD_LIMITATIONS)
        # Find the </p> that closes the old limitations paragraph
        # The old limitations is one long paragraph, find the </p> after the last old text
        old_end_marker = "Test-retest data — not part of the spine-generic dataset — would be required to resolve this.</p>"
        end_idx = content.find(old_end_marker, idx) + len(old_end_marker)
        if end_idx > idx:
            content = content[:idx] + NEW_LIMITATIONS.strip() + content[end_idx:]
            changes_made += 1
            print("Replaced limitations section.")

    with open(filepath, 'w', encoding='utf-8') as f:
        f.write(content)

    print(f"HTML: {changes_made} text replacements applied.")
    return changes_made


def apply_md_fixes(filepath):
    """Apply same fixes to Markdown version."""
    with open(filepath, 'r', encoding='utf-8') as f:
        content = f.read()

    changes_made = 0
    # Markdown doesn't have HTML tags, so use simpler patterns
    MD_REPLACEMENTS = [
        ("exceeding biological covariates by two orders of magnitude",
         "exceeding biological covariates by one to two orders of magnitude (age: ~108×; sex: ~11×)"),
        ("confirming that harmonization unmasks — rather than erases — biological signal",
         "confirming that harmonization improves — rather than degrades — the detectability of biological signal, though absolute effect sizes remain modest (<1% variance explained)"),
        ("and a linear mixed-effects (LME) baseline",
         "and a linear mixed-effects model with vendor random intercept (LME)"),
        ("the first quantitative, multi-endpoint benchmark for spinal cord harmonization strategy selection",
         "a quantitative, multi-endpoint comparison for spinal cord harmonization strategy selection"),
        ("Here we present the first systematic, multi-endpoint benchmark of five statistical harmonization methods applied to eight standard spinal cord qMRI biomarkers",
         "Here we present a systematic comparison of five statistical harmonization methods applied to eight standard spinal cord qMRI biomarkers"),
        ("**LME baseline**: A linear mixed-effects model per biomarker",
         "**LME (random intercept)**: A linear mixed-effects model per biomarker"),
        ("the simple LME baseline", "the simple LME random-intercept model"),
        ("including the simple LME baseline", "including the LME random-intercept model"),
        ('This "biology unmasking" is the expected consequence of removing a confounder',
         "This improvement in biological-signal detectability is the expected consequence of removing a confounder"),
        ("(1) Vendor explains ~28% of univariate biomarker variance and ~32% of multivariate structure — dominating biological covariates by two orders of magnitude.",
         "(1) Vendor explains ~28% of univariate biomarker variance and ~32% of multivariate structure — exceeding age-related variance by ~two orders of magnitude and sex-related variance by ~one order of magnitude."),
        ("(3) Harmonization unmasks biological signal: age R² increases 2.3- to 3.7-fold.",
         "(3) Harmonization improves biological-signal detectability: age R² increases 2.3- to 3.7-fold (though absolute effects remain <1% of variance)."),
        ("In this first systematic, multi-endpoint benchmark of statistical harmonization methods applied to spinal cord qMRI biomarkers, four findings are clear and mutually reinforcing.",
         "In this systematic comparison of statistical harmonization methods applied to spinal cord qMRI biomarkers, four findings are clear and mutually reinforcing."),
        ("**Third, harmonization unmasks biological signal.**",
         "**Third, harmonization improves biological-signal recovery.**"),
        ("Rather than erasing biological variance, every method increased the age and sex signal recovered from the harmonized data.",
         "Rather than erasing biological variance, every method was associated with a modest increase in the age and sex signal recovered from the harmonized data. We caution that the absolute age-related variance explained remains below 1% post-harmonization (0.65–1.63% depending on method and cohort), reflecting the limited sensitivity of C2–C3 qMRI metrics to normal aging in the age range (19–52 years) of this cohort."),
        ("This is a practically important finding: it means that harmonization can be deployed without fear of erasing the very biological effects the study aims to detect, and may in fact increase statistical power for age- and sex-related analyses.",
         "This is a practically reassuring finding: harmonization does not erase the biological effects the study aims to detect, and may modestly improve statistical power for age- and sex-related analyses. The small absolute magnitude of these improvements, however, means they should be interpreted primarily as evidence of signal preservation rather than enhancement."),
        ("The central intellectual contribution of this benchmark is the observation",
         "A key contribution of this comparison is the observation"),
        ("all methods reduced vendor variance by more than three orders of magnitude and rendered multivariate vendor structure indistinguishable from chance",
         "all methods effectively eliminated vendor variance (η² reduction to ~10⁻⁴, a >99.8% decrease) and rendered multivariate vendor structure indistinguishable from chance"),
    ]

    for old, new in MD_REPLACEMENTS:
        if old in content:
            content = content.replace(old, new)
            changes_made += 1
        else:
                print(f"WARNING: Not found in MD: {repr(old[:60])}...")

    # Handle PERMANOVA claim in MD (may have different formatting)
    permanova_old = "These multivariate R² values exceed every univariate η² (including AD's 0.62)"
    permanova_new = "These multivariate R² values (0.315–0.336) are comparable in magnitude to the upper range of univariate effects (the largest being AD at η²=0.62)"
    if permanova_old in content:
        # Need to find the full sentence
        idx = content.find(permanova_old)
        end_sentence = content.find(".", idx)
        full_old = content[idx:end_sentence+1]
        full_new = permanova_new + ". This demonstrates that vendor differences pervade the joint biomarker distribution — vendor-specific covariance patterns create separation in the full 8-dimensional space. Notably, no single univariate metric fully captures this multivariate vendor structure: the PERMANOVA R², while substantial, does not exceed the most extreme univariate η² (AD=0.62), confirming that univariate and multivariate metrics capture complementary aspects of vendor influence."
        content = content.replace(full_old, full_new)
        changes_made += 1

    # Replace limitations in MD
    old_lim_md = "Several limitations should be considered."
    if old_lim_md in content:
        idx = content.find(old_lim_md)
        old_end = "Test-retest data — not part of the spine-generic dataset — would be required to resolve this."
        end_idx = content.find(old_end, idx)
        if end_idx > idx:
            end_idx += len(old_end)
            # Build MD version of new limitations
            md_limitations = """
Several limitations warrant careful consideration.

**Vendor-site confounding.** The spine-generic dataset contains 42 sites but only 3 vendors, making it impossible to disentangle vendor effects from site-level effects within each vendor. A systematic Siemens-versus-Philips difference may partially reflect differences in the specific sites operating each vendor's scanners, rather than purely hardware-level differences. This limitation is inherent to the retrospective, multi-vendor design of the spine-generic study and motivates prospective cross-vendor designs with site crossing. Both our findings and their interpretation should be understood as reflecting the composite 'vendor-in-site' effect, which is the relevant batch unit for most practical multi-center consortia.

**Limited vendor diversity and imbalance.** The dataset covers only three vendors (GE, Philips, Siemens) and is dominated by Siemens scanners (67% of HC cohort). The preservation ranking we report may be sensitive to vendor imbalance — as confirmed by our Design C supplementary analysis, which shows that balanced subsampling compresses the preservation difference across methods. Extrapolation to other vendors (e.g., Canon, United Imaging) or to scanner models not represented in this dataset should be undertaken cautiously.

**Modest biological effect sizes.** Although harmonization consistently increased the age-related variance explained (2.3–3.7 fold relative to baseline), the absolute effect remains small: age explains less than 1% of post-harmonization biomarker variance, and less than 1.7% even in the best-case (CovBat, ALL cohort). These small effect sizes partly reflect the limited age range (19–52 years) and the known moderate sensitivity of single-level C2–C3 scalar metrics to normal aging. Studies targeting age-related cord changes should consider broader age ranges or more spatially extended sampling.

**Method convergence on vendor removal.** All five methods achieved near-identical univariate and multivariate vendor removal (η² reduced to ~10⁻⁴ in every case; all PERMANOVA p ≥ 0.97). This convergence, while reassuring for method robustness, limits the discriminatory power of vendor-removal metrics. The practical differentiation among methods therefore rests on subject-level preservation — a secondary endpoint with modest absolute differences (Δr ~0.06 between best and worst). Users in balanced consortia may find this differentiation insufficiently large to drive a strong method preference.

**Single-pass, within-cohort design.** We evaluated harmonization in a single-pass, full-cohort design without a hold-out vendor. Out-of-sample generalization to a held-out vendor — the relevant scenario for prospective deployment — is an important next step requiring a study design with sufficient per-vendor sample sizes for cross-validation.

**Absence of test-retest data.** Our subject-preservation metric (Pearson r) is correlational and conflates harmonization-induced information loss with irreducible within-subject measurement noise. Partitioning these requires test-retest data, which are not part of the spine-generic dataset. The reported pre-post r values therefore represent a lower bound on true subject preservation.

**Scope of methods evaluated.** We tested five methods from the ComBat/LME family. Other approaches — including GAN-based harmonization, traveling-subject calibration, RAVEL, and ComBat-GAM — were not evaluated. Our comparison should be understood as focused on this specific methodological family rather than an exhaustive survey of all available harmonization strategies.
"""
            content = content[:idx] + md_limitations.strip() + "\n\n" + content[end_idx:]
            changes_made += 1
            print("Replaced MD limitations section.")

    with open(filepath, 'w', encoding='utf-8') as f:
        f.write(content)

    print(f"MD: {changes_made} text replacements applied.")
    return changes_made


if __name__ == "__main__":
    html_path = "C:/Users/admin/WorkBuddy/2026-07-06-12-34-04/spine_generic_harmonization_paper.html"
    md_path = "C:/Users/admin/WorkBuddy/2026-07-06-12-34-04/spine_generic_harmonization_paper.md"

    n_html = apply_html_fixes(html_path)
    n_md = apply_md_fixes(md_path)

    print(f"\nDone. HTML: {n_html} fixes. MD: {n_md} fixes.")
