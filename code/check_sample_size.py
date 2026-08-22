"""
Check actual sample size of 7-biomarker analysis.
"""
import pandas as pd
import numpy as np
import os, glob

# ── Load data ──
perlevel = pd.read_csv("E:/boshi/spine-generic-multi-subject/perlevel_metrics.csv")
participants = pd.read_csv("E:/boshi/spine-generic-multi-subject/participants.tsv", sep="\t")

print("=" * 70)
print("RAW DATA")
print("=" * 70)
print(f"perlevel:  {perlevel.shape[0]} rows, {perlevel['subject'].nunique()} subjects")
print(f"participants: {participants.shape[0]} subjects")
print(f"levels: {sorted(perlevel['level'].unique())}")

# ── Filter C2-C3 ──
c23 = perlevel[perlevel["level"].isin(["C2", "C3"])].copy()
print(f"\nC2-C3: {c23.shape[0]} rows, {c23['subject'].nunique()} subjects")

# ── Aggregate mean per subject at C2-C3 ──
BIOMARKERS_7 = ["CSA", "FA", "MD", "AD", "RD", "MTR", "MTsat"]
agg = c23.groupby("subject")[BIOMARKERS_7].mean()
print(f"\nC2-C3 mean aggregated: {agg.shape[0]} subjects")

# NaN check per biomarker
print("\nNaN counts per biomarker:")
for bm in BIOMARKERS_7:
    n = int(agg[bm].isna().sum())
    print(f"  {bm}: {n} NaN  →  {263 - n} valid")

# Complete cases (all 7)
complete_7 = agg.dropna()
print(f"\nComplete 7 biomarkers: {complete_7.shape[0]} subjects")

# ── Add site info from perlevel ──
site_info = c23.groupby("subject")["site"].first()
complete_7_w_site = complete_7.copy()
complete_7_w_site["site"] = site_info

# ── Vendor mapping ──
VENDOR_MAP = {
    "amu": "Siemens", "balgrist": "Siemens", "barcelona": "Siemens",
    "beijingGE": "GE", "beijingPrisma": "Siemens", "beijingVerio": "Siemens",
    "brnoCeitec": "Siemens", "brnoUhb": "GE", "cardiff": "Siemens",
    "cmrra": "Siemens", "cmrrb": "Siemens", "dresden": "Siemens",
    "fslAchieva": "Philips", "fslPrisma": "Siemens", "geneva": "Siemens",
    "hamburg": "Siemens", "juntendo750w": "GE", "mgh": "Siemens",
    "milan": "Siemens", "mniPilot": "Siemens", "mniS": "Siemens",
    "mountSinai": "Siemens", "mpicbs": "Siemens", "nottwil": "Philips",
    "nwu": "Siemens", "oxfordFmrib": "Siemens", "oxfordOhba": "Siemens",
    "pavia": "Siemens", "perform": "GE", "queensland": "Siemens",
    "sherbrooke": "Philips", "stanford": "GE", "strasbourg": "Siemens",
    "tehranS": "Siemens", "tokyo750w": "GE", "tokyoIngenia": "Philips",
    "tokyoSkyra": "Siemens", "ubc": "Philips", "ucdavis": "Siemens",
    "ucl": "Philips", "unf": "Siemens", "vallHebron": "Siemens",
    "vuiisAchieva": "Philips", "vuiisIngenia": "Philips",
}

c7 = complete_7_w_site.reset_index()
c7["vendor"] = c7["site"].map(VENDOR_MAP)
unknown = c7[c7["vendor"].isna()]
if len(unknown) > 0:
    print(f"\n⚠ Unknown vendors: {unknown['site'].unique()}")

# ── Merge with pathology ──
p_info = participants[["participant_id", "pathology"]].copy()
merged = c7.merge(p_info, left_on="subject", right_on="participant_id", how="left")
hc = merged[merged["pathology"] == "HC"]
non_hc = merged[merged["pathology"] != "HC"]

print("\n" + "=" * 70)
print("7-BIOMARKER RESULTS")
print("=" * 70)
print(f"ALL subjects:  {len(c7)}")
print(f"  HC:          {len(hc)}")
print(f"  Non-HC:      {len(non_hc)}")
print(f"Sites:         {c7['site'].nunique()}")
print(f"\nVendor distribution (ALL):")
for v in ["Siemens", "GE", "Philips", "Unknown"]:
    n = int((c7["vendor"] == v).sum())
    if n > 0:
        print(f"  {v}: {n}")
print(f"\nVendor distribution (HC):")
for v in ["Siemens", "GE", "Philips", "Unknown"]:
    n = int((hc["vendor"] == v).sum())
    if n > 0:
        print(f"  {v}: {n}")

sites_per_v = c7.groupby("vendor")["site"].nunique()
print(f"\nSites per vendor:")
for v, n in sites_per_v.items():
    print(f"  {v}: {n} sites")

# ── 8-biomarker (with T2*) ──
t2s_dir = "E:/qm_harmonization_paper/derivatives/csa_t2star"
t2s_files = glob.glob(os.path.join(t2s_dir, "*_csa_t2star.csv"))
t2s_csa = {}
for f in t2s_files:
    df = pd.read_csv(f)
    if "GM_CSA_mm2" in df.columns and len(df) > 0:
        sub = str(df["participant_id"].iloc[0])
        t2s_csa[sub] = float(df["GM_CSA_mm2"].iloc[0])
print(f"\nT2* GM CSA loaded: {len(t2s_csa)} subjects")

c7_df = c7.copy()
c7_df["T2s_GM_CSA"] = c7_df["subject"].map(t2s_csa)
c8 = c7_df.dropna(subset=["T2s_GM_CSA"])
c8_hc = c8.merge(p_info, left_on="subject", right_on="participant_id", how="left")
c8_hc = c8_hc[c8_hc["pathology"] == "HC"]

print("\n" + "=" * 70)
print("8-BIOMARKER RESULTS (with T2*)")
print("=" * 70)
print(f"ALL subjects:  {len(c8)}")
print(f"  HC:          {len(c8_hc)}")
print(f"Sites:         {c8['site'].nunique()}")
print(f"\nVendor distribution (ALL):")
for v in ["Siemens", "GE", "Philips"]:
    n = int((c8["vendor"] == v).sum())
    if n > 0:
        print(f"  {v}: {n}")
print(f"\nVendor distribution (HC):")
for v in ["Siemens", "GE", "Philips"]:
    n = int((c8_hc["vendor"] == v).sum())
    if n > 0:
        print(f"  {v}: {n}")

# ── Missing sites ──
sites_in_data = set(c7["site"].unique())
all_vendor_sites = set(VENDOR_MAP.keys())
sites_missing = all_vendor_sites - sites_in_data
print(f"\nSites in VENDOR_MAP but NOT in data: {sorted(sites_missing) if sites_missing else 'NONE'}")

# ── SUMMARY ──
print("\n" + "=" * 70)
print("PAPER vs ACTUAL")
print("=" * 70)
print(f"Paper claims:         267 subjects, 44 sites")
print(f"7-bio (C2-C3, ALL):   {len(c7)} subjects, {c7['site'].nunique()} sites")
print(f"7-bio (C2-C3, HC):    {len(hc)} subjects")
print(f"8-bio (C2-C3, ALL):   {len(c8)} subjects, {c8['site'].nunique()} sites")
print(f"8-bio (C2-C3, HC):    {len(c8_hc)} subjects")

# perlevel covers 4 levels (C2-C5), let's check if using all levels changes count
c_all_levels = perlevel.groupby("subject")[BIOMARKERS_7].mean().dropna()
print(f"\n7-bio (ALL levels C2-C5, complete): {c_all_levels.shape[0]} subjects")
