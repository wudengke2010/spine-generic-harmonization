#!/usr/bin/env python3
"""Cross-check: Paper VENDOR_MAP vs Official Dataset vendor assignments."""

from collections import Counter

# Official manufacturer mapping from participants.tsv
OFFICIAL_MFR = {
    'amu': ('Siemens', 'Verio'),
    'balgrist': ('Siemens', 'Prisma'),
    'barcelona': ('Siemens', 'Prisma-fit'),
    'beijingGE': ('GE', 'MR750'),
    'beijingPrisma': ('Siemens', 'Prisma'),
    'beijingVerio': ('Siemens', 'Verio'),
    'brnoCeitec': ('Siemens', 'Prisma'),
    'brnoUhb': ('GE', 'Signa-PETMR'),
    'cardiff': ('Siemens', 'Prisma'),
    'cmrra': ('Siemens', 'Prisma-fit'),
    'cmrrb': ('Siemens', 'Prisma-fit'),
    'dresden': ('Siemens', 'Prisma'),
    'fslAchieva': ('Philips', 'Achieva'),
    'fslPrisma': ('Siemens', 'Prisma'),
    'geneva': ('Siemens', 'Prisma'),
    'hamburg': ('Siemens', 'Prisma-fit'),
    'juntendo': ('GE', 'MR750w'),
    'mgh': ('Siemens', 'Skyra'),
    'milan': ('Siemens', 'Prisma'),
    'mniPilot': ('Siemens', 'Prisma-fit'),
    'mniS': ('Siemens', 'Prisma-fit'),
    'mountSinai': ('Siemens', 'Skyra'),
    'mpicbs': ('Siemens', 'Prisma-fit'),
    'nottwil': ('Philips', 'Achieva-dStream'),
    'nwu': ('Siemens', 'Prisma-fit'),
    'oxfordFmrib': ('Siemens', 'Prisma'),
    'oxfordOhba': ('Siemens', 'Prisma'),
    'pavia': ('Siemens', 'Skyra'),
    'perform': ('GE', 'MR750'),
    'queensland': ('Siemens', 'Prisma-fit'),
    'sherbrooke': ('Philips', 'Ingenia'),
    'stanford': ('GE', 'MR750'),
    'strasbourg': ('Siemens', 'Verio'),
    'tehranS': ('Siemens', 'Prisma-fit'),
    'tokyo': ('GE', 'MR750w'),
    'tokyoIngenia': ('Philips', 'Ingenia'),
    'tokyoSkyra': ('Siemens', 'Skyra'),
    'ubc': ('Philips', 'Ingenia-ElitionX'),
    'ucdavis': ('Siemens', 'Vida'),
    'ucl': ('Philips', 'Achieva-dStream'),
    'unf': ('Siemens', 'Prisma-fit'),
    'vallHebron': ('Siemens', 'TrioTim'),
    'vuiisAchieva': ('Philips', 'Achieva-dStream'),
    'vuiisIngenia': ('Philips', 'Ingenia-ElitionX'),
}

# Paper VENDOR_MAP from analyze_spine_generic.py
PAPER_VENDOR = {
    'amu': 'Siemens', 'balgrist': 'Siemens', 'barcelona': 'Siemens',
    'beijingPrisma': 'Siemens', 'beijingVerio': 'Siemens', 'fslPrisma': 'Siemens',
    'tokyoSkyra': 'Siemens',
    'beijingGE': 'GE', 'juntendo': 'GE', 'tokyo': 'GE',
    'fslAchieva': 'Philips', 'tokyoIngenia': 'Philips', 'vuiisAchieva': 'Philips',
    'vuiisIngenia': 'Philips',
    'brnoCeitec': 'Siemens', 'brnoUhb': 'Siemens', 'cardiff': 'Siemens',
    'cmrra': 'Siemens', 'cmrrb': 'Siemens', 'dresden': 'Siemens',
    'geneva': 'Siemens', 'hamburg': 'Siemens', 'mgh': 'Siemens',
    'milan': 'Siemens', 'mniPilot': 'Siemens', 'mniS': 'Siemens',
    'mountSinai': 'Siemens', 'mpicbs': 'Siemens', 'nottwil': 'Siemens',
    'nwu': 'Siemens', 'oxfordFmrib': 'Siemens', 'oxfordOhba': 'Siemens',
    'pavia': 'Siemens', 'perform': 'Siemens', 'queensland': 'Siemens',
    'sherbrooke': 'Siemens', 'stanford': 'Siemens', 'strasbourg': 'Siemens',
    'tehranS': 'Siemens', 'ucdavis': 'Siemens', 'ucl': 'Siemens',
    'unf': 'Siemens', 'vallHebron': 'Siemens', 'ubc': 'Siemens',
}

SUBJECT_COUNTS = {
    'amu':5,'balgrist':6,'barcelona':6,'beijingGE':4,'beijingPrisma':5,'beijingVerio':4,
    'brnoCeitec':6,'brnoUhb':8,'cardiff':6,'cmrra':6,'cmrrb':7,'dresden':2,
    'fslAchieva':6,'fslPrisma':6,'geneva':6,'hamburg':6,'juntendo':6,
    'mgh':6,'milan':7,'mniPilot':1,'mniS':9,'mountSinai':6,'mpicbs':6,
    'nottwil':6,'nwu':6,'oxfordFmrib':11,'oxfordOhba':5,'pavia':6,'perform':6,
    'queensland':6,'sherbrooke':7,'stanford':6,'strasbourg':6,'tehranS':6,
    'tokyo':7,'tokyoIngenia':7,'tokyoSkyra':7,'ubc':6,'ucdavis':7,'ucl':6,
    'unf':7,'vallHebron':7,'vuiisAchieva':6,'vuiisIngenia':6,
}

print("=" * 70)
print("VENDOR MAPPING CROSS-CHECK: Paper VENDOR_MAP vs Official Dataset")
print("=" * 70)
print()

errors = []
paper_vendor_total = Counter()
official_vendor_total = Counter()

for site, official in OFFICIAL_MFR.items():
    paper = PAPER_VENDOR.get(site, 'MISSING')
    count = SUBJECT_COUNTS.get(site, 0)
    
    official_vendor_total[official[0]] += count
    paper_vendor_total[paper] += count
    
    if paper != official[0]:
        errors.append((site, paper, official[0], official[1], count))

print(f"Correct mappings: {len(OFFICIAL_MFR) - len(errors)}/{len(OFFICIAL_MFR)} sites")
print(f"INCORRECT mappings: {len(errors)}/{len(OFFICIAL_MFR)} sites")
print()

if errors:
    print("*** ERRORS FOUND IN analyze_spine_generic.py VENDOR_MAP ***")
    print(f"  {'Site':<20} {'Script says':<12} {'Actually':<12} {'Model':<22} {'N subj':>6}")
    print("  " + "-" * 72)
    for site, paper, actual, model, count in errors:
        print(f"  {site:<20} {paper:<12} {actual:<12} {model:<22} {count:>6}")
    print()

print("--- Vendor totals (ALL 267 subjects from local data) ---")
print(f"  {'Vendor':<12} {'Script (WRONG)':>16} {'Official (CORRECT)':>20}")
print("  " + "-" * 50)
for v in ['Siemens', 'GE', 'Philips']:
    print(f"  {v:<12} {paper_vendor_total.get(v,0):>16} {official_vendor_total.get(v,0):>20}")

print()
print("--- Cross-check with paper claims ---")
print(f"  Paper HC cohort claims:  GE=28,  Philips=36,  Siemens=139  (Total 203)")
print(f"  Paper ALL complete-case: GE=33,  Philips=47,  Siemens=166  (Total 246)")
print(f"  Valosek 2024 HC:         GE=28,  Philips=36,  Siemens=139  (Total 203)")
print()

# Compute what HC vendor distribution should be if paper numbers are correct
# Our paper says 64 pathology out of 267 total = 203 HC
# From official data: Siemens 180, GE 37, Philips 50 total
# We need to verify how many HC per vendor

# Count pathology per vendor from the participants.tsv data (using the official pathology labels)
# This requires cross-referencing with the datase
print("*** NOTE: Pathology assignment in paper may differ from dataset 'pathology' column.")
print("*** The dataset 'pathology' column labels subjects as HC/MildCompression/DCM")
print("*** But the paper reclassifies based on presence/absence of spinal cord compression.")
print("*** Verification of exact HC counts per vendor requires the actual analysis data.")
