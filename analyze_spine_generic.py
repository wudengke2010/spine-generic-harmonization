#!/usr/bin/env python3
"""
Spine-Generic QM Harmonization Paper: Data Collection & Statistical Analysis
============================================================================
Reads CSA (Cross-Sectional Area) derivatives from spine-generic public database,
performs multi-site harmonization analysis, and generates publication figures.
"""

# ============================================================
# DEPRECATED — 2026-07-08
# This script predates the 8-biomarker harmonization analysis
# (spine_generic_harmonization_paper.tex). It performs basic
# T2*/T2w descriptive statistics only.
#
# VENDOR_MAP was corrected on 2026-07-08 to fix 7 site mappings
# that were misclassified as Siemens:
#   brnoUhb:   Siemens → GE      (Signa-PETMR)
#   nottwil:   Siemens → Philips (Achieva-dStream)
#   perform:   Siemens → GE      (MR750)
#   sherbrooke: Siemens → Philips (Ingenia)
#   stanford:  Siemens → GE      (MR750)
#   ubc:       GE      → Philips (Ingenia-ElitionX)
#   ucl:       Siemens → Philips (Achieva-dStream)
#
# All results from this script are superseded by the main
# harmonization paper. Do not use for publication figures.
# ============================================================

import os
import re
import glob
import json
import warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
from scipy import stats
from scipy.stats import shapiro, kruskal, mannwhitneyu, levene, f_oneway
from statsmodels.stats.multicomp import pairwise_tukeyhsd
import statsmodels.api as sm
from statsmodels.formula.api import ols

warnings.filterwarnings('ignore')

# ============================================================
# Configuration
# ============================================================
DERIVATIVES_DIR = "E:/qm_harmonization_paper/derivatives"
OUTPUT_DIR = "C:/Users/admin/WorkBuddy/2026-07-06-12-34-04/results"
FIGURES_DIR = os.path.join(OUTPUT_DIR, "figures")
DATA_DIR = os.path.join(OUTPUT_DIR, "data")
os.makedirs(FIGURES_DIR, exist_ok=True)
os.makedirs(DATA_DIR, exist_ok=True)

# Vendor mapping based on spine-generic dataset
# Sites with scanner model in name
VENDOR_MAP = {
    # Siemens sites
    'beijingPrisma': 'Siemens', 'beijingVerio': 'Siemens', 'fslPrisma': 'Siemens',
    'tokyoSkyra': 'Siemens',
    # GE sites
    'beijingGE': 'GE', 'juntendo750w': 'GE', 'tokyo750w': 'GE',
    # Philips sites
    'fslAchieva': 'Philips', 'tokyoIngenia': 'Philips', 'vuiisAchieva': 'Philips',
    'vuiisIngenia': 'Philips',
    # Siemens-only sites (from spine-generic publication)
    'amu': 'Siemens', 'balgrist': 'Siemens', 'barcelona': 'Siemens',
    'brnoCeitec': 'Siemens', 'cardiff': 'Siemens',
    'cmrra': 'Siemens', 'cmrrb': 'Siemens', 'dresden': 'Siemens',
    'geneva': 'Siemens', 'hamburg': 'Siemens', 'mgh': 'Siemens',
    'milan': 'Siemens', 'mniPilot': 'Siemens', 'mniS': 'Siemens',
    'mountSinai': 'Siemens', 'mpicbs': 'Siemens',
    'nwu': 'Siemens', 'oxfordFmrib': 'Siemens', 'oxfordOhba': 'Siemens',
    'pavia': 'Siemens', 'queensland': 'Siemens',
    'strasbourg': 'Siemens',
    'tehranS': 'Siemens', 'ucdavis': 'Siemens',
    'unf': 'Siemens', 'vallHebron': 'Siemens',
    # GE sites (additional, not in name)
    'brnoUhb': 'GE', 'perform': 'GE', 'stanford': 'GE',
    # Philips sites (additional, not in name)
    'nottwil': 'Philips', 'sherbrooke': 'Philips', 'ucl': 'Philips',
    'ubc': 'Philips',
}

# Field strength mapping (from spine-generic: most are 3T, some 1.5T)
FIELD_STRENGTH_MAP = {
    'balgrist': '3T', 'brnoUhb': '1.5T', 'tehranS': '1.5T',
    # Default 3T for all others
}

# Region mapping (from spine-generic: C2-C3, C4-C5, C6-C7)
# We'll analyze by vertebral level

plt.rcParams.update({
    'font.size': 11,
    'axes.titlesize': 13,
    'axes.labelsize': 12,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'legend.fontsize': 10,
    'figure.dpi': 150,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'font.family': 'sans-serif',
})

# ============================================================
# 1. Data Collection
# ============================================================
def parse_subject_id(sub_id):
    """Parse subject ID to extract site and vendor."""
    # Remove 'sub-' prefix
    name = sub_id.replace('sub-', '')
    # Try to match site name (letters only, before numbers)
    match = re.match(r'^([a-zA-Z]+)', name)
    if match:
        site_key = match.group(1)
    else:
        site_key = name
    
    # Check if site_key is in vendor map, if not try full name without numbers
    if site_key not in VENDOR_MAP:
        # Try removing trailing digits more carefully
        site_full = re.sub(r'\d+$', '', name)
        site_key = site_full if site_full in VENDOR_MAP else site_key
    
    vendor = VENDOR_MAP.get(site_key, 'Unknown')
    field_strength = FIELD_STRENGTH_MAP.get(site_key, '3T')
    
    return site_key, vendor, field_strength


def collect_t2star_data():
    """Collect all T2* CSA summary data."""
    files = sorted(glob.glob(os.path.join(DERIVATIVES_DIR, "csa_t2star", "*_csa_t2star.csv")))
    records = []
    
    for f in files:
        df = pd.read_csv(f)
        if df.empty:
            continue
        row = df.iloc[0]
        sub_id = row['participant_id']
        site, vendor, field = parse_subject_id(sub_id)
        
        records.append({
            'subject': sub_id,
            'site': site,
            'vendor': vendor,
            'field_strength': field,
            'n_slices_t2star': row.get('n_valid_slices', np.nan),
            'SC_CSA_t2star': row.get('SC_CSA_mm2', np.nan),
            'GM_CSA_t2star': row.get('GM_CSA_mm2', np.nan),
            'WM_CSA_t2star': row.get('WM_CSA_mm2', np.nan),
            'GM_WM_ratio_t2star': row.get('GM_WM_ratio', np.nan),
            'GM_SC_frac_t2star': row.get('GM_SC_frac', np.nan),
            'voxel_dx': row.get('voxel_dx', np.nan),
            'voxel_dy': row.get('voxel_dy', np.nan),
            'voxel_dz': row.get('voxel_dz', np.nan),
        })
    
    return pd.DataFrame(records)


def collect_t2w_data():
    """Collect all T2w CSA data."""
    files = sorted(glob.glob(os.path.join(DERIVATIVES_DIR, "csa_t2w", "*_csa.csv")))
    records = []
    
    for f in files:
        fname = os.path.basename(f)
        sub_id = fname.replace('_csa.csv', '')
        site, vendor, field = parse_subject_id(sub_id)
        
        df = pd.read_csv(f)
        if df.empty:
            continue
        
        # The T2w CSA file has detailed per-slice data
        # Extract mean CSA across valid slices
        area_col = 'MEAN(area)'
        if area_col not in df.columns:
            continue
        
        areas = pd.to_numeric(df[area_col], errors='coerce').dropna()
        if areas.empty:
            continue
        
        # Parse vertebral levels if available
        vert_levels = df.get('VertLevel', pd.Series(dtype=str)).dropna()
        
        records.append({
            'subject': sub_id,
            'site': site,
            'vendor': vendor,
            'field_strength': field,
            'SC_CSA_t2w': areas.mean(),
            'SC_CSA_std_t2w': areas.std(),
            'SC_CSA_min_t2w': areas.min(),
            'SC_CSA_max_t2w': areas.max(),
            'n_slices_t2w': len(areas),
            'n_vert_levels': len(vert_levels.unique()) if not vert_levels.empty else 0,
        })
    
    return pd.DataFrame(records)


def collect_t2w_by_vertebral_level():
    """Collect T2w CSA by vertebral level."""
    files = sorted(glob.glob(os.path.join(DERIVATIVES_DIR, "csa_t2w", "*_csa.csv")))
    records = []
    
    for f in files:
        fname = os.path.basename(f)
        sub_id = fname.replace('_csa.csv', '')
        site, vendor, field = parse_subject_id(sub_id)
        
        df = pd.read_csv(f)
        if df.empty or 'VertLevel' not in df.columns:
            continue
        
        df['area'] = pd.to_numeric(df['MEAN(area)'], errors='coerce')
        
        for vl, group in df.groupby('VertLevel'):
            if pd.isna(vl) or vl == '':
                continue
            areas = group['area'].dropna()
            if areas.empty:
                continue
            records.append({
                'subject': sub_id,
                'site': site,
                'vendor': vendor,
                'field_strength': field,
                'vert_level': vl,
                'SC_CSA': areas.mean(),
            })
    
    return pd.DataFrame(records)


def collect_t2star_detailed():
    """Collect detailed T2* per-component data (GM, SC, WM separately)."""
    files = sorted(glob.glob(os.path.join(DERIVATIVES_DIR, "csa_t2star", "*_csa_gm.csv")))
    records = []
    
    for f in files:
        fname = os.path.basename(f)
        sub_id = fname.replace('_csa_gm.csv', '')
        site, vendor, field = parse_subject_id(sub_id)
        
        df = pd.read_csv(f)
        if df.empty:
            continue
        
        area_col = 'MEAN(area)'
        if area_col not in df.columns:
            continue
        
        areas = pd.to_numeric(df[area_col], errors='coerce').dropna()
        if areas.empty:
            continue
        
        records.append({
            'subject': sub_id,
            'site': site,
            'vendor': vendor,
            'field_strength': field,
            'GM_CSA_detailed': areas.mean(),
        })
    
    return pd.DataFrame(records)


# ============================================================
# 2. Statistical Analysis
# ============================================================
def descriptive_stats(df, columns, group_col=None):
    """Compute descriptive statistics."""
    if group_col:
        stats_df = df.groupby(group_col)[columns].agg(['count', 'mean', 'std', 'min', 'max']).round(2)
    else:
        stats_list = []
        for col in columns:
            if col in df.columns:
                vals = df[col].dropna()
                stats_list.append({
                    'variable': col,
                    'n': len(vals),
                    'mean': round(vals.mean(), 2),
                    'std': round(vals.std(), 2),
                    'min': round(vals.min(), 2),
                    'max': round(vals.max(), 2),
                    'median': round(vals.median(), 2),
                    'IQR': round(vals.quantile(0.75) - vals.quantile(0.25), 2),
                })
        stats_df = pd.DataFrame(stats_list)
    return stats_df


def coefficient_of_variation(df, value_col, group_col):
    """Calculate CV within each group and overall."""
    cv_results = []
    for group, data in df.groupby(group_col):
        vals = data[value_col].dropna()
        if len(vals) >= 2:
            cv = (vals.std() / vals.mean()) * 100
            cv_results.append({
                group_col: group,
                'n': len(vals),
                'mean': round(vals.mean(), 2),
                'std': round(vals.std(), 2),
                'CV(%)': round(cv, 2),
            })
    return pd.DataFrame(cv_results).sort_values('CV(%)', ascending=False)


def perform_kruskal_wallis(df, value_col, group_col):
    """Kruskal-Wallis test across groups."""
    groups = []
    group_names = []
    for name, group in df.groupby(group_col):
        vals = group[value_col].dropna()
        if len(vals) >= 3:
            groups.append(vals.values)
            group_names.append(name)
    
    if len(groups) >= 2:
        stat, p = kruskal(*groups)
        return {'test': 'Kruskal-Wallis', 'statistic': round(stat, 4), 'p_value': p, 'n_groups': len(groups)}
    return None


def perform_anova(df, value_col, group_col):
    """One-way ANOVA across groups."""
    groups = []
    for name, group in df.groupby(group_col):
        vals = group[value_col].dropna()
        if len(vals) >= 3:
            groups.append(vals.values)
    
    if len(groups) >= 2:
        stat, p = f_oneway(*groups)
        return {'test': 'One-way ANOVA', 'statistic': round(stat, 4), 'p_value': p, 'n_groups': len(groups)}
    return None


def calculate_icc(df, value_col, subject_col, rater_col):
    """Calculate ICC(2,1) for inter-rater reliability."""
    # Pivot to wide format
    pivot = df.pivot_table(index=subject_col, columns=rater_col, values=value_col)
    pivot = pivot.dropna()
    
    if pivot.shape[1] < 2 or pivot.shape[0] < 3:
        return None
    
    n = pivot.shape[0]
    k = pivot.shape[1]
    
    # Two-way random, single measure ICC(2,1)
    grand_mean = pivot.values.mean()
    ss_between = n * ((pivot.mean(axis=1) - grand_mean) ** 2).sum()
    ss_within = ((pivot - pivot.mean(axis=0)) ** 2).values.sum()
    ss_raters = n * ((pivot.mean(axis=0) - grand_mean) ** 2).sum()
    
    ms_between = ss_between / (n - 1)
    ms_within = ss_within / ((n - 1) * (k - 1))
    ms_raters = ss_raters / (k - 1)
    
    icc = (ms_between - ms_within) / (ms_between + (k - 1) * ms_within + k * (ms_raters - ms_within) / n)
    
    return {'ICC(2,1)': round(icc, 4), 'n_subjects': n, 'n_raters': k}


# ============================================================
# 3. Figure Generation
# ============================================================
def plot_site_variability_t2star(df, savepath):
    """Box plot of T2* CSA metrics across sites."""
    fig, axes = plt.subplots(2, 2, figsize=(18, 14))
    
    metrics = [
        ('SC_CSA_t2star', 'Spinal Cord CSA (mm²)', 'A'),
        ('GM_CSA_t2star', 'Gray Matter CSA (mm²)', 'B'),
        ('WM_CSA_t2star', 'White Matter CSA (mm²)', 'C'),
        ('GM_SC_frac_t2star', 'GM/SC Fraction', 'D'),
    ]
    
    for ax, (col, label, panel) in zip(axes.flat, metrics):
        if col not in df.columns:
            ax.set_visible(False)
            continue
        
        plot_df = df[['site', col]].dropna()
        site_order = plot_df.groupby('site')[col].median().sort_values().index
        
        colors = []
        for s in site_order:
            v = VENDOR_MAP.get(s, 'Unknown')
            colors.append({'Siemens': '#1f77b4', 'GE': '#d62728', 'Philips': '#2ca02c'}.get(v, '#7f7f7f'))
        
        sns.boxplot(data=plot_df, x='site', y=col, order=site_order, ax=ax,
                    palette=colors, width=0.6, fliersize=3)
        ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right', fontsize=8)
        ax.set_xlabel('')
        ax.set_ylabel(label)
        ax.set_title(f'{panel}. {label} by Site', fontweight='bold', loc='left')
        
        # Add overall mean line
        overall_mean = plot_df[col].mean()
        ax.axhline(y=overall_mean, color='red', linestyle='--', alpha=0.5, linewidth=1)
    
    # Add legend
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='#1f77b4', label='Siemens'),
        Patch(facecolor='#d62728', label='GE'),
        Patch(facecolor='#2ca02c', label='Philips'),
    ]
    fig.legend(handles=legend_elements, loc='lower center', ncol=3, fontsize=11,
               bbox_to_anchor=(0.5, -0.02))
    
    plt.suptitle('T2* CSA Metrics Across Acquisition Sites', fontsize=15, fontweight='bold', y=1.01)
    plt.tight_layout()
    plt.savefig(savepath, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {savepath}")


def plot_vendor_comparison(df, savepath):
    """Compare CSA metrics across vendors."""
    fig, axes = plt.subplots(1, 3, figsize=(16, 6))
    
    metrics = [
        ('SC_CSA_t2star', 'Spinal Cord CSA (mm²)'),
        ('GM_CSA_t2star', 'Gray Matter CSA (mm²)'),
        ('GM_SC_frac_t2star', 'GM/SC Fraction'),
    ]
    
    for ax, (col, label) in zip(axes, metrics):
        if col not in df.columns:
            ax.set_visible(False)
            continue
        plot_df = df[['vendor', col]].dropna()
        plot_df = plot_df[plot_df['vendor'] != 'Unknown']
        
        sns.violinplot(data=plot_df, x='vendor', y=col, ax=ax,
                       palette={'Siemens': '#1f77b4', 'GE': '#d62728', 'Philips': '#2ca02c'},
                       inner='box', cut=0)
        sns.stripplot(data=plot_df, x='vendor', y=col, ax=ax,
                      color='black', alpha=0.3, size=3, jitter=True)
        
        ax.set_xlabel('Vendor')
        ax.set_ylabel(label)
        ax.set_title(label, fontweight='bold')
    
    plt.suptitle('T2* CSA Metrics: Vendor Comparison', fontsize=14, fontweight='bold', y=1.02)
    plt.tight_layout()
    plt.savefig(savepath, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {savepath}")


def plot_t2w_vs_t2star(merged_df, savepath):
    """Bland-Altman and correlation plot: T2w vs T2* SC CSA."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    plot_df = merged_df[['SC_CSA_t2w', 'SC_CSA_t2star', 'vendor']].dropna()
    plot_df = plot_df[plot_df['vendor'] != 'Unknown']
    
    # Correlation plot
    ax = axes[0]
    for vendor, color in [('Siemens', '#1f77b4'), ('GE', '#d62728'), ('Philips', '#2ca02c')]:
        subset = plot_df[plot_df['vendor'] == vendor]
        ax.scatter(subset['SC_CSA_t2w'], subset['SC_CSA_t2star'], 
                   c=color, label=vendor, alpha=0.7, s=40, edgecolors='white', linewidths=0.5)
    
    # Add diagonal line
    lims = [min(plot_df['SC_CSA_t2w'].min(), plot_df['SC_CSA_t2star'].min()) - 2,
            max(plot_df['SC_CSA_t2w'].max(), plot_df['SC_CSA_t2star'].max()) + 2]
    ax.plot(lims, lims, 'k--', alpha=0.5, linewidth=1)
    
    r, p = stats.pearsonr(plot_df['SC_CSA_t2w'], plot_df['SC_CSA_t2star'])
    ax.set_xlabel('T2w SC CSA (mm²)')
    ax.set_ylabel('T2* SC CSA (mm²)')
    ax.set_title(f'A. T2w vs T2* SC CSA\nr = {r:.3f}, p < 0.001' if p < 0.001 else f'A. T2w vs T2* SC CSA\nr = {r:.3f}, p = {p:.3f}',
                 fontweight='bold', loc='left')
    ax.legend()
    
    # Bland-Altman plot
    ax = axes[1]
    plot_df['mean_csa'] = (plot_df['SC_CSA_t2w'] + plot_df['SC_CSA_t2star']) / 2
    plot_df['diff'] = plot_df['SC_CSA_t2w'] - plot_df['SC_CSA_t2star']
    
    for vendor, color in [('Siemens', '#1f77b4'), ('GE', '#d62728'), ('Philips', '#2ca02c')]:
        subset = plot_df[plot_df['vendor'] == vendor]
        ax.scatter(subset['mean_csa'], subset['diff'], 
                   c=color, label=vendor, alpha=0.7, s=40, edgecolors='white', linewidths=0.5)
    
    mean_diff = plot_df['diff'].mean()
    sd_diff = plot_df['diff'].std()
    ax.axhline(y=mean_diff, color='red', linestyle='-', linewidth=1.5, label=f'Mean bias: {mean_diff:.2f}')
    ax.axhline(y=mean_diff + 1.96*sd_diff, color='gray', linestyle='--', linewidth=1, 
               label=f'±1.96 SD: [{mean_diff+1.96*sd_diff:.2f}, {mean_diff-1.96*sd_diff:.2f}]')
    ax.axhline(y=mean_diff - 1.96*sd_diff, color='gray', linestyle='--', linewidth=1)
    
    ax.set_xlabel('Mean SC CSA (mm²)')
    ax.set_ylabel('T2w - T2* Difference (mm²)')
    ax.set_title('B. Bland-Altman: T2w vs T2*', fontweight='bold', loc='left')
    ax.legend(fontsize=9)
    
    plt.tight_layout()
    plt.savefig(savepath, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {savepath}")


def plot_cv_by_site(df, savepath):
    """Bar plot of coefficient of variation by site."""
    fig, axes = plt.subplots(1, 2, figsize=(16, 7))
    
    # T2* SC CSA CV
    ax = axes[0]
    cv_data = coefficient_of_variation(df, 'SC_CSA_t2star', 'site')
    if not cv_data.empty:
        colors = [{'Siemens': '#1f77b4', 'GE': '#d62728', 'Philips': '#2ca02c'}.get(VENDOR_MAP.get(s, 'Unknown'), '#7f7f7f') 
                  for s in cv_data['site']]
        bars = ax.barh(range(len(cv_data)), cv_data['CV(%)'], color=colors, edgecolor='white', linewidth=0.5)
        ax.set_yticks(range(len(cv_data)))
        ax.set_yticklabels(cv_data['site'], fontsize=9)
        ax.set_xlabel('Coefficient of Variation (%)')
        ax.set_title('A. T2* SC CSA: Within-Site CV', fontweight='bold', loc='left')
        ax.axvline(x=cv_data['CV(%)'].mean(), color='red', linestyle='--', linewidth=1, 
                   label=f'Mean CV: {cv_data["CV(%)"].mean():.2f}%')
        ax.legend()
    
    # T2* GM CSA CV
    ax = axes[1]
    cv_data_gm = coefficient_of_variation(df, 'GM_CSA_t2star', 'site')
    if not cv_data_gm.empty:
        colors = [{'Siemens': '#1f77b4', 'GE': '#d62728', 'Philips': '#2ca02c'}.get(VENDOR_MAP.get(s, 'Unknown'), '#7f7f7f') 
                  for s in cv_data_gm['site']]
        bars = ax.barh(range(len(cv_data_gm)), cv_data_gm['CV(%)'], color=colors, edgecolor='white', linewidth=0.5)
        ax.set_yticks(range(len(cv_data_gm)))
        ax.set_yticklabels(cv_data_gm['site'], fontsize=9)
        ax.set_xlabel('Coefficient of Variation (%)')
        ax.set_title('B. T2* GM CSA: Within-Site CV', fontweight='bold', loc='left')
        ax.axvline(x=cv_data_gm['CV(%)'].mean(), color='red', linestyle='--', linewidth=1,
                   label=f'Mean CV: {cv_data_gm["CV(%)"].mean():.2f}%')
        ax.legend()
    
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='#1f77b4', label='Siemens'),
        Patch(facecolor='#d62728', label='GE'),
        Patch(facecolor='#2ca02c', label='Philips'),
    ]
    fig.legend(handles=legend_elements, loc='lower center', ncol=3, fontsize=10,
               bbox_to_anchor=(0.5, -0.02))
    
    plt.suptitle('Within-Site Variability (Coefficient of Variation)', fontsize=14, fontweight='bold', y=1.01)
    plt.tight_layout()
    plt.savefig(savepath, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {savepath}")


def plot_vertebral_level_distribution(df_vl, savepath):
    """CSA distribution by vertebral level."""
    fig, ax = plt.subplots(figsize=(12, 7))
    
    if df_vl.empty:
        plt.close()
        return
    
    # Sort vertebral levels
    vl_order = sorted(df_vl['vert_level'].unique())
    
    sns.boxplot(data=df_vl, x='vert_level', y='SC_CSA', order=vl_order, ax=ax,
                palette='Set2', width=0.6, fliersize=3)
    sns.stripplot(data=df_vl, x='vert_level', y='SC_CSA', order=vl_order, ax=ax,
                  color='black', alpha=0.2, size=2, jitter=True)
    
    ax.set_xlabel('Vertebral Level')
    ax.set_ylabel('Spinal Cord CSA (mm²)')
    ax.set_title('T2w SC CSA Distribution by Vertebral Level', fontweight='bold')
    
    plt.tight_layout()
    plt.savefig(savepath, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {savepath}")


def plot_correlation_matrix(df, savepath):
    """Correlation matrix of all CSA metrics."""
    cols = ['SC_CSA_t2star', 'GM_CSA_t2star', 'WM_CSA_t2star', 
            'GM_WM_ratio_t2star', 'GM_SC_frac_t2star', 'SC_CSA_t2w']
    cols_present = [c for c in cols if c in df.columns]
    
    if len(cols_present) < 2:
        return
    
    corr_df = df[cols_present].dropna()
    if corr_df.empty:
        return
    
    corr = corr_df.corr()
    
    fig, ax = plt.subplots(figsize=(10, 8))
    mask = np.triu(np.ones_like(corr, dtype=bool), k=1)
    sns.heatmap(corr, mask=mask, annot=True, fmt='.2f', cmap='RdBu_r', 
                center=0, vmin=-1, vmax=1, square=True, ax=ax,
                linewidths=0.5, cbar_kws={'shrink': 0.8})
    
    labels = [c.replace('_t2star', ' (T2*)').replace('_t2w', ' (T2w)').replace('_', ' ') for c in cols_present]
    ax.set_xticklabels(labels, rotation=45, ha='right', fontsize=9)
    ax.set_yticklabels(labels, rotation=0, fontsize=9)
    ax.set_title('Correlation Matrix of CSA Metrics', fontweight='bold')
    
    plt.tight_layout()
    plt.savefig(savepath, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {savepath}")


def plot_site_means_forest(df, savepath):
    """Forest plot of site-level mean SC CSA with 95% CI."""
    fig, ax = plt.subplots(figsize=(10, 10))
    
    site_stats = []
    for site, group in df.groupby('site'):
        vals = group['SC_CSA_t2star'].dropna()
        if len(vals) >= 2:
            mean = vals.mean()
            sem = vals.std() / np.sqrt(len(vals))
            ci_low = mean - 1.96 * sem
            ci_high = mean + 1.96 * sem
            vendor = VENDOR_MAP.get(site, 'Unknown')
            site_stats.append({
                'site': site, 'mean': mean, 'ci_low': ci_low, 'ci_high': ci_high,
                'n': len(vals), 'vendor': vendor
            })
    
    site_stats = pd.DataFrame(site_stats).sort_values('mean')
    
    colors = {'Siemens': '#1f77b4', 'GE': '#d62728', 'Philips': '#2ca02c', 'Unknown': '#7f7f7f'}
    
    for i, row in site_stats.iterrows():
        y = list(site_stats.index).index(i)
        color = colors.get(row['vendor'], '#7f7f7f')
        ax.errorbar(row['mean'], y, xerr=[[row['mean'] - row['ci_low']], [row['ci_high'] - row['mean']]],
                    fmt='o', color=color, capsize=4, markersize=8, linewidth=1.5)
        ax.text(row['ci_high'] + 0.5, y, f"n={row['n']}", va='center', fontsize=8)
    
    ax.set_yticks(range(len(site_stats)))
    ax.set_yticklabels(site_stats['site'], fontsize=9)
    
    grand_mean = df['SC_CSA_t2star'].mean()
    ax.axvline(x=grand_mean, color='red', linestyle='--', linewidth=1.5, label=f'Grand mean: {grand_mean:.2f} mm²')
    ax.legend()
    
    ax.set_xlabel('T2* SC CSA (mm²) [Mean ± 95% CI]')
    ax.set_title('Forest Plot: Site-Level Mean T2* SC CSA', fontweight='bold')
    
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='#1f77b4', label='Siemens'),
        Patch(facecolor='#d62728', label='GE'),
        Patch(facecolor='#2ca02c', label='Philips'),
    ]
    ax.legend(handles=legend_elements + [plt.Line2D([0], [0], color='red', linestyle='--', label=f'Grand mean: {grand_mean:.2f}')],
              loc='lower right')
    
    plt.tight_layout()
    plt.savefig(savepath, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {savepath}")


# ============================================================
# 4. Main Execution
# ============================================================
def main():
    print("=" * 70)
    print("Spine-Generic QM Harmonization Analysis")
    print("=" * 70)
    
    # --- Collect Data ---
    print("\n[1] Collecting T2* CSA data...")
    df_t2star = collect_t2star_data()
    print(f"  Collected {len(df_t2star)} subjects from {df_t2star['site'].nunique()} sites")
    print(f"  Vendors: {df_t2star['vendor'].value_counts().to_dict()}")
    
    print("\n[2] Collecting T2w CSA data...")
    df_t2w = collect_t2w_data()
    print(f"  Collected {len(df_t2w)} subjects from {df_t2w['site'].nunique()} sites")
    
    print("\n[3] Collecting T2w CSA by vertebral level...")
    df_vl = collect_t2w_by_vertebral_level()
    print(f"  Collected {len(df_vl)} records across {df_vl['vert_level'].nunique() if not df_vl.empty else 0} vertebral levels")
    
    # --- Merge T2w and T2* data ---
    print("\n[4] Merging T2w and T2* data...")
    merged = pd.merge(df_t2star, df_t2w[['subject', 'SC_CSA_t2w', 'n_slices_t2w']], 
                      on='subject', how='outer')
    print(f"  Merged dataset: {len(merged)} subjects")
    overlap = merged.dropna(subset=['SC_CSA_t2star', 'SC_CSA_t2w'])
    print(f"  Subjects with both T2* and T2w: {len(overlap)}")
    
    # --- Save merged data ---
    merged.to_csv(os.path.join(DATA_DIR, 'merged_csa_data.csv'), index=False)
    df_t2star.to_csv(os.path.join(DATA_DIR, 't2star_csa_data.csv'), index=False)
    df_t2w.to_csv(os.path.join(DATA_DIR, 't2w_csa_data.csv'), index=False)
    if not df_vl.empty:
        df_vl.to_csv(os.path.join(DATA_DIR, 't2w_by_vertebral_level.csv'), index=False)
    
    # --- Descriptive Statistics ---
    print("\n[5] Computing descriptive statistics...")
    
    # Overall stats for T2*
    t2star_cols = ['SC_CSA_t2star', 'GM_CSA_t2star', 'WM_CSA_t2star', 
                   'GM_WM_ratio_t2star', 'GM_SC_frac_t2star']
    desc_t2star = descriptive_stats(df_t2star, t2star_cols)
    desc_t2star.to_csv(os.path.join(DATA_DIR, 'descriptive_stats_t2star.csv'), index=False)
    print("\n  T2* Descriptive Statistics:")
    print(desc_t2star.to_string(index=False))
    
    # T2w stats
    t2w_cols = ['SC_CSA_t2w']
    desc_t2w = descriptive_stats(df_t2w, t2w_cols)
    desc_t2w.to_csv(os.path.join(DATA_DIR, 'descriptive_stats_t2w.csv'), index=False)
    print("\n  T2w Descriptive Statistics:")
    print(desc_t2w.to_string(index=False))
    
    # By vendor
    desc_vendor = descriptive_stats(df_t2star, t2star_cols, group_col='vendor')
    desc_vendor.to_csv(os.path.join(DATA_DIR, 'descriptive_stats_by_vendor.csv'))
    print("\n  T2* Stats by Vendor:")
    print(desc_vendor.to_string())
    
    # By site
    desc_site = descriptive_stats(df_t2star, ['SC_CSA_t2star', 'GM_CSA_t2star'], group_col='site')
    desc_site.to_csv(os.path.join(DATA_DIR, 'descriptive_stats_by_site.csv'))
    print("\n  T2* Stats by Site (first 10):")
    print(desc_site.head(10).to_string())
    
    # --- CV Analysis ---
    print("\n[6] Coefficient of Variation analysis...")
    cv_sc = coefficient_of_variation(df_t2star, 'SC_CSA_t2star', 'site')
    cv_gm = coefficient_of_variation(df_t2star, 'GM_CSA_t2star', 'site')
    cv_sc.to_csv(os.path.join(DATA_DIR, 'cv_sc_by_site.csv'), index=False)
    cv_gm.to_csv(os.path.join(DATA_DIR, 'cv_gm_by_site.csv'), index=False)
    
    # Overall CV (across all subjects)
    overall_cv_sc = (df_t2star['SC_CSA_t2star'].std() / df_t2star['SC_CSA_t2star'].mean()) * 100
    overall_cv_gm = (df_t2star['GM_CSA_t2star'].std() / df_t2star['GM_CSA_t2star'].mean()) * 100
    overall_cv_wm = (df_t2star['WM_CSA_t2star'].std() / df_t2star['WM_CSA_t2star'].mean()) * 100
    print(f"  Overall CV - SC: {overall_cv_sc:.2f}%, GM: {overall_cv_gm:.2f}%, WM: {overall_cv_wm:.2f}%")
    print(f"  Mean within-site CV - SC: {cv_sc['CV(%)'].mean():.2f}%, GM: {cv_gm['CV(%)'].mean():.2f}%")
    
    # Pooled within-site CV
    pooled_cv_sc = np.sqrt(np.mean(cv_sc['CV(%)']**2))
    pooled_cv_gm = np.sqrt(np.mean(cv_gm['CV(%)']**2))
    print(f"  Pooled within-site CV - SC: {pooled_cv_sc:.2f}%, GM: {pooled_cv_gm:.2f}%")
    
    # Between-site CV
    site_means = df_t2star.groupby('site')['SC_CSA_t2star'].mean()
    between_cv_sc = (site_means.std() / site_means.mean()) * 100
    site_means_gm = df_t2star.groupby('site')['GM_CSA_t2star'].mean()
    between_cv_gm = (site_means_gm.std() / site_means_gm.mean()) * 100
    print(f"  Between-site CV - SC: {between_cv_sc:.2f}%, GM: {between_cv_gm:.2f}%")
    
    # --- Statistical Tests ---
    print("\n[7] Statistical tests...")
    
    # Normality test
    for col in ['SC_CSA_t2star', 'GM_CSA_t2star']:
        vals = df_t2star[col].dropna()
        if len(vals) >= 3:
            stat, p = shapiro(vals)
            print(f"  Shapiro-Wilk {col}: W={stat:.4f}, p={p:.4f} {'(normal)' if p > 0.05 else '(non-normal)'}")
    
    # Kruskal-Wallis across sites
    kw_sc = perform_kruskal_wallis(df_t2star, 'SC_CSA_t2star', 'site')
    kw_gm = perform_kruskal_wallis(df_t2star, 'GM_CSA_t2star', 'site')
    print(f"  Kruskal-Wallis SC CSA by site: H={kw_sc['statistic']}, p={kw_sc['p_value']:.6e}")
    print(f"  Kruskal-Wallis GM CSA by site: H={kw_gm['statistic']}, p={kw_gm['p_value']:.6e}")
    
    # Kruskal-Wallis across vendors
    kw_vendor_sc = perform_kruskal_wallis(df_t2star, 'SC_CSA_t2star', 'vendor')
    kw_vendor_gm = perform_kruskal_wallis(df_t2star, 'GM_CSA_t2star', 'vendor')
    if kw_vendor_sc:
        print(f"  Kruskal-Wallis SC CSA by vendor: H={kw_vendor_sc['statistic']}, p={kw_vendor_sc['p_value']:.6e}")
    if kw_vendor_gm:
        print(f"  Kruskal-Wallis GM CSA by vendor: H={kw_vendor_gm['statistic']}, p={kw_vendor_gm['p_value']:.6e}")
    
    # Pairwise vendor comparison (Mann-Whitney U)
    vendors_present = df_t2star['vendor'].unique()
    vendors_present = [v for v in vendors_present if v != 'Unknown']
    print(f"\n  Pairwise vendor comparisons (Mann-Whitney U):")
    pairwise_results = []
    for i in range(len(vendors_present)):
        for j in range(i+1, len(vendors_present)):
            v1, v2 = vendors_present[i], vendors_present[j]
            d1 = df_t2star[df_t2star['vendor'] == v1]['SC_CSA_t2star'].dropna()
            d2 = df_t2star[df_t2star['vendor'] == v2]['SC_CSA_t2star'].dropna()
            if len(d1) >= 3 and len(d2) >= 3:
                stat, p = mannwhitneyu(d1, d2, alternative='two-sided')
                pairwise_results.append({'comparison': f'{v1} vs {v2}', 'metric': 'SC_CSA', 
                                         'U': stat, 'p_value': p, 'n1': len(d1), 'n2': len(d2)})
                print(f"    {v1} vs {v2} (SC CSA): U={stat:.1f}, p={p:.6e}")
            
            d1 = df_t2star[df_t2star['vendor'] == v1]['GM_CSA_t2star'].dropna()
            d2 = df_t2star[df_t2star['vendor'] == v2]['GM_CSA_t2star'].dropna()
            if len(d1) >= 3 and len(d2) >= 3:
                stat, p = mannwhitneyu(d1, d2, alternative='two-sided')
                pairwise_results.append({'comparison': f'{v1} vs {v2}', 'metric': 'GM_CSA',
                                         'U': stat, 'p_value': p, 'n1': len(d1), 'n2': len(d2)})
                print(f"    {v1} vs {v2} (GM CSA): U={stat:.1f}, p={p:.6e}")
    
    pd.DataFrame(pairwise_results).to_csv(os.path.join(DATA_DIR, 'pairwise_vendor_comparisons.csv'), index=False)
    
    # Levene's test for variance homogeneity
    vendors_data = [df_t2star[df_t2star['vendor'] == v]['SC_CSA_t2star'].dropna().values 
                    for v in vendors_present if len(df_t2star[df_t2star['vendor'] == v]) >= 3]
    if len(vendors_data) >= 2:
        stat, p = levene(*vendors_data)
        print(f"\n  Levene's test (SC CSA variance across vendors): W={stat:.4f}, p={p:.4f}")
    
    # T2w vs T2* comparison
    overlap = merged.dropna(subset=['SC_CSA_t2star', 'SC_CSA_t2w'])
    if len(overlap) >= 3:
        r, p = stats.pearsonr(overlap['SC_CSA_t2w'], overlap['SC_CSA_t2star'])
        bias = (overlap['SC_CSA_t2w'] - overlap['SC_CSA_t2star']).mean()
        bias_sd = (overlap['SC_CSA_t2w'] - overlap['SC_CSA_t2star']).std()
        loa_low = bias - 1.96 * bias_sd
        loa_high = bias + 1.96 * bias_sd
        print(f"\n  T2w vs T2* SC CSA comparison (n={len(overlap)}):")
        print(f"    Pearson r = {r:.4f}, p = {p:.6e}")
        print(f"    Mean bias (T2w - T2*) = {bias:.2f} mm²")
        print(f"    95% LoA: [{loa_low:.2f}, {loa_high:.2f}] mm²")
        
        # Paired test
        stat_w, p_w = stats.wilcoxon(overlap['SC_CSA_t2w'], overlap['SC_CSA_t2star'])
        print(f"    Wilcoxon signed-rank: W={stat_w:.1f}, p={p_w:.6e}")
    
    # --- Linear Mixed Model / Regression ---
    print("\n[8] Linear regression: CSA ~ vendor + field_strength...")
    try:
        reg_df = df_t2star.dropna(subset=['SC_CSA_t2star', 'vendor']).copy()
        reg_df = reg_df[reg_df['vendor'] != 'Unknown']
        model = ols('SC_CSA_t2star ~ C(vendor)', data=reg_df).fit()
        print(f"  R² = {model.rsquared:.4f}")
        print(f"  F-statistic = {model.fvalue:.4f}, p = {model.f_pvalue:.6e}")
        with open(os.path.join(DATA_DIR, 'regression_results.txt'), 'w') as f:
            f.write(str(model.summary()))
    except Exception as e:
        print(f"  Regression error: {e}")
    
    # --- Voxel resolution analysis ---
    print("\n[9] Voxel resolution analysis...")
    voxel_stats = df_t2star.groupby('vendor')[['voxel_dx', 'voxel_dy', 'voxel_dz']].agg(['mean', 'std']).round(3)
    print(voxel_stats)
    voxel_stats.to_csv(os.path.join(DATA_DIR, 'voxel_resolution_by_vendor.csv'))
    
    # --- Generate Figures ---
    print("\n[10] Generating figures...")
    plot_site_variability_t2star(df_t2star, os.path.join(FIGURES_DIR, 'fig1_site_variability_t2star.png'))
    plot_vendor_comparison(df_t2star, os.path.join(FIGURES_DIR, 'fig2_vendor_comparison.png'))
    plot_t2w_vs_t2star(merged, os.path.join(FIGURES_DIR, 'fig3_t2w_vs_t2star.png'))
    plot_cv_by_site(df_t2star, os.path.join(FIGURES_DIR, 'fig4_cv_by_site.png'))
    if not df_vl.empty:
        plot_vertebral_level_distribution(df_vl, os.path.join(FIGURES_DIR, 'fig5_vertebral_level.png'))
    plot_correlation_matrix(merged, os.path.join(FIGURES_DIR, 'fig6_correlation_matrix.png'))
    plot_site_means_forest(df_t2star, os.path.join(FIGURES_DIR, 'fig7_forest_plot.png'))
    
    # --- Summary Report ---
    print("\n" + "=" * 70)
    print("ANALYSIS COMPLETE")
    print("=" * 70)
    print(f"\nData saved to: {DATA_DIR}")
    print(f"Figures saved to: {FIGURES_DIR}")
    
    # Return key stats for paper writing
    return {
        'n_subjects_t2star': len(df_t2star),
        'n_sites_t2star': df_t2star['site'].nunique(),
        'n_subjects_t2w': len(df_t2w),
        'n_sites_t2w': df_t2w['site'].nunique(),
        'n_overlap': len(overlap),
        'vendor_counts': df_t2star['vendor'].value_counts().to_dict(),
        'overall_cv_sc': overall_cv_sc,
        'overall_cv_gm': overall_cv_gm,
        'overall_cv_wm': overall_cv_wm,
        'mean_within_cv_sc': cv_sc['CV(%)'].mean(),
        'mean_within_cv_gm': cv_gm['CV(%)'].mean(),
        'pooled_within_cv_sc': pooled_cv_sc,
        'pooled_within_cv_gm': pooled_cv_gm,
        'between_cv_sc': between_cv_sc,
        'between_cv_gm': between_cv_gm,
        'kw_sc': kw_sc,
        'kw_gm': kw_gm,
        'kw_vendor_sc': kw_vendor_sc,
        'kw_vendor_gm': kw_vendor_gm,
        'desc_t2star': desc_t2star,
        'desc_t2w': desc_t2w,
    }


if __name__ == '__main__':
    results = main()
