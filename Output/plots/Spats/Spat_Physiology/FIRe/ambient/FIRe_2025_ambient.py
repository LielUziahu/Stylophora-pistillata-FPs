import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patches as mpatches
import matplotlib.ticker as ticker
from scipy import stats

# ==========================================
# 1. GLOBAL SETUP & DATA LOADING
# ==========================================
filename = "FIRe_2025_ambient.csv"
df = pd.read_csv(filename)

# Clean columns
df.columns = df.columns.str.strip()
if 'Site' in df.columns:
    df = df.rename(columns={'Site': 'morph'})
df['morph'] = df['morph'].str.strip()

# --- DATA CLEANING (Global Filters)---
# 1. Connectivity (p) < 0.8
if 'p' in df.columns:
    df = df[df['p'] <= 0.8].copy()

# 2. Pmax < 450
if 'Pmax.e.s' in df.columns:
    df = df[df['Pmax.e.s'] <= 450].copy()

# 3. Sigma < 850
if 'Sigma' in df.columns:
    df = df[df['Sigma'] <= 850].copy()

# 4. Fv/Fm > 0 (Sanity Check)
if 'fv_fm' in df.columns:
    df = df[df['fv_fm'] > 0].copy()

# Style Settings
plt.rcParams.update({
    "font.family": "serif",
    "axes.edgecolor": "black", "axes.linewidth": 1,
    "xtick.color": "black", "ytick.color": "black",
    "text.color": "black", "axes.labelcolor": "black",
    "legend.frameon": False,
    "ytick.labelsize": 11,
    "xtick.labelsize": 11
})

variant_colors = {"HF": "#2ca02c", "NF": "#d62728"}
morph_order = ["HF", "NF"]

# Helper to calculate stats
def get_pval(dataframe, col):
    hf = dataframe[dataframe['morph'] == 'HF'][col]
    nf = dataframe[dataframe['morph'] == 'NF'][col]
    if len(hf) > 1 and len(nf) > 1:
        _, p = stats.ttest_ind(hf, nf)
        return p
    return 1.0

def generate_stats_table(dataframe, col, morph_order):
    stats_data = {}
    for morph in morph_order:
        data = dataframe[dataframe['morph'] == morph][col].dropna()
        stats_data[f'{morph} Mean'] = data.mean()
        stats_data[f'{morph} Std'] = data.std()
        stats_data[f'{morph} Count'] = data.count()

    p_val = get_pval(dataframe, col)

    table_data = {
        'Variable': [col],
        'HF Mean': [f"{stats_data['HF Mean']:.3f}"],
        'HF Std': [f"{stats_data['HF Std']:.3f}"],
        'HF Count': [int(stats_data['HF Count'])],
        'NF Mean': [f"{stats_data['NF Mean']:.3f}"],
        'NF Std': [f"{stats_data['NF Std']:.3f}"],
        'NF Count': [int(stats_data['NF Count'])],
        'P-value': [f"{p_val:.4f}" if p_val < 0.05 else f"{p_val:.3f} (ns)"]
    }
    return pd.DataFrame(table_data)

# ==========================================
# PLOT 1: Fm (Maximum Fluorescence)
# ==========================================
col = "Fm"
ylabel = "Fm (a.u.)"
title = "Maximum Fluorescence (Fm) under ambient conditions"
tick_spacing = 10.0  # <--- CHANGE TICK SPACING HERE

if col in df.columns:
    sub = df.dropna(subset=[col, 'morph'])
    sub = sub[sub['morph'].isin(morph_order)]
    p_val = get_pval(sub, col)

    print()
    # Display statistical test name
    print("Statistical Test: Independent t-test")
    # Display stats table
    display(generate_stats_table(sub, col, morph_order))

    fig, ax = plt.subplots(figsize=(6, 6))
    
    # Boxplot & Jitter
    sns.boxplot(data=sub, x='morph', y=col, hue='morph', order=morph_order, width=0.5, showfliers=False,
                palette=variant_colors, ax=ax, boxprops=dict(edgecolor='black', alpha=0.6),
                whiskerprops=dict(color='black'), capprops=dict(color='black'), medianprops=dict(color='black'), legend=False)
    sns.stripplot(data=sub, x='morph', y=col, hue='morph', order=morph_order, palette=variant_colors, size=5, 
                  jitter=0.15, edgecolor='black', linewidth=1, alpha=0.7, zorder=3, ax=ax, legend=False)
    
    # Mean Diamond
    means = sub.groupby('morph')[col].mean().reindex(morph_order)
    ax.scatter(range(len(morph_order)), means, c=[variant_colors[m] for m in morph_order], 
               marker='D', s=60, edgecolors='black', linewidths=1.5, alpha=1.0, zorder=4)

    # Labels & Annotation
    ax.set_ylabel(ylabel, fontsize=13)
    ax.set_xlabel("")
    ax.set_title(title, fontsize=14, pad=15)
    if p_val < 0.05:
        p_text = f"p = {p_val:.4f}"
        ax.text(0.95, 0.02, p_text, transform=ax.transAxes, ha='right', va='bottom', fontsize=10, style='italic')
    
    # Ticks & Limits
    ax.yaxis.set_major_locator(ticker.MultipleLocator(tick_spacing))
    ax.set_ylim(0, sub[col].max() * 1.15)
    
    sns.despine(ax=ax)
    plt.tight_layout()
    #plt.savefig("FIRe_Fm.png", dpi=600)
    plt.show()

# ==========================================
# PLOT 2: Fv/Fm (Efficiency)
# ==========================================
col = "fv_fm"
ylabel = "Fv/Fm"
title = "Photosynthetic Efficiency (Fv/Fm)\nunder ambient conditions"
tick_spacing = 0.3  # <--- CHANGE TICK SPACING HERE

if col in df.columns:
    sub = df.dropna(subset=[col, 'morph'])
    sub = sub[sub['morph'].isin(morph_order)]
    p_val = get_pval(sub, col)

    print()
    # Display statistical test name
    print("Statistical Test: Independent t-test")
    # Display stats table
    display(generate_stats_table(sub, col, morph_order))

    fig, ax = plt.subplots(figsize=(6, 6))
    sns.boxplot(data=sub, x='morph', y=col, hue='morph', order=morph_order, width=0.5, showfliers=False,
                palette=variant_colors, ax=ax, boxprops=dict(edgecolor='black', alpha=0.6),
                whiskerprops=dict(color='black'), capprops=dict(color='black'), medianprops=dict(color='black'), legend=False)
    sns.stripplot(data=sub, x='morph', y=col, hue='morph', order=morph_order, palette=variant_colors, size=5, 
                  jitter=0.15, edgecolor='black', linewidth=1, alpha=0.7, zorder=3, ax=ax, legend=False)
    means = sub.groupby('morph')[col].mean().reindex(morph_order)
    ax.scatter(range(len(morph_order)), means, c=[variant_colors[m] for m in morph_order], 
               marker='D', s=60, edgecolors='black', linewidths=1.5, alpha=1.0, zorder=4)

    ax.set_ylabel(ylabel, fontsize=13)
    ax.set_xlabel("")
    ax.set_title(title, fontsize=14, pad=15)
    if p_val < 0.05:
        p_text = f"p = {p_val:.4f}"
        ax.text(0.95, 0.02, p_text, transform=ax.transAxes, ha='right', va='bottom', fontsize=10, style='italic')

    ax.yaxis.set_major_locator(ticker.MultipleLocator(tick_spacing))
    ax.set_ylim(0, sub[col].max() * 1.15)
    
    sns.despine(ax=ax)
    plt.tight_layout()
    #plt.savefig("FIRe_FvFm.png", dpi=600)
    plt.show()

# ==========================================
# PLOT 3: Pmax (Max Photosynthesis)
# ==========================================
col = "Pmax.e.s"
ylabel = "Pmax ($e^- s^{-1}$)"
title = "Max Photosynthesis (Pmax)\nunder ambient conditions"
tick_spacing = 100  # <--- CHANGE TICK SPACING HERE

if col in df.columns:
    sub = df.dropna(subset=[col, 'morph'])
    sub = sub[sub['morph'].isin(morph_order)]
    p_val = get_pval(sub, col)

    print()
    # Display statistical test name
    print("Statistical Test: Independent t-test")
    # Display stats table
    display(generate_stats_table(sub, col, morph_order))

    fig, ax = plt.subplots(figsize=(6, 6))
    sns.boxplot(data=sub, x='morph', y=col, hue='morph', order=morph_order, width=0.5, showfliers=False,
                palette=variant_colors, ax=ax, boxprops=dict(edgecolor='black', alpha=0.6),
                whiskerprops=dict(color='black'), capprops=dict(color='black'), medianprops=dict(color='black'), legend=False)
    sns.stripplot(data=sub, x='morph', y=col, hue='morph', order=morph_order, palette=variant_colors, size=5, 
                  jitter=0.15, edgecolor='black', linewidth=1, alpha=0.7, zorder=3, ax=ax, legend=False)
    means = sub.groupby('morph')[col].mean().reindex(morph_order)
    ax.scatter(range(len(morph_order)), means, c=[variant_colors[m] for m in morph_order], 
               marker='D', s=60, edgecolors='black', linewidths=1.5, alpha=1.0, zorder=4)

    ax.set_ylabel(ylabel, fontsize=13)
    ax.set_xlabel("")
    ax.set_title(title, fontsize=14, pad=15)
    if p_val < 0.05:
        p_text = f"p = {p_val:.4f}"
        ax.text(0.95, 0.02, p_text, transform=ax.transAxes, ha='right', va='bottom', fontsize=10, style='italic')

    ax.yaxis.set_major_locator(ticker.MultipleLocator(tick_spacing))
    ax.set_ylim(0, sub[col].max() * 1.15)
    
    sns.despine(ax=ax)
    plt.tight_layout()
    #plt.savefig("FIRe_Pmax.png", dpi=600)
    plt.show()

# ==========================================
# PLOT 4: Connectivity (p)
# ==========================================
col = "p"
ylabel = "p (connectivity)"
title = "Connectivity (p)\nunder ambient conditions"
tick_spacing = 0.3  # <--- CHANGE TICK SPACING HERE

if col in df.columns:
    sub = df.dropna(subset=[col, 'morph'])
    sub = sub[sub['morph'].isin(morph_order)]
    p_val = get_pval(sub, col)

    print()
    # Display statistical test name
    print("Statistical Test: Independent t-test")
    # Display stats table
    display(generate_stats_table(sub, col, morph_order))

    fig, ax = plt.subplots(figsize=(6, 6))
    sns.boxplot(data=sub, x='morph', y=col, hue='morph', order=morph_order, width=0.5, showfliers=False,
                palette=variant_colors, ax=ax, boxprops=dict(edgecolor='black', alpha=0.6),
                whiskerprops=dict(color='black'), capprops=dict(color='black'), medianprops=dict(color='black'), legend=False)
    sns.stripplot(data=sub, x='morph', y=col, hue='morph', order=morph_order, palette=variant_colors, size=5, 
                  jitter=0.15, edgecolor='black', linewidth=1, alpha=0.7, zorder=3, ax=ax, legend=False)
    means = sub.groupby('morph')[col].mean().reindex(morph_order)
    ax.scatter(range(len(morph_order)), means, c=[variant_colors[m] for m in morph_order], 
               marker='D', s=60, edgecolors='black', linewidths=1.5, alpha=1.0, zorder=4)

    ax.set_ylabel(ylabel, fontsize=13)
    ax.set_xlabel("")
    ax.set_title(title, fontsize=14, pad=15)
    if p_val < 0.05:
        p_text = f"p = {p_val:.4f}"
        ax.text(0.95, 0.02, p_text, transform=ax.transAxes, ha='right', va='bottom', fontsize=10, style='italic')

    ax.yaxis.set_major_locator(ticker.MultipleLocator(tick_spacing))
    ax.set_ylim(0, sub[col].max() * 1.15)
    
    sns.despine(ax=ax)
    plt.tight_layout()
    #plt.savefig("FIRe_p.png", dpi=600)
    plt.show()

# ==========================================
# PLOT 5: Sigma (Cross-section)
# ==========================================
col = "Sigma"
ylabel = "Sigma ($\\AA^2$)"
title = "Effective Cross-section ($\\sigma_{PSII}$) \nunder ambient conditions"
tick_spacing = 250.0  # <--- CHANGE TICK SPACING HERE

if col in df.columns:
    sub = df.dropna(subset=[col, 'morph'])
    sub = sub[sub['morph'].isin(morph_order)]
    p_val = get_pval(sub, col)

    print()
    # Display statistical test name
    print("Statistical Test: Independent t-test")
    # Display stats table
    display(generate_stats_table(sub, col, morph_order))

    fig, ax = plt.subplots(figsize=(6, 6))
    sns.boxplot(data=sub, x='morph', y=col, hue='morph', order=morph_order, width=0.5, showfliers=False,
                palette=variant_colors, ax=ax, boxprops=dict(edgecolor='black', alpha=0.6),
                whiskerprops=dict(color='black'), capprops=dict(color='black'), medianprops=dict(color='black'), legend=False)
    sns.stripplot(data=sub, x='morph', y=col, hue='morph', order=morph_order, palette=variant_colors, size=5, 
                  jitter=0.15, edgecolor='black', linewidth=1, alpha=0.7, zorder=3, ax=ax, legend=False)
    means = sub.groupby('morph')[col].mean().reindex(morph_order)
    ax.scatter(range(len(morph_order)), means, c=[variant_colors[m] for m in morph_order], 
               marker='D', s=60, edgecolors='black', linewidths=1.5, alpha=1.0, zorder=4)

    ax.set_ylabel(ylabel, fontsize=13)
    ax.set_xlabel("")
    ax.set_title(title, fontsize=14, pad=15)
    if p_val < 0.05:
        p_text = f"p = {p_val:.4f}"
        ax.text(0.95, 0.02, p_text, transform=ax.transAxes, ha='right', va='bottom', fontsize=10, style='italic')

    ax.yaxis.set_major_locator(ticker.MultipleLocator(tick_spacing))
    ax.set_ylim(0, sub[col].max() * 1.15)
    
    sns.despine(ax=ax)
    plt.tight_layout()
    #plt.savefig("FIRe_Sigma.png", dpi=600)
    plt.show()