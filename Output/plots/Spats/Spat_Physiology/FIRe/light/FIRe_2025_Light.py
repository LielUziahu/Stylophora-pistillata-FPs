import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patches as mpatches
import matplotlib.ticker as ticker
from statsmodels.formula.api import ols
import statsmodels.api as sm
from statsmodels.stats.multicomp import pairwise_tukeyhsd # Import for Tukey's HSD

# ==========================================
# 1. GLOBAL SETUP & DATA LOADING
# ==========================================
filename = "FIRe_2025_Light.csv"
df = pd.read_csv(filename)

# Clean columns
df.columns = df.columns.str.strip()
if 'Site' in df.columns:
    df = df.rename(columns={'Site': 'morph'})
df['morph'] = df['morph'].str.strip()

# Clean Light column and set order
df['Light'] = df['Light'].astype(str).str.strip()
light_order = ["10m", "50m", "Dark"]
df = df[df['Light'].isin(light_order)] # Keep only valid light levels

# --- DATA CLEANING (Global Filters) ---
# 1. Connectivity (p) < 0.8
if 'p' in df.columns:
    df = df[df['p'] <= 0.8].copy()

# 2. Pmax < 10000 
# Note: Ambient used < 450. Here, High Light samples are naturally higher. 
# We filter only obvious errors (> 10,000).
if 'Pmax.e.s' in df.columns:
    df = df[df['Pmax.e.s'] <= 450].copy()

# 3. Sigma < 1200
# Note: Ambient used < 850. Adjusted slightly for Light dataset range (max was ~1300).
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

def generate_stats_table_2way(dataframe, col):
    """Calculates Two-Way ANOVA stats."""
    try:
        model = ols(f"Q('{col}') ~ C(morph) + C(Light) + C(morph):C(Light)", data=dataframe).fit()
        anova = sm.stats.anova_lm(model, typ=2)
        
        # Extract P-values
        p_morph = anova.loc['C(morph)', 'PR(>F)']
        p_light = anova.loc['C(Light)', 'PR(>F)']
        p_inter = anova.loc['C(morph):C(Light)', 'PR(>F)']
        
        return pd.DataFrame({
            'Factor': ['Morph', 'Light', 'Interaction'],
            'P-value': [p_morph, p_light, p_inter],
            'Significance': [
                "ns" if p >= 0.05 else ("*" if p >= 0.01 else ("**" if p >= 0.001 else "***")) 
                for p in [p_morph, p_light, p_inter]
            ]
        })
    except Exception as e:
        return pd.DataFrame({'Error': [str(e)]})

def generate_tukey_hsd_table(dataframe, col_name):
    """Calculates Tukey's HSD post-hoc test."""
    try:
        # Create a combined group for pairwise comparisons
        df_tukey = dataframe.copy()
        df_tukey['group'] = df_tukey['morph'].astype(str) + '_' + df_tukey['Light'].astype(str)
        
        # Ensure there's enough data for each group
        if df_tukey['group'].nunique() > 1:
            tukey_results = pairwise_tukeyhsd(endog=df_tukey[col_name], groups=df_tukey['group'], alpha=0.05)
            tukey_df = pd.DataFrame(data=tukey_results._results_table.data[1:], columns=tukey_results._results_table.data[0])
            tukey_df['reject'] = tukey_df['reject'].apply(lambda x: 'Yes' if x else 'No') # Convert boolean to 'Yes'/'No'
            return tukey_df
        else:
            return pd.DataFrame({'Note': ['Not enough groups for Tukey HSD.']})
    except Exception as e:
        return pd.DataFrame({'Error': [str(e)]})

# ==========================================
# PLOT 1: Fm (Maximum Fluorescence)
# ==========================================
col = "Fm"
ylabel = "Fm (a.u.)"
title = "Maximum Fluorescence (Fm)\nunder different light conditions"
tick_spacing = 10.0  # <--- CHANGE TICK SPACING HERE

if col in df.columns:
    sub = df.dropna(subset=[col, 'morph', 'Light'])
    sub = sub[sub['morph'].isin(morph_order)]

    print(f"\n\nTwo-Way ANOVA Statistical Analysis for {col}:")
    stats_df = generate_stats_table_2way(sub, col)
    print(stats_df)
    
    print(f"\n\nTukey's HSD Post-Hoc Test for {col}:")
    tukey_hsd_df = generate_tukey_hsd_table(sub, col)
    print(tukey_hsd_df)

    fig, ax = plt.subplots(figsize=(7, 6))

    # Calculate positions
    box_width = 0.6
    x_pos = np.arange(len(light_order))
    offset = box_width / 4

    hf_data = [sub[(sub['Light']==l) & (sub['morph']=='HF')][col].values for l in light_order]
    bp1 = ax.boxplot(hf_data, positions=x_pos - offset, widths=box_width/2, patch_artist=True, showfliers=False)
    for patch in bp1['boxes']:
        patch.set_facecolor(variant_colors['HF'])
        patch.set_alpha(0.6)
        patch.set_edgecolor('black')
    for element in ['whiskers', 'caps', 'medians']:
        plt.setp(bp1[element], color='black')

    nf_data = [sub[(sub['Light']==l) & (sub['morph']=='NF')][col].values for l in light_order]
    bp2 = ax.boxplot(nf_data, positions=x_pos + offset, widths=box_width/2, patch_artist=True, showfliers=False)
    for patch in bp2['boxes']:
        patch.set_facecolor(variant_colors['NF'])
        patch.set_alpha(0.6)
        patch.set_edgecolor('black')
    for element in ['whiskers', 'caps', 'medians']:
        plt.setp(bp2[element], color='black')

    for i, light in enumerate(light_order):
        y_hf = sub[(sub['Light']==light) & (sub['morph']=='HF')][col]
        x_hf = np.random.normal(i - offset, 0.04, size=len(y_hf))
        ax.scatter(x_hf, y_hf, s=20, color=variant_colors['HF'], edgecolors='black', linewidth=0.5, alpha=0.7, zorder=3)
        
        y_nf = sub[(sub['Light']==light) & (sub['morph']=='NF')][col]
        x_nf = np.random.normal(i + offset, 0.04, size=len(y_nf))
        ax.scatter(x_nf, y_nf, s=20, color=variant_colors['NF'], edgecolors='black', linewidth=0.5, alpha=0.7, zorder=3)

    # Mean Diamonds
    means_hf = sub[sub['morph']=='HF'].groupby('Light')[col].mean().reindex(light_order)
    ax.scatter(x_pos - offset, means_hf, marker='D', s=50, color=variant_colors['HF'], edgecolors='black', zorder=4)
    means_nf = sub[sub['morph']=='NF'].groupby('Light')[col].mean().reindex(light_order)
    ax.scatter(x_pos + offset, means_nf, marker='D', s=50, color=variant_colors['NF'], edgecolors='black', zorder=4)

    # Labels
    ax.set_ylabel(ylabel, fontsize=13)
    ax.set_xlabel("Light/Depth Condition", fontsize=13)
    ax.set_title(title, fontsize=14, pad=15)
    ax.set_xticks(x_pos)
    ax.set_xticklabels(light_order)
    
    # Annotate P-value (Morph)
    p_morph = stats_df.loc[stats_df['Factor']=='Morph', 'P-value'].values[0] if 'P-value' in stats_df.columns else 1.0
    if p_morph < 0.05: # Only annotate if significant
        p_text = f"p(morph) = {p_morph:.3f}" if p_morph >= 0.001 else "p(morph) < 0.001"
        ax.text(0.95, 0.05, p_text, transform=ax.transAxes, ha='right', va='bottom', fontsize=10, style='italic')

    # Ticks & Limits
    ax.yaxis.set_major_locator(ticker.MultipleLocator(tick_spacing))
    min_val, max_val = sub[col].min(), sub[col].max()
    padding = (max_val - min_val) * 0.15
    ax.set_ylim(max(0, min_val - padding), max_val + padding)

    # Legend
    handles = [mpatches.Patch(color=variant_colors[m], label=m, alpha=0.6, ec='black') for m in morph_order]
    ax.legend(handles=handles, title=None, frameon=False, loc="upper right")

    sns.despine(ax=ax)
    plt.tight_layout()
    plt.savefig(f"FIRe_Light_{col}.tiff", dpi=600)
    plt.savefig(f"FIRe_Light_{col}.png", dpi=600)
    plt.savefig(f"FIRe_Light_{col}.pdf")
    plt.savefig(f"FIRe_Light_{col}_transparent.svg", transparent=True)
    plt.savefig(f"FIRe_Light_{col}_white.svg", facecolor='white', transparent=False)
    plt.show()

# ==========================================
# PLOT 2: Fv/Fm (Efficiency)
# ==========================================
col = "fv_fm"
ylabel = "Fv/Fm"
title = "Photosynthetic Efficiency (Fv/Fm)\nunder different light conditions"
tick_spacing = 0.3  # <--- CHANGE TICK SPACING HERE

if col in df.columns:
    sub = df.dropna(subset=[col, 'morph', 'Light'])
    sub = sub[sub['morph'].isin(morph_order)]

    print(f"\n\nTwo-Way ANOVA Statistical Analysis for {col}:")
    stats_df = generate_stats_table_2way(sub, col)
    print(stats_df)

    print(f"\n\nTukey's HSD Post-Hoc Test for {col}:")
    tukey_hsd_df = generate_tukey_hsd_table(sub, col)
    print(tukey_hsd_df)
    
    fig, ax = plt.subplots(figsize=(7, 6))
    box_width = 0.6
    x_pos = np.arange(len(light_order))
    offset = box_width / 4

    hf_data = [sub[(sub['Light']==l) & (sub['morph']=='HF')][col].values for l in light_order]
    bp1 = ax.boxplot(hf_data, positions=x_pos - offset, widths=box_width/2, patch_artist=True, showfliers=False)
    for patch in bp1['boxes']:
        patch.set_facecolor(variant_colors['HF'])
        patch.set_alpha(0.6)
        patch.set_edgecolor('black')
    for element in ['whiskers', 'caps', 'medians']:
        plt.setp(bp1[element], color='black')

    nf_data = [sub[(sub['Light']==l) & (sub['morph']=='NF')][col].values for l in light_order]
    bp2 = ax.boxplot(nf_data, positions=x_pos + offset, widths=box_width/2, patch_artist=True, showfliers=False)
    for patch in bp2['boxes']:
        patch.set_facecolor(variant_colors['NF'])
        patch.set_alpha(0.6)
        patch.set_edgecolor('black')
    for element in ['whiskers', 'caps', 'medians']:
        plt.setp(bp2[element], color='black')

    for i, light in enumerate(light_order):
        y_hf = sub[(sub['Light']==light) & (sub['morph']=='HF')][col]
        x_hf = np.random.normal(i - offset, 0.04, size=len(y_hf))
        ax.scatter(x_hf, y_hf, s=20, color=variant_colors['HF'], edgecolors='black', linewidth=0.5, alpha=0.7, zorder=3)
        y_nf = sub[(sub['Light']==light) & (sub['morph']=='NF')][col]
        x_nf = np.random.normal(i + offset, 0.04, size=len(y_nf))
        ax.scatter(x_nf, y_nf, s=20, color=variant_colors['NF'], edgecolors='black', linewidth=0.5, alpha=0.7, zorder=3)

    means_hf = sub[sub['morph']=='HF'].groupby('Light')[col].mean().reindex(light_order)
    ax.scatter(x_pos - offset, means_hf, marker='D', s=50, color=variant_colors['HF'], edgecolors='black', zorder=4)
    means_nf = sub[sub['morph']=='NF'].groupby('Light')[col].mean().reindex(light_order)
    ax.scatter(x_pos + offset, means_nf, marker='D', s=50, color=variant_colors['NF'], edgecolors='black', zorder=4)

    ax.set_ylabel(ylabel, fontsize=13)
    ax.set_xlabel("Light/Depth Condition", fontsize=13)
    ax.set_title(title, fontsize=14, pad=15)
    ax.set_xticks(x_pos)
    ax.set_xticklabels(light_order)
    
    p_morph = stats_df.loc[stats_df['Factor']=='Morph', 'P-value'].values[0] if 'P-value' in stats_df.columns else 1.0
    if p_morph < 0.05: # Only annotate if significant
        p_text = f"p(morph) = {p_morph:.3f}" if p_morph >= 0.001 else "p(morph) < 0.001"
        ax.text(0.95, 0.05, p_text, transform=ax.transAxes, ha='right', va='bottom', fontsize=10, style='italic')

    ax.yaxis.set_major_locator(ticker.MultipleLocator(tick_spacing))
    min_val, max_val = sub[col].min(), sub[col].max()
    padding = (max_val - min_val) * 0.15
    ax.set_ylim(max(0, min_val - padding), max_val + padding)

    handles = [mpatches.Patch(color=variant_colors[m], label=m, alpha=0.6, ec='black') for m in morph_order]
    ax.legend(handles=handles, title=None, frameon=False, loc="upper right")

    sns.despine(ax=ax)
    plt.tight_layout()
    plt.savefig(f"FIRe_Light_{col}.tiff", dpi=600)
    plt.savefig(f"FIRe_Light_{col}.png", dpi=600)
    plt.savefig(f"FIRe_Light_{col}.pdf")
    plt.savefig(f"FIRe_Light_{col}_transparent.svg", transparent=True)
    plt.savefig(f"FIRe_Light_{col}_white.svg", facecolor='white', transparent=False)
    plt.show()

# ==========================================
# PLOT 3: Pmax (Max Photosynthesis)
# ==========================================
col = "Pmax.e.s"
ylabel = "Pmax ($e^- s^{-1}$)"
title = "Max Photosynthesis (Pmax) \nunder different light conditions"
tick_spacing = 100.0  # <--- CHANGE TICK SPACING HERE

if col in df.columns:
    sub = df.dropna(subset=[col, 'morph', 'Light'])
    sub = sub[sub['morph'].isin(morph_order)]
    
    print(f"\n\nTwo-Way ANOVA Statistical Analysis for {col}:")
    stats_df = generate_stats_table_2way(sub, col)
    print(stats_df)

    print(f"\n\nTukey's HSD Post-Hoc Test for {col}:")
    tukey_hsd_df = generate_tukey_hsd_table(sub, col)
    print(tukey_hsd_df)

    fig, ax = plt.subplots(figsize=(7, 6))
    box_width = 0.6
    x_pos = np.arange(len(light_order))
    offset = box_width / 4

    hf_data = [sub[(sub['Light']==l) & (sub['morph']=='HF')][col].values for l in light_order]
    bp1 = ax.boxplot(hf_data, positions=x_pos - offset, widths=box_width/2, patch_artist=True, showfliers=False)
    for patch in bp1['boxes']:
        patch.set_facecolor(variant_colors['HF'])
        patch.set_alpha(0.6)
        patch.set_edgecolor('black')
    for element in ['whiskers', 'caps', 'medians']:
        plt.setp(bp1[element], color='black')

    nf_data = [sub[(sub['Light']==l) & (sub['morph']=='NF')][col].values for l in light_order]
    bp2 = ax.boxplot(nf_data, positions=x_pos + offset, widths=box_width/2, patch_artist=True, showfliers=False)
    for patch in bp2['boxes']:
        patch.set_facecolor(variant_colors['NF'])
        patch.set_alpha(0.6)
        patch.set_edgecolor('black')
    for element in ['whiskers', 'caps', 'medians']:
        plt.setp(bp2[element], color='black')

    for i, light in enumerate(light_order):
        y_hf = sub[(sub['Light']==light) & (sub['morph']=='HF')][col]
        x_hf = np.random.normal(i - offset, 0.04, size=len(y_hf))
        ax.scatter(x_hf, y_hf, s=20, color=variant_colors['HF'], edgecolors='black', linewidth=0.5, alpha=0.7, zorder=3)
        y_nf = sub[(sub['Light']==light) & (sub['morph']=='NF')][col]
        x_nf = np.random.normal(i + offset, 0.04, size=len(y_nf))
        ax.scatter(x_nf, y_nf, s=20, color=variant_colors['NF'], edgecolors='black', linewidth=0.5, alpha=0.7, zorder=3)

    means_hf = sub[sub['morph']=='HF'].groupby('Light')[col].mean().reindex(light_order)
    ax.scatter(x_pos - offset, means_hf, marker='D', s=50, color=variant_colors['HF'], edgecolors='black', zorder=4)
    means_nf = sub[sub['morph']=='NF'].groupby('Light')[col].mean().reindex(light_order)
    ax.scatter(x_pos + offset, means_nf, marker='D', s=50, color=variant_colors['NF'], edgecolors='black', zorder=4)

    ax.set_ylabel(ylabel, fontsize=13)
    ax.set_xlabel("Light/Depth Condition", fontsize=13)
    ax.set_title(title, fontsize=14, pad=15)
    ax.set_xticks(x_pos)
    ax.set_xticklabels(light_order)
    
    p_morph = stats_df.loc[stats_df['Factor']=='Morph', 'P-value'].values[0] if 'P-value' in stats_df.columns else 1.0
    if p_morph < 0.05: # Only annotate if significant
        p_text = f"p(morph) = {p_morph:.3f}" if p_morph >= 0.001 else "p(morph) < 0.001"
        ax.text(0.95, 0.05, p_text, transform=ax.transAxes, ha='right', va='bottom', fontsize=10, style='italic')

    ax.yaxis.set_major_locator(ticker.MultipleLocator(tick_spacing))
    min_val, max_val = sub[col].min(), sub[col].max()
    padding = (max_val - min_val) * 0.15
    ax.set_ylim(max(0, min_val - padding), max_val + padding)

    handles = [mpatches.Patch(color=variant_colors[m], label=m, alpha=0.6, ec='black') for m in morph_order]
    ax.legend(handles=handles, title=None, frameon=False, loc="upper right")

    sns.despine(ax=ax)
    plt.tight_layout()
    plt.savefig(f"FIRe_Light_{col}.tiff", dpi=600)
    plt.savefig(f"FIRe_Light_{col}.png", dpi=600)
    plt.savefig(f"FIRe_Light_{col}.pdf")
    plt.savefig(f"FIRe_Light_{col}_transparent.svg", transparent=True)
    plt.savefig(f"FIRe_Light_{col}_white.svg", facecolor='white', transparent=False)
    plt.show()

# ==========================================
# PLOT 4: Connectivity (p)
# ==========================================
col = "p"
ylabel = "p (connectivity)"
title = "Connectivity (p)\nunder different light conditions"
tick_spacing = 0.3  # <--- CHANGE TICK SPACING HERE

if col in df.columns:
    sub = df.dropna(subset=[col, 'morph', 'Light'])
    sub = sub[sub['morph'].isin(morph_order)]
    
    print(f"\n\nTwo-Way ANOVA Statistical Analysis for {col}:")
    stats_df = generate_stats_table_2way(sub, col)
    print(stats_df)

    print(f"\n\nTukey's HSD Post-Hoc Test for {col}:")
    tukey_hsd_df = generate_tukey_hsd_table(sub, col)
    print(tukey_hsd_df)

    fig, ax = plt.subplots(figsize=(7, 6))
    box_width = 0.6
    x_pos = np.arange(len(light_order))
    offset = box_width / 4

    hf_data = [sub[(sub['Light']==l) & (sub['morph']=='HF')][col].values for l in light_order]
    bp1 = ax.boxplot(hf_data, positions=x_pos - offset, widths=box_width/2, patch_artist=True, showfliers=False)
    for patch in bp1['boxes']:
        patch.set_facecolor(variant_colors['HF'])
        patch.set_alpha(0.6)
        patch.set_edgecolor('black')
    for element in ['whiskers', 'caps', 'medians']:
        plt.setp(bp1[element], color='black')

    nf_data = [sub[(sub['Light']==l) & (sub['morph']=='NF')][col].values for l in light_order]
    bp2 = ax.boxplot(nf_data, positions=x_pos + offset, widths=box_width/2, patch_artist=True, showfliers=False)
    for patch in bp2['boxes']:
        patch.set_facecolor(variant_colors['NF'])
        patch.set_alpha(0.6)
        patch.set_edgecolor('black')
    for element in ['whiskers', 'caps', 'medians']:
        plt.setp(bp2[element], color='black')

    for i, light in enumerate(light_order):
        y_hf = sub[(sub['Light']==light) & (sub['morph']=='HF')][col]
        x_hf = np.random.normal(i - offset, 0.04, size=len(y_hf))
        ax.scatter(x_hf, y_hf, s=20, color=variant_colors['HF'], edgecolors='black', linewidth=0.5, alpha=0.7, zorder=3)
        y_nf = sub[(sub['Light']==light) & (sub['morph']=='NF')][col]
        x_nf = np.random.normal(i + offset, 0.04, size=len(y_nf))
        ax.scatter(x_nf, y_nf, s=20, color=variant_colors['NF'], edgecolors='black', linewidth=0.5, alpha=0.7, zorder=3)

    means_hf = sub[sub['morph']=='HF'].groupby('Light')[col].mean().reindex(light_order)
    ax.scatter(x_pos - offset, means_hf, marker='D', s=50, color=variant_colors['HF'], edgecolors='black', zorder=4)
    means_nf = sub[sub['morph']=='NF'].groupby('Light')[col].mean().reindex(light_order)
    ax.scatter(x_pos + offset, means_nf, marker='D', s=50, color=variant_colors['NF'], edgecolors='black', zorder=4)

    ax.set_ylabel(ylabel, fontsize=13)
    ax.set_xlabel("Light/Depth Condition", fontsize=13)
    ax.set_title(title, fontsize=14, pad=15)
    ax.set_xticks(x_pos)
    ax.set_xticklabels(light_order)
    
    p_morph = stats_df.loc[stats_df['Factor']=='Morph', 'P-value'].values[0] if 'P-value' in stats_df.columns else 1.0
    if p_morph < 0.05: # Only annotate if significant
        p_text = f"p(morph) = {p_morph:.3f}" if p_morph >= 0.001 else "p(morph) < 0.001"
        ax.text(0.95, 0.05, p_text, transform=ax.transAxes, ha='right', va='bottom', fontsize=10, style='italic')

    ax.yaxis.set_major_locator(ticker.MultipleLocator(tick_spacing))
    min_val, max_val = sub[col].min(), sub[col].max()
    padding = (max_val - min_val) * 0.15
    ax.set_ylim(max(0, min_val - padding), max_val + padding)

    handles = [mpatches.Patch(color=variant_colors[m], label=m, alpha=0.6, ec='black') for m in morph_order]
    ax.legend(handles=handles, title=None, frameon=False, loc="upper right")

    sns.despine(ax=ax)
    plt.tight_layout()
    plt.savefig(f"FIRe_Light_{col}.tiff", dpi=600)
    plt.savefig(f"FIRe_Light_{col}.png", dpi=600)
    plt.savefig(f"FIRe_Light_{col}.pdf")
    plt.savefig(f"FIRe_Light_{col}_transparent.svg", transparent=True)
    plt.savefig(f"FIRe_Light_{col}_white.svg", facecolor='white', transparent=False)
    plt.show()

# ==========================================
# PLOT 5: Sigma (Cross-section)
# ==========================================
col = "Sigma"
ylabel ="Sigma ($\AA^2$)"
title = "Effective Cross-section ($\sigma_{PSII}$)\nunder different light conditions"
tick_spacing = 250.0  # <--- CHANGE TICK SPACING HERE

if col in df.columns:
    sub = df.dropna(subset=[col, 'morph', 'Light'])
    sub = sub[sub['morph'].isin(morph_order)]
    
    print(f"\n\nTwo-Way ANOVA Statistical Analysis for {col}:")
    stats_df = generate_stats_table_2way(sub, col)
    print(stats_df)

    print(f"\n\nTukey's HSD Post-Hoc Test for {col}:")
    tukey_hsd_df = generate_tukey_hsd_table(sub, col)
    print(tukey_hsd_df)

    fig, ax = plt.subplots(figsize=(7, 6))
    box_width = 0.6
    x_pos = np.arange(len(light_order))
    offset = box_width / 4

    hf_data = [sub[(sub['Light']==l) & (sub['morph']=='HF')][col].values for l in light_order]
    bp1 = ax.boxplot(hf_data, positions=x_pos - offset, widths=box_width/2, patch_artist=True, showfliers=False)
    for patch in bp1['boxes']:
        patch.set_facecolor(variant_colors['HF'])
        patch.set_alpha(0.6)
        patch.set_edgecolor('black')
    for element in ['whiskers', 'caps', 'medians']:
        plt.setp(bp1[element], color='black')

    nf_data = [sub[(sub['Light']==l) & (sub['morph']=='NF')][col].values for l in light_order]
    bp2 = ax.boxplot(nf_data, positions=x_pos + offset, widths=box_width/2, patch_artist=True, showfliers=False)
    for patch in bp2['boxes']:
        patch.set_facecolor(variant_colors['NF'])
        patch.set_alpha(0.6)
        patch.set_edgecolor('black')
    for element in ['whiskers', 'caps', 'medians']:
        plt.setp(bp2[element], color='black')

    for i, light in enumerate(light_order):
        y_hf = sub[(sub['Light']==light) & (sub['morph']=='HF')][col]
        x_hf = np.random.normal(i - offset, 0.04, size=len(y_hf))
        ax.scatter(x_hf, y_hf, s=20, color=variant_colors['HF'], edgecolors='black', linewidth=0.5, alpha=0.7, zorder=3)
        y_nf = sub[(sub['Light']==light) & (sub['morph']=='NF')][col]
        x_nf = np.random.normal(i + offset, 0.04, size=len(y_nf))
        ax.scatter(x_nf, y_nf, s=20, color=variant_colors['NF'], edgecolors='black', linewidth=0.5, alpha=0.7, zorder=3)

    means_hf = sub[sub['morph']=='HF'].groupby('Light')[col].mean().reindex(light_order)
    ax.scatter(x_pos - offset, means_hf, marker='D', s=50, color=variant_colors['HF'], edgecolors='black', zorder=4)
    means_nf = sub[sub['morph']=='NF'].groupby('Light')[col].mean().reindex(light_order)
    ax.scatter(x_pos + offset, means_nf, marker='D', s=50, color=variant_colors['NF'], edgecolors='black', zorder=4)

    ax.set_ylabel(ylabel, fontsize=13)
    ax.set_xlabel("Light/Depth Condition", fontsize=13)
    ax.set_title(title, fontsize=14, pad=15)
    ax.set_xticks(x_pos)
    ax.set_xticklabels(light_order)
    
    p_morph = stats_df.loc[stats_df['Factor']=='Morph', 'P-value'].values[0] if 'P-value' in stats_df.columns else 1.0
    if p_morph < 0.05: # Only annotate if significant
        p_text = f"p(morph) = {p_morph:.3f}" if p_morph >= 0.001 else "p(morph) < 0.001"
        ax.text(0.95, 0.05, p_text, transform=ax.transAxes, ha='right', va='bottom', fontsize=10, style='italic')

    ax.yaxis.set_major_locator(ticker.MultipleLocator(tick_spacing))
    min_val, max_val = sub[col].min(), sub[col].max()
    padding = (max_val - min_val) * 0.15
    ax.set_ylim(max(0, min_val - padding), max_val + padding)

    handles = [mpatches.Patch(color=variant_colors[m], label=m, alpha=0.6, ec='black') for m in morph_order]
    ax.legend(handles=handles, title=None, frameon=False, loc="upper right")

    sns.despine(ax=ax)
    plt.tight_layout()
    plt.savefig(f"FIRe_Light_{col}.tiff", dpi=600)
    plt.savefig(f"FIRe_Light_{col}.png", dpi=600)
    plt.savefig(f"FIRe_Light_{col}.pdf")
    plt.savefig(f"FIRe_Light_{col}_transparent.svg", transparent=True)
    plt.savefig(f"FIRe_Light_{col}_white.svg", facecolor='white', transparent=False)
    plt.show()