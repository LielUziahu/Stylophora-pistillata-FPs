import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

# ==========================================
# 1. SETUP & DATA LOADING
# ==========================================
# Make sure to upload this file to Colab first!
filename = 'GFP_raw_sum_normalized.csv'
df = pd.read_csv(filename)

# ==========================================
# 2. CLEANING & STATS
# ==========================================
df['Variant'] = df['Variant'].str.strip()
df = df.dropna(subset=['Intensity_per_µm2', 'Variant'])

# Run T-test
hf_vals = df[df['Variant'] == 'HF']['Intensity_per_µm2']
nf_vals = df[df['Variant'] == 'NF']['Intensity_per_µm2']
t_stat, p_val = stats.ttest_ind(hf_vals, nf_vals)

# ==========================================
# 3. STYLING
# ==========================================
plt.rcParams.update({
    'font.family': 'serif',
    'axes.edgecolor': 'black', 'axes.linewidth': 1,
    'xtick.color': 'black', 'ytick.color': 'black',
    'text.color': 'black', 'axes.labelcolor': 'black',
    'legend.frameon': False,
    'ytick.labelsize': 11,
    'xtick.labelsize': 11
})

variant_colors = {'HF': '#2ca02c', 'NF': '#d62728'}
morph_order = ['HF', 'NF']

plt.figure(figsize=(6, 6))

# ==========================================
# 4. PLOT CONSTRUCTION
# ==========================================

# A. Box Plot
ax = sns.boxplot(
    data=df, x='Variant', y='Intensity_per_µm2', order=morph_order, palette=variant_colors,
    width=0.6, showfliers=False, linewidth=1,
    boxprops=dict(edgecolor='black', alpha=0.5),
    whiskerprops=dict(color='black'), capprops=dict(color='black'), medianprops=dict(color='black')
)

# B. Jitter Points
sns.stripplot(
    data=df, x='Variant', y='Intensity_per_µm2', order=morph_order, hue='Variant', palette=variant_colors,
    jitter=0.15, size=5, linewidth=1, edgecolor='black',
    alpha=0.6, legend=False, zorder=3
)

# C. Mean Diamond
means = df.groupby('Variant')['Intensity_per_µm2'].mean().reindex(morph_order)
mean_colors = [variant_colors[m] for m in morph_order]
plt.scatter(range(2), means, c=mean_colors, marker='D', s=60, edgecolors='black', linewidths=2, alpha=0.7, zorder=4)

# ==========================================
# 5. ANNOTATIONS
# ==========================================
# Letters 'a' and 'b' (Since P < 0.0001, groups are different)
letters_map = {'HF': 'a', 'NF': 'b'}
y_offset_percent = 0.05
max_val_global = df['Intensity_per_µm2'].max()

for i, m in enumerate(morph_order):
    subset = df[df['Variant'] == m]['Intensity_per_µm2']
    highest_point = subset.max()
    pos_y = highest_point + (max_val_global * y_offset_percent)

    letter = letters_map[m]
    plt.text(x=i, y=pos_y, s=letter, ha='center', va='bottom', size=12, weight='bold')

# P-value (Bottom Right, under axis)
# Positioned closer to axis as requested
p_text = f"p = {p_val:.4f}" if p_val >= 0.0001 else "p < 0.0001"
plt.text(1.0, -0.12, s=p_text, transform=ax.transAxes,
         ha='right', va='top', size=10, style='italic', color='black')

# Labels & Title
plt.ylabel("Intensity per µm$^2$", size=13)
plt.xlabel("")
plt.title("GFP Intensity by Morph", size=14, pad=15)
plt.ylim(0, max_val_global * 1.15)
plt.yticks(np.arange(0, (max_val_global * 1.15) + 1, 100))

sns.despine()
plt.tight_layout()
plt.savefig("GFP Intensity by Morph.pdf", dpi=600)
plt.savefig("GFP Intensity by Morph.tiff", dpi=600)
plt.savefig("GFP Intensity by Morph.png", dpi=600)
plt.show()