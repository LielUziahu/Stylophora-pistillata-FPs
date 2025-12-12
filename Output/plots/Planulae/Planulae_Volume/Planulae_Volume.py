import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

# ==========================================
# 1. SETUP & DATA LOADING
# ==========================================
# Make sure to upload this file to Colab first!
filename = 'Planulae size.xlsx'
df = pd.read_excel('Planulae size.xlsx')

# 2. FILTERING
samples_to_exclude = [135, 148]
df = df[~df['no'].isin(samples_to_exclude)].copy()

# 3. CLEANING
df['morph'] = df['morph'].str.strip()
df = df.dropna(subset=['Volume mm³', 'morph'])

# Calculate P-value for text
hf_vol = df[df['morph'] == 'HF']['Volume mm³']
nf_vol = df[df['morph'] == 'NF']['Volume mm³']
t_stat, p_val = stats.ttest_ind(hf_vol, nf_vol)

# 4. STYLING
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

# Box Plot
ax = sns.boxplot(
    data=df, x='morph', y='Volume mm³', order=morph_order, hue='morph', palette=variant_colors,
    width=0.8, showfliers=False, linewidth=1, legend=False,
    boxprops=dict(edgecolor='black', alpha=0.5),
    whiskerprops=dict(color='black'), capprops=dict(color='black'), medianprops=dict(color='black')
)

# Jitter Points
sns.stripplot(
    data=df, x='morph', y='Volume mm³', order=morph_order, hue='morph', palette=variant_colors,
    jitter=0.15, size=5, linewidth=1, edgecolor='black',
    alpha=0.6, legend=False, zorder=3
)

# Mean Diamond
means = df.groupby('morph')['Volume mm³'].mean().reindex(morph_order)
mean_colors = [variant_colors[m] for m in morph_order]
plt.scatter(range(2), means, c=mean_colors, marker='D', s=60, edgecolors='black', linewidths=2, alpha=0.7, zorder=4)

# 5. ANNOTATIONS (Letters & P-value)
letters_map = {'HF': 'a', 'NF': 'b'}
y_offset_percent = 0.05

max_val_global = df['Volume mm³'].max()

for i, m in enumerate(morph_order):
    subset = df[df['morph'] == m]['Volume mm³']
    highest_point = subset.max()

    # Add offset
    pos_y = highest_point + (max_val_global * y_offset_percent)

    letter = letters_map[m]
    plt.text(x=i, y=pos_y, s=letter, ha='center', va='bottom', size=12, weight='bold')

# P-value (Bottom Right, under axis)
# Positioned at (1.0, -0.15) relative to axes coordinates (bottom-right corner, below axis)
p_text = f"p = {p_val:.4f}" if p_val >= 0.0001 else "p < 0.0001"
plt.text(1.0, -0.1, s=p_text, transform=ax.transAxes, # Moved p-value closer to x-axis
         ha='right', va='top', size=10, style='italic', color='black')

# 6. LABELS & TITLE
plt.ylabel("Volume ($mm^3$)", size=13)
plt.xlabel("")
plt.title("Volume distribution of different planulae morphs", size=14, pad=15)

# Adjust Y limit to fit annotations
plt.ylim(0.1, max_val_global * 1.15)
plt.yticks(np.arange(0.1, 0.8, 0.3)) # Adjusted y-axis ticks to start from 0.2

sns.despine()
plt.tight_layout()
plt.savefig("Volume_by_morph.pdf", dpi=600)
plt.savefig("Volume_by_morph.tiff", dpi=600)
plt.savefig("Volume_by_morph.png", dpi=600)
plt.show()