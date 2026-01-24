import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

# ==========================================
# 1. SETUP & DATA LOADING
# ==========================================
filename = 'Ambient_size.csv'
# Using read_csv since the provided file is a CSV
df = pd.read_csv(filename)

# 2. FILTERING
# (Optional) List samples to exclude if necessary, currently commented out
# samples_to_exclude = ['sample_id_1', 'sample_id_2']
# df = df[~df['SampleCode'].isin(samples_to_exclude)].copy()

# 3. CLEANING
# Clean whitespace from string columns
df['morph'] = df['morph'].str.strip()
df['SampleCode'] = df['SampleCode'].str.strip()

# Drop rows where Area or Morph is missing
df = df.dropna(subset=['area mm^2', 'morph'])

# ==========================================
# 4. STATISTICS
# ==========================================
# Extract groups
hf_area = df[df['morph'] == 'HF']['area mm^2']
nf_area = df[df['morph'] == 'NF']['area mm^2']
pf_area = df[df['morph'] == 'PF']['area mm^2']

# Calculate P-value (comparing HF vs NF to match previous guideline logic)
# If you want an ANOVA across all 3, use: stats.f_oneway(hf_area, nf_area, pf_area)
t_stat, p_val = stats.ttest_ind(hf_area, nf_area)

# 5. STYLING
plt.rcParams.update(
    {
        'font.family': 'serif',
        'axes.edgecolor': 'black',
        'axes.linewidth': 1,
        'xtick.color': 'black',
        'ytick.color': 'black',
        'text.color': 'black',
        'axes.labelcolor': 'black',
        'legend.frameon': False,
        'ytick.labelsize': 11,
        'xtick.labelsize': 11,
    }
)

# Define colors (Added specific color for PF)
# HF=Green, NF=Red, PF=Blue/Orange (using Orange for contrast)
variant_colors = {'HF': '#2ca02c', 'NF': '#d62728', 'PF': '#ff7f0e'}
morph_order = ['HF', 'PF', 'NF'] # Logical order: High, Partial, None

plt.figure(figsize=(6, 6))

# Box Plot
ax = sns.boxplot(
    data=df,
    x='morph',
    y='area mm^2',
    order=morph_order,
    hue='morph',
    palette=variant_colors,
    width=0.8,
    showfliers=False,
    linewidth=1,
    legend=False,
    boxprops=dict(edgecolor='black', alpha=0.5),
    whiskerprops=dict(color='black'),
    capprops=dict(color='black'),
    medianprops=dict(color='black'),
)

# Jitter Points
sns.stripplot(
    data=df,
    x='morph',
    y='area mm^2',
    order=morph_order,
    hue='morph',
    palette=variant_colors,
    jitter=0.15,
    size=5,
    linewidth=1,
    edgecolor='black',
    alpha=0.6,
    legend=False,
    zorder=3,
)

# Mean Diamond
means = df.groupby('morph')['area mm^2'].mean().reindex(morph_order)
mean_colors = [variant_colors[m] for m in morph_order]
# scatter x-range must match the number of morphs (0, 1, 2)
plt.scatter(
    range(len(morph_order)),
    means,
    c=mean_colors,
    marker='D',
    s=60,
    edgecolors='black',
    linewidths=2,
    alpha=0.7,
    zorder=4,
)

# 6. ANNOTATIONS (P-value only)
# UPDATE THESE manually based on your actual post-hoc stats results
# letters_map = {'HF': 'a', 'PF': 'ab', 'NF': 'b'}
# y_offset_percent = 0.05

max_val_global = df['area mm^2'].max()

# for i, m in enumerate(morph_order):
#     if m in df['morph'].unique(): # Only plot if morph exists in data
#         subset = df[df['morph'] == m]['area mm^2']
#         highest_point = subset.max()

#         # Add offset
#         pos_y = highest_point + (max_val_global * y_offset_percent)

#         letter = letters_map.get(m, "")
#         plt.text(x=i, y=pos_y, s=letter, ha='center', va='bottom', size=12, weight='bold')

# P-value (Bottom Right)
# Note: This P-value reflects the HF vs NF comparison calculated above
# p_text = f"p(HF vs NF) = {p_val:.4f}" if p_val >= 0.0001 else "p(HF vs NF) < 0.0001"
# plt.text(
#     0.95,
#     0.02,
#     s=p_text,
#     transform=ax.transAxes,
#     ha='right',
#     va='bottom',
#     size=10,
#     style='italic',
#     color='black',
# )

# 7. LABELS & TITLE
plt.ylabel("Spat Area ($mm^2$)", size=13)
plt.xlabel("") # Hide X label as categories are self-explanatory
plt.title("Size distribution of settled spats from difrrent morphs", size=14, pad=15)

# Adjust Y limit to fit annotations
plt.ylim(0, max_val_global * 1.15)
# Optional: Set specific ticks if needed
plt.yticks(np.arange(0, max_val_global, 2))

sns.despine()
plt.tight_layout()

# 8. SAVING
plt.savefig("Spat_Size_by_morph.pdf", dpi=600)
plt.savefig("Spat_Size_by_morph.tiff", dpi=600)
plt.savefig("Spat_Size_by_morph_transparent.svg", transparent=True)
plt.savefig("Spat_Size_by_morph_white_bg.svg", facecolor='white')
plt.savefig("Spat_Size_by_morph.png", dpi=600)

plt.show()