import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

# 1. DATA (8 points per group - Matches JMP Output)
data = {
    'Colony': [
        'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', # HF Group
        'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', # PF Group
        'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H'  # NF Group
    ],
    'Morph': [
        'HF', 'HF', 'HF', 'HF', 'HF', 'HF', 'HF', 'HF',
        'PF', 'PF', 'PF', 'PF', 'PF', 'PF', 'PF', 'PF',
        'NF', 'NF', 'NF', 'NF', 'NF', 'NF', 'NF', 'NF'
    ],
    'Total Release': [
        199, 23, 621, 57, 5, 16, 9, 5,  # HF Values
        164, 22, 950, 39, 4, 10, 2, 4,  # PF Values
        30, 13, 505, 15, 6, 8, 4, 2     # NF Values
    ]
}
df = pd.DataFrame(data)
df['ColonyTotal'] = df.groupby('Colony')['Total Release'].transform('sum')
df['Percent'] = (df['Total Release'] / df['ColonyTotal']) * 100

# 2. CALC WHISKERS (To ensure annotations don't overlap data)
def get_upper_whisker(series):
    Q1 = series.quantile(0.25)
    Q3 = series.quantile(0.75)
    IQR = Q3 - Q1
    upper_bound = Q3 + 1.5 * IQR
    whisker = series[series <= upper_bound].max()
    return whisker if not pd.isna(whisker) else Q3

# 3. STYLING
plt.rcParams.update({
    'font.family': 'serif',
    'axes.edgecolor': 'black', 'axes.linewidth': 1,
    'xtick.color': 'black', 'ytick.color': 'black',
    'text.color': 'black', 'axes.labelcolor': 'black',
    'legend.frameon': False,
    'ytick.labelsize': 11,
    'xtick.labelsize': 11
})

variant_colors = {'HF': '#2ca02c', 'PF': '#ff7f0e', 'NF': '#d62728'}
morph_order = ['HF', 'PF', 'NF']

plt.figure(figsize=(5, 6))

# 4. PLOT LAYERS
ax = sns.boxplot(
    data=df, x='Morph', y='Percent', order=morph_order, hue='Morph', palette=variant_colors,
    width=0.7, showfliers=False, linewidth=1, legend=False,
    boxprops=dict(edgecolor='black', alpha=0.5),
    whiskerprops=dict(color='black'), capprops=dict(color='black'), medianprops=dict(color='black')
)

sns.stripplot(
    data=df, x='Morph', y='Percent', order=morph_order, hue='Morph', palette=variant_colors,
    jitter=0.15, size=6, linewidth=1, edgecolor='black', alpha=0.7, legend=False, zorder=3
)

means = df.groupby('Morph')['Percent'].mean().reindex(morph_order)
mean_colors = [variant_colors[m] for m in morph_order]
plt.scatter(range(3), means, c=mean_colors, marker='D', s=40, edgecolors='black', linewidths=2, alpha=0.7, zorder=4)

# 5. ANNOTATIONS (Corrected based on JMP Text)
# HF is A -> 'a'
# PF is A and B -> 'ab'
# NF is B -> 'b'
letters_map = {'HF': 'a', 'PF': 'ab', 'NF': 'b'}
y_offset = 5

for i, m in enumerate(morph_order):
    subset = df[df['Morph'] == m]['Percent']
    # Calculate highest point (dot or whisker) for clearance
    highest_point = subset.max()
    pos_y = highest_point + y_offset

    letter = letters_map[m]
    plt.text(x=i, y=pos_y, s=letter, ha='center', va='bottom', size=12, weight='bold')

# P-value (Exact from JMP ANOVA table)
plt.text(x=1.0, y=-0.12, s="p = 0.0028", transform=ax.transAxes,
         ha='right', va='top', size=10, style='italic')

# 6. LABELS
plt.ylabel("Release (%)", size=13)
plt.xlabel("")
plt.title("Planulae released by morph\n(Ambient Conditions)", size=14, pad=15)
plt.ylim(0, 85)
plt.yticks(np.arange(0, 81, 20))

sns.despine()
plt.tight_layout()
plt.savefig("Planuae_Release_by_morph.pdf", dpi=600)
plt.savefig("Planuae_Release_by_morph.tiff", dpi=600)
plt.savefig("Planuae_Release_by_morph.png", dpi=600)
plt.show()