import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
import matplotlib.ticker as ticker
from statsmodels.stats.multicomp import pairwise_tukeyhsd # Import for Tukey's HSD

# ==========================================
# 1. SETUP & DATA LOADING
# ==========================================
filename = 'Ambient_Assay.csv'
df_resp = pd.read_csv(filename)

# Filter for the relevant morphs (just in case)
target_morphs = ['HF', 'PF', 'NF']
df_resp = df_resp[df_resp['morph'].isin(target_morphs)]

# --- NEW: Filter samples above 2.09 and print removed samples ---
removed_ambient_samples = df_resp[df_resp['OxyRate nmol/mm2/min'] > 2.09]
if not removed_ambient_samples.empty:
    print("--- Ambient Samples Removed (OxyRate > 2.09) ---")
    for index, row in removed_ambient_samples.iterrows():
        print(f"Morph: {row['morph']}, SampleCode: {row['SampleCode']}, OxyRate: {row['OxyRate nmol/mm2/min']:.2f}")
    print("---------------------------------------------------")
df_resp = df_resp[df_resp['OxyRate nmol/mm2/min'] <= 2.09]
# ----------------------------------------------------------------

# ==========================================
# 2. STATISTICS (One-Way ANOVA)
# ==========================================
# We use ANOVA because we now have 3 groups (HF, PF, NF)
hf_vals = df_resp[df_resp['morph'] == 'HF']['OxyRate nmol/mm2/min']
pf_vals = df_resp[df_resp['morph'] == 'PF']['OxyRate nmol/mm2/min']
nf_vals = df_resp[df_resp['morph'] == 'NF']['OxyRate nmol/mm2/min']

f_stat, p_val_resp = stats.f_oneway(hf_vals, pf_vals, nf_vals)

# Present ANOVA results in a table
anova_table_data = {'Statistic': ['F-statistic', 'P-value'],
                    'Value': [f_stat, p_val_resp]}
anova_df = pd.DataFrame(anova_table_data)
print("\n=== One-Way ANOVA Results ===")
display(anova_df)

# --- NEW: Tukey's HSD Post-hoc Test ---
# Create a 'group' column for Tukey's HSD
df_resp['group'] = df_resp['morph']

tukey_result = pairwise_tukeyhsd(endog=df_resp['OxyRate nmol/mm2/min'], groups=df_resp['group'], alpha=0.05)

print("\n=== Tukey's HSD Post-hoc Test Results ===")
print(tukey_result)
# -------------------------------------------

# Letter assignment logic
# If p > 0.05, no difference (all 'a'). If p < 0.05, we would need post-hoc tests.
# For now, if non-significant, all get 'a'.
if p_val_resp < 0.05:
    # Based on Tukey's HSD: HF is different from NF and PF, but NF and PF are not different.
    letters_map = {'HF': 'a', 'PF': 'b', 'NF': 'b'}
else:
    letters_map = {'HF': 'a', 'PF': 'a', 'NF': 'a'}

# ==========================================
# 3. STYLING SETUP
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

# Define Colors and Order
variant_colors = {'HF': '#2ca02c', 'PF': '#ff7f0e', 'NF': '#d62728'}
morph_order = ['HF', 'PF', 'NF']

plt.figure(figsize=(6, 6))

# ==========================================
# 4. PLOT CONSTRUCTION
# ==========================================

# A. Box Plot
ax = sns.boxplot(
    data=df_resp, x='morph', y='OxyRate nmol/mm2/min', hue='morph', order=morph_order, palette=variant_colors, legend=False,
    width=0.7, showfliers=False, linewidth=1,
    boxprops=dict(edgecolor='black', alpha=0.5),
    whiskerprops=dict(color='black'), capprops=dict(color='black'), medianprops=dict(color='black')
)

# B. Jitter Points
sns.stripplot(
    data=df_resp, x='morph', y='OxyRate nmol/mm2/min', order=morph_order, hue='morph', palette=variant_colors,
    jitter=0.15, size=6, linewidth=1, edgecolor='black',
    alpha=0.7, legend=False, zorder=3
)

# C. Mean Diamond
means_resp = df_resp.groupby('morph')['OxyRate nmol/mm2/min'].mean().reindex(morph_order)
mean_colors_resp = [variant_colors[m] for m in morph_order]
plt.scatter(range(len(morph_order)), means_resp, c=mean_colors_resp, marker='D', s=70, edgecolors='black', linewidths=2, alpha=0.9, zorder=4)

# ==========================================
# 5. ANNOTATIONS
# ==========================================
y_offset_percent = 0.08
max_val_global_resp = df_resp['OxyRate nmol/mm2/min'].max()

for i, m in enumerate(morph_order):
    subset_resp = df_resp[df_resp['morph'] == m]['OxyRate nmol/mm2/min']
    highest_point_resp = subset_resp.max()

    # Add offset
    pos_y_resp = highest_point_resp + (max_val_global_resp * y_offset_percent)

    letter = letters_map[m]
    plt.text(x=i, y=pos_y_resp, s=letter, ha='center', va='bottom', size=14, weight='bold')

# P-value (Bottom Right)
p_text_resp = f"p = {p_val_resp:.3f}" if p_val_resp >= 0.0001 else "p < 0.0001"
plt.text(0.95, -0.1, s=p_text_resp, transform=ax.transAxes,
         ha='right', va='bottom', size=11, style='italic', color='black')

# ==========================================
# 6. LABELS & TITLE
# ==========================================
# Updated Y-label to match your CSV unit (nmol/mm2/min)
plt.ylabel("Respiration rate (nmol/mm$^2$/min)", size=13)
plt.xlabel("")
plt.title("Respiration rate of settled spats by morph \nunder ambient conditions", size=14, pad=15)

# Adjust Y limit and set ticks
plt.ylim(0, df_resp['OxyRate nmol/mm2/min'].max() * 1.25)
ax.yaxis.set_major_locator(ticker.MultipleLocator(1))

sns.despine()
plt.tight_layout()

# Save        ---> remove # to get a downloadedble plod
plt.savefig("Respiration_Rate_Spats_Ambient.png", dpi=600)
#plt.savefig("Respiration_Rate_Spats_Ambient.tiff", dpi=600)
#plt.savefig("Respiration_Rate_Spats_Ambient_white.svg", dpi=600, facecolor='white')
#plt.savefig("Respiration_Rate_Spats_Ambient_transparent.svg", dpi=600, transparent=True)
#plt.savefig("Respiration_Rate_Spats_Ambient.pdf", dpi=600)

plt.show()