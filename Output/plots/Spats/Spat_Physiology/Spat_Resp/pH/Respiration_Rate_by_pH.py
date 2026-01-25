import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
import statsmodels.api as sm
from statsmodels.formula.api import ols
import matplotlib.ticker as ticker
import matplotlib.patches as mpatches
from statsmodels.stats.multicomp import pairwise_tukeyhsd

# ==========================================
# 1. SETUP & DATA LOADING
# ==========================================
filename = "pH_Assay.csv"
df = pd.read_csv(filename)

# Normalize / rename the response column
df = df.rename(columns={"OxyRate nmol/mm2/min": "OxyRate"})

# Keep only HF / NF
target_morphs = ["HF", "NF"]
df = df[df["morph"].isin(target_morphs)].copy()

# --- MODIFIED: Specific Outlier Removal ---
# 1. Remove pH 8.2 NF Technical Errors (Plate D failures > 2.0)
# (Note: pH is likely float 8.2 here)
df = df[~((df['pH'] == 8.2) & (df['morph'] == 'NF') & (df['OxyRate'] > 2.0))]

# 2. Remove pH 7.8 HF Statistical Outliers (> 1.5)
# This removes both 2.36 AND 1.89 (which the previous >2.09 filter missed)
df = df[~((df['pH'] == 7.8) & (df['morph'] == 'HF') & (df['OxyRate'] > 1.5))]
# ------------------------------------------

# Ensure pH is clean and ordered
df["pH"] = df["pH"].astype(str).str.strip()
ph_order = ["8.2", "7.8", "7.6"]
df["pH"] = pd.Categorical(df["pH"], categories=ph_order, ordered=True)

# Optional: drop rows with unexpected pH values
df = df[df["pH"].isin(ph_order)].copy()

# ==========================================
# 2. STATISTICS (Two-Way ANOVA & Tukey's HSD)
# ==========================================

# 2.1 Calculate descriptive statistics (N, Mean, Std Dev) grouped by 'morph' and 'pH'
descriptive_stats_ph = df.groupby(['morph', 'pH'], observed=False)['OxyRate'].agg(
    N='count',
    Mean='mean',
    Std_Dev_Col='std'
).reset_index()
descriptive_stats_ph = descriptive_stats_ph.rename(columns={'Std_Dev_Col': 'Std Dev'})

# 2.2 Format the descriptive statistics table to 3 decimal places and print it
descriptive_stats_ph_formatted = descriptive_stats_ph.round(3)
print("\n=== Descriptive Statistics (pH Assay) ===")
print(descriptive_stats_ph_formatted.to_string(index=False))

model = ols("OxyRate ~ C(morph) + C(pH) + C(morph):C(pH)", data=df).fit()
anova_table = sm.stats.anova_lm(model, typ=2)

# 2.3 Extract F-value, df, and P-value for main effects and interaction
anova_results_data_ph = {
    'Source': ['C(morph)', 'C(pH)', 'C(morph):C(pH)']
}

for stat in ['df', 'F', 'PR(>F)']:
    col_name = 'F-value' if stat == 'F' else ('P-value' if stat == 'PR(>F)' else stat)
    anova_results_data_ph[col_name] = [
        anova_table.loc['C(morph)', stat] if 'C(morph)' in anova_table.index else np.nan,
        anova_table.loc['C(pH)', stat] if 'C(pH)' in anova_table.index else np.nan,
        anova_table.loc['C(morph):C(pH)', stat] if 'C(morph):C(pH)' in anova_table.index else np.nan
    ]

anova_results_df_ph = pd.DataFrame(anova_results_data_ph)

# 2.4 Format the ANOVA results table to 3 decimal places and print it
anova_results_df_ph_formatted = anova_results_df_ph.round(3)
print("\n=== Two-Way ANOVA Results (pH Assay) ===")
print(anova_results_df_ph_formatted.to_string(index=False))

# 2.5 Tukey's HSD Post-hoc Test
df['group'] = df['morph'].astype(str) + "_" + df['pH'].astype(str)
tukey_result = pairwise_tukeyhsd(endog=df['OxyRate'], groups=df['group'], alpha=0.05)

print("\n=== Tukey's HSD Post-hoc Test Results ===")
print(tukey_result)

# ==========================================
# 3. STYLING SETUP
# ==========================================
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

fig, ax = plt.subplots(figsize=(7, 6))

# ==========================================
# 4. PLOT CONSTRUCTION (MANUAL ALIGNMENT)
# ==========================================
box_width = 0.6
offset = box_width / 4

x_positions = {
    "HF": -offset,
    "NF": +offset
}

for morph in morph_order:
    sub = df[df["morph"] == morph]

    # A. Boxplots
    sns.boxplot(
        data=sub, x="pH", y="OxyRate",
        order=ph_order,
        width=box_width / 2,
        color=variant_colors[morph],
        showfliers=False, linewidth=1,
        boxprops=dict(edgecolor="black", alpha=0.6),
        whiskerprops=dict(color="black"),
        capprops=dict(color="black"),
        medianprops=dict(color="black"),
        ax=ax,
        positions=np.arange(len(ph_order)) + x_positions[morph]
    )

    # B. Jittered points
    for i, ph in enumerate(ph_order):
        y = sub[sub["pH"] == ph]["OxyRate"]
        positions_for_jitter = np.arange(len(ph_order)) + x_positions[morph]
        x = np.random.normal(positions_for_jitter[i], 0.04, size=len(y))

        ax.scatter(
            x, y,
            s=40, c=variant_colors[morph],
            edgecolors="black", linewidths=1,
            alpha=0.7, zorder=3
        )

    # C. Mean diamonds
    means = sub.groupby("pH", observed=False)["OxyRate"].mean().reindex(ph_order)

    ax.plot(
        np.arange(len(ph_order)) + x_positions[morph],
        means,
        marker="D", linestyle="--",
        markersize=8,
        markerfacecolor=variant_colors[morph],
        markeredgecolor="black",
        markeredgewidth=1.5,
        linewidth=2,
        color=variant_colors[morph],
        zorder=4
    )

# ==========================================
# 6. LABELS & TITLE
# ==========================================
ax.set_ylabel("Respiration rate (nmol/mm$^2$/min)", fontsize=13)
ax.set_xlabel("pH Level", fontsize=13)
ax.set_title("Respiration Rate of settled spats by pH and morph", fontsize=14, pad=15)

# Custom Legend
legend_handles = []
for morph in morph_order:
    patch = mpatches.Patch(color=variant_colors[morph], label=morph, edgecolor='black', linewidth=1, alpha=0.6)
    legend_handles.append(patch)
ax.legend(handles=legend_handles, labels=morph_order, title=None, frameon=False, loc="upper right")

# Y-axis formatting
max_val = df["OxyRate"].max()
ax.set_ylim(0, max_val * 1.2)
ax.yaxis.set_major_locator(ticker.MultipleLocator(1))

sns.despine(ax=ax)
fig.tight_layout()

plt.savefig("Respiration_Rate_by_pH.png", dpi=600)
#plt.savefig("Respiration_Rate_by_pH.pdf", dpi=600)
#plt.savefig("Respiration_Rate_by_pH.tiff", dpi=600)
#plt.savefig("Respiration_Rate_by_pH_transparent.svg", dpi=600, transparent=True)
#plt.savefig("Respiration_Rate_by_pH_white.svg", dpi=600, facecolor='white')

plt.show()