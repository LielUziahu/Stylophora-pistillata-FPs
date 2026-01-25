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
filename = "Light_Assay.csv"
df = pd.read_csv(filename)

# --- Handle 'morph' and 'color' columns consistently ---
# The prompt states 'color' contains HF/NF info, so we use that as the canonical 'morph' column.
if 'color' in df.columns:
    # Save the desired 'morph' data from 'color' before dropping anything
    desired_morph_data = df['color'].copy()
    # Drop the original 'color' column
    df = df.drop(columns=['color'])
    # If an old 'morph' column exists, drop it to prevent duplicates
    if 'morph' in df.columns:
        df = df.drop(columns=['morph'])
    # Add the desired 'morph' data as the canonical 'morph' column
    df['morph'] = desired_morph_data
# else: If 'color' column doesn't exist, assume 'morph' column is already correctly named and contains HF/NF.

# Rename 'OxyRate nmol/mm2/min' to 'OxyRate'
df = df.rename(columns={
    "OxyRate nmol/mm2/min": "OxyRate"
})

# --- NEW: Filter samples above 2.09 and print removed samples ---
removed_light_samples = df[df['OxyRate'] > 2.09]
if not removed_light_samples.empty:
    print("--- Light Samples Removed (OxyRate > 2.09) ---")
    for index, row in removed_light_samples.iterrows():
        print(f"Morph: {row['morph']}, SampleCode: {row['SampleCode']}, OxyRate: {row['OxyRate']:.2f}")
    print("--------------------------------------------------")
df = df[df['OxyRate'] <= 2.09]
# ----------------------------------------------------------------

# Clean and Order the LightDepth column
# Logic: High Light (10m) -> Low Light (50m) -> No Light (Dark)
df["LightDepth"] = df["LightDepth"].astype(str).str.strip()
light_order = ["10m", "50m", "Dark"]
df["LightDepth"] = pd.Categorical(df["LightDepth"], categories=light_order, ordered=True)

# Drop rows with unexpected LightDepth values
df = df[df["LightDepth"].isin(light_order)].copy()
df = df.reset_index(drop=True) # Ensure unique index after filtering

# ==========================================
# 2. STATISTICS (Two-Way ANOVA & Tukey's HSD)
# ==========================================

# 2.1 Calculate descriptive statistics (N, Mean, Std Dev) grouped by 'morph' and 'LightDepth'
descriptive_stats_light = df.groupby(['morph', 'LightDepth'], observed=False)['OxyRate'].agg(
    N='count',
    Mean='mean',
    Std_Dev_Col='std'
).reset_index()
descriptive_stats_light = descriptive_stats_light.rename(columns={'Std_Dev_Col': 'Std Dev'})

# 2.2 Format the descriptive statistics table to 3 decimal places and print it
descriptive_stats_light_formatted = descriptive_stats_light.round(3)
print("\n=== Descriptive Statistics (Light Assay) ===")
print(descriptive_stats_light_formatted.to_string(index=False))

# Model: OxyRate depends on morph and LightDepth
model = ols("OxyRate ~ C(morph) + C(LightDepth) + C(morph):C(LightDepth)", data=df).fit()
anova_table = sm.stats.anova_lm(model, typ=2)

# 2.3 Extract F-value, df, and P-value for main effects and interaction
anova_results_data = {
    'Source': ['C(morph)', 'C(LightDepth)', 'C(morph):C(LightDepth)']
}

for stat_name in ['df', 'F', 'PR(>F)']:
    col_name = 'F-value' if stat_name == 'F' else ('P-value' if stat_name == 'PR(>F)' else stat_name)
    anova_results_data[col_name] = [
        anova_table.loc['C(morph)', stat_name] if 'C(morph)' in anova_table.index else np.nan,
        anova_table.loc['C(LightDepth)', stat_name] if 'C(LightDepth)' in anova_table.index else np.nan,
        anova_table.loc['C(morph):C(LightDepth)', stat_name] if 'C(morph):C(LightDepth)' in anova_table.index else np.nan
    ]

anova_results_df_light = pd.DataFrame(anova_results_data)

# 2.4 Format the ANOVA results table to 3 decimal places and print it
anova_results_df_light_formatted = anova_results_df_light.round(3)
print("\n=== Two-Way ANOVA Results (Light Assay) ===")
print(anova_results_df_light_formatted.to_string(index=False))

# 2.5 Tukey's HSD Post-hoc Test
df['group'] = df['morph'].astype(str) + "_" + df['LightDepth'].astype(str)
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
# 4. PLOT CONSTRUCTION
# ==========================================
box_width = 0.6
offset = box_width / 4

x_positions = {
    "HF": -offset,
    "NF": +offset
}

# Variable to plot on X-axis
x_col = "LightDepth"
x_order = light_order

for morph in morph_order:
    sub = df[df["morph"] == morph]

    # A. Boxplots
    sns.boxplot(
        data=sub, x=x_col, y="OxyRate",
        order=x_order,
        width=box_width / 2,
        color=variant_colors[morph],
        showfliers=False, linewidth=1,
        boxprops=dict(edgecolor="black", alpha=0.6),
        whiskerprops=dict(color="black"),
        capprops=dict(color="black"),
        medianprops=dict(color="black"),
        ax=ax,
        positions=np.arange(len(x_order)) + x_positions[morph]
    )

    # B. Jittered points
    for i, level in enumerate(x_order):
        y = sub[sub[x_col] == level]["OxyRate"]
        # Add random jitter to x
        x = np.random.normal(i + x_positions[morph], 0.04, size=len(y))

        ax.scatter(
            x, y,
            s=40, c=variant_colors[morph],
            edgecolors="black", linewidths=1,
            alpha=0.7, zorder=3
        )

    # C. Mean diamonds (changed from median)
    means = sub.groupby(x_col)["OxyRate"].mean().reindex(x_order)

    ax.plot(
        np.arange(len(x_order)) + x_positions[morph],
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
# 5. LABELS, LEGEND & ANNOTATION
# ==========================================
ax.set_ylabel("Respiration rate (nmol/mm$^2$/min)", fontsize=13)
ax.set_xlabel("Light/Depth Conditions", fontsize=13)
ax.set_title("Respiration Rate of settled spats by Light and morph", fontsize=14, pad=15)

# Custom Legend
legend_handles = []
for morph in morph_order:
    patch = mpatches.Patch(color=variant_colors[morph], label=morph, edgecolor='black', linewidth=1, alpha=0.6)
    legend_handles.append(patch)
ax.legend(handles=legend_handles, labels=morph_order, title=None, frameon=False, loc="upper right")

# Fix X-axis ticks
ax.set_xticks(np.arange(len(x_order)))
ax.set_xticklabels(x_order)

# Y-axis formatting
max_val = df["OxyRate"].max()
ax.set_ylim(0, max_val * 1.2)
ax.yaxis.set_major_locator(ticker.MultipleLocator(1)) # Adjusted tick frequency to 1

sns.despine(ax=ax)
fig.tight_layout()

# Save options   ---> remove # to get a downloadedble plot
plt.savefig("Respiration_Rate_by_Light.png", dpi=600)
#plt.savefig("Respiration_Rate_by_Light.pdf", dpi=600)
#plt.savefig("Respiration_Rate_by_Light.tiff", dpi=600)
#plt.savefig("Respiration_Rate_by_Light_transparent.svg", dpi=600, transparent=True)
#plt.savefig("Respiration_Rate_by_Light_white.svg", dpi=600, facecolor='white')
plt.show()