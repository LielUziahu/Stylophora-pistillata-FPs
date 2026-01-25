import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import statsmodels.api as sm
from statsmodels.formula.api import ols
import matplotlib.ticker as ticker
import matplotlib.patches as mpatches
from statsmodels.stats.multicomp import pairwise_tukeyhsd

# ==========================================
# 1. SETUP & DATA LOADING
# ==========================================
filename = "Temp_Assay.csv"
df = pd.read_csv(filename)

# Clean column names
df = df.rename(columns={"OxyRate nmol/mm2/min": "OxyRate"})
if 'color' in df.columns:
    df = df.rename(columns={"color": "morph"})

# --- NEW: Specific Defect Sample Removal ---
specific_point_to_remove = df[
    (df['SampleCode'] == 'plateD_NF_8.2_2') &
    (df['OxyRate'] == 2.098348468) # Using the precise value from Light_Assay.csv output
].copy()

if not specific_point_to_remove.empty:
    print("--- Temperature Samples Removed (Specific Defect) ---")
    for index, row in specific_point_to_remove.iterrows():
        print(f"Morph: {row['morph']}, SampleCode: {row['SampleCode']}, Temp: {row['temp']}, OxyRate: {row['OxyRate']:.2f}")
    print("---------------------------------------------------")
    # Remove the identified row(s) from the main DataFrame
    df = df.drop(specific_point_to_remove.index)
# ------------------------------------------

# --- DATA CLEANING (Existing) ---
# 1. Grouped 1.5 IQR Outlier Removal
# Initialize a list to store removed samples
removed_outliers_list = []

# Iterate through each unique combination of 'morph' and 'temp'
for (morph_group, temp_group), group_df in df.groupby(['morph', 'temp'], observed=False):
    Q1 = group_df['OxyRate'].quantile(0.25)
    Q3 = group_df['OxyRate'].quantile(0.75)
    IQR = Q3 - Q1

    lower_bound = Q1 - 1.5 * IQR
    upper_bound = Q3 + 1.5 * IQR

    # Identify outliers in the current group
    outliers = group_df[(group_df['OxyRate'] < lower_bound) | (group_df['OxyRate'] > upper_bound)]

    # Append outliers to the list
    if not outliers.empty:
        removed_outliers_list.append(outliers)

# Concatenate all removed outliers into a single DataFrame
removed_temp_samples_iqr = pd.DataFrame(columns=df.columns) # Initialize with columns to preserve dtypes
if removed_outliers_list:
    removed_temp_samples_iqr = pd.concat(removed_outliers_list)

# Filter the main DataFrame `df` to remove all rows that are present in `removed_temp_samples_iqr`
# We use `isin` on the index to efficiently remove rows.
df = df[~df.index.isin(removed_temp_samples_iqr.index)]

# Print removed samples if any
if not removed_temp_samples_iqr.empty:
    print("--- Temperature Samples Removed (Grouped 1.5 IQR) ---")
    for index, row in removed_temp_samples_iqr.iterrows():
        print(f"Morph: {row['morph']}, SampleCode: {row['SampleCode']}, Temp: {row['temp']}, OxyRate: {row['OxyRate']:.2f}")
    print("-----------------------------------------------------------")

# --- NEW: Filter out 'PF' morph (as observed in last execution) ---
target_morphs_temp = ['HF', 'NF']
df = df[df['morph'].isin(target_morphs_temp)].copy()
# ------------------------------------------------------------------

# 2. Reset Index
df = df.reset_index(drop=True)

# 3. Ensure Temp is Ordered Categorical
temp_order = [23, 27, 32]
df["temp"] = pd.Categorical(df["temp"], categories=temp_order, ordered=True)

# ==========================================
# 2. STATISTICS (Two-Way ANOVA)
# ==========================================

# 2.1 Calculate descriptive statistics (N, Mean, Std Dev) grouped by 'morph' and 'temp'
descriptive_stats_temp = df.groupby(['morph', 'temp'], observed=False)['OxyRate'].agg(
    N='count',
    Mean='mean',
    Std_Dev_Col='std'
).reset_index()
descriptive_stats_temp = descriptive_stats_temp.rename(columns={'Std_Dev_Col': 'Std Dev'})

# 2.2 Format the descriptive statistics table to 3 decimal places and print it
descriptive_stats_temp_formatted = descriptive_stats_temp.round(3)
print("\n=== Descriptive Statistics (Temperature Assay) ===")
print(descriptive_stats_temp_formatted.to_string(index=False))

# Model: OxyRate depends on morph and temp
model = ols("OxyRate ~ C(morph) + C(temp) + C(morph):C(temp)", data=df).fit()
anova_table = sm.stats.anova_lm(model, typ=2)

# 2.3 Extract F-value, df, and P-value for main effects and interaction
anova_results_data_temp = {
    'Source': ['C(morph)', 'C(temp)', 'C(morph):C(temp)']
}

for stat in ['df', 'F', 'PR(>F)']:
    col_name = 'F-value' if stat == 'F' else ('P-value' if stat == 'PR(>F)' else stat)
    anova_results_data_temp[col_name] = [
        anova_table.loc['C(morph)', stat] if 'C(morph)' in anova_table.index else np.nan,
        anova_table.loc['C(temp)', stat] if 'C(temp)' in anova_table.index else np.nan,
        anova_table.loc['C(morph):C(temp)', stat] if 'C(morph):C(temp)' in anova_table.index else np.nan
    ]

anova_results_df_temp = pd.DataFrame(anova_results_data_temp)

# 2.4 Format the ANOVA results table to 3 decimal places and print it
anova_results_df_temp_formatted = anova_results_df_temp.round(3)
print("\n=== Two-Way ANOVA Results (Temperature Assay) ===")
print(anova_results_df_temp_formatted.to_string(index=False))

# Post-hoc Test
df['group'] = df['morph'].astype(str) + "_" + df['temp'].astype(str)
tukey_result = pairwise_tukeyhsd(endog=df['OxyRate'], groups=df['group'], alpha=0.05)
print("\n=== Tukey HSD ===")
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

    # Prepare data for manual boxplot (list of arrays)
    plot_data = [sub[sub["temp"] == t]["OxyRate"].values for t in temp_order]

    positions = np.arange(len(temp_order)) + x_positions[morph]

    # A. Boxplots
    bp = ax.boxplot(
        plot_data,
        positions=positions,
        widths=box_width / 2,
        patch_artist=True,
        showfliers=False,
        zorder=2
    )

    # Style
    for patch in bp['boxes']:
        patch.set_facecolor(variant_colors[morph])
        patch.set_alpha(0.6)
        patch.set_edgecolor('black')
        patch.set_linewidth(1)
    for element in ['whiskers', 'caps', 'medians']:
        plt.setp(bp[element], color='black', linewidth=1)

    # B. Jittered points
    for i, t in enumerate(temp_order):
        y = sub[sub["temp"] == t]["OxyRate"]
        # Jitter X
        x = np.random.normal(positions[i], 0.04, size=len(y))

        ax.scatter(
            x, y,
            s=40, c=variant_colors[morph],
            edgecolors="black", linewidths=1,
            alpha=0.7, zorder=3
        )

    # C. Mean Diamonds (changed from median)
    means = sub.groupby("temp", observed=False)["OxyRate"].mean().reindex(temp_order)

    ax.plot(
        positions,
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
# 5. LABELS & LEGEND
# ==========================================
ax.set_ylabel("Respiration rate (nmol/mm$^2$/min)", fontsize=13)
ax.set_xlabel("Temperature (°C)", fontsize=13)
ax.set_title("Respiration Rate of settled spats by Temperature and morph", fontsize=14, pad=15)

# Fix X-axis
ax.set_xticks(np.arange(len(temp_order)))
ax.set_xticklabels(temp_order)

# Legend
legend_handles = []
for morph in morph_order:
    patch = mpatches.Patch(facecolor=variant_colors[morph], label=morph, edgecolor='black', linewidth=1, alpha=0.6)
    legend_handles.append(patch)
ax.legend(handles=legend_handles, labels=morph_order, title=None, frameon=False, loc="upper right")

# Y-axis scaling
max_val = df["OxyRate"].max()
ax.set_ylim(0, max_val * 1.2)
ax.yaxis.set_major_locator(ticker.MultipleLocator(1))

sns.despine(ax=ax)
fig.tight_layout()

# Save   ---> remove # to get a downloadedble plot
plt.savefig("Respiration_Rate_by_Temp.png", dpi=600)
#plt.savefig("Respiration_Rate_by_Temp.pdf", dpi=600)
#plt.savefig("Respiration_Rate_by_Temp.tiff", dpi=600)
#plt.savefig("Respiration_Rate_by_Temp_transparent.svg", dpi=600, transparent=True)
#plt.savefig("Respiration_Rate_by_Temp_white.svg", dpi=600, facecolor='white')

plt.show()
