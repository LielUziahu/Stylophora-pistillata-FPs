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
filename = "Light_size.csv"
df = pd.read_csv(filename)

# --- Handle 'morph' and 'color' columns consistently ---
# The file has 'color' column for morphs (HF/NF/PF)
if 'color' in df.columns:
    # Rename 'color' to 'morph'
    df = df.rename(columns={'color': 'morph'})

# Drop rows with missing values in critical columns
df = df.dropna(subset=['area mm^2', 'morph', 'LightDepth'])

# --- FILTERING ---
# 1. Remove PF (Partial Fluorescence)
df = df[df['morph'] != 'PF'].copy()

# 2. Rename Target Variable
# Rename 'area mm^2' to 'Size' to match the script logic
df = df.rename(columns={"area mm^2": "Size"})

# Clean 'morph' strings
df['morph'] = df['morph'].str.strip()

# Clean and Order the LightDepth column
# Logic: High Light (10m) -> Low Light (50m) -> No Light (Dark)
df["LightDepth"] = df["LightDepth"].astype(str).str.strip()
light_order = ["10m", "50m", "Dark"]
df["LightDepth"] = pd.Categorical(df["LightDepth"], categories=light_order, ordered=True)

# Drop rows with unexpected LightDepth values
df = df[df["LightDepth"].isin(light_order)].copy()
df = df.reset_index(drop=True)

# ==========================================
# 2. STATISTICS (Two-Way ANOVA & Tukey's HSD)
# ==========================================
print("=== Two-Way ANOVA Results (Size) ===")
# Model: Size depends on morph and LightDepth
model = ols("Size ~ C(morph) + C(LightDepth) + C(morph):C(LightDepth)", data=df).fit()
anova_table = sm.stats.anova_lm(model, typ=2)
print(anova_table)

# 2.1 Tukey's HSD Post-hoc Test
df['group'] = df['morph'].astype(str) + "_" + df['LightDepth'].astype(str)
tukey_result = pairwise_tukeyhsd(endog=df['Size'], groups=df['group'], alpha=0.05)

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

# Colors matching the previous file
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
    if sub.empty: continue

    # A. Boxplots
    sns.boxplot(
        data=sub, x=x_col, y="Size",
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
        y = sub[sub[x_col] == level]["Size"]
        if len(y) > 0:
            # Add random jitter to x
            x_vals = np.random.normal(i + x_positions[morph], 0.04, size=len(y))

            ax.scatter(
                x_vals, y,
                s=40, c=variant_colors[morph],
                edgecolors="black", linewidths=1,
                alpha=0.7, zorder=3
            )

    # C. Median diamonds
    medians = sub.groupby(x_col)["Size"].median().reindex(x_order)

    # Handle potential missing data for lines
    valid_indices = ~medians.isna()
    if valid_indices.any():
        ax.plot(
            np.arange(len(x_order))[valid_indices] + x_positions[morph],
            medians[valid_indices],
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
ax.set_ylabel("Spat Area ($mm^2$)", fontsize=13)
ax.set_xlabel("Light/Depth Conditions", fontsize=13)
ax.set_title("Size distribution of settled spats by Light and morph", fontsize=14, pad=15)

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
max_val = df["Size"].max()
ax.set_ylim(0, max_val * 1.2)
# Adjusted tick frequency to 0.5 for size data (since values are smaller than respiration)
ax.yaxis.set_major_locator(ticker.MultipleLocator(3)) # Changed from 0.5 to 3

sns.despine(ax=ax)
fig.tight_layout()

# Save options
plt.savefig("Size_by_Light_Morph_Style.png", dpi=600)
plt.savefig("Size_by_Light_Morph_Style.pdf", dpi=600)
plt.savefig("Size_by_Light_Morph_Style.tiff", dpi=600) # Added TIFF save option
plt.savefig("Size_by_Light_Morph_Style_transparent.svg", transparent=True) # Added SVG transparent save option
plt.savefig("Size_by_Light_Morph_Style_white_bg.svg", facecolor='white') # Added SVG white background save option
plt.show()