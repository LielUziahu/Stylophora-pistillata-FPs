import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patches as mpatches
import matplotlib.ticker as ticker
from scipy import stats
import statsmodels.api as sm
from statsmodels.formula.api import ols
from statsmodels.stats.multicomp import pairwise_tukeyhsd

# ==========================================
# 1. SETUP & DATA LOADING
# ==========================================
filename = "pH_size.csv"
df = pd.read_csv(filename)

# Standardize columns
df.columns = df.columns.str.strip()
if 'color' in df.columns: df = df.rename(columns={'color': 'morph'})
if 'area mm^2' in df.columns: df = df.rename(columns={'area mm^2': 'Size'})
else:
    cols = [c for c in df.columns if 'area' in c.lower()]
    if cols: df = df.rename(columns={cols[0]: 'Size'})

df['morph'] = df['morph'].str.strip()
df['pH'] = df['pH'].astype(str).str.strip()

# Drop NaNs
df = df.dropna(subset=['Size', 'morph', 'pH'])

# ==========================================
# 2. FILTER & PREP FOR PLOT
# ==========================================
# Remove PF
df_plot = df[df['morph'] != 'PF'].copy()

# Order pH: 8.2 (Control) -> 7.8 -> 7.6 (Stress)
ph_order = ["8.2", "7.8", "7.6"]
# Filter to keep only these levels if they exist
df_plot = df_plot[df_plot["pH"].isin(ph_order)].copy()
df_plot["pH"] = pd.Categorical(df_plot["pH"], categories=ph_order, ordered=True)

# ==========================================
# 2.1 STATISTICS (Two-Way ANOVA & Tukey's HSD) - ADDED HERE
# ==========================================
print("=== Two-Way ANOVA Results (Size) ===")
# Model: Size depends on morph and pH
model = ols("Size ~ C(morph) + C(pH) + C(morph):C(pH)", data=df_plot).fit()
anova_table = sm.stats.anova_lm(model, typ=2)
display(anova_table) # Use display for DataFrame

# Tukey's HSD Post-hoc Test
df_plot['group'] = df_plot['morph'].astype(str) + "_" + df_plot['pH'].astype(str)
tukey_result = pairwise_tukeyhsd(endog=df_plot['Size'], groups=df_plot['group'], alpha=0.05)

print("\n=== Tukey's HSD Post-hoc Test Results ===")
print(tukey_result) # print for summary table

# ==========================================
# 3. PLOTTING
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

# COLORS FROM Light_size.py
variant_colors = {"HF": "#2ca02c", "NF": "#d62728"}
morph_order = ["HF", "NF"]

fig, ax = plt.subplots(figsize=(7, 6))

box_width = 0.8 # Changed to make boxes wider
offset = box_width / 4
x_positions = {"HF": -offset, "NF": +offset}

for morph in morph_order:
    sub = df_plot[df_plot["morph"] == morph]
    if sub.empty: continue

    positions = np.arange(len(ph_order)) + x_positions[morph]

    # A. Boxplots
    # Using ax.boxplot to mimic sns.boxplot style from Light_size.py
    data_to_plot = [sub[sub["pH"] == level]["Size"].values for level in ph_order]
    
    # Handle empty data arrays to avoid errors
    data_to_plot = [d if len(d) > 0 else np.array([]) for d in data_to_plot]
    
    bp = ax.boxplot(data_to_plot, positions=positions, widths=box_width/2, 
                    patch_artist=True, showfliers=False)
    
    # Style match: alpha=0.6, black edges
    for patch in bp['boxes']:
        patch.set_facecolor(variant_colors[morph])
        patch.set_alpha(0.6) # Matched to Light_size.py alpha
        patch.set_edgecolor('black')
        patch.set_linewidth(1)
        
    for element in ['whiskers', 'caps', 'medians']:
        plt.setp(bp[element], color='black', linewidth=1)

    # B. Jittered Points
    for i, level in enumerate(ph_order):
        y = sub[sub["pH"] == level]["Size"]
        if len(y) > 0:
            # Jitter
            x_vals = np.random.normal(i + x_positions[morph], 0.04, size=len(y))
            
            ax.scatter(
                x_vals, y,
                s=40, c=variant_colors[morph],
                edgecolors="black", linewidths=1,
                alpha=0.7, zorder=3
            )

    # C. Mean Diamonds (Changed from Median Diamonds)
    means = sub.groupby("pH")["Size"].mean().reindex(ph_order) # Changed to mean()
    valid = ~means.isna()
    if valid.any():
        ax.plot(
            np.arange(len(ph_order))[valid] + x_positions[morph],
            means[valid], # Changed to means
            marker="D", linestyle="--",
            markersize=8,
            markerfacecolor=variant_colors[morph],
            markeredgecolor="black",
            markeredgewidth=1.5,
            linewidth=2,
            color=variant_colors[morph],
            zorder=4
        )

# Labels & Legend
ax.set_ylabel("Spat Area ($mm^2$)", fontsize=13)
ax.set_xlabel("pH Conditions", fontsize=13)
ax.set_title("Size distribution of settled spats by pH and morph", fontsize=14, pad=15)

# Legend
handles = [mpatches.Patch(color=variant_colors[m], label=m, alpha=0.6, ec='black') for m in morph_order]
ax.legend(handles=handles, title=None, frameon=False, loc="upper right")

# Ticks
ax.set_xticks(np.arange(len(ph_order)))
ax.set_xticklabels(ph_order)

# Y-axis formatting
max_val = df_plot["Size"].max()
ax.set_ylim(0, max_val * 1.2)
ax.yaxis.set_major_locator(ticker.MultipleLocator(2.0)) # Changed to tick every 2 units

sns.despine(ax=ax)
fig.tight_layout()

plt.savefig("Size_by_pH_Morph_Style.png", dpi=600)
plt.savefig("Size_by_pH_Morph_Style.pdf", dpi=600)
plt.savefig("Size_by_pH_Morph_Style.tiff", dpi=600) # Added TIFF save option
plt.savefig("Size_by_pH_Morph_Style_transparent.svg", transparent=True) # Added SVG transparent save option
plt.savefig("Size_by_pH_Morph_Style_white_bg.svg", facecolor='white') # Added SVG white background save option
plt.show()