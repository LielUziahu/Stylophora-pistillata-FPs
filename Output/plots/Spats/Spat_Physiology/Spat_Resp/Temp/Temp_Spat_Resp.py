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
filename = "Temp_Spat_Resp.csv"
df = pd.read_csv(filename)

# Clean column names
df = df.rename(columns={"OxyRate nmol/mm2/min": "OxyRate"})
if 'color' in df.columns:
    df = df.rename(columns={"color": "morph"})

# --- DATA CLEANING ---
# 1. Filter Outliers (0.1 - 2.07)
# (Your file is already clean, but this is a safety line)
df = df[(df['OxyRate'] >= 0.1) & (df['OxyRate'] <= 2.07)]

# 2. Reset Index
df = df.reset_index(drop=True)

# 3. Ensure Temp is Ordered Categorical
temp_order = [23, 27, 32]
df["temp"] = pd.Categorical(df["temp"], categories=temp_order, ordered=True)

# ==========================================
# 2. STATISTICS (Two-Way ANOVA)
# ==========================================
model = ols("OxyRate ~ C(morph) + C(temp) + C(morph):C(temp)", data=df).fit()
anova_table = sm.stats.anova_lm(model, typ=2)

print("=== Two-Way ANOVA Results ===")
print(anova_table)

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
    means = sub.groupby("temp")["OxyRate"].mean().reindex(temp_order)

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
    patch = mpatches.Patch(color=variant_colors[morph], label=morph, edgecolor='black', linewidth=1, alpha=0.6)
    legend_handles.append(patch)
ax.legend(handles=legend_handles, labels=morph_order, title=None, frameon=False, loc="upper right")

# Y-axis scaling
ax.set_ylim(0, df["OxyRate"].max() * 1.2)
ax.yaxis.set_major_locator(ticker.MultipleLocator(1))

sns.despine(ax=ax)
fig.tight_layout()

# Save
plt.savefig("Respiration_Rate_by_Temp.png", dpi=600)
plt.savefig("Respiration_Rate_by_Temp.pdf", dpi=600)
plt.savefig("Respiration_Rate_by_Temp.tiff", dpi=600)
plt.savefig("Respiration_Rate_by_Temp_transparent.svg", dpi=600, transparent=True)
plt.savefig("Respiration_Rate_by_Temp_white.svg", dpi=600, facecolor='white')

plt.show()