import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
import statsmodels.api as sm
from statsmodels.formula.api import ols
import matplotlib.ticker as ticker
import matplotlib.patches as mpatches # Import for creating legend patches
from statsmodels.stats.multicomp import pairwise_tukeyhsd # Import for Tukey's HSD

# ==========================================
# 1. SETUP & DATA LOADING
# ==========================================
filename = "pH_Spat_Resp.csv"
df = pd.read_csv(filename)

# Normalize / rename the response column once (safer than repeating the long name everywhere)
df = df.rename(columns={"OxyRate nmol/mm2/min": "OxyRate"})

# Keep only HF / NF
target_morphs = ["HF", "NF"]
df = df[df["morph"].isin(target_morphs)].copy()

# Ensure pH is clean and ordered (robust to 7.60 vs 7.6)
df["pH"] = df["pH"].astype(str).str.strip()
ph_order = ["8.2", "7.8", "7.6"] # Reversed order
df["pH"] = pd.Categorical(df["pH"], categories=ph_order, ordered=True)

# Optional: drop rows with unexpected pH values (prevents empty categories / odd plotting)
df = df[df["pH"].isin(ph_order)].copy()

# ==========================================
# 2. STATISTICS (Two-Way ANOVA & Tukey's HSD)
# ==========================================
model = ols("OxyRate ~ C(morph) + C(pH) + C(morph):C(pH)", data=df).fit()
anova_table = sm.stats.anova_lm(model, typ=2)

print("=== Two-Way ANOVA Results ===")
print(anova_table)

p_morph = anova_table.loc["C(morph)", "PR(>F)"]
p_ph = anova_table.loc["C(pH)", "PR(>F)"]
p_inter = anova_table.loc["C(morph):C(pH)", "PR(>F)"]

def format_pvalue(p):
    return "< 0.0001" if p < 0.0001 else f"= {p:.3f}"

# 2.1 Tukey's HSD Post-hoc Test
# Create a "group" column for interaction comparisons
df['group'] = df['morph'].astype(str) + "_" + df['pH'].astype(str)

# Run Tukey's HSD
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

# Consistent dodge so layers align
dodge_width = 0.35

fig, ax = plt.subplots(figsize=(7, 6))

# ==========================================
# 4. PLOT CONSTRUCTION (MANUAL ALIGNMENT)
# ==========================================
box_width = 0.6
offset = box_width / 4   # half of half-width → perfect centering

x_positions = {
    "HF": -offset,
    "NF": +offset
}

for morph in morph_order:
    sub = df[df["morph"] == morph]

    # A. Boxplots (no hue)
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

    # B. Jittered points (exact same positions)
    for i, ph in enumerate(ph_order):
        y = sub[sub["pH"] == ph]["OxyRate"]
        x = np.random.normal(i + x_positions[morph], 0.04, size=len(y))
        ax.scatter(
            x, y,
            s=40, c=variant_colors[morph],
            edgecolors="black", linewidths=1,
            alpha=0.7, zorder=3
        )

    # C. Median diamonds (PERFECTLY centered)
    medians = sub.groupby("pH")["OxyRate"].median().reindex(ph_order)

    ax.plot(
        np.arange(len(ph_order)) + x_positions[morph],
        medians,
        marker="D", linestyle="--", 
        markersize=8,
        markerfacecolor=variant_colors[morph],
        markeredgecolor="black",   # ✅ black outline
        markeredgewidth=1.5,
        linewidth=2,
        color=variant_colors[morph],
        zorder=4
    )

# ==========================================
# 5. LABELS, LEGEND & ANNOTATION
# ==========================================
ax.set_ylabel("Respiration rate (nmol/mm$^2$/min)", fontsize=13)
ax.set_xlabel("pH Level", fontsize=13)
ax.set_title("Respiration Rate of seettled spats by pH and morph", fontsize=14, pad=15)

# Custom Legend (Top Right) for box plots
legend_handles = []
for morph in morph_order:
    # Create a colored patch resembling a box, with matching alpha
    patch = mpatches.Patch(color=variant_colors[morph], label=morph, edgecolor='black', linewidth=1, alpha=0.6)
    legend_handles.append(patch)
ax.legend(handles=legend_handles, labels=morph_order, title=None, frameon=False, loc="upper right")

# Removed the ANOVA text box as requested.

# Y-axis formatting
max_val = df["OxyRate"].max()
ax.set_ylim(0, max_val * 1.2)
ax.yaxis.set_major_locator(ticker.MultipleLocator(1))

sns.despine(ax=ax)
fig.tight_layout()

# Save options
plt.savefig("Respiration_Rate_by_pH.png", dpi=600)
plt.savefig("Respiration_Rate_by_pH.tiff", dpi=600)
plt.savefig("Respiration_Rate_by_pH_white.svg", dpi=600, facecolor='white')
plt.savefig("Respiration_Rate_by_pH_transparent.svg", dpi=600, transparent=True)
plt.savefig("Respiration_Rate_by_pH.pdf", dpi=600)

plt.show()