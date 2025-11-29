plt.rcParams.update({
    'font.family': 'serif',
    'axes.edgecolor': 'black', 'axes.linewidth': 1,
    'xtick.color': 'black', 'ytick.color': 'black',
    'text.color': 'black', 'axes.labelcolor': 'black',
    'legend.frameon': False,
    'ytick.labelsize': 11,
    'xtick.labelsize': 11
})

variant_colors = {'HF': '#2ca02c', 'NF': '#d62728'}
morph_order = ['HF', 'NF']

plt.figure(figsize=(6, 6))

# Box Plot
ax = sns.boxplot(
    data=df_resp, x='morph', y='rate umol/mm3/min', order=morph_order, hue='morph', palette=variant_colors,
    width=0.8, showfliers=False, linewidth=1, legend=False,
    boxprops=dict(edgecolor='black', alpha=0.5),
    whiskerprops=dict(color='black'), capprops=dict(color='black'), medianprops=dict(color='black')
)

# Jitter Points
sns.stripplot(
    data=df_resp, x='morph', y='rate umol/mm3/min', order=morph_order, hue='morph', palette=variant_colors,
    jitter=0.15, size=5, linewidth=1, edgecolor='black',
    alpha=0.6, legend=False, zorder=3
)

# Mean Diamond
means_resp = df_resp.groupby('morph')['rate umol/mm3/min'].mean().reindex(morph_order)
mean_colors_resp = [variant_colors[m] for m in morph_order]
plt.scatter(range(len(morph_order)), means_resp, c=mean_colors_resp, marker='D', s=60, edgecolors='black', linewidths=2, alpha=0.7, zorder=4)

# 5. ANNOTATIONS (Letters & P-value)
letters_map = {'HF': 'a', 'NF': 'a'} # Assuming no significant difference based on p_val_resp > 0.05, so same letter
y_offset_percent = 0.05

max_val_global_resp = df_resp['rate umol/mm3/min'].max()

for i, m in enumerate(morph_order):
    subset_resp = df_resp[df_resp['morph'] == m]['rate umol/mm3/min']
    highest_point_resp = subset_resp.max()

    # Add offset
    pos_y_resp = highest_point_resp + (max_val_global_resp * y_offset_percent)

    letter = letters_map[m]
    plt.text(x=i, y=pos_y_resp, s=letter, ha='center', va='bottom', size=12, weight='bold')

# P-value (Bottom Right, under axis)
p_text_resp = f"p = {p_val_resp:.4f}" if p_val_resp >= 0.0001 else "p < 0.0001"
plt.text(1.0, -0.1, s=p_text_resp, transform=ax.transAxes,
         ha='right', va='top', size=10, style='italic', color='black')

# 6. LABELS & TITLE
plt.ylabel("Respiration rate ($µ$mol/mm$^3$/min)", size=13)
plt.xlabel("")
plt.title("Respiration rate distribution of different planulae morphs", size=14, pad=15)

# Adjust Y limit to fit annotations and data
plt.ylim(0, max_val_global_resp * 1.15)
plt.yticks(np.arange(0, round(max_val_global_resp * 1.15) + 0.2, 1.0))

sns.despine()
plt.tight_layout()
plt.savefig("Respiration_Rate_by_Morph.pdf", dpi=600)
plt.savefig("Respiration_Rate_by_Morph.tiff", dpi=600)
plt.savefig("Respiration_Rate_by_Morph.png", dpi=600)
plt.show()