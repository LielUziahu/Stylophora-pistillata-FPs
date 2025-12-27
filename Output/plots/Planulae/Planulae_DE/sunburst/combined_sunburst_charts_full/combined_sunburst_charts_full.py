#-------------------------------------------------------------#
#----------------Sunburst (WITH GENES)------------------------#
#----Synced Rotation HF/NF + No Black Wedges + EXPORTS--------#
#   HTML + SVG (transparent/white) + PNG/TIFF 600dpi + PDF     #
#-------------------------------------------------------------#

import os
import copy
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import plotly.io as pio
from plotly.subplots import make_subplots
from PIL import Image

# -----------------------------
# Settings
# -----------------------------
file_path = "Categorized_Gene_Table_HF_NF_with_LFC_no_unknown.csv"

ROTATION_DEG = 0          # same start angle for HF & NF
INSIDE_FONT_SIZE = 22

W, H = 1920, 1080
DPI = 600

BASE_OUT = "combined_sunburst_charts_WITH_GENE"

HTML_OUT            = f"{BASE_OUT}.html"
SVG_TRANSPARENT_OUT = f"{BASE_OUT}_transparent.svg"
SVG_WHITE_OUT       = f"{BASE_OUT}_white.svg"
PNG_600_OUT         = f"{BASE_OUT}_600dpi.png"
TIFF_600_OUT        = f"{BASE_OUT}_600dpi.tiff"
PDF_OUT             = f"{BASE_OUT}.pdf"

# Kaleido must be active for static export
if pio.kaleido.scope is None:
    raise RuntimeError("Kaleido is not active (pio.kaleido.scope is None).")

# -----------------------------
# 1) Data Loading and Cleaning
# -----------------------------
try:
    df = pd.read_csv(file_path)
except FileNotFoundError:
    print("File not found. Using dummy data.")
    df = pd.DataFrame({
        "Morph": ["HF", "HF", "NF", "NF", "HF"],
        "Category": ["Adhesion", "Adhesion", "Immunity", "Immunity", "Signaling"],
        "Sub-Category": ["General", "Receptors", "Cytokines", "Receptors", "Kinases"],
        "Gene Symbol": ["GENE1", "GENE2", "GENE3", "GENE4", "GENE5"],
        "LFC": [2.5, -1.8, 3.1, -2.0, 1.5],
        "Description": ["Desc1", "Desc2", "Desc3", "Desc4", "Desc5"],
    })

df["LFC"] = pd.to_numeric(df["LFC"], errors="coerce")
df = df.dropna(subset=["Morph", "Category", "Sub-Category", "Gene Symbol", "LFC", "Description"]).copy()

# Normalize strings (strip whitespace)
df["Morph"] = df["Morph"].astype(str).str.strip()
df["Category"] = df["Category"].astype(str).str.strip()
df["Sub-Category"] = df["Sub-Category"].astype(str).str.strip()
df["Gene Symbol"] = df["Gene Symbol"].astype(str).str.strip()
df["Description"] = df["Description"].astype(str).str.strip()

print("Data loaded and cleaned.")

# -----------------------------
# 2) Prepare Data for go.Sunburst
#    (Important fix for rotation comparability)
#    -> enforce consistent ordering of categories/subcats/genes
# -----------------------------
# Global ordering ensures HF and NF build wedges in the same sequence
category_order = sorted(df["Category"].unique())
subcat_order = (
    df[["Category", "Sub-Category"]].drop_duplicates()
      .sort_values(["Category", "Sub-Category"])
)
subcat_order_map = {c: [] for c in category_order}
for c in category_order:
    subcat_order_map[c] = subcat_order[subcat_order["Category"] == c]["Sub-Category"].tolist()

# Build the full sunburst arrays
sunburst_ids = []
sunburst_labels = []
sunburst_parents = []
sunburst_values = []
sunburst_hover_lfc_text = []
sunburst_hover_desc_text = []
sunburst_node_category = []

# Use a stable morph order (HF then NF if present)
morph_order = sorted(df["Morph"].unique().tolist())
if "HF" in morph_order and "NF" in morph_order:
    morph_order = ["HF", "NF"] + [m for m in morph_order if m not in ("HF", "NF")]

for morph in morph_order:
    morph_id = morph
    sunburst_ids.append(morph_id)
    sunburst_labels.append(morph)
    sunburst_parents.append("")
    sunburst_values.append(0)  # allowed; keeps center visible
    sunburst_hover_lfc_text.append("")
    sunburst_hover_desc_text.append("")
    sunburst_node_category.append(None)

    df_morph = df[df["Morph"] == morph]

    # enforce global category order
    for category in category_order:
        df_category = df_morph[df_morph["Category"] == category]
        if df_category.empty:
            continue

        category_id = f"{morph_id}_{category}"
        sunburst_ids.append(category_id)
        sunburst_labels.append(category)
        sunburst_parents.append(morph_id)
        sunburst_values.append(0)
        sunburst_hover_lfc_text.append("")
        sunburst_hover_desc_text.append("")
        sunburst_node_category.append(category)

        # enforce per-category subcat order
        for subcategory in subcat_order_map.get(category, []):
            df_subcategory = df_category[df_category["Sub-Category"] == subcategory]
            if df_subcategory.empty:
                continue

            subcategory_id = f"{category_id}_{subcategory}"
            sunburst_ids.append(subcategory_id)
            sunburst_labels.append(subcategory)
            sunburst_parents.append(category_id)
            sunburst_values.append(0)
            sunburst_hover_lfc_text.append("")
            sunburst_hover_desc_text.append("")
            sunburst_node_category.append(category)

            # enforce stable gene order within subcategory (by abs(LFC) desc, then name)
            df_subcategory = df_subcategory.assign(absLFC=df_subcategory["LFC"].abs())
            df_subcategory = df_subcategory.sort_values(["absLFC", "Gene Symbol"], ascending=[False, True])

            for _, row in df_subcategory.iterrows():
                gene_symbol = row["Gene Symbol"]
                gene_id = f"{subcategory_id}_{gene_symbol}"
                sunburst_ids.append(gene_id)
                sunburst_labels.append(gene_symbol)
                sunburst_parents.append(subcategory_id)
                sunburst_values.append(float(abs(row["LFC"])))
                sunburst_hover_lfc_text.append(f"LFC: {row['LFC']:.2f}")
                sunburst_hover_desc_text.append(f"Description: {row['Description']}")
                sunburst_node_category.append(category)

sunburst_data_df = pd.DataFrame({
    "ids": sunburst_ids,
    "labels": sunburst_labels,
    "parents": sunburst_parents,
    "values": sunburst_values,
    "hover_lfc_text": sunburst_hover_lfc_text,
    "hover_desc_text": sunburst_hover_desc_text,
    "node_category": sunburst_node_category
})
print("Sunburst data DataFrame created.")

# -----------------------------
# 3) Generate Individual go.Sunburst Charts (one per morph)
# -----------------------------
morph_go_charts = {}

unique_morphs_from_data = sunburst_data_df["ids"][sunburst_data_df["parents"] == ""].unique()

all_unique_categories = category_order
palette = px.colors.qualitative.Plotly
category_color_map = {c: palette[i % len(palette)] for i, c in enumerate(all_unique_categories)}

for morph_name in unique_morphs_from_data:
    df_filtered_for_morph = sunburst_data_df[
        (sunburst_data_df["ids"] == morph_name) |
        (sunburst_data_df["ids"].str.startswith(f"{morph_name}_"))
    ].copy()

    node_colors = []
    for category_val in df_filtered_for_morph["node_category"]:
        if category_val is None:
            node_colors.append("#cccccc")
        else:
            node_colors.append(category_color_map.get(str(category_val).strip(), "#cccccc"))

    customdata_list = df_filtered_for_morph[["hover_lfc_text", "hover_desc_text"]].values.tolist()

    hovertemplate = (
        "<b>%{label}</b><br>"
        "%{customdata[0]}<br>"
        "%{customdata[1]}<extra></extra>"
    )

    sunburst_trace = go.Sunburst(
        ids=df_filtered_for_morph["ids"],
        labels=df_filtered_for_morph["labels"],
        parents=df_filtered_for_morph["parents"],
        values=df_filtered_for_morph["values"],
        customdata=customdata_list,
        hovertemplate=hovertemplate,
        marker=dict(colors=node_colors),
        insidetextfont=dict(size=INSIDE_FONT_SIZE),
        rotation=ROTATION_DEG,   # <-- synced rotation
        sort=False               # <-- critical: prevent Plotly re-sorting wedges
    )

    fig_m = go.Figure(sunburst_trace)
    fig_m.update_layout(title_text=f"{morph_name} Gene Expression Sunburst Chart")
    morph_go_charts[morph_name] = fig_m

print(f"Generated {len(morph_go_charts)} individual go.Sunburst charts.")

# -----------------------------
# 4) Combine and Display Sunburst Charts
# -----------------------------
num_charts = len(morph_go_charts)
fig = make_subplots(
    rows=1,
    cols=num_charts,
    specs=[[{"type": "domain"}] * num_charts],
    subplot_titles=[""] * num_charts
)

# Keep subplot order stable using morph_order
ordered_keys = [m for m in morph_order if m in morph_go_charts]
for i, morph_name in enumerate(ordered_keys):
    fig.add_trace(morph_go_charts[morph_name].data[0], row=1, col=i + 1)

fig.update_layout(
    title_text="Gene Expression in <i>S. pistillata</i> planulae morphs",
    title_x=0.5,
    title_font_size=36,
    width=W,
    height=H,
    font_family="serif",
    margin=dict(t=120, l=20, r=20, b=20)
)

fig.show()
print("Combined Sunburst charts displayed.")

# -----------------------------
# 5) Exports (same set as before)
# -----------------------------
fig.write_html(HTML_OUT)
print(f"Saved HTML: {HTML_OUT}")

def with_bg(fig_in, transparent: bool):
    f = copy.deepcopy(fig_in)
    if transparent:
        f.update_layout(paper_bgcolor="rgba(0,0,0,0)", plot_bgcolor="rgba(0,0,0,0)")
    else:
        f.update_layout(paper_bgcolor="white", plot_bgcolor="white")
    return f

# SVG: transparent + white
with_bg(fig, True).write_image(SVG_TRANSPARENT_OUT, format="svg", width=W, height=H, scale=1)
print(f"Saved SVG (transparent): {SVG_TRANSPARENT_OUT}")

with_bg(fig, False).write_image(SVG_WHITE_OUT, format="svg", width=W, height=H, scale=1)
print(f"Saved SVG (white): {SVG_WHITE_OUT}")

# PNG: export then embed 600 dpi metadata
tmp_png = "_tmp_export.png"
with_bg(fig, False).write_image(tmp_png, format="png", width=W, height=H, scale=1)
im = Image.open(tmp_png)
im.save(PNG_600_OUT, dpi=(DPI, DPI))
im.close()
os.remove(tmp_png)
print(f"Saved PNG ({DPI} dpi metadata): {PNG_600_OUT}")

# TIFF: export then convert to TIFF with 600 dpi metadata
tmp_png = "_tmp_export.png"
with_bg(fig, False).write_image(tmp_png, format="png", width=W, height=H, scale=1)
im = Image.open(tmp_png).convert("RGB")
im.save(TIFF_600_OUT, dpi=(DPI, DPI), compression="tiff_lzw")
im.close()
os.remove(tmp_png)
print(f"Saved TIFF ({DPI} dpi metadata): {TIFF_600_OUT}")

# PDF: vector (dpi not meaningful)
with_bg(fig, False).write_image(PDF_OUT, format="pdf", width=W, height=H, scale=1)
print(f"Saved PDF: {PDF_OUT}")
