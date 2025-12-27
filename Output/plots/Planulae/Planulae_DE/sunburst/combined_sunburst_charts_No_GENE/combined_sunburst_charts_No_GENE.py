#-------------------------------------------------------------#
#----------------Sunburst Without Gene Circle-----------------#
#----Grey Root + Synced Rotation + No Black Wedges + EXPORTS--#
#   SVG (transparent) + SVG (white) + PNG/TIFF 600dpi + PDF    #
#-------------------------------------------------------------#

import os
import copy
import pandas as pd
import plotly.express as px
import plotly.io as pio
from plotly.subplots import make_subplots
from PIL import Image

# -----------------------------
# Settings
# -----------------------------
file_path = "Categorized_Gene_Table_HF_NF_with_LFC_no_unknown.csv"

ROOT_COLOR = "#BFBFBF"        # grey center
ROTATION_DEG = 0              # same rotation for HF and NF
INSIDE_FONT_SIZE = 22

W, H = 1920, 1080
DPI = 600

# One source of truth for ALL output names:
BASE_OUT = "combined_sunburst_charts_No_GENE"

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
# 1) Load + clean
# -----------------------------
print(f"[1/5] Loading CSV: {file_path}")
df = pd.read_csv(file_path)
print(f"     Loaded shape: {df.shape}")

df["LFC"] = pd.to_numeric(df["LFC"], errors="coerce")
before_drop = df.shape[0]
df = df.dropna(subset=["Morph", "Category", "Sub-Category", "LFC"]).copy()
after_drop = df.shape[0]
print(f"     After dropna: {after_drop} rows (dropped {before_drop - after_drop})")

df["Morph"] = df["Morph"].astype(str).str.strip()
df["Category"] = df["Category"].astype(str).str.strip()
df["Sub-Category"] = df["Sub-Category"].astype(str).str.strip()
df["Gene Symbol"] = df["Gene Symbol"].astype(str).str.strip()

df["absLFC"] = df["LFC"].abs()
print(f"     Morphs: {sorted(df['Morph'].unique().tolist())}")
print(f"     Categories: {df['Category'].nunique()} | Sub-categories: {df['Sub-Category'].nunique()}")
print(f"     Σ absLFC (sanity): {df['absLFC'].sum():.3f}")

# -----------------------------
# 2) Aggregate (removes gene level)
# -----------------------------
print("[2/5] Aggregating to Morph -> Category -> Sub-Category")
agg = (
    df.groupby(["Morph", "Category", "Sub-Category"], as_index=False)
      .agg(
          absLFC=("absLFC", "sum"),
          n_genes=("Gene Symbol", "nunique"),
      )
)
print(f"     Aggregated rows: {agg.shape[0]}")

print("     Computing top genes (top 5 by absLFC) per sub-category")
top5 = (
    df.sort_values(["Morph", "Category", "Sub-Category", "absLFC"],
                   ascending=[True, True, True, False])
      .groupby(["Morph", "Category", "Sub-Category"])["Gene Symbol"]
      .apply(lambda s: ", ".join(s.head(5).astype(str)))
      .reset_index(name="top_genes")
)

agg = agg.merge(top5, on=["Morph", "Category", "Sub-Category"], how="left")
agg["Morph"] = agg["Morph"].astype(str).str.strip()
agg["Category"] = agg["Category"].astype(str).str.strip()
agg["Sub-Category"] = agg["Sub-Category"].astype(str).str.strip()
agg["top_genes"] = agg["top_genes"].fillna("")
print(f"     Σ absLFC over agg (sanity): {agg['absLFC'].sum():.3f}")

# -----------------------------
# 3) Colors
# -----------------------------
print("[3/5] Building category color map")
all_categories = sorted(agg["Category"].unique())
palette = px.colors.qualitative.Plotly
category_color_map = {c: palette[i % len(palette)] for i, c in enumerate(all_categories)}
print(f"     Mapped categories: {len(category_color_map)}")

# -----------------------------
# 4) Build combined subplot figure (one per morph)
# -----------------------------
morphs = sorted(agg["Morph"].unique())
print(f"[4/5] Building sunbursts for morphs: {morphs}")

fig = make_subplots(
    rows=1,
    cols=len(morphs),
    specs=[[{"type": "domain"}] * len(morphs)],
    subplot_titles=[""] * len(morphs),
)

for col_i, morph in enumerate(morphs, start=1):
    agg_m = agg[agg["Morph"] == morph].copy()
    agg_m["Root"] = morph  # Option B: center circle

    print(f"     -> Morph '{morph}': subcat rows={agg_m.shape[0]} | Σ absLFC={agg_m['absLFC'].sum():.3f}")

    fig_m = px.sunburst(
        agg_m,
        path=["Root", "Category", "Sub-Category"],
        values="absLFC",
        custom_data=["absLFC", "n_genes", "top_genes"],
    )

    trace = fig_m.data[0]
    ids = list(trace.ids)
    print(f"        Nodes rendered (root+categories+subcats): {len(ids)}")

    marker_colors = []
    for node_id in ids:
        parts = [p.strip() for p in str(node_id).split("/")]

        # Root
        if len(parts) == 1:
            marker_colors.append(ROOT_COLOR)
        # Category / Sub-Category -> inherit category color
        else:
            cat = parts[1]
            marker_colors.append(category_color_map.get(cat, ROOT_COLOR))

    fig_m.update_traces(
        marker=dict(colors=marker_colors),
        hovertemplate=(
            "<b>%{label}</b><br>"
            "Σ|LFC|: %{customdata[0]:.2f}<br>"
            "n genes: %{customdata[1]}<br>"
            "Top genes: %{customdata[2]}"
            "<extra></extra>"
        ),
        insidetextfont=dict(size=INSIDE_FONT_SIZE),
        insidetextorientation="radial",
        rotation=ROTATION_DEG,
        sort=False
    )

    fig.add_trace(fig_m.data[0], row=1, col=col_i)

# -----------------------------
# 5) Layout + export
# -----------------------------
print("[5/5] Rendering + exporting")
fig.update_layout(
    title_text="Gene Expression in <i>S. pistillata</i> planulae morphs",
    title_x=0.5,
    title_font_size=36,
    width=W,
    height=H,
    font_family="serif",
    margin=dict(t=120, l=20, r=20, b=20),
)

fig.show()

# HTML export (interactive)
fig.write_html(HTML_OUT)
print(f"Saved HTML: {HTML_OUT}")

def with_bg(fig_in, transparent: bool):
    f = copy.deepcopy(fig_in)
    if transparent:
        f.update_layout(paper_bgcolor="rgba(0,0,0,0)", plot_bgcolor="rgba(0,0,0,0)")
    else:
        f.update_layout(paper_bgcolor="white", plot_bgcolor="white")
    return f

# SVG exports
with_bg(fig, True).write_image(SVG_TRANSPARENT_OUT, format="svg", width=W, height=H, scale=1)
print(f"Saved SVG (transparent): {SVG_TRANSPARENT_OUT}")

with_bg(fig, False).write_image(SVG_WHITE_OUT, format="svg", width=W, height=H, scale=1)
print(f"Saved SVG (white): {SVG_WHITE_OUT}")

# PNG export + 600dpi metadata
tmp_png = "_tmp_export.png"
with_bg(fig, False).write_image(tmp_png, format="png", width=W, height=H, scale=1)
im = Image.open(tmp_png)
im.save(PNG_600_OUT, dpi=(DPI, DPI))
im.close()
os.remove(tmp_png)
print(f"Saved PNG ({DPI} dpi metadata): {PNG_600_OUT}")

# TIFF export + 600dpi metadata
tmp_png = "_tmp_export.png"
with_bg(fig, False).write_image(tmp_png, format="png", width=W, height=H, scale=1)
im = Image.open(tmp_png).convert("RGB")
im.save(TIFF_600_OUT, dpi=(DPI, DPI), compression="tiff_lzw")
im.close()
os.remove(tmp_png)
print(f"Saved TIFF ({DPI} dpi metadata): {TIFF_600_OUT}")

# PDF export (vector; dpi not meaningful)
with_bg(fig, False).write_image(PDF_OUT, format="pdf", width=W, height=H, scale=1)
print(f"Saved PDF: {PDF_OUT}")
