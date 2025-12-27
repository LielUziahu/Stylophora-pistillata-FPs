#-------------------------------------------------------------#
#----------Sunburst WITHOUT Gene Circle AND WITHOUT-----------#
#----------------------Sub-Categories-------------------------#
#----Grey Root + Synced Rotation + Multi-format Export--------#
#   SVG (transparent) + SVG (white) + PNG/TIFF 600dpi + PDF    #
#-------------------------------------------------------------#

import os
import copy
import pandas as pd
import plotly.express as px
import plotly.io as pio
import kaleido
kaleido.get_chrome_sync()
from plotly.subplots import make_subplots
from PIL import Image

# -----------------------------
# Settings
# -----------------------------
file_path = "Categorized_Gene_Table_HF_NF_with_LFC_no_unknown.csv"

ROOT_COLOR = "#BFBFBF"
ROTATION_DEG = 0
INSIDE_FONT_SIZE = 22

W, H = 1920, 1080
DPI = 600

OUT_HTML = "combined_sunburst_charts_No_GENE_No_SUBCAT.html"

OUT_SVG_TRANSPARENT = "combined_sunburst_charts_No_GENE_No_SUBCAT_transparent.svg"
OUT_SVG_WHITE       = "combined_sunburst_charts_No_GENE_No_SUBCAT_white.svg"
OUT_PNG_600         = "combined_sunburst_charts_No_GENE_No_SUBCAT_600dpi.png"
OUT_TIFF_600        = "combined_sunburst_charts_No_GENE_No_SUBCAT_600dpi.tiff"
OUT_PDF             = "combined_sunburst_charts_No_GENE_No_SUBCAT.pdf"

# -----------------------------
# Ensure Chrome for Kaleido v1+
# -----------------------------
def ensure_chrome_for_kaleido():
    """
    Kaleido v1+ needs Chrome/Chromium. If it's missing, install via kaleido helper.
    """
    try:
        import kaleido
        # Will raise if Chrome isn't available
        from kaleido.errors import ChromeNotFoundError
        try:
            # Quick probe: ask kaleido to launch a minimal session indirectly by checking a calc
            # Plotly will do it on first write_image; we just pre-install if needed.
            return True
        except ChromeNotFoundError:
            pass
    except Exception:
        # If kaleido isn't importable, plotly would fail anyway
        raise RuntimeError("Kaleido is not available. Install: pip install -U kaleido")

    # Install Chrome using kaleido's official helper
    try:
        import kaleido
        # Prefer sync in notebooks
        kaleido.get_chrome_sync()
        return True
    except Exception as e:
        raise RuntimeError(
            "Chrome is required for Kaleido image export but could not be installed automatically.\n"
            "Try running in a terminal:\n"
            "  plotly_get_chrome\n"
            "or in a notebook cell:\n"
            "  !plotly_get_chrome\n"
            f"\nDetails: {repr(e)}"
        )

# -----------------------------
# 1) Load + clean
# -----------------------------
df = pd.read_csv(file_path)

df["LFC"] = pd.to_numeric(df["LFC"], errors="coerce")
df = df.dropna(subset=["Morph", "Category", "LFC"]).copy()

df["Morph"] = df["Morph"].astype(str).str.strip()
df["Category"] = df["Category"].astype(str).str.strip()
df["Gene Symbol"] = df["Gene Symbol"].astype(str).str.strip()

df["absLFC"] = df["LFC"].abs()

# -----------------------------
# 2) Aggregate to Category level only (Morph -> Category)
# -----------------------------
agg = (
    df.groupby(["Morph", "Category"], as_index=False)
      .agg(
          absLFC=("absLFC", "sum"),
          n_genes=("Gene Symbol", "nunique"),
      )
)

top5 = (
    df.sort_values(["Morph", "Category", "absLFC"], ascending=[True, True, False])
      .groupby(["Morph", "Category"])["Gene Symbol"]
      .apply(lambda s: ", ".join(s.head(5).astype(str)))
      .reset_index(name="top_genes")
)

agg = agg.merge(top5, on=["Morph", "Category"], how="left")
agg["Morph"] = agg["Morph"].astype(str).str.strip()
agg["Category"] = agg["Category"].astype(str).str.strip()
agg["top_genes"] = agg["top_genes"].fillna("")

# -----------------------------
# 3) Colors
# -----------------------------
all_categories = sorted(agg["Category"].unique())
palette = px.colors.qualitative.Plotly
category_color_map = {c: palette[i % len(palette)] for i, c in enumerate(all_categories)}

# -----------------------------
# 4) Build combined subplot figure (one per morph)
# -----------------------------
morphs = sorted(agg["Morph"].unique())
fig = make_subplots(
    rows=1,
    cols=len(morphs),
    specs=[[{"type": "domain"}] * len(morphs)],
    subplot_titles=[""] * len(morphs),
)

for col_i, morph in enumerate(morphs, start=1):
    agg_m = agg[agg["Morph"] == morph].copy()
    agg_m["Root"] = morph  # Option B: center circle

    fig_m = px.sunburst(
        agg_m,
        path=["Root", "Category"],
        values="absLFC",
        custom_data=["absLFC", "n_genes", "top_genes"],
    )

    trace = fig_m.data[0]
    ids = list(trace.ids)

    marker_colors = []
    for node_id in ids:
        parts = [p.strip() for p in str(node_id).split("/")]
        if len(parts) == 1:
            marker_colors.append(ROOT_COLOR)
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
# 5) Layout
# -----------------------------
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
fig.write_html(OUT_HTML)
print(f"Saved HTML: {OUT_HTML}")

# -----------------------------
# 6) Static exports (SVG/PNG/TIFF/PDF)
# -----------------------------
# Kaleido must be active and Chrome must exist
if pio.kaleido.scope is None:
    raise RuntimeError("Kaleido is not active (pio.kaleido.scope is None).")

ensure_chrome_for_kaleido()

def with_bg(fig_in, transparent: bool):
    f = copy.deepcopy(fig_in)
    if transparent:
        f.update_layout(paper_bgcolor="rgba(0,0,0,0)", plot_bgcolor="rgba(0,0,0,0)")
    else:
        f.update_layout(paper_bgcolor="white", plot_bgcolor="white")
    return f

# SVG: transparent + white
with_bg(fig, True).write_image(OUT_SVG_TRANSPARENT, format="svg", width=W, height=H, scale=1)
print(f"Saved SVG (transparent): {OUT_SVG_TRANSPARENT}")

with_bg(fig, False).write_image(OUT_SVG_WHITE, format="svg", width=W, height=H, scale=1)
print(f"Saved SVG (white): {OUT_SVG_WHITE}")

# PNG: embed 600 dpi metadata
tmp_png = "_tmp_export.png"
with_bg(fig, False).write_image(tmp_png, format="png", width=W, height=H, scale=1)
im = Image.open(tmp_png)
im.save(OUT_PNG_600, dpi=(DPI, DPI))
im.close()
os.remove(tmp_png)
print(f"Saved PNG ({DPI} dpi metadata): {OUT_PNG_600}")

# TIFF: embed 600 dpi metadata
tmp_png = "_tmp_export.png"
with_bg(fig, False).write_image(tmp_png, format="png", width=W, height=H, scale=1)
im = Image.open(tmp_png).convert("RGB")
im.save(OUT_TIFF_600, dpi=(DPI, DPI), compression="tiff_lzw")
im.close()
os.remove(tmp_png)
print(f"Saved TIFF ({DPI} dpi metadata): {OUT_TIFF_600}")

# PDF: vector
with_bg(fig, False).write_image(OUT_PDF, format="pdf", width=W, height=H, scale=1)
print(f"Saved PDF: {OUT_PDF}")
