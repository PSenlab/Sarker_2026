#!/usr/bin/env python
# coding: utf-8

# In[9]:


#!/usr/bin/env python
# coding: utf-8

# ============================================================
# GREAT REACTOME (REGION-BASED) DOT MATRIX PLOT
# RELAXED VERSION (corCutOff = 0.25)
#
# Filter: p_adjust <= 0.05, fold_enrichment >= 1, 
#         observed_region_hits >= 5
#
# Plot:
#   - Rows = Pathways
#   - Columns = Male | Female
#   - Background color = fold_enrichment
#   - Dot color = Stability class
#   - Dot size = observed_gene_hits
#   - Sorted by = observed_gene_hits
# ============================================================

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Rectangle
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
import warnings
warnings.filterwarnings("ignore")

# ============================================================
# CONFIGURATION
# ============================================================

# RELAXED VERSION
great_outdir = "/data/sarkern2/multiome_liver/Seurat/epigenome/downstream_stability/relaxed/peak2gene_stability_analysis/GREAT_analysis"

STABILITY_CLASSES = [
    "Stable_Active",
    "Stable_Repressive",
    "Monotonic_A_to_R",
    "Monotonic_R_to_A",
    "Non_Monotonic"
]

SEXES = ["male", "female"]

# Stability class colors
STABILITY_COLS = {
    "Stable_Active": "#FF6B6B",
    "Stable_Repressive": "#4E79A7",
    "Monotonic_A_to_R": "darkred",
    "Monotonic_R_to_A": "darkblue",
    "Non_Monotonic": "#708238"
}

STABILITY_LABELS = {
    "Stable_Active": "Stable Active",
    "Stable_Repressive": "Stable Repressive",
    "Monotonic_A_to_R": "A → R",
    "Monotonic_R_to_A": "R → A",
    "Non_Monotonic": "Non-Monotonic"
}

# ============================================================
# LOAD AND COMBINE ALL REACTOME RESULTS
# ============================================================

print("\n" + "="*60)
print("LOADING GREAT REACTOME RESULTS (RELAXED)")
print("="*60)

# Find all Reactome files
reactome_files = [f for f in os.listdir(great_outdir) 
                  if f.startswith("GREAT_") and f.endswith("_Reactome.tsv")]

print(f"\n>>> Found {len(reactome_files)} Reactome files")

# Load and combine
all_dfs = []
for f in reactome_files:
    filepath = os.path.join(great_outdir, f)
    df = pd.read_csv(filepath, sep="\t")
    all_dfs.append(df)
    print(f"  ✓ Loaded: {f} ({len(df)} terms)")

dt = pd.concat(all_dfs, ignore_index=True)
print(f"\n>>> Combined: {len(dt)} total terms")

# ============================================================
# PREPARE DATA
# ============================================================

print("\n>>> Preparing data...")

# Extract Sex and Stability from Analysis column
if 'Sex' not in dt.columns:
    dt['Sex'] = dt['Analysis'].str.extract(r'^(male|female)_')[0]
if 'Stability' not in dt.columns:
    dt['Stability'] = dt['Analysis'].str.replace(r'^(male|female)_', '', regex=True)

# Clean pathway names
dt['Pathway_Clean'] = dt['Pathway_Name'].str.replace(r'^Mus musculus: ', '', regex=True)

# Ensure numeric
for col in ['observed_region_hits', 'fold_enrichment', 'p_adjust', 
            'observed_gene_hits', 'gene_set_size']:
    if col in dt.columns:
        dt[col] = pd.to_numeric(dt[col], errors='coerce')

print(f"  Sex values: {dt['Sex'].unique().tolist()}")
print(f"  Stability values: {dt['Stability'].unique().tolist()}")

# ============================================================
# FILTER: REGION-BASED (BINOMIAL) TEST
# ============================================================

print("\n>>> Filtering (region-based)...")

dt_filtered = dt[
    (dt['p_adjust'] <= 0.05) &
    (dt['fold_enrichment'] >= 1) &
    (dt['observed_region_hits'] >= 5) 
].copy()

print(f"  Significant pathways: {len(dt_filtered)} (from {len(dt)})")

# ============================================================
# SELECT TOP TERMS PER ANALYSIS (SORTED BY observed_gene_hits)
# ============================================================

print("\n>>> Selecting top terms per Analysis...")

topN = 10

selected_list = []
for analysis in dt_filtered['Analysis'].unique():
    subset = dt_filtered[dt_filtered['Analysis'] == analysis].copy()
    subset = subset.sort_values(['observed_gene_hits', 'p_adjust'], 
                                 ascending=[False, True])
    selected_list.append(subset.head(topN))

plot_df = pd.concat(selected_list, ignore_index=True)

print(f"  Selected: {len(plot_df)} pathways (top {topN} per analysis)")

# Truncate pathway names
plot_df['Term'] = plot_df['Pathway_Clean'].apply(
    lambda x: x[:42] + '...' if len(x) > 45 else x
)

# ============================================================
# SUMMARY
# ============================================================

print("\n>>> Summary per Analysis:")
summary = plot_df.groupby(['Sex', 'Stability']).agg({
    'Term': 'count',
    'fold_enrichment': 'max',
    'observed_gene_hits': 'max'
}).reset_index()
summary.columns = ['Sex', 'Stability', 'N_pathways', 'Max_FE', 'Max_GeneHits']
print(summary.to_string(index=False))

# ============================================================
# DOT MATRIX PLOT FUNCTION
# ============================================================

def plot_dot_matrix_great(plot_df, stability_cols, outfile_prefix, 
                           title="GREAT Reactome Pathways", 
                           figsize_width=8):
    """
    Dot matrix plot for GREAT region-based results.
    
    - Rows = Pathways
    - Columns = Male | Female
    - Background color = fold_enrichment (max per cell)
    - Dot color = Stability class (multiple dots if multiple classes)
    - Dot size = observed_gene_hits (normalized scaling)
    - Sorted by = observed_gene_hits
    """
    
    df = plot_df.copy()
    
    # Get all unique terms, sorted by max observed_gene_hits (descending)
    term_order = (df.groupby('Term')['observed_gene_hits']
                   .max()
                   .sort_values(ascending=False)
                   .index.tolist())
    
    term_to_y = {term: i for i, term in enumerate(term_order)}
    sex_to_x = {'male': 0, 'female': 1}
    
    # Figure sizing
    n_terms = len(term_order)
    fig_height = max(4, 0.20 * n_terms)
    fig, ax = plt.subplots(figsize=(figsize_width, fig_height))
    
    # Color normalization for fold_enrichment background
    fe_max = df['fold_enrichment'].max()
    norm = Normalize(vmin=0, vmax=fe_max)
    cmap = plt.cm.Oranges
    
    # Calculate max fold_enrichment per (Term, Sex) for background
    max_fe_per_cell = df.groupby(['Term', 'Sex'])['fold_enrichment'].max().to_dict()
    
    # Draw background rectangles based on MAX fold_enrichment in each cell
    for term in term_order:
        for sex in ['male', 'female']:
            x = sex_to_x[sex]
            y = term_to_y[term]
            
            fe_val = max_fe_per_cell.get((term, sex), 0)
            
            if fe_val > 0:
                facecolor = cmap(norm(fe_val))
            else:
                facecolor = '#f5f5f5'
            
            rect = Rectangle(
                (x - 0.48, y - 0.48), 0.96, 0.96,
                facecolor=facecolor,
                edgecolor='white',
                linewidth=1.5,
                zorder=1
            )
            ax.add_patch(rect)
    
    # Count how many stability classes per (Term, Sex)
    stab_counts = df.groupby(['Term', 'Sex']).size().to_dict()
    
    # Define offsets for multiple dots in same cell
    def get_dot_offsets(n_dots):
        """Return x,y offsets for n_dots within a cell"""
        if n_dots == 1:
            return [(0, 0)]
        elif n_dots == 2:
            return [(-0.15, 0), (0.15, 0)]
        elif n_dots == 3:
            return [(-0.2, 0), (0, 0), (0.2, 0)]
        elif n_dots == 4:
            return [(-0.15, 0.12), (0.15, 0.12), (-0.15, -0.12), (0.15, -0.12)]
        else:  # 5+
            return [(-0.2, 0.12), (0, 0.12), (0.2, 0.12), (-0.15, -0.12), (0.15, -0.12)]
    
    # Size parameters
    min_size = 30
    max_single_size = 400
    
    # Get min/max for normalization
    min_hits_global = df['observed_gene_hits'].min()
    max_hits_global = df['observed_gene_hits'].max()
    
    # Size mapping function
    def get_dot_size(hits, n_dots=1):
        # Normalize to 0-1 range, then scale to size range
        if max_hits_global == min_hits_global:
            base_size = (min_size + max_single_size) / 2
        else:
            normalized = (hits - min_hits_global) / (max_hits_global - min_hits_global)
            base_size = min_size + normalized * (max_single_size - min_size)
        
        # Reduce size when multiple dots in cell
        if n_dots == 1:
            return base_size
        elif n_dots == 2:
            return base_size * 0.7
        else:
            return base_size * 0.5
    
    # Draw dots - handle multiple stability classes per cell
    for (term, sex), group in df.groupby(['Term', 'Sex']):
        x_base = sex_to_x[sex]
        y_base = term_to_y[term]
        
        n_dots = len(group)
        offsets = get_dot_offsets(n_dots)
        
        # Sort by observed_gene_hits so largest dot is drawn last (on top)
        group_sorted = group.sort_values('observed_gene_hits', ascending=True)
        
        for i, (_, row) in enumerate(group_sorted.iterrows()):
            stab = row['Stability']
            gene_hits = row['observed_gene_hits']
            
            # Get offset
            if i < len(offsets):
                x_off, y_off = offsets[i]
            else:
                x_off, y_off = 0, 0
            
            # Calculate size using normalized scaling
            size = get_dot_size(gene_hits, n_dots)
            
            ax.scatter(
                x_base + x_off, 
                y_base + y_off, 
                s=size,
                c=stability_cols.get(stab, 'gray'),
                edgecolor='black',
                linewidth=0.8,
                alpha=0.95,
                zorder=5,
                clip_on=False
            )
    
    # Axis setup
    ax.set_xlim(-0.6, 1.6)
    ax.set_ylim(-0.6, n_terms - 0.4)
    
    # X-axis (Sex labels at top)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(['MALE', 'FEMALE'], fontsize=11, fontweight='bold')
    ax.xaxis.set_ticks_position('top')
    ax.xaxis.set_label_position('top')
    
    # Y-axis (Pathway terms)
    ax.set_yticks(range(n_terms))
    ax.set_yticklabels(term_order, fontsize=8)
    
    # Title
    ax.set_title(f'{title}\n', fontsize=10, fontweight='bold', pad=5)
    
    # Invert y-axis (highest gene_hits at top)
    ax.invert_yaxis()
    
    # Remove spines
    for spine in ax.spines.values():
        spine.set_visible(False)
    
    ax.tick_params(axis='both', length=0)
    
    # ===== LEGENDS =====
    
    plt.subplots_adjust(right=0.70)
    
    # 1. Colorbar for fold_enrichment (background)
    sm = ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar_ax = fig.add_axes([0.73, 0.60, 0.02, 0.25])
    cbar = plt.colorbar(sm, cax=cbar_ax)
    cbar.set_label('Fold Enrichment', fontsize=9, fontweight='bold')
    cbar.ax.tick_params(labelsize=8)
    cbar.outline.set_linewidth(0.5)
    
    # 2. Stability class legend (dot color)
    unique_stab = df['Stability'].unique()
    color_patches = [
        mpatches.Patch(
            color=stability_cols.get(s, 'gray'),
            label=STABILITY_LABELS.get(s, s),
            edgecolor='black',
            linewidth=0.5
        )
        for s in STABILITY_COLS.keys() if s in unique_stab
    ]
    
    legend1 = fig.legend(
        handles=color_patches,
        title='Stability',
        loc='center left',
        bbox_to_anchor=(0.72, 0.38),
        framealpha=0.95,
        fontsize=8,
        title_fontsize=9,
        edgecolor='gray'
    )
    
    # 3. Size legend (observed_gene_hits) - Simple scatter legend
    min_hits = df['observed_gene_hits'].min()
    max_hits = df['observed_gene_hits'].max()
    
    # Choose 3 representative legend values
    if max_hits - min_hits < 10:
        legend_vals = [int(min_hits), int(max_hits)]
    else:
        mid_hits = (min_hits + max_hits) / 2
        legend_vals = [int(min_hits), int(mid_hits), int(max_hits)]
    legend_vals = sorted(set(legend_vals))
    
    # Create size mapping function (same as used in plotting)
    def get_legend_size(hits):
        if max_hits == min_hits:
            return (min_size + max_single_size) / 2
        normalized = (hits - min_hits) / (max_hits - min_hits)
        return min_size + normalized * (max_single_size - min_size)
    
    size_elements = [
        plt.scatter(
            [], [], 
            s=get_legend_size(v),
            c='white',
            edgecolor='black',
            linewidth=0.8,
            label=f'{v}'
        )
        for v in legend_vals
    ]
    
    legend2 = fig.legend(
        handles=size_elements,
        title='Gene Hits',
        loc='center left',
        bbox_to_anchor=(0.72, 0.12),
        framealpha=0.95,
        fontsize=8,
        title_fontsize=9,
        edgecolor='gray'
    )
    
    # Save
    pdf_file = os.path.join(great_outdir, f"{outfile_prefix}_dotmatrix.pdf")
    png_file = os.path.join(great_outdir, f"{outfile_prefix}_dotmatrix.png")
    
    plt.savefig(pdf_file, dpi=300, bbox_inches='tight', facecolor='white')
    plt.savefig(png_file, dpi=300, bbox_inches='tight', facecolor='white')
    plt.show()
    
    print(f"\n  ✓ Saved: {pdf_file}")
    print(f"  ✓ Saved: {png_file}")
    print(f"  ✓ Terms: {n_terms}")
    print(f"  ✓ Fold Enrichment range: 0 - {fe_max:.1f}")
    print(f"  ✓ Gene Hits range: {min_hits} - {max_hits}")
    
    # Report multi-dot cells
    multi_dot_cells = {k: v for k, v in stab_counts.items() if v > 1}
    if multi_dot_cells:
        print(f"\n  ✓ Cells with multiple dots: {len(multi_dot_cells)}")
        for (term, sex), count in sorted(multi_dot_cells.items(), key=lambda x: -x[1]):
            print(f"      {sex}: {term} ({count} stability classes)")

# ============================================================
# GENERATE FINAL PLOT
# ============================================================

print("\n" + "="*60)
print("GENERATING GREAT REACTOME DOT MATRIX PLOT (RELAXED)")
print("="*60)

plot_dot_matrix_great(
    plot_df=plot_df,
    stability_cols=STABILITY_COLS,
    outfile_prefix="GREAT_Reactome_region_relaxed",
    title="GREAT Reactome (Region-Based) - RELAXED\np_adjust < 0.05, FE ≥ 1, Region Hits ≥ 5",
    figsize_width=3
)

print("\n✓ Final plot generated!")

# ============================================================
# SAVE FILTERED DATA
# ============================================================

print("\n>>> Saving filtered data...")

plot_df.to_csv(os.path.join(great_outdir, "GREAT_Reactome_region_plot_data_relaxed.tsv"), 
               sep="\t", index=False)
print(f"  ✓ Saved: GREAT_Reactome_region_plot_data_relaxed.tsv")

print("\n" + "="*60)
print("✓ DONE!")
print("="*60)


# In[ ]:




