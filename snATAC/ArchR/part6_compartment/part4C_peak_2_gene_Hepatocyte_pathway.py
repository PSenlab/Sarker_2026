#!/usr/bin/env python
# coding: utf-8

# In[8]:


#!/usr/bin/env python
# coding: utf-8

# In[58]:


# ============================================================
# PATHWAY ENRICHMENT WITH PROPER BACKGROUND FROM ADATA
# ============================================================

get_ipython().run_line_magic('load_ext', 'autoreload')
get_ipython().run_line_magic('autoreload', '2')

import pandas as pd
import numpy as np
import os
import glob
import gseapy as gp
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
import mygene

adata = sc.read_h5ad("/data/sarkern2/multiome_liver/final_object/final_rna_wnn.h5ad")
print(f"  Full adata: {adata.shape[0]} cells x {adata.shape[1]} genes")


# In[59]:


print("\n>>> Unique celltype labels:")
print(sorted(adata.obs["celltype"].unique()))


# In[60]:


# Subset to Hepatocytes
adata = adata[adata.obs["celltype"] == "Hepatocyte"].copy()
print(f"  Hepatocyte adata: {adata.shape[0]} cells x {adata.shape[1]} genes")

# Get expressed genes (background)
expr = adata.X
gene_names = adata.var_names
avg_expr = np.asarray(expr.mean(axis=0)).flatten()

# Set threshold for expressed genes
threshold = 0.05
expressed_genes = gene_names[avg_expr > threshold].tolist()

print(f"\n✓ Background genes (expressed in Hepatocytes): {len(expressed_genes)}")
print(f"  Threshold: mean expression > {threshold}")
print(f"  Sample genes: {expressed_genes[:10]}")


# In[67]:


# ============================================================
# LOAD GENES FROM BOTH DEFAULT AND RELAXED
# ============================================================
import os

print("\n>>> Loading all stability class genes from DEFAULT and RELAXED...")

# Base directory
base_outdir = "/data/sarkern2/multiome_liver/Seurat/epigenome/downstream_stability"

# Two P2G cutoff versions
P2G_VERSIONS = {
    "default": os.path.join(base_outdir, "default", "peak2gene_stability_analysis"),
    "relaxed": os.path.join(base_outdir, "relaxed", "peak2gene_stability_analysis")
}

STABILITY_CLASSES = ["Stable_Active", "Stable_Repressive", "Monotonic_A_to_R", "Monotonic_R_to_A", "Non_Monotonic"]
SEXES = ["male", "female"]

# Store genes: genes_by_version[version][sex][stab_class]
genes_by_version = {}

for version, p2g_outdir in P2G_VERSIONS.items():
    print(f"\n{'='*60}")
    print(f"Loading: {version.upper()} (corCutOff = {'0.45' if version == 'default' else '0.25'})")
    print(f"{'='*60}")
    print(f"  Path: {p2g_outdir}")
    
    genes_by_version[version] = {}
    
    for sex in SEXES:
        genes_by_version[version][sex] = {}
        for stab_class in STABILITY_CLASSES:
            gene_file = os.path.join(p2g_outdir, f"genes_{sex}_{stab_class}.txt")
            if os.path.exists(gene_file):
                with open(gene_file, 'r') as f:
                    genes_by_version[version][sex][stab_class] = [l.strip() for l in f if l.strip()]
                print(f"  {sex:6} | {stab_class:20}: {len(genes_by_version[version][sex][stab_class]):5} genes")
            else:
                genes_by_version[version][sex][stab_class] = []
                print(f"  {sex:6} | {stab_class:20}: FILE NOT FOUND")

print("\n✓ Genes loaded into: genes_by_version[version][sex][stab_class]")
print("  Example: genes_by_version['default']['male']['Stable_Active']")
print("  Example: genes_by_version['relaxed']['female']['Monotonic_A_to_R']")


# In[ ]:


# ============================================================
# SELECT WHICH VERSION TO USE FOR ANALYSIS
# ============================================================

# Choose version: "default" (corCutOff=0.45) or "relaxed" (corCutOff=0.25)
SELECTED_VERSION = "relaxed"  # <-- CHANGE THIS TO SWITCH BETWEEN VERSIONS

print(f"\n>>> Using {SELECTED_VERSION.upper()} version for downstream analysis")

# Set genes_by_sex_class to selected version
genes_by_sex_class = genes_by_version[SELECTED_VERSION]

# Summary
print(f"\nGene counts for {SELECTED_VERSION.upper()}:")
for sex in SEXES:
    print(f"\n  {sex.upper()}:")
    for stab_class in STABILITY_CLASSES:
        n_genes = len(genes_by_sex_class[sex][stab_class])
        print(f"    {stab_class:20}: {n_genes:5} genes")


# In[ ]:


# ============================================================
# STEP 2: DEFINE CONVERSION FUNCTION
# ============================================================

def convert_mouse_to_human(mouse_genes):
    """Convert mouse gene symbols to human ortholog symbols"""
    
    mg = mygene.MyGeneInfo()
    
    # Ensure unique, clean input
    mouse_genes = list(set([g for g in mouse_genes if pd.notna(g) and g]))
    
    if len(mouse_genes) == 0:
        return []
    
    # Query mouse genes for homologs
    gene_conversion = mg.querymany(
        mouse_genes,
        scopes="symbol",
        fields="homologene",
        species="mouse",
        as_dataframe=False,
        verbose=False
    )
    
    # Collect human Entrez IDs
    human_entrez_ids = set()
    for entry in gene_conversion:
        if isinstance(entry, dict) and "homologene" in entry:
            genes = entry["homologene"].get("genes", [])
            for tax_id, entrez_id in genes:
                if tax_id == 9606:  # Human
                    human_entrez_ids.add(entrez_id)
    
    if len(human_entrez_ids) == 0:
        return []
    
    # Convert Entrez IDs → human symbols
    human_conversion = mg.querymany(
        list(human_entrez_ids),
        scopes="entrezgene",
        fields="symbol",
        species="human",
        as_dataframe=False,
        verbose=False
    )
    
    human_symbols = sorted({
        entry["symbol"].upper()
        for entry in human_conversion
        if isinstance(entry, dict) and "symbol" in entry
    })
    
    return human_symbols

# ============================================================
# STEP 3: CONVERT ALL STABILITY CLASSES
# ============================================================
print("\n>>> Step 2: Converting mouse genes to human orthologs...")

human_genes_by_sex_class = {}

for sex in SEXES:
    print(f"\n{sex.upper()}:")
    human_genes_by_sex_class[sex] = {}
    
    for stab_class in STABILITY_CLASSES:
        mouse_genes = genes_by_sex_class[sex][stab_class]
        
        if len(mouse_genes) == 0:
            human_genes_by_sex_class[sex][stab_class] = []
            print(f"    {stab_class:20}: 0 mouse → 0 human (skipped)")
            continue
        
        human_genes = convert_mouse_to_human(mouse_genes)
        human_genes_by_sex_class[sex][stab_class] = human_genes
        
        conversion_rate = 100 * len(human_genes) / len(mouse_genes) if len(mouse_genes) > 0 else 0
        print(f"    {stab_class:20}: {len(mouse_genes):5} mouse → {len(human_genes):5} human ({conversion_rate:.1f}%)")

# ============================================================
# STEP 4: SUMMARY
# ============================================================
print("\n" + "="*60)
print("CONVERSION SUMMARY")
print("="*60)

for sex in SEXES:
    print(f"\n{sex.upper()}:")
    for stab_class in STABILITY_CLASSES:
        n_mouse = len(genes_by_sex_class[sex][stab_class])
        n_human = len(human_genes_by_sex_class[sex][stab_class])
        print(f"    {stab_class:20}: {n_mouse:5} mouse → {n_human:5} human")

# ============================================================
# ACCESS EXAMPLES
# ============================================================
# human_genes_by_sex_class["male"]["Monotonic_A_to_R"]
# human_genes_by_sex_class["female"]["Non_Monotonic"]
# human_genes_by_sex_class["male"]["Stable_Active"]

print("\n✓ Human genes stored in: human_genes_by_sex_class[sex][stability_class]")


# In[ ]:


# ============================================================
# ENRICHR ANALYSIS FUNCTION (ROBUST, NO FILTERING)
# ============================================================
def run_enrichr_analysis(gene_list, background, gene_set, title, outfile_prefix):
    """Run Enrichr with explicit background, no post-filtering"""
    # Standardize gene symbols
    gene_list = sorted(set(g.upper() for g in gene_list if isinstance(g, str)))
    background = sorted(set(g.upper() for g in background if isinstance(g, str)))
    # Ensure gene list is a subset of background
    gene_list = sorted(set(gene_list) & set(background))
    if len(gene_list) < 3:
        print(f"  ⚠️ Too few genes ({len(gene_list)}) for {title}")
        return pd.DataFrame()
    try:
        enr = gp.enrichr(
            gene_list=gene_list,
            gene_sets=gene_set,
            organism="Human",
            background=background,
            outdir=None,
            no_plot=True
        )
        if enr is None or not hasattr(enr, "results"):
            print(f"  ⚠️ Enrichr returned no results for {title}")
            return pd.DataFrame()
        df = enr.results.copy()
        if df is None or df.empty:
            print(f"  ⚠️ No enrichment terms for {title}")
            return pd.DataFrame()
        # Recount overlap with input genes
        input_genes = set(gene_list)
        def recount_overlap(gene_str):
            if not isinstance(gene_str, str):
                return 0
            genes = {g.strip().upper() for g in gene_str.split(";") if g.strip()}
            return len(genes & input_genes)
        df["Gene_Count"] = df["Genes"].apply(recount_overlap)
        df["-log10(pval)"] = -np.log10(df["P-value"].clip(lower=1e-300))
        df["-log10(padj)"] = -np.log10(df["Adjusted P-value"].clip(lower=1e-300))
        # Minimal provenance metadata
        df["GeneSet"] = gene_set
        df["N_input_genes"] = len(gene_list)
        df["N_background_genes"] = len(background)
        # Save full results
        out_csv = f"{outfile_prefix}_results.csv"
        df.to_csv(out_csv, index=False)
        print(f"  ✓ {title}: {len(df)} enriched terms")
        print(f"    • Input genes: {len(gene_list)}")
        print(f"    • Background genes: {len(background)}")
        return df
    except Exception as e:
        print(f"  ❌ Error in {title}: {e}")
        return pd.DataFrame()

# ============================================================
# RUN PATHWAY ANALYSIS (PER SEX × STABILITY CLASS)
# ============================================================
print("\n" + "="*60)
print(f"RUNNING PATHWAY ANALYSIS (PER SEX × STABILITY CLASS) - {SELECTED_VERSION.upper()}")
print("="*60)

# Gene set libraries to query
gene_sets = [
    "Reactome_Pathways_2024"
]

# Background = all P2G human genes (union across all stability classes)
background_human = set()
for sex in SEXES:
    for stab_class in STABILITY_CLASSES:
        background_human.update(human_genes_by_sex_class[sex][stab_class])
background_human = sorted(background_human)

print(f"\nBackground genes (human): {len(background_human)}")

# Store results
all_results = {}

for sex in SEXES:
    all_results[sex] = {}
    
    for stab_class in STABILITY_CLASSES:
        
        # Get human genes for this Sex × Stability
        p2g_genes_human = human_genes_by_sex_class[sex][stab_class]
        
        print(f"\n{'='*50}")
        print(f"{sex.upper()} - {stab_class}")
        print(f"{'='*50}")
        print(f"Foreground genes (human): {len(p2g_genes_human)}")
        print(f"Background genes (human): {len(background_human)}")
        
        if len(p2g_genes_human) < 3:
            print(f"  ⚠️ Skipping (too few genes)")
            all_results[sex][stab_class] = {}
            continue
        
        all_results[sex][stab_class] = {}
        
        for gs in gene_sets:
            gs_short = gs.split("_")[0]
            
            result = run_enrichr_analysis(
                gene_list=p2g_genes_human,
                background=background_human,
                gene_set=gs,
                title=f"{sex}_{stab_class} - {gs_short}",
                outfile_prefix=f"{SELECTED_VERSION}_{sex}_{stab_class}_{gs}"
            )
            
            if result is not None and not result.empty:
                all_results[sex][stab_class][gs] = result

print("\n" + "="*60)
print("✓ Pathway analysis complete")
print("="*60)


# In[ ]:


# ============================================================
# COMBINED DOT PLOT - ALL MODULES
# ============================================================

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib.patches as mpatches

print("="*60)
print("CREATING COMBINED DOT PLOT")
print("="*60)



# In[ ]:


# ============================================================
# FILTER AND PREPARE DATA (PER SEX × STABILITY CLASS)
# Using P-value (not Adjusted P-value) for consistency
# ============================================================
import numpy as np
import pandas as pd

print("\n" + "="*60)
print("FILTERING AND PREPARING DATA FOR PLOTTING")
print("="*60)

# Store filtered results separately (preserve original all_results)
filtered_results = {}
combined_data = []

for sex in SEXES:
    filtered_results[sex] = {}
    
    for stab_class in STABILITY_CLASSES:
        filtered_results[sex][stab_class] = {}
        
        if stab_class not in all_results[sex]:
            continue
        
        for gs, df in all_results[sex][stab_class].items():
            if df is None or df.empty:
                continue
            
            # Filter criteria - using P-value
            df_filtered = df[
                (df['P-value'] < 0.05) &
                (df['Gene_Count'] > 3) &
                (np.isfinite(df['Odds Ratio']))
            ].copy()
            
            # Store filtered version
            filtered_results[sex][stab_class][gs] = df_filtered
            
            if df_filtered.empty:
                print(f"  {sex} | {stab_class} | {gs}: No terms pass filter")
                continue
            
            # Sort by Odds Ratio and take top 10
            df_top = (
                df_filtered
                .sort_values('Odds Ratio', ascending=False)
                .head(10)
            )
            
            # Add metadata columns
            df_top = df_top.copy()
            df_top['Sex'] = sex
            df_top['Stability'] = stab_class
            df_top['GeneSet'] = gs
            
            combined_data.append(df_top)
            print(f"  {sex:6} | {stab_class:20} | {gs}: Selected {len(df_top)} terms")

# Combine all
if combined_data:
    plot_df = pd.concat(combined_data, ignore_index=True)
    print(f"\n✓ Total terms for plotting: {len(plot_df)}")
else:
    print("❌ No data to plot!")
    plot_df = pd.DataFrame()

# ============================================================
# SUMMARY BY SEX × STABILITY
# ============================================================
if not plot_df.empty:
    print("\n>>> Terms per Sex × Stability:")
    summary = plot_df.groupby(['Sex', 'Stability']).size().reset_index(name='N_terms')
    print(summary.to_string(index=False))


# ============================================================
# FINALIZED DOT MATRIX PLOT - USING P-VALUE (NOT ADJUSTED)
# ============================================================
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Rectangle
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable

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

def plot_dot_matrix_final(plot_df, stability_cols, outfile_prefix, 
                          title="Reactome Pathways", 
                          figsize_width=8,
                          size_scale=30,
                          min_size=20):
    """
    Finalized dot matrix plot:
    - Rows = Pathways
    - Columns = Male | Female
    - Background color = Odds Ratio
    - Dot color = Stability class
    - Dot size = -log10(P-value)  # Changed from Adjusted P-value
    """
    
    df = plot_df.copy()
    
    # Use P-value (not Adjusted P-value) for consistency with filter
    df['neglog10_pval'] = -np.log10(df['P-value'].clip(lower=1e-300))
    
    # Get all unique terms, sorted by max Odds Ratio (descending)
    term_order = (df.groupby('Term')['Odds Ratio']
                   .max()
                   .sort_values(ascending=False)
                   .index.tolist())
    
    term_to_y = {term: i for i, term in enumerate(term_order)}
    sex_to_x = {'male': 0, 'female': 1}
    
    # Figure sizing
    n_terms = len(term_order)
    fig_height = max(2, 0.20 * n_terms)
    fig, ax = plt.subplots(figsize=(figsize_width, fig_height))
    
    # Color normalization for Odds Ratio background
    or_max = df['Odds Ratio'].max()
    norm = Normalize(vmin=0, vmax=or_max)
    cmap = plt.cm.Blues
    
    # Draw background rectangles for each cell
    for _, row in df.iterrows():
        x = sex_to_x.get(row['Sex'], 0)
        y = term_to_y.get(row['Term'], 0)
        or_val = row['Odds Ratio']
        
        rect = Rectangle(
            (x - 0.48, y - 0.48), 0.96, 0.96,
            facecolor=cmap(norm(or_val)),
            edgecolor='white',
            linewidth=1.5,
            zorder=1
        )
        ax.add_patch(rect)
    
    # Draw empty cells (gray) for missing combinations
    all_sexes = ['male', 'female']
    existing_combos = set(zip(df['Term'], df['Sex']))
    
    for term in term_order:
        for sex in all_sexes:
            if (term, sex) not in existing_combos:
                x = sex_to_x[sex]
                y = term_to_y[term]
                rect = Rectangle(
                    (x - 0.48, y - 0.48), 0.96, 0.96,
                    facecolor='#f5f5f5',
                    edgecolor='white',
                    linewidth=1.5,
                    zorder=1
                )
                ax.add_patch(rect)
    
    # Draw dots on top - using P-value
    for _, row in df.iterrows():
        x = sex_to_x.get(row['Sex'], 0)
        y = term_to_y.get(row['Term'], 0)
        stab = row['Stability']
        pval = row['neglog10_pval']  # Changed from neglog10_adjP
        
        # Cap size to fit within cell
        size = min(pval * size_scale + min_size, 400)
        
        ax.scatter(
            x, y, 
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
    ax.set_title(f'{title}\n', fontsize=9, fontweight='bold', pad=5)
    
    # Invert y-axis (highest OR at top)
    ax.invert_yaxis()
    
    # Remove spines
    for spine in ax.spines.values():
        spine.set_visible(False)
    
    ax.tick_params(axis='both', length=0)
    
    # ===== LEGENDS - REPOSITIONED =====
    
    # Adjust plot to make room for legends
    plt.subplots_adjust(right=0.72)
    
    # 1. Colorbar for Odds Ratio (background)
    sm = ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar_ax = fig.add_axes([0.75, 0.60, 0.02, 0.25])
    cbar = plt.colorbar(sm, cax=cbar_ax)
    cbar.set_label('Odds Ratio', fontsize=9, fontweight='bold')
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
        bbox_to_anchor=(0.74, 0.38),
        framealpha=0.95,
        fontsize=8,
        title_fontsize=9,
        edgecolor='gray'
    )
    
    # 3. Size legend (-log10 P-value)
    min_pval = df['neglog10_pval'].min()
    max_pval = df['neglog10_pval'].max()
    
    # Create nice legend values
    if max_pval - min_pval < 2:
        legend_vals = [round(min_pval, 1), round(max_pval, 1)]
    else:
        mid_pval = (min_pval + max_pval) / 2
        legend_vals = [round(min_pval, 1), round(mid_pval, 1), round(max_pval, 1)]
    legend_vals = sorted(set(legend_vals))
    
    size_elements = [
        plt.scatter(
            [], [], 
            s=min(v * size_scale + min_size, 400),
            c='white',
            edgecolor='black',
            linewidth=0.8,
            label=f'{v:.1f}'
        )
        for v in legend_vals
    ]
    
    legend2 = fig.legend(
        handles=size_elements,
        title='-log₁₀(P)',  # Changed label to reflect P-value
        loc='center left',
        bbox_to_anchor=(0.74, 0.12),
        framealpha=0.95,
        fontsize=8,
        title_fontsize=10,
        edgecolor='gray'
    )
    
    # Save
    pdf_file = f"{outfile_prefix}_dotmatrix.pdf"
    png_file = f"{outfile_prefix}_dotmatrix.png"
    
    plt.savefig(pdf_file, dpi=300, bbox_inches='tight', facecolor='white')
    plt.savefig(png_file, dpi=300, bbox_inches='tight', facecolor='white')
    plt.show()
    
    print(f"\n  ✓ Saved: {pdf_file}")
    print(f"  ✓ Saved: {png_file}")
    print(f"  ✓ Terms: {n_terms}")
    print(f"  ✓ Odds Ratio range: 0 - {or_max:.1f}")
    print(f"  ✓ -log10(P) range: {min_pval:.1f} - {max_pval:.1f}")


# ============================================================
# GENERATE FINAL PLOT
# ============================================================
print("\n" + "="*60)
print(f"GENERATING FINAL DOT MATRIX PLOT - {SELECTED_VERSION.upper()}")
print("="*60)

plot_dot_matrix_final(
    plot_df=plot_df,
    stability_cols=STABILITY_COLS,
    outfile_prefix=f"Reactome_Hepatocyte_Compartment_{SELECTED_VERSION}",
    title=f"Reactome Pathways ({SELECTED_VERSION.upper()}), P-value < 0.05, Gene_Count > 3",
    figsize_width=3,
    size_scale=30,
    min_size=20
)

print("\n✓ Final plot generated!")


# In[10]:


# ============================================================
# FILTER AND PREPARE DATA (PER SEX × STABILITY CLASS)
# Using P-value (not Adjusted P-value) for consistency
# ============================================================
import numpy as np
import pandas as pd

print("\n" + "="*60)
print("FILTERING AND PREPARING DATA FOR PLOTTING")
print("="*60)

# Store filtered results separately (preserve original all_results)
filtered_results = {}
combined_data = []

for sex in SEXES:
    filtered_results[sex] = {}
    
    for stab_class in STABILITY_CLASSES:
        filtered_results[sex][stab_class] = {}
        
        if stab_class not in all_results[sex]:
            continue
        
        for gs, df in all_results[sex][stab_class].items():
            if df is None or df.empty:
                continue
            
            # Filter criteria - using P-value
            df_filtered = df[
                (df['P-value'] < 0.05) &
                (df['Gene_Count'] > 2) &
                (np.isfinite(df['Odds Ratio']))
            ].copy()
            
            # Store filtered version
            filtered_results[sex][stab_class][gs] = df_filtered
            
            if df_filtered.empty:
                print(f"  {sex} | {stab_class} | {gs}: No terms pass filter")
                continue
            
            # Sort by Odds Ratio and take top 10
            df_top = (
                df_filtered
                .sort_values('Odds Ratio', ascending=False)
                .head(10)
            )
            
            # Add metadata columns
            df_top = df_top.copy()
            df_top['Sex'] = sex
            df_top['Stability'] = stab_class
            df_top['GeneSet'] = gs
            
            combined_data.append(df_top)
            print(f"  {sex:6} | {stab_class:20} | {gs}: Selected {len(df_top)} terms")

# Combine all
if combined_data:
    plot_df = pd.concat(combined_data, ignore_index=True)
    print(f"\n✓ Total terms for plotting: {len(plot_df)}")
else:
    print("❌ No data to plot!")
    plot_df = pd.DataFrame()

# ============================================================
# SUMMARY BY SEX × STABILITY
# ============================================================
if not plot_df.empty:
    print("\n>>> Terms per Sex × Stability:")
    summary = plot_df.groupby(['Sex', 'Stability']).size().reset_index(name='N_terms')
    print(summary.to_string(index=False))
    
    # Check for duplicate pathways (same term in multiple stability classes)
    print("\n>>> Pathways appearing in multiple stability classes:")
    term_counts = plot_df.groupby(['Term', 'Sex']).size().reset_index(name='N_stab_classes')
    duplicates = term_counts[term_counts['N_stab_classes'] > 1]
    if not duplicates.empty:
        print(duplicates.to_string(index=False))
    else:
        print("  None found")


# ============================================================
# FINALIZED DOT MATRIX PLOT - MULTIPLE DOTS PER CELL
# ============================================================
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Rectangle
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable

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

def plot_dot_matrix_multidot(plot_df, stability_cols, outfile_prefix, 
                              title="Reactome Pathways", 
                              figsize_width=8,
                              size_scale=30,
                              min_size=20):
    """
    Dot matrix plot with MULTIPLE DOTS per cell when pathway appears
    in multiple stability classes.
    
    - Rows = Pathways
    - Columns = Male | Female
    - Background color = Max Odds Ratio for that cell
    - Dot color = Stability class (multiple dots if multiple classes)
    - Dot size = -log10(P-value)
    """
    
    df = plot_df.copy()
    
    # Use P-value for consistency with filter
    df['neglog10_pval'] = -np.log10(df['P-value'].clip(lower=1e-300))
    
    # Get all unique terms, sorted by max Odds Ratio (descending)
    term_order = (df.groupby('Term')['Odds Ratio']
                   .max()
                   .sort_values(ascending=False)
                   .index.tolist())
    
    term_to_y = {term: i for i, term in enumerate(term_order)}
    sex_to_x = {'male': 0, 'female': 1}
    
    # Figure sizing
    n_terms = len(term_order)
    fig_height = max(2, 0.25 * n_terms)
    fig, ax = plt.subplots(figsize=(figsize_width, fig_height))
    
    # Color normalization for Odds Ratio background
    or_max = df['Odds Ratio'].max()
    norm = Normalize(vmin=0, vmax=or_max)
    cmap = plt.cm.Blues
    
    # Calculate max OR per (Term, Sex) for background
    max_or_per_cell = df.groupby(['Term', 'Sex'])['Odds Ratio'].max().to_dict()
    
    # Draw background rectangles based on MAX Odds Ratio in each cell
    drawn_cells = set()
    for term in term_order:
        for sex in ['male', 'female']:
            x = sex_to_x[sex]
            y = term_to_y[term]
            
            or_val = max_or_per_cell.get((term, sex), 0)
            
            if or_val > 0:
                facecolor = cmap(norm(or_val))
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
            drawn_cells.add((term, sex))
    
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
    
    # Draw dots - handle multiple stability classes per cell
    for (term, sex), group in df.groupby(['Term', 'Sex']):
        x_base = sex_to_x[sex]
        y_base = term_to_y[term]
        
        n_dots = len(group)
        offsets = get_dot_offsets(n_dots)
        
        # Sort by Odds Ratio so largest dot is drawn last (on top)
        group_sorted = group.sort_values('Odds Ratio', ascending=True)
        
        for i, (_, row) in enumerate(group_sorted.iterrows()):
            stab = row['Stability']
            pval = row['neglog10_pval']
            
            # Get offset
            if i < len(offsets):
                x_off, y_off = offsets[i]
            else:
                x_off, y_off = 0, 0
            
            # Cap size to fit within cell (smaller when multiple dots)
            max_size = 300 if n_dots == 1 else (200 if n_dots <= 2 else 150)
            size = min(pval * size_scale + min_size, max_size)
            
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
    ax.set_title(f'{title}\n', fontsize=9, fontweight='bold', pad=5)
    
    # Invert y-axis (highest OR at top)
    ax.invert_yaxis()
    
    # Remove spines
    for spine in ax.spines.values():
        spine.set_visible(False)
    
    ax.tick_params(axis='both', length=0)
    
    # ===== LEGENDS =====
    
    plt.subplots_adjust(right=0.72)
    
    # 1. Colorbar for Odds Ratio (background)
    sm = ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar_ax = fig.add_axes([0.75, 0.60, 0.02, 0.25])
    cbar = plt.colorbar(sm, cax=cbar_ax)
    cbar.set_label('Odds Ratio', fontsize=9, fontweight='bold')
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
        bbox_to_anchor=(0.74, 0.38),
        framealpha=0.95,
        fontsize=8,
        title_fontsize=9,
        edgecolor='gray'
    )
    
    # 3. Size legend (-log10 P-value)
    min_pval = df['neglog10_pval'].min()
    max_pval = df['neglog10_pval'].max()
    
    if max_pval - min_pval < 2:
        legend_vals = [round(min_pval, 1), round(max_pval, 1)]
    else:
        mid_pval = (min_pval + max_pval) / 2
        legend_vals = [round(min_pval, 1), round(mid_pval, 1), round(max_pval, 1)]
    legend_vals = sorted(set(legend_vals))
    
    size_elements = [
        plt.scatter(
            [], [], 
            s=min(v * size_scale + min_size, 300),
            c='white',
            edgecolor='black',
            linewidth=0.8,
            label=f'{v:.1f}'
        )
        for v in legend_vals
    ]
    
    legend2 = fig.legend(
        handles=size_elements,
        title='-log₁₀(P)',
        loc='center left',
        bbox_to_anchor=(0.74, 0.12),
        framealpha=0.95,
        fontsize=8,
        title_fontsize=10,
        edgecolor='gray'
    )
    
    # Save
    pdf_file = f"{outfile_prefix}_dotmatrix.pdf"
    png_file = f"{outfile_prefix}_dotmatrix.png"
    
    plt.savefig(pdf_file, dpi=300, bbox_inches='tight', facecolor='white')
    plt.savefig(png_file, dpi=300, bbox_inches='tight', facecolor='white')
    plt.show()
    
    print(f"\n  ✓ Saved: {pdf_file}")
    print(f"  ✓ Saved: {png_file}")
    print(f"  ✓ Terms: {n_terms}")
    print(f"  ✓ Odds Ratio range: 0 - {or_max:.1f}")
    print(f"  ✓ -log10(P) range: {min_pval:.1f} - {max_pval:.1f}")
    
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
print(f"GENERATING FINAL DOT MATRIX PLOT - {SELECTED_VERSION.upper()}")
print("="*60)

plot_dot_matrix_multidot(
    plot_df=plot_df,
    stability_cols=STABILITY_COLS,
    outfile_prefix=f"Reactome_Hepatocyte_Compartment_{SELECTED_VERSION}",
    title=f"Reactome Pathways ({SELECTED_VERSION.upper()}), P-value < 0.05, Gene_Count > 2",
    figsize_width=4,
    size_scale=30,
    min_size=20
)

print("\n✓ Final plot generated!")


# In[ ]:




