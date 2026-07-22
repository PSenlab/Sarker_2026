#!/usr/bin/env python
# coding: utf-8
# =============================================================================
# SCENIC+ Complete Analysis Pipeline
# - RSS (Regulon Specificity Scores) by Age AND Age_Sex
# - Correlation Diamond Plot
# - Network Plotting with 4 TF Modules
# =============================================================================

import mudata
import os
import scanpy as sc
import anndata
import matplotlib
import matplotlib.pyplot as plt
import adjustText
import numpy as np
import pandas as pd
import re
import seaborn as sns
from scipy.cluster.hierarchy import linkage, leaves_list
from scenicplus.plotting.dotplot import heatmap_dotplot
from scenicplus.RSS import regulon_specificity_scores, plot_rss
from scenicplus.scenicplus_class import mudata_to_scenicplus
from scenicplus.networks import create_nx_tables, create_nx_graph, plot_networkx

# =============================================================================
# SECTION 1: Load Data
# =============================================================================
print("="*70)
print("LOADING DATA")
print("="*70)

scplus_mdata = mudata.read("scplusmdata.h5mu")

# Filter eRegulon metadata for direct +/+ only
eregulon_metadata = scplus_mdata.uns["direct_e_regulon_metadata"]
filtered_eregulon_metadata = eregulon_metadata[
    eregulon_metadata['eRegulon_name'].str.endswith('direct_+/+')
]
scplus_mdata.uns["direct_e_regulon_metadata"] = filtered_eregulon_metadata

print("[OK] Data loaded successfully")

# =============================================================================
# SECTION 2: Heatmap Dotplots
# =============================================================================
print("\n" + "="*70)
print("GENERATING HEATMAP DOTPLOTS")
print("="*70)

# ---------------------------------------------------------------
# 2A: Dotplot by Age
# ---------------------------------------------------------------
heatmap_dotplot(
    scplus_mudata=scplus_mdata,
    color_modality="direct_gene_based_AUC",
    size_modality="direct_region_based_AUC",
    group_variable="scRNA_counts:age",
    group_variable_order=['young', 'mid_age', 'old', 'pre_geriatric', 'geriatric'],
    eRegulon_metadata_key="direct_e_regulon_metadata",
    color_feature_key="Gene_signature_name",
    size_feature_key="Region_signature_name",
    feature_name_key="eRegulon_name",
    sort_data_by="direct_gene_based_AUC",
    orientation="vertical",
    figsize=(6, 17),
    save="Heps_age.pdf"
)
print("[OK] Saved: Heps_age.pdf")

# ---------------------------------------------------------------
# 2B: Dotplot by Sex
# ---------------------------------------------------------------
heatmap_dotplot(
    scplus_mudata=scplus_mdata,
    color_modality="direct_gene_based_AUC",
    size_modality="direct_region_based_AUC",
    group_variable="scRNA_counts:sex",
    eRegulon_metadata_key="direct_e_regulon_metadata",
    color_feature_key="Gene_signature_name",
    size_feature_key="Region_signature_name",
    feature_name_key="eRegulon_name",
    sort_data_by="direct_gene_based_AUC",
    orientation="horizontal",
    figsize=(13, 5),
    save="Heps_sex_horizontal.pdf"
)
print("[OK] Saved: Heps_sex_horizontal.pdf")

# ---------------------------------------------------------------
# 2C: Dotplot by Cell Type
# ---------------------------------------------------------------
heatmap_dotplot(
    scplus_mudata=scplus_mdata,
    color_modality="direct_gene_based_AUC",
    size_modality="direct_region_based_AUC",
    group_variable="scRNA_counts:celltype2",
    eRegulon_metadata_key="direct_e_regulon_metadata",
    color_feature_key="Gene_signature_name",
    size_feature_key="Region_signature_name",
    feature_name_key="eRegulon_name",
    sort_data_by="direct_gene_based_AUC",
    orientation="vertical",
    figsize=(6, 25),
    save="celltype2.pdf"
)
print("[OK] Saved: celltype2.pdf")

# ---------------------------------------------------------------
# 2D: Dotplot by Age_Sex (Combined)
# ---------------------------------------------------------------
# Create combined age_sex column
scplus_mdata.obs['age_sex'] = (
    scplus_mdata.obs['scRNA_counts:age'].astype(str) + "_" +
    scplus_mdata.obs['scRNA_counts:sex'].astype(str)
)

# Define desired order: females first, then males
desired_order_age_sex = [
    'young_female', 'mid_age_female', 'old_female',
    'pre_geriatric_female', 'geriatric_female',
    'young_male', 'mid_age_male', 'old_male',
    'pre_geriatric_male', 'geriatric_male'
]

# Convert to categorical with specified order
scplus_mdata.obs['age_sex'] = pd.Categorical(
    scplus_mdata.obs['age_sex'],
    categories=desired_order_age_sex,
    ordered=True
)

# Filter to only existing values
actual_values = scplus_mdata.obs['age_sex'].unique()
filtered_order_age_sex = [g for g in desired_order_age_sex if g in actual_values]

heatmap_dotplot(
    scplus_mudata=scplus_mdata,
    color_modality="direct_gene_based_AUC",
    size_modality="direct_region_based_AUC",
    group_variable="age_sex",
    group_variable_order=filtered_order_age_sex,
    eRegulon_metadata_key="direct_e_regulon_metadata",
    color_feature_key="Gene_signature_name",
    size_feature_key="Region_signature_name",
    feature_name_key="eRegulon_name",
    sort_data_by="direct_gene_based_AUC",
    orientation="vertical",
    figsize=(10, 20),
    save="age_sex_dotplot.pdf"
)
print("[OK] Saved: age_sex_dotplot.pdf")

# =============================================================================
# SECTION 3: RSS Calculations - AGE ONLY
# =============================================================================
print("\n" + "="*70)
print("COMPUTING RSS BY AGE (Gene-Based)")
print("="*70)

# ---------------------------------------------------------------
# 3A: Gene-Based RSS by Age
# ---------------------------------------------------------------
scplus_mdata.obs['age'] = scplus_mdata.obs['scRNA_counts:age'].astype(str)
scplus_mdata.obs['age'] = pd.Categorical(
    scplus_mdata.obs['age'],
    categories=['young', 'mid_age', 'old', 'pre_geriatric', 'geriatric'],
    ordered=True
)

rss_gene_age = regulon_specificity_scores(
    scplus_mudata=scplus_mdata,
    variable="age",
    modalities=["direct_gene_based_AUC"]
)

filtered_cols = [
    col for col in rss_gene_age.columns
    if re.search(r'_direct_\+/\+_\(\d+g\)$', col)
]
rss_gene_age = rss_gene_age[filtered_cols]

desired_order_age = ['young', 'mid_age', 'old', 'pre_geriatric', 'geriatric']
rss_gene_age.index = rss_gene_age.index.str.strip()
filtered_order_age = [g for g in desired_order_age if g in rss_gene_age.index]
rss_gene_age = rss_gene_age.loc[filtered_order_age]
rss_gene_age.index = pd.Index([f"{i:02d}_{g}" for i, g in enumerate(filtered_order_age)])

plot_rss(
    data_matrix=rss_gene_age,
    top_n=30,
    num_columns=5,
    save="age_rss_gene_top30.pdf"
)
rss_gene_age.to_csv("age_rss_matrix_gene.csv")
print("[OK] Saved: age_rss_gene_top30.pdf")
print("[OK] Saved: age_rss_matrix_gene.csv")

# ---------------------------------------------------------------
# 3B: Region-Based RSS by Age
# ---------------------------------------------------------------
print("\n" + "="*70)
print("COMPUTING RSS BY AGE (Region-Based)")
print("="*70)

rss_region_age = regulon_specificity_scores(
    scplus_mudata=scplus_mdata,
    variable="age",
    modalities=["direct_region_based_AUC"]
)

filtered_cols = [
    col for col in rss_region_age.columns
    if re.search(r'_direct_\+/\+_\(\d+r\)$', col)
]
rss_region_age = rss_region_age[filtered_cols]

rss_region_age.index = rss_region_age.index.str.strip()
filtered_order_age = [r for r in desired_order_age if r in rss_region_age.index]
rss_region_age = rss_region_age.loc[filtered_order_age]
rss_region_age.index = pd.Index([f"{i:02d}_{r}" for i, r in enumerate(filtered_order_age)])

plot_rss(
    data_matrix=rss_region_age,
    top_n=25,
    num_columns=5,
    save="age_rss_region_top25.pdf"
)
rss_region_age.to_csv("age_rss_matrix_region.csv")
print("[OK] Saved: age_rss_region_top25.pdf")
print("[OK] Saved: age_rss_matrix_region.csv")

# =============================================================================
# SECTION 4: RSS Calculations - AGE_SEX
# =============================================================================
print("\n" + "="*70)
print("COMPUTING RSS BY AGE_SEX (Gene-Based)")
print("="*70)

# ---------------------------------------------------------------
# 4A: Gene-Based RSS by Age_Sex
# ---------------------------------------------------------------
scplus_mdata.obs['age_sex'] = (
    scplus_mdata.obs['scRNA_counts:age'].astype(str) + "_" +
    scplus_mdata.obs['scRNA_counts:sex'].astype(str)
)
scplus_mdata.obs['age_sex'] = pd.Categorical(
    scplus_mdata.obs['age_sex'],
    categories=desired_order_age_sex,
    ordered=True
)

rss_gene_age_sex = regulon_specificity_scores(
    scplus_mudata=scplus_mdata,
    variable="age_sex",
    modalities=["direct_gene_based_AUC"]
)

filtered_cols = [
    col for col in rss_gene_age_sex.columns
    if re.search(r'_direct_\+/\+_\(\d+g\)$', col)
]
rss_gene_age_sex = rss_gene_age_sex[filtered_cols]

rss_gene_age_sex.index = rss_gene_age_sex.index.str.strip()
filtered_order = [g for g in desired_order_age_sex if g in rss_gene_age_sex.index]
rss_gene_age_sex = rss_gene_age_sex.loc[filtered_order]
rss_gene_age_sex.index = pd.Index([f"{i:02d}_{g}" for i, g in enumerate(filtered_order)])

plot_rss(
    data_matrix=rss_gene_age_sex,
    top_n=30,
    num_columns=5,
    save="age_sex_rss_gene_top30.pdf"
)
rss_gene_age_sex.to_csv("age_sex_rss_matrix_gene.csv")
print("[OK] Saved: age_sex_rss_gene_top30.pdf")
print("[OK] Saved: age_sex_rss_matrix_gene.csv")

# ---------------------------------------------------------------
# 4B: Region-Based RSS by Age_Sex
# ---------------------------------------------------------------
print("\n" + "="*70)
print("COMPUTING RSS BY AGE_SEX (Region-Based)")
print("="*70)

rss_region_age_sex = regulon_specificity_scores(
    scplus_mudata=scplus_mdata,
    variable="age_sex",
    modalities=["direct_region_based_AUC"]
)

filtered_cols = [
    col for col in rss_region_age_sex.columns
    if re.search(r'_direct_\+/\+_\(\d+r\)$', col)
]
rss_region_age_sex = rss_region_age_sex[filtered_cols]

rss_region_age_sex.index = rss_region_age_sex.index.str.strip()
filtered_order = [r for r in desired_order_age_sex if r in rss_region_age_sex.index]
rss_region_age_sex = rss_region_age_sex.loc[filtered_order]
rss_region_age_sex.index = pd.Index([f"{i:02d}_{r}" for i, r in enumerate(filtered_order)])

plot_rss(
    data_matrix=rss_region_age_sex,
    top_n=25,
    num_columns=5,
    save="age_sex_rss_region_top25.pdf"
)
rss_region_age_sex.to_csv("age_sex_rss_matrix_region.csv")
print("[OK] Saved: age_sex_rss_region_top25.pdf")
print("[OK] Saved: age_sex_rss_matrix_region.csv")

# =============================================================================
# SECTION 5: Convert to SCENIC+ Object for Networks
# =============================================================================
print("\n" + "="*70)
print("CONVERTING TO SCENICPLUS OBJECT")
print("="*70)

scplus_obj = mudata_to_scenicplus(
    mdata=scplus_mdata,
    path_to_cistarget_h5="ctx_results.hdf5",
    path_to_dem_h5="dem_results.hdf5"
)
print("[OK] Conversion complete")

# =============================================================================
# SECTION 6: Correlation Diamond Plot (Gene vs Region AUC)
# =============================================================================
print("\n" + "="*70)
print("GENERATING CORRELATION DIAMOND PLOT")
print("="*70)

def core_name(name):
    """Remove trailing numeric + modality suffix like (123g) or (123r)"""
    return re.sub(r'_\(\d+[gr]\)$', '', name)

# Use age-based RSS for correlation (can also use age_sex)
rss_gene_corr = rss_gene_age.copy()
rss_region_corr = rss_region_age.copy()
rss_gene_corr.columns = [core_name(c) for c in rss_gene_corr.columns]
rss_region_corr.columns = [core_name(c) for c in rss_region_corr.columns]

# Extract AUC matrices
auc_gene = scplus_obj.uns['eRegulon_AUC']['Gene_based'].copy()
auc_region = scplus_obj.uns['eRegulon_AUC']['Region_based'].copy()
auc_gene.columns = [core_name(c) for c in auc_gene.columns]
auc_region.columns = [core_name(c) for c in auc_region.columns]

print(f"RSS gene-based regulons: {len(rss_gene_corr.columns)}")
print(f"RSS region-based regulons: {len(rss_region_corr.columns)}")
print(f"AUC gene-based regulons: {len(auc_gene.columns)}")
print(f"AUC region-based regulons: {len(auc_region.columns)}")

# Find overlapping regulons
common = list(
    set(auc_gene.columns) & set(auc_region.columns) &
    set(rss_gene_corr.columns) & set(rss_region_corr.columns)
)
print(f"[OK] {len(common)} regulons overlap between RSS + AUC (gene & region)")

if len(common) > 0:
    auc_gene = auc_gene[common]
    auc_region = auc_region[common]
    
    corr_gene = auc_gene.corr(method='pearson')
    corr_region = auc_region.corr(method='pearson')
    
    combined = pd.DataFrame(
        np.full_like(corr_gene, np.nan),
        index=corr_gene.index,
        columns=corr_gene.columns
    )
    for i in range(len(combined)):
        for j in range(len(combined)):
            combined.iat[i, j] = corr_gene.iat[i, j] if i >= j else corr_region.iat[i, j]
    
    link = linkage(combined.fillna(0), method='average')
    order = leaves_list(link)
    combined = combined.iloc[order, order]
    
    plt.figure(figsize=(24, 24))
    sns.heatmap(
        combined,
        cmap='plasma',
        vmin=-1,
        vmax=1,
        square=True,
        cbar_kws={'label': 'Pearson correlation'}
    )
    plt.plot([0, len(combined)], [0, len(combined)], color='black', lw=1)
    plt.text(len(combined)*0.05, len(combined)*0.93, 'Region-based',
             color='white', weight='bold', fontsize=14)
    plt.text(len(combined)*0.70, len(combined)*0.08, 'Gene-based',
             color='white', weight='bold', fontsize=14)
    plt.title('Joint correlation: Region-based (upper) vs Gene-based (lower)', fontsize=16)
    plt.tight_layout()
    plt.savefig('combined_correlation_gene_vs_region_diamond.pdf', bbox_inches='tight', dpi=300)
    plt.show()
    print("[OK] Saved: combined_correlation_gene_vs_region_diamond.pdf")
else:
    print("[WARN] No overlapping regulons found for correlation plot")

# =============================================================================
# SECTION 7: Highly Variable Features
# =============================================================================
print("\n" + "="*70)
print("FINDING HIGHLY VARIABLE FEATURES")
print("="*70)

from pycisTopic.diff_features import find_highly_variable_features

hvr = find_highly_variable_features(
    scplus_obj.to_df('ACC').loc[list(set(scplus_obj.uns['eRegulon_metadata']['Region']))],
    n_top_features=1000,
    plot=False
)
hvg = find_highly_variable_features(
    scplus_obj.to_df('EXP')[list(set(scplus_obj.uns['eRegulon_metadata']['Gene']))].T,
    n_top_features=1000,
    plot=False
)
print(f"[OK] Found {len(hvr)} highly variable regions")
print(f"[OK] Found {len(hvg)} highly variable genes")

# =============================================================================
# SECTION 8: Define 4 TF Modules
# =============================================================================
print("\n" + "="*70)
print("DEFINING TF MODULES")
print("="*70)

# Module 1: Hepatic Function (Green)
hepatic_function_tfs = [
    'Klf13', 'Tef', 'Trim28', 'Nr1h4', 'Nfib', 'Klf12', 'Tcf12', 'Chd2',
    'Pparg', 'Pbx1', 'Nfia', 'Nfe2l2', 'Ppara', 'Foxo3', 'Cux2', 'Thrb',
    'Nr1i3', 'Srebf1', 'Bhlhe40', 'Nr6a1', 'Foxo1', 'Nr1i2', 'Esr1', 'Ar', 'Hnf4a'
]

# Module 2: Immediate Early Genes (Red)
immediate_early_tfs = [
    'Bach2', 'Egr1', 'Nfil3', 'Crem', 'Rora', 'Bach1', 'Fosl2', 'Maff',
    'Fosb', 'Atf3', 'Jun', 'Fos', 'Jund'
]

# Module 3: Inflammation-Related (Cyan)
inflammation_tfs = [
    'Klf10', 'Stat3', 'Junb', 'Ets2', 'Foxa3', 'Nfkb1', 'Onecut1', 'Cbfb',
    'Creb3l2', 'Irf6', 'Atf7', 'Stat1', 'Irf1', 'Relb'
]

# Module 4: Fibro-Balance Module (Blue)
fibro_balance_tfs = [
    'Rxra', 'Nr3c1', 'Phf21a', 'Nr2c1', 'Onecut2', 'Ascl1', 'Nfix', 'Nr2f2',
    'Nfyc', 'Bcl6', 'Smad3', 'Stat5b', 'Ahctf1', 'Irf2', 'Zfp148'
]

MODULE_COLORS = {
    'hepatic_function': '#2ca02c',   # Green
    'immediate_early': '#d62728',     # Red
    'inflammation': '#17becf',        # Cyan
    'fibro_balance': '#1f77b4'        # Blue
}

MODULES = {
    'hepatic_function': hepatic_function_tfs,
    'immediate_early': immediate_early_tfs,
    'inflammation': inflammation_tfs,
    'fibro_balance': fibro_balance_tfs
}

for name, tfs in MODULES.items():
    print(f"  {name}: {len(tfs)} TFs")

# =============================================================================
# SECTION 9: Network Plotting Functions
# =============================================================================

def plot_module_network(
    scplus_obj,
    module_name,
    keep_tfs,
    module_color,
    age_contrast='old',
    figsize=(12, 12),
    save_prefix=None
):
    """Plot network for a single TF module."""
    
    scplus_obj.eRegulon_metadata = scplus_obj.uns["eRegulon_metadata"]
    eregulon_meta = scplus_obj.eRegulon_metadata
    
    direct_pos = eregulon_meta[
        (eregulon_meta['is_extended'] == False) &
        (eregulon_meta['regulation'] == 1)
    ]
    unique_direct_pos_genes = direct_pos['Gene'].unique().tolist()
    
    # FIXED: Create list of eRegulon names (variable, not string!)
    subset_eRegulons = [f"{tf}_direct_+/+" for tf in keep_tfs]
    existing_eRegulons = eregulon_meta['eRegulon_name'].unique()
    subset_eRegulons = [e for e in subset_eRegulons if e in existing_eRegulons]
    
    if len(subset_eRegulons) == 0:
        print(f"[WARN] No eRegulons found for module: {module_name}")
        return None
    
    print(f"\n{'-'*50}")
    print(f"Module: {module_name}")
    print(f"TFs requested: {len(keep_tfs)} | eRegulons found: {len(subset_eRegulons)}")
    print(f"{'-'*50}")
    
    category_color = {tf: module_color for tf in keep_tfs}
    
    try:
        nx_tables = create_nx_tables(
            scplus_obj,
            eRegulon_metadata_key='eRegulon_metadata',
            subset_eRegulons=subset_eRegulons,  # [OK] FIXED: variable, not "subset_eRegulons"
            subset_regions=None,
            subset_genes=unique_direct_pos_genes,
            add_differential_gene_expression=True,
            add_differential_region_accessibility=True,
            differential_variable=['age']
        )
    except Exception as e:
        print(f"Error creating nx_tables: {e}")
        return None
    
    # Filter Edge tables
    nx_tables['Edge']['TF2R'] = nx_tables['Edge']['TF2R'][
        nx_tables['Edge']['TF2R']['TF'].isin(keep_tfs)
    ]
    nx_tables['Edge']['TF2G'] = nx_tables['Edge']['TF2G'][
        nx_tables['Edge']['TF2G']['TF'].isin(keep_tfs)
    ]
    
    regions = set(nx_tables['Edge']['TF2R']['Region'])
    nx_tables['Edge']['R2G'] = nx_tables['Edge']['R2G'][
        nx_tables['Edge']['R2G']['Region'].isin(regions)
    ]
    genes = set(nx_tables['Edge']['R2G']['Gene'])
    genes |= set(nx_tables['Edge']['TF2G']['Gene'])
    
    nx_tables['Edge']['TF2G'] = nx_tables['Edge']['TF2G'][
        nx_tables['Edge']['TF2G']['Gene'].isin(genes)
    ]
    
    # Filter Node tables
    nx_tables['Node']['TF'] = nx_tables['Node']['TF'][
        nx_tables['Node']['TF']['TF'].isin(keep_tfs)
    ].drop_duplicates(subset='TF', keep='first').reset_index(drop=True)
    
    nx_tables['Node']['Region'] = nx_tables['Node']['Region'][
        nx_tables['Node']['Region']['Region'].isin(regions)
    ].drop_duplicates(subset='Region', keep='first').reset_index(drop=True)
    
    nx_tables['Node']['Gene'] = nx_tables['Node']['Gene'][
        nx_tables['Node']['Gene']['Gene'].isin(genes)
    ].drop_duplicates(subset='Gene', keep='first').reset_index(drop=True)
    
    print(f"  TFs: {nx_tables['Node']['TF']['TF'].nunique()}")
    print(f"  Regions: {nx_tables['Node']['Region']['Region'].nunique()}")
    print(f"  Genes: {nx_tables['Node']['Gene']['Gene'].nunique()}")
    
    if nx_tables['Node']['TF'].empty:
        print(f"[WARN] No TF nodes found for module: {module_name}")
        return None
    
    log2fc_col = f'age_Log2FC_{age_contrast}'
    
    G_c, pos_c, edge_tables_c, node_tables_c = create_nx_graph(
        nx_tables,
        use_edge_tables=['TF2R', 'R2G'],
        color_edge_by={
            'TF2R': {'variable': 'TF', 'category_color': category_color},
            'R2G': {'variable': 'rho_R2G', 'continuous_color': 'viridis', 'v_min': -1, 'v_max': 1}
        },
        transparency_edge_by={
            'R2G': {'variable': 'importance_R2G', 'min_alpha': 0.1, 'v_min': 0}
        },
        width_edge_by={
            'R2G': {'variable': 'importance_R2G', 'max_size': 1.5, 'min_size': 1}
        },
        color_node_by={
            'TF': {'variable': 'TF', 'category_color': category_color},
            'Gene': {'variable': log2fc_col, 'continuous_color': 'Oranges'},
            'Region': {'variable': log2fc_col, 'continuous_color': 'Purples'}
        },
        transparency_node_by={
            'Region': {'variable': log2fc_col, 'min_alpha': 0.1},
            'Gene': {'variable': log2fc_col, 'min_alpha': 0.1}
        },
        size_node_by={
            'TF': {'variable': 'fixed_size', 'fixed_size': 12},
            'Gene': {'variable': 'fixed_size', 'fixed_size': 5},
            'Region': {'variable': 'fixed_size', 'fixed_size': 5}
        },
        shape_node_by={
            'TF': {'variable': 'fixed_shape', 'fixed_shape': 'ellipse'},
            'Gene': {'variable': 'fixed_shape', 'fixed_shape': 'ellipse'},
            'Region': {'variable': 'fixed_shape', 'fixed_shape': 'diamond'}
        },
        label_size_by={
            'TF': {'variable': 'fixed_label_size', 'fixed_label_size': 14.0},
            'Gene': {'variable': 'fixed_label_size', 'fixed_label_size': 0.0},
            'Region': {'variable': 'fixed_label_size', 'fixed_label_size': 0.0}
        },
        layout='spring_layout',
        scale_position_by=250
    )
    
    fig, ax = plt.subplots(figsize=figsize)
    plot_networkx(G_c, pos_c)
    plt.title(f'{module_name.replace("_", " ").title()} Network', fontsize=16, fontweight='bold')
    
    if save_prefix:
        save_path = f'{save_prefix}_{module_name}_network.pdf'
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"[OK] Saved: {save_path}")
    
    plt.show()
    
    return G_c, pos_c, nx_tables


def plot_all_modules(scplus_obj, save_prefix='module'):
    """Plot networks for all 4 TF modules."""
    results = {}
    
    for module_name, tfs in MODULES.items():
        color = MODULE_COLORS[module_name]
        result = plot_module_network(
            scplus_obj=scplus_obj,
            module_name=module_name,
            keep_tfs=tfs,
            module_color=color,
            age_contrast='old',
            figsize=(12, 12),
            save_prefix=save_prefix
        )
        if result is not None:
            results[module_name] = result
    
    return results


def plot_combined_network(scplus_obj, figsize=(20, 20), save_path='combined_all_modules_network.pdf'):
    """Plot a single network with all modules, each TF colored by its module."""
    
    scplus_obj.eRegulon_metadata = scplus_obj.uns["eRegulon_metadata"]
    eregulon_meta = scplus_obj.eRegulon_metadata
    
    direct_pos = eregulon_meta[
        (eregulon_meta['is_extended'] == False) &
        (eregulon_meta['regulation'] == 1)
    ]
    unique_direct_pos_genes = direct_pos['Gene'].unique().tolist()
    
    all_tfs = []
    category_color = {}
    
    for module_name, tfs in MODULES.items():
        color = MODULE_COLORS[module_name]
        for tf in tfs:
            all_tfs.append(tf)
            category_color[tf] = color
    
    subset_eRegulons = [f"{tf}_direct_+/+" for tf in all_tfs]
    existing_eRegulons = eregulon_meta['eRegulon_name'].unique()
    subset_eRegulons = [e for e in subset_eRegulons if e in existing_eRegulons]
    
    print(f"\n{'='*60}")
    print(f"Combined Network: All 4 Modules")
    print(f"Total TFs: {len(all_tfs)} | eRegulons found: {len(subset_eRegulons)}")
    print(f"{'='*60}")
    
    nx_tables = create_nx_tables(
        scplus_obj,
        eRegulon_metadata_key='eRegulon_metadata',
        subset_eRegulons=subset_eRegulons,
        subset_regions=None,
        subset_genes=unique_direct_pos_genes,
        add_differential_gene_expression=True,
        add_differential_region_accessibility=True,
        differential_variable=['age']
    )
    
    nx_tables['Edge']['TF2R'] = nx_tables['Edge']['TF2R'][
        nx_tables['Edge']['TF2R']['TF'].isin(all_tfs)
    ]
    nx_tables['Edge']['TF2G'] = nx_tables['Edge']['TF2G'][
        nx_tables['Edge']['TF2G']['TF'].isin(all_tfs)
    ]
    
    regions = set(nx_tables['Edge']['TF2R']['Region'])
    nx_tables['Edge']['R2G'] = nx_tables['Edge']['R2G'][
        nx_tables['Edge']['R2G']['Region'].isin(regions)
    ]
    genes = set(nx_tables['Edge']['R2G']['Gene'])
    genes |= set(nx_tables['Edge']['TF2G']['Gene'])
    
    nx_tables['Edge']['TF2G'] = nx_tables['Edge']['TF2G'][
        nx_tables['Edge']['TF2G']['Gene'].isin(genes)
    ]
    
    nx_tables['Node']['TF'] = nx_tables['Node']['TF'][
        nx_tables['Node']['TF']['TF'].isin(all_tfs)
    ].drop_duplicates(subset='TF', keep='first').reset_index(drop=True)
    
    nx_tables['Node']['Region'] = nx_tables['Node']['Region'][
        nx_tables['Node']['Region']['Region'].isin(regions)
    ].drop_duplicates(subset='Region', keep='first').reset_index(drop=True)
    
    nx_tables['Node']['Gene'] = nx_tables['Node']['Gene'][
        nx_tables['Node']['Gene']['Gene'].isin(genes)
    ].drop_duplicates(subset='Gene', keep='first').reset_index(drop=True)
    
    G_c, pos_c, edge_tables_c, node_tables_c = create_nx_graph(
        nx_tables,
        use_edge_tables=['TF2R', 'R2G'],
        color_edge_by={
            'TF2R': {'variable': 'TF', 'category_color': category_color},
            'R2G': {'variable': 'rho_R2G', 'continuous_color': 'viridis', 'v_min': -1, 'v_max': 1}
        },
        transparency_edge_by={
            'R2G': {'variable': 'importance_R2G', 'min_alpha': 0.1, 'v_min': 0}
        },
        width_edge_by={
            'R2G': {'variable': 'importance_R2G', 'max_size': 1.5, 'min_size': 1}
        },
        color_node_by={
            'TF': {'variable': 'TF', 'category_color': category_color},
            'Gene': {'variable': 'age_Log2FC_old', 'continuous_color': 'Greys'},
            'Region': {'variable': 'age_Log2FC_old', 'continuous_color': 'Greys'}
        },
        transparency_node_by={
            'Region': {'variable': 'age_Log2FC_old', 'min_alpha': 0.1},
            'Gene': {'variable': 'age_Log2FC_old', 'min_alpha': 0.1}
        },
        size_node_by={
            'TF': {'variable': 'fixed_size', 'fixed_size': 15},
            'Gene': {'variable': 'fixed_size', 'fixed_size': 3},
            'Region': {'variable': 'fixed_size', 'fixed_size': 3}
        },
        shape_node_by={
            'TF': {'variable': 'fixed_shape', 'fixed_shape': 'ellipse'},
            'Gene': {'variable': 'fixed_shape', 'fixed_shape': 'ellipse'},
            'Region': {'variable': 'fixed_shape', 'fixed_shape': 'diamond'}
        },
        label_size_by={
            'TF': {'variable': 'fixed_label_size', 'fixed_label_size': 10.0},
            'Gene': {'variable': 'fixed_label_size', 'fixed_label_size': 0.0},
            'Region': {'variable': 'fixed_label_size', 'fixed_label_size': 0.0}
        },
        layout='spring_layout',
        scale_position_by=350
    )
    
    fig, ax = plt.subplots(figsize=figsize)
    plot_networkx(G_c, pos_c)
    
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor=MODULE_COLORS['hepatic_function'], label='Hepatic Function'),
        Patch(facecolor=MODULE_COLORS['immediate_early'], label='Immediate Early Genes'),
        Patch(facecolor=MODULE_COLORS['inflammation'], label='Inflammation-Related'),
        Patch(facecolor=MODULE_COLORS['fibro_balance'], label='Fibro-Balance Module')
    ]
    ax.legend(handles=legend_elements, loc='upper left', fontsize=12)
    
    plt.title('Combined eRegulon Network - All Modules', fontsize=18, fontweight='bold')
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    print(f"[OK] Saved: {save_path}")
    plt.show()
    
    return G_c, pos_c, nx_tables


# =============================================================================
# SECTION 10: Execute Network Plotting
# =============================================================================
print("\n" + "="*70)
print("PLOTTING INDIVIDUAL MODULE NETWORKS")
print("="*70)

results = plot_all_modules(scplus_obj, save_prefix='module')

print("\n" + "="*70)
print("PLOTTING COMBINED NETWORK")
print("="*70)

plot_combined_network(scplus_obj)

print("\n" + "="*70)
print("ANALYSIS COMPLETE!")
print("="*70)
print("\nOutput files generated:")
print("  Dotplots:")
print("    - Heps_age.pdf")
print("    - Heps_sex_horizontal.pdf")
print("    - celltype2.pdf")
print("    - age_sex_dotplot.pdf")
print("  RSS (Age only):")
print("    - age_rss_gene_top30.pdf / .csv")
print("    - age_rss_region_top25.pdf / .csv")
print("  RSS (Age_Sex):")
print("    - age_sex_rss_gene_top30.pdf / .csv")
print("    - age_sex_rss_region_top25.pdf / .csv")
print("  Correlation:")
print("    - combined_correlation_gene_vs_region_diamond.pdf")
print("  Networks:")
print("    - module_hepatic_function_network.pdf")
print("    - module_immediate_early_network.pdf")
print("    - module_inflammation_network.pdf")
print("    - module_fibro_balance_network.pdf")
print("    - combined_all_modules_network.pdf")
