#!/usr/bin/env python3
#===============================================================================
# Multiome Integration: RNA and ATAC Label Transfer
#===============================================================================
# Description: Integrates single-nucleus RNA-seq and ATAC-seq data from 
#              multiome experiment by transferring cell type annotations from
#              RNA to ATAC based on shared cell barcodes
#
# Input:       - RNA: Annotated snRNA-seq AnnData (from scVI integration)
#              - ATAC: SnapATAC2 processed AnnDataSet
# Output:      - Annotated RNA and ATAC objects with matched barcodes
#              - Shared cell subset for downstream multiome analysis
#
# Steps:       1. Standardize barcode formats between modalities
#              2. Transfer cell type labels from RNA to ATAC
#              3. Generate UMAP visualizations
#              4. Export matched datasets
#===============================================================================

import os
import logging
import numpy as np
import pandas as pd
import scanpy as sc
import snapatac2 as snap
import warnings

warnings.filterwarnings("ignore")
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)

#-------------------------------------------------------------------------------
# Configuration
#-------------------------------------------------------------------------------
RNA_INPUT = "rna_annotated.h5ad"
ATAC_INPUT = "combined_atac.h5ads"

RNA_OUTPUT = "rna_shared_atac.h5ad"
ATAC_OUTPUT = "atac_shared_rna.h5ad"

# Cell type color mapping
CELLTYPE_COLORS = {
    "Hepatocyte": "#17becf",
    "Endothelial_01": "#a6cee3",
    "Endothelial_02": "#1f78b4",
    "Stellate": "#fb9a99",
    "Kupffer": "#b2df8a",
    "lymp_T": "#fdbf6f",
    "MoMFs": "#33a02c",
    "Cholangiocyte_01": "#e31a1c",
    "Cholangiocyte_02": "#cab2d6",
    "lymp_B": "#e377c2",
    "Unassigned": "#d3d3d3",
    "not_annotated": "#696969",
}

#-------------------------------------------------------------------------------
# Barcode Standardization Functions
#-------------------------------------------------------------------------------
def standardize_suffix(name):
    """Standardize sample suffixes to match between RNA and ATAC."""
    name = name.replace("_0", "_")
    if "pre_geriatric" in name:
        return name.replace("pre_geriatric_pg_", "pre_ger_")
    elif "middle_age" in name:
        return name.replace("middle_age_ma_", "mid_age_")
    elif "young" in name:
        return name.replace("young_y", "young_")
    elif "old" in name:
        return name.replace("old_o", "old_")
    elif "geriatric" in name:
        return name.replace("geriatric_g", "geriatric_")
    return name

def standardize_atac_suffix(atac_suffix):
    """Transform ATAC suffixes to match RNA naming convention."""
    if "Pre_Geriatric" in atac_suffix:
        return atac_suffix.replace("Pre_Geriatric_PG_", "pre_ger_").lower()
    elif "Middle_age" in atac_suffix:
        return atac_suffix.replace("Middle_age_MA_", "mid_age_").lower()
    elif "Young" in atac_suffix:
        return atac_suffix.replace("Young_Y", "young_").lower()
    elif "Old" in atac_suffix:
        return atac_suffix.replace("Old_O", "old_").lower()
    elif "Geriatric" in atac_suffix:
        return atac_suffix.replace("Geriatric_G", "geriatric_").lower()
    return atac_suffix.lower()

def standardize_obs_names(obs_names, suffix_standardizer):
    """Standardize observation names with given suffix standardizer."""
    standardized_names = []
    for name in obs_names:
        parts = name.split('-')
        barcode = '-'.join(parts[:-1])
        suffix = parts[-1]
        standardized_suffix = suffix_standardizer(suffix)
        standardized_names.append(f"{barcode}-{standardized_suffix}")
    return pd.Index(standardized_names)

#-------------------------------------------------------------------------------
# Main Pipeline
#-------------------------------------------------------------------------------
def main():
    # Load data
    logger.info("Loading RNA data...")
    rna = sc.read_h5ad(RNA_INPUT)
    logger.info(f"RNA shape: {rna.shape}")
    
    logger.info("Loading ATAC data...")
    data = snap.read_dataset(ATAC_INPUT)
    atac = data.to_adata()
    logger.info(f"ATAC shape: {atac.shape}")
    
    # Create unique cell IDs for ATAC
    logger.info("Standardizing ATAC barcodes...")
    unique_cell_ids = [f"{sa}:{bc}" for sa, bc in zip(atac.obs['sample'], atac.obs_names)]
    atac.obs_names = unique_cell_ids
    assert atac.n_obs == np.unique(atac.obs_names).size
    
    # Restructure ATAC indices
    atac.obs_names = [
        f"{name.split(':')[1]}-{name.split(':')[0].split('.')[0]}" 
        for name in atac.obs_names
    ]
    
    # Standardize observation names
    logger.info("Standardizing observation names...")
    rna.obs_names = standardize_obs_names(rna.obs_names, standardize_atac_suffix)
    atac.obs_names = standardize_obs_names(atac.obs_names, standardize_atac_suffix)
    
    # Final suffix standardization
    rna.obs_names = rna.obs_names.map(
        lambda name: name.split('-')[0] + '-' + standardize_suffix(name.split('-')[-1])
    )
    atac.obs_names = [
        name.split('-')[0] + '-' + standardize_suffix(name.split('-')[-1])
        for name in atac.obs_names
    ]
    
    # Find shared barcodes
    rna_obs_set = set(rna.obs_names)
    atac_obs_set = set(atac.obs_names)
    shared_obs = rna_obs_set.intersection(atac_obs_set)
    
    logger.info(f"Shared barcodes: {len(shared_obs):,}")
    logger.info(f"RNA only: {len(rna_obs_set - atac_obs_set):,}")
    logger.info(f"ATAC only: {len(atac_obs_set - rna_obs_set):,}")
    
    # Transfer cell type labels from RNA to ATAC
    logger.info("Transferring cell type labels...")
    celltype_map = dict(zip(rna.obs_names, rna.obs['celltype']))
    atac.obs['celltype'] = pd.Categorical(
        [celltype_map.get(x, 'not_annotated') for x in atac.obs_names]
    )
    
    # Assign colors
    atac.uns["celltype_colors"] = [
        CELLTYPE_COLORS[ct] for ct in atac.obs["celltype"].cat.categories
    ]
    
    # Summary
    logger.info(f"Cell type distribution in ATAC:\n{atac.obs['celltype'].value_counts()}")
    
    # Compute UMAP for ATAC
    logger.info("Computing UMAP for ATAC...")
    snap.tl.umap(atac, use_rep="X_spectral_harmony")
    
    # Plot ATAC UMAP
    logger.info("Generating UMAP plots...")
    sc.pl.umap(
        atac,
        color="celltype",
        title=f"ATAC - Cell Types (Total: {atac.n_obs:,} cells)",
        size=3,
        legend_loc="right margin",
        save="_ATAC_celltypes.pdf"
    )
    
    # Save full annotated objects
    logger.info("Saving annotated objects...")
    rna.write_h5ad("rna_annotated.h5ad")
    atac.write_h5ad("atac_annotated.h5ad")
    
    # Filter to shared barcodes only
    logger.info("Filtering to shared barcodes...")
    rna_shared = rna[rna.obs_names.isin(shared_obs)].copy()
    atac_shared = atac[atac.obs_names.isin(shared_obs)].copy()
    
    # Sort by barcode for alignment
    rna_shared = rna_shared[rna_shared.obs_names.sort_values()]
    atac_shared = atac_shared[atac_shared.obs_names.sort_values()]
    
    # Save shared objects
    rna_shared.write_h5ad(RNA_OUTPUT)
    atac_shared.write_h5ad(ATAC_OUTPUT)
    
    logger.info(f"Saved RNA (shared): {rna_shared.n_obs:,} cells")
    logger.info(f"Saved ATAC (shared): {atac_shared.n_obs:,} cells")
    logger.info("Pipeline complete!")

if __name__ == "__main__":
    main()
