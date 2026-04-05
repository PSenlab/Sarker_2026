#!/usr/bin/env python3
#===============================================================================
# scVI Integration for single-nucleus RNA-seq Data
#===============================================================================
# Description: Batch integration using scVI (single-cell Variational Inference)
#              with Leiden clustering at multiple resolutions
#
# Input:       Filtered AnnData (from snrna_preprocessing.py)
# Output:      Integrated AnnData with scVI latent space, UMAP, and clusters
#
# Reference:   Lopez R et al. (2018). Deep generative modeling for single-cell
#              transcriptomics. Nature Methods, 15(12):1053-1058
#===============================================================================

import logging
import scanpy as sc
import scvi
import torch

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)

# GPU configuration
torch.set_float32_matmul_precision("high")
device = "cuda" if torch.cuda.is_available() else "cpu"
logger.info(f"Using device: {device}")

#-------------------------------------------------------------------------------
# Configuration
#-------------------------------------------------------------------------------
INPUT_FILE = "filtered_singlets.h5ad"
OUTPUT_FILE = "integrated_scvi.h5ad"

# HVG selection
N_TOP_GENES = 3000

# Clustering resolutions
RESOLUTIONS = [0.1, 0.2, 0.3, 0.4, 0.5]

#-------------------------------------------------------------------------------
# Main Pipeline
#-------------------------------------------------------------------------------
def main():
    # Load data
    logger.info(f"Loading {INPUT_FILE}...")
    adata = sc.read_h5ad(INPUT_FILE)
    logger.info(f"Data shape: {adata.shape}")
    
    # Preserve raw counts
    adata.layers["counts"] = adata.X.copy()
    adata.raw = adata
    
    # Normalize and log-transform
    logger.info("Normalizing and log-transforming...")
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    
    # Highly variable genes
    logger.info("Selecting highly variable genes...")
    sc.pp.highly_variable_genes(
        adata,
        flavor="seurat_v3",
        n_top_genes=N_TOP_GENES,
        layer="counts",
        batch_key="sample",
        subset=True
    )
    logger.info(f"HVGs selected: {adata.n_vars}")
    
    # scVI setup and training
    logger.info("Setting up scVI model...")
    scvi.model.SCVI.setup_anndata(adata, layer="counts", batch_key="sample")
    vae = scvi.model.SCVI(adata)
    
    logger.info("Training scVI model...")
    vae.train()
    
    # Extract latent representation
    logger.info("Extracting latent representation...")
    adata.obsm["X_scVI"] = vae.get_latent_representation()
    
    # Neighbors and initial clustering
    logger.info("Computing neighbors...")
    sc.pp.neighbors(adata, use_rep="X_scVI")
    sc.tl.leiden(adata)
    
    # MDE embedding
    logger.info("Computing MDE embedding...")
    adata.obsm["X_scVI_MDE"] = scvi.model.utils.mde(adata.obsm["X_scVI"])
    
    # Leiden clustering at multiple resolutions
    logger.info("Running Leiden clustering at multiple resolutions...")
    for res in RESOLUTIONS:
        key = f"leiden_{int(res * 10)}"
        sc.tl.leiden(adata, resolution=res, key_added=key)
        logger.info(f"  {key}: {adata.obs[key].nunique()} clusters")
    
    # UMAP
    logger.info("Computing UMAP...")
    sc.tl.umap(adata)
    
    # Save
    adata.write_h5ad(OUTPUT_FILE)
    logger.info(f"Saved to {OUTPUT_FILE}")

if __name__ == "__main__":
    main()
