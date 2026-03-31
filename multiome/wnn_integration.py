#!/usr/bin/env python3
#===============================================================================
# Weighted Nearest Neighbors (WNN) Integration for Multiome Data
#===============================================================================
# Description: Integrates snRNA-seq and snATAC-seq modalities using Weighted
#              Nearest Neighbors (WNN) approach via MuOn for joint embedding
#
# Input:       - RNA: Filtered AnnData with scVI embedding (X_scVI_MDE)
#              - ATAC: Filtered AnnData with spectral harmony embedding
# Output:      - Combined MuData object with WNN graph
#              - Individual modalities with shared WNN UMAP
#
# Reference:   Hao Y et al. (2021). Integrated analysis of multimodal single-cell
#              data. Cell, 184(13):3573-3587.e29
#===============================================================================

import logging
import numpy as np
import scanpy as sc
import muon as mu
from muon import atac as ac

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)

#-------------------------------------------------------------------------------
# Configuration
#-------------------------------------------------------------------------------
RNA_INPUT = "rna_shared_atac.h5ad"
ATAC_INPUT = "atac_shared_rna.h5ad"

MUDATA_OUTPUT = "multiome_wnn.h5mu"
RNA_OUTPUT = "rna_wnn.h5ad"
ATAC_OUTPUT = "atac_wnn.h5ad"

# Embedding keys
RNA_REP = "X_scVI_MDE"
ATAC_REP = "X_spectral_harmony"

#-------------------------------------------------------------------------------
# Main Pipeline
#-------------------------------------------------------------------------------
def main():
    # Load data
    logger.info("Loading RNA data...")
    rna = sc.read_h5ad(RNA_INPUT)
    logger.info(f"RNA shape: {rna.shape}")
    
    logger.info("Loading ATAC data...")
    atac = sc.read_h5ad(ATAC_INPUT)
    logger.info(f"ATAC shape: {atac.shape}")
    
    # Create MuData object
    logger.info("Creating MuData object...")
    mdata = mu.MuData({"rna": rna, "atac": atac})
    
    # Compute neighbors for each modality
    logger.info(f"Computing RNA neighbors using {RNA_REP}...")
    sc.pp.neighbors(mdata["rna"], use_rep=RNA_REP)
    
    logger.info(f"Computing ATAC neighbors using {ATAC_REP}...")
    mdata["atac"].X = mdata["atac"].X.astype(float)
    sc.pp.neighbors(mdata["atac"], use_rep=ATAC_REP)
    
    # Compute WNN graph
    logger.info("Computing WNN graph...")
    mu.pp.neighbors(mdata, key_added="wnn")
    logger.info(f"WNN params: {mdata.uns['wnn']['params']['use_rep']}")
    
    # Compute joint UMAP
    logger.info("Computing WNN UMAP...")
    mu.tl.umap(mdata, neighbors_key="wnn", random_state=42)
    mdata.obsm["X_wnn_umap"] = mdata.obsm["X_umap"]
    
    # Transfer WNN embedding to individual modalities
    logger.info("Transferring WNN embedding to modalities...")
    wnn_embedding = mdata.obsm["X_wnn_umap"]
    mdata.mod["rna"].obsm["X_wnn"] = wnn_embedding
    mdata.mod["atac"].obsm["X_wnn"] = wnn_embedding
    
    # Save MuData
    logger.info(f"Saving MuData to {MUDATA_OUTPUT}...")
    mdata.write(MUDATA_OUTPUT)
    
    # Save individual modalities
    logger.info("Saving individual modalities...")
    mdata.mod["rna"].write(RNA_OUTPUT)
    mdata.mod["atac"].write(ATAC_OUTPUT)
    
    logger.info(f"Saved RNA: {mdata.mod['rna'].n_obs:,} cells")
    logger.info(f"Saved ATAC: {mdata.mod['atac'].n_obs:,} cells")
    logger.info("Pipeline complete!")

if __name__ == "__main__":
    main()
