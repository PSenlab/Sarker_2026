#!/usr/bin/env python3
#===============================================================================
# snRNA-seq Preprocessing Pipeline
#===============================================================================
# Description: Complete preprocessing pipeline for single-nucleus RNA-seq data
#              including loading, QC filtering, and doublet removal
#
# Steps:       1. Load SoupX-corrected count matrices
#              2. Calculate QC metrics (mt, ribo, hb)
#              3. Filter cells by QC thresholds
#              4. Detect and remove doublets using Scrublet
#
# Input:       SoupX-corrected 10x matrices (one folder per sample)
# Output:      Filtered AnnData with singlets only (h5ad format)
#
# References:  
#   - SoupX: Young MD, Behjati S (2020). GigaScience, 9(12):giaa151
#   - Scrublet: Wolock SL et al. (2019). Cell Systems, 8(4):281-291.e9
#===============================================================================

import os
import glob
import logging
import numpy as np
import pandas as pd
import scanpy as sc
import scrublet as scr
import matplotlib.pyplot as plt

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)

#-------------------------------------------------------------------------------
# Configuration
#-------------------------------------------------------------------------------
INPUT_DIR = "path/to/soupx/corrected"
OUTPUT_FILE = "filtered_singlets.h5ad"

# Sample prefixes by age group
AGE_GROUPS = {
    "young": [f"Y{i}" for i in range(1, 9)],
    "mid_age": [f"MA{i}" for i in range(1, 9)],
    "old": [f"O{i}" for i in range(1, 9)],
    "pre_ger": [f"PG{i}" for i in range(1, 9)],
    "geriatric": [f"G{i}" for i in range(1, 9)]
}

# QC thresholds
QC_PARAMS = {
    "min_genes": 200,
    "max_genes": 10000,
    "min_counts": 500,
    "max_counts": 100000,
    "max_pct_hb": 10,
    "max_pct_ribo": 30,
    "min_cells_per_gene": 5
}

#-------------------------------------------------------------------------------
# Step 1: Load and Concatenate Samples
#-------------------------------------------------------------------------------
def load_samples(input_dir):
    """Load all samples and return concatenated AnnData."""
    adata_list = []
    sample_names = []
    
    for age_group, samples in AGE_GROUPS.items():
        for sample_id in samples:
            sample_path = os.path.join(input_dir, sample_id)
            
            if not os.path.exists(sample_path):
                logger.warning(f"Sample not found: {sample_id}")
                continue
            
            adata = sc.read_10x_mtx(sample_path, cache=False)
            adata.var_names_make_unique()
            adata.obs["sample"] = sample_id
            adata.obs["age_group"] = age_group
            
            adata_list.append(adata)
            sample_names.append(sample_id)
            logger.info(f"Loaded: {sample_id} ({adata.n_obs} cells)")
    
    # Concatenate
    if len(adata_list) > 1:
        adata = adata_list[0].concatenate(
            *adata_list[1:],
            batch_key="sample",
            batch_categories=sample_names
        )
    else:
        adata = adata_list[0]
    
    logger.info(f"Total cells loaded: {adata.n_obs}")
    return adata

#-------------------------------------------------------------------------------
# Step 2: Calculate QC Metrics
#-------------------------------------------------------------------------------
def calculate_qc_metrics(adata):
    """Annotate genes and calculate QC metrics."""
    
    # Annotate gene categories (mouse gene symbols)
    adata.var["mt"] = adata.var_names.str.startswith("mt-")
    adata.var["ribo"] = adata.var_names.str.startswith(("Rps", "Rpl"))
    adata.var["hb"] = adata.var_names.str.contains("^Hb[bag]")
    
    # Calculate QC metrics
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt", "ribo", "hb"],
        percent_top=None, log1p=False, inplace=True
    )
    
    return adata

#-------------------------------------------------------------------------------
# Step 3: QC Filtering
#-------------------------------------------------------------------------------
def filter_cells(adata, params):
    """Apply QC filters."""
    
    cells_before = adata.obs.groupby("sample").size()
    logger.info(f"Cells before filtering: {adata.n_obs}")
    
    # Filter genes
    sc.pp.filter_genes(adata, min_cells=params["min_cells_per_gene"])
    
    # Filter cells by QC metrics
    adata = adata[adata.obs["pct_counts_hb"] < params["max_pct_hb"]]
    adata = adata[adata.obs["pct_counts_ribo"] < params["max_pct_ribo"]]
    adata = adata[
        (adata.obs["n_genes_by_counts"] > params["min_genes"]) &
        (adata.obs["n_genes_by_counts"] < params["max_genes"])
    ]
    adata = adata[
        (adata.obs["total_counts"] > params["min_counts"]) &
        (adata.obs["total_counts"] < params["max_counts"])
    ].copy()
    
    cells_after = adata.obs.groupby("sample").size()
    logger.info(f"Cells after filtering: {adata.n_obs}")
    
    # Summary
    summary = pd.DataFrame({
        "Before Filtering": cells_before,
        "After Filtering": cells_after
    })
    summary["Change"] = summary["Before Filtering"] - summary["After Filtering"]
    logger.info(f"\n{summary}")
    
    return adata

#-------------------------------------------------------------------------------
# Step 4: Doublet Detection and Removal
#-------------------------------------------------------------------------------
def remove_doublets(adata):
    """Detect and remove doublets using Scrublet."""
    
    logger.info("Running Scrublet doublet detection...")
    
    # Get raw counts
    raw_counts = adata.X.toarray() if hasattr(adata.X, "toarray") else adata.X
    
    # Run Scrublet
    scrub = scr.Scrublet(raw_counts)
    doublet_scores, predicted_doublets = scrub.scrub_doublets()
    
    # Add results to AnnData
    adata.obs["doublet_scores"] = doublet_scores
    adata.obs["predicted_doublets"] = predicted_doublets
    
    n_doublets = predicted_doublets.sum()
    logger.info(f"Detected {n_doublets} doublets ({100*n_doublets/len(predicted_doublets):.1f}%)")
    
    # Filter singlets
    cells_before = adata.obs.groupby("sample").size()
    adata_filtered = adata[~adata.obs["predicted_doublets"]].copy()
    cells_after = adata_filtered.obs.groupby("sample").size()
    
    # Summary
    summary = pd.DataFrame({
        "Before Filtering": cells_before,
        "After Filtering": cells_after
    })
    summary["Change"] = summary["Before Filtering"] - summary["After Filtering"]
    logger.info(f"Cells after doublet removal: {adata_filtered.n_obs}")
    logger.info(f"\n{summary}")
    
    return adata_filtered

#-------------------------------------------------------------------------------
# Main Pipeline
#-------------------------------------------------------------------------------
def main():
    # Load samples
    adata = load_samples(INPUT_DIR)
    
    # Calculate QC metrics
    adata = calculate_qc_metrics(adata)
    
    # Filter cells
    adata = filter_cells(adata, QC_PARAMS)
    
    # Remove doublets
    adata = remove_doublets(adata)
    
    # Save
    adata.write_h5ad(OUTPUT_FILE)
    logger.info(f"Saved to {OUTPUT_FILE}")
    logger.info(f"Final cells per sample:\n{adata.obs['sample'].value_counts()}")

if __name__ == "__main__":
    main()
