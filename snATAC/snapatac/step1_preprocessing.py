#!/usr/bin/env python3
#===============================================================================
# SnapATAC2 Preprocessing Pipeline for single-nucleus ATAC-seq data
#===============================================================================
# Description: Processes snATAC-seq data from multiome experiment using 
#              SnapATAC2 including fragment import, QC filtering, doublet 
#              removal, and batch integration with Harmony
#
# Input:       CellRanger Arc output (atac_fragments.tsv.gz files)
# Output:      Combined AnnDataSet with integrated ATAC data
#
# Reference:   Zhang K et al. (2024). SnapATAC2: A fast, scalable and versatile
#              tool for analysis of single-cell omics data. Nature Methods
#===============================================================================

import os
import numpy as np
import multiprocessing as mp
import snapatac2 as snap

#-------------------------------------------------------------------------------
# Configuration
#-------------------------------------------------------------------------------
INPUT_DIR = "path/to/cellranger_arc/outputs"
OUTPUT_FILE = "combined_atac.h5ads"

# Sample prefixes by age group
AGE_GROUPS = {
    "Young": [f"Y{i}" for i in range(1, 9)],
    "Mid_age": [f"MA{i}" for i in range(1, 9)],
    "Old": [f"O{i}" for i in range(1, 9)],
    "Pre_Geriatric": [f"PG{i}" for i in range(1, 9)],
    "Geriatric": [f"G{i}" for i in range(1, 9)],
}

# QC parameters
QC_PARAMS = {
    "min_num_fragments": 500,
    "min_tsse": 4,
    "n_features": 400000,  # for combined dataset feature selection
}

#-------------------------------------------------------------------------------
# Step 1: Import and Preprocess Individual Samples
#-------------------------------------------------------------------------------
def process_samples(input_dir, age_groups, qc_params):
    """Import fragments and preprocess each sample."""
    
    files = []
    
    # Collect fragment file paths
    for age_group, samples in age_groups.items():
        for sample_id in samples:
            fragment_path = os.path.join(
                input_dir, sample_id, "outs", "atac_fragments.tsv.gz"
            )
            if os.path.exists(fragment_path):
                sample_key = f"{age_group}_{sample_id}"
                files.append((sample_key, fragment_path))
                print(f"Found: {sample_key}")
            else:
                print(f"Not found: {fragment_path}")
    
    # Process each sample
    adatas_list = []
    for name, fragment_file in files:
        print(f"Processing: {name}")
        
        # Import fragments
        adatas = snap.pp.import_data(
            [fragment_file],
            file=[f"{name}.h5ad"],
            chrom_sizes=snap.genome.mm10,
            min_num_fragments=qc_params["min_num_fragments"],
            sorted_by_barcode=False,
        )
        
        # QC and preprocessing
        snap.metrics.tsse(adatas, snap.genome.mm10)
        snap.pp.filter_cells(adatas, min_tsse=qc_params["min_tsse"])
        snap.pp.add_tile_matrix(adatas)
        snap.pp.select_features(adatas)
        
        # Doublet detection and removal
        snap.pp.scrublet(adatas)
        snap.pp.filter_doublets(adatas)
        
        adatas_list.extend(adatas)
    
    return files, adatas_list

#-------------------------------------------------------------------------------
# Step 2: Combine and Integrate Samples
#-------------------------------------------------------------------------------
def combine_and_integrate(files, adatas_list, output_file, n_features):
    """Combine samples and perform batch integration."""
    
    # Create combined dataset
    sample_names = [f[0] for f in files]
    data = snap.AnnDataSet(
        adatas=[(name, adata) for name, adata in zip(sample_names, adatas_list)],
        filename=output_file
    )
    
    print(f"Combined cells: {data.n_obs}")
    print(f"Unique barcodes: {np.unique(data.obs_names).size}")
    
    # Make cell IDs unique
    unique_cell_ids = [f"{sa}:{bc}" for sa, bc in zip(data.obs["sample"], data.obs_names)]
    data.obs_names = unique_cell_ids
    assert data.n_obs == np.unique(data.obs_names).size
    
    # Feature selection
    snap.pp.select_features(data, n_features=n_features)
    
    # Dimensionality reduction
    snap.tl.spectral(data)
    snap.tl.umap(data)
    
    # Batch correction
    snap.pp.mnc_correct(data, batch="sample")
    snap.pp.harmony(data, batch="sample", max_iter_harmony=20)
    
    # Clustering
    snap.pp.knn(data, use_rep="X_spectral_harmony")
    snap.tl.leiden(data)
    
    # Create gene activity matrix
    import scanpy as sc
    gene_matrix = snap.pp.make_gene_matrix(data, snap.genome.mm10)
    
    # Process gene matrix
    sc.pp.filter_genes(gene_matrix, min_cells=5)
    sc.pp.normalize_total(gene_matrix)
    sc.pp.log1p(gene_matrix)
    sc.external.pp.magic(gene_matrix, solver="approximate")
    
    # Transfer UMAP coordinates and save
    gene_matrix.obsm["X_umap"] = data.obsm["X_umap"]
    gene_matrix.write("gene_activity_matrix.h5ad", compression="gzip")
    
    data.close()
    print(f"Saved to {output_file}")

#-------------------------------------------------------------------------------
# Main Pipeline
#-------------------------------------------------------------------------------
def main():
    # Process samples
    files, adatas_list = process_samples(INPUT_DIR, AGE_GROUPS, QC_PARAMS)
    
    # Combine and integrate
    combine_and_integrate(
        files, adatas_list, OUTPUT_FILE, QC_PARAMS["n_features"]
    )

if __name__ == "__main__":
    mp.set_start_method("fork", force=True)
    main()
