#!/usr/bin/env python3
"""
RNA-ATAC Data Integration Preparation for SCENIC+
This script prepares matched RNA (GEX) and ATAC (cisTopic) data for SCENIC+ analysis
by aligning barcodes, filtering to common cells, and ensuring consistent metadata.

"""

# =============================================================================
# IMPORTS
# =============================================================================

import os
import re
import pickle
import warnings
from typing import Optional, List, Tuple

import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc

# Suppress warnings for cleaner output
warnings.filterwarnings('ignore')


# =============================================================================
# CONFIGURATION
# =============================================================================

# Input file paths
RNA_ANNDATA_PATH = "integrated_scvi.h5ad"
CISTOPIC_OBJ_PATH = "cisTopicObject_lda_hep_all_cluster.pkl"

# Output file paths
OUTPUT_DIR = "scenicplus_input"
OUTPUT_GEX_PATH = os.path.join(OUTPUT_DIR, "GEX_anndata.h5ad")
OUTPUT_CISTOPIC_PATH = os.path.join(OUTPUT_DIR, "cisTopic_obj.pkl")

# Cell type to extract (set to None to use all cells)
CELL_TYPE_FILTER = "Hepatocyte"

# Age order for categorical variable
DESIRED_AGE_ORDER = ['young', 'mid_age', 'old', 'pre_geriatric', 'geriatric']

# Normalization parameters
TARGET_SUM = 1e4


# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

def transform_barcode_index(index: str) -> str:
    """
    Transform RNA barcode index to match cisTopic format.
    
    Input format:  ACGTACGT-geriatric_3
    Output format: ACGTACGT-1-geriatric_03___geriatric_03
    
    Parameters
    ----------
    index : str
        Original barcode index from RNA AnnData
        
    Returns
    -------
    str
        Transformed barcode index matching cisTopic format
    """
    # Split the index into barcode and sample group (split at last hyphen)
    parts = index.rsplit('-', maxsplit=1)
    if len(parts) != 2:
        raise ValueError(f"Index format is invalid: {index}")
    
    barcode, sample = parts
    
    # Zero-pad the sample identifier (e.g., geriatric_3 -> geriatric_03)
    sample_padded = re.sub(r'_(\d+)$', lambda x: f"_{int(x.group(1)):02d}", sample)
    
    # Add -1- and the '___sample' suffix
    transformed_index = f"{barcode}-1-{sample_padded}___{sample_padded}"
    
    return transformed_index


def load_rna_anndata(filepath: str, cell_type: Optional[str] = None) -> ad.AnnData:
    """
    Load RNA AnnData object and optionally filter to specific cell type.
    
    Parameters
    ----------
    filepath : str
        Path to the h5ad file
    cell_type : str, optional
        Cell type to filter for (uses 'celltype' column)
        
    Returns
    -------
    ad.AnnData
        Loaded (and optionally filtered) AnnData object
    """
    print(f" Loading RNA AnnData from {filepath}...")
    adata = sc.read(filepath)
    print(f"   [OK] Loaded {adata.n_obs} cells x {adata.n_vars} genes")
    
    if cell_type is not None:
        if 'celltype' not in adata.obs.columns:
            raise ValueError("'celltype' column not found in AnnData.obs")
        
        print(f"\n Filtering to {cell_type} cells...")
        adata = adata[adata.obs['celltype'] == cell_type].copy()
        print(f"   [OK] Filtered to {adata.n_obs} {cell_type} cells")
    
    return adata


def load_cistopic_object(filepath: str):
    """
    Load cisTopic object from pickle file.
    
    Parameters
    ----------
    filepath : str
        Path to the pickle file
        
    Returns
    -------
    CistopicObject
        Loaded cisTopic object
    """
    print(f"\n Loading cisTopic object from {filepath}...")
    with open(filepath, "rb") as f:
        cistopic_obj = pickle.load(f)
    
    n_cells = len(cistopic_obj.cell_names)
    n_regions = len(cistopic_obj.region_names) if hasattr(cistopic_obj, 'region_names') else 'N/A'
    print(f"   [OK] Loaded {n_cells} cells x {n_regions} regions")
    
    return cistopic_obj


def preprocess_rna(adata: ad.AnnData, target_sum: float = 1e4) -> ad.AnnData:
    """
    Preprocess RNA data: store raw, normalize, and log-transform.
    
    Parameters
    ----------
    adata : ad.AnnData
        Raw RNA AnnData object
    target_sum : float
        Target sum for normalization
        
    Returns
    -------
    ad.AnnData
        Preprocessed AnnData object
    """
    print("\n Preprocessing RNA data...")
    
    # Store raw counts
    adata.raw = adata.copy()
    print("   [OK] Stored raw counts in adata.raw")
    
    # Normalize
    sc.pp.normalize_total(adata, target_sum=target_sum)
    print(f"   [OK] Normalized to target sum {target_sum:.0e}")
    
    # Log transform
    sc.pp.log1p(adata)
    print("   [OK] Log-transformed (log1p)")
    
    return adata


def transform_rna_indices(adata: ad.AnnData) -> ad.AnnData:
    """
    Transform RNA AnnData indices to match cisTopic barcode format.
    
    Parameters
    ----------
    adata : ad.AnnData
        RNA AnnData with original indices
        
    Returns
    -------
    ad.AnnData
        AnnData with transformed indices
    """
    print("\n Transforming RNA barcode indices...")
    
    original_indices = adata.obs_names.tolist()
    
    # Transform indices
    transformed_indices = []
    failed_indices = []
    
    for idx in original_indices:
        try:
            transformed_indices.append(transform_barcode_index(idx))
        except ValueError as e:
            failed_indices.append(idx)
            transformed_indices.append(idx)  # Keep original if transformation fails
    
    if failed_indices:
        print(f"   [WARN] {len(failed_indices)} indices could not be transformed")
        print(f"      Examples: {failed_indices[:3]}")
    
    adata.obs_names = pd.Index(transformed_indices)
    print(f"   [OK] Transformed {len(transformed_indices)} indices")
    print(f"      Example: {original_indices[0]} -> {transformed_indices[0]}")
    
    return adata


def remove_duplicates(adata: ad.AnnData, cistopic_obj) -> Tuple[ad.AnnData, object]:
    """
    Remove duplicate indices from both RNA and ATAC data.
    
    Parameters
    ----------
    adata : ad.AnnData
        RNA AnnData object
    cistopic_obj : CistopicObject
        cisTopic object
        
    Returns
    -------
    Tuple[ad.AnnData, CistopicObject]
        Objects with duplicates removed
    """
    print("\n Removing duplicate indices...")
    
    # Check for duplicates in RNA
    rna_dups = adata.obs_names.duplicated()
    n_rna_dups = rna_dups.sum()
    if n_rna_dups > 0:
        print(f"   [WARN] Found {n_rna_dups} duplicate indices in RNA data")
        adata = adata[~rna_dups].copy()
        print(f"   [OK] Removed RNA duplicates: {adata.n_obs} cells remaining")
    
    # Check for duplicates in cisTopic
    cistopic_dups = cistopic_obj.cell_data.index.duplicated()
    n_cistopic_dups = cistopic_dups.sum()
    if n_cistopic_dups > 0:
        print(f"   [WARN] Found {n_cistopic_dups} duplicate indices in cisTopic data")
        cistopic_obj.cell_data = cistopic_obj.cell_data[~cistopic_dups]
        print(f"   [OK] Removed cisTopic duplicates: {len(cistopic_obj.cell_data)} cells remaining")
    
    if n_rna_dups == 0 and n_cistopic_dups == 0:
        print("   [OK] No duplicates found")
    
    return adata, cistopic_obj


def find_matching_cells(adata: ad.AnnData, cistopic_obj) -> Tuple[ad.AnnData, object, set]:
    """
    Find cells present in both RNA and ATAC data.
    
    Parameters
    ----------
    adata : ad.AnnData
        RNA AnnData object
    cistopic_obj : CistopicObject
        cisTopic object
        
    Returns
    -------
    Tuple[ad.AnnData, CistopicObject, set]
        Filtered objects and set of shared indices
    """
    print("\n Finding matching cells between RNA and ATAC...")
    
    rna_indices = set(adata.obs_names)
    atac_indices = set(cistopic_obj.cell_data.index)
    
    shared_indices = rna_indices.intersection(atac_indices)
    missing_from_atac = rna_indices - atac_indices
    missing_from_rna = atac_indices - rna_indices
    
    print(f"   [STATS] RNA cells: {len(rna_indices)}")
    print(f"   [STATS] ATAC cells: {len(atac_indices)}")
    print(f"   [OK] Shared cells: {len(shared_indices)}")
    print(f"   [ERROR] In RNA only: {len(missing_from_atac)}")
    print(f"   [ERROR] In ATAC only: {len(missing_from_rna)}")
    
    if len(shared_indices) == 0:
        raise ValueError("No matching cells found between RNA and ATAC data!")
    
    # Show examples of mismatches if any
    if missing_from_atac:
        print(f"\n   Examples missing from ATAC: {list(missing_from_atac)[:3]}")
    if missing_from_rna:
        print(f"   Examples missing from RNA: {list(missing_from_rna)[:3]}")
    
    return shared_indices


def subset_to_shared_cells(adata: ad.AnnData, cistopic_obj, shared_indices: set) -> Tuple[ad.AnnData, object]:
    """
    Subset both RNA and ATAC data to shared cells and ensure alignment.
    
    Parameters
    ----------
    adata : ad.AnnData
        RNA AnnData object
    cistopic_obj : CistopicObject
        cisTopic object
    shared_indices : set
        Set of shared cell indices
        
    Returns
    -------
    Tuple[ad.AnnData, CistopicObject]
        Subset and aligned objects
    """
    print("\n Subsetting to shared cells and aligning indices...")
    
    shared_indices_list = list(shared_indices)
    
    # Subset RNA
    adata = adata[shared_indices_list, :].copy()
    print(f"   [OK] RNA subset: {adata.n_obs} cells")
    
    # Subset cisTopic cell_data
    cistopic_obj.cell_data = cistopic_obj.cell_data.loc[shared_indices_list]
    print(f"   [OK] cisTopic subset: {len(cistopic_obj.cell_data)} cells")
    
    # Reorder cisTopic to match RNA order
    cistopic_obj.cell_data = cistopic_obj.cell_data.loc[adata.obs_names]
    
    # Verify alignment
    assert list(adata.obs_names) == list(cistopic_obj.cell_data.index), \
        "Index alignment failed!"
    print("   [OK] Indices aligned and verified")
    
    return adata, cistopic_obj


def set_categorical_age(adata: ad.AnnData, cistopic_obj, age_order: List[str]) -> Tuple[ad.AnnData, object]:
    """
    Set age column as ordered categorical in both objects.
    
    Parameters
    ----------
    adata : ad.AnnData
        RNA AnnData object
    cistopic_obj : CistopicObject
        cisTopic object
    age_order : List[str]
        Desired order of age categories
        
    Returns
    -------
    Tuple[ad.AnnData, CistopicObject]
        Objects with categorical age
    """
    print("\n Setting age as ordered categorical...")
    
    # Set in RNA AnnData
    if 'age' in adata.obs.columns:
        # Get unique ages in data
        unique_ages = set(adata.obs['age'].unique())
        valid_order = [a for a in age_order if a in unique_ages]
        
        adata.obs['age'] = pd.Categorical(
            adata.obs['age'],
            categories=valid_order,
            ordered=True
        )
        print(f"   [OK] RNA age categories: {adata.obs['age'].cat.categories.tolist()}")
    else:
        print("   [WARN] 'age' column not found in RNA data")
    
    # Set in cisTopic
    if 'age' in cistopic_obj.cell_data.columns:
        unique_ages = set(cistopic_obj.cell_data['age'].unique())
        valid_order = [a for a in age_order if a in unique_ages]
        
        cistopic_obj.cell_data['age'] = pd.Categorical(
            cistopic_obj.cell_data['age'],
            categories=valid_order,
            ordered=True
        )
        print(f"   [OK] cisTopic age categories: {cistopic_obj.cell_data['age'].cat.categories.tolist()}")
    else:
        print("   [WARN] 'age' column not found in cisTopic data")
    
    return adata, cistopic_obj


def save_outputs(adata: ad.AnnData, cistopic_obj, gex_path: str, cistopic_path: str):
    """
    Save processed RNA AnnData and cisTopic objects.
    
    Parameters
    ----------
    adata : ad.AnnData
        Processed RNA AnnData object
    cistopic_obj : CistopicObject
        Processed cisTopic object
    gex_path : str
        Output path for GEX AnnData
    cistopic_path : str
        Output path for cisTopic object
    """
    print("\n[SAVE] Saving outputs...")
    
    # Create output directory if needed
    output_dir = os.path.dirname(gex_path)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)
        print(f"    Created output directory: {output_dir}")
    
    # Save GEX AnnData
    adata.write(gex_path)
    print(f"   [OK] Saved GEX AnnData: {gex_path}")
    
    # Save cisTopic object
    with open(cistopic_path, "wb") as f:
        pickle.dump(cistopic_obj, f)
    print(f"   [OK] Saved cisTopic object: {cistopic_path}")


def print_final_summary(adata: ad.AnnData, cistopic_obj):
    """Print final summary of processed data."""
    print("\n" + "="*60)
    print("[STATS] FINAL SUMMARY")
    print("="*60)
    
    print(f"\n GEX AnnData:")
    print(f"   Cells: {adata.n_obs}")
    print(f"   Genes: {adata.n_vars}")
    print(f"   Obs columns: {adata.obs.columns.tolist()}")
    if 'age' in adata.obs.columns:
        print(f"   Age distribution:")
        for age, count in adata.obs['age'].value_counts().items():
            print(f"      {age}: {count}")
    
    print(f"\n cisTopic Object:")
    print(f"   Cells: {len(cistopic_obj.cell_data)}")
    if hasattr(cistopic_obj, 'region_names'):
        print(f"   Regions: {len(cistopic_obj.region_names)}")
    if hasattr(cistopic_obj, 'selected_model') and cistopic_obj.selected_model is not None:
        print(f"   Topics: {cistopic_obj.selected_model.n_topic}")
    print(f"   Cell data columns: {cistopic_obj.cell_data.columns.tolist()}")
    
    # Verify alignment
    print(f"\n[OK] Index alignment verified: {list(adata.obs_names[:3])} ...")


# =============================================================================
# MAIN PIPELINE
# =============================================================================

def main(
    rna_path: str = RNA_ANNDATA_PATH,
    cistopic_path: str = CISTOPIC_OBJ_PATH,
    cell_type: Optional[str] = CELL_TYPE_FILTER,
    output_gex: str = OUTPUT_GEX_PATH,
    output_cistopic: str = OUTPUT_CISTOPIC_PATH,
    age_order: List[str] = DESIRED_AGE_ORDER
):
    """
    Run the complete RNA-ATAC integration preparation pipeline.
    
    Parameters
    ----------
    rna_path : str
        Path to input RNA AnnData (.h5ad)
    cistopic_path : str
        Path to input cisTopic object (.pkl)
    cell_type : str, optional
        Cell type to filter for (None for all cells)
    output_gex : str
        Output path for processed GEX AnnData
    output_cistopic : str
        Output path for processed cisTopic object
    age_order : List[str]
        Desired order of age categories
    """
    print("\n" + "="*60)
    print("RNA-ATAC INTEGRATION PREPARATION FOR SCENIC+")
    print("="*60)
    
    # Step 1: Load RNA data
    print("\n" + "-"*40)
    print("STEP 1: Load RNA Data")
    print("-"*40)
    adata = load_rna_anndata(rna_path, cell_type=cell_type)
    
    # Step 2: Preprocess RNA
    print("\n" + "-"*40)
    print("STEP 2: Preprocess RNA")
    print("-"*40)
    adata = preprocess_rna(adata, target_sum=TARGET_SUM)
    
    # Step 3: Transform RNA indices
    print("\n" + "-"*40)
    print("STEP 3: Transform RNA Indices")
    print("-"*40)
    adata = transform_rna_indices(adata)
    
    # Step 4: Load cisTopic object
    print("\n" + "-"*40)
    print("STEP 4: Load cisTopic Object")
    print("-"*40)
    cistopic_obj = load_cistopic_object(cistopic_path)
    
    # Step 5: Remove duplicates
    print("\n" + "-"*40)
    print("STEP 5: Remove Duplicates")
    print("-"*40)
    adata, cistopic_obj = remove_duplicates(adata, cistopic_obj)
    
    # Step 6: Find matching cells
    print("\n" + "-"*40)
    print("STEP 6: Find Matching Cells")
    print("-"*40)
    shared_indices = find_matching_cells(adata, cistopic_obj)
    
    # Step 7: Subset to shared cells
    print("\n" + "-"*40)
    print("STEP 7: Subset to Shared Cells")
    print("-"*40)
    adata, cistopic_obj = subset_to_shared_cells(adata, cistopic_obj, shared_indices)
    
    # Step 8: Set categorical age
    print("\n" + "-"*40)
    print("STEP 8: Set Categorical Age")
    print("-"*40)
    adata, cistopic_obj = set_categorical_age(adata, cistopic_obj, age_order)
    
    # Step 9: Save outputs
    print("\n" + "-"*40)
    print("STEP 9: Save Outputs")
    print("-"*40)
    save_outputs(adata, cistopic_obj, output_gex, output_cistopic)
    
    # Print summary
    print_final_summary(adata, cistopic_obj)
    
    print("\n" + "="*60)
    print("[OK] PIPELINE COMPLETED SUCCESSFULLY")
    print("="*60)
    
    return adata, cistopic_obj


# =============================================================================
# ENTRY POINT
# =============================================================================

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(
        description='Prepare RNA and ATAC data for SCENIC+ integration'
    )
    parser.add_argument('--rna', type=str, default=RNA_ANNDATA_PATH,
                        help='Path to RNA AnnData (.h5ad)')
    parser.add_argument('--atac', type=str, default=CISTOPIC_OBJ_PATH,
                        help='Path to cisTopic object (.pkl)')
    parser.add_argument('--celltype', type=str, default=CELL_TYPE_FILTER,
                        help='Cell type to filter (use "all" for no filter)')
    parser.add_argument('--output-dir', type=str, default=OUTPUT_DIR,
                        help='Output directory')
    parser.add_argument('--output-prefix', type=str, default='',
                        help='Prefix for output files')
    
    args = parser.parse_args()
    
    # Handle cell type argument
    cell_type = None if args.celltype.lower() == 'all' else args.celltype
    
    # Set output paths
    prefix = f"{args.output_prefix}_" if args.output_prefix else ""
    output_gex = os.path.join(args.output_dir, f"{prefix}GEX_anndata.h5ad")
    output_cistopic = os.path.join(args.output_dir, f"{prefix}cisTopic_obj.pkl")
    
    # Run pipeline
    adata, cistopic_obj = main(
        rna_path=args.rna,
        cistopic_path=args.atac,
        cell_type=cell_type,
        output_gex=output_gex,
        output_cistopic=output_cistopic
    )
