#!/usr/bin/env python3
"""
pycisTopic Preprocessing Pipeline for SCENIC+
Complete version with QC visualization, per-sample doublet detection,
cell annotation, and consistent path handling.

Author: Nishat Sarker
Date: 2025
"""

# === Standard Library ===
import os
import re
import json
import pickle
import subprocess
import time
from collections import Counter

# === Scientific Computing ===
import numpy as np
from numpy import array
import pandas as pd
import polars as pl

# === Plotting ===
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.pyplot import rc_context
import seaborn as sns

# === Scanpy / AnnData ===
import scanpy as sc
import anndata as ad

# === Graph/Clustering ===
import leidenalg

# === pycisTopic Modules ===
import pycisTopic
from pycisTopic.qc import get_barcodes_passing_qc_for_sample
from pycisTopic.cistopic_class import create_cistopic_object_from_fragments, merge
from pycisTopic.pseudobulk_peak_calling import export_pseudobulk, peak_calling
from pycisTopic.iterative_peak_calling import get_consensus_peaks
from pycisTopic.plotting.qc_plot import plot_sample_stats, plot_barcode_stats

# === External Tools / Genomics ===
import pyranges as pr
import requests
import scrublet as scr
import packaging
import pybiomart as pbm

# === Parallelism ===
import ray


# =============================================================================
# CONFIGURATION
# =============================================================================

# Set directories
outDir = "outs_trial"
tmpDir = os.path.join(outDir, "temp_sbatch")

# Ensure directories exist
os.makedirs(tmpDir, exist_ok=True)
os.makedirs(os.path.join(outDir, 'consensus_peak_calling', 'pseudobulk_bed_files'), exist_ok=True)
os.makedirs(os.path.join(outDir, 'consensus_peak_calling', 'pseudobulk_bw_files'), exist_ok=True)
os.makedirs(os.path.join(outDir, 'qc'), exist_ok=True)
os.makedirs(os.path.join(outDir, 'qc', 'plots'), exist_ok=True)
os.makedirs(os.path.join(outDir, 'individual_cistopic_objects'), exist_ok=True)

# Define paths for intermediate results
fragments_dict_path = os.path.join(tmpDir, 'fragments_dict.pkl')
narrow_peak_dict_path = os.path.join(tmpDir, 'narrow_peak_dict.pkl')
filtered_fragments_dict_path = os.path.join(tmpDir, "filtered_fragments_dict.pkl")

# Paths for QC and regions
regions_bed_filename = os.path.join(outDir, "consensus_peak_calling", "consensus_regions.bed")
tss_bed_filename = os.path.join(outDir, "qc", "tss.bed")
path_to_blacklist = "/data/sarkern2/scenicplus1/blacklist/mm10-blacklist.v2.bed"

# Cell annotation file path (for hepatocytes or other cell types)
CELL_DATA_ANNOTATION_PATH = "cell_data_Heps.tsv"

# pycistopic binary path
PYCISTOPIC_BIN = "/data/sarkern2/conda/envs/scenicplus3/bin/pycistopic"

# Processing parameters
N_CPU = 15
DOUBLET_RATE = 0.1
DOUBLET_THRESHOLD = 0.22
PEAK_HALF_WIDTH = 250

# Split pattern - IMPORTANT: Use consistent pattern throughout
# '___' for pseudobulk export (sample___barcode format)
# '-' for cisTopic object creation (barcode-1 format from CellRanger)
SPLIT_PATTERN_PSEUDOBULK = '___'
SPLIT_PATTERN_CISTOPIC = '-'


# =============================================================================
# CATEGORY MAPPINGS
# =============================================================================

categories = {
    "young": {
        "folder_prefix": "Y",
        "sample_prefix": "young_",
        "folder_format": "{}",
        "folder_name": "Young"
    },
    "mid_age": {
        "folder_prefix": "MA_",
        "sample_prefix": "mid_age_",
        "folder_format": "{:02d}",
        "folder_name": "Middle_age"
    },
    "old": {
        "folder_prefix": "O",
        "sample_prefix": "old_",
        "folder_format": "{}",
        "folder_name": "Old"
    },
    "pre_ger": {
        "folder_prefix": "PG_",
        "sample_prefix": "pre_ger_",
        "folder_format": "{:02d}",
        "folder_name": "Pre_Geriatric"
    },
    "geriatric": {
        "folder_prefix": "G",
        "sample_prefix": "geriatric_",
        "folder_format": "{}",
        "folder_name": "Geriatric"
    }
}

sample_counts = 8  # Number of samples per category


# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

def save_paths_to_tsv(paths_dict, filepath):
    """Save dictionary of paths to TSV file."""
    with open(filepath, 'w') as f:
        for key, path in paths_dict.items():
            f.write(f"{key}\t{path}\n")


def load_paths_from_tsv(filepath):
    """Load dictionary of paths from TSV file."""
    paths_dict = {}
    with open(filepath, 'r') as f:
        for line in f:
            v, p = line.strip().split('\t')
            paths_dict[v] = p
    return paths_dict


def save_cistopic_object(cistopic_obj, sample_id, output_dir):
    """Saves a CistopicObject to a file."""
    output_file = os.path.join(output_dir, f"{sample_id}_cistopic_obj.pkl")
    with open(output_file, "wb") as f:
        pickle.dump(cistopic_obj, f)
    print(f"✅ Saved CistopicObject for {sample_id} to {output_file}")


def run_scrublet_per_sample(cistopic_obj, expected_doublet_rate=0.1, threshold=0.22):
    """
    Run Scrublet doublet detection on a single sample's cisTopic object.
    
    Parameters
    ----------
    cistopic_obj : CistopicObject
        Single sample cisTopic object
    expected_doublet_rate : float
        Expected doublet rate for Scrublet
    threshold : float
        Threshold for calling doublets
        
    Returns
    -------
    CistopicObject
        Filtered cisTopic object with only singlets
    int
        Number of doublets detected
    """
    scrub = scr.Scrublet(cistopic_obj.fragment_matrix.T, expected_doublet_rate=expected_doublet_rate)
    doublet_scores, predicted_doublets = scrub.scrub_doublets()
    predicted_doublets = scrub.call_doublets(threshold=threshold)
    
    scrublet_df = pd.DataFrame({
        'Doublet_scores_fragments': scrub.doublet_scores_obs_,
        'Predicted_doublets_fragments': scrub.predicted_doublets_
    }, index=cistopic_obj.cell_names)
    
    cistopic_obj.add_cell_data(scrublet_df, split_pattern=SPLIT_PATTERN_CISTOPIC)
    
    n_doublets = sum(scrublet_df['Predicted_doublets_fragments'])
    singlets = cistopic_obj.cell_data[~cistopic_obj.cell_data.Predicted_doublets_fragments].index.tolist()
    
    if len(singlets) > 0:
        cistopic_obj = cistopic_obj.subset(singlets, copy=True, split_pattern=SPLIT_PATTERN_CISTOPIC)
    
    return cistopic_obj, n_doublets


def regenerate_qc_metrics(sample_id, fragments_file, regions_bed, tss_bed, output_dir):
    """Regenerate QC metrics for a sample using pycistopic CLI."""
    os.makedirs(output_dir, exist_ok=True)
    command = [
        PYCISTOPIC_BIN, "qc",
        "--fragments", fragments_file,
        "--regions", regions_bed,
        "--tss", tss_bed,
        "--output", output_dir
    ]
    try:
        subprocess.run(command, check=True, capture_output=True)
        print(f"🔁 QC metrics regenerated for sample {sample_id}.")
        return True
    except subprocess.CalledProcessError as e:
        print(f"❌ Error regenerating QC for sample {sample_id}: {e}")
        return False


def fix_pyranges_head():
    """Override PyRanges head method to use np.bool_ instead of deprecated np.bool."""
    def head_override(self, n=8):
        """Return the n first rows."""
        subsetter = np.zeros(len(self), dtype=np.bool_)
        subsetter[:n] = True
        return self[subsetter]
    pr.PyRanges.head = head_override


# =============================================================================
# STEP 1: GENERATE FRAGMENTS DICTIONARY
# =============================================================================

def generate_fragments_dict():
    """Generate dictionary mapping sample IDs to fragment file paths."""
    print("\n" + "="*60)
    print("STEP 1: Generating fragments dictionary")
    print("="*60)
    
    fragments_dict = {}
    
    for category, prefixes in categories.items():
        folder_prefix = prefixes["folder_prefix"]
        sample_prefix = prefixes["sample_prefix"]
        folder_format = prefixes["folder_format"]
        folder_name = prefixes["folder_name"]
        
        for i in range(1, sample_counts + 1):
            folder = folder_prefix + folder_format.format(i)
            sample = sample_prefix + f"{i:02d}"
            
            frag_file = os.path.join(
                "/data/sarkern2/multiome_liver/CellRanger",
                folder_name, folder, "outs", "atac_fragments.tsv.gz"
            )
            
            if os.path.exists(frag_file):
                fragments_dict[sample] = frag_file
                print(f"✅ Found: {sample}")
            else:
                print(f"⚠️ Missing: {sample} at {frag_file}")
    
    # Save fragments_dict
    pd.to_pickle(fragments_dict, fragments_dict_path)
    print(f"\n📁 Generated fragments_dict with {len(fragments_dict)} entries")
    
    return fragments_dict


# =============================================================================
# STEP 2: GET CHROMOSOME SIZES
# =============================================================================

def get_chromsizes():
    """Download and format mm10 chromosome sizes."""
    print("\n" + "="*60)
    print("STEP 2: Getting chromosome sizes for mm10")
    print("="*60)
    
    target_url = 'http://hgdownload.cse.ucsc.edu/goldenPath/mm10/bigZips/mm10.chrom.sizes'
    chromsizes = pd.read_csv(target_url, sep='\t', header=None)
    chromsizes.columns = ['Chromosome', 'End']
    chromsizes['Start'] = [0] * chromsizes.shape[0]
    chromsizes = chromsizes.loc[:, ['Chromosome', 'Start', 'End']]
    
    # Adjust Chromosome names to match CellRangerARC annotations
    chromsizes['Chromosome'] = [
        chromsizes['Chromosome'][x].replace('v', '.') 
        for x in range(len(chromsizes['Chromosome']))
    ]
    chromsizes['Chromosome'] = [
        chromsizes['Chromosome'][x].split('_')[1] 
        if len(chromsizes['Chromosome'][x].split('_')) > 1 
        else chromsizes['Chromosome'][x]
        for x in range(len(chromsizes['Chromosome']))
    ]
    
    chromsizes = pr.PyRanges(chromsizes)
    print(f"✅ Loaded chromosome sizes: {len(chromsizes)} entries")
    
    return chromsizes


# =============================================================================
# STEP 3: LOAD AND PROCESS CELL DATA
# =============================================================================

def load_and_process_cell_data(fragments_dict):
    """Load cell metadata and filter to valid samples."""
    print("\n" + "="*60)
    print("STEP 3: Loading and processing cell data")
    print("="*60)
    
    cell_data_path = "cell_data.tsv"
    cell_data = pd.read_csv(cell_data_path, sep="\t", index_col=0)
    print(f"📊 Loaded cell_data with shape: {cell_data.shape}")
    print(f"   Columns: {cell_data.columns.tolist()}")
    
    # Move 'sample' from the index to a regular column
    cell_data = cell_data.reset_index()
    
    # Rename 'sample' to 'sample_id'
    cell_data.rename(columns={"sample": "sample_id"}, inplace=True)
    
    # Zero-pad sample IDs: e.g., geriatric_1 -> geriatric_01
    cell_data['sample_id'] = cell_data['sample_id'].apply(
        lambda x: re.sub(r'_(\d)$', r'_0\1', x)
    )
    
    # Drop rows with NaN values in the 'sample_id' column
    before_drop = len(cell_data)
    cell_data = cell_data.dropna(subset=['sample_id'])
    print(f"⚠️ Dropped {before_drop - len(cell_data)} rows with NaN sample_id")
    
    # Clean the barcode by removing extra suffixes
    cell_data['barcode'] = cell_data['barcode'].astype(str)
    cell_data['barcode'] = cell_data['barcode'].apply(
        lambda x: re.sub(r'-.*$', '-1', x.split('___')[0])
    )
    
    # Set 'barcode' as the new index
    cell_data = cell_data.set_index('barcode')
    
    # Add barcode column for downstream analysis
    cell_data['barcode'] = [str(x).split('___')[0] for x in cell_data.index.tolist()]
    
    # Filter out (sample, age) pairs with no cells
    pair_counts = cell_data.groupby(['sample_id', 'age']).size().reset_index(name='count')
    valid_pairs = pair_counts[pair_counts['count'] > 0][['sample_id', 'age']]
    
    before_merge = len(cell_data)
    cell_data = cell_data.merge(valid_pairs, on=['sample_id', 'age'], how='inner')
    print(f"⚠️ Cells dropped during valid_pairs merge: {before_merge - len(cell_data)}")
    
    # Get valid sample IDs
    filtered_pairs = cell_data.groupby(['sample_id', 'age']).size().reset_index(name='count')
    filtered_pairs = filtered_pairs[filtered_pairs['count'] > 0]
    valid_sample_ids = set(filtered_pairs['sample_id'].unique())
    
    # Final filter
    cell_data = cell_data[cell_data['sample_id'].isin(valid_sample_ids)]
    print(f"✅ Final filtered sample_ids: {len(valid_sample_ids)}")
    
    # Build filtered fragments dictionary
    filtered_fragments_dict = {
        sample: path for sample, path in fragments_dict.items()
        if sample in valid_sample_ids
    }
    
    # Save filtered fragments dict
    pd.to_pickle(filtered_fragments_dict, filtered_fragments_dict_path)
    print(f"✅ Filtered fragments_dict includes {len(filtered_fragments_dict)} samples")
    
    if not filtered_fragments_dict:
        raise ValueError("❌ Filtered fragments_dict is empty. Check your filtering logic.")
    
    return cell_data, filtered_fragments_dict


# =============================================================================
# STEP 4: PSEUDOBULK EXPORT
# =============================================================================

def run_pseudobulk_export(cell_data, filtered_fragments_dict, chromsizes):
    """Export pseudobulk profiles for peak calling."""
    print("\n" + "="*60)
    print("STEP 4: Pseudobulk export")
    print("="*60)
    
    try:
        bw_paths, bed_paths = export_pseudobulk(
            input_data=cell_data,
            variable='age',
            sample_id_col='sample_id',
            chromsizes=chromsizes,
            bed_path=os.path.join(outDir, 'consensus_peak_calling', 'pseudobulk_bed_files'),
            bigwig_path=os.path.join(outDir, 'consensus_peak_calling', 'pseudobulk_bw_files'),
            path_to_fragments=filtered_fragments_dict,
            n_cpu=N_CPU,
            normalize_bigwig=True,
            temp_dir=os.path.join(tmpDir, 'ray_spill'),
            split_pattern=SPLIT_PATTERN_PSEUDOBULK
        )
        print("✅ Pseudobulk export completed successfully.")
        
        # Save paths
        save_paths_to_tsv(bw_paths, os.path.join(outDir, 'consensus_peak_calling', 'bw_paths.tsv'))
        save_paths_to_tsv(bed_paths, os.path.join(outDir, 'consensus_peak_calling', 'bed_paths.tsv'))
        
        return bw_paths, bed_paths
        
    except Exception as e:
        raise RuntimeError(f"❌ Error during export_pseudobulk: {e}")


# =============================================================================
# STEP 5: PEAK CALLING
# =============================================================================

def run_peak_calling(bed_paths):
    """Call peaks using MACS2."""
    print("\n" + "="*60)
    print("STEP 5: Peak calling with MACS2")
    print("="*60)
    
    macs_path = 'macs2'
    macs_outdir = os.path.join(outDir, 'consensus_peak_calling', 'MACS')
    os.makedirs(macs_outdir, exist_ok=True)
    
    narrow_peaks_dict = peak_calling(
        macs_path,
        bed_paths,
        macs_outdir,
        genome_size='mm',
        n_cpu=N_CPU,
        input_format='BEDPE',
        shift=73,
        ext_size=146,
        keep_dup='all',
        q_value=0.05,
        _temp_dir=None
    )
    
    # Save narrow peaks dict
    with open(os.path.join(macs_outdir, 'narrow_peaks_dict.pkl'), 'wb') as f:
        pickle.dump(narrow_peaks_dict, f)
    pd.to_pickle(narrow_peaks_dict, narrow_peak_dict_path)
    
    print(f"✅ Peak calling completed. Found peaks for {len(narrow_peaks_dict)} groups.")
    
    return narrow_peaks_dict


# =============================================================================
# STEP 6: CONSENSUS PEAKS
# =============================================================================

def get_consensus_peaks_wrapper(narrow_peaks_dict, chromsizes):
    """Get consensus peaks from narrow peaks."""
    print("\n" + "="*60)
    print("STEP 6: Getting consensus peaks")
    print("="*60)
    
    # Fix PyRanges head method for NumPy 2.x compatibility
    fix_pyranges_head()
    
    # Convert keys to strings if they are integers
    narrow_peaks_dict = {str(k): v for k, v in narrow_peaks_dict.items()}
    
    # Ensure each PyRanges object has the required columns
    required_columns = ['Chromosome', 'Start', 'End', 'Summit', 'Name', 'Score']
    for key, pyrange in narrow_peaks_dict.items():
        df = pyrange.df
        for col in required_columns:
            if col not in df.columns:
                df[col] = 0
        narrow_peaks_dict[key] = pr.PyRanges(df)
    
    try:
        consensus_peaks = get_consensus_peaks(
            narrow_peaks_dict=narrow_peaks_dict,
            peak_half_width=PEAK_HALF_WIDTH,
            chromsizes=chromsizes,
            path_to_blacklist=path_to_blacklist
        )
        print(f"✅ Consensus peaks obtained: {len(consensus_peaks)} regions")
        
        # Save consensus peaks
        consensus_peaks.to_bed(
            path=regions_bed_filename,
            keep=True,
            compression='infer',
            chain=False
        )
        print(f"✅ Saved consensus peaks to {regions_bed_filename}")
        
        return consensus_peaks
        
    except Exception as e:
        print(f"❌ Error getting consensus peaks: {e}")
        raise


# =============================================================================
# STEP 7: GET TSS ANNOTATIONS
# =============================================================================

def get_tss_annotations():
    """Get TSS annotations for mouse genome."""
    print("\n" + "="*60)
    print("STEP 7: Getting TSS annotations")
    print("="*60)
    
    try:
        subprocess.run(
            [
                PYCISTOPIC_BIN, "tss", "get_tss",
                "--output", tss_bed_filename,
                "--name", "mmusculus_gene_ensembl",
                "--to-chrom-source", "ucsc",
                "--ucsc", "mm10"
            ],
            check=True
        )
        print(f"✅ TSS annotations saved to {tss_bed_filename}")
    except subprocess.CalledProcessError as e:
        print(f"❌ Error running pycistopic tss get_tss: {e}")
        raise


# =============================================================================
# STEP 8: RUN QC
# =============================================================================

def run_qc(filtered_fragments_dict):
    """Run QC on all samples."""
    print("\n" + "="*60)
    print("STEP 8: Running QC")
    print("="*60)
    
    pycistopic_qc_commands_filename = "pycistopic_qc_commands.txt"
    
    # Create command file
    with open(pycistopic_qc_commands_filename, "w") as fh:
        for sample, fragment_filename in filtered_fragments_dict.items():
            sample_output = os.path.join(outDir, 'qc', sample)
            print(
                f"{PYCISTOPIC_BIN} qc",
                f"--fragments {fragment_filename}",
                f"--regions {regions_bed_filename}",
                f"--tss {tss_bed_filename}",
                f"--output {sample_output}",
                sep=" ",
                file=fh,
            )
    
    os.chmod(pycistopic_qc_commands_filename, 0o755)
    
    # Run QC commands in parallel (FIXED: proper shell command)
    try:
        with open("qc_commands_output.log", "w") as output_log:
            subprocess.run(
                f"cat {pycistopic_qc_commands_filename} | parallel -j 4",
                shell=True,
                check=True,
                stdout=output_log,
                stderr=output_log
            )
        print("✅ QC commands executed successfully.")
    except subprocess.CalledProcessError as e:
        print(f"⚠️ Error running pycistopic QC commands in parallel: {e}")
        print("   Attempting sequential execution...")
        
        # Fallback: run sequentially
        with open(pycistopic_qc_commands_filename, 'r') as f:
            for line in f:
                cmd = line.strip()
                if cmd:
                    try:
                        subprocess.run(cmd, shell=True, check=True)
                    except subprocess.CalledProcessError as e2:
                        print(f"⚠️ Failed: {cmd}")
    
    # Wait for filesystem to sync
    time.sleep(5)
    
    # Verify and regenerate missing QC files
    print("\n🔍 Verifying QC outputs...")
    for sample_id in filtered_fragments_dict:
        sample_output_dir = os.path.join(outDir, "qc", sample_id)
        sample_metrics_path = os.path.join(sample_output_dir, f"{sample_id}.fragments_stats_per_cb.parquet")
        
        if not os.path.exists(sample_metrics_path):
            print(f"⚠️ Missing metrics for {sample_id} — regenerating...")
            regenerate_qc_metrics(
                sample_id,
                filtered_fragments_dict[sample_id],
                regions_bed_filename,
                tss_bed_filename,
                sample_output_dir
            )
    
    # Final verification
    for sample_id in filtered_fragments_dict:
        sample_metrics_path = os.path.join(outDir, "qc", sample_id, f"{sample_id}.fragments_stats_per_cb.parquet")
        if not os.path.exists(sample_metrics_path):
            print(f"❌ File still missing: {sample_metrics_path}")
        else:
            print(f"✅ Verified: {sample_id}")


# =============================================================================
# STEP 9: GET BARCODES PASSING QC
# =============================================================================

def get_barcodes_passing_qc(filtered_fragments_dict):
    """Get barcodes passing QC for each sample."""
    print("\n" + "="*60)
    print("STEP 9: Getting barcodes passing QC")
    print("="*60)
    
    pycistopic_qc_output_dir = os.path.join(outDir, "qc")
    sample_id_to_barcodes_passing_filters = {}
    sample_id_to_thresholds = {}
    
    for sample_id in filtered_fragments_dict:
        try:
            barcodes, thresholds = get_barcodes_passing_qc_for_sample(
                sample_id=sample_id,
                pycistopic_qc_output_dir=pycistopic_qc_output_dir,
                unique_fragments_threshold=None,
                tss_enrichment_threshold=None,
                frip_threshold=0,
                use_automatic_thresholds=True,
            )
            sample_id_to_barcodes_passing_filters[sample_id] = barcodes
            sample_id_to_thresholds[sample_id] = thresholds
            print(f"✅ {sample_id}: {len(barcodes)} barcodes passing QC")
        except FileNotFoundError as e:
            print(f"⚠️ File not found for {sample_id}: {e}")
        except Exception as e:
            print(f"⚠️ Error for {sample_id}: {e}")
    
    return sample_id_to_barcodes_passing_filters, sample_id_to_thresholds


# =============================================================================
# STEP 9B: QC VISUALIZATION
# =============================================================================

def generate_qc_plots(filtered_fragments_dict, sample_id_to_barcodes_passing_filters, sample_id_to_thresholds):
    """
    Generate QC visualization plots for all samples.
    
    Parameters
    ----------
    filtered_fragments_dict : dict
        Dictionary of sample IDs to fragment file paths
    sample_id_to_barcodes_passing_filters : dict
        Dictionary of sample IDs to barcodes passing QC
    sample_id_to_thresholds : dict
        Dictionary of sample IDs to QC thresholds
    """
    print("\n" + "="*60)
    print("STEP 9B: Generating QC visualization plots")
    print("="*60)
    
    pycistopic_qc_output_dir = os.path.join(outDir, "qc")
    qc_plots_dir = os.path.join(outDir, "qc", "plots")
    os.makedirs(qc_plots_dir, exist_ok=True)
    
    # =========================================================================
    # 1. Summary statistics across all samples
    # =========================================================================
    print("\n📊 Generating sample-level summary statistics...")
    
    # Collect summary stats for all samples
    summary_stats = []
    for sample_id in filtered_fragments_dict:
        if sample_id in sample_id_to_barcodes_passing_filters:
            n_cells = len(sample_id_to_barcodes_passing_filters[sample_id])
            thresholds = sample_id_to_thresholds.get(sample_id, {})
            summary_stats.append({
                'sample_id': sample_id,
                'n_cells_passing_qc': n_cells,
                'unique_fragments_threshold': thresholds.get('unique_fragments_threshold', 'N/A'),
                'tss_enrichment_threshold': thresholds.get('tss_enrichment_threshold', 'N/A'),
                'frip_threshold': thresholds.get('frip_threshold', 'N/A')
            })
    
    summary_df = pd.DataFrame(summary_stats)
    summary_df.to_csv(os.path.join(qc_plots_dir, 'qc_summary_stats.csv'), index=False)
    print(f"✅ Saved QC summary stats to {os.path.join(qc_plots_dir, 'qc_summary_stats.csv')}")
    
    # Plot cells passing QC per sample
    fig, ax = plt.subplots(figsize=(14, 6))
    bars = ax.bar(summary_df['sample_id'], summary_df['n_cells_passing_qc'], 
                  color='steelblue', edgecolor='black', alpha=0.8)
    ax.set_xlabel('Sample ID', fontsize=12)
    ax.set_ylabel('Number of Cells Passing QC', fontsize=12)
    ax.set_title('Cells Passing QC per Sample', fontsize=14, fontweight='bold')
    ax.tick_params(axis='x', rotation=45)
    plt.xticks(ha='right')
    
    # Add value labels on bars
    for bar, val in zip(bars, summary_df['n_cells_passing_qc']):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 50, 
                str(val), ha='center', va='bottom', fontsize=8)
    
    plt.tight_layout()
    fig.savefig(os.path.join(qc_plots_dir, 'cells_passing_qc_per_sample.png'), dpi=150, bbox_inches='tight')
    fig.savefig(os.path.join(qc_plots_dir, 'cells_passing_qc_per_sample.pdf'), bbox_inches='tight')
    plt.close(fig)
    print(f"✅ Saved cells per sample plot")
    
    # =========================================================================
    # 2. Per-sample barcode statistics plots
    # =========================================================================
    print("\n📊 Generating per-sample barcode statistics plots...")
    
    for sample_id in filtered_fragments_dict:
        try:
            # Check if the required parquet file exists
            sample_metrics_path = os.path.join(
                pycistopic_qc_output_dir, sample_id, 
                f"{sample_id}.fragments_stats_per_cb.parquet"
            )
            
            if not os.path.exists(sample_metrics_path):
                print(f"⚠️ Skipping {sample_id} - metrics file not found")
                continue
            
            # Load metrics
            sample_metrics = pl.read_parquet(sample_metrics_path).to_pandas()
            
            # Create multi-panel QC plot for this sample
            fig, axes = plt.subplots(2, 3, figsize=(15, 10))
            fig.suptitle(f'QC Metrics for {sample_id}', fontsize=14, fontweight='bold')
            
            # 1. Total fragments distribution
            if 'total_fragments_count' in sample_metrics.columns:
                ax = axes[0, 0]
                ax.hist(np.log10(sample_metrics['total_fragments_count'] + 1), 
                       bins=50, color='steelblue', edgecolor='black', alpha=0.7)
                ax.set_xlabel('log10(Total Fragments + 1)')
                ax.set_ylabel('Count')
                ax.set_title('Total Fragments Distribution')
                ax.axvline(x=np.log10(1000), color='red', linestyle='--', label='1000 threshold')
                ax.legend()
            
            # 2. Unique fragments distribution
            if 'unique_fragments_count' in sample_metrics.columns:
                ax = axes[0, 1]
                ax.hist(np.log10(sample_metrics['unique_fragments_count'] + 1), 
                       bins=50, color='forestgreen', edgecolor='black', alpha=0.7)
                ax.set_xlabel('log10(Unique Fragments + 1)')
                ax.set_ylabel('Count')
                ax.set_title('Unique Fragments Distribution')
            
            # 3. TSS enrichment distribution
            if 'tss_enrichment' in sample_metrics.columns:
                ax = axes[0, 2]
                ax.hist(sample_metrics['tss_enrichment'], 
                       bins=50, color='darkorange', edgecolor='black', alpha=0.7)
                ax.set_xlabel('TSS Enrichment Score')
                ax.set_ylabel('Count')
                ax.set_title('TSS Enrichment Distribution')
                ax.axvline(x=2, color='red', linestyle='--', label='2.0 threshold')
                ax.legend()
            
            # 4. FRIP (Fraction of Reads in Peaks)
            if 'fraction_of_fragments_in_peaks' in sample_metrics.columns:
                ax = axes[1, 0]
                ax.hist(sample_metrics['fraction_of_fragments_in_peaks'], 
                       bins=50, color='purple', edgecolor='black', alpha=0.7)
                ax.set_xlabel('Fraction of Fragments in Peaks')
                ax.set_ylabel('Count')
                ax.set_title('FRIP Distribution')
            
            # 5. Duplication ratio
            if 'duplication_ratio' in sample_metrics.columns:
                ax = axes[1, 1]
                ax.hist(sample_metrics['duplication_ratio'], 
                       bins=50, color='crimson', edgecolor='black', alpha=0.7)
                ax.set_xlabel('Duplication Ratio')
                ax.set_ylabel('Count')
                ax.set_title('Duplication Ratio Distribution')
            
            # 6. Nucleosome signal (if available)
            if 'nucleosome_signal' in sample_metrics.columns:
                ax = axes[1, 2]
                ax.hist(sample_metrics['nucleosome_signal'], 
                       bins=50, color='teal', edgecolor='black', alpha=0.7)
                ax.set_xlabel('Nucleosome Signal')
                ax.set_ylabel('Count')
                ax.set_title('Nucleosome Signal Distribution')
            else:
                axes[1, 2].text(0.5, 0.5, 'Nucleosome Signal\nNot Available', 
                              ha='center', va='center', fontsize=12)
                axes[1, 2].set_axis_off()
            
            plt.tight_layout()
            fig.savefig(os.path.join(qc_plots_dir, f'{sample_id}_qc_metrics.png'), dpi=150, bbox_inches='tight')
            plt.close(fig)
            print(f"✅ Saved QC plot for {sample_id}")
            
        except Exception as e:
            print(f"⚠️ Error generating QC plot for {sample_id}: {e}")
    
    # =========================================================================
    # 3. Aggregate QC metrics across all samples
    # =========================================================================
    print("\n📊 Generating aggregate QC plots...")
    
    all_metrics = []
    for sample_id in filtered_fragments_dict:
        sample_metrics_path = os.path.join(
            pycistopic_qc_output_dir, sample_id, 
            f"{sample_id}.fragments_stats_per_cb.parquet"
        )
        if os.path.exists(sample_metrics_path):
            df = pl.read_parquet(sample_metrics_path).to_pandas()
            df['sample_id'] = sample_id
            # Extract age group from sample_id
            age_group = sample_id.rsplit('_', 1)[0]
            df['age_group'] = age_group
            all_metrics.append(df)
    
    if all_metrics:
        combined_metrics = pd.concat(all_metrics, ignore_index=True)
        
        # Plot TSS enrichment by age group
        if 'tss_enrichment' in combined_metrics.columns:
            fig, ax = plt.subplots(figsize=(12, 6))
            age_groups = combined_metrics['age_group'].unique()
            positions = range(len(age_groups))
            
            data_to_plot = [combined_metrics[combined_metrics['age_group'] == ag]['tss_enrichment'].dropna() 
                           for ag in age_groups]
            
            bp = ax.boxplot(data_to_plot, positions=positions, patch_artist=True)
            colors = plt.cm.viridis(np.linspace(0, 1, len(age_groups)))
            for patch, color in zip(bp['boxes'], colors):
                patch.set_facecolor(color)
                patch.set_alpha(0.7)
            
            ax.set_xticks(positions)
            ax.set_xticklabels(age_groups, rotation=45, ha='right')
            ax.set_xlabel('Age Group', fontsize=12)
            ax.set_ylabel('TSS Enrichment', fontsize=12)
            ax.set_title('TSS Enrichment by Age Group', fontsize=14, fontweight='bold')
            ax.axhline(y=2, color='red', linestyle='--', alpha=0.7, label='Threshold')
            ax.legend()
            
            plt.tight_layout()
            fig.savefig(os.path.join(qc_plots_dir, 'tss_enrichment_by_age_group.png'), dpi=150, bbox_inches='tight')
            fig.savefig(os.path.join(qc_plots_dir, 'tss_enrichment_by_age_group.pdf'), bbox_inches='tight')
            plt.close(fig)
            print(f"✅ Saved TSS enrichment by age group plot")
        
        # Plot FRIP by age group
        if 'fraction_of_fragments_in_peaks' in combined_metrics.columns:
            fig, ax = plt.subplots(figsize=(12, 6))
            
            data_to_plot = [combined_metrics[combined_metrics['age_group'] == ag]['fraction_of_fragments_in_peaks'].dropna() 
                           for ag in age_groups]
            
            bp = ax.boxplot(data_to_plot, positions=positions, patch_artist=True)
            for patch, color in zip(bp['boxes'], colors):
                patch.set_facecolor(color)
                patch.set_alpha(0.7)
            
            ax.set_xticks(positions)
            ax.set_xticklabels(age_groups, rotation=45, ha='right')
            ax.set_xlabel('Age Group', fontsize=12)
            ax.set_ylabel('Fraction of Fragments in Peaks (FRIP)', fontsize=12)
            ax.set_title('FRIP by Age Group', fontsize=14, fontweight='bold')
            
            plt.tight_layout()
            fig.savefig(os.path.join(qc_plots_dir, 'frip_by_age_group.png'), dpi=150, bbox_inches='tight')
            fig.savefig(os.path.join(qc_plots_dir, 'frip_by_age_group.pdf'), bbox_inches='tight')
            plt.close(fig)
            print(f"✅ Saved FRIP by age group plot")
        
        # Scatter plot: TSS enrichment vs log10(unique fragments)
        if 'tss_enrichment' in combined_metrics.columns and 'unique_fragments_count' in combined_metrics.columns:
            fig, ax = plt.subplots(figsize=(10, 8))
            
            # Sample for visualization (if too many points)
            plot_data = combined_metrics.sample(n=min(50000, len(combined_metrics)), random_state=42)
            
            scatter = ax.scatter(
                np.log10(plot_data['unique_fragments_count'] + 1),
                plot_data['tss_enrichment'],
                c=plot_data['age_group'].astype('category').cat.codes,
                cmap='viridis',
                alpha=0.3,
                s=5
            )
            
            ax.set_xlabel('log10(Unique Fragments + 1)', fontsize=12)
            ax.set_ylabel('TSS Enrichment', fontsize=12)
            ax.set_title('TSS Enrichment vs Unique Fragments', fontsize=14, fontweight='bold')
            ax.axhline(y=2, color='red', linestyle='--', alpha=0.7, label='TSS threshold')
            ax.axvline(x=np.log10(1000), color='red', linestyle='--', alpha=0.7, label='Fragments threshold')
            
            # Add colorbar
            cbar = plt.colorbar(scatter, ax=ax)
            cbar.set_label('Age Group')
            
            plt.tight_layout()
            fig.savefig(os.path.join(qc_plots_dir, 'tss_vs_fragments_scatter.png'), dpi=150, bbox_inches='tight')
            plt.close(fig)
            print(f"✅ Saved TSS vs fragments scatter plot")
    
    print(f"\n✅ All QC plots saved to {qc_plots_dir}")


# =============================================================================
# STEP 10: CREATE CISTOPIC OBJECTS (WITH PER-SAMPLE DOUBLET DETECTION)
# =============================================================================

def create_cistopic_objects(filtered_fragments_dict, sample_id_to_barcodes_passing_filters):
    """Create cisTopic objects for each sample with per-sample doublet detection."""
    print("\n" + "="*60)
    print("STEP 10: Creating cisTopic objects (with per-sample doublet detection)")
    print("="*60)
    
    individual_cistopic_objects_dir = os.path.join(outDir, "individual_cistopic_objects")
    pycistopic_qc_output_dir = os.path.join(outDir, "qc")
    
    cistopic_obj_list = []
    total_doublets = 0
    total_cells_before = 0
    
    for sample_id in filtered_fragments_dict:
        cistopic_obj_file = os.path.join(individual_cistopic_objects_dir, f"{sample_id}_cistopic_obj.pkl")
        
        if os.path.exists(cistopic_obj_file):
            # Load existing object
            with open(cistopic_obj_file, "rb") as f:
                cistopic_obj = pickle.load(f)
            print(f"📂 Loaded existing CistopicObject for {sample_id}")
            cistopic_obj_list.append(cistopic_obj)
            continue
        
        # Check for metrics file
        sample_metrics_path = os.path.join(
            pycistopic_qc_output_dir, sample_id, 
            f"{sample_id}.fragments_stats_per_cb.parquet"
        )
        
        if not os.path.exists(sample_metrics_path):
            print(f"⚠️ Metrics file missing for {sample_id}: {sample_metrics_path}")
            continue
        
        if sample_id not in sample_id_to_barcodes_passing_filters:
            print(f"⚠️ No barcodes passing QC for {sample_id}. Skipping...")
            continue
        
        # Load metrics
        sample_metrics = pl.read_parquet(sample_metrics_path).to_pandas()
        valid_barcodes = sample_id_to_barcodes_passing_filters[sample_id]
        sample_metrics = sample_metrics.set_index("CB").loc[valid_barcodes]
        
        print(f"\n🔬 Creating CistopicObject for {sample_id}...")
        
        # Create cisTopic object
        cistopic_obj = create_cistopic_object_from_fragments(
            path_to_fragments=filtered_fragments_dict[sample_id],
            path_to_regions=regions_bed_filename,
            path_to_blacklist=path_to_blacklist,
            metrics=sample_metrics,
            valid_bc=valid_barcodes,
            n_cpu=N_CPU,
            project=sample_id,
            split_pattern=SPLIT_PATTERN_CISTOPIC
        )
        
        n_cells_before = len(cistopic_obj.cell_names)
        total_cells_before += n_cells_before
        
        # Run Scrublet per sample (IMPROVED: per-sample doublet detection)
        cistopic_obj, n_doublets = run_scrublet_per_sample(
            cistopic_obj, 
            expected_doublet_rate=DOUBLET_RATE,
            threshold=DOUBLET_THRESHOLD
        )
        
        total_doublets += n_doublets
        n_cells_after = len(cistopic_obj.cell_names)
        
        print(f"   📊 {sample_id}: {n_cells_before} cells → {n_cells_after} singlets ({n_doublets} doublets removed)")
        
        cistopic_obj_list.append(cistopic_obj)
        save_cistopic_object(cistopic_obj, sample_id, individual_cistopic_objects_dir)
    
    print(f"\n📈 Summary: {total_cells_before} total cells, {total_doublets} doublets removed")
    
    return cistopic_obj_list


# =============================================================================
# STEP 11: MERGE CISTOPIC OBJECTS
# =============================================================================

def merge_cistopic_objects(cistopic_obj_list):
    """Merge all cisTopic objects into one."""
    print("\n" + "="*60)
    print("STEP 11: Merging cisTopic objects")
    print("="*60)
    
    if len(cistopic_obj_list) == 0:
        raise ValueError("❌ No cisTopic objects to merge!")
    
    cistopic_obj = merge(cistopic_obj_list)
    print(f"✅ Merged {len(cistopic_obj_list)} cisTopic objects")
    print(f"   📊 Final shape: {cistopic_obj.fragment_matrix.shape}")
    
    # Clean the cell_data by removing rows with any NaN values
    nan_count = cistopic_obj.cell_data.isna().sum().sum()
    print(f"   ⚠️ Found {nan_count} total NaN entries in cell metadata.")
    
    if nan_count > 0:
        print("   🧹 Dropping rows with NaNs from cell_data.")
        cistopic_obj.cell_data = cistopic_obj.cell_data.dropna()
    
    # Save merged object
    merged_path = os.path.join(outDir, 'cisTopicObject_merged_dbl_filtered.pkl')
    with open(merged_path, 'wb') as f:
        pickle.dump(cistopic_obj, f)
    print(f"✅ Saved merged cisTopic object to {merged_path}")
    
    return cistopic_obj


# =============================================================================
# STEP 12: ADD CELL TYPE ANNOTATIONS
# =============================================================================

def add_cell_annotations(cistopic_obj, annotation_file_path=None):
    """
    Add cell type annotations from external file to cisTopic object.
    
    Parameters
    ----------
    cistopic_obj : CistopicObject
        The merged cisTopic object
    annotation_file_path : str, optional
        Path to the annotation TSV file. If None, uses CELL_DATA_ANNOTATION_PATH.
        
    Returns
    -------
    CistopicObject
        Annotated cisTopic object
    """
    print("\n" + "="*60)
    print("STEP 12: Adding cell type annotations")
    print("="*60)
    
    if annotation_file_path is None:
        annotation_file_path = CELL_DATA_ANNOTATION_PATH
    
    if not os.path.exists(annotation_file_path):
        print(f"⚠️ Annotation file not found: {annotation_file_path}")
        print("   Skipping cell annotation step.")
        return cistopic_obj
    
    # Load annotation data
    cell_data = pd.read_csv(annotation_file_path, sep="\t", index_col=0)
    print(f"📊 Loaded annotation data with shape: {cell_data.shape}")
    
    # Check for barcodes ending in something other than -1
    non_one_barcodes = cistopic_obj.cell_data.index[
        cistopic_obj.cell_data.index.str.extract(r'-([0-9]+)')[0] != '1'
    ]
    if len(non_one_barcodes) > 0:
        print(f"⚠️ Found {len(non_one_barcodes)} barcodes with suffix other than -1")
    
    # Transform cell_data index to match cisTopic format
    cell_data = cell_data.copy()
    cell_data.reset_index(inplace=True)
    
    # Extract barcode and sample_id from old index
    # Expected format: BARCODE-sample_id (e.g., ACGTACGT-geriatric_3)
    cell_data[['barcode_part', 'sample_id']] = cell_data['barcode'].str.extract(
        r'^([A-Z0-9]+)-([a-zA-Z0-9_]+)$'
    )
    
    # Zero-pad sample_id (e.g., geriatric_3 → geriatric_03)
    cell_data['sample_id'] = cell_data['sample_id'].apply(
        lambda x: re.sub(r'_(\d)$', r'_0\1', str(x)) if pd.notna(x) else x
    )
    
    # Build cisTopic-style index
    # Format: BARCODE-1-sample_id___sample_id
    suffix_number = '1'
    cell_data['cisTopic_id'] = cell_data.apply(
        lambda row: f"{row['barcode_part']}-{suffix_number}-{row['sample_id']}___{row['sample_id']}" 
        if pd.notna(row['barcode_part']) and pd.notna(row['sample_id']) 
        else None,
        axis=1
    )
    
    # Drop rows with failed index transformation
    before_drop = len(cell_data)
    cell_data = cell_data.dropna(subset=['cisTopic_id'])
    print(f"⚠️ Dropped {before_drop - len(cell_data)} rows with invalid index transformation")
    
    # Set new index
    cell_data = cell_data.set_index('cisTopic_id')
    
    # Remove helper columns
    cols_to_drop = ['barcode', 'barcode_part']
    cols_to_drop = [c for c in cols_to_drop if c in cell_data.columns]
    if cols_to_drop:
        cell_data.drop(columns=cols_to_drop, inplace=True)
    
    print(f"✅ Transformed annotation index. Sample entries:\n   {cell_data.index[:3].tolist()}")
    
    # Check overlap with cisTopic object
    cis_barcodes = set(cistopic_obj.cell_names)
    transformed_barcodes = set(cell_data.index)
    intersect = cis_barcodes.intersection(transformed_barcodes)
    missing_from_cistopic = transformed_barcodes - cis_barcodes
    missing_from_annotation = cis_barcodes - transformed_barcodes
    
    print(f"\n📊 Barcode matching statistics:")
    print(f"   cisTopic cells: {len(cis_barcodes)}")
    print(f"   Annotation cells: {len(transformed_barcodes)}")
    print(f"   ✅ Matched: {len(intersect)}")
    print(f"   ❌ In annotation but not in cisTopic: {len(missing_from_cistopic)}")
    print(f"   ❌ In cisTopic but not in annotation: {len(missing_from_annotation)}")
    
    if len(intersect) == 0:
        print("❌ No matching barcodes found! Check barcode format.")
        return cistopic_obj
    
    # Remove duplicates in cell_data
    duplicates_in_cell_data = cell_data.index[cell_data.index.duplicated()]
    if len(duplicates_in_cell_data) > 0:
        print(f"⚠️ Found {len(duplicates_in_cell_data)} duplicate indices in annotation data")
        cell_data = cell_data[~cell_data.index.duplicated(keep='first')]
    
    # Remove duplicates in CistopicObject cell data
    duplicates_in_cistopic = cistopic_obj.cell_data.index[cistopic_obj.cell_data.index.duplicated()]
    if len(duplicates_in_cistopic) > 0:
        print(f"⚠️ Found {len(duplicates_in_cistopic)} duplicate indices in cisTopic object")
        cistopic_obj.cell_data = cistopic_obj.cell_data[~cistopic_obj.cell_data.index.duplicated(keep='first')]
    
    # Add cell data to cisTopic object
    try:
        cistopic_obj.add_cell_data(cell_data)
        print("✅ Successfully added annotation data to cisTopic object")
    except Exception as e:
        print(f"❌ Error adding cell data: {e}")
        return cistopic_obj
    
    # Clean up NaN values
    nan_count = cistopic_obj.cell_data.isna().sum().sum()
    if nan_count > 0:
        print(f"⚠️ Found {nan_count} NaN values. Cleaning...")
        cistopic_obj.cell_data = cistopic_obj.cell_data.dropna()
    
    print(f"✅ Final annotated cell_data shape: {cistopic_obj.cell_data.shape}")
    
    return cistopic_obj


# =============================================================================
# STEP 13: FILTER AND FINALIZE
# =============================================================================

def filter_and_finalize(cistopic_obj, required_columns=None):
    """
    Filter cells based on required columns and finalize the object.
    
    Parameters
    ----------
    cistopic_obj : CistopicObject
        The annotated cisTopic object
    required_columns : list, optional
        List of columns that must have non-NaN values
        
    Returns
    -------
    CistopicObject
        Filtered and finalized cisTopic object
    """
    print("\n" + "="*60)
    print("STEP 13: Filtering and finalizing")
    print("="*60)
    
    if required_columns is None:
        # Default required columns for SCENIC+
        required_columns = [
            'sample_id', 'celltype', 'age', 'sex'
        ]
    
    # Check which required columns exist
    existing_required = [col for col in required_columns if col in cistopic_obj.cell_data.columns]
    missing_required = [col for col in required_columns if col not in cistopic_obj.cell_data.columns]
    
    if missing_required:
        print(f"⚠️ Missing expected columns: {missing_required}")
    
    if existing_required:
        print(f"📋 Filtering on columns: {existing_required}")
        before_filter = len(cistopic_obj.cell_data)
        cistopic_obj.cell_data = cistopic_obj.cell_data.dropna(subset=existing_required)
        after_filter = len(cistopic_obj.cell_data)
        print(f"   Cells before: {before_filter}")
        print(f"   Cells after: {after_filter}")
        print(f"   Dropped: {before_filter - after_filter}")
    
    # Print sample distribution
    print("\n🧾 Cell counts per sample:")
    sample_names = [x.split('-')[2].split('___')[0] for x in cistopic_obj.cell_names]
    sample_counts_dict = Counter(sample_names)
    for sample, count in sorted(sample_counts_dict.items()):
        print(f"   {sample}: {count} cells")
    
    print(f"\n📊 Total unique samples: {len(sample_counts_dict)}")
    print(f"📊 Total cells: {sum(sample_counts_dict.values())}")
    
    # Final NaN check
    any_nans = cistopic_obj.cell_data.isna().any().any()
    print(f"\n🔍 Any remaining NaN values: {any_nans}")
    
    return cistopic_obj


# =============================================================================
# STEP 14: SAVE FINAL OBJECT
# =============================================================================

def save_final_object(cistopic_obj, suffix="Heps"):
    """Save the final annotated and filtered cisTopic object."""
    print("\n" + "="*60)
    print("STEP 14: Saving final object")
    print("="*60)
    
    output_path = os.path.join(outDir, f'cisTopicObject_filtered_annotated_{suffix}.pkl')
    
    with open(output_path, 'wb') as f:
        pickle.dump(cistopic_obj, f)
    
    print(f"✅ Saved final cisTopic object to {output_path}")
    
    return output_path


# =============================================================================
# STEP 15: FINAL SUMMARY
# =============================================================================

def print_final_summary(cistopic_obj):
    """Print final summary statistics."""
    print("\n" + "="*60)
    print("📊 FINAL OBJECT SUMMARY FOR SCENIC+")
    print("="*60)
    
    print(f"Fragment matrix shape: {cistopic_obj.fragment_matrix.shape}")
    print(f"Number of cells: {len(cistopic_obj.cell_names)}")
    print(f"Number of regions: {len(cistopic_obj.region_names)}")
    print(f"Cell metadata columns: {cistopic_obj.cell_data.columns.tolist()}")
    
    if 'sample_id' in cistopic_obj.cell_data.columns:
        print(f"Samples represented: {cistopic_obj.cell_data['sample_id'].nunique()}")
        
        # Summary by sample
        summary = cistopic_obj.cell_data['sample_id'].value_counts().sort_index()
        print("\n🧾 Cell counts per sample (final):")
        print(summary)
        
        # Save summary
        summary.to_csv(os.path.join(outDir, "final_cell_counts_per_sample.csv"), header=["cell_count"])
    
    # Cell type distribution if available
    if 'celltype' in cistopic_obj.cell_data.columns:
        print("\n🧬 Cell type distribution:")
        celltype_counts = cistopic_obj.cell_data['celltype'].value_counts()
        print(celltype_counts)
    
    # Age distribution if available
    if 'age' in cistopic_obj.cell_data.columns:
        print("\n📅 Age group distribution:")
        age_counts = cistopic_obj.cell_data['age'].value_counts()
        print(age_counts)
    
    # Count cells with complete metadata
    complete_meta = cistopic_obj.cell_data.dropna()
    print(f"\n✅ Final number of cells with complete metadata: {complete_meta.shape[0]}")
    
    # Check for required columns for SCENIC+
    required_cols = ['sample_id', 'celltype', 'age', 'sex']
    missing_cols = [col for col in required_cols if col not in cistopic_obj.cell_data.columns]
    if missing_cols:
        print(f"⚠️ Warning: Missing recommended columns for SCENIC+: {missing_cols}")
    else:
        print("✅ All recommended metadata columns present for SCENIC+")


# =============================================================================
# MAIN EXECUTION
# =============================================================================

def main(run_preprocessing=True, run_annotation=True, annotation_suffix="Heps", generate_plots=True):
    """
    Main execution pipeline.
    
    Parameters
    ----------
    run_preprocessing : bool
        Whether to run the preprocessing steps (1-11)
    run_annotation : bool
        Whether to run the annotation steps (12-14)
    annotation_suffix : str
        Suffix for the output file (e.g., "Heps" for hepatocytes)
    generate_plots : bool
        Whether to generate QC visualization plots
    """
    print("\n" + "="*60)
    print("pycisTopic PREPROCESSING PIPELINE FOR SCENIC+")
    print("="*60)
    
    # =========================================================================
    # PREPROCESSING STEPS (1-11)
    # =========================================================================
    
    if run_preprocessing:
        # Step 1: Generate fragments dictionary
        if os.path.exists(fragments_dict_path):
            print(f"📂 Loading existing fragments_dict from {fragments_dict_path}")
            fragments_dict = pd.read_pickle(fragments_dict_path)
        else:
            fragments_dict = generate_fragments_dict()
        
        # Step 2: Get chromosome sizes
        chromsizes = get_chromsizes()
        
        # Step 3: Load and process cell data
        cell_data, filtered_fragments_dict = load_and_process_cell_data(fragments_dict)
        
        # Step 4: Pseudobulk export (skip if bed files exist)
        bed_paths_file = os.path.join(outDir, 'consensus_peak_calling', 'bed_paths.tsv')
        if os.path.exists(bed_paths_file):
            print(f"📂 Loading existing bed paths from {bed_paths_file}")
            bed_paths = load_paths_from_tsv(bed_paths_file)
            bw_paths = load_paths_from_tsv(os.path.join(outDir, 'consensus_peak_calling', 'bw_paths.tsv'))
        else:
            bw_paths, bed_paths = run_pseudobulk_export(cell_data, filtered_fragments_dict, chromsizes)
        
        # Step 5: Peak calling (skip if narrow peaks exist)
        narrow_peaks_file = os.path.join(outDir, 'consensus_peak_calling', 'MACS', 'narrow_peaks_dict.pkl')
        if os.path.exists(narrow_peaks_file):
            print(f"📂 Loading existing narrow peaks from {narrow_peaks_file}")
            with open(narrow_peaks_file, 'rb') as f:
                narrow_peaks_dict = pickle.load(f)
        else:
            narrow_peaks_dict = run_peak_calling(bed_paths)
        
        # Step 6: Consensus peaks (skip if consensus regions exist)
        if os.path.exists(regions_bed_filename):
            print(f"📂 Consensus regions already exist at {regions_bed_filename}")
        else:
            consensus_peaks = get_consensus_peaks_wrapper(narrow_peaks_dict, chromsizes)
        
        # Step 7: TSS annotations (skip if exists)
        if os.path.exists(tss_bed_filename):
            print(f"📂 TSS annotations already exist at {tss_bed_filename}")
        else:
            get_tss_annotations()
        
        # Step 8: Run QC
        run_qc(filtered_fragments_dict)
        
        # Step 9: Get barcodes passing QC
        sample_id_to_barcodes, sample_id_to_thresholds = get_barcodes_passing_qc(filtered_fragments_dict)
        
        # Step 9B: QC Visualization
        if generate_plots:
            generate_qc_plots(filtered_fragments_dict, sample_id_to_barcodes, sample_id_to_thresholds)
        
        # Step 10: Create cisTopic objects with per-sample doublet detection
        cistopic_obj_list = create_cistopic_objects(filtered_fragments_dict, sample_id_to_barcodes)
        
        # Step 11: Merge objects
        cistopic_obj = merge_cistopic_objects(cistopic_obj_list)
    
    # =========================================================================
    # ANNOTATION STEPS (12-14)
    # =========================================================================
    
    if run_annotation:
        # Load merged object if not running preprocessing
        if not run_preprocessing:
            merged_path = os.path.join(outDir, 'cisTopicObject_merged_dbl_filtered.pkl')
            print(f"📂 Loading merged cisTopic object from {merged_path}")
            with open(merged_path, 'rb') as infile:
                cistopic_obj = pickle.load(infile)
            print("✅ cisTopic object loaded successfully")
            print(f"   meta_data shape: {cistopic_obj.cell_data.shape}")
        
        # Step 12: Add cell annotations
        cistopic_obj = add_cell_annotations(cistopic_obj)
        
        # Step 13: Filter and finalize
        cistopic_obj = filter_and_finalize(cistopic_obj)
        
        # Step 14: Save final object
        save_final_object(cistopic_obj, suffix=annotation_suffix)
    
    # =========================================================================
    # FINAL SUMMARY
    # =========================================================================
    
    # Step 15: Final summary
    print_final_summary(cistopic_obj)
    
    print("\n" + "="*60)
    print("✅ PIPELINE COMPLETED SUCCESSFULLY")
    print("="*60)
    
    return cistopic_obj


# =============================================================================
# ENTRY POINT
# =============================================================================

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description='pycisTopic preprocessing pipeline for SCENIC+')
    parser.add_argument('--skip-preprocessing', action='store_true',
                        help='Skip preprocessing steps (1-11) and only run annotation')
    parser.add_argument('--skip-annotation', action='store_true',
                        help='Skip annotation steps (12-14)')
    parser.add_argument('--skip-plots', action='store_true',
                        help='Skip QC visualization plots')
    parser.add_argument('--suffix', type=str, default='Heps',
                        help='Suffix for output file (e.g., Heps for hepatocytes)')
    parser.add_argument('--annotation-file', type=str, default=None,
                        help='Path to annotation TSV file')
    
    args = parser.parse_args()
    
    # Update annotation file path if provided
    if args.annotation_file:
        CELL_DATA_ANNOTATION_PATH = args.annotation_file
    
    cistopic_obj = main(
        run_preprocessing=not args.skip_preprocessing,
        run_annotation=not args.skip_annotation,
        annotation_suffix=args.suffix,
        generate_plots=not args.skip_plots
    )
