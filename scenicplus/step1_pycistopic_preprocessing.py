#!/usr/bin/env python3
"""
pycisTopic Preprocessing Pipeline for SCENIC+ — ALL COMPARTMENTS
=================================================================
One script covering the four compartment runs of the liver aging multiome
atlas: Hepatocyte, endothelial_Kupffer02, myeloid, and T_ILC.

IMPORTANT — why each compartment is a FULL independent run
----------------------------------------------------------
The compartment subset is applied at Step 3, and it propagates into the peaks:

    Step 3  loads a compartment-specific cell_data.tsv (that compartment's
            barcodes only)
    Step 4  passes that frame to export_pseudobulk  -> pseudobulk BEDs contain
            only that compartment's fragments
    Step 5  MACS2 on those BEDs                     -> compartment-specific peaks
    Step 6  consensus regions                       -> compartment-specific
    Step 8  QC / FRIP computed against those regions
    Step 10 cisTopic objects built on those regions

So peaks are called WITHIN each compartment and are NOT shared or comparable
across compartments. Steps 1-11 therefore cannot be run once and reused; every
compartment gets its own outDir, its own consensus_regions.bed, and its own
merged object. Only chromosome sizes (Step 2) are genuinely shared, so they are
fetched once and passed into the loop.

A second, different subsetting happens at Step 13: Step 10 admits every barcode
passing ATAC QC in the fragments file (not only compartment cells), so the
annotation-based cistopic_obj.subset() in Step 13 is what actually trims the
object down to the compartment.

Per-compartment inputs
----------------------
Each entry in COMPARTMENTS needs a cell_data TSV exported from that
compartment's .h5ad, with:
    index   = barcode (adata.obs_names)
    columns : sample_id (or 'sample'), celltype, age, sex

Usage
-----
    # all four compartments, end to end
    python pycistopic_preprocessing_all.py

    # one compartment
    python pycistopic_preprocessing_all.py --run T_ILC

    # annotation only (merged objects already built)
    python pycistopic_preprocessing_all.py --skip-preprocessing

    # reuse existing QC parquets
    python pycistopic_preprocessing_all.py --run myeloid --skip-qc

Author: Nishat Sarker
"""

# === Standard Library ===
import os
import re
import pickle
import subprocess
import time
import traceback
from collections import Counter

# === Scientific Computing ===
import numpy as np
import pandas as pd
import polars as pl

# === Plotting ===
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# === pycisTopic Modules ===
import pycisTopic
from pycisTopic.qc import get_barcodes_passing_qc_for_sample
from pycisTopic.cistopic_class import create_cistopic_object_from_fragments, merge
from pycisTopic.pseudobulk_peak_calling import export_pseudobulk, peak_calling
from pycisTopic.iterative_peak_calling import get_consensus_peaks

# === External Tools / Genomics ===
import pyranges as pr
import scrublet as scr


# =============================================================================
# SHARED CONFIGURATION (identical for every compartment)
# =============================================================================

path_to_blacklist = "/data/sarkern2/scenicplus1/blacklist/mm10-blacklist.v2.bed"
PYCISTOPIC_BIN = "/data/sarkern2/conda/envs/scenicplus5/bin/pycistopic"
CELLRANGER_BASE_DIR = "/data/sarkern2/multiome_liver/CellRanger"

N_CPU = 15
DOUBLET_RATE = 0.1
DOUBLET_THRESHOLD = 0.22
PEAK_HALF_WIDTH = 250

SPLIT_PATTERN_PSEUDOBULK = '___'   # sample___barcode for pseudobulk export
SPLIT_PATTERN_CISTOPIC = '-'       # barcode-1 for cisTopic object creation

# Variable used to build pseudobulks for peak calling. 'age' pools all cells of
# a compartment per age group (5 robust pseudobulks). Switch to 'celltype' only
# if every subtype has enough cells for MACS2.
PSEUDOBULK_VARIABLE = 'age'


# =============================================================================
# CATEGORY MAPPINGS (40 samples; study layout, shared)
# =============================================================================

categories = {
    "young":     {"folder_prefix": "Y",   "sample_prefix": "young_",
                  "folder_format": "{}",     "folder_name": "Young"},
    "mid_age":   {"folder_prefix": "MA_", "sample_prefix": "mid_age_",
                  "folder_format": "{:02d}", "folder_name": "Middle_age"},
    "old":       {"folder_prefix": "O",   "sample_prefix": "old_",
                  "folder_format": "{}",     "folder_name": "Old"},
    "pre_ger":   {"folder_prefix": "PG_", "sample_prefix": "pre_ger_",
                  "folder_format": "{:02d}", "folder_name": "Pre_Geriatric"},
    "geriatric": {"folder_prefix": "G",   "sample_prefix": "geriatric_",
                  "folder_format": "{}",     "folder_name": "Geriatric"},
}

sample_counts = 8  # samples per category


# =============================================================================
# PER-COMPARTMENT REGISTRY
#   out_dir      : isolated output tree (own peaks, QC, objects)
#   cell_data    : TSV exported from that compartment's .h5ad
#   suffix       : final object name -> cisTopicObject_filtered_annotated_<suffix>.pkl
#   celltype_col : name of the annotation column IN THAT TSV. Compartment
#                  exports use 'celltype'; the all-cell-types atlas export uses
#                  'cell_type'. Whatever it is, it is renamed to 'celltype'
#                  internally by normalize_celltype_col(), so every downstream
#                  step and the final object use one consistent name.
#   pseudobulk_variable : optional per-compartment override of the global
#                  PSEUDOBULK_VARIABLE used to build peak-calling pseudobulks.
# =============================================================================

COMPARTMENTS = {

    # ---- full atlas: every annotated cell type, global peak set --------------
    "all_celltypes": dict(
        out_dir="outs_all_celltypes",
        cell_data="cell_data_all.tsv",
        suffix="all_celltypes",
        celltype_col="cell_type",     # <-- atlas export names it 'cell_type'
    ),

    "Hepatocyte": dict(
        out_dir="outs_Heps",
        cell_data="cell_data_Heps.tsv",
        suffix="Heps",
        celltype_col="celltype",
    ),

    "endothelial_Kupffer02": dict(
        out_dir="outs_endothelial_Kupffer02",
        cell_data="cell_data_endothelial_Kupffer02.tsv",
        suffix="endothelial_Kupffer02",
        celltype_col="celltype",
    ),

    "myeloid": dict(
        out_dir="outs_myeloid",
        cell_data="cell_data_myeloid.tsv",
        suffix="myeloid",
        celltype_col="celltype",
    ),

    "T_ILC": dict(
        out_dir="outs_T_ILC",
        cell_data="cell_data_T_ILC.tsv",
        suffix="T_ILC",
        celltype_col="celltype",
    ),
}


# =============================================================================
# COMPARTMENT CONFIG OBJECT
# =============================================================================

class Cfg:
    """All per-compartment paths, derived from out_dir. Created once per run."""

    def __init__(self, name, out_dir, cell_data, suffix,
                 celltype_col="celltype", pseudobulk_variable=None):
        self.name = name
        self.outDir = out_dir
        self.suffix = suffix
        # name of the annotation column as it appears in the input TSV
        self.celltype_col = celltype_col
        # optional per-compartment override of the pseudobulk grouping variable
        self.pseudobulk_variable = pseudobulk_variable or PSEUDOBULK_VARIABLE
        self.tmpDir = os.path.join(out_dir, "temp_sbatch")

        # cell_data feeds BOTH Step 3 (which cells define the peaks) and Step 12
        self.cell_data_path = cell_data
        self.annotation_path = cell_data

        # intermediates
        self.fragments_dict_path = os.path.join(self.tmpDir, "fragments_dict.pkl")
        self.narrow_peak_dict_path = os.path.join(self.tmpDir, "narrow_peak_dict.pkl")
        self.filtered_fragments_dict_path = os.path.join(
            self.tmpDir, "filtered_fragments_dict.pkl")

        # regions / QC
        self.regions_bed = os.path.join(out_dir, "consensus_peak_calling",
                                        "consensus_regions.bed")
        self.tss_bed = os.path.join(out_dir, "qc", "tss.bed")
        self.qc_dir = os.path.join(out_dir, "qc")
        self.qc_plots_dir = os.path.join(out_dir, "qc", "plots")
        self.individual_objs_dir = os.path.join(out_dir, "individual_cistopic_objects")
        self.merged_path = os.path.join(out_dir, "cisTopicObject_merged_dbl_filtered.pkl")

        # per-compartment command / log filenames (so parallel runs don't collide)
        self.qc_commands_file = f"pycistopic_qc_commands_{suffix}.txt"
        self.qc_log_file = f"qc_commands_output_{suffix}.log"

    def makedirs(self):
        for d in (self.tmpDir,
                  os.path.join(self.outDir, "consensus_peak_calling", "pseudobulk_bed_files"),
                  os.path.join(self.outDir, "consensus_peak_calling", "pseudobulk_bw_files"),
                  self.qc_dir, self.qc_plots_dir, self.individual_objs_dir):
            os.makedirs(d, exist_ok=True)


# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

def save_paths_to_tsv(paths_dict, filepath):
    with open(filepath, 'w') as f:
        for key, path in paths_dict.items():
            f.write(f"{key}\t{path}\n")


def load_paths_from_tsv(filepath):
    paths_dict = {}
    with open(filepath, 'r') as f:
        for line in f:
            v, p = line.strip().split('\t')
            paths_dict[v] = p
    return paths_dict


def save_cistopic_object(cistopic_obj, sample_id, output_dir):
    output_file = os.path.join(output_dir, f"{sample_id}_cistopic_obj.pkl")
    with open(output_file, "wb") as f:
        pickle.dump(cistopic_obj, f)
    print(f"[OK] Saved CistopicObject for {sample_id} to {output_file}")


def run_scrublet_per_sample(cistopic_obj, expected_doublet_rate=0.1, threshold=0.22):
    """Run Scrublet doublet detection on a single sample's cisTopic object."""
    scrub = scr.Scrublet(cistopic_obj.fragment_matrix.T,
                         expected_doublet_rate=expected_doublet_rate)
    doublet_scores, predicted_doublets = scrub.scrub_doublets()
    predicted_doublets = scrub.call_doublets(threshold=threshold)

    scrublet_df = pd.DataFrame({
        'Doublet_scores_fragments': scrub.doublet_scores_obs_,
        'Predicted_doublets_fragments': scrub.predicted_doublets_
    }, index=cistopic_obj.cell_names)

    cistopic_obj.add_cell_data(scrublet_df, split_pattern=SPLIT_PATTERN_CISTOPIC)

    n_doublets = sum(scrublet_df['Predicted_doublets_fragments'])
    singlets = cistopic_obj.cell_data[
        ~cistopic_obj.cell_data.Predicted_doublets_fragments].index.tolist()

    if len(singlets) > 0:
        cistopic_obj = cistopic_obj.subset(singlets, copy=True,
                                           split_pattern=SPLIT_PATTERN_CISTOPIC)

    return cistopic_obj, n_doublets


def regenerate_qc_metrics(sample_id, fragments_file, regions_bed, tss_bed, output_prefix):
    """Regenerate QC metrics via the pycistopic CLI.

    output_prefix is a path PREFIX (e.g. qc/<sample>); pycistopic writes
    <prefix>.fragments_stats_per_cb.parquet, so only the parent dir is created.
    """
    os.makedirs(os.path.dirname(output_prefix), exist_ok=True)
    command = [PYCISTOPIC_BIN, "qc",
               "--fragments", fragments_file,
               "--regions", regions_bed,
               "--tss", tss_bed,
               "--output", output_prefix]
    try:
        subprocess.run(command, check=True, capture_output=True)
        print(f" QC metrics regenerated for sample {sample_id}.")
        return True
    except subprocess.CalledProcessError as e:
        print(f"[ERROR] Error regenerating QC for sample {sample_id}: {e}")
        return False


def fix_pyranges_head():
    """Override PyRanges.head to use np.bool_ (np.bool removed in NumPy 2.x)."""
    def head_override(self, n=8):
        subsetter = np.zeros(len(self), dtype=np.bool_)
        subsetter[:n] = True
        return self[subsetter]
    pr.PyRanges.head = head_override


def clean_to_fragment_barcode(bc):
    """Normalize any obs_name to the 10x fragment barcode form 'BARCODE-1'.

    Handles 'BARCODE-1-sample', 'BARCODE-sample', 'BARCODE-1', 'BARCODE___sample'.
    """
    bc = str(bc).split('___')[0]
    dna = bc.split('-')[0]
    return f"{dna}-1"


def normalize_celltype_col(df, celltype_col, context=""):
    """Rename the annotation column to 'celltype' so downstream code is uniform.

    Compartment exports use 'celltype'; the all-cell-types atlas export uses
    'cell_type'. If neither the configured name nor a known alias is present,
    the frame is returned unchanged and a warning is printed -- Step 13 will then
    simply not filter on celltype.
    """
    if celltype_col in df.columns:
        if celltype_col != 'celltype':
            df = df.rename(columns={celltype_col: 'celltype'})
            print(f"[OK] {context}renamed annotation column "
                  f"'{celltype_col}' -> 'celltype'")
        return df

    # configured name absent -> try the usual aliases before giving up
    for alias in ('celltype', 'cell_type', 'celltype2', 'cell_types'):
        if alias in df.columns:
            if alias != 'celltype':
                df = df.rename(columns={alias: 'celltype'})
            print(f"[WARN] {context}configured celltype_col='{celltype_col}' not "
                  f"found; using '{alias}' instead")
            return df

    print(f"[WARN] {context}no celltype column found "
          f"(looked for '{celltype_col}', 'celltype', 'cell_type'). "
          f"Available: {df.columns.tolist()}")
    return df


# =============================================================================
# STEP 1: GENERATE FRAGMENTS DICTIONARY
# =============================================================================

def generate_fragments_dict(cfg):
    print("\n" + "=" * 60)
    print("STEP 1: Generating fragments dictionary")
    print("=" * 60)

    fragments_dict = {}
    for category, prefixes in categories.items():
        for i in range(1, sample_counts + 1):
            folder = prefixes["folder_prefix"] + prefixes["folder_format"].format(i)
            sample = prefixes["sample_prefix"] + f"{i:02d}"
            frag_file = os.path.join(CELLRANGER_BASE_DIR, prefixes["folder_name"],
                                     folder, "outs", "atac_fragments.tsv.gz")
            if os.path.exists(frag_file):
                fragments_dict[sample] = frag_file
                print(f"[OK] Found: {sample}")
            else:
                print(f"[WARN] Missing: {sample} at {frag_file}")

    pd.to_pickle(fragments_dict, cfg.fragments_dict_path)
    print(f"\n Generated fragments_dict with {len(fragments_dict)} entries")
    return fragments_dict


# =============================================================================
# STEP 2: GET CHROMOSOME SIZES  (shared across compartments)
# =============================================================================

def get_chromsizes():
    print("\n" + "=" * 60)
    print("STEP 2: Getting chromosome sizes for mm10 (shared)")
    print("=" * 60)

    target_url = 'http://hgdownload.cse.ucsc.edu/goldenPath/mm10/bigZips/mm10.chrom.sizes'
    chromsizes = pd.read_csv(target_url, sep='\t', header=None)
    chromsizes.columns = ['Chromosome', 'End']
    chromsizes['Start'] = [0] * chromsizes.shape[0]
    chromsizes = chromsizes.loc[:, ['Chromosome', 'Start', 'End']]

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
    print(f"[OK] Loaded chromosome sizes: {len(chromsizes)} entries")
    return chromsizes


# =============================================================================
# STEP 3: LOAD AND PROCESS CELL DATA  (defines the compartment -> defines peaks)
# =============================================================================

def load_and_process_cell_data(cfg, fragments_dict):
    """Load compartment metadata, normalize barcodes/sample_id, filter samples.

    The frame returned here is what Step 4 turns into pseudobulks, so THIS is
    where the compartment restriction enters the peak set.
    """
    print("\n" + "=" * 60)
    print(f"STEP 3: Loading and processing cell_data ({cfg.name})")
    print("=" * 60)

    cell_data = pd.read_csv(cfg.cell_data_path, sep="\t", index_col=0)
    print(f"[STATS] Loaded {cfg.cell_data_path} with shape: {cell_data.shape}")
    print(f"   Columns: {cell_data.columns.tolist()}")

    cell_data = cell_data.reset_index()
    first_col = cell_data.columns[0]
    if first_col != 'barcode':
        cell_data.rename(columns={first_col: 'barcode'}, inplace=True)

    # Ensure a sample_id column (accept 'sample_id' or 'sample')
    if 'sample_id' not in cell_data.columns:
        if 'sample' in cell_data.columns:
            cell_data.rename(columns={'sample': 'sample_id'}, inplace=True)
        else:
            raise ValueError("[ERROR] cell_data needs a 'sample_id' or 'sample' column.")

    # Normalize the annotation column ('cell_type' -> 'celltype' for the atlas run)
    cell_data = normalize_celltype_col(cell_data, cfg.celltype_col)

    # Zero-pad sample IDs: geriatric_1 -> geriatric_01
    cell_data['sample_id'] = cell_data['sample_id'].astype(str).apply(
        lambda x: re.sub(r'_(\d)$', r'_0\1', x))
    print(f"[STATS] Unique sample_id in cell_data: "
          f"{sorted(cell_data['sample_id'].unique())}")

    before_drop = len(cell_data)
    cell_data = cell_data.dropna(subset=['sample_id'])
    print(f"[WARN] Dropped {before_drop - len(cell_data)} rows with NaN sample_id")

    # Normalize barcode to the 10x fragment form 'BARCODE-1'
    cell_data['barcode'] = cell_data['barcode'].apply(clean_to_fragment_barcode)
    cell_data = cell_data.set_index('barcode')
    cell_data['barcode'] = [str(x).split('___')[0] for x in cell_data.index.tolist()]

    # Diagnostic: sample_id overlap with fragments_dict
    cd_samples = set(cell_data['sample_id'].unique())
    frag_samples = set(fragments_dict.keys())
    print(f"[STATS] sample_id matched to fragments_dict: "
          f"{len(cd_samples & frag_samples)} / {len(cd_samples)}")
    only_cd = cd_samples - frag_samples
    if only_cd:
        print(f"[WARN] sample_id in cell_data but NOT in fragments_dict "
              f"(will be dropped): {sorted(only_cd)}")

    # Drop (sample, age) pairs with no cells
    pair_counts = cell_data.groupby(['sample_id', 'age']).size().reset_index(name='count')
    valid_pairs = pair_counts[pair_counts['count'] > 0][['sample_id', 'age']]
    before_merge = len(cell_data)
    cell_data = cell_data.merge(valid_pairs, on=['sample_id', 'age'], how='inner')
    print(f"[WARN] Cells dropped during valid_pairs merge: {before_merge - len(cell_data)}")

    # Keep only samples present in fragments_dict
    cell_data = cell_data[cell_data['sample_id'].isin(frag_samples)]
    valid_sample_ids = set(cell_data['sample_id'].unique())
    print(f"[OK] Final valid sample_ids: {len(valid_sample_ids)}")

    filtered_fragments_dict = {s: p for s, p in fragments_dict.items()
                               if s in valid_sample_ids}
    pd.to_pickle(filtered_fragments_dict, cfg.filtered_fragments_dict_path)
    print(f"[OK] Filtered fragments_dict includes {len(filtered_fragments_dict)} samples")

    if not filtered_fragments_dict:
        raise ValueError("[ERROR] Filtered fragments_dict is empty. Check sample_id naming.")

    print(f"[STATS] {cfg.name} cells retained: {len(cell_data)}")
    if 'celltype' in cell_data.columns:
        print(f"[STATS] {cfg.name} subtype counts:")
        print(cell_data['celltype'].value_counts())

    return cell_data, filtered_fragments_dict


# =============================================================================
# STEP 4: PSEUDOBULK EXPORT  (compartment cells only)
# =============================================================================

def run_pseudobulk_export(cfg, cell_data, filtered_fragments_dict, chromsizes):
    print("\n" + "=" * 60)
    print(f"STEP 4: Pseudobulk export (variable='{cfg.pseudobulk_variable}')")
    print("=" * 60)

    try:
        bw_paths, bed_paths = export_pseudobulk(
            input_data=cell_data,
            variable=cfg.pseudobulk_variable,
            sample_id_col='sample_id',
            chromsizes=chromsizes,
            bed_path=os.path.join(cfg.outDir, 'consensus_peak_calling',
                                  'pseudobulk_bed_files'),
            bigwig_path=os.path.join(cfg.outDir, 'consensus_peak_calling',
                                     'pseudobulk_bw_files'),
            path_to_fragments=filtered_fragments_dict,
            n_cpu=N_CPU,
            normalize_bigwig=True,
            temp_dir=os.path.join(cfg.tmpDir, 'ray_spill'),
            split_pattern=SPLIT_PATTERN_PSEUDOBULK
        )
        print("[OK] Pseudobulk export completed successfully.")
        save_paths_to_tsv(bw_paths, os.path.join(cfg.outDir, 'consensus_peak_calling',
                                                 'bw_paths.tsv'))
        save_paths_to_tsv(bed_paths, os.path.join(cfg.outDir, 'consensus_peak_calling',
                                                  'bed_paths.tsv'))
        return bw_paths, bed_paths
    except Exception as e:
        raise RuntimeError(f"[ERROR] Error during export_pseudobulk: {e}")


# =============================================================================
# STEP 5: PEAK CALLING
# =============================================================================

def run_peak_calling(cfg, bed_paths):
    print("\n" + "=" * 60)
    print("STEP 5: Peak calling with MACS2")
    print("=" * 60)

    macs_outdir = os.path.join(cfg.outDir, 'consensus_peak_calling', 'MACS')
    os.makedirs(macs_outdir, exist_ok=True)

    narrow_peaks_dict = peak_calling(
        'macs2', bed_paths, macs_outdir,
        genome_size='mm', n_cpu=N_CPU, input_format='BEDPE',
        shift=73, ext_size=146, keep_dup='all', q_value=0.05, _temp_dir=None
    )

    with open(os.path.join(macs_outdir, 'narrow_peaks_dict.pkl'), 'wb') as f:
        pickle.dump(narrow_peaks_dict, f)
    pd.to_pickle(narrow_peaks_dict, cfg.narrow_peak_dict_path)

    print(f"[OK] Peak calling completed. Found peaks for {len(narrow_peaks_dict)} groups.")
    return narrow_peaks_dict


# =============================================================================
# STEP 6: CONSENSUS PEAKS
# =============================================================================

def get_consensus_peaks_wrapper(cfg, narrow_peaks_dict, chromsizes):
    print("\n" + "=" * 60)
    print("STEP 6: Getting consensus peaks")
    print("=" * 60)

    fix_pyranges_head()
    narrow_peaks_dict = {str(k): v for k, v in narrow_peaks_dict.items()}

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
        print(f"[OK] Consensus peaks obtained: {len(consensus_peaks)} regions")
        consensus_peaks.to_bed(path=cfg.regions_bed, keep=True,
                               compression='infer', chain=False)
        print(f"[OK] Saved consensus peaks to {cfg.regions_bed}")
        return consensus_peaks
    except Exception as e:
        print(f"[ERROR] Error getting consensus peaks: {e}")
        raise


# =============================================================================
# STEP 7: GET TSS ANNOTATIONS
# =============================================================================

def get_tss_annotations(cfg):
    print("\n" + "=" * 60)
    print("STEP 7: Getting TSS annotations")
    print("=" * 60)
    try:
        subprocess.run([PYCISTOPIC_BIN, "tss", "get_tss",
                        "--output", cfg.tss_bed,
                        "--name", "mmusculus_gene_ensembl",
                        "--to-chrom-source", "ucsc",
                        "--ucsc", "mm10"], check=True)
        print(f"[OK] TSS annotations saved to {cfg.tss_bed}")
    except subprocess.CalledProcessError as e:
        print(f"[ERROR] Error running pycistopic tss get_tss: {e}")
        raise


# =============================================================================
# STEP 8: RUN QC
# =============================================================================

def run_qc(cfg, filtered_fragments_dict):
    print("\n" + "=" * 60)
    print("STEP 8: Running QC")
    print("=" * 60)

    with open(cfg.qc_commands_file, "w") as fh:
        for sample, fragment_filename in filtered_fragments_dict.items():
            sample_output = os.path.join(cfg.qc_dir, sample)
            print(f"{PYCISTOPIC_BIN} qc",
                  f"--fragments {fragment_filename}",
                  f"--regions {cfg.regions_bed}",
                  f"--tss {cfg.tss_bed}",
                  f"--output {sample_output}",
                  sep=" ", file=fh)

    os.chmod(cfg.qc_commands_file, 0o755)

    try:
        with open(cfg.qc_log_file, "w") as output_log:
            subprocess.run(f"cat {cfg.qc_commands_file} | parallel -j 4",
                           shell=True, check=True,
                           stdout=output_log, stderr=output_log)
        print("[OK] QC commands executed successfully.")
    except subprocess.CalledProcessError as e:
        print(f"[WARN] Error running pycistopic QC commands in parallel: {e}")
        print("   Attempting sequential execution...")
        with open(cfg.qc_commands_file, 'r') as f:
            for line in f:
                cmd = line.strip()
                if cmd:
                    try:
                        subprocess.run(cmd, shell=True, check=True)
                    except subprocess.CalledProcessError:
                        print(f"[WARN] Failed: {cmd}")

    time.sleep(5)

    print("\n Verifying QC outputs...")
    for sample_id in filtered_fragments_dict:
        p = os.path.join(cfg.qc_dir, f"{sample_id}.fragments_stats_per_cb.parquet")
        if not os.path.exists(p):
            print(f"[WARN] Missing metrics for {sample_id} - regenerating...")
            regenerate_qc_metrics(sample_id, filtered_fragments_dict[sample_id],
                                  cfg.regions_bed, cfg.tss_bed,
                                  os.path.join(cfg.qc_dir, sample_id))

    for sample_id in filtered_fragments_dict:
        p = os.path.join(cfg.qc_dir, f"{sample_id}.fragments_stats_per_cb.parquet")
        print(f"[ERROR] File still missing: {p}" if not os.path.exists(p)
              else f"[OK] Verified: {sample_id}")


# =============================================================================
# STEP 9: GET BARCODES PASSING QC
# =============================================================================

def get_barcodes_passing_qc(cfg, filtered_fragments_dict):
    print("\n" + "=" * 60)
    print("STEP 9: Getting barcodes passing QC")
    print("=" * 60)

    sample_id_to_barcodes_passing_filters = {}
    sample_id_to_thresholds = {}

    for sample_id in filtered_fragments_dict:
        try:
            barcodes, thresholds = get_barcodes_passing_qc_for_sample(
                sample_id=sample_id,
                pycistopic_qc_output_dir=cfg.qc_dir,
                unique_fragments_threshold=None,
                tss_enrichment_threshold=None,
                frip_threshold=0,
                use_automatic_thresholds=True,
            )
            sample_id_to_barcodes_passing_filters[sample_id] = barcodes
            sample_id_to_thresholds[sample_id] = thresholds
            print(f"[OK] {sample_id}: {len(barcodes)} barcodes passing QC")
        except FileNotFoundError as e:
            print(f"[WARN] File not found for {sample_id}: {e}")
        except Exception as e:
            print(f"[WARN] Error for {sample_id}: {e}")

    return sample_id_to_barcodes_passing_filters, sample_id_to_thresholds


# =============================================================================
# STEP 9B: QC VISUALIZATION
# =============================================================================

def generate_qc_plots(cfg, filtered_fragments_dict,
                      sample_id_to_barcodes_passing_filters, sample_id_to_thresholds):
    print("\n" + "=" * 60)
    print("STEP 9B: Generating QC visualization plots")
    print("=" * 60)

    os.makedirs(cfg.qc_plots_dir, exist_ok=True)

    # --- 1. sample-level summary ---
    print("\n[STATS] Generating sample-level summary statistics...")
    summary_stats = []
    for sample_id in filtered_fragments_dict:
        if sample_id in sample_id_to_barcodes_passing_filters:
            thresholds = sample_id_to_thresholds.get(sample_id, {})
            summary_stats.append({
                'sample_id': sample_id,
                'n_cells_passing_qc': len(sample_id_to_barcodes_passing_filters[sample_id]),
                'unique_fragments_threshold': thresholds.get('unique_fragments_threshold', 'N/A'),
                'tss_enrichment_threshold': thresholds.get('tss_enrichment_threshold', 'N/A'),
                'frip_threshold': thresholds.get('frip_threshold', 'N/A')
            })

    summary_df = pd.DataFrame(summary_stats)
    summary_df.to_csv(os.path.join(cfg.qc_plots_dir, 'qc_summary_stats.csv'), index=False)
    print("[OK] Saved QC summary stats")

    if not summary_df.empty:
        fig, ax = plt.subplots(figsize=(14, 6))
        bars = ax.bar(summary_df['sample_id'], summary_df['n_cells_passing_qc'],
                      color='steelblue', edgecolor='black', alpha=0.8)
        ax.set_xlabel('Sample ID', fontsize=12)
        ax.set_ylabel('Number of Cells Passing QC', fontsize=12)
        ax.set_title(f'{cfg.name} cells passing QC per sample',
                     fontsize=14, fontweight='bold')
        ax.tick_params(axis='x', rotation=45)
        plt.xticks(ha='right')
        for bar, val in zip(bars, summary_df['n_cells_passing_qc']):
            ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 50,
                    str(val), ha='center', va='bottom', fontsize=8)
        plt.tight_layout()
        fig.savefig(os.path.join(cfg.qc_plots_dir, 'cells_passing_qc_per_sample.png'),
                    dpi=150, bbox_inches='tight')
        fig.savefig(os.path.join(cfg.qc_plots_dir, 'cells_passing_qc_per_sample.pdf'),
                    bbox_inches='tight')
        plt.close(fig)
        print("[OK] Saved cells per sample plot")

    # --- 2. per-sample barcode statistics ---
    print("\n[STATS] Generating per-sample barcode statistics plots...")
    for sample_id in filtered_fragments_dict:
        try:
            p = os.path.join(cfg.qc_dir, f"{sample_id}.fragments_stats_per_cb.parquet")
            if not os.path.exists(p):
                print(f"[WARN] Skipping {sample_id} - metrics file not found")
                continue

            sample_metrics = pl.read_parquet(p).to_pandas()
            fig, axes = plt.subplots(2, 3, figsize=(15, 10))
            fig.suptitle(f'QC Metrics for {sample_id} ({cfg.name})',
                         fontsize=14, fontweight='bold')

            if 'total_fragments_count' in sample_metrics.columns:
                ax = axes[0, 0]
                ax.hist(np.log10(sample_metrics['total_fragments_count'] + 1),
                        bins=50, color='steelblue', edgecolor='black', alpha=0.7)
                ax.set_xlabel('log10(Total Fragments + 1)'); ax.set_ylabel('Count')
                ax.set_title('Total Fragments Distribution')
                ax.axvline(x=np.log10(1000), color='red', linestyle='--',
                           label='1000 threshold')
                ax.legend()

            if 'unique_fragments_count' in sample_metrics.columns:
                ax = axes[0, 1]
                ax.hist(np.log10(sample_metrics['unique_fragments_count'] + 1),
                        bins=50, color='forestgreen', edgecolor='black', alpha=0.7)
                ax.set_xlabel('log10(Unique Fragments + 1)'); ax.set_ylabel('Count')
                ax.set_title('Unique Fragments Distribution')

            if 'tss_enrichment' in sample_metrics.columns:
                ax = axes[0, 2]
                ax.hist(sample_metrics['tss_enrichment'],
                        bins=50, color='darkorange', edgecolor='black', alpha=0.7)
                ax.set_xlabel('TSS Enrichment Score'); ax.set_ylabel('Count')
                ax.set_title('TSS Enrichment Distribution')
                ax.axvline(x=2, color='red', linestyle='--', label='2.0 threshold')
                ax.legend()

            if 'fraction_of_fragments_in_peaks' in sample_metrics.columns:
                ax = axes[1, 0]
                ax.hist(sample_metrics['fraction_of_fragments_in_peaks'],
                        bins=50, color='purple', edgecolor='black', alpha=0.7)
                ax.set_xlabel('Fraction of Fragments in Peaks'); ax.set_ylabel('Count')
                ax.set_title('FRIP Distribution')

            if 'duplication_ratio' in sample_metrics.columns:
                ax = axes[1, 1]
                ax.hist(sample_metrics['duplication_ratio'],
                        bins=50, color='crimson', edgecolor='black', alpha=0.7)
                ax.set_xlabel('Duplication Ratio'); ax.set_ylabel('Count')
                ax.set_title('Duplication Ratio Distribution')

            if 'nucleosome_signal' in sample_metrics.columns:
                ax = axes[1, 2]
                ax.hist(sample_metrics['nucleosome_signal'],
                        bins=50, color='teal', edgecolor='black', alpha=0.7)
                ax.set_xlabel('Nucleosome Signal'); ax.set_ylabel('Count')
                ax.set_title('Nucleosome Signal Distribution')
            else:
                axes[1, 2].text(0.5, 0.5, 'Nucleosome Signal\nNot Available',
                                ha='center', va='center', fontsize=12)
                axes[1, 2].set_axis_off()

            plt.tight_layout()
            fig.savefig(os.path.join(cfg.qc_plots_dir, f'{sample_id}_qc_metrics.png'),
                        dpi=150, bbox_inches='tight')
            plt.close(fig)
            print(f"[OK] Saved QC plot for {sample_id}")
        except Exception as e:
            print(f"[WARN] Error generating QC plot for {sample_id}: {e}")

    # --- 3. aggregate across samples ---
    print("\n[STATS] Generating aggregate QC plots...")
    all_metrics = []
    for sample_id in filtered_fragments_dict:
        p = os.path.join(cfg.qc_dir, f"{sample_id}.fragments_stats_per_cb.parquet")
        if os.path.exists(p):
            df = pl.read_parquet(p).to_pandas()
            df['sample_id'] = sample_id
            df['age_group'] = sample_id.rsplit('_', 1)[0]
            all_metrics.append(df)

    if all_metrics:
        combined_metrics = pd.concat(all_metrics, ignore_index=True)
        age_groups = combined_metrics['age_group'].unique()
        positions = range(len(age_groups))
        colors = plt.cm.viridis(np.linspace(0, 1, len(age_groups)))

        if 'tss_enrichment' in combined_metrics.columns:
            fig, ax = plt.subplots(figsize=(12, 6))
            data_to_plot = [combined_metrics[combined_metrics['age_group'] == ag]
                            ['tss_enrichment'].dropna() for ag in age_groups]
            bp = ax.boxplot(data_to_plot, positions=positions, patch_artist=True)
            for patch, color in zip(bp['boxes'], colors):
                patch.set_facecolor(color); patch.set_alpha(0.7)
            ax.set_xticks(positions)
            ax.set_xticklabels(age_groups, rotation=45, ha='right')
            ax.set_xlabel('Age Group', fontsize=12)
            ax.set_ylabel('TSS Enrichment', fontsize=12)
            ax.set_title(f'TSS Enrichment by Age Group ({cfg.name})',
                         fontsize=14, fontweight='bold')
            ax.axhline(y=2, color='red', linestyle='--', alpha=0.7, label='Threshold')
            ax.legend()
            plt.tight_layout()
            fig.savefig(os.path.join(cfg.qc_plots_dir, 'tss_enrichment_by_age_group.png'),
                        dpi=150, bbox_inches='tight')
            fig.savefig(os.path.join(cfg.qc_plots_dir, 'tss_enrichment_by_age_group.pdf'),
                        bbox_inches='tight')
            plt.close(fig)
            print("[OK] Saved TSS enrichment by age group plot")

        if 'fraction_of_fragments_in_peaks' in combined_metrics.columns:
            fig, ax = plt.subplots(figsize=(12, 6))
            data_to_plot = [combined_metrics[combined_metrics['age_group'] == ag]
                            ['fraction_of_fragments_in_peaks'].dropna() for ag in age_groups]
            bp = ax.boxplot(data_to_plot, positions=positions, patch_artist=True)
            for patch, color in zip(bp['boxes'], colors):
                patch.set_facecolor(color); patch.set_alpha(0.7)
            ax.set_xticks(positions)
            ax.set_xticklabels(age_groups, rotation=45, ha='right')
            ax.set_xlabel('Age Group', fontsize=12)
            ax.set_ylabel('Fraction of Fragments in Peaks (FRIP)', fontsize=12)
            ax.set_title(f'FRIP by Age Group ({cfg.name})', fontsize=14, fontweight='bold')
            plt.tight_layout()
            fig.savefig(os.path.join(cfg.qc_plots_dir, 'frip_by_age_group.png'),
                        dpi=150, bbox_inches='tight')
            fig.savefig(os.path.join(cfg.qc_plots_dir, 'frip_by_age_group.pdf'),
                        bbox_inches='tight')
            plt.close(fig)
            print("[OK] Saved FRIP by age group plot")

        if ('tss_enrichment' in combined_metrics.columns
                and 'unique_fragments_count' in combined_metrics.columns):
            fig, ax = plt.subplots(figsize=(10, 8))
            plot_data = combined_metrics.sample(
                n=min(50000, len(combined_metrics)), random_state=42)
            scatter = ax.scatter(np.log10(plot_data['unique_fragments_count'] + 1),
                                 plot_data['tss_enrichment'],
                                 c=plot_data['age_group'].astype('category').cat.codes,
                                 cmap='viridis', alpha=0.3, s=5)
            ax.set_xlabel('log10(Unique Fragments + 1)', fontsize=12)
            ax.set_ylabel('TSS Enrichment', fontsize=12)
            ax.set_title(f'TSS Enrichment vs Unique Fragments ({cfg.name})',
                         fontsize=14, fontweight='bold')
            ax.axhline(y=2, color='red', linestyle='--', alpha=0.7, label='TSS threshold')
            ax.axvline(x=np.log10(1000), color='red', linestyle='--', alpha=0.7,
                       label='Fragments threshold')
            cbar = plt.colorbar(scatter, ax=ax)
            cbar.set_label('Age Group')
            plt.tight_layout()
            fig.savefig(os.path.join(cfg.qc_plots_dir, 'tss_vs_fragments_scatter.png'),
                        dpi=150, bbox_inches='tight')
            plt.close(fig)
            print("[OK] Saved TSS vs fragments scatter plot")

    print(f"\n[OK] All QC plots saved to {cfg.qc_plots_dir}")


# =============================================================================
# STEP 10: CREATE CISTOPIC OBJECTS (WITH PER-SAMPLE DOUBLET DETECTION)
# =============================================================================

def create_cistopic_objects(cfg, filtered_fragments_dict,
                            sample_id_to_barcodes_passing_filters):
    """Create cisTopic objects per sample.

    NOTE: valid_bc here is every barcode passing ATAC QC in that sample's
    fragments file, NOT only compartment cells. The compartment restriction is
    applied later, by the annotation-driven subset() in Step 13.
    """
    print("\n" + "=" * 60)
    print("STEP 10: Creating cisTopic objects (with per-sample doublet detection)")
    print("=" * 60)

    cistopic_obj_list = []
    total_doublets = 0
    total_cells_before = 0

    for sample_id in filtered_fragments_dict:
        obj_file = os.path.join(cfg.individual_objs_dir, f"{sample_id}_cistopic_obj.pkl")
        if os.path.exists(obj_file):
            with open(obj_file, "rb") as f:
                cistopic_obj = pickle.load(f)
            print(f" Loaded existing CistopicObject for {sample_id}")
            cistopic_obj_list.append(cistopic_obj)
            continue

        p = os.path.join(cfg.qc_dir, f"{sample_id}.fragments_stats_per_cb.parquet")
        if not os.path.exists(p):
            print(f"[WARN] Metrics file missing for {sample_id}: {p}")
            continue
        if sample_id not in sample_id_to_barcodes_passing_filters:
            print(f"[WARN] No barcodes passing QC for {sample_id}. Skipping...")
            continue

        sample_metrics = pl.read_parquet(p).to_pandas()
        valid_barcodes = sample_id_to_barcodes_passing_filters[sample_id]
        sample_metrics = sample_metrics.set_index("CB").loc[valid_barcodes]

        print(f"\n Creating CistopicObject for {sample_id}...")
        cistopic_obj = create_cistopic_object_from_fragments(
            path_to_fragments=filtered_fragments_dict[sample_id],
            path_to_regions=cfg.regions_bed,
            path_to_blacklist=path_to_blacklist,
            metrics=sample_metrics,
            valid_bc=valid_barcodes,
            n_cpu=N_CPU,
            project=sample_id,
            split_pattern=SPLIT_PATTERN_CISTOPIC
        )

        n_cells_before = len(cistopic_obj.cell_names)
        total_cells_before += n_cells_before

        cistopic_obj, n_doublets = run_scrublet_per_sample(
            cistopic_obj, expected_doublet_rate=DOUBLET_RATE,
            threshold=DOUBLET_THRESHOLD)

        total_doublets += n_doublets
        print(f"   [STATS] {sample_id}: {n_cells_before} cells -> "
              f"{len(cistopic_obj.cell_names)} singlets ({n_doublets} doublets removed)")

        cistopic_obj_list.append(cistopic_obj)
        save_cistopic_object(cistopic_obj, sample_id, cfg.individual_objs_dir)

    print(f"\n Summary: {total_cells_before} total cells, {total_doublets} doublets removed")
    return cistopic_obj_list


# =============================================================================
# STEP 11: MERGE CISTOPIC OBJECTS
# =============================================================================

def merge_cistopic_objects(cfg, cistopic_obj_list):
    print("\n" + "=" * 60)
    print("STEP 11: Merging cisTopic objects")
    print("=" * 60)

    if len(cistopic_obj_list) == 0:
        raise ValueError("[ERROR] No cisTopic objects to merge!")

    cistopic_obj = merge(cistopic_obj_list)
    print(f"[OK] Merged {len(cistopic_obj_list)} cisTopic objects")
    print(f"   [STATS] Final shape: {cistopic_obj.fragment_matrix.shape}")

    nan_count = cistopic_obj.cell_data.isna().sum().sum()
    print(f"   [WARN] Found {nan_count} total NaN entries in cell metadata.")
    if nan_count > 0:
        print("    Dropping rows with NaNs from cell_data.")
        cistopic_obj.cell_data = cistopic_obj.cell_data.dropna()

    with open(cfg.merged_path, 'wb') as f:
        pickle.dump(cistopic_obj, f)
    print(f"[OK] Saved merged cisTopic object to {cfg.merged_path}")

    return cistopic_obj


# =============================================================================
# STEP 12: ADD CELL TYPE ANNOTATIONS  (robust barcode reconstruction)
# =============================================================================

def add_cell_annotations(cfg, cistopic_obj, annotation_file_path=None):
    """Add compartment annotations to the cisTopic object.

    Builds the cisTopic-style index 'BARCODE-1-sample___sample' from the DNA
    barcode + sample_id columns (robust to BARCODE-1-sample or BARCODE-sample
    obs_names), rather than a fragile single-regex extraction.
    """
    print("\n" + "=" * 60)
    print(f"STEP 12: Adding {cfg.name} annotations")
    print("=" * 60)

    if annotation_file_path is None:
        annotation_file_path = cfg.annotation_path

    if not os.path.exists(annotation_file_path):
        print(f"[WARN] Annotation file not found: {annotation_file_path}. Skipping.")
        return cistopic_obj

    cell_data = pd.read_csv(annotation_file_path, sep="\t", index_col=0)
    print(f"[STATS] Loaded annotation data with shape: {cell_data.shape}")

    cell_data = cell_data.reset_index()
    first_col = cell_data.columns[0]
    if first_col != 'barcode':
        cell_data.rename(columns={first_col: 'barcode'}, inplace=True)

    # DNA barcode = segment before the first '-'
    cell_data['dna_barcode'] = cell_data['barcode'].astype(str).str.split('-').str[0]

    if 'sample_id' not in cell_data.columns:
        if 'sample' in cell_data.columns:
            cell_data.rename(columns={'sample': 'sample_id'}, inplace=True)
        else:
            print("[ERROR] No sample_id/sample column in annotation file.")
            return cistopic_obj
    cell_data['sample_id'] = cell_data['sample_id'].astype(str).apply(
        lambda x: re.sub(r'_(\d)$', r'_0\1', x))

    # Normalize the annotation column so keep_cols/Step 13 always see 'celltype'
    cell_data = normalize_celltype_col(cell_data, cfg.celltype_col)

    # cisTopic-style index: BARCODE-1-sample___sample
    cell_data['cisTopic_id'] = (cell_data['dna_barcode'] + '-1-'
                                + cell_data['sample_id'] + '___'
                                + cell_data['sample_id'])

    before_drop = len(cell_data)
    cell_data = cell_data.dropna(subset=['cisTopic_id'])
    print(f"[WARN] Dropped {before_drop - len(cell_data)} rows with invalid index")

    cell_data = cell_data.set_index('cisTopic_id')

    keep_cols = [c for c in ['sample_id', 'celltype', 'age', 'sex']
                 if c in cell_data.columns]
    cell_data = cell_data[keep_cols]
    print(f"[OK] Annotation columns kept: {keep_cols}")
    print(f"[OK] Example transformed index: {cell_data.index[:3].tolist()}")

    cis_barcodes = set(cistopic_obj.cell_names)
    transformed = set(cell_data.index)
    intersect = cis_barcodes & transformed

    print("\n[STATS] Barcode matching statistics:")
    print(f"   cisTopic cells: {len(cis_barcodes)}")
    print(f"   Annotation cells: {len(transformed)}")
    print(f"   [OK] Matched: {len(intersect)}")
    print(f"   [INFO] In annotation but not cisTopic: {len(transformed - cis_barcodes)}")
    print(f"   [INFO] In cisTopic but not annotation: {len(cis_barcodes - transformed)}")

    if len(intersect) == 0:
        print("[ERROR] No matching barcodes found! Check obs_name / cisTopic naming.")
        print("        Sample cisTopic cell_names:", list(cis_barcodes)[:3])
        return cistopic_obj

    if cell_data.index.duplicated().any():
        n = cell_data.index.duplicated().sum()
        print(f"[WARN] {n} duplicate indices in annotation; keeping first")
        cell_data = cell_data[~cell_data.index.duplicated(keep='first')]
    if cistopic_obj.cell_data.index.duplicated().any():
        n = cistopic_obj.cell_data.index.duplicated().sum()
        print(f"[WARN] {n} duplicate indices in cisTopic; keeping first")
        cistopic_obj.cell_data = cistopic_obj.cell_data[
            ~cistopic_obj.cell_data.index.duplicated(keep='first')]

    try:
        cistopic_obj.add_cell_data(cell_data)
        print("[OK] Successfully added annotation data to cisTopic object")
    except Exception as e:
        print(f"[ERROR] Error adding cell data: {e}")
        return cistopic_obj

    # IMPORTANT: do NOT dropna cell_data here. The object still holds every
    # QC-passing cell; trimming cell_data now desyncs it from cell_names and makes
    # the Step 13 subset fail (IndexError: positional indexers out-of-bounds).
    # Non-compartment cells (NaN annotation) are removed properly -- together with
    # cell_names and the fragment matrix -- by cistopic_obj.subset() in Step 13.
    n_annotated = (cistopic_obj.cell_data.dropna(subset=keep_cols).shape[0]
                   if keep_cols else 0)
    print(f"[OK] cell_data kept aligned to object: {cistopic_obj.cell_data.shape} "
          f"({n_annotated} cells carry full compartment annotation)")
    return cistopic_obj


# =============================================================================
# STEP 13: FILTER AND FINALIZE  (the actual compartment trim of the object)
# =============================================================================

def filter_and_finalize(cfg, cistopic_obj, required_columns=None):
    print("\n" + "=" * 60)
    print("STEP 13: Filtering and finalizing")
    print("=" * 60)

    if required_columns is None:
        required_columns = ['sample_id', 'celltype', 'age', 'sex']

    existing = [c for c in required_columns if c in cistopic_obj.cell_data.columns]
    missing = [c for c in required_columns if c not in cistopic_obj.cell_data.columns]
    if missing:
        print(f"[WARN] Missing expected columns: {missing}")

    if existing:
        print(f"[INFO] Restricting object to cells with complete annotation: {existing}")
        before_n = len(cistopic_obj.cell_names)
        # Cells belonging to this compartment carry non-NaN celltype/age/sex after
        # Step 12; everything else entered via ATAC QC and must be removed.
        cells_to_keep = cistopic_obj.cell_data.dropna(subset=existing).index.tolist()
        print(f"   Cells in merged object (all QC-passing): {before_n}")
        print(f"   Cells with compartment annotation:       {len(cells_to_keep)}")
        if len(cells_to_keep) == 0:
            raise ValueError("[ERROR] No annotated cells to keep -- check Step 12 matching.")
        # CRITICAL: subset the OBJECT (fragment_matrix + cell_names + cell_data),
        # not just the metadata table.
        cistopic_obj = cistopic_obj.subset(cells_to_keep, copy=True,
                                           split_pattern=SPLIT_PATTERN_CISTOPIC)
        print(f"   Cells after subset (compartment only):   {len(cistopic_obj.cell_names)}")

    print("\n Cell counts per sample:")
    try:
        sample_names = [x.split('-')[2].split('___')[0] for x in cistopic_obj.cell_names]
        counts = Counter(sample_names)
        for sample, count in sorted(counts.items()):
            print(f"   {sample}: {count} cells")
        print(f"\n[STATS] Total unique samples: {len(counts)}")
        print(f"[STATS] Total cells: {sum(counts.values())}")
    except Exception as e:
        print(f"[WARN] Could not parse sample names from cell_names: {e}")

    print(f"\n Any remaining NaN values: {cistopic_obj.cell_data.isna().any().any()}")
    return cistopic_obj


# =============================================================================
# STEP 14: SAVE FINAL OBJECT
# =============================================================================

def save_final_object(cfg, cistopic_obj):
    print("\n" + "=" * 60)
    print("STEP 14: Saving final object")
    print("=" * 60)

    output_path = os.path.join(cfg.outDir,
                               f'cisTopicObject_filtered_annotated_{cfg.suffix}.pkl')
    with open(output_path, 'wb') as f:
        pickle.dump(cistopic_obj, f)
    print(f"[OK] Saved final cisTopic object to {output_path}")
    return output_path


# =============================================================================
# STEP 15: FINAL SUMMARY
# =============================================================================

def print_final_summary(cfg, cistopic_obj):
    print("\n" + "=" * 60)
    print(f"[STATS] FINAL OBJECT SUMMARY FOR SCENIC+ ({cfg.name})")
    print("=" * 60)

    print(f"Fragment matrix shape: {cistopic_obj.fragment_matrix.shape}")
    print(f"Number of cells: {len(cistopic_obj.cell_names)}")
    print(f"Number of regions: {len(cistopic_obj.region_names)}")
    print(f"Cell metadata columns: {cistopic_obj.cell_data.columns.tolist()}")

    if 'sample_id' in cistopic_obj.cell_data.columns:
        print(f"Samples represented: {cistopic_obj.cell_data['sample_id'].nunique()}")
        summary = cistopic_obj.cell_data['sample_id'].value_counts().sort_index()
        print("\n Cell counts per sample (final):")
        print(summary)
        summary.to_csv(os.path.join(cfg.outDir, "final_cell_counts_per_sample.csv"),
                       header=["cell_count"])

    if 'celltype' in cistopic_obj.cell_data.columns:
        print(f"\n {cfg.name} subtype distribution:")
        print(cistopic_obj.cell_data['celltype'].value_counts())

    if 'age' in cistopic_obj.cell_data.columns:
        print("\n Age group distribution:")
        print(cistopic_obj.cell_data['age'].value_counts())

    print(f"\n[OK] Final cells with complete metadata: "
          f"{cistopic_obj.cell_data.dropna().shape[0]}")

    required_cols = ['sample_id', 'celltype', 'age', 'sex']
    missing = [c for c in required_cols if c not in cistopic_obj.cell_data.columns]
    if missing:
        print(f"[WARN] Missing recommended columns for SCENIC+: {missing}")
    else:
        print("[OK] All recommended metadata columns present for SCENIC+")


# =============================================================================
# PER-COMPARTMENT DRIVER
# =============================================================================

def run_compartment(cfg, chromsizes, run_preprocessing=True, run_annotation=True,
                    generate_plots=True, skip_qc=False):
    """Run the full pipeline for ONE compartment (its own peaks, own outputs)."""
    print("\n" + "#" * 70)
    print(f"#  COMPARTMENT: {cfg.name}   ->  {cfg.outDir}")
    print("#" * 70)

    cfg.makedirs()
    cistopic_obj = None

    if run_preprocessing:
        # Step 1
        if os.path.exists(cfg.fragments_dict_path):
            print(f" Loading existing fragments_dict from {cfg.fragments_dict_path}")
            fragments_dict = pd.read_pickle(cfg.fragments_dict_path)
        else:
            fragments_dict = generate_fragments_dict(cfg)

        # Step 3 — compartment restriction enters here
        cell_data, filtered_fragments_dict = load_and_process_cell_data(cfg, fragments_dict)

        # Step 4
        bed_paths_file = os.path.join(cfg.outDir, 'consensus_peak_calling', 'bed_paths.tsv')
        if os.path.exists(bed_paths_file):
            print(f" Loading existing bed paths from {bed_paths_file}")
            bed_paths = load_paths_from_tsv(bed_paths_file)
        else:
            _, bed_paths = run_pseudobulk_export(cfg, cell_data,
                                                 filtered_fragments_dict, chromsizes)

        # Step 5
        narrow_peaks_file = os.path.join(cfg.outDir, 'consensus_peak_calling',
                                         'MACS', 'narrow_peaks_dict.pkl')
        if os.path.exists(narrow_peaks_file):
            print(f" Loading existing narrow peaks from {narrow_peaks_file}")
            with open(narrow_peaks_file, 'rb') as f:
                narrow_peaks_dict = pickle.load(f)
        else:
            narrow_peaks_dict = run_peak_calling(cfg, bed_paths)

        # Step 6
        if os.path.exists(cfg.regions_bed):
            print(f" Consensus regions already exist at {cfg.regions_bed}")
        else:
            get_consensus_peaks_wrapper(cfg, narrow_peaks_dict, chromsizes)

        # Step 7
        if os.path.exists(cfg.tss_bed):
            print(f" TSS annotations already exist at {cfg.tss_bed}")
        else:
            get_tss_annotations(cfg)

        # Step 8
        if skip_qc:
            print(f"[INFO] --skip-qc set: reusing existing QC parquets in {cfg.qc_dir}")
        else:
            run_qc(cfg, filtered_fragments_dict)

        # Step 9
        sample_id_to_barcodes, sample_id_to_thresholds = get_barcodes_passing_qc(
            cfg, filtered_fragments_dict)

        # Step 9B
        if generate_plots:
            generate_qc_plots(cfg, filtered_fragments_dict,
                              sample_id_to_barcodes, sample_id_to_thresholds)

        # Step 10
        cistopic_obj_list = create_cistopic_objects(cfg, filtered_fragments_dict,
                                                    sample_id_to_barcodes)

        # Step 11
        cistopic_obj = merge_cistopic_objects(cfg, cistopic_obj_list)

    if run_annotation:
        if cistopic_obj is None:
            print(f" Loading merged cisTopic object from {cfg.merged_path}")
            with open(cfg.merged_path, 'rb') as infile:
                cistopic_obj = pickle.load(infile)
            print("[OK] cisTopic object loaded successfully")
            print(f"   meta_data shape: {cistopic_obj.cell_data.shape}")

        cistopic_obj = add_cell_annotations(cfg, cistopic_obj)   # Step 12
        cistopic_obj = filter_and_finalize(cfg, cistopic_obj)    # Step 13
        save_final_object(cfg, cistopic_obj)                     # Step 14

    if cistopic_obj is not None:
        print_final_summary(cfg, cistopic_obj)                   # Step 15

    print("\n" + "=" * 60)
    print(f"[OK] {cfg.name} PIPELINE COMPLETED")
    print("=" * 60)
    return cistopic_obj


# =============================================================================
# MAIN
# =============================================================================

def main(run=None, run_preprocessing=True, run_annotation=True,
         generate_plots=True, skip_qc=False, continue_on_error=True):
    print("\n" + "=" * 60)
    print("pycisTopic PREPROCESSING PIPELINE FOR SCENIC+  (ALL COMPARTMENTS)")
    print("=" * 60)

    to_run = run or list(COMPARTMENTS)
    print(f"Compartments to run: {to_run}")

    # Step 2 is the only genuinely shared step -> fetch once
    chromsizes = get_chromsizes() if run_preprocessing else None

    results = {}
    for name in to_run:
        spec = COMPARTMENTS[name]
        cfg = Cfg(name, spec["out_dir"], spec["cell_data"], spec["suffix"],
                  celltype_col=spec.get("celltype_col", "celltype"),
                  pseudobulk_variable=spec.get("pseudobulk_variable"))
        try:
            results[name] = run_compartment(
                cfg, chromsizes,
                run_preprocessing=run_preprocessing,
                run_annotation=run_annotation,
                generate_plots=generate_plots,
                skip_qc=skip_qc,
            )
        except Exception as e:
            print(f"\n[ERROR] {name} failed: {e}")
            traceback.print_exc()
            if not continue_on_error:
                raise
            print(f"[INFO] Continuing to next compartment...\n")
            results[name] = None

    print("\n" + "=" * 60)
    print("ALL COMPARTMENTS FINISHED")
    for name, obj in results.items():
        status = (f"{len(obj.cell_names)} cells" if obj is not None else "FAILED")
        print(f"   {name}: {status}")
    print("=" * 60)
    return results


# =============================================================================
# ENTRY POINT
# =============================================================================

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description='pycisTopic preprocessing pipeline for SCENIC+ (all compartments)')
    parser.add_argument('--run', nargs='+', default=None, choices=sorted(COMPARTMENTS),
                        help='Which compartments to run (default: all)')
    parser.add_argument('--skip-preprocessing', action='store_true',
                        help='Skip steps 1-11 and only run annotation')
    parser.add_argument('--skip-annotation', action='store_true',
                        help='Skip annotation steps (12-14)')
    parser.add_argument('--skip-plots', action='store_true',
                        help='Skip QC visualization plots')
    parser.add_argument('--skip-qc', action='store_true',
                        help='Skip Step 8 QC (reuse existing qc/*.parquet files)')
    parser.add_argument('--stop-on-error', action='store_true',
                        help='Abort on the first compartment failure '
                             '(default: continue to the next)')

    args = parser.parse_args()

    main(
        run=args.run,
        run_preprocessing=not args.skip_preprocessing,
        run_annotation=not args.skip_annotation,
        generate_plots=not args.skip_plots,
        skip_qc=args.skip_qc,
        continue_on_error=not args.stop_on_error,
    )
