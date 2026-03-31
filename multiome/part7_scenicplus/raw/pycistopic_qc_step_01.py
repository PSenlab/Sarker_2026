# === Standard Library ===
import os
import re
import json
import pickle
import subprocess

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
from pycisTopic.cistopic_class import create_cistopic_object_from_fragments
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


# Set directories
outDir = "outs_trial"
tmpDir = os.path.join(outDir, "temp_sbatch")

# Ensure directories exist
os.makedirs(tmpDir, exist_ok=True)
os.makedirs(os.path.join(outDir, 'consensus_peak_calling', 'pseudobulk_bed_files'), exist_ok=True)
os.makedirs(os.path.join(outDir, 'consensus_peak_calling', 'pseudobulk_bw_files'), exist_ok=True)

# Define paths for intermediate results
fragments_dict_path = os.path.join(tmpDir, 'fragments_dict.pkl')
narrow_peak_dict_path = os.path.join(tmpDir, 'narrow_peak_dict.pkl')


# Define category mappings
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

# Initialize the fragments dictionary
fragments_dict = {}

# Generate fragments_dict
for category, prefixes in categories.items():
    folder_prefix = prefixes["folder_prefix"]
    sample_prefix = prefixes["sample_prefix"]
    folder_format = prefixes["folder_format"]
    folder_name = prefixes["folder_name"]

    for i in range(1, sample_counts + 1):
        # Construct folder and sample names
        folder = folder_prefix + folder_format.format(i)
        sample = sample_prefix + f"{i:02d}"

        # Construct fragment file path
        frag_file = os.path.join(
            "/data/sarkern2/multiome_liver/CellRanger",
            folder_name,
            folder,
            "outs",
            "atac_fragments.tsv.gz"
        )

        print(f"Category: {category}, folder: {folder}, sample: {sample}")
        print(f"Checking path: {frag_file}")
        if os.path.exists(frag_file):
            print("Exists: True")
            fragments_dict[sample] = frag_file
        else:
            print("Exists: False")

# Save fragments_dict
pd.to_pickle(fragments_dict, fragments_dict_path)
print(f"Generated fragments_dict with {len(fragments_dict)} entries:")
for sample, path in fragments_dict.items():
    print(f"{sample}: {path}")


# Get chromosome sizes for mm10
target_url = 'http://hgdownload.cse.ucsc.edu/goldenPath/mm10/bigZips/mm10.chrom.sizes'
chromsizes = pd.read_csv(target_url, sep='\t', header=None)
chromsizes.columns = ['Chromosome', 'End']
chromsizes['Start'] = [0]*chromsizes.shape[0]
chromsizes = chromsizes.loc[:, ['Chromosome', 'Start', 'End']]
# Adjust Chromosome names to match CellRangerARC annotations
chromsizes['Chromosome'] = [chromsizes['Chromosome'][x].replace('v', '.') for x in range(len(chromsizes['Chromosome']))]
chromsizes['Chromosome'] = [
    chromsizes['Chromosome'][x].split('_')[1] if len(chromsizes['Chromosome'][x].split('_')) > 1
    else chromsizes['Chromosome'][x] for x in range(len(chromsizes['Chromosome']))
]
chromsizes = pr.PyRanges(chromsizes)


import re
import pandas as pd

# Load cell data
cell_data_path = "cell_data.tsv"
cell_data = pd.read_csv(cell_data_path, sep="\t", index_col=0)
print(cell_data.head())
print("Columns in cell_data:", cell_data.columns)

# Move 'sample' from the index to a regular column
cell_data = cell_data.reset_index()
print("Columns before rename:", cell_data.columns)

# Rename 'sample' to 'sample_id'
cell_data.rename(columns={"sample": "sample_id"}, inplace=True)

# Zero-pad sample IDs: e.g., geriatric_1 -> geriatric_01
cell_data['sample_id'] = cell_data['sample_id'].apply(
    lambda x: re.sub(r'_(\d)$', r'_0\1', x)
)
print("Unique sample_id in cell_data:", cell_data['sample_id'].unique())


# Drop rows with NaN values in the 'sample_id' column after mapping
cell_data = cell_data.dropna(subset=['sample_id'])


# Clean the barcode by removing extra suffixes:
cell_data['barcode'] = cell_data['barcode'].astype(str)
cell_data['barcode'] = cell_data['barcode'].apply(lambda x: re.sub(r'-.*$', '-1', x.split('___')[0]))



# Set 'barcode' as the new index
cell_data = cell_data.set_index('barcode')
print("New index is barcode:\n", cell_data.index[:5])

# Add barcode column (if needed for downstream analysis)
cell_data['barcode'] = [str(x).split('___')[0] for x in cell_data.index.tolist()]

# Filter out (sample, age) pairs with no cells
pair_counts = cell_data.groupby(['sample_id', 'age']).size().reset_index(name='count')
valid_pairs = pair_counts[pair_counts['count'] > 0][['sample_id', 'age']]
cell_data = cell_data.merge(valid_pairs, on=['sample_id', 'age'], how='inner')



# Load fragments_dict
fragments_dict = pd.read_pickle(fragments_dict_path)


# Recompute (sample_id, age) pairs in cell_data after all filtering
filtered_pairs = cell_data.groupby(['sample_id', 'age']).size().reset_index(name='count')
filtered_pairs = filtered_pairs[filtered_pairs['count'] > 0]

# Convert to set of sample_ids with valid barcodes and valid age assignment
valid_sample_ids = set(filtered_pairs['sample_id'].unique())

# Final filter: retain only rows with valid sample_ids
cell_data = cell_data[cell_data['sample_id'].isin(valid_sample_ids)]
print("Final filtered sample_id values:", sorted(cell_data['sample_id'].unique()))



# Build a fragments_dict that only includes filtered sample_ids
filtered_fragments_dict = {
    sample: path
    for sample, path in fragments_dict.items()
    if sample in valid_sample_ids
}

# Save to disk
filtered_fragments_dict_path = "filtered_fragments_dict.pkl"
pd.to_pickle(filtered_fragments_dict, filtered_fragments_dict_path)

print(f"✅ Final filtered_fragments_dict includes {len(filtered_fragments_dict)} samples:")
for s, p in filtered_fragments_dict.items():
    print(f"  {s}: {p}")




# Ensure the filtered fragments_dict has entries
if not filtered_fragments_dict:
    raise ValueError("Filtered fragments_dict is empty. Ensure proper filtering of fragments.")



# Attempt pseudobulk export
try:
    bw_paths, bed_paths = export_pseudobulk(
        input_data=cell_data,
        variable='age',
        sample_id_col='sample_id',
        chromsizes=chromsizes,
        bed_path=os.path.join(outDir, 'consensus_peak_calling', 'pseudobulk_bed_files'),
        bigwig_path=os.path.join(outDir, 'consensus_peak_calling', 'pseudobulk_bw_files'),
        path_to_fragments=filtered_fragments_dict,  # Use filtered_fragments_dict here
        n_cpu=15,
        normalize_bigwig=True,
        temp_dir=os.path.join(tmpDir, 'ray_spill'),
        split_pattern='___'
    )
    print("Pseudobulk export completed successfully.")
except Exception as e:
    raise RuntimeError(f"Error during export_pseudobulk: {e}")
# Save paths to TSV files
def save_paths_to_tsv(paths_dict, filepath):
    with open(filepath, 'w') as f:
        for key, path in paths_dict.items():
            f.write(f"{key}\t{path}\n")

save_paths_to_tsv(bw_paths, os.path.join(outDir, 'consensus_peak_calling', 'bw_paths.tsv'))
save_paths_to_tsv(bed_paths, os.path.join(outDir, 'consensus_peak_calling', 'bed_paths.tsv'))



# Load paths from TSV files
def load_paths_from_tsv(filepath):
    paths_dict = {}
    with open(filepath, 'r') as f:
        for line in f:
            v, p = line.strip().split('\t')
            paths_dict[v] = p
    return paths_dict

bw_paths = load_paths_from_tsv(os.path.join(outDir, 'consensus_peak_calling', 'bw_paths.tsv'))
bed_paths = load_paths_from_tsv(os.path.join(outDir, 'consensus_peak_calling', 'bed_paths.tsv'))




# Step2: Inferring consensus peaks
from pycisTopic.pseudobulk_peak_calling import *


macs_path = 'macs2'
macs_outdir = os.path.join(outDir, 'consensus_peak_calling', 'MACS')
os.makedirs(macs_outdir, exist_ok=True)
peak_half_width = 250
path_to_blacklist = "/data/sarkern2/scenicplus1/blacklist/mm10-blacklist.v2.bed"  # Update to mouse blacklist


narrow_peaks_dict = peak_calling(
    macs_path,
    bed_paths,
    macs_outdir,
    genome_size='mm',
    n_cpu=15,
    input_format='BEDPE',
    shift=73,
    ext_size=146,
    keep_dup='all',
    q_value=0.05,
    _temp_dir=None
)

with open(os.path.join(macs_outdir, 'narrow_peaks_dict.pkl'), 'wb') as f:
    pickle.dump(narrow_peaks_dict, f)
pd.to_pickle(narrow_peaks_dict, narrow_peak_dict_path)


# Load the narrow_peaks_dict.pkl file
with open('outs_trial/consensus_peak_calling/MACS/narrow_peaks_dict.pkl', 'rb') as f:
    narrow_peaks_dict = pickle.load(f)


# Override the head method to use np.bool_
def head_override(self, n=8):
    """Return the n first rows."""
    subsetter = np.zeros(len(self), dtype=np.bool_)  # Use np.bool_ instead of np.bool
    subsetter[:n] = True
    return self[subsetter]

# Apply the override to the PyRanges class
pr.PyRanges.head = head_override


# Convert keys to strings if they are integers
narrow_peaks_dict = {str(k): v for k, v in narrow_peaks_dict.items()}

# Ensure each PyRanges object has the required columns
required_columns = ['Chromosome', 'Start', 'End', 'Summit', 'Name', 'Score']

for key, pyrange in narrow_peaks_dict.items():
    df = pyrange.df
    
    # Add any missing required columns with default values
    for col in required_columns:
        if col not in df.columns:
            df[col] = 0  # or another appropriate default value
    
    # Convert DataFrame back to PyRanges
    narrow_peaks_dict[key] = pr.PyRanges(df)

# Verify the structure again using the overridden method
for key, pyrange in narrow_peaks_dict.items():
    print(f"Key: {key}")
    print(pyrange.head())

from pycisTopic.iterative_peak_calling import get_consensus_peaks


# Get consensus peaks
try:
    consensus_peaks = get_consensus_peaks(
        narrow_peaks_dict=narrow_peaks_dict,
        peak_half_width=peak_half_width,
        chromsizes=chromsizes,
        path_to_blacklist=path_to_blacklist
    )
    print("Consensus peaks obtained successfully.")
except AttributeError as e:
    print(f"An error occurred: {e}")


consensus_peaks.to_bed(
    path = os.path.join(outDir, "consensus_peak_calling/consensus_regions.bed"),
    keep =True,
    compression = 'infer',
    chain = False)


# Step 2: QC

# Create QC directory if it doesn't exist
os.makedirs(os.path.join(outDir, "qc"), exist_ok=True)
out_dir = outDir

# Get TSS annotations

# Run pycistopic TSS command to get TSS annotations for mouse
try:
    subprocess.run(
        [
            "/data/sarkern2/conda/envs/scenicplus3/bin/pycistopic", "tss", "get_tss",
            "--output", os.path.join(outDir, "qc", "tss.bed"),
            "--name", "mmusculus_gene_ensembl",
            "--to-chrom-source", "ucsc",
            "--ucsc", "mm10"
        ],
        check=True
    )
    print("TSS annotations for mouse obtained successfully.")
except subprocess.CalledProcessError as e:
    print(f"Error running pycistopic tss get_tss: {e}")

# Preview the TSS file
try:
    subprocess.run(
        ["head", os.path.join(outDir, "qc", "tss.bed")],
        check=True
    )
except subprocess.CalledProcessError as e:
    print(f"Error previewing TSS file: {e}")

# Define paths for QC
regions_bed_filename = os.path.join(outDir, "consensus_peak_calling", "consensus_regions.bed")
tss_bed_filename = os.path.join(outDir, "qc", "tss.bed")
pycistopic_qc_commands_filename = "pycistopic_qc_commands.txt"

# Create text file with all pycistopic QC command lines
with open(pycistopic_qc_commands_filename, "w") as fh:
    for sample, fragment_filename in filtered_fragments_dict.items():
        print(
            "/data/sarkern2/conda/envs/scenicplus3/bin/pycistopic qc",
            f"--fragments {fragment_filename}",
            f"--regions {regions_bed_filename}",
            f"--tss {tss_bed_filename}",
            f"--output {os.path.join(outDir, 'qc', sample)}",
            sep=" ",
            file=fh,
        )

# Ensure the command script is executable
os.chmod(pycistopic_qc_commands_filename, 0o755)



# Run the QC commands in parallel
try:
    with open("qc_commands_output.log", "w") as output_log:
        subprocess.run(["cat", pycistopic_qc_commands_filename, "|", "parallel", "-j", "4", "{}"], shell=True, check=True, stdout=output_log, stderr=output_log)
    print("QC commands executed successfully.")
except subprocess.CalledProcessError as e:
    print(f"Error running pycistopic QC commands in parallel: {e}")


# === Get barcodes passing QC for each sample ===
sample_id_to_barcodes_passing_filters = {}
sample_id_to_thresholds = {}

for sample_id in filtered_fragments_dict:
    try:
        barcodes, thresholds = get_barcodes_passing_qc_for_sample(
            sample_id=sample_id,
            pycistopic_qc_output_dir=os.path.join(out_dir, "qc"),
            unique_fragments_threshold=None,
            tss_enrichment_threshold=None,
            frip_threshold=0,
            use_automatic_thresholds=True,
        )
        sample_id_to_barcodes_passing_filters[sample_id] = barcodes
        sample_id_to_thresholds[sample_id] = thresholds
    except FileNotFoundError as e:
        print(f"⚠️ File not found for sample {sample_id}: {e}")
    except Exception as e:
        print(f"⚠️ Error getting barcodes for {sample_id}: {e}")

# === Summary: barcodes passing QC ===
for sample_id, barcodes in sample_id_to_barcodes_passing_filters.items():
    print(f"🧬 Sample {sample_id} has {len(barcodes)} barcodes passing QC.")

# === Function to regenerate QC metrics for a sample ===
def regenerate_qc_metrics(sample_id, fragments_file, regions_bed_filename, tss_bed_filename, output_dir):
    command = [
        "/data/sarkern2/conda/envs/scenicplus3/bin/pycistopic", "qc",
        "--fragments", fragments_file,
        "--regions", regions_bed_filename,
        "--tss", tss_bed_filename,
        "--output", output_dir
    ]
    try:
        subprocess.run(command, check=True)
        print(f"🔁 QC metrics regenerated for sample {sample_id}.")
    except subprocess.CalledProcessError as e:
        print(f"❌ Error regenerating QC for sample {sample_id}: {e}")

# === Regenerate QC metrics for missing files ===
for sample_id in filtered_fragments_dict:
    sample_output_dir = os.path.join(out_dir, "qc", sample_id)
    sample_metrics_path = os.path.join(sample_output_dir, f"{sample_id}.fragments_stats_per_cb.parquet")

    if not os.path.exists(sample_metrics_path):
        print(f"🔍 Missing metrics for {sample_id} — regenerating.")
        regenerate_qc_metrics(
            sample_id,
            filtered_fragments_dict[sample_id],
            regions_bed_filename,
            tss_bed_filename,
            sample_output_dir
        )

# === Final check: verify regenerated files ===
for sample_id in filtered_fragments_dict:
    sample_metrics_path = os.path.join(out_dir, "qc", sample_id, f"{sample_id}.fragments_stats_per_cb.parquet")
    if not os.path.exists(sample_metrics_path):
        print(f"❌ File still missing after regeneration: {sample_metrics_path}")
    elif not os.access(sample_metrics_path, os.R_OK):
        print(f"⚠️ File exists but is not readable: {sample_metrics_path}")
    else:
        print(f"✅ Verified metrics file for {sample_id}")




# Define the save_cistopic_object function
def save_cistopic_object(cistopic_obj, sample_id, output_dir):
    """Saves a CistopicObject to a file."""
    output_file = os.path.join(output_dir, f"{sample_id}_cistopic_obj.pkl")
    with open(output_file, "wb") as f:
        pickle.dump(cistopic_obj, f)
    print(f"Saved CistopicObject for {sample_id} to {output_file}")



# Paths and directories
out_dir = "outs_trial"
individual_cistopic_objects_dir = os.path.join(out_dir, "individual_cistopic_objects")
os.makedirs(individual_cistopic_objects_dir, exist_ok=True)
pycistopic_qc_output_dir = "outs_trial/qc"

# Sample IDs (only filtered ones!)
sample_ids = list(filtered_fragments_dict.keys())

# Load barcodes passing QC
sample_id_to_barcodes_passing_filters = {}
sample_id_to_thresholds = {}
for sample_id in filtered_fragments_dict:
    try:
        (
            sample_id_to_barcodes_passing_filters[sample_id],
            sample_id_to_thresholds[sample_id]
        ) = get_barcodes_passing_qc_for_sample(
                sample_id=sample_id,
                pycistopic_qc_output_dir=pycistopic_qc_output_dir,
                unique_fragments_threshold=None,
                tss_enrichment_threshold=None,
                frip_threshold=0,
                use_automatic_thresholds=True,
        )
    except FileNotFoundError as e:
        print(f"File not found: {e}")
    except Exception as e:
        print(f"An error occurred while getting barcodes passing QC for sample {sample_id}: {e}")



# Load or create cistopic objects
cistopic_obj_list = []
for sample_id in sample_ids:
    cistopic_obj_file = os.path.join(individual_cistopic_objects_dir, f"{sample_id}_cistopic_obj.pkl")
    if os.path.exists(cistopic_obj_file):
        with open(cistopic_obj_file, "rb") as f:
            cistopic_obj = pickle.load(f)
            cistopic_obj_list.append(cistopic_obj)
            print(f"Loaded CistopicObject for {sample_id} from {cistopic_obj_file}")
    else:
        sample_metrics_path = os.path.join(pycistopic_qc_output_dir, f"{sample_id}.fragments_stats_per_cb.parquet")

        if not os.path.exists(sample_metrics_path):
            print(f"Metrics file missing for sample {sample_id}: {sample_metrics_path}")
            continue

        sample_metrics = pl.read_parquet(sample_metrics_path).to_pandas()

        if sample_id not in sample_id_to_barcodes_passing_filters:
            print(f"No barcodes passing QC for sample {sample_id}. Skipping...")
            continue

        sample_metrics = sample_metrics.set_index("CB").loc[sample_id_to_barcodes_passing_filters[sample_id]]

        cistopic_obj = create_cistopic_object_from_fragments(
            path_to_fragments=filtered_fragments_dict[sample_id],  # ✅ FILTERED
            path_to_regions=regions_bed_filename,
            path_to_blacklist=path_to_blacklist,
            metrics=sample_metrics,
            valid_bc=sample_id_to_barcodes_passing_filters[sample_id],
            n_cpu=20,
            project=sample_id,
            split_pattern='-'
        )
        cistopic_obj_list.append(cistopic_obj)

        save_cistopic_object(cistopic_obj, sample_id, individual_cistopic_objects_dir)









# Load individual cistopic objects
individual_cistopic_objects_dir = "outs_trial/individual_cistopic_objects"

cistopic_obj_list = []

for filename in os.listdir(individual_cistopic_objects_dir):
    if filename.endswith("_cistopic_obj.pkl"):
        with open(os.path.join(individual_cistopic_objects_dir, filename), "rb") as f:
            cistopic_obj = pickle.load(f)
            cistopic_obj_list.append(cistopic_obj)
            print(f"Loaded CistopicObject from {filename}")
from pycisTopic.cistopic_class import *
cistopic_obj = merge(cistopic_obj_list)

# Save
with open(os.path.join(outDir, 'cisTopicObject_merged.pkl'), 'wb') as f:
  pickle.dump(cistopic_obj, f)



#doublet removal
scrub = scr.Scrublet(cistopic_obj.fragment_matrix.T, expected_doublet_rate=0.1)
doublet_scores, predicted_doublets = scrub.scrub_doublets()
predicted_doublets = scrub.call_doublets(threshold=0.22)  # capture output
scrub.plot_histogram();
scrublet = pd.DataFrame([scrub.doublet_scores_obs_, scrub.predicted_doublets_], columns=cistopic_obj.cell_names, index=['Doublet_scores_fragments', 'Predicted_doublets_fragments']).T

cistopic_obj.add_cell_data(scrublet, split_pattern = '-')
sum(cistopic_obj.cell_data.Predicted_doublets_fragments == True)

# Remove doublets
singlets = cistopic_obj.cell_data[cistopic_obj.cell_data.Predicted_doublets_fragments == False].index.tolist()
# Subset cisTopic object
cistopic_obj = cistopic_obj.subset(singlets, copy=True, split_pattern='-')
print(cistopic_obj)

print(f"🔍 Total cells before doublet filtering: {len(cistopic_obj.cell_names)}")
print(f"🧪 Detected doublets: {sum(scrublet['Predicted_doublets_fragments'])}")
print(f"✅ Remaining singlets: {len(singlets)}")


# Clean the cell_data by removing rows with any NaN values
nan_count = cistopic_obj.cell_data.isna().sum().sum()
print(f"Found {nan_count} total NaN entries in cell metadata.")

if nan_count > 0:
    print("⚠️ Dropping rows with NaNs from cell_data.")
    cistopic_obj.cell_data = cistopic_obj.cell_data.dropna()


# Save final doublet filtered merged cisTopic object
with open(os.path.join(outDir, 'cisTopicObject_merged_dbl_filtered.pkl'), 'wb') as f:
    pickle.dump(cistopic_obj, f)


# === Final Sanity Checks ===

print(f"🧬 Final cisTopic object shape: {cistopic_obj.fragment_matrix.shape}")
print(f"📋 Metadata columns in cell_data:\n{cistopic_obj.cell_data.columns.tolist()}")

# Count cells with complete metadata (no NaNs)
complete_meta = cistopic_obj.cell_data.dropna()
print(f"✅ Final number of singlets with complete metadata: {complete_meta.shape[0]}")

# Optional: summary by sample
summary = complete_meta['sample_id'].value_counts().sort_index()
print("\n🧾 Cell counts per sample (after doublet removal and metadata merge):")
print(summary)

# Optional: save summary to CSV
summary.to_csv(os.path.join(outDir, "final_cell_counts_per_sample.csv"), header=["cell_count"])

