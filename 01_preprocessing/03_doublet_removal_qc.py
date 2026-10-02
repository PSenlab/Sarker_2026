#!/usr/bin/env python3
#===============================================================================
# Step 03: snRNA-seq doublet removal and QC filtering
#===============================================================================
# Description: Loads the SoupX-corrected snRNA-seq counts for every sample,
#              removes doublets with Scrublet (run separately on each sample,
#              since doublets only form within one 10x channel), then computes
#              QC metrics and filters low-quality nuclei.
#
# Steps:       1. Load SoupX-corrected count matrices (Step 02 output)
#              2. Detect doublets with Scrublet, per sample
#              3. Combine samples and remove doublets
#              4. Calculate QC metrics (mitochondrial, ribosomal, hemoglobin)
#              5. Filter genes and nuclei by QC thresholds
#
# Input:       <input_dir>/<sample_id>/   SoupX-corrected 10x v3 matrices
#
# Output:      <output_dir>/filtered_singlets.h5ad       singlets passing QC,
#                                                        raw counts in X
#              <output_dir>/cell_counts_per_sample.csv   nuclei kept at each step
#              <output_dir>/qc/<sample_id>_scrublet.pdf  doublet score histograms
#              <output_dir>/versions.txt                 Python package versions
#
# Usage:       python 03_doublet_removal_qc.py
#
# References:
#   - SoupX: Young MD, Behjati S (2020). GigaScience, 9(12):giaa151
#   - Scrublet: Wolock SL, Lopez R, Klein AM (2019). Cell Systems,
#     8(4):281-291.e9
#===============================================================================

import os
import sys
import logging
from importlib.metadata import version, PackageNotFoundError

import matplotlib
matplotlib.use("Agg")   # plots are written to files; no display needed
import matplotlib.pyplot as plt

import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import scrublet as scr

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)

#-------------------------------------------------------------------------------
# Configuration: update paths for your environment
#-------------------------------------------------------------------------------
INPUT_DIR = "path/to/soupx/corrected"     # Step 02 output folder
OUTPUT_DIR = "path/to/qc_filtered"
OUTPUT_FILE = "filtered_singlets.h5ad"

# Sample sheet: same sample IDs as Steps 01 and 02 (5 age groups x 8 animals).
# These age_group labels are used by the downstream scripts.
AGE_GROUPS = {
    "young":     [f"Y{i}" for i in range(1, 9)],
    "mid_age":   [f"MA{i}" for i in range(1, 9)],
    "old":       [f"O{i}" for i in range(1, 9)],
    "pre_ger":   [f"PG{i}" for i in range(1, 9)],
    "geriatric": [f"G{i}" for i in range(1, 9)],
}

# Scrublet settings (package defaults, written out for reproducibility)
SCRUBLET_PARAMS = {
    "expected_doublet_rate": 0.1,
    "random_state": 0,
}

# Manual doublet score thresholds, only for samples where Scrublet cannot
# set a threshold automatically. Example: {"O3": 0.25}
SCRUBLET_THRESHOLDS = {}

# QC thresholds
QC_PARAMS = {
    "min_genes": 200,
    "max_genes": 10000,
    "min_counts": 500,
    "max_counts": 100000,
    "max_pct_hb": 10,
    "max_pct_ribo": 30,
    "min_cells_per_gene": 5,
}

#-------------------------------------------------------------------------------
# Steps 1 and 2: load each sample and score doublets with Scrublet
#-------------------------------------------------------------------------------
def load_sample(sample_id, age_group, input_dir):
    """Read one SoupX-corrected 10x matrix."""
    adata = sc.read_10x_mtx(os.path.join(input_dir, sample_id), cache=False)
    adata.var_names_make_unique()
    adata.obs["age_group"] = age_group
    return adata


def score_doublets(adata, sample_id, qc_dir):
    """Run Scrublet on one sample. Returns False if no threshold could be set."""
    scrub = scr.Scrublet(adata.X, **SCRUBLET_PARAMS)
    scores, calls = scrub.scrub_doublets(verbose=False)

    if sample_id in SCRUBLET_THRESHOLDS:
        calls = scrub.call_doublets(threshold=SCRUBLET_THRESHOLDS[sample_id],
                                    verbose=False)

    if calls is None:
        logger.error(f"{sample_id}: Scrublet could not set a doublet "
                     f"threshold automatically")
        return False

    adata.obs["doublet_scores"] = scores
    adata.obs["predicted_doublets"] = calls

    fig, _ = scrub.plot_histogram()
    fig.savefig(os.path.join(qc_dir, f"{sample_id}_scrublet.pdf"))
    plt.close(fig)

    logger.info(f"{sample_id}: {int(calls.sum())} doublets of {adata.n_obs} "
                f"nuclei ({100 * calls.mean():.1f}%), "
                f"threshold {scrub.threshold_:.3f}")
    return True


def load_and_score_samples(input_dir, qc_dir):
    """Load all samples, score doublets per sample, and combine them."""
    adata_list, sample_ids = [], []
    missing, no_threshold = [], []

    for age_group, samples in AGE_GROUPS.items():
        for sample_id in samples:
            if not os.path.isdir(os.path.join(input_dir, sample_id)):
                missing.append(sample_id)
                continue

            adata = load_sample(sample_id, age_group, input_dir)
            logger.info(f"Loaded: {sample_id} ({adata.n_obs} nuclei)")

            if not score_doublets(adata, sample_id, qc_dir):
                no_threshold.append(sample_id)
                continue

            adata_list.append(adata)
            sample_ids.append(sample_id)

    # Stop before writing anything if any sample is missing or unscored
    if missing:
        raise FileNotFoundError(
            f"{len(missing)} sample folder(s) not found in {input_dir}: "
            f"{', '.join(missing)}")
    if no_threshold:
        raise RuntimeError(
            f"Scrublet could not set a threshold for: {', '.join(no_threshold)}."
            f" Check their histograms and add manual thresholds to "
            f"SCRUBLET_THRESHOLDS.")

    # Cell names become <barcode>-<sample_id>, e.g. AAACAGCCAAACATAG-1-Y1
    adata = ad.concat(adata_list, label="sample", keys=sample_ids,
                      index_unique="-", merge="same")
    logger.info(f"Total nuclei loaded: {adata.n_obs}")
    return adata

#-------------------------------------------------------------------------------
# Step 3: remove doublets
#-------------------------------------------------------------------------------
def remove_doublets(adata):
    """Keep predicted singlets only."""
    adata = adata[~adata.obs["predicted_doublets"].astype(bool).values].copy()
    logger.info(f"Nuclei after doublet removal: {adata.n_obs}")
    return adata

#-------------------------------------------------------------------------------
# Step 4: QC metrics
#-------------------------------------------------------------------------------
def calculate_qc_metrics(adata):
    """Annotate gene categories (mouse gene symbols) and compute QC metrics."""
    adata.var["mt"] = adata.var_names.str.startswith("mt-")
    adata.var["ribo"] = adata.var_names.str.startswith(("Rps", "Rpl"))
    adata.var["hb"] = adata.var_names.str.contains("^Hb[bag]")

    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt", "ribo", "hb"],
        percent_top=None, log1p=False, inplace=True
    )
    return adata

#-------------------------------------------------------------------------------
# Step 5: QC filtering
#-------------------------------------------------------------------------------
def filter_cells(adata, params):
    """Filter genes, then nuclei, by the thresholds in QC_PARAMS.

    pct_counts_mt is recorded in obs but is not used as a filter.
    """
    sc.pp.filter_genes(adata, min_cells=params["min_cells_per_gene"])

    obs = adata.obs
    keep = (
        (obs["pct_counts_hb"] < params["max_pct_hb"]) &
        (obs["pct_counts_ribo"] < params["max_pct_ribo"]) &
        (obs["n_genes_by_counts"] > params["min_genes"]) &
        (obs["n_genes_by_counts"] < params["max_genes"]) &
        (obs["total_counts"] > params["min_counts"]) &
        (obs["total_counts"] < params["max_counts"])
    )
    adata = adata[keep.values].copy()
    logger.info(f"Nuclei after QC filtering: {adata.n_obs}")
    return adata

#-------------------------------------------------------------------------------
# Records
#-------------------------------------------------------------------------------
def count_per_sample(adata, sample_order):
    return adata.obs["sample"].value_counts().reindex(sample_order, fill_value=0)


def write_versions(path):
    packages = ["scanpy", "anndata", "scrublet", "numpy", "pandas",
                "scipy", "scikit-learn"]
    with open(path, "w") as f:
        f.write(f"python {sys.version.split()[0]}\n")
        for pkg in packages:
            try:
                f.write(f"{pkg} {version(pkg)}\n")
            except PackageNotFoundError:
                f.write(f"{pkg} not installed\n")

#-------------------------------------------------------------------------------
# Main pipeline
#-------------------------------------------------------------------------------
def main():
    qc_dir = os.path.join(OUTPUT_DIR, "qc")
    os.makedirs(qc_dir, exist_ok=True)
    write_versions(os.path.join(OUTPUT_DIR, "versions.txt"))

    sample_order = [s for samples in AGE_GROUPS.values() for s in samples]
    sample_to_group = {s: g for g, samples in AGE_GROUPS.items() for s in samples}

    # Steps 1 and 2: load and score doublets per sample
    adata = load_and_score_samples(INPUT_DIR, qc_dir)
    n_loaded = count_per_sample(adata, sample_order)

    # Step 3: remove doublets
    adata = remove_doublets(adata)
    n_singlets = count_per_sample(adata, sample_order)

    # Steps 4 and 5: QC metrics and filtering
    adata = calculate_qc_metrics(adata)
    adata = filter_cells(adata, QC_PARAMS)
    n_final = count_per_sample(adata, sample_order)

    # Per-sample record of nuclei kept at each step
    summary = pd.DataFrame({
        "sample": sample_order,
        "age_group": [sample_to_group[s] for s in sample_order],
        "loaded": n_loaded.values,
        "doublets_removed": (n_loaded - n_singlets).values,
        "doublet_pct": (100 * (n_loaded - n_singlets) / n_loaded).round(2).values,
        "after_doublet_removal": n_singlets.values,
        "removed_by_qc": (n_singlets - n_final).values,
        "final": n_final.values,
    })
    summary.to_csv(os.path.join(OUTPUT_DIR, "cell_counts_per_sample.csv"),
                   index=False)
    logger.info(f"\n{summary.to_string(index=False)}")

    # Save
    output_path = os.path.join(OUTPUT_DIR, OUTPUT_FILE)
    adata.write_h5ad(output_path)
    logger.info(f"Saved {adata.n_obs} nuclei to {output_path}")


if __name__ == "__main__":
    main()
