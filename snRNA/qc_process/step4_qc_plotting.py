#!/usr/bin/env python3
#===============================================================================
# QC visualization for single-nucleus RNA-seq data
#===============================================================================
# Description: Generates QC violin plots (RNA metrics by age group) and
#              UMAP overlays (mito%, ribo%, Hb%) for quality assessment
#
# Input:       Integrated AnnData with scVI embedding (from scvi_integration.py)
# Output:      RNA QC violin plots (n_counts, n_genes)
#              UMAP grid (age x QC metric: mito, ribo, Hb)
# #===============================================================================

import logging
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.gridspec import GridSpec

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)

# Global plot settings
mpl.rcParams["font.family"] = "Arial"
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42
mpl.rcParams["axes.linewidth"] = 0.8
mpl.rcParams["xtick.major.width"] = 0.8
mpl.rcParams["ytick.major.width"] = 0.8

#-------------------------------------------------------------------------------
# Configuration
#-------------------------------------------------------------------------------
INPUT_FILE = "integrated_scvi.h5ad"
OUTPUT_DIR = "figures"

# Age group order and colors
AGE_ORDER = ["young", "mid-age", "old", "pre-geriatric", "geriatric"]
AGE_COLORS = {
    "young": "#a6cee3",
    "mid-age": "#fdbf6f",
    "old": "#b2df8a",
    "pre-geriatric": "#cab2d6",
    "geriatric": "#fb9a99",
}

# QC metrics to plot
RNA_METRICS = {
    "n_counts_RNA": "n counts RNA",
    "n_genes_by_counts": "n genes by counts",
}

UMAP_METRICS = {
    "pct_counts_mito": "% counts mito",
    "pct_counts_ribo": "% counts ribo",
    "pct_counts_Hb": "% counts Hb",
}

# UMAP embedding key
UMAP_BASIS = "X_scVI_MDE"

# Age group column in adata.obs
AGE_COL = "age"


#-------------------------------------------------------------------------------
# Plotting Functions
#-------------------------------------------------------------------------------
def plot_qc_violins(adata, metrics, panel_label, filename):
    """
    Generate paired violin plots for QC metrics split by age group.
    Matches panels c/d in the figure.
    """
    fig, axes = plt.subplots(1, len(metrics), figsize=(5 * len(metrics), 4))
    if len(metrics) == 1:
        axes = [axes]

    for ax, (metric, label) in zip(axes, metrics.items()):
        parts = ax.violinplot(
            [adata.obs.loc[adata.obs[AGE_COL] == age, metric].dropna().values
             for age in AGE_ORDER],
            positions=range(len(AGE_ORDER)),
            showmeans=False,
            showmedians=True,
            showextrema=False,
        )

        # Style violin bodies
        for i, body in enumerate(parts["bodies"]):
            body.set_facecolor(AGE_COLORS[AGE_ORDER[i]])
            body.set_edgecolor("black")
            body.set_linewidth(0.5)
            body.set_alpha(0.85)

        # Style median lines
        parts["cmedians"].set_edgecolor("black")
        parts["cmedians"].set_linewidth(0.8)

        ax.set_xticks(range(len(AGE_ORDER)))
        ax.set_xticklabels(AGE_ORDER, rotation=45, ha="right", fontsize=9)
        ax.set_ylabel(label, fontsize=10)
        ax.tick_params(axis="y", labelsize=9)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    fig.text(0.01, 0.98, panel_label, fontsize=14, fontweight="bold", va="top")
    plt.tight_layout()
    plt.savefig(f"{OUTPUT_DIR}/{filename}", bbox_inches="tight", dpi=300)
    plt.close()
    logger.info(f"Saved {filename}")


def plot_umap_qc_grid(adata, metrics, panel_label, filename):
    """
    Generate age x metric UMAP grid overlay.
    Matches panel e in the figure: columns = age groups, rows = QC metrics.
    """
    n_rows = len(metrics)
    n_cols = len(AGE_ORDER)

    fig, axes = plt.subplots(
        n_rows, n_cols,
        figsize=(3 * n_cols, 3 * n_rows),
        gridspec_kw={"wspace": 0.05, "hspace": 0.15},
    )

    # Get UMAP coordinates
    umap_key = UMAP_BASIS
    coords = adata.obsm[umap_key]

    for row_idx, (metric, label) in enumerate(metrics.items()):
        # Compute consistent color range across all age groups
        vals = adata.obs[metric].values
        vmin = np.nanpercentile(vals, 1)
        vmax = np.nanpercentile(vals, 99)

        for col_idx, age in enumerate(AGE_ORDER):
            ax = axes[row_idx, col_idx] if n_rows > 1 else axes[col_idx]
            mask = adata.obs[AGE_COL] == age

            # Background: all cells in light gray
            ax.scatter(
                coords[:, 0], coords[:, 1],
                s=0.5, c="lightgray", alpha=0.3, rasterized=True,
            )

            # Overlay: age-specific cells colored by metric
            sc_plot = ax.scatter(
                coords[mask, 0], coords[mask, 1],
                s=0.5,
                c=adata.obs.loc[mask, metric].values,
                cmap="Reds",
                vmin=vmin,
                vmax=vmax,
                alpha=0.8,
                rasterized=True,
            )

            ax.set_xticks([])
            ax.set_yticks([])
            ax.axis("off")

            # Column titles (top row only)
            if row_idx == 0:
                ax.set_title(age, fontsize=10, fontweight="bold")

        # Colorbar for each row
        cbar = fig.colorbar(
            sc_plot, ax=axes[row_idx, :] if n_rows > 1 else axes,
            shrink=0.6, aspect=15, pad=0.02,
        )
        cbar.set_label(label, fontsize=9)
        cbar.ax.tick_params(labelsize=8)

    fig.text(0.01, 0.98, panel_label, fontsize=14, fontweight="bold", va="top")
    plt.savefig(f"{OUTPUT_DIR}/{filename}", bbox_inches="tight", dpi=300)
    plt.close()
    logger.info(f"Saved {filename}")


#-------------------------------------------------------------------------------
# Main Pipeline
#-------------------------------------------------------------------------------
def main():
    import os
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # Load data
    logger.info(f"Loading {INPUT_FILE}...")
    adata = sc.read_h5ad(INPUT_FILE)
    logger.info(f"Data shape: {adata.shape}")

    # Ensure age group is categorical with correct order
    adata.obs[AGE_COL] = pd.Categorical(
        adata.obs[AGE_COL], categories=AGE_ORDER, ordered=True
    )

    # -------------------------------------------------------------------------
    # Panel c: RNA QC violins
    # -------------------------------------------------------------------------
    logger.info("Generating Panel c: RNA QC violin plots...")
    plot_qc_violins(adata, RNA_METRICS, panel_label="c", filename="panel_c_rna_qc_violins.pdf")

    # -------------------------------------------------------------------------
    # Panel d: UMAP QC grid (age x metric)
    # -------------------------------------------------------------------------
    logger.info("Generating Panel d: UMAP QC metric grid...")
    plot_umap_qc_grid(adata, UMAP_METRICS, panel_label="d", filename="panel_d_umap_qc_grid.pdf")

    logger.info("QC visualization complete.")


if __name__ == "__main__":
    main()
