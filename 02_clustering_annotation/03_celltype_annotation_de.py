#!/usr/bin/env python3
#===============================================================================
# Cell Type Annotation & Differential Expression for Single-Nucleus Multi-Ome
#===============================================================================
# Description: Annotates cell types based on Leiden clustering, validates
#              with canonical marker gene expression, performs Wilcoxon DE,
#              and generates dotplots + UMAP marker overlays per cell type
#
# Input:       Integrated AnnData with clustering (from 01_scvi_integration_clustering.py)
# Output:      Annotated AnnData with cell type labels + DE results
#
# Cell Types:  hepatocyte, cholangiocyte 01, cholangiocyte 02, endothelial,
#              Kupffer 01, Kupffer 02, non-resident myeloid, mesenchymal,
#              T/ILC cells, B cells
#===============================================================================

import logging
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)

# Global plot settings
mpl.rcParams["font.family"] = "Arial"
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42

#-------------------------------------------------------------------------------
# Configuration
#-------------------------------------------------------------------------------
INPUT_FILE = "integrated_scvi.h5ad"
OUTPUT_FILE = "annotated.h5ad"

# Marker genes for cell type validation
MARKER_GENES = {
    "cholangiocyte": ["Kcnma1", "Ctnnd2", "Krt19", "Epcam"],
    "endothelial": ["Flt1", "Pecam1"],
    "hepatocyte": ["Cps1", "Cyp2e1", "Cyp7b1"],
    "Kupffer": ["Clec4f", "Adgre1"],
    "non-resident myeloid": ["Cx3cr1", "Ccr2"],
    "mesenchymal": ["Reln", "Lrat"],
    "B cells": ["Ms4a1", "Pax5"],
    "T/ILC cells": ["Itk", "Camk4"],
}

# Cluster to cell type mapping (based on leiden_5 resolution)
CLUSTER_ANNOTATION = {
    "0": "hepatocyte",
    "1": "hepatocyte",
    "2": "endothelial",
    "3": "hepatocyte",
    "4": "Kupffer 01",
    "5": "mesenchymal",
    "6": "hepatocyte",
    "7": "hepatocyte",
    "8": "T/ILC cells",
    "9": "hepatocyte",
    "10": "hepatocyte",
    "11": "non-resident myeloid",
    "12": "Kupffer 02",
    "13": "B cells",
    "14": "cholangiocyte 01",
    "15": "cholangiocyte 02",
    "17": "Unassigned",
    "18": "Unassigned",
    "19": "Unassigned",
}

# Cell type color mapping
CELLTYPE_COLORS = {
    "hepatocyte": "#17becf",
    "endothelial": "#a6cee3",
    "Kupffer 02": "#1f78b4",
    "mesenchymal": "#fb9a99",
    "Kupffer 01": "#b2df8a",
    "T/ILC cells": "#fdbf6f",
    "non-resident myeloid": "#33a02c",
    "cholangiocyte 01": "#e31a1c",
    "cholangiocyte 02": "#cab2d6",
    "B cells": "#e377c2",
    "Unassigned": "#d3d3d3",
}

#-------------------------------------------------------------------------------
# Main Pipeline
#-------------------------------------------------------------------------------
def main():
    # -------------------------------------------------------------------------
    # Step 1: Load data
    # -------------------------------------------------------------------------
    logger.info(f"Loading {INPUT_FILE}...")
    adata = sc.read_h5ad(INPUT_FILE)
    logger.info(f"Data shape: {adata.shape}")

    # -------------------------------------------------------------------------
    # Step 2: Annotate cell types based on leiden_5 clustering
    # -------------------------------------------------------------------------
    logger.info("Annotating cell types...")
    adata.obs["celltype"] = adata.obs["leiden_5"].map(CLUSTER_ANNOTATION).astype("category")

    # Assign colors
    adata.uns["celltype_colors"] = [
        CELLTYPE_COLORS[ct] for ct in adata.obs["celltype"].cat.categories
    ]

    logger.info(f"Cell type distribution:\n{adata.obs['celltype'].value_counts()}")

    # -------------------------------------------------------------------------
    # Step 3: Generate canonical marker gene dotplot
    # -------------------------------------------------------------------------
    logger.info("Generating canonical marker gene dotplot...")
    sc.pl.dotplot(adata, MARKER_GENES, groupby="celltype", save="_canonical_markers.pdf")

    # -------------------------------------------------------------------------
    # Step 4: Perform differential expression analysis (Wilcoxon) per cell type
    # -------------------------------------------------------------------------
    logger.info("Running Wilcoxon rank-sum DE across cell types...")
    sc.tl.rank_genes_groups(adata, groupby="celltype", method="wilcoxon")

    # Extract DE results per cell type
    result = adata.uns["rank_genes_groups"]
    groups = result["names"].dtype.names
    de_results = {}

    for group in groups:
        df = pd.DataFrame({
            "Gene": result["names"][group],
            "Score": result["scores"][group],
            "LogFC": result["logfoldchanges"][group],
            "pvals": result["pvals"][group],
            "pvals_adj": result["pvals_adj"][group],
        })
        de_results[group] = df
        logger.info(f"  {group}: {(df['pvals_adj'] < 0.05).sum()} significant genes (FDR < 0.05)")

    # Save DE results to Excel (one sheet per cell type)
    with pd.ExcelWriter("de_results_wilcoxon.xlsx", engine="openpyxl") as writer:
        for group, df in de_results.items():
            sheet_name = group[:31]  # Excel sheet name limit
            df.to_excel(writer, sheet_name=sheet_name, index=False)
    logger.info("Saved DE results to de_results_wilcoxon.xlsx")

    # -------------------------------------------------------------------------
    # Step 5: Ranked genes dotplot (all cell types)
    # -------------------------------------------------------------------------
    logger.info("Generating ranked genes dotplot (all cell types)...")
    sc.pl.rank_genes_groups_dotplot(
        adata,
        n_genes=5,
        values_to_plot="logfoldchanges",
        min_logfoldchange=6,
        vmax=8,
        vmin=-8,
        cmap="bwr",
        show=False,
    )
    plt.savefig("rank_genes_groups_dotplot_all.pdf", bbox_inches="tight", dpi=300)
    plt.close()
    logger.info("Saved rank_genes_groups_dotplot_all.pdf")

    # -------------------------------------------------------------------------
    # Step 6: UMAP overlay of top marker genes per cell type
    # -------------------------------------------------------------------------
    logger.info("Generating UMAP marker overlays per cell type...")

    # Exclude Unassigned from UMAP overlays
    celltypes_to_plot = [ct for ct in adata.obs["celltype"].cat.categories if ct != "Unassigned"]

    for celltype in celltypes_to_plot:
        top_genes = sc.get.rank_genes_groups_df(adata, group=celltype).head(9)["names"]

        sc.pl.embedding(
            adata,
            basis="X_wnn",
            color=[*top_genes, "celltype"],
            legend_loc="on data",
            frameon=False,
            ncols=3,
            show=False,
        )

        safe_name = celltype.replace("/", "_").replace(" ", "_")
        filename = f"wnn_umap_{safe_name}.pdf"
        plt.savefig(filename, bbox_inches="tight", dpi=300)
        plt.close()
        logger.info(f"  Saved {filename}")

    # -------------------------------------------------------------------------
    # Step 7: Save annotated AnnData
    # -------------------------------------------------------------------------
    adata.write_h5ad(OUTPUT_FILE)
    logger.info(f"Saved annotated AnnData to {OUTPUT_FILE}")


if __name__ == "__main__":
    main()