#!/usr/bin/env python3
#===============================================================================
# Cell Type Annotation for Single-Nucleus RNA-seq Data
#===============================================================================
# Description: Annotates cell types based on Leiden clustering and validates
#              with canonical marker gene expression
#
# Input:       Integrated AnnData with clustering (from scvi_integration.py)
# Output:      Annotated AnnData with cell type labels
#
# Cell Types:  Hepatocytes, Endothelial, Kupffer, Stellate, MoMFs,
#              Cholangiocytes, Lymphoid (T and B cells)
#===============================================================================

import logging
import scanpy as sc

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)

#-------------------------------------------------------------------------------
# Configuration
#-------------------------------------------------------------------------------
INPUT_FILE = "integrated_scvi.h5ad"
OUTPUT_FILE = "annotated.h5ad"

# Marker genes for cell type validation
MARKER_GENES = {
    "Cholangiocyte": ["Kcnma1", "Ctnnd2", "Krt19", "Epcam"],
    "Endothelial": ["Flt1", "Pecam1"],
    "Hepatocyte": ["Cps1", "Cyp2e1", "Cyp7b1"],
    "Kupffer": ["Clec4f", "Adgre1"],
    "MoMFs": ["Cx3cr1", "Ccr2"],
    "Stellate": ["Reln", "Lrat"],
    "Lymphoid_B": ["Ms4a1", "Pax5"],
    "Lymphoid_T": ["Itk", "Camk4"],
}

# Cluster to cell type mapping (based on leiden_5 resolution)
CLUSTER_ANNOTATION = {
    "0": "Hep-01/Hepatocyte",
    "1": "Hep-02/Hepatocyte",
    "2": "Endothelial-01/Endothelial",
    "3": "Hep-03/Hepatocyte",
    "4": "Kupffer",
    "5": "Stellate",
    "6": "Hep-04/Hepatocyte",
    "7": "Hep-05/Hepatocyte",
    "8": "Lymp_T",
    "9": "Hep-06/Hepatocyte",
    "10": "Hep-07/Hepatocyte",
    "11": "MoMFs",
    "12": "Endothelial-02/Endothelial",
    "13": "Lymp_B",
    "14": "Cholangiocyte-01/Cholangiocyte",
    "15": "Cholangiocyte-02/Cholangiocyte",
    "17": "Unassigned",
    "18": "Unassigned",
    "19": "Unassigned",
}

# Cell type color mapping
CELLTYPE_COLORS = {
    "Hep-01/Hepatocyte": "#17becf",
    "Hep-02/Hepatocyte": "#17becf",
    "Hep-03/Hepatocyte": "#17becf",
    "Hep-04/Hepatocyte": "#17becf",
    "Hep-05/Hepatocyte": "#17becf",
    "Hep-06/Hepatocyte": "#17becf",
    "Hep-07/Hepatocyte": "#17becf",
    "Endothelial-01/Endothelial": "#a6cee3",
    "Endothelial-02/Endothelial": "#1f78b4",
    "Stellate": "#fb9a99",
    "Kupffer": "#b2df8a",
    "Lymp_T": "#fdbf6f",
    "MoMFs": "#33a02c",
    "Cholangiocyte-01/Cholangiocyte": "#e31a1c",
    "Cholangiocyte-02/Cholangiocyte": "#cab2d6",
    "Lymp_B": "#e377c2",
    "Unassigned": "#d3d3d3",
}

#-------------------------------------------------------------------------------
# Main Pipeline
#-------------------------------------------------------------------------------
def main():
    # Load data
    logger.info(f"Loading {INPUT_FILE}...")
    adata = sc.read_h5ad(INPUT_FILE)
    logger.info(f"Data shape: {adata.shape}")
    
    # Annotate cell types
    logger.info("Annotating cell types based on leiden_5 clustering...")
    adata.obs["celltype"] = adata.obs["leiden_5"].map(CLUSTER_ANNOTATION).astype("category")
    
    # Assign colors
    adata.uns["celltype_colors"] = [
        CELLTYPE_COLORS[ct] for ct in adata.obs["celltype"].cat.categories
    ]
    
    # Summary
    logger.info(f"Cell type distribution:\n{adata.obs['celltype'].value_counts()}")
    
    # Generate marker gene dotplot
    logger.info("Generating marker gene dotplot...")
    sc.pl.dotplot(adata, MARKER_GENES, groupby="celltype", save="_marker_genes")
    
    # Save
    adata.write_h5ad(OUTPUT_FILE)
    logger.info(f"Saved to {OUTPUT_FILE}")

if __name__ == "__main__":
    main()
