#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Lymphoid-T subclustering — liver aging multiome (RNA).

Subsets the T/NK (lymp_T) compartment from the full annotated atlas, re-integrates
it with scVI on subset-specific HVGs, and runs a Leiden resolution sweep with a
matched UMAP on the scVI latent.

Why re-integrate the subset rather than reuse the global embedding:
the subset inherits a neighbors graph and HVG set shaped by the global,
hepatocyte-dominated atlas. We therefore (a) clear the stale graph and
(b) re-derive HVGs + scVI on the T cells alone, so the clustering and the UMAP
reflect T-cell variation rather than global structure. Clustering and UMAP are
both computed on the same scVI latent, so they are matched.

Expected input AnnData:
    .layers['counts']  raw counts
    .X                 log-normalized expression
    .obs[CELLTYPE_KEY] major-type annotation containing CELLTYPE_VAL
    .obs[BATCH_KEY]    batch / sample id for integration

Dependencies: scanpy, scvi-tools, torch, pymde, igraph + leidenalg (Leiden backend).

Usage:
    python tcell_subcluster_scvi.py
"""

import scanpy as sc
import scvi
import torch
import pymde

# --------------------------------------------------------------------------- #
# Config
# --------------------------------------------------------------------------- #
INPUT_H5AD   = "/data/sarkern2/multiome_liver/final_object/final_rna_wnn.h5ad"
OUTPUT_H5AD  = "/data/sarkern2/multiome_liver/final_object/lymp_T_scvi_subclustered.h5ad"

CELLTYPE_KEY = "celltype"      # obs column with the major-type label
CELLTYPE_VAL = "lymp_T"        # compartment to subset (change to re-use for other lineages)
BATCH_KEY    = "sample"        # batch key for scVI / HVG
COUNTS_LAYER = "counts"        # raw counts layer

N_TOP_GENES  = 3000
RESOLUTIONS  = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
SEED         = 0

sc.settings.verbosity = 1
scvi.settings.seed = SEED
torch.set_float32_matmul_precision("high")


def main():
    # --- 1. load + subset to the compartment of interest ---
    adata = sc.read_h5ad(INPUT_H5AD)
    adata = adata[adata.obs[CELLTYPE_KEY] == CELLTYPE_VAL].copy()
    print(f"{adata.n_obs} {CELLTYPE_VAL} cells")

    # --- 2. clear the stale (global-atlas) graph so it isn't reused downstream ---
    for k in ("connectivities", "distances"):
        adata.obsp.pop(k, None)
    adata.uns.pop("neighbors", None)

    # --- 3. drop genes detected in no cell of the subset; keep log-norm for DE ---
    sc.pp.filter_genes(adata, min_cells=1)
    adata.layers["lognorm"] = adata.X.copy()
    print(f"{adata.n_vars} genes kept")

    # --- 4. subset-specific HVGs from raw counts (batch-aware; keep all genes) ---
    sc.pp.highly_variable_genes(
        adata, flavor="seurat_v3", n_top_genes=N_TOP_GENES,
        layer=COUNTS_LAYER, batch_key=BATCH_KEY, subset=False,
    )

    # --- 5. train scVI on the HVG subset (raw counts, batch = sample) ---
    adata_hvg = adata[:, adata.var.highly_variable].copy()
    scvi.model.SCVI.setup_anndata(adata_hvg, layer=COUNTS_LAYER, batch_key=BATCH_KEY)
    vae = scvi.model.SCVI(adata_hvg)
    vae.train()
    adata.obsm["X_scVI"] = vae.get_latent_representation()   # cells align with full object

    # --- 5b. MDE embedding of the scVI latent (alternative 2D view to UMAP) ---
    adata.obsm["X_scVI_MDE"] = (
        pymde.preserve_neighbors(adata.obsm["X_scVI"], embedding_dim=2)
        .embed().cpu().numpy()
    )

    # --- 6. neighbors -> Leiden sweep -> UMAP, all on the scVI latent ---
    sc.pp.neighbors(adata, use_rep="X_scVI", random_state=SEED)
    for res in RESOLUTIONS:
        key = f"leiden_scvi_{int(round(res * 10))}"
        sc.tl.leiden(
            adata, resolution=res, key_added=key,
            flavor="igraph", n_iterations=2, directed=False, random_state=SEED,
        )
        print(f"{key}: {adata.obs[key].nunique()} clusters")
    sc.tl.umap(adata, random_state=SEED)   # UMAP on the same scVI neighbors graph

    # --- 7. save ---
    adata.write_h5ad(OUTPUT_H5AD)
    print(f"saved -> {OUTPUT_H5AD}")


if __name__ == "__main__":
    main()
