#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Compartment subclustering — liver aging multiome (RNA modality).

Generic scVI subclustering for any major compartment of the annotated atlas.
Subsets a compartment (one or more coarse labels), re-integrates it with scVI on
subset-specific highly variable genes (HVGs), and runs a Leiden resolution sweep
with a matched UMAP computed on the scVI latent.

The pipeline is compartment-agnostic. Everything compartment-specific — the
coarse label(s) to subset, the marker panel, the cluster->celltype map, the
group order, and the output obs column — lives in the PANELS registry below and
is selected with ``--compartment``. Adding a new compartment = adding one entry
to PANELS; the pipeline body never changes.

Why re-integrate the subset rather than reuse the global embedding
------------------------------------------------------------------
The subset inherits a neighbors graph and HVG set shaped by the global,
hepatocyte-dominated atlas. We therefore (a) clear the stale graph and
(b) re-derive HVGs + scVI on the subset alone, so that clustering and the UMAP
reflect within-compartment variation rather than global structure. Clustering
and UMAP are both computed on the same scVI latent and are therefore matched.

Expected input AnnData
----------------------
    .layers['counts']   raw counts
    .X                  log-normalized expression
    .obs[celltype_key]  major-type annotation containing the compartment label(s)
    .obs[batch_key]     batch / sample id for integration

Reproducibility
---------------
Seeds are fixed for scVI, PyTorch, and every downstream Scanpy step. The exact
package versions used are logged at run time (``--print-versions``); the versions
this analysis was validated against are pinned in ``requirements.txt``.

Usage
-----
    # T / NK compartment
    python subcluster_scvi.py --compartment lymph_T \
        --input  data/final_rna_wnn.h5ad \
        --output results/lymp_T_scvi_subclustered.h5ad

    # Myeloid compartment (Kupffer + MoMFs)
    python subcluster_scvi.py --compartment myeloid \
        --input  data/final_rna_wnn.h5ad \
        --output results/myeloid_scvi_subclustered.h5ad

Run ``python subcluster_scvi.py --help`` for all options.

Citation
--------
If you use this code, please cite:
    Sarker N. et al. "<paper title>." <journal>, <year>. DOI: <doi>

License: MIT (see LICENSE).
"""

import argparse
import sys

import scanpy as sc              # noqa: E402
import scvi                      # noqa: E402
import torch                     # noqa: E402
import pymde                     # noqa: E402


# --------------------------------------------------------------------------- #
# Compartment registry (dataset-specific)
#
# One entry per compartment. Each holds everything the generic pipeline needs to
# subset, annotate, and validate that lineage:
#
#   subset_labels : coarse label(s) in obs[celltype_key] defining the compartment
#   annotation    : cluster (str) -> fine celltype, for the resolution in
#                   annotate_key. NOTE: re-integrating a subset yields NEW cluster
#                   numbers each time — run once, inspect the clusters at
#                   annotate_key, then edit these. Unmapped clusters fail loudly.
#   annotate_key  : leiden key whose clusters `annotation` maps from
#   group_order   : left-to-right final category order (== legend order)
#   drop_labels   : fine labels removed after annotation (junk / low-quality)
#   celltype_out  : obs column the fine annotation is written to
# --------------------------------------------------------------------------- #
PANELS = {
    "lymph_T": {
        "subset_labels": ["T/ILC cells"],
        "annotate_key": "leiden_scvi_8",
        "annotation": {
            "0":  "CD4T",
            "1":  "CD8T",
            "2":  "CD8T",
            "3":  "iNKT",
            "4":  "Treg",
            "5":  "neutrophil",
            "6":  "NK",
            "7":  "ILC1",
            "8":  "CD8T",
            "9":  "CD8T",
            "10": "gdT",
            "11": "Low-quality",
        },
        "group_order": ["CD4T", "CD8T", "Treg", "iNKT", "gdT", "NK", "ILC1", "neutrophil"],
        "drop_labels": ["Low-quality"],
        "celltype_out": "T/ILC cells",
    },

    "myeloid": {
        "subset_labels": ["Kupffer 01", "non-resident myeloid"],
        "annotate_key": "leiden_scvi_8",
        "annotation": {
            # PLACEHOLDER seeded from the earlier leiden_scvi_3 calls — will NOT
            # match this fresh re-integration. Inspect the clusters, then edit.
            "0":  "Kupffer",
            "1":  "Kupffer",
            "2":  "LAM",
            "3":  "Kupffer_cycling",
            "4":  "neutrophil",
            "5":  "cDC1",
            "6":  "pDC",
            "7":  "MoMF",
            "8":  "MoMF",
        },
        "group_order": ["Kupffer", "Kupffer_cycling", "LAM", "MoMF", "cDC1", "pDC", "neutrophil"],
        "drop_labels": [],
        "celltype_out": "myeloid",
    },

    "endo": {
        "subset_labels": ["endothelial", "Kupffer 02"],
        "annotate_key": "leiden_scvi_8",
        "annotation": {
            # PLACEHOLDER — inspect the clusters at annotate_key, then edit.
            "0":  "LSEC",
            "1":  "LSEC",
            "2":  "MV portal",
            "3":  "MV central",
            "4":  "LSEC cycling",
            "5":  "Kupffer-like",
        },
        "group_order": ["LSEC", "LSEC cycling", "MV portal", "MV central", "Kupffer-like"],
        "drop_labels": [],
        "celltype_out": "endothelial_Kupffer02",
    },

    # Kupffer 01 + Kupffer 02 combined (resident Kupffer + the reassigned
    # Endothelial-02 / Kupffer-like population). Re-integrating them together with
    # scVI checks whether the two co-cluster; if they do, the combined identity is
    # simply "Kupffer". Run WITHOUT --annotate (annotation left empty on purpose);
    # inspect the Leiden sweep, then fill in later.
    "kupffer": {
        "subset_labels": ["Kupffer 01", "Kupffer 02"],
        "annotate_key": "leiden_scvi_8",
        "annotation": {},        # fill after inspecting the sweep, if annotating
        "group_order": [],
        "drop_labels": [],
        "celltype_out": "Kupffer",
    },
}


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #
def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="Generic scVI subclustering of a chosen compartment.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    # compartment selection
    p.add_argument("--compartment", required=True, choices=sorted(PANELS),
                   help="Which compartment panel to use (see PANELS registry).")

    # I/O
    p.add_argument("--input", required=True,
                   help="Path to the annotated atlas .h5ad (RNA).")
    p.add_argument("--output", required=True,
                   help="Path to write the subclustered .h5ad.")

    # obs / layer keys
    p.add_argument("--celltype-key", default="celltype",
                   help="obs column holding the major-type annotation.")
    p.add_argument("--celltype-val", nargs="+", default=None,
                   help="Override the compartment's coarse label(s) to subset. "
                        "Defaults to PANELS[compartment]['subset_labels'].")
    p.add_argument("--batch-key", default="sample",
                   help="obs column used as batch for scVI / HVG.")
    p.add_argument("--counts-layer", default="counts",
                   help="Layer holding raw counts.")

    # analysis params
    p.add_argument("--n-top-genes", type=int, default=3000,
                   help="Number of subset-specific HVGs.")
    p.add_argument("--resolutions", type=float, nargs="+",
                   default=[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0],
                   help="Leiden resolutions to sweep.")
    p.add_argument("--max-epochs", type=int, default=None,
                   help="scVI training epochs (None = scvi-tools heuristic).")
    p.add_argument("--seed", type=int, default=0,
                   help="Global random seed.")

    # annotation
    p.add_argument("--annotate", action="store_true",
                   help="Apply the manual annotation map after the sweep.")
    p.add_argument("--annotate-key", default=None,
                   help="Override the compartment's annotate_key.")

    p.add_argument("--print-versions", action="store_true",
                   help="Log key package versions and exit environment info.")
    return p.parse_args(argv)


# --------------------------------------------------------------------------- #
# Helpers
# --------------------------------------------------------------------------- #
def log_versions():
    """Record versions of the packages that affect numerical results."""
    import anndata
    import numpy
    import leidenalg
    print("=== environment ===")
    print(f"python      {sys.version.split()[0]}")
    print(f"scanpy      {sc.__version__}")
    print(f"anndata     {anndata.__version__}")
    print(f"scvi-tools  {scvi.__version__}")
    print(f"torch       {torch.__version__}")
    print(f"numpy       {numpy.__version__}")
    print(f"leidenalg   {leidenalg.version}")
    print(f"pymde       {pymde.__version__}")
    print(f"cuda        {torch.cuda.is_available()}")
    print("===================")


def leiden_key(res):
    """Key name for a Leiden resolution, e.g. 0.8 -> 'leiden_scvi_8'.

    Matches the naming used throughout the analysis and the annotation maps.
    Unique for the 0.1-step sweep used here (resolutions 0.1-1.0 map to keys
    1-10); if you sweep finer than 0.1 steps, switch to f-string
    'leiden_scvi_{res:g}' and update the panel annotate_key accordingly.
    """
    return f"leiden_scvi_{int(round(res * 10))}"


def annotate_and_validate(adata, panel, annotate_key):
    """Map clusters -> cell types and drop junk labels.

    All compartment-specific info comes from `panel` (an entry of PANELS), so the
    function itself is generic.

    Returns a (possibly smaller) AnnData with:
        .obs[panel['celltype_out']]  categorical annotation, drop_labels removed
    """
    annotation   = panel["annotation"]
    group_order  = panel["group_order"]
    drop_labels  = panel["drop_labels"]
    celltype_out = panel["celltype_out"]

    if annotate_key not in adata.obs:
        raise KeyError(
            f"{annotate_key!r} not in obs. Available leiden keys: "
            f"{[c for c in adata.obs if c.startswith('leiden')]}"
        )

    # unmapped clusters would become NaN and silently vanish -> fail loudly
    clusters = set(adata.obs[annotate_key].astype(str).unique())
    unmapped = clusters - set(annotation)
    if unmapped:
        raise ValueError(
            f"Clusters in {annotate_key} with no annotation entry: "
            f"{sorted(unmapped)}"
        )

    adata.obs[celltype_out] = (
        adata.obs[annotate_key].astype(str).map(annotation).astype("category")
    )
    print(f"\n{celltype_out} (from {annotate_key}):")
    print(adata.obs[celltype_out].value_counts())

    # --- drop junk labels, tidy categories ---
    n_before = adata.n_obs
    adata = adata[~adata.obs[celltype_out].isin(drop_labels)].copy()
    adata.obs[celltype_out] = adata.obs[celltype_out].cat.remove_unused_categories()
    order = [g for g in group_order if g in adata.obs[celltype_out].cat.categories]
    adata.obs[celltype_out] = adata.obs[celltype_out].cat.reorder_categories(order)
    print(f"dropped {drop_labels}: {n_before} -> {adata.n_obs} cells")
    return adata


# --------------------------------------------------------------------------- #
# Pipeline
# --------------------------------------------------------------------------- #
def main(argv=None):
    args = parse_args(argv)

    panel = PANELS[args.compartment]
    subset_labels = args.celltype_val or panel["subset_labels"]
    annotate_key  = args.annotate_key or panel["annotate_key"]

    sc.settings.verbosity = 1
    scvi.settings.seed = args.seed
    torch.set_float32_matmul_precision("high")

    log_versions()
    if args.print_versions:
        return

    # --- 1. load + subset to the compartment of interest (one or more labels) ---
    adata = sc.read_h5ad(args.input)
    n_before = adata.n_obs
    adata = adata[adata.obs[args.celltype_key].isin(subset_labels)].copy()
    print(f"subset {args.compartment} {subset_labels}: {adata.n_obs} / {n_before} cells")
    if adata.n_obs == 0:
        raise ValueError(
            f"No cells with {args.celltype_key} in {subset_labels}. "
            f"Available: {sorted(adata.obs[args.celltype_key].unique())}"
        )
    # drop coarse-label categories no longer present so scVI batch/HVG stay clean
    if hasattr(adata.obs[args.batch_key], "cat"):
        adata.obs[args.batch_key] = adata.obs[args.batch_key].cat.remove_unused_categories()

    # --- 2. clear the stale (global-atlas) graph so it isn't reused downstream ---
    for k in ("connectivities", "distances"):
        adata.obsp.pop(k, None)
    adata.uns.pop("neighbors", None)

    # --- 3. drop genes detected in no cell of the subset; keep log-norm for DE ---
    sc.pp.filter_genes(adata, min_cells=1)
    adata.layers["lognorm"] = adata.X.copy()
    print(f"{adata.n_vars} genes kept after subset filtering")

    # --- 4. subset-specific HVGs from raw counts (batch-aware; keep all genes) ---
    sc.pp.highly_variable_genes(
        adata, flavor="seurat_v3", n_top_genes=args.n_top_genes,
        layer=args.counts_layer, batch_key=args.batch_key, subset=False,
    )

    # --- 5. train scVI on the HVG subset (raw counts, batch = sample) ---
    adata_hvg = adata[:, adata.var.highly_variable].copy()
    scvi.model.SCVI.setup_anndata(
        adata_hvg, layer=args.counts_layer, batch_key=args.batch_key
    )
    vae = scvi.model.SCVI(adata_hvg)
    vae.train(max_epochs=args.max_epochs)
    # cells in adata_hvg align with adata, so the latent maps back directly
    adata.obsm["X_scVI"] = vae.get_latent_representation()

    # --- 5b. MDE embedding of the scVI latent (alternative 2D view to UMAP) ---
    adata.obsm["X_scVI_MDE"] = (
        pymde.preserve_neighbors(adata.obsm["X_scVI"], embedding_dim=2)
        .embed().cpu().numpy()
    )

    # --- 6. neighbors -> Leiden sweep -> UMAP, all on the scVI latent ---
    sc.pp.neighbors(adata, use_rep="X_scVI", random_state=args.seed)
    for res in args.resolutions:
        key = leiden_key(res)
        sc.tl.leiden(
            adata, resolution=res, key_added=key,
            flavor="igraph", n_iterations=2, directed=False,
            random_state=args.seed,
        )
        print(f"{key}: {adata.obs[key].nunique()} clusters")
    sc.tl.umap(adata, random_state=args.seed)   # UMAP on the same scVI neighbors graph

    # --- 7. manual annotation (optional) ---
    if args.annotate:
        adata = annotate_and_validate(
            adata, panel=panel, annotate_key=annotate_key,
        )

    # --- 8. save ---
    adata.write_h5ad(args.output)
    print(f"saved -> {args.output}")


if __name__ == "__main__":
    main()
