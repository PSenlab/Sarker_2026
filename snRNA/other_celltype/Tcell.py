#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Lymphoid-T subclustering — liver aging multiome (RNA modality).

Subsets the T/NK (``lymp_T``) compartment from the full annotated atlas,
re-integrates it with scVI on subset-specific highly variable genes (HVGs), and
runs a Leiden resolution sweep with a matched UMAP computed on the scVI latent.

Why re-integrate the subset rather than reuse the global embedding
------------------------------------------------------------------
The subset inherits a neighbors graph and HVG set shaped by the global,
hepatocyte-dominated atlas. We therefore (a) clear the stale graph and
(b) re-derive HVGs + scVI on the T cells alone, so that clustering and the UMAP
reflect T-cell variation rather than global structure. Clustering and UMAP are
both computed on the same scVI latent and are therefore matched.

Expected input AnnData
----------------------
    .layers['counts']   raw counts
    .X                  log-normalized expression
    .obs[celltype_key]  major-type annotation containing ``celltype_val``
    .obs[batch_key]     batch / sample id for integration

Reproducibility
---------------
Seeds are fixed for scVI, PyTorch, and every downstream Scanpy step. The exact
package versions used are logged at run time (``--print-versions``); the versions
this analysis was validated against are pinned in ``requirements.txt``.

Usage
-----
    python tcell_subcluster_scvi.py \
        --input  data/final_rna_wnn.h5ad \
        --output results/lymp_T_scvi_subclustered.h5ad

Run ``python tcell_subcluster_scvi.py --help`` for all options.

Citation
--------
If you use this code, please cite:
    Sarker N. et al. "<paper title>." <journal>, <year>. DOI: <doi>

License: MIT (see LICENSE).
"""

import argparse
import sys

import matplotlib
matplotlib.use("Agg")            # headless: write figures without a display
import matplotlib.pyplot as plt  # noqa: E402
import scanpy as sc              # noqa: E402
import scvi                      # noqa: E402
import torch                     # noqa: E402
import pymde                     # noqa: E402


# --------------------------------------------------------------------------- #
# Manual annotation (dataset-specific)
#
# Cluster -> cell-type map for the chosen Leiden resolution (ANNOTATE_KEY).
# Assignments were made by inspecting the marker dot plot below; the rationale
# for non-obvious calls is noted inline. Editing these three blocks is all that
# is needed to re-annotate at a different resolution.
# --------------------------------------------------------------------------- #
ANNOTATE_KEY = "leiden_scvi_8"   # resolution 0.8 clustering used for annotation

ANNOTATION = {
    "0":  "CD4T",
    "1":  "CD8T",
    "2":  "CD8T",
    "3":  "iNKT",         # confirmed iNKT (Zbtb16+ Klrb1c+), not MAIT
    "4":  "Treg",
    "5":  "neutrophil",   # Cd3+ Csf3r+ granulocyte signal; retained, not merged
    "6":  "NK",
    "7":  "ILC1",
    "8":  "CD8T",
    "9":  "CD8T",
    "10": "gdT",
    "11": "Low-quality",  # all-marker-zero; removed before saving
}

# Marker panel used to validate the annotation (mouse symbols).
MARKERS = {
    "T core":     ["Cd3d", "Cd3e", "Cd3g"],
    "CD8T":       ["Cd8a", "Cd8b1", "Gzmk", "Gzmb"],
    "CD4T":       ["Cd4", "Tcf7", "Lef1"],
    "Treg":       ["Foxp3", "Ikzf2", "Ctla4"],
    "iNKT":       ["Zbtb16", "Klrb1c"],
    "ILC1":       ["Tnfsf10", "Cd200r1", "Tmem176b", "Tbx21"],
    "NK":         ["Klra8", "Klra4", "Zeb2", "Sell", "Ncr1", "Nkg7"],
    "gdT":        ["Trdc"],
    "neutrophil": ["S100a8", "S100a9", "Retnlg", "Csf3r"],
}

# Left-to-right group order for the dot plot / final categories.
GROUP_ORDER = ["CD4T", "CD8T", "Treg", "iNKT", "gdT", "NK", "ILC1", "neutrophil"]

# Labels removed after annotation (kept out of the saved object).
DROP_LABELS = ["Low-quality"]

CELLTYPE_OUT = "celltype_T"      # obs column the annotation is written to


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #
def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="scVI subclustering of the lymphoid-T compartment.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    # I/O
    p.add_argument("--input", required=True,
                   help="Path to the annotated atlas .h5ad (RNA).")
    p.add_argument("--output", required=True,
                   help="Path to write the subclustered .h5ad.")

    # obs / layer keys
    p.add_argument("--celltype-key", default="celltype",
                   help="obs column holding the major-type annotation.")
    p.add_argument("--celltype-val", default="lymp_T",
                   help="Compartment label to subset (reuse for other lineages).")
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
                   help="Apply the manual ANNOTATION map after the sweep.")
    p.add_argument("--annotate-key", default=ANNOTATE_KEY,
                   help="Leiden key whose clusters ANNOTATION maps from.")
    p.add_argument("--dotplot-out", default=None,
                   help="If set, write the marker-validation dot plot here (e.g. .pdf).")

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

    Matches the naming used throughout the analysis and the annotation map
    below. Unique for the 0.1-step sweep used here (resolutions 0.1-1.0 map
    to keys 1-10); if you sweep finer than 0.1 steps, switch to f-string
    'leiden_scvi_{res:g}' and update ANNOTATE_KEY accordingly.
    """
    return f"leiden_scvi_{int(round(res * 10))}"


def annotate_and_validate(adata, annotate_key, dotplot_out=None):
    """Map clusters -> cell types, optionally write a marker dot plot, drop junk.

    Returns a (possibly smaller) AnnData with:
        .obs[CELLTYPE_OUT]  categorical annotation, DROP_LABELS removed
    """
    if annotate_key not in adata.obs:
        raise KeyError(
            f"{annotate_key!r} not in obs. Available leiden keys: "
            f"{[c for c in adata.obs if c.startswith('leiden')]}"
        )

    # unmapped clusters would become NaN and silently vanish -> fail loudly
    clusters = set(adata.obs[annotate_key].astype(str).unique())
    unmapped = clusters - set(ANNOTATION)
    if unmapped:
        raise ValueError(
            f"Clusters in {annotate_key} with no ANNOTATION entry: "
            f"{sorted(unmapped)}"
        )

    adata.obs[CELLTYPE_OUT] = (
        adata.obs[annotate_key].astype(str).map(ANNOTATION).astype("category")
    )
    print(f"\n{CELLTYPE_OUT} (from {annotate_key}):")
    print(adata.obs[CELLTYPE_OUT].value_counts())

    # --- marker validation dot plot on log-norm expression ---
    if dotplot_out is not None:
        present = {k: [g for g in v if g in adata.var_names]
                   for k, v in MARKERS.items()}
        missing = {k: [g for g in v if g not in adata.var_names]
                   for k, v in MARKERS.items()}
        missing = {k: v for k, v in missing.items() if v}
        if missing:
            print(f"[dotplot] markers absent from var_names: {missing}")
        present = {k: v for k, v in present.items() if v}

        order = [g for g in GROUP_ORDER if g in adata.obs[CELLTYPE_OUT].cat.categories]
        sc.pl.dotplot(
            adata, present, groupby=CELLTYPE_OUT, categories_order=order,
            standard_scale="var", dendrogram=False, show=False,
        )
        plt.savefig(dotplot_out, bbox_inches="tight")
        plt.close()
        print(f"[dotplot] saved -> {dotplot_out}")

    # --- drop junk labels, tidy categories ---
    n_before = adata.n_obs
    adata = adata[~adata.obs[CELLTYPE_OUT].isin(DROP_LABELS)].copy()
    adata.obs[CELLTYPE_OUT] = adata.obs[CELLTYPE_OUT].cat.remove_unused_categories()
    order = [g for g in GROUP_ORDER if g in adata.obs[CELLTYPE_OUT].cat.categories]
    adata.obs[CELLTYPE_OUT] = adata.obs[CELLTYPE_OUT].cat.reorder_categories(order)
    print(f"dropped {DROP_LABELS}: {n_before} -> {adata.n_obs} cells")
    return adata


# --------------------------------------------------------------------------- #
# Pipeline
# --------------------------------------------------------------------------- #
def main(argv=None):
    args = parse_args(argv)

    sc.settings.verbosity = 1
    scvi.settings.seed = args.seed
    torch.set_float32_matmul_precision("high")

    log_versions()
    if args.print_versions:
        return

    # --- 1. load + subset to the compartment of interest ---
    adata = sc.read_h5ad(args.input)
    n_before = adata.n_obs
    adata = adata[adata.obs[args.celltype_key] == args.celltype_val].copy()
    print(f"subset {args.celltype_val}: {adata.n_obs} / {n_before} cells")
    if adata.n_obs == 0:
        raise ValueError(
            f"No cells with {args.celltype_key} == {args.celltype_val!r}. "
            f"Available: {sorted(adata.obs[args.celltype_key].unique())}"
        )

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

    # --- 7. manual annotation + marker validation (optional) ---
    if args.annotate:
        adata = annotate_and_validate(
            adata, annotate_key=args.annotate_key, dotplot_out=args.dotplot_out,
        )

    # --- 8. save ---
    adata.write_h5ad(args.output)
    print(f"saved -> {args.output}")


if __name__ == "__main__":
    main()
