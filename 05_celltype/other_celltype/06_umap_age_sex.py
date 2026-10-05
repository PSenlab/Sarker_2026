#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Compartment UMAPs split by age and sex, liver aging multiome (RNA modality).

For each compartment it draws:
  1. Atlas inset: the full-atlas UMAP in grey with the compartment's coarse
     cell types highlighted and labelled                 -> <prefix>_atlas_inset.pdf
  2. Age x sex grid: the compartment UMAP split into male (top row) and
     female (bottom row) x five age groups, cells colored by subcluster
                                                         -> <prefix>_umap_age_sex.pdf

Uses the same COMPARTMENTS keys, h5ad names, label columns and output folders
as 02_composition.py.

Each compartment object is expected to carry:
  .obsm[basis]       2D embedding of the compartment (UMAP from 01_subcluster.py)
  .obs[group_col]    subcluster annotation
  .obs['age'], .obs['sex']
  .uns[f'{group_col}_colors'] (optional) subcluster colors

Usage
-----
    python 06_umap_age_sex.py --parent liver_atlas.h5ad
    python 06_umap_age_sex.py --parent liver_atlas.h5ad --run T_ILC myeloid
    python 06_umap_age_sex.py --skip-inset --run Kupffer
"""

import argparse
import os

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


# --------------------------------------------------------------------------- #
# Dataset conventions
# --------------------------------------------------------------------------- #
AGE_GROUPS = ["young", "mid_age", "old", "pre_geriatric", "geriatric"]
AGE_DISP = {"young": "young", "mid_age": "mid-age", "old": "old",
            "pre_geriatric": "pre-geriatric", "geriatric": "geriatric"}
SEX_ROWS = ["male", "female"]          # male on top, as in the figures
ATLAS_CELLTYPE_COL = "cell_type"

# Coarse cell type colors (same palette as the atlas figures)
CELLTYPE_COLS = {
    "hepatocyte":           "#17becf",
    "endothelial":          "#a6cee3",
    "Kupffer 02":           "#1f78b4",
    "cholangiocyte 01":     "#e31a1c",
    "cholangiocyte 02":     "#cab2d6",
    "Kupffer 01":           "#b2df8a",
    "non-resident myeloid": "#33a02c",
    "mesenchymal":          "#fb9a99",
    "T/ILC cells":          "#fdbf6f",
    "B cells":              "#e377c2",
}


# --------------------------------------------------------------------------- #
# Per-compartment registry (keys, files and labels as in 02_composition.py)
#   coarse : atlas cell types highlighted in the inset
# --------------------------------------------------------------------------- #
COMPARTMENTS = {
    "myeloid": dict(
        h5ad="myeloid.h5ad",
        group_col="celltype_myeloid",
        outdir="path/to/myeloid_DE_csvs",
        prefix="myeloid",
        labels=["Kupffer", "Kupffer cycling", "LAM", "MoMF",
                "cDC1", "pDC", "Neutrophil"],
        coarse=["Kupffer 01", "non-resident myeloid"],
    ),
    "Kupffer": dict(
        h5ad="Kupffer.h5ad",
        group_col="celltype_sub",
        outdir="path/to/Kupffer_DE_csvs",
        prefix="Kupffer",
        labels=["Kupffer", "Kupffer cycling", "LAM", "LSEC-like"],
        coarse=["Kupffer 01", "Kupffer 02"],
    ),
    "endothelial_Kupffer02": dict(
        h5ad="endothelial_Kupffer02.h5ad",
        group_col="endo_subcluster",
        outdir="path/to/endo_DE_csvs",
        prefix="endothelial_Kupffer02",
        labels=["LSEC", "LSEC cycling", "MV portal", "MV central", "Kupffer-like"],
        coarse=["endothelial", "Kupffer 02"],
    ),
    "T_ILC": dict(
        h5ad="T_ILC.h5ad",
        group_col="celltype_T",
        outdir="path/to/T_DE_csvs",
        prefix="T_ILC",
        labels=["CD4T", "CD8T", "ILC1", "NK", "neutrophil", "Treg", "gdT", "iNKT"],
        coarse=["T/ILC cells"],
    ),
}


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #
def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="Compartment UMAPs split by age x sex, plus atlas inset.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--parent", default=None,
                   help="Full liver atlas .h5ad (for the inset). "
                        "Required unless --skip-inset.")
    p.add_argument("--run", nargs="+", default=None, choices=sorted(COMPARTMENTS),
                   help="Which compartments to run (default: all).")
    p.add_argument("--basis", default="X_umap",
                   help="obsm key of the compartment embedding.")
    p.add_argument("--atlas-basis", default=None,
                   help="obsm key of the atlas embedding "
                        "(default: first of X_wnn, X_umap present).")
    p.add_argument("--point-size", type=float, default=3.0)
    p.add_argument("--skip-inset", action="store_true")
    return p.parse_args(argv)


# --------------------------------------------------------------------------- #
# Helpers
# --------------------------------------------------------------------------- #
def subcluster_colors(sub, group_col, labels):
    """Colors stored in the object if they match its categories, else tab20."""
    key = f"{group_col}_colors"
    col = sub.obs[group_col]
    if key in sub.uns and hasattr(col, "cat") and \
            len(sub.uns[key]) == len(col.cat.categories):
        return dict(zip(col.cat.categories.astype(str), sub.uns[key]))
    cmap = plt.get_cmap("tab20")
    return {lab: matplotlib.colors.to_hex(cmap(i % 20)) for i, lab in enumerate(labels)}


def clean_axis(ax):
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)


# --------------------------------------------------------------------------- #
# Plots
# --------------------------------------------------------------------------- #
def plot_atlas_inset(atlas_xy, atlas_ct, cfg, point_size):
    fig, ax = plt.subplots(figsize=(3, 3))
    ax.scatter(atlas_xy[:, 0], atlas_xy[:, 1], s=0.2, c="#d9d9d9",
               linewidths=0, rasterized=True)
    for ct in cfg["coarse"]:
        m = (atlas_ct == ct)
        if not m.any():
            print(f"  [WARN] '{ct}' not found in atlas {ATLAS_CELLTYPE_COL}")
            continue
        ax.scatter(atlas_xy[m, 0], atlas_xy[m, 1], s=1.0,
                   c=CELLTYPE_COLS.get(ct, "#ff7f00"), linewidths=0, rasterized=True)
        cx, cy = np.median(atlas_xy[m, 0]), np.median(atlas_xy[m, 1])
        ax.text(cx, cy, ct, fontsize=8, ha="right", va="center")
    clean_axis(ax)
    fn = os.path.join(cfg["outdir"], f"{cfg['prefix']}_atlas_inset.pdf")
    fig.savefig(fn, bbox_inches="tight", dpi=300)
    plt.close(fig)
    print(f"  saved {os.path.basename(fn)}")


def plot_age_sex_grid(sub, cfg, basis, point_size):
    gc = cfg["group_col"]
    xy = sub.obsm[basis][:, :2]
    colors = subcluster_colors(sub, gc, cfg["labels"])
    lab = sub.obs[gc].astype(str).values
    age = sub.obs["age"].astype(str).values
    sex = sub.obs["sex"].astype(str).values

    pad_x = 0.05 * np.ptp(xy[:, 0]); pad_y = 0.05 * np.ptp(xy[:, 1])
    xlim = (xy[:, 0].min() - pad_x, xy[:, 0].max() + pad_x)
    ylim = (xy[:, 1].min() - pad_y, xy[:, 1].max() + pad_y)

    fig, axes = plt.subplots(len(SEX_ROWS), len(AGE_GROUPS),
                             figsize=(2.2 * len(AGE_GROUPS) + 1.8, 2.2 * len(SEX_ROWS)),
                             squeeze=False)
    for r, sx in enumerate(SEX_ROWS):
        for c, ag in enumerate(AGE_GROUPS):
            ax = axes[r, c]
            m = (sex == sx) & (age == ag)
            for l in cfg["labels"]:           # draw in label order
                mm = m & (lab == l)
                if mm.any():
                    ax.scatter(xy[mm, 0], xy[mm, 1], s=point_size, c=colors.get(l, "grey"),
                               linewidths=0, rasterized=True)
            ax.set_xlim(xlim); ax.set_ylim(ylim)
            clean_axis(ax)
            if r == 0:
                ax.set_title(AGE_DISP[ag], fontsize=10, fontweight="bold")
            if c == 0:
                ax.set_ylabel(sx, fontsize=10, fontweight="bold")

    handles = [Line2D([0], [0], marker="o", ls="", markersize=6,
                      markerfacecolor=colors.get(l, "grey"), markeredgewidth=0, label=l)
               for l in cfg["labels"]]
    fig.legend(handles=handles, loc="center left", bbox_to_anchor=(1.0, 0.5),
               frameon=False, fontsize=9)
    plt.tight_layout()
    fn = os.path.join(cfg["outdir"], f"{cfg['prefix']}_umap_age_sex.pdf")
    fig.savefig(fn, bbox_inches="tight", dpi=300)
    plt.close(fig)
    print(f"  saved {os.path.basename(fn)}")


# --------------------------------------------------------------------------- #
# Main
# --------------------------------------------------------------------------- #
def main(argv=None):
    args = parse_args(argv)
    to_run = args.run or list(COMPARTMENTS)

    atlas_xy = atlas_ct = None
    if not args.skip_inset:
        if args.parent is None:
            raise ValueError("--parent is required (or pass --skip-inset).")
        print(f"loading parent atlas (obs + embedding only): {args.parent}")
        atlas = sc.read_h5ad(args.parent, backed="r")
        basis = args.atlas_basis or next(
            (k for k in ("X_wnn", "X_umap") if k in atlas.obsm), None)
        if basis is None:
            raise KeyError(f"no atlas embedding found (obsm: {list(atlas.obsm.keys())})")
        atlas_xy = np.asarray(atlas.obsm[basis])[:, :2]
        atlas_ct = atlas.obs[ATLAS_CELLTYPE_COL].astype(str).values
        print(f"  {atlas.n_obs} cells, embedding '{basis}'")
        atlas.file.close()

    for name in to_run:
        cfg = COMPARTMENTS[name]
        print(f"\n===== {name} =====")
        os.makedirs(cfg["outdir"], exist_ok=True)

        if not args.skip_inset:
            plot_atlas_inset(atlas_xy, atlas_ct, cfg, args.point_size)

        sub = sc.read_h5ad(cfg["h5ad"])
        for need in (cfg["group_col"], "age", "sex"):
            if need not in sub.obs:
                raise KeyError(f"{need!r} not in {cfg['h5ad']} .obs")
        if args.basis not in sub.obsm:
            raise KeyError(f"{args.basis!r} not in {cfg['h5ad']} .obsm "
                           f"(have: {list(sub.obsm.keys())})")
        missing = sorted(set(sub.obs[cfg["group_col"]].dropna().astype(str))
                         - set(cfg["labels"]))
        if missing:
            print(f"  [WARN] labels in object but not in config (not drawn): {missing}")
        plot_age_sex_grid(sub, cfg, args.basis, args.point_size)
        del sub

    print("\ndone.")


if __name__ == "__main__":
    main()
