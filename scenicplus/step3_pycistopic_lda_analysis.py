#!/usr/bin/env python3
"""
SCENIC+ Step 4: RNA-ATAC Integration Prep — ALL COMPARTMENTS
=============================================================
Aligns matched RNA (GEX) and ATAC (cisTopic) data per compartment: reconstructs
barcodes, intersects cells, orders both objects identically, and writes the
SCENIC+ input bundle.

Consumes
--------
  <out_dir>/cisTopicObject_lda_complete.pkl     (Step 3)
  <out_dir>/region_sets/                        (Step 3)
  <rna_h5ad>                                    (compartment RNA object)

Produces
--------
  <out_dir>/scenicplus_input/<suffix>_GEX_anndata.h5ad
  <out_dir>/scenicplus_input/<suffix>_cisTopic_obj.pkl
  <out_dir>/scenicplus_input/region_sets/

NOTE: this duplicates step 9 of step3_lda_all.py. Keep ONE as canonical -- the
recommended split is to run Step 3 with `--skip 9` and use this script, so the
integration can be redone (different RNA object, different filter) without
repeating LDA.

Usage
-----
  python step4_rna_atac_prep_all.py                  # all compartments
  python step4_rna_atac_prep_all.py --run T_ILC
  python step4_rna_atac_prep_all.py --run myeloid --rna other.h5ad
  python step4_rna_atac_prep_all.py --dry-run        # report overlap, write nothing

Author: Nishat Sarker
License: MIT
"""

import argparse
import os
import pickle
import re
import shutil
import traceback
import warnings
from typing import List, Optional, Tuple

import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc

warnings.filterwarnings("ignore")


# =============================================================================
# SHARED CONFIGURATION
# =============================================================================

DESIRED_AGE_ORDER = ["young", "mid_age", "old", "pre_geriatric", "geriatric"]
TARGET_SUM = 1e4

# RNA obs_names spell some sample groups differently from the cisTopic cell
# names (RNA 'pre_geriatric_3' vs cisTopic 'pre_ger_03', because Step 1
# abbreviated it when building sample ids). Applied in transform_barcode_index().
# If one age group shows zero shared cells, the fix belongs here.
SAMPLE_RENAME = {"pre_geriatric": "pre_ger"}


# =============================================================================
# PER-COMPARTMENT REGISTRY  (out_dir / suffix must match Steps 1-3)
#   celltype_col     : annotation column ('cell_type' for the full-atlas run)
#   cell_type_filter : None = no pre-filter; the RNA<->ATAC intersection already
#                      defines the compartment when rna_h5ad is the subset object
# =============================================================================

COMPARTMENTS = {

    "all_celltypes": dict(
        out_dir="outs_all_celltypes",
        suffix="all_celltypes",
        rna_h5ad="liver_atlas.h5ad",
        celltype_col="cell_type",
        cell_type_filter=None,
    ),

    "Hepatocyte": dict(
        out_dir="outs_Heps",
        suffix="Heps",
        rna_h5ad="Hepatocyte.h5ad",
        celltype_col="celltype",
        cell_type_filter=None,
    ),

    "endothelial_Kupffer02": dict(
        out_dir="outs_endothelial_Kupffer02",
        suffix="endothelial_Kupffer02",
        rna_h5ad="endothelial_Kupffer02.h5ad",
        celltype_col="celltype",
        cell_type_filter=None,
    ),

    "myeloid": dict(
        out_dir="outs_myeloid",
        suffix="myeloid",
        rna_h5ad="myeloid.h5ad",
        celltype_col="celltype",
        cell_type_filter=None,
    ),

    "T_ILC": dict(
        out_dir="outs_T_ILC",
        suffix="T_ILC",
        rna_h5ad="T_ILC.h5ad",
        celltype_col="celltype",
        cell_type_filter=None,
    ),
}


# =============================================================================
# CONFIG OBJECT
# =============================================================================

class Cfg:
    def __init__(self, name, spec, rna_override=None):
        self.name = name
        self.outDir = spec["out_dir"]
        self.suffix = spec["suffix"]
        self.rna_h5ad = rna_override or spec["rna_h5ad"]
        self.celltype_col = spec.get("celltype_col", "celltype")
        self.cell_type_filter = spec.get("cell_type_filter")

        self.cistopic_pkl = os.path.join(self.outDir, "cisTopicObject_lda_complete.pkl")
        self.region_sets_dir = os.path.join(self.outDir, "region_sets")
        self.scenicplus_dir = os.path.join(self.outDir, "scenicplus_input")
        self.gex_out = os.path.join(self.scenicplus_dir,
                                    f"{self.suffix}_GEX_anndata.h5ad")
        self.ctx_out = os.path.join(self.scenicplus_dir,
                                    f"{self.suffix}_cisTopic_obj.pkl")


# =============================================================================
# HELPERS
# =============================================================================

def transform_barcode_index(index: str) -> str:
    """RNA obs_name -> cisTopic cell name.

    Handles BOTH 'ACGT-geriatric_3' and 'ACGT-1-geriatric_3'. The previous
    rsplit('-', 1) implementation turned the second form into
    'ACGT-1-1-geriatric_03___geriatric_03', which matches nothing -- the failure
    is silent because it just lowers the shared-cell count.

        ACGT-geriatric_3        -> ACGT-1-geriatric_03___geriatric_03
        ACGT-1-geriatric_3      -> ACGT-1-geriatric_03___geriatric_03
        ACGT-pre_geriatric_3    -> ACGT-1-pre_ger_03___pre_ger_03
    """
    s = str(index).split("___")[0]
    parts = s.split("-")
    dna = parts[0]
    rest = [p for p in parts[1:] if p != "1"]
    if not rest:
        raise ValueError(f"cannot parse sample from index: {index}")
    sample = "-".join(rest)
    for old, new in SAMPLE_RENAME.items():
        sample = re.sub(rf"^{old}", new, sample)
    sample = re.sub(r"_(\d+)$", lambda m: f"_{int(m.group(1)):02d}", sample)
    return f"{dna}-1-{sample}___{sample}"


def model_topic_count(model):
    for attr in ("n_topic", "topic_no", "n_topics"):
        if hasattr(model, attr):
            return int(getattr(model, attr))
    return None


def load_rna_anndata(filepath: str, celltype_col: str = "celltype",
                     cell_type: Optional[str] = None) -> ad.AnnData:
    """Load RNA object, preferring RAW COUNTS from layers['counts'].

    In the subclustered objects .X is ALREADY log-normalized (raw counts live in
    layers['counts']). Reading .X and then running normalize_total + log1p would
    double-transform, and adata.raw would hold lognorm despite the 'raw counts'
    label. Reads via h5py so a .uns entry an older anndata cannot decode (e.g. a
    stored dendrogram) does not break the load.
    """
    import h5py
    from anndata import AnnData
    try:
        from anndata.experimental import read_elem
    except ImportError:
        from anndata._io.specs import read_elem

    print(f" Loading RNA AnnData from {filepath}...")
    with h5py.File(filepath, "r") as f:
        has_counts = "layers" in f and "counts" in f["layers"]
        X = read_elem(f["layers"]["counts"]) if has_counts else read_elem(f["X"])
        adata = AnnData(X=X, obs=read_elem(f["obs"]), var=read_elem(f["var"]))

    src = "layers['counts'] (raw)" if has_counts else ".X"
    print(f"   [OK] {adata.n_obs} cells x {adata.n_vars} genes  [matrix: {src}]")
    if not has_counts:
        print("   [WARN] no layers['counts'] -- if .X is already log-normalized, "
              "the normalize/log1p step below will DOUBLE-TRANSFORM it. Verify.")

    if cell_type is not None:
        col = celltype_col if celltype_col in adata.obs.columns else "celltype"
        if col not in adata.obs.columns:
            raise ValueError(f"'{celltype_col}' not found in .obs "
                             f"(have: {adata.obs.columns.tolist()})")
        print(f"\n Filtering to {cell_type} cells (column '{col}')...")
        adata = adata[adata.obs[col] == cell_type].copy()
        print(f"   [OK] {adata.n_obs} cells retained")

    return adata


def load_cistopic_object(filepath: str):
    print(f"\n Loading cisTopic object from {filepath}...")
    if not os.path.exists(filepath):
        raise FileNotFoundError(
            f"[ERROR] Step 3 output missing: {filepath}\n"
            f"        Run step3_lda_all.py for this compartment first.")
    with open(filepath, "rb") as f:
        obj = pickle.load(f)
    n_regions = len(obj.region_names) if hasattr(obj, "region_names") else "N/A"
    print(f"   [OK] {len(obj.cell_names)} cells x {n_regions} regions")
    return obj


def preprocess_rna(adata: ad.AnnData, target_sum: float = 1e4) -> ad.AnnData:
    print("\n Preprocessing RNA data...")
    adata.raw = adata.copy()
    print("   [OK] raw counts stored in adata.raw")
    sc.pp.normalize_total(adata, target_sum=target_sum)
    print(f"   [OK] normalized to target sum {target_sum:.0e}")
    sc.pp.log1p(adata)
    print("   [OK] log1p")
    return adata


def transform_rna_indices(adata: ad.AnnData) -> ad.AnnData:
    print("\n Transforming RNA barcode indices...")
    original = adata.obs_names.tolist()
    transformed, failed = [], []
    for idx in original:
        try:
            transformed.append(transform_barcode_index(idx))
        except ValueError:
            failed.append(idx)
            transformed.append(idx)
    if failed:
        print(f"   [WARN] {len(failed)} indices could not be transformed")
        print(f"      examples: {failed[:3]}")
    adata.obs_names = pd.Index(transformed)
    print(f"   [OK] {len(transformed)} indices")
    print(f"      example: {original[0]} -> {transformed[0]}")
    return adata


def remove_duplicates(adata, cistopic_obj):
    print("\n Removing duplicate indices...")
    rna_dups = adata.obs_names.duplicated()
    if rna_dups.sum():
        print(f"   [WARN] {rna_dups.sum()} duplicate RNA indices")
        adata = adata[~rna_dups].copy()
    ctx_dups = cistopic_obj.cell_data.index.duplicated()
    if ctx_dups.sum():
        print(f"   [WARN] {ctx_dups.sum()} duplicate cisTopic indices")
        cistopic_obj.cell_data = cistopic_obj.cell_data[~ctx_dups]
    if not rna_dups.sum() and not ctx_dups.sum():
        print("   [OK] none found")
    return adata, cistopic_obj


def find_matching_cells(adata, cistopic_obj) -> List[str]:
    """Return SORTED shared indices.

    Sorted, not set-ordered: iterating a set gives an arbitrary (and between-run
    variable) order, which would make the written cell order non-reproducible
    even though the RNA/ATAC alignment itself stayed correct.
    """
    print("\n Finding matching cells between RNA and ATAC...")
    rna_idx = set(adata.obs_names)
    atac_idx = set(cistopic_obj.cell_data.index)
    shared = sorted(rna_idx & atac_idx)

    print(f"   [STATS] RNA cells:    {len(rna_idx)}")
    print(f"   [STATS] ATAC cells:   {len(atac_idx)}")
    print(f"   [OK]    shared:       {len(shared)}")
    print(f"           RNA only:     {len(rna_idx - atac_idx)}")
    print(f"           ATAC only:    {len(atac_idx - rna_idx)}")

    if len(shared) == 0:
        raise ValueError(
            "No matching cells between RNA and ATAC -- barcode format mismatch.\n"
            f"  RNA example:  {sorted(rna_idx)[:1]}\n"
            f"  ATAC example: {sorted(atac_idx)[:1]}\n"
            "  Check transform_barcode_index() and SAMPLE_RENAME.")

    # Per-sample-group breakdown. A global count can look healthy while ONE age
    # group silently contributes zero (the classic pre_geriatric/pre_ger case).
    def group_of(x):
        return x.split("___")[-1].rsplit("_", 1)[0]

    shared_grp = pd.Series([group_of(s) for s in shared]).value_counts()
    atac_grp = pd.Series([group_of(s) for s in atac_idx]).value_counts()
    print("\n   shared cells per sample group (shared / ATAC):")
    for g in sorted(atac_grp.index):
        s = int(shared_grp.get(g, 0))
        flag = "   <-- ZERO, check SAMPLE_RENAME" if s == 0 else ""
        print(f"      {g:16s} {s:7d} / {int(atac_grp[g]):7d}{flag}")

    if len(rna_idx - atac_idx):
        print(f"\n   examples missing from ATAC: {sorted(rna_idx - atac_idx)[:2]}")
    if len(atac_idx - rna_idx):
        print(f"   examples missing from RNA:  {sorted(atac_idx - rna_idx)[:2]}")

    return shared


def subset_to_shared_cells(adata, cistopic_obj, shared: List[str]):
    print("\n Subsetting to shared cells and aligning...")
    adata = adata[shared, :].copy()
    cistopic_obj.cell_data = cistopic_obj.cell_data.loc[adata.obs_names]
    assert list(adata.obs_names) == list(cistopic_obj.cell_data.index), \
        "Index alignment failed!"
    print(f"   [OK] aligned: {adata.n_obs} cells")
    return adata, cistopic_obj


def set_categorical_age(adata, cistopic_obj, age_order: List[str]):
    print("\n Setting age as ordered categorical...")
    if "age" in adata.obs.columns:
        valid = [a for a in age_order if a in set(adata.obs["age"].unique())]
        adata.obs["age"] = pd.Categorical(adata.obs["age"], categories=valid,
                                          ordered=True)
        print(f"   [OK] RNA: {valid}")
    else:
        print("   [WARN] 'age' not in RNA .obs")

    if "age" in cistopic_obj.cell_data.columns:
        valid = [a for a in age_order
                 if a in set(cistopic_obj.cell_data["age"].unique())]
        cistopic_obj.cell_data["age"] = pd.Categorical(
            cistopic_obj.cell_data["age"], categories=valid, ordered=True)
        print(f"   [OK] cisTopic: {valid}")
    else:
        print("   [WARN] 'age' not in cisTopic cell_data")
    return adata, cistopic_obj


def copy_region_sets(cfg):
    """SCENIC+ reads the region sets from the input bundle, so copy them in."""
    if not os.path.isdir(cfg.region_sets_dir):
        print(f"   [WARN] region_sets not found at {cfg.region_sets_dir}")
        return
    dst = os.path.join(cfg.scenicplus_dir, "region_sets")
    if os.path.exists(dst):
        shutil.rmtree(dst)
    shutil.copytree(cfg.region_sets_dir, dst)
    print(f"   [OK] region sets -> {dst}")
    for sub in sorted(os.listdir(dst)):
        p = os.path.join(dst, sub)
        if os.path.isdir(p):
            n = len([f for f in os.listdir(p) if f.endswith(".bed")])
            print(f"      {sub:20s} {n} bed")


def save_outputs(cfg, adata, cistopic_obj):
    print("\n[SAVE] Saving SCENIC+ input bundle...")
    os.makedirs(cfg.scenicplus_dir, exist_ok=True)
    adata.write(cfg.gex_out)
    print(f"   [OK] {cfg.gex_out}")
    with open(cfg.ctx_out, "wb") as f:
        pickle.dump(cistopic_obj, f)
    print(f"   [OK] {cfg.ctx_out}")
    copy_region_sets(cfg)


def print_final_summary(cfg, adata, cistopic_obj):
    print("\n" + "=" * 60)
    print(f"[STATS] FINAL SUMMARY — {cfg.name}")
    print("=" * 60)
    print(f"\n GEX AnnData: {adata.n_obs} cells x {adata.n_vars} genes")
    print(f"   obs columns: {adata.obs.columns.tolist()}")
    if "age" in adata.obs.columns:
        print("   age distribution:")
        for age, count in adata.obs["age"].value_counts().sort_index().items():
            print(f"      {age}: {count}")
    print(f"\n cisTopic: {len(cistopic_obj.cell_data)} cells")
    if hasattr(cistopic_obj, "region_names"):
        print(f"   regions: {len(cistopic_obj.region_names)}")
    if getattr(cistopic_obj, "selected_model", None) is not None:
        print(f"   topics: {model_topic_count(cistopic_obj.selected_model)}")
    print(f"   cell_data columns: {cistopic_obj.cell_data.columns.tolist()}")
    print(f"\n[OK] aligned index example: {list(adata.obs_names[:2])}")


# =============================================================================
# PER-COMPARTMENT DRIVER
# =============================================================================

def run_compartment(cfg: Cfg, dry_run: bool = False):
    print("\n" + "#" * 70)
    print(f"#  COMPARTMENT: {cfg.name}  ->  {cfg.scenicplus_dir}")
    print("#" * 70)

    if not os.path.exists(cfg.rna_h5ad):
        raise FileNotFoundError(f"[ERROR] RNA object not found: {cfg.rna_h5ad}")

    adata = load_rna_anndata(cfg.rna_h5ad, celltype_col=cfg.celltype_col,
                             cell_type=cfg.cell_type_filter)
    adata = preprocess_rna(adata, target_sum=TARGET_SUM)
    adata = transform_rna_indices(adata)

    cistopic_obj = load_cistopic_object(cfg.cistopic_pkl)
    adata, cistopic_obj = remove_duplicates(adata, cistopic_obj)

    shared = find_matching_cells(adata, cistopic_obj)

    if dry_run:
        print("\n[DRY RUN] overlap reported; nothing written.")
        return None, None

    adata, cistopic_obj = subset_to_shared_cells(adata, cistopic_obj, shared)
    adata, cistopic_obj = set_categorical_age(adata, cistopic_obj, DESIRED_AGE_ORDER)
    save_outputs(cfg, adata, cistopic_obj)
    print_final_summary(cfg, adata, cistopic_obj)
    return adata, cistopic_obj


# =============================================================================
# MAIN
# =============================================================================

def main():
    p = argparse.ArgumentParser(
        description="Prepare RNA + ATAC data for SCENIC+ (all compartments)")
    p.add_argument("--run", nargs="+", default=None, choices=sorted(COMPARTMENTS),
                   help="Which compartments to prepare (default: all)")
    p.add_argument("--rna", type=str, default=None,
                   help="Override the RNA .h5ad (single compartment only)")
    p.add_argument("--dry-run", action="store_true",
                   help="Report RNA/ATAC overlap and exit without writing")
    p.add_argument("--stop-on-error", action="store_true",
                   help="Abort on first failure (default: continue)")
    args = p.parse_args()

    names = args.run or list(COMPARTMENTS)
    if args.rna and len(names) > 1:
        p.error("--rna applies to a single compartment; "
                "set rna_h5ad in COMPARTMENTS for multi-compartment runs")

    print("\n" + "=" * 60)
    print("RNA-ATAC INTEGRATION PREPARATION FOR SCENIC+")
    print(f"Compartments: {names}")
    print("=" * 60)

    results = {}
    for name in names:
        try:
            run_compartment(Cfg(name, COMPARTMENTS[name], rna_override=args.rna),
                            dry_run=args.dry_run)
            results[name] = "OK"
        except Exception as e:
            print(f"\n[ERROR] {name} failed: {e}")
            traceback.print_exc()
            if args.stop_on_error:
                raise
            print("[INFO] Continuing to next compartment...\n")
            results[name] = "FAILED"

    print("\n" + "=" * 60)
    print("STEP 4 FINISHED")
    for k, v in results.items():
        print(f"   {k}: {v}")
    print("=" * 60)


if __name__ == "__main__":
    main()
