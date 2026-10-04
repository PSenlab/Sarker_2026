#!/usr/bin/env python
# coding: utf-8
# =============================================================================
# SCENIC+ Downstream: eRegulon Dot Plots for the Atlas and Non-Hepatocyte Compartments
# =============================================================================
#
# Description:
#   For each compartment, draws two SCENIC+ heatmap dot plots of the direct
#   (+/+) eRegulons:
#     1. by subcluster
#     2. by age (young -> geriatric)
#   Color = direct gene-based AUC, size = direct region-based AUC (same
#   settings as 07_downstream_hepatocyte.py).
#
#   Compartments: all_celltypes, T_ILC, endothelial_Kupffer02, myeloid
#   (the hepatocyte compartment is handled by 07_downstream_hepatocyte.py)
#
# Input (per compartment, from Step 6):
#   <out_dir>/scenicplus/outs/scplusmdata.h5mu
#
# Output (per compartment, in <out_dir>/scenicplus/figures/):
#   <suffix>_subcluster_dotplot.pdf
#   <suffix>_age_dotplot.pdf
#
# Usage:
#   python 09_downstream_other_compartments.py
#   python 09_downstream_other_compartments.py --run all_celltypes T_ILC
#
# =============================================================================

import argparse
import os
import traceback

import mudata
from scenicplus.plotting.dotplot import heatmap_dotplot


# =============================================================================
# CONFIGURATION
# =============================================================================

# Base directory holding the outs_* trees (Steps 1-6). "." if you run in place.
BASE_DIR = "."

AGE_ORDER = ["young", "mid_age", "old", "pre_geriatric", "geriatric"]

# out_dir / suffix must match Steps 1-5.
# group_col: subcluster column in the compartment RNA object (Step 4 input);
#            it appears in the SCENIC+ mudata as "scRNA_counts:<group_col>".
# order:     column order of the subcluster dot plot (None = alphabetical)
COMPARTMENTS = {
    "all_celltypes": dict(
        out_dir="outs_all_celltypes",
        suffix="all_celltypes",
        group_col="cell_type",
        order=["cholangiocyte 01", "cholangiocyte 02", "endothelial", "Kupffer 02",
               "hepatocyte", "Kupffer 01", "non-resident myeloid", "mesenchymal",
               "B cells", "T/ILC cells"],
    ),
    "T_ILC": dict(
        out_dir="outs_T_ILC",
        suffix="T_ILC",
        group_col="celltype_T",
    ),
    "endothelial_Kupffer02": dict(
        out_dir="outs_endothelial_Kupffer02",
        suffix="endothelial_Kupffer02",
        group_col="endo_subcluster",
    ),
    "myeloid": dict(
        out_dir="outs_myeloid",
        suffix="myeloid",
        group_col="celltype_myeloid",
    ),
}


# =============================================================================
# HELPERS
# =============================================================================

def banner(text):
    print("\n" + "=" * 70)
    print(text)
    print("=" * 70)


def find_group_column(obs_columns, group_col):
    """Subcluster column in the mudata obs: RNA label first, ATAC 'celltype'
    (from the Step 1 cell_data TSV) as fallback."""
    for col in (f"scRNA_counts:{group_col}",
                "scRNA_counts:celltype",
                "scATAC_counts:celltype",
                "scATAC_counts:cell_type"):
        if col in obs_columns:
            return col
    raise KeyError(f"no subcluster column found (looked for "
                   f"'scRNA_counts:{group_col}', 'scRNA_counts:celltype', "
                   f"'scATAC_counts:celltype')")


def dotplot_size(n_eregulons, n_groups):
    """Figure size that scales with the number of eRegulons and groups."""
    return (max(5, 2.5 + 0.45 * n_groups), max(6, 0.28 * n_eregulons))


def draw(scplus_mdata, group_variable, n_eregulons, out_pdf, order=None):
    groups = scplus_mdata.obs[group_variable].dropna().unique()
    kwargs = dict(
        scplus_mudata=scplus_mdata,
        color_modality="direct_gene_based_AUC",
        size_modality="direct_region_based_AUC",
        group_variable=group_variable,
        eRegulon_metadata_key="direct_e_regulon_metadata",
        color_feature_key="Gene_signature_name",
        size_feature_key="Region_signature_name",
        feature_name_key="eRegulon_name",
        sort_data_by="direct_gene_based_AUC",
        orientation="vertical",
        figsize=dotplot_size(n_eregulons, len(groups)),
        save=out_pdf,
    )
    if order is not None:
        extra = sorted(set(groups) - set(order))
        if extra:
            print(f"  [WARN] groups not in the configured order (appended): {extra}")
        kwargs["group_variable_order"] = [g for g in order if g in set(groups)] + extra
    heatmap_dotplot(**kwargs)
    print(f"  [OK] {out_pdf}")


# =============================================================================
# PER-COMPARTMENT RUNNER
# =============================================================================

def run_compartment(name, spec, base_dir):
    banner(f"COMPARTMENT: {name}")

    sp_dir = os.path.join(base_dir, spec["out_dir"], "scenicplus")
    mdata_path = os.path.join(sp_dir, "outs", "scplusmdata.h5mu")
    fig_dir = os.path.join(sp_dir, "figures")
    os.makedirs(fig_dir, exist_ok=True)

    if not os.path.exists(mdata_path):
        raise FileNotFoundError(f"[ERROR] SCENIC+ output not found: {mdata_path}")

    scplus_mdata = mudata.read(mdata_path)
    print(f"  [OK] loaded {mdata_path}: {scplus_mdata.n_obs} cells")

    # direct +/+ eRegulons only (as in 07_downstream_hepatocyte.py)
    meta = scplus_mdata.uns["direct_e_regulon_metadata"]
    meta = meta[meta["eRegulon_name"].str.endswith("direct_+/+")]
    scplus_mdata.uns["direct_e_regulon_metadata"] = meta
    n_eregulons = meta["eRegulon_name"].nunique()
    print(f"  [OK] direct +/+ eRegulons: {n_eregulons}")

    # 1. by subcluster
    group_variable = find_group_column(scplus_mdata.obs.columns, spec["group_col"])
    print(f"  subcluster column: {group_variable}")
    draw(scplus_mdata, group_variable, n_eregulons,
         os.path.join(fig_dir, f"{spec['suffix']}_subcluster_dotplot.pdf"),
         order=spec.get("order"))

    # 2. by age
    draw(scplus_mdata, "scRNA_counts:age", n_eregulons,
         os.path.join(fig_dir, f"{spec['suffix']}_age_dotplot.pdf"),
         order=AGE_ORDER)


# =============================================================================
# MAIN
# =============================================================================

def main():
    p = argparse.ArgumentParser(
        description="SCENIC+ eRegulon dot plots for the non-hepatocyte compartments")
    p.add_argument("--run", nargs="+", default=None, choices=sorted(COMPARTMENTS),
                   help="Which compartments (default: all)")
    p.add_argument("--base-dir", default=BASE_DIR,
                   help="Directory containing the outs_* trees")
    args = p.parse_args()

    results = {}
    for name in args.run or list(COMPARTMENTS):
        try:
            run_compartment(name, COMPARTMENTS[name], args.base_dir)
            results[name] = "OK"
        except Exception as e:
            print(f"\n[ERROR] {name} failed: {e}")
            traceback.print_exc()
            results[name] = "FAILED"

    banner("DONE")
    for k, v in results.items():
        print(f"   {k}: {v}")


if __name__ == "__main__":
    main()
