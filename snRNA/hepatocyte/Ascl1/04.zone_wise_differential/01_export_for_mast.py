#!/usr/bin/env python3
"""
01_export_for_mast.py
============================================================
Export log-normalized expression + Ascl1 status for MAST DE (run in R).

Pipeline:
  load adata -> subset female hepatocytes -> map zonation
  -> binarize Ascl1 (detection) -> export global + each zone

Outputs (per subset {allfemhep, periportal, midlobular, pericentral}):
  mast_export/<subset>_expr.mtx     genes x cells, log-normalized
  mast_export/<subset>_genes.csv    gene order
  mast_export/<subset>_meta.csv     barcode, Ascl1_status

Note: .X must be the LOG-NORMALIZED layer (same as used with use_raw=False).
"""
import os
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.io as sio
import scipy.sparse as sp

# ------------------------- config -------------------------
ADATA_PATH   = "adata_female_hepatocytes.h5ad"   # <-- set to your object
SEX_COL      = "sex"
CELLTYPE_COL = "celltype"
CT2_COL      = "celltype2"        # hepatocyte sub-clusters Hep-01..Hep-07
ZONE_COL     = "Zonation"
GENE         = "Ascl1"
EXPORT_DIR   = "mast_export"
MIN_FRAC     = 0.01               # keep genes detected in >=1% of cells per subset

HEP_ZONATION_MAP = {
    "Hep-01": "Periportal", "Hep-02": "Midlobular", "Hep-03": "Pericentral",
    "Hep-04": "Pericentral", "Hep-05": "Periportal", "Hep-06": "Midlobular",
    "Hep-07": "Midlobular",
}
ZONE_ORDER = ["Periportal", "Midlobular", "Pericentral"]
# ----------------------------------------------------------

os.makedirs(EXPORT_DIR, exist_ok=True)
adata = sc.read_h5ad(ADATA_PATH)

# female hepatocytes (no-op if the object is already subset)
mask = (adata.obs[SEX_COL] == "female") & (adata.obs[CELLTYPE_COL] == "Hepatocyte")
adata_fh = adata[mask].copy() if mask.any() else adata.copy()

# zonation
if ZONE_COL not in adata_fh.obs.columns:
    adata_fh.obs[ZONE_COL] = adata_fh.obs[CT2_COL].map(HEP_ZONATION_MAP)

# Ascl1 status: detection-based (> 0), identical to the submitted analysis
a = adata_fh[:, GENE].X
a = a.toarray().ravel() if sp.issparse(a) else np.asarray(a).ravel()
adata_fh.obs["Ascl1_status"] = np.where(a > 0, "Ascl1_pos", "Ascl1_neg")


def export_subset(sub, tag):
    X = sub.X
    X = X.tocsr() if sp.issparse(X) else sp.csr_matrix(X)
    if MIN_FRAC > 0:
        det  = np.asarray((X > 0).sum(axis=0)).ravel()
        keep = det >= (MIN_FRAC * X.shape[0])
        X, genes = X[:, keep], sub.var_names[keep]
    else:
        genes = sub.var_names
    sio.mmwrite(f"{EXPORT_DIR}/{tag}_expr.mtx", X.T.tocsr())     # genes x cells
    pd.DataFrame({"gene": genes}).to_csv(f"{EXPORT_DIR}/{tag}_genes.csv", index=False)
    pd.DataFrame({"barcode": sub.obs_names,
                  "Ascl1_status": sub.obs["Ascl1_status"].values}
                 ).to_csv(f"{EXPORT_DIR}/{tag}_meta.csv", index=False)
    npos = int((sub.obs["Ascl1_status"] == "Ascl1_pos").sum())
    print(f"{tag:12s}: {sub.n_obs:>7,} cells x {len(genes):>6,} genes  (Ascl1+: {npos:,})")


# global (all female hepatocytes)
export_subset(adata_fh, "allfemhep")

# zone-wise
for zone in ZONE_ORDER:
    sub = adata_fh[adata_fh.obs[ZONE_COL] == zone]
    export_subset(sub, zone.lower())

print("\nDone ->", os.path.abspath(EXPORT_DIR))
