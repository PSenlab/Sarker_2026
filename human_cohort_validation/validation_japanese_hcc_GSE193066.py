#!/usr/bin/env python3
# ==============================================================================
# Japanese MASLD-HCC Risk Validation - GSE193066
# ==============================================================================
#
# Description:
#   Validates Ascl1 sex dimorphism in the Japanese MASLD-HCC progression
#   cohort (serial biopsies). Expression is distributed as a pre-normalized
#   GCT file with gene symbols in the index. Sex is annotated in the
#   metadata as "male" / "female".
#
# Dataset-specific quirks:
#   - GCT format (not series matrix, not per-sample counts)
#   - Pre-normalized values (log2(raw + 0.1) applied for ASCL1)
#   - Sample matching via Title column (NOT GSM)
#   - Sex annotated as "male" / "female" (lowercase)
#   - Pseudocount 0.1 instead of 1
#   - Both biopsies combined (n=164)
#
# Inputs:
#   GSE193066_raw/GSE193066_NAFLD.HUn164.gct
#   GSE193066_raw/GSE193066_series_matrix.txt.gz
#
# Output:
#   verify_GSE193066_ASCL1_barplot.{png,pdf}
#
# ==============================================================================

import os
import numpy as np
import pandas as pd

from validation_utils import (
    banner, parse_series_matrix, mann_whitney_sex, barplot_male_vs_female,
)

# ==============================================================================
# CONFIGURATION - UPDATE THESE PATHS
# ==============================================================================
DATA_DIR    = "GSE193066_raw"
GCT_FILE    = os.path.join(DATA_DIR, "GSE193066_NAFLD.HUn164.gct")
SERIES_FILE = os.path.join(DATA_DIR, "GSE193066_series_matrix.txt.gz")


# ==============================================================================
# STEP 1: Load GCT expression file
# ==============================================================================
banner("STEP 1: Load GCT expression file")

with open(GCT_FILE) as f:
    version = f.readline().strip()
    dims = f.readline().strip().split("\t")
    print(f"  GCT version: {version}")
    print(f"  Dimensions: {dims[0]} genes x {dims[1]} samples")

expr = pd.read_csv(GCT_FILE, sep="\t", skiprows=2, index_col=0)
if "Description" in expr.columns:
    expr = expr.drop(columns="Description")
non_numeric = [c for c in expr.columns if not np.issubdtype(expr[c].dtype, np.number)]
if non_numeric:
    expr = expr.drop(columns=non_numeric)
expr = expr.apply(pd.to_numeric, errors="coerce")
print(f"  Expression: {expr.shape}")
print(f"  Value range: {expr.min().min():.2f} - {expr.max().max():.2f}")


# ==============================================================================
# STEP 2: Metadata
# ==============================================================================
banner("STEP 2: Metadata extraction")

_, meta, _ = parse_series_matrix(SERIES_FILE, include_expr=False)
print(f"  Samples in metadata: {len(meta)}")
print(f"  Columns: {[c for c in meta.columns if c not in ['GSM', 'Title']]}")
if "Sex" in meta.columns:
    print(f"  [OK] Sex annotated: {meta['Sex'].value_counts().to_dict()}")


# ==============================================================================
# STEP 3: Sample matching via Title (not GSM)
# ==============================================================================
banner("STEP 3: Match samples (via Title column)")

overlap_gsm   = len(set(expr.columns) & set(meta["GSM"]))
overlap_title = len(set(expr.columns) & set(meta["Title"]))
print(f"  Overlap with GSM:   {overlap_gsm}")
print(f"  Overlap with Title: {overlap_title}")
meta = meta[meta["Title"].isin(expr.columns)].copy()
print(f"  Matched samples: {len(meta)}")


# ==============================================================================
# STEP 4: ASCL1 extraction (pseudocount 0.1)
# ==============================================================================
banner("STEP 4: ASCL1 extraction (log2 with pseudocount 0.1)")

if "ASCL1" not in expr.index:
    raise KeyError("ASCL1 not found in GCT expression matrix")

meta["ASCL1_raw"]  = expr.loc["ASCL1", meta["Title"].values].values
meta["ASCL1_log2"] = np.log2(meta["ASCL1_raw"].astype(float) + 0.1)
print(f"  Raw range: {meta['ASCL1_raw'].min():.2f} - {meta['ASCL1_raw'].max():.2f}")
print(f"  Non-zero: {(meta['ASCL1_raw'] > 0).sum()} / {len(meta)}")
print(f"  log2(raw + 0.1) mean: {meta['ASCL1_log2'].mean():.4f}")


# ==============================================================================
# STEP 5: Mann-Whitney M vs F + bar plot
# ==============================================================================
banner("STEP 5: Mann-Whitney U test (Male vs Female)")

mv, fv, u_stat, p_val, sig = mann_whitney_sex(
    meta, ascl1_col="ASCL1_log2", sex_col="Sex",
    male_label="male", female_label="female",
)

banner("STEP 6: Generating barplot")

barplot_male_vs_female(
    mv, fv, p_val,
    title="GSE193066 - Japanese MASLD-HCC\nASCL1: Male vs Female",
    y_label="ASCL1 expression\n(log2 raw + 0.1)",
    out_prefix="verify_GSE193066_ASCL1_barplot",
    male_label="male", female_label="female",
)


# ==============================================================================
# SUMMARY
# ==============================================================================
banner("PIPELINE SUMMARY")
n_1st = int((meta.get("biopsy") == "1st biopsy").sum()) if "biopsy" in meta.columns else 0
n_2nd = int((meta.get("biopsy") == "2nd biopsy").sum()) if "biopsy" in meta.columns else 0

print(f"  Platform:       RNA-seq (pre-normalized GCT)")
print(f"  Gene IDs:       Gene symbols")
print(f"  Normalization:  log2(raw + 0.1)")
print(f"  Sex:            ANNOTATED (male / female)")
print(f"  Sample match:   Title column (not GSM)")
print(f"  Samples:        {len(meta)} total")
if n_1st or n_2nd:
    print(f"  Biopsies:       {n_1st} 1st + {n_2nd} 2nd")
print(f"  ASCL1 result:   Mann-Whitney U, p = {p_val:.4e} {sig}")
