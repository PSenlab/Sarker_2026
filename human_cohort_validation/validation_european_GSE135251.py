#!/usr/bin/env python3
# ==============================================================================
# European MASLD Validation - GSE135251
# ==============================================================================
#
# Description:
#   Validates Ascl1 sex dimorphism in the European MASLD RNA-seq cohort.
#   Raw counts are distributed as per-sample TSV files with versioned
#   Ensembl IDs. Sex is not annotated, so it is inferred via KMeans on
#   XIST plus the mean of 7 Y-chromosome genes (more robust than a single
#   Y-gene because of individual variability in Y-gene expression).
#
# Dataset-specific quirks:
#   - Per-sample count files (one TSV per GSM)
#   - Ensembl IDs with version numbers (e.g. ENSG00000139352.5) -> strip
#   - Gene matching by Ensembl ID, not symbol
#   - CPM + log2(x+1) normalization
#   - Sex INFERRED via KMeans (XIST + mean of 7 Y-chromosome genes)
#
# Inputs:
#   GSE135251_raw/GSE135251_series_matrix.txt.gz
#   GSE135251_raw/*counts*txt (per-sample)
#
# Output:
#   verify_GSE135251_ASCL1_barplot.{png,pdf}
#
# ==============================================================================

import os
import glob
import numpy as np
import pandas as pd

from validation_utils import (
    banner, parse_series_matrix, infer_sex_kmeans,
    mann_whitney_sex, barplot_male_vs_female,
)

# ==============================================================================
# CONFIGURATION - UPDATE THESE PATHS
# ==============================================================================
DATA_DIR    = "GSE135251_raw"
SERIES_FILE = os.path.join(DATA_DIR, "GSE135251_series_matrix.txt.gz")
COUNTS_GLOB = os.path.join(DATA_DIR, "*counts*txt")

ASCL1_ENSG = "ENSG00000139352"
XIST_ENSG  = "ENSG00000229807"
Y_GENES = {
    "ENSG00000129824": "RPS4Y1",
    "ENSG00000067048": "DDX3Y",
    "ENSG00000198692": "EIF1AY",
    "ENSG00000012817": "KDM5D",
    "ENSG00000183878": "UTY",
    "ENSG00000067646": "ZFY",
    "ENSG00000114374": "USP9Y",
}


# ==============================================================================
# STEP 1: Metadata
# ==============================================================================
banner("STEP 1: Metadata extraction")

_, meta, _ = parse_series_matrix(SERIES_FILE, include_expr=False)
print(f"  Samples: {len(meta)}")


# ==============================================================================
# STEP 2: Load per-sample count files
# ==============================================================================
banner("STEP 2: Load per-sample count files")

count_files = sorted([f for f in glob.glob(COUNTS_GLOB) if not f.endswith(".gz")])
print(f"  Count files: {len(count_files)}")

expr_list = []
for f in count_files:
    gsm = os.path.basename(f).split("_")[0]
    df = pd.read_csv(f, sep="\t", header=None, index_col=0)
    df.columns = [gsm]
    expr_list.append(df)

expr = pd.concat(expr_list, axis=1)
print(f"  Expression: {expr.shape}")
print(f"  [OK] Example Ensembl ID (versioned): {expr.index[0]}")

# Strip Ensembl version numbers
expr.index = [g.split(".")[0] for g in expr.index]
print(f"  [OK] After version strip: {expr.index[0]}")


# ==============================================================================
# STEP 3: Match samples
# ==============================================================================
banner("STEP 3: Match samples to metadata")

common = sorted(set(expr.columns) & set(meta["GSM"]))
meta = meta[meta["GSM"].isin(common)].copy()
print(f"  Common samples: {len(common)} / {len(meta) + len(set(meta['GSM']) - set(expr.columns))}")


# ==============================================================================
# STEP 4: CPM + log2 normalization
# ==============================================================================
banner("STEP 4: CPM + log2(x+1) normalization")

expr_sub = expr[meta["GSM"]]
lib_sizes = expr_sub.sum(axis=0)
cpm = expr_sub.div(lib_sizes, axis=1) * 1e6
log2cpm = np.log2(cpm + 1)

print(f"  Library size: {lib_sizes.min():.0f} - {lib_sizes.max():.0f} "
      f"(mean {lib_sizes.mean():.0f})")


# ==============================================================================
# STEP 5: Sex inference
# ==============================================================================
banner("STEP 5: Sex inference (KMeans: XIST + 7 Y-genes)")

meta["Sex"] = infer_sex_kmeans(
    log2_expr   = log2cpm,
    sample_ids  = meta["GSM"].tolist(),
    xist_key    = XIST_ENSG,
    y_gene_keys = list(Y_GENES.keys()),
)
print(f"  Counts: {meta['Sex'].value_counts().to_dict()}")


# ==============================================================================
# STEP 6: ASCL1 extraction
# ==============================================================================
banner("STEP 6: ASCL1 extraction")

if ASCL1_ENSG not in expr.index:
    raise KeyError(f"ASCL1 ({ASCL1_ENSG}) not found in expression matrix")

meta["ASCL1_raw"]     = expr.loc[ASCL1_ENSG, meta["GSM"]].values
meta["ASCL1_log2cpm"] = log2cpm.loc[ASCL1_ENSG, meta["GSM"]].values
print(f"  Raw count range: {meta['ASCL1_raw'].min():.0f} - {meta['ASCL1_raw'].max():.0f}")
print(f"  Non-zero: {(meta['ASCL1_raw'] > 0).sum()} / {len(meta)}")
print(f"  log2(CPM+1) mean: {meta['ASCL1_log2cpm'].mean():.4f}")


# ==============================================================================
# STEP 7: Mann-Whitney M vs F + bar plot
# ==============================================================================
banner("STEP 7: Mann-Whitney U test (Male vs Female)")

mv, fv, u_stat, p_val, sig = mann_whitney_sex(
    meta, ascl1_col="ASCL1_log2cpm", sex_col="Sex",
    male_label="Male", female_label="Female",
)

banner("STEP 8: Generating barplot")

barplot_male_vs_female(
    mv, fv, p_val,
    title="GSE135251 - European MASLD\nASCL1: Male vs Female",
    y_label="ASCL1 expression\n(log2 CPM+1)",
    out_prefix="verify_GSE135251_ASCL1_barplot",
    male_label="Male", female_label="Female",
)


# ==============================================================================
# SUMMARY
# ==============================================================================
banner("PIPELINE SUMMARY")
print(f"  Platform:       RNA-seq (per-sample count files)")
print(f"  Gene IDs:       Ensembl IDs (version-stripped)")
print(f"  Normalization:  log2(CPM + 1)")
print(f"  Sex:            INFERRED (KMeans on XIST + 7 Y-genes)")
print(f"  Samples:        {len(meta)} total "
      f"({(meta['Sex']=='Male').sum()}M / {(meta['Sex']=='Female').sum()}F)")
print(f"  ASCL1 result:   Mann-Whitney U, p = {p_val:.4e} {sig}")
