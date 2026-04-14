#!/usr/bin/env python3
# ==============================================================================
# Japanese MASLD Validation - GSE167523
# ==============================================================================
#
# Description:
#   Validates Ascl1 sex dimorphism in the Japanese MASLD RNA-seq cohort.
#   Raw gene counts are distributed as a single matrix with gene symbols
#   as the row index. Sex (M/F) is annotated directly in the series matrix,
#   so no KMeans inference is needed.
#
# Dataset-specific quirks:
#   - Single raw counts matrix (gene symbols in index)
#   - Sex directly annotated as "M" / "F"
#   - Sample matching via Title extraction (NAFLD### pattern) with
#     positional fallback
#   - CPM + log2(x+1) normalization
#
# Inputs:
#   GSE167523_raw/GSE167523_series_matrix.txt.gz
#   GSE167523_raw/GSE167523_Raw_gene_counts_matrix.txt.gz
#
# Output:
#   verify_GSE167523_ASCL1_barplot.{png,pdf}
#
# ==============================================================================

import os
import gzip
import numpy as np
import pandas as pd

from validation_utils import (
    banner, mann_whitney_sex, barplot_male_vs_female,
)

# ==============================================================================
# CONFIGURATION - UPDATE THESE PATHS
# ==============================================================================
DATA_DIR    = "GSE167523_raw"
SERIES_FILE = os.path.join(DATA_DIR, "GSE167523_series_matrix.txt.gz")
COUNTS_FILE = os.path.join(DATA_DIR, "GSE167523_Raw_gene_counts_matrix.txt.gz")


# ==============================================================================
# STEP 1: Metadata (sex is directly annotated)
# ==============================================================================
banner("STEP 1: Metadata extraction (sex annotated)")

metadata = {}
with gzip.open(SERIES_FILE, "rt") as f:
    for line in f:
        if line.startswith("!Sample_geo_accession"):
            metadata["GSM"] = [x.strip('"') for x in line.strip().split("\t")[1:]]
        if line.startswith("!Sample_title"):
            metadata["Title"] = [x.strip('"') for x in line.strip().split("\t")[1:]]
        if "disease subtype:" in line:
            metadata["Subtype"] = [
                x.split(": ")[1].strip('"') for x in line.strip().split("\t")[1:]
            ]
        if "age:" in line:
            metadata["Age"] = [
                int(x.split(": ")[1].strip('"')) for x in line.strip().split("\t")[1:]
            ]
        if "gender:" in line:
            metadata["Sex"] = [
                x.split(": ")[1].strip('"') for x in line.strip().split("\t")[1:]
            ]

meta = pd.DataFrame(metadata)
print(f"  Samples: {len(meta)}")
print(f"  Sex counts: {meta['Sex'].value_counts().to_dict()} (M/F pre-annotated)")
print(f"  Subtypes: {meta['Subtype'].value_counts().to_dict()}")
print(f"  Age: {meta['Age'].min()}-{meta['Age'].max()}")


# ==============================================================================
# STEP 2: Raw counts matrix
# ==============================================================================
banner("STEP 2: Load raw counts matrix")

counts = pd.read_csv(COUNTS_FILE, sep="\t", index_col=0)
print(f"  Shape: {counts.shape}")
print(f"  Gene index (first 5): {list(counts.index[:5])}")
if "ASCL1" not in counts.index:
    raise KeyError("ASCL1 not in counts index")
print(f"  ASCL1 raw count range: {counts.loc['ASCL1'].min():.0f} - "
      f"{counts.loc['ASCL1'].max():.0f}")


# ==============================================================================
# STEP 3: Sample matching
# ==============================================================================
banner("STEP 3: Sample matching")

meta["sample_id"] = meta["Title"].str.extract(r"(NAFLD\d+)", expand=False)
if meta["sample_id"].isna().any():
    if counts.shape[1] == len(meta):
        print(f"  Falling back to positional matching ({len(meta)} samples)")
        meta["sample_id"] = counts.columns.tolist()
    else:
        raise ValueError(f"counts cols ({counts.shape[1]}) != meta rows ({len(meta)})")

missing = [s for s in meta["sample_id"] if s not in counts.columns]
if missing:
    raise ValueError(f"{len(missing)} sample IDs not in counts columns")
print(f"  [OK] All {len(meta)} sample_ids found in counts columns")


# ==============================================================================
# STEP 4: CPM + log2 normalization
# ==============================================================================
banner("STEP 4: CPM + log2(x+1) normalization")

lib_sizes = counts.sum(axis=0)
cpm = counts.div(lib_sizes, axis=1) * 1e6
log2cpm = np.log2(cpm + 1)
print(f"  Library size: {lib_sizes.min():.0f} - {lib_sizes.max():.0f} "
      f"(mean {lib_sizes.mean():.0f})")


# ==============================================================================
# STEP 5: ASCL1 extraction
# ==============================================================================
banner("STEP 5: ASCL1 extraction")

meta["ASCL1_log2cpm"] = log2cpm.loc["ASCL1", meta["sample_id"]].values
meta["ASCL1_raw"]     = counts.loc["ASCL1", meta["sample_id"]].values
print(f"  log2(CPM+1): mean={meta['ASCL1_log2cpm'].mean():.4f}, "
      f"std={meta['ASCL1_log2cpm'].std():.4f}")
print(f"  Zero counts: {int((meta['ASCL1_raw'] == 0).sum())} / {len(meta)}")


# ==============================================================================
# STEP 6: Mann-Whitney M vs F + bar plot
# ==============================================================================
banner("STEP 6: Mann-Whitney U test (Male vs Female)")

mv, fv, u_stat, p_val, sig = mann_whitney_sex(
    meta, ascl1_col="ASCL1_log2cpm", sex_col="Sex",
    male_label="M", female_label="F",
)

banner("STEP 7: Generating barplot")

barplot_male_vs_female(
    mv, fv, p_val,
    title="GSE167523 - Japanese MASLD\nASCL1: Male vs Female",
    y_label="ASCL1 expression\n(log2 CPM+1)",
    out_prefix="verify_GSE167523_ASCL1_barplot",
    male_label="M", female_label="F",
)


# ==============================================================================
# SUMMARY
# ==============================================================================
banner("PIPELINE SUMMARY")
print(f"  Platform:       RNA-seq (single counts matrix)")
print(f"  Gene IDs:       Gene symbols")
print(f"  Normalization:  log2(CPM + 1)")
print(f"  Sex:            ANNOTATED in metadata (no inference)")
print(f"  Samples:        {len(meta)} total "
      f"({(meta['Sex']=='M').sum()}M / {(meta['Sex']=='F').sum()}F)")
print(f"  Subtypes:       "
      f"{int((meta['Subtype']=='NAFL').sum())} NAFL + "
      f"{int((meta['Subtype']=='NASH').sum())} NASH")
print(f"  ASCL1 result:   Mann-Whitney U, p = {p_val:.4e} {sig}")
