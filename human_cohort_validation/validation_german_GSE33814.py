#!/usr/bin/env python3
# ==============================================================================
# German MASLD Validation - GSE33814
# ==============================================================================
#
# Description:
#   Validates Ascl1 sex dimorphism in the German MASLD microarray cohort
#   (Illumina HumanWG-6 v3, GPL6884). Data are already quantile-normalized
#   and log2-scale. Sex is not annotated, so it is inferred via KMeans on
#   XIST + RPS4Y1 expression.
#
# Dataset-specific quirks:
#   - Illumina BeadArray (ILMN_ probe IDs)
#   - Probe -> gene symbol via GPL6884 BGX annotation file
#   - Multi-probe genes collapsed by MAX (not mean)
#   - Sex inferred (XIST + RPS4Y1, single Y-gene)
#   - ASCL1 extracted by probe ID (ILMN_1701653) - no log transform
#
# Inputs:
#   GSE33814_raw/GSE33814_series_matrix.txt.gz
#   GSE33814_raw/GPL6884_HumanWG-6_V3_0_R0_11282955_A.bgx.gz
#
# Output:
#   verify_GSE33814_ASCL1_barplot.{png,pdf}
#
# Reference:
#   Sarker et al. (2026) Cell Metabolism
# ==============================================================================

import os
import gzip
from io import StringIO
import pandas as pd

from validation_utils import (
    banner, parse_series_matrix, infer_sex_kmeans,
    mann_whitney_sex, barplot_male_vs_female,
)

# ==============================================================================
# CONFIGURATION - UPDATE THESE PATHS
# ==============================================================================
DATA_DIR      = "GSE33814_raw"
SERIES_FILE   = os.path.join(DATA_DIR, "GSE33814_series_matrix.txt.gz")
BGX_FILE      = os.path.join(DATA_DIR, "GPL6884_HumanWG-6_V3_0_R0_11282955_A.bgx.gz")

ASCL1_PROBE   = "ILMN_1701653"


# ==============================================================================
# STEP 1: Load probe-level expression + metadata
# ==============================================================================
banner("STEP 1: Load expression + metadata (series matrix)")

expr_raw, meta, platform_id = parse_series_matrix(SERIES_FILE)
print(f"  Platform: {platform_id}")
print(f"  Expression: {expr_raw.shape[0]} probes x {expr_raw.shape[1]} samples")
print(f"  Samples: {len(meta)}")


# ==============================================================================
# STEP 2: Probe-to-gene mapping via BGX annotation
# ==============================================================================
banner("STEP 2: Probe-to-gene mapping (GPL6884 BGX)")

probe_lines = []
reading_probes = False
with gzip.open(BGX_FILE, "rt", errors="replace") as f:
    for line in f:
        if line.strip() == "[Probes]":
            reading_probes = True
            continue
        if reading_probes and line.startswith("["):
            break
        if reading_probes:
            probe_lines.append(line.strip())

annot = pd.read_csv(StringIO("\n".join(probe_lines)), sep="\t", low_memory=False)
probe_to_gene = annot.set_index("Probe_Id")["Symbol"].dropna().to_dict()

print(f"  Probes with gene symbol: {len(probe_to_gene)}")
print(f"  Unique genes: {len(set(probe_to_gene.values()))}")

# Map and collapse multi-probe genes by MAX
expr_genes = expr_raw.copy()
expr_genes.index = expr_genes.index.map(lambda x: probe_to_gene.get(x, x))
expr_genes = expr_genes[~expr_genes.index.str.startswith("ILMN_")]
expr_genes = expr_genes.groupby(expr_genes.index).max()
print(f"  [OK] Gene-level expression: {expr_genes.shape[0]} genes "
      f"(multi-probe collapsed by MAX)")


# ==============================================================================
# STEP 3: Sex inference
# ==============================================================================
banner("STEP 3: Sex inference (KMeans: XIST + RPS4Y1)")

meta["Sex"] = infer_sex_kmeans(
    log2_expr  = expr_genes,
    sample_ids = meta["GSM"].tolist(),
    xist_key   = "XIST",
    y_gene_keys = ["RPS4Y1"],
)
print(f"  Counts: {meta['Sex'].value_counts().to_dict()}")


# ==============================================================================
# STEP 4: ASCL1 extraction (probe-level)
# ==============================================================================
banner("STEP 4: ASCL1 extraction (pre-normalized log2)")

meta["ASCL1_log2"] = expr_raw.loc[ASCL1_PROBE, meta["GSM"]].values
print(f"  ASCL1 probe: {ASCL1_PROBE}")
print(f"  Range: {meta['ASCL1_log2'].min():.4f} - {meta['ASCL1_log2'].max():.4f}")
print(f"  Mean:  {meta['ASCL1_log2'].mean():.4f} +/- {meta['ASCL1_log2'].std():.4f}")


# ==============================================================================
# STEP 5: Mann-Whitney M vs F + bar plot
# ==============================================================================
banner("STEP 5: Mann-Whitney U test (Male vs Female)")

mv, fv, u_stat, p_val, sig = mann_whitney_sex(
    meta, ascl1_col="ASCL1_log2", sex_col="Sex",
    male_label="Male", female_label="Female",
)

banner("STEP 6: Generating barplot")

barplot_male_vs_female(
    mv, fv, p_val,
    title="GSE33814 - German MASLD\nASCL1: Male vs Female",
    y_label="ASCL1 expression\n(log2)",
    out_prefix="verify_GSE33814_ASCL1_barplot",
    male_label="Male", female_label="Female",
)


# ==============================================================================
# SUMMARY
# ==============================================================================
banner("PIPELINE SUMMARY")
print(f"  Platform:       Illumina HumanWG-6 v3 ({platform_id})")
print(f"  Probe mapping:  GPL6884 BGX annotation")
print(f"  Gene collapse:  MAX across probes")
print(f"  Normalization:  Quantile-normalized log2 (no additional transform)")
print(f"  Sex:            INFERRED (KMeans on XIST + RPS4Y1)")
print(f"  Samples:        {len(meta)} total "
      f"({(meta['Sex']=='Male').sum()}M / {(meta['Sex']=='Female').sum()}F)")
print(f"  ASCL1 result:   Mann-Whitney U, p = {p_val:.4e} {sig}")
