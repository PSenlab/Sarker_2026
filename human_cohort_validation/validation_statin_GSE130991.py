#!/usr/bin/env python3
# ==============================================================================
# Statin Treatment MASLD Validation - GSE130991
# ==============================================================================
#
# Description:
#   Validates Ascl1 sex dimorphism in the statin treatment MASLD cohort
#   (Affymetrix HuGene 2.0 ST, GPL20265). Data are RMA-normalized
#   (log2-scale). Sex is annotated in the metadata ("gender" column).
#
# Dataset-specific quirks:
#   - Affymetrix HuGene 2.0 ST (numeric probe IDs)
#   - Probe -> gene via GPL20265 SOFT file "gene_assignment" field
#     (format: "NM_xxx // GENE_SYMBOL // description // cytoband // entrez_id")
#   - RMA-normalized; no additional log2 transform
#   - Multi-probe genes averaged
#   - Sex annotated as "Male" / "Female" in the "gender" field
#   - Clinical variables: statin treatment, BMI, age
#
# Inputs:
#   GSE130991_raw/GSE130991_series_matrix.txt.gz
#   GSE130991_raw/GPL20265_family.soft.gz
#
# Output:
#   verify_GSE130991_ASCL1_barplot.{png,pdf}
#
# Reference:
#   Sarker et al. (2026) Cell Metabolism
# ==============================================================================

import os
import gzip
import numpy as np
import pandas as pd

from validation_utils import (
    banner, parse_series_matrix, mann_whitney_sex, barplot_male_vs_female,
)

# ==============================================================================
# CONFIGURATION - UPDATE THESE PATHS
# ==============================================================================
DATA_DIR    = "GSE130991_raw"
SERIES_FILE = os.path.join(DATA_DIR, "GSE130991_series_matrix.txt.gz")
GPL_FILE    = os.path.join(DATA_DIR, "GPL20265_family.soft.gz")


# ==============================================================================
# STEP 1: Load probe-level expression + metadata
# ==============================================================================
banner("STEP 1: Load expression + metadata (series matrix)")

expr_raw, meta, platform_id = parse_series_matrix(SERIES_FILE)
print(f"  Platform: {platform_id}")
print(f"  Expression: {expr_raw.shape[0]} probes x {expr_raw.shape[1]} samples")
print(f"  Samples: {len(meta)}")


# ==============================================================================
# STEP 2: Parse GPL20265 SOFT file for probe -> gene
# ==============================================================================
banner("STEP 2: Probe-to-gene mapping (GPL20265 SOFT)")

gpl_header = None
gpl_rows = []
in_table = False

with gzip.open(GPL_FILE, "rt", errors="replace") as f:
    for line in f:
        if line.startswith("!platform_table_begin"):
            in_table = True
            continue
        if line.startswith("!platform_table_end"):
            break
        if in_table:
            if gpl_header is None:
                gpl_header = line.strip().split("\t")
            else:
                gpl_rows.append(line.strip().split("\t"))

gpl = pd.DataFrame(
    gpl_rows,
    columns=gpl_header[:len(gpl_rows[0])] if gpl_rows else gpl_header,
)
print(f"  GPL annotation rows: {len(gpl)}")

def _extract_symbol(ga):
    """
    gene_assignment format: "NM_xxx // GENE_SYMBOL // description // cytoband // entrez_id"
    We want the second '//' delimited field.
    """
    if pd.isna(ga) or ga.startswith("---"):
        return None
    parts = [p.strip() for p in ga.split("//")]
    if len(parts) >= 2:
        symbol = parts[1].strip()
        if symbol and symbol != "---":
            return symbol
    return None

gpl["SYMBOL"] = gpl["gene_assignment"].apply(_extract_symbol)
gpl_mapped = gpl[gpl["SYMBOL"].notna()][["ID", "SYMBOL"]].copy()
gpl_mapped["ID"] = gpl_mapped["ID"].astype(str)
print(f"  Mapped probes: {len(gpl_mapped)} -> "
      f"{gpl_mapped['SYMBOL'].nunique()} unique genes")

# Map and collapse (mean) to gene level
expr_raw.index = expr_raw.index.astype(str)
common_probes = expr_raw.index.intersection(gpl_mapped["ID"])
mapped = expr_raw.loc[common_probes].copy()
probe_dedup = gpl_mapped.drop_duplicates(subset="ID").set_index("ID")
mapped["SYMBOL"] = probe_dedup.loc[mapped.index, "SYMBOL"].values
expr = mapped.groupby("SYMBOL").mean()
print(f"  [OK] Gene-level expression: {expr.shape[0]} genes "
      f"(multi-probe collapsed by mean)")


# ==============================================================================
# STEP 3: Metadata cleanup
# ==============================================================================
banner("STEP 3: Metadata cleanup")

meta["Sex"] = meta["gender"].str.strip()
if "age" in meta.columns:
    meta["Age"] = pd.to_numeric(meta["age"], errors="coerce")
if "bmi" in meta.columns:
    meta["BMI"] = pd.to_numeric(meta["bmi"], errors="coerce")
if "statin treatment" in meta.columns:
    meta["Statin"] = meta["statin treatment"].str.strip()

print(f"  [OK] Sex annotated: {meta['Sex'].value_counts().to_dict()}")
if "Statin" in meta.columns:
    print(f"  Statin: {meta['Statin'].value_counts().to_dict()}")


# ==============================================================================
# STEP 4: ASCL1 extraction (RMA, no transform)
# ==============================================================================
banner("STEP 4: ASCL1 extraction (RMA values direct)")

if "ASCL1" not in expr.index:
    raise KeyError("ASCL1 not found in gene-level expression")

meta["ASCL1_expr"] = expr.loc["ASCL1", meta["GSM"].values].values
print(f"  RMA range: {meta['ASCL1_expr'].min():.4f} - "
      f"{meta['ASCL1_expr'].max():.4f}")
print(f"  Mean:      {meta['ASCL1_expr'].mean():.4f} +/- "
      f"{meta['ASCL1_expr'].std():.4f}")


# ==============================================================================
# STEP 5: Mann-Whitney M vs F + bar plot
# ==============================================================================
banner("STEP 5: Mann-Whitney U test (Male vs Female)")

mv, fv, u_stat, p_val, sig = mann_whitney_sex(
    meta, ascl1_col="ASCL1_expr", sex_col="Sex",
    male_label="Male", female_label="Female",
)

banner("STEP 6: Generating barplot")

barplot_male_vs_female(
    mv, fv, p_val,
    title="GSE130991 - Statin Treatment MASLD\nASCL1: Male vs Female",
    y_label="ASCL1 expression\n(RMA)",
    out_prefix="verify_GSE130991_ASCL1_barplot",
    male_label="Male", female_label="Female",
)


# ==============================================================================
# SUMMARY
# ==============================================================================
banner("PIPELINE SUMMARY")
print(f"  Platform:       Affymetrix HuGene 2.0 ST ({platform_id})")
print(f"  Probe mapping:  GPL20265 SOFT 'gene_assignment' field")
print(f"  Gene collapse:  mean across probes")
print(f"  Normalization:  RMA (no additional transform)")
print(f"  Sex:            ANNOTATED ('gender' column)")
print(f"  Samples:        {len(meta)} total "
      f"({(meta['Sex']=='Male').sum()}M / {(meta['Sex']=='Female').sum()}F)")
print(f"  ASCL1 result:   Mann-Whitney U, p = {p_val:.4e} {sig}")
