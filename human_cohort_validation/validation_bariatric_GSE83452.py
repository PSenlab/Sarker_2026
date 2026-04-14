#!/usr/bin/env python3
# ==============================================================================
# Bariatric Surgery MASLD Validation - GSE83452
# ==============================================================================
#
# Description:
#   Validates Ascl1 sex dimorphism in the bariatric surgery MASLD cohort
#   (Affymetrix HuGene 2.0 ST, GPL16686). Data are RMA-normalized. Sex is
#   annotated in the metadata ("gender" column). Longitudinal design
#   (baseline + follow-up).
#
# Dataset-specific quirks:
#   - Affymetrix HuGene 2.0 ST (numeric probe IDs)
#   - Probe -> gene via pre-computed GPL16686_probe_to_gene.csv
#     (PROBEID, SYMBOL columns)
#   - RMA-normalized; no additional log2 transform
#   - Multi-probe genes averaged
#   - Sex annotated as "male" / "female" (lowercase)
#   - Longitudinal: baseline + follow-up biopsies
#   - Clinical variables: NASH status, intervention type, time
#
# Inputs:
#   GSE83452_raw/GSE83452_series_matrix.txt.gz
#   GSE83452_raw/GPL16686_probe_to_gene.csv
#
# Output:
#   verify_GSE83452_ASCL1_barplot.{png,pdf}
#
# Reference:
#   Sarker et al. (2026) Cell Metabolism
# ==============================================================================

import os
import pandas as pd

from validation_utils import (
    banner, parse_series_matrix, mann_whitney_sex, barplot_male_vs_female,
)

# ==============================================================================
# CONFIGURATION - UPDATE THESE PATHS
# ==============================================================================
DATA_DIR      = "GSE83452_raw"
SERIES_FILE   = os.path.join(DATA_DIR, "GSE83452_series_matrix.txt.gz")
PROBE_MAP_CSV = os.path.join(DATA_DIR, "GPL16686_probe_to_gene.csv")


# ==============================================================================
# STEP 1: Load expression + metadata
# ==============================================================================
banner("STEP 1: Load expression + metadata (series matrix)")

expr_raw, meta, platform_id = parse_series_matrix(SERIES_FILE)
print(f"  Platform: {platform_id}")
print(f"  Expression: {expr_raw.shape[0]} probes x {expr_raw.shape[1]} samples")
print(f"  Samples: {len(meta)}")


# ==============================================================================
# STEP 2: Probe-to-gene mapping (pre-computed CSV)
# ==============================================================================
banner("STEP 2: Probe-to-gene mapping (GPL16686 CSV)")

probe_map = pd.read_csv(PROBE_MAP_CSV).dropna(subset=["SYMBOL"])
probe_map["PROBEID"] = probe_map["PROBEID"].astype(str)
expr_raw.index = expr_raw.index.astype(str)

print(f"  Mapped probes: {len(probe_map)}")
print(f"  Unique genes:  {probe_map['SYMBOL'].nunique()}")

common_probes = expr_raw.index.intersection(probe_map["PROBEID"])
mapped = expr_raw.loc[common_probes].copy()
probe_dedup = probe_map.drop_duplicates(subset="PROBEID").set_index("PROBEID")
mapped["SYMBOL"] = probe_dedup.loc[mapped.index, "SYMBOL"].values
expr = mapped.groupby("SYMBOL").mean()
print(f"  [OK] Gene-level expression: {expr.shape[0]} genes "
      f"(multi-probe collapsed by mean)")


# ==============================================================================
# STEP 3: Metadata cleanup
# ==============================================================================
banner("STEP 3: Metadata cleanup")

meta["Sex"] = meta["gender"].str.strip().str.lower()
if "age" in meta.columns:
    meta["Age"] = pd.to_numeric(meta["age"], errors="coerce")
if "liver status" in meta.columns:
    meta["NASH_status"] = meta["liver status"].str.strip()
if "type of intervention" in meta.columns:
    meta["Intervention"] = meta["type of intervention"].str.strip()
if "time" in meta.columns:
    meta["Time"] = meta["time"].str.strip()

print(f"  [OK] Sex annotated: {meta['Sex'].value_counts().to_dict()}")
if "NASH_status" in meta.columns:
    print(f"  NASH status: {meta['NASH_status'].value_counts().to_dict()}")
if "Time" in meta.columns:
    print(f"  Time: {meta['Time'].value_counts().to_dict()}")


# ==============================================================================
# STEP 4: ASCL1 extraction (RMA direct)
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
    male_label="male", female_label="female",
)

banner("STEP 6: Generating barplot")

barplot_male_vs_female(
    mv, fv, p_val,
    title="GSE83452 - Bariatric Surgery MASLD\nASCL1: Male vs Female",
    y_label="ASCL1 expression\n(RMA)",
    out_prefix="verify_GSE83452_ASCL1_barplot",
    male_label="male", female_label="female",
)


# ==============================================================================
# SUMMARY
# ==============================================================================
banner("PIPELINE SUMMARY")
print(f"  Platform:       Affymetrix HuGene 2.0 ST ({platform_id})")
print(f"  Probe mapping:  GPL16686_probe_to_gene.csv (pre-computed)")
print(f"  Gene collapse:  mean across probes")
print(f"  Normalization:  RMA (no additional transform)")
print(f"  Sex:            ANNOTATED ('gender' column)")
print(f"  Samples:        {len(meta)} total "
      f"({(meta['Sex']=='male').sum()}M / {(meta['Sex']=='female').sum()}F)")
print(f"  ASCL1 result:   Mann-Whitney U, p = {p_val:.4e} {sig}")
