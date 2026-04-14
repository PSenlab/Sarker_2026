#!/usr/bin/env python3
# ==============================================================================
# TCGA-LIHC Hepatocellular Carcinoma Validation
# ==============================================================================
#
# Description:
#   Validates Ascl1 sex dimorphism in TCGA-LIHC (RNA-seq, Illumina HiSeq,
#   RSEM normalized). Pan-cancer expression TSV is filtered to LIHC using
#   the TCGA-CDR clinical supplemental table. Sex is annotated in the
#   clinical data ("gender" column).
#
# Dataset-specific quirks:
#   - Pan-cancer expression TSV - filtered to LIHC via TCGA-CDR
#   - Index format is "GENE|ENTREZ" - split on "|" to extract symbol
#   - RSEM normalized - requires log2(x + 1) transform
#   - TCGA barcode matching: -01 = tumor, -11 = normal
#   - Clinical: gender, age, stage, vital status
#   - Reports BOTH all-samples and tumor-only comparisons
#
# Inputs:
#   EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv
#   TCGA-CDR-SupplementalTableS1.xlsx
#
# Output:
#   verify_TCGA_LIHC_ASCL1_barplot.{png,pdf}
#
# ==============================================================================

import os
import numpy as np
import pandas as pd

from validation_utils import (
    banner, mann_whitney_sex, barplot_male_vs_female,
)

# ==============================================================================
# CONFIGURATION - UPDATE THESE PATHS
# ==============================================================================
DATA_DIR  = "path/to/tcga_raw"
EXPR_FILE = os.path.join(
    DATA_DIR, "EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv"
)
CLIN_FILE = os.path.join(DATA_DIR, "TCGA-CDR-SupplementalTableS1.xlsx")


# ==============================================================================
# STEP 1: Pan-cancer expression
# ==============================================================================
banner("STEP 1: Load pan-cancer expression TSV")

expr_raw = pd.read_csv(EXPR_FILE, sep="\t", index_col=0)
print(f"  Pan-cancer expression: {expr_raw.shape[0]} rows x "
      f"{expr_raw.shape[1]} samples")
print(f"  Index format: {expr_raw.index[0]}  (GENE|ENTREZ)")


# ==============================================================================
# STEP 2: Parse GENE|ENTREZ -> symbol, collapse
# ==============================================================================
banner("STEP 2: GENE|ENTREZ parsing + gene-level collapse")

gene_symbol_map = {}
for idx in expr_raw.index:
    parts = str(idx).split("|")
    symbol = parts[0]
    if symbol != "?" and symbol not in gene_symbol_map:
        gene_symbol_map[idx] = symbol

expr_mapped = expr_raw.loc[expr_raw.index.isin(gene_symbol_map.keys())].copy()
expr_mapped["SYMBOL"] = [gene_symbol_map[idx] for idx in expr_mapped.index]
expr = expr_mapped.groupby("SYMBOL").mean()
expr = expr.loc[:, ~expr.columns.duplicated()]

print(f"  Unique genes: {len(set(gene_symbol_map.values()))}")
print(f"  [OK] Gene-level expression: {expr.shape}")


# ==============================================================================
# STEP 3: Clinical metadata - LIHC subset
# ==============================================================================
banner("STEP 3: Clinical metadata (TCGA-CDR -> LIHC)")

clinical = pd.read_excel(CLIN_FILE, sheet_name="TCGA-CDR")
lihc_clin = clinical[clinical["type"] == "LIHC"].copy()
lihc_barcodes = lihc_clin["bcr_patient_barcode"].tolist()
print(f"  LIHC patients in TCGA-CDR: {len(lihc_clin)}")

# Match tumor (-01) and normal (-11) columns
lihc_tumor_cols = [
    c for c in expr.columns
    if any(c.startswith(b) for b in lihc_barcodes) and "-01" in c
]
lihc_normal_cols = [
    c for c in expr.columns
    if any(c.startswith(b) for b in lihc_barcodes) and "-11" in c
]
print(f"  Tumor samples (-01):  {len(lihc_tumor_cols)}")
print(f"  Normal samples (-11): {len(lihc_normal_cols)}")

def _barcode_to_patient(col):
    return "-".join(col.split("-")[:3])

CLIN_COLS = [
    "bcr_patient_barcode", "gender",
    "age_at_initial_pathologic_diagnosis",
    "ajcc_pathologic_tumor_stage", "vital_status",
]

meta_t = pd.DataFrame({"Sample": lihc_tumor_cols})
meta_t["Patient"] = meta_t["Sample"].apply(_barcode_to_patient)
meta_t = meta_t.merge(
    lihc_clin[CLIN_COLS], left_on="Patient",
    right_on="bcr_patient_barcode", how="left",
)
meta_t["Tissue"] = "Tumor"

meta_n = pd.DataFrame({"Sample": lihc_normal_cols})
meta_n["Patient"] = meta_n["Sample"].apply(_barcode_to_patient)
meta_n = meta_n.merge(
    lihc_clin[CLIN_COLS], left_on="Patient",
    right_on="bcr_patient_barcode", how="left",
)
meta_n["Tissue"] = "Normal"

meta = pd.concat([meta_t, meta_n], ignore_index=True)
meta["Sex"] = meta["gender"].str.upper()
meta["Age"] = pd.to_numeric(
    meta["age_at_initial_pathologic_diagnosis"], errors="coerce"
)

def _simplify_stage(s):
    if pd.isna(s) or str(s).startswith("["):
        return np.nan
    s = str(s).replace("Stage ", "")
    if s.startswith("IV"):  return "IV"
    if s.startswith("III"): return "III"
    if s.startswith("II"):  return "II"
    if s.startswith("I"):   return "I"
    return np.nan

meta["Stage"] = meta["ajcc_pathologic_tumor_stage"].apply(_simplify_stage)

n_tumor  = int((meta["Tissue"] == "Tumor").sum())
n_normal = int((meta["Tissue"] == "Normal").sum())
n_m = int((meta["Sex"] == "MALE").sum())
n_f = int((meta["Sex"] == "FEMALE").sum())
print(f"  [OK] Samples: {len(meta)} ({n_tumor} tumor, {n_normal} normal)")
print(f"  Sex annotated: MALE={n_m}, FEMALE={n_f}")


# ==============================================================================
# STEP 4: ASCL1 extraction (log2(RSEM + 1))
# ==============================================================================
banner("STEP 4: ASCL1 extraction (log2(RSEM + 1))")

if "ASCL1" not in expr.index:
    raise KeyError("ASCL1 not found in gene-level expression")

ascl1_series = expr.loc["ASCL1"]
meta["ASCL1_expr"] = meta["Sample"].map(ascl1_series).astype(float)
meta["ASCL1_log2"] = np.log2(meta["ASCL1_expr"] + 1)

print(f"  RSEM range:       {meta['ASCL1_expr'].min():.4f} - "
      f"{meta['ASCL1_expr'].max():.4f}")
print(f"  log2(RSEM+1):     mean={meta['ASCL1_log2'].mean():.4f} +/- "
      f"{meta['ASCL1_log2'].std():.4f}")
print(f"  Non-zero:         {(meta['ASCL1_expr'] > 0).sum()} / {len(meta)}")


# ==============================================================================
# STEP 5: Mann-Whitney tests (all samples + tumor-only)
# ==============================================================================
banner("STEP 5: Mann-Whitney U test (all samples)")

mv, fv, u_stat, p_val, sig = mann_whitney_sex(
    meta, ascl1_col="ASCL1_log2", sex_col="Sex",
    male_label="MALE", female_label="FEMALE",
)

banner("STEP 5b: Mann-Whitney U test (tumor only)")

mv_t, fv_t, u_t, p_t, sig_t = mann_whitney_sex(
    meta[meta["Tissue"] == "Tumor"],
    ascl1_col="ASCL1_log2", sex_col="Sex",
    male_label="MALE", female_label="FEMALE",
)


# ==============================================================================
# STEP 6: Bar plot (all samples)
# ==============================================================================
banner("STEP 6: Generating barplot")

barplot_male_vs_female(
    mv, fv, p_val,
    title="TCGA-LIHC - Hepatocellular Carcinoma\nASCL1: Male vs Female (all samples)",
    y_label="ASCL1 expression\n(log2(RSEM + 1))",
    out_prefix="verify_TCGA_LIHC_ASCL1_barplot",
    male_label="Male", female_label="Female",
)


# ==============================================================================
# SUMMARY
# ==============================================================================
banner("PIPELINE SUMMARY")
print(f"  Platform:       RNA-seq (Illumina HiSeq, RSEM normalized)")
print(f"  Gene IDs:       GENE|ENTREZ -> symbol (first field)")
print(f"  Gene collapse:  mean")
print(f"  Normalization:  log2(RSEM + 1)")
print(f"  Sex:            ANNOTATED (TCGA-CDR 'gender')")
print(f"  Sample match:   TCGA barcode (-01 tumor / -11 normal)")
print(f"  Samples:        {len(meta)} total "
      f"({n_tumor} tumor, {n_normal} normal; {n_m}M / {n_f}F)")
print(f"  ASCL1 all:      Mann-Whitney U, p = {p_val:.4e} {sig}")
print(f"  ASCL1 tumor:    Mann-Whitney U, p = {p_t:.4e} {sig_t}")
