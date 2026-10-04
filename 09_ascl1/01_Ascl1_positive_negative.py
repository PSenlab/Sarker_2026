#!/usr/bin/env python3
# ==============================================================================
# Ascl1 Expression Analysis in Liver Aging (Female Hepatocytes)
# ==============================================================================
#
# Description:
#   Profiles Ascl1 expression across the full single-nucleus multi-ome liver
#   aging dataset, then focuses on female hepatocytes where Ascl1 shows a
#   sex-specific aging signal. Ascl1+ vs Ascl1- differential expression and
#   GSEA are done per lobular zone in 04_Ascl1_zonation.py.
#
# Pipeline:
#   1. Load data
#   2. Global Ascl1 profiling (sex / age / cell type)
#   3. Female-only Ascl1 analysis
#   4. Female hepatocyte Ascl1 analysis (age / subcluster)
#   5. Count-based dot plots
#   6. Ascl1+ subset characterization
#
# Input:
#   - rna_wnn.h5ad  (annotated AnnData with cell_type, subcluster, sex, age)
#
# Output:
#   - Figures in figures/
#   - Tables in results/
#
#
# ==============================================================================

import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns

warnings.filterwarnings("ignore")

# Global plot settings
mpl.rcParams["font.family"] = "Arial"
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42


# ==============================================================================
# CONFIGURATION - UPDATE THIS PATH
# ==============================================================================

H5AD_PATH = "rna_wnn.h5ad"

FIG_DIR = Path("figures")
RES_DIR = Path("results")
FIG_DIR.mkdir(exist_ok=True)
RES_DIR.mkdir(exist_ok=True)

# Column names
CELLTYPE_COL  = "cell_type"
CELLTYPE2_COL = "subcluster"
SEX_COL       = "sex"
AGE_COL       = "age"

# Ordering
AGE_ORDER = ["young", "mid_age", "old", "pre_geriatric", "geriatric"]
AGE_ALIASES = {"midage": "mid_age", "pregeriatric": "pre_geriatric"}

GENE = "Ascl1"


# ==============================================================================
# HELPERS
# ==============================================================================

def banner(text):
    line = "=" * 70
    print()
    print(line)
    print(text)
    print(line)


def harmonize_age(adata):
    """Harmonize age labels and set ordered categorical."""
    age = adata.obs[AGE_COL].astype(str).str.strip().str.lower().replace(AGE_ALIASES)
    present = [a for a in AGE_ORDER if a in age.unique()]
    adata.obs[AGE_COL] = pd.Categorical(age, categories=present, ordered=True)
    return adata


def extract_gene_expression(adata, gene):
    """Return 1-D expression vector for a gene (handles sparse/dense)."""
    expr = adata[:, gene].X
    if hasattr(expr, "toarray"):
        expr = expr.toarray()
    return np.asarray(expr).ravel()


def summarize_by(adata, by_cols):
    """Grouped summary of Ascl1 positivity + mean expression."""
    if isinstance(by_cols, str):
        by_cols = [by_cols]
    df = adata.obs.groupby(by_cols, observed=True).agg(
        total=("Ascl1_positive", "size"),
        positive=("Ascl1_positive", "sum"),
        mean_expr=("Ascl1_expr", "mean"),
    ).reset_index()
    df["pct_positive"] = 100 * df["positive"] / df["total"]
    return df


# ==============================================================================
# 1. LOAD DATA
# ==============================================================================

banner("STEP 1: Load data")

print(f"  Loading {H5AD_PATH}")
adata = sc.read_h5ad(H5AD_PATH)
print(f"  [OK] {adata.n_obs:,} cells x {adata.n_vars:,} genes")

adata = harmonize_age(adata)

if GENE not in adata.var_names:
    raise ValueError(f"{GENE} not found in var_names")

adata.obs["Ascl1_expr"] = extract_gene_expression(adata, GENE)
adata.obs["Ascl1_positive"] = adata.obs["Ascl1_expr"] > 0


# ==============================================================================
# 2. GLOBAL ASCL1 PROFILING
# ==============================================================================

banner("STEP 2: Global Ascl1 profiling")

total_cells = adata.n_obs
ascl1_pos = int(adata.obs["Ascl1_positive"].sum())
print(f"\n  Overall:")
print(f"    Total cells:  {total_cells:,}")
print(f"    Ascl1+ cells: {ascl1_pos:,} ({100 * ascl1_pos / total_cells:.1f}%)")

sex_summary     = summarize_by(adata, SEX_COL)
age_summary     = summarize_by(adata, AGE_COL)
ct_summary      = summarize_by(adata, CELLTYPE_COL).sort_values("pct_positive", ascending=False)
sex_age_summary = summarize_by(adata, [SEX_COL, AGE_COL])

print(f"\n  By Sex:\n{sex_summary.to_string(index=False)}")
print(f"\n  By Age:\n{age_summary.to_string(index=False)}")
print(f"\n  By Celltype:\n{ct_summary.to_string(index=False)}")
print(f"\n  By Sex x Age:\n{sex_age_summary.to_string(index=False)}")

sex_summary.to_csv(RES_DIR / "Ascl1_summary_by_sex.csv", index=False)
age_summary.to_csv(RES_DIR / "Ascl1_summary_by_age.csv", index=False)
ct_summary.to_csv(RES_DIR / "Ascl1_summary_by_celltype.csv", index=False)
sex_age_summary.to_csv(RES_DIR / "Ascl1_summary_by_sex_age.csv", index=False)

fig, axes = plt.subplots(1, 3, figsize=(13, 4))
sc.pl.dotplot(adata, var_names=[GENE], groupby=SEX_COL,
              ax=axes[0], show=False, cmap="Reds")
axes[0].set_title("Ascl1 by Sex")
sc.pl.dotplot(adata, var_names=[GENE], groupby=AGE_COL,
              ax=axes[1], show=False, cmap="Reds")
axes[1].set_title("Ascl1 by Age")
sc.pl.dotplot(adata, var_names=[GENE], groupby=CELLTYPE_COL,
              ax=axes[2], show=False, cmap="Reds")
axes[2].set_title("Ascl1 by Cell Type")
plt.tight_layout()
fig.savefig(FIG_DIR / "Ascl1_dotplots_all_sex_age_celltype.pdf",
            dpi=300, bbox_inches="tight")
plt.close(fig)
print(f"\n  [OK] {FIG_DIR}/Ascl1_dotplots_all_sex_age_celltype.pdf")


# ==============================================================================
# 3. FEMALE-ONLY ASCL1 ANALYSIS
# ==============================================================================

banner("STEP 3: Female-only Ascl1 analysis")

adata_f = adata[adata.obs[SEX_COL] == "female"].copy()
adata_f.obs["Ascl1_expr"] = extract_gene_expression(adata_f, GENE)
adata_f.obs["Ascl1_positive"] = adata_f.obs["Ascl1_expr"] > 0
adata_f = harmonize_age(adata_f)

print(f"  Cells: {adata_f.n_obs:,}")
total = adata_f.n_obs
pos = int(adata_f.obs["Ascl1_positive"].sum())
print(f"  Ascl1+: {pos:,} ({100 * pos / total:.1f}%)")

age_summary_f = summarize_by(adata_f, AGE_COL)
ct_summary_f  = summarize_by(adata_f, CELLTYPE_COL).sort_values("pct_positive", ascending=False)

print(f"\n  By Age (Female):\n{age_summary_f.to_string(index=False)}")
print(f"\n  By Celltype (Female):\n{ct_summary_f.to_string(index=False)}")

age_summary_f.to_csv(RES_DIR / "Ascl1_summary_female_by_age.csv", index=False)
ct_summary_f.to_csv(RES_DIR / "Ascl1_summary_female_by_celltype.csv", index=False)

fig, axes = plt.subplots(1, 2, figsize=(10, 4))
sc.pl.dotplot(adata_f, var_names=[GENE], groupby=AGE_COL,
              ax=axes[0], show=False, cmap="Reds")
axes[0].set_title("Ascl1 by Age (Female)")
sc.pl.dotplot(adata_f, var_names=[GENE], groupby=CELLTYPE_COL,
              ax=axes[1], show=False, cmap="Reds")
axes[1].set_title("Ascl1 by Cell Type (Female)")
plt.tight_layout()
fig.savefig(FIG_DIR / "Ascl1_dotplots_female_only_age_celltype.pdf",
            dpi=300, bbox_inches="tight")
plt.close(fig)
print(f"\n  [OK] {FIG_DIR}/Ascl1_dotplots_female_only_age_celltype.pdf")


# ==============================================================================
# 4. FEMALE HEPATOCYTE ASCL1 ANALYSIS
# ==============================================================================

banner("STEP 4: Female hepatocyte Ascl1 analysis (age + subcluster)")

adata_fh = adata[
    (adata.obs[SEX_COL] == "female")
    & (adata.obs[CELLTYPE_COL] == "hepatocyte")
].copy()

if CELLTYPE2_COL not in adata_fh.obs.columns:
    raise ValueError(f"adata.obs does not contain '{CELLTYPE2_COL}' column")

adata_fh.obs["Ascl1_expr"] = extract_gene_expression(adata_fh, GENE)
adata_fh.obs["Ascl1_positive"] = adata_fh.obs["Ascl1_expr"] > 0
adata_fh = harmonize_age(adata_fh)

print(f"  Cells: {adata_fh.n_obs:,}")
total = adata_fh.n_obs
pos = int(adata_fh.obs["Ascl1_positive"].sum())
print(f"  Ascl1+: {pos:,} ({100 * pos / total:.1f}%)")

age_summary_fh = summarize_by(adata_fh, AGE_COL)
ct2_summary    = summarize_by(adata_fh, CELLTYPE2_COL).sort_values("pct_positive", ascending=False)

print(f"\n  By Age (Female hepatocytes):\n{age_summary_fh.to_string(index=False)}")
print(f"\n  By subcluster:\n{ct2_summary.to_string(index=False)}")

age_summary_fh.to_csv(RES_DIR / "Ascl1_summary_female_hep_by_age.csv", index=False)
ct2_summary.to_csv(RES_DIR / "Ascl1_summary_female_hep_by_celltype2.csv", index=False)

fig, axes = plt.subplots(1, 2, figsize=(10, 4))
sc.pl.dotplot(adata_fh, var_names=[GENE], groupby=AGE_COL,
              ax=axes[0], show=False, cmap="Reds")
axes[0].set_title("Ascl1 by Age (Female hepatocytes)")
sc.pl.dotplot(adata_fh, var_names=[GENE], groupby=CELLTYPE2_COL,
              ax=axes[1], show=False, cmap="Reds")
axes[1].set_title("Ascl1 by subcluster (Female hepatocytes)")
plt.tight_layout()
fig.savefig(FIG_DIR / "Ascl1_dotplots_female_hepatocyte_age_celltype2.pdf",
            dpi=300, bbox_inches="tight")
plt.close(fig)
print(f"\n  [OK] {FIG_DIR}/Ascl1_dotplots_female_hepatocyte_age_celltype2.pdf")


# ==============================================================================
# 5. COUNT-BASED DOT PLOTS
# ==============================================================================

banner("STEP 5: Count-based dot plots (Ascl1+ cell counts)")

age_rev = list(reversed(AGE_ORDER))
age_rev_present = [a for a in age_rev if a in adata_fh.obs[AGE_COL].unique()]
adata_fh_rev = adata_fh.copy()
adata_fh_rev.obs[AGE_COL] = pd.Categorical(
    adata_fh_rev.obs[AGE_COL].astype(str),
    categories=age_rev_present, ordered=True,
)

age_counts = (
    adata_fh_rev.obs.groupby(AGE_COL, observed=True)["Ascl1_positive"]
    .sum().reset_index(name="Ascl1_count")
)
ct2_counts = (
    adata_fh_rev.obs.groupby(CELLTYPE2_COL, observed=True)["Ascl1_positive"]
    .sum().reset_index(name="Ascl1_count")
)

fig, axes = plt.subplots(1, 2, figsize=(8, 4))

sns.scatterplot(
    data=age_counts, x=[GENE] * len(age_counts), y=AGE_COL,
    size="Ascl1_count", hue="Ascl1_count",
    palette="Reds", sizes=(50, 600),
    edgecolor="black", linewidth=0.6, ax=axes[0],
)
axes[0].set_title("Ascl1+ count by Age (Female Hep)", weight="bold", fontsize=12)
axes[0].set_xlabel("")
axes[0].set_ylabel("Age")
axes[0].invert_yaxis()
axes[0].legend(title="Ascl1+ count", loc="center left",
               bbox_to_anchor=(1.15, 0.5), frameon=True)

sns.scatterplot(
    data=ct2_counts, x=[GENE] * len(ct2_counts), y=CELLTYPE2_COL,
    size="Ascl1_count", hue="Ascl1_count",
    palette="Reds", sizes=(50, 600),
    edgecolor="black", linewidth=0.6, ax=axes[1],
)
axes[1].set_title("Ascl1+ count by subcluster", weight="bold", fontsize=12)
axes[1].set_xlabel("")
axes[1].set_ylabel("Subcluster")
axes[1].legend(title="Ascl1+ count", loc="center left",
               bbox_to_anchor=(1.15, 0.5), frameon=True)

plt.tight_layout()
fig.savefig(FIG_DIR / "Ascl1_count_dotplots_female_hepatocyte.pdf",
            dpi=300, bbox_inches="tight")
plt.close(fig)
print(f"  [OK] {FIG_DIR}/Ascl1_count_dotplots_female_hepatocyte.pdf")


# ==============================================================================
# 6. ASCL1+ SUBSET CHARACTERIZATION
# ==============================================================================

banner("STEP 6: Ascl1+ subset characterization")

adata_fh_pos = adata_fh[adata_fh.obs["Ascl1_positive"]].copy()

print(f"  Female hepatocytes total:  {adata_fh.n_obs:,}")
print(f"  Ascl1+ female hepatocytes: {adata_fh_pos.n_obs:,}")
print(f"  Fraction:                  {adata_fh_pos.n_obs / adata_fh.n_obs:.3f}")

print(f"\n  Ascl1 expression in Ascl1+ cells:")
print(f"    Mean:   {adata_fh_pos.obs['Ascl1_expr'].mean():.3f}")
print(f"    Median: {adata_fh_pos.obs['Ascl1_expr'].median():.3f}")
print(f"    Std:    {adata_fh_pos.obs['Ascl1_expr'].std():.3f}")


# ==============================================================================
# DONE
# ==============================================================================

banner("ALL STEPS COMPLETE")
print(f"  Figures: {FIG_DIR}/")
print(f"  Tables:  {RES_DIR}/")
