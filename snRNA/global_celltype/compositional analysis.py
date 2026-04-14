#!/usr/bin/env python3
# ==============================================================================
# Global Cell Type Analysis Across Aging
# ==============================================================================
#
# Description:
#   Global (all cell types) visualization and compositional analysis of the
#   liver aging multi-ome dataset, stratified by sex. Complements
#   hepatocyte_subcluster_analysis.py, which handles the hepatocyte
#   sub-cluster deep dive (celltype2 + zonation).
#
# Input:
#   - Annotated AnnData (.h5ad) with celltype, sex, age, sample labels and
#     WNN UMAP coordinates
#
# Output:
#   Figures:
#     - umap_WNN_age.pdf           (global UMAP colored by age)
#     - umap_WNN_celltype.pdf      (global UMAP colored by cell type)
#     - boxplot_celltype_by_age_male_female.pdf
#         (panel-style figure: male top, female bottom, one-way ANOVA + BH FDR,
#          * p<0.05, ** p<0.01, *** p<0.001, n.s. non-significant,
#          red dots = outliers)
#     - stacked_bar_celltype_pct_{sex}.pdf
#     - stacked_bar_celltype_counts_{sex}.pdf
#
#   Tables:
#     - celltype_percentages_per_sample.csv
#     - celltype_age_ANOVA_by_sex.csv
#     - {sex}_celltype_percentages.csv
#     - {sex}_celltype_anova.csv
#     - stacked_celltype_proportions_{sex}.csv
#
# Pipeline:
#   0. Load data
#   1. Global WNN UMAP by age
#   2. Global WNN UMAP by cell type
#   3. Cell type proportions and ANOVA by sex
#   4. Boxplots of cell type proportions by age
#   5. Stacked barplots of cell type proportions by age (percentage)
#   6. Stacked barplots of cell type counts by age
#
#
# ==============================================================================

import os
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.patches import Patch
from scipy.stats import f_oneway
from statsmodels.stats.multitest import multipletests

# Global plot settings
mpl.rcParams["font.family"] = "Arial"
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42


# ==============================================================================
# CONFIGURATION - UPDATE THIS PATH
# ==============================================================================

DATA_PATH = "rna_wnn.h5ad"
FIGURES_DIR = "figures"
RESULTS_DIR = "results"

os.makedirs(FIGURES_DIR, exist_ok=True)
os.makedirs(RESULTS_DIR, exist_ok=True)

# Column names
SAMPLE_COL = "sample"
CELLTYPE_COL = "celltype"
SEX_COL = "sex"
AGE_COL = "age"
WNN_BASIS = "X_wnn"

# Ordering
AGE_ORDER = ["young", "mid_age", "old", "pre_geriatric", "geriatric"]
AGE_ALIASES = {"midage": "mid_age", "pregeriatric": "pre_geriatric"}
SEX_ORDER = ["male", "female"]

# Palettes
AGE_PALETTE = {
    "young":         "#1ABC9C",
    "mid_age":       "#F1C40F",
    "old":           "#C39BD3",
    "pre_geriatric": "#2980B9",
    "geriatric":     "#E84393",
}


# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================

def pval_to_star(p):
    if pd.isna(p):
        return ""
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    return ""


# ==============================================================================
# 0. LOAD DATA
# ==============================================================================

print()
print("=" * 70)
print("ANALYSIS 0: Load data")
print("=" * 70)

adata = ad.read_h5ad(DATA_PATH)

age = adata.obs[AGE_COL].astype(str).str.strip().str.lower().replace(AGE_ALIASES)
adata.obs[AGE_COL] = pd.Categorical(age, categories=AGE_ORDER, ordered=True)

print(f"  Full dataset: {adata.n_obs:,} cells x {adata.n_vars:,} genes")
print(f"  Cell types: {sorted(adata.obs[CELLTYPE_COL].unique().tolist())}")
print(f"  Sex: {dict(adata.obs[SEX_COL].value_counts())}")
print(f"  Age: {dict(adata.obs[AGE_COL].value_counts())}")

if WNN_BASIS not in adata.obsm:
    raise KeyError(f"{WNN_BASIS} not found in adata.obsm")

present_ages = [a for a in AGE_ORDER if a in adata.obs[AGE_COL].cat.categories]


# ==============================================================================
# 1. GLOBAL WNN UMAP BY AGE
# ==============================================================================

print()
print("=" * 70)
print("ANALYSIS 1: Global WNN UMAP by age")
print("=" * 70)

adata.uns[f"{AGE_COL}_colors"] = [AGE_PALETTE[a] for a in present_ages]
adata.obsm["X_umap"] = adata.obsm[WNN_BASIS]

sc.pl.umap(
    adata,
    color=AGE_COL,
    title="Age (WNN UMAP)",
    legend_loc=None,
    frameon=False,
    size=10,
    show=False,
)
handles = [Patch(color=AGE_PALETTE[a], label=a) for a in present_ages]
plt.legend(
    handles=handles, title="Age Group",
    loc="center left", bbox_to_anchor=(1.02, 0.5),
    frameon=False, borderaxespad=0.0,
)
plt.tight_layout()
fname = f"{FIGURES_DIR}/umap_WNN_age.pdf"
plt.savefig(fname, bbox_inches="tight", dpi=300)
plt.close()
print(f"  [OK] {fname}")


# ==============================================================================
# 2. GLOBAL WNN UMAP BY CELL TYPE
# ==============================================================================

print()
print("=" * 70)
print("ANALYSIS 2: Global WNN UMAP by cell type")
print("=" * 70)

if adata.obs[CELLTYPE_COL].dtype.name != "category":
    adata.obs[CELLTYPE_COL] = adata.obs[CELLTYPE_COL].astype("category")

sc.pl.embedding(
    adata, basis=WNN_BASIS, color=CELLTYPE_COL,
    title="WNN UMAP - Cell Types",
    size=3, legend_loc="on data", frameon=False, show=False,
)
plt.tight_layout()
fname = f"{FIGURES_DIR}/umap_WNN_celltype.pdf"
plt.savefig(fname, bbox_inches="tight", dpi=300)
plt.close()
print(f"  [OK] {fname}")


# ==============================================================================
# 3. CELL TYPE PROPORTIONS AND ANOVA BY SEX
# ==============================================================================

print()
print("=" * 70)
print("ANALYSIS 3: Cell type proportions and ANOVA")
print("=" * 70)

# Global percentages per sample
counts_all = (
    adata.obs.groupby([SAMPLE_COL, CELLTYPE_COL])
    .size().reset_index(name="cell_count")
)
counts_all["total"] = counts_all.groupby(SAMPLE_COL)["cell_count"].transform("sum")
counts_all["percentage"] = (counts_all["cell_count"] / counts_all["total"]) * 100

pivot = counts_all.pivot_table(
    index=SAMPLE_COL, columns=CELLTYPE_COL, values="percentage", fill_value=0,
)
out = f"{RESULTS_DIR}/celltype_percentages_per_sample.csv"
pivot.to_csv(out)
print(f"  [OK] {out}")

# ANOVA by sex
meta = adata.obs[[SAMPLE_COL, AGE_COL, SEX_COL]].drop_duplicates().set_index(SAMPLE_COL)

anova_frames = []
for sex in SEX_ORDER:
    if sex not in adata.obs[SEX_COL].unique():
        continue
    sub = adata[adata.obs[SEX_COL] == sex].copy()
    ct_counts = (
        sub.obs.groupby([SAMPLE_COL, CELLTYPE_COL])
        .size().reset_index(name="cell_count")
    )
    ct_counts["total"] = ct_counts.groupby(SAMPLE_COL)["cell_count"].transform("sum")
    ct_counts["percentage"] = (ct_counts["cell_count"] / ct_counts["total"]) * 100
    ct_counts[AGE_COL] = ct_counts[SAMPLE_COL].map(meta[AGE_COL])

    results = []
    for ct in ct_counts[CELLTYPE_COL].unique():
        ct_sub = ct_counts[ct_counts[CELLTYPE_COL] == ct]
        groups = [g["percentage"].values for _, g in ct_sub.groupby(AGE_COL) if len(g) >= 2]
        if len(groups) >= 2:
            stat, p = f_oneway(*groups)
        else:
            stat, p = np.nan, np.nan
        results.append({CELLTYPE_COL: ct, "anova_stat": stat, "p_value": p})

    df_res = pd.DataFrame(results)
    df_res["adj_p"] = multipletests(df_res["p_value"].fillna(1), method="fdr_bh")[1]
    df_res["significance"] = df_res["adj_p"].apply(pval_to_star)
    df_res["sex"] = sex
    anova_frames.append(df_res)

anova_combined = pd.concat(anova_frames, ignore_index=True)
out = f"{RESULTS_DIR}/celltype_age_ANOVA_by_sex.csv"
anova_combined.to_csv(out, index=False)
print(f"  [OK] {out}")
print(anova_combined.sort_values("adj_p").head(20).to_string(index=False))


# ==============================================================================
# 4. BOXPLOTS - CELL TYPE PROPORTIONS BY AGE
# ==============================================================================

print()
print("=" * 70)
print("ANALYSIS 4: Boxplots (celltype proportions by age)")
print("=" * 70)

#  Per-sex counts + ANOVA, then a single combined figure
#  (male top, female bottom) matching the paper's panel (a) legend.

sex_results = {}

for sex in SEX_ORDER:
    if sex not in adata.obs[SEX_COL].unique():
        continue

    sub = adata[adata.obs[SEX_COL] == sex].copy()
    ct_counts = (
        sub.obs.groupby([SAMPLE_COL, CELLTYPE_COL])
        .size().reset_index(name="cell_count")
    )
    ct_counts["total"] = ct_counts.groupby(SAMPLE_COL)["cell_count"].transform("sum")
    ct_counts["percentage"] = (ct_counts["cell_count"] / ct_counts["total"]) * 100
    ct_counts[AGE_COL] = ct_counts[SAMPLE_COL].map(meta[AGE_COL])
    ct_counts[SEX_COL] = sex

    # Per-celltype one-way ANOVA, FDR-BH corrected
    anova_rows = []
    for ct in ct_counts[CELLTYPE_COL].unique():
        ct_sub = ct_counts[ct_counts[CELLTYPE_COL] == ct]
        groups = [g["percentage"].values for _, g in ct_sub.groupby(AGE_COL) if len(g) >= 2]
        if len(groups) >= 2:
            stat, p = f_oneway(*groups)
        else:
            stat, p = np.nan, np.nan
        anova_rows.append({CELLTYPE_COL: ct, "anova_stat": stat, "p_value": p})
    anova_df = pd.DataFrame(anova_rows)
    anova_df["adj_p"] = multipletests(anova_df["p_value"].fillna(1), method="fdr_bh")[1]
    anova_df["significance"] = anova_df["adj_p"].apply(pval_to_star)

    ct_counts = ct_counts.merge(
        anova_df[[CELLTYPE_COL, "significance"]], on=CELLTYPE_COL, how="left",
    )
    ct_counts["significance"] = ct_counts["significance"].fillna("")

    sex_results[sex] = {"counts": ct_counts, "anova": anova_df}

    # Per-sex tables
    ct_counts.to_csv(f"{RESULTS_DIR}/{sex}_celltype_percentages.csv", index=False)
    anova_df.to_csv(f"{RESULTS_DIR}/{sex}_celltype_anova.csv", index=False)

# Combined figure - male top, female bottom, shared cell-type order
ct_order = sorted(
    set().union(*[sex_results[s]["counts"][CELLTYPE_COL].unique()
                  for s in sex_results])
)
global_present_ages = [
    a for a in AGE_ORDER
    if any(a in sex_results[s]["counts"][AGE_COL].dropna().unique()
           for s in sex_results)
]

# Red-dot outlier styling per the figure legend
flierprops = dict(
    marker="o", markerfacecolor="red", markeredgecolor="red",
    markersize=4, linestyle="none",
)

panel_sexes = [s for s in ["male", "female"] if s in sex_results]
fig, axes = plt.subplots(
    len(panel_sexes), 1, figsize=(14, 5 * len(panel_sexes)),
    sharex=True, sharey=False,
)
if len(panel_sexes) == 1:
    axes = [axes]

ymax_global = max(
    sex_results[s]["counts"]["percentage"].max() for s in panel_sexes
)

for ax, sex in zip(axes, panel_sexes):
    ct_counts = sex_results[sex]["counts"]
    sex_present_ages = [
        a for a in AGE_ORDER if a in ct_counts[AGE_COL].dropna().unique()
    ]

    sns.boxplot(
        data=ct_counts, x=CELLTYPE_COL, y="percentage", hue=AGE_COL,
        order=ct_order, hue_order=sex_present_ages,
        palette={k: AGE_PALETTE[k] for k in sex_present_ages},
        showcaps=True, flierprops=flierprops, dodge=True, ax=ax,
    )
    ax.set_ylim(0, ymax_global * 1.2)

    # Significance annotation per celltype
    for i, ct in enumerate(ct_order):
        sub_ct = ct_counts[ct_counts[CELLTYPE_COL] == ct]
        if sub_ct.empty:
            continue
        star = sub_ct["significance"].iloc[0]
        label, color = (star, "red") if star else ("n.s.", "darkblue")
        ax.text(
            i, sub_ct["percentage"].max() * 1.08, label,
            ha="center", va="bottom", fontsize=16, color=color,
        )

    # Legend only on the top panel
    if sex == panel_sexes[0]:
        handles = [
            Patch(facecolor=AGE_PALETTE[a], edgecolor="black", label=a)
            for a in global_present_ages
        ]
        ax.legend(
            handles=handles, title="Age Group",
            bbox_to_anchor=(1.01, 1), loc="upper left", frameon=False,
        )
    else:
        leg = ax.get_legend()
        if leg is not None:
            leg.remove()

    ax.set_title(sex.capitalize(), fontsize=14, fontweight="bold", loc="left")
    ax.set_ylabel("% of cells per sample", fontsize=12)
    ax.set_xlabel("")

axes[-1].set_xlabel("Cell Type", fontsize=12)
plt.setp(axes[-1].get_xticklabels(), rotation=45, ha="right")
plt.tight_layout()

fname = f"{FIGURES_DIR}/boxplot_celltype_by_age_male_female.pdf"
fig.savefig(fname, bbox_inches="tight", dpi=300)
plt.close(fig)
print(f"  [OK] {fname}")


# ==============================================================================
# 5 & 6. STACKED BARPLOTS - CELL TYPE BY AGE (PERCENTAGE AND COUNTS)
# ==============================================================================

print()
print("=" * 70)
print("ANALYSIS 5 & 6: Stacked barplots (celltype by age)")
print("=" * 70)

data = adata.obs[[AGE_COL, SEX_COL, CELLTYPE_COL]].copy()
data[AGE_COL] = pd.Categorical(data[AGE_COL], categories=AGE_ORDER, ordered=True)


def prepare_stacked(df, sex, value_type="percentage"):
    """Return a (age x celltype) pivot of percentages or raw counts."""
    sub = df[df[SEX_COL] == sex]
    counts = sub.groupby([AGE_COL, CELLTYPE_COL]).size().reset_index(name="count")
    if value_type == "percentage":
        counts["value"] = (
            counts["count"]
            / counts.groupby(AGE_COL)["count"].transform("sum") * 100
        )
    else:
        counts["value"] = counts["count"]
    counts = counts.sort_values(AGE_COL)
    stacked = counts.pivot(index=AGE_COL, columns=CELLTYPE_COL, values="value").fillna(0)
    return stacked


def plot_stacked(stacked, sex, value_type, out_pdf):
    fig, ax = plt.subplots(figsize=(12, 8))
    stacked.plot(
        kind="bar", stacked=True,
        colormap="tab20", edgecolor="none", width=1.0, ax=ax,
    )
    ylabel = "Percentage (%)" if value_type == "percentage" else "Cell Count"
    title_prefix = ("Distribution of Cell Types by Age"
                    if value_type == "percentage"
                    else "Cell Type Counts by Age")
    ax.set_title(f"{title_prefix} ({sex.capitalize()}s)", fontsize=16)
    ax.set_ylabel(ylabel, fontsize=14)
    ax.set_xlabel("Age Groups", fontsize=14)
    plt.xticks(rotation=45, ha="right", fontsize=12)
    ax.legend(title="Cell Type",
              bbox_to_anchor=(1.05, 1), loc="upper left", fontsize=10)
    plt.tight_layout()
    fig.savefig(out_pdf, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  [OK] {out_pdf}")


for sex in SEX_ORDER:
    if sex not in adata.obs[SEX_COL].unique():
        continue

    # Percentage
    stacked_pct = prepare_stacked(data, sex, value_type="percentage")
    plot_stacked(
        stacked_pct, sex, "percentage",
        f"{FIGURES_DIR}/stacked_bar_celltype_pct_{sex}.pdf",
    )
    stacked_pct.to_csv(f"{RESULTS_DIR}/stacked_celltype_proportions_{sex}.csv")

    # Counts
    stacked_cnt = prepare_stacked(data, sex, value_type="count")
    plot_stacked(
        stacked_cnt, sex, "count",
        f"{FIGURES_DIR}/stacked_bar_celltype_counts_{sex}.pdf",
    )


# ==============================================================================
# SUMMARY
# ==============================================================================

print()
print("=" * 70)
print("ALL ANALYSES COMPLETE")
print("=" * 70)
print(f"  Figures: {FIGURES_DIR}/")
print(f"  Tables:  {RESULTS_DIR}/")
print("=" * 70)
