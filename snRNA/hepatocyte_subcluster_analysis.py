#!/usr/bin/env python3
# ==============================================================================
# Hepatocyte Sub-Cluster and Zonation Analysis Across Aging
# ==============================================================================
#
# Description:
#   Analyzes hepatocyte sub-cluster composition and zonation patterns across
#   aging in mouse liver, stratified by sex. Includes differential expression,
#   GSEA prerank (Reactome), compositional analysis with ANOVA, gene module
#   scoring for zonation assignment, and cell-level + sample-level zonation
#   comparisons across age groups.
#
# Input:
#   - Annotated AnnData (.h5ad) with celltype, celltype2 (sub-clusters),
#     sex, age, sample labels, and WNN UMAP coordinates
#
# Output:
#   Figures:
#     - WNN-UMAP by sub-cluster + per-cluster highlight grids
#     - DE dotplot (Wilcoxon)
#     - GSEA Reactome bubble plot
#     - Sub-cluster proportion boxplots and stacked barplots by age
#     - Zonation module score heatmap
#     - Zonation UMAP grid, boxplots, stacked barplots
#     - Cell-level zonation violin plots with ANOVA
#     - Sample-level zonation boxplots with ANOVA
#
#   Tables:
#     - DE results (all genes, per cluster)
#     - GSEA results per cluster
#     - Sub-cluster proportions and ANOVA
#     - Zonation proportions and ANOVA
#     - Sample-level zonation scores
#
# Pipeline:
#   0. Load and subset to Hepatocytes
#   1. WNN-UMAP by celltype2 + per-cluster highlight grids
#   2. Differential expression (Wilcoxon) + dotplot
#   3. GSEA prerank (Reactome) + bubble plot
#   4. Celltype2 proportions and ANOVA by sex
#   5. Boxplots of celltype2 proportions by age
#   5b. Stacked barplots of celltype2 proportions by age
#   6. Gene module scoring and heatmap (global)
#   7. Hepatocyte zonation: assignment, UMAP grid, boxplots, stacked bars
#   8. Per-subset zonation violin plots (cell-level) + ANOVA
#   9. Per-subset zonation boxplots (sample-level means) + ANOVA
#
# Reference:
#   Sarker et al. (2026) Cell Metabolism
#
# ==============================================================================

import os
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import gseapy as gp
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from scipy.stats import f_oneway
from statsmodels.stats.multitest import multipletests

# Global plot settings
mpl.rcParams["font.family"] = "Arial"
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42


# ==============================================================================
# CONFIGURATION - UPDATE THESE PATHS
# ==============================================================================

DATA_PATH = "integrated_scvi.h5ad"
FIGURES_DIR = "figures"
RESULTS_DIR = "results"

os.makedirs(FIGURES_DIR, exist_ok=True)
os.makedirs(RESULTS_DIR, exist_ok=True)

# Column names
SAMPLE_COL = "sample"
CELLTYPE_COL = "celltype"
CELLTYPE2_COL = "celltype2"
SEX_COL = "sex"
AGE_COL = "age"
ZONATION_COL = "hepatocyte_zonation"

# Ordering
AGE_ORDER = ["young", "mid_age", "old", "pre_geriatric", "geriatric"]
AGE_ALIASES = {"midage": "mid_age", "pregeriatric": "pre_geriatric"}
SEX_ORDER = ["male", "female"]
HEP_ORDER = [f"Hep-0{i}" for i in range(1, 8)]
ZONE_ORDER = ["Periportal", "Midlobular", "Pericentral"]

# Zonation mapping (sub-cluster to zone)
HEPATOCYTE_ZONATION_MAP = {
    "Hep-01": "Periportal",
    "Hep-02": "Midlobular",
    "Hep-03": "Pericentral",
    "Hep-04": "Pericentral",
    "Hep-05": "Periportal",
    "Hep-06": "Midlobular",
    "Hep-07": "Midlobular",
}

# Color palettes
AGE_PALETTE = {
    "young": "#1ABC9C",
    "mid_age": "#F1C40F",
    "old": "#C39BD3",
    "pre_geriatric": "#2980B9",
    "geriatric": "#E84393",
}

ZONATION_PALETTE = {
    "Periportal": "#66C2A5",
    "Midlobular": "#E7298A",
    "Pericentral": "#B2D33C",
}

CELLTYPE2_PALETTE = {
    "Hep-01": "#66C2A5",
    "Hep-02": "#FC8D62",
    "Hep-03": "#5E4FA2",
    "Hep-04": "#B2D33C",
    "Hep-05": "#8B0000",
    "Hep-06": "#00BFFF",
    "Hep-07": "#E7298A",
}

# Zonation gene list
ZONATION_GENE_LIST = [
    "Gm42418", "Glul", "Cyp2f2", "Cyp2e1", "Aldh1a1", "Cyp1a2", "Cyp2a5", "Ephx1",
    "Gsta1", "Gsta2", "Gstm1", "Gstm2", "Gstm3", "Gstm6", "Gstt1",
    "Abcc3", "Abcc4", "Abcc6", "Abcc9", "Abcg2", "Abcg5", "Abcg8",
    "Hpd", "Atp1a1", "Atp5g3", "Ndufa2", "Ndufa3", "Ndufa9",
    "Ndufb1", "Ndufb7", "Ndufs1",
    "Ldha", "Ldhb", "Ldhc", "Ldhf", "Ldhx",
    "Lpl", "Lipa", "Lipe", "Lpin1", "Lpin2", "Lpin3",
    "Hmgcs2", "Hmgcr", "Hmgcl", "Mvd", "Mvk",
    "Me1", "Me2", "Me3", "Mdh1", "Mdh2", "Mdh1b",
    "Mlxip", "Mlxipl", "Mlx", "Malat1", "Mat1a", "Mat2a", "Mat2b",
    "Mettl7a1", "Mettl7a2", "Mettl7b",
    "Otc", "Oat", "Oaz1", "Oaz2", "Oaz3", "Odc1",
    "Pah", "Psat1", "Phgdh", "Phka2", "Phkb",
    "Pklr", "Pkm", "Pfkl", "Pfkp", "Pcx", "Pdha1", "Pdha2", "Pdhb",
    "Pdk1", "Pdk2", "Pdk3", "Pdk4",
    "Ppp1r3b", "Pparg", "Ppara", "Ppard", "Ppp1r3c", "Ppp1r3a",
    "Rgn", "Rnr1", "Rnr2", "Rnr3", "Rnmt",
    "Rbm3", "Rbm4", "Rbm5", "Rbm6", "Rbm7", "Rbmx",
    "Serpina1a", "Serpina1b", "Serpina1c", "Serpina1d", "Serpina1e",
    "Slc1a1", "Slc1a2", "Slc1a3", "Slc1a4", "Slc1a5", "Slc1a6", "Slc1a7",
    "Slc2a1", "Slc2a2", "Slc2a3", "Slc2a4", "Slc2a5", "Slc2a6", "Slc2a7",
    "Slc2a8", "Slc2a9", "Slc2a10", "Slc2a11", "Slc2a12", "Slc2a13",
    "Slc2a14", "Slc2a15", "Slc2a16", "Slc2a17", "Slc2a18", "Slc2a19",
    "Slc2a20", "Slc2a21", "Slc2a22", "Slc2a23", "Slc2a24", "Slc2a25",
    "Slc2a26", "Slc2a27", "Slc2a28", "Slc2a29", "Slc2a30",
    "Ugt1a1", "Ugt1a6", "Ugt2b1", "Ugt2b34", "Ugt2b5", "Ugt2b7",
    "Ugt2b38", "Ugt3a1", "Ugt3a2", "Ugt8a",
    "Ugt6a1", "Ugt7a1", "Ugt7a2", "Ugt7a3",
    "Vim", "Vnn1", "Vnn2", "Vnn3",
    "Zfp36l1", "Zfp36l2", "Zfp36l3", "Zfp36",
    "Zhx1", "Zhx2", "Zhx3", "Uhrf1",
]


# ==============================================================================
# 0. LOAD DATA AND SUBSET TO HEPATOCYTES
# ==============================================================================

print()
print("=" * 70)
print("ANALYSIS 0: Load data and subset to Hepatocytes")
print("=" * 70)

adata = ad.read_h5ad(DATA_PATH)

# Harmonize age labels
age = adata.obs[AGE_COL].astype(str).str.strip().str.lower().replace(AGE_ALIASES)
adata.obs[AGE_COL] = pd.Categorical(age, categories=AGE_ORDER, ordered=True)

print(f"  Full dataset: {adata.n_obs:,} cells x {adata.n_vars:,} genes")

adata = adata[adata.obs[CELLTYPE_COL] == "Hepatocyte"].copy()
print(f"  Hepatocytes: {adata.n_obs:,} cells")
print(f"  Sub-clusters: {sorted(adata.obs[CELLTYPE2_COL].unique().tolist())}")
print(f"  Sex: {dict(adata.obs[SEX_COL].value_counts())}")
print(f"  Age: {dict(adata.obs[AGE_COL].value_counts())}")


# ==============================================================================
# 1. WNN-UMAP BY CELLTYPE2 + PER-CLUSTER HIGHLIGHT GRIDS
# ==============================================================================

print()
print("=" * 70)
print("ANALYSIS 1: WNN-UMAP by celltype2 + highlight grids")
print("=" * 70)

# 1a. Full celltype2 UMAP
if not adata.obs[CELLTYPE2_COL].dtype.name == "category":
    adata.obs[CELLTYPE2_COL] = adata.obs[CELLTYPE2_COL].astype("category")
adata.uns[f"{CELLTYPE2_COL}_colors"] = [
    CELLTYPE2_PALETTE[ct] for ct in adata.obs[CELLTYPE2_COL].cat.categories
]

sc.pl.embedding(
    adata, basis="X_wnn", color=CELLTYPE2_COL,
    title="WNN UMAP - Hepatocyte Sub-clusters",
    size=3, legend_loc="on data", show=False,
)
plt.tight_layout()
fname = f"{FIGURES_DIR}/umap_WNN_celltype2.pdf"
plt.savefig(fname, bbox_inches="tight", dpi=300)
plt.close()
print(f"  [OK] {fname}")

# 1b. Highlight grid for each sub-cluster
x_coords = adata.obsm["X_wnn"][:, 0]
y_coords = adata.obsm["X_wnn"][:, 1]
xlim_global = (x_coords.min(), x_coords.max())
ylim_global = (y_coords.min(), y_coords.max())

sorted_ages = [a for a in AGE_ORDER if a in adata.obs[AGE_COL].unique()]
sexes = adata.obs[SEX_COL].unique()

for highlight_ct in sorted(adata.obs[CELLTYPE2_COL].unique()):
    nrow, ncol = len(sexes), len(sorted_ages)
    fig, axes = plt.subplots(
        nrow, ncol, figsize=(5 * ncol, 5 * nrow), sharex=True, sharey=True,
    )
    if nrow == 1:
        axes = axes[np.newaxis, :]

    for col_idx, age_val in enumerate(sorted_ages):
        for row_idx, sex_val in enumerate(sexes):
            ax = axes[row_idx, col_idx]
            mask = (adata.obs[AGE_COL] == age_val) & (adata.obs[SEX_COL] == sex_val)
            subset = adata[mask].copy()

            subset.obs["highlight"] = pd.Categorical(
                np.where(
                    subset.obs[CELLTYPE2_COL] == highlight_ct,
                    highlight_ct, "Other",
                ),
                categories=[highlight_ct, "Other"],
            )

            sc.pl.embedding(
                subset, basis="X_wnn", color="highlight",
                palette=["#FF0000", "#F0F0F0"],
                title=f"{sex_val.capitalize()} - {age_val}",
                size=20, legend_loc=None, ax=ax, show=False,
            )
            ax.set_xlim(xlim_global)
            ax.set_ylim(ylim_global)

    plt.suptitle(f"Highlighted: {highlight_ct}", fontsize=18, y=1.01)
    plt.tight_layout()
    fname = f"{FIGURES_DIR}/umap_highlight_{highlight_ct}.pdf"
    fig.savefig(fname, bbox_inches="tight", dpi=300)
    plt.close(fig)
    print(f"  [OK] {fname}")


# ==============================================================================
# 2. DIFFERENTIAL EXPRESSION (CELLTYPE2) + DOTPLOT
# ==============================================================================

print()
print("=" * 70)
print("ANALYSIS 2: Differential expression (celltype2)")
print("=" * 70)

sc.tl.rank_genes_groups(adata, groupby=CELLTYPE2_COL, method="wilcoxon")

df_deg = sc.get.rank_genes_groups_df(adata, group=None, key="rank_genes_groups")
df_deg["names"] = df_deg["names"].astype(str)
df_deg["gene_excel_safe"] = "'" + df_deg["names"]

out = f"{RESULTS_DIR}/rank_genes_groups_celltype2_all.csv"
df_deg.to_csv(out, index=False)
print(f"  [OK] {out}")

sc.pl.rank_genes_groups_dotplot(
    adata, n_genes=5, values_to_plot="logfoldchanges",
    min_logfoldchange=1, vmax=4, vmin=-4, cmap="bwr", show=False,
)
fname = f"{FIGURES_DIR}/rank_genes_dotplot_celltype2.pdf"
plt.savefig(fname, bbox_inches="tight", dpi=300)
plt.savefig(fname.replace(".pdf", ".png"), bbox_inches="tight", dpi=300)
plt.close()
print(f"  [OK] {fname}")


# ==============================================================================
# 3. GSEA PRERANK (REACTOME) + BUBBLE PLOT
# ==============================================================================

print()
print("=" * 70)
print("ANALYSIS 3: GSEA prerank (Reactome) + bubble plot")
print("=" * 70)

# 3a. Save per-cluster DEG files
GSEA_DIR = f"{RESULTS_DIR}/GSEA"
os.makedirs(GSEA_DIR, exist_ok=True)

for cluster in sorted(df_deg["group"].unique()):
    cluster_df = df_deg[df_deg["group"] == cluster][
        ["names", "scores", "pvals", "pvals_adj", "logfoldchanges"]
    ].copy()
    cluster_df.columns = ["gene", "score", "pval", "pval_adj", "logfoldchange"]
    out_path = f"{GSEA_DIR}/differential_expression_{cluster}.csv"
    cluster_df.to_csv(out_path, index=False)
    print(f"  [OK] {out_path}")

# 3b. Run GSEA prerank for Hep-01 through Hep-07
for i in range(1, 8):
    hep = f"Hep-0{i}"
    file_path = f"{GSEA_DIR}/differential_expression_{hep}.csv"
    if not os.path.exists(file_path):
        print(f"  [SKIP] DEG file not found for {hep}")
        continue

    print(f"  Running GSEA prerank for {hep}...")
    df_gsea = pd.read_csv(file_path)

    ranked_genes = df_gsea[["gene", "score"]].dropna()
    ranked_genes.columns = ["Gene", "Score"]
    ranked_genes["Gene"] = ranked_genes["Gene"].str.upper()
    ranked_genes["Score"] += np.random.normal(0, 1e-6, len(ranked_genes))
    ranked_genes = ranked_genes.sort_values("Score", ascending=False)

    gsea_results = gp.prerank(
        rnk=ranked_genes,
        gene_sets="Reactome_Pathways_2024",
        organism="Human",
        permutation_num=1000,
        min_size=5,
        max_size=2000,
        outdir=f"{GSEA_DIR}/GSEA_results_0{i}",
        seed=42,
    )

    gsea_results_df = gsea_results.res2d

    # Extract Overlap_N and numeric Tag %
    if "Tag %" in gsea_results_df.columns:
        tag_raw = gsea_results_df["Tag %"].astype(str)
        if tag_raw.str.contains("/").any():
            parts = tag_raw.str.extract(r"(\d+)/(\d+)")
            gsea_results_df["Overlap_N"] = parts[0].astype(float)
            gsea_results_df["Tag %_numeric"] = (
                parts[0].astype(float) / parts[1].astype(float) * 100
            )
        else:
            gsea_results_df["Tag %_numeric"] = (
                tag_raw.str.replace("%", "", regex=False).astype(float)
            )

    # Convert Gene % to numeric
    if "Gene %" in gsea_results_df.columns:
        gene_raw = gsea_results_df["Gene %"].astype(str)
        if gene_raw.str.contains("/").any():
            parts = gene_raw.str.extract(r"(\d+)/(\d+)")
            gsea_results_df["Gene %_numeric"] = (
                parts[0].astype(float) / parts[1].astype(float) * 100
            )
        else:
            gsea_results_df["Gene %_numeric"] = (
                gene_raw.str.replace("%", "", regex=False).astype(float)
            )

    out_path = f"{GSEA_DIR}/GSEA_Reactome_2025_hep_0{i}.csv"
    gsea_results_df.to_csv(out_path, index=False)
    print(f"  [OK] {out_path}")

# 3c. Bubble plot - top 10 upregulated pathways per cluster
gsea_combined = []
for i in range(1, 8):
    cluster_id = f"Hep-0{i}"
    file_path = f"{GSEA_DIR}/GSEA_Reactome_2025_hep_0{i}.csv"
    if not os.path.exists(file_path):
        continue

    df_gp = pd.read_csv(file_path)

    if "Overlap_N" not in df_gp.columns and "Tag %" in df_gp.columns:
        tag_raw = df_gp["Tag %"].astype(str)
        if tag_raw.str.contains("/").any():
            df_gp["Overlap_N"] = tag_raw.str.extract(r"(\d+)/")[0].astype(float)
        else:
            df_gp["Overlap_N"] = np.nan

    if "Gene %_numeric" in df_gp.columns:
        df_gp["Gene_pct"] = df_gp["Gene %_numeric"]
    elif "Gene %" in df_gp.columns:
        gene_raw = df_gp["Gene %"].astype(str)
        if gene_raw.str.contains("/").any():
            parts = gene_raw.str.extract(r"(\d+)/(\d+)")
            df_gp["Gene_pct"] = parts[0].astype(float) / parts[1].astype(float) * 100
        else:
            df_gp["Gene_pct"] = gene_raw.str.replace("%", "", regex=False).astype(float)
    else:
        df_gp["Gene_pct"] = np.nan

    # Filter: NOM p-val < 0.05, Overlap_N > 4, NES > 0
    df_gp = df_gp[
        (df_gp["NOM p-val"] < 0.05)
        & (df_gp["Overlap_N"] > 4)
        & (df_gp["NES"] > 0)
    ]
    if df_gp.empty:
        print(f"  [SKIP] No pathways passed filters for {cluster_id}")
        continue

    df_gp = df_gp.sort_values("NES", ascending=False).head(10)
    df_gp["Cluster"] = cluster_id
    gsea_combined.append(df_gp)

if gsea_combined:
    plot_df = pd.concat(gsea_combined, ignore_index=True)

    MAX_CHARS = 60
    plot_df["y_label"] = plot_df["Term"].apply(
        lambda x: x[:MAX_CHARS] + "..." if len(str(x)) > MAX_CHARS else x
    )

    # Deduplicate identical pathway names across clusters
    seen = {}
    unique_labels = []
    for _, row in plot_df.iterrows():
        label = row["y_label"]
        if label in seen:
            seen[label] += 1
            label = label + " " * seen[label]
        else:
            seen[label] = 0
        unique_labels.append(label)
    plot_df["y_label"] = unique_labels

    plot_df["Cluster_rank"] = plot_df["Cluster"].map(
        {c: idx for idx, c in enumerate(HEP_ORDER)}
    )
    plot_df = plot_df.sort_values(
        ["Cluster_rank", "NES"], ascending=[True, False],
    ).reset_index(drop=True)

    SIZE_SCALE = 30
    plot_df["bubble_size"] = plot_df["Gene_pct"] * SIZE_SCALE

    fig, ax = plt.subplots(figsize=(6.5, max(6, 0.25 * len(plot_df))))

    ax.scatter(
        plot_df["NES"], plot_df["y_label"],
        s=plot_df["bubble_size"],
        c=plot_df["Cluster"].map(CELLTYPE2_PALETTE),
        edgecolors="black", linewidths=0.5, alpha=0.85,
        zorder=3, clip_on=False,
    )

    for _, row in plot_df.iterrows():
        ax.hlines(
            y=row["y_label"], xmin=0, xmax=row["NES"],
            color="gray", alpha=0.3, linestyle="--", linewidth=0.6, zorder=1,
        )

    prev_cluster = None
    for idx_row, row in plot_df.iterrows():
        if prev_cluster is not None and row["Cluster"] != prev_cluster:
            ax.axhline(
                y=idx_row - 0.5, color="lightgray", linewidth=1.2,
                linestyle="-", zorder=0,
            )
        prev_cluster = row["Cluster"]

    nes_min = plot_df["NES"].min()
    nes_max = plot_df["NES"].max()
    nes_range = nes_max - nes_min if nes_max > nes_min else 1.0
    ax.set_xlim(nes_min - 0.08 * nes_range, nes_max + 0.08 * nes_range)
    ax.set_ylim(len(plot_df) - 0.5, -0.5)
    ax.margins(y=0.02)

    ax.axvline(x=0, color="red", linestyle="--", linewidth=1, alpha=0.6)
    ax.set_xlabel("Normalized Enrichment Score (NES)", fontsize=12)
    ax.set_ylabel("")
    ax.set_title(
        "Top 10 Upregulated Reactome Pathways per Hepatocyte Sub-cluster\n"
        "(NOM p-val < 0.05, Overlap > 4, NES > 0)",
        fontsize=13, fontweight="bold",
    )
    ax.tick_params(axis="y", labelsize=8)
    sns.despine(ax=ax, left=True)

    plotted_clusters = plot_df["Cluster"].unique()
    cluster_handles = [
        Line2D(
            [0], [0], marker="o", color="w", label=c,
            markerfacecolor=CELLTYPE2_PALETTE[c], markersize=10,
            markeredgecolor="black", markeredgewidth=0.5,
        )
        for c in plotted_clusters
    ]

    gene_pct_values = [5, 10, 20]
    size_handles = [
        Line2D(
            [0], [0], marker="o", color="w", label=f"{v}%",
            markerfacecolor="gray", markeredgecolor="black",
            markeredgewidth=0.5, markersize=np.sqrt(v * SIZE_SCALE) * 0.8,
        )
        for v in gene_pct_values
    ]

    first_legend = ax.legend(
        handles=cluster_handles, title="Cluster",
        bbox_to_anchor=(1.02, 1), loc="upper left", frameon=False,
    )
    ax.add_artist(first_legend)
    ax.legend(
        handles=size_handles, title="Gene %",
        bbox_to_anchor=(1.02, 0.55), loc="upper left", frameon=False,
    )

    plt.tight_layout()
    fname = f"{FIGURES_DIR}/gsea_reactome_bubble_plot.pdf"
    fig.savefig(fname, bbox_inches="tight", dpi=300)
    fig.savefig(fname.replace(".pdf", ".png"), bbox_inches="tight", dpi=300)
    plt.close(fig)
    print(f"  [OK] {fname}")
else:
    print("  [SKIP] No clusters had significant pathways")


# ==============================================================================
# 4. CELLTYPE2 PROPORTIONS AND ANOVA BY SEX
# ==============================================================================

print()
print("=" * 70)
print("ANALYSIS 4: Celltype2 proportions and ANOVA")
print("=" * 70)

counts_all = (
    adata.obs.groupby([SAMPLE_COL, CELLTYPE2_COL])
    .size().reset_index(name="cell_count")
)
counts_all["total"] = counts_all.groupby(SAMPLE_COL)["cell_count"].transform("sum")
counts_all["percentage"] = (counts_all["cell_count"] / counts_all["total"]) * 100

pivot = counts_all.pivot_table(
    index=SAMPLE_COL, columns=CELLTYPE2_COL, values="percentage", fill_value=0,
)
out = f"{RESULTS_DIR}/celltype2_percentages_per_sample.csv"
pivot.to_csv(out)
print(f"  [OK] {out}")

# ANOVA by sex
meta = adata.obs[[SAMPLE_COL, AGE_COL, SEX_COL]].drop_duplicates().set_index(SAMPLE_COL)

anova_frames = []
for sex in adata.obs[SEX_COL].unique():
    sub = adata[adata.obs[SEX_COL] == sex].copy()
    ct_counts = (
        sub.obs.groupby([SAMPLE_COL, CELLTYPE2_COL])
        .size().reset_index(name="cell_count")
    )
    ct_counts["total"] = ct_counts.groupby(SAMPLE_COL)["cell_count"].transform("sum")
    ct_counts["percentage"] = (ct_counts["cell_count"] / ct_counts["total"]) * 100
    ct_counts[AGE_COL] = ct_counts[SAMPLE_COL].map(meta[AGE_COL])

    results = []
    for ct in ct_counts[CELLTYPE2_COL].unique():
        ct_sub = ct_counts[ct_counts[CELLTYPE2_COL] == ct]
        groups = [g["percentage"].values for _, g in ct_sub.groupby(AGE_COL) if len(g) >= 2]
        if len(groups) >= 2:
            stat, p = f_oneway(*groups)
        else:
            stat, p = np.nan, np.nan
        results.append({CELLTYPE2_COL: ct, "anova_stat": stat, "p_value": p})

    df_res = pd.DataFrame(results)
    df_res["adj_p"] = multipletests(df_res["p_value"].fillna(1), method="fdr_bh")[1]
    df_res["significance"] = df_res["adj_p"].apply(
        lambda p: "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else ""
    )
    df_res["sex"] = sex
    anova_frames.append(df_res)

anova_combined = pd.concat(anova_frames, ignore_index=True)
out = f"{RESULTS_DIR}/celltype2_age_ANOVA_by_sex.csv"
anova_combined.to_csv(out, index=False)
print(f"  [OK] {out}")
print(anova_combined.sort_values("adj_p").head(20).to_string(index=False))


# ==============================================================================
# 5. BOXPLOTS - CELLTYPE2 PROPORTIONS BY AGE
# ==============================================================================

print()
print("=" * 70)
print("ANALYSIS 5: Boxplots (celltype2 proportions by age)")
print("=" * 70)

for sex in ["female", "male"]:
    sub = adata[adata.obs[SEX_COL] == sex].copy()
    ct_counts = (
        sub.obs.groupby([SAMPLE_COL, CELLTYPE2_COL])
        .size().reset_index(name="cell_count")
    )
    ct_counts["total"] = ct_counts.groupby(SAMPLE_COL)["cell_count"].transform("sum")
    ct_counts["percentage"] = (ct_counts["cell_count"] / ct_counts["total"]) * 100
    ct_counts[AGE_COL] = ct_counts[SAMPLE_COL].map(meta[AGE_COL])
    ct_counts[SEX_COL] = sex

    anova_rows = []
    for ct in ct_counts[CELLTYPE2_COL].unique():
        ct_sub = ct_counts[ct_counts[CELLTYPE2_COL] == ct]
        groups = [g["percentage"].values for _, g in ct_sub.groupby(AGE_COL) if len(g) >= 2]
        if len(groups) >= 2:
            stat, p = f_oneway(*groups)
        else:
            stat, p = np.nan, np.nan
        anova_rows.append({CELLTYPE2_COL: ct, "anova_stat": stat, "p_value": p})
    anova_df = pd.DataFrame(anova_rows)
    anova_df["adj_p"] = multipletests(anova_df["p_value"].fillna(1), method="fdr_bh")[1]
    anova_df["significance"] = anova_df["adj_p"].apply(
        lambda p: "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else ""
    )

    ct_counts = ct_counts.merge(
        anova_df[[CELLTYPE2_COL, "significance"]], on=CELLTYPE2_COL, how="left",
    )
    ct_counts["significance"] = ct_counts["significance"].fillna("")

    present_ages = [a for a in AGE_ORDER if a in ct_counts[AGE_COL].dropna().unique()]
    ct_order = ct_counts[CELLTYPE2_COL].unique()

    fig, ax = plt.subplots(figsize=(14, 8))
    sns.boxplot(
        data=ct_counts, x=CELLTYPE2_COL, y="percentage", hue=AGE_COL,
        order=ct_order, hue_order=present_ages,
        palette={k: AGE_PALETTE[k] for k in present_ages},
        showcaps=True, fliersize=3, dodge=True, ax=ax,
    )
    ymax = ct_counts["percentage"].max()
    ax.set_ylim(0, ymax * 1.2)

    for i, ct in enumerate(ct_order):
        sub_ct = ct_counts[ct_counts[CELLTYPE2_COL] == ct]
        if sub_ct.empty:
            continue
        star = sub_ct["significance"].iloc[0]
        label, color = (star, "red") if star else ("ns", "darkblue")
        ax.text(i, sub_ct["percentage"].max() * 1.1, label,
                ha="center", va="bottom", fontsize=18, color=color)

    handles = [Patch(facecolor=AGE_PALETTE[a], edgecolor="black", label=a) for a in present_ages]
    ax.legend(
        handles=handles, title="Age Group",
        bbox_to_anchor=(1.01, 1), loc="upper left", frameon=False,
    )
    ax.set_title(f"{sex.capitalize()} Hepatocytes: Sub-cluster Proportions by Age", fontsize=16)
    ax.set_ylabel("Percentage of Cells per Sample", fontsize=12)
    ax.set_xlabel("Hepatocyte Sub-cluster", fontsize=12)
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()

    fname = f"{FIGURES_DIR}/{sex}_boxplot_celltype2_by_age.pdf"
    fig.savefig(fname, bbox_inches="tight")
    plt.close(fig)
    print(f"  [OK] {fname}")

    ct_counts.to_csv(f"{RESULTS_DIR}/{sex}_celltype2_percentages.csv", index=False)
    anova_df.to_csv(f"{RESULTS_DIR}/{sex}_celltype2_anova.csv", index=False)


# ==============================================================================
# 5b. STACKED BARPLOTS - CELLTYPE2 PROPORTIONS BY AGE
# ==============================================================================

print()
print("=" * 70)
print("ANALYSIS 5b: Stacked barplots (celltype2 by age)")
print("=" * 70)

data = adata.obs[[AGE_COL, SEX_COL, CELLTYPE2_COL]].copy()
data[AGE_COL] = pd.Categorical(data[AGE_COL], categories=AGE_ORDER, ordered=True)


def prepare_stacked_data(df, sex):
    subset = df[df[SEX_COL] == sex]
    counts = subset.groupby([AGE_COL, CELLTYPE2_COL]).size().reset_index(name="count")
    counts["percentage"] = (
        counts["count"] / counts.groupby(AGE_COL)["count"].transform("sum") * 100
    )
    counts = counts.sort_values(AGE_COL)
    stacked = counts.pivot(index=AGE_COL, columns=CELLTYPE2_COL, values="percentage").fillna(0)
    stacked = stacked[[ct for ct in CELLTYPE2_PALETTE if ct in stacked.columns]]
    return stacked


for sex in ["male", "female"]:
    stacked_data = prepare_stacked_data(data, sex)

    fig, ax = plt.subplots(figsize=(12, 8))
    stacked_data.plot(
        kind="bar", stacked=True,
        color=[CELLTYPE2_PALETTE[ct] for ct in stacked_data.columns],
        edgecolor="none", width=1.0, ax=ax,
    )
    ax.set_title(f"Distribution of Hepatocyte Sub-clusters by Age ({sex.capitalize()}s)",
                 fontsize=16)
    ax.set_ylabel("Percentage (%)", fontsize=14)
    ax.set_xlabel("Age Groups", fontsize=14)
    plt.xticks(rotation=45, ha="right", fontsize=12)
    ax.legend(title="Sub-cluster", bbox_to_anchor=(1.05, 1), loc="upper left", fontsize=10)
    plt.tight_layout()

    fname = f"{FIGURES_DIR}/stacked_bar_celltype2_by_age_{sex}.pdf"
    fig.savefig(fname, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  [OK] {fname}")

    stacked_data.to_csv(f"{RESULTS_DIR}/stacked_celltype2_proportions_{sex}.csv")


# ==============================================================================
# 6. GENE MODULE SCORING AND HEATMAP (GLOBAL)
# ==============================================================================

print()
print("=" * 70)
print("ANALYSIS 6: Gene module scoring (global, informs zonation)")
print("=" * 70)

score_name_global = "Zonation_score_ModuleScore"

present_genes = [g for g in ZONATION_GENE_LIST if g in adata.var_names]
use_raw = False
if not present_genes and adata.raw is not None:
    present_genes = [g for g in ZONATION_GENE_LIST if g in adata.raw.var_names]
    use_raw = bool(present_genes)
print(f"  Zonation genes found: {len(present_genes)} / {len(ZONATION_GENE_LIST)}")

sc.tl.score_genes(
    adata, gene_list=present_genes,
    score_name=score_name_global, use_raw=use_raw,
)

mean_df = (
    adata.obs.groupby(CELLTYPE2_COL)[score_name_global]
    .mean().to_frame(name="Mean")
)

fig, ax = plt.subplots(figsize=(6, max(4, 0.4 * len(mean_df))))
sns.heatmap(
    mean_df, annot=True, cmap="coolwarm", linewidths=0.5,
    cbar_kws={"label": "Avg Module Score"}, ax=ax,
)
ax.set_title(f"{score_name_global} (Mean) by Sub-cluster")
ax.set_xlabel("Statistic")
ax.set_ylabel("Hepatocyte Sub-cluster")
plt.tight_layout()

fname = f"{FIGURES_DIR}/zonation_module_score_heatmap.pdf"
fig.savefig(fname, dpi=300, bbox_inches="tight")
plt.close(fig)
print(f"  [OK] {fname}")


# ==============================================================================
# 7. HEPATOCYTE ZONATION: ASSIGNMENT, UMAP GRID, BOXPLOTS, STACKED BARS
# ==============================================================================

print()
print("=" * 70)
print("ANALYSIS 7: Hepatocyte zonation (assignment, UMAP, stats)")
print("=" * 70)

# 7a. Assign zonation labels
adata.obs[ZONATION_COL] = (
    adata.obs[CELLTYPE2_COL].replace(HEPATOCYTE_ZONATION_MAP).astype("category")
)
adata.obs[ZONATION_COL] = adata.obs[ZONATION_COL].cat.set_categories(
    ZONE_ORDER, ordered=True,
)
adata.uns[f"{ZONATION_COL}_colors"] = [ZONATION_PALETTE[z] for z in ZONE_ORDER]
print(adata.obs[ZONATION_COL].value_counts().to_string())

# 7b. Zonation UMAP grid (age x sex)
x_coords = adata.obsm["X_wnn"][:, 0]
y_coords = adata.obsm["X_wnn"][:, 1]
xlim_global = (x_coords.min(), x_coords.max())
ylim_global = (y_coords.min(), y_coords.max())

sorted_ages = [a for a in AGE_ORDER if a in adata.obs[AGE_COL].unique()]
sexes = adata.obs[SEX_COL].unique()
nrow, ncol_grid = len(sexes), len(sorted_ages)

fig, axes = plt.subplots(nrow, ncol_grid, figsize=(7 * ncol_grid, 5 * nrow))
if nrow == 1:
    axes = axes[np.newaxis, :]

for i, age_val in enumerate(sorted_ages):
    for j, sex_val in enumerate(sexes):
        ax = axes[j, i]
        mask = (adata.obs[AGE_COL] == age_val) & (adata.obs[SEX_COL] == sex_val)
        subset = adata[mask]
        legend_loc = "right margin" if (i == 0 and j == 0) else None

        sc.pl.embedding(
            subset, basis="X_wnn", color=ZONATION_COL,
            title=f"{age_val}, {sex_val}", palette=ZONATION_PALETTE,
            size=25, ax=ax, show=False, legend_loc=legend_loc,
        )
        ax.set_xlim(xlim_global)
        ax.set_ylim(ylim_global)

for ax in fig.axes:
    leg = ax.get_legend()
    if leg:
        for t in leg.get_texts():
            t.set_fontsize(18)
        leg.get_title().set_fontsize(20)

plt.tight_layout()
fname = f"{FIGURES_DIR}/umap_hepatocyte_zonation_grid.pdf"
fig.savefig(fname, bbox_inches="tight")
plt.close(fig)
print(f"  [OK] {fname}")

# 7c. Zonation proportions and ANOVA
zon_counts = (
    adata.obs.groupby([SAMPLE_COL, ZONATION_COL])
    .size().reset_index(name="cell_count")
)
zon_counts["total"] = zon_counts.groupby(SAMPLE_COL)["cell_count"].transform("sum")
zon_counts["percentage"] = (zon_counts["cell_count"] / zon_counts["total"]) * 100
zon_counts[AGE_COL] = zon_counts[SAMPLE_COL].map(meta[AGE_COL])
zon_counts[SEX_COL] = zon_counts[SAMPLE_COL].map(meta[SEX_COL])
zon_counts[AGE_COL] = pd.Categorical(zon_counts[AGE_COL], categories=AGE_ORDER, ordered=True)

zon_anova_rows = []
for sex in zon_counts[SEX_COL].dropna().unique():
    df_sex = zon_counts[zon_counts[SEX_COL] == sex]
    for zone in ZONE_ORDER:
        df_zone = df_sex[df_sex[ZONATION_COL] == zone]
        groups = [g["percentage"].values for _, g in df_zone.groupby(AGE_COL) if len(g) >= 2]
        if len(groups) >= 2:
            stat, p = f_oneway(*groups)
        else:
            stat, p = np.nan, np.nan
        zon_anova_rows.append({"sex": sex, "zonation": zone, "F_stat": stat, "p_value": p})

zon_anova_df = pd.DataFrame(zon_anova_rows)
zon_anova_df["adj_p"] = multipletests(zon_anova_df["p_value"].fillna(1), method="fdr_bh")[1]
zon_anova_df["significance"] = zon_anova_df["adj_p"].apply(
    lambda p: "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else ""
)
out = f"{RESULTS_DIR}/zonation_anova_by_sex.csv"
zon_anova_df.to_csv(out, index=False)
print(f"  [OK] {out}")

# 7d. Zonation boxplots
for sex in zon_counts[SEX_COL].dropna().unique():
    df_sex = zon_counts[zon_counts[SEX_COL] == sex]
    present_ages = [a for a in AGE_ORDER if a in df_sex[AGE_COL].dropna().unique()]

    fig, ax = plt.subplots(figsize=(12, 6))
    sns.boxplot(
        data=df_sex, x=ZONATION_COL, y="percentage", hue=AGE_COL,
        order=ZONE_ORDER, hue_order=present_ages,
        palette={a: AGE_PALETTE[a] for a in present_ages},
        showcaps=True, fliersize=3, dodge=True, ax=ax,
    )

    sub_anova = zon_anova_df[zon_anova_df["sex"] == sex]
    for i, zone in enumerate(ZONE_ORDER):
        row = sub_anova[sub_anova["zonation"] == zone]
        if row.empty:
            continue
        star = row["significance"].iloc[0]
        label, color = (star, "red") if star else ("ns", "darkblue")
        max_y = df_sex[df_sex[ZONATION_COL] == zone]["percentage"].max()
        ax.text(i, max_y + 1.5, label, ha="center", fontsize=14, color=color)

    ax.set_title(f"{sex.capitalize()} Samples: Zonation Percentage by Age")
    ax.set_ylabel("Percentage")
    ax.set_xlabel("Hepatocyte Zonation")
    ax.legend(title="Age", bbox_to_anchor=(1.01, 1), loc="upper left")
    plt.tight_layout()

    fname = f"{FIGURES_DIR}/{sex}_zonation_boxplot.pdf"
    fig.savefig(fname, bbox_inches="tight")
    plt.close(fig)
    print(f"  [OK] {fname}")

# 7e. Zonation stacked barplots
for sex in zon_counts[SEX_COL].dropna().unique():
    sub = zon_counts[zon_counts[SEX_COL] == sex]
    for value_col, ylabel, tag in [
        ("cell_count", "Cell Count", "stacked_counts"),
        ("percentage", "Percentage", "stacked_pct"),
    ]:
        pivot_zon = sub.pivot_table(
            index=[AGE_COL, SAMPLE_COL], columns=ZONATION_COL,
            values=value_col, fill_value=0,
        ).reset_index()
        pivot_zon = pivot_zon.sort_values(
            AGE_COL, key=lambda s: pd.Categorical(s, categories=AGE_ORDER, ordered=True),
        )

        fig, ax = plt.subplots(figsize=(12, 6))
        bottom = None
        x_labels = pivot_zon[AGE_COL].astype(str) + "_" + pivot_zon[SAMPLE_COL].astype(str)
        for zone in ZONE_ORDER:
            if zone not in pivot_zon.columns:
                continue
            ax.bar(
                x_labels, pivot_zon[zone], bottom=bottom, label=zone,
                color=ZONATION_PALETTE[zone], width=1.0, linewidth=0,
            )
            bottom = pivot_zon[zone] if bottom is None else bottom + pivot_zon[zone]

        ax.set_ylabel(ylabel)
        ax.set_title(f"{sex.capitalize()}: Zonation {ylabel}")
        ax.set_xticklabels(x_labels, rotation=45, ha="right", fontsize=7)
        ax.legend(title="Hepatocyte Zonation")
        plt.tight_layout()

        fname = f"{FIGURES_DIR}/{sex}_zonation_{tag}.pdf"
        fig.savefig(fname, bbox_inches="tight")
        plt.close(fig)
        print(f"  [OK] {fname}")


# ==============================================================================
# 8. PER-SUBSET ZONATION: CELL-LEVEL VIOLIN PLOTS + ANOVA
# ==============================================================================

print()
print("=" * 70)
print("ANALYSIS 8: Cell-level zonation violin plots + ANOVA")
print("=" * 70)

# 8a. Score per (age x sex) subset
score_col = "Zonation_score"
cell_parts = []

for sx in SEX_ORDER:
    for ag in AGE_ORDER:
        sub = adata[(adata.obs[AGE_COL] == ag) & (adata.obs[SEX_COL] == sx)].copy()
        if sub.n_obs == 0:
            continue

        present = [g for g in ZONATION_GENE_LIST if g in sub.var_names]
        use_raw_sub = False
        if not present and sub.raw is not None:
            present = [g for g in ZONATION_GENE_LIST if g in sub.raw.var_names]
            use_raw_sub = bool(present)
        if not present:
            print(f"  [SKIP] {sx} x {ag}: no zonation genes present")
            continue

        sc.tl.score_genes(sub, gene_list=present, score_name=score_col, use_raw=use_raw_sub)
        out_df = sub.obs[[score_col, AGE_COL, SEX_COL, CELLTYPE2_COL, SAMPLE_COL]].copy()
        cell_parts.append(out_df)

df_cells = pd.concat(cell_parts, ignore_index=True)
df_cells.rename(columns={score_col: "zonation_score"}, inplace=True)
df_cells[AGE_COL] = pd.Categorical(df_cells[AGE_COL], categories=AGE_ORDER, ordered=True)
df_cells[SEX_COL] = pd.Categorical(df_cells[SEX_COL], categories=SEX_ORDER, ordered=True)

present_ct = [ct for ct in HEP_ORDER if ct in df_cells[CELLTYPE2_COL].unique()]

# 8b. Sample-level means and ANOVA
df_sample = (
    df_cells
    .groupby([SAMPLE_COL, CELLTYPE2_COL, AGE_COL, SEX_COL])["zonation_score"]
    .mean().reset_index()
    .rename(columns={"zonation_score": "score_sample_mean"})
)

zonation_anova_rows = []
for sx in SEX_ORDER:
    sx_df = df_sample[df_sample[SEX_COL] == sx]
    for ct in present_ct:
        ct_df = sx_df[sx_df[CELLTYPE2_COL] == ct]
        groups = [g["score_sample_mean"].values
                  for _, g in ct_df.groupby(AGE_COL) if len(g) >= 2]
        if len(groups) >= 2:
            stat, p = f_oneway(*groups)
        else:
            stat, p = np.nan, np.nan
        zonation_anova_rows.append({
            "sex": sx, CELLTYPE2_COL: ct, "F_stat": stat, "p_value": p,
        })

zonation_anova_df = pd.DataFrame(zonation_anova_rows)
zonation_anova_df["adj_p"] = multipletests(
    zonation_anova_df["p_value"].fillna(1), method="fdr_bh",
)[1]
zonation_anova_df["significance"] = zonation_anova_df["adj_p"].apply(
    lambda p: "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "ns"
)

anova_csv_out = f"{RESULTS_DIR}/zonation_score_ANOVA_by_sex.csv"
zonation_anova_df.to_csv(anova_csv_out, index=False)
print(f"  [OK] {anova_csv_out}")

anova_xlsx_out = f"{RESULTS_DIR}/zonation_score_ANOVA_by_sex.xlsx"
zonation_anova_df.to_excel(anova_xlsx_out, index=False, sheet_name="ANOVA_results")
print(f"  [OK] {anova_xlsx_out}")

print(zonation_anova_df.to_string(index=False))

# 8c. Violin plot with ANOVA stars
xmin = np.nanmin(df_cells["zonation_score"])
xmax = np.nanmax(df_cells["zonation_score"])
xrng = xmax - xmin if np.isfinite(xmax - xmin) else 1.0
xmin_plot = xmin - 0.05 * xrng
xmax_plot = xmax + 0.08 * xrng

sns.set_style("whitegrid")
ncols_v = 7
nrows_block = int(np.ceil(len(present_ct) / ncols_v))

fig, axes = plt.subplots(
    2 * nrows_block, ncols_v,
    figsize=(ncols_v * 4.2, 2 * nrows_block * 3.8),
    squeeze=False,
)

for block_idx, (sex_label, sex_df) in enumerate([
    ("male", df_cells[df_cells[SEX_COL] == "male"]),
    ("female", df_cells[df_cells[SEX_COL] == "female"]),
]):
    row_offset = block_idx * nrows_block
    ax_block = axes[row_offset: row_offset + nrows_block, :]

    if sex_df.empty:
        for r in range(nrows_block):
            for c_ in range(ncols_v):
                ax_block[r, c_].axis("off")
        continue

    for i, ct in enumerate(present_ct):
        r, c = divmod(i, ncols_v)
        ax = ax_block[r, c]
        dat = sex_df[sex_df[CELLTYPE2_COL] == ct]
        if dat.empty:
            ax.axis("off")
            continue

        sns.violinplot(
            data=dat, y=AGE_COL, x="zonation_score",
            order=AGE_ORDER, palette=AGE_PALETTE,
            inner="box", cut=0, scale="width",
            linewidth=0.9, saturation=1.0, ax=ax,
        )
        ax.set_xlim(xmin_plot, xmax_plot)
        ax.set_title(ct, fontsize=11, fontweight="bold")
        ax.set_xlabel("Zonation score", fontsize=10)
        ax.set_ylabel("")
        ax.tick_params(axis="y", labelsize=9)
        sns.despine(ax=ax, offset=6, trim=True)

        row = zonation_anova_df[
            (zonation_anova_df["sex"] == sex_label)
            & (zonation_anova_df[CELLTYPE2_COL] == ct)
        ]
        if not row.empty:
            star = row["significance"].iloc[0]
            star_color = "red" if star != "ns" else "darkblue"
            ax.text(
                0.97, 0.95, star, transform=ax.transAxes, fontsize=18,
                fontweight="bold", color=star_color, ha="right", va="top",
            )

    for k in range(len(present_ct), nrows_block * ncols_v):
        r, c = divmod(k, ncols_v)
        ax_block[r, c].axis("off")

    ax_block[0, 0].set_title(
        f"{sex_label.capitalize()}", loc="left", fontsize=12, fontweight="bold",
    )

plt.suptitle(
    "Zonation (cell-level) by Age - Hep-01..Hep-07 (male top, female bottom)\n"
    "One-way ANOVA on sample means (score ~ age, FDR-corrected)",
    y=1.05, fontsize=16, fontweight="bold",
)
handles = [
    plt.Line2D(
        [0], [0], marker="s", markersize=18, color=col,
        linestyle="", markerfacecolor=col,
    )
    for col in AGE_PALETTE.values()
]
fig.legend(
    handles, list(AGE_PALETTE.keys()),
    title="Age group", title_fontsize=13,
    loc="upper center", bbox_to_anchor=(0.5, 1.015),
    ncol=len(AGE_PALETTE), fontsize=12, frameon=False,
)

plt.tight_layout()
plt.subplots_adjust(top=0.90)
fname = f"{FIGURES_DIR}/zonation_violin_cellLevel_Hep01to07.pdf"
fig.savefig(fname, dpi=450, bbox_inches="tight")
plt.close(fig)
print(f"  [OK] {fname}")


# ==============================================================================
# 9. PER-SUBSET ZONATION: SAMPLE-LEVEL BOXPLOTS + ANOVA STARS
# ==============================================================================

print()
print("=" * 70)
print("ANALYSIS 9: Sample-level zonation boxplots + ANOVA stars")
print("=" * 70)

out_csv = f"{RESULTS_DIR}/zonation_score_sample_means.csv"
df_sample.to_csv(out_csv, index=False)
print(f"  [OK] {out_csv}")

out_xlsx = f"{RESULTS_DIR}/zonation_score_sample_means.xlsx"
df_sample.to_excel(out_xlsx, index=False, sheet_name="sample_means")
print(f"  [OK] {out_xlsx}")

present_ct_s = [ct for ct in HEP_ORDER if ct in df_sample[CELLTYPE2_COL].unique()]

xmin_s = np.nanmin(df_sample["score_sample_mean"])
xmax_s = np.nanmax(df_sample["score_sample_mean"])
xrng_s = xmax_s - xmin_s if np.isfinite(xmax_s - xmin_s) else 1.0
xmin_plot_s = xmin_s - 0.05 * xrng_s
xmax_plot_s = xmax_s + 0.08 * xrng_s

ncols_b = 7
nrows_block_s = int(np.ceil(len(present_ct_s) / ncols_b))

fig, axes = plt.subplots(
    2 * nrows_block_s, ncols_b,
    figsize=(ncols_b * 4.2, 2 * nrows_block_s * 4.0),
    squeeze=False,
)

for block_idx, (sex_label, sex_df) in enumerate([
    ("male", df_sample[df_sample[SEX_COL] == "male"]),
    ("female", df_sample[df_sample[SEX_COL] == "female"]),
]):
    row_offset = block_idx * nrows_block_s
    ax_block = axes[row_offset: row_offset + nrows_block_s, :]

    if sex_df.empty:
        for r in range(nrows_block_s):
            for c_ in range(ncols_b):
                ax_block[r, c_].axis("off")
        continue

    for i, ct in enumerate(present_ct_s):
        r, c = divmod(i, ncols_b)
        ax = ax_block[r, c]
        dat = sex_df[sex_df[CELLTYPE2_COL] == ct]
        if dat.empty:
            ax.axis("off")
            continue

        sns.boxplot(
            data=dat, y=AGE_COL, x="score_sample_mean",
            order=AGE_ORDER, palette=AGE_PALETTE,
            linewidth=0.9, fliersize=1.8, ax=ax,
        )
        sns.stripplot(
            data=dat, y=AGE_COL, x="score_sample_mean",
            order=AGE_ORDER, color="black", alpha=0.85, size=2.8,
            jitter=0.25, dodge=False, ax=ax,
        )
        ax.set_xlim(xmin_plot_s, xmax_plot_s)
        ax.set_title(ct, fontsize=11, fontweight="bold")
        ax.set_xlabel("Zonation score", fontsize=10)
        ax.set_ylabel("")
        ax.tick_params(axis="y", labelsize=9)
        sns.despine(ax=ax, offset=6, trim=True)

        row = zonation_anova_df[
            (zonation_anova_df["sex"] == sex_label)
            & (zonation_anova_df[CELLTYPE2_COL] == ct)
        ]
        if not row.empty:
            star = row["significance"].iloc[0]
            star_color = "red" if star != "ns" else "darkblue"
            ax.text(
                0.97, 0.95, star, transform=ax.transAxes, fontsize=18,
                fontweight="bold", color=star_color, ha="right", va="top",
            )

    for k in range(len(present_ct_s), nrows_block_s * ncols_b):
        r, c = divmod(k, ncols_b)
        ax_block[r, c].axis("off")

    ax_block[0, 0].set_title(
        f"{sex_label.capitalize()}", loc="left", fontsize=12, fontweight="bold",
    )

plt.suptitle(
    "Zonation by Age - Hep-01..Hep-07 (male top, female bottom)\n"
    "One-way ANOVA (score ~ age, FDR-corrected)",
    y=1.05, fontsize=16, fontweight="bold",
)
handles = [
    plt.Line2D(
        [0], [0], marker="s", markersize=18, color=col,
        linestyle="", markerfacecolor=col,
    )
    for col in AGE_PALETTE.values()
]
fig.legend(
    handles, list(AGE_PALETTE.keys()),
    title="Age group", title_fontsize=13,
    loc="upper center", bbox_to_anchor=(0.5, 1.015),
    ncol=len(AGE_PALETTE), fontsize=12, frameon=False,
)

plt.tight_layout()
plt.subplots_adjust(top=0.90)
fname = f"{FIGURES_DIR}/zonation_boxplot_sampleLevel_Hep01to07.pdf"
fig.savefig(fname, dpi=450, bbox_inches="tight")
plt.close(fig)
print(f"  [OK] {fname}")


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
