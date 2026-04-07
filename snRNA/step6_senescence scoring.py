#!/usr/bin/env python3
# ==============================================================================
# Senescence Scoring in Hepatocyte Sub-Clusters Across Aging
# ==============================================================================
#
# Description:
#   Scores hepatocyte cells for senescence using four complementary gene sets:
#     1. Classical senescence markers (per-gene scoring)
#     2. SHGS  - Senescent Hepatocyte Gene Signature (Du et al. 2025 Nat Commun,
#                100-gene palbociclib / NRAS(G12V) mouse ortholog set)
#     3. SenMayo - 125-gene SASP/senescence panel (Saul et al. 2022 Nat Commun)
#     4. SenNet  - NIH SenNet consortium core biomarkers (2024)
#
#   For each module, produces side-by-side male/female heatmaps of mean module
#   score per sub-cluster x age, and a combined summary heatmap across all
#   three signatures.
#
# Input:
#   - Annotated AnnData (.h5ad) with celltype, celltype2 (sub-clusters),
#     sex, age labels
#
# Output:
#   Figures:
#     - panel_g_senescence_module_scores_combined.pdf
#         (single figure: 3 signatures x 2 sexes, each cell = celltype2 x age)
#     - classical_senescence_heatmap_male_female.pdf
#         (panel h: z-scored classical markers, male top / female bottom)
#
#   Tables:
#     - classical_senescence_zscore_{sex}.csv
#     - module_score_{signature}_{sex}.csv
#     - senescence_module_scores_summary.csv
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

# Global plot settings
mpl.rcParams["font.family"] = "Arial"
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42


# ==============================================================================
# CONFIGURATION - UPDATE THIS PATH
# ==============================================================================

DATA_PATH = "integrated_scvi.h5ad"
FIGURES_DIR = "figures"
RESULTS_DIR = "results"

os.makedirs(FIGURES_DIR, exist_ok=True)
os.makedirs(RESULTS_DIR, exist_ok=True)

# Column names
CELLTYPE_COL = "celltype"
CELLTYPE2_COL = "celltype2"
SEX_COL = "sex"
AGE_COL = "age"

# Ordering
AGE_ORDER = ["young", "mid_age", "old", "pre_geriatric", "geriatric"]
AGE_ALIASES = {"midage": "mid_age", "pregeriatric": "pre_geriatric"}
SEX_ORDER = ["male", "female"]
HEP_ORDER = [f"Hep-0{i}" for i in range(1, 8)]


# ==============================================================================
# SENESCENCE GENE SETS
# ==============================================================================

# Classical senescence markers (per-gene scoring)
CLASSICAL_SENESCENCE_GENES = [
    "Cdkn1a", "Cdkn2a", "Trp53", "Rb1", "Mdm2",
    "Il6", "Ccl2", "Cxcl1", "Serpine1", "Igfbp7",
]

# SHGS: Senescent Hepatocyte Gene Signature
# Du et al. 2025, Nat Commun - 100-gene overlap of in vitro palbociclib-treated
# Huh7 + in vivo NRAS(G12V) hepatocytes, mouse orthologs from the Zenodo
# senescence_signatures.rds
SHGS_GENES = [
    "Gm29685", "Gm19510", "Cxcl10", "4930546K05Rik", "Tox2",
    "Gm15418", "Glis3", "Rsad2", "Ly6m", "Osmr",
    "Gm48682", "Socs1", "Nid1", "Cib4", "Tgtp2",
    "Gbp2", "Hcar2", "Cd14", "Gm33782", "Anxa2",
    "Dusp5", "Gbp8", "Serpine1", "Gm30692", "Gm16685",
    "Tnfrsf12a", "Cdkn1a", "Hilpda", "Il7", "S1pr3",
    "2310034O05Rik", "Rasal1", "Tmem154", "Ralgds", "Bco1",
    "Relb", "Gldn", "Tgtp1", "Gm42765", "Cmpk2",
    "Cxcl9", "Cacng4", "Ier3", "Sorcs2", "Ifit3b",
    "Ubd", "Numbl", "Shc2", "Dchs1", "Gadd45b",
    "Nrg1", "Sult1e1", "Col4a1", "Gm29282", "Adcy1",
    "Clcf1", "Ifit2", "Gm48341", "Gbp9", "Fhl3",
    "Tgm2", "Gm37359", "Prr5l", "Tubb6", "Ly6a",
    "Gm17746", "Gm12185", "Cxcl13", "Oasl2", "Slc7a1",
    "2010003K11Rik", "Zfp462", "Olfr56", "Gm39348", "H2-Q1",
    "Igdcc4", "4933401D09Rik", "Casp12", "Srxn1", "Gucy2c",
    "Igfbp1", "Gldnos", "Scn2a", "Gpc6", "Tlr2",
    "Mt2", "Gipc2", "Rragd", "Mab21l3", "5330417C22Rik",
    "Casp4", "Arid5a", "Rims4", "Coq8b", "9930111J21Rik2",
    "Gm45774", "Ifit3", "BC016579", "Atp6v0d2", "Ccdc120",
    "Lpl", "Cacna1b", "Ly6d", "Spdef", "Ier5",
    "Mt1", "Rnf144a", "Alpk1", "Gm11454", "Slc4a9",
    "Slc1a4", "9030624G23Rik", "Olfr16", "Akr1b7", "Ddr1",
    "2700038G22Rik", "Gm15567", "Lepr", "Inhbb", "Lcn2",
    "Gm26788", "Syt12", "Dusp8", "Slc37a1", "Gm33699",
    "Asb11", "Tnfaip3", "Ifit1", "Epha4", "Serpina3n",
    "Gm36037", "Cerk", "Rp1", "Hbegf", "Gm50383",
    "Spats2l", "Gbp6", "Scn8a", "Npas2", "Klf6",
    "Ifi44", "Sox9", "Ffar4", "Gm20125", "Tagln2",
    "Trim6", "Gm26714", "Cd274", "Ccdc180", "Dnah1",
    "1200007C13Rik", "Sytl5", "Crybb3", "Ets2", "Chka",
    "Pde6h", "Irak3", "Gm41609", "Il1rn", "Prag1",
    "Arhgef2", "Tnfrsf23", "Itpkc", "Nlrp9c",
]

# SenMayo: 125-gene SASP/senescence panel
# Saul et al. 2022, Nat Commun - mouse gene names
SENMAYO_GENES = [
    "Acvr1b", "Ang", "Angpt1", "Angptl4", "Areg", "Axl", "Bex3",
    "Bmp2", "Bmp6", "C3", "Ccl1", "Ccl2", "Ccl20", "Ccl24",
    "Ccl26", "Ccl3", "Ccl4", "Ccl5", "Ccl7", "Ccl8", "Cd55",
    "Cd9", "Csf1", "Csf2", "Csf2rb", "Ctnnb1", "Ctsb", "Cxcl1",
    "Cxcl10", "Cxcl12", "Cxcl16", "Cxcl2", "Cxcl3", "Cxcr2",
    "Dkk1", "Edn1", "Egf", "Egfr", "Ereg", "Esm1", "Ets2",
    "Fas", "Fgf1", "Fgf2", "Fgf7", "Gdf15", "Gem", "Gmfg",
    "Hgf", "Hmgb1", "Icam1", "Icam3", "Igf1", "Igfbp1", "Igfbp2",
    "Igfbp3", "Igfbp4", "Igfbp5", "Igfbp6", "Igfbp7", "Il10",
    "Il13", "Il15", "Il18", "Il1a", "Il1b", "Il2", "Il6", "Il6st",
    "Il7", "Inha", "Iqgap2", "Itga2", "Itpka", "Jun", "Kitl",
    "Lcp1", "Mif", "Mmp10", "Mmp12", "Mmp13", "Mmp14", "Mmp2",
    "Mmp3", "Mmp9", "Nap1l4", "Nrg1", "Pappa", "Pecam1", "Pgf",
    "Pigf", "Plat", "Plau", "Plaur", "Ptbp1", "Ptger2", "Ptges",
    "Rps6ka5", "Scamp4", "Selplg", "Sema3f", "Serpinb4", "Serpine1",
    "Serpine2", "Spp1", "Spx", "Timp2", "Tnf", "Tnfrsf11b",
    "Tnfrsf1a", "Tnfrsf1b", "Tubgcp2", "Vegfa", "Vegfc", "Vgf",
    "Wnt16", "Wnt2",
]

# SenNet: Core senescence biomarkers (NIH SenNet consortium 2024)
SENNET_GENES = [
    # Cell Cycle Arrest
    "Cdkn2a", "Cdkn1a", "Trp53", "Rb1",
    # DNA Damage Response
    "H2afx", "Trp53bp1", "Atm", "Atr",
    # Nuclear Changes
    "Lmnb1",
    # Anti-Apoptotic (pro-survival)
    "Bcl2", "Bcl2l1", "Bcl2l2", "Serpinb2",
    # SASP Core
    "Il1a", "Il1b", "Il6", "Cxcl1", "Ccl2", "Mmp3", "Mmp9",
    # Surface Markers
    "Dpp4", "Cd36", "Icam1", "Plaur", "Tnfrsf10d",
    # Lysosomal
    "Glb1",
]


# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================

def score_module(adata, gene_list, score_name, title_label):
    """
    Score cells with sc.tl.score_genes and return a dict of per-sex pivot
    tables (celltype2 rows x age columns, mean score).
    """
    present = [g for g in gene_list if g in adata.var_names]
    use_raw = False
    if not present and adata.raw is not None:
        present = [g for g in gene_list if g in adata.raw.var_names]
        use_raw = bool(present)
    print(f"  {title_label}: {len(present)}/{len(gene_list)} genes found")

    if not present:
        print(f"  [SKIP] No genes found for {title_label}")
        return None

    sc.tl.score_genes(adata, gene_list=present, score_name=score_name, use_raw=use_raw)

    dfs = {}
    for sex in SEX_ORDER:
        if sex not in adata.obs[SEX_COL].unique():
            continue
        sub = adata.obs[adata.obs[SEX_COL] == sex]
        pivot = sub.pivot_table(
            index=CELLTYPE2_COL, columns=AGE_COL,
            values=score_name, aggfunc="mean",
        )
        row_order = [ct for ct in HEP_ORDER if ct in pivot.index]
        col_order = [a for a in AGE_ORDER if a in pivot.columns]
        pivot = pivot.loc[row_order, col_order]
        dfs[sex] = pivot

        # Save underlying data
        out = f"{RESULTS_DIR}/module_score_{title_label}_{sex}.csv"
        pivot.to_csv(out)
        print(f"  [OK] {out}")

    return dfs


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


# ==============================================================================
# 1. CLASSICAL SENESCENCE MARKERS - PER-GENE SCORING + Z-SCORED HEATMAP
# ==============================================================================

print()
print("=" * 70)
print("ANALYSIS 1: Classical senescence markers (per-gene scoring)")
print("=" * 70)

# Score each gene individually
for gene in CLASSICAL_SENESCENCE_GENES:
    score_name = f"{gene}_score"
    if gene in adata.var_names:
        sc.tl.score_genes(adata, gene_list=[gene], score_name=score_name)
    elif adata.raw is not None and gene in adata.raw.var_names:
        sc.tl.score_genes(adata, gene_list=[gene], score_name=score_name, use_raw=True)
    else:
        print(f"  [SKIP] {gene} not found")

score_cols = [
    f"{g}_score" for g in CLASSICAL_SENESCENCE_GENES
    if f"{g}_score" in adata.obs.columns
]
scored_genes = [col.replace("_score", "") for col in score_cols]

# Per-sex z-scored (gene x celltype2) matrices
zscore_per_sex = {}
for sex in SEX_ORDER:
    if sex not in adata.obs[SEX_COL].unique():
        continue
    df_sex = adata.obs[adata.obs[SEX_COL] == sex]

    # Mean score per gene across celltype2
    df = df_sex.groupby(CELLTYPE2_COL)[score_cols].mean().T  # genes x celltypes

    # Reorder celltype2 columns by HEP_ORDER
    col_order = [ct for ct in HEP_ORDER if ct in df.columns]
    df = df[col_order]

    # Row-wise z-score (per gene, across cell types)
    df_zscore = df.sub(df.mean(axis=1), axis=0).div(df.std(axis=1), axis=0)
    df_zscore["Gene"] = df_zscore.index.str.replace("_score", "", regex=False)
    df_zscore = df_zscore.set_index("Gene").loc[scored_genes]

    zscore_per_sex[sex] = df_zscore

    # Save underlying data
    out_csv = f"{RESULTS_DIR}/classical_senescence_zscore_{sex}.csv"
    df_zscore.to_csv(out_csv)
    print(f"  [OK] {out_csv}")

# Combined figure - male on top, female on bottom (matching panel h)
panel_sexes = [s for s in SEX_ORDER if s in zscore_per_sex]

if panel_sexes:
    # Shared color range across both sexes
    all_vals = pd.concat(
        [zscore_per_sex[s].stack() for s in panel_sexes], ignore_index=True,
    ).dropna()
    vmax = float(max(abs(all_vals.min()), abs(all_vals.max())))
    vmin = -vmax

    n_rows = len(scored_genes)
    fig, axes = plt.subplots(
        len(panel_sexes), 1,
        figsize=(8, max(6, 0.45 * n_rows * len(panel_sexes))),
        sharex=True,
    )
    if len(panel_sexes) == 1:
        axes = [axes]

    for ax, sex in zip(axes, panel_sexes):
        df_zscore = zscore_per_sex[sex]
        sns.heatmap(
            df_zscore, ax=ax,
            cmap="coolwarm", center=0,
            vmin=vmin, vmax=vmax,
            linewidths=0.05,
            xticklabels=True, yticklabels=True,
            cbar=True,
            cbar_kws={"label": "Z-score"},
        )
        ax.set_title(sex.capitalize(), fontsize=14, fontweight="bold", loc="left")
        ax.set_ylabel("")
        if sex == panel_sexes[-1]:
            ax.set_xlabel("classical senescence mRNAs", fontsize=12)
        else:
            ax.set_xlabel("")
        for label in ax.get_xticklabels():
            label.set_rotation(0)
            label.set_horizontalalignment("center")

    plt.tight_layout()

    fname = f"{FIGURES_DIR}/classical_senescence_heatmap_male_female.pdf"
    fig.savefig(fname, dpi=300, bbox_inches="tight")
    fig.savefig(fname.replace(".pdf", ".png"), dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  [OK] {fname}")


# ==============================================================================
# 2. PANEL (g): COMBINED MODULE SCORE HEATMAP
#    3 rows (SenMayo, SenNet, SHGS) x 2 columns (male, female)
#    Each cell = celltype2 (rows) x age (columns), mean module score
# ==============================================================================

print()
print("=" * 70)
print("ANALYSIS 2: Combined senescence module score heatmap (panel g)")
print("=" * 70)

signature_specs = [
    ("SenMayo", SENMAYO_GENES, "SenMayo_score"),
    ("SenNet",  SENNET_GENES,  "SenNet_score"),
    ("SHGS",    SHGS_GENES,    "SHGS_score"),
]

signature_dfs = {}
for label, gene_list, score_col in signature_specs:
    print(f"\n  Scoring {label}...")
    dfs = score_module(adata, gene_list, score_col, label)
    if dfs is not None:
        signature_dfs[label] = dfs

# Build the combined 3-row x 2-col figure
panel_sexes = [s for s in SEX_ORDER if s in adata.obs[SEX_COL].unique()]
n_sigs = len(signature_dfs)

if n_sigs and panel_sexes:
    fig, axes = plt.subplots(
        n_sigs, len(panel_sexes),
        figsize=(6 * len(panel_sexes), 3 * n_sigs),
        squeeze=False,
    )

    for row_idx, (label, _, _) in enumerate(
        [s for s in signature_specs if s[0] in signature_dfs]
    ):
        dfs = signature_dfs[label]

        # Independent colorbar range per signature row
        row_vals = pd.concat(
            [dfs[s].stack() for s in panel_sexes if s in dfs],
            ignore_index=True,
        ).dropna()
        vmin = float(row_vals.min())
        vmax = float(row_vals.max())

        for col_idx, sex in enumerate(panel_sexes):
            ax = axes[row_idx, col_idx]
            if sex not in dfs or dfs[sex].empty:
                ax.axis("off")
                continue

            df_plot = dfs[sex]
            show_cbar = (col_idx == len(panel_sexes) - 1)

            sns.heatmap(
                df_plot, ax=ax, annot=True, fmt=".3f",
                cmap="coolwarm", vmin=vmin, vmax=vmax,
                linewidths=0.05,
                cbar=show_cbar,
                cbar_kws={"label": "mean module score"} if show_cbar else None,
            )

            # Column titles (male/female) only on the top row
            if row_idx == 0:
                ax.set_title(sex, fontsize=13, fontweight="bold")
            else:
                ax.set_title("")

            # X-axis label only on bottom row
            if row_idx == n_sigs - 1:
                ax.set_xlabel("")
                for lab in ax.get_xticklabels():
                    lab.set_rotation(0)
                    lab.set_horizontalalignment("center")
            else:
                ax.set_xlabel("")
                ax.set_xticklabels([])

            # Y-axis label only on left column; signature label on right col
            if col_idx == 0:
                ax.set_ylabel("")
            else:
                ax.set_ylabel("")
                # Put signature label on the right side
                ax.text(
                    1.25, 0.5, label,
                    transform=ax.transAxes,
                    fontsize=13, fontweight="bold",
                    ha="left", va="center", rotation=270,
                )

    plt.tight_layout()
    fname = f"{FIGURES_DIR}/panel_g_senescence_module_scores_combined.pdf"
    fig.savefig(fname, dpi=300, bbox_inches="tight")
    fig.savefig(fname.replace(".pdf", ".png"), dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"\n  [OK] {fname}")


# ==============================================================================
# 3. COMBINED SENESCENCE SUMMARY TABLE
# ==============================================================================

print()
print("=" * 70)
print("ANALYSIS 3: Combined senescence summary table")
print("=" * 70)

summary_rows = []
for sig_name, score_col in [
    ("SHGS", "SHGS_score"),
    ("SenMayo", "SenMayo_score"),
    ("SenNet", "SenNet_score"),
]:
    if score_col not in adata.obs.columns:
        continue
    for sex in SEX_ORDER:
        if sex not in adata.obs[SEX_COL].unique():
            continue
        sub = adata.obs[adata.obs[SEX_COL] == sex]
        means = sub.groupby(CELLTYPE2_COL)[score_col].mean()
        for ct, val in means.items():
            summary_rows.append({
                "signature": sig_name,
                "sex": sex,
                CELLTYPE2_COL: ct,
                "mean_score": val,
            })

if summary_rows:
    df_summary = pd.DataFrame(summary_rows)
    out = f"{RESULTS_DIR}/senescence_module_scores_summary.csv"
    df_summary.to_csv(out, index=False)
    print(f"  [OK] {out}")


# ==============================================================================
# SUMMARY
# ==============================================================================

print()
print("=" * 70)
print("ALL SENESCENCE ANALYSES COMPLETE")
print("=" * 70)
print(f"  Figures: {FIGURES_DIR}/")
print(f"  Tables:  {RESULTS_DIR}/")
print("=" * 70)
