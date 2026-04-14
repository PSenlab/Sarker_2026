#!/usr/bin/env python3
# ==============================================================================
# Ascl1 GMM High/Low Classification in Female Hepatocytes
# ==============================================================================
#
# Description:
#   Classifies Ascl1+ female hepatocytes into Ascl1-High vs Ascl1-Low
#   populations using a 2-component Gaussian Mixture Model (GMM) fit to
#   the non-zero expression distribution, validates bimodality via AIC/BIC
#   model selection, then performs downstream differential expression
#   (Wilcoxon) and Reactome GSEA between the two populations.
#
# Pipeline:
#   1. Load data and subset to female hepatocytes
#   2. GMM classification (Ascl1-High vs Ascl1-Low) among Ascl1+ cells
#   3. GMM High/Low count dot plots for selected celltype2 sub-clusters
#   4. Wilcoxon DE: Ascl1-High vs Ascl1-Low
#   5. GSEA Reactome Pathways 2024
#   6. GSEA dot plots (top 30 positive NES, top 30 negative NES)
#
# Input:
#   - integrated_scvi.h5ad  (annotated AnnData with celltype, celltype2,
#                            sex, age labels)
#
# Output:
#   - Figures in figures/
#   - Tables in results/
#
#
# ==============================================================================

import re
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
import gseapy as gp
from sklearn.mixture import GaussianMixture

warnings.filterwarnings("ignore")

# Global plot settings
mpl.rcParams["font.family"] = "Arial"
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42


# ==============================================================================
# CONFIGURATION - UPDATE THIS PATH
# ==============================================================================

H5AD_PATH = "integrated_scvi.h5ad"

FIG_DIR = Path("figures")
RES_DIR = Path("results")
FIG_DIR.mkdir(exist_ok=True)
RES_DIR.mkdir(exist_ok=True)

# Column names
CELLTYPE_COL  = "celltype"
CELLTYPE2_COL = "celltype2"
SEX_COL       = "sex"
AGE_COL       = "age"

GENE = "Ascl1"

# Selected hepatocyte sub-clusters (panel d celltype2 subtypes with
# sufficient Ascl1+ cells)
KEEP_CT2 = ["Hep-07", "Hep-06", "Hep-04", "Hep-03", "Hep-02"]

# Age order
AGE_ORDER     = ["young", "mid_age", "old", "pre_geriatric", "geriatric"]
AGE_ORDER_REV = list(reversed(AGE_ORDER))
AGE_ALIASES   = {"midage": "mid_age", "pregeriatric": "pre_geriatric"}

# DE thresholds
PVAL_THRESHOLD       = 0.05
LOGFC_STRONG_CUTOFF  = 0.5

# GSEA filters
NOM_P_CUTOFF  = 0.05
TOP_N_GSEA    = 30


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


def plot_gsea_dotplot(df, direction, outfile):
    """
    Generic dot plot for the top N Reactome pathways of one NES direction.
    df must contain columns: term_clean, NES, Gene %, NOM p-val.
    direction: "positive" or "negative" (drives title wording only).
    """
    if df.empty:
        print(f"  [SKIP] No pathways to plot for {direction} NES")
        return

    df = df.copy()
    df["pval_size"] = -np.log10(df["NOM p-val"].replace(0, 1e-10))

    min_s, max_s = 60, 400
    pmin, pmax = df["pval_size"].min(), df["pval_size"].max()
    if pmax > pmin:
        df["size_scaled"] = (
            min_s + (df["pval_size"] - pmin) / (pmax - pmin) * (max_s - min_s)
        )
    else:
        df["size_scaled"] = (min_s + max_s) / 2

    # Sort so most positive NES is at top (or least negative at top for negative)
    df = df.sort_values("NES", ascending=False)
    df["term_clean"] = pd.Categorical(
        df["term_clean"], categories=df["term_clean"], ordered=True,
    )
    nes_min, nes_max = df["NES"].min(), df["NES"].max()

    label_map = {
        "positive": "Enriched in Ascl1-High (NES > 0)",
        "negative": "Enriched in Ascl1-Low (NES < 0)",
    }

    plt.figure(figsize=(14, 7))
    ax = sns.scatterplot(
        data=df, x="NES", y="term_clean",
        size="pval_size", hue="Gene %",
        palette="RdYlBu_r", sizes=(min_s, max_s),
        edgecolor="black", linewidth=0.6,
    )
    ax.set_xlim(nes_min - 0.05, nes_max + 0.05)
    plt.axvline(0, color="gray", lw=1)
    plt.xlabel("Normalized enrichment score (NES)", fontsize=12, weight="bold")
    plt.ylabel("Reactome pathway", fontsize=12, weight="bold")
    plt.title(
        f"Top {TOP_N_GSEA} Reactome pathways - {label_map[direction]}\n"
        f"Female hepatocyte CT2 subset (GMM High vs Low)   "
        f"NES range: {nes_min:.2f} to {nes_max:.2f}",
        fontsize=14, weight="bold",
    )

    norm = plt.Normalize(df["Gene %"].min(), df["Gene %"].max())
    sm = plt.cm.ScalarMappable(cmap="RdYlBu_r", norm=norm)
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax)
    cbar.set_label("Gene %", fontsize=11, weight="bold")

    ax.legend(
        title="-log10(NOM p-val)",
        bbox_to_anchor=(1.45, 1), loc="upper left",
        frameon=True, fontsize=9, title_fontsize=10,
    )

    plt.tight_layout()
    plt.savefig(outfile, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"  [OK] {outfile}")


# ==============================================================================
# 1. LOAD DATA
# ==============================================================================

banner("STEP 1: Load data and subset to female hepatocytes")

print(f"  Loading {H5AD_PATH}")
adata = sc.read_h5ad(H5AD_PATH)
print(f"  [OK] {adata.n_obs:,} cells x {adata.n_vars:,} genes")

adata = harmonize_age(adata)

adata_fh = adata[
    (adata.obs[SEX_COL].astype(str) == "female")
    & (adata.obs[CELLTYPE_COL].astype(str) == "Hepatocyte")
].copy()

if adata_fh.n_obs == 0:
    raise ValueError("No female hepatocytes found")
if GENE not in adata_fh.var_names:
    raise ValueError(f"{GENE} not found in var_names")

print(f"  Female hepatocytes: {adata_fh.n_obs:,}")

adata_fh.obs["Ascl1_expr"] = extract_gene_expression(adata_fh, GENE)
adata_fh.obs["Ascl1_positive"] = adata_fh.obs["Ascl1_expr"] > 0

adata_pos = adata_fh[adata_fh.obs["Ascl1_positive"]].copy()
print(f"  Ascl1+ female hepatocytes: {adata_pos.n_obs:,} "
      f"({100 * adata_pos.n_obs / adata_fh.n_obs:.1f}%)")


# ==============================================================================
# 2. GMM CLASSIFICATION (Ascl1-High vs Ascl1-Low)
# ==============================================================================

banner("STEP 2: GMM classification - Ascl1-High vs Ascl1-Low")

# Fit 2-component GMM to non-zero Ascl1 expression
expr_vals = adata_pos.obs["Ascl1_expr"].values.reshape(-1, 1)

gmm = GaussianMixture(n_components=2, random_state=0).fit(expr_vals)
clusters = gmm.predict(expr_vals)
means = gmm.means_.flatten()

# The cluster with the higher mean is the "High" group
high_cluster = int(np.argmax(means))
adata_pos.obs["Ascl1_group"] = np.where(clusters == high_cluster, "High", "Low")

# Posterior probabilities
probs = gmm.predict_proba(expr_vals)
adata_pos.obs["prob_high"] = probs[:, high_cluster]
adata_pos.obs["prob_low"]  = probs[:, 1 - high_cluster]

# Model selection: compare 1/2/3-component models
gmm_1 = GaussianMixture(n_components=1, random_state=0).fit(expr_vals)
gmm_3 = GaussianMixture(n_components=3, random_state=0).fit(expr_vals)

low_count  = int((adata_pos.obs["Ascl1_group"] == "Low").sum())
high_count = int((adata_pos.obs["Ascl1_group"] == "High").sum())

print(f"\n  GMM classification:")
print(f"    Low:  {low_count:,} cells ({100 * low_count / adata_pos.n_obs:.1f}%)")
print(f"    High: {high_count:,} cells ({100 * high_count / adata_pos.n_obs:.1f}%)")

print(f"\n  Cluster parameters:")
print(f"    Low mean:  {means[1 - high_cluster]:.3f}")
print(f"    High mean: {means[high_cluster]:.3f}")
print(f"    Fold diff: {means[high_cluster] / means[1 - high_cluster]:.2f}x")

print(f"\n  Expression statistics by group:")
for group in ["Low", "High"]:
    sub = adata_pos.obs[adata_pos.obs["Ascl1_group"] == group]["Ascl1_expr"]
    print(f"    {group}: mean={sub.mean():.3f}, median={sub.median():.3f}, "
          f"std={sub.std():.3f}, range=[{sub.min():.3f}, {sub.max():.3f}]")

print(f"\n  Model selection (AIC / BIC):")
print(f"    {'Model':<15} {'AIC':>12} {'BIC':>12}")
print(f"    {'-' * 15} {'-' * 12} {'-' * 12}")
print(f"    {'1-component':<15} {gmm_1.aic(expr_vals):>12.1f} {gmm_1.bic(expr_vals):>12.1f}")
print(f"    {'2-component':<15} {gmm.aic(expr_vals):>12.1f} {gmm.bic(expr_vals):>12.1f}")
print(f"    {'3-component':<15} {gmm_3.aic(expr_vals):>12.1f} {gmm_3.bic(expr_vals):>12.1f}")

delta_aic = gmm_1.aic(expr_vals) - gmm.aic(expr_vals)
delta_bic = gmm_1.bic(expr_vals) - gmm.bic(expr_vals)
print(f"\n    dAIC (1 vs 2): {delta_aic:.1f} (2-component better)")
print(f"    dBIC (1 vs 2): {delta_bic:.1f} (2-component better)")

high_conf = int(((adata_pos.obs["prob_high"] > 0.9)
                 | (adata_pos.obs["prob_low"] > 0.9)).sum())
print(f"\n  High-confidence (prob > 0.9): {high_conf:,} cells "
      f"({100 * high_conf / adata_pos.n_obs:.1f}%)")

print(f"\n  GMM groups by age:")
print(pd.crosstab(adata_pos.obs[AGE_COL],
                  adata_pos.obs["Ascl1_group"], margins=True).to_string())

print(f"\n  GMM groups by celltype2:")
print(pd.crosstab(adata_pos.obs[CELLTYPE2_COL],
                  adata_pos.obs["Ascl1_group"], margins=True).to_string())

# Save results
out_df = adata_pos.obs[["Ascl1_expr", "Ascl1_group", "prob_low", "prob_high"]]
out_df.to_csv(RES_DIR / "Ascl1_GMM_HighLow_female_hepatocytes.csv", index=True)
print(f"\n  [OK] {RES_DIR}/Ascl1_GMM_HighLow_female_hepatocytes.csv")

age_gmm_summary = adata_pos.obs.groupby(
    [AGE_COL, "Ascl1_group"], observed=True
).agg(
    count=("Ascl1_expr", "size"),
    mean_expr=("Ascl1_expr", "mean"),
).reset_index()
age_gmm_summary.to_csv(RES_DIR / "Ascl1_GMM_summary_by_age.csv", index=False)

ct2_gmm_summary = adata_pos.obs.groupby(
    [CELLTYPE2_COL, "Ascl1_group"], observed=True
).agg(
    count=("Ascl1_expr", "size"),
    mean_expr=("Ascl1_expr", "mean"),
).reset_index()
ct2_gmm_summary.to_csv(RES_DIR / "Ascl1_GMM_summary_by_celltype2.csv", index=False)
print(f"  [OK] {RES_DIR}/Ascl1_GMM_summary_by_age.csv")
print(f"  [OK] {RES_DIR}/Ascl1_GMM_summary_by_celltype2.csv")


# ==============================================================================
# 3. GMM HIGH/LOW COUNT DOT PLOTS (selected celltype2 subtypes)
# ==============================================================================

banner("STEP 3: GMM High/Low count dot plots (selected celltype2)")

adata_hep_f_ct2_pos = adata_pos[
    adata_pos.obs[CELLTYPE2_COL].isin(KEEP_CT2)
].copy()

print(f"  Ascl1+ cells in selected subtypes: {adata_hep_f_ct2_pos.n_obs:,}")
print(f"\n  GMM group distribution:")
print(adata_hep_f_ct2_pos.obs["Ascl1_group"].value_counts().to_string())

# Ordered categories
age_order_rev_present = [
    a for a in AGE_ORDER_REV if a in adata_hep_f_ct2_pos.obs[AGE_COL].unique()
]
adata_hep_f_ct2_pos.obs[AGE_COL] = pd.Categorical(
    adata_hep_f_ct2_pos.obs[AGE_COL].astype(str),
    categories=age_order_rev_present, ordered=True,
)
adata_hep_f_ct2_pos.obs[CELLTYPE2_COL] = pd.Categorical(
    adata_hep_f_ct2_pos.obs[CELLTYPE2_COL].astype(str),
    categories=KEEP_CT2, ordered=True,
)

# Count tables
age_ct = (
    adata_hep_f_ct2_pos.obs
    .groupby([AGE_COL, "Ascl1_group"], observed=True)
    .size().reset_index(name="cell_count")
)
ct2_ct = (
    adata_hep_f_ct2_pos.obs
    .groupby([CELLTYPE2_COL, "Ascl1_group"], observed=True)
    .size().reset_index(name="cell_count")
)

# Summary
total = adata_hep_f_ct2_pos.n_obs
n_high = int((adata_hep_f_ct2_pos.obs["Ascl1_group"] == "High").sum())
n_low  = int((adata_hep_f_ct2_pos.obs["Ascl1_group"] == "Low").sum())
print(f"\n  Overall:")
print(f"    Total Ascl1+: {total:,}")
print(f"    High: {n_high:,} ({100 * n_high / total:.1f}%)")
print(f"    Low:  {n_low:,} ({100 * n_low / total:.1f}%)")

age_summary = adata_hep_f_ct2_pos.obs.groupby(AGE_COL, observed=True).agg(
    total=("Ascl1_group", "size"),
    high=("Ascl1_group", lambda x: int((x == "High").sum())),
    low=("Ascl1_group", lambda x: int((x == "Low").sum())),
).reset_index()
age_summary["pct_high"] = 100 * age_summary["high"] / age_summary["total"]
age_summary["pct_low"]  = 100 * age_summary["low"]  / age_summary["total"]
print(f"\n  GMM groups by age:\n{age_summary.to_string(index=False)}")

ct2_summary = adata_hep_f_ct2_pos.obs.groupby(CELLTYPE2_COL, observed=True).agg(
    total=("Ascl1_group", "size"),
    high=("Ascl1_group", lambda x: int((x == "High").sum())),
    low=("Ascl1_group", lambda x: int((x == "Low").sum())),
).reset_index()
ct2_summary["pct_high"] = 100 * ct2_summary["high"] / ct2_summary["total"]
ct2_summary["pct_low"]  = 100 * ct2_summary["low"]  / ct2_summary["total"]
ct2_summary = ct2_summary.sort_values("pct_high", ascending=False)
print(f"\n  GMM groups by celltype2:\n{ct2_summary.to_string(index=False)}")

age_summary.to_csv(RES_DIR / "Ascl1_GMM_dotplot_summary_by_age.csv", index=False)
ct2_summary.to_csv(RES_DIR / "Ascl1_GMM_dotplot_summary_by_celltype2.csv", index=False)

# Dodge Low/High positions on x-axis
dodge_map = {"Low": -0.15, "High": 0.15}
group_x_map = {"Low": 0, "High": 1}

age_ct["x_pos"] = (
    age_ct["Ascl1_group"].map(group_x_map)
    + age_ct["Ascl1_group"].map(dodge_map)
)
ct2_ct["x_pos"] = (
    ct2_ct["Ascl1_group"].map(group_x_map)
    + ct2_ct["Ascl1_group"].map(dodge_map)
)

fig, axes = plt.subplots(1, 2, figsize=(10, 5))

sns.scatterplot(
    data=age_ct, x="x_pos", y=AGE_COL,
    size="cell_count", hue="cell_count",
    palette="Reds", sizes=(60, 550),
    edgecolor="black", linewidth=0.6, ax=axes[0],
)
axes[0].set_title("Ascl1 GMM High/Low - Age (Female Hep)", weight="bold")
axes[0].set_xticks([0, 1])
axes[0].set_xticklabels(["Low", "High"])
axes[0].set_xlabel("Ascl1 group")
axes[0].set_ylabel("Age")
axes[0].set_xlim(-0.5, 1.5)
axes[0].set_ylim(-0.4, len(age_order_rev_present) - 0.6)
axes[0].set_yticks(range(len(age_order_rev_present)))
axes[0].set_yticklabels(age_order_rev_present)
axes[0].legend(title="Cell count", bbox_to_anchor=(1.18, 0.5),
               loc="center left", frameon=True, fontsize=9)

sns.scatterplot(
    data=ct2_ct, x="x_pos", y=CELLTYPE2_COL,
    size="cell_count", hue="cell_count",
    palette="Reds", sizes=(60, 550),
    edgecolor="black", linewidth=0.6, ax=axes[1],
)
axes[1].set_title("Ascl1 GMM High/Low - Celltype2 (Female Hep)", weight="bold")
axes[1].set_xticks([0, 1])
axes[1].set_xticklabels(["Low", "High"])
axes[1].set_xlabel("Ascl1 group")
axes[1].set_ylabel("Celltype2")
axes[1].set_xlim(-0.5, 1.5)
axes[1].set_ylim(-0.4, len(KEEP_CT2) - 0.6)
axes[1].set_yticks(range(len(KEEP_CT2)))
axes[1].set_yticklabels(KEEP_CT2)
axes[1].legend(title="Cell count", bbox_to_anchor=(1.18, 0.5),
               loc="center left", frameon=True, fontsize=9)

plt.tight_layout()
fig.savefig(
    FIG_DIR / "Ascl1_GMM_HighLow_count_dotplots_female_hep_age_ct2.pdf",
    dpi=300, bbox_inches="tight",
)
plt.close(fig)
print(f"\n  [OK] {FIG_DIR}/Ascl1_GMM_HighLow_count_dotplots_female_hep_age_ct2.pdf")


# ==============================================================================
# 4. DIFFERENTIAL EXPRESSION: Ascl1-High vs Ascl1-Low (Wilcoxon)
# ==============================================================================

banner("STEP 4: Wilcoxon DE (Ascl1-High vs Ascl1-Low)")

# Fresh copy for DE
adata_hep_f_ct2_pos = adata_pos[
    adata_pos.obs[CELLTYPE2_COL].isin(KEEP_CT2)
].copy()

print(f"  Cells: {adata_hep_f_ct2_pos.n_obs:,}")
print(f"  Subtypes: {', '.join(KEEP_CT2)}")
print(f"\n  Group distribution:")
print(adata_hep_f_ct2_pos.obs["Ascl1_group"].value_counts().to_string())

# Rename High/Low for scanpy compatibility
adata_hep_f_ct2_pos.obs["Ascl1_group"] = (
    adata_hep_f_ct2_pos.obs["Ascl1_group"]
    .astype(str)
    .replace({"High": "Ascl1_high", "Low": "Ascl1_low"})
)

print("\n  Running Wilcoxon DE test (Ascl1_high vs Ascl1_low)...")
sc.tl.rank_genes_groups(
    adata_hep_f_ct2_pos,
    groupby="Ascl1_group",
    groups=["Ascl1_high"],
    reference="Ascl1_low",
    method="wilcoxon",
    use_raw=False,
    n_genes=adata_hep_f_ct2_pos.shape[1],
)
print("  [OK] DE complete")

de_res = sc.get.rank_genes_groups_df(adata_hep_f_ct2_pos, group="Ascl1_high")

sig_genes  = de_res[de_res["pvals_adj"] < PVAL_THRESHOLD]
up_genes   = sig_genes[sig_genes["logfoldchanges"] > 0]
down_genes = sig_genes[sig_genes["logfoldchanges"] < 0]

up_strong   = up_genes[up_genes["logfoldchanges"] > LOGFC_STRONG_CUTOFF]
down_strong = down_genes[down_genes["logfoldchanges"] < -LOGFC_STRONG_CUTOFF]

print(f"\n  Total genes tested: {len(de_res):,}")
print(f"  Significant (adj p < {PVAL_THRESHOLD}): {len(sig_genes):,} "
      f"({100 * len(sig_genes) / max(len(de_res), 1):.1f}%)")
print(f"    Up in Ascl1-High:   {len(up_genes):,}")
print(f"    Down in Ascl1-High: {len(down_genes):,}")
print(f"\n  Strong effects (|logFC| > {LOGFC_STRONG_CUTOFF}):")
print(f"    Upregulated:   {len(up_strong):,}")
print(f"    Downregulated: {len(down_strong):,}")

print(f"\n  Top 15 upregulated:")
print(up_genes.nlargest(15, "scores")[["names", "scores", "logfoldchanges", "pvals_adj"]]
      .rename(columns={"names": "Gene", "scores": "Score",
                       "logfoldchanges": "logFC", "pvals_adj": "adj_pval"})
      .to_string(index=False))

print(f"\n  Top 15 downregulated:")
print(down_genes.nsmallest(15, "scores")[["names", "scores", "logfoldchanges", "pvals_adj"]]
      .rename(columns={"names": "Gene", "scores": "Score",
                       "logfoldchanges": "logFC", "pvals_adj": "adj_pval"})
      .to_string(index=False))

de_res.to_csv(
    RES_DIR / "DE_Ascl1_GMM_high_vs_low_female_hep_CT2subset.csv", index=False,
)
sig_genes.to_csv(
    RES_DIR / "DE_Ascl1_GMM_high_vs_low_female_hep_CT2subset_significant.csv",
    index=False,
)
print(f"\n  [OK] {RES_DIR}/DE_Ascl1_GMM_high_vs_low_female_hep_CT2subset.csv")


# ==============================================================================
# 5. GSEA: REACTOME PATHWAYS 2024
# ==============================================================================

banner("STEP 5: GSEA Reactome Pathways 2024")

csv_path = RES_DIR / "DE_Ascl1_GMM_high_vs_low_female_hep_CT2subset.csv"
results_df = pd.read_csv(csv_path)
results_df = results_df.rename(columns={"names": "gene"})
results_df["gene"] = results_df["gene"].astype(str)

# Rank ALL genes by Wilcoxon z-score (no pre-filtering by significance)
results_df["Score"] = pd.to_numeric(results_df["scores"], errors="coerce")

rnk = (
    results_df[["gene", "Score"]]
    .dropna()
    .sort_values("Score", ascending=False)
    .copy()
)
rnk["gene"] = rnk["gene"].str.upper()  # Reactome uses human (uppercase) symbols

print(f"  Ranked genes: {len(rnk):,}")
print(f"  Score range: [{rnk['Score'].min():.2f}, {rnk['Score'].max():.2f}]")
print(f"  Positive scores (up in High): {int((rnk['Score'] > 0).sum()):,}")
print(f"  Negative scores (down in High): {int((rnk['Score'] < 0).sum()):,}")

print("\n  Running GSEA prerank (Reactome_Pathways_2024)...")
gsea = gp.prerank(
    rnk=rnk,
    gene_sets="Reactome_Pathways_2024",
    permutation_num=1000,
    min_size=2,
    max_size=2000,
    seed=42,
    verbose=False,
    outdir=None,
)
print("  [OK] GSEA complete")

gsea_df = gsea.res2d.copy()
gsea_df["FDR q-val"] = pd.to_numeric(gsea_df["FDR q-val"], errors="coerce")
gsea_df["NES"] = pd.to_numeric(gsea_df["NES"], errors="coerce")
gsea_df["NOM p-val"] = pd.to_numeric(gsea_df["NOM p-val"], errors="coerce")

if "Gene %" in gsea_df.columns:
    gsea_df["Gene %"] = (
        gsea_df["Gene %"].astype(str)
        .str.replace("%", "", regex=False)
        .astype(float)
    )
else:
    gsea_df["Gene %"] = np.nan

total_pathways = len(gsea_df)
sig_nom = gsea_df[gsea_df["NOM p-val"] < NOM_P_CUTOFF]
sig_fdr = gsea_df[gsea_df["FDR q-val"] < 0.25]
pos_nes = gsea_df[gsea_df["NES"] > 0]
neg_nes = gsea_df[gsea_df["NES"] < 0]
pos_sig = sig_nom[sig_nom["NES"] > 0]
neg_sig = sig_nom[sig_nom["NES"] < 0]

print(f"\n  Total pathways tested: {total_pathways:,}")
print(f"    Positive NES (enriched in Ascl1-High): {len(pos_nes):,}")
print(f"    Negative NES (enriched in Ascl1-Low):  {len(neg_nes):,}")
print(f"\n  Significant:")
print(f"    NOM p < {NOM_P_CUTOFF}: {len(sig_nom):,} "
      f"({100 * len(sig_nom) / total_pathways:.1f}%)")
print(f"      Positive NES: {len(pos_sig):,}")
print(f"      Negative NES: {len(neg_sig):,}")
print(f"    FDR q < 0.25: {len(sig_fdr):,} "
      f"({100 * len(sig_fdr) / total_pathways:.1f}%)")

print(f"\n  Top 10 positive NES:")
print(pos_nes.nlargest(10, "NES")[["Term", "NES", "NOM p-val", "FDR q-val"]]
      .to_string(index=False))

print(f"\n  Top 10 negative NES:")
print(neg_nes.nsmallest(10, "NES")[["Term", "NES", "NOM p-val", "FDR q-val"]]
      .to_string(index=False))

out_file = RES_DIR / "GSEA_Ascl1_GMM_high_vs_low_female_hep_CT2subset_Reactome2024.csv"
gsea_df.to_csv(out_file, index=False)
sig_nom.to_csv(
    RES_DIR / "GSEA_Ascl1_GMM_significant_NOMp05_Reactome2024.csv",
    index=False,
)
print(f"\n  [OK] {out_file}")
print(f"  [OK] {RES_DIR}/GSEA_Ascl1_GMM_significant_NOMp05_Reactome2024.csv")


# ==============================================================================
# 6. GSEA DOT PLOTS (TOP 30 POSITIVE NES, TOP 30 NEGATIVE NES)
# ==============================================================================

banner("STEP 6: GSEA dot plots (top 30 positive + top 30 negative NES)")

gsea = pd.read_csv(out_file)
gsea.columns = gsea.columns.str.strip()

if "Gene %" in gsea.columns:
    gsea["Gene %"] = (
        gsea["Gene %"].astype(str)
        .str.replace("%", "", regex=False)
        .astype(float)
    )
else:
    gsea["Gene %"] = np.nan

gsea["term_clean"] = gsea["Term"].apply(
    lambda t: re.sub(r"\(GO:\d+\)", "", str(t)).strip()
)

# Top 30 positive NES
pos_filt = gsea[(gsea["NES"] > 0) & (gsea["NOM p-val"] < NOM_P_CUTOFF)].copy()
top30_pos = pos_filt.sort_values("Gene %", ascending=False).head(TOP_N_GSEA)
print(f"  Positive NES pathways (NES > 0, NOM p < {NOM_P_CUTOFF}): {len(top30_pos)}")
plot_gsea_dotplot(
    top30_pos, direction="positive",
    outfile=FIG_DIR / "GSEA_GMM_top30_positiveNES_female_hep_CT2subset.pdf",
)

# Top 30 negative NES
neg_filt = gsea[(gsea["NES"] < 0) & (gsea["NOM p-val"] < NOM_P_CUTOFF)].copy()
top30_neg = neg_filt.sort_values("Gene %", ascending=False).head(TOP_N_GSEA)
print(f"  Negative NES pathways (NES < 0, NOM p < {NOM_P_CUTOFF}): {len(top30_neg)}")
plot_gsea_dotplot(
    top30_neg, direction="negative",
    outfile=FIG_DIR / "GSEA_GMM_top30_negativeNES_female_hep_CT2subset.pdf",
)


# ==============================================================================
# DONE
# ==============================================================================

banner("ALL STEPS COMPLETE")
print(f"  Figures: {FIG_DIR}/")
print(f"  Tables:  {RES_DIR}/")
