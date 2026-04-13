#!/usr/bin/env python3
# ==============================================================================
# Ascl1 Expression Analysis in Liver Aging (Female Hepatocytes)
# ==============================================================================
#
# Description:
#   Profiles Ascl1 expression across the full single-nucleus multi-ome liver
#   aging dataset, then focuses on female hepatocytes where Ascl1 shows a
#   sex-specific aging signal. Performs differential expression between
#   Ascl1+ and Ascl1- female hepatocytes (Wilcoxon), generates a volcano
#   plot with functional category annotations, and runs GSEA against
#   Reactome Pathways 2024.
#
# Pipeline:
#   1. Load data
#   2. Global Ascl1 profiling (sex / age / celltype)
#   3. Female-only Ascl1 analysis
#   4. Female hepatocyte Ascl1 analysis (age / celltype2)
#   5. Count-based dot plots
#   6. Ascl1+ subset characterization
#   7. Wilcoxon DE: Ascl1+ vs Ascl1-
#   8. Volcano plot with functional categories
#   9. GSEA Reactome Pathways 2024
#  10. GSEA dot plot (top 30 by Gene %)
#
# Input:
#   - integrated_scvi.h5ad  (annotated AnnData with celltype, celltype2, sex, age)
#
# Output:
#   - Figures in figures/
#   - Tables in results/
#
# Reference:
#   Sarker et al. (2026) Cell Metabolism
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
from matplotlib.patches import Patch

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

# Ordering
AGE_ORDER = ["young", "mid_age", "old", "pre_geriatric", "geriatric"]
AGE_ALIASES = {"midage": "mid_age", "pregeriatric": "pre_geriatric"}

GENE = "Ascl1"

# DE thresholds
PVAL_THRESHOLD  = 0.05
LOGFC_THRESHOLD = 0  # adj-p drives significance; any direction counts


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

banner("STEP 4: Female hepatocyte Ascl1 analysis (age + celltype2)")

adata_fh = adata[
    (adata.obs[SEX_COL] == "female")
    & (adata.obs[CELLTYPE_COL] == "Hepatocyte")
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
print(f"\n  By Celltype2:\n{ct2_summary.to_string(index=False)}")

age_summary_fh.to_csv(RES_DIR / "Ascl1_summary_female_hep_by_age.csv", index=False)
ct2_summary.to_csv(RES_DIR / "Ascl1_summary_female_hep_by_celltype2.csv", index=False)

fig, axes = plt.subplots(1, 2, figsize=(10, 4))
sc.pl.dotplot(adata_fh, var_names=[GENE], groupby=AGE_COL,
              ax=axes[0], show=False, cmap="Reds")
axes[0].set_title("Ascl1 by Age (Female hepatocytes)")
sc.pl.dotplot(adata_fh, var_names=[GENE], groupby=CELLTYPE2_COL,
              ax=axes[1], show=False, cmap="Reds")
axes[1].set_title("Ascl1 by Celltype2 (Female hepatocytes)")
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
axes[1].set_title("Ascl1+ count by Celltype2", weight="bold", fontsize=12)
axes[1].set_xlabel("")
axes[1].set_ylabel("Celltype2")
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
# 7. DIFFERENTIAL EXPRESSION: Ascl1+ vs Ascl1- (Wilcoxon)
# ==============================================================================

banner("STEP 7: Wilcoxon DE (Ascl1+ vs Ascl1-) in female hepatocytes")

adata_fh.obs["Ascl1_status"] = pd.Categorical(
    np.where(adata_fh.obs["Ascl1_expr"] > 0, "Ascl1_pos", "Ascl1_neg"),
    categories=["Ascl1_pos", "Ascl1_neg"], ordered=True,
)

print("  Group sizes:")
print(adata_fh.obs["Ascl1_status"].value_counts().to_string())

print("\n  Running Wilcoxon rank-sum test...")
sc.tl.rank_genes_groups(
    adata_fh,
    groupby="Ascl1_status",
    groups=["Ascl1_pos"],
    reference="Ascl1_neg",
    method="wilcoxon",
    use_raw=False,
)
print("  [OK] DE complete")

de_res = sc.get.rank_genes_groups_df(adata_fh, group="Ascl1_pos")

sig_genes  = de_res[de_res["pvals_adj"] < PVAL_THRESHOLD]
up_genes   = sig_genes[sig_genes["logfoldchanges"] > 0]
down_genes = sig_genes[sig_genes["logfoldchanges"] < 0]

print(f"\n  Total genes tested: {len(de_res):,}")
print(f"  Significant (adj p < {PVAL_THRESHOLD}): {len(sig_genes):,}")
print(f"    Up in Ascl1+:   {len(up_genes):,}")
print(f"    Down in Ascl1+: {len(down_genes):,}")

de_res.to_csv(RES_DIR / "wilcoxon_Ascl1_pos_vs_neg_female_hepatocytes.csv", index=False)
sig_genes.to_csv(
    RES_DIR / "wilcoxon_Ascl1_pos_vs_neg_female_hepatocytes_significant.csv",
    index=False,
)
print(f"\n  [OK] {RES_DIR}/wilcoxon_Ascl1_pos_vs_neg_female_hepatocytes.csv")


# ==============================================================================
# 8. VOLCANO PLOT
# ==============================================================================

banner("STEP 8: Volcano plot with functional gene categories")

de_res = pd.read_csv(RES_DIR / "wilcoxon_Ascl1_pos_vs_neg_female_hepatocytes.csv")

epsilon = 1e-300
de_res["-log10_pval_adj"] = -np.log10(de_res["pvals_adj"] + epsilon)

de_res["sig"] = "Not Sig"
de_res.loc[
    (de_res["pvals_adj"] < PVAL_THRESHOLD)
    & (de_res["logfoldchanges"] > LOGFC_THRESHOLD),
    "sig",
] = "Up in Ascl1+"
de_res.loc[
    (de_res["pvals_adj"] < PVAL_THRESHOLD)
    & (de_res["logfoldchanges"] < -LOGFC_THRESHOLD),
    "sig",
] = "Down in Ascl1+"

print(de_res["sig"].value_counts().to_string())

sig_colors = {
    "Not Sig":        "lightgray",
    "Up in Ascl1+":   "lightpink",
    "Down in Ascl1+": "#008080",
}

# Functional gene categories to label on the volcano
genes_to_label = {
    "Cell Cycle (Up)": {
        "genes": ["Mki67", "Ccne2", "Ccnd1", "Ccnd2", "Ccng1", "Ccnf"],
        "color": "#8B0000",
    },
    "Translation (Up)": {
        "genes": ["Mrps14", "Mrps36", "Rps13", "Rpl5"],
        "color": "#4B0082",
    },
    "mTOR Pathway (Up)": {
        "genes": ["Lamtor5", "Lamtor4", "Lamtor3", "Lamtor2", "Slc38a9"],
        "color": "darkgreen",
    },
    "Immediate Early (Down)": {
        "genes": ["Jun", "Fos", "Egr1", "Jund", "Atf3"],
        "color": "#DC143C",
    },
    "Cell Cycle (Down)": {
        "genes": ["Wee1", "Plk3"],
        "color": "#FF1493",
    },
    "Hepatic Function (Down)": {
        "genes": ["Aldob", "G6pc", "Cps1", "Gys2", "Cyp2c37",
                  "Cyp2c67", "Abcc3", "Scd1", "Fasn", "Hmgcr"],
        "color": "#FF7F50",
    },
    "Others": {
        "genes": ["Lifr", "Esr1", "Angpt1", "Cd36", "Vwf",
                  "Fgl1", "Hif1a", "Hamp", "Hamp2"],
        "color": "#8B4513",
    },
}

fig, ax = plt.subplots(figsize=(20, 10))

for category, color in sig_colors.items():
    sub = de_res[de_res["sig"] == category]
    ax.scatter(
        sub["logfoldchanges"], sub["-log10_pval_adj"],
        c=color, label=category, s=24, alpha=0.7, edgecolors="none",
    )

ax.axhline(-np.log10(PVAL_THRESHOLD), color="black", ls="--", lw=0.8, alpha=0.6)
ax.axvline(0, color="black", ls="--", lw=0.8, alpha=0.6)

gene_positions = []
for category, info in genes_to_label.items():
    for g in info["genes"]:
        row = de_res[de_res["names"] == g]
        if not row.empty:
            r = row.iloc[0]
            gene_positions.append({
                "gene": r["names"],
                "x": r["logfoldchanges"],
                "y": r["-log10_pval_adj"],
                "color": info["color"],
            })

try:
    from adjustText import adjust_text
    texts = []
    for pos in gene_positions:
        t = ax.text(
            pos["x"], pos["y"], pos["gene"],
            fontsize=11, fontweight="bold", color=pos["color"],
            ha="center", va="center",
            bbox=dict(boxstyle="round,pad=0.3", facecolor="white",
                      edgecolor="none", alpha=0.85),
            zorder=10,
        )
        texts.append(t)
    adjust_text(
        texts,
        arrowprops=dict(arrowstyle="-", color="gray", lw=0.5, alpha=0.6),
        expand_points=(1.5, 1.5), expand_text=(1.2, 1.2),
        force_points=(0.5, 0.5), force_text=(0.5, 0.5),
        ax=ax,
    )
except ImportError:
    print("  [INFO] adjustText not installed; using offset placement")
    for pos in gene_positions:
        ax.annotate(
            pos["gene"], (pos["x"], pos["y"]),
            fontsize=11, fontweight="bold", color=pos["color"],
            textcoords="offset points", xytext=(0, 10), ha="center",
            bbox=dict(boxstyle="round,pad=0.3", facecolor="white",
                      edgecolor="none", alpha=0.85),
        )

ax.set_xlim(-5, 5)
ax.set_xlabel("Log2 fold change (Ascl1+ vs Ascl1-)", fontsize=13, fontweight="bold")
ax.set_ylabel("-log10(adjusted p-value)", fontsize=13, fontweight="bold")
ax.set_title(
    "Ascl1+ vs Ascl1- female hepatocytes (Wilcoxon rank-sum)",
    fontsize=15, fontweight="bold",
)

sig_legend = ax.legend(
    frameon=True, loc="lower right", bbox_to_anchor=(1.15, 0),
    title="Significance", fontsize=12, title_fontsize=13, markerscale=2.5,
)

functional_patches = [
    Patch(facecolor=info["color"], label=category)
    for category, info in genes_to_label.items()
]
ax.legend(
    handles=functional_patches, frameon=True, framealpha=0.9,
    loc="upper right", bbox_to_anchor=(1.18, 1),
    title="Functional categories", fontsize=12, title_fontsize=13,
)
ax.add_artist(sig_legend)

plt.tight_layout()
fig.savefig(FIG_DIR / "volcano_Ascl1_hepatocytes.pdf", dpi=300, bbox_inches="tight")
fig.savefig(FIG_DIR / "volcano_Ascl1_hepatocytes.png", dpi=300, bbox_inches="tight")
plt.close(fig)
print(f"  [OK] {FIG_DIR}/volcano_Ascl1_hepatocytes.pdf")


# ==============================================================================
# 9. GSEA (REACTOME PATHWAYS 2024)
# ==============================================================================

banner("STEP 9: GSEA Reactome Pathways 2024")

de_res = pd.read_csv(RES_DIR / "wilcoxon_Ascl1_pos_vs_neg_female_hepatocytes.csv")
de_res = de_res.rename(columns={"names": "gene"})
de_res["gene"] = de_res["gene"].astype(str)

# Small jitter (seeded) to break ties required by prerank
rng = np.random.default_rng(42)
de_res["Score"] = de_res["scores"] + rng.normal(0, 1e-6, len(de_res))

rnk = de_res[["gene", "Score"]].dropna().sort_values("Score", ascending=False)
rnk["gene"] = rnk["gene"].str.upper()  # Reactome uses human (uppercase) symbols
print(f"  Ranked {len(rnk):,} genes for GSEA")

gsea = gp.prerank(
    rnk=rnk,
    gene_sets="Reactome_Pathways_2024",
    permutation_num=1000,
    min_size=2,
    max_size=2000,
    seed=42,
    outdir=None,
    verbose=False,
)

gsea_df = gsea.res2d.copy()
gsea_df["FDR q-val"] = pd.to_numeric(gsea_df["FDR q-val"], errors="coerce")
if "Gene %" in gsea_df.columns:
    gsea_df["Gene %"] = (
        gsea_df["Gene %"].astype(str).str.replace("%", "", regex=False).astype(float)
    )
else:
    gsea_df["Gene %"] = 1.0

out_file = RES_DIR / "GSEA_Ascl1_pos_neg_female_hepatocytes_Reactome2024.csv"
gsea_df.to_csv(out_file, index=False)
print(f"  [OK] {out_file} ({len(gsea_df)} pathways)")


# ==============================================================================
# 10. GSEA DOT PLOT (TOP 30 BY GENE %)
# ==============================================================================

banner("STEP 10: GSEA dot plot (top 30 by Gene %)")

gsea = pd.read_csv(out_file)
gsea.columns = gsea.columns.str.strip()

if "Gene %" in gsea.columns:
    gsea["Gene %"] = (
        gsea["Gene %"].astype(str).str.replace("%", "", regex=False).astype(float)
    )
else:
    gsea["Gene %"] = np.nan

# Filter: NOM p < 0.05, top 30 by Gene %
filtered = gsea[gsea["NOM p-val"] < 0.05].copy()
top = filtered.sort_values("Gene %", ascending=False).head(30)
print(f"  Selected {len(top)} Reactome terms (NOM p < 0.05, top 30 by Gene %)")

top["term_clean"] = top["Term"].apply(
    lambda t: re.sub(r"\(GO:\d+\)", "", str(t)).strip()
)

top["pval_size"] = -np.log10(top["NOM p-val"].replace(0, 1e-10))

min_size, max_size = 60, 400
pmin, pmax = top["pval_size"].min(), top["pval_size"].max()
if pmax > pmin:
    top["size_scaled"] = (
        min_size + (top["pval_size"] - pmin) / (pmax - pmin) * (max_size - min_size)
    )
else:
    top["size_scaled"] = (min_size + max_size) / 2

top = top.sort_values("NES", ascending=False)
top["term_clean"] = pd.Categorical(
    top["term_clean"], categories=top["term_clean"], ordered=True,
)

nes_min, nes_max = top["NES"].min(), top["NES"].max()
nes_range = f"NES range: {nes_min:.2f} to {nes_max:.2f}"

plt.figure(figsize=(7, 9))
ax = sns.scatterplot(
    data=top, x="NES", y="term_clean",
    size="size_scaled", hue="Gene %",
    palette="RdYlBu_r", sizes=(min_size, max_size),
    edgecolor="black", linewidth=0.6,
)

plt.xlim(nes_min - 0.1, nes_max + 0.1)
plt.axvline(0, color="gray", lw=1)
plt.xlabel("Normalized enrichment score (NES)", fontsize=12, weight="bold")
plt.ylabel("Reactome pathway", fontsize=12, weight="bold")
plt.title(
    f"Top 30 Reactome pathways by Gene %\n"
    f"Ascl1+ vs Ascl1- (Female hepatocytes)\n{nes_range}",
    fontsize=14, weight="bold",
)

norm = plt.Normalize(top["Gene %"].min(), top["Gene %"].max())
sm = plt.cm.ScalarMappable(cmap="RdYlBu_r", norm=norm)
sm.set_array([])
cbar = plt.colorbar(sm, ax=ax)
cbar.set_label("Gene %", fontsize=11, weight="bold")

ax.legend(
    title="-log10(NOM p-val)",
    bbox_to_anchor=(1.42, 1), loc="upper left",
    frameon=True, fontsize=9, title_fontsize=10,
)

plt.tight_layout()
plt.savefig(
    FIG_DIR / "GSEA_Ascl1_pos_neg_female_hepatocytes_Reactome2024.pdf",
    dpi=300, bbox_inches="tight",
)
plt.close()
print(f"  [OK] {FIG_DIR}/GSEA_Ascl1_pos_neg_female_hepatocytes_Reactome2024.pdf")


# ==============================================================================
# DONE
# ==============================================================================

banner("ALL STEPS COMPLETE")
print(f"  Figures: {FIG_DIR}/")
print(f"  Tables:  {RES_DIR}/")
