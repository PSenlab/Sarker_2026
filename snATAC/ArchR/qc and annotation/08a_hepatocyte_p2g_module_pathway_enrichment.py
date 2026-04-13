#!/usr/bin/env python3
# ==============================================================================
# Hepatocyte P2G Cluster Pathway Enrichment
# ==============================================================================
#
# Description:
#   Runs Reactome pathway enrichment on hepatocyte P2G k-means clusters
#   produced by archr_p2g_analysis.R (Step 12). For hepatocytes the P2G
#   heatmap is clustered with k=4. Mouse genes from each cluster are
#   converted to human orthologs via mygene, then run through Enrichr
#   with a hepatocyte-expressed background.
#
# Input:
#   - AnnData (.h5ad) with celltype labels - used to derive the
#     hepatocyte-expressed background gene set
#   - Per-cluster gene CSVs from Step 12.4:
#       P2G_C1_genes.csv, P2G_C2_genes.csv, P2G_C3_genes.csv, P2G_C4_genes.csv
#
# Output:
#   - Per-cluster Reactome results CSV
#     (P2G_Module{N}_Reactome_Pathways_2024_results.csv)
#   - Combined dot plot across all modules
#     (p2g_module_pathway_dotplot.pdf/png)
#
# Reference:
#   Sarker et al. (2026) Cell Metabolism
#
# ==============================================================================

import os
import glob
import numpy as np
import pandas as pd
import scanpy as sc
import gseapy as gp
import mygene
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.lines import Line2D
import warnings

warnings.filterwarnings("ignore")

# Global plot settings
mpl.rcParams["font.family"] = "Arial"
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42


# ==============================================================================
# CONFIGURATION - UPDATE THESE PATHS
# ==============================================================================

# AnnData for background gene extraction
ADATA_PATH = "integrated_scvi.h5ad"

# Directory containing the per-cluster CSVs produced by archr_p2g_analysis.R
# (Step 12.4 writes P2G_C1_genes.csv .. P2G_CK_genes.csv into its OUTPUT_DIR)
P2G_RESULTS_DIR = "path/to/Step12_P2G_Analysis"

# Output directory for this pathway enrichment step
OUTPUT_DIR = "p2g_module_pathway_enrichment"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Analysis parameters
EXPRESSION_THRESHOLD = 0.05
GENE_SET_LIBRARY = "Reactome_Pathways_2024"

# Filters for plotting
PVAL_THRESHOLD = 0.01
MIN_GENE_COUNT = 5
TOP_N_PER_MODULE = 10

# Module colors (k=4 for hepatocyte)
MODULE_COLORS = {
    1: "#E41A1C",  # Red
    2: "#377EB8",  # Blue
    3: "#4DAF4A",  # Green
    4: "#984EA3",  # Purple
}


# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================

def convert_mouse_to_human(mouse_genes):
    """Convert mouse gene symbols to human orthologs via HomoloGene."""
    mouse_genes = list(set(g for g in mouse_genes if pd.notna(g) and g))
    if not mouse_genes:
        return []

    mg = mygene.MyGeneInfo()
    gene_conversion = mg.querymany(
        mouse_genes, scopes="symbol", fields="homologene",
        species="mouse", as_dataframe=False, verbose=False,
    )

    human_entrez_ids = set()
    for entry in gene_conversion:
        if isinstance(entry, dict) and "homologene" in entry:
            for tax_id, entrez_id in entry["homologene"].get("genes", []):
                if tax_id == 9606:
                    human_entrez_ids.add(entrez_id)

    if not human_entrez_ids:
        return []

    human_conversion = mg.querymany(
        list(human_entrez_ids), scopes="entrezgene", fields="symbol",
        species="human", as_dataframe=False, verbose=False,
    )

    human_symbols = sorted({
        entry["symbol"].upper()
        for entry in human_conversion
        if isinstance(entry, dict) and "symbol" in entry
    })
    return human_symbols


def run_enrichr_analysis(gene_list, background, gene_set, title, outfile_prefix):
    """Run Enrichr with explicit background, return full results DataFrame."""
    if len(gene_list) < 5:
        print(f"  [SKIP] Too few genes ({len(gene_list)}) for {title}")
        return pd.DataFrame()

    try:
        enr = gp.enrichr(
            gene_list=gene_list,
            gene_sets=gene_set,
            organism="Human",
            background=background,
            outdir=None,
            no_plot=True,
        )

        if enr is None or not hasattr(enr, "results"):
            print(f"  [SKIP] No Enrichr results for {title}")
            return pd.DataFrame()

        df = enr.results.copy()
        if df is None or df.empty:
            print(f"  [SKIP] No enrichment terms for {title}")
            return pd.DataFrame()

        # Recount gene overlap against input gene list
        input_set = {g.upper() for g in gene_list}

        def recount_overlap(gene_str):
            if not isinstance(gene_str, str) or gene_str.strip() == "":
                return 0
            genes = {g.strip().upper() for g in gene_str.split(";") if g.strip()}
            return len(genes & input_set)

        df["Gene_Count"] = df["Genes"].apply(recount_overlap)
        df["-log10(pval)"] = -np.log10(df["P-value"].clip(lower=1e-300))
        df["-log10(padj)"] = -np.log10(df["Adjusted P-value"].clip(lower=1e-300))
        df["GeneSet"] = gene_set

        out_csv = f"{outfile_prefix}_results.csv"
        df.to_csv(out_csv, index=False)
        print(f"  [OK] {title}: {len(df)} terms -> {os.path.basename(out_csv)}")
        return df

    except Exception as e:
        print(f"  [ERROR] {title}: {e}")
        return pd.DataFrame()


# ==============================================================================
# STEP 1: LOAD ADATA AND DERIVE BACKGROUND
# ==============================================================================

print()
print("=" * 70)
print("STEP 1: Load AnnData and derive hepatocyte background")
print("=" * 70)

adata = sc.read_h5ad(ADATA_PATH)
print(f"  Full dataset: {adata.n_obs:,} cells x {adata.n_vars:,} genes")

adata = adata[adata.obs["celltype"] == "Hepatocyte"].copy()
print(f"  Hepatocytes: {adata.n_obs:,} cells")

expr = adata.X
gene_names = adata.var_names
avg_expr = np.asarray(expr.mean(axis=0)).flatten()
expressed_genes = gene_names[avg_expr > EXPRESSION_THRESHOLD].tolist()
print(f"  Background genes (mean expr > {EXPRESSION_THRESHOLD}): {len(expressed_genes)}")


# ==============================================================================
# STEP 2: LOAD P2G CLUSTER GENE FILES
# ==============================================================================

print()
print("=" * 70)
print("STEP 2: Load P2G cluster gene files")
print("=" * 70)

# Matches archr_p2g_analysis.R output: P2G_C1_genes.csv .. P2G_CK_genes.csv
cluster_files = sorted(glob.glob(os.path.join(P2G_RESULTS_DIR, "P2G_C*_genes.csv")))
print(f"  Found {len(cluster_files)} cluster files in {P2G_RESULTS_DIR}")

if not cluster_files:
    raise FileNotFoundError(
        f"No P2G_C*_genes.csv files found in {P2G_RESULTS_DIR}. "
        "Run archr_p2g_analysis.R (Step 12) first."
    )

cluster_genes = {}
for f in cluster_files:
    # Extract cluster number from "P2G_C1_genes.csv"
    basename = os.path.basename(f)
    cluster_num = int(basename.split("_C")[1].split("_")[0])
    df = pd.read_csv(f)
    genes = df["gene"].dropna().astype(str).tolist()
    cluster_genes[cluster_num] = genes
    print(f"  Cluster {cluster_num}: {len(genes)} genes")

K = len(cluster_genes)
print(f"  K = {K} clusters loaded")

if K != 4:
    print(f"  [WARN] Expected K=4 for hepatocyte, found K={K}")


# ==============================================================================
# STEP 3: CONVERT MOUSE TO HUMAN
# ==============================================================================

print()
print("=" * 70)
print("STEP 3: Convert mouse genes to human orthologs")
print("=" * 70)

module_genes_human = {}
for m in sorted(cluster_genes.keys()):
    human_genes = convert_mouse_to_human(cluster_genes[m])
    module_genes_human[m] = human_genes
    print(f"  Module {m}: {len(cluster_genes[m])} mouse -> {len(human_genes)} human")

print(f"\n  Converting background genes...")
background_human = convert_mouse_to_human(expressed_genes)
print(f"  Background: {len(expressed_genes)} mouse -> {len(background_human)} human")


# ==============================================================================
# STEP 4: RUN REACTOME ENRICHMENT PER MODULE
# ==============================================================================

print()
print("=" * 70)
print("STEP 4: Reactome enrichment per module")
print("=" * 70)

all_results = {}
for m in sorted(cluster_genes.keys()):
    print(f"\n  Module {m} ({len(module_genes_human[m])} human genes)")
    result = run_enrichr_analysis(
        gene_list=module_genes_human[m],
        background=background_human,
        gene_set=GENE_SET_LIBRARY,
        title=f"Module {m}",
        outfile_prefix=os.path.join(
            OUTPUT_DIR, f"P2G_Module{m}_{GENE_SET_LIBRARY}"
        ),
    )
    if not result.empty:
        all_results[m] = result


# ==============================================================================
# STEP 5: COMBINED DOT PLOT
# ==============================================================================

print()
print("=" * 70)
print("STEP 5: Combined dot plot across modules")
print("=" * 70)

combined_data = []
for m, df in all_results.items():
    df_filtered = df[
        (df["P-value"] < PVAL_THRESHOLD)
        & (df["Gene_Count"] > MIN_GENE_COUNT)
        & (np.isfinite(df["Odds Ratio"]))
    ].copy()

    if df_filtered.empty:
        print(f"  Module {m}: no terms passed filter "
              f"(P < {PVAL_THRESHOLD}, Gene_Count > {MIN_GENE_COUNT})")
        continue

    df_top = df_filtered.nlargest(TOP_N_PER_MODULE, "Odds Ratio").copy()
    df_top["Module"] = m
    combined_data.append(df_top)
    print(f"  Module {m}: {len(df_top)} top terms selected")

if not combined_data:
    print("  [SKIP] No modules had significant terms - no dot plot generated")
else:
    plot_df = pd.concat(combined_data, ignore_index=True)
    plot_df = plot_df[np.isfinite(plot_df["Odds Ratio"])].copy()
    print(f"\n  Total terms for plotting: {len(plot_df)}")

    # Sort terms by max odds ratio (high to low)
    term_max_or = plot_df.groupby("Term")["Odds Ratio"].max().sort_values(ascending=False)
    all_terms_sorted = term_max_or.index.tolist()

    # Truncate long term names
    term_labels = [
        (t[:55] + "...") if len(t) > 55 else t
        for t in all_terms_sorted
    ]

    max_gene_count = plot_df["Gene_Count"].max()
    or_min = plot_df["Odds Ratio"].min()
    or_max = plot_df["Odds Ratio"].max()

    fig, ax = plt.subplots(figsize=(7, max(7, len(all_terms_sorted) * 0.25)))

    for i, term in enumerate(all_terms_sorted):
        term_data = plot_df[plot_df["Term"] == term]
        for _, row in term_data.iterrows():
            size = (row["Gene_Count"] / max_gene_count) * 500 + 50
            ax.scatter(
                row["Odds Ratio"], i,
                s=size,
                c=MODULE_COLORS.get(row["Module"], "gray"),
                alpha=0.85,
                edgecolors="darkgray", linewidths=0.3,
            )

    ax.set_xlim(0, or_max * 1.1)
    ax.set_xlabel("Odds Ratio", fontsize=12, fontweight="bold")
    ax.set_yticks(range(len(all_terms_sorted)))
    ax.set_yticklabels(term_labels, fontsize=9)
    ax.set_title(
        f"P2G Module Pathway Enrichment\n"
        f"(P < {PVAL_THRESHOLD}, Gene Count > {MIN_GENE_COUNT}, "
        f"top {TOP_N_PER_MODULE} per module by Odds Ratio)",
        fontsize=12, fontweight="bold", pad=10,
    )

    # Light horizontal guides
    for i in range(len(all_terms_sorted)):
        ax.axhline(y=i, color="lightgray", linestyle="-", linewidth=0.3, zorder=0)

    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)

    # Module color legend
    present_modules = sorted(plot_df["Module"].unique())
    module_handles = [
        Line2D(
            [0], [0], marker="o", color="w",
            markerfacecolor=MODULE_COLORS.get(m, "gray"),
            markersize=10, markeredgecolor="gray", markeredgewidth=0.3,
            label=f"{m}",
        )
        for m in present_modules
    ]
    legend1 = ax.legend(
        handles=module_handles, title="Module",
        loc="upper left", bbox_to_anchor=(1.02, 1), framealpha=0.95,
    )
    ax.add_artist(legend1)

    # Gene count size legend
    size_values = [15, 30, 45, 60]
    size_handles = [
        Line2D(
            [0], [0], marker="o", color="w", markerfacecolor="lightgray",
            markersize=np.sqrt((v / max_gene_count) * 500 + 50) / 2,
            markeredgecolor="gray", markeredgewidth=0.3, label=f"{v}",
        )
        for v in size_values
    ]
    ax.legend(
        handles=size_handles, title="Gene Count",
        loc="upper left", bbox_to_anchor=(1.02, 0.6), framealpha=0.95,
    )

    plt.tight_layout()
    pdf_out = os.path.join(OUTPUT_DIR, "p2g_module_pathway_dotplot.pdf")
    png_out = os.path.join(OUTPUT_DIR, "p2g_module_pathway_dotplot.png")
    plt.savefig(pdf_out, dpi=300, bbox_inches="tight")
    plt.savefig(png_out, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"  [OK] {pdf_out}")
    print(f"  [OK] {png_out}")


# ==============================================================================
# SUMMARY
# ==============================================================================

print()
print("=" * 70)
print("COMPLETE")
print("=" * 70)
print(f"  Modules processed: {K}")
print(f"  Filters: P < {PVAL_THRESHOLD}, Gene_Count > {MIN_GENE_COUNT}, "
      f"top {TOP_N_PER_MODULE} per module by Odds Ratio")
print(f"  Output: {OUTPUT_DIR}/")
print("=" * 70)
