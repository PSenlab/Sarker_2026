#!/usr/bin/env python3
# ==============================================================================
# Hepatocyte Sub-Cluster Reactome Pathway Enrichment (Panels g/h/i)
# ==============================================================================
#
# Description:
#   For each hepatocyte sub-cluster (hep_01, hep_02, hep_04), loads the
#   per-k-means-cluster gene CSVs produced by STEP 4 of
#   hepatocyte_pseudobulk_deseq2_pipeline.R
#   (e.g. hep_01_cluster_1_genes_cluster.csv, hep_01_cluster_2_genes_cluster.csv),
#   converts mouse genes to human orthologs via mygene, and runs Enrichr
#   Reactome enrichment with a hepatocyte-expressed background.
#
#   Pathways are filtered by adjusted P < 0.05 and sorted by odds ratio.
#   One side-by-side dot plot per sub-cluster (cluster 1 left, cluster 2 right)
#   matching figure panels (g) hep_01, (h) hep_02, (i) hep_04.
#
# Input:
#   - integrated_scvi.h5ad - canonical AnnData (for hepatocyte background)
#   - Per-cluster gene CSVs from hepatocyte_pseudobulk_deseq2_pipeline.R STEP 4:
#       hep_01_cluster_1_genes_cluster.csv
#       hep_01_cluster_2_genes_cluster.csv
#       hep_02_cluster_1_genes_cluster.csv
#       hep_02_cluster_2_genes_cluster.csv
#       hep_04_cluster_1_genes_cluster.csv
#       hep_04_cluster_2_genes_cluster.csv
#
# Output:
#   - hep_<NN>_cluster_<k>_Reactome_results.csv  (full Enrichr output)
#   - panel_g_hep_01_reactome_dotplot.pdf/png
#   - panel_h_hep_02_reactome_dotplot.pdf/png
#   - panel_i_hep_04_reactome_dotplot.pdf/png
#
# Reference:
#   Sarker et al. (2026) Cell Metabolism
#
# ==============================================================================

import os
import numpy as np
import pandas as pd
import scanpy as sc
import gseapy as gp
import mygene
import matplotlib.pyplot as plt
import matplotlib as mpl
import warnings

warnings.filterwarnings("ignore")

# Global plot settings
mpl.rcParams["font.family"] = "Arial"
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42


# ==============================================================================
# CONFIGURATION - UPDATE THESE PATHS
# ==============================================================================

ADATA_PATH    = "integrated_scvi.h5ad"

# Directory containing the per-cluster CSVs from STEP 4 of the R pipeline
CSV_INPUT_DIR = "."

# Output directory
OUTPUT_DIR    = "hepatocyte_cluster_pathway_enrichment"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Sub-clusters to process (panels g, h, i)
PANEL_MAP = {
    "hep_01": "g",
    "hep_02": "h",
    "hep_04": "i",
}

# K-means clusters per sub-cluster (KMEANS_K = 2 in the R pipeline)
K_CLUSTERS = [1, 2]

# Cluster colors
CLUSTER_COLORS = {
    1: "#2ECC71",  # green
    2: "#E74C3C",  # red
}

# Analysis parameters
EXPRESSION_THRESHOLD = 0.05
GENE_SET_LIBRARY     = "Reactome_Pathways_2024"

# Pathway filters (panels g/h/i)
ADJ_P_THRESHOLD = 0.05
TOP_N_PATHWAYS  = 10


# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================

def convert_mouse_to_human(mouse_genes):
    """Convert mouse gene symbols to human orthologs via HomoloGene."""
    mouse_genes = sorted({g for g in mouse_genes if pd.notna(g) and g})
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

    return sorted({
        entry["symbol"].upper()
        for entry in human_conversion
        if isinstance(entry, dict) and "symbol" in entry
    })


def run_enrichr_analysis(gene_list, background, gene_set, title):
    """Run Enrichr with explicit background; return full results DataFrame."""
    gene_list = sorted({g.upper() for g in gene_list if isinstance(g, str)})
    background = sorted({g.upper() for g in background if isinstance(g, str)})
    gene_list = sorted(set(gene_list) & set(background))

    if len(gene_list) < 3:
        print(f"    [SKIP] Too few genes ({len(gene_list)}) for {title}")
        return pd.DataFrame()

    try:
        enr = gp.enrichr(
            gene_list=gene_list, gene_sets=gene_set,
            organism="Human", background=background,
            outdir=None, no_plot=True,
        )
        if enr is None or not hasattr(enr, "results"):
            return pd.DataFrame()

        df = enr.results.copy()
        if df is None or df.empty:
            return pd.DataFrame()

        input_set = set(gene_list)

        def recount_overlap(gene_str):
            if not isinstance(gene_str, str):
                return 0
            genes = {g.strip().upper() for g in gene_str.split(";") if g.strip()}
            return len(genes & input_set)

        df["Gene_Count"] = df["Genes"].apply(recount_overlap)
        df["-log10(padj)"] = -np.log10(df["Adjusted P-value"].clip(lower=1e-300))
        df["N_input_genes"] = len(gene_list)
        df["N_background_genes"] = len(background)

        print(f"    [OK] {title}: {len(df)} terms "
              f"({len(gene_list)} input / {len(background)} bg)")
        return df

    except Exception as e:
        print(f"    [ERROR] {title}: {e}")
        return pd.DataFrame()


def plot_panel(cluster_results, subcluster, panel_letter, outfile_prefix):
    """
    Side-by-side dot plot for one hepatocyte sub-cluster:
    cluster 1 (left) and cluster 2 (right). Filtered by adj P < 0.05,
    sorted by odds ratio, top N retained.
    """
    plot_data = {}
    for k in K_CLUSTERS:
        df = cluster_results.get(k)
        if df is None or df.empty:
            continue
        df_filt = df[
            (df["Adjusted P-value"] < ADJ_P_THRESHOLD)
            & (np.isfinite(df["Odds Ratio"]))
        ].copy()
        if df_filt.empty:
            print(f"    [SKIP] {subcluster} cluster {k}: "
                  f"no terms pass adj P < {ADJ_P_THRESHOLD}")
            continue
        df_top = df_filt.sort_values("Odds Ratio", ascending=False).head(TOP_N_PATHWAYS)
        plot_data[k] = df_top
        print(f"    {subcluster} cluster {k}: {len(df_top)} top terms")

    if not plot_data:
        print(f"    [SKIP] {subcluster}: no clusters had significant terms")
        return

    n_panels = len(plot_data)
    fig, axes = plt.subplots(
        1, n_panels, figsize=(5.5 * n_panels, 5.5),
        sharey=False,
    )
    if n_panels == 1:
        axes = [axes]

    for ax, (k, df) in zip(axes, sorted(plot_data.items())):
        df = df.copy()
        # Truncate long term names
        df["Term_short"] = df["Term"].apply(
            lambda t: (t[:55] + "...") if len(t) > 55 else t
        )
        df = df.sort_values("Odds Ratio", ascending=True)  # bottom-up

        max_count = df["Gene_Count"].max() if len(df) else 1

        ax.scatter(
            df["Odds Ratio"], df["Term_short"],
            s=(df["Gene_Count"] / max_count) * 400 + 60,
            c=CLUSTER_COLORS.get(k, "gray"),
            alpha=0.85, edgecolors="black", linewidths=0.4,
        )
        ax.set_xlabel("Odds Ratio", fontsize=11, fontweight="bold")
        ax.set_title(
            f"{subcluster.replace('_', ' ').title()} - cluster {k}",
            fontsize=12, fontweight="bold",
        )
        ax.tick_params(axis="y", labelsize=9)
        ax.tick_params(axis="x", labelsize=9)
        ax.grid(axis="x", linestyle="--", alpha=0.3)
        ax.set_axisbelow(True)
        for spine in ["top", "right"]:
            ax.spines[spine].set_visible(False)
        ax.margins(y=0.05)

    fig.suptitle(
        f"({panel_letter}) Top Reactome pathways - "
        f"{subcluster.replace('_', ' ').title()} "
        f"(adj P < {ADJ_P_THRESHOLD}, sorted by odds ratio)",
        fontsize=12, fontweight="bold", y=1.02,
    )
    plt.tight_layout()

    pdf_out = f"{outfile_prefix}.pdf"
    png_out = f"{outfile_prefix}.png"
    plt.savefig(pdf_out, dpi=300, bbox_inches="tight")
    plt.savefig(png_out, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"    [OK] {pdf_out}")


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
print(f"  Hepatocytes:  {adata.n_obs:,} cells")

avg_expr = np.asarray(adata.X.mean(axis=0)).flatten()
expressed_genes = adata.var_names[avg_expr > EXPRESSION_THRESHOLD].tolist()
print(f"  Background mouse genes (mean expr > {EXPRESSION_THRESHOLD}): "
      f"{len(expressed_genes):,}")

print("\n  Converting background to human orthologs...")
background_human = convert_mouse_to_human(expressed_genes)
print(f"  Background human genes: {len(background_human):,}")


# ==============================================================================
# STEP 2: PROCESS EACH SUB-CLUSTER (panels g, h, i)
# ==============================================================================

print()
print("=" * 70)
print("STEP 2: Reactome enrichment per hepatocyte sub-cluster")
print("=" * 70)

for subcluster, panel_letter in PANEL_MAP.items():

    print(f"\n--- {subcluster.upper()} (panel {panel_letter}) ---")

    cluster_results = {}

    for k in K_CLUSTERS:

        csv_name = f"{subcluster}_cluster_{k}_genes_cluster.csv"
        csv_path = os.path.join(CSV_INPUT_DIR, csv_name)

        if not os.path.exists(csv_path):
            print(f"    [SKIP] {csv_name} not found")
            continue

        df_genes = pd.read_csv(csv_path)
        if "gene" not in df_genes.columns:
            print(f"    [SKIP] {csv_name}: 'gene' column missing")
            continue

        mouse_genes = df_genes["gene"].dropna().astype(str).tolist()
        print(f"    Cluster {k}: {len(mouse_genes)} mouse genes")

        human_genes = convert_mouse_to_human(mouse_genes)
        print(f"    Cluster {k}: {len(human_genes)} human orthologs")

        result = run_enrichr_analysis(
            gene_list=human_genes,
            background=background_human,
            gene_set=GENE_SET_LIBRARY,
            title=f"{subcluster} cluster {k}",
        )

        if not result.empty:
            out_csv = os.path.join(
                OUTPUT_DIR,
                f"{subcluster}_cluster_{k}_Reactome_results.csv",
            )
            result.to_csv(out_csv, index=False)
            print(f"    [OK] {os.path.basename(out_csv)}")
            cluster_results[k] = result

    # Build the panel figure
    plot_panel(
        cluster_results=cluster_results,
        subcluster=subcluster,
        panel_letter=panel_letter,
        outfile_prefix=os.path.join(
            OUTPUT_DIR,
            f"panel_{panel_letter}_{subcluster}_reactome_dotplot",
        ),
    )


# ==============================================================================
# SUMMARY
# ==============================================================================

print()
print("=" * 70)
print("COMPLETE")
print("=" * 70)
print(f"  Panels: {', '.join(f'({v}) {k}' for k, v in PANEL_MAP.items())}")
print(f"  Filter: adjusted P < {ADJ_P_THRESHOLD}, "
      f"sorted by Odds Ratio, top {TOP_N_PATHWAYS} per cluster")
print(f"  Output: {OUTPUT_DIR}/")
print("=" * 70)
