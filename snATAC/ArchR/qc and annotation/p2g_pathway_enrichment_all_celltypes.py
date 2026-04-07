#!/usr/bin/env python3
# ==============================================================================
# P2G-Linked Gene Pathway Enrichment for Non-Hepatocyte Cell Types
# ==============================================================================
#
# Description:
#   Runs Reactome pathway enrichment on ArchR P2G-linked genes for each
#   non-hepatocyte cell type (Endothelial, Stellate, Cholangiocyte, Kupffer,
#   MoMFs, T cells, B cells). For each cell type, uses cell type-specific
#   expressed genes as background. Mouse genes are converted to human
#   orthologs via mygene before running Enrichr.
#
# Input:
#   - Annotated AnnData (.h5ad) with celltype labels
#   - Per-celltype ArchR_P2G_AllClusters_GENES.xlsx files
#     (from compartment_04_p2g_stability_overlap.R / P2G extraction step)
#
# Output (per cell type):
#   - {celltype}_Reactome_Pathways_2024_results.csv (full enrichment results)
#   - {celltype}_dotplot_top5_pathways.pdf
#
# Reference:
#   Sarker et al. (2026) Cell Metabolism
#
# ==============================================================================

import os
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import matplotlib as mpl
import gseapy as gp
import mygene
import warnings

warnings.filterwarnings("ignore")

# Global plot settings
mpl.rcParams["font.family"] = "Arial"
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42


# ==============================================================================
# CONFIGURATION - UPDATE THESE PATHS
# ==============================================================================

# AnnData with celltype annotations
ADATA_PATH = "path/to/final_rna_wnn.h5ad"

# Base directory containing per-celltype P2G output folders
# Expected structure: P2G_BASE_DIR/{celltype}/ArchR_P2G_AllClusters_GENES.xlsx
P2G_BASE_DIR = "path/to/RNA_ATAC_archR/allcelltype"

# Output directory
OUTPUT_DIR = "p2g_pathway_enrichment"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Cell types to process (Hepatocyte excluded - analyzed separately)
CELL_TYPES = [
    "Endothelial_01",
    "Endothelial_02",
    "Stellate",
    "Cholangiocyte_01",
    "Cholangiocyte_02",
    "Kupffer",
    "MoMFs",
    "Tcells",
    "Bcells",
]

# Map P2G folder name to AnnData celltype label (if they differ)
# Update if your celltype labels in adata.obs use different naming
CELLTYPE_TO_ADATA_LABEL = {
    "Endothelial_01": "Endothelial.01",
    "Endothelial_02": "Endothelial.02",
    "Stellate": "Stellate",
    "Cholangiocyte_01": "Cholangiocyte.01",
    "Cholangiocyte_02": "Cholangiocyte.02",
    "Kupffer": "Kupffer",
    "MoMFs": "MoMFs",
    "Tcells": "lymp_T",
    "Bcells": "lymp_B",
}

# Analysis parameters
EXPRESSION_THRESHOLD = 0.05
GENE_SET_LIBRARY = "Reactome_Pathways_2024"

# Pathway filters (for plotting)
PVAL_THRESHOLD = 0.05
MIN_GENE_COUNT = 4
TOP_N_PATHWAYS = 5


# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================

def get_background_genes(adata, celltype_label):
    """Compute expressed background genes for a given celltype."""
    sub = adata[adata.obs["celltype"] == celltype_label].copy()
    if sub.n_obs == 0:
        return []

    expr = sub.X
    avg_expr = np.asarray(expr.mean(axis=0)).flatten()
    expressed = sub.var_names[avg_expr > EXPRESSION_THRESHOLD].tolist()

    print(f"    Cells: {sub.n_obs:,}, Background genes: {len(expressed):,}")
    return expressed


def convert_mouse_to_human(mouse_genes):
    """Convert mouse gene symbols to human orthologs via HomoloGene."""
    mg = mygene.MyGeneInfo()
    mouse_genes = list(set(g for g in mouse_genes if pd.notna(g) and g))

    if len(mouse_genes) == 0:
        return []

    gene_conversion = mg.querymany(
        mouse_genes,
        scopes="symbol",
        fields="homologene",
        species="mouse",
        as_dataframe=False,
        verbose=False,
    )

    human_entrez_ids = set()
    for entry in gene_conversion:
        if isinstance(entry, dict) and "homologene" in entry:
            for tax_id, entrez_id in entry["homologene"].get("genes", []):
                if tax_id == 9606:
                    human_entrez_ids.add(entrez_id)

    if len(human_entrez_ids) == 0:
        return []

    human_conversion = mg.querymany(
        list(human_entrez_ids),
        scopes="entrezgene",
        fields="symbol",
        species="human",
        as_dataframe=False,
        verbose=False,
    )

    human_symbols = sorted({
        entry["symbol"].upper()
        for entry in human_conversion
        if isinstance(entry, dict) and "symbol" in entry
    })

    return human_symbols


def run_enrichr_analysis(gene_list, background, gene_set, title):
    """Run Enrichr with explicit background. Returns full results DataFrame."""
    gene_list = sorted(set(g.upper() for g in gene_list if isinstance(g, str)))
    background = sorted(set(g.upper() for g in background if isinstance(g, str)))
    gene_list = sorted(set(gene_list) & set(background))

    if len(gene_list) < 3:
        print(f"    [SKIP] Too few genes ({len(gene_list)}) for {title}")
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
            print(f"    [SKIP] No Enrichr results for {title}")
            return pd.DataFrame()

        df = enr.results.copy()
        if df is None or df.empty:
            print(f"    [SKIP] No enrichment terms for {title}")
            return pd.DataFrame()

        # Recount overlap with input genes
        input_genes = set(gene_list)

        def recount_overlap(gene_str):
            if not isinstance(gene_str, str):
                return 0
            genes = {g.strip().upper() for g in gene_str.split(";") if g.strip()}
            return len(genes & input_genes)

        df["Gene_Count"] = df["Genes"].apply(recount_overlap)
        df["-log10(pval)"] = -np.log10(df["P-value"].clip(lower=1e-300))
        df["-log10(padj)"] = -np.log10(df["Adjusted P-value"].clip(lower=1e-300))
        df["GeneSet"] = gene_set
        df["N_input_genes"] = len(gene_list)
        df["N_background_genes"] = len(background)

        print(f"    [OK] {title}: {len(df)} terms ({len(gene_list)} input / {len(background)} background)")
        return df

    except Exception as e:
        print(f"    [ERROR] {title}: {e}")
        return pd.DataFrame()


def plot_top_pathways(plot_df, celltype, outfile_prefix):
    """Single-cell-type dot plot of top pathways."""
    if plot_df.empty:
        return

    plot_df = plot_df.copy()
    plot_df["neglog10_adjP"] = -np.log10(plot_df["Adjusted P-value"].clip(lower=1e-300))
    plot_df["Term"] = plot_df["Term"].astype(str)
    plot_df = plot_df.sort_values("Odds Ratio", ascending=False).reset_index(drop=True)

    n_terms = len(plot_df)
    fig_height = max(1.5, 0.20 * n_terms)
    fig, ax = plt.subplots(figsize=(3, fig_height))

    size_scale = 10

    scatter = ax.scatter(
        plot_df["Odds Ratio"], plot_df["Term"],
        s=plot_df["Gene_Count"] * size_scale,
        c=plot_df["neglog10_adjP"],
        cmap="Spectral_r",
        alpha=0.85, edgecolor="black", linewidth=0.5,
        clip_on=False,
    )

    ax.set_xlabel("Odds Ratio", fontsize=9)
    ax.set_ylabel("")
    ax.set_title(f"Top {TOP_N_PATHWAYS} Pathways - {celltype}",
                 fontsize=10, fontweight="bold")
    ax.tick_params(axis="y", labelsize=8)
    ax.tick_params(axis="x", labelsize=8)

    x_min = plot_df["Odds Ratio"].min()
    x_max = plot_df["Odds Ratio"].max()
    x_padding = (x_max - x_min) * 0.15 if x_max > x_min else 0.5
    ax.set_xlim(x_min - x_padding, x_max + x_padding)

    ax.grid(axis="x", linestyle="--", alpha=0.3)
    ax.set_axisbelow(True)
    ax.margins(y=0.1)
    plt.subplots_adjust(right=0.68)

    # Colorbar
    cax = fig.add_axes([0.72, 0.35, 0.025, 0.45])
    cbar = plt.colorbar(scatter, cax=cax)
    cbar.set_label("-log10(Adj. P)", fontsize=8)
    cbar.ax.tick_params(labelsize=7)

    # Size legend
    min_genes = int(plot_df["Gene_Count"].min())
    max_genes = int(plot_df["Gene_Count"].max())

    if max_genes <= 10:
        legend_sizes = [min_genes, max_genes]
    elif max_genes <= 20:
        legend_sizes = [min_genes, 10, max_genes]
    else:
        legend_sizes = [min_genes, (min_genes + max_genes) // 2, max_genes]
    legend_sizes = sorted(set(legend_sizes))

    legend_elements = [
        plt.scatter(
            [], [], s=s * size_scale, c="gray", alpha=0.6,
            edgecolor="black", linewidth=0.5, label=f"{s}",
        )
        for s in legend_sizes
    ]
    ax.legend(
        handles=legend_elements, title="Gene Count",
        loc="upper left", bbox_to_anchor=(1.0, 0.35),
        framealpha=0.9, fontsize=7, title_fontsize=8,
    )

    pdf_file = f"{outfile_prefix}_dotplot_top{TOP_N_PATHWAYS}_pathways.pdf"
    png_file = f"{outfile_prefix}_dotplot_top{TOP_N_PATHWAYS}_pathways.png"
    plt.savefig(pdf_file, dpi=300, bbox_inches="tight")
    plt.savefig(png_file, dpi=300, bbox_inches="tight")
    plt.close(fig)

    print(f"    [OK] Saved: {os.path.basename(pdf_file)}")


# ==============================================================================
# MAIN PIPELINE
# ==============================================================================

print()
print("=" * 70)
print("Loading AnnData")
print("=" * 70)

adata = sc.read_h5ad(ADATA_PATH)
print(f"  Full dataset: {adata.n_obs:,} cells x {adata.n_vars:,} genes")
print(f"  Cell types: {sorted(adata.obs['celltype'].unique().tolist())}")


# Process each cell type
for celltype in CELL_TYPES:

    print()
    print("=" * 70)
    print(f"PROCESSING: {celltype}")
    print("=" * 70)

    adata_label = CELLTYPE_TO_ADATA_LABEL.get(celltype, celltype)

    # ------------------------------------------------------------------
    # Step 1: Background genes
    # ------------------------------------------------------------------
    print(f"\n  Step 1: Computing background genes ({adata_label})")
    expressed_genes = get_background_genes(adata, adata_label)

    if len(expressed_genes) == 0:
        print(f"    [SKIP] No cells found for {adata_label}")
        continue

    # ------------------------------------------------------------------
    # Step 2: Load P2G genes
    # ------------------------------------------------------------------
    print(f"\n  Step 2: Loading P2G genes")
    p2g_file = os.path.join(P2G_BASE_DIR, celltype, "ArchR_P2G_AllClusters_GENES.xlsx")

    if not os.path.exists(p2g_file):
        print(f"    [SKIP] P2G file not found: {p2g_file}")
        continue

    df_p2g = pd.read_excel(p2g_file)
    all_genes = df_p2g["gene"].dropna().astype(str).unique().tolist()
    print(f"    Loaded {len(all_genes):,} unique mouse genes")

    # ------------------------------------------------------------------
    # Step 3: Convert mouse to human
    # ------------------------------------------------------------------
    print(f"\n  Step 3: Converting mouse to human orthologs")
    p2g_genes_human = convert_mouse_to_human(all_genes)
    background_human = convert_mouse_to_human(expressed_genes)
    print(f"    P2G genes:        {len(all_genes):,} mouse -> {len(p2g_genes_human):,} human")
    print(f"    Background genes: {len(expressed_genes):,} mouse -> {len(background_human):,} human")

    # ------------------------------------------------------------------
    # Step 4: Run Enrichr
    # ------------------------------------------------------------------
    print(f"\n  Step 4: Running Reactome enrichment")
    result = run_enrichr_analysis(
        gene_list=p2g_genes_human,
        background=background_human,
        gene_set=GENE_SET_LIBRARY,
        title=celltype,
    )

    if result.empty:
        continue

    # Save full results
    out_csv = os.path.join(OUTPUT_DIR, f"{celltype}_{GENE_SET_LIBRARY}_results.csv")
    result.to_csv(out_csv, index=False)
    print(f"    [OK] Saved: {os.path.basename(out_csv)}")

    # ------------------------------------------------------------------
    # Step 5: Filter and plot
    # ------------------------------------------------------------------
    print(f"\n  Step 5: Filtering and plotting top pathways")

    # Filter: Adjusted P < 0.05, Gene_Count > 4
    filtered = result[
        (result["Adjusted P-value"] < PVAL_THRESHOLD)
        & (result["Gene_Count"] > MIN_GENE_COUNT)
        & (np.isfinite(result["Odds Ratio"]))
    ].copy()

    if filtered.empty:
        print(f"    [SKIP] No terms pass filter (Adj P < {PVAL_THRESHOLD}, Gene_Count > {MIN_GENE_COUNT})")
        continue

    # Top N by Odds Ratio
    top_pathways = filtered.sort_values("Odds Ratio", ascending=False).head(TOP_N_PATHWAYS)
    print(f"    Selected top {len(top_pathways)} pathways")

    plot_top_pathways(
        top_pathways,
        celltype=celltype,
        outfile_prefix=os.path.join(OUTPUT_DIR, celltype),
    )


# ==============================================================================
# SUMMARY
# ==============================================================================

print()
print("=" * 70)
print("  COMPLETE")
print("=" * 70)
print(f"  Cell types processed: {len(CELL_TYPES)}")
print(f"  Filters: Adj P < {PVAL_THRESHOLD}, Gene_Count > {MIN_GENE_COUNT}, top {TOP_N_PATHWAYS} by Odds Ratio")
print(f"  Output: {OUTPUT_DIR}/")
print("    {celltype}_Reactome_Pathways_2024_results.csv (full results)")
print(f"    {{celltype}}_dotplot_top{TOP_N_PATHWAYS}_pathways.pdf/png")
print("=" * 70)
