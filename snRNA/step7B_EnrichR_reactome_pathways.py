#!/usr/bin/env python
"""
reactome_from_heatmap_excel.py
===============================
Reads Heatmap_genes_stat40_with_clusters.xlsx (each tab = celltype,
each row has gene + cluster column), converts mouse → human,
runs Enrichr (Reactome_Pathways_2024) per celltype × cluster,
and plots cluster-wise bubble plots.

Input:  Heatmap_genes_stat40_with_clusters.xlsx
Output: Per-celltype bubble plot PDFs + CSVs
"""

import os
import numpy as np
import pandas as pd
import scanpy as sc
import gseapy as gp
import mygene
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# ============================================================
# Configuration
# ============================================================
EXCEL_PATH   = "Heatmap_genes_stat40_with_clusters.xlsx"
ADATA_PATH   = "adata_with_resolutions2.h5ad"
FIGURES_DIR  = "figures"
RESULTS_DIR  = "results/enrichr"

os.makedirs(FIGURES_DIR, exist_ok=True)
os.makedirs(RESULTS_DIR, exist_ok=True)

PVAL_CUTOFF  = 0.05
TOP_N        = 20
EXPR_THRESH  = 0.05

CLUSTER_PALETTE = {
    1: "#2C3E50",   # dark
    2: "#95A5A6",   # gray
    3: "#E74C3C",   # red
    4: "#3498DB",   # blue
}


# ============================================================
# 1. Build background gene list from adata
# ============================================================
print("=" * 60)
print("STEP 1: Building background gene list")
print("=" * 60)

adata = sc.read_h5ad(ADATA_PATH)
adata = adata[adata.obs["celltype"] == "Hepatocyte"].copy()

avg_expr = np.asarray(adata.X.mean(axis=0)).flatten()
expressed_genes = [g.upper() for g in adata.var_names[avg_expr > EXPR_THRESH]]
print(f"  Background genes (avg > {EXPR_THRESH}): {len(expressed_genes)}")


# ============================================================
# 2. Mouse → Human converter
# ============================================================
print("\n" + "=" * 60)
print("STEP 2: Setting up mouse → human converter")
print("=" * 60)

mg = mygene.MyGeneInfo()


def mouse_to_human(mouse_gene_list):
    if not mouse_gene_list:
        return []
    results = mg.querymany(
        mouse_gene_list, scopes="symbol", fields="homologene",
        species="mouse", verbose=False
    )
    human_entrez = []
    for entry in results:
        if "homologene" in entry and "genes" in entry["homologene"]:
            for homolog in entry["homologene"]["genes"]:
                if homolog[0] == 9606:
                    human_entrez.append(homolog[1])
    if not human_entrez:
        return []
    sym_results = mg.querymany(
        human_entrez, scopes="entrezgene", fields="symbol",
        species="human", verbose=False
    )
    return list(set(e["symbol"].upper() for e in sym_results if "symbol" in e))


# ============================================================
# 3. Read Excel + run Enrichr per celltype × cluster
# ============================================================
print("\n" + "=" * 60)
print("STEP 3: Reading Excel + running Enrichr")
print("=" * 60)

xl = pd.ExcelFile(EXCEL_PATH)
all_enrich = {}

for sheet in xl.sheet_names:

    ct = sheet  # e.g., "Hep-01"
    ct_safe = ct.replace("-", "_")
    df = pd.read_excel(xl, sheet_name=sheet)

    if "gene" not in df.columns or "cluster" not in df.columns:
        print(f"\n  ⚠️ Skipping {ct}: missing 'gene' or 'cluster' column")
        continue

    clusters = sorted(df["cluster"].dropna().unique())
    print(f"\n📊 {ct}: {len(df)} genes, clusters = {clusters}")

    ct_results = {}

    for cluster_id in clusters:
        cluster_id = int(cluster_id)
        cluster_df = df[df["cluster"] == cluster_id]
        mouse_genes = cluster_df["gene"].dropna().tolist()
        print(f"    Cluster {cluster_id}: {len(mouse_genes)} mouse genes")

        if len(mouse_genes) < 3:
            print(f"    ⚠️ Too few genes, skipping")
            continue

        # Mouse → Human
        human_genes = mouse_to_human(mouse_genes)
        print(f"      → {len(human_genes)} human orthologs")

        if len(human_genes) < 3:
            print(f"    ⚠️ Too few orthologs, skipping")
            continue

        # Enrichr
        try:
            enrich = gp.enrichr(
                gene_list=human_genes,
                gene_sets="Reactome_Pathways_2024",
                organism="Human",
                background=expressed_genes,
                outdir=None,
            )
        except Exception as e:
            print(f"    ❌ Enrichr failed: {e}")
            continue

        enrich_df = pd.DataFrame(enrich.results)
        if enrich_df.empty:
            print(f"    ⚠️ No results")
            continue

        # Recount gene overlap
        input_set = set(human_genes)

        def recount_overlap(gene_str):
            if not isinstance(gene_str, str) or gene_str.strip() == "":
                return 0
            genes = {g.strip().upper() for g in gene_str.split(";") if g.strip()}
            return len(genes & input_set)

        enrich_df["Gene Count"] = enrich_df["Genes"].apply(recount_overlap)
        enrich_df["-log10(p-value)"] = -np.log10(enrich_df["P-value"].clip(lower=1e-300))
        enrich_df["cluster"] = cluster_id
        enrich_df["celltype"] = ct

        # Filter significant
        enrich_sig = enrich_df[enrich_df["Adjusted P-value"] < PVAL_CUTOFF].copy()

        # Save CSV
        out_csv = f"{RESULTS_DIR}/{ct_safe}_cluster_{cluster_id}_Reactome.csv"
        enrich_df.to_csv(out_csv, index=False)
        print(f"      Saved: {out_csv} ({len(enrich_sig)} sig pathways)")

        ct_results[cluster_id] = enrich_sig

    if ct_results:
        all_enrich[ct] = ct_results


# ============================================================
# 4. Per-celltype bubble plots
# ============================================================
print("\n" + "=" * 60)
print("STEP 4: Generating bubble plots")
print("=" * 60)

for ct, cluster_dict in all_enrich.items():

    ct_safe = ct.replace("-", "_")

    # Combine top pathways per cluster
    plot_frames = []
    for cid, edf in sorted(cluster_dict.items()):
        if edf.empty:
            continue
        top = edf.sort_values("Adjusted P-value").head(TOP_N).copy()
        top["cluster"] = cid
        plot_frames.append(top)

    if not plot_frames:
        print(f"  ⚠️ No significant pathways for {ct}")
        continue

    plot_df = pd.concat(plot_frames, ignore_index=True)

    # Truncate long names
    MAX_CHARS = 55
    plot_df["Term_short"] = plot_df["Term"].apply(
        lambda x: x[:MAX_CHARS] + "…" if len(str(x)) > MAX_CHARS else x
    )

    # Deduplicate labels
    seen = {}
    unique_labels = []
    for _, row in plot_df.iterrows():
        label = row["Term_short"]
        if label in seen:
            seen[label] += 1
            label = label + " " * seen[label]
        else:
            seen[label] = 0
        unique_labels.append(label)
    plot_df["y_label"] = unique_labels

    # X-axis column
    x_col = "Odds Ratio" if "Odds Ratio" in plot_df.columns else (
        "Combined Score" if "Combined Score" in plot_df.columns else "-log10(p-value)")

    plot_df = plot_df.sort_values(
        ["cluster", x_col], ascending=[True, False]
    ).reset_index(drop=True)

    # ---- Plot ----
    fig, ax = plt.subplots(figsize=(10, max(6, 0.35 * len(plot_df))))

    # Lollipop stems
    for _, row in plot_df.iterrows():
        ax.hlines(
            y=row["y_label"], xmin=0, xmax=row[x_col],
            color="gray", linewidth=1.2, linestyle="dotted", zorder=1,
        )

    # Bubbles per cluster
    clusters_sorted = sorted(plot_df["cluster"].unique())
    for cid in clusters_sorted:
        sub = plot_df[plot_df["cluster"] == cid]
        color = CLUSTER_PALETTE.get(cid, "#333333")
        ax.scatter(
            sub[x_col], sub["y_label"],
            s=sub["Gene Count"] * 20,
            c=color, edgecolors="black", linewidths=0.5,
            alpha=0.85, zorder=3,
            label=f"cluster {cid}",
        )

    # Cluster separator line
    if len(clusters_sorted) > 1:
        for i in range(len(clusters_sorted) - 1):
            last_idx = plot_df[plot_df["cluster"] == clusters_sorted[i]].index.max()
            if not pd.isna(last_idx):
                ax.axhline(
                    y=last_idx + 0.5, color="lightgray",
                    linewidth=1.2, linestyle="-", zorder=0,
                )

    # Formatting
    ax.set_xlabel(x_col, fontsize=12)
    ax.set_ylabel("")
    ax.set_title(
        f"{ct}\nREACTOME pathways padj < {PVAL_CUTOFF}",
        fontsize=14, fontweight="bold",
    )
    ax.tick_params(axis="y", labelsize=9)
    ax.set_ylim(len(plot_df) - 0.5, -0.5)
    sns.despine(ax=ax, left=True)

    # --- Legends ---
    # Cluster legend
    cluster_handles = [
        Line2D([0], [0], marker="o", color="w",
               markerfacecolor=CLUSTER_PALETTE.get(cid, "#333"),
               markersize=10, markeredgecolor="black", markeredgewidth=0.5,
               label=f"cluster {cid}")
        for cid in clusters_sorted
    ]

    # Gene count size legend
    gc_values = [1, 5, 10, 15]
    size_handles = [
        Line2D([0], [0], marker="o", color="w",
               markerfacecolor="gray", markeredgecolor="black",
               markeredgewidth=0.5, markersize=np.sqrt(v * 20) * 0.8,
               label=str(v))
        for v in gc_values
    ]

    leg1 = ax.legend(
        handles=cluster_handles, title="",
        bbox_to_anchor=(1.02, 1), loc="upper left", frameon=False,
    )
    ax.add_artist(leg1)
    ax.legend(
        handles=size_handles, title="gene count",
        bbox_to_anchor=(1.02, 0.6), loc="upper left", frameon=False,
    )

    plt.tight_layout()

    for ext in ["pdf", "png"]:
        fname = f"{FIGURES_DIR}/{ct_safe}_Reactome_clusters.{ext}"
        fig.savefig(fname, bbox_inches="tight", dpi=300)
        print(f"  ✅ {fname}")
    plt.close(fig)


# ============================================================
# Done
# ============================================================
print("\n" + "=" * 60)
print("REACTOME ENRICHMENT COMPLETE")
print(f"  Input:  {EXCEL_PATH}")
print(f"  CSVs:   {RESULTS_DIR}/")
print(f"  Plots:  {FIGURES_DIR}/*_Reactome_clusters.pdf")
print("=" * 60)