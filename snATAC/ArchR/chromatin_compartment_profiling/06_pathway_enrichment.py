#!/usr/bin/env python3
# ==============================================================================
# Pathway Enrichment of P2G-linked genes by compartment stability class
# ==============================================================================
#
# Description:
#   Runs Reactome pathway enrichment (via Enrichr) on P2G-linked genes
#   stratified by sex and compartment stability class. Uses Hepatocyte-
#   expressed genes as background. Mouse genes are converted to human
#   orthologs via mygene before enrichment. Results are visualized as
#   dot matrix plots (rows = pathways, columns = male/female, dot color
#   = stability class, dot size = -log10 P-value, background = odds ratio).
#
# Input:
#   - genes_[sex]_[stability].txt files
#     (from compartment_04_p2g_stability_overlap.R)
#   - Final annotated AnnData (.h5ad) with celltype labels
#
#
# ==============================================================================

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib as mpl
from matplotlib.patches import Rectangle
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
import scanpy as sc
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

# Base directory for stability analysis outputs
BASE_OUTDIR = "path/to/output/downstream_stability"

# P2G cutoff versions
P2G_VERSIONS = {
    "default": os.path.join(BASE_OUTDIR, "default", "peak2gene_stability_analysis"),
    "relaxed": os.path.join(BASE_OUTDIR, "relaxed", "peak2gene_stability_analysis"),
}

# Which version to run enrichment on
SELECTED_VERSION = "relaxed"

# Analysis parameters
STABILITY_CLASSES = [
    "Stable_Active",
    "Stable_Repressive",
    "Monotonic_A_to_R",
    "Monotonic_R_to_A",
    "Non_Monotonic",
]
SEXES = ["male", "female"]

# Enrichment parameters
GENE_SET_LIBRARIES = ["Reactome_Pathways_2024"]
PVALUE_THRESHOLD = 0.05
MIN_GENE_COUNT = 2
TOP_N_TERMS = 10
EXPRESSION_THRESHOLD = 0.05

# Visualization constants
STABILITY_COLS = {
    "Stable_Active": "#FF6B6B",
    "Stable_Repressive": "#4E79A7",
    "Monotonic_A_to_R": "darkred",
    "Monotonic_R_to_A": "darkblue",
    "Non_Monotonic": "#708238",
}

STABILITY_LABELS = {
    "Stable_Active": "Stable Active",
    "Stable_Repressive": "Stable Repressive",
    "Monotonic_A_to_R": "A to R",
    "Monotonic_R_to_A": "R to A",
    "Non_Monotonic": "Non-Monotonic",
}


# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================

def convert_mouse_to_human(mouse_genes):
    """Convert mouse gene symbols to human orthologs via mygene."""
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
            genes = entry["homologene"].get("genes", [])
            for tax_id, entrez_id in genes:
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


def run_enrichr_analysis(gene_list, background, gene_set, title, outfile_prefix):
    """Run Enrichr with explicit background. Returns full results DataFrame."""
    gene_list = sorted(set(g.upper() for g in gene_list if isinstance(g, str)))
    background = sorted(set(g.upper() for g in background if isinstance(g, str)))
    gene_list = sorted(set(gene_list) & set(background))

    if len(gene_list) < 3:
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
            print(f"  [SKIP] Enrichr returned no results for {title}")
            return pd.DataFrame()

        df = enr.results.copy()
        if df is None or df.empty:
            print(f"  [SKIP] No enrichment terms for {title}")
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

        out_csv = f"{outfile_prefix}_results.csv"
        df.to_csv(out_csv, index=False)

        print(f"  [OK] {title}: {len(df)} terms | {len(gene_list)} input | {len(background)} background")
        return df

    except Exception as e:
        print(f"  [ERROR] {title}: {e}")
        return pd.DataFrame()


def plot_dot_matrix(plot_df, stability_cols, outfile_prefix,
                    title="Reactome Pathways", figsize_width=4,
                    size_scale=30, min_size=20):
    """
    Dot matrix plot with multiple dots per cell when a pathway appears
    in multiple stability classes.

    Rows = pathways, columns = male/female.
    Background color = max odds ratio for that cell.
    Dot color = stability class, dot size = -log10(P-value).
    """
    df = plot_df.copy()
    df["neglog10_pval"] = -np.log10(df["P-value"].clip(lower=1e-300))

    term_order = (
        df.groupby("Term")["Odds Ratio"]
        .max()
        .sort_values(ascending=False)
        .index.tolist()
    )

    term_to_y = {term: i for i, term in enumerate(term_order)}
    sex_to_x = {"male": 0, "female": 1}

    n_terms = len(term_order)
    fig_height = max(2, 0.25 * n_terms)
    fig, ax = plt.subplots(figsize=(figsize_width, fig_height))

    or_max = df["Odds Ratio"].max()
    norm = Normalize(vmin=0, vmax=or_max)
    cmap = plt.cm.Blues

    max_or_per_cell = df.groupby(["Term", "Sex"])["Odds Ratio"].max().to_dict()
    stab_counts = df.groupby(["Term", "Sex"]).size().to_dict()

    # Background rectangles
    for term in term_order:
        for sex in ["male", "female"]:
            x = sex_to_x[sex]
            y = term_to_y[term]
            or_val = max_or_per_cell.get((term, sex), 0)
            facecolor = cmap(norm(or_val)) if or_val > 0 else "#f5f5f5"

            rect = Rectangle(
                (x - 0.48, y - 0.48), 0.96, 0.96,
                facecolor=facecolor, edgecolor="white",
                linewidth=1.5, zorder=1,
            )
            ax.add_patch(rect)

    # Dot offsets for multiple stability classes per cell
    def get_dot_offsets(n):
        if n == 1:
            return [(0, 0)]
        elif n == 2:
            return [(-0.15, 0), (0.15, 0)]
        elif n == 3:
            return [(-0.2, 0), (0, 0), (0.2, 0)]
        elif n == 4:
            return [(-0.15, 0.12), (0.15, 0.12), (-0.15, -0.12), (0.15, -0.12)]
        else:
            return [(-0.2, 0.12), (0, 0.12), (0.2, 0.12),
                    (-0.15, -0.12), (0.15, -0.12)]

    # Draw dots
    for (term, sex), group in df.groupby(["Term", "Sex"]):
        x_base = sex_to_x[sex]
        y_base = term_to_y[term]
        n_dots = len(group)
        offsets = get_dot_offsets(n_dots)
        group_sorted = group.sort_values("Odds Ratio", ascending=True)

        for i, (_, row) in enumerate(group_sorted.iterrows()):
            stab = row["Stability"]
            pval = row["neglog10_pval"]
            x_off, y_off = offsets[i] if i < len(offsets) else (0, 0)
            max_sz = 300 if n_dots == 1 else (200 if n_dots <= 2 else 150)
            size = min(pval * size_scale + min_size, max_sz)

            ax.scatter(
                x_base + x_off, y_base + y_off,
                s=size, c=stability_cols.get(stab, "gray"),
                edgecolor="black", linewidth=0.8,
                alpha=0.95, zorder=5, clip_on=False,
            )

    # Axes
    ax.set_xlim(-0.6, 1.6)
    ax.set_ylim(-0.6, n_terms - 0.4)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["MALE", "FEMALE"], fontsize=11, fontweight="bold")
    ax.xaxis.set_ticks_position("top")
    ax.xaxis.set_label_position("top")
    ax.set_yticks(range(n_terms))
    ax.set_yticklabels(term_order, fontsize=8)
    ax.set_title(f"{title}\n", fontsize=9, fontweight="bold", pad=5)
    ax.invert_yaxis()

    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.tick_params(axis="both", length=0)

    plt.subplots_adjust(right=0.72)

    # Colorbar (odds ratio)
    sm = ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar_ax = fig.add_axes([0.75, 0.60, 0.02, 0.25])
    cbar = plt.colorbar(sm, cax=cbar_ax)
    cbar.set_label("Odds Ratio", fontsize=9, fontweight="bold")
    cbar.ax.tick_params(labelsize=8)
    cbar.outline.set_linewidth(0.5)

    # Stability legend
    unique_stab = df["Stability"].unique()
    color_patches = [
        mpatches.Patch(
            color=stability_cols.get(s, "gray"),
            label=STABILITY_LABELS.get(s, s),
            edgecolor="black", linewidth=0.5,
        )
        for s in STABILITY_COLS if s in unique_stab
    ]
    fig.legend(
        handles=color_patches, title="Stability",
        loc="center left", bbox_to_anchor=(0.74, 0.38),
        framealpha=0.95, fontsize=8, title_fontsize=9, edgecolor="gray",
    )

    # Size legend
    min_pval = df["neglog10_pval"].min()
    max_pval = df["neglog10_pval"].max()
    if max_pval - min_pval < 2:
        legend_vals = [round(min_pval, 1), round(max_pval, 1)]
    else:
        mid_pval = (min_pval + max_pval) / 2
        legend_vals = [round(min_pval, 1), round(mid_pval, 1), round(max_pval, 1)]
    legend_vals = sorted(set(legend_vals))

    size_elements = [
        plt.scatter(
            [], [], s=min(v * size_scale + min_size, 300),
            c="white", edgecolor="black", linewidth=0.8, label=f"{v:.1f}",
        )
        for v in legend_vals
    ]
    fig.legend(
        handles=size_elements, title="-log10(P)",
        loc="center left", bbox_to_anchor=(0.74, 0.12),
        framealpha=0.95, fontsize=8, title_fontsize=10, edgecolor="gray",
    )

    # Save
    pdf_file = f"{outfile_prefix}_dotmatrix.pdf"
    png_file = f"{outfile_prefix}_dotmatrix.png"
    plt.savefig(pdf_file, dpi=300, bbox_inches="tight", facecolor="white")
    plt.savefig(png_file, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)

    print(f"  [OK] Saved: {os.path.basename(pdf_file)}")
    print(f"  [OK] Terms: {n_terms}")
    print(f"  [OK] Odds Ratio range: 0 - {or_max:.1f}")
    print(f"  [OK] -log10(P) range: {min_pval:.1f} - {max_pval:.1f}")

    # Report multi-dot cells
    multi_dot_cells = {k: v for k, v in stab_counts.items() if v > 1}
    if multi_dot_cells:
        print(f"  [OK] Cells with multiple dots: {len(multi_dot_cells)}")
        for (term, sex), count in sorted(multi_dot_cells.items(), key=lambda x: -x[1]):
            print(f"       {sex}: {term} ({count} stability classes)")


# ==============================================================================
# STEP 1: DEFINE BACKGROUND GENES FROM HEPATOCYTE EXPRESSION
# ==============================================================================

print()
print("=" * 70)
print("STEP 1: Define background genes from Hepatocyte expression")
print("=" * 70)

adata = sc.read_h5ad(ADATA_PATH)
print(f"  Full adata: {adata.shape[0]} cells x {adata.shape[1]} genes")

adata = adata[adata.obs["celltype"] == "Hepatocyte"].copy()
print(f"  Hepatocyte adata: {adata.shape[0]} cells x {adata.shape[1]} genes")

expr = adata.X
gene_names = adata.var_names
avg_expr = np.asarray(expr.mean(axis=0)).flatten()

expressed_genes = gene_names[avg_expr > EXPRESSION_THRESHOLD].tolist()
print(f"  Background genes (mean expr > {EXPRESSION_THRESHOLD}): {len(expressed_genes)}")

# Free memory
del adata, expr


# ==============================================================================
# STEP 2: LOAD P2G GENE LISTS
# ==============================================================================

print()
print("=" * 70)
print("STEP 2: Load P2G gene lists")
print("=" * 70)

genes_by_version = {}

for version, p2g_outdir in P2G_VERSIONS.items():
    cutoff = "0.45" if version == "default" else "0.25"
    print(f"\n  {version.upper()} (corCutOff = {cutoff})")

    genes_by_version[version] = {}

    for sex in SEXES:
        genes_by_version[version][sex] = {}
        for stab_class in STABILITY_CLASSES:
            gene_file = os.path.join(p2g_outdir, f"genes_{sex}_{stab_class}.txt")
            if os.path.exists(gene_file):
                with open(gene_file, "r") as f:
                    genes = [l.strip() for l in f if l.strip()]
                genes_by_version[version][sex][stab_class] = genes
                print(f"    {sex:6} | {stab_class:20}: {len(genes):5} genes")
            else:
                genes_by_version[version][sex][stab_class] = []
                print(f"    {sex:6} | {stab_class:20}: FILE NOT FOUND")

genes_by_sex_class = genes_by_version[SELECTED_VERSION]
print(f"\n  [OK] Using {SELECTED_VERSION.upper()} version for enrichment")


# ==============================================================================
# STEP 3: CONVERT MOUSE GENES TO HUMAN ORTHOLOGS
# ==============================================================================

print()
print("=" * 70)
print("STEP 3: Convert mouse genes to human orthologs")
print("=" * 70)

human_genes_by_sex_class = {}

for sex in SEXES:
    print(f"\n  {sex.upper()}:")
    human_genes_by_sex_class[sex] = {}

    for stab_class in STABILITY_CLASSES:
        mouse_genes = genes_by_sex_class[sex][stab_class]

        if len(mouse_genes) == 0:
            human_genes_by_sex_class[sex][stab_class] = []
            print(f"    {stab_class:20}: 0 mouse -> 0 human (skipped)")
            continue

        human_genes = convert_mouse_to_human(mouse_genes)
        human_genes_by_sex_class[sex][stab_class] = human_genes

        rate = 100 * len(human_genes) / len(mouse_genes) if mouse_genes else 0
        print(f"    {stab_class:20}: {len(mouse_genes):5} mouse -> {len(human_genes):5} human ({rate:.1f}%)")


# ==============================================================================
# STEP 4: RUN ENRICHR PATHWAY ANALYSIS
# ==============================================================================

print()
print("=" * 70)
print(f"STEP 4: Run pathway analysis ({SELECTED_VERSION.upper()})")
print("=" * 70)

# Background = union of all P2G human genes across stability classes
background_human = set()
for sex in SEXES:
    for stab_class in STABILITY_CLASSES:
        background_human.update(human_genes_by_sex_class[sex][stab_class])
background_human = sorted(background_human)

print(f"  Background genes (human): {len(background_human)}")

all_results = {}

for sex in SEXES:
    all_results[sex] = {}

    for stab_class in STABILITY_CLASSES:
        p2g_genes_human = human_genes_by_sex_class[sex][stab_class]

        print(f"\n  {sex.upper()} - {stab_class} ({len(p2g_genes_human)} genes)")

        if len(p2g_genes_human) < 3:
            print(f"  [SKIP] Too few genes")
            all_results[sex][stab_class] = {}
            continue

        all_results[sex][stab_class] = {}

        for gs in GENE_SET_LIBRARIES:
            result = run_enrichr_analysis(
                gene_list=p2g_genes_human,
                background=background_human,
                gene_set=gs,
                title=f"{sex}_{stab_class}",
                outfile_prefix=f"{SELECTED_VERSION}_{sex}_{stab_class}_{gs}",
            )

            if result is not None and not result.empty:
                all_results[sex][stab_class][gs] = result


# ==============================================================================
# STEP 5: FILTER AND PREPARE DATA FOR PLOTTING
# ==============================================================================

print()
print("=" * 70)
print("STEP 5: Filter and prepare data for plotting")
print("=" * 70)

combined_data = []

for sex in SEXES:
    for stab_class in STABILITY_CLASSES:
        if stab_class not in all_results[sex]:
            continue

        for gs, df in all_results[sex][stab_class].items():
            if df is None or df.empty:
                continue

            df_filtered = df[
                (df["P-value"] < PVALUE_THRESHOLD)
                & (df["Gene_Count"] > MIN_GENE_COUNT)
                & (np.isfinite(df["Odds Ratio"]))
            ].copy()

            if df_filtered.empty:
                print(f"  {sex:6} | {stab_class:20}: No terms pass filter")
                continue

            df_top = (
                df_filtered
                .sort_values("Odds Ratio", ascending=False)
                .head(TOP_N_TERMS)
                .copy()
            )

            df_top["Sex"] = sex
            df_top["Stability"] = stab_class
            df_top["GeneSet"] = gs

            combined_data.append(df_top)
            print(f"  {sex:6} | {stab_class:20}: {len(df_top)} terms")

if combined_data:
    plot_df = pd.concat(combined_data, ignore_index=True)
    print(f"\n  [OK] Total terms for plotting: {len(plot_df)}")
else:
    print("  [SKIP] No data to plot")
    plot_df = pd.DataFrame()


# ==============================================================================
# STEP 6: GENERATE DOT MATRIX PLOT
# ==============================================================================

if not plot_df.empty:
    print()
    print("=" * 70)
    print(f"STEP 6: Generate dot matrix plot ({SELECTED_VERSION.upper()})")
    print("=" * 70)

    print("\n  Terms per Sex x Stability:")
    summary = plot_df.groupby(["Sex", "Stability"]).size().reset_index(name="N_terms")
    print(summary.to_string(index=False))

    plot_dot_matrix(
        plot_df=plot_df,
        stability_cols=STABILITY_COLS,
        outfile_prefix=f"Reactome_Hepatocyte_Compartment_{SELECTED_VERSION}",
        title=(
            f"Reactome Pathways ({SELECTED_VERSION.upper()}), "
            f"P < {PVALUE_THRESHOLD}, Gene_Count > {MIN_GENE_COUNT}"
        ),
        figsize_width=4,
        size_scale=30,
        min_size=20,
    )


# ==============================================================================
# SUMMARY
# ==============================================================================

print()
print("=" * 70)
print("  COMPLETE")
print("=" * 70)
print()
print(f"  Version: {SELECTED_VERSION.upper()} (corCutOff = "
      f"{'0.45' if SELECTED_VERSION == 'default' else '0.25'})")
print(f"  Background: {len(background_human)} human genes")
print(f"  Filter: P < {PVALUE_THRESHOLD}, Gene_Count > {MIN_GENE_COUNT}")
print(f"  Top terms per group: {TOP_N_TERMS}")
print()
print("  Output files:")
print(f"    Reactome_Hepatocyte_Compartment_{SELECTED_VERSION}_dotmatrix.pdf/png")
print(f"    {SELECTED_VERSION}_[sex]_[stability]_Reactome_Pathways_2024_results.csv")
print("=" * 70)
