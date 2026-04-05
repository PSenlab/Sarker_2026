#!/usr/bin/env python3
# ==============================================================================
# UpSet plots for P2G-linked gene overlaps across stability classes
# ==============================================================================
#
# Description:
#   Generates UpSet plots showing gene overlaps across the 5 compartment
#   stability classes for Hepatocyte P2G-linked genes. Runs for both
#   default (corCutOff=0.45) and relaxed (corCutOff=0.25) P2G cutoffs,
#   and produces a comparison table between the two.
#
# Input:
#   - genes_[sex]_[stability].txt files
#     (from compartment_04_p2g_stability_overlap.R)
#
# Output (per cutoff):
#   - UpSet_Male_Stability_[version].pdf/png
#   - UpSet_Female_Stability_[version].pdf/png
#   - UpSet_MaleFemale_Switching_[version].pdf/png
#   - UpSet_All_SexStability_[version].pdf/png
#   - male_female_overlap_summary_[version].tsv
#   - comparison_default_vs_relaxed.tsv
#
#
# ==============================================================================

import os
import re
import inspect
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from upsetplot import UpSet, from_contents
import warnings

warnings.filterwarnings("ignore")

# Global plot settings
mpl.rcParams["font.family"] = "Arial"
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42


# ==============================================================================
# CONFIGURATION - UPDATE THIS PATH
# ==============================================================================

base_outdir = "path/to/output/downstream_stability"

P2G_VERSIONS = {
    "default": os.path.join(base_outdir, "default", "peak2gene_stability_analysis"),
    "relaxed": os.path.join(base_outdir, "relaxed", "peak2gene_stability_analysis"),
}

STABILITY_CLASSES = [
    "Stable_Active",
    "Stable_Repressive",
    "Monotonic_A_to_R",
    "Monotonic_R_to_A",
    "Non_Monotonic",
]

SWITCHING_CLASSES = ["Monotonic_A_to_R", "Monotonic_R_to_A", "Non_Monotonic"]

SEXES = ["female", "male"]

STAB_LABEL = {
    "Stable_Active": "Stable Active",
    "Stable_Repressive": "Stable Repressive",
    "Monotonic_A_to_R": "Monotonic A to R",
    "Monotonic_R_to_A": "Monotonic R to A",
    "Non_Monotonic": "Non-Monotonic",
}

SHORT_LABEL = {
    ("female", "Stable_Active"): "F_SA",
    ("female", "Stable_Repressive"): "F_SR",
    ("female", "Monotonic_A_to_R"): "F_AtoR",
    ("female", "Monotonic_R_to_A"): "F_RtoA",
    ("female", "Non_Monotonic"): "F_NM",
    ("male", "Stable_Active"): "M_SA",
    ("male", "Stable_Repressive"): "M_SR",
    ("male", "Monotonic_A_to_R"): "M_AtoR",
    ("male", "Monotonic_R_to_A"): "M_RtoA",
    ("male", "Non_Monotonic"): "M_NM",
}

DEFAULT_MIN_SUBSET_SIZE = 1
DEFAULT_MAX_INTERSECTIONS = 30


# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================

def read_gene_file(path):
    """Read a gene list file (one gene per line)."""
    genes = set()
    with open(path, "r") as f:
        for line in f:
            g = line.strip()
            if not g:
                continue
            g = re.sub(r"\s+", "", g)
            genes.add(g)
    return genes


def create_upset_plot(gene_dict, title, outfile_prefix, output_dir,
                      min_subset_size=1, max_intersections=30):
    """
    Create UpSet plot from dictionary {label: set(genes)}.
    Compatible across upsetplot versions.
    """
    gene_dict_filtered = {
        k: v for k, v in gene_dict.items() if isinstance(v, set) and len(v) > 0
    }

    if len(gene_dict_filtered) < 2:
        print(f"  [SKIP] Not enough non-empty sets for: {title}")
        return

    upset_data = from_contents(gene_dict_filtered)
    fig = plt.figure(figsize=(13, 8))

    kwargs = dict(
        subset_size="count",
        show_counts=True,
        show_percentages=False,
        sort_by="cardinality",
        sort_categories_by="cardinality",
        min_subset_size=min_subset_size,
        element_size=38,
    )

    # Handle version differences in upsetplot API
    sig = inspect.signature(UpSet.__init__).parameters
    if "max_intersections" in sig:
        kwargs["max_intersections"] = max_intersections
    elif "max_rank" in sig:
        kwargs["max_rank"] = max_intersections

    upset = UpSet(upset_data, **kwargs)
    upset.plot(fig=fig)
    plt.suptitle(title, fontsize=14, fontweight="bold", y=1.02)

    pdf_file = os.path.join(output_dir, f"{outfile_prefix}_upset.pdf")
    png_file = os.path.join(output_dir, f"{outfile_prefix}_upset.png")
    plt.savefig(pdf_file, dpi=300, bbox_inches="tight")
    plt.savefig(png_file, dpi=300, bbox_inches="tight")
    plt.close(fig)

    print(f"  [OK] Saved: {os.path.basename(pdf_file)}")
    print(f"  [OK] Saved: {os.path.basename(png_file)}")


# ==============================================================================
# LOAD GENE LISTS
# ==============================================================================

print()
print("=" * 70)
print("LOADING GENE LISTS FROM DEFAULT AND RELAXED VERSIONS")
print("=" * 70)

genes_by_version = {}

for version, p2g_outdir in P2G_VERSIONS.items():
    cutoff = "0.45" if version == "default" else "0.25"
    print(f"\n  {version.upper()} (corCutOff = {cutoff})")
    print(f"  Path: {p2g_outdir}")

    os.makedirs(p2g_outdir, exist_ok=True)
    genes_by_version[version] = {sex: {} for sex in SEXES}
    missing = []

    for sex in SEXES:
        for stab in STABILITY_CLASSES:
            gene_file = os.path.join(p2g_outdir, f"genes_{sex}_{stab}.txt")
            if os.path.exists(gene_file):
                genes = read_gene_file(gene_file)
                genes_by_version[version][sex][stab] = genes
                print(f"    {sex:6} | {stab:20}: {len(genes):5} genes")
            else:
                genes_by_version[version][sex][stab] = set()
                missing.append(gene_file)
                print(f"    {sex:6} | {stab:20}: FILE NOT FOUND")

    if missing:
        print(f"\n  [NOTE] {len(missing)} files not found for {version}")


# ==============================================================================
# GENERATE UPSET PLOTS
# ==============================================================================

for version, p2g_outdir in P2G_VERSIONS.items():

    print()
    print("=" * 70)
    print(f"GENERATING UPSET PLOTS: {version.upper()}")
    print("=" * 70)

    genes_by_sex_class = genes_by_version[version]

    # Plot 1: Male - all stability classes
    print(f"\n  Plot 1: Male ({version})")
    male_genes = {
        STAB_LABEL[stab]: genes_by_sex_class["male"][stab]
        for stab in STABILITY_CLASSES
    }
    create_upset_plot(
        male_genes,
        f"Gene overlaps across stability classes (Male) - {version.upper()}",
        f"UpSet_Male_Stability_{version}",
        p2g_outdir,
    )

    # Plot 2: Female - all stability classes
    print(f"\n  Plot 2: Female ({version})")
    female_genes = {
        STAB_LABEL[stab]: genes_by_sex_class["female"][stab]
        for stab in STABILITY_CLASSES
    }
    create_upset_plot(
        female_genes,
        f"Gene overlaps across stability classes (Female) - {version.upper()}",
        f"UpSet_Female_Stability_{version}",
        p2g_outdir,
    )

    # Plot 3: Switching classes - male vs female
    print(f"\n  Plot 3: Switching Male vs Female ({version})")
    switching_genes = {}
    for sex in SEXES:
        for stab in SWITCHING_CLASSES:
            label = f"{sex.capitalize()} {STAB_LABEL[stab]}"
            switching_genes[label] = genes_by_sex_class[sex][stab]

    create_upset_plot(
        switching_genes,
        f"Gene overlaps: Male vs Female (Switching classes) - {version.upper()}",
        f"UpSet_MaleFemale_Switching_{version}",
        p2g_outdir,
    )

    # Plot 4: All sex x stability combinations
    print(f"\n  Plot 4: All Sex x Stability ({version})")
    all_genes = {
        SHORT_LABEL[(sex, stab)]: genes_by_sex_class[sex][stab]
        for sex in SEXES
        for stab in STABILITY_CLASSES
    }
    create_upset_plot(
        all_genes,
        f"Gene overlaps: All Sex x Stability - {version.upper()}",
        f"UpSet_All_SexStability_{version}",
        p2g_outdir,
        min_subset_size=DEFAULT_MIN_SUBSET_SIZE,
        max_intersections=DEFAULT_MAX_INTERSECTIONS,
    )

    # Summary table: male vs female overlap per stability class
    print(f"\n  Summary table ({version})")
    rows = []
    for stab in STABILITY_CLASSES:
        m = genes_by_sex_class["male"][stab]
        f = genes_by_sex_class["female"][stab]
        rows.append({
            "Stability": stab,
            "Male": len(m),
            "Female": len(f),
            "Overlap": len(m & f),
            "Male_only": len(m - f),
            "Female_only": len(f - m),
        })
    summary_df = pd.DataFrame(rows)
    summary_path = os.path.join(
        p2g_outdir, f"male_female_overlap_summary_{version}.tsv"
    )
    summary_df.to_csv(summary_path, sep="\t", index=False)

    print(f"  [OK] Saved: {os.path.basename(summary_path)}")
    print(summary_df.to_string(index=False))


# ==============================================================================
# COMPARISON: DEFAULT VS RELAXED
# ==============================================================================

print()
print("=" * 70)
print("COMPARISON: DEFAULT vs RELAXED")
print("=" * 70)

comparison_rows = []
for stab in STABILITY_CLASSES:
    for sex in SEXES:
        n_default = len(genes_by_version["default"][sex][stab])
        n_relaxed = len(genes_by_version["relaxed"][sex][stab])
        diff = n_relaxed - n_default
        pct_increase = 100 * diff / n_default if n_default > 0 else 0

        comparison_rows.append({
            "Sex": sex,
            "Stability": stab,
            "Default_0.45": n_default,
            "Relaxed_0.25": n_relaxed,
            "Difference": diff,
            "Pct_Increase": f"{pct_increase:.1f}%",
        })

comparison_df = pd.DataFrame(comparison_rows)
comparison_path = os.path.join(base_outdir, "comparison_default_vs_relaxed.tsv")
comparison_df.to_csv(comparison_path, sep="\t", index=False)

print(f"\n  [OK] Saved: {comparison_path}")
print(comparison_df.to_string(index=False))


# ==============================================================================
# SUMMARY
# ==============================================================================

print()
print("=" * 70)
print("  COMPLETE - ALL UPSET PLOTS GENERATED")
print("=" * 70)
print()
print("  Per cutoff:")
print("    UpSet_Male_Stability_[version].pdf/png")
print("    UpSet_Female_Stability_[version].pdf/png")
print("    UpSet_MaleFemale_Switching_[version].pdf/png")
print("    UpSet_All_SexStability_[version].pdf/png")
print("    male_female_overlap_summary_[version].tsv")
print()
print("  Comparison:")
print("    comparison_default_vs_relaxed.tsv")
print("=" * 70)
