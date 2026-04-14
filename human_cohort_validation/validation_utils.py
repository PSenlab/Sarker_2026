#!/usr/bin/env python3
# ==============================================================================
# Shared Utilities for Human MASLD Cohort ASCL1 Validation
# ==============================================================================
#
# Description:
#   Common helper functions used by all human validation scripts in this
#   folder (Duke/TCGA-LIHC, German GSE33814, European GSE135251, Japanese
#   GSE167523, Japanese GSE193066, Statin GSE130991, Bariatric GSE83452).
#   Factoring these out keeps each per-cohort script focused on just the
#   dataset-specific preprocessing quirks.
#
# Provided helpers:
#   - banner()               : consistent section divider
#   - sig_label()            : p-value -> "*" / "**" / "***" / "n.s."
#   - parse_series_matrix()  : GEO series matrix -> (expr_df, meta_df)
#   - infer_sex_kmeans()     : KMeans sex inference from XIST + Y-gene proxy
#   - extract_gene()         : robust per-gene extraction with fallback
#   - mann_whitney_sex()     : male vs female Mann-Whitney U test
#   - barplot_male_vs_female(): publication bar + jittered points + bracket
#
# Reference:
#   Sarker et al. (2026) Cell Metabolism
#
# ==============================================================================

import gzip
import os
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy import stats
from sklearn.cluster import KMeans

# Global plot settings
mpl.rcParams["font.family"] = "Arial"
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42


# ==============================================================================
# SHARED CONFIGURATION
# ==============================================================================

SEX_COLORS = {
    "male":   "#3498db",
    "Male":   "#3498db",
    "MALE":   "#3498db",
    "M":      "#3498db",
    "female": "#e74c3c",
    "Female": "#e74c3c",
    "FEMALE": "#e74c3c",
    "F":      "#e74c3c",
}

JITTER_SEED = 42


# ==============================================================================
# HELPERS
# ==============================================================================

def banner(text):
    line = "=" * 70
    print()
    print(line)
    print(text)
    print(line)


def sig_label(p):
    """Return n.s. / * / ** / *** for a p-value."""
    if p < 0.001: return "***"
    if p < 0.01:  return "**"
    if p < 0.05:  return "*"
    return "n.s."


def parse_series_matrix(series_file, include_expr=True):
    """
    Parse a GEO series matrix (.txt.gz) into (expr_df, meta_df, platform_id).

    expr_df:  probe/gene x sample DataFrame (or None if include_expr=False)
    meta_df:  sample metadata (GSM, Title, and per-characteristic columns)
    platform_id: GPL accession string (e.g. "GPL6884")
    """
    metadata = {}
    char_rows = []
    expr_data = []
    reading_expr = False
    platform_id = None

    with gzip.open(series_file, "rt", errors="replace") as f:
        for line in f:
            if "!Series_platform_id" in line:
                platform_id = line.strip().split("\t")[1].strip('"')

            if line.startswith("!Sample_geo_accession"):
                metadata["GSM"] = [x.strip('"') for x in line.strip().split("\t")[1:]]
            elif line.startswith("!Sample_title"):
                metadata["Title"] = [x.strip('"') for x in line.strip().split("\t")[1:]]
            elif line.startswith("!Sample_characteristics_ch1"):
                vals = [x.strip('"') for x in line.strip().split("\t")[1:]]
                char_rows.append(vals)
                if ":" in vals[0]:
                    key = vals[0].split(":")[0].strip()
                    values = [
                        x.split(": ", 1)[1].strip('"') if ": " in x else ""
                        for x in vals
                    ]
                    if key not in metadata:
                        metadata[key] = values

            if line.startswith("!series_matrix_table_begin"):
                reading_expr = True
                continue
            if line.startswith("!series_matrix_table_end"):
                reading_expr = False
                continue
            if reading_expr and include_expr:
                expr_data.append(line.strip().split("\t"))

    meta = pd.DataFrame(metadata)

    expr_df = None
    if include_expr and expr_data:
        expr_df = pd.DataFrame(expr_data[1:], columns=expr_data[0])
        expr_df = expr_df.set_index(expr_df.columns[0])
        expr_df.index = [x.strip('"') for x in expr_df.index]
        expr_df.columns = [x.strip('"') for x in expr_df.columns]
        expr_df = expr_df.apply(pd.to_numeric, errors="coerce")

    return expr_df, meta, platform_id


def infer_sex_kmeans(log2_expr, sample_ids, xist_key, y_gene_keys,
                     female_label="Female", male_label="Male", random_state=42):
    """
    KMeans sex inference from XIST vs. Y-chromosome genes.

    log2_expr:    gene x sample DataFrame (log2-scale expression)
    sample_ids:   ordered sample IDs to classify
    xist_key:     index label for XIST (gene symbol or Ensembl ID)
    y_gene_keys:  list of possible Y-chromosome gene keys
    Returns:      list of sex labels aligned with sample_ids
    """
    if xist_key not in log2_expr.index:
        raise KeyError(f"XIST key '{xist_key}' not in expression index")

    y_available = [k for k in y_gene_keys if k in log2_expr.index]
    if not y_available:
        raise KeyError(f"No Y-chromosome genes in expression index "
                       f"(looked for {y_gene_keys})")

    xist_vals = log2_expr.loc[xist_key, sample_ids].values.astype(float)

    if len(y_available) == 1:
        y_proxy = log2_expr.loc[y_available[0], sample_ids].values.astype(float)
    else:
        y_proxy = log2_expr.loc[y_available, sample_ids].mean(axis=0).values.astype(float)

    features = np.column_stack([xist_vals, y_proxy])
    km = KMeans(n_clusters=2, random_state=random_state, n_init=10)
    labels = km.fit_predict(features)

    cluster_xist_mean = pd.Series(xist_vals).groupby(labels).mean()
    if cluster_xist_mean[0] > cluster_xist_mean[1]:
        sex_map = {0: female_label, 1: male_label}
    else:
        sex_map = {1: female_label, 0: male_label}

    print(f"  [OK] Sex inferred via KMeans (XIST + {len(y_available)} Y-gene(s))")
    print(f"    Cluster 0 XIST mean: {cluster_xist_mean[0]:.2f} -> {sex_map[0]}")
    print(f"    Cluster 1 XIST mean: {cluster_xist_mean[1]:.2f} -> {sex_map[1]}")

    return [sex_map[l] for l in labels]


def mann_whitney_sex(meta, ascl1_col, sex_col, male_label, female_label):
    """
    Run a two-sided Mann-Whitney U test of ASCL1 expression between sexes.
    Returns (mv_array, fv_array, u_stat, p_val, sig_str).
    """
    mv = meta.loc[meta[sex_col] == male_label, ascl1_col].dropna().values
    fv = meta.loc[meta[sex_col] == female_label, ascl1_col].dropna().values

    u_stat, p_val = stats.mannwhitneyu(mv, fv)
    sig = sig_label(p_val)

    print(f"  Male:   n={len(mv)}, mean={np.mean(mv):.4f} +/- {np.std(mv):.4f}")
    print(f"  Female: n={len(fv)}, mean={np.mean(fv):.4f} +/- {np.std(fv):.4f}")
    print(f"  U = {u_stat:.1f}, p = {p_val:.4e} {sig}")

    return mv, fv, u_stat, p_val, sig


def _add_sig_bracket(ax, x1, x2, y, p, h=0.02, lw=1.2, fontsize=11):
    """Draw a significance bracket with stars / n.s."""
    sig_text = sig_label(p)
    ax.plot([x1, x1, x2, x2], [y, y + h, y + h, y], lw=lw, color="black")
    ax.text(
        (x1 + x2) / 2, y + h, sig_text,
        ha="center", va="bottom", fontsize=fontsize, fontweight="bold",
    )


def barplot_male_vs_female(
    mv, fv, p_val,
    title, y_label, out_prefix,
    male_label="Male", female_label="Female",
    figsize=(5, 5.5),
):
    """
    Publication bar plot: mean +/- SEM with jittered points and a
    significance bracket. Writes both PNG and PDF.
    """
    fig, ax = plt.subplots(figsize=figsize)

    means = [np.mean(mv), np.mean(fv)]
    sems  = [stats.sem(mv), stats.sem(fv)]

    ax.bar(
        [0, 1], means, yerr=sems, width=0.55,
        color=[SEX_COLORS.get(male_label, "#3498db"),
               SEX_COLORS.get(female_label, "#e74c3c")],
        edgecolor="black", linewidth=1.2,
        capsize=5, error_kw={"lw": 1.5},
    )

    rng = np.random.default_rng(JITTER_SEED)
    for i, vals in enumerate([mv, fv]):
        jitter = rng.normal(0, 0.06, len(vals))
        ax.scatter(
            np.full(len(vals), i) + jitter, vals,
            color="black", alpha=0.35, s=15, zorder=3, linewidths=0,
        )

    y_max = max(np.max(mv), np.max(fv))
    bracket_y = y_max + 0.04 * abs(y_max)
    _add_sig_bracket(ax, 0, 1, bracket_y, p_val, h=0.015 * abs(y_max))

    ax.set_xticks([0, 1])
    ax.set_xticklabels(
        [f"{male_label}\n(n={len(mv)})", f"{female_label}\n(n={len(fv)})"],
        fontsize=11,
    )
    ax.set_ylabel(y_label, fontsize=12)
    ax.set_title(title, fontweight="bold", fontsize=13)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    plt.tight_layout()
    plt.savefig(f"{out_prefix}.png", dpi=300, bbox_inches="tight")
    plt.savefig(f"{out_prefix}.pdf", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  [OK] {out_prefix}.png / .pdf")
