#!/usr/bin/env python3
# ==============================================================================
# Ascl1+ vs Ascl1- female hepatocytes, by lobular zone  -  4-STEP PIPELINE
# ==============================================================================
#
#   STEP 1  load h5ad -> subset to female hepatocytes        -> hep_female.h5ad
#   STEP 2  define lobular zone from subcluster              -> (in-memory obs)
#   STEP 3  Wilcoxon DE per zone (Ascl1_pos vs Ascl1_neg)    -> results/*.csv
#   STEP 4  volcano per zone + combined 1x3 row              -> figures/
#
#   Every step is CHECKPOINTED. Step 3 skips work whose output csv already
#   exists (set REUSE_EXISTING = False to force a recompute), and step 4
#   reads those csvs off disk. So you can rerun just the figures:
#
#       python 04_Ascl1_zonation.py --steps 4
#
#   Step 4 does not import scanpy (imported lazily inside steps 1 and 3), so
#   figure-only reruns work in a light environment.
#
#   SCOPE: the three lobular zones. The pooled all-female-hepatocyte
#   comparison is in 01_Ascl1_positive_negative.py.
# ==============================================================================

import argparse
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

warnings.filterwarnings("ignore")

mpl.rcParams["font.family"] = "Arial"
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42
# Editable SVG: keeps text as <text> elements instead of outlining every glyph.
mpl.rcParams["svg.fonttype"] = "none"


# ==============================================================================
# CONFIGURATION
# ==============================================================================

ADATA_PATH = "rna_wnn.h5ad"
HEP_CACHE  = Path("hep_female.h5ad")     # STEP 1 checkpoint

SEX_COL, CELLTYPE_COL, CT2_COL = "sex", "cell_type", "subcluster"
SEX_VALUE, CELLTYPE_VALUE = "female", "hepatocyte"
ZONE_COL = "Zonation"
GENE     = "Ascl1"

HEP_ZONATION_MAP = {
    "hep_01": "Periportal",  "hep_02": "Midlobular",  "hep_03": "Pericentral",
    "hep_04": "Pericentral", "hep_05": "Periportal",  "hep_06": "Midlobular",
    "hep_07": "Midlobular",
}
ZONE_ORDER = ["Periportal", "Midlobular", "Pericentral"]

MIN_CELLS_PER_GROUP = 20     # skip a zone if either Ascl1 group is smaller
PVAL_THRESHOLD      = 0.05
LOGFC_THRESHOLD     = 0      # log2 units; 0 = adjusted p alone defines significance

REUSE_EXISTING = True        # skip STEP 3 where the output csv already exists
FORMATS = ("pdf", "png", "svg")

FIG_DIR = Path("figures")
RES_DIR = Path("results")
FIG_DIR.mkdir(exist_ok=True)
RES_DIR.mkdir(exist_ok=True)

DE_STEM   = "wilcoxon_Ascl1_pos_vs_neg_female_hep_{zone}"

# --- volcano styling ----------------------------------------------------------
XLIM = (-5, 5)               # None to autoscale
RASTERIZE_POINTS = False     # True shrinks SVG/PDF; labels stay vector

COLORS = {"up": "lightpink", "down": "#008080", "ns": "lightgray"}
# lightpink is unreadable as bold text on white, so the count footer darkens it
TEXT_COLORS = {"up": "#c9647f", "down": "#008080"}

GENES_TO_LABEL = {
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
    "Acute Phase Proteins": {
        "genes": ["Crp", "Saa1", "Saa2", "Lcn2"],
        "color": "#A6D854",
    },
    "Others": {
        "genes": ["Lifr", "Esr1", "Angpt1", "Cd36", "Vwf",
                  "Fgl1", "Hif1a", "Hamp", "Hamp2"],
        "color": "#8B4513",
    },
}


# ==============================================================================
# HELPERS
# ==============================================================================

def banner(text):
    line = "=" * 70
    print(f"\n{line}\n{text}\n{line}")


def de_path(zone):
    return RES_DIR / f"{DE_STEM.format(zone=zone)}.csv"


def savefig(fig, directory, stem):
    for ext in FORMATS:
        fig.savefig(Path(directory) / f"{stem}.{ext}", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"    [OK] {directory}/{stem}.{{{','.join(FORMATS)}}}")


# ==============================================================================
# STEP 1 - LOAD + SUBSET TO FEMALE HEPATOCYTES
# ==============================================================================

def step1_load_subset():
    """Read the full object, keep female hepatocytes, attach Ascl1 expression."""
    banner("STEP 1: load h5ad + subset to female hepatocytes")
    import scanpy as sc

    if REUSE_EXISTING and HEP_CACHE.exists():
        print(f"  Reusing checkpoint {HEP_CACHE}")
        hep_female = sc.read_h5ad(HEP_CACHE)
        print(f"  [OK] {hep_female.n_obs:,} cells x {hep_female.n_vars:,} genes")
        return hep_female

    print(f"  Loading {ADATA_PATH}")
    adata = sc.read_h5ad(ADATA_PATH)
    print(f"  [OK] full object: {adata.n_obs:,} cells x {adata.n_vars:,} genes")

    for col in (CELLTYPE_COL, SEX_COL, CT2_COL):
        if col not in adata.obs.columns:
            raise ValueError(f"obs has no '{col}' column")
    if GENE not in adata.var_names:
        raise ValueError(f"{GENE} not found in var_names")

    hep_female = adata[
        (adata.obs[CELLTYPE_COL] == CELLTYPE_VALUE)
        & (adata.obs[SEX_COL] == SEX_VALUE)
    ].copy()
    del adata
    print(f"  Subset to {SEX_VALUE} {CELLTYPE_VALUE}: "
          f"{hep_female.n_obs:,} cells x {hep_female.n_vars:,} genes")

    # Ascl1 expression as a plain 1-D vector (handles sparse or dense X)
    expr = hep_female[:, GENE].X
    if hasattr(expr, "toarray"):
        expr = expr.toarray()
    hep_female.obs["Ascl1_expr"] = np.asarray(expr).ravel()
    hep_female.obs["Ascl1_positive"] = hep_female.obs["Ascl1_expr"] > 0

    n_pos = int(hep_female.obs["Ascl1_positive"].sum())
    pct = 100 * n_pos / hep_female.n_obs if hep_female.n_obs else 0
    print(f"  Ascl1+: {n_pos:,} / {hep_female.n_obs:,} ({pct:.2f}%)")

    hep_female.write(HEP_CACHE)
    print(f"  [OK] checkpoint written -> {HEP_CACHE}")
    return hep_female


# ==============================================================================
# STEP 2 - DEFINE LOBULAR ZONE
# ==============================================================================

def step2_define_zones(hep_female):
    """Map subcluster (hep_01..hep_07) to Periportal / Midlobular / Pericentral."""
    banner("STEP 2: define lobular zone from subcluster")

    ct2 = hep_female.obs[CT2_COL].astype(str)
    zones = ct2.map(HEP_ZONATION_MAP)

    unmapped = sorted(set(ct2[zones.isna()].unique()))
    if unmapped:
        print(f"  [WARN] subcluster labels with no zone mapping (dropped from all "
              f"downstream steps): {unmapped}")

    hep_female.obs[ZONE_COL] = pd.Categorical(
        zones, categories=ZONE_ORDER, ordered=True)

    counts = hep_female.obs[ZONE_COL].value_counts().reindex(ZONE_ORDER)
    print("  Cells per zone:")
    for z in ZONE_ORDER:
        n = int(counts.get(z, 0))
        n_pos = int(((hep_female.obs[ZONE_COL] == z)
                     & hep_female.obs["Ascl1_positive"]).sum())
        print(f"    {z:<12} {n:>8,} cells   Ascl1+ {n_pos:>6,}")
    return hep_female


# ==============================================================================
# STEP 3 - WILCOXON DE PER ZONE
# ==============================================================================

def step3_wilcoxon(hep_female):
    """Wilcoxon rank-sum Ascl1_pos vs Ascl1_neg inside each zone."""
    banner("STEP 3: Wilcoxon DE per zone")
    import scanpy as sc

    summary_rows = []
    for zone in ZONE_ORDER:
        print(f"\n  --- {zone} ---")
        out = de_path(zone)

        adata_zone = hep_female[hep_female.obs[ZONE_COL] == zone].copy()
        if adata_zone.n_obs == 0:
            print("    [SKIP] no cells in this zone")
            continue

        adata_zone.obs["Ascl1_status"] = pd.Categorical(
            np.where(adata_zone.obs["Ascl1_expr"] > 0, "Ascl1_pos", "Ascl1_neg"),
            categories=["Ascl1_pos", "Ascl1_neg"], ordered=True)
        n_pos = int((adata_zone.obs["Ascl1_status"] == "Ascl1_pos").sum())
        n_neg = int((adata_zone.obs["Ascl1_status"] == "Ascl1_neg").sum())
        print(f"    cells {adata_zone.n_obs:,} | Ascl1_pos {n_pos:,} | "
              f"Ascl1_neg {n_neg:,}")

        if n_pos < MIN_CELLS_PER_GROUP or n_neg < MIN_CELLS_PER_GROUP:
            print(f"    [SKIP] a group is below "
                  f"MIN_CELLS_PER_GROUP={MIN_CELLS_PER_GROUP}")
            summary_rows.append({"zone": zone, "n_cells": adata_zone.n_obs,
                                 "n_Ascl1_pos": n_pos, "n_Ascl1_neg": n_neg,
                                 "n_sig": 0, "n_up": 0, "n_down": 0,
                                 "status": "skipped_small_group"})
            continue

        if REUSE_EXISTING and out.exists():
            print(f"    Reusing {out}")
            de_res = pd.read_csv(out)
        else:
            sc.tl.rank_genes_groups(
                adata_zone, groupby="Ascl1_status",
                groups=["Ascl1_pos"], reference="Ascl1_neg",
                method="wilcoxon", use_raw=False)
            de_res = sc.get.rank_genes_groups_df(adata_zone, group="Ascl1_pos")
            de_res.to_csv(out, index=False)
            print(f"    [OK] {out}")

        sig = de_res[de_res["pvals_adj"] < PVAL_THRESHOLD]
        n_up   = int((sig["logfoldchanges"] >  LOGFC_THRESHOLD).sum())
        n_down = int((sig["logfoldchanges"] < -LOGFC_THRESHOLD).sum())
        print(f"    significant (adj p < {PVAL_THRESHOLD}): {len(sig):,} "
              f"(up {n_up:,}, down {n_down:,})")

        summary_rows.append({"zone": zone, "n_cells": adata_zone.n_obs,
                             "n_Ascl1_pos": n_pos, "n_Ascl1_neg": n_neg,
                             "n_sig": len(sig), "n_up": n_up, "n_down": n_down,
                             "status": "ok"})

    summary = pd.DataFrame(summary_rows)
    print("\n  Cross-zone summary:")
    print(summary.to_string(index=False))
    summary.to_csv(RES_DIR / "zone_summary.csv", index=False)
    print(f"\n  [OK] {RES_DIR}/zone_summary.csv")
    return summary


# ==============================================================================
# STEP 4 - VOLCANO
# ==============================================================================

def _load_de(zone):
    """Read a DE csv and add -log10(adj p) plus the up/down/ns category."""
    path = de_path(zone)
    if not path.exists():
        print(f"    [skip] missing {path}")
        return None

    df = pd.read_csv(path).dropna(subset=["logfoldchanges", "pvals_adj"]).copy()
    df["names"] = df["names"].astype(str)

    nonzero = df.loc[df["pvals_adj"] > 0, "pvals_adj"]
    floor = nonzero.min() if len(nonzero) else 1e-300
    df["neglog10_fdr"] = -np.log10(df["pvals_adj"].clip(lower=floor))

    df["cat"] = "ns"
    up   = (df["pvals_adj"] < PVAL_THRESHOLD) & (df["logfoldchanges"] >  LOGFC_THRESHOLD)
    down = (df["pvals_adj"] < PVAL_THRESHOLD) & (df["logfoldchanges"] < -LOGFC_THRESHOLD)
    df.loc[up, "cat"], df.loc[down, "cat"] = "up", "down"
    return df


def _draw_volcano(ax, df, title):
    for c in ("ns", "down", "up"):
        d = df[df["cat"] == c]
        ax.scatter(d["logfoldchanges"], d["neglog10_fdr"],
                   c=COLORS[c], edgecolor="none",
                   s=10 if c == "ns" else 16,
                   alpha=0.5 if c == "ns" else 0.9,
                   rasterized=RASTERIZE_POINTS,
                   label={"up": "Up in Ascl1+", "down": "Down in Ascl1+",
                          "ns": "Not significant"}[c])

    ax.axhline(-np.log10(PVAL_THRESHOLD), color="gray", ls="--", lw=0.8)
    if LOGFC_THRESHOLD > 0:
        ax.axvline( LOGFC_THRESHOLD, color="gray", ls="--", lw=0.8)
        ax.axvline(-LOGFC_THRESHOLD, color="gray", ls="--", lw=0.8)
    else:
        ax.axvline(0, color="gray", ls="--", lw=0.8)

    if XLIM is not None:
        ax.set_xlim(*XLIM)
    ax.set_xlabel("Log2 fold change (Ascl1+ vs Ascl1-)", fontsize=11, fontweight="bold")
    ax.set_ylabel("-log10(adjusted p-value)", fontsize=11, fontweight="bold")
    ax.set_title(title, fontsize=12, fontweight="bold")
    ax.spines[["top", "right"]].set_visible(False)

    lookup = df.set_index("names")
    positions, present_cats = [], []
    for cat, info in GENES_TO_LABEL.items():
        hits = [g for g in info["genes"] if g in lookup.index]
        if hits:
            present_cats.append(cat)
        for g in hits:
            r = lookup.loc[g]
            if isinstance(r, pd.DataFrame):
                r = r.iloc[0]
            positions.append({"gene": g, "x": float(r["logfoldchanges"]),
                              "y": float(r["neglog10_fdr"]), "color": info["color"]})
    return positions, present_cats


def _label_genes(ax, positions, fontsize=9):
    texts = [ax.text(p["x"], p["y"], p["gene"], fontsize=fontsize,
                     fontweight="bold", color=p["color"], ha="center",
                     va="center", zorder=10,
                     bbox=dict(boxstyle="round,pad=0.25", facecolor="white",
                               edgecolor="none", alpha=0.85))
             for p in positions]
    try:
        from adjustText import adjust_text
        adjust_text(texts, ax=ax,
                    arrowprops=dict(arrowstyle="-", color="gray", lw=0.5, alpha=0.6),
                    expand_points=(1.5, 1.5), expand_text=(1.2, 1.2),
                    force_points=(0.5, 0.5), force_text=(0.5, 0.5))
    except ImportError:
        print("    [info] adjustText not installed - labels drawn unadjusted")


def _annotate_counts(ax, df, fontsize=9):
    n_up = int((df["cat"] == "up").sum())
    n_dn = int((df["cat"] == "down").sum())
    box = dict(boxstyle="round,pad=0.3", facecolor="white",
               edgecolor="none", alpha=0.75)
    ax.text(0.02, 0.02, f"\u2193 {n_dn:,} down", transform=ax.transAxes,
            ha="left", va="bottom", fontsize=fontsize, fontweight="bold",
            color=TEXT_COLORS["down"], bbox=box, zorder=12)
    ax.text(0.98, 0.02, f"\u2191 {n_up:,} up", transform=ax.transAxes,
            ha="right", va="bottom", fontsize=fontsize, fontweight="bold",
            color=TEXT_COLORS["up"], bbox=box, zorder=12)


def step4_volcano():
    banner("STEP 4: volcano plots")
    loaded = []
    for zone in ZONE_ORDER:
        df = _load_de(zone)
        if df is not None:
            loaded.append((zone, df))
            print(f"    {zone}: {len(df):,} genes | "
                  f"up {(df['cat'] == 'up').sum():,} | "
                  f"down {(df['cat'] == 'down').sum():,}")
    if not loaded:
        print("    nothing to plot - run STEP 3 first")
        return

    for zone, df in loaded:
        fig, ax = plt.subplots(figsize=(7.5, 6))
        positions, cats = _draw_volcano(ax, df, f"Ascl1+ vs Ascl1- - {zone}")
        _label_genes(ax, positions)
        _annotate_counts(ax, df)
        sig_leg = ax.legend(frameon=True, loc="upper left", fontsize=8,
                            markerscale=1.8, title="Significance", title_fontsize=9)
        ax.legend(handles=[Patch(facecolor=GENES_TO_LABEL[c]["color"], label=c)
                           for c in cats],
                  frameon=True, framealpha=0.9, loc="upper left",
                  bbox_to_anchor=(1.02, 1.0), title="Functional categories",
                  fontsize=8, title_fontsize=9)
        ax.add_artist(sig_leg)
        fig.tight_layout()
        savefig(fig, FIG_DIR, f"volcano_Ascl1_{zone}")

    if len(loaded) > 1:
        n = len(loaded)
        fig, axes = plt.subplots(1, n, figsize=(5.2 * n, 5.4))
        axes = np.atleast_1d(axes)
        all_cats = []
        for ax, (zone, df) in zip(axes, loaded):
            positions, cats = _draw_volcano(ax, df, zone)
            _label_genes(ax, positions, fontsize=8)
            _annotate_counts(ax, df, fontsize=8)
            ax.legend(frameon=False, loc="upper left", fontsize=7, markerscale=1.5)
            all_cats += [c for c in cats if c not in all_cats]
        fig.legend(handles=[Patch(facecolor=GENES_TO_LABEL[c]["color"], label=c)
                            for c in all_cats],
                   loc="center left", bbox_to_anchor=(1.0, 0.5), frameon=True,
                   title="Functional categories", fontsize=9, title_fontsize=10)
        fig.suptitle("Ascl1+ vs Ascl1- female hepatocytes by lobular zone",
                     fontsize=14, fontweight="bold")
        fig.tight_layout(rect=[0, 0, 1, 0.95])
        savefig(fig, FIG_DIR, "volcano_Ascl1_COMBINED_zones")


# ==============================================================================
# MAIN
# ==============================================================================

def main():
    ap = argparse.ArgumentParser(
        description="Ascl1 zone pipeline. Steps are checkpointed; rerun a "
                    "subset with --steps (e.g. --steps 4 for figures only).")
    ap.add_argument("--steps", default="1,2,3,4",
                    help="comma-separated step numbers to run (default: all)")
    ap.add_argument("--force", action="store_true",
                    help="recompute even where output csv already exists")
    args = ap.parse_args()

    global REUSE_EXISTING
    if args.force:
        REUSE_EXISTING = False

    steps = {int(s) for s in args.steps.split(",") if s.strip()}
    print(f"Running steps: {sorted(steps)}   "
          f"(reuse existing outputs: {REUSE_EXISTING})")

    hep_female = None
    if 1 in steps:
        hep_female = step1_load_subset()
    if 2 in steps:
        if hep_female is None:
            hep_female = step1_load_subset()
        hep_female = step2_define_zones(hep_female)
    if 3 in steps:
        if hep_female is None or ZONE_COL not in hep_female.obs.columns:
            hep_female = step2_define_zones(step1_load_subset())
        step3_wilcoxon(hep_female)
    if 4 in steps:
        step4_volcano()

    banner("DONE")
    print(f"  Figures: {FIG_DIR}/    Tables: {RES_DIR}/")


if __name__ == "__main__":
    sys.exit(main())
