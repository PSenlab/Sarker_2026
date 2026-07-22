#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Multi-compartment downstream analysis — liver aging multiome (RNA modality).

One script for every annotated compartment. For each, it runs:
  1. Marker dot plot (variance-scaled) supporting the annotation      -> PDF
  2. Per-subcluster differential expression (Wilcoxon)                -> CSV per subcluster
  3. Per-subcluster expression stats (pct expressed + mean expression)-> TSV
  4. ATAC label export (rna_barcode -> subcluster) for ArchR transfer -> CSV
  5. Composition vs. age/sex (% of ALL liver cells): counts + % +
     per-sex ANOVA (BH across subclusters)          -> multi-block .xlsx + CSVs
     plus a subcluster x age boxplot faceted by sex -> PDF

The full liver atlas is loaded ONCE and reused as the composition denominator
for every compartment, so all panels sit on a common baseline.

Pairs with the R scripts (compartment_marker_bubbles.R, browser tracks, module
scores): `prefix` here determines the filenames they read, so the `prefix` and
`labels` in COMPARTMENTS below must match their configs exactly.

Label vocabulary
----------------
`labels` are the VERBATIM subcluster names as stored in obs (e.g. "Kupffer
cycling", "LSEC-like", "MV portal"). They are used verbatim in every table,
plot axis, and the ATAC label CSV. ONLY the per-subcluster DE filenames are
sanitized (space/hyphen -> underscore) via safe_name(), which the R side
reproduces with its sanitize().

Each subset object is expected to carry:
  .X                 log-normalized expression (DE + expression stats)
  .obs[group_col]    subcluster annotation (categorical)
  .obs[sample_key]   per-sample id
  .obs['age']        young / mid_age / old / pre_geriatric / geriatric
  .obs['sex']        'male' / 'female'

Usage
-----
    python compartment_downstream.py --parent liver_atlas.h5ad
    python compartment_downstream.py --parent liver_atlas.h5ad --run Kupffer T_ILC
    python compartment_downstream.py --run myeloid --skip-composition

Edit COMPARTMENTS below to point at your .h5ad paths / output dirs.

Citation
--------
If you use this code, please cite:
    Sarker N. et al. "<paper title>." <journal>, <year>. DOI: <doi>

License: MIT (see LICENSE).
"""

import argparse
import os
import sys

import numpy as np
import pandas as pd
import scanpy as sc


# --------------------------------------------------------------------------- #
# Dataset conventions
# --------------------------------------------------------------------------- #
AGE_GROUPS = ["young", "mid_age", "old", "pre_geriatric", "geriatric"]
SEX_ORDER = ["female", "male"]
AGE_DISP = {"young": "young", "mid_age": "mid-age", "old": "old",
            "pre_geriatric": "pre-geriatric", "geriatric": "geriatric"}
AGE_PALETTE = {"young": "#1ABC9C", "mid_age": "#F1C40F", "old": "#C39BD3",
               "pre_geriatric": "#2980B9", "geriatric": "#E84393"}
# alphabetical age order for the pooled counts table rows (matches paper layout)
AGE_ALPHA = ["geriatric", "mid_age", "old", "pre_geriatric", "young"]
SAMPLE_KEY = "sample"


# --------------------------------------------------------------------------- #
# Per-compartment registry
#   h5ad      : annotated compartment subset
#   group_col : obs column holding the subcluster labels
#   labels    : VERBATIM labels, in publication order
#   markers   : dot plot panel, one entry per population
#   outdir    : where all outputs for this compartment are written
#   prefix    : filename prefix (must match the R scripts' `prefix`)
# --------------------------------------------------------------------------- #
COMPARTMENTS = {

    "myeloid": dict(
        h5ad="myeloid.h5ad",
        group_col="celltype_myeloid",
        outdir="/data/sarkern2/multiome_liver/myeloid_DE_csvs",
        prefix="myeloid",
        labels=["Kupffer", "Kupffer cycling", "LAM", "MoMF",
                "cDC1", "pDC", "Neutrophil"],
        markers={
            "Kupffer":         ["Clec4f", "Vsig4", "Timd4", "Cd5l", "Marco", "C1qa"],
            "Kupffer cycling": ["Mki67", "Top2a", "Cenpf"],
            "LAM":             ["Gpnmb", "Trem2", "Cd9", "Arhgap22", "Sirpb1a",
                                "Cd74", "Pparg"],
            "MoMF":            ["Ccr2", "Cx3cr1", "Ly6c2", "Chil3", "Plac8", "Fn1"],
            "cDC1":            ["Xcr1", "Clec9a", "Batf3", "Cadm1"],
            "pDC":             ["Siglech", "Bst2", "Il3ra", "Tcf4"],
            "Neutrophil":      ["S100a8", "S100a9", "Retnlg", "Csf3r"],
        },
    ),

    "Kupffer": dict(
        h5ad="Kupffer.h5ad",
        group_col="celltype_sub",
        outdir="/data/sarkern2/multiome_liver/Kupffer_DE_csvs",
        prefix="Kupffer",
        labels=["Kupffer", "Kupffer cycling", "LAM", "LSEC-like"],
        markers={
            "Kupffer":         ["Clec4f", "Vsig4", "Timd4", "Cd5l", "Marco", "C1qa"],
            "Kupffer cycling": ["Mki67", "Top2a", "Cenpf"],
            "LAM":             ["Gpnmb", "Trem2", "Cd9", "Arhgap22", "Sirpb1a",
                                "Cd74", "Pparg"],
            "LSEC-like":       ["Stab2", "Clec4g", "Dnase1l3", "Oit3", "Kdr",
                                "Gata4", "Ptprb", "Pecam1"],
        },
    ),

    "endothelial_Kupffer02": dict(
        h5ad="endothelial_Kupffer02.h5ad",
        group_col="endo_subcluster",
        outdir="/data/sarkern2/multiome_liver/endo_DE_csvs",
        prefix="endothelial_Kupffer02",
        labels=["LSEC", "LSEC cycling", "MV portal", "MV central", "Kupffer-like"],
        markers={
            "LSEC":         ["Stab2", "Clec4g", "Gata4", "Oit3", "Dnase1l3", "Mrc1"],
            "LSEC cycling": ["Mki67", "Top2a", "Cenpf"],
            "MV portal":    ["Vwf", "Sox17", "Efnb2"],
            "MV central":   ["Rspo3", "Wnt9b", "Wnt2", "Thbd"],
            "Kupffer-like": ["Clec4f", "Vsig4", "Timd4", "Cd5l", "C1qa",
                             "Marco", "Csf1r"],
        },
    ),

    "T_ILC": dict(
        h5ad="T_ILC.h5ad",
        group_col="celltype_T",
        outdir="/data/sarkern2/multiome_liver/T_DE_csvs",
        prefix="T_ILC",
        labels=["CD4T", "Treg", "CD8T", "gdT", "iNKT", "NK", "ILC1", "neutrophil"],
        markers={
            "CD4T":       ["Cd3e", "Cd4"],
            "Treg":       ["Foxp3", "Ikzf2"],
            "CD8T":       ["Cd8a", "Cd8b1", "Gzmk"],
            "gdT":        ["Trdc"],
            "iNKT":       ["Zbtb16"],
            "NK":         ["Ncr1", "Klrb1c", "Gzmb"],
            "ILC1":       ["Tbx21", "Eomes"],
            "neutrophil": ["S100a8", "S100a9", "Retnlg"],
        },
    ),
}


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #
def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="Multi-compartment downstream: dot plot + DE + expression "
                    "stats + ATAC label export + composition.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--parent", default=None,
                   help="Full liver atlas .h5ad (composition denominator). "
                        "Required unless --skip-composition.")
    p.add_argument("--run", nargs="+", default=None, choices=sorted(COMPARTMENTS),
                   help="Which compartments to run (default: all).")
    p.add_argument("--sample-key", default=SAMPLE_KEY,
                   help="obs column holding per-sample id.")
    p.add_argument("--skip-dotplot", action="store_true")
    p.add_argument("--skip-de", action="store_true")
    p.add_argument("--skip-expr", action="store_true")
    p.add_argument("--skip-labels", action="store_true")
    p.add_argument("--skip-composition", action="store_true")
    return p.parse_args(argv)


# --------------------------------------------------------------------------- #
# Helpers
# --------------------------------------------------------------------------- #
def log_versions():
    import anndata
    print("=== environment ===")
    print(f"python   {sys.version.split()[0]}")
    print(f"scanpy   {sc.__version__}")
    print(f"anndata  {anndata.__version__}")
    print(f"numpy    {np.__version__}")
    print(f"pandas   {pd.__version__}")
    print("===================")


def safe_name(label):
    """Filename-safe label (the label itself is never changed)."""
    return label.replace(" ", "_").replace("-", "_")


def prepare_subset(sub, group_col, labels):
    """Validate labels against the object and fix the categorical order."""
    col = sub.obs[group_col]
    present = set(col.dropna().unique())
    if present != set(labels):
        raise ValueError(
            f"label mismatch in {group_col!r}: object has {sorted(present)}, "
            f"config has {sorted(labels)}"
        )
    sub.obs[group_col] = pd.Categorical(col, categories=labels, ordered=True)
    sub.obs["age"] = pd.Categorical(sub.obs["age"], AGE_GROUPS, ordered=True)
    return sub


# --------------------------------------------------------------------------- #
# Steps
# --------------------------------------------------------------------------- #
def run_dotplot(sub, cfg):
    """Variance-scaled marker dot plot supporting the annotation."""
    import matplotlib.pyplot as plt
    present = {k: [g for g in v if g in sub.var_names] for k, v in cfg["markers"].items()}
    dropped = {g for v in cfg["markers"].values() for g in v} - set(sub.var_names)
    if dropped:
        print("  dropped (not in var_names):", ", ".join(sorted(dropped)))
    present = {k: v for k, v in present.items() if v}
    sc.pl.dotplot(sub, present, groupby=cfg["group_col"],
                  standard_scale="var", dendrogram=False, show=False)
    fn = os.path.join(cfg["outdir"], f"{cfg['prefix']}_RNA_dotplot_stdscale.pdf")
    plt.savefig(fn, bbox_inches="tight")
    plt.close("all")
    print(f"  saved {os.path.basename(fn)}")


def run_de(sub, cfg):
    """Per-subcluster Wilcoxon DE, one CSV per subcluster (Seurat-style cols)."""
    sc.tl.rank_genes_groups(sub, groupby=cfg["group_col"], method="wilcoxon", pts=True)
    for grp in cfg["labels"]:
        df = sc.get.rank_genes_groups_df(sub, group=grp, key="rank_genes_groups")
        df = df.rename(columns={"names": "gene", "logfoldchanges": "avg_log2FC",
                                "pvals_adj": "p_val_adj", "pvals": "p_val",
                                "scores": "score"}).set_index("gene")
        fn = os.path.join(cfg["outdir"],
                          f"{cfg['prefix']}_DE_Cluster_{safe_name(grp)}.csv")
        df.to_csv(fn)
        print(f"  saved {os.path.basename(fn)}  ({df.shape[0]} genes)")


def run_expr_stats(sub, cfg):
    """Per-subcluster pct-expressed and mean-expression, long format TSV."""
    X = sub.X.toarray() if hasattr(sub.X, "toarray") else np.asarray(sub.X)
    expr = pd.DataFrame(X, columns=sub.var_names, index=sub.obs_names)
    grp = sub.obs[cfg["group_col"]]
    pct = expr.gt(0).groupby(grp, observed=True).mean().T * 100
    mean = expr.groupby(grp, observed=True).mean().T
    pct_long = (pct.reset_index()
                .melt(id_vars="index", var_name="celltype", value_name="pct_expressed")
                .rename(columns={"index": "gene"}))
    mean_long = (mean.reset_index()
                 .melt(id_vars="index", var_name="celltype", value_name="avg_expression")
                 .rename(columns={"index": "gene"}))
    stats = pct_long.merge(mean_long, on=["gene", "celltype"])
    fn = os.path.join(cfg["outdir"], f"{cfg['prefix']}_expr_stats.tsv")
    stats.to_csv(fn, sep="\t", index=False)
    print(f"  saved {os.path.basename(fn)}  {stats.shape}")


def run_label_export(sub, cfg):
    """Export rna_barcode -> subcluster label (verbatim) for ArchR transfer."""
    gc = cfg["group_col"]
    lab = pd.DataFrame({"rna_barcode": sub.obs_names,
                        gc: sub.obs[gc].astype(str).values})
    fn = os.path.join(cfg["outdir"], f"{cfg['prefix']}_sub_labels.csv")
    lab.to_csv(fn, index=False)
    print(f"  saved {os.path.basename(fn)}  ({lab.shape[0]} cells)")
    print(lab[gc].value_counts().to_string())


def run_composition(sub, parent, cfg, sample_key):
    """% of ALL liver cells per sample; counts + % + per-sex ANOVA; boxplot."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from scipy.stats import f_oneway
    from statsmodels.stats.multitest import multipletests
    plt.rcParams["font.family"] = "Arial"

    group_col, order = cfg["group_col"], cfg["labels"]
    outdir, prefix = cfg["outdir"], cfg["prefix"]

    # denominator: all liver cells per sample (parent)
    total = parent.obs.groupby(sample_key, observed=True).size().rename("total")

    # numerator: subcluster counts per sample, aligned to every parent sample
    counts = (sub.obs.groupby([sample_key, group_col], observed=True).size()
              .unstack(fill_value=0).reindex(columns=order, fill_value=0)
              .reindex(index=total.index, fill_value=0))

    meta = (parent.obs[[sample_key, "age", "sex"]].drop_duplicates()
            .set_index(sample_key).reindex(total.index))
    cwm = counts.join(meta)

    pct = (100 * counts.div(total, axis=0)).join(meta)
    pct["age"] = pd.Categorical(pct["age"], AGE_GROUPS, ordered=True)

    # --- pooled counts table (sample groups as rows) ---
    rows = []
    for ag in AGE_ALPHA:
        for sx in ["female", "male"]:
            members = sorted(cwm.index[(cwm["sex"] == sx) & (cwm["age"] == ag)].tolist())
            if not members:
                continue
            n = counts.loc[members, order].sum(0)
            tot = int(total.loc[members].sum())
            rows.append(pd.Series({**n.to_dict(), "total cell count": tot},
                                  name=f"{AGE_DISP[ag]}_{sx} ({' + '.join(members)})"))
    counts_tbl = pd.DataFrame(rows).astype(int)
    counts_tbl.loc["total cell count"] = counts_tbl.sum(0)
    pct_tbl = (100 * counts_tbl[order].div(counts_tbl["total cell count"], axis=0)).round(2)

    # --- per-sex ANOVA across age, BH across subclusters ---
    stat_rows = []
    for sx in SEX_ORDER:
        ssub = pct[pct["sex"] == sx]
        F, P = {}, {}
        for ct in order:
            grp = [ssub.loc[ssub["age"] == a, ct].values for a in AGE_GROUPS]
            grp = [g for g in grp if len(g) >= 2]
            F[ct], P[ct] = f_oneway(*grp) if len(grp) >= 2 else (np.nan, np.nan)
        valid = [c for c in order if not np.isnan(P[c])]
        padj = (dict(zip(valid, multipletests([P[c] for c in valid],
                                              method="fdr_bh")[1])) if valid else {})
        for ct in order:
            stat_rows.append({"sex": sx, "subcluster": ct, "F_stat": F[ct],
                              "pval": P[ct], "padj": padj.get(ct, np.nan)})
    stats = pd.DataFrame(stat_rows)

    # --- xlsx (3 side-by-side blocks) + CSVs ---
    counts_out = counts_tbl.reset_index().rename(columns={"index": "sample group"})
    xlsx = os.path.join(outdir, f"{prefix}_subcluster_composition.xlsx")
    with pd.ExcelWriter(xlsx, engine="openpyxl") as xl:
        counts_out.to_excel(xl, "composition", startrow=1, startcol=0, index=False)
        c1 = counts_out.shape[1] + 1
        pct_tbl.reset_index(drop=True).to_excel(xl, "composition", startrow=1,
                                                startcol=c1, index=False)
        c2 = c1 + pct_tbl.shape[1] + 1
        stats.to_excel(xl, "composition", startrow=1, startcol=c2, index=False)
        ws = xl.sheets["composition"]
        ws.cell(row=1, column=1, value="counts")
        ws.cell(row=1, column=c1 + 1, value="percentage (% of all liver cells)")
        ws.cell(row=1, column=c2 + 1, value="statistics")

    counts_tbl.to_csv(os.path.join(outdir, f"{prefix}_counts_per_group.csv"))
    pct_tbl.to_csv(os.path.join(outdir, f"{prefix}_percentage_of_all_cells_per_group.csv"))
    stats.to_csv(os.path.join(outdir, f"{prefix}_stats_anova_age.csv"), index=False)
    print(f"  saved {os.path.basename(xlsx)} + 3 CSVs")

    # --- boxplot: subcluster x age, faceted by sex, FDR stars ---
    comp = (counts.reset_index()
            .melt(id_vars=sample_key, var_name=group_col, value_name="n")
            .merge(total, on=sample_key))
    comp["pct"] = 100 * comp["n"] / comp["total"]
    comp = comp.merge(meta.reset_index(), on=sample_key)
    comp["age"] = pd.Categorical(comp["age"], AGE_GROUPS, ordered=True)
    comp.to_csv(os.path.join(outdir, f"{prefix}_composition_per_sample.csv"), index=False)

    stars = lambda q: ("***" if q < 1e-3 else "**" if q < 1e-2
                       else "*" if q < 5e-2 else "n.s.")
    n_age = len(AGE_GROUPS)
    width = max(4, 0.55 * len(order) + 2)
    fig, axes = plt.subplots(1, len(SEX_ORDER), figsize=(width * len(SEX_ORDER), 3.5))
    for ax, sx in zip(np.atleast_1d(axes), SEX_ORDER):
        s = comp[comp.sex == sx]
        pv = {}
        for ct in order:
            g = [s[(s[group_col] == ct) & (s.age == a)]["pct"].values for a in AGE_GROUPS]
            g = [x for x in g if len(x) >= 2]
            pv[ct] = f_oneway(*g)[1] if len(g) >= 2 else np.nan
        cts = [c for c in order if not np.isnan(pv[c])]
        qmap = (dict(zip(cts, multipletests([pv[c] for c in cts],
                                            method="fdr_bh")[1])) if cts else {})
        for ci, ct in enumerate(order):
            for ai, ag in enumerate(AGE_GROUPS):
                vals = s[(s[group_col] == ct) & (s.age == ag)]["pct"].values
                pos = ci * (n_age + 1) + ai
                bp = ax.boxplot(vals, positions=[pos], widths=0.7, patch_artist=True,
                                flierprops=dict(marker="o", mfc="red", ms=3, mec="none"))
                bp["boxes"][0].set_facecolor(AGE_PALETTE[ag])
                bp["boxes"][0].set_alpha(0.85)
                bp["medians"][0].set_color("black")
            q = qmap.get(ct, np.nan)
            lab = stars(q) if not np.isnan(q) else "n.s."
            center = ci * (n_age + 1) + (n_age - 1) / 2
            ymax = s[s[group_col] == ct]["pct"].max()
            ax.text(center, (ymax * 1.05 + 1) if np.isfinite(ymax) else 1, lab,
                    ha="center", color="red" if lab != "n.s." else "darkblue",
                    fontsize=11, fontweight="bold")
        ax.set_xticks([ci * (n_age + 1) + (n_age - 1) / 2 for ci in range(len(order))])
        ax.set_xticklabels(order, rotation=45, ha="right")
        ax.set_title(sx)
        ax.set_ylabel("% of all cells")
        ax.set_facecolor("#f7f7f7")
    fig.legend([plt.Rectangle((0, 0), 1, 1, fc=AGE_PALETTE[a], alpha=0.85)
                for a in AGE_GROUPS], AGE_GROUPS, loc="center left",
               bbox_to_anchor=(1.0, 0.5), frameon=False)
    plt.tight_layout()
    pdf = os.path.join(outdir, f"{prefix}_composition_box.pdf")
    plt.savefig(pdf, bbox_inches="tight")
    plt.close(fig)
    print(f"  saved {os.path.basename(pdf)}")

    pooled = int(counts_tbl.loc["total cell count", "total cell count"])
    print(f"  parent: {parent.n_obs} | subset: {sub.n_obs} | pooled denominator: {pooled}"
          + ("  OK" if pooled == parent.n_obs else "  <-- MISMATCH, check sample ids"))


# --------------------------------------------------------------------------- #
# Main
# --------------------------------------------------------------------------- #
def main(argv=None):
    args = parse_args(argv)
    to_run = args.run or list(COMPARTMENTS)
    log_versions()

    # parent atlas loaded ONCE, shared across every compartment
    parent = None
    if not args.skip_composition:
        if args.parent is None:
            raise ValueError("--parent is required (or pass --skip-composition).")
        print(f"loading parent atlas: {args.parent}")
        parent = sc.read_h5ad(args.parent)
        print(f"  {parent.n_obs} cells, "
              f"{parent.obs[args.sample_key].nunique()} samples")

    for name in to_run:
        cfg = COMPARTMENTS[name]
        print(f"\n===== {name} =====")
        os.makedirs(cfg["outdir"], exist_ok=True)

        sub = sc.read_h5ad(cfg["h5ad"])
        if cfg["group_col"] not in sub.obs:
            raise KeyError(f"{cfg['group_col']!r} not in {cfg['h5ad']} .obs "
                           f"(have: {list(sub.obs.columns)})")
        sub = prepare_subset(sub, cfg["group_col"], cfg["labels"])
        print(f"loaded {cfg['h5ad']}: {sub.n_obs} cells, "
              f"{len(cfg['labels'])} subclusters")

        if not args.skip_dotplot:
            print("[1/5] marker dot plot")
            run_dotplot(sub, cfg)
        if not args.skip_de:
            print("[2/5] differential expression")
            run_de(sub, cfg)
        if not args.skip_expr:
            print("[3/5] expression stats")
            run_expr_stats(sub, cfg)
        if not args.skip_labels:
            print("[4/5] ATAC label export")
            run_label_export(sub, cfg)
        if not args.skip_composition:
            print("[5/5] composition + boxplot")
            run_composition(sub, parent, cfg, args.sample_key)

        del sub   # free memory before the next compartment

    print("\ndone.")


if __name__ == "__main__":
    main()
