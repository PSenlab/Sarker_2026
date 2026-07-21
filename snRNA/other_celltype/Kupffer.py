#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Compartment downstream analysis — liver aging multiome (RNA modality).

Runs, for a single annotated compartment subset:
  1. Per-subcluster differential expression (Wilcoxon) -> one CSV per subcluster
  2. Per-subcluster expression stats (pct expressed + mean expression) -> TSV
  3. ATAC label export (rna_barcode -> subcluster) for ArchR transfer -> CSV
  4. Composition vs. age/sex (% of all liver cells), counts + % + per-sex ANOVA
     -> multi-block .xlsx + CSVs, and a subcluster x age boxplot faceted by sex

Inputs
------
  --subset   annotated compartment .h5ad  (obs[group_col] = subcluster labels;
             obs_names = RNA barcodes e.g. 'AAACAGCCAAACCCTA-geriatric_5')
  --parent   full liver atlas .h5ad, used ONLY as the composition denominator
             (all liver cells per sample). Required for step 4.

The subset object is expected to carry:
  .X                 log-normalized expression (used for DE + expression stats)
  .obs[group_col]    subcluster annotation (categorical)
  .obs[sample_key]   per-sample id
  .obs['age']        age group  (young / mid_age / old / pre_geriatric / geriatric)
  .obs['sex']        'male' / 'female'

Subcluster labels are used VERBATIM in table/plot/label outputs; only the DE
filenames are sanitized (space/hyphen -> underscore) so they glob cleanly.

Usage
-----
    python compartment_downstream.py \
        --subset  Kupffer.h5ad \
        --parent  liver_atlas.h5ad \
        --group-col celltype_sub \
        --outdir  results/Kupffer \
        --prefix  Kupffer \
        --order "Kupffer" "Kupffer cycling" "LAM" "LSEC-like"

Omit --order to use the object's stored categorical order. Steps can be
toggled with --skip-de / --skip-expr / --skip-labels / --skip-composition.

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
# Constants (dataset conventions)
# --------------------------------------------------------------------------- #
AGE_GROUPS = ["young", "mid_age", "old", "pre_geriatric", "geriatric"]
SEX_ORDER = ["female", "male"]
AGE_DISP = {"young": "young", "mid_age": "mid-age", "old": "old",
            "pre_geriatric": "pre-geriatric", "geriatric": "geriatric"}
AGE_PALETTE = {"young": "#1ABC9C", "mid_age": "#F1C40F", "old": "#C39BD3",
               "pre_geriatric": "#2980B9", "geriatric": "#E84393"}
# alphabetical age order for the pooled counts table rows (matches paper layout)
AGE_ALPHA = ["geriatric", "mid_age", "old", "pre_geriatric", "young"]


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #
def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="Downstream DE + expression + label export + composition "
                    "for one annotated compartment subset.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--subset", required=True,
                   help="Annotated compartment .h5ad (subcluster labels in obs).")
    p.add_argument("--parent", default=None,
                   help="Full liver atlas .h5ad (composition denominator). "
                        "Required unless --skip-composition.")
    p.add_argument("--group-col", default="celltype_sub",
                   help="obs column holding subcluster labels.")
    p.add_argument("--sample-key", default="sample",
                   help="obs column holding per-sample id.")
    p.add_argument("--outdir", required=True, help="Output directory.")
    p.add_argument("--prefix", default=None,
                   help="Filename prefix for outputs (default: --group-col).")
    p.add_argument("--order", nargs="+", default=None,
                   help="Explicit subcluster order (verbatim labels). "
                        "Default: object's stored categorical order.")

    # step toggles
    p.add_argument("--skip-de", action="store_true", help="Skip Wilcoxon DE.")
    p.add_argument("--skip-expr", action="store_true", help="Skip expression stats.")
    p.add_argument("--skip-labels", action="store_true", help="Skip ATAC label export.")
    p.add_argument("--skip-composition", action="store_true",
                   help="Skip composition table + boxplot (needs --parent).")
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


def resolve_order(adata, group_col, order):
    """Return the subcluster order, validating against the object's labels."""
    col = adata.obs[group_col]
    stored = list(col.cat.categories) if hasattr(col, "cat") \
        else sorted(col.dropna().unique())
    if order is None:
        return stored
    missing = set(order) - set(stored)
    extra = set(stored) - set(order)
    if missing or extra:
        raise ValueError(
            f"--order mismatch. object has {sorted(stored)}, "
            f"--order has {sorted(order)} (missing={sorted(missing)}, "
            f"extra={sorted(extra)})"
        )
    return list(order)


def safe_name(label):
    """Filename-safe version of a subcluster label (label itself unchanged)."""
    return label.replace(" ", "_").replace("-", "_")


# --------------------------------------------------------------------------- #
# Steps
# --------------------------------------------------------------------------- #
def run_de(adata, group_col, order, outdir, prefix):
    """Per-subcluster Wilcoxon DE, one CSV per subcluster (Seurat-style cols)."""
    sc.tl.rank_genes_groups(adata, groupby=group_col, method="wilcoxon", pts=True)
    for grp in order:
        df = sc.get.rank_genes_groups_df(adata, group=grp, key="rank_genes_groups")
        df = df.rename(columns={"names": "gene", "logfoldchanges": "avg_log2FC",
                                "pvals_adj": "p_val_adj", "pvals": "p_val",
                                "scores": "score"}).set_index("gene")
        fn = os.path.join(outdir, f"{prefix}_DE_Cluster_{safe_name(grp)}.csv")
        df.to_csv(fn)
        print(f"  saved {os.path.basename(fn)}  ({df.shape[0]} genes)")


def run_expr_stats(adata, group_col, outdir, prefix):
    """Per-subcluster pct-expressed and mean-expression, long format TSV."""
    X = adata.X.toarray() if hasattr(adata.X, "toarray") else np.asarray(adata.X)
    expr = pd.DataFrame(X, columns=adata.var_names, index=adata.obs_names)
    grp = adata.obs[group_col]
    pct = expr.gt(0).groupby(grp, observed=True).mean().T * 100
    mean = expr.groupby(grp, observed=True).mean().T
    pct_long = (pct.reset_index()
                .melt(id_vars="index", var_name="celltype", value_name="pct_expressed")
                .rename(columns={"index": "gene"}))
    mean_long = (mean.reset_index()
                 .melt(id_vars="index", var_name="celltype", value_name="avg_expression")
                 .rename(columns={"index": "gene"}))
    stats = pct_long.merge(mean_long, on=["gene", "celltype"])
    fn = os.path.join(outdir, f"{prefix}_expr_stats.tsv")
    stats.to_csv(fn, sep="\t", index=False)
    print(f"  saved {os.path.basename(fn)}  {stats.shape}")


def run_label_export(adata, group_col, outdir, prefix):
    """Export rna_barcode -> subcluster label (verbatim) for ArchR transfer."""
    lab = pd.DataFrame({
        "rna_barcode": adata.obs_names,
        group_col: adata.obs[group_col].astype(str).values,
    })
    fn = os.path.join(outdir, f"{prefix}_sub_labels.csv")
    lab.to_csv(fn, index=False)
    print(f"  saved {os.path.basename(fn)}  ({lab.shape[0]} cells)")
    print(lab[group_col].value_counts().to_string())


def run_composition(sub, parent, group_col, sample_key, order, outdir, prefix):
    """% of ALL liver cells per sample; counts + % + per-sex ANOVA; boxplot."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from scipy.stats import f_oneway
    from statsmodels.stats.multitest import multipletests
    plt.rcParams["font.family"] = "Arial"

    # denominator: all liver cells per sample (parent)
    total = parent.obs.groupby(sample_key, observed=True).size().rename("total")

    # numerator: subcluster counts per sample, aligned to every parent sample
    counts = (sub.obs.groupby([sample_key, group_col], observed=True).size()
              .unstack(fill_value=0).reindex(columns=order, fill_value=0)
              .reindex(index=total.index, fill_value=0))

    meta = (parent.obs[[sample_key, "age", "sex"]].drop_duplicates()
            .set_index(sample_key).reindex(total.index))
    cwm = counts.join(meta)

    # per-sample % -> feeds ANOVA and boxplot
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

    # --- write xlsx (3 side-by-side blocks) + CSVs ---
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
    fig, axes = plt.subplots(1, len(SEX_ORDER), figsize=(4 * len(SEX_ORDER), 3.5))
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

    # sanity
    pooled = int(counts_tbl.loc["total cell count", "total cell count"])
    print(f"  parent cells: {parent.n_obs} | subset cells: {sub.n_obs} | "
          f"pooled denominator: {pooled}"
          + ("  OK" if pooled == parent.n_obs else "  <-- MISMATCH, check sample ids"))


# --------------------------------------------------------------------------- #
# Main
# --------------------------------------------------------------------------- #
def main(argv=None):
    args = parse_args(argv)
    prefix = args.prefix or args.group_col
    os.makedirs(args.outdir, exist_ok=True)
    log_versions()

    print(f"loading subset: {args.subset}")
    sub = sc.read_h5ad(args.subset)
    if args.group_col not in sub.obs:
        raise KeyError(f"{args.group_col!r} not in subset .obs "
                       f"(have: {list(sub.obs.columns)})")
    order = resolve_order(sub, args.group_col, args.order)
    print(f"subcluster order: {order}")

    if not args.skip_de:
        print("[1/4] differential expression")
        run_de(sub, args.group_col, order, args.outdir, prefix)
    if not args.skip_expr:
        print("[2/4] expression stats")
        run_expr_stats(sub, args.group_col, args.outdir, prefix)
    if not args.skip_labels:
        print("[3/4] ATAC label export")
        run_label_export(sub, args.group_col, args.outdir, prefix)
    if not args.skip_composition:
        if args.parent is None:
            raise ValueError("--parent is required for composition "
                             "(or pass --skip-composition).")
        print("[4/4] composition + boxplot")
        parent = sc.read_h5ad(args.parent)
        run_composition(sub, parent, args.group_col, args.sample_key,
                        order, args.outdir, prefix)

    print("done.")


if __name__ == "__main__":
    main()
