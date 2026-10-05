#!/usr/bin/env python
# coding: utf-8
# ============================================================================
# Overlap of Ascl1+ vs Ascl1- DE gene lists with ASCL1 peak-annotated genes
# Figure 7, panel e (proportional Venn, all peaks)
#
#   query        = gene list read from Excel/CSV (no expression filtering)
#   peak filters = allPeaks, promoterOnly (annotation startswith "Promoter")
#
#   Table (per filter x sex): overlap counts and gene lists
#   Venn  (per filter): query  x  male peaks  x  female peaks
#   UpSet (all)       : query + {allPeaks,promoter} x {male,female}
#
# UpSet uses the upsetplot package, with a self-contained matplotlib fallback
# if the import/plot fails (upsetplot breaks on some matplotlib/pandas combos).
# All mouse-symbol space, no ortholog map.
#
# Input : all_hep_up.csv, all_hep_down.csv    genes up / down in Ascl1+ vs
#         Ascl1- female hepatocytes (results/ of 01_Ascl1_positive_negative.py)
#         chipseeker_out/all_peaks.csv        (07_chipseeker_annotation.R)
# Output (per gene list <tag>):
#         genelist_peak_overlap_<tag>.xlsx    summary + overlap gene lists
#         venn_proportional_<tag>_<filter>.pdf
#         upset_genelist_<tag>.pdf
# ============================================================================
import os
import numpy as np
import pandas as pd
from collections import Counter
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib_venn import venn3

# ---- style: Arial, editable text in PDF ----
matplotlib.rcParams["pdf.fonttype"] = 42
_avail = {f.name for f in font_manager.fontManager.ttflist}
for _f in ["Arial", "Liberation Sans", "Nimbus Sans", "DejaVu Sans"]:
    if _f in _avail:
        plt.rcParams["font.family"] = _f
        break

# ================================================================ config =====
GENE_LISTS = {                     # gene list file -> circle / row label
    "all_hep_up.csv":   "Ascl1+ up",
    "all_hep_down.csv": "Ascl1+ down",
}
GENES_SHEET = 0                    # sheet name or index
GENES_COL   = None                 # None = auto-detect, or e.g. "gene"
PEAKS_CSV   = "chipseeker_out/all_peaks.csv"
OUTDIR      = "."

PROPORTIONAL = True                # False -> equal-size circles (unweighted)
COLORS       = ["#4E79A7", "#59A14F", "#E15759"]   # query, male, female

# ================================================= load input gene list ======
def load_gene_list(GENES_XLSX):
    if str(GENES_XLSX).lower().endswith((".csv", ".tsv", ".txt")):
        sep = "\t" if str(GENES_XLSX).lower().endswith((".tsv", ".txt")) else ","
        gdf = pd.read_csv(GENES_XLSX, sep=sep)
    else:
        gdf = pd.read_excel(GENES_XLSX, sheet_name=GENES_SHEET)

    col = GENES_COL
    if col is None:
        # auto-detect: first column whose name looks like a gene column, else col 0
        cand = [c for c in gdf.columns
                if str(c).strip().lower() in
                {"gene", "genes", "symbol", "gene_symbol", "gene_name", "geneid"}]
        col = cand[0] if cand else gdf.columns[0]
        print(f"gene column: '{col}'  (auto-detected from {list(gdf.columns)})")

    input_genes = (gdf[col].dropna().astype(str).str.strip())
    input_genes = input_genes[input_genes != ""]
    query = set(input_genes)
    n_dup = len(input_genes) - len(query)
    print(f"input genes: {len(query)} unique"
          + (f"  ({n_dup} duplicate rows dropped)" if n_dup else ""))
    if not query:
        raise RuntimeError(f"no genes read from {GENES_XLSX} / column '{col}'")
    return query

# ================================================= peak gene sets ============
peaks = pd.read_csv(PEAKS_CSV)

def syms(df):
    return set(df["SYMBOL"].dropna().astype(str).str.strip())

def peak_gene_sets(peak_filter):
    d = peaks
    if peak_filter == "promoterOnly":
        d = d[d["annotation"].astype(str).str.startswith("Promoter")]
    return {"all":    syms(d),
            "male":   syms(d[d["sex"] == "male"]),
            "female": syms(d[d["sex"] == "female"])}

def peak_sets(pf):
    d = peaks[peaks["annotation"].astype(str).str.startswith("Promoter")] \
        if pf == "promoterOnly" else peaks
    return syms(d[d["sex"] == "male"]), syms(d[d["sex"] == "female"])

if not PROPORTIONAL:
    from matplotlib_venn.layout.venn3 import DefaultLayoutAlgorithm
    _LAYOUT = {"layout_algorithm":
               DefaultLayoutAlgorithm(fixed_subset_sizes=(1, 1, 1, 1, 1, 1, 1))}
else:
    _LAYOUT = {}

def upset_fallback(sets, outfile, min_size=1, title=""):
    """Self-contained UpSet, used if the upsetplot package misbehaves."""
    names = list(sets)
    universe = set().union(*sets.values())
    memb = {g: tuple(n for n in names if g in sets[n]) for g in universe}
    patt = Counter(c for c in memb.values() if c)
    rows = sorted([(c, s) for c, s in patt.items() if s >= min_size],
                  key=lambda x: -x[1])
    combos = [c for c, _ in rows]
    sizes = [s for _, s in rows]
    set_tot = {n: len(sets[n]) for n in names}
    n_c, n_s = len(combos), len(names)

    fig = plt.figure(figsize=(max(6, 0.7 * n_c + 3), 0.5 * n_s + 4))
    gs = fig.add_gridspec(2, 2, width_ratios=[1, 3.2], height_ratios=[3, 1.1],
                          wspace=0.05, hspace=0.05)
    ax_bar = fig.add_subplot(gs[0, 1])
    ax_mat = fig.add_subplot(gs[1, 1], sharex=ax_bar)
    ax_tot = fig.add_subplot(gs[1, 0])

    x = np.arange(n_c)
    ax_bar.bar(x, sizes, color="#4C4C4C", width=0.6)
    for xi, s in zip(x, sizes):
        ax_bar.text(xi, s, str(s), ha="center", va="bottom", fontsize=8)
    ax_bar.set_ylabel("Intersection size", fontsize=11)
    ax_bar.spines[["top", "right"]].set_visible(False)
    ax_bar.set_xticks([])
    ax_bar.margins(y=0.15)
    if title:
        ax_bar.set_title(title, fontsize=12, weight="bold")

    y = np.arange(n_s)[::-1]
    for xi, combo in enumerate(combos):
        ax_mat.scatter([xi] * n_s, y, color="#D9D9D9", s=90, zorder=1)
        yy = [y[names.index(nm)] for nm in combo]
        ax_mat.scatter([xi] * len(yy), yy, color="#4C4C4C", s=90, zorder=2)
        if len(yy) > 1:
            ax_mat.plot([xi, xi], [min(yy), max(yy)], color="#4C4C4C", lw=2, zorder=2)
    ax_mat.set_xlim(-0.5, n_c - 0.5)
    ax_mat.set_ylim(-0.6, n_s - 0.4)
    ax_mat.set_yticks(y)
    ax_mat.set_yticklabels(names, fontsize=9)
    ax_mat.set_xticks([])
    ax_mat.spines[:].set_visible(False)

    ax_tot.barh(y, [set_tot[n] for n in names], color="#9E9E9E", height=0.5)
    ax_tot.set_yticks(y)
    ax_tot.set_yticklabels([])
    ax_tot.invert_xaxis()
    ax_tot.set_xlabel("Set size", fontsize=11)
    ax_tot.set_ylim(-0.6, n_s - 0.4)
    ax_tot.spines[["top", "left", "right"]].set_visible(False)
    fig.savefig(outfile, bbox_inches="tight", dpi=300)
    plt.close(fig)

# ================================================= run per gene list =========
for GENES_XLSX, QUERY_LABEL in GENE_LISTS.items():
    tag = os.path.splitext(os.path.basename(GENES_XLSX))[0]
    print(f"\n######## {GENES_XLSX}  ({QUERY_LABEL}) ########")
    query = load_gene_list(GENES_XLSX)

    # ---- overlap table -------------------------------------------------------
    print("\n=== input gene list  ∩  peak-annotated genes ===")
    print(f"{'peak_filter':<13} {'sex':<7} {'n_peak_genes':>12} {'overlap':>8} "
          f"{'% of input':>11} {'% of peaks':>11}")

    rows = []
    overlap_lists = {}                                # (pf, sex) -> sorted gene list
    for pf in ["allPeaks", "promoterOnly"]:
        for sex, pg in peak_gene_sets(pf).items():
            ov   = sorted(query & pg)
            n_ov = len(ov)
            pct_in  = round(100 * n_ov / len(query), 1)
            pct_pk  = round(100 * n_ov / len(pg), 1) if pg else np.nan
            print(f"{pf:<13} {sex:<7} {len(pg):>12} {n_ov:>8} "
                  f"{pct_in:>11} {pct_pk:>11}")
            rows.append({"peak_filter": pf, "sex": sex,
                         "n_input_genes": len(query), "n_peak_genes": len(pg),
                         "n_overlap": n_ov,
                         "pct_of_input": pct_in, "pct_of_peaks": pct_pk})
            overlap_lists[(pf, sex)] = ov

    summary = pd.DataFrame(rows)

    # input genes with no peak at all (allPeaks/all)
    no_peak = sorted(query - peak_gene_sets("allPeaks")["all"])
    print(f"\ninput genes with no annotated peak: {len(no_peak)}")

    OUT_XLSX = f"{OUTDIR}/genelist_peak_overlap_{tag}.xlsx"
    with pd.ExcelWriter(OUT_XLSX, engine="openpyxl") as xw:
        summary.to_excel(xw, sheet_name="summary", index=False)
        for (pf, sex), genes in overlap_lists.items():
            sheet = f"{pf}_{sex}"[:31]                 # Excel sheet-name cap
            pd.DataFrame({"gene": genes}).to_excel(xw, sheet_name=sheet, index=False)
        pd.DataFrame({"gene": no_peak}).to_excel(xw, sheet_name="input_no_peak",
                                                 index=False)
    print(f"\n-> {OUT_XLSX}  (summary + {len(overlap_lists)} overlap sheets + input_no_peak)")

    # ---- Venn per peak filter ------------------------------------------------
    for pf in ["allPeaks", "promoterOnly"]:
        male, female = peak_sets(pf)
        fig, ax = plt.subplots(figsize=(7, 7))
        venn = venn3(
            [query, male, female],
            set_labels=[QUERY_LABEL, "Male peaks", "Female peaks"],
            set_colors=COLORS,
            alpha=0.65,
            ax=ax,
            **_LAYOUT,
        )
        # ---- style ----
        for patch in venn.patches:
            if patch is not None:
                patch.set_edgecolor("black")
                patch.set_linewidth(1.8)
        for txt in (venn.set_labels or []):
            if txt:
                txt.set_fontsize(14)
                txt.set_fontweight("bold")
        for txt in (venn.subset_labels or []):
            if txt:
                txt.set_fontsize(12)
                txt.set_fontweight("bold")
        ax.set_title(f"{pf} ({'proportional' if PROPORTIONAL else 'unweighted'})",
                     fontsize=16, fontweight="bold")
        plt.tight_layout()
        out = f"{OUTDIR}/venn_proportional_{tag}_{pf}.pdf"
        plt.savefig(out, dpi=600, bbox_inches="tight")
        plt.close(fig)
        print(f"-> {out}   query={len(query)}  male={len(male)}  female={len(female)}  "
              f"| query&male={len(query & male)}  query&female={len(query & female)}  "
              f"query&both={len(query & male & female)}  "
              f"query only={len(query - male - female)}")

    # ---- UpSet across all sets -----------------------------------------------
    am, af = peak_sets("allPeaks")
    pm, pf_ = peak_sets("promoterOnly")
    contents = {
        QUERY_LABEL:       query,
        "allPeaks male":   am,
        "allPeaks female": af,
        "Promoter male":   pm,
        "Promoter female": pf_,
    }

    UPSET_OUT = f"{OUTDIR}/upset_genelist_{tag}.pdf"
    try:
        from upsetplot import from_contents, UpSet
        data = from_contents(contents)
        up = UpSet(
            data,
            subset_size="count",
            show_counts=True,
            sort_by="cardinality",
            sort_categories_by=None,
        )
        fig = plt.figure(figsize=(12, 7))
        up.plot(fig=fig)
        plt.savefig(UPSET_OUT, dpi=300, bbox_inches="tight")
        plt.close(fig)
        print(f"-> {UPSET_OUT}  (upsetplot)")
    except Exception as e:
        print(f"upsetplot failed ({type(e).__name__}: {e}) -- matplotlib fallback")
        upset_fallback(contents, UPSET_OUT, title=f"{QUERY_LABEL} vs peak sets")
        print(f"-> {UPSET_OUT}  (fallback)")
