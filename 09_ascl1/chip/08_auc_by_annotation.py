#!/usr/bin/env python
# coding: utf-8
# ============================================================================
# ASCL1 AUC by peak-annotation category, per sex (violins)
# Figure 7, panels b (male) and c (female), AUC sum on a log axis
#
# Input : Ascl1_F/M_narrowpeaks_q0.001_AUC.txt  (from 05_bwtool_rpkm.sh)
#         chipseeker_out/all_peaks.csv          (from 07_chipseeker_annotation.R)
# Output: chipseeker_out/auc_violin_<sex>_<sum|mean>.pdf
#         chipseeker_out/auc_violin_combined_<sum|mean>.pdf
#         chipseeker_out/auc_by_annotation_summary.csv
#
# Scores:
#   sum  -- per-base signal summed over the peak (shown in the figure)
#   mean -- sum / size, width-normalized
#
# Join: the AUC files are 0-based half-open (BED, '#chrom' column);
# all_peaks.csv is 1-based (ChIPseeker GRanges, BED read with
# starts.in.df.are.0based = TRUE), so AUC start + 1 = annotation start.
#
# Violins are sorted by median, with jittered points and an inner quartile
# box; categories with < MIN_VIOLIN peaks are drawn as points only. No
# statistical test is applied; medians and quartiles go to the summary CSV.
# ============================================================================

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager

for _f in ["Arial", "Liberation Sans", "Nimbus Sans", "DejaVu Sans"]:
    if _f in {f.name for f in font_manager.fontManager.ttflist}:
        plt.rcParams["font.family"] = _f
        break

INDIR = "chipseeker_out"
SEX_COLS = {"male": "#89d9e1", "female": "#fb9a99"}
GREY = "#4d4d4d"

AUC_FILES = {
    "female": "Ascl1_F_narrowpeaks_q0.001_AUC.txt",
    "male":   "Ascl1_M_narrowpeaks_q0.001_AUC.txt",
}

AUC_CHR, AUC_START, AUC_END = "#chrom", "start", "end"
AUC_IS_0BASED = True
MIN_VIOLIN = 15          # below this, draw points only, no violin
LOG_Y      = True        # AUC sum spans ~400 to ~6000; log reads better

# ---------------------------------------------------------------- read/join --
def read_auc(path):
    return pd.read_csv(path, sep=None, engine="python")

def coordkey(df, c, s, e, shift=0):
    return (df[c].astype(str) + ":" +
            (df[s].astype(int) + shift).astype(str) + "-" +
            df[e].astype(int).astype(str))

anno = pd.read_csv(os.path.join(INDIR, "all_peaks.csv"))

frames = []
for sex, path in AUC_FILES.items():
    d = read_auc(path)
    # keep both score columns; derive mean if the file lacks it
    if "mean" not in d.columns and {"sum", "size"}.issubset(d.columns):
        d["mean"] = d["sum"] / d["size"]
    d = d.rename(columns={"sum": "AUC_sum", "mean": "AUC_mean"})
    a = anno[anno["sex"] == sex].copy()

    d["_key"] = coordkey(d, AUC_CHR, AUC_START, AUC_END,
                         shift=1 if AUC_IS_0BASED else 0)
    a["_key"] = coordkey(a, "seqnames", "start", "end")

    keep_cols = ["_key", "AUC_sum"] + (["AUC_mean"] if "AUC_mean" in d else [])
    m = a.merge(d[keep_cols], on="_key", how="inner")
    print(f"[{sex}] {len(a)} annotated, {len(d)} AUC, {len(m)} matched "
          f"({100*len(m)/max(len(a),1):.0f}%)")
    m["sex"] = sex
    frames.append(m)

df = pd.concat(frames, ignore_index=True)

def broad(a):
    a = str(a)
    for pref, lab in [("Promoter", "Promoter"), ("5' UTR", "5' UTR"),
                      ("3' UTR", "3' UTR"), ("Exon", "Exon"),
                      ("Intron", "Intron"), ("Downstream", "Downstream")]:
        if a.startswith(pref):
            return lab
    if "Distal Intergenic" in a:
        return "Distal Intergenic"
    return a

df["region"] = df["annotation"].map(broad)

# ---------------------------------------------------------------- plotting ---
def violin_panel(sub, value_col, title, ylabel, outfile, fill):
    d = sub[["region", value_col]].dropna().rename(columns={value_col: "v"})
    if LOG_Y:
        d = d[d["v"] > 0]        # log needs positive
    # order categories by median, high to low -- makes the contrast read L->R
    med = d.groupby("region")["v"].median().sort_values(ascending=False)
    cats = list(med.index)
    if not cats:
        print(f"[{title}] nothing to plot"); return

    fig, ax = plt.subplots(figsize=(max(7, 1.1 * len(cats) + 2), 5))
    rng = np.random.default_rng(0)

    for i, c in enumerate(cats, 1):
        vals = d.loc[d["region"] == c, "v"].values
        n = len(vals)

        if n >= MIN_VIOLIN:
            parts = ax.violinplot([vals], positions=[i], widths=0.8,
                                   showmedians=False, showextrema=False)
            for pc in parts["bodies"]:
                pc.set_facecolor(fill)
                pc.set_edgecolor(GREY)
                pc.set_alpha(0.55)
                pc.set_linewidth(0.8)
            # inner quartile box
            q1, md, q3 = np.percentile(vals, [25, 50, 75])
            ax.add_patch(plt.Rectangle((i - 0.05, q1), 0.10, q3 - q1,
                                       facecolor=GREY, edgecolor="none",
                                       zorder=4))
            ax.hlines(md, i - 0.11, i + 0.11, color="white",
                      linewidth=1.6, zorder=5)

        # jittered points: all if few, a capped sample if many (keeps PDF light)
        show = vals if n <= 400 else rng.choice(vals, 400, replace=False)
        jit = rng.uniform(-0.18, 0.18, len(show))
        ax.scatter(np.full(len(show), i) + jit, show, s=6,
                   color=GREY, alpha=0.35, linewidth=0, zorder=3)

        # median marker in the sex colour
        ax.scatter([i], [np.median(vals)], s=45, color=fill,
                   edgecolor="black", linewidth=0.8, zorder=6)

    ax.set_xticks(range(1, len(cats) + 1))
    ax.set_xticklabels(cats, rotation=30, ha="right")
    ax.set_ylabel(ylabel, fontsize=12, weight="bold")
    ax.set_xlabel("Peak annotation category", fontsize=12, weight="bold")
    ax.set_title(title, fontsize=11, weight="bold")
    if LOG_Y:
        ax.set_yscale("log")
    ax.spines[["top", "right"]].set_visible(False)

    # n per category, above each violin
    ytop = ax.get_ylim()[1]
    for i, c in enumerate(cats, 1):
        n = int((d["region"] == c).sum())
        ax.text(i, ytop, f"n={n}", ha="center", va="bottom", fontsize=7.5,
                color=GREY)

    fig.tight_layout()
    fig.savefig(outfile, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"saved -> {outfile}")

# --- draw: sum and mean (width-normalized) for each sex ----------------------
score_specs = [
    ("AUC_sum",  "Ascl1 AUC (sum)",  "sum"),
    ("AUC_mean", "Ascl1 AUC (mean, width-normalized)", "mean"),
]

for sex in ["female", "male"]:
    sub = df[df["sex"] == sex]
    if not len(sub):
        continue
    for col, ylabel, tag in score_specs:
        if col not in sub.columns:
            continue
        violin_panel(sub, col,
                     f"Ascl1 AUC by annotation ({sex}, {tag})",
                     ylabel,
                     os.path.join(INDIR, f"auc_violin_{sex}_{tag}.pdf"),
                     SEX_COLS[sex])

# combined
for col, ylabel, tag in score_specs:
    if col in df.columns:
        violin_panel(df, col, f"Ascl1 AUC by annotation (both sexes, {tag})",
                     ylabel, os.path.join(INDIR, f"auc_violin_combined_{tag}.pdf"),
                     "#7fb3bd")

# ---------------------------------------------------------------- summary ----
agg = {"n": ("AUC_sum", "count"),
       "sum_median": ("AUC_sum", "median"),
       "sum_q25": ("AUC_sum", lambda x: x.quantile(.25)),
       "sum_q75": ("AUC_sum", lambda x: x.quantile(.75))}
if "AUC_mean" in df.columns:
    agg.update({"mean_median": ("AUC_mean", "median"),
                "mean_q25": ("AUC_mean", lambda x: x.quantile(.25)),
                "mean_q75": ("AUC_mean", lambda x: x.quantile(.75))})
summ = (df.groupby(["sex", "region"]).agg(**agg).reset_index()
          .sort_values(["sex", "sum_median"], ascending=[True, False]))
summ.to_csv(os.path.join(INDIR, "auc_by_annotation_summary.csv"), index=False)
print("\n" + summ.to_string(index=False))
