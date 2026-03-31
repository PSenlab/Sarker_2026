# ============================================================
# 10. Senescence scoring — Classical markers, SHGS, SenMayo,
#     SenNet module scores + heatmaps
# ============================================================
print("\n" + "=" * 60)
print("ANALYSIS 10: Senescence scoring (SHGS, SenMayo, SenNet)")
print("=" * 60)

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

# ------------------------------------------------------------------
# 10a. Define senescence gene sets
# ------------------------------------------------------------------

# --- Classical senescence markers (individual gene scoring) ---
CLASSICAL_SENESCENCE_GENES = [
    "Cdkn1a", "Cdkn2a", "Trp53", "Rb1", "Mdm2",
    "Il6", "Ccl2", "Cxcl1", "Serpine1", "Igfbp7",
]

# --- SHGS: Senescent Hepatocyte Gene Signature ---
# (Du et al. 2025, Nat Commun — 100-gene overlap of in vitro
#  palbociclib-treated Huh7 + in vivo NRASG12V hepatocytes,
#  mouse orthologs from the Zenodo senescence_signatures.rds)
SHGS_GENES = [
    "Gm29685", "Gm19510", "Cxcl10", "4930546K05Rik", "Tox2",
    "Gm15418", "Glis3", "Rsad2", "Ly6m", "Osmr",
    "Gm48682", "Socs1", "Nid1", "Cib4", "Tgtp2",
    "Gbp2", "Hcar2", "Cd14", "Gm33782", "Anxa2",
    "Dusp5", "Gbp8", "Serpine1", "Gm30692", "Gm16685",
    "Tnfrsf12a", "Cdkn1a", "Hilpda", "Il7", "S1pr3",
    "2310034O05Rik", "Rasal1", "Tmem154", "Ralgds", "Bco1",
    "Relb", "Gldn", "Tgtp1", "Gm42765", "Cmpk2",
    "Cxcl9", "Cacng4", "Ier3", "Sorcs2", "Ifit3b",
    "Ubd", "Numbl", "Shc2", "Dchs1", "Gadd45b",
    "Nrg1", "Sult1e1", "Col4a1", "Gm29282", "Adcy1",
    "Clcf1", "Ifit2", "Gm48341", "Gbp9", "Fhl3",
    "Tgm2", "Gm37359", "Prr5l", "Tubb6", "Ly6a",
    "Gm17746", "Gm12185", "Cxcl13", "Oasl2", "Slc7a1",
    "2010003K11Rik", "Zfp462", "Olfr56", "Gm39348", "H2-Q1",
    "Igdcc4", "4933401D09Rik", "Casp12", "Srxn1", "Gucy2c",
    "Igfbp1", "Gldnos", "Scn2a", "Gpc6", "Tlr2",
    "Mt2", "Gipc2", "Rragd", "Mab21l3", "5330417C22Rik",
    "Casp4", "Arid5a", "Rims4", "Coq8b", "9930111J21Rik2",
    "Gm45774", "Ifit3", "BC016579", "Atp6v0d2", "Ccdc120",
    "Lpl", "Cacna1b", "Ly6d", "Spdef", "Ier5",
    "Mt1", "Rnf144a", "Alpk1", "Gm11454", "Slc4a9",
    "Slc1a4", "9030624G23Rik", "Olfr16", "Akr1b7", "Ddr1",
    "2700038G22Rik", "Gm15567", "Lepr", "Inhbb", "Lcn2",
    "Gm26788", "Syt12", "Dusp8", "Slc37a1", "Gm33699",
    "Asb11", "Tnfaip3", "Ifit1", "Epha4", "Serpina3n",
    "Gm36037", "Cerk", "Rp1", "Hbegf", "Gm50383",
    "Spats2l", "Gbp6", "Scn8a", "Npas2", "Klf6",
    "Ifi44", "Sox9", "Ffar4", "Gm20125", "Tagln2",
    "Trim6", "Gm26714", "Cd274", "Ccdc180", "Dnah1",
    "1200007C13Rik", "Sytl5", "Crybb3", "Ets2", "Chka",
    "Pde6h", "Irak3", "Gm41609", "Il1rn", "Prag1",
    "Arhgef2", "Tnfrsf23", "Itpkc", "Nlrp9c",
]

# --- SenMayo: 125-gene SASP/senescence panel ---
# (Saul et al. 2022, Nat Commun — mouse gene names)
SENMAYO_GENES = [
    "Acvr1b", "Ang", "Angpt1", "Angptl4", "Areg", "Axl", "Bex3",
    "Bmp2", "Bmp6", "C3", "Ccl1", "Ccl2", "Ccl20", "Ccl24",
    "Ccl26", "Ccl3", "Ccl4", "Ccl5", "Ccl7", "Ccl8", "Cd55",
    "Cd9", "Csf1", "Csf2", "Csf2rb", "Ctnnb1", "Ctsb", "Cxcl1",
    "Cxcl10", "Cxcl12", "Cxcl16", "Cxcl2", "Cxcl3", "Cxcr2",
    "Dkk1", "Edn1", "Egf", "Egfr", "Ereg", "Esm1", "Ets2",
    "Fas", "Fgf1", "Fgf2", "Fgf7", "Gdf15", "Gem", "Gmfg",
    "Hgf", "Hmgb1", "Icam1", "Icam3", "Igf1", "Igfbp1", "Igfbp2",
    "Igfbp3", "Igfbp4", "Igfbp5", "Igfbp6", "Igfbp7", "Il10",
    "Il13", "Il15", "Il18", "Il1a", "Il1b", "Il2", "Il6", "Il6st",
    "Il7", "Inha", "Iqgap2", "Itga2", "Itpka", "Jun", "Kitl",
    "Lcp1", "Mif", "Mmp10", "Mmp12", "Mmp13", "Mmp14", "Mmp2",
    "Mmp3", "Mmp9", "Nap1l4", "Nrg1", "Pappa", "Pecam1", "Pgf",
    "Pigf", "Plat", "Plau", "Plaur", "Ptbp1", "Ptger2", "Ptges",
    "Rps6ka5", "Scamp4", "Selplg", "Sema3f", "Serpinb4", "Serpine1",
    "Serpine2", "Spp1", "Spx", "Timp2", "Tnf", "Tnfrsf11b",
    "Tnfrsf1a", "Tnfrsf1b", "Tubgcp2", "Vegfa", "Vegfc", "Vgf",
    "Wnt16", "Wnt2",
]

# --- SenNet: Core senescence biomarkers ---
# (NIH SenNet consortium 2024 recommendations)
SENNET_GENES = [
    # Cell Cycle Arrest
    "Cdkn2a", "Cdkn1a", "Trp53", "Rb1",
    # DNA Damage Response
    "H2afx", "Trp53bp1", "Atm", "Atr",
    # Nuclear Changes
    "Lmnb1",
    # Anti-Apoptotic (pro-survival)
    "Bcl2", "Bcl2l1", "Bcl2l2", "Serpinb2",
    # SASP Core
    "Il1a", "Il1b", "Il6", "Cxcl1", "Ccl2", "Mmp3", "Mmp9",
    # Surface Markers
    "Dpp4", "Cd36", "Icam1", "Plaur", "Tnfrsf10d",
    # Lysosomal
    "Glb1",
]

# ------------------------------------------------------------------
# 10b. Classical senescence markers — per-gene scoring + z-scored
#      heatmap per sex (celltype2 as columns)
# ------------------------------------------------------------------
print("\n--- 10b: Classical senescence marker heatmaps ---")

# Score each gene individually
for gene in CLASSICAL_SENESCENCE_GENES:
    score_name = f"{gene}_score"
    if gene in adata.var_names:
        sc.tl.score_genes(adata, gene_list=[gene], score_name=score_name)
    elif adata.raw is not None and gene in adata.raw.var_names:
        sc.tl.score_genes(adata, gene_list=[gene], score_name=score_name, use_raw=True)
    else:
        print(f"  ⚠️ {gene} not found — skipping")

# Build list of score columns that exist
valid_genes = [
    "Cdkn1a", "Cdkn2a", "Trp53", "Rb1", "Mdm2",
    "Il6", "Ccl2", "Cxcl1", "Serpine1", "Igfbp7",
]
score_cols = [f"{g}_score" for g in valid_genes if f"{g}_score" in adata.obs.columns]
scored_genes = [col.replace("_score", "") for col in score_cols]

for sex in adata.obs["sex"].unique():
    print(f"  Processing: {sex}")
    df_sex = adata.obs[adata.obs["sex"] == sex]

    # Mean score per gene across celltype2
    df = df_sex.groupby("celltype2")[score_cols].mean().T  # genes × celltypes
    # Z-score row-wise
    df_zscore = df.sub(df.mean(axis=1), axis=0).div(df.std(axis=1), axis=0)

    # Reindex to gene order
    df_zscore["Gene"] = df_zscore.index.str.replace("_score", "", regex=False)
    df_zscore = df_zscore.set_index("Gene").loc[scored_genes]

    # Plot
    plt.figure(figsize=(10, 0.4 * len(df_zscore)))
    sns.heatmap(
        df_zscore, cmap="coolwarm", center=0, linewidths=0.05,
        xticklabels=True, yticklabels=True,
        cbar_kws={"label": "Z-score"},
    )
    plt.title(f"Classical Senescence Markers – {sex}")
    plt.xlabel("Cell Type")
    plt.ylabel("Gene")
    plt.tight_layout()

    fname = f"{FIGURES_DIR}/senescence_score_heatmap_{sex}.pdf"
    plt.savefig(fname, dpi=300)
    plt.close()
    print(f"  Saved: {fname}")


# ------------------------------------------------------------------
# 10c. Helper — score a gene set & produce side-by-side male/female
#      heatmap (celltype2 × age, mean module score)
# ------------------------------------------------------------------
def score_and_plot_module(adata, gene_list, score_name, title_label):
    """
    Score cells with sc.tl.score_genes, then create a side-by-side
    male / female heatmap of mean module score (celltype2 × age).
    Returns the score column name added to adata.obs.
    """
    # Filter to genes present
    present = [g for g in gene_list if g in adata.var_names]
    use_raw = False
    if not present and adata.raw is not None:
        present = [g for g in gene_list if g in adata.raw.var_names]
        use_raw = bool(present)
    print(f"  {title_label}: {len(present)}/{len(gene_list)} genes found")

    if not present:
        print(f"  ⚠️ No genes found for {title_label} — skipping")
        return None

    sc.tl.score_genes(adata, gene_list=present, score_name=score_name, use_raw=use_raw)

    # Pivot: celltype2 (rows) × age (columns), mean score — per sex
    dfs = {}
    for sex in SEX_ORDER:
        sub = adata.obs[adata.obs[SEX_COL] == sex]
        pivot = sub.pivot_table(
            index=CELLTYPE2_COL, columns=AGE_COL,
            values=score_name, aggfunc="mean",
        )
        # Reorder
        row_order = [ct for ct in HEP_ORDER if ct in pivot.index]
        col_order = [a for a in AGE_ORDER if a in pivot.columns]
        pivot = pivot.loc[row_order, col_order]
        dfs[sex] = pivot

    df_male = dfs.get("male", pd.DataFrame())
    df_female = dfs.get("female", pd.DataFrame())

    if df_male.empty and df_female.empty:
        print(f"  ⚠️ No data for {title_label} — skipping plot")
        return score_name

    # Determine shared color range
    all_vals = pd.concat(
        [df_male.stack(), df_female.stack()], ignore_index=True
    ).dropna()
    vmin = float(all_vals.min())
    vmax = float(all_vals.max())

    n_rows = max(len(df_male), len(df_female), 1)
    fig, axes = plt.subplots(
        1, 2, figsize=(18, max(4, 0.4 * n_rows)), sharey=True
    )

    # Male
    if not df_male.empty:
        sns.heatmap(
            df_male, ax=axes[0], annot=True, fmt=".3f",
            cmap="coolwarm", vmin=vmin, vmax=vmax,
            cbar_kws={"label": "Avg Module Score"},
        )
    axes[0].set_title("Male")
    axes[0].set_xlabel("Age Group")
    axes[0].set_ylabel("Cell Type")

    # Female
    if not df_female.empty:
        sns.heatmap(
            df_female, ax=axes[1], annot=True, fmt=".3f",
            cmap="coolwarm", vmin=vmin, vmax=vmax,
            cbar_kws={"label": "Avg Module Score"},
        )
    axes[1].set_title("Female")
    axes[1].set_xlabel("Age Group")
    axes[1].set_ylabel("")

    plt.suptitle(f"module_score_{title_label}", fontsize=14)
    plt.tight_layout(rect=[0, 0, 1, 0.95])

    fname_pdf = f"{FIGURES_DIR}/module_score_{title_label}.pdf"
    fname_png = f"{FIGURES_DIR}/module_score_{title_label}.png"
    fig.savefig(fname_pdf, dpi=300, bbox_inches="tight")
    fig.savefig(fname_png, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {fname_pdf}")
    print(f"  Saved: {fname_png}")

    # Save underlying data
    for sex, df in dfs.items():
        out = f"{RESULTS_DIR}/module_score_{title_label}_{sex}.csv"
        df.to_csv(out)
        print(f"  Saved: {out}")

    return score_name


# ------------------------------------------------------------------
# 10d. Run module scoring for SHGS, SenMayo, SenNet
# ------------------------------------------------------------------
print("\n--- 10d: SHGS module score ---")
score_and_plot_module(adata, SHGS_GENES, "SHGS_score", "SHGS")

print("\n--- 10e: SenMayo module score ---")
score_and_plot_module(adata, SENMAYO_GENES, "SenMayo_score", "SenMayo")

print("\n--- 10f: SenNet module score ---")
score_and_plot_module(adata, SENNET_GENES, "SenNet_score", "SenNet")


# ------------------------------------------------------------------
# 10g. Combined summary heatmap — all three signatures, male vs
#      female, mean score per celltype2 (collapsed across age)
# ------------------------------------------------------------------
print("\n--- 10g: Combined senescence summary heatmap ---")

summary_rows = []
for sig_name, score_col in [
    ("SHGS", "SHGS_score"),
    ("SenMayo", "SenMayo_score"),
    ("SenNet", "SenNet_score"),
]:
    if score_col not in adata.obs.columns:
        continue
    for sex in SEX_ORDER:
        sub = adata.obs[adata.obs[SEX_COL] == sex]
        means = sub.groupby(CELLTYPE2_COL)[score_col].mean()
        for ct, val in means.items():
            summary_rows.append({
                "signature": sig_name,
                "sex": sex,
                CELLTYPE2_COL: ct,
                "mean_score": val,
            })

if summary_rows:
    df_summary = pd.DataFrame(summary_rows)
    out = f"{RESULTS_DIR}/senescence_module_scores_summary.csv"
    df_summary.to_csv(out, index=False)
    print(f"  Saved: {out}")

print("\n" + "=" * 60)
print("ANALYSIS 10 COMPLETE — Senescence scoring done")
print(f"  Figures → {FIGURES_DIR}/")
print(f"  Tables  → {RESULTS_DIR}/")
print("=" * 60)