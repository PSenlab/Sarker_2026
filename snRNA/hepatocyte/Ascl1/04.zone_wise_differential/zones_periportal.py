#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import pandas as pd
 
# ------------------------------------------------------------
# 1. Load both
# ------------------------------------------------------------
mast = pd.read_csv("Ascl1_MAST_periportal_up.csv")


# In[3]:


import gseapy as gp

# detect the gene column ONCE, up front
gcol = "gene" if "gene" in mast.columns else "primerid"
print("Using gene column:", gcol)

mast_up = mast[(mast["logFC"] > 0.1) & (mast["p_value"] < 0.05)]
up_mast = set(mast_up[gcol])

gene_list = [str(g).upper() for g in sorted(up_mast)]
print(f"Input: {len(gene_list)} genes")


# In[5]:


enr = gp.enrichr(
    gene_list=gene_list,
    gene_sets="Reactome_Pathways_2024",
    organism="human",
    outdir=None,
)

res = enr.results.copy()
ov = res["Overlap"].str.split("/", expand=True).astype(int)
res["n_overlap"] = ov[0]
res["set_size"]  = ov[1]
res["Gene %"]    = (100 * res["n_overlap"] / res["set_size"]).round(1)

cols = ["Term", "n_overlap", "set_size", "Gene %",
        "P-value", "Adjusted P-value", "Odds Ratio", "Combined Score", "Genes"]
cols = [c for c in cols if c in res.columns]
res = res[cols].sort_values("Adjusted P-value")

res.to_csv("enrichr_reactome_mast_up_periportal.csv", index=False)
print(res.drop(columns="Genes").head(25).to_string(index=False))


# In[6]:


#!/usr/bin/env python3
# ============================================================
# Enrichr Reactome dotplot - UP in Ascl1+ (Female Hepatocytes)
#   Filter:  n_overlap > 5  AND  Adjusted P-value < 0.05
#   Sort:    Odds Ratio (desc), plot top 30, highest OR at TOP
#   x = Odds Ratio | size = -log10(Adj P) | hue = Gene %
# ============================================================

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

plt.rcParams["font.family"] = "Arial"

PATH    = "enrichr_reactome_mast_up_periportal.csv"
OUTFILE = "enrichr_reactome_mast_up_periportal_dotplot.pdf"
TOP_N   = 30
MIN_SIZE, MAX_SIZE = 60, 400

# ------------------------------------------------------------
# load + clean
# ------------------------------------------------------------
df = pd.read_csv(PATH)
df.columns = df.columns.str.strip()

# rebuild n_overlap / set_size / Gene % if not already present
if "n_overlap" not in df.columns:
    ov = df["Overlap"].str.split("/", expand=True).astype(int)
    df["n_overlap"], df["set_size"] = ov[0], ov[1]
if "Gene %" not in df.columns:
    df["Gene %"] = (100 * df["n_overlap"] / df["set_size"]).round(1)

for c in ["Adjusted P-value", "Odds Ratio", "Gene %"]:
    df[c] = pd.to_numeric(df[c], errors="coerce")

# ------------------------------------------------------------
# filter + select
# ------------------------------------------------------------
sub = df[(df["n_overlap"] > 5) & (df["Adjusted P-value"] < 0.05)].copy()
sub = sub.sort_values("Odds Ratio", ascending=True).reset_index(drop=True)
sub["Term"] = pd.Categorical(sub["Term"], categories=sub["Term"], ordered=True)

# dot size = -log10(Adjusted P-value)
sub["pval_size"] = -np.log10(sub["Adjusted P-value"].replace(0, 1e-300))
pmin, pmax = sub["pval_size"].min(), sub["pval_size"].max()
rng = (pmax - pmin) if pmax > pmin else 1.0
sub["size_scaled"] = MIN_SIZE + (sub["pval_size"] - pmin) / rng * (MAX_SIZE - MIN_SIZE)

# y-axis order: highest odds ratio FIRST in dataframe; invert_yaxis() puts it on top
sub = sub.sort_values("Odds Ratio", ascending=False).reset_index(drop=True)
sub["Term"] = pd.Categorical(sub["Term"], categories=sub["Term"], ordered=True)

# ------------------------------------------------------------
# plot
# ------------------------------------------------------------
fig, ax = plt.subplots(figsize=(3, 9))
sns.scatterplot(
    data=sub, x="Odds Ratio", y="Term",
    size="size_scaled", sizes=(MIN_SIZE, MAX_SIZE),
    hue="Gene %", palette="RdYlBu_r",
    edgecolor="black", linewidth=0.6, ax=ax, legend="brief",
)
ax.invert_yaxis()   # highest odds ratio at the TOP

ax.set_xlabel("Odds Ratio", fontsize=12, weight="bold")
ax.set_ylabel("Reactome pathway", fontsize=12, weight="bold")
ax.set_title("Reactome enrichment - Up in Ascl1+ (Periportal)\n"
             "gene count > 5, adj p < 0.05, top 30 by odds ratio",
             fontsize=13, weight="bold")
ax.margins(x=0.08)

norm = plt.Normalize(sub["Gene %"].min(), sub["Gene %"].max())
sm = plt.cm.ScalarMappable(cmap="RdYlBu_r", norm=norm); sm.set_array([])
cbar = fig.colorbar(sm, ax=ax, pad=0.02)
cbar.set_label("Gene %", fontsize=11, weight="bold")

ax.legend(title="-log10(adj p)", bbox_to_anchor=(1.32, 1),
          loc="upper left", frameon=True, fontsize=9, title_fontsize=10)

fig.tight_layout()
fig.savefig(OUTFILE, dpi=300, bbox_inches="tight")
print(f"saved -> {OUTFILE}")
plt.show()


# In[7]:


import numpy as np
import pandas as pd
 
# ------------------------------------------------------------
# 1. Load both
# ------------------------------------------------------------
mast = pd.read_csv("Ascl1_MAST_periportal_down.csv")


# In[8]:


import gseapy as gp

# detect the gene column ONCE, up front
gcol = "gene" if "gene" in mast.columns else "primerid"
print("Using gene column:", gcol)

mast_up = mast[(mast["logFC"] < -0.1) & (mast["p_value"] < 0.05)]
up_mast = set(mast_up[gcol])

gene_list = [str(g).upper() for g in sorted(up_mast)]
print(f"Input: {len(gene_list)} genes")


# In[9]:


enr = gp.enrichr(
    gene_list=gene_list,
    gene_sets="Reactome_Pathways_2024",
    organism="human",
    outdir=None,
)

res = enr.results.copy()
ov = res["Overlap"].str.split("/", expand=True).astype(int)
res["n_overlap"] = ov[0]
res["set_size"]  = ov[1]
res["Gene %"]    = (100 * res["n_overlap"] / res["set_size"]).round(1)

cols = ["Term", "n_overlap", "set_size", "Gene %",
        "P-value", "Adjusted P-value", "Odds Ratio", "Combined Score", "Genes"]
cols = [c for c in cols if c in res.columns]
res = res[cols].sort_values("Adjusted P-value")

res.to_csv("enrichr_reactome_mast_down_periportal.csv", index=False)
print(res.drop(columns="Genes").head(25).to_string(index=False))


# In[ ]:




