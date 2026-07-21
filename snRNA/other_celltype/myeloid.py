import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import f_oneway
from statsmodels.stats.multitest import multipletests
plt.rcParams['font.family'] = 'Arial'

# =========================================================================
# 0. objects / keys
#   `adata`     = myeloid subset (Kupffer + MoMFs), scVI-subclustered + annotated
#   `parent`    = full 295,364-cell liver atlas (composition denominator)
# =========================================================================
parent      = adata_parent          # <-- full liver atlas (all cells)
GROUP_COL   = 'celltype_myeloid'
SAMPLE_KEY  = 'sample'
age_groups  = ['young', 'mid_age', 'old', 'pre_geriatric', 'geriatric']
sex_order   = ['female', 'male']
disp        = {'young':'young','mid_age':'mid-age','old':'old',
               'pre_geriatric':'pre-geriatric','geriatric':'geriatric'}

# clean label + fixed publication order
adata.obs[GROUP_COL] = (adata.obs[GROUP_COL].astype(str)
                        .replace({'Kupffer cycling': 'Kupffer_cycling'}).astype('category'))
ORDER = ['Kupffer', 'Kupffer_cycling', 'LAM', 'MoMF', 'cDC1', 'pDC', 'Neutrophil']
assert set(adata.obs[GROUP_COL].unique()) == set(ORDER), \
    f"label mismatch: {set(adata.obs[GROUP_COL].unique())} vs {set(ORDER)}"
adata.obs[GROUP_COL] = adata.obs[GROUP_COL].cat.reorder_categories(ORDER)

# =========================================================================
# 1. marker dot plot (annotation support, variance-scaled)
# =========================================================================
marker_genes = {
    "Kupffer":         ["Clec4f","Vsig4","Timd4","Cd5l","Marco","C1qa"],
    "Kupffer_cycling": ["Mki67","Top2a","Cenpf"],
    "LAM":             ["Gpnmb","Trem2","Cd9","Arhgap22","Sirpb1a","Cd74","Pparg"],
    "MoMF":            ["Ccr2","Cx3cr1","Ly6c2","Chil3","Plac8","Fn1"],
    "cDC1":            ["Xcr1","Clec9a","Batf3","Cadm1"],
    "pDC":             ["Siglech","Bst2","Il3ra","Tcf4"],
    "Neutrophil":      ["S100a8","S100a9","Retnlg","Csf3r"],
}
present = {k:[g for g in v if g in adata.var_names] for k,v in marker_genes.items()}
present = {k:v for k,v in present.items() if v}
sc.pl.dotplot(adata, present, groupby=GROUP_COL,
              standard_scale='var', dendrogram=False, show=False)
plt.savefig("myeloid_RNA_dotplot_stdscale.pdf", bbox_inches="tight")
plt.show()

# =========================================================================
# 2. composition:  % of ALL liver cells per sample  (denominator = parent)
# =========================================================================
total_per_sample = parent.obs.groupby(SAMPLE_KEY, observed=True).size().rename('total')

counts_wide = (adata.obs.groupby([SAMPLE_KEY, GROUP_COL], observed=True)
                 .size().unstack(fill_value=0).reindex(columns=ORDER, fill_value=0)
                 .reindex(index=total_per_sample.index, fill_value=0))
comp = counts_wide.reset_index().melt(id_vars=SAMPLE_KEY, var_name=GROUP_COL, value_name='n')
comp = comp.merge(total_per_sample, on=SAMPLE_KEY)
comp['pct'] = 100 * comp['n'] / comp['total']
comp = comp.merge(parent.obs[[SAMPLE_KEY,'age','sex']].drop_duplicates(), on=SAMPLE_KEY)
comp['age'] = pd.Categorical(comp['age'], age_groups, ordered=True)

# =========================================================================
# 3. boxplot: subcluster x age, faceted by sex, with FDR stars
# =========================================================================
AGE_PALETTE = {"young":"#1ABC9C","mid_age":"#F1C40F","old":"#C39BD3",
               "pre_geriatric":"#2980B9","geriatric":"#E84393"}
stars = lambda p: '***' if p<1e-3 else '**' if p<1e-2 else '*' if p<5e-2 else 'n.s.'

fig, axes = plt.subplots(1, len(sex_order), figsize=(4*len(sex_order), 3.5))
for ax, sx in zip(np.atleast_1d(axes), sex_order):
    sub = comp[comp.sex == sx]
    pvals = {}
    for ct in ORDER:
        g = [sub[(sub[GROUP_COL]==ct)&(sub.age==a)]['pct'].values for a in age_groups]
        g = [x for x in g if len(x) >= 2]
        pvals[ct] = f_oneway(*g)[1] if len(g) >= 2 else np.nan
    cts  = [c for c in ORDER if not np.isnan(pvals[c])]
    qmap = dict(zip(cts, multipletests([pvals[c] for c in cts], method='fdr_bh')[1])) if cts else {}
    for ci, ct in enumerate(ORDER):
        for ai, ag in enumerate(age_groups):
            vals = sub[(sub[GROUP_COL]==ct)&(sub.age==ag)]['pct'].values
            pos  = ci*(len(age_groups)+1) + ai
            bp = ax.boxplot(vals, positions=[pos], widths=0.7, patch_artist=True,
                            flierprops=dict(marker='o', mfc='red', ms=3, mec='none'))
            bp['boxes'][0].set_facecolor(AGE_PALETTE[ag]); bp['boxes'][0].set_alpha(0.85)
            bp['medians'][0].set_color('black')
        q = qmap.get(ct, np.nan); s = stars(q) if not np.isnan(q) else 'n.s.'
        center = ci*(len(age_groups)+1) + (len(age_groups)-1)/2
        ymax = sub[sub[GROUP_COL]==ct]['pct'].max()
        ax.text(center, (ymax*1.05+1) if np.isfinite(ymax) else 1, s, ha='center',
                color='red' if s!='n.s.' else 'darkblue', fontsize=11, fontweight='bold')
    ax.set_xticks([ci*(len(age_groups)+1)+(len(age_groups)-1)/2 for ci in range(len(ORDER))])
    ax.set_xticklabels(ORDER, rotation=45, ha='right'); ax.set_title(sx)
    ax.set_ylabel('% of all cells'); ax.set_facecolor('#f7f7f7')
fig.legend([plt.Rectangle((0,0),1,1, fc=AGE_PALETTE[a], alpha=0.85) for a in age_groups],
           age_groups, loc='center left', bbox_to_anchor=(1.0,0.5), frameon=False)
plt.tight_layout()
plt.savefig('myeloid_composition_box.pdf', bbox_inches='tight'); plt.show()
comp.to_csv('myeloid_composition_per_sample.csv', index=False)

# =========================================================================
# 4. supplementary table: pooled counts + % + per-sex ANOVA (BH across subclusters)
# =========================================================================
meta = parent.obs[[SAMPLE_KEY,'age','sex']].drop_duplicates().set_index(SAMPLE_KEY).reindex(total_per_sample.index)
cwm  = counts_wide.join(meta)
pct_mat = (100 * counts_wide.div(total_per_sample, axis=0)).join(meta)
pct_mat['age'] = pd.Categorical(pct_mat['age'], age_groups, ordered=True)

age_alpha = ['geriatric','mid_age','old','pre_geriatric','young']
rows = []
for ag in age_alpha:
    for sx in ['female','male']:
        members = sorted(cwm.index[(cwm['sex']==sx)&(cwm['age']==ag)].tolist())
        if not members: continue
        n   = counts_wide.loc[members, ORDER].sum(0)
        tot = int(total_per_sample.loc[members].sum())
        rows.append(pd.Series({**n.to_dict(), 'total cell count': tot},
                              name=f"{disp[ag]}_{sx} ({' + '.join(members)})"))
counts_tbl = pd.DataFrame(rows).astype(int)
counts_tbl.loc['total cell count'] = counts_tbl.sum(0)
pct_tbl = (100 * counts_tbl[ORDER].div(counts_tbl['total cell count'], axis=0)).round(2)

stat_rows = []
for sx in sex_order:
    ssub = pct_mat[pct_mat['sex']==sx]
    F, P = {}, {}
    for ct in ORDER:
        grp = [ssub.loc[ssub['age']==a, ct].values for a in age_groups]
        grp = [g for g in grp if len(g) >= 2]
        F[ct], P[ct] = f_oneway(*grp) if len(grp) >= 2 else (np.nan, np.nan)
    valid = [c for c in ORDER if not np.isnan(P[c])]
    padj  = dict(zip(valid, multipletests([P[c] for c in valid], method='fdr_bh')[1])) if valid else {}
    for ct in ORDER:
        stat_rows.append({'sex':sx, 'subcluster':ct,
                          'F_stat':F[ct], 'pval':P[ct], 'padj':padj.get(ct, np.nan)})
stats = pd.DataFrame(stat_rows)

counts_out = counts_tbl.reset_index().rename(columns={'index':'sample group'})
with pd.ExcelWriter('myeloid_subcluster_composition.xlsx', engine='openpyxl') as xl:
    counts_out.to_excel(xl, 'composition', startrow=1, startcol=0, index=False)
    c1 = counts_out.shape[1] + 1
    pct_tbl.reset_index(drop=True).to_excel(xl, 'composition', startrow=1, startcol=c1, index=False)
    c2 = c1 + pct_tbl.shape[1] + 1
    stats.to_excel(xl, 'composition', startrow=1, startcol=c2, index=False)
    ws = xl.sheets['composition']
    ws.cell(row=1, column=1,    value='counts')
    ws.cell(row=1, column=c1+1, value='percentage (% of all liver cells)')
    ws.cell(row=1, column=c2+1, value='statistics')

counts_tbl.to_csv('myeloid_counts_per_group.csv')
pct_tbl.to_csv('myeloid_percentage_of_all_cells_per_group.csv')
stats.to_csv('myeloid_stats_anova_age.csv', index=False)

print('parent cells :', parent.n_obs, '| subset cells:', adata.n_obs)
print('pooled total (should equal parent):', int(counts_tbl.loc["total cell count","total cell count"]))
