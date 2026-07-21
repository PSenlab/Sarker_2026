library(schard)
library(scCustomize)
library(ggplot2)

# --- load; unify object name so downstream lines match ---
data <- schard::h5ad2seurat("myeloid.h5ad")   # <-- myeloid subset .h5ad
DefaultAssay(data) <- "RNA"
reduction_used <- "Xumap_"

# --- myeloid marker panel (per-lineage, one module score per gene) ---
gene_list <- c(
  # Kupffer cells
  "Clec4f", "Vsig4", "Timd4", "Cd5l", "Marco", "C1qa",
  # Cycling Kupffer cells
  "Mki67", "Top2a", "Cenpf",
  # Lipid-associated macrophages (LAM)
  "Gpnmb", "Trem2", "Cd9", "Arhgap22", "Sirpb1a", "Cd74", "Pparg",
  # Monocyte-derived macrophages (MoMF)
  "Ccr2", "Cx3cr1", "Ly6c2", "Chil3", "Plac8", "Fn1",
  # Conventional dendritic cells type 1 (cDC1)
  "Xcr1", "Clec9a", "Batf3", "Cadm1",
  # Plasmacytoid dendritic cells (pDC)
  "Siglech", "Bst2", "Il3ra", "Tcf4",
  # Neutrophils
  "S100a8", "S100a9", "Retnlg", "Csf3r"
)

# --- keep only genes actually present, so AddModuleScore never errors ---
missing   <- setdiff(gene_list, rownames(data))
if (length(missing)) message("skipped (not in RNA assay): ",
                             paste(missing, collapse = ", "))
gene_list <- gene_list[gene_list %in% rownames(data)]

# --- per-gene module score (seed BEFORE the loop; AddModuleScore is stochastic) ---
set.seed(1)
for (gene in gene_list) {
  data <- AddModuleScore(
    data,
    features = list(gene),
    name     = paste0(gene, "_ModuleScore")
  )
}

# --- one FeaturePlot per gene, clean gene-name title ---
score_cols <- paste0(gene_list, "_ModuleScore1")   # AddModuleScore appends "1"
for (i in seq_along(gene_list)) {
  gene <- gene_list[i]
  feat <- score_cols[i]
  p <- FeaturePlot_scCustom(
    data,
    reduction = reduction_used,
    features  = feat,
    order     = TRUE
  ) + ggtitle(gene)
  ggsave(
    filename = paste0("ModuleScore_", gene, ".pdf"),
    plot     = p,
    width    = 5,
    height   = 4
  )
}
