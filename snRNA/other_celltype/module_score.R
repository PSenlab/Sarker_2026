library(schard)
library(scCustomize)
library(ggplot2)

# =========================================================================
# per-compartment config: h5ad, marker panel, output dir / prefix
# =========================================================================
COMPARTMENTS <- list(

  myeloid = list(
    h5ad   = "myeloid.h5ad",
    outdir = "/data/sarkern2/multiome_liver/myeloid_DE_csvs/module_scores",
    prefix = "myeloid",
    genes  = c(
      "Clec4f","Vsig4","Timd4","Cd5l","Marco","C1qa",             # Kupffer
      "Mki67","Top2a","Cenpf",                                    # cycling
      "Gpnmb","Trem2","Cd9","Arhgap22","Sirpb1a","Cd74","Pparg",  # LAM
      "Ccr2","Cx3cr1","Ly6c2","Chil3","Plac8","Fn1",              # MoMF
      "Xcr1","Clec9a","Batf3","Cadm1",                            # cDC1
      "Siglech","Bst2","Il3ra","Tcf4",                            # pDC
      "S100a8","S100a9","Retnlg","Csf3r"                          # Neutrophil
    )
  ),

  Kupffer = list(
    h5ad   = "Kupffer.h5ad",
    outdir = "/data/sarkern2/multiome_liver/Kupffer_DE_csvs/module_scores",
    prefix = "Kupffer",
    genes  = c(
      "Clec4f","Vsig4","Timd4","Cd5l","Marco","C1qa",             # Kupffer
      "Mki67","Top2a","Cenpf",                                    # cycling
      "Gpnmb","Trem2","Cd9","Arhgap22","Sirpb1a","Cd74","Pparg",  # LAM
      "Stab2","Clec4g","Dnase1l3","Oit3","Kdr","Gata4"            # LSEC-like
    )
  ),

  endothelial_Kupffer02 = list(
    h5ad   = "endothelial_Kupffer02.h5ad",
    outdir = "/data/sarkern2/multiome_liver/endo_DE_csvs/module_scores",
    prefix = "endothelial_Kupffer02",
    genes  = c(
      "Stab2","Clec4g","Gata4","Oit3","Dnase1l3","Mrc1",          # pan-LSEC
      "Mki67","Top2a","Cenpf",                                    # cycling
      "Vwf","Sox17","Efnb2",                                      # MV portal
      "Rspo3","Wnt9b","Wnt2","Thbd",                              # MV central
      "Clec4f","Vsig4","Timd4","Cd5l","C1qa","Marco","Csf1r"      # Kupffer-like / KC2
    )
  ),

  T_ILC = list(
    h5ad   = "T_ILC.h5ad",
    outdir = "/data/sarkern2/multiome_liver/T_DE_csvs/module_scores",
    prefix = "T_ILC",
    genes  = c(
      "Cd3e","Cd4","Foxp3","Ikzf2",                               # CD4T / Treg
      "Cd8a","Cd8b1","Gzmk",                                      # CD8T
      "Trdc","Zbtb16",                                            # gdT / iNKT
      "Ncr1","Klrb1c","Gzmb","Eomes","Tbx21",                     # NK / ILC1
      "S100a8","S100a9","Retnlg"                                  # neutrophil
    )
  )
)

REDUCTION <- "Xumap_"   # schard-converted UMAP key (check Reductions(data) if blank)
RUN <- NULL             # NULL = all; or e.g. c("Kupffer","endothelial_Kupffer02")

# =========================================================================
# per-compartment runner
# =========================================================================
run_module_scores <- function(cfg) {
  cat("\n===== ", cfg$prefix, " =====\n", sep = "")
  dir.create(cfg$outdir, showWarnings = FALSE, recursive = TRUE)

  data <- schard::h5ad2seurat(cfg$h5ad)
  DefaultAssay(data) <- "RNA"

  # keep only genes present so AddModuleScore never errors
  genes   <- cfg$genes
  missing <- setdiff(genes, rownames(data))
  if (length(missing)) message("  skipped (not in RNA assay): ",
                               paste(missing, collapse = ", "))
  genes <- genes[genes %in% rownames(data)]

  # per-gene module score (seed before the loop; AddModuleScore is stochastic)
  set.seed(1)
  for (gene in genes) {
    data <- AddModuleScore(data, features = list(gene),
                           name = paste0(gene, "_ModuleScore"))
  }

  # one FeaturePlot per gene, clean gene-name title, prefixed filename
  score_cols <- paste0(genes, "_ModuleScore1")   # AddModuleScore appends "1"
  for (i in seq_along(genes)) {
    gene <- genes[i]
    p <- FeaturePlot_scCustom(
      data, reduction = REDUCTION, features = score_cols[i], order = TRUE
    ) + ggtitle(gene)
    ggsave(
      filename = file.path(cfg$outdir,
                           paste0(cfg$prefix, "_ModuleScore_", gene, ".pdf")),
      plot = p, width = 5, height = 4
    )
  }
  cat("  saved ", length(genes), " module-score plots to ", cfg$outdir, "\n", sep = "")
}

# =========================================================================
# run
# =========================================================================
to_run <- if (is.null(RUN)) names(COMPARTMENTS) else RUN
for (name in to_run) {
  run_module_scores(COMPARTMENTS[[name]])
}
cat("\ndone.\n")
