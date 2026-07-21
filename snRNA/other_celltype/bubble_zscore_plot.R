#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Multi-compartment marker bubble plots (RNA + ATAC) — liver aging multiome.
#
# One script for every compartment. Pairs with compartment_downstream.py and
# reads the files it writes, per compartment:
#   {PREFIX}_DE_Cluster_{sanitized}.csv   (per-subcluster Wilcoxon DE, RNA)
#   {PREFIX}_expr_stats.tsv               (pct_expressed + avg_expression, RNA)
#   {PREFIX}_sub_labels.csv               (rna_barcode -> subcluster, ATAC transfer)
#
# For each compartment it produces two z-scored bubble plots (size = % cells
# 20-100, color = Z-score across groups):
#   PART A (RNA)  : python DE + expr-stats exports
#   PART B (ATAC) : ArchR GeneScoreMatrix marker features on transferred labels
#
# Label vocabulary — the trap this handles:
#   `labels` in each COMPARTMENTS entry are the VERBATIM subcluster names as
#   stored in obs / the CSVs (e.g. "Kupffer cycling", "LSEC-like", "MV portal").
#   Joins + ArchR grouping use them verbatim; DE FILENAMES are sanitized
#   (space/hyphen -> underscore) to match python's safe_name(), via sanitize().
#   Set `labels` to whatever your DE files + expr_stats 'celltype' column use.
#
# Usage: edit COMPARTMENTS + the RUN switches, then
#        `Rscript compartment_marker_bubbles.R`
#        PART B needs the ArchR project + a cluster session; set DO_ATAC <- FALSE
#        to run only the RNA panels. RUN restricts which compartments execute.
#
# License: MIT. Cite: Sarker N. et al. "<paper title>." <journal>, <year>.
# ---------------------------------------------------------------------------

## ===================== GLOBAL SWITCHES =====================
DO_RNA  <- TRUE
DO_ATAC <- TRUE
ARCHR_PROJECT <- "/data/sarkern2/multiome_liver/Seurat/archR/ArchR_Projects/Step6_Xwnn_UMAP"
RUN <- NULL   # NULL = all; or e.g. c("Kupffer","endo") to run a subset
## ===========================================================

## ===================== PER-COMPARTMENT CONFIG =====================
COMPARTMENTS <- list(

  Kupffer = list(
    out_dir   = "/data/sarkern2/multiome_liver/Kupffer_DE_csvs",
    prefix    = "Kupffer",
    label_col = "celltype_sub",
    labels    = c("Kupffer", "Kupffer cycling", "LAM", "LSEC-like"),
    genes     = c("Clec4f","Vsig4","Timd4","Cd5l","Marco","C1qa",
                  "Mki67","Top2a","Cenpf",
                  "Gpnmb","Trem2","Cd9","Arhgap22","Sirpb1a","Cd74","Pparg",
                  "Stab2","Clec4g","Dnase1l3","Oit3","Kdr","Gata4",
                  "Ptprb","Egfl7","Pecam1","Cdh5"),
    fig_w = 9, fig_h = 3.5
  ),

  myeloid = list(
    out_dir   = "/data/sarkern2/multiome_liver/myeloid_DE_csvs",
    prefix    = "myeloid",
    label_col = "celltype_myeloid",
    labels    = c("Kupffer", "Kupffer cycling", "LAM", "MoMF",
                  "cDC1", "pDC", "Neutrophil"),
    genes     = c("Clec4f","Vsig4","Timd4","Cd5l","Marco","C1qa",
                  "Mki67","Top2a","Cenpf",
                  "Gpnmb","Trem2","Cd9","Arhgap22","Sirpb1a","Cd74","Pparg",
                  "Ccr2","Cx3cr1","Ly6c2","Chil3","Plac8","Fn1",
                  "Xcr1","Clec9a","Batf3","Cadm1",
                  "Siglech","Bst2","Il3ra","Tcf4",
                  "S100a8","S100a9","Retnlg","Csf3r"),
    fig_w = 11, fig_h = 4
  ),

  endo = list(
    out_dir   = "/data/sarkern2/multiome_liver/endo_DE_csvs",
    prefix    = "endo",
    label_col = "endo_subcluster",
    labels    = c("LSEC", "LSEC cycling", "MV portal", "MV central", "Kupffer-like"),
    genes     = c("Stab2","Clec4g","Gata4","Oit3","Dnase1l3","Mrc1",
                  "Mki67","Top2a","Cenpf",
                  "Vwf","Sox17","Efnb2",
                  "Rspo3","Wnt9b","Wnt2","Thbd",
                  "Clec4f","Vsig4","Timd4","Cd5l","C1qa","Marco","Csf1r"),
    fig_w = 11, fig_h = 4
  ),

  T_ILC = list(
    out_dir   = "/data/sarkern2/multiome_liver/T_DE_csvs",
    prefix    = "T",
    label_col = "celltype_T",
    labels    = c("CD4T", "Treg", "CD8T", "gdT", "iNKT", "NK", "ILC1", "neutrophil"),
    genes     = c("Cd3e","Cd4","Foxp3","Ikzf2",
                  "Cd8a","Cd8b1","Gzmk",
                  "Trdc","Zbtb16",
                  "Ncr1","Klrb1c","Gzmb","Eomes","Tbx21",
                  "S100a8","S100a9","Retnlg"),
    fig_w = 10, fig_h = 4
  )
)
## ==================================================================

suppressPackageStartupMessages({ library(ggplot2); library(dplyr) })

## ---- helpers -------------------------------------------------------------
sanitize <- function(x) gsub("-", "_", gsub(" ", "_", x))   # == python safe_name

clamp20_100 <- function(x) pmin(pmax(x, 20), 100)
size_scale_20_100 <- scale_size_continuous(
  range = c(2, 7), limits = c(20, 100), breaks = seq(20, 100, by = 20),
  name  = "fraction of cells in group (%)")

bubble_plot <- function(df, title, labels) {
  ggplot(df, aes(x = gene, y = celltype)) +
    geom_point(aes(size = PctPos_plot, color = Zscore)) +
    size_scale_20_100 +
    scale_color_gradient2(low = "blue", mid = "white", high = "red",
                          midpoint = 0, name = "Z-score") +
    scale_y_discrete(limits = rev(labels)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9, face = "italic"),
          axis.text.y = element_text(size = 9),
          axis.title  = element_blank()) +
    ggtitle(title)
}

zscore_per_gene <- function(df, value_col) {
  df %>% group_by(gene) %>%
    mutate(Zscore = { v <- .data[[value_col]]
                      as.numeric((v - mean(v)) / (sd(v) + 1e-9)) }) %>%
    ungroup()
}

usable_genes <- function(gene_order, present) {
  drop <- setdiff(gene_order, present)
  if (length(drop)) cat("  dropped (not in matrix):", paste(drop, collapse = ", "), "\n")
  gene_order[gene_order %in% present]
}

# RNA obs_name -> ArchR cell name (validated barcode reconstruction)
reconstruct_key <- function(rna_barcode) {
  parts   <- strsplit(rna_barcode, "-", fixed = TRUE)
  barcode <- vapply(parts, `[`, "", 1)
  sample  <- vapply(parts, function(z) paste(z[-1], collapse = "-"), "")
  sample  <- sub("^pre_geriatric", "pre_ger", sample)
  sample  <- sub("_([0-9]+)$", "", sample) |>
             paste0("_", sprintf("%02d", as.integer(sub(".*_", "", sample))))
  paste0(sample, "#", barcode, "-1")
}

## ---- per-compartment RNA panel ------------------------------------------
run_rna <- function(cfg) {
  cat("[RNA] ", cfg$prefix, "\n", sep = "")
  de <- bind_rows(lapply(cfg$labels, function(lab) {
    fn <- file.path(cfg$out_dir, paste0(cfg$prefix, "_DE_Cluster_", sanitize(lab), ".csv"))
    if (!file.exists(fn)) { cat("  MISSING:", basename(fn), "\n"); return(NULL) }
    d <- read.csv(fn); d$celltype <- lab; d          # verbatim label
  }))
  expr_stats <- read.delim(file.path(cfg$out_dir, paste0(cfg$prefix, "_expr_stats.tsv")))
  gene_use <- usable_genes(cfg$genes, unique(de$gene))

  plot_df <- de %>%
    select(gene, celltype, LogFC = avg_log2FC) %>%
    left_join(expr_stats, by = c("gene", "celltype")) %>%
    filter(gene %in% gene_use) %>%
    rename(PctPos = pct_expressed) %>%
    zscore_per_gene("LogFC") %>%
    mutate(celltype    = factor(celltype, levels = cfg$labels),
           gene        = factor(gene,     levels = gene_use),
           PctPos_plot = clamp20_100(PctPos))

  p <- bubble_plot(plot_df,
    sprintf("%s subclusters (RNA) — size: %% cells (20-100), color: Z-score", cfg$prefix),
    cfg$labels)
  fn <- file.path(cfg$out_dir, sprintf("PanelC_%s_RNA_bubble.pdf", cfg$prefix))
  ggsave(fn, p, width = cfg$fig_w, height = cfg$fig_h, limitsize = FALSE)
  cat("  saved", basename(fn), "\n")
}

## ---- per-compartment ATAC panel -----------------------------------------
run_atac <- function(cfg, proj) {
  cat("[ATAC] ", cfg$prefix, "\n", sep = "")
  lab <- read.csv(file.path(cfg$out_dir, paste0(cfg$prefix, "_sub_labels.csv")),
                  stringsAsFactors = FALSE)
  colnames(lab)[1] <- "rna_barcode"
  lab$archr_key <- reconstruct_key(lab$rna_barcode)

  idx <- match(getCellNames(proj), lab$archr_key)
  cat("  matched:", sum(!is.na(idx)), "/", nCells(proj), "\n")

  GRP <- paste0(cfg$prefix, "_sub")
  proj[[GRP]] <- lab[[cfg$label_col]][idx]           # verbatim labels
  print(table(proj[[GRP]], useNA = "ifany"))

  proj_c <- proj[!is.na(proj[[GRP]]), ]
  cat("  compartment cells:", nCells(proj_c), "\n")

  markersGS <- getMarkerFeatures(
    ArchRProj = proj_c, useMatrix = "GeneScoreMatrix", groupBy = GRP,
    bias = c("TSSEnrichment", "log10(nFrags)"), testMethod = "wilcoxon")
  saveRDS(markersGS, file.path(cfg$out_dir, sprintf("markersGS_%s_sub.rds", cfg$prefix)))

  L2FC <- assays(markersGS)[["Log2FC"]]; rn <- rowData(markersGS)$name
  rownames(L2FC) <- rn; colnames(L2FC) <- colnames(markersGS)

  gene_use  <- usable_genes(cfg$genes, rn)
  group_use <- cfg$labels[cfg$labels %in% colnames(L2FC)]

  plot_df <- expand.grid(gene = gene_use, celltype = group_use, stringsAsFactors = FALSE)
  plot_df$Log2FC <- mapply(function(g, c) L2FC[g, c], plot_df$gene, plot_df$celltype)

  gsm <- getMatrixFromProject(proj_c, useMatrix = "GeneScoreMatrix", binarize = FALSE)
  mat <- assay(gsm); rownames(mat) <- rowData(gsm)$name
  grp <- getCellColData(proj_c, GRP)[, 1]
  pct <- sapply(group_use, function(g) {
    cols <- which(grp == g)
    Matrix::rowMeans(mat[gene_use, cols, drop = FALSE] > 0) * 100
  })
  plot_df$PctPos <- mapply(function(g, c) pct[g, c], plot_df$gene, plot_df$celltype)
  write.csv(plot_df, file.path(cfg$out_dir, sprintf("plot_df_%s_ATAC.csv", cfg$prefix)),
            row.names = FALSE)

  plot_df <- plot_df %>%
    zscore_per_gene("Log2FC") %>%
    mutate(celltype    = factor(celltype, levels = cfg$labels),
           gene        = factor(gene,     levels = gene_use),
           PctPos_plot = clamp20_100(PctPos))

  p <- bubble_plot(plot_df,
    sprintf("%s subclusters (ATAC GeneScore) — size: %% cells (20-100), color: Z-score",
            cfg$prefix),
    cfg$labels)
  fn <- file.path(cfg$out_dir, sprintf("PanelC_%s_ATAC_bubble.pdf", cfg$prefix))
  ggsave(fn, p, width = cfg$fig_w, height = cfg$fig_h, limitsize = FALSE)
  cat("  saved", basename(fn), "\n")
}

## =========================================================================
## run
## =========================================================================
to_run <- if (is.null(RUN)) names(COMPARTMENTS) else RUN

# RNA panels (no ArchR needed)
if (DO_RNA) {
  for (name in to_run) {
    cat("\n===== ", name, " (RNA) =====\n", sep = "")
    run_rna(COMPARTMENTS[[name]])
  }
}

# ATAC panels (ArchR loaded once, reused across compartments)
if (DO_ATAC) {
  suppressPackageStartupMessages({ library(ArchR); library(Matrix) })
  addArchRThreads(16); addArchRGenome("mm10")
  proj <- loadArchRProject(ARCHR_PROJECT)
  for (name in to_run) {
    cat("\n===== ", name, " (ATAC) =====\n", sep = "")
    run_atac(COMPARTMENTS[[name]], proj)
  }
}

cat("\ndone.\n")
