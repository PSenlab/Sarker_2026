#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Compartment marker bubble plots (RNA + ATAC) — liver aging multiome.
#
# Pairs with compartment_downstream.py. Reads the files that script writes:
#   {PREFIX}_DE_Cluster_{sanitized}.csv   (per-subcluster Wilcoxon DE, RNA)
#   {PREFIX}_expr_stats.tsv               (pct_expressed + avg_expression, RNA)
#   {PREFIX}_sub_labels.csv               (rna_barcode -> subcluster, ATAC transfer)
#
# Produces two z-scored bubble plots (size = % cells 20-100, color = Z-score):
#   PART A (RNA)  : from the python DE + expr-stats exports
#   PART B (ATAC) : ArchR GeneScoreMatrix marker features on the transferred labels
#
# Label vocabulary — the single trap this script exists to handle:
#   LABELS below are the VERBATIM subcluster names as stored in obs / the CSVs
#   (e.g. "Kupffer cycling", "LSEC-like"). Joins and ArchR grouping use these.
#   DE filenames are SANITIZED (space/hyphen -> underscore) to match python's
#   safe_name(); sanitize() reproduces that mapping. Do not hand-edit one
#   without the other.
#
# Usage: edit the CONFIG block, then `Rscript compartment_marker_bubbles.R`
#        or source it. PART B needs the ArchR project + a cluster session; set
#        DO_ATAC <- FALSE to run the RNA panel alone.
#
# License: MIT. Cite: Sarker N. et al. "<paper title>." <journal>, <year>.
# ---------------------------------------------------------------------------

## ===================== CONFIG (edit per compartment) =====================
OUT_DIR <- "/data/sarkern2/multiome_liver/Kupffer_DE_csvs"
PREFIX  <- "Kupffer"          # must match --prefix used in compartment_downstream.py

# VERBATIM subcluster labels, in the display order you want on the y-axis.
LABELS <- c("Kupffer", "Kupffer cycling", "LAM", "LSEC-like")

# Marker panel (x-axis order). Genes absent from a matrix are dropped w/ a note.
GENE_ORDER <- c("Clec4f","Vsig4","Timd4","Cd5l","Marco","C1qa",     # Kupffer
                "Mki67","Top2a","Cenpf",                            # cycling
                "Gpnmb","Trem2","Cd9","Arhgap22","Sirpb1a","Cd74","Pparg", # LAM
                "Stab2","Clec4g","Dnase1l3","Oit3","Kdr","Gata4",   # LSEC
                "Ptprb","Egfl7","Pecam1","Cdh5")

ARCHR_PROJECT <- "/data/sarkern2/multiome_liver/Seurat/archR/ArchR_Projects/Step6_Xwnn_UMAP"
LABEL_COL     <- "celltype_sub"   # column name inside {PREFIX}_sub_labels.csv

DO_RNA  <- TRUE
DO_ATAC <- TRUE
FIG_W <- 9; FIG_H <- 3.5          # shared figure size for both panels
## =========================================================================

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr)
})

## ---- helpers -------------------------------------------------------------
sanitize <- function(x) gsub("-", "_", gsub(" ", "_", x))   # == python safe_name

clamp20_100 <- function(x) pmin(pmax(x, 20), 100)
size_scale_20_100 <- scale_size_continuous(
  range = c(2, 7), limits = c(20, 100), breaks = seq(20, 100, by = 20),
  name  = "fraction of cells in group (%)")

bubble_plot <- function(df, title) {
  # df needs columns: gene, celltype (verbatim), Zscore, PctPos_plot
  ggplot(df, aes(x = gene, y = celltype)) +
    geom_point(aes(size = PctPos_plot, color = Zscore)) +
    size_scale_20_100 +
    scale_color_gradient2(low = "blue", mid = "white", high = "red",
                          midpoint = 0, name = "Z-score") +
    scale_y_discrete(limits = rev(LABELS)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9, face = "italic"),
          axis.text.y = element_text(size = 9),
          axis.title  = element_blank()) +
    ggtitle(title)
}

zscore_per_gene <- function(df, value_col) {
  df %>% group_by(gene) %>%
    mutate(Zscore = {
      v <- .data[[value_col]]
      as.numeric((v - mean(v)) / (sd(v) + 1e-9))   # z across the groups present
    }) %>% ungroup()
}

# genes actually plottable (order preserved), with a dropped-gene note
usable_genes <- function(present) {
  drop <- setdiff(GENE_ORDER, present)
  if (length(drop)) cat("  dropped (not in matrix):", paste(drop, collapse = ", "), "\n")
  GENE_ORDER[GENE_ORDER %in% present]
}

## =========================================================================
## PART A — RNA bubble (from python DE + expr-stats exports)
## =========================================================================
if (DO_RNA) {
  cat("[RNA] loading python exports\n")

  # DE: one file per label, filename sanitized; tag celltype with the VERBATIM
  # label so it joins to expr_stats (which stores verbatim celltype).
  de <- bind_rows(lapply(LABELS, function(lab) {
    fn <- file.path(OUT_DIR, paste0(PREFIX, "_DE_Cluster_", sanitize(lab), ".csv"))
    if (!file.exists(fn)) { cat("  MISSING:", basename(fn), "\n"); return(NULL) }
    d <- read.csv(fn); d$celltype <- lab; d       # verbatim label
  }))
  expr_stats <- read.delim(file.path(OUT_DIR, paste0(PREFIX, "_expr_stats.tsv")))

  gene_use <- usable_genes(unique(de$gene))

  plot_df <- de %>%
    select(gene, celltype, LogFC = avg_log2FC) %>%
    left_join(expr_stats, by = c("gene", "celltype")) %>%   # verbatim <-> verbatim
    filter(gene %in% gene_use) %>%
    rename(PctPos = pct_expressed)

  plot_df <- zscore_per_gene(plot_df, "LogFC") %>%
    mutate(celltype    = factor(celltype, levels = LABELS),
           gene        = factor(gene,     levels = gene_use),
           PctPos_plot = clamp20_100(PctPos))

  p_rna <- bubble_plot(plot_df,
    sprintf("%s subclusters (RNA) — size: %% cells (20-100), color: Z-score", PREFIX))
  fn <- file.path(OUT_DIR, sprintf("PanelC_%s_RNA_bubble.pdf", PREFIX))
  ggsave(fn, p_rna, width = FIG_W, height = FIG_H, limitsize = FALSE)
  cat("  saved", basename(fn), "\n")
}

## =========================================================================
## PART B — ATAC bubble (ArchR GeneScoreMatrix on transferred labels)
## =========================================================================
if (DO_ATAC) {
  suppressPackageStartupMessages({
    library(ArchR); library(Matrix)
  })
  addArchRThreads(16); addArchRGenome("mm10")
  cat("[ATAC] loading ArchR project\n")

  proj <- loadArchRProject(ARCHR_PROJECT)
  lab  <- read.csv(file.path(OUT_DIR, paste0(PREFIX, "_sub_labels.csv")),
                   stringsAsFactors = FALSE)
  colnames(lab)[1] <- "rna_barcode"

  # --- barcode reconstruction: RNA obs_name -> ArchR cell name (validated) ---
  parts   <- strsplit(lab$rna_barcode, "-", fixed = TRUE)
  barcode <- vapply(parts, `[`, "", 1)
  sample  <- vapply(parts, function(z) paste(z[-1], collapse = "-"), "")
  sample  <- sub("^pre_geriatric", "pre_ger", sample)
  sample  <- sub("_([0-9]+)$", "", sample) |>
             paste0("_", sprintf("%02d", as.integer(sub(".*_", "", sample))))
  lab$archr_key <- paste0(sample, "#", barcode, "-1")

  cat("  ArchR names :", paste(head(getCellNames(proj), 2), collapse = " | "), "\n")
  cat("  archr_key   :", paste(head(lab$archr_key, 2),     collapse = " | "), "\n")

  idx <- match(getCellNames(proj), lab$archr_key)
  cat("  matched:", sum(!is.na(idx)), "/", nCells(proj), "\n")

  # transfer VERBATIM labels; group + factor on the same vocabulary as RNA
  GRP <- paste0(PREFIX, "_sub")
  proj[[GRP]] <- lab[[LABEL_COL]][idx]
  print(table(proj[[GRP]], useNA = "ifany"))

  proj_c <- proj[!is.na(proj[[GRP]]), ]
  cat("  compartment cells:", nCells(proj_c), "\n")

  # 1. bias-corrected marker features
  markersGS <- getMarkerFeatures(
    ArchRProj = proj_c, useMatrix = "GeneScoreMatrix", groupBy = GRP,
    bias = c("TSSEnrichment", "log10(nFrags)"), testMethod = "wilcoxon")
  saveRDS(markersGS, file.path(OUT_DIR, sprintf("markersGS_%s_sub.rds", PREFIX)))

  # 2. Log2FC + Mean -> named matrices
  L2FC <- assays(markersGS)[["Log2FC"]]; MEAN <- assays(markersGS)[["Mean"]]
  rn <- rowData(markersGS)$name
  rownames(L2FC) <- rownames(MEAN) <- rn
  colnames(L2FC) <- colnames(MEAN) <- colnames(markersGS)

  gene_use  <- usable_genes(rn)
  group_use <- LABELS[LABELS %in% colnames(L2FC)]   # verbatim, order preserved

  plot_df <- expand.grid(gene = gene_use, celltype = group_use,
                         stringsAsFactors = FALSE)
  plot_df$Log2FC <- mapply(function(g, c) L2FC[g, c], plot_df$gene, plot_df$celltype)

  # 3. PctPos from raw gene-score matrix
  gsm <- getMatrixFromProject(proj_c, useMatrix = "GeneScoreMatrix", binarize = FALSE)
  mat <- assay(gsm); rownames(mat) <- rowData(gsm)$name
  grp <- getCellColData(proj_c, GRP)[, 1]
  pct <- sapply(group_use, function(g) {
    cols <- which(grp == g)
    Matrix::rowMeans(mat[gene_use, cols, drop = FALSE] > 0) * 100
  })
  plot_df$PctPos <- mapply(function(g, c) pct[g, c], plot_df$gene, plot_df$celltype)
  write.csv(plot_df, file.path(OUT_DIR, sprintf("plot_df_%s_ATAC.csv", PREFIX)),
            row.names = FALSE)

  # 4. z-score + plot (same helper, same vocabulary as RNA panel)
  plot_df <- zscore_per_gene(plot_df, "Log2FC") %>%
    mutate(celltype    = factor(celltype, levels = LABELS),
           gene        = factor(gene,     levels = gene_use),
           PctPos_plot = clamp20_100(PctPos))

  p_atac <- bubble_plot(plot_df,
    sprintf("%s subclusters (ATAC GeneScore) — size: %% cells (20-100), color: Z-score",
            PREFIX))
  fn <- file.path(OUT_DIR, sprintf("PanelC_%s_ATAC_bubble.pdf", PREFIX))
  ggsave(fn, p_atac, width = FIG_W, height = FIG_H, limitsize = FALSE)
  cat("  saved", basename(fn), "\n")
}

cat("done.\n")
