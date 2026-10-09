#!/usr/bin/env Rscript
# ==============================================================================
# ArchR browser-track coverage plots - liver aging multiome
# ==============================================================================
# Generates ReadsInTSS-normalized coverage tracks over marker/zonation genes
# for two groupings of the same ArchR project:
#   (1) all cell types      (groupBy = "cell_type")
#   (2) hepatocyte zonation (groupBy = "subcluster", ordered hep_01 -> hep_07)
#
# RNA-derived cell labels (from scRNA/Seurat) are matched back onto the ArchR
# cell names by reconstructing the ArchR barcode key from each RNA barcode.
#
# Usage:
#   Rscript browser_tracks_liver_multiome.R
#
# Dependencies: ArchR, Matrix, dplyr, ggplot2, grid
# ==============================================================================

suppressPackageStartupMessages({
  library(ArchR)
  library(Matrix)
  library(dplyr)
  library(ggplot2)
})

addArchRThreads(16)
addArchRGenome("mm10")

# ------------------------------------------------------------------------------
# Config
# ------------------------------------------------------------------------------
ARCHR_PROJECT <- "path/to/ArchR_Projects/Step7_Xwnn_UMAP"
OUTDIR        <- "Plots"

LABELS_ALL <- "celltype_all.csv"   # provides column `cell_type`
LABELS_HEP <- "celltype_hep.csv"   # provides column `subcluster`

# base R pdf() has no "Arial"; "sans" maps to the Helvetica/Arial-equivalent.
# For embedded true Arial, use showtext/extrafont instead.
FONT_FAMILY <- "sans"

UPSTREAM   <- 1e4
DOWNSTREAM <- 1e4

# ------------------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------------------

#' Rebuild the ArchR cell key ("<sample>#<barcode>-1") from an RNA barcode.
#' Mirrors the reconstruction validated on the endothelial compartment.
reconstruct_archr_keys <- function(rna_barcode) {
  parts   <- strsplit(rna_barcode, "-", fixed = TRUE)
  barcode <- vapply(parts, `[`, "", 1)
  sample  <- vapply(parts, function(z) paste(z[-1], collapse = "-"), "")
  sample  <- sub("^pre_geriatric", "pre_ger", sample)
  sample  <- sub("_([0-9]+)$", "", sample) |>
             paste0("_", sprintf("%02d", as.integer(sub(".*_", "", sample))))
  paste0(sample, "#", barcode, "-1")
}

#' Read a label CSV and return the requested label column aligned to proj cells
#' (NA for cells absent from the CSV, e.g. other compartments).
label_for_project <- function(proj, csv, label_col) {
  lab <- read.csv(csv, stringsAsFactors = FALSE)
  colnames(lab)[1] <- "rna_barcode"
  lab$archr_key <- reconstruct_archr_keys(lab$rna_barcode)
  idx <- match(getCellNames(proj), lab$archr_key)
  message(sprintf("[%s] matched %d / %d cells",
                  label_col, sum(!is.na(idx)), nCells(proj)))
  lab[[label_col]][idx]
}

#' Plot browser tracks for a set of genes and write one page per gene to a PDF.
#' `use_groups` (optional) both subsets and orders the tracks: the first entry
#' becomes the top panel, filling downward. Passing only the groups of interest
#' also drops the <NA> track (cells outside this labelling).
browser_track_pdf <- function(proj, group_by, markers, out_pdf,
                              use_groups = NULL,
                              width = 5, height = 5,
                              upstream = UPSTREAM, downstream = DOWNSTREAM,
                              font_family = FONT_FAMILY) {
  p <- plotBrowserTrack(
    ArchRProj  = proj,
    groupBy    = group_by,
    useGroups  = use_groups,
    geneSymbol = markers,
    upstream   = upstream,
    downstream = downstream,
    normMethod = "ReadsInTSS"
  )

  while (dev.cur() > 1) dev.off()          # clear any dangling device
  pdf(out_pdf, width = width, height = height, onefile = TRUE)
  on.exit(dev.off(), add = TRUE)

  for (g in names(p)) {
    grid::grid.newpage()
    grid::grid.draw(p[[g]])
    grid::grid.text(
      g, x = 0.5, y = 0.985,
      gp = grid::gpar(fontsize = 12, fontface = "italic", fontfamily = font_family)
    )
  }

  message(sprintf("wrote %s (%d genes)", out_pdf, length(p)))
  invisible(p)
}

# ------------------------------------------------------------------------------
# Load project
# ------------------------------------------------------------------------------
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
proj <- loadArchRProject(ARCHR_PROJECT)

# ------------------------------------------------------------------------------
# (1) All cell types
# ------------------------------------------------------------------------------
marker_panel <- list(
  # Parenchymal
  hepatocyte             = "Abcc2",
  cholangiocyte          = "Spp1",
  # Mesenchymal
  mesenchymal            = "Dcn",
  # Endothelial
  endothelial            = "Ptprb",
  # Myeloid
  Kupffer                = "Cd5l",
  `non-resident myeloid` = "Mctp1",
  # Lymphoid
  `T/ILC cells`          = "Ms4a4b",
  `B cells`              = "Ebf1"
)
markers_all <- unique(unlist(marker_panel, use.names = FALSE))

# Track order (top -> bottom). These strings must match the labels in the
# `cell_type` column EXACTLY (case + spacing); verify against the table below.
celltype_order <- c(
  "hepatocyte",
  "cholangiocyte 01",
  "cholangiocyte 02",
  "Kupffer 01",
  "Kupffer 02",
  "endothelial",
  "mesenchymal",
  "T/ILC cells",
  "B cells",
  "non-resident myeloid"
)

proj$cell_type <- label_for_project(proj, LABELS_ALL, "cell_type")
print(table(proj$cell_type, useNA = "ifany"))

# Sanity check: any name here that isn't in the column is silently ignored;
# any cell type not listed here is dropped from the plot.
missing <- setdiff(celltype_order, unique(stats::na.omit(proj$cell_type)))
if (length(missing)) warning("celltype_order not found in data: ",
                             paste(missing, collapse = ", "))

browser_track_pdf(
  proj, group_by = "cell_type", markers = markers_all,
  use_groups = celltype_order,   # top -> bottom; drops <NA>
  out_pdf = file.path(OUTDIR, "BrowserTrack_allcelltype_10_labeled.pdf")
)

# ------------------------------------------------------------------------------
# (2) Hepatocyte zonation (hep_01 periportal -> hep_07 pericentral)
# ------------------------------------------------------------------------------
zonation_markers <- list(
  periportal  = c("Cyp2f2", "Cdh1", "Hal"),   # zone 1
  midlobular  = "Hamp2",                       # zone 2
  pericentral = c("Cyp2e1", "Glul", "Cyp7a1")  # zone 3
)
markers_hep <- unique(unlist(zonation_markers, use.names = FALSE))

hep_order <- paste0("hep_", sprintf("%02d", 1:7))   # top -> bottom

proj$subcluster <- label_for_project(proj, LABELS_HEP, "subcluster")
print(table(proj$subcluster, useNA = "ifany"))

browser_track_pdf(
  proj, group_by = "subcluster", markers = markers_hep,
  use_groups = hep_order,   # hep_01 on top, hep_07 at bottom; drops <NA>
  out_pdf = file.path(OUTDIR, "BrowserTrack_Heps_zone_10_labeled.pdf")
)
