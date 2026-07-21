library(ArchR); library(Matrix); library(dplyr); library(ggplot2)
addArchRThreads(16); addArchRGenome("mm10")

proj <- loadArchRProject(
    "/data/sarkern2/multiome_liver/Seurat/archR/ArchR_Projects/Step6_Xwnn_UMAP")

# =========================================================================
# per-compartment config: CSV, label column, marker panel, output name
# =========================================================================
COMPARTMENTS <- list(
  myeloid = list(
    csv       = "celltype_myeloid.csv",
    col       = "celltype_myeloid",
    outdir    = "/data/sarkern2/multiome_liver/scanpy_subcluster/subcluster/Kupffer/coverage_plots/Plots",
    pdf       = "BrowserTrack_myeloid_subclusters_labeled.pdf",
    panel     = list(
      Kupffer          = c("Clec4f", "Vsig4", "Timd4", "Cd5l"),
      Kupffer_cycling  = c("Mki67", "Top2a"),
      LAM              = c("Gpnmb", "Trem2", "Cd9"),
      MoMF             = c("Ccr2", "Ly6c2", "Plac8"),
      cDC1             = c("Xcr1", "Clec9a", "Batf3"),
      pDC              = c("Siglech", "Bst2", "Tcf4"),
      Neutrophil       = c("S100a8", "S100a9", "Retnlg"),
      Hepatocyte       = c("Alb"),          # cross-lineage negative control
      Endothelial      = c("Ptprb")         # cross-lineage negative control
    )
  ),
  Kupffer = list(
    csv       = "Kupffer_sub_labels.csv",
    col       = "celltype_sub",
    outdir    = "/data/sarkern2/multiome_liver/scanpy_subcluster/subcluster/Kupffer/coverage_plots/Plots",
    pdf       = "BrowserTrack_Kupffer_subclusters_labeled.pdf",
    panel     = list(
      Kupffer          = c("Clec4f", "Vsig4", "Timd4", "Cd5l"),
      `Kupffer cycling`= c("Mki67", "Top2a"),
      LAM              = c("Gpnmb", "Trem2", "Cd9"),
      `LSEC-like`      = c("Stab2", "Clec4g", "Kdr", "Gata4"),  # endothelial-like / KC2
      Hepatocyte       = c("Alb")           # cross-lineage negative control
    )
  )
)

# =========================================================================
# shared helpers
# =========================================================================
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

run_browser_tracks <- function(proj, cfg, name) {
  cat("\n===== ", name, " =====\n", sep = "")

  # --- barcode reconstruction + label transfer ---
  lab <- read.csv(cfg$csv, stringsAsFactors = FALSE)
  colnames(lab)[1] <- "rna_barcode"
  lab$archr_key <- reconstruct_key(lab$rna_barcode)

  idx <- match(getCellNames(proj), lab$archr_key)
  cat("matched:", sum(!is.na(idx)), "/", nCells(proj), "\n")

  proj[[cfg$col]] <- lab[[cfg$col]][idx]
  print(table(proj[[cfg$col]], useNA = "ifany"))

  # --- marker vector + gene-annotation guard ---
  markers <- unique(unlist(cfg$panel, use.names = FALSE))
  avail   <- getGenes(proj)$symbol
  missing <- setdiff(markers, avail)
  if (length(missing)) cat("skipped (not in mm10 geneAnnotation):",
                           paste(missing, collapse = ", "), "\n")
  markers <- markers[markers %in% avail]

  # --- tracks ---
  p <- plotBrowserTrack(
    ArchRProj  = proj,
    groupBy    = cfg$col,
    geneSymbol = markers,
    upstream   = 1e4, downstream = 1e4,
    normMethod = "ReadsInTSS"
  )

  # --- one labeled page per gene ---
  while (dev.cur() > 1) dev.off()          # clear any dangling device
  dir.create(cfg$outdir, showWarnings = FALSE, recursive = TRUE)
  pdf(file.path(cfg$outdir, cfg$pdf), width = 5, height = 5, onefile = TRUE)
  for (g in names(p)) {
    grid::grid.newpage()
    grid::grid.draw(p[[g]])
    grid::grid.text(g, x = 0.5, y = 0.985,
                    gp = grid::gpar(fontsize = 12, fontface = "italic",
                                    fontfamily = "sans"))   # base PDF: 'sans', not 'Arial'
  }
  dev.off()
  cat("saved", cfg$pdf, "\n")
}

# =========================================================================
# run both compartments
# =========================================================================
for (name in names(COMPARTMENTS)) {
  run_browser_tracks(proj, COMPARTMENTS[[name]], name)
}
cat("\ndone.\n")
