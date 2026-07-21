library(ArchR); library(Matrix); library(dplyr); library(ggplot2)
addArchRThreads(16); addArchRGenome("mm10")
proj <- loadArchRProject(
    "/data/sarkern2/multiome_liver/Seurat/archR/ArchR_Projects/Step6_Xwnn_UMAP")

# --- barcode reconstruction (identical validated logic) ---
lab <- read.csv("celltype_myeloid.csv", stringsAsFactors = FALSE)   # <-- myeloid label CSV
colnames(lab)[1] <- "rna_barcode"
parts   <- strsplit(lab$rna_barcode, "-", fixed = TRUE)
barcode <- vapply(parts, `[`, "", 1)
sample  <- vapply(parts, function(z) paste(z[-1], collapse = "-"), "")
sample  <- sub("^pre_geriatric", "pre_ger", sample)
sample  <- sub("_([0-9]+)$", "", sample) |>
           paste0("_", sprintf("%02d", as.integer(sub(".*_", "", sample))))
lab$archr_key <- paste0(sample, "#", barcode, "-1")
idx <- match(getCellNames(proj), lab$archr_key)
cat("matched:", sum(!is.na(idx)), "/", nCells(proj), "\n")   # expect the myeloid count

# --- myeloid subcluster marker panel (lineage-defining, one+ per label) ---
marker_panel <- list(
  Kupffer         = c("Clec4f", "Vsig4", "Timd4", "Cd5l"),
  Kupffer_cycling = c("Mki67", "Top2a"),
  LAM             = c("Gpnmb", "Trem2", "Cd9"),
  MoMF            = c("Ccr2", "Ly6c2", "Plac8"),
  cDC1            = c("Xcr1", "Clec9a", "Batf3"),
  pDC             = c("Siglech", "Bst2", "Tcf4"),
  Neutrophil      = c("S100a8", "S100a9", "Retnlg"),
  # cross-lineage purity checks (should be closed in myeloid groups)
  Hepatocyte      = c("Alb"),
  Endothelial     = c("Ptprb")
)
markers <- unique(unlist(marker_panel, use.names = FALSE))

# --- assign labels (NA = non-myeloid cells kept visible in the table) ---
COL <- "celltype_myeloid"
proj[[COL]] <- lab[[COL]][idx]
table(proj[[COL]], useNA = "ifany")

# guard: markers must exist in the gene annotation before plotting
avail   <- getGenes(proj)$symbol
missing <- setdiff(markers, avail)
if (length(missing)) cat("skipped (not in mm10 geneAnnotation):",
                         paste(missing, collapse = ", "), "\n")
markers <- markers[markers %in% avail]

p <- plotBrowserTrack(
    ArchRProj  = proj,
    groupBy    = COL,
    geneSymbol = markers,
    upstream   = 1e4, downstream = 1e4,
    normMethod = "ReadsInTSS"
)

# --- write one labeled page per gene ---
while (dev.cur() > 1) dev.off()          # clear any dangling device
outdir <- "/data/sarkern2/multiome_liver/scanpy_subcluster/subcluster/Kupffer/coverage_plots/Plots"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
pdf(file.path(outdir, "BrowserTrack_myeloid_subclusters_labeled.pdf"),
    width = 5, height = 5, onefile = TRUE)
for (g in names(p)) {
    grid::grid.newpage()
    grid::grid.draw(p[[g]])
    grid::grid.text(g, x = 0.5, y = 0.985,
                    gp = grid::gpar(fontsize = 12, fontface = "italic",
                                    fontfamily = "sans"))   # base PDF: 'sans', not 'Arial'
}
dev.off()
