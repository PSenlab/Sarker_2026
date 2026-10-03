#!/usr/bin/env Rscript
# ==============================================================================
# ArchR downstream analysis: Gene scores and marker discovery (Step 8)
# ==============================================================================
#
# Description:
#   This script performs downstream analysis on the preprocessed ArchR project
#   including imputation, gene score extraction, and marker gene discovery.
#
#
# Input:
#   - Preprocessed ArchR project 
#
# Output:
#   - Gene score matrix (full and grouped by cell type)
#   - Marker genes per cell type
#   - ATAC gene score marker bubble plot (z-scored Log2FC, % cells > 0)
#
#
# Requirements:
#   - R >= 4.0
#   - ArchR >= 1.0.2
#   - BSgenome.Mmusculus.UCSC.mm10
#
# ==============================================================================

# ==============================================================================
# SETUP AND CONFIGURATION
# ==============================================================================

# Load Required Libraries
suppressPackageStartupMessages({
    library(ArchR)
    library(GenomicRanges)
    library(stringr)
    library(dplyr)
    library(ggplot2)
    library(BSgenome.Mmusculus.UCSC.mm10)
})

# Set Global Parameters
addArchRThreads(threads = 60)
addArchRGenome("mm10")

# Define Paths (modify according to your directory structure)
# Load project from Step 7 of 03_archr_preprocessing.R
STEP7_PROJ_PATH <- "path/to/ArchR_Projects/Step7_Xwnn_UMAP"
OUTPUT_DIR <- "path/to/output/directory"

# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================

#' Print a formatted banner for pipeline steps
#' @param text Text to display in banner
banner <- function(text) {
    line <- paste(rep("=", 80), collapse = "")
    message("\n", line)
    message(text)
    message(line, "\n")
}

#' Ensure directory exists
#' @param dir_path Path to directory
#' @return Invisibly returns the directory path
ensure_dir <- function(dir_path) {
    if (!dir.exists(dir_path)) {
        dir.create(dir_path, recursive = TRUE)
    }
    invisible(dir_path)
}

# Ensure output directory exists
ensure_dir(OUTPUT_DIR)

# ==============================================================================
# 8.1: LOAD ARCHR PROJECT
# ==============================================================================

banner("8.1: Load ArchR Project from Step 7")

proj <- loadArchRProject(path = STEP7_PROJ_PATH)
message(sprintf("[8.1] Loaded project from Step 7: %s", STEP7_PROJ_PATH))
message(sprintf("[8.1] Total cells: %d", nCells(proj)))

# ==============================================================================
# 8.2: ADD IMPUTATION WEIGHTS
# ==============================================================================

banner("8.2: Add Imputation Weights")

step8_dir <- file.path(dirname(STEP7_PROJ_PATH), "Step8_Imputed")
ensure_dir(step8_dir)

if (file.exists(file.path(step8_dir, "ArchRProject.rds"))) {
    proj <- loadArchRProject(step8_dir)
    message(sprintf("[8.2] Loaded existing project with imputation weights"))
} else {
    proj <- addImputeWeights(proj)
    proj <- saveArchRProject(proj, outputDirectory = step8_dir, load = TRUE)
    message(sprintf("[8.2] Imputation weights added"))
    message(sprintf("[8.2] Saved: %s", step8_dir))
}

# ==============================================================================
# 8.3: EXPORT FULL GENE SCORE MATRIX
# ==============================================================================

banner("8.3: Export Full Gene Score Matrix")

genescore_mat <- getMatrixFromProject(
    ArchRProj = proj,
    useMatrix = "GeneScoreMatrix",
    verbose   = TRUE,
    binarize  = FALSE,
    threads   = getArchRThreads()
)

saveRDS(genescore_mat, file = file.path(OUTPUT_DIR, "proj.GeneScoreMatrix.rds"))
message(sprintf("[8.3] Saved: proj.GeneScoreMatrix.rds"))

# ==============================================================================
# 8.4: EXPORT GROUPED GENE SCORES BY CELL TYPE
# ==============================================================================

banner("8.4: Export Grouped Gene Scores by Cell Type")

gene_se <- getGroupSE(
    ArchRProj = proj,
    useMatrix = "GeneScoreMatrix",
    groupBy   = "cell_type",
    divideN   = TRUE,
    verbose   = TRUE
)

gene_scores <- assay(gene_se, "GeneScoreMatrix")
rownames(gene_scores) <- rowData(gene_se)$name

write.table(
    data.frame(gene = rownames(gene_scores), gene_scores),
    file      = file.path(OUTPUT_DIR, "GeneMat_by_Cell_Subtype.tsv"),
    sep       = "\t",
    quote     = FALSE,
    row.names = FALSE
)

message(sprintf("[8.4] Saved: GeneMat_by_Cell_Subtype.tsv"))

# ==============================================================================
# 8.5: MARKER GENE DISCOVERY
# ==============================================================================

banner("8.5: Marker Gene Discovery")

markers <- getMarkerFeatures(
    ArchRProj  = proj,
    useMatrix  = "GeneScoreMatrix",
    groupBy    = "cell_type",
    bias       = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
)

# Save markers to project directory
saveRDS(markers, file = file.path(step8_dir, "MarkerFeatures_celltype.rds"))
message(sprintf("[8.5] Saved markers to project: %s", file.path(step8_dir, "MarkerFeatures_celltype.rds")))

# Also export to output directory
saveRDS(markers, file = file.path(OUTPUT_DIR, "proj.Cell_Subtype.MarkerGenes.rds"))
message(sprintf("[8.5] Exported: proj.Cell_Subtype.MarkerGenes.rds"))

# Extract significant markers
marker_list <- getMarkers(markers, cutOff = "FDR <= 0.05 & Log2FC >= 1")

# Save each marker list per cell type
for (subtype in names(marker_list)) {
    output_file <- file.path(OUTPUT_DIR, paste0("MarkerGenes_", subtype, ".csv"))
    write.csv(
        marker_list[[subtype]],
        file      = output_file,
        row.names = FALSE,
        quote     = FALSE
    )
    message(sprintf("[8.5] Saved: MarkerGenes_%s.csv", subtype))
}

# ==============================================================================
# 8.6: ATAC GENE SCORE MARKER BUBBLE PLOT
# ==============================================================================

banner("8.6: ATAC Gene Score Marker Bubble Plot")

# Cell type order (top to bottom) and marker genes
celltype_order <- c("mesenchymal", "Kupffer 02", "cholangiocyte 02",
                    "cholangiocyte 01", "B cells", "non-resident myeloid",
                    "T/ILC cells", "Kupffer 01", "endothelial", "hepatocyte")
gene_order <- c(
    "Dpyd", "Adk", "Errfi1", "Ptprd", "Abcc2", "Kel", "Cldn13", "Gm44204",
    "Traj3", "Trav8-2", "Slc8a1", "Myo9a", "Hdac9", "Fyb", "Cd5l", "Skap1",
    "Ptprc", "Dock2", "Arhgap15", "Ms4a4b", "Lyz2", "Gpr141", "Mctp1", "Ccr2",
    "Itgam", "Ebf1", "Bank1", "Atf3", "Cd74", "Pkhd1", "Bicc1", "Spp1",
    "Glis3", "2610307P16Rik", "Dmbt1", "Chrm3", "Thsd4", "Gm609", "Ptprb",
    "Stab2", "Meis2", "Fbxl7", "Sema6a", "Rbms3", "Reln", "Ank3", "Dcn",
    "Prkg1"
)

# Color: bias-corrected Log2FC from getMarkerFeatures (8.5)
L2FC <- assays(markers)[["Log2FC"]]
rownames(L2FC) <- rowData(markers)$name
colnames(L2FC) <- colnames(markers)

missing_genes  <- setdiff(gene_order, rownames(L2FC))
gene_order     <- gene_order[gene_order %in% rownames(L2FC)]
celltype_order <- celltype_order[celltype_order %in% colnames(L2FC)]
if (length(missing_genes) > 0) {
    message(sprintf("[8.6] Not in GeneScoreMatrix (dropped): %s",
                    paste(missing_genes, collapse = ", ")))
}

plot_df <- expand.grid(gene = gene_order, celltype = celltype_order,
                       stringsAsFactors = FALSE)
plot_df$Log2FC <- mapply(function(g, c) L2FC[g, c], plot_df$gene, plot_df$celltype)

# Size: % of cells with raw gene score > 0 (GeneScoreMatrix from 8.3)
gs_mat <- assay(genescore_mat)
rownames(gs_mat) <- rowData(genescore_mat)$name
grp <- colData(genescore_mat)$cell_type

pct <- sapply(celltype_order, function(g) {
    cols <- which(grp == g)
    Matrix::rowMeans(gs_mat[gene_order, cols, drop = FALSE] > 0) * 100
})
plot_df$PctPos <- mapply(function(g, c) pct[g, c], plot_df$gene, plot_df$celltype)

# Z-score of Log2FC per gene across cell types
plot_df <- plot_df %>%
    group_by(gene) %>%
    mutate(Zscore = as.numeric(scale(Log2FC))) %>%
    ungroup() %>%
    mutate(celltype = factor(celltype, levels = celltype_order),
           gene     = factor(gene,     levels = gene_order))

write.csv(plot_df, file.path(OUTPUT_DIR, "ATAC_GeneScore_bubble_data.csv"),
          row.names = FALSE)

# Size scale clamped to 20-100 %
clamp20_100 <- function(x) pmin(pmax(x, 20), 100)
size_scale_20_100 <- scale_size_continuous(
    range = c(2, 7), limits = c(20, 100), breaks = seq(20, 100, by = 20),
    name = "fraction of cells\nin group (%)")
plot_df$PctPos_plot <- clamp20_100(plot_df$PctPos)

p_atac <- ggplot(plot_df, aes(x = gene, y = celltype)) +
    geom_point(aes(size = PctPos_plot, color = Zscore)) +
    size_scale_20_100 +
    scale_color_gradient2(low = "blue", mid = "white", high = "red",
                          midpoint = 0, name = "z-score") +
    scale_y_discrete(limits = rev(celltype_order)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9, face = "italic"),
          axis.text.y = element_text(size = 9),
          axis.title  = element_blank(),
          plot.title  = element_text(hjust = 0.5)) +
    ggtitle("ATAC")

ggsave(file.path(OUTPUT_DIR, "ATAC_GeneScore_marker_bubble.pdf"),
       p_atac, width = 15, height = 4, limitsize = FALSE)
message("[8.6] Saved: ATAC_GeneScore_marker_bubble.pdf")

# ==============================================================================
# PIPELINE SUMMARY
# ==============================================================================

banner("STEP 8 COMPLETE: DOWNSTREAM ANALYSIS")

message("Summary:")
message(sprintf("  Cells analyzed: %d", nCells(proj)))
message(sprintf("  Cell types: %s", paste(unique(proj$cell_type), collapse = ", ")))
message(sprintf("  Marker gene cutoff: FDR <= 0.05 & Log2FC >= 1"))
message(sprintf("  Output directory: %s", OUTPUT_DIR))

# Save session info for reproducibility
sessionInfo()
