#!/usr/bin/env Rscript
# ==============================================================================
# ArchR Downstream Analysis: Gene Scores and Marker Discovery (Step 8)
# ==============================================================================
#
# Description:
#   This script performs downstream analysis on the preprocessed ArchR project
#   including imputation, gene score extraction, and marker gene discovery.
#
# Prerequisites:
#   - Run archr_arrow_creation.R (Step 1)
#   - Run archr_preprocessing.R (Steps 2-7)
#
# Input:
#   - Preprocessed ArchR project (from archr_preprocessing.R Step 7)
#
# Output:
#   - Gene score matrix (full and grouped by cell type)
#   - Marker genes per cell type
#
# Pipeline Overview:
#   8.1: Load ArchR Project from Step 7
#   8.2: Add Imputation Weights
#   8.3: Export Full Gene Score Matrix
#   8.4: Export Grouped Gene Scores by Cell Type
#   8.5: Marker Gene Discovery
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
    library(BSgenome.Mmusculus.UCSC.mm10)
})

# Set Global Parameters
addArchRThreads(threads = 60)
addArchRGenome("mm10")

# Define Paths (modify according to your directory structure)
# Load project from Step 7 of archr_preprocessing.R
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
    groupBy   = "celltype",
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
    groupBy    = "celltype",
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
# PIPELINE SUMMARY
# ==============================================================================

banner("STEP 8 COMPLETE: DOWNSTREAM ANALYSIS")

message("Summary:")
message(sprintf("  Cells analyzed: %d", nCells(proj)))
message(sprintf("  Cell types: %s", paste(unique(proj$celltype), collapse = ", ")))
message(sprintf("  Marker gene cutoff: FDR <= 0.05 & Log2FC >= 1"))
message(sprintf("  Output directory: %s", OUTPUT_DIR))

# Save session info for reproducibility
sessionInfo()
