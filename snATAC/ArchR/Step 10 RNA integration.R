#!/usr/bin/env Rscript
# ==============================================================================
# ArchR RNA Integration: Add Gene Expression Matrix (Step 10)
# ==============================================================================
#
# Description:
#   This script integrates RNA expression data from a Scanpy h5ad object into
#   the ArchR project, enabling multi-omic analysis (ATAC + RNA).
#
# Prerequisites:
#   - Run archr_arrow_creation.R (Step 1)
#   - Run archr_preprocessing.R (Steps 2-7)
#   - Run archr_downstream_analysis.R (Step 8)
#   - Run archr_bigwig_generation.R (Step 9) 
#
# Input:
#   - ArchR project from Step 8
#   - RNA h5ad file (from multiome_integration.py)
#
# Output:
#   - ArchR project with GeneExpressionMatrix added
#
# Pipeline Overview:
#   10.1: Load ArchR Project from Step 8
#   10.2: Load RNA h5ad and Convert to SCE
#   10.3: Align RNA Barcodes to ArchR Format
#   10.4: Align Gene Universe
#   10.5: Add rowRanges from Gene Annotation
#   10.6: Build Clean SummarizedExperiment
#   10.7: Add RNA to ArchR Project
#   10.8: Save Project
#
# Requirements:
#   - R >= 4.0
#   - ArchR >= 1.0.2
#   - schard (for h5ad conversion)
#   - SummarizedExperiment
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
    library(SummarizedExperiment)
    library(schard)
    library(BSgenome.Mmusculus.UCSC.mm10)
})

# Set Global Parameters
addArchRThreads(threads = 60)
addArchRGenome("mm10")

# Define Paths (modify according to your directory structure)
STEP8_PROJ_PATH <- "path/to/ArchR_Projects/Step8_Imputed"
RNA_H5AD_PATH <- "path/to/rna_wnn.h5ad"
OUTPUT_BASE <- "path/to/ArchR_Projects"

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

# ==============================================================================
# 10.1: LOAD ARCHR PROJECT
# ==============================================================================

banner("10.1: Load ArchR Project from Step 8")

if (!dir.exists(STEP8_PROJ_PATH)) {
    stop(sprintf("[FATAL] Project path does not exist: %s", STEP8_PROJ_PATH))
}

proj <- loadArchRProject(path = STEP8_PROJ_PATH)
message(sprintf("[10.1] Loaded project from Step 8: %s", STEP8_PROJ_PATH))
message(sprintf("[10.1] Total cells: %d", nCells(proj)))

# ==============================================================================
# 10.2: LOAD RNA H5AD AND CONVERT TO SCE
# ==============================================================================

banner("10.2: Load RNA h5ad and Convert to SCE")

if (!file.exists(RNA_H5AD_PATH)) {
    stop(sprintf("[FATAL] RNA h5ad file not found: %s", RNA_H5AD_PATH))
}

# Convert h5ad to SingleCellExperiment using schard
sce <- schard::h5ad2sce(RNA_H5AD_PATH)

# Fix rownames if missing
if (is.null(rownames(sce)) || length(rownames(sce)) == 0) {
    rd <- rowData(sce)
    if ("symbol" %in% colnames(rd)) {
        rownames(sce) <- rd$symbol
    } else if ("gene_name" %in% colnames(rd)) {
        rownames(sce) <- rd$gene_name
    } else if ("index" %in% colnames(rd)) {
        rownames(sce) <- rd$index
    } else {
        stop("[FATAL] Gene names not found in rowData.")
    }
}

message(sprintf("[10.2] RNA loaded: %d genes x %d cells", nrow(sce), ncol(sce)))

# ==============================================================================
# 10.3: ALIGN RNA BARCODES TO ARCHR FORMAT
# ==============================================================================

banner("10.3: Align RNA Barcodes to ArchR Format")

rna_cells <- colnames(sce)

# Extract barcode and sample from RNA cell names
rna_barcode <- sub("^([A-Z0-9]+)-.*$", "\\1", rna_cells)
rna_sample <- sub("^[A-Z0-9]+-(.*)$", "\\1", rna_cells)

# Convert sample names to ArchR format (add leading zero)
rna_sample_archr <- sub("_([0-9])$", "_0\\1", rna_sample)

# Create ArchR-compatible cell names
rna_archr_names <- paste0(rna_sample_archr, "#", rna_barcode, "-1")

# Get ArchR cell names
archr_cells <- rownames(proj@cellColData)

# Filter to matching cells
keep <- rna_archr_names %in% archr_cells
sce <- sce[, keep, drop = FALSE]
colnames(sce) <- rna_archr_names[keep]

# Reorder strictly to match ArchR
common_cells <- intersect(archr_cells, colnames(sce))
sce <- sce[, common_cells]

stopifnot(identical(colnames(sce), common_cells))
message(sprintf("[10.3] RNA-ATAC cell alignment: %d cells", ncol(sce)))

# ==============================================================================
# 10.4: ALIGN GENE UNIVERSE
# ==============================================================================

banner("10.4: Align Gene Universe")

genes_archr <- getFeatures(proj, "GeneScoreMatrix")
common_genes <- intersect(genes_archr, rownames(sce))
sce <- sce[common_genes, ]

stopifnot(!anyDuplicated(rownames(sce)))
message(sprintf("[10.4] Gene alignment: %d common genes", nrow(sce)))

# ==============================================================================
# 10.5: ADD ROWRANGES FROM GENE ANNOTATION
# ==============================================================================
 
banner("10.5: Add rowRanges from Gene Annotation")

geneAnnotation <- getGeneAnnotation(proj)
genes_gr <- geneAnnotation$genes

# Match genes to annotation
idx <- match(rownames(sce), genes_gr$symbol)
valid <- !is.na(idx)

sce <- sce[valid, ]
matched_gr <- genes_gr[idx[valid]]
names(matched_gr) <- rownames(sce)

# Filter to valid chromosomes
chromSizes <- getChromSizes(proj)
valid_chr <- as.character(seqnames(matched_gr)) %in% as.character(seqnames(chromSizes))

sce <- sce[valid_chr, ]
matched_gr <- matched_gr[valid_chr]

rowRanges(sce) <- matched_gr
message(sprintf("[10.5] rowRanges added: %d genes with valid coordinates", nrow(sce)))

# ==============================================================================
# 10.6: BUILD CLEAN SUMMARIZEDEXPERIMENT
# ==============================================================================

banner("10.6: Build Clean SummarizedExperiment")

# Determine assay name
assay_name <- if ("counts" %in% assayNames(sce)) {
    "counts"
} else {
    assayNames(sce)[1]
}

counts_mat <- assay(sce, assay_name)
rownames(counts_mat) <- rownames(sce)
colnames(counts_mat) <- colnames(sce)

seRNA <- SummarizedExperiment(
    assays = list(counts = counts_mat),
    rowRanges = rowRanges(sce)
)

# Validate
stopifnot(
    identical(colnames(seRNA), common_cells),
    length(rowRanges(seRNA)) == nrow(seRNA)
)

message(sprintf("[10.6] SummarizedExperiment built: %d genes x %d cells", 
                nrow(seRNA), ncol(seRNA)))

# ==============================================================================
# 10.7: ADD RNA TO ARCHR PROJECT
# ==============================================================================

banner("10.7: Add RNA to ArchR Project")

proj <- addGeneExpressionMatrix(
    input       = proj,
    seRNA       = seRNA,
    strictMatch = TRUE,
    force       = TRUE
)

message(sprintf("[10.7] GeneExpressionMatrix added to project"))

# ==============================================================================
# 10.8: SAVE PROJECT
# ==============================================================================

banner("10.8: Save Project")

step10_dir <- file.path(OUTPUT_BASE, "Step10_RNA_Integrated")
ensure_dir(step10_dir)

proj <- saveArchRProject(
    ArchRProj = proj,
    outputDirectory = step10_dir,
    load = TRUE
)

message(sprintf("[10.8] Saved: %s", step10_dir))

# ==============================================================================
# PIPELINE SUMMARY
# ==============================================================================

banner("STEP 10 COMPLETE: RNA INTEGRATION")

message("Summary:")
message(sprintf("  Total cells: %d", nCells(proj)))
message(sprintf("  Genes in GeneExpressionMatrix: %d", nrow(seRNA)))
message(sprintf("  Available matrices: %s", paste(getAvailableMatrices(proj), collapse = ", ")))
message(sprintf("  Project location: %s", step10_dir))

message("\nNext steps:")
message("  - Peak-to-gene linkage analysis")
message("  - Regulatory network inference")
message("  - Multi-omic integration analysis")

# Save session info for reproducibility
sessionInfo()
