#!/usr/bin/env Rscript
# ==============================================================================
# ArchR Preprocessing Pipeline for Single-Cell ATAC-seq (Multiome)
# ==============================================================================
#
# Description:
#   This script performs preprocessing of single-cell ATAC-seq data from a 
#   10x Multiome experiment using ArchR. The pipeline includes quality control,
#   doublet removal, dimensionality reduction, batch correction, and integration
#   with Seurat-derived cell type annotations and embeddings.
#
# Prerequisites:
#   - Run archr_arrow_creation.R (Step 1) first to generate Arrow files
#
# Study Design:
#   - Species: Mouse (mm10)
#   - Tissue: Liver
#   - Age Groups: Young, Middle-age, Old, Pre-geriatric, Geriatric
#   - Replicates: 8 biological replicates per age group (n=40 total)
#
# Pipeline Overview:
#   Step 1: Arrow File Creation (archr_arrow_creation.R)
#   Step 2: ArchR Project Creation from Arrow files
#   Step 3: Doublet Detection and Filtering
#   Step 4: Iterative Latent Semantic Indexing (LSI)
#   Step 5: Harmony Batch Correction
#   Step 6: Seurat Metadata Integration
#   Step 7: WNN UMAP Embedding Transfer
#
# Requirements:
#   - R >= 4.0
#   - ArchR >= 1.0.2
#   - Seurat >= 4.0
#   - schard (for h5ad conversion)
#   - harmony
#   - Pre-generated Arrow files (from archr_arrow_creation.R)
#   - Integrated h5ad object with cell type annotations (rna_wnn.h5ad)
#
# Citation:
#   - ArchR: Granja et al., Nature Genetics 2021
#   - Harmony: Korsunsky et al., Nature Methods 2019
#   - Seurat: Hao et al., Cell 2021
#
# ==============================================================================

# ==============================================================================
# SETUP AND CONFIGURATION
# ==============================================================================

# Load Required Libraries
suppressPackageStartupMessages({
    library(ArchR)
    library(Seurat)
    library(schard)
    library(dplyr)
    library(stringr)
})

# Set Global Parameters
set.seed(10918)
addArchRThreads(threads = 60)
addArchRGenome("mm10")

# Define Paths (modify according to your directory structure)
WORK_DIR <- "path/to/working/directory"
OUTPUT_BASE <- file.path(WORK_DIR, "ArchR_Projects")
H5AD_PATH <- "path/to/rna_wnn.h5ad"

setwd(WORK_DIR)

if (!dir.exists(OUTPUT_BASE)) {
    dir.create(OUTPUT_BASE, recursive = TRUE)
}

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

#' Get expected Arrow file path
#' @param samp Sample name
#' @return Full path to Arrow file
expected_arrow <- function(samp) {
    file.path(WORK_DIR, paste0(samp, ".arrow"))
}

#' Convert Seurat barcodes to ArchR format
#' @param barcode Seurat-style barcode
#' @return ArchR-compatible cell name
convert_barcode_to_archr <- function(barcode) {
    raw_barcode <- str_extract(barcode, "^[^-]+")
    sample_raw <- str_extract(barcode, "(?<=-).+$")
    
    sample_archr <- case_when(
        str_detect(sample_raw, "^young_")     ~ str_replace(sample_raw, "^young_",     "young_0"),
        str_detect(sample_raw, "^mid_age_")   ~ str_replace(sample_raw, "^mid_age_",   "mid_age_0"),
        str_detect(sample_raw, "^old_")       ~ str_replace(sample_raw, "^old_",       "old_0"),
        str_detect(sample_raw, "^pre_ger_")   ~ str_replace(sample_raw, "^pre_ger_",   "pre_ger_0"),
        str_detect(sample_raw, "^geriatric_") ~ str_replace(sample_raw, "^geriatric_", "geriatric_0"),
        TRUE ~ sample_raw
    )
    
    paste0(sample_archr, "#", raw_barcode, "-1")
}

# ==============================================================================
# SAMPLE METADATA DEFINITION
# ==============================================================================

# Folder names correspond to original 10x output directories
folders <- c(
    paste0("Y", 1:8),
    paste0("MA_", sprintf("%02d", 1:8)),
    paste0("O", 1:8),
    paste0("PG_", sprintf("%02d", 1:8)),
    paste0("G", 1:8)
)

# Sample names for ArchR Arrow files
samples <- c(
    paste0("young_",     sprintf("%02d", 1:8)),
    paste0("mid_age_",   sprintf("%02d", 1:8)),
    paste0("old_",       sprintf("%02d", 1:8)),
    paste0("pre_ger_",   sprintf("%02d", 1:8)),
    paste0("geriatric_", sprintf("%02d", 1:8))
)

# Age group assignments
dataset_group <- c(
    rep("Young",         8),
    rep("Middle_age",    8),
    rep("Old",           8),
    rep("Pre_Geriatric", 8),
    rep("Geriatric",     8)
)

# Validate metadata consistency
stopifnot(
    length(folders) == length(samples),
    length(samples) == length(dataset_group)
)

# Create metadata data frame
meta_df <- data.frame(
    SampleName   = samples,
    DatasetGroup = dataset_group,
    Folder       = folders,
    stringsAsFactors = FALSE
)

banner("SAMPLE METADATA SUMMARY")
print(table(meta_df$DatasetGroup))

# ==============================================================================
# 2: ARCHR PROJECT CREATION
# ==============================================================================

banner("2: ArchR Project Creation")

step2_dir <- file.path(OUTPUT_BASE, "Step2_Project")
ensure_dir(step2_dir)

# Collect existing Arrow files (from archr_arrow_creation.R Step 1)
arrow_files <- sapply(meta_df$SampleName, expected_arrow)
arrow_files <- arrow_files[file.exists(arrow_files)]

if (length(arrow_files) == 0L) {
    stop("No Arrow files found. Run archr_arrow_creation.R (Step 1) first.")
}

message(sprintf("[2] Found %d Arrow files from Step 1", length(arrow_files)))

# Create or load project
if (file.exists(file.path(step2_dir, "ArchRProject.rds"))) {
    proj <- loadArchRProject(step2_dir)
    message(sprintf("[2] Loaded existing project"))
} else {
    message(sprintf("[2] Creating new ArchR project..."))
    
    proj <- ArchRProject(
        ArrowFiles      = arrow_files,
        outputDirectory = step2_dir,
        copyArrows      = TRUE
    )
    
    # Add age group metadata
    sample_to_group <- setNames(meta_df$DatasetGroup, meta_df$SampleName)
    proj$DatasetGroup <- ArchR::mapLabels(
        x         = proj$Sample,
        newLabels = sample_to_group,
        oldLabels = names(sample_to_group)
    )
    
    proj <- saveArchRProject(proj, outputDirectory = step2_dir, load = TRUE)
    message(sprintf("[2] Saved: %s", step2_dir))
}

# Ensure DatasetGroup exists
if (!"DatasetGroup" %in% colnames(proj@cellColData)) {
    sample_to_group <- setNames(meta_df$DatasetGroup, meta_df$SampleName)
    proj$DatasetGroup <- ArchR::mapLabels(
        x         = proj$Sample,
        newLabels = sample_to_group,
        oldLabels = names(sample_to_group)
    )
    proj <- saveArchRProject(proj, step2_dir, load = TRUE)
}

message(sprintf("[2] Total cells: %d", nCells(proj)))

# ==============================================================================
# 3: DOUBLET DETECTION AND FILTERING
# ==============================================================================

banner("3: Doublet Detection and Filtering")

step3_dir <- file.path(OUTPUT_BASE, "Step3_DoubletsFiltered")
ensure_dir(step3_dir)

if (file.exists(file.path(step3_dir, "ArchRProject.rds"))) {
    proj <- loadArchRProject(step3_dir)
    message(sprintf("[3] Loaded existing project"))
} else {
    message(sprintf("[3] Running doublet detection..."))
    
    # Add doublet scores if not present
    if (is.null(proj$DoubletEnrichment)) {
        proj <- addDoubletScores(
            input     = proj,
            k         = 10,
            knnMethod = "UMAP",
            LSIMethod = 1
        )
    }
    
    # Filter doublets
    proj <- filterDoublets(proj)
    
    proj <- saveArchRProject(proj, outputDirectory = step3_dir, load = TRUE)
    message(sprintf("[3] Saved: %s", step3_dir))
}

message(sprintf("[3] Cells after doublet filtering: %d", nCells(proj)))

# ==============================================================================
# 4: ITERATIVE LSI DIMENSIONALITY REDUCTION
# ==============================================================================

banner("4: Iterative LSI Dimensionality Reduction")

step4_dir <- file.path(OUTPUT_BASE, "Step4_LSI")
ensure_dir(step4_dir)

if (file.exists(file.path(step4_dir, "ArchRProject.rds"))) {
    proj <- loadArchRProject(step4_dir)
    message(sprintf("[4] Loaded existing project"))
} else {
    message(sprintf("[4] Running Iterative LSI..."))
    
    proj <- addIterativeLSI(
        ArchRProj     = proj,
        useMatrix     = "TileMatrix",
        name          = "IterativeLSI",
        iterations    = 2,
        clusterParams = list(
            resolution  = c(0.1),
            sampleCells = 10000,
            n.start     = 10
        ),
        varFeatures   = 25000,
        dimsToUse     = 1:30
    )
    
    proj <- saveArchRProject(proj, outputDirectory = step4_dir, load = TRUE)
    message(sprintf("[4] Saved: %s", step4_dir))
}

# ==============================================================================
# 5: HARMONY BATCH CORRECTION
# ==============================================================================

banner("5: Harmony Batch Correction")

step5_dir <- file.path(OUTPUT_BASE, "Step5_Harmony")
ensure_dir(step5_dir)

if (file.exists(file.path(step5_dir, "ArchRProject.rds"))) {
    proj <- loadArchRProject(step5_dir)
    message(sprintf("[5] Loaded existing project"))
} else {
    message(sprintf("[5] Running Harmony batch correction..."))
    
    # Increase memory limit for large datasets
    options(
        future.globals.maxSize = 2000 * 1024^3
    )
    
    proj <- addHarmony(
        ArchRProj   = proj,
        reducedDims = "IterativeLSI",
        name        = "Harmony",
        groupBy     = "DatasetGroup"
    )
    
    proj <- saveArchRProject(proj, outputDirectory = step5_dir, load = TRUE)
    message(sprintf("[5] Saved: %s", step5_dir))
}

# ==============================================================================
# 6: SEURAT METADATA INTEGRATION
# ==============================================================================

banner("6: Seurat Metadata Integration")

step6_dir <- file.path(OUTPUT_BASE, "Step6_MetadataAdded")
ensure_dir(step6_dir)

if (file.exists(file.path(step6_dir, "ArchRProject.rds"))) {
    proj <- loadArchRProject(step6_dir)
    message(sprintf("[6] Loaded existing project"))
} else {
    message(sprintf("[6] Integrating Seurat metadata..."))
    
    # Convert h5ad to Seurat using schard
    # Install: remotes::install_github("cellgeni/schard")
    seurat_obj <- schard::h5ad2seurat(H5AD_PATH)
    seurat_meta <- seurat_obj@meta.data
    seurat_meta$barcode <- rownames(seurat_meta)
    
    # Convert barcodes to ArchR format
    seurat_meta <- seurat_meta %>%
        mutate(archr_cell = sapply(barcode, convert_barcode_to_archr))
    
    # Find matching cells
    matched_cells <- intersect(seurat_meta$archr_cell, getCellNames(proj))
    message(sprintf("[6] Matched cells: %d", length(matched_cells)))
    
    # Subset ArchR project to matched cells
    proj <- subsetArchRProject(
        ArchRProj       = proj,
        cells           = matched_cells,
        outputDirectory = file.path(step6_dir, "ArchRProject_Subset"),
        force           = TRUE
    )
    
    # Prepare metadata for transfer
    seurat_sub <- seurat_meta %>% filter(archr_cell %in% matched_cells)
    rownames(seurat_sub) <- seurat_sub$archr_cell
    
    # Transfer metadata columns
    cols_to_add <- c("celltype", "sample", "celltype2", "age", "sex")
    
    for (colname in cols_to_add) {
        proj <- addCellColData(
            ArchRProj = proj,
            data      = seurat_sub[[colname]],
            cells     = rownames(seurat_sub),
            name      = colname,
            force     = TRUE
        )
        message(sprintf("[6] Added metadata: %s", colname))
    }
    
    proj <- saveArchRProject(proj, outputDirectory = step6_dir, load = TRUE)
    message(sprintf("[6] Saved: %s", step6_dir))
    
    # Clean up
    rm(seurat_obj, seurat_meta, seurat_sub)
    gc()
}

# ==============================================================================
# 7: WNN UMAP EMBEDDING TRANSFER
# ==============================================================================

banner("7: WNN UMAP Embedding Transfer")

step7_dir <- file.path(OUTPUT_BASE, "Step7_Xwnn_UMAP")
ensure_dir(step7_dir)

if (file.exists(file.path(step7_dir, "ArchRProject.rds"))) {
    proj <- loadArchRProject(step7_dir)
    message(sprintf("[7] Loaded existing project"))
} else {
    message(sprintf("[7] Transferring WNN UMAP embedding..."))
    
    # Convert h5ad to Seurat using schard
    seurat_obj <- schard::h5ad2seurat(H5AD_PATH)
    
    # Extract WNN UMAP embedding
    seurat_embed <- Embeddings(seurat_obj, reduction = "Xwnn_") %>%
        as.data.frame() %>%
        setNames(c("UMAP_1", "UMAP_2"))
    
    seurat_embed$barcode <- rownames(seurat_embed)
    
    # Convert barcodes to ArchR format
    seurat_embed <- seurat_embed %>%
        mutate(archr_cell = sapply(barcode, convert_barcode_to_archr))
    
    # Filter for matched cells
    matched_df <- seurat_embed %>% filter(archr_cell %in% getCellNames(proj))
    
    # Transfer UMAP coordinates
    for (dim in c("UMAP_1", "UMAP_2")) {
        vec <- matched_df[[dim]]
        names(vec) <- matched_df$archr_cell
        
        proj <- addCellColData(
            ArchRProj = proj,
            data      = vec,
            cells     = names(vec),
            name      = paste0("Xwnn_", dim)
        )
        message(sprintf("[7] Transferred: %s", dim))
    }
    
    proj <- saveArchRProject(proj, outputDirectory = step7_dir, load = TRUE)
    message(sprintf("[7] Saved: %s", step7_dir))
    
    # Clean up
    rm(seurat_obj, seurat_embed, matched_df)
    gc()
}

# ==============================================================================
# PIPELINE SUMMARY
# ==============================================================================

banner("PREPROCESSING PIPELINE COMPLETE")

message("Summary:")
message(sprintf("  Final cell count: %d", nCells(proj)))
message(sprintf("  Samples: %d", length(unique(proj$Sample))))
message(sprintf("  Age groups: %s", paste(unique(proj$DatasetGroup), collapse = ", ")))

if ("celltype" %in% colnames(proj@cellColData)) {
    message("\n  Cell types:")
    print(table(proj$celltype))
}

message(sprintf("\n  Final project location: %s", step7_dir))

message("\nNext step: Run archr_downstream_analysis.R or archr_bigwig_generation.R")

# Save session info for reproducibility
sessionInfo()
