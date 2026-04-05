#!/usr/bin/env Rscript
# ==============================================================================
# ArchR Peak Calling, Co-Accessibility, and Peak-to-Gene Linkage (Step 11)
# ==============================================================================
#
# Description:
#   This script performs MACS2 peak calling, co-accessibility analysis, and
#   peak-to-gene linkage for each cell type. Processes all cell types in a
#   loop, creating separate projects per cell type.
#
# Prerequisites:
#   - Run archr_arrow_creation.R (Step 1)
#   - Run archr_preprocessing.R (Steps 2-7)
#   - Run archr_downstream_analysis.R (Step 8)
#   - Run archr_bigwig_generation.R (Step 9) [optional]
#   - Run archr_rna_integration.R (Step 10)
#
# Input:
#   - ArchR project with RNA integration (from Step 10)
#
# Output:
#   - Per-celltype ArchR projects with:
#     - MACS2 peak calls (grouped by sex_age)
#     - Peak matrix
#     - Co-accessibility network
#     - Peak-to-gene linkages
#
# Pipeline Overview:
#   11.1: Load ArchR Project from Step 10
#   11.2: Define Cell Types to Process
#   For each cell type:
#     11.3: Subset Cells by Cell Type
#     11.4: Create sex_age Grouping Variable
#     11.5: Add Group Coverages
#     11.6: Call Peaks with MACS2
#     11.7: Add Peak Matrix
#     11.8: Save Intermediate Project (Step11a)
#     11.9: Add Co-Accessibility
#     11.10: Add Peak-to-Gene Links
#     11.11: Save Final Project (Step11b)
#
# Requirements:
#   - R >= 4.0
#   - ArchR >= 1.0.2
#   - MACS2 installed and accessible
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
    library(dplyr)
    library(stringr)
    library(igraph)
    library(Matrix)
    library(matrixStats)
    library(SummarizedExperiment)
    library(BSgenome.Mmusculus.UCSC.mm10)
})

# Set Global Parameters
addArchRThreads(threads = 60)
addArchRGenome("mm10")

# Define Paths (modify according to your directory structure)
STEP10_PROJ_PATH <- "path/to/ArchR_Projects/Step10_RNA_Integrated"
OUTPUT_BASE <- "path/to/ArchR_Projects"
PATH_TO_MACS2 <- "path/to/macs2"  # e.g., /path/to/conda/envs/env_name/bin/macs2

# Processing Parameters
MIN_CELLS_PER_CELLTYPE <- 100    # Skip cell types with fewer cells
CELLTYPES_TO_EXCLUDE <- c()      # Add cell types to skip, e.g., c("Unassigned")

# Peak Calling Parameters (MACS2)
PEAK_CUTOFF <- 0.01              # FDR cutoff for peak calling
MAX_PEAKS <- 500000              # Maximum peaks per group

# Group Coverage Parameters
MIN_CELLS_COVERAGE <- 40         # Minimum cells per group for coverage
MAX_CELLS_COVERAGE <- 10000      # Maximum cells per group for coverage
MAX_FRAGMENTS <- 25e6            # Maximum fragments per group
MIN_REPLICATES <- 2              # Minimum replicates per group
MAX_REPLICATES <- 40             # Maximum replicates per group

# Co-Accessibility Parameters
COA_COR_CUTOFF <- 0.4            # Correlation cutoff for co-accessibility
COA_K <- 100                     # Number of nearest neighbors
COA_MAX_DIST <- 1000000          # Maximum distance (1 Mb)

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

#' Clean cell type name for use in directory names
#' @param ct Cell type name
#' @return Cleaned name safe for file paths
clean_celltype_name <- function(ct) {
    gsub("[^A-Za-z0-9_]", "_", ct)
}

# ==============================================================================
# 11.1: LOAD ARCHR PROJECT
# ==============================================================================

banner("11.1: Load ArchR Project from Step 10")

if (!dir.exists(STEP10_PROJ_PATH)) {
    stop(sprintf("[FATAL] Project path does not exist: %s", STEP10_PROJ_PATH))
}

proj_main <- loadArchRProject(path = STEP10_PROJ_PATH)
message(sprintf("[11.1] Loaded project from Step 10: %s", STEP10_PROJ_PATH))
message(sprintf("[11.1] Total cells: %d", nCells(proj_main)))

# Validate required metadata
stopifnot("celltype" %in% colnames(proj_main@cellColData))
stopifnot("sex" %in% colnames(proj_main@cellColData))
stopifnot("age" %in% colnames(proj_main@cellColData))

# Validate GeneExpressionMatrix exists (required for P2G)
available_matrices <- getAvailableMatrices(proj_main)
if (!"GeneExpressionMatrix" %in% available_matrices) {
    stop("[FATAL] GeneExpressionMatrix not found. Run archr_rna_integration.R (Step 10) first.")
}

message(sprintf("[11.1] Available matrices: %s", paste(available_matrices, collapse = ", ")))

# ==============================================================================
# 11.2: DEFINE CELL TYPES TO PROCESS
# ==============================================================================

banner("11.2: Define Cell Types to Process")

all_celltypes <- unique(proj_main$celltype)
celltypes_to_process <- setdiff(all_celltypes, CELLTYPES_TO_EXCLUDE)

message(sprintf("[11.2] Total cell types found: %d", length(all_celltypes)))
message(sprintf("[11.2] Cell types to process: %d", length(celltypes_to_process)))

# Print cell counts per cell type
celltype_counts <- table(proj_main$celltype)
for (ct in celltypes_to_process) {
    n_cells <- celltype_counts[ct]
    status <- if (n_cells >= MIN_CELLS_PER_CELLTYPE) "OK" else "SKIP"
    message(sprintf("  - %s: %d cells [%s]", ct, n_cells, status))
}

# ==============================================================================
# 11.3-11.11: PROCESS EACH CELL TYPE
# ==============================================================================

banner("11.3-11.11: Process Each Cell Type")

results_summary <- data.frame(
    celltype = character(),
    n_cells = integer(),
    n_peaks = integer(),
    n_p2g_links = integer(),
    status = character(),
    stringsAsFactors = FALSE
)

for (ct in celltypes_to_process) {
    
    message("\n", strrep("=", 70))
    message(sprintf(">>> PROCESSING CELL TYPE: %s", ct))
    message(strrep("=", 70))
    
    ct_clean <- clean_celltype_name(ct)
    
    tryCatch({
        
        # -----------------------------------------------------------------
        # 11.3: SUBSET CELLS BY CELL TYPE
        # -----------------------------------------------------------------
        message(sprintf("[11.3] Subsetting %s cells", ct))
        
        ct_cells <- rownames(proj_main@cellColData)[proj_main$celltype == ct]
        
        if (length(ct_cells) < MIN_CELLS_PER_CELLTYPE) {
            message(sprintf("[11.3] SKIPPING %s: only %d cells (need >= %d)", 
                            ct, length(ct_cells), MIN_CELLS_PER_CELLTYPE))
            
            results_summary <- rbind(results_summary, data.frame(
                celltype = ct,
                n_cells = length(ct_cells),
                n_peaks = NA,
                n_p2g_links = NA,
                status = "SKIPPED_LOW_CELLS"
            ))
            next
        }
        
        step11a_dir <- file.path(OUTPUT_BASE, paste0("Step11a_", ct_clean, "_Peaks"))
        ensure_dir(step11a_dir)
        
        proj <- subsetArchRProject(
            ArchRProj = proj_main,
            cells = ct_cells,
            outputDirectory = step11a_dir,
            dropCells = TRUE,
            force = TRUE
        )
        
        message(sprintf("[11.3] Subset project: %d cells", nCells(proj)))
        
        # -----------------------------------------------------------------
        # 11.4: CREATE SEX_AGE GROUPING VARIABLE
        # -----------------------------------------------------------------
        message(sprintf("[11.4] Creating sex_age grouping variable"))
        
        proj$sex_age <- paste0(proj$sex, "_", proj$age)
        
        message(sprintf("[11.4] sex_age groups:"))
        print(table(proj$sex_age))
        
        # -----------------------------------------------------------------
        # 11.5: ADD GROUP COVERAGES
        # -----------------------------------------------------------------
        message(sprintf("[11.5] Adding group coverages"))
        
        proj <- addGroupCoverages(
            ArchRProj     = proj,
            groupBy       = "sex_age",
            useLabels     = TRUE,
            minCells      = MIN_CELLS_COVERAGE,
            maxCells      = MAX_CELLS_COVERAGE,
            maxFragments  = MAX_FRAGMENTS,
            minReplicates = MIN_REPLICATES,
            maxReplicates = MAX_REPLICATES,
            sampleRatio   = 0.8,
            kmerLength    = 6,
            threads       = getArchRThreads(),
            force         = TRUE
        )
        
        message(sprintf("[11.5] Group coverages added"))
        
        # -----------------------------------------------------------------
        # 11.6: CALL PEAKS WITH MACS2
        # -----------------------------------------------------------------
        message(sprintf("[11.6] Calling peaks with MACS2"))
        
        proj <- addReproduciblePeakSet(
            ArchRProj   = proj,
            groupBy     = "sex_age",
            pathToMacs2 = PATH_TO_MACS2,
            cutOff      = PEAK_CUTOFF,
            maxPeaks    = MAX_PEAKS,
            plot        = FALSE,
            force       = TRUE
        )
        
        peak_set <- getPeakSet(proj)
        n_peaks <- length(peak_set)
        message(sprintf("[11.6] Peaks called: %d", n_peaks))
        
        # -----------------------------------------------------------------
        # 11.7: ADD PEAK MATRIX
        # -----------------------------------------------------------------
        message(sprintf("[11.7] Adding peak matrix"))
        
        proj <- addPeakMatrix(proj)
        
        message(sprintf("[11.7] Peak matrix added"))
        
        # -----------------------------------------------------------------
        # 11.8: SAVE INTERMEDIATE PROJECT (STEP 11a)
        # -----------------------------------------------------------------
        message(sprintf("[11.8] Saving intermediate project"))
        
        proj <- saveArchRProject(
            ArchRProj = proj,
            outputDirectory = step11a_dir,
            load = TRUE
        )
        
        message(sprintf("[11.8] Saved: %s", step11a_dir))
        
        # Reload to fix Arrow indexing issues
        proj <- loadArchRProject(step11a_dir)
        
        # -----------------------------------------------------------------
        # 11.9: ADD CO-ACCESSIBILITY
        # -----------------------------------------------------------------
        message(sprintf("[11.9] Adding co-accessibility"))
        
        proj <- addCoAccessibility(
            ArchRProj   = proj,
            reducedDims = "IterativeLSI",
            corCutOff   = COA_COR_CUTOFF,
            k           = COA_K,
            maxDist     = COA_MAX_DIST
        )
        
        # Validate co-accessibility was added
        coa <- getCoAccessibility(proj)
        if (is.null(coa)) {
            warning(sprintf("[11.9] Co-accessibility not added for %s", ct))
        } else {
            message(sprintf("[11.9] Co-accessibility added: %d links", length(coa)))
        }
        
        # -----------------------------------------------------------------
        # 11.10: ADD PEAK-TO-GENE LINKS
        # -----------------------------------------------------------------
        message(sprintf("[11.10] Adding peak-to-gene links"))
        
        proj <- addPeak2GeneLinks(
            ArchRProj   = proj,
            reducedDims = "IterativeLSI",
            useMatrix   = "GeneExpressionMatrix"
        )
        
        # Get P2G link count
        p2g <- getPeak2GeneLinks(proj)
        n_p2g <- if (!is.null(p2g)) length(p2g) else 0
        message(sprintf("[11.10] Peak-to-gene links added: %d", n_p2g))
        
        # -----------------------------------------------------------------
        # 11.11: SAVE FINAL PROJECT (STEP 11b)
        # -----------------------------------------------------------------
        message(sprintf("[11.11] Saving final project"))
        
        step11b_dir <- file.path(OUTPUT_BASE, paste0("Step11b_", ct_clean, "_P2G"))
        ensure_dir(step11b_dir)
        
        proj <- saveArchRProject(
            ArchRProj = proj,
            outputDirectory = step11b_dir,
            load = TRUE
        )
        
        message(sprintf("[11.11] Saved: %s", step11b_dir))
        
        # Record results
        results_summary <- rbind(results_summary, data.frame(
            celltype = ct,
            n_cells = nCells(proj),
            n_peaks = n_peaks,
            n_p2g_links = n_p2g,
            status = "SUCCESS"
        ))
        
        message(sprintf("[11.11] COMPLETED %s: %d cells, %d peaks, %d P2G links", 
                        ct, nCells(proj), n_peaks, n_p2g))
        
        # Clean up memory
        rm(proj)
        gc()
        
    }, error = function(e) {
        message(sprintf("[ERROR] Processing %s failed: %s", ct, e$message))
        
        results_summary <<- rbind(results_summary, data.frame(
            celltype = ct,
            n_cells = NA,
            n_peaks = NA,
            n_p2g_links = NA,
            status = paste0("ERROR: ", e$message)
        ))
    })
}

# ==============================================================================
# PIPELINE SUMMARY
# ==============================================================================

banner("STEP 11 COMPLETE: PEAK CALLING AND P2G LINKAGE")

message("Summary by Cell Type:")
print(results_summary)

# Save summary
summary_file <- file.path(OUTPUT_BASE, "Step11_summary.tsv")
write.table(
    results_summary,
    file = summary_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)
message(sprintf("\nSummary saved: %s", summary_file))

# Overall statistics
n_success <- sum(results_summary$status == "SUCCESS", na.rm = TRUE)
n_skipped <- sum(grepl("SKIPPED", results_summary$status), na.rm = TRUE)
n_error <- sum(grepl("ERROR", results_summary$status), na.rm = TRUE)

message(sprintf("\nOverall Results:"))
message(sprintf("  Successful: %d", n_success))
message(sprintf("  Skipped: %d", n_skipped))
message(sprintf("  Errors: %d", n_error))

message("\nOutput directories:")
message(sprintf("  Peaks (11a): %s/Step11a_*_Peaks/", OUTPUT_BASE))
message(sprintf("  P2G (11b): %s/Step11b_*_P2G/", OUTPUT_BASE))

# Save session info for reproducibility
sessionInfo()
