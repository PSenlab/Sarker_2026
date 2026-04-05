#!/usr/bin/env Rscript
# ==============================================================================
# ArchR Arrow File Creation for single-nucleus ATAC-seq data
# ==============================================================================
#
#
# Input:
#   - CellRanger Arc output (atac_fragments.tsv.gz files)
#
# Output:
#   - Arrow files (one per sample)
#   - QC plots organized by age group
#
#
# Citation:
#   - ArchR: Granja et al., Nature Genetics 2021
#
# ==============================================================================

# ==============================================================================
# SETUP AND CONFIGURATION
# ==============================================================================

# Load Required Libraries
suppressPackageStartupMessages({
    library(ArchR)
})

# Set Global Parameters
set.seed(10918)
addArchRThreads(threads = 60)
addArchRGenome("mm10")

# Define Paths (modify according to your directory structure)
WORK_DIR <- "path/to/working/directory"
CELLRANGER_BASE <- "path/to/cellranger_arc/outputs"

setwd(WORK_DIR)

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
# SAMPLE METADATA DEFINITION
# ==============================================================================

# Base directories for fragment files (per age group)
BASES <- list(
    young         = file.path(CELLRANGER_BASE, "Young"),
    mid_age       = file.path(CELLRANGER_BASE, "Middle_age"),
    old           = file.path(CELLRANGER_BASE, "Old"),
    pre_geriatric = file.path(CELLRANGER_BASE, "Pre_Geriatric"),
    geriatric     = file.path(CELLRANGER_BASE, "Geriatric")
)

# Folder names inside each base (Cell Ranger output folders)
FOLDERS <- list(
    young         = paste0("Y", 1:8),
    mid_age       = paste0("MA_", sprintf("%02d", 1:8)),
    old           = paste0("O", 1:8),
    pre_geriatric = paste0("PG_", sprintf("%02d", 1:8)),
    geriatric     = paste0("G", 1:8)
)

# Sample names for ArchR
SAMPLES <- list(
    young         = paste0("young_",     sprintf("%02d", 1:8)),
    mid_age       = paste0("mid_age_",   sprintf("%02d", 1:8)),
    old           = paste0("old_",       sprintf("%02d", 1:8)),
    pre_geriatric = paste0("pre_ger_",   sprintf("%02d", 1:8)),
    geriatric     = paste0("geriatric_", sprintf("%02d", 1:8))
)

# Age group order for processing
GROUP_ORDER <- c("young", "mid_age", "old", "pre_geriatric", "geriatric")

# Arrow file QC parameters
QC_PARAMS <- list(
    minTSS   = 2,
    minFrags = 1000
)

# ==============================================================================
# 1: BUILD INPUT FILES
# ==============================================================================

banner("1: Build Input File List")

inputFiles <- character(0)

for (grp in names(BASES)) {
    base_dir <- BASES[[grp]]
    folders  <- FOLDERS[[grp]]
    samples  <- SAMPLES[[grp]]
    
    stopifnot(length(folders) == length(samples))
    
    frag_paths <- file.path(base_dir, folders, "outs", "atac_fragments.tsv.gz")
    
    # Check for existing files
    exists_vec <- file.exists(frag_paths)
    if (any(!exists_vec)) {
        missing <- paste(samples[!exists_vec], collapse = ", ")
        warning(sprintf("[WARN] Missing fragments for group '%s': %s", grp, missing))
    }
    
    frag_paths <- frag_paths[exists_vec]
    samp_names <- samples[exists_vec]
    
    if (length(frag_paths)) {
        names(frag_paths) <- samp_names
        inputFiles <- c(inputFiles, frag_paths)
    }
}

if (length(inputFiles) == 0L) {
    stop("[FATAL] No fragment files found. Check CELLRANGER_BASE path.")
}

message(sprintf("[1] Found %d fragment files", length(inputFiles)))

# ==============================================================================
# 2: CREATE ARROW FILES BY AGE GROUP
# ==============================================================================

banner("2: Create Arrow Files")

ArrowFiles <- character(0)

for (grp in GROUP_ORDER) {
    # Match samples for this group by sample name prefix
    if (grp == "pre_geriatric") {
        grp_samples <- grep("^pre_ger_", names(inputFiles), value = TRUE)
    } else {
        grp_samples <- grep(paste0("^", grp, "_"), names(inputFiles), value = TRUE)
    }
    
    if (length(grp_samples) == 0L) {
        message(sprintf("[2] No samples for group '%s' - skipping", grp))
        next
    }
    
    files_this_grp   <- inputFiles[grp_samples]
    samples_this_grp <- names(files_this_grp)
    
    # QC directory for this group
    qc_dir <- file.path(WORK_DIR, "QualityControl", grp)
    ensure_dir(qc_dir)
    
    message(sprintf("[2] Processing group: %-14s | Samples: %d | QC: %s",
                    grp, length(files_this_grp), qc_dir))
    
    # Create Arrow files
    af <- tryCatch({
        createArrowFiles(
            inputFiles      = files_this_grp,
            sampleNames     = samples_this_grp,
            minTSS          = QC_PARAMS$minTSS,
            minFrags        = QC_PARAMS$minFrags,
            addTileMat      = TRUE,
            addGeneScoreMat = TRUE,
            QCDir           = qc_dir,
            threads         = getArchRThreads()
        )
    }, error = function(e) {
        warning(sprintf("[ERROR] createArrowFiles failed for group '%s': %s", grp, e$message))
        return(NULL)
    })
    
    if (!is.null(af)) {
        ArrowFiles <- c(ArrowFiles, af)
    }
}

ArrowFiles <- unique(ArrowFiles)

if (length(ArrowFiles) == 0L) {
    stop("[FATAL] No Arrow files were created. Check logs above.")
}

# ==============================================================================
# PIPELINE SUMMARY
# ==============================================================================

banner("ARROW FILE CREATION COMPLETE")

message("Summary:")
message(sprintf("  Total Arrow files: %d", length(ArrowFiles)))

# Count by group
by_grp_counts <- sapply(GROUP_ORDER, function(g) {
    if (g == "pre_geriatric") {
        sum(grepl("^pre_ger_", basename(ArrowFiles)))
    } else {
        sum(grepl(paste0("^", g, "_"), basename(ArrowFiles)))
    }
})
names(by_grp_counts) <- GROUP_ORDER

message("\n  Arrow files per age group:")
for (grp in GROUP_ORDER) {
    message(sprintf("    %s: %d", grp, by_grp_counts[grp]))
}

message(sprintf("\n  QC output: %s", file.path(WORK_DIR, "QualityControl")))
message(sprintf("  Arrow files: %s", WORK_DIR))

message("\nNext step: Run archr_preprocessing.R")

# Close any open graphics devices
while (dev.cur() > 1) dev.off()

# Save session info
sessionInfo()
