#!/usr/bin/env Rscript
# ==============================================================================
# Generate Group-Level BigWigs (Age × Sex × Celltype) using 80 kb Bins (Step 9)
# ==============================================================================
#
# Description:
#   This script generates normalized BigWig files for chromatin accessibility
#   grouped by Age, Sex, and Cell Type combinations. Uses 80 kb bins for
#   downstream compartment switching analysis.
#
# Prerequisites:
#   - Run archr_arrow_creation.R (Step 1)
#   - Run archr_preprocessing.R (Steps 2-7)
#   - Run archr_downstream_analysis.R (Step 8)
#
# Input:
#   - Preprocessed ArchR project with imputation (from Step 8)
#
# Output:
#   - BigWig files per Age × Sex × Celltype group
#   - Group counts table
#   - BigWig manifest file
#
# Pipeline Overview:
#   9.1: Load ArchR Project from Step 8
#   9.2: Validate Metadata
#   9.3: Add Group Labels (Age__Sex__Celltype)
#   9.4: Filter Groups by Minimum Cell Count
#   9.5: Generate Group Count Tables
#   9.6: Generate BigWigs (80 kb Tiles)
#   9.7: Write Output Files
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
    library(dplyr)
    library(stringr)
    library(tidyr)
    library(tibble)
    library(BSgenome.Mmusculus.UCSC.mm10)
})

# Define Paths (modify according to your directory structure)
# Load project from Step 8 of archr_downstream_analysis.R
STEP8_PROJ_PATH <- "path/to/ArchR_Projects/Step8_Imputed"

# BigWig Generation Parameters
TILE_SIZE   <- 80000L      # Fixed at 80 kb for compartment analysis
NORM_METHOD <- "nFrags"    # Normalize by total fragments per cell
MAX_CELLS   <- 50000       # Maximum cells to sample per group
CEILING_VAL <- 4           # Clip extreme coverage values
THREADS     <- 60          # Parallel threads
MIN_CELLS   <- 1           # Skip groups with fewer cells

# Naming Conventions
GROUP_DIR_ROOT <- "AgeSexCelltype"
LOG_BASE       <- "getGroupBW_AgeSexCelltype"
MANIFEST_BASE  <- "AgeSexCelltype.bigwig_manifest_tile80k.tsv"
COUNTS_BASE    <- "AgeSexCelltype.group_counts_tile80k.tsv"

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

#' Get project save path with fallback
#' @param proj ArchR project
#' @param fallback Fallback path if save path not found
.get_save_path <- function(proj, fallback = STEP8_PROJ_PATH) {
    sp <- tryCatch(
        proj@projectMetadata$SavePath,
        error = function(e) NA_character_
    )
    if (is.null(sp) || is.na(sp) || !nzchar(sp)) fallback else sp
}

#' Extract group label from BigWig filename
#' @param fn Filename
#' @param tile_size Tile size used
.label_from_filename <- function(fn, tile_size) {
    base <- basename(fn)
    base <- sub("\\.bw$", "", base)
    sub(paste0("-TileSize-", tile_size, "-.*$"), "", base)
}

#' Format tile size for display
#' @param x Tile size in bp
.pretty_tile <- function(x) {
    if (x %% 1000L == 0L) paste0(x %/% 1000L, "k") else as.character(x)
}

#' Get cell names from project
#' @param proj ArchR project
.get_cells <- function(proj) {
    rownames(getCellColData(proj))
}

# ==============================================================================
# 9.1: LOAD ARCHR PROJECT
# ==============================================================================

banner("9.1: Load ArchR Project from Step 8")

addArchRThreads(THREADS)
addArchRGenome("mm10")

if (!dir.exists(STEP8_PROJ_PATH)) {
    stop(sprintf("[FATAL] Project path does not exist: %s", STEP8_PROJ_PATH))
}

proj <- loadArchRProject(STEP8_PROJ_PATH)
proj_dir <- .get_save_path(proj)

message(sprintf("[9.1] Loaded project from Step 8: %s", STEP8_PROJ_PATH))
message(sprintf("[9.1] Total cells: %d", nCells(proj)))

# ==============================================================================
# 9.2: VALIDATE METADATA
# ==============================================================================

banner("9.2: Validate Metadata")

meta <- as.data.frame(getCellColData(proj))
need_cols <- c("age", "sex", "celltype")
missing <- setdiff(need_cols, colnames(meta))

if (length(missing)) {
    stop(sprintf("[FATAL] Missing metadata columns: %s", paste(missing, collapse = ", ")))
}

meta <- meta %>%
    mutate(across(all_of(need_cols), as.character))

message(sprintf("[9.2] Unique cell types (n=%d):", length(unique(meta$celltype))))
for (ct in sort(unique(meta$celltype))) {
    message(sprintf('  - "%s"', ct))
}

# ==============================================================================
# 9.3: ADD GROUP LABELS
# ==============================================================================

banner("9.3: Add Group Labels (Age__Sex__Celltype)")

cells0 <- .get_cells(proj)
canon_group_vec <- paste(meta$age, meta$sex, meta$celltype, sep = "__")
stopifnot(length(canon_group_vec) == length(cells0))

proj <- addCellColData(
    ArchRProj = proj,
    data      = canon_group_vec,
    name      = "AgeSexCelltype",
    cells     = cells0,
    force     = TRUE
)

message(sprintf("[9.3] Added AgeSexCelltype labels to %d cells", length(cells0)))

# ==============================================================================
# 9.4: FILTER BY MINIMUM CELLS (OPTIONAL)
# ==============================================================================

banner("9.4: Filter Groups by Minimum Cell Count")

grp_vals0 <- getCellColData(proj, select = "AgeSexCelltype")[, 1]
grp_tbl_pre <- sort(table(grp_vals0), decreasing = TRUE)

if (MIN_CELLS > 0) {
    keep_groups <- names(grp_tbl_pre)[grp_tbl_pre >= MIN_CELLS]
    ac_df <- getCellColData(proj, select = "AgeSexCelltype")
    cells_keep <- rownames(ac_df)[!is.na(ac_df$AgeSexCelltype) & 
                                   ac_df$AgeSexCelltype %in% keep_groups]
    
    if (length(cells_keep) > 0 && length(cells_keep) < nrow(ac_df)) {
        proj <- subsetArchRProject(
            proj,
            cells           = cells_keep,
            outputDirectory = paste0("subset_", GROUP_DIR_ROOT),
            force           = TRUE
        )
        proj_dir <- .get_save_path(
            proj,
            fallback = file.path(getwd(), paste0("subset_", GROUP_DIR_ROOT))
        )
        
        # Re-add group labels after subsetting
        meta2 <- as.data.frame(getCellColData(proj))
        meta2 <- meta2 %>% mutate(across(all_of(need_cols), as.character))
        cells_sub <- .get_cells(proj)
        canon_group_vec2 <- paste(meta2$age, meta2$sex, meta2$celltype, sep = "__")
        
        proj <- addCellColData(
            ArchRProj = proj,
            data      = canon_group_vec2,
            name      = "AgeSexCelltype",
            cells     = cells_sub,
            force     = TRUE
        )
        
        message(sprintf("[9.4] Filtered to %d cells (min_cells=%d)", 
                        length(cells_keep), MIN_CELLS))
    }
} else {
    message(sprintf("[9.4] No filtering applied (min_cells=%d)", MIN_CELLS))
}

# ==============================================================================
# 9.5: GENERATE COUNT TABLES
# ==============================================================================

banner("9.5: Generate Group Count Tables")

grp_vals_used <- getCellColData(proj, select = "AgeSexCelltype")[, 1]
grp_tbl_used <- sort(table(grp_vals_used), decreasing = TRUE)

split_counts_used <- tibble(
    group   = names(grp_tbl_used),
    n_cells = as.integer(grp_tbl_used)
) %>%
    separate(group, into = c("age", "sex", "celltype"), sep = "__", remove = FALSE) %>%
    select(age, sex, celltype, n_cells, group)

message(sprintf("[9.5] Generated counts for %d groups", nrow(split_counts_used)))

# ==============================================================================
# 9.6: GENERATE BIGWIGS (80 KB TILES)
# ==============================================================================

banner("9.6: Generate BigWigs (80 kb Tiles)")

ts <- TILE_SIZE
ts_tag <- .pretty_tile(ts)
group_col <- paste0(GROUP_DIR_ROOT, "__tile", ts_tag)
log_name <- paste0(LOG_BASE, "_tile", ts_tag, ".log")

archr_out_dir <- file.path(proj_dir, "GroupBigWigs", group_col)
dir.create(archr_out_dir, showWarnings = FALSE, recursive = TRUE)

cells_now <- .get_cells(proj)
agesex_vec <- as.character(getCellColData(proj, select = "AgeSexCelltype")[, 1])
stopifnot(length(agesex_vec) == length(cells_now))

proj <- addCellColData(
    ArchRProj = proj,
    data      = agesex_vec,
    name      = group_col,
    cells     = cells_now,
    force     = TRUE
)

message(sprintf("[9.6] Generating 80 kb BigWigs (groupBy=%s)", group_col))

bw_files <- getGroupBW(
    ArchRProj  = proj,
    groupBy    = group_col,
    normMethod = NORM_METHOD,
    tileSize   = ts,
    maxCells   = MAX_CELLS,
    ceiling    = CEILING_VAL,
    threads    = getArchRThreads(),
    logFile    = createLogFile(log_name)
)

# Collect BigWig file paths
bw_vec <- tryCatch(unlist(bw_files), error = function(e) character(0))
if (length(bw_vec) == 0 && dir.exists(archr_out_dir)) {
    pat <- paste0("-TileSize-", ts, "-.*ArchR\\.bw$")
    bw_vec <- list.files(archr_out_dir, pattern = pat, full.names = TRUE)
}
if (length(bw_vec) == 0) {
    stop(sprintf("[FATAL] No BigWigs found for %s in %s", group_col, archr_out_dir))
}

grp_labels <- vapply(bw_vec, .label_from_filename, character(1), tile_size = ts)

message(sprintf("[9.6] Generated %d BigWig files", length(bw_vec)))

# ==============================================================================
# 9.7: WRITE OUTPUT FILES
# ==============================================================================

banner("9.7: Write Output Files")

# Write counts table
counts_file <- file.path(archr_out_dir, COUNTS_BASE)
write.table(
    split_counts_used,
    file      = counts_file,
    sep       = "\t",
    quote     = FALSE,
    row.names = FALSE
)
message(sprintf("[9.7] Saved: %s", counts_file))

# Write manifest
manifest <- tibble(
    group  = grp_labels,
    bigwig = bw_vec
) %>%
    separate(group, into = c("age", "sex", "celltype"), sep = "__", remove = FALSE) %>%
    select(age, sex, celltype, group, bigwig)

manifest_file <- file.path(archr_out_dir, MANIFEST_BASE)
write.table(
    manifest,
    file      = manifest_file,
    sep       = "\t",
    quote     = FALSE,
    row.names = FALSE
)
message(sprintf("[9.7] Saved: %s", manifest_file))

# ==============================================================================
# PIPELINE SUMMARY
# ==============================================================================

banner("STEP 9 COMPLETE: BIGWIG GENERATION")

message("Summary:")
message(sprintf("  Tile size: %s", ts_tag))
message(sprintf("  Groups processed: %d", length(bw_vec)))
message(sprintf("  Output directory: %s", archr_out_dir))
message(sprintf("  Manifest: %s", manifest_file))
message(sprintf("  Counts: %s", counts_file))

# Save session info for reproducibility
sessionInfo()
