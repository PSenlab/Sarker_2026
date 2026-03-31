#!/usr/bin/env Rscript
## =========================================================
## COMPARTMENT ANALYSIS - SHARED CONFIGURATION
## 
## This script defines all shared parameters, color palettes,
## and helper functions used across the compartment analysis
## pipeline.
##
## Source this file at the beginning of each analysis script:
##   source("compartment_00_config.R")
##
## =========================================================

## =========================================================
## SECTION 1: PATHS AND PARAMETERS
## =========================================================

# Project paths (modify for your environment)
CONFIG <- list(
    # Input paths
    proj_path       = "path/to/ArchR_Projects/Step6_Xwnn_UMAP",
    mm10_blacklist  = "path/to/mm10-blacklist.v2.bed.gz",
    chromhmm_bed    = "path/to/mm10_100_segments_segments.bed.gz",
    
    # Output directory
    outdir          = "path/to/compartment_output",
    
    # Analysis parameters
    tile_size       = 80000L,
    non_gap_thresh  = 0.75,
    k_bins          = 25L,
    threads         = 60L,
    
    # Stability thresholds
    t_on            = 0.6,   # Fraction for "Active" classification
    t_off           = 0.4,   # Fraction for "Repressive" classification
    
    # Random seed
    seed            = 10918
)

# Derive additional paths
CONFIG$group_dir <- file.path(CONFIG$proj_path, "GroupBigWigs", "AgeSexCelltype__tile80k")
CONFIG$manifest_path <- file.path(CONFIG$group_dir, "AgeSexCelltype.bigwig_manifest_tile80k.tsv")


## =========================================================
## SECTION 2: FACTOR LEVELS
## =========================================================

# Age levels (in biological order)
AGE_LEVELS <- c("young", "mid_age", "old", "pre_geriatric", "geriatric")

# Sex levels
SEX_LEVELS <- c("male", "female", "unknown")

# Stability classes (5-class model)
STABILITY_LEVELS <- c(
    "Stable_Active",
    "Stable_Repressive", 
    "Monotonic_A_to_R",
    "Monotonic_R_to_A",
    "Non_Monotonic"
)

# Switching classes only
SWITCHING_CLASSES <- c("Monotonic_A_to_R", "Monotonic_R_to_A", "Non_Monotonic")

# Non-monotonic pattern types
NONMONO_PATTERNS <- c("Transient_Return", "Two_Switch", "Triple_Switch", "Highly_Dynamic", "Other")


## =========================================================
## SECTION 3: COLOR PALETTES
## =========================================================

# Age colors
AGE_COLS <- c(
    "young"         = "#1ABC9C",
    "mid_age"       = "#F1C40F",
    "old"           = "#C39BD3",
    "pre_geriatric" = "#2980B9",
    "geriatric"     = "#E84393"
)

# Sex colors
SEX_COLS <- c(
    "male"    = "#2E86C1",
    "female"  = "#EC7063",
    "unknown" = "#B2BABB"
)

# Cell type colors
CELLTYPE_COLS <- c(
    "Hepatocyte"       = "#17becf",
    "Endothelial.01"   = "#a6cee3",
    "Endothelial.02"   = "#1f78b4",
    "Cholangiocyte.01" = "#e31a1c",
    "Cholangiocyte.02" = "#cab2d6",
    "Kupffer"          = "#b2df8a",
    "MoMFs"            = "#33a02c",
    "Stellate"         = "#fb9a99",
    "lymp_T"           = "#fdbf6f",
    "lymp_B"           = "#e377c2"
)

# Compartment state colors (for alluvial strata)
STATE_COLS <- c(
    "Active"     = "#FF6B6B",
    "Repressive" = "#4E79A7"
)

# 5-class stability colors
STABILITY_COLS <- c(
    "Stable_Active"     = "#FF6B6B",
    "Stable_Repressive" = "#4E79A7",
    "Monotonic_A_to_R"  = "darkred",
    "Monotonic_R_to_A"  = "darkblue",
    "Non_Monotonic"     = "#708238",
    "Missing"           = "grey80",
    "Insufficient_Data" = "grey60"
)

# Switching-only colors (subset)
SWITCHING_COLS <- c(
    "Monotonic_A_to_R" = "darkred",
    "Monotonic_R_to_A" = "darkblue",
    "Non_Monotonic"    = "#708238"
)

# Non-monotonic pattern colors
NONMONO_PATTERN_COLS <- c(
    "Transient_Return" = "#2E8B57",
    "Two_Switch"       = "#E67E22",
    "Triple_Switch"    = "#CD5C5C",
    "Highly_Dynamic"   = "#DDA0DD",
    "Other"            = "#95A5A6"
)

# Compartment colors (for heatmap)
COMPARTMENT_COLS <- c(
    "0" = "#0A2A43",   # Repressive = dark navy
    "1" = "#F2D28B"    # Active = tan yellow
)

# Chromosome colors
chr_order <- paste0("chr", c(1:19, "X", "Y"))
n_chr <- length(chr_order)
hues <- seq(0, 360 - 360/n_chr, length.out = n_chr)
CHR_COLS <- setNames(grDevices::hcl(h = hues, c = 90, l = 60), chr_order)


## =========================================================
## SECTION 4: HELPER FUNCTIONS
## =========================================================

#' Print a formatted banner for pipeline steps
#' @param text Text to display in banner
#' @param char Character to use for border (default "=")
#' @param width Width of banner (default 70)
banner <- function(text, char = "=", width = 70) {
    line <- paste(rep(char, width), collapse = "")
    message("\n", line)
    message(text)
    message(line, "\n")
}

#' Print a fancy box banner
#' @param text Text to display
box_banner <- function(text) {
    cat("\n")
    cat("╔══════════════════════════════════════════════════════════════════╗\n")
    cat(sprintf("║ %-66s ║\n", text))
    cat("╚══════════════════════════════════════════════════════════════════╝\n\n")
}

#' Ensure directory exists
#' @param dir_path Path to directory
#' @return Invisibly returns the directory path
ensure_dir <- function(dir_path) {
    if (!dir.exists(dir_path)) {
        dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
    }
    invisible(dir_path)
}

#' Clean cell type name for use in file/directory names
#' @param ct Cell type name
#' @return Cleaned name safe for file paths
clean_name <- function(ct) {
    gsub("[^A-Za-z0-9_]", "_", ct)
}

#' Compute fraction active per category
#' @param comp_bin Compartment binary matrix (bins × samples)
#' @param meta_vec Named vector of metadata (e.g., age_vec)
#' @param levels_use Levels to compute (e.g., AGE_LEVELS)
#' @return Matrix of fraction active (bins × levels)
frac_active_by <- function(comp_bin, meta_vec, levels_use) {
    stopifnot(all(colnames(comp_bin) %in% names(meta_vec)))
    
    out <- lapply(levels_use, function(lv) {
        idx <- which(meta_vec[colnames(comp_bin)] == lv)
        if (length(idx) == 0) return(rep(NA_real_, nrow(comp_bin)))
        rowMeans(comp_bin[, idx, drop = FALSE] == 1L, na.rm = TRUE)
    })
    
    M <- do.call(cbind, out)
    colnames(M) <- levels_use
    rownames(M) <- rownames(comp_bin)
    return(M)
}

#' Write BED file from data.table
#' @param dt Data table with chr, start, end columns
#' @param filename Output file path
#' @param name_col Column to use for BED name field
#' @param score_col Column to use for BED score field (optional)
write_bed <- function(dt, filename, name_col = "bin_id", score_col = NULL) {
    bed <- data.table::data.table(
        chr   = dt$chr,
        start = dt$start - 1,  # BED is 0-based
        end   = dt$end,
        name  = dt[[name_col]]
    )
    if (!is.null(score_col) && score_col %in% names(dt)) {
        bed[, score := dt[[score_col]]]
    } else {
        bed[, score := 0]
    }
    bed[, strand := "."]
    
    data.table::fwrite(bed, filename, sep = "\t", col.names = FALSE)
    message(sprintf("    → %s (%d regions)", basename(filename), nrow(bed)))
}


## =========================================================
## SECTION 5: VALIDATION FUNCTIONS
## =========================================================

#' Check if required objects exist in environment
#' @param required_objects Character vector of object names
#' @param stop_on_missing If TRUE, stop with error; if FALSE, return missing names
check_required_objects <- function(required_objects, stop_on_missing = TRUE) {
    missing <- required_objects[!sapply(required_objects, exists, envir = .GlobalEnv)]
    
    if (length(missing) > 0) {
        msg <- sprintf("Missing required objects: %s", paste(missing, collapse = ", "))
        if (stop_on_missing) {
            stop(msg)
        } else {
            warning(msg)
            return(missing)
        }
    }
    return(character(0))
}

#' Validate compartment binary matrix
#' @param comp_bin Compartment binary matrix
#' @return Logical indicating if valid
validate_comp_bin <- function(comp_bin) {
    checks <- list(
        "is matrix" = is.matrix(comp_bin),
        "has rownames" = !is.null(rownames(comp_bin)),
        "has colnames" = !is.null(colnames(comp_bin)),
        "values are 0/1/NA" = all(comp_bin %in% c(0L, 1L, NA), na.rm = TRUE)
    )
    
    failed <- names(checks)[!unlist(checks)]
    if (length(failed) > 0) {
        warning(sprintf("comp_bin validation failed: %s", paste(failed, collapse = ", ")))
        return(FALSE)
    }
    return(TRUE)
}


## =========================================================
## SECTION 6: INITIALIZE
## =========================================================

# Set seed for reproducibility
set.seed(CONFIG$seed)

# Create output directory
ensure_dir(CONFIG$outdir)

# Print configuration
message("\n[CONFIG] Compartment Analysis Configuration Loaded")
message(sprintf("  Output directory: %s", CONFIG$outdir))
message(sprintf("  Tile size: %d bp", CONFIG$tile_size))
message(sprintf("  K-means clusters: %d", CONFIG$k_bins))
message(sprintf("  Threads: %d", CONFIG$threads))
