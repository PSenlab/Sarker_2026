#!/usr/bin/env Rscript
# ==============================================================================
# Peak-to-Gene Linkage Overlap with Compartment Stability Classes
# ==============================================================================
#
# Description:
#   Overlaps Hepatocyte P2G-linked peaks with compartment stability bins
#   (Stable_Active, Stable_Repressive, Monotonic_A_to_R, Monotonic_R_to_A,
#   Non_Monotonic) for each sex. Runs at two correlation cutoffs (default
#   0.45, relaxed 0.25) to assess sensitivity.
#
# Input:
#   - ArchR project with Hepatocyte P2G linkages (from peak calling step)
#   - compartment_stability_5class.tsv (from compartment_01_main_analysis.R)
#
# Output (per cutoff):
#   - Gene lists per sex x stability class (.txt)
#   - Peak info per sex x stability class (.tsv)
#   - Stability bins per sex x stability class (.rds, .bed)
#   - Overlap statistics summary (.tsv)
#   - R objects for downstream analysis (.rds)
#
# Reference:
#   Sarker et al. (2026) Cell Metabolism
#
# ==============================================================================


# ==============================================================================
# SETUP
# ==============================================================================

suppressPackageStartupMessages({
    library(ArchR)
    library(data.table)
    library(GenomicRanges)
    library(IRanges)
    library(S4Vectors)
    library(dplyr)
})

# ------------------------------------------------------------------------------
# Configuration - Update these paths
# ------------------------------------------------------------------------------

# Compartment analysis output directory
compartment_outdir <- "path/to/output/downstream_stability"

# Hepatocyte ArchR project with P2G linkages
archr_proj_path <- "path/to/ArchR_Projects/Step8_Hepatocyte_CCAN_P2G"

# Stability classification file
stability_file <- "path/to/output/compartment_stability_5class.tsv"

# Threads
addArchRThreads(threads = 16)

# ------------------------------------------------------------------------------
# Analysis parameters
# ------------------------------------------------------------------------------

P2G_CUTOFFS <- list(
    default = 0.45,
    relaxed = 0.25
)

STABILITY_CLASSES <- c(
    "Stable_Active",
    "Stable_Repressive",
    "Monotonic_A_to_R",
    "Monotonic_R_to_A",
    "Non_Monotonic"
)

SEXES <- c("male", "female")


# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================

get_count <- function(tab, key) {
    if (!is.null(tab) && key %in% names(tab)) as.integer(tab[[key]]) else 0L
}

split_genes <- function(x) {
    x <- as.character(x)
    x <- x[!is.na(x) & nzchar(x)]
    if (length(x) == 0) return(character(0))
    g <- unlist(strsplit(x, ";", fixed = TRUE), use.names = FALSE)
    g <- trimws(g)
    g <- g[nzchar(g)]
    unique(g)
}

parse_bin_id <- function(bin_ids) {
    bin_ids <- as.character(bin_ids)
    bad <- is.na(bin_ids) | !grepl("^[^:]+:[0-9]+-[0-9]+$", bin_ids)
    if (any(bad)) {
        stop("Invalid bin_id format detected. Examples: ",
             paste(head(bin_ids[bad], 5), collapse = ", "))
    }
    parts <- strsplit(bin_ids, ":", fixed = TRUE)
    chr <- vapply(parts, `[`, character(1), 1)
    coords <- vapply(parts, `[`, character(1), 2)
    coord_parts <- strsplit(coords, "-", fixed = TRUE)
    start_pos <- as.integer(vapply(coord_parts, `[`, character(1), 1))
    end_pos   <- as.integer(vapply(coord_parts, `[`, character(1), 2))
    data.table(chr = chr, start = start_pos, end = end_pos)
}

banner <- function(text) {
    line <- paste(rep("=", 70), collapse = "")
    message("\n", line)
    message(text)
    message(line, "\n")
}


# ==============================================================================
# STEP 1: LOAD ARCHR PROJECT
# ==============================================================================

banner("STEP 1: Load ArchR project")

proj <- loadArchRProject(archr_proj_path)
message(sprintf("  [OK] Project loaded: %s", proj@projectMetadata$outputDirectory))
message(sprintf("  [OK] Cells: %d", nCells(proj)))

projPeakSet <- getPeakSet(proj)


# ==============================================================================
# STEP 2: LOAD COMPARTMENT STABILITY DATA
# ==============================================================================

banner("STEP 2: Load compartment stability data")

if (!file.exists(stability_file)) {
    stop("Stability file not found: ", stability_file)
}

stability_dt <- fread(stability_file)
setnames(stability_dt, trimws(names(stability_dt)))

stability_dt[, Sex := tolower(trimws(as.character(Sex)))]
stability_dt[, Celltype := trimws(as.character(Celltype))]
stability_dt[, Stability := trimws(as.character(Stability))]
stability_dt[, bin_id := trimws(as.character(bin_id))]

# Filter for Hepatocyte
hep_dt <- stability_dt[Celltype == "Hepatocyte"]

message(sprintf("  [OK] Hepatocyte bins: %s", format(nrow(hep_dt), big.mark = ",")))

message("\n  Stability distribution:")
print(hep_dt[, .N, by = Stability][order(-N)])

message("\n  By sex:")
print(hep_dt[, .N, by = .(Sex, Stability)][order(Sex, -N)])


# ==============================================================================
# STEP 3: CREATE STABILITY GRANGES
# ==============================================================================

banner("STEP 3: Create stability GRanges")

stability_gr_list <- list()

for (sx in SEXES) {
    message(sprintf("\n  %s:", toupper(sx)))
    stability_gr_list[[sx]] <- list()

    for (stab_class in STABILITY_CLASSES) {
        bins_dt <- hep_dt[Sex == sx & Stability == stab_class]

        if (nrow(bins_dt) == 0) {
            message(sprintf("    %s: 0 bins (skipping)", stab_class))
            next
        }

        coords <- parse_bin_id(bins_dt$bin_id)

        gr <- GRanges(
            seqnames  = coords$chr,
            ranges    = IRanges(start = coords$start, end = coords$end),
            bin_id    = bins_dt$bin_id,
            Sex       = sx,
            Stability = stab_class
        )

        stability_gr_list[[sx]][[stab_class]] <- gr
        message(sprintf("    %s: %d bins", stab_class, length(gr)))
    }
}


# ==============================================================================
# STEP 4: OVERLAP P2G PEAKS WITH STABILITY BINS (BOTH CUTOFFS)
# ==============================================================================

for (cutoff_name in names(P2G_CUTOFFS)) {

    corCutOff_value <- P2G_CUTOFFS[[cutoff_name]]

    banner(sprintf("PROCESSING: %s (corCutOff = %.2f)",
                   toupper(cutoff_name), corCutOff_value))

    # Output directory for this cutoff
    p2g_outdir <- file.path(compartment_outdir, cutoff_name,
                            "peak2gene_stability_analysis")
    dir.create(p2g_outdir, recursive = TRUE, showWarnings = FALSE)
    message(sprintf("  Output: %s", p2g_outdir))

    # ------------------------------------------------------------------
    # 4a: Extract P2G links
    # ------------------------------------------------------------------

    message("\n  Extracting P2G links...")

    p2g <- getPeak2GeneLinks(
        ArchRProj  = proj,
        corCutOff  = corCutOff_value,
        resolution = 1,
        returnLoops = FALSE
    )

    p2g$geneName <- mcols(metadata(p2g)$geneSet)$name[p2g$idxRNA]
    p2g$peakName <- (metadata(p2g)$peakSet %>% {
        paste0(seqnames(.), "_", start(.), "_", end(.))
    })[p2g$idxATAC]

    p2g_dt <- as.data.table(as.data.frame(p2g))

    # Parse peak coordinates
    p2g_dt[, peak_chr   := gsub("_.*", "", peakName)]
    p2g_dt[, peak_start := as.integer(gsub("^[^_]+_([0-9]+)_.*", "\\1", peakName))]
    p2g_dt[, peak_end   := as.integer(gsub("^.*_", "", peakName))]
    p2g_dt[, gene := geneName]
    p2g_dt[, peak := paste0(peak_chr, ":", peak_start, "-", peak_end)]

    # Add peak annotations
    p2g_dt[, peakType    := projPeakSet$peakType[idxATAC]]
    p2g_dt[, nearestGene := projPeakSet$nearestGene[idxATAC]]
    p2g_dt[, distToTSS   := projPeakSet$distToTSS[idxATAC]]

    message(sprintf("  [OK] P2G links: %s", format(nrow(p2g_dt), big.mark = ",")))
    message(sprintf("  [OK] Unique peaks: %s", format(length(unique(p2g_dt$peak)), big.mark = ",")))
    message(sprintf("  [OK] Unique genes: %s", format(length(unique(p2g_dt$gene)), big.mark = ",")))

    saveRDS(p2g_dt, file.path(p2g_outdir, "hepatocyte_p2g_links.rds"))
    fwrite(p2g_dt, file.path(p2g_outdir, "hepatocyte_p2g_links.tsv"), sep = "\t")

    # ------------------------------------------------------------------
    # 4b: Create peak GRanges
    # ------------------------------------------------------------------

    message("\n  Creating peak GRanges...")

    peaks_unique <- p2g_dt[, .(
        genes    = paste(unique(gene), collapse = ";"),
        n_genes  = length(unique(gene)),
        mean_cor = mean(Correlation),
        peakType = peakType[1]
    ), by = .(peak, peak_chr, peak_start, peak_end)]

    peak_gr <- GRanges(
        seqnames = peaks_unique$peak_chr,
        ranges   = IRanges(start = peaks_unique$peak_start, end = peaks_unique$peak_end),
        peak_id  = peaks_unique$peak,
        genes    = peaks_unique$genes,
        n_genes  = peaks_unique$n_genes,
        mean_cor = peaks_unique$mean_cor,
        peakType = peaks_unique$peakType
    )

    message(sprintf("  [OK] Peak GRanges: %d peaks", length(peak_gr)))
    saveRDS(peak_gr, file.path(p2g_outdir, "hepatocyte_p2g_peaks_gr.rds"))

    # ------------------------------------------------------------------
    # 4c: Overlap peaks with stability bins
    # ------------------------------------------------------------------

    message("\n  Overlapping peaks with stability bins...")

    genes_by_sex_class <- list()
    peaks_by_sex_class <- list()
    overlap_stats_list <- list()

    for (sx in SEXES) {

        message(sprintf("\n  %s:", toupper(sx)))

        genes_by_sex_class[[sx]] <- list()
        peaks_by_sex_class[[sx]] <- list()

        for (stab_class in STABILITY_CLASSES) {

            bins_gr <- stability_gr_list[[sx]][[stab_class]]

            # Handle empty bins
            if (is.null(bins_gr) || length(bins_gr) == 0) {
                genes_by_sex_class[[sx]][[stab_class]] <- character(0)
                peaks_by_sex_class[[sx]][[stab_class]] <- character(0)
                overlap_stats_list[[paste(sx, stab_class, sep = "_")]] <- data.table(
                    Sex = sx, Stability = stab_class,
                    n_bins = 0L, n_peaks_overlap = 0L, n_genes_linked = 0L,
                    n_promoter = 0L, n_distal = 0L, n_intronic = 0L, n_exonic = 0L
                )
                next
            }

            olaps <- findOverlaps(peak_gr, bins_gr, ignore.strand = TRUE)

            # Handle no overlaps
            if (length(olaps) == 0) {
                genes_by_sex_class[[sx]][[stab_class]] <- character(0)
                peaks_by_sex_class[[sx]][[stab_class]] <- character(0)
                overlap_stats_list[[paste(sx, stab_class, sep = "_")]] <- data.table(
                    Sex = sx, Stability = stab_class,
                    n_bins = as.integer(length(bins_gr)),
                    n_peaks_overlap = 0L, n_genes_linked = 0L,
                    n_promoter = 0L, n_distal = 0L, n_intronic = 0L, n_exonic = 0L
                )
                next
            }

            peak_hits <- unique(queryHits(olaps))
            overlapping_peaks <- peak_gr[peak_hits]
            linked_genes <- split_genes(overlapping_peaks$genes)

            genes_by_sex_class[[sx]][[stab_class]] <- linked_genes
            peaks_by_sex_class[[sx]][[stab_class]] <- as.character(overlapping_peaks$peak_id)

            # Peak type breakdown
            pt <- as.character(overlapping_peaks$peakType)
            pt <- ifelse(is.na(pt) | !nzchar(pt), "Unknown", pt)
            peak_types <- table(pt)

            nP <- get_count(peak_types, "Promoter")
            nD <- get_count(peak_types, "Distal")
            nI <- get_count(peak_types, "Intronic")
            nE <- get_count(peak_types, "Exonic")

            message(sprintf("    %s: Bins=%d, Peaks=%d, Genes=%d",
                            stab_class, length(bins_gr),
                            length(peak_hits), length(linked_genes)))

            # Save gene list
            writeLines(linked_genes, file.path(
                p2g_outdir, sprintf("genes_%s_%s.txt", sx, stab_class)
            ))

            # Save peak info
            peak_info <- data.table(
                peak     = as.character(overlapping_peaks$peak_id),
                genes    = as.character(overlapping_peaks$genes),
                n_genes  = as.integer(overlapping_peaks$n_genes),
                mean_cor = round(as.numeric(overlapping_peaks$mean_cor), 4),
                peakType = as.character(overlapping_peaks$peakType)
            )
            fwrite(peak_info, file.path(
                p2g_outdir, sprintf("peaks_%s_%s.tsv", sx, stab_class)
            ), sep = "\t")

            # Save bins (RDS + BED)
            saveRDS(bins_gr, file.path(
                p2g_outdir, sprintf("bins_%s_%s.rds", sx, stab_class)
            ))

            bed_dt <- data.table(
                chr    = as.character(seqnames(bins_gr)),
                start  = start(bins_gr) - 1L,
                end    = end(bins_gr),
                name   = bins_gr$bin_id,
                score  = 0,
                strand = "."
            )
            fwrite(bed_dt, file.path(
                p2g_outdir, sprintf("bins_%s_%s.bed", sx, stab_class)
            ), sep = "\t", col.names = FALSE)

            overlap_stats_list[[paste(sx, stab_class, sep = "_")]] <- data.table(
                Sex = sx, Stability = stab_class,
                n_bins          = as.integer(length(bins_gr)),
                n_peaks_overlap = as.integer(length(peak_hits)),
                n_genes_linked  = as.integer(length(linked_genes)),
                n_promoter = nP, n_distal = nD, n_intronic = nI, n_exonic = nE
            )
        }
    }

    # Combine stats
    overlap_stats <- rbindlist(overlap_stats_list, use.names = TRUE, fill = TRUE)
    den <- fifelse(overlap_stats$n_peaks_overlap > 0L, overlap_stats$n_peaks_overlap, 1L)
    overlap_stats[, pct_promoter := round(100 * n_promoter / den, 1)]
    overlap_stats[, pct_distal   := round(100 * n_distal / den, 1)]

    fwrite(overlap_stats, file.path(p2g_outdir, "overlap_statistics_all.tsv"), sep = "\t")

    # ------------------------------------------------------------------
    # 4d: Save R objects
    # ------------------------------------------------------------------

    message("\n  Saving R objects...")

    saveRDS(genes_by_sex_class, file.path(p2g_outdir, "genes_by_sex_class.rds"))
    saveRDS(peaks_by_sex_class, file.path(p2g_outdir, "peaks_by_sex_class.rds"))
    saveRDS(stability_gr_list,  file.path(p2g_outdir, "stability_gr_list.rds"))
    saveRDS(overlap_stats,      file.path(p2g_outdir, "overlap_stats.rds"))
    saveRDS(p2g_dt,             file.path(p2g_outdir, "p2g_dt.rds"))
    saveRDS(peak_gr,            file.path(p2g_outdir, "peak_gr.rds"))

    # ------------------------------------------------------------------
    # Summary for this cutoff
    # ------------------------------------------------------------------

    banner(sprintf("%s (corCutOff = %.2f) SUMMARY", toupper(cutoff_name), corCutOff_value))

    message(sprintf("  P2G links: %s", format(nrow(p2g_dt), big.mark = ",")))
    message(sprintf("  Unique peaks: %s", format(length(unique(p2g_dt$peak)), big.mark = ",")))
    message(sprintf("  Unique genes: %s", format(length(unique(p2g_dt$gene)), big.mark = ",")))
    message(sprintf("  Output: %s", p2g_outdir))

    message("\n  Gene counts per Sex x Stability:")
    for (sx in SEXES) {
        message(sprintf("  %s:", toupper(sx)))
        for (stab_class in STABILITY_CLASSES) {
            n_genes <- length(genes_by_sex_class[[sx]][[stab_class]])
            message(sprintf("    %-20s: %5d genes", stab_class, n_genes))
        }
    }
}


# ==============================================================================
# SESSION INFO
# ==============================================================================

writeLines(
    capture.output(sessionInfo()),
    file.path(compartment_outdir, "sessionInfo_p2g_stability.txt")
)


# ==============================================================================
# FINAL SUMMARY
# ==============================================================================

cat("\n")
cat("========================================================================\n")
cat("  COMPLETE - BOTH CUTOFFS PROCESSED\n")
cat("========================================================================\n")
cat("\n")
cat("  Output directories:\n")
cat(sprintf("    DEFAULT (0.45): %s/default/peak2gene_stability_analysis\n", compartment_outdir))
cat(sprintf("    RELAXED (0.25): %s/relaxed/peak2gene_stability_analysis\n", compartment_outdir))
cat("\n")
cat("  Files per cutoff:\n")
cat("    genes_[sex]_[stability].txt\n")
cat("    peaks_[sex]_[stability].tsv\n")
cat("    bins_[sex]_[stability].rds/bed\n")
cat("    hepatocyte_p2g_links.rds/tsv\n")
cat("    hepatocyte_p2g_peaks_gr.rds\n")
cat("    genes_by_sex_class.rds\n")
cat("    peaks_by_sex_class.rds\n")
cat("    overlap_stats.rds\n")
cat("    overlap_statistics_all.tsv\n")
cat("\n")
cat("  Reproducibility:\n")
cat("    sessionInfo_p2g_stability.txt\n")
cat("========================================================================\n")
