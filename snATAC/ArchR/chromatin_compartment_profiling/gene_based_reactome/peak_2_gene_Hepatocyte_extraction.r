#!/usr/bin/env Rscript
## =========================================================
## OVERLAP ALL STABILITY CLASSES WITH HEPATOCYTE P2G
## TWO VERSIONS: DEFAULT (corCutOff=0.45) and RELAXED (corCutOff=0.25)
## =========================================================

suppressPackageStartupMessages({
    library(ArchR)
    library(data.table)
    library(GenomicRanges)
    library(IRanges)
    library(S4Vectors)
    library(dplyr)
})

## =========================================================
## PATHS - UPDATE THESE
## =========================================================

# Compartment analysis output
compartment_outdir <- "/data/sarkern2/multiome_liver/Seurat/epigenome/downstream_stability"

# Hepatocyte ArchR project with P2G
archr_proj_path <- "/data/sarkern2/multiome_liver/Seurat/archR/ArchR_Projects/Step8_Hepatocyte_CCAN_P2G"

# Stability file
stability_file <- "/data/sarkern2/multiome_liver/Seurat/epigenome/downstream_stability/newfolder/compartment_stability_5class.tsv"

# Threads
addArchRThreads(threads = 16)

## =========================================================
## DEFINE P2G CUTOFFS
## =========================================================

P2G_CUTOFFS <- list(
    default = 0.45,
    relaxed = 0.25
)

## =========================================================
## DEFINE ALL STABILITY CLASSES
## =========================================================

STABILITY_CLASSES <- c(
    "Stable_Active",
    "Stable_Repressive",
    "Monotonic_A_to_R",
    "Monotonic_R_to_A",
    "Non_Monotonic"
)

SEXES <- c("male", "female")

## =========================================================
## HELPER FUNCTIONS
## =========================================================

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
    start <- as.integer(vapply(coord_parts, `[`, character(1), 1))
    end   <- as.integer(vapply(coord_parts, `[`, character(1), 2))
    data.table(chr = chr, start = start, end = end)
}

## =========================================================
## STEP 1: LOAD ARCHR PROJECT
## =========================================================

message("\n")
message("╔══════════════════════════════════════════════════════════════════╗")
message("║  STEP 1: LOAD ARCHR PROJECT                                      ║")
message("╚══════════════════════════════════════════════════════════════════╝\n")

proj <- loadArchRProject(archr_proj_path)
message(sprintf("  ✓ Project loaded: %s", proj@projectMetadata$outputDirectory))
message(sprintf("  ✓ Cells: %d", nCells(proj)))

# Get peak annotations from project
projPeakSet <- getPeakSet(proj)

## =========================================================
## STEP 2: LOAD COMPARTMENT STABILITY DATA
## =========================================================

message("\n")
message("╔══════════════════════════════════════════════════════════════════╗")
message("║  STEP 2: LOAD COMPARTMENT STABILITY DATA                         ║")
message("╚══════════════════════════════════════════════════════════════════╝\n")

if (!file.exists(stability_file)) {
    stop("Stability file not found: ", stability_file)
}

stability_dt <- fread(stability_file)
setnames(stability_dt, trimws(names(stability_dt)))

# Normalize columns
stability_dt[, Sex := tolower(trimws(as.character(Sex)))]
stability_dt[, Celltype := trimws(as.character(Celltype))]
stability_dt[, Stability := trimws(as.character(Stability))]
stability_dt[, bin_id := trimws(as.character(bin_id))]

# Filter for Hepatocyte
hep_dt <- stability_dt[Celltype == "Hepatocyte"]

message(sprintf("  ✓ Total Hepatocyte bins: %s", format(nrow(hep_dt), big.mark = ",")))

message("\n  Stability distribution:")
print(hep_dt[, .N, by = Stability][order(-N)])

message("\n  By Sex:")
print(hep_dt[, .N, by = .(Sex, Stability)][order(Sex, -N)])

## =========================================================
## STEP 3: CREATE STABILITY GRANGES (SHARED)
## =========================================================

message("\n")
message("╔══════════════════════════════════════════════════════════════════╗")
message("║  STEP 3: CREATE STABILITY GRANGES                                ║")
message("╚══════════════════════════════════════════════════════════════════╝\n")

stability_gr_list <- list()

for (sx in SEXES) {
    message(sprintf("\n>>> Processing: %s", toupper(sx)))
    stability_gr_list[[sx]] <- list()
    
    for (stab_class in STABILITY_CLASSES) {
        bins_dt <- hep_dt[Sex == sx & Stability == stab_class]
        
        if (nrow(bins_dt) == 0) {
            message(sprintf("    %s: 0 bins (skipping)", stab_class))
            next
        }
        
        coords <- parse_bin_id(bins_dt$bin_id)
        
        gr <- GRanges(
            seqnames = coords$chr,
            ranges = IRanges(start = coords$start, end = coords$end),
            bin_id = bins_dt$bin_id,
            Sex = sx,
            Stability = stab_class
        )
        
        stability_gr_list[[sx]][[stab_class]] <- gr
        message(sprintf("    %s: %d bins", stab_class, length(gr)))
    }
}

## =========================================================
## STEP 4: LOOP OVER BOTH CUTOFFS
## =========================================================

for (cutoff_name in names(P2G_CUTOFFS)) {
    
    corCutOff_value <- P2G_CUTOFFS[[cutoff_name]]
    
    message("\n")
    message("╔══════════════════════════════════════════════════════════════════╗")
    message(sprintf("║  PROCESSING: %s (corCutOff = %.2f)                            ║", 
                    toupper(cutoff_name), corCutOff_value))
    message("╚══════════════════════════════════════════════════════════════════╝\n")
    
    # Output directory for this cutoff
    p2g_outdir <- file.path(compartment_outdir, cutoff_name, "peak2gene_stability_analysis")
    dir.create(p2g_outdir, recursive = TRUE, showWarnings = FALSE)
    
    message(sprintf("  Output directory: %s", p2g_outdir))
    
    ## =========================================================
    ## EXTRACT P2G DATA
    ## =========================================================
    
    message("\n>>> Extracting P2G links...")
    
    p2g <- getPeak2GeneLinks(
        ArchRProj = proj,
        corCutOff = corCutOff_value,
        resolution = 1,
        returnLoops = FALSE
    )
    
    p2g$geneName <- mcols(metadata(p2g)$geneSet)$name[p2g$idxRNA]
    p2g$peakName <- (metadata(p2g)$peakSet %>% {paste0(seqnames(.), "_", start(.), "_", end(.))})[p2g$idxATAC]
    
    # Build P2G data table
    p2g_df <- as.data.frame(p2g)
    p2g_dt <- as.data.table(p2g_df)
    
    # Parse peak coordinates
    p2g_dt[, peak_chr := gsub("_.*", "", peakName)]
    p2g_dt[, peak_start := as.integer(gsub("^[^_]+_([0-9]+)_.*", "\\1", peakName))]
    p2g_dt[, peak_end := as.integer(gsub("^.*_", "", peakName))]
    p2g_dt[, gene := geneName]
    p2g_dt[, peak := paste0(peak_chr, ":", peak_start, "-", peak_end)]
    
    # Add peak annotations
    p2g_dt[, peakType := projPeakSet$peakType[idxATAC]]
    p2g_dt[, nearestGene := projPeakSet$nearestGene[idxATAC]]
    p2g_dt[, distToTSS := projPeakSet$distToTSS[idxATAC]]
    
    message(sprintf("  ✓ P2G links: %s", format(nrow(p2g_dt), big.mark = ",")))
    message(sprintf("  ✓ Unique peaks: %s", format(length(unique(p2g_dt$peak)), big.mark = ",")))
    message(sprintf("  ✓ Unique genes: %s", format(length(unique(p2g_dt$gene)), big.mark = ",")))
    
    # Save P2G data
    saveRDS(p2g_dt, file.path(p2g_outdir, "hepatocyte_p2g_links.rds"))
    fwrite(p2g_dt, file.path(p2g_outdir, "hepatocyte_p2g_links.tsv"), sep = "\t")
    
    ## =========================================================
    ## CREATE PEAK GRANGES
    ## =========================================================
    
    message("\n>>> Creating peak GRanges...")
    
    peaks_unique <- p2g_dt[, .(
        genes = paste(unique(gene), collapse = ";"),
        n_genes = length(unique(gene)),
        mean_cor = mean(Correlation),
        peakType = peakType[1]
    ), by = .(peak, peak_chr, peak_start, peak_end)]
    
    peak_gr <- GRanges(
        seqnames = peaks_unique$peak_chr,
        ranges = IRanges(start = peaks_unique$peak_start, end = peaks_unique$peak_end),
        peak_id = peaks_unique$peak,
        genes = peaks_unique$genes,
        n_genes = peaks_unique$n_genes,
        mean_cor = peaks_unique$mean_cor,
        peakType = peaks_unique$peakType
    )
    
    message(sprintf("  ✓ Peak GRanges: %d peaks", length(peak_gr)))
    
    saveRDS(peak_gr, file.path(p2g_outdir, "hepatocyte_p2g_peaks_gr.rds"))
    
    ## =========================================================
    ## OVERLAP PEAKS WITH STABILITY BINS
    ## =========================================================
    
    message("\n>>> Overlapping peaks with stability bins...")
    
    genes_by_sex_class <- list()
    peaks_by_sex_class <- list()
    overlap_stats_list <- list()
    
    for (sx in SEXES) {
        
        message(sprintf("\n  %s:", toupper(sx)))
        
        genes_by_sex_class[[sx]] <- list()
        peaks_by_sex_class[[sx]] <- list()
        
        for (stab_class in STABILITY_CLASSES) {
            
            bins_gr <- stability_gr_list[[sx]][[stab_class]]
            
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
                            stab_class, length(bins_gr), length(peak_hits), length(linked_genes)))
            
            # Save gene list
            gene_file <- file.path(p2g_outdir, sprintf("genes_%s_%s.txt", sx, stab_class))
            writeLines(linked_genes, gene_file)
            
            # Save peak info
            peak_info <- data.table(
                peak = as.character(overlapping_peaks$peak_id),
                genes = as.character(overlapping_peaks$genes),
                n_genes = as.integer(overlapping_peaks$n_genes),
                mean_cor = round(as.numeric(overlapping_peaks$mean_cor), 4),
                peakType = as.character(overlapping_peaks$peakType)
            )
            peak_file <- file.path(p2g_outdir, sprintf("peaks_%s_%s.tsv", sx, stab_class))
            fwrite(peak_info, peak_file, sep = "\t")
            
            # Save bins
            rds_file <- file.path(p2g_outdir, sprintf("bins_%s_%s.rds", sx, stab_class))
            saveRDS(bins_gr, rds_file)
            
            bed_dt <- data.table(
                chr = as.character(seqnames(bins_gr)),
                start = start(bins_gr) - 1L,
                end = end(bins_gr),
                name = bins_gr$bin_id,
                score = 0,
                strand = "."
            )
            bed_file <- file.path(p2g_outdir, sprintf("bins_%s_%s.bed", sx, stab_class))
            fwrite(bed_dt, bed_file, sep = "\t", col.names = FALSE)
            
            overlap_stats_list[[paste(sx, stab_class, sep = "_")]] <- data.table(
                Sex = sx, Stability = stab_class,
                n_bins = as.integer(length(bins_gr)),
                n_peaks_overlap = as.integer(length(peak_hits)),
                n_genes_linked = as.integer(length(linked_genes)),
                n_promoter = nP, n_distal = nD, n_intronic = nI, n_exonic = nE
            )
        }
    }
    
    # Combine and save stats
    overlap_stats <- rbindlist(overlap_stats_list, use.names = TRUE, fill = TRUE)
    den <- fifelse(overlap_stats$n_peaks_overlap > 0L, overlap_stats$n_peaks_overlap, 1L)
    overlap_stats[, pct_promoter := round(100 * n_promoter / den, 1)]
    overlap_stats[, pct_distal := round(100 * n_distal / den, 1)]
    
    fwrite(overlap_stats, file.path(p2g_outdir, "overlap_statistics_all.tsv"), sep = "\t")
    
    ## =========================================================
    ## SAVE R OBJECTS
    ## =========================================================
    
    message("\n>>> Saving R objects...")
    
    saveRDS(genes_by_sex_class, file.path(p2g_outdir, "genes_by_sex_class.rds"))
    saveRDS(peaks_by_sex_class, file.path(p2g_outdir, "peaks_by_sex_class.rds"))
    saveRDS(stability_gr_list, file.path(p2g_outdir, "stability_gr_list.rds"))
    saveRDS(overlap_stats, file.path(p2g_outdir, "overlap_stats.rds"))
    saveRDS(p2g_dt, file.path(p2g_outdir, "p2g_dt.rds"))
    saveRDS(peak_gr, file.path(p2g_outdir, "peak_gr.rds"))
    
    ## =========================================================
    ## SUMMARY FOR THIS CUTOFF
    ## =========================================================
    
    message("\n")
    message("═══════════════════════════════════════════════════════════════════")
    message(sprintf("  %s (corCutOff = %.2f) SUMMARY", toupper(cutoff_name), corCutOff_value))
    message("═══════════════════════════════════════════════════════════════════")
    message(sprintf("  P2G links: %s", format(nrow(p2g_dt), big.mark = ",")))
    message(sprintf("  Unique peaks: %s", format(length(unique(p2g_dt$peak)), big.mark = ",")))
    message(sprintf("  Unique genes: %s", format(length(unique(p2g_dt$gene)), big.mark = ",")))
    message(sprintf("  Output: %s", p2g_outdir))
    
    message("\n  Gene counts per Sex × Stability:")
    for (sx in SEXES) {
        message(sprintf("  %s:", toupper(sx)))
        for (stab_class in STABILITY_CLASSES) {
            n_genes <- length(genes_by_sex_class[[sx]][[stab_class]])
            message(sprintf("    %-20s: %5d genes", stab_class, n_genes))
        }
    }
}

## =========================================================
## FINAL SUMMARY
## =========================================================

cat("\n")
cat("╔══════════════════════════════════════════════════════════════════╗\n")
cat("║              BOTH CUTOFFS COMPLETE                               ║\n")
cat("╠══════════════════════════════════════════════════════════════════╣\n")
cat("║                                                                  ║\n")
cat("║  Output directories:                                             ║\n")
cat(sprintf("║    DEFAULT (0.45): %s/default/peak2gene_stability_analysis\n", compartment_outdir))
cat(sprintf("║    RELAXED (0.25): %s/relaxed/peak2gene_stability_analysis\n", compartment_outdir))
cat("║                                                                  ║\n")
cat("║  Files per cutoff:                                               ║\n")
cat("║    • genes_[sex]_[stability].txt                                 ║\n")
cat("║    • peaks_[sex]_[stability].tsv                                 ║\n")
cat("║    • bins_[sex]_[stability].rds/bed                              ║\n")
cat("║    • genes_by_sex_class.rds                                      ║\n")
cat("║    • overlap_stats.rds                                           ║\n")
cat("║                                                                  ║\n")
cat("╚══════════════════════════════════════════════════════════════════╝\n")

message("\n[DONE] Both default and relaxed P2G analyses complete!")



