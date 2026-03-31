#!/usr/bin/env Rscript
## =========================================================
## DOWNLOAD REPEATMASKER ANNOTATIONS FOR mm10
## 
## Downloads RepeatMasker data from UCSC or AnnotationHub
## and prepares it for stability class analysis
## =========================================================

cat("\n")
cat("╔══════════════════════════════════════════════════════════════════╗\n")
cat("║       DOWNLOAD REPEATMASKER ANNOTATIONS - mm10                   ║\n")
cat("╚══════════════════════════════════════════════════════════════════╝\n\n")

## =========================================================
## CONFIGURATION
## =========================================================

# Output directory
output_dir <- "/data/sarkern2/multiome_liver/Seurat/epigenome/repeatmasker"

if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    message(sprintf("Created output directory: %s", output_dir))
}

## =========================================================
## LOAD PACKAGES
## =========================================================

message("[SETUP] Loading packages...")

suppressPackageStartupMessages({
    library(AnnotationHub)
    library(GenomicRanges)
    library(rtracklayer)
    library(data.table)
})

message("  ✅ Packages loaded")


## =========================================================
## METHOD 1: DOWNLOAD FROM ANNOTATIONHUB
## =========================================================

message("\n[METHOD 1] Trying AnnotationHub...")

tryCatch({
    
    ah <- AnnotationHub()
    
    # Query for mm10 RepeatMasker
    rmsk_query <- query(ah, c("RepeatMasker", "mm10", "UCSC"))
    
    message(sprintf("  Found %d RepeatMasker entries:", length(rmsk_query)))
    print(rmsk_query)
    
    if (length(rmsk_query) > 0) {
        
        message("\n  Downloading RepeatMasker from AnnotationHub...")
        message("  (This may take a few minutes...)")
        
        # Download the first matching entry
        rmsk_gr <- ah[[names(rmsk_query)[1]]]
        
        message(sprintf("  ✅ Downloaded: %d repeat elements", length(rmsk_gr)))
        
        # Check metadata columns
        message("\n  Metadata columns:")
        print(head(mcols(rmsk_gr)))
        
        # Save as RDS
        saveRDS(rmsk_gr, file.path(output_dir, "mm10_repeatmasker_AnnotationHub.rds"))
        message(sprintf("  ✅ Saved: mm10_repeatmasker_AnnotationHub.rds"))
        
        rmsk_source <- "AnnotationHub"
        
    } else {
        message("  ⚠️ No RepeatMasker found in AnnotationHub, trying UCSC...")
        rmsk_source <- NULL
    }
    
}, error = function(e) {
    message(sprintf("  ❌ AnnotationHub failed: %s", e$message))
    rmsk_source <- NULL
})


## =========================================================
## METHOD 2: DOWNLOAD DIRECTLY FROM UCSC
## =========================================================

if (!exists("rmsk_gr") || is.null(rmsk_gr)) {
    
    message("\n[METHOD 2] Downloading from UCSC Table Browser...")
    
    # UCSC RepeatMasker URL
    ucsc_url <- "https://hgdownload.soe.ucsc.edu/goldenPath/mm10/database/rmsk.txt.gz"
    
    local_file <- file.path(output_dir, "rmsk.txt.gz")
    
    tryCatch({
        
        message(sprintf("  Downloading: %s", ucsc_url))
        download.file(ucsc_url, local_file, mode = "wb", quiet = FALSE)
        
        message("  Parsing RepeatMasker file...")
        
        # Read the file (UCSC format)
        # Columns: bin, swScore, milliDiv, milliDel, milliIns, genoName, genoStart, genoEnd, 
        #          genoLeft, strand, repName, repClass, repFamily, repStart, repEnd, repLeft, id
        
        rmsk_dt <- fread(
            cmd = sprintf("zcat %s", local_file),
            header = FALSE,
            col.names = c("bin", "swScore", "milliDiv", "milliDel", "milliIns",
                          "genoName", "genoStart", "genoEnd", "genoLeft", "strand",
                          "repName", "repClass", "repFamily", "repStart", "repEnd", 
                          "repLeft", "id")
        )
        
        message(sprintf("  ✅ Loaded: %d repeat elements", nrow(rmsk_dt)))
        
        # Convert to GRanges
        rmsk_gr <- GRanges(
            seqnames = rmsk_dt$genoName,
            ranges = IRanges(start = rmsk_dt$genoStart + 1, end = rmsk_dt$genoEnd),
            strand = rmsk_dt$strand,
            repName = rmsk_dt$repName,
            repClass = rmsk_dt$repClass,
            repFamily = rmsk_dt$repFamily,
            swScore = rmsk_dt$swScore
        )
        
        # Keep standard chromosomes only
        rmsk_gr <- keepStandardChromosomes(rmsk_gr, pruning.mode = "coarse")
        
        message(sprintf("  ✅ GRanges created: %d elements (standard chromosomes)", length(rmsk_gr)))
        
        # Save
        saveRDS(rmsk_gr, file.path(output_dir, "mm10_repeatmasker_UCSC.rds"))
        fwrite(rmsk_dt, file.path(output_dir, "mm10_repeatmasker_UCSC.tsv.gz"), sep = "\t")
        
        message("  ✅ Saved: mm10_repeatmasker_UCSC.rds")
        
        rmsk_source <- "UCSC"
        
    }, error = function(e) {
        message(sprintf("  ❌ UCSC download failed: %s", e$message))
    })
}


## =========================================================
## SUMMARIZE REPEAT CLASSES
## =========================================================

if (exists("rmsk_gr") && !is.null(rmsk_gr)) {
    
    message("\n[SUMMARY] Repeat element statistics...")
    
    # Get repeat class distribution
    class_counts <- as.data.table(table(mcols(rmsk_gr)$repClass))
    setnames(class_counts, c("repClass", "N"))
    class_counts <- class_counts[order(-N)]
    
    message("\n  Top repeat classes:")
    print(head(class_counts, 15))
    
    # Get repeat family distribution
    family_counts <- as.data.table(table(mcols(rmsk_gr)$repFamily))
    setnames(family_counts, c("repFamily", "N"))
    family_counts <- family_counts[order(-N)]
    
    message("\n  Top repeat families:")
    print(head(family_counts, 15))
    
    # Save summaries
    fwrite(class_counts, file.path(output_dir, "repeat_class_counts.tsv"), sep = "\t")
    fwrite(family_counts, file.path(output_dir, "repeat_family_counts.tsv"), sep = "\t")
    
    
    ## =========================================================
    ## CREATE SIMPLIFIED REPEAT CATEGORIES
    ## =========================================================
    
    message("\n[SIMPLIFY] Creating major repeat categories...")
    
    # Define major categories
    rmsk_dt_summary <- data.table(
        repClass = as.character(mcols(rmsk_gr)$repClass),
        repFamily = as.character(mcols(rmsk_gr)$repFamily)
    )
    
    # Create major categories
    rmsk_dt_summary[, Major_Class := fcase(
        repClass == "LINE", "LINE",
        repClass == "SINE", "SINE",
        repClass == "LTR", "LTR",
        repClass == "DNA", "DNA_Transposon",
        repClass == "Satellite", "Satellite",
        repClass == "Simple_repeat", "Simple_repeat",
        repClass == "Low_complexity", "Low_complexity",
        repClass == "RC", "Rolling_Circle",
        repClass == "RNA", "RNA",
        repClass == "srpRNA", "RNA",
        repClass == "scRNA", "RNA",
        repClass == "snRNA", "RNA",
        repClass == "tRNA", "RNA",
        repClass == "rRNA", "RNA",
        default = "Other"
    )]
    
    # Add to GRanges
    mcols(rmsk_gr)$Major_Class <- rmsk_dt_summary$Major_Class
    
    # Summary by major class
    major_counts <- rmsk_dt_summary[, .N, by = Major_Class][order(-N)]
    
    message("\n  Major repeat categories:")
    print(major_counts)
    
    fwrite(major_counts, file.path(output_dir, "repeat_major_class_counts.tsv"), sep = "\t")
    
    
    ## =========================================================
    ## SAVE FINAL PROCESSED DATA
    ## =========================================================
    
    message("\n[SAVE] Saving processed RepeatMasker data...")
    
    # Save final GRanges with Major_Class
    saveRDS(rmsk_gr, file.path(output_dir, "mm10_repeatmasker_processed.rds"))
    
    # Also save as BED for external tools
    rmsk_bed <- data.table(
        chr = as.character(seqnames(rmsk_gr)),
        start = start(rmsk_gr) - 1,  # BED is 0-based
        end = end(rmsk_gr),
        name = mcols(rmsk_gr)$repName,
        score = mcols(rmsk_gr)$swScore,
        strand = as.character(strand(rmsk_gr)),
        repClass = mcols(rmsk_gr)$repClass,
        repFamily = mcols(rmsk_gr)$repFamily,
        Major_Class = mcols(rmsk_gr)$Major_Class
    )
    
    fwrite(rmsk_bed, file.path(output_dir, "mm10_repeatmasker.bed.gz"), sep = "\t", col.names = FALSE)
    
    message("  ✅ Saved: mm10_repeatmasker_processed.rds")
    message("  ✅ Saved: mm10_repeatmasker.bed.gz")
    
}


## =========================================================
## FINAL SUMMARY
## =========================================================

cat("\n")
cat("╔══════════════════════════════════════════════════════════════════╗\n")
cat("║                    DOWNLOAD COMPLETE                             ║\n")
cat("╠══════════════════════════════════════════════════════════════════╣\n")
cat(sprintf("║  Source: %-55s║\n", rmsk_source))
cat(sprintf("║  Total elements: %-47s║\n", format(length(rmsk_gr), big.mark = ",")))
cat(sprintf("║  Output: %-55s║\n", substr(output_dir, 1, 55)))
cat("║                                                                  ║\n")
cat("║  Files created:                                                  ║\n")
cat("║    • mm10_repeatmasker_processed.rds (GRanges with Major_Class)  ║\n")
cat("║    • mm10_repeatmasker.bed.gz (BED format)                       ║\n")
cat("║    • repeat_class_counts.tsv                                     ║\n")
cat("║    • repeat_family_counts.tsv                                    ║\n")
cat("║    • repeat_major_class_counts.tsv                               ║\n")
cat("║                                                                  ║\n")
cat("║  Major categories for analysis:                                  ║\n")
cat("║    • LINE (heterochromatin-associated)                           ║\n")
cat("║    • SINE (euchromatin-associated)                               ║\n")
cat("║    • LTR (ERVs - age-related)                                    ║\n")
cat("║    • DNA_Transposon                                              ║\n")
cat("║    • Satellite (constitutive heterochromatin)                    ║\n")
cat("║    • Simple_repeat, Low_complexity, Other                        ║\n")
cat("╚══════════════════════════════════════════════════════════════════╝\n")

message("\n[DONE] RepeatMasker ready for stability class analysis!")
message("\nNext step: Run the mapping script to analyze repeat enrichment per stability class")


