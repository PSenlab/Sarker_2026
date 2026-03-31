#!/usr/bin/env Rscript
## =========================================================
## LOCAL GREAT (REACTOME ONLY) – RELAXED P2G (corCutOff = 0.25)
## =========================================================
suppressPackageStartupMessages({
  library(rGREAT)
  library(reactome.db)
  library(GenomicRanges)
  library(IRanges)
  library(S4Vectors)
  library(TxDb.Mmusculus.UCSC.mm10.knownGene)
  library(data.table)
})

## =========================================================
## PATHS
## =========================================================
base_outdir <- "/data/sarkern2/multiome_liver/Seurat/epigenome/downstream_stability"
p2g_outdir  <- file.path(base_outdir, "relaxed", "peak2gene_stability_analysis")
great_outdir <- file.path(p2g_outdir, "GREAT_analysis")
dir.create(great_outdir, recursive = TRUE, showWarnings = FALSE)

message("\n")
message("╔══════════════════════════════════════════════════════════════════╗")
message("║  LOCAL GREAT ANALYSIS: RELAXED P2G (corCutOff = 0.25)            ║")
message("╚══════════════════════════════════════════════════════════════════╝\n")

message(sprintf("  Input: %s", p2g_outdir))
message(sprintf("  Output: %s", great_outdir))

## =========================================================
## LOAD INPUT OBJECTS (ROBUST)
## =========================================================
message("\n>>> Loading input data...")

peak_gr <- readRDS(
  file.path(p2g_outdir, "hepatocyte_p2g_peaks_gr.rds")
)
stability_gr_list <- readRDS(
  file.path(p2g_outdir, "stability_gr_list.rds")
)

SEXES <- names(stability_gr_list)
STABILITY_CLASSES <- unique(unlist(lapply(stability_gr_list, names)))

message(sprintf("  ✓ Loaded %d peaks", length(peak_gr)))
message(sprintf("  ✓ SEXES: %s", paste(SEXES, collapse=", ")))
message(sprintf("  ✓ STABILITY_CLASSES: %s", paste(STABILITY_CLASSES, collapse=", ")))

## =========================================================
## BUILD REACTOME GENE SETS (MOUSE)
## =========================================================
message("\n>>> Building Reactome gene sets for mouse...")

reactome_all <- as.list(reactomePATHID2EXTID)
reactome_mouse <- reactome_all[grep("^R-MMU-", names(reactome_all))]
pathway_names <- as.list(reactomePATHID2NAME)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

message(sprintf("  ✓ Mouse Reactome pathways: %d", length(reactome_mouse)))

## =========================================================
## FUNCTION: RUN LOCAL GREAT (REACTOME)
## =========================================================
run_great_reactome <- function(gr, name) {
  if (length(gr) < 5) {
    message(sprintf("    %s: too few regions (%d), skip", name, length(gr)))
    return(NULL)
  }
  
  message(sprintf("    Running GREAT: %s (%d regions)", name, length(gr)))
  
  job <- great(
    gr,
    gene_sets = reactome_mouse,
    tss_source = "txdb:mm10",
    biomart_dataset = NULL,
    min_gene_set_size = 2,
    mode = "basalPlusExt",
    basal_upstream = 5000,
    basal_downstream = 1000,
    extension = 1000000
  )
  
  tb <- getEnrichmentTable(job)
  if (nrow(tb) == 0) return(NULL)
  
  dt <- as.data.table(tb)
  dt[, Analysis := name]
  dt[, Pathway_Name := pathway_names[id]]
  dt[, neglog10_pval := -log10(p_adjust + 1e-300)]
  
  fwrite(dt, file.path(great_outdir, sprintf("GREAT_%s_Reactome.tsv", name)), sep="\t")
  
  # Count significant
  n_sig <- sum(dt$p_adjust < 0.05, na.rm = TRUE)
  message(sprintf("      ✓ %d total terms, %d significant (p < 0.05)", nrow(dt), n_sig))
  
  dt
}

## =========================================================
## MAIN LOOP
## =========================================================
message("\n")
message("═══════════════════════════════════════════════════════════════════")
message("  RUNNING LOCAL GREAT (REACTOME) ON RELAXED P2G PEAKS")
message("═══════════════════════════════════════════════════════════════════")

great_results <- list()

for (sx in SEXES) {
  great_results[[sx]] <- list()
  
  message(sprintf("\n>>> %s:", toupper(sx)))
  message(paste(rep("-", 50), collapse = ""))
  
  for (stab in STABILITY_CLASSES) {
    bins_gr <- stability_gr_list[[sx]][[stab]]
    
    if (is.null(bins_gr) || length(bins_gr) == 0) {
      message(sprintf("    %s: No bins, skipping", stab))
      great_results[[sx]][[stab]] <- NULL
      next
    }
    
    olaps <- findOverlaps(peak_gr, bins_gr, ignore.strand = TRUE)
    
    if (length(olaps) == 0) {
      message(sprintf("    %s: %d bins → 0 peaks, skipping", stab, length(bins_gr)))
      great_results[[sx]][[stab]] <- NULL
      next
    }
    
    peak_hits <- unique(S4Vectors::queryHits(olaps))
    peaks_use <- peak_gr[peak_hits]
    
    message(sprintf("\n    %s: %d bins → %d peaks", stab, length(bins_gr), length(peaks_use)))
    
    # Save BED
    bed_file <- file.path(great_outdir, sprintf("peaks_%s_%s.bed", sx, stab))
    fwrite(
      data.table(
        chr = as.character(seqnames(peaks_use)),
        start = start(peaks_use) - 1L,
        end = end(peaks_use),
        name = peaks_use$peak_id
      ),
      bed_file,
      sep = "\t",
      col.names = FALSE
    )
    message(sprintf("      Saved: %s", basename(bed_file)))
    
    nm <- sprintf("%s_%s", sx, stab)
    great_results[[sx]][[stab]] <- run_great_reactome(peaks_use, nm)
  }
}

## =========================================================
## SAVE COMBINED RESULTS
## =========================================================
message("\n>>> Saving R objects...")

saveRDS(great_results, file.path(great_outdir, "great_results_reactome_relaxed.rds"))
message("  ✓ great_results_reactome_relaxed.rds")

## =========================================================
## CREATE SUMMARY TABLE
## =========================================================
message("\n>>> Creating summary table...")

summary_rows <- list()

for (sx in SEXES) {
  for (stab in STABILITY_CLASSES) {
    res <- great_results[[sx]][[stab]]
    
    n_total <- if (!is.null(res)) nrow(res) else 0
    n_sig <- if (!is.null(res)) sum(res$p_adjust < 0.05, na.rm = TRUE) else 0
    
    summary_rows[[paste(sx, stab, sep = "_")]] <- data.table(
      Sex = sx,
      Stability = stab,
      Reactome_Total = n_total,
      Reactome_Sig = n_sig
    )
  }
}

summary_dt <- rbindlist(summary_rows)
summary_file <- file.path(great_outdir, "GREAT_summary_reactome_relaxed.tsv")
fwrite(summary_dt, summary_file, sep = "\t")

message(sprintf("\n  ✓ %s", basename(summary_file)))
print(summary_dt)

## =========================================================
## COMBINE ALL SIGNIFICANT RESULTS
## =========================================================
message("\n>>> Combining significant results...")

all_sig_results <- list()

for (sx in SEXES) {
  for (stab in STABILITY_CLASSES) {
    res <- great_results[[sx]][[stab]]
    
    if (!is.null(res) && nrow(res) > 0) {
      sig_res <- res[p_adjust < 0.05]
      if (nrow(sig_res) > 0) {
        sig_res[, Sex := sx]
        sig_res[, Stability := stab]
        all_sig_results[[paste(sx, stab, sep = "_")]] <- sig_res
      }
    }
  }
}

if (length(all_sig_results) > 0) {
  combined_sig <- rbindlist(all_sig_results, fill = TRUE)
  combined_file <- file.path(great_outdir, "GREAT_Reactome_significant_combined_relaxed.tsv")
  fwrite(combined_sig, combined_file, sep = "\t")
  message(sprintf("  ✓ Combined significant results: %d terms → %s", 
                  nrow(combined_sig), basename(combined_file)))
} else {
  message("  ⚠ No significant results to combine")
}

## =========================================================
## FINAL SUMMARY
## =========================================================
cat("\n")
cat("╔══════════════════════════════════════════════════════════════════╗\n")
cat("║        LOCAL GREAT ANALYSIS COMPLETE (RELAXED P2G)               ║\n")
cat("╠══════════════════════════════════════════════════════════════════╣\n")
cat("║                                                                  ║\n")
cat("║  Version: RELAXED (corCutOff = 0.25)                             ║\n")
cat("║                                                                  ║\n")
cat(sprintf("║  Mouse Reactome pathways: %d                                ║\n", length(reactome_mouse)))
cat("║                                                                  ║\n")
cat(sprintf("║  Output: %s\n", great_outdir))
cat("║                                                                  ║\n")
cat("║  Files generated:                                                ║\n")
cat("║    • peaks_[sex]_[stability].bed                                 ║\n")
cat("║    • GREAT_[sex]_[stability]_Reactome.tsv                        ║\n")
cat("║    • GREAT_summary_reactome_relaxed.tsv                          ║\n")
cat("║    • GREAT_Reactome_significant_combined_relaxed.tsv             ║\n")
cat("║    • great_results_reactome_relaxed.rds                          ║\n")
cat("║                                                                  ║\n")
cat("╚══════════════════════════════════════════════════════════════════╝\n")

message("\n[DONE] Local GREAT Reactome analysis complete for RELAXED P2G!")


