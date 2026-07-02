#!/usr/bin/env Rscript
#===============================================================================
# DEFINE AGE-DYNAMIC PEAKS  (gaining / losing / stable) + P2G ANNOTATION
#-------------------------------------------------------------------------------
# Classifies every ATAC peak by how its accessibility changes with age, and
# flags which peaks are Peak-to-Gene (P2G) linked. This is the classification
# used to seed the age-directional GWAS-SCAVENGE analysis.
#
# Method:
#   1. Pull the ArchR PeakMatrix (raw counts) and pseudobulk by biological
#      replicate (sample).
#   2. Fit a sex-adjusted limma-voom trend against an ordinal age rank
#      (young=1 ... geriatric=5); the age_rank coefficient is the per-stage
#      change in accessibility.
#   3. Classify each peak by the sign of its FDR-significant age slope:
#        gaining = FDR < DA_FDR  &  slope > +DA_LFC   (opens with age)
#        losing  = FDR < DA_FDR  &  slope < -DA_LFC   (closes with age)
#        stable  = otherwise
#   4. Annotate which peaks are P2G-linked (getPeak2GeneLinks). The subset used
#      to seed SCAVENGE is (progressive != "stable") & is_p2g.
#
# Outputs (in OUT_DIR):
#   age_dynamic_peaks_classification.csv  full table (all peaks)
#   age_da_results.rds                    list(res, params, counts)
#   peaks_gaining.txt / peaks_losing.txt              genome-wide peak ids
#   peaks_gaining_p2g.txt / peaks_losing_p2g.txt      P2G-linked peak ids
#
# Usage:
#   Rscript Define_AgeDynamic_Peaks.R
#===============================================================================

suppressPackageStartupMessages({
  library(ArchR)
  library(GenomicRanges)
  library(GenomeInfoDb)
  library(Matrix)
  library(SummarizedExperiment)
  library(edgeR)
  library(limma)
  library(data.table)
})

addArchRThreads(60)
addArchRGenome("mm10")
set.seed(42)

#===============================================================================
# PARAMETERS  (edit these)
#===============================================================================
ARCHR_PROJECT <- "/data/sarkern2/multiome_liver/Seurat/archR/ArchR_Projects/Step8_Hepatocyte_CCAN_P2G"
OUT_DIR       <- "/data/sarkern2/multiome_liver/Seurat/epigenome/gwas_peak_gene/age_dynamic_peaks_02"

# Age ordering (young -> geriatric), coded as an ordinal 1..5 rank
age_levels <- c("young", "mid_age", "old", "pre_geriatric", "geriatric")
SAMPLE_COL <- "Sample"   # colData column: biological replicate id
AGE_COL    <- "age"      # colData column: age stage
SEX_COL    <- "sex"      # colData column: sex (adjusted for in the model)

DA_FDR <- 0.05           # FDR threshold to call a peak age-dynamic
DA_LFC <- 0              # min |age slope| (0 = classify on direction only)

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
log_file <- file.path(OUT_DIR, paste0("log_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt"))
write_log <- function(msg) {
  cat(msg, "\n")
  cat(paste0("[", Sys.time(), "] ", msg, "\n"), file = log_file, append = TRUE)
}

write_log("===============================================================================")
write_log("DEFINE AGE-DYNAMIC PEAKS (gaining / losing / stable) + P2G annotation")
write_log("===============================================================================")
write_log(paste("Project   :", ARCHR_PROJECT))
write_log(paste("Age levels:", paste(age_levels, collapse = " < ")))
write_log(paste("DA_FDR    :", DA_FDR, "| DA_LFC:", DA_LFC))

#===============================================================================
# STEP 1: Load ArchR project
#===============================================================================
write_log("\n====== STEP 1: Load ArchR Project ======")
proj <- tryCatch(
  loadArchRProject(ARCHR_PROJECT),
  error = function(e) { write_log(paste("CRITICAL ERROR:", e$message)); stop("Cannot load project") })
write_log(paste("Loaded:", length(getCellNames(proj)), "cells"))

#===============================================================================
# STEP 2: Pull PeakMatrix (RAW counts) and pseudobulk by sample
#===============================================================================
write_log("\n====== STEP 2: PeakMatrix -> pseudobulk by sample ======")

peakMat <- getMatrixFromProject(ArchRProj = proj, useMatrix = "PeakMatrix", binarize = FALSE)
write_log(paste("Full PeakMatrix:", nrow(peakMat), "peaks x", ncol(peakMat), "cells"))
pk_all <- rowRanges(peakMat)

cdat <- as.data.frame(colData(peakMat))
samp <- factor(cdat[[SAMPLE_COL]])
pb   <- as.matrix(assay(peakMat) %*% t(Matrix::fac2sparse(samp)))   # peaks x samples
colnames(pb) <- levels(samp)

meta <- cdat[!duplicated(samp), c(SAMPLE_COL, SEX_COL, AGE_COL)]
meta <- meta[match(colnames(pb), meta[[SAMPLE_COL]]), ]
meta$age_rank <- as.integer(factor(meta[[AGE_COL]], levels = age_levels, ordered = TRUE))
meta$sex      <- factor(meta[[SEX_COL]])
stopifnot(identical(colnames(pb), meta[[SAMPLE_COL]]))
write_log(paste("Pseudobulk:", nrow(pb), "peaks x", ncol(pb), "samples"))
write_log(paste("Samples per age stage:",
                paste(names(table(meta[[AGE_COL]])), table(meta[[AGE_COL]]),
                      sep = "=", collapse = ", ")))

#===============================================================================
# STEP 3: Sex-adjusted limma-voom age trend
#===============================================================================
write_log("\n====== STEP 3: limma-voom age trend (~ sex + age_rank) ======")

design <- model.matrix(~ sex + age_rank, data = meta)
dge  <- DGEList(pb)
keep <- filterByExpr(dge, design)
dge  <- dge[keep, , keep.lib.sizes = FALSE]
pkf  <- pk_all[keep]
write_log(paste("Peaks after filterByExpr:", sum(keep), "of", length(keep)))

dge <- calcNormFactors(dge)
v   <- voom(dge, design)
fit <- eBayes(lmFit(v, design))
tt  <- topTable(fit, coef = "age_rank", number = Inf, sort.by = "none")

#===============================================================================
# STEP 4: Classify peaks by sign of (FDR-significant) age slope
#===============================================================================
write_log("\n====== STEP 4: Classify gaining / losing / stable ======")

age_res <- data.frame(
  peak_id     = paste0(seqnames(pkf), ":", start(pkf), "-", end(pkf)),
  seqnames    = as.character(seqnames(pkf)),
  start       = start(pkf),
  end         = end(pkf),
  logFC       = tt$logFC,        # per-age-step change in accessibility
  AveExpr     = tt$AveExpr,
  t           = tt$t,
  P.Value     = tt$P.Value,
  adj.P.Val   = tt$adj.P.Val,
  progressive = ifelse(tt$adj.P.Val < DA_FDR & tt$logFC >  DA_LFC, "gaining",
                ifelse(tt$adj.P.Val < DA_FDR & tt$logFC < -DA_LFC, "losing", "stable")),
  stringsAsFactors = FALSE)

write_log(paste("  gaining:", sum(age_res$progressive == "gaining")))
write_log(paste("  losing :", sum(age_res$progressive == "losing")))
write_log(paste("  stable :", sum(age_res$progressive == "stable")))
write_log(paste("  total  :", nrow(age_res)))

#===============================================================================
# STEP 5: Annotate P2G-linked peaks
#   The SCAVENGE seeding uses only peaks that are BOTH age-dynamic AND P2G.
#===============================================================================
write_log("\n====== STEP 5: Annotate P2G-linked peaks ======")

p2g <- getPeak2GeneLinks(ArchRProj = proj, returnLoops = FALSE)
p2g_peak_idx <- unique(p2g$idxATAC)
p2g_gr <- pk_all[p2g_peak_idx]                       # P2G peaks in full-peak coords
p2g_ids <- paste0(seqnames(p2g_gr), ":", start(p2g_gr), "-", end(p2g_gr))
write_log(paste("P2G links:", nrow(p2g), "| unique P2G peaks:", length(p2g_peak_idx)))

age_res$is_p2g <- age_res$peak_id %in% p2g_ids

g_all <- sum(age_res$progressive == "gaining")
l_all <- sum(age_res$progressive == "losing")
g_p2g <- sum(age_res$progressive == "gaining" & age_res$is_p2g)
l_p2g <- sum(age_res$progressive == "losing"  & age_res$is_p2g)
write_log(paste("  gaining -- genome-wide:", g_all, "| P2G-linked:", g_p2g))
write_log(paste("  losing  -- genome-wide:", l_all, "| P2G-linked:", l_p2g))

#===============================================================================
# STEP 6: Save
#===============================================================================
write_log("\n====== STEP 6: Save ======")

fwrite(age_res, file.path(OUT_DIR, "age_dynamic_peaks_classification.csv"))
saveRDS(list(res = age_res,
             params = list(age_levels = age_levels, DA_FDR = DA_FDR, DA_LFC = DA_LFC),
             counts = list(gaining_all = g_all, losing_all = l_all,
                           gaining_p2g = g_p2g, losing_p2g = l_p2g)),
        file.path(OUT_DIR, "age_da_results.rds"))

writeLines(age_res$peak_id[age_res$progressive == "gaining"],
           file.path(OUT_DIR, "peaks_gaining.txt"))
writeLines(age_res$peak_id[age_res$progressive == "losing"],
           file.path(OUT_DIR, "peaks_losing.txt"))
writeLines(age_res$peak_id[age_res$progressive == "gaining" & age_res$is_p2g],
           file.path(OUT_DIR, "peaks_gaining_p2g.txt"))
writeLines(age_res$peak_id[age_res$progressive == "losing"  & age_res$is_p2g],
           file.path(OUT_DIR, "peaks_losing_p2g.txt"))

write_log(paste("Saved:", file.path(OUT_DIR, "age_dynamic_peaks_classification.csv")))
write_log(paste("Saved:", file.path(OUT_DIR, "age_da_results.rds")))
write_log("Saved: peaks_gaining.txt, peaks_losing.txt, peaks_gaining_p2g.txt, peaks_losing_p2g.txt")

write_log("\n===============================================================================")
write_log("DONE")
write_log(paste("gaining (all/P2G):", g_all, "/", g_p2g,
                "| losing (all/P2G):", l_all, "/", l_p2g))
write_log("===============================================================================")
write_log(paste("Finished:", Sys.time()))
