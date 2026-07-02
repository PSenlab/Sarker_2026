#!/usr/bin/env Rscript
#===============================================================================
# DEFINE AGE-DYNAMIC PEAKS  (gaining / losing / stable)
#
# Standalone classification step: pseudobulk the ATAC PeakMatrix by sample,
# fit a sex-adjusted limma-voom age trend, and label every peak by the sign of
# its (FDR-significant) age slope.
#
#   gaining = FDR < DA_FDR  &  age slope > +DA_LFC   (opens with age)
#   losing  = FDR < DA_FDR  &  age slope < -DA_LFC   (closes with age)
#   stable  = everything else
#
# Output: a table of peak_id, logFC (per-age-step slope), adj.P.Val, progressive.
# This is the exact classification the SCAVENGE age-seeding scripts use to build
# TARGET_PEAK_IDS; running it standalone lets you inspect the peaks first.
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
# PARAMETERS
#===============================================================================
ARCHR_PROJECT <- "/data/sarkern2/multiome_liver/Seurat/archR/ArchR_Projects/Step8_Hepatocyte_CCAN_P2G"
OUT_DIR       <- "/data/sarkern2/multiome_liver/Seurat/epigenome/gwas_peak_gene/age_dynamic_peaks"

# Age ordering (young -> geriatric) coded as an ordinal 1..5 rank
age_levels <- c("young", "mid_age", "old", "pre_geriatric", "geriatric")
SAMPLE_COL <- "Sample"   # colData column holding the biological replicate id
AGE_COL    <- "age"      # colData column holding the age stage
SEX_COL    <- "sex"      # colData column holding sex (adjusted for)

DA_FDR <- 0.05           # FDR threshold to call a peak age-dynamic
DA_LFC <- 0              # min |age slope| (0 = classify on direction only)

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
log_file <- file.path(OUT_DIR, paste0("log_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt"))
write_log <- function(msg) {
  cat(msg, "\n")
  cat(paste0("[", Sys.time(), "] ", msg, "\n"), file = log_file, append = TRUE)
}

write_log("===============================================================================")
write_log("DEFINE AGE-DYNAMIC PEAKS (gaining / losing / stable)")
write_log("===============================================================================")
write_log(paste("Project:", ARCHR_PROJECT))
write_log(paste("Age levels:", paste(age_levels, collapse = " < ")))
write_log(paste("DA_FDR:", DA_FDR, "| DA_LFC:", DA_LFC))

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

n_gain <- sum(age_res$progressive == "gaining")
n_lose <- sum(age_res$progressive == "losing")
n_stab <- sum(age_res$progressive == "stable")
write_log(paste("  gaining:", n_gain))
write_log(paste("  losing :", n_lose))
write_log(paste("  stable :", n_stab))
write_log(paste("  total  :", nrow(age_res)))

#===============================================================================
# STEP 5: Save
#===============================================================================
write_log("\n====== STEP 5: Save classification ======")

fwrite(age_res, file.path(OUT_DIR, "age_dynamic_peaks_classification.csv"))
saveRDS(list(res = age_res, params = list(age_levels = age_levels,
             DA_FDR = DA_FDR, DA_LFC = DA_LFC)),
        file.path(OUT_DIR, "age_da_results.rds"))

# convenience: bare peak-id lists per direction
writeLines(age_res$peak_id[age_res$progressive == "gaining"],
           file.path(OUT_DIR, "peaks_gaining.txt"))
writeLines(age_res$peak_id[age_res$progressive == "losing"],
           file.path(OUT_DIR, "peaks_losing.txt"))

write_log(paste("Saved:", file.path(OUT_DIR, "age_dynamic_peaks_classification.csv")))
write_log(paste("Saved:", file.path(OUT_DIR, "age_da_results.rds")))
write_log(paste("Saved:", file.path(OUT_DIR, "peaks_gaining.txt"), "and peaks_losing.txt"))

write_log("\n===============================================================================")
write_log("DONE")
write_log(paste("gaining:", n_gain, "| losing:", n_lose, "| stable:", n_stab))
write_log("===============================================================================")
write_log(paste("Finished:", Sys.time()))
