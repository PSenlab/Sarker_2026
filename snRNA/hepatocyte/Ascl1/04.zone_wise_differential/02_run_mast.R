#!/usr/bin/env Rscript
# ============================================================
# 02_run_mast.R  <subset>
# CDR-controlled MAST DE: Ascl1_pos vs Ascl1_neg (female hepatocytes).
# One script for global and every zone - no copy-paste.
#
#   subset in: allfemhep | periportal | midlobular | pericentral
#
# Run (serial; give RAM headroom - global ~90k cells needs the most):
#   sinteractive --mem=128g --cpus-per-task=4     # allfemhep
#   sinteractive --mem=96g  --cpus-per-task=4     # a single zone
#   Rscript 02_run_mast.R allfemhep
#   Rscript 02_run_mast.R periportal
# ============================================================

suppressPackageStartupMessages({
  library(MAST)
  library(Matrix)
  library(data.table)
  library(BiocParallel)
})

# force serial execution: avoids fork-death + matrix-copy memory blowup
options(mc.cores = 1)
register(SerialParam())

# ---- argument ----
args   <- commandArgs(trailingOnly = TRUE)
SUBSET <- if (length(args) >= 1) args[1] else "allfemhep"
VALID  <- c("allfemhep", "periportal", "midlobular", "pericentral")
if (!SUBSET %in% VALID)
  stop("subset must be one of: ", paste(VALID, collapse = ", "))

EXPORT_DIR <- "mast_export"
OUTDIR     <- "ascl1_de_by_zone_MAST"
FDR_CUTOFF <- 0.05
dir.create(OUTDIR, showWarnings = FALSE)

cat(strrep("=", 60), "\n")
cat("SUBSET:", SUBSET, "- MAST Ascl1_pos vs Ascl1_neg (CDR-controlled)\n")
cat(strrep("=", 60), "\n")

# ---- load ----
expr  <- readMM(file.path(EXPORT_DIR, paste0(SUBSET, "_expr.mtx")))   # genes x cells
genes <- fread(file.path(EXPORT_DIR, paste0(SUBSET, "_genes.csv")))$gene
meta  <- fread(file.path(EXPORT_DIR, paste0(SUBSET, "_meta.csv")))

expr <- as(expr, "CsparseMatrix")
rownames(expr) <- genes
colnames(expr) <- meta$barcode

# ---- MAST object on log-normalized matrix ----
sca <- FromMatrix(as.matrix(expr),
                  cData = data.frame(meta),
                  fData = data.frame(primerid = genes))
colData(sca)$Ascl1_status <- factor(colData(sca)$Ascl1_status,
                                    levels = c("Ascl1_neg", "Ascl1_pos"))

# ---- depth control: scaled cellular detection rate ----
cdr <- colSums(assay(sca) > 0)
colData(sca)$cngeneson <- scale(cdr)[, 1]

npos <- sum(colData(sca)$Ascl1_status == "Ascl1_pos")
nneg <- sum(colData(sca)$Ascl1_status == "Ascl1_neg")
cat(sprintf("   Cells: %d | Ascl1+: %d  Ascl1-: %d\n", ncol(sca), npos, nneg))

# ---- hurdle model with CDR covariate ----
cat("   Fitting MAST hurdle model (slow step)...\n")
zfit <- zlm(~ Ascl1_status + cngeneson, sca, parallel = FALSE)

coefName <- "Ascl1_statusAscl1_pos"
summ <- summary(zfit, doLRT = coefName)
dt   <- summ$datatable

res <- merge(
  dt[contrast == coefName & component == "H",
     .(primerid, p_value = `Pr(>Chisq)`)],
  dt[contrast == coefName & component == "logFC",
     .(primerid, logFC = coef, ci_lo = ci.lo, ci_hi = ci.hi)],
  by = "primerid"
)
res[, fdr := p.adjust(p_value, "BH")]
setorder(res, fdr, -logFC)

sig  <- res[fdr < FDR_CUTOFF & !is.na(logFC)]
up   <- sig[logFC > 0][order(-logFC)]
down <- sig[logFC < 0][order(logFC)]

# output filename: allfemhep -> female_hepatocytes; zones keep their name
tag <- if (SUBSET == "allfemhep") "female_hepatocytes" else SUBSET
fwrite(res,  file.path(OUTDIR, sprintf("Ascl1_MAST_%s_all.csv",         tag)))
fwrite(sig,  file.path(OUTDIR, sprintf("Ascl1_MAST_%s_significant.csv", tag)))
fwrite(up,   file.path(OUTDIR, sprintf("Ascl1_MAST_%s_up.csv",          tag)))
fwrite(down, file.path(OUTDIR, sprintf("Ascl1_MAST_%s_down.csv",        tag)))

cat(sprintf("   Significant: %d (up %d / down %d)\n", nrow(sig), nrow(up), nrow(down)))
cat(sprintf("   Saved -> %s/Ascl1_MAST_%s_*.csv\n", OUTDIR, tag))
cat("Done:", SUBSET, "\n")
