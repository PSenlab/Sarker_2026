#!/usr/bin/env Rscript
#===============================================================================
# GWAS-SCAVENGE: LIVER TRAITS - P2G PEAKS, SEEDED FROM AGE-DYNAMIC PEAKS
#
# Identical to the all-P2G hepatocyte pipeline, EXCEPT the trait seeding is
# restricted to peaks classified as age-GAINING (or age-LOSING) by the limma
# age-trend analysis. The cell-cell network is still built on the full P2G
# peak matrix (so LSI/kNN are unaffected); only which peaks carry trait signal
# into computeWeightedDeviations is masked.
#
# Run TWICE:
#   DIRECTION <- "gaining"   ->  results land in .../hepatocyte_gaining
#   DIRECTION <- "losing"    ->  results land in .../hepatocyte_losing
# then contrast median TRS / %FDR by age between the two.
#
# Usage:
#   Rscript Run_GWAS_Liver_P2G_ageSeeded.R
#===============================================================================

suppressPackageStartupMessages({
  library(ArchR)
  library(GenomicRanges)
  library(GenomeInfoDb)
  library(dplyr)
  library(tidyr)
  library(Matrix)
  library(SummarizedExperiment)
  library(BSgenome.Mmusculus.UCSC.mm10)
  library(edgeR)
  library(limma)
  library(SCAVENGE)
  library(chromVAR)
  library(gchromVAR)
  library(parallel)
  library(data.table)
  library(rtracklayer)
})

addArchRThreads(60)
addArchRGenome("mm10")
set.seed(42)

#===============================================================================
# PARAMETERS
#===============================================================================
CATEGORY <- "Liver_P2G"
PEAK_EXTENSION <- 500
LD_WINDOW <- 25000
MIN_REGIONS <- 2
MIN_SEED_CELLS <- 5
MIN_VALID_CELLS <- 10
SCAVENGE_CORES <- 55
FDR_THRESHOLD <- 0.05

#-------------------------------------------------------------------------------
# AGE-SEEDING CONFIG  -- the only new block vs the all-P2G script
#-------------------------------------------------------------------------------
# Which age-dynamic class to seed from: "gaining" or "losing"
DIRECTION <- "losing"           # LOSING run

# Age-trend is computed INLINE from the ArchR project (no external file needed).
# Peaks are classified gaining / losing by a sample-level limma age trend,
# identical to the standalone age analysis (pseudobulk by Sample, sex-adjusted).
age_levels <- c("young", "mid_age", "old", "pre_geriatric", "geriatric")
SAMPLE_COL <- "Sample"; AGE_COL <- "age"; SEX_COL <- "sex"
DA_FDR <- 0.05      # FDR to call a peak age-dynamic
DA_LFC <- 0         # min |age slope| (0 = direction only)

#===============================================================================
# INCLUSION PATTERNS - Comprehensive liver biology
#===============================================================================
LIVER_PATTERNS <- c(
  "\\bliver\\b", "hepat",
  "NAFLD", "NASH", "steatosis", "cirrhosis",
  "\\bALT\\b", "\\bAST\\b", "\\bGGT\\b", "\\bALP\\b",
  "alanine aminotransferase", "aspartate aminotransferase",
  "gamma glutamyl transferase", "gamma glutamyltransferase",
  "gamma-glutamyl transpeptidase", "alkaline phosphatase", "bilirubin",
  "albumin", "fibrinogen", "prothrombin", "coagulation factor",
  "complement C", "haptoglobin", "transferrin", "ceruloplasmin",
  "alpha.1.antitrypsin", "angiotensinogen",
  "apolipoprotein", "\\bAPO[ABCDEFHLM]\\b", "\\bVLDL\\b",
  "ferritin", "hemochromatosis", "hepcidin", "\\bliver iron\\b",
  "hepatitis B", "hepatitis C", "\\bHBV\\b", "\\bHCV\\b",
  "autoimmune hepatitis", "primary biliary", "primary sclerosing",
  "\\bPBC\\b", "\\bPSC\\b", "drug.induced liver", "hepatotoxicity",
  "cholangiocarcinoma", "biliary", "bile duct", "gallbladder"
)

#===============================================================================
# EXCLUSION PATTERNS - Non-liver traits to remove
#===============================================================================
EXCLUSION_PATTERNS <- c(
  "putamen iron", "pallidum iron", "caudate iron", "substantia nigra iron",
  "thalamus iron", "hippocampus iron", "accumbens iron", "amygdala iron",
  "\\bR2\\* MRI", "quantitative susceptibility mapping",
  "environmental pollutants", "environmental tobacco", "environmental stress",
  "environmental factor", "environmentalism", "environmental exposures",
  "occupational environmental",
  "gamma-glutamylglutamate", "gamma-glutamylglutamine", "gamma-glutamylleucine",
  "gamma-glutamylvaline", "gamma-glutamyltyrosine", "gamma-glutamylphenylalanine",
  "gamma-glutamylthreonine", "gamma-glutamylisoleucine", "gamma-glutamylmethionine",
  "gamma-glutamylhistidine", "gamma-glutamylalanine", "gamma-glutamyltryptophan",
  "gamma-glutamylglycine", "gamma-glutamyl-alpha-lysine", "gamma-glutamyl-epsilon-lysine",
  "gamma-glutamyl-2-aminobutyrate", "gamma-glutamylcitrulline",
  "gamma-glutamylaminecyclotransferase", "gamma-glutamylcyclotransferase",
  "gamma-glutamyl hydrolase", "Glutathione-specific gamma-glutamylcyclotransferase",
  "Protein-glutamine gamma-glutamyltransferase", "Inactive gamma-glutamyltranspeptidase",
  "gamma-glutamyltransferase 5",
  "very large VLDL", "very small VLDL", "large VLDL", "medium VLDL", "small VLDL",
  "very very large VLDL", "chylomicrons and extremely large VLDL",
  "chylomicron and extremely large VLDL", "Fasting.*VLDL", "Postprandial.*VLDL",
  "VLDL particle concentration", "VLDL particles", "diameter.*VLDL",
  "to Total Lipids.*VLDL percentage", "to total lipids ratio.*VLDL", "UKB data field 23",
  "Alzheimer", "APOE e4", "APOEe4", "APOE E4", "dementia", "amyloid", "cerebral", "cerebrospinal",
  "apoptosis", "apoptotic", "APOBEC",
  "cystic fibrosis", "pulmonary fibrosis", "idiopathic pulmonary",
  "lung disease", "lung function", "FEV1", "FVC", "meconium ileus", "CFTR",
  "intestinal.type alkaline phosphatase", "placental.type alkaline phosphatase",
  "placental alkaline phosphatase", "alkaline phosphatase, placental",
  "alkaline phosphatase, tissue-nonspecific",
  "urinary albumin", "microalbuminuria", "macroalbuminuria", "albuminuria",
  "albumin-to-creatinine ratio", "albumin excretion", "microalbumin",
  "Chronic kidney disease", "chronic kidney disease", "renal", "kidney disease",
  "end.stage.*kidney", "end stage.*kidney",
  "interstitial cystitis", "multiple sclerosis", "rheumatoid arthritis",
  "inflammatory bowel disease", "Crohn's disease", "ulcerative colitis",
  "ankylosing spondylitis", "psoriasis", "colorectal cancer", "colorectal carcinoma",
  "acute lymphoblastic leukemia", "B-cell lymphoblastic leukemia",
  "elite athletes", "mercapturic acid", "vaginal microbiome", "urea cycle",
  "Biological.*Liver Condition", "Liver Liking", "pancreas iron", "spleen iron"
)

#===============================================================================
# DIRECTORIES  (direction-suffixed so gaining/losing never clobber)
#===============================================================================
out_dir <- file.path("/data/sarkern2/multiome_liver/Seurat/epigenome/gwas_peak_gene",
                     paste0("hepatocyte_", DIRECTION))

gwas_dir <- file.path(out_dir, "gwas_p2g/Liver")
checkpoint_dir <- file.path(gwas_dir, "checkpoints")
bed_dir <- file.path(gwas_dir, "trait_beds")
results_dir <- file.path(gwas_dir, "results")
raw_dir <- file.path(results_dir, "raw")
filtered_dir <- file.path(results_dir, "filtered")

for (d in c(gwas_dir, checkpoint_dir, bed_dir, results_dir, raw_dir, filtered_dir))
  dir.create(d, showWarnings = FALSE, recursive = TRUE)

log_file <- file.path(results_dir, paste0("log_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt"))

#===============================================================================
# Helper Functions
#===============================================================================
write_log <- function(msg) {
  cat(msg, "\n")
  cat(paste0("[", Sys.time(), "] ", msg, "\n"), file = log_file, append = TRUE)
}
save_checkpoint <- function(obj, name) {
  saveRDS(obj, file.path(checkpoint_dir, paste0(name, ".rds")))
  write_log(paste("CHECKPOINT SAVED:", name))
}
load_checkpoint <- function(name) {
  f <- file.path(checkpoint_dir, paste0(name, ".rds"))
  if (file.exists(f)) { write_log(paste("CHECKPOINT LOADED:", name)); return(readRDS(f)) }
  NULL
}
checkpoint_exists <- function(name) file.exists(file.path(checkpoint_dir, paste0(name, ".rds")))
clean_trait_name <- function(trait) {
  clean <- gsub("[^a-zA-Z0-9]", "_", trait)
  clean <- gsub("_+", "_", clean)
  substr(gsub("^_|_$", "", clean), 1, 100)
}

write_log("===============================================================================")
write_log(paste("GWAS-SCAVENGE: LIVER TRAITS - P2G PEAKS, SEEDED FROM AGE-", toupper(DIRECTION), "PEAKS"))
write_log("===============================================================================")
write_log(paste("Direction:", DIRECTION))
write_log(paste("Inclusion patterns:", length(LIVER_PATTERNS), "| Exclusion patterns:", length(EXCLUSION_PATTERNS)))
write_log(paste("FDR threshold:", FDR_THRESHOLD, "| Output:", results_dir))

#===============================================================================
# STEP 1: Load ArchR Project
#===============================================================================
write_log("\n====== STEP 1: Load ArchR Project ======")
proj <- tryCatch(
  loadArchRProject("/data/sarkern2/multiome_liver/Seurat/archR/ArchR_Projects/Step8_Hepatocyte_CCAN_P2G"),
  error = function(e) { write_log(paste("CRITICAL ERROR:", e$message)); stop("Cannot continue without ArchR project") })
write_log(paste("Loaded:", length(getCellNames(proj)), "cells"))

#===============================================================================
# STEP 2: Pull PeakMatrix ONCE -> (a) inline age trend, (b) P2G SCAVENGE SE
#   One non-binarized pull serves both: the age model needs raw counts
#   pseudobulked by Sample; SCAVENGE needs the binarized P2G subset.
#===============================================================================
write_log("\n====== STEP 2: PeakMatrix -> age trend + P2G SE ======")

if (checkpoint_exists("step2_SE_Data") && checkpoint_exists("step2_peaks_p2g") &&
    checkpoint_exists("step2_p2g_stats") && checkpoint_exists("step2_age_res")) {
  SE_Data   <- load_checkpoint("step2_SE_Data")
  peaks_p2g <- load_checkpoint("step2_peaks_p2g")
  p2g_stats <- load_checkpoint("step2_p2g_stats")
  age_res   <- load_checkpoint("step2_age_res")
  write_log(paste("Loaded from checkpoint:", nrow(SE_Data), "peaks x", ncol(SE_Data), "cells"))
} else {
  p2g <- getPeak2GeneLinks(ArchRProj = proj, returnLoops = FALSE)
  p2g_peak_idx <- unique(p2g$idxATAC)
  p2g_stats <- list(n_links = nrow(p2g), n_peaks = length(p2g_peak_idx),
                    n_genes = length(unique(p2g$idxRNA)))
  write_log(paste("P2G links:", p2g_stats$n_links, "| peaks:", p2g_stats$n_peaks,
                  "| genes:", p2g_stats$n_genes))

  # ---- one pull, RAW counts (needed for the age model) ----
  peakMat <- getMatrixFromProject(ArchRProj = proj, useMatrix = "PeakMatrix", binarize = FALSE)
  write_log(paste("Full PeakMatrix:", nrow(peakMat), "peaks x", ncol(peakMat), "cells"))
  pk_all <- rowRanges(peakMat)

  # ---------------------------------------------------------------------------
  # STEP 2a: INLINE AGE TREND (pseudobulk by Sample, sex-adjusted limma-voom)
  # ---------------------------------------------------------------------------
  write_log("---- Step 2a: age trend (limma-voom) ----")
  cdat <- as.data.frame(colData(peakMat))
  samp <- factor(cdat[[SAMPLE_COL]])
  pb   <- as.matrix(assay(peakMat) %*% t(Matrix::fac2sparse(samp)))   # peaks x samples
  colnames(pb) <- levels(samp)

  meta <- cdat[!duplicated(samp), c(SAMPLE_COL, SEX_COL, AGE_COL)]
  meta <- meta[match(colnames(pb), meta[[SAMPLE_COL]]), ]
  meta$age_rank <- as.integer(factor(meta[[AGE_COL]], levels = age_levels, ordered = TRUE))
  meta$sex      <- factor(meta[[SEX_COL]])
  stopifnot(identical(colnames(pb), meta[[SAMPLE_COL]]))
  write_log(paste("  pseudobulk:", nrow(pb), "peaks x", ncol(pb), "samples"))

  design <- model.matrix(~ sex + age_rank, data = meta)
  dge  <- DGEList(pb)
  keep <- filterByExpr(dge, design)
  dge  <- dge[keep, , keep.lib.sizes = FALSE]
  pkf  <- pk_all[keep]
  dge  <- calcNormFactors(dge)
  v    <- voom(dge, design)
  fit  <- eBayes(lmFit(v, design))
  tt   <- topTable(fit, coef = "age_rank", number = Inf, sort.by = "none")

  age_res <- data.frame(
    peak_id     = paste0(seqnames(pkf), ":", start(pkf), "-", end(pkf)),
    logFC       = tt$logFC,
    adj.P.Val   = tt$adj.P.Val,
    progressive = ifelse(tt$adj.P.Val < DA_FDR & tt$logFC >  DA_LFC, "gaining",
                  ifelse(tt$adj.P.Val < DA_FDR & tt$logFC < -DA_LFC, "losing", "stable")),
    stringsAsFactors = FALSE)
  write_log(paste("  age peaks -> gaining:", sum(age_res$progressive == "gaining"),
                  "| losing:", sum(age_res$progressive == "losing"),
                  "| stable:", sum(age_res$progressive == "stable")))

  # ---------------------------------------------------------------------------
  # STEP 2: binarized P2G SE for SCAVENGE (binarize in-memory, no second pull)
  # ---------------------------------------------------------------------------
  SE_Data <- peakMat[p2g_peak_idx, ]
  assay(SE_Data) <- as((assay(SE_Data) > 0) * 1, "CsparseMatrix")
  if (!"counts" %in% assayNames(SE_Data)) assayNames(SE_Data)[1] <- "counts"
  peaks_p2g <- rowRanges(SE_Data)

  save_checkpoint(SE_Data, "step2_SE_Data")
  save_checkpoint(peaks_p2g, "step2_peaks_p2g")
  save_checkpoint(p2g_stats, "step2_p2g_stats")
  save_checkpoint(age_res, "step2_age_res")
  rm(peakMat, pb, dge, v, fit); gc()
  write_log(paste("SE_Data (P2G subset, binarized):", nrow(SE_Data), "peaks x", ncol(SE_Data), "cells"))
}

# ---- TARGET peak ids for this DIRECTION (from the inline age trend) ----
TARGET_PEAK_IDS <- unique(age_res$peak_id[age_res$progressive == DIRECTION])
write_log(paste("Age-", DIRECTION, "peaks:", length(TARGET_PEAK_IDS)))

# ---- align TARGET_PEAK_IDS to the SE peak order (logical mask) ----
se_peak_ids <- paste0(seqnames(peaks_p2g), ":", start(peaks_p2g), "-", end(peaks_p2g))
TARGET_MASK <- se_peak_ids %in% TARGET_PEAK_IDS

# --- mask verification (peak-ID format must match between age-DA and SE) ---
write_log("---- TARGET_MASK verification ----")
write_log(paste("  example age-DA peak_id :", paste(head(TARGET_PEAK_IDS, 3), collapse = " | ")))
write_log(paste("  example SE peak id     :", paste(head(se_peak_ids, 3), collapse = " | ")))
write_log(paste("  TARGET_PEAK_IDS (total):", length(TARGET_PEAK_IDS)))
write_log(paste("  TARGET_MASK TRUE       :", sum(TARGET_MASK), "of", length(se_peak_ids)))
write_log(paste("  matched / available    :",
                round(100 * sum(TARGET_MASK) / length(TARGET_PEAK_IDS), 1),
                "% of age-", DIRECTION, "peaks are P2G-linked (most are not, expected)"))
# Hard stop on total failure. A low % of genome-wide age peaks being P2G-linked
# is EXPECTED (only ~10k of ~260k peaks are P2G), so only warn if NONE matched.
if (sum(TARGET_MASK) == 0)
  stop("TARGET_MASK is empty -- peak-ID format mismatch between the age trend and the SE peak set.")

#===============================================================================
# STEP 2b: GC bias + background peaks  (computed ONCE, reused for every trait)
#   The peak matrix is identical across all traits, so addGCBias and
#   getBackgroundPeaks only need to run once -- not inside the per-trait loop.
#===============================================================================
write_log("\n====== STEP 2b: GC bias + background peaks (one-time) ======")

if (checkpoint_exists("step2b_SE_gc") && checkpoint_exists("step2b_SE_bg")) {
  SE_gc <- load_checkpoint("step2b_SE_gc")
  SE_bg <- load_checkpoint("step2b_SE_bg")
  write_log("Loaded GC-biased SE + background peaks from checkpoint")
} else {
  SE_gc <- addGCBias(SE_Data, genome = BSgenome.Mmusculus.UCSC.mm10)
  SE_bg <- getBackgroundPeaks(SE_gc, niterations = 200)
  if (!is.matrix(SE_bg) || nrow(SE_bg) != nrow(SE_gc))
    stop("Background peaks validation failed (nrow mismatch)")
  save_checkpoint(SE_gc, "step2b_SE_gc")
  save_checkpoint(SE_bg, "step2b_SE_bg")
  write_log(paste("Computed background peaks:", nrow(SE_bg), "x", ncol(SE_bg)))
}

#===============================================================================
# STEP 3: Load GWAS Catalog & Filter Liver Traits
#===============================================================================
write_log("\n====== STEP 3: Load GWAS Catalog & Filter Liver Traits ======")

if (checkpoint_exists("step3_gwas_filtered") && checkpoint_exists("step3_liver_traits")) {
  gwas_filtered <- load_checkpoint("step3_gwas_filtered")
  liver_traits <- load_checkpoint("step3_liver_traits")
} else {
  gwas_file <- file.path("/data/sarkern2/multiome_liver/Seurat/epigenome/gwas_peak_gene/hepatocyte",
                         "gwas_catalog_v1.0.2-associations.tsv")
  if (!file.exists(gwas_file)) stop("GWAS Catalog file not found: ", gwas_file)
  write_log(paste("Loading GWAS Catalog from:", gwas_file))
  gwas <- fread(gwas_file, quote = "")
  unique_traits <- unique(gwas$`DISEASE/TRAIT`)
  write_log(paste("Total unique traits:", length(unique_traits)))

  find_traits <- function(patterns, traits)
    unique(unlist(lapply(patterns, function(p) grep(p, traits, value = TRUE, ignore.case = TRUE))))
  liver_traits_raw <- find_traits(LIVER_PATTERNS, unique_traits)
  exclusion_regex <- paste(EXCLUSION_PATTERNS, collapse = "|")
  liver_traits <- liver_traits_raw[!grepl(exclusion_regex, liver_traits_raw, ignore.case = TRUE)]
  write_log(paste("Liver traits:", length(liver_traits_raw), "->", length(liver_traits), "after exclusion"))

  gwas_filtered <- gwas[`DISEASE/TRAIT` %in% liver_traits]
  save_checkpoint(gwas_filtered, "step3_gwas_filtered")
  save_checkpoint(liver_traits, "step3_liver_traits")
  rm(gwas); gc()
}
write_log(paste("GWAS entries for liver traits:", nrow(gwas_filtered)))

#===============================================================================
# STEP 4: LiftOver hg38 -> mm10
#===============================================================================
write_log("\n====== STEP 4: LiftOver hg38 -> mm10 ======")

if (checkpoint_exists("step4_gwas_mm10")) {
  gwas_mm10 <- load_checkpoint("step4_gwas_mm10")
} else {
  gwas_clean <- gwas_filtered[!is.na(CHR_ID) & !is.na(CHR_POS) & CHR_ID != "" & CHR_POS != ""]
  gwas_clean <- gwas_clean[!grepl(";|x", CHR_ID, ignore.case = TRUE)]
  gwas_clean[, CHR_POS := as.numeric(CHR_POS)]
  gwas_clean <- gwas_clean[!is.na(CHR_POS)]
  gwas_gr <- GRanges(seqnames = paste0("chr", gwas_clean$CHR_ID),
                     ranges = IRanges(start = gwas_clean$CHR_POS, width = 1),
                     snp = gwas_clean$SNPS, trait = gwas_clean$`DISEASE/TRAIT`,
                     pvalue = gwas_clean$`P-VALUE`, gene = gwas_clean$MAPPED_GENE)
  gwas_gr <- keepStandardChromosomes(gwas_gr, pruning.mode = "coarse")

  chain_file <- file.path("/data/sarkern2/multiome_liver/Seurat/epigenome/hepatocyte", "hg38ToMm10.over.chain")
  if (!file.exists(chain_file))
    chain_file <- file.path("/data/sarkern2/multiome_liver/Seurat/epigenome/gwas_peak_gene/hepatocyte",
                            "hg38ToMm10.over.chain")
  chain <- import.chain(chain_file)
  gwas_mm10 <- unlist(liftOver(gwas_gr, chain))
  write_log(paste("Mapped to mm10:", length(gwas_mm10)))
  save_checkpoint(gwas_mm10, "step4_gwas_mm10")
}

#===============================================================================
# STEP 5: LD Expansion
#===============================================================================
write_log("\n====== STEP 5: LD Expansion ======")

if (checkpoint_exists("step5_gwas_mm10_ld")) {
  gwas_mm10_ld <- load_checkpoint("step5_gwas_mm10_ld")
} else {
  gwas_mm10_expanded <- gwas_mm10
  start(gwas_mm10_expanded) <- pmax(1, start(gwas_mm10) - LD_WINDOW)
  end(gwas_mm10_expanded) <- end(gwas_mm10) + LD_WINDOW
  unique_traits <- unique(gwas_mm10_expanded$trait)
  merged_list <- vector("list", length(unique_traits))
  for (i in seq_along(unique_traits)) {
    trait <- unique_traits[i]
    trait_gr <- gwas_mm10_expanded[gwas_mm10_expanded$trait == trait]
    trait_reduced <- reduce(trait_gr)
    if (length(trait_reduced) > 0) {
      mcols(trait_reduced)$trait <- trait
      mcols(trait_reduced)$n_snps_merged <- countOverlaps(trait_reduced, trait_gr)
      merged_list[[i]] <- trait_reduced
    }
  }
  merged_list <- merged_list[!sapply(merged_list, is.null)]
  gwas_mm10_ld <- unlist(GRangesList(merged_list), use.names = FALSE)
  write_log(paste("Merged LD regions:", length(gwas_mm10_ld)))
  save_checkpoint(gwas_mm10_ld, "step5_gwas_mm10_ld")
}

#===============================================================================
# STEP 6: Find GWAS Overlaps with P2G Peaks
#===============================================================================
write_log("\n====== STEP 6: Find GWAS Overlaps ======")

if (checkpoint_exists("step6_trait_counts")) {
  trait_counts <- load_checkpoint("step6_trait_counts")
} else {
  peaks_extended <- peaks_p2g
  start(peaks_extended) <- pmax(1, start(peaks_p2g) - PEAK_EXTENSION)
  end(peaks_extended) <- end(peaks_p2g) + PEAK_EXTENSION
  overlaps <- findOverlaps(peaks_extended, gwas_mm10_ld, ignore.strand = TRUE)
  if (length(overlaps) > 0) {
    overlap_dt <- data.table(peak_idx = queryHits(overlaps), gwas_idx = subjectHits(overlaps),
                             trait = gwas_mm10_ld$trait[subjectHits(overlaps)])
    trait_counts <- overlap_dt[, .(n_regions = uniqueN(gwas_idx), n_peaks = uniqueN(peak_idx)), by = trait]
    trait_counts[, passes_filter := n_regions >= MIN_REGIONS]
  } else {
    trait_counts <- data.table(trait = character(0), n_regions = integer(0),
                               n_peaks = integer(0), passes_filter = logical(0))
  }
  save_checkpoint(trait_counts, "step6_trait_counts")
}
passing_traits <- trait_counts[passes_filter == TRUE, trait]
write_log(paste("Traits passing filter (>=", MIN_REGIONS, "regions):", length(passing_traits)))

#===============================================================================
# STEP 7: Create BED Files
#===============================================================================
write_log("\n====== STEP 7: Create BED Files ======")

if (checkpoint_exists("step7_bed_files")) {
  trait_bed_files <- load_checkpoint("step7_bed_files")
} else {
  create_bed <- function(trait_name, gwas_gr, out_dir) {
    trait_idx <- which(gwas_gr$trait == trait_name)
    if (length(trait_idx) == 0) return(NULL)
    trait_gr <- gwas_gr[trait_idx]
    raw_score <- if (!is.null(mcols(trait_gr)$n_snps_merged)) mcols(trait_gr)$n_snps_merged else rep(1, length(trait_gr))
    score <- if (length(unique(raw_score)) > 1 && max(raw_score) > 1) log1p(raw_score) else raw_score
    bed8 <- data.frame(chr = as.character(seqnames(trait_gr)), start = start(trait_gr), end = end(trait_gr),
                       name = paste0(gsub("[^a-zA-Z0-9]", "_", trait_name), "_", seq_along(trait_gr)),
                       score = score, strand = ".", thickStart = start(trait_gr), thickEnd = end(trait_gr))
    clean_name <- substr(gsub("_+", "_", gsub("[^a-zA-Z0-9]", "_", trait_name)), 1, 100)
    bed_file <- file.path(out_dir, paste0(clean_name, ".bed"))
    write.table(bed8, bed_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
    bed_file
  }
  trait_bed_files <- list()
  for (i in seq_along(passing_traits)) {
    bed_file <- tryCatch(create_bed(passing_traits[i], gwas_mm10_ld, bed_dir), error = function(e) NULL)
    if (!is.null(bed_file)) trait_bed_files[[passing_traits[i]]] <- bed_file
  }
  write_log(paste("Total BED files:", length(trait_bed_files)))
  save_checkpoint(trait_bed_files, "step7_bed_files")
}

#===============================================================================
# STEP 8: Create Run List
#===============================================================================
run_list <- data.frame(trait = names(trait_bed_files), stringsAsFactors = FALSE)
write_log(paste("\nTotal traits to run:", nrow(run_list)))

#===============================================================================
# STEP 9: Define SCAVENGE Function  (seeding restricted to TARGET_MASK)
#===============================================================================
Run_SCAVENGE_p2g <- function(SE_subset, SE_bg, trait_file, mycores = 55,
                             min_seed = 2, min_cells = 10) {
  result <- tryCatch({
    # SE_subset is already GC-biased; SE_bg is the precomputed background (Step 2b)
    if (!is.matrix(SE_bg) || nrow(SE_bg) != nrow(SE_subset))
      return(list(success = FALSE, error = "Background peaks validation failed"))

    trait_df <- read.table(trait_file, sep = "\t", stringsAsFactors = FALSE)
    if (nrow(trait_df) == 0) return(list(success = FALSE, error = "Empty BED"))
    colnames(trait_df) <- c("chr", "start", "end", "name", "score", "strand", "thickStart", "thickEnd")
    trait_gr <- GRanges(seqnames = trait_df$chr, ranges = IRanges(start = trait_df$start, end = trait_df$end), score = trait_df$score)
    peaks <- rowRanges(SE_subset)
    overlaps <- findOverlaps(peaks, trait_gr)
    if (length(overlaps) == 0) return(list(success = FALSE, error = "No overlaps"))

    trait_import <- rep(0, length(peaks))
    scores_by_peak <- tapply(trait_gr$score[subjectHits(overlaps)], queryHits(overlaps), sum)
    trait_import[as.integer(names(scores_by_peak))] <- scores_by_peak

    #-----------------------------------------------------------------------
    # AGE SEEDING MASK: keep trait signal only on age-DIRECTION peaks.
    # TARGET_MASK is aligned to rowRanges(SE_Data) == peaks here.
    #-----------------------------------------------------------------------
    trait_import[!TARGET_MASK] <- 0

    n_peaks_with_signal <- sum(trait_import > 0)
    if (n_peaks_with_signal < 1)
      return(list(success = FALSE, error = "No age-direction peaks with signal"))

    trait_import <- Matrix(trait_import, ncol = 1, sparse = TRUE)
    colnames(trait_import) <- "trait"

    SE_DEV <- computeWeightedDeviations(SE_subset, trait_import, background_peaks = SE_bg)
    z <- assays(SE_DEV)[["z"]]
    if (sum(is.na(z)) / length(z) * 100 > 90)
      return(list(success = FALSE, error = "NA rate > 90%"))

    z_score_mat <- data.frame(cell = colnames(z), z_score = as.numeric(z),
                              colData(SE_subset)[colnames(z), , drop = FALSE], row.names = colnames(z))
    z_score_mat <- z_score_mat[complete.cases(z_score_mat), ]
    if (nrow(z_score_mat) < min_cells) return(list(success = FALSE, error = paste("Only", nrow(z_score_mat), "cells")))

    seed_idx <- seedindex(z_score_mat$z_score, 0.05)
    if (sum(seed_idx) < min_seed) return(list(success = FALSE, error = paste("Only", sum(seed_idx), "seeds")))
    scale_factor <- cal_scalefactor(z_score = z_score_mat$z_score, 0.01)

    peak_by_cell_mat <- assay(SE_subset)[, rownames(z_score_mat)]
    tfidf_mat <- tfidf(bmat = peak_by_cell_mat, mat_binary = TRUE, TF = TRUE, log_TF = TRUE)
    lsi_mat <- do_lsi(tfidf_mat, dims = 30)
    mutualknn30 <- getmutualknn(lsi_mat, 30)

    seed_cells <- intersect(rownames(z_score_mat)[seed_idx], rownames(mutualknn30))
    if (length(seed_cells) < min_seed) return(list(success = FALSE, error = paste("Only", length(seed_cells), "valid seeds")))

    np_score <- randomWalk_sparse(intM = mutualknn30, queryCells = seed_cells, gamma = 0.05)
    omit_idx <- np_score == 0
    if (sum(!omit_idx) < min_cells) return(list(success = FALSE, error = paste("Only", sum(!omit_idx), "after walk")))

    mutualknn30 <- mutualknn30[!omit_idx, !omit_idx]
    np_score <- np_score[!omit_idx]
    TRS <- np_score %>% capOutlierQuantile(., 0.95) %>% max_min_scale
    TRS <- TRS * scale_factor

    median_TRS <- median(TRS, na.rm = TRUE)
    if (median_TRS < 0.001) return(list(success = FALSE, error = paste("Random walk collapse - median TRS:", round(median_TRS, 6))))

    filtered_cells <- names(np_score)
    mono_mat <- data.frame(z_score_mat[filtered_cells, ], seed_idx = filtered_cells %in% seed_cells,
                           np_score = np_score, TRS = TRS)

    get_sigcell_custom <- function(knn_sparse_mat, seed_idx, topseed_npscore,
                                   permutation_times = 1000, true_cell_significance = 0.05,
                                   mycores = 55, rw_gamma = 0.05) {
      cell_mat <- data.frame(cell = 1:nrow(knn_sparse_mat), degree = colSums(knn_sparse_mat))
      seed_mat_top <- data.frame(seed = which(seed_idx),
                                 degree = colSums(knn_sparse_mat[, which(seed_idx == TRUE), drop = FALSE]))
      seed_table_top <- data.frame(table(seed_mat_top$degree))
      if (nrow(seed_table_top) == 0) return(data.frame(seed_idx = seed_idx, true_cell_top_idx = FALSE, empirical_pval = 1))
      xx_top <- tapply(cell_mat[, 1], cell_mat[, 2], list)
      xx2_top <- xx_top[names(xx_top) %in% seed_table_top$Var1]
      if (length(xx2_top) == 0) return(data.frame(seed_idx = seed_idx, true_cell_top_idx = FALSE, empirical_pval = 1))
      permutation_score_top <- mclapply(1:permutation_times, mc.cores = mycores, function(i) {
        tryCatch({
          sampled_cellid <- xx2_top %>% mapply(sample, ., seed_table_top$Freq) %>% unlist() %>% sort()
          randomWalk_sparse(intM = knn_sparse_mat, queryCells = rownames(knn_sparse_mat)[as.numeric(sampled_cellid)], gamma = rw_gamma)
        }, error = function(e) rep(NA, nrow(knn_sparse_mat)))
      })
      permutation_score_top <- data.frame(sapply(permutation_score_top, c))
      valid_perms <- colSums(is.na(permutation_score_top)) == 0
      if (sum(valid_perms) == 0) return(data.frame(seed_idx = seed_idx, true_cell_top_idx = FALSE, empirical_pval = 1))
      permutation_score_top <- permutation_score_top[, valid_perms, drop = FALSE]
      permutation_df_top <- apply(permutation_score_top, 2, function(x) x > topseed_npscore)
      if (is.null(dim(permutation_df_top))) permutation_df_top <- matrix(permutation_df_top, ncol = 1)
      empirical_pval <- rowSums(permutation_df_top) / ncol(permutation_df_top)
      data.frame(seed_idx, true_cell_top_idx = empirical_pval <= true_cell_significance, empirical_pval)
    }

    mono_permu <- get_sigcell_custom(mutualknn30, mono_mat$seed_idx, mono_mat$np_score, 1000, 0.05, mycores, 0.05)
    mono_mat2 <- data.frame(mono_mat, mono_permu)
    mono_mat2$fdr <- p.adjust(mono_mat2$empirical_pval, method = "BH")

    list(success = TRUE, result = mono_mat2, n_significant = sum(mono_mat2$true_cell_top_idx),
         n_cells = nrow(mono_mat2), n_peaks = nrow(SE_subset),
         n_peaks_with_signal = n_peaks_with_signal, n_seed = length(seed_cells), median_TRS = median_TRS)
  }, error = function(e) list(success = FALSE, error = e$message))
  result
}

#===============================================================================
# STEP 10: Run SCAVENGE
#===============================================================================
write_log("\n====== STEP 10: Running SCAVENGE ======")

results_file <- file.path(results_dir, "results_Liver_P2G_ALL.rds")
summary_file <- file.path(results_dir, "summary_Liver_P2G.csv")
skipped_file <- file.path(results_dir, "skipped_Liver_P2G.csv")

if (file.exists(results_file)) {
  scavenge_results <- readRDS(results_file)
  scavenge_summary <- read.csv(summary_file, stringsAsFactors = FALSE)
  skipped_runs <- if (file.exists(skipped_file)) read.csv(skipped_file, stringsAsFactors = FALSE) else data.frame()
  write_log(paste("Resuming:", nrow(scavenge_summary), "completed"))
} else {
  scavenge_results <- list(); scavenge_summary <- data.frame(); skipped_runs <- data.frame()
}

remaining <- run_list[!run_list$trait %in% scavenge_summary$trait, , drop = FALSE]
write_log(paste("Completed:", nrow(scavenge_summary), "| Remaining:", nrow(remaining)))

if (nrow(remaining) > 0) {
  start_time <- Sys.time()
  for (i in 1:nrow(remaining)) {
    trait <- remaining$trait[i]
    bed_file <- trait_bed_files[[trait]]
    trait_info <- trait_counts[trait == ..trait]
    n_regions <- if (nrow(trait_info) > 0) trait_info$n_regions[1] else 0
    write_log(paste0("\n[", i, "/", nrow(remaining), "] ", trait))

    run_result <- tryCatch(Run_SCAVENGE_p2g(SE_gc, SE_bg, bed_file, SCAVENGE_CORES, MIN_SEED_CELLS, MIN_VALID_CELLS),
                           error = function(e) list(success = FALSE, error = e$message))

    if (run_result$success) {
      result <- run_result$result
      scavenge_results[[trait]] <- result
      n_fdr_sig <- sum(result$fdr < FDR_THRESHOLD, na.rm = TRUE)
      n_true_cell <- sum(result$true_cell_top_idx, na.rm = TRUE)
      median_TRS <- if (!is.null(run_result$median_TRS)) run_result$median_TRS else median(result$TRS, na.rm = TRUE)

      scavenge_summary <- bind_rows(scavenge_summary, data.frame(
        trait = trait, category = CATEGORY, direction = DIRECTION, n_regions = n_regions,
        n_peaks = run_result$n_peaks, n_peaks_with_signal = run_result$n_peaks_with_signal,
        n_cells = run_result$n_cells, n_seed = run_result$n_seed,
        mean_zscore = mean(result$z_score, na.rm = TRUE),
        n_significant = run_result$n_significant,
        pct_significant = run_result$n_significant / run_result$n_cells * 100,
        n_fdr_sig = n_fdr_sig, pct_fdr_sig = n_fdr_sig / run_result$n_cells * 100,
        n_true_cell = n_true_cell, pct_true_cell = n_true_cell / run_result$n_cells * 100,
        median_TRS = median_TRS, mean_TRS = mean(result$TRS, na.rm = TRUE),
        status = "success", stringsAsFactors = FALSE))
      write_log(paste("  SUCCESS:", n_true_cell, "sig | FDR<0.05:", n_fdr_sig, "| medTRS:", round(median_TRS, 4)))

      clean_name <- clean_trait_name(trait)
      result$trait <- trait
      fwrite(result, file.path(raw_dir, paste0(clean_name, "_raw.csv")))
      result_filtered <- result[result$fdr < FDR_THRESHOLD, ]
      if (nrow(result_filtered) > 0) fwrite(result_filtered, file.path(filtered_dir, paste0(clean_name, "_filtered.csv")))
    } else {
      write_log(paste("  SKIPPED:", run_result$error))
      skipped_runs <- bind_rows(skipped_runs, data.frame(
        trait = trait, category = CATEGORY, direction = DIRECTION,
        reason = run_result$error, stringsAsFactors = FALSE))
    }

    if (i %% 5 == 0 || i == nrow(remaining)) {
      saveRDS(scavenge_results, results_file)
      write.csv(scavenge_summary, summary_file, row.names = FALSE)
      write.csv(skipped_runs, skipped_file, row.names = FALSE)
      write_log(paste("  [SAVED] Elapsed:", round(as.numeric(difftime(Sys.time(), start_time, units = "hours")), 2), "h"))
    }
  }
}

saveRDS(scavenge_results, results_file)
write.csv(scavenge_summary, summary_file, row.names = FALSE)
write.csv(skipped_runs, skipped_file, row.names = FALSE)

#===============================================================================
# STEP 11: Create Combined Output Files
#===============================================================================
write_log("\n====== STEP 11: Create Combined Output Files ======")

raw_files <- list.files(raw_dir, pattern = "_raw.csv$", full.names = TRUE)
if (length(raw_files) > 0) {
  all_raw <- rbindlist(lapply(raw_files, fread), fill = TRUE)
  fwrite(all_raw, file.path(results_dir, "combined_ALL_raw.csv"))
  write_log(paste("  Combined raw:", nrow(all_raw), "rows from", length(raw_files), "traits"))
}
filtered_files <- list.files(filtered_dir, pattern = "_filtered.csv$", full.names = TRUE)
if (length(filtered_files) > 0) {
  all_filtered <- rbindlist(lapply(filtered_files, fread), fill = TRUE)
  fwrite(all_filtered, file.path(results_dir, "combined_ALL_filtered.csv"))
  write_log(paste("  Combined filtered:", nrow(all_filtered), "rows from", length(filtered_files), "traits"))

  filtered_summary <- all_filtered[, .(n_sig_cells = .N, mean_TRS = mean(TRS, na.rm = TRUE),
                                       median_TRS = median(TRS, na.rm = TRUE),
                                       mean_fdr = mean(fdr, na.rm = TRUE)), by = trait][order(-n_sig_cells)]
  fwrite(filtered_summary, file.path(results_dir, "summary_filtered_by_trait.csv"))
}

write_log("\n===============================================================================")
write_log(paste("LIVER P2G GWAS-SCAVENGE (age-", DIRECTION, ") COMPLETE"))
write_log(paste("Successful:", nrow(scavenge_summary), "| Skipped:", nrow(skipped_runs)))
write_log("===============================================================================")
gc()
write_log(paste("Finished:", Sys.time()))
