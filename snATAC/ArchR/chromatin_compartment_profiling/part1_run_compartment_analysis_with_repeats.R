#!/usr/bin/env Rscript
# ==============================================================================
# Chromatin Compartment Analysis with Repeat Element Integration
# ==============================================================================
#
# Description:
#   Computes A/B compartment scores across aging timepoints using HMM
#   segmentation of binned bigWig signal, classifies compartment switching
#   (Stable_Active, Stable_Repressive, Monotonic_A→R, Monotonic_R→A,
#   Non_Monotonic), and integrates ChromHMM and RepeatMasker annotations
#   for enrichment visualization.
#
# Input:
#   - ArchR project with GroupBigWigs (AgeSexCelltype, 80 kb tiles)
#   - mm10 blacklist BED (ENCODE)
#   - ChromHMM segmentation BED (Full-Stack)
#   - RepeatMasker processed RDS (from repeat_mask_download.R)
#
# Output:
#   - Compartment binary matrix (Active=1, Repressive=0)
#   - ComplexHeatmap (main + ChromHMM + RepeatMasker panels)
#   - Stability classification (5-class model)
#   - Repeat enrichment by compartment and stability class
#   - Alluvial plots of compartment transitions
#
# Pipeline Steps:
#   1.  Load manifest
#   2.  Create genomic bins (80 kb, blacklist/gap-filtered)
#   3.  Extract signal from bigWigs
#   4.  HMM segmentation (BaumWelch per chromosome)
#   5.  ChromHMM enrichment (log2 fold-enrichment per bin)
#   6.  RepeatMasker enrichment (4 classes: SINE, DNA, LTR, LINE)
#   7.  CpG density and AT content
#   8.  K-means clustering of compartment patterns
#   9.  Bin ordering and ECG relabeling
#   10. Build annotations (row + column)
#   11. Build ChromHMM heatmap panel
#   12. Build RepeatMasker heatmap panel
#   13. Assemble and save heatmap
#   14. Stability analysis (5-class model)
#   15. Repeat enrichment by stability class
#   16. Stability barplots
#   17. Alluvial plots
#   18. Session info
#
# Reference:
#   Sarker et al. (2026) Cell Metabolism
#
# ==============================================================================


# ==============================================================================
# SETUP AND CONFIGURATION
# ==============================================================================

cat("\n")
cat("========================================================================\n")
cat("  COMPARTMENT ANALYSIS PIPELINE - FULL RUN\n")
cat("  Mouse Liver Aging Single-Nucleus Multi-Ome Study\n")
cat("  WITH REPEATMASKER ENRICHMENT\n")
cat("========================================================================\n\n")

# ------------------------------------------------------------------------------
# User Configuration - Modify These Paths
# ------------------------------------------------------------------------------

# Project paths
proj_path      <- "path/to/ArchR_Projects/Step6_Xwnn_UMAP"
group_dir      <- file.path(proj_path, "GroupBigWigs", "AgeSexCelltype__tile80k")
manifest_path  <- file.path(group_dir, "AgeSexCelltype.bigwig_manifest_tile80k.tsv")
outdir         <- "path/to/output/compartment_main_heatmap"

# Reference files
mm10_blacklist_bed <- "path/to/mm10-blacklist.v2.bed.gz"
fs_bed             <- file.path(proj_path, "mm10_100_segments_segments.bed.gz")

# RepeatMasker file (from repeat_mask_download.R)
repeatmasker_rds   <- "path/to/repeatmasker/mm10_repeatmasker_processed.rds"

# Analysis parameters
tile_size      <- 80000L
non_gap_thresh <- 0.75
k_bins         <- 25L
threads        <- 60L

# Stability thresholds
t_on           <- 0.6
t_off          <- 0.4

# Random seed
seed           <- 10918


# ------------------------------------------------------------------------------
# Build CONFIG list
# ------------------------------------------------------------------------------

CONFIG <- list(
    proj_path        = proj_path,
    group_dir        = group_dir,
    manifest_path    = manifest_path,
    outdir           = outdir,
    mm10_blacklist   = mm10_blacklist_bed,
    chromhmm_bed     = fs_bed,
    repeatmasker_rds = repeatmasker_rds,
    tile_size        = tile_size,
    non_gap_thresh   = non_gap_thresh,
    k_bins           = k_bins,
    threads          = threads,
    t_on             = t_on,
    t_off            = t_off,
    seed             = seed
)


# ------------------------------------------------------------------------------
# Validate Paths
# ------------------------------------------------------------------------------

message("[SETUP] Validating paths...")

paths_to_check <- list(
    "Project directory"  = proj_path,
    "BigWig directory"   = group_dir,
    "Manifest file"      = manifest_path,
    "Blacklist BED"      = mm10_blacklist_bed,
    "ChromHMM BED"       = fs_bed,
    "RepeatMasker RDS"   = repeatmasker_rds
)

all_exist <- TRUE
for (name in names(paths_to_check)) {
    path <- paths_to_check[[name]]
    exists_flag <- file.exists(path)
    status <- if (exists_flag) "[OK]" else "[MISSING]"
    message(sprintf("  %s %s: %s", status, name, path))
    if (!exists_flag) all_exist <- FALSE
}

if (!all_exist) {
    stop("\n[ERROR] Some required paths do not exist. Please check paths above.\n",
         "       If RepeatMasker is missing, run repeat_mask_download.R first.")
}

if (!dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE)
    message(sprintf("  [OK] Created output directory: %s", outdir))
} else {
    message(sprintf("  [OK] Output directory exists: %s", outdir))
}


# ------------------------------------------------------------------------------
# Load Packages
# ------------------------------------------------------------------------------

message("\n[SETUP] Loading required packages...")

required_packages <- c(
    "BSgenome.Mmusculus.UCSC.mm10",
    "GenomicRanges",
    "IRanges",
    "GenomeInfoDb",
    "rtracklayer",
    "data.table",
    "stringr",
    "circlize",
    "ComplexHeatmap",
    "HMMt",
    "Biostrings",
    "ggplot2",
    "ggalluvial",
    "cowplot",
    "ArchR"
)

missing_packages <- c()
for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        missing_packages <- c(missing_packages, pkg)
    }
}

if (length(missing_packages) > 0) {
    message("\n[WARNING] Missing packages:")
    for (pkg in missing_packages) {
        message(sprintf("  - %s", pkg))
    }
    message("\nInstall with:")
    message('  BiocManager::install(c("', paste(missing_packages, collapse = '", "'), '"))')
    stop("\n[ERROR] Please install missing packages before running.")
}

suppressPackageStartupMessages({
    library(BSgenome.Mmusculus.UCSC.mm10)
    library(GenomicRanges)
    library(IRanges)
    library(GenomeInfoDb)
    library(rtracklayer)
    library(data.table)
    library(stringr)
    library(circlize)
    library(ComplexHeatmap)
    library(HMMt)
    library(Biostrings)
    library(ggplot2)
    library(ggalluvial)
    library(cowplot)
    library(ArchR)
})

addArchRThreads(threads = CONFIG$threads)
message("  [OK] All packages loaded")


# ==============================================================================
# SHARED CONSTANTS
# ==============================================================================

message("\n[SETUP] Setting up shared constants...")

# Age levels
AGE_LEVELS <- c("young", "mid_age", "old", "pre_geriatric", "geriatric")

# Sex levels
SEX_LEVELS <- c("male", "female", "unknown")

# Stability classes
STABILITY_LEVELS <- c(
    "Stable_Active",
    "Stable_Repressive",
    "Monotonic_A_to_R",
    "Monotonic_R_to_A",
    "Non_Monotonic"
)

SWITCHING_CLASSES <- c("Monotonic_A_to_R", "Monotonic_R_to_A", "Non_Monotonic")

NONMONO_PATTERNS <- c(
    "Transient_Return", "Two_Switch", "Triple_Switch",
    "Highly_Dynamic", "Other"
)

# Color palettes
AGE_COLS <- c(
    "young"         = "#1ABC9C",
    "mid_age"       = "#F1C40F",
    "old"           = "#C39BD3",
    "pre_geriatric" = "#2980B9",
    "geriatric"     = "#E84393"
)

SEX_COLS <- c(
    "male"    = "#2E86C1",
    "female"  = "#EC7063",
    "unknown" = "#B2BABB"
)

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

STATE_COLS <- c(
    "Active"     = "#FF6B6B",
    "Repressive" = "#4E79A7"
)

STABILITY_COLS <- c(
    "Stable_Active"     = "#FF6B6B",
    "Stable_Repressive" = "#4E79A7",
    "Monotonic_A_to_R"  = "darkred",
    "Monotonic_R_to_A"  = "darkblue",
    "Non_Monotonic"     = "#708238",
    "Missing"           = "grey80",
    "Insufficient_Data" = "grey60"
)

SWITCHING_COLS <- c(
    "Monotonic_A_to_R" = "darkred",
    "Monotonic_R_to_A" = "darkblue",
    "Non_Monotonic"    = "#708238"
)

NONMONO_PATTERN_COLS <- c(
    "Transient_Return" = "#2E8B57",
    "Two_Switch"       = "#E67E22",
    "Triple_Switch"    = "#CD5C5C",
    "Highly_Dynamic"   = "#DDA0DD",
    "Other"            = "#95A5A6"
)

COMPARTMENT_COLS <- c(
    "0" = "#0A2A43",
    "1" = "#F2D28B"
)

# Repeat element colors (euchromatin -> heterochromatin)
# Satellite removed (centromeric artifact)
REPEAT_COLS <- c(
    "SINE"           = "#00CED1",
    "DNA_Transposon" = "#FF6347",
    "LTR"            = "#9370DB",
    "LINE"           = "#20B2AA"
)

REPEAT_CLASSES_KEEP <- c("SINE", "DNA_Transposon", "LTR", "LINE")

# Chromosome colors
chr_order <- paste0("chr", c(1:19, "X", "Y"))
n_chr     <- length(chr_order)
hues      <- seq(0, 360 - 360 / n_chr, length.out = n_chr)
CHR_COLS  <- setNames(grDevices::hcl(h = hues, c = 90, l = 60), chr_order)

message("  [OK] Constants defined")


# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================

banner <- function(text, char = "=", width = 70) {
    line <- paste(rep(char, width), collapse = "")
    message("\n", line)
    message(text)
    message(line, "\n")
}

ensure_dir <- function(dir_path) {
    if (!dir.exists(dir_path)) {
        dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
    }
    invisible(dir_path)
}

clean_name <- function(ct) {
    gsub("[^A-Za-z0-9_]", "_", ct)
}

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

write_bed <- function(dt, filename, name_col = "bin_id", score_col = NULL) {
    bed <- data.table::data.table(
        chr   = dt$chr,
        start = dt$start - 1,
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
    message(sprintf("    -> %s (%d regions)", basename(filename), nrow(bed)))
}

per_bin_active_strict <- function(comp_bin, meta_vec, categories, ord_bins_ids,
                                   unknown_label = "Other/NA") {
    stopifnot(all(colnames(comp_bin) %in% names(meta_vec)))
    meta <- as.character(meta_vec[colnames(comp_bin)])
    meta[is.na(meta) | meta == ""] <- unknown_label
    cats_use <- unique(c(categories, setdiff(unique(meta), categories)))
    G <- model.matrix(~ 0 + factor(meta, levels = cats_use))
    colnames(G) <- cats_use
    A  <- as.matrix(comp_bin) %*% G
    rs <- rowSums(comp_bin)
    F_mat <- matrix(0, nrow = nrow(A), ncol = ncol(A),
                    dimnames = list(rownames(A), colnames(A)))
    nz <- rs > 0
    if (any(nz)) F_mat[nz, ] <- sweep(A[nz, , drop = FALSE], 1, rs[nz], "/")
    F_mat <- pmin(pmax(F_mat, 0), 1)
    F_mat <- F_mat[ord_bins_ids, , drop = FALSE]
    t(F_mat)
}

message("  [OK] Helper functions defined")


# ==============================================================================
# SET SEED
# ==============================================================================

set.seed(CONFIG$seed)
message(sprintf("\n[CONFIG] Random seed: %d", CONFIG$seed))
message(sprintf("[CONFIG] Threads: %d", CONFIG$threads))

# Track timing
start_time <- Sys.time()


# ==============================================================================
# STEP 1: LOAD MANIFEST
# ==============================================================================

banner("STEP 1: Loading manifest")

man <- fread(CONFIG$manifest_path)
setnames(man, tolower(names(man)))

if (!all(c("age", "sex", "celltype") %in% names(man))) {
    clean <- sub("-TileSize.*$", "", man$group)
    parts <- strsplit(clean, "__", fixed = TRUE)
    man[, age      := vapply(parts, `[`, character(1), 1)]
    man[, sex      := vapply(parts, `[`, character(1), 2)]
    man[, celltype := vapply(parts, `[`, character(1), 3)]
}

man[, sex := {
    s <- tolower(as.character(sex))
    s[!s %in% c("male", "female")] <- "unknown"
    s
}]
man[, age := factor(as.character(age), levels = AGE_LEVELS)]
man[, celltype := sub("-TileSize.*$", "", celltype)]

fwrite(man, file.path(CONFIG$outdir, "manifest.normalized.tsv"), sep = "\t")

bw_vec    <- setNames(man$bigwig, man$group)
age_vec   <- setNames(as.character(man$age), man$group)
sex_vec   <- setNames(as.character(man$sex), man$group)
ctype_vec <- setNames(as.character(man$celltype), man$group)

message(sprintf("  [OK] %d samples loaded", nrow(man)))
message(sprintf("  [OK] Ages: %s", paste(unique(age_vec), collapse = ", ")))
message(sprintf("  [OK] Sexes: %s", paste(unique(sex_vec), collapse = ", ")))
message(sprintf("  [OK] Cell types: %d unique", length(unique(ctype_vec))))


# ==============================================================================
# STEP 2: CREATE GENOMIC BINS
# ==============================================================================

banner("STEP 2: Creating genomic bins")

seqs <- seqlengths(BSgenome.Mmusculus.UCSC.mm10)
seqs <- seqs[grepl("^chr([0-9]{1,2}|X|Y)$", names(seqs))]
bins100k <- tileGenome(seqs, tilewidth = CONFIG$tile_size,
                       cut.last.tile.in.chrom = TRUE)

mm10_blacklist <- import(CONFIG$mm10_blacklist)
bins_clean <- subsetByOverlaps(bins100k, mm10_blacklist, invert = TRUE)

# Load ChromHMM segmentation
bed_dt <- tryCatch(
    fread(CONFIG$chromhmm_bed, header = FALSE),
    error = function(e) fread(cmd = paste("zcat -f", shQuote(CONFIG$chromhmm_bed)))
)
setnames(bed_dt, c("chr", "start0", "end", "label"))
bed_dt[, start := start0 + 1L]
bed_dt[, mnemonic := sub("[0-9]+$", "", sub("^[^_]*_", "", label))]

fs_gr <- GRanges(
    seqnames = bed_dt$chr,
    ranges   = IRanges(start = bed_dt$start, end = bed_dt$end),
    mnemonic = bed_dt$mnemonic
)
seqlevelsStyle(fs_gr) <- "UCSC"
fs_gr <- keepStandardChromosomes(fs_gr, pruning.mode = "coarse")
fs_gr <- sort(fs_gr)

# Filter bins by non-gap fraction
gap_states <- grep("^mGapArtf", unique(mcols(fs_gr)$mnemonic), value = TRUE)
ol <- findOverlaps(bins_clean, fs_gr, ignore.strand = TRUE)
cov_dt <- data.table(
    bin = queryHits(ol),
    mn  = as.character(mcols(fs_gr)$mnemonic[subjectHits(ol)]),
    w   = width(pintersect(bins_clean[queryHits(ol)], fs_gr[subjectHits(ol)]))
)

bw_all <- width(bins_clean)
cov_sum  <- cov_dt[, .(covered = sum(w)), by = .(bin, mn)]
gap_cov  <- cov_sum[mn %in% gap_states, .(gap_bp = sum(covered)), by = bin]

# Full join: ensure all bins are represented
non_gap <- data.table(bin = seq_along(bins_clean))
non_gap <- gap_cov[non_gap, on = "bin"]
non_gap[is.na(gap_bp), gap_bp := 0L]
non_gap[, non_gap_frac := (bw_all[bin] - gap_bp) / bw_all[bin]]

bins_use <- bins_clean[non_gap$non_gap_frac >= CONFIG$non_gap_thresh]
bins_use$bin_id <- paste0(seqnames(bins_use), ":", start(bins_use), "-", end(bins_use))

bins_df <- data.frame(
    chr   = as.character(seqnames(bins_use)),
    start = start(bins_use),
    end   = end(bins_use),
    row   = seq_len(length(bins_use))
)
bins_df$bin_id <- bins_use$bin_id
rownames(bins_df) <- bins_df$bin_id

export(bins_use, file.path(
    CONFIG$outdir,
    sprintf("bins100k.mm10.nonGap%02d.bed", CONFIG$non_gap_thresh * 100)
))
saveRDS(bins_use, file.path(
    CONFIG$outdir,
    sprintf("bins100k.mm10.nonGap%02d.rds", CONFIG$non_gap_thresh * 100)
))
fwrite(bins_df, file.path(CONFIG$outdir, "bins_table.tsv"), sep = "\t")

message(sprintf("  [OK] %d bins created after filtering", length(bins_use)))


# ==============================================================================
# STEP 3: EXTRACT SIGNAL FROM BIGWIGS
# ==============================================================================

banner("STEP 3: Extracting signal from bigWigs")

bw_means_in_bins <- function(bw, bins) {
    cov <- import(BigWigFile(bw), as = "RleList")
    common <- intersect(names(cov), seqlevels(bins))
    bins2 <- keepSeqlevels(bins, common, pruning.mode = "coarse")
    cov2  <- cov[common]
    out   <- binnedAverage(bins2, cov2, "score")
    as.numeric(mcols(out)$score)
}

mat <- sapply(seq_along(bw_vec), function(i) {
    if (i %% 10 == 0) message(sprintf("  Processing sample %d/%d", i, length(bw_vec)))
    bw_means_in_bins(bw_vec[[i]], bins_use)
})
colnames(mat) <- names(bw_vec)
rownames(mat) <- bins_df$bin_id

mat_z <- scale(mat)
mat_z[is.na(mat_z)] <- 0

saveRDS(mat_z, file.path(CONFIG$outdir, "matrix_signal_z.rds"))

message(sprintf("  [OK] Signal matrix: %d bins x %d samples", nrow(mat_z), ncol(mat_z)))


# ==============================================================================
# STEP 4: HMM SEGMENTATION
# ==============================================================================

banner("STEP 4: Running HMM segmentation")

joint_hmm_one_chr <- function(x_mat_chr) {
    n_bins   <- nrow(x_mat_chr)
    n_groups <- ncol(x_mat_chr)
    x_mat_chr[!is.finite(x_mat_chr)] <- 0
    if (all(x_mat_chr == 0)) {
        return(matrix(0L, nrow = n_bins, ncol = n_groups,
                      dimnames = dimnames(x_mat_chr)))
    }
    fit <- HMMt::BaumWelchT(
        x = as.numeric(x_mat_chr),
        series.length = rep(n_bins, n_groups),
        maxiter = 100
    )
    v      <- fit$ViterbiPath
    active <- which.max(fit$mu)
    matrix(ifelse(v == active, 1L, 0L),
           nrow = n_bins, ncol = n_groups,
           byrow = FALSE, dimnames = dimnames(x_mat_chr))
}

chr_levels <- unique(bins_df$chr)
comp_list <- lapply(chr_levels, function(cc) {
    message(sprintf("  Processing %s", cc))
    idx <- which(bins_df$chr == cc)
    joint_hmm_one_chr(mat_z[idx, , drop = FALSE])
})
comp_bin <- do.call(rbind, comp_list)
colnames(comp_bin) <- colnames(mat_z)
rownames(comp_bin) <- rownames(mat_z)

saveRDS(comp_bin, file.path(CONFIG$outdir, "compartments_binary_A1_R0.rds"))

frac_active <- mean(comp_bin == 1L, na.rm = TRUE)
message(sprintf("  [OK] Compartments: %d bins x %d samples", nrow(comp_bin), ncol(comp_bin)))
message(sprintf("  [OK] Active fraction: %.1f%%", 100 * frac_active))


# ==============================================================================
# STEP 5: CHROMHMM ENRICHMENT
# ==============================================================================

banner("STEP 5: Computing ChromHMM enrichments")

ol2 <- findOverlaps(bins_use, fs_gr, ignore.strand = TRUE)
q2  <- queryHits(ol2)
s2  <- subjectHits(ol2)
ov2 <- width(pintersect(bins_use[q2], fs_gr[s2]))

cov_dt_mn    <- data.table(bin = q2, mn = as.character(mcols(fs_gr)$mnemonic[s2]), w = ov2)
cov_wide_mn  <- dcast(cov_dt_mn, bin ~ mn, value.var = "w", fill = 0L, fun.aggregate = sum)
bw_use       <- width(bins_use)
cov_wide_mn[, total := bw_use[bin]]

mn_cols <- setdiff(names(cov_wide_mn), c("bin", "total"))
for (nm in mn_cols) set(cov_wide_mn, j = nm, value = cov_wide_mn[[nm]] / cov_wide_mn$total)

state_len <- tapply(width(fs_gr), as.character(mcols(fs_gr)$mnemonic), sum, na.rm = TRUE)
state_len <- state_len[!grepl("^mGapArtf", names(state_len))]
state_pct <- state_len / sum(state_len)

keep_states <- intersect(mn_cols, names(state_pct))
fe_raw  <- as.matrix(cov_wide_mn[, ..keep_states])
row_tot <- rowSums(fe_raw)
fe_norm <- sweep(fe_raw, 1, pmax(row_tot, 1e-12), "/")

eps <- 1e-12
log2fe_bin <- log2(sweep(fe_norm + eps, 2, state_pct[colnames(fe_norm)] + eps, "/"))

bin_ids <- bins_use$bin_id
rownames(log2fe_bin) <- bin_ids[cov_wide_mn$bin]
log2fe_bin <- log2fe_bin[bin_ids, , drop = FALSE]

saveRDS(log2fe_bin, file.path(
    CONFIG$outdir,
    sprintf("log2fe_bin.nonGap%02d.rds", CONFIG$non_gap_thresh * 100)
))

message("  [OK] ChromHMM enrichments computed")


# ==============================================================================
# STEP 6: REPEATMASKER ENRICHMENT
# ==============================================================================

banner("STEP 6: Computing RepeatMasker enrichments")

message("  Loading RepeatMasker annotations...")
rmsk_gr <- readRDS(CONFIG$repeatmasker_rds)
message(sprintf("  [OK] Loaded %s repeat elements", format(length(rmsk_gr), big.mark = ",")))

# Add Major_Class if missing
if (!"Major_Class" %in% names(mcols(rmsk_gr))) {
    message("  Adding Major_Class annotations...")
    rmsk_dt_tmp <- data.table(repClass = as.character(mcols(rmsk_gr)$repClass))
    rmsk_dt_tmp[, Major_Class := fcase(
        repClass == "LINE", "LINE",
        repClass == "SINE", "SINE",
        repClass == "LTR",  "LTR",
        repClass == "DNA",  "DNA_Transposon",
        default = "Other"
    )]
    mcols(rmsk_gr)$Major_Class <- rmsk_dt_tmp$Major_Class
}

# Filter to biologically relevant classes
message(sprintf("  Filtering to %d repeat classes: %s",
                length(REPEAT_CLASSES_KEEP),
                paste(REPEAT_CLASSES_KEEP, collapse = ", ")))
rmsk_gr <- rmsk_gr[mcols(rmsk_gr)$Major_Class %in% REPEAT_CLASSES_KEEP]
rmsk_gr <- keepStandardChromosomes(rmsk_gr, pruning.mode = "coarse")
message(sprintf("  [OK] Retained %s repeat elements", format(length(rmsk_gr), big.mark = ",")))

# Overlap with genomic bins
message("  Computing overlap with genomic bins...")
ol_rmsk <- findOverlaps(bins_use, rmsk_gr, ignore.strand = TRUE)
q_rmsk  <- queryHits(ol_rmsk)
s_rmsk  <- subjectHits(ol_rmsk)
ov_rmsk <- width(pintersect(bins_use[q_rmsk], rmsk_gr[s_rmsk]))

cov_dt_rmsk  <- data.table(
    bin = q_rmsk,
    repeat_class = as.character(mcols(rmsk_gr)$Major_Class[s_rmsk]),
    w = ov_rmsk
)
cov_wide_rmsk <- dcast(cov_dt_rmsk, bin ~ repeat_class,
                       value.var = "w", fill = 0L, fun.aggregate = sum)
cov_wide_rmsk[, total := bw_use[bin]]

rmsk_cols <- setdiff(names(cov_wide_rmsk), c("bin", "total"))
for (nm in rmsk_cols) {
    set(cov_wide_rmsk, j = nm, value = cov_wide_rmsk[[nm]] / cov_wide_rmsk$total)
}

# Genome-wide expected fractions
repeat_len <- tapply(width(rmsk_gr), as.character(mcols(rmsk_gr)$Major_Class),
                     sum, na.rm = TRUE)
genome_len <- sum(as.numeric(seqlengths(BSgenome.Mmusculus.UCSC.mm10)[chr_order]))
repeat_pct <- repeat_len / genome_len

message("\n  Genome-wide repeat fractions:")
for (rclass in names(sort(repeat_pct, decreasing = TRUE))) {
    message(sprintf("    %s: %.2f%%", rclass, 100 * repeat_pct[rclass]))
}

# Log2 fold-enrichment (simple: observed / expected)
keep_repeats <- intersect(rmsk_cols, names(repeat_pct))
fe_raw_rmsk  <- as.matrix(cov_wide_rmsk[, ..keep_repeats])

log2fe_rmsk_simple <- log2(sweep(fe_raw_rmsk + eps, 2,
                                  repeat_pct[colnames(fe_raw_rmsk)] + eps, "/"))

rownames(log2fe_rmsk_simple) <- bin_ids[cov_wide_rmsk$bin]

# Pad missing bins with 0
missing_bins <- setdiff(bin_ids, rownames(log2fe_rmsk_simple))
if (length(missing_bins) > 0) {
    pad_rmsk <- matrix(0, nrow = length(missing_bins), ncol = ncol(log2fe_rmsk_simple),
                       dimnames = list(missing_bins, colnames(log2fe_rmsk_simple)))
    log2fe_rmsk_simple <- rbind(log2fe_rmsk_simple, pad_rmsk)
}
log2fe_rmsk_simple <- log2fe_rmsk_simple[bin_ids, , drop = FALSE]

saveRDS(log2fe_rmsk_simple, file.path(
    CONFIG$outdir,
    sprintf("log2fe_repeatmasker_simple.nonGap%02d.rds", CONFIG$non_gap_thresh * 100)
))

message(sprintf("  [OK] RepeatMasker enrichments: %d bins x %d repeat classes",
                nrow(log2fe_rmsk_simple), ncol(log2fe_rmsk_simple)))


# ==============================================================================
# STEP 7: CPG AND AT CONTENT
# ==============================================================================

banner("STEP 7: Computing sequence features")

bin_seqs   <- getSeq(BSgenome.Mmusculus.UCSC.mm10, bins_use)
cpg_count  <- vcountPattern("CG", bin_seqs, fixed = TRUE)
cpg_density <- cpg_count / width(bins_use)
at_count   <- letterFrequency(bin_seqs, letters = c("A", "T"))
at_content <- rowSums(at_count) / width(bins_use)

gc_dt <- data.frame(
    bin_id      = bins_use$bin_id,
    CpG_density = cpg_density,
    AT_content  = at_content,
    stringsAsFactors = FALSE
)
rownames(gc_dt) <- gc_dt$bin_id

saveRDS(gc_dt, file.path(CONFIG$outdir, "bins_CpG_AT.rds"))

message("  [OK] Sequence features computed")


# ==============================================================================
# STEP 8: K-MEANS CLUSTERING
# ==============================================================================

banner("STEP 8: Running k-means clustering")

km <- kmeans(comp_bin, centers = CONFIG$k_bins, nstart = 25)
groups_kept <- km$cluster
names(groups_kept) <- rownames(comp_bin)

cluster_assignments <- data.frame(
    bin_id      = rownames(comp_bin),
    ecg_cluster = groups_kept,
    stringsAsFactors = FALSE
)
fwrite(cluster_assignments,
       file.path(CONFIG$outdir, "kmeans_cluster_assignments.tsv"), sep = "\t")

message(sprintf("  [OK] %d k-means clusters created", CONFIG$k_bins))


# ==============================================================================
# STEP 9: BIN ORDERING AND ECG RELABELING
# ==============================================================================

banner("STEP 9: Ordering bins and relabeling ECG")

frac_active_all <- rowMeans(comp_bin == 1L, na.rm = TRUE)

ord_idx      <- order(groups_kept, -frac_active_all, na.last = TRUE)
ord_bins_ids <- rownames(comp_bin)[ord_idx]
ecg_split    <- factor(groups_kept[ord_bins_ids])
ecg_split    <- droplevels(ecg_split)

cluster_activity  <- tapply(frac_active_all[ord_bins_ids], ecg_split, mean)
cluster_order     <- order(cluster_activity, decreasing = TRUE)
ecg_levels_ordered <- levels(ecg_split)[cluster_order]
ecg_split <- factor(ecg_split, levels = ecg_levels_ordered)

ord_idx2     <- order(as.integer(ecg_split), -frac_active_all[ord_bins_ids])
ord_bins_ids <- ord_bins_ids[ord_idx2]
ecg_split    <- ecg_split[ord_idx2]

# Relabel consecutively
present_levels <- levels(ecg_split)
n_clusters     <- length(present_levels)
new_labels     <- as.character(1:n_clusters)
ecg_remap      <- setNames(new_labels, present_levels)
ecg_split      <- factor(ecg_remap[as.character(ecg_split)], levels = new_labels)

# Build rotated matrix
comp_ord    <- comp_bin[ord_bins_ids, , drop = FALSE]
mat_chr_rot <- t(apply(comp_ord, 2, function(x) ifelse(x == 1L, "1", "0")))
colnames(mat_chr_rot) <- ord_bins_ids

message(sprintf("  [OK] Bins ordered, ECG clusters relabeled (1-%d)", n_clusters))


# ==============================================================================
# STEP 10: BUILD ANNOTATIONS
# ==============================================================================

banner("STEP 10: Building annotations")

# Row annotations
annot_df <- data.frame(
    Sex      = sex_vec[colnames(comp_bin)],
    Age      = age_vec[colnames(comp_bin)],
    Celltype = ctype_vec[colnames(comp_bin)]
)
row_ha <- rowAnnotation(
    Age      = factor(annot_df$Age, levels = AGE_LEVELS),
    Sex      = factor(annot_df$Sex, levels = SEX_LEVELS),
    Celltype = factor(annot_df$Celltype, levels = names(CELLTYPE_COLS)),
    col = list(Age = AGE_COLS, Sex = SEX_COLS, Celltype = CELLTYPE_COLS)
)

# Column annotation tracks
age_track      <- per_bin_active_strict(comp_bin, age_vec, AGE_LEVELS, ord_bins_ids)
sex_track      <- per_bin_active_strict(comp_bin, sex_vec, SEX_LEVELS, ord_bins_ids)
celltype_track <- per_bin_active_strict(comp_bin, ctype_vec, names(CELLTYPE_COLS), ord_bins_ids)

age_mat_for_anno      <- t(age_track)
sex_mat_for_anno      <- t(sex_track)
celltype_mat_for_anno <- t(celltype_track)

ecg_blocks <- HeatmapAnnotation(
    ECG = anno_block(
        labels = new_labels,
        gp = gpar(fill = rep(c("#F2F2F2", "#FFFFFF"), length.out = n_clusters)),
        labels_gp = gpar(fontsize = 14, fontface = "bold"),
        labels_rot = 0
    ),
    which = "column",
    height = unit(6, "mm")
)

ha_age <- HeatmapAnnotation(
    Age = anno_barplot(
        age_mat_for_anno, which = "column", stack = TRUE,
        gp = gpar(fill = unname(AGE_COLS[colnames(age_mat_for_anno)]), col = NA),
        border = FALSE, axis = FALSE, ylim = c(0, 1), bar_width = 2,
        height = unit(22, "mm")
    ),
    which = "column"
)
ha_sex <- HeatmapAnnotation(
    Sex = anno_barplot(
        sex_mat_for_anno, which = "column", stack = TRUE,
        gp = gpar(fill = unname(SEX_COLS[colnames(sex_mat_for_anno)]), col = NA),
        border = FALSE, axis = FALSE, ylim = c(0, 1), bar_width = 2,
        height = unit(16, "mm")
    ),
    which = "column"
)
ha_cell <- HeatmapAnnotation(
    Celltype = anno_barplot(
        celltype_mat_for_anno, which = "column", stack = TRUE,
        gp = gpar(fill = unname(CELLTYPE_COLS[colnames(celltype_mat_for_anno)]), col = NA),
        border = FALSE, axis = FALSE, ylim = c(0, 1), bar_width = 2,
        height = unit(30, "mm")
    ),
    which = "column"
)
top_anno_all <- c(ecg_blocks, ha_age, ha_sex, ha_cell)

# Chromosome track
chr_vec <- bins_df$chr[match(ord_bins_ids, bins_df$bin_id)]
chr_mat <- matrix(chr_vec, nrow = 1, dimnames = list("Chromosome", ord_bins_ids))
ht_chr  <- Heatmap(
    chr_mat, name = "Chromosome", col = CHR_COLS,
    cluster_rows = FALSE, cluster_columns = FALSE,
    show_row_names = TRUE, show_column_names = FALSE,
    row_names_gp = gpar(fontsize = 11, fontface = "bold"),
    height = unit(0.8, "cm")
)

# CpG and AT tracks
cpg_vec <- gc_dt[ord_bins_ids, "CpG_density", drop = TRUE]
at_vec  <- gc_dt[ord_bins_ids, "AT_content", drop = TRUE]
cpg_brk <- quantile(cpg_vec, c(0.05, 0.50, 0.95), na.rm = TRUE)
at_brk  <- quantile(at_vec, c(0.05, 0.50, 0.95), na.rm = TRUE)

cpg_mat <- matrix(cpg_vec, nrow = 1, dimnames = list("CpG_density", ord_bins_ids))
at_mat  <- matrix(at_vec, nrow = 1, dimnames = list("AT_content", ord_bins_ids))

ht_cpg <- Heatmap(
    cpg_mat, name = "CpG density",
    col = circlize::colorRamp2(cpg_brk, c("white", "lightpink", "red")),
    cluster_rows = FALSE, cluster_columns = FALSE,
    show_row_names = TRUE, na_col = "white",
    row_names_gp = gpar(fontsize = 11, fontface = "bold"),
    height = unit(0.8, "cm")
)
ht_at <- Heatmap(
    at_mat, name = "AT content",
    col = circlize::colorRamp2(at_brk, c("white", "lightblue", "purple")),
    cluster_rows = FALSE, cluster_columns = FALSE,
    show_row_names = TRUE, na_col = "white",
    row_names_gp = gpar(fontsize = 11, fontface = "bold"),
    height = unit(0.8, "cm")
)

message("  [OK] All annotations built")


# ==============================================================================
# STEP 11: CHROMHMM HEATMAP PANEL
# ==============================================================================

banner("STEP 11: Building ChromHMM heatmap panel")

# Pad missing bins
missing_cols <- setdiff(ord_bins_ids, rownames(log2fe_bin))
if (length(missing_cols)) {
    pad <- matrix(NA_real_, nrow = length(missing_cols), ncol = ncol(log2fe_bin),
                  dimnames = list(missing_cols, colnames(log2fe_bin)))
    log2fe_bin <- rbind(log2fe_bin, pad)
}
log2fe_bin <- log2fe_bin[ord_bins_ids, , drop = FALSE]

state_mat   <- t(log2fe_bin)
state_mat_z <- t(scale(t(state_mat)))
state_mat_z[!is.finite(state_mat_z)] <- NA_real_

state_order <- c(
    "mTSS", "mTx", "mTxEx", "mTxWk", "mTxEnh",
    "mEnhA", "mEnhWk", "mBivProm",
    "mPromF", "mOpenC",
    "mReprPC", "mReprPC_openC",
    "mQuies", "mHET", "mZNF"
)
state_order <- intersect(state_order, rownames(state_mat_z))
state_mat_z <- state_mat_z[state_order, , drop = FALSE]

cap_z <- 2.0
state_mat_z_clip <- pmax(pmin(state_mat_z, cap_z), -cap_z)
col_fun_chromhmm_z <- circlize::colorRamp2(
    c(-cap_z, 0, cap_z),
    c("#2166ac", "white", "#b2182b")
)

ht_states <- Heatmap(
    state_mat_z_clip,
    name = "ChromHMM z",
    col  = col_fun_chromhmm_z,
    cluster_rows = FALSE, cluster_columns = FALSE,
    show_row_names = TRUE, show_column_names = FALSE,
    row_names_gp = gpar(fontsize = 14, fontface = "bold"),
    row_title = "ChromHMM",
    row_title_gp = gpar(fontsize = 14, fontface = "bold"),
    height = unit(6, "cm")
)

message("  [OK] ChromHMM panel built")


# ==============================================================================
# STEP 12: REPEATMASKER HEATMAP PANEL
# ==============================================================================

banner("STEP 12: Building RepeatMasker heatmap panel")

log2fe_rmsk_ord <- log2fe_rmsk_simple[ord_bins_ids, , drop = FALSE]
rmsk_mat        <- t(log2fe_rmsk_ord)
rmsk_mat_z      <- t(scale(t(rmsk_mat)))
rmsk_mat_z[!is.finite(rmsk_mat_z)] <- NA_real_

repeat_order <- intersect(REPEAT_CLASSES_KEEP, rownames(rmsk_mat_z))
rmsk_mat_z   <- rmsk_mat_z[repeat_order, , drop = FALSE]

cap_z_rmsk       <- 2.0
rmsk_mat_z_clip  <- pmax(pmin(rmsk_mat_z, cap_z_rmsk), -cap_z_rmsk)

# Teal-White-Coral (distinct from ChromHMM blue-white-red)
col_fun_rmsk_z <- circlize::colorRamp2(
    c(-cap_z_rmsk, 0, cap_z_rmsk),
    c("#008080", "white", "#FF7F50")
)

ht_repeats <- Heatmap(
    rmsk_mat_z_clip,
    name = "Repeat z",
    col  = col_fun_rmsk_z,
    cluster_rows = FALSE, cluster_columns = FALSE,
    show_row_names = TRUE, show_column_names = FALSE,
    row_names_gp = gpar(fontsize = 14, fontface = "bold"),
    row_title = "Repeats",
    row_title_gp = gpar(fontsize = 14, fontface = "bold"),
    row_names_side = "left",
    height = unit(2.5, "cm")
)

message(sprintf("  [OK] RepeatMasker panel: %d classes x %d bins",
                nrow(rmsk_mat_z_clip), ncol(rmsk_mat_z_clip)))
message(sprintf("  [OK] Order: %s",
                paste(rownames(rmsk_mat_z_clip), collapse = " -> ")))


# ==============================================================================
# ROW ORDERING (within-age k-means)
# ==============================================================================

comp_ord_numeric <- comp_bin[ord_bins_ids, , drop = FALSE]
mat_rowsamples   <- t(comp_ord_numeric)
col_w <- rep(1, ncol(mat_rowsamples))
mat_rowsamples_w <- sweep(mat_rowsamples, 2, col_w, "*")

row_age <- factor(as.character(age_vec[rownames(mat_chr_rot)]), levels = AGE_LEVELS)

k_max_per_age <- 6
order_rows_in_slice <- function(rids) {
    if (length(rids) <= 2) return(rids)
    X  <- mat_rowsamples_w[rids, , drop = FALSE]
    k  <- min(k_max_per_age, max(2, floor(sqrt(nrow(X)))))
    km_inner <- kmeans(X, centers = k, nstart = 25)
    cl      <- km_inner$cluster
    centers <- km_inner$centers
    sw      <- col_w / sum(col_w)
    cl_score <- as.numeric(centers %*% sw)
    cl_order <- order(cl_score, decreasing = TRUE)
    out <- character(0)
    for (cid in cl_order) {
        members <- names(cl)[cl == cid]
        if (length(members) == 1) { out <- c(out, members); next }
        d2 <- as.numeric(rowSums((X[members, , drop = FALSE] - centers[cid, ])^2))
        out <- c(out, members[order(d2, decreasing = TRUE)])
    }
    out
}

row_order_ids <- unlist(lapply(levels(row_age), function(ag) {
    ids <- rownames(mat_chr_rot)[row_age == ag]
    if (!length(ids)) return(character(0))
    order_rows_in_slice(ids)
}), use.names = FALSE)

if (length(row_order_ids) < nrow(mat_chr_rot)) {
    leftovers <- setdiff(rownames(mat_chr_rot), row_order_ids)
    row_order_ids <- c(row_order_ids, leftovers)
}

row_order_idx <- match(row_order_ids, rownames(mat_chr_rot))
row_split_age <- row_age


# ==============================================================================
# STEP 13: ASSEMBLE AND SAVE HEATMAP
# ==============================================================================

banner("STEP 13: Creating and saving heatmap")

ComplexHeatmap::ht_opt(
    legend_title_gp  = gpar(fontsize = 20, fontface = "bold"),
    legend_labels_gp = gpar(fontsize = 16)
)

ht_main <- Heatmap(
    mat_chr_rot,
    name = "Compartment",
    col  = COMPARTMENT_COLS,
    row_split        = row_split_age,
    row_order        = row_order_idx,
    cluster_rows     = FALSE,
    cluster_columns  = FALSE,
    show_row_names   = FALSE,
    left_annotation  = row_ha,
    row_title        = "Groups (Age fixed; k-means within Age)",
    use_raster       = FALSE,
    na_col           = "white",
    heatmap_legend_param = list(
        at = c("0", "1"),
        labels = c("Repressive", "Active"),
        labels_gp = gpar(fontsize = 14)
    ),
    column_split     = ecg_split,
    column_gap       = unit(0, "mm"),
    gap              = unit(0.5, "mm"),
    show_column_dend = FALSE,
    show_column_names = FALSE,
    top_annotation   = top_anno_all
)

# With RepeatMasker
ht_stack <- ht_main %v% ht_chr %v% ht_cpg %v% ht_at %v% ht_states %v% ht_repeats

message("  Saving heatmap with RepeatMasker...")
pdf(file.path(CONFIG$outdir, "compartment_heatmap_full_with_repeats.pdf"),
    width = 24, height = 18)
draw(ht_stack,
     heatmap_legend_side = "right",
     annotation_legend_side = "right",
     merge_legend = TRUE)
dev.off()

png(file.path(CONFIG$outdir, "compartment_heatmap_full_with_repeats.png"),
    width = 4800, height = 3600, res = 200)
draw(ht_stack,
     heatmap_legend_side = "right",
     annotation_legend_side = "right",
     merge_legend = TRUE)
dev.off()

# Without RepeatMasker
ht_stack_no_rmsk <- ht_main %v% ht_chr %v% ht_cpg %v% ht_at %v% ht_states

message("  Saving heatmap without RepeatMasker...")
pdf(file.path(CONFIG$outdir, "compartment_heatmap_full.pdf"),
    width = 24, height = 16)
draw(ht_stack_no_rmsk,
     heatmap_legend_side = "right",
     annotation_legend_side = "right",
     merge_legend = TRUE)
dev.off()

png(file.path(CONFIG$outdir, "compartment_heatmap_full.png"),
    width = 4800, height = 3200, res = 200)
draw(ht_stack_no_rmsk,
     heatmap_legend_side = "right",
     annotation_legend_side = "right",
     merge_legend = TRUE)
dev.off()

message("  [OK] Both heatmap versions saved")


# ==============================================================================
# STEP 14: REPEAT ENRICHMENT SUMMARY
# ==============================================================================

banner("STEP 14: Repeat enrichment by compartment state")

active_bins     <- names(which(rowMeans(comp_bin == 1L, na.rm = TRUE) >= 0.8))
repressive_bins <- names(which(rowMeans(comp_bin == 0L, na.rm = TRUE) >= 0.8))

if (length(active_bins) > 0 && length(repressive_bins) > 0) {

    mean_active     <- colMeans(log2fe_rmsk_simple[active_bins, , drop = FALSE], na.rm = TRUE)
    mean_repressive <- colMeans(log2fe_rmsk_simple[repressive_bins, , drop = FALSE], na.rm = TRUE)

    repeat_summary <- data.table(
        Repeat_Class           = names(mean_active),
        Mean_log2FE_Active     = mean_active,
        Mean_log2FE_Repressive = mean_repressive,
        Difference             = mean_active - mean_repressive
    )
    repeat_summary[, Repeat_Class := factor(Repeat_Class, levels = REPEAT_CLASSES_KEEP)]
    repeat_summary <- repeat_summary[order(Repeat_Class)]

    message("\n  Repeat enrichment by compartment:")
    message("  (Positive difference = enriched in Active compartment)")
    print(repeat_summary)

    fwrite(repeat_summary,
           file.path(CONFIG$outdir, "repeat_enrichment_by_compartment.tsv"), sep = "\t")

    # Barplot
    repeat_summary_long <- melt(
        repeat_summary[, .(Repeat_Class, Mean_log2FE_Active, Mean_log2FE_Repressive)],
        id.vars = "Repeat_Class",
        variable.name = "Compartment",
        value.name = "Mean_log2FE"
    )
    repeat_summary_long[, Compartment := gsub("Mean_log2FE_", "", Compartment)]
    repeat_summary_long[, Repeat_Class := factor(Repeat_Class, levels = REPEAT_CLASSES_KEEP)]

    p_repeat <- ggplot(repeat_summary_long,
                       aes(x = Repeat_Class, y = Mean_log2FE, fill = Compartment)) +
        geom_col(position = "dodge", width = 0.7) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "gray40", linewidth = 1) +
        scale_fill_manual(values = c("Active" = "#FF7F50", "Repressive" = "#008080")) +
        labs(
            title = "Repeat Element Enrichment by Compartment",
            subtitle = "Mean log2 fold-enrichment in stable Active vs Repressive bins",
            x = "Repeat Class",
            y = "Mean log2(Observed/Expected)",
            fill = "Compartment"
        ) +
        theme_bw(base_size = 14) +
        theme(
            text = element_text(family = "Arial"),
            axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 12),
            plot.title = element_text(face = "bold", size = 16),
            plot.subtitle = element_text(size = 11, color = "gray30"),
            legend.position = "top"
        )

    ggsave(file.path(CONFIG$outdir, "Barplot_Repeat_Enrichment_by_Compartment.pdf"),
           p_repeat, width = 10, height = 7, useDingbats = FALSE)

    message("  [OK] Repeat enrichment summary saved")
}


# ==============================================================================
# STEP 15: STABILITY ANALYSIS (5-CLASS MODEL)
# ==============================================================================

banner("STEP 15: Stability analysis")

check_monotonicity <- function(states) {
    states <- states[!is.na(states)]
    if (length(states) < 2) {
        return(list(is_monotonic = NA, direction = "Insufficient", n_transitions = NA))
    }
    transitions <- states[-1] != states[-length(states)]
    n_trans     <- sum(transitions)
    if (n_trans == 0) {
        return(list(is_monotonic = TRUE, direction = "Stable", n_transitions = 0))
    }
    numeric_states <- ifelse(states == "Active", 1L, 0L)
    is_monotonic_AR <- all(diff(numeric_states) <= 0) & any(diff(numeric_states) < 0)
    is_monotonic_RA <- all(diff(numeric_states) >= 0) & any(diff(numeric_states) > 0)
    if (is_monotonic_AR) {
        return(list(is_monotonic = TRUE, direction = "A_to_R", n_transitions = n_trans))
    } else if (is_monotonic_RA) {
        return(list(is_monotonic = TRUE, direction = "R_to_A", n_transitions = n_trans))
    } else {
        return(list(is_monotonic = FALSE, direction = "Non_Monotonic", n_transitions = n_trans))
    }
}

classify_stability_extended <- function(fracs, t_on = 0.6, t_off = 0.4) {
    if (all(is.na(fracs))) return("Missing")
    n_valid <- sum(!is.na(fracs))
    if (n_valid < 3) return("Insufficient_Data")

    states <- rep("Intermediate", length(fracs))
    states[fracs >= t_on]  <- "Active"
    states[fracs <= t_off] <- "Repressive"
    states[is.na(fracs)]   <- NA

    valid_states <- states[states %in% c("Active", "Repressive")]
    if (length(valid_states) < 3) return("Insufficient_Data")

    n_active     <- sum(valid_states == "Active")
    n_repressive <- sum(valid_states == "Repressive")
    n_total      <- length(valid_states)

    if (n_active / n_total >= 0.8) return("Stable_Active")
    if (n_repressive / n_total >= 0.8) return("Stable_Repressive")

    mono_check <- check_monotonicity(valid_states)
    if (is.na(mono_check$is_monotonic)) return("Insufficient_Data")
    if (mono_check$is_monotonic) {
        if (mono_check$direction == "A_to_R") return("Monotonic_A_to_R")
        if (mono_check$direction == "R_to_A") return("Monotonic_R_to_A")
        return("Stable_Mixed")
    }
    return("Non_Monotonic")
}

compute_stability_extended <- function(comp_bin, age_vec, sex_vec, ctype_vec,
                                        age_levels, t_on = 0.6, t_off = 0.4) {
    res      <- list()
    sexes    <- intersect(c("male", "female"), unique(as.character(sex_vec)))
    celltypes <- sort(unique(as.character(ctype_vec)))

    for (sx in sexes) {
        use_sx <- names(sex_vec)[as.character(sex_vec) == sx]
        for (ct in celltypes) {
            use <- intersect(use_sx, names(ctype_vec)[as.character(ctype_vec) == ct])
            if (length(use) < 2) next
            message(sprintf("  [STABILITY] %s - %s (%d samples)", sx, ct, length(use)))
            M <- comp_bin[, use, drop = FALSE]
            F_mat <- frac_active_by(M, age_vec[use], age_levels)
            if (!is.matrix(F_mat) || nrow(F_mat) == 0 || ncol(F_mat) < 2) next
            stability <- apply(F_mat, 1, classify_stability_extended,
                               t_on = t_on, t_off = t_off)
            res[[paste(sx, ct, sep = "|")]] <- data.table(
                bin_id    = rownames(M),
                Sex       = sx,
                Celltype  = ct,
                Stability = stability
            )
        }
    }
    return(rbindlist(res, use.names = TRUE, fill = TRUE))
}

stability_ext_dt <- compute_stability_extended(
    comp_bin, age_vec, sex_vec, ctype_vec, AGE_LEVELS,
    t_on = CONFIG$t_on, t_off = CONFIG$t_off
)

fwrite(stability_ext_dt,
       file.path(CONFIG$outdir, "compartment_stability_5class.tsv"), sep = "\t")

message("\n  Stability distribution:")
print(stability_ext_dt[, .N, by = Stability][order(-N)])


# ==============================================================================
# STEP 16: REPEAT ENRICHMENT BY STABILITY CLASS
# ==============================================================================

banner("STEP 16: Repeat enrichment by stability class")

stability_with_repeats <- merge(
    stability_ext_dt,
    data.table(bin_id = rownames(log2fe_rmsk_simple), log2fe_rmsk_simple),
    by = "bin_id",
    all.x = TRUE
)

repeat_cols_use <- colnames(log2fe_rmsk_simple)
stability_repeat_summary <- stability_with_repeats[
    !Stability %in% c("Missing", "Insufficient_Data"),
    lapply(.SD, mean, na.rm = TRUE),
    by = .(Stability, Sex),
    .SDcols = repeat_cols_use
]

stability_repeat_long <- melt(
    stability_repeat_summary,
    id.vars = c("Stability", "Sex"),
    variable.name = "Repeat_Class",
    value.name = "Mean_log2FE"
)
stability_repeat_long[, Stability := factor(Stability, levels = STABILITY_LEVELS)]
stability_repeat_long[, Repeat_Class := factor(Repeat_Class, levels = repeat_order)]

for (sx in c("male", "female")) {

    plot_dt <- stability_repeat_long[Sex == sx]
    if (nrow(plot_dt) == 0) next

    mat_stab_rep   <- dcast(plot_dt, Repeat_Class ~ Stability, value.var = "Mean_log2FE")
    mat_stab_rep_m <- as.matrix(mat_stab_rep[, -1, with = FALSE])
    rownames(mat_stab_rep_m) <- mat_stab_rep$Repeat_Class

    mat_stab_rep_m <- mat_stab_rep_m[
        intersect(REPEAT_CLASSES_KEEP, rownames(mat_stab_rep_m)), , drop = FALSE
    ]
    mat_stab_rep_m <- mat_stab_rep_m[
        , intersect(STABILITY_LEVELS, colnames(mat_stab_rep_m)), drop = FALSE
    ]

    cap_val     <- 1.5
    mat_clipped <- pmax(pmin(mat_stab_rep_m, cap_val), -cap_val)

    col_fun_stab <- circlize::colorRamp2(
        c(-cap_val, 0, cap_val),
        c("#008080", "white", "#FF7F50")
    )

    ht_stab_repeat <- Heatmap(
        mat_clipped,
        name = "Mean log2FE",
        col  = col_fun_stab,
        cluster_rows = FALSE, cluster_columns = FALSE,
        show_row_names = TRUE, show_column_names = TRUE,
        column_names_rot = 45,
        row_names_gp    = gpar(fontsize = 12, fontface = "bold"),
        column_names_gp = gpar(fontsize = 10),
        cell_fun = function(j, i, x, y, width, height, fill) {
            grid.text(sprintf("%.2f", mat_stab_rep_m[i, j]), x, y,
                      gp = gpar(fontsize = 9, fontface = "bold"))
        },
        column_title    = sprintf("Repeat Enrichment by Stability Class - %s",
                                  tools::toTitleCase(sx)),
        column_title_gp = gpar(fontsize = 14, fontface = "bold"),
        row_title       = "Euchromatin -> Heterochromatin",
        row_title_gp    = gpar(fontsize = 10, fontface = "italic")
    )

    pdf(file.path(CONFIG$outdir,
                  sprintf("Heatmap_Repeat_by_Stability_%s.pdf", tools::toTitleCase(sx))),
        width = 10, height = 8)
    draw(ht_stab_repeat)
    dev.off()

    message(sprintf("  [OK] Saved: Heatmap_Repeat_by_Stability_%s.pdf",
                    tools::toTitleCase(sx)))
}

fwrite(stability_repeat_summary,
       file.path(CONFIG$outdir, "repeat_enrichment_by_stability_class.tsv"), sep = "\t")
message("  [OK] Summary table saved")


# ==============================================================================
# STEP 17: STABILITY BARPLOTS AND ALLUVIAL PLOTS
# ==============================================================================

banner("STEP 17: Stability barplots and alluvial plots")

# --- State-per-age data ---
state_per_age_list <- list()
for (sample_name in colnames(comp_bin)) {
    state_per_age_list[[sample_name]] <- data.table(
        bin_id   = rownames(comp_bin),
        State    = ifelse(comp_bin[, sample_name] == 1L, "Active", "Repressive"),
        Age      = age_vec[sample_name],
        Sex      = sex_vec[sample_name],
        Celltype = ctype_vec[sample_name]
    )
}
state_per_age <- rbindlist(state_per_age_list)
state_per_age[, Age := factor(Age, levels = AGE_LEVELS)]

state_with_stability <- merge(
    state_per_age,
    stability_ext_dt[, .(bin_id, Sex, Celltype, Stability)],
    by = c("bin_id", "Sex", "Celltype"),
    all.x = TRUE
)

message(sprintf("  [OK] State data: %d rows", nrow(state_with_stability)))

# --- Stability barplot ---
plot_dt <- stability_ext_dt[!Stability %in% c("Missing", "Insufficient_Data")]

p_5class <- ggplot(
    plot_dt[, .N, by = .(Sex, Celltype, Stability)],
    aes(x = Celltype, y = N, fill = Stability)
) +
    geom_col(position = "stack") +
    facet_wrap(~ Sex, ncol = 1) +
    scale_fill_manual(values = STABILITY_COLS) +
    labs(
        title = "Compartment Stability - 5-Class Model",
        y = "Number of 80 kb bins",
        x = "Cell Type"
    ) +
    theme_bw(base_size = 12) +
    theme(
        text = element_text(family = "Arial"),
        axis.text.x = element_text(angle = 45, hjust = 1)
    )

ggsave(file.path(CONFIG$outdir, "Barplot_Stability_5class.pdf"),
       p_5class, width = 12, height = 8, useDingbats = FALSE)

message("  [OK] Stability barplot saved")

# --- Alluvial plots ---
message("  Generating alluvial plots...")

switching_data <- state_with_stability[Stability %in% SWITCHING_CLASSES]

# Variant 1: All switching bins, faceted by sex
if (nrow(switching_data) > 0) {

    alluvial_summary <- switching_data[
        , .N, by = .(Age, State, Stability, Sex)
    ]

    p_alluvial <- ggplot(
        alluvial_summary,
        aes(x = Age, stratum = State, alluvium = interaction(Stability, State),
            y = N, fill = State)
    ) +
        geom_flow(alpha = 0.4) +
        geom_stratum(width = 0.3) +
        scale_fill_manual(values = STATE_COLS) +
        facet_grid(Stability ~ Sex, scales = "free_y") +
        labs(
            title = "Compartment Switching Across Aging",
            y = "Number of 80 kb bins",
            x = "Age Group"
        ) +
        theme_bw(base_size = 12) +
        theme(
            text = element_text(family = "Arial"),
            axis.text.x = element_text(angle = 45, hjust = 1),
            strip.text = element_text(face = "bold")
        )

    ggsave(file.path(CONFIG$outdir, "Alluvial_switching_by_sex.pdf"),
           p_alluvial, width = 14, height = 10, useDingbats = FALSE)

    message("  [OK] Alluvial plot saved")

    # Variant 2: Per cell type
    for (ct in unique(switching_data$Celltype)) {
        ct_data <- switching_data[Celltype == ct]
        if (nrow(ct_data) == 0) next

        ct_summary <- ct_data[, .N, by = .(Age, State, Stability, Sex)]

        p_ct <- ggplot(
            ct_summary,
            aes(x = Age, stratum = State, alluvium = interaction(Stability, State),
                y = N, fill = State)
        ) +
            geom_flow(alpha = 0.4) +
            geom_stratum(width = 0.3) +
            scale_fill_manual(values = STATE_COLS) +
            facet_grid(Stability ~ Sex, scales = "free_y") +
            labs(
                title = sprintf("Compartment Switching - %s", ct),
                y = "Number of 80 kb bins",
                x = "Age Group"
            ) +
            theme_bw(base_size = 12) +
            theme(
                text = element_text(family = "Arial"),
                axis.text.x = element_text(angle = 45, hjust = 1),
                strip.text = element_text(face = "bold")
            )

        safe_ct <- clean_name(ct)
        ggsave(file.path(CONFIG$outdir, sprintf("Alluvial_switching_%s.pdf", safe_ct)),
               p_ct, width = 14, height = 10, useDingbats = FALSE)
    }

    message(sprintf("  [OK] Per-celltype alluvial plots saved (%d cell types)",
                    length(unique(switching_data$Celltype))))

    # Variant 3: Switching class colored (not state colored)
    switch_summary <- switching_data[, .N, by = .(Age, Stability, Sex)]

    p_switch_class <- ggplot(
        switch_summary,
        aes(x = Age, stratum = Stability, alluvium = Stability,
            y = N, fill = Stability)
    ) +
        geom_flow(alpha = 0.4) +
        geom_stratum(width = 0.3) +
        scale_fill_manual(values = SWITCHING_COLS) +
        facet_wrap(~ Sex, ncol = 2) +
        labs(
            title = "Switching Class Distribution Across Aging",
            y = "Number of 80 kb bins",
            x = "Age Group"
        ) +
        theme_bw(base_size = 12) +
        theme(
            text = element_text(family = "Arial"),
            axis.text.x = element_text(angle = 45, hjust = 1)
        )

    ggsave(file.path(CONFIG$outdir, "Alluvial_switching_class_by_sex.pdf"),
           p_switch_class, width = 12, height = 6, useDingbats = FALSE)

    message("  [OK] Switching class alluvial saved")
}


# ==============================================================================
# STEP 18: SESSION INFO
# ==============================================================================

banner("STEP 18: Session info")

writeLines(
    capture.output(sessionInfo()),
    file.path(CONFIG$outdir, "sessionInfo.txt")
)

message("  [OK] sessionInfo.txt saved")


# ==============================================================================
# FINAL SUMMARY
# ==============================================================================

end_time <- Sys.time()
duration <- difftime(end_time, start_time, units = "mins")

cat("\n")
cat("========================================================================\n")
cat("  PIPELINE COMPLETE\n")
cat("========================================================================\n")
cat(sprintf("  Duration: %.1f minutes\n", as.numeric(duration)))
cat(sprintf("  Output:   %s\n", CONFIG$outdir))
cat("\n")
cat("  Heatmaps:\n")
cat("    compartment_heatmap_full_with_repeats.pdf/png\n")
cat("    compartment_heatmap_full.pdf/png\n")
cat("\n")
cat("  Core data:\n")
cat("    compartments_binary_A1_R0.rds\n")
cat("    compartment_stability_5class.tsv\n")
cat("    kmeans_cluster_assignments.tsv\n")
cat("    matrix_signal_z.rds\n")
cat("    bins_CpG_AT.rds\n")
cat("\n")
cat("  Enrichment:\n")
cat(sprintf("    log2fe_bin.nonGap%02d.rds (ChromHMM)\n", CONFIG$non_gap_thresh * 100))
cat(sprintf("    log2fe_repeatmasker_simple.nonGap%02d.rds\n", CONFIG$non_gap_thresh * 100))
cat("    repeat_enrichment_by_compartment.tsv\n")
cat("    repeat_enrichment_by_stability_class.tsv\n")
cat("\n")
cat("  Plots:\n")
cat("    Barplot_Stability_5class.pdf\n")
cat("    Barplot_Repeat_Enrichment_by_Compartment.pdf\n")
cat("    Heatmap_Repeat_by_Stability_Male/Female.pdf\n")
cat("    Alluvial_switching_by_sex.pdf\n")
cat("    Alluvial_switching_[celltype].pdf\n")
cat("    Alluvial_switching_class_by_sex.pdf\n")
cat("\n")
cat("  Repeat classes (Euchromatin -> Heterochromatin):\n")
cat("    SINE -> DNA_Transposon -> LTR -> LINE\n")
cat("\n")
cat("  Reproducibility:\n")
cat("    sessionInfo.txt\n")
cat("    manifest.normalized.tsv\n")
cat("========================================================================\n")
