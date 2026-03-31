#!/usr/bin/env Rscript
## =========================================================
## COMPARTMENT ANALYSIS - MASTER RUN SCRIPT
## 
## This script runs the full compartment analysis pipeline
## with your specific paths and parameters.
##
## Usage:
##   Rscript run_compartment_analysis.R
##
## Or in R:
##   source("run_compartment_analysis.R")
##
## =========================================================

cat("\n")
cat("╔══════════════════════════════════════════════════════════════════╗\n")
cat("║       COMPARTMENT ANALYSIS PIPELINE - FULL RUN                   ║\n")
cat("║       Mouse Liver Aging Multiome Study                           ║\n")
cat("╚══════════════════════════════════════════════════════════════════╝\n\n")

## =========================================================
## USER-SPECIFIC PATHS AND PARAMETERS
## =========================================================

# Project paths
proj_path      <- "/data/sarkern2/multiome_liver/Seurat/archR/ArchR_Projects/Step6_Xwnn_UMAP"
group_dir      <- file.path(proj_path, "GroupBigWigs", "AgeSexCelltype__tile80k")
manifest_path  <- file.path(group_dir, "AgeSexCelltype.bigwig_manifest_tile80k.tsv")
outdir         <- "/data/sarkern2/multiome_liver/Seurat/epigenome/compartment_main_heatmap"

# Reference files
mm10_blacklist_bed <- "/data/sarkern2/multiome_liver/Seurat/archR/mm10-blacklist.v2.bed.gz"
fs_bed             <- file.path(proj_path, "mm10_100_segments_segments.bed.gz")

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


## =========================================================
## OVERRIDE CONFIG WITH USER PATHS
## =========================================================

# This will be used by the config script
CONFIG <- list(
    proj_path       = proj_path,
    group_dir       = group_dir,
    manifest_path   = manifest_path,
    outdir          = outdir,
    mm10_blacklist  = mm10_blacklist_bed,
    chromhmm_bed    = fs_bed,
    tile_size       = tile_size,
    non_gap_thresh  = non_gap_thresh,
    k_bins          = k_bins,
    threads         = threads,
    t_on            = t_on,
    t_off           = t_off,
    seed            = seed
)


## =========================================================
## VALIDATE PATHS
## =========================================================

message("[SETUP] Validating paths...")

paths_to_check <- list(
    "Project directory" = proj_path,
    "BigWig directory"  = group_dir,
    "Manifest file"     = manifest_path,
    "Blacklist BED"     = mm10_blacklist_bed,
    "ChromHMM BED"      = fs_bed
)

all_exist <- TRUE
for (name in names(paths_to_check)) {
    path <- paths_to_check[[name]]
    exists <- file.exists(path)
    status <- if (exists) "✅" else "❌"
    message(sprintf("  %s %s: %s", status, name, path))
    if (!exists) all_exist <- FALSE
}

if (!all_exist) {
    stop("\n[ERROR] Some required paths do not exist. Please check paths above.")
}

# Create output directory
if (!dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE)
    message(sprintf("  📁 Created output directory: %s", outdir))
} else {
    message(sprintf("  📁 Output directory exists: %s", outdir))
}


## =========================================================
## LOAD REQUIRED PACKAGES
## =========================================================

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
        message(sprintf("  • %s", pkg))
    }
    message("\nInstall with:")
    message('  BiocManager::install(c("', paste(missing_packages, collapse = '", "'), '"))')
    stop("\n[ERROR] Please install missing packages before running.")
}

message("  ✅ All required packages available")


## =========================================================
## DEFINE SHARED CONSTANTS (FROM CONFIG)
## =========================================================

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
NONMONO_PATTERNS <- c("Transient_Return", "Two_Switch", "Triple_Switch", "Highly_Dynamic", "Other")

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

chr_order <- paste0("chr", c(1:19, "X", "Y"))
n_chr <- length(chr_order)
hues <- seq(0, 360 - 360/n_chr, length.out = n_chr)
CHR_COLS <- setNames(grDevices::hcl(h = hues, c = 90, l = 60), chr_order)


## =========================================================
## HELPER FUNCTIONS
## =========================================================

banner <- function(text, char = "=", width = 70) {
    line <- paste(rep(char, width), collapse = "")
    message("\n", line)
    message(text)
    message(line, "\n")
}

box_banner <- function(text) {
    cat("\n")
    cat("╔══════════════════════════════════════════════════════════════════╗\n")
    cat(sprintf("║ %-66s ║\n", text))
    cat("╚══════════════════════════════════════════════════════════════════╝\n\n")
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
    message(sprintf("    → %s (%d regions)", basename(filename), nrow(bed)))
}

check_required_objects <- function(required_objects, stop_on_missing = TRUE) {
    missing <- required_objects[!sapply(required_objects, exists, envir = .GlobalEnv)]
    if (length(missing) > 0) {
        msg <- sprintf("Missing required objects: %s", paste(missing, collapse = ", "))
        if (stop_on_missing) stop(msg) else warning(msg)
        return(missing)
    }
    return(character(0))
}

message("  ✅ Helper functions defined")


## =========================================================
## SET SEED AND THREADS
## =========================================================

set.seed(CONFIG$seed)
message(sprintf("\n[CONFIG] Random seed: %d", CONFIG$seed))
message(sprintf("[CONFIG] Threads: %d", CONFIG$threads))


## =========================================================
## LOAD LIBRARIES
## =========================================================

message("\n[SETUP] Loading libraries...")

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

message("  ✅ All libraries loaded")


## =========================================================
## PIPELINE EXECUTION
## =========================================================

# Track timing
start_time <- Sys.time()


## -----------------------------------------------------------
## STEP 1: LOAD MANIFEST
## -----------------------------------------------------------

box_banner("STEP 1/15: Loading manifest")

man <- fread(CONFIG$manifest_path)
setnames(man, tolower(names(man)))

if (!all(c("age","sex","celltype") %in% names(man))) {
    clean <- sub("-TileSize.*$", "", man$group)
    parts <- strsplit(clean, "__", fixed = TRUE)
    man[, age      := vapply(parts, `[`, character(1), 1)]
    man[, sex      := vapply(parts, `[`, character(1), 2)]
    man[, celltype := vapply(parts, `[`, character(1), 3)]
}

man[, sex := {s <- tolower(as.character(sex)); s[!s %in% c("male","female")] <- "unknown"; s}]
man[, age := factor(as.character(age), levels = AGE_LEVELS)]
man[, celltype := sub("-TileSize.*$", "", celltype)]

fwrite(man, file.path(CONFIG$outdir, "manifest.normalized.tsv"), sep = "\t")

bw_vec    <- setNames(man$bigwig,   man$group)
age_vec   <- setNames(as.character(man$age),      man$group)
sex_vec   <- setNames(as.character(man$sex),      man$group)
ctype_vec <- setNames(as.character(man$celltype), man$group)

message(sprintf("  ✅ %d samples loaded", nrow(man)))
message(sprintf("  ✅ Ages: %s", paste(unique(age_vec), collapse = ", ")))
message(sprintf("  ✅ Sexes: %s", paste(unique(sex_vec), collapse = ", ")))
message(sprintf("  ✅ Celltypes: %d unique", length(unique(ctype_vec))))


## -----------------------------------------------------------
## STEP 2: CREATE GENOMIC BINS
## -----------------------------------------------------------

box_banner("STEP 2/15: Creating genomic bins")

seqs <- seqlengths(BSgenome.Mmusculus.UCSC.mm10)
seqs <- seqs[grepl("^chr([0-9]{1,2}|X|Y)$", names(seqs))]
bins100k <- tileGenome(seqs, tilewidth = CONFIG$tile_size, cut.last.tile.in.chrom = TRUE)

mm10_blacklist <- import(CONFIG$mm10_blacklist)
bins_clean <- subsetByOverlaps(bins100k, mm10_blacklist, invert = TRUE)

bed_dt <- tryCatch(fread(CONFIG$chromhmm_bed, header = FALSE),
                   error = function(e) fread(cmd = paste("zcat -f", shQuote(CONFIG$chromhmm_bed))))
setnames(bed_dt, c("chr","start0","end","label"))
bed_dt[, start := start0 + 1L]
bed_dt[, mnemonic := sub("[0-9]+$", "", sub("^[^_]*_", "", label))]

fs_gr <- GRanges(seqnames = bed_dt$chr,
                 ranges   = IRanges(start = bed_dt$start, end = bed_dt$end),
                 mnemonic = bed_dt$mnemonic)
seqlevelsStyle(fs_gr) <- "UCSC"
fs_gr <- keepStandardChromosomes(fs_gr, pruning.mode = "coarse")
fs_gr <- sort(fs_gr)

gap_states <- grep("^mGapArtf", unique(mcols(fs_gr)$mnemonic), value = TRUE)
ol <- findOverlaps(bins_clean, fs_gr, ignore.strand = TRUE)
cov_dt <- data.table(bin = queryHits(ol),
                     mn  = as.character(mcols(fs_gr)$mnemonic[subjectHits(ol)]),
                     w   = width(pintersect(bins_clean[queryHits(ol)], fs_gr[subjectHits(ol)])))
bw <- width(bins_clean)
cov_sum <- cov_dt[, .(covered = sum(w)), by = .(bin, mn)]
gap_cov <- cov_sum[mn %in% gap_states, .(gap_bp = sum(covered)), by = bin]
non_gap <- data.table(bin = seq_along(bins_clean))[gap_cov, on = "bin"]
non_gap[, gap_bp := fifelse(is.na(gap_bp), 0L, gap_bp)]
non_gap[, non_gap_frac := (bw[bin] - gap_bp)/bw[bin]]
bins_use <- bins_clean[non_gap$non_gap_frac >= CONFIG$non_gap_thresh]
bins_use$bin_id <- paste0(seqnames(bins_use), ":", start(bins_use), "-", end(bins_use))

bins_df <- data.frame(chr = as.character(seqnames(bins_use)),
                      start = start(bins_use), end = end(bins_use),
                      row = seq_len(length(bins_use)))
bins_df$bin_id <- bins_use$bin_id
rownames(bins_df) <- bins_df$bin_id

export(bins_use, file.path(CONFIG$outdir, sprintf("bins100k.mm10.nonGap%02d.bed", CONFIG$non_gap_thresh*100)))
saveRDS(bins_use, file.path(CONFIG$outdir, sprintf("bins100k.mm10.nonGap%02d.rds", CONFIG$non_gap_thresh*100)))
fwrite(bins_df, file.path(CONFIG$outdir, "bins_table.tsv"), sep = "\t")

message(sprintf("  ✅ %d bins created after filtering", length(bins_use)))


## -----------------------------------------------------------
## STEP 3: EXTRACT SIGNAL FROM BIGWIGS
## -----------------------------------------------------------

box_banner("STEP 3/15: Extracting signal from bigWigs")

bw_means_in_bins <- function(bw, bins){
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

message(sprintf("  ✅ Signal matrix: %d bins × %d samples", nrow(mat_z), ncol(mat_z)))


## -----------------------------------------------------------
## STEP 4: HMM SEGMENTATION
## -----------------------------------------------------------

box_banner("STEP 4/15: Running HMM segmentation")

joint_hmm_one_chr <- function(x_mat_chr) {
    n_bins   <- nrow(x_mat_chr)
    n_groups <- ncol(x_mat_chr)
    x_mat_chr[!is.finite(x_mat_chr)] <- 0
    if (all(x_mat_chr == 0))
        return(matrix(0L, nrow = n_bins, ncol = n_groups, dimnames = dimnames(x_mat_chr)))
    fit <- HMMt::BaumWelchT(x = as.numeric(x_mat_chr),
                            series.length = rep(n_bins, n_groups),
                            maxiter = 100)
    v   <- fit$ViterbiPath
    active <- which.max(fit$mu)
    matrix(ifelse(v == active, 1L, 0L), nrow = n_bins, ncol = n_groups,
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
message(sprintf("  ✅ Compartments: %d bins × %d samples", nrow(comp_bin), ncol(comp_bin)))
message(sprintf("  ✅ Active fraction: %.1f%%", 100 * frac_active))


## -----------------------------------------------------------
## STEP 5: CHROMHMM ENRICHMENT
## -----------------------------------------------------------

box_banner("STEP 5/15: Computing ChromHMM enrichments")

ol2 <- findOverlaps(bins_use, fs_gr, ignore.strand = TRUE)
q2  <- queryHits(ol2)
s2  <- subjectHits(ol2)
ov2 <- width(pintersect(bins_use[q2], fs_gr[s2]))

cov_dt_mn <- data.table(bin = q2, mn = as.character(mcols(fs_gr)$mnemonic[s2]), w = ov2)
cov_wide_mn <- dcast(cov_dt_mn, bin ~ mn, value.var = "w", fill = 0L, fun.aggregate = sum)
bw_use <- width(bins_use)
cov_wide_mn[, total := bw_use[bin]]

mn_cols <- setdiff(names(cov_wide_mn), c("bin","total"))
for (nm in mn_cols) set(cov_wide_mn, j = nm, value = cov_wide_mn[[nm]]/cov_wide_mn$total)

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

saveRDS(log2fe_bin, file.path(CONFIG$outdir, sprintf("log2fe_bin.nonGap%02d.rds", CONFIG$non_gap_thresh * 100)))

message("  ✅ ChromHMM enrichments computed")


## -----------------------------------------------------------
## STEP 6: CPG AND AT CONTENT
## -----------------------------------------------------------

box_banner("STEP 6/15: Computing sequence features")

bin_seqs <- getSeq(BSgenome.Mmusculus.UCSC.mm10, bins_use)
cpg_count <- vcountPattern("CG", bin_seqs, fixed = TRUE)
cpg_density <- cpg_count / width(bins_use)
at_count <- letterFrequency(bin_seqs, letters = c("A","T"))
at_content <- rowSums(at_count) / width(bins_use)

gc_dt <- data.frame(
    bin_id = bins_use$bin_id,
    CpG_density = cpg_density,
    AT_content  = at_content,
    stringsAsFactors = FALSE
)
rownames(gc_dt) <- gc_dt$bin_id

saveRDS(gc_dt, file.path(CONFIG$outdir, "bins_CpG_AT.rds"))

message("  ✅ Sequence features computed")


## -----------------------------------------------------------
## STEP 7: K-MEANS CLUSTERING
## -----------------------------------------------------------

box_banner("STEP 7/15: Running k-means clustering")

km <- kmeans(comp_bin, centers = CONFIG$k_bins, nstart = 25)
groups_kept <- km$cluster
names(groups_kept) <- rownames(comp_bin)

cluster_assignments <- data.frame(
    bin_id = rownames(comp_bin),
    ecg_cluster = groups_kept,
    stringsAsFactors = FALSE
)
fwrite(cluster_assignments, file.path(CONFIG$outdir, "kmeans_cluster_assignments.tsv"), sep = "\t")

message(sprintf("  ✅ %d k-means clusters created", CONFIG$k_bins))


## -----------------------------------------------------------
## STEP 8-9: BIN ORDERING AND ECG RELABELING
## -----------------------------------------------------------

box_banner("STEP 8-9/15: Ordering bins and relabeling ECG")

frac_active_all <- rowMeans(comp_bin == 1L, na.rm = TRUE)

ord_idx <- order(groups_kept, -frac_active_all, na.last = TRUE)
ord_bins_ids <- rownames(comp_bin)[ord_idx]
ecg_split <- factor(groups_kept[ord_bins_ids])
ecg_split <- droplevels(ecg_split)

cluster_activity <- tapply(frac_active_all[ord_bins_ids], ecg_split, mean)
cluster_order <- order(cluster_activity, decreasing = TRUE)
ecg_levels_ordered <- levels(ecg_split)[cluster_order]
ecg_split <- factor(ecg_split, levels = ecg_levels_ordered)

ord_idx2 <- order(as.integer(ecg_split), -frac_active_all[ord_bins_ids])
ord_bins_ids <- ord_bins_ids[ord_idx2]
ecg_split <- ecg_split[ord_idx2]

# Relabel consecutively
present_levels <- levels(ecg_split)
n_clusters <- length(present_levels)
new_labels <- as.character(1:n_clusters)
ecg_remap <- setNames(new_labels, present_levels)
ecg_split <- factor(ecg_remap[as.character(ecg_split)], levels = new_labels)

# Build rotated matrix
comp_ord <- comp_bin[ord_bins_ids, , drop = FALSE]
mat_chr_rot <- t(apply(comp_ord, 2, function(x) ifelse(x == 1L, "1", "0")))
colnames(mat_chr_rot) <- ord_bins_ids

message(sprintf("  ✅ Bins ordered, ECG clusters relabeled (1-%d)", n_clusters))


## -----------------------------------------------------------
## STEP 10-13: BUILD ANNOTATIONS AND ROW ORDERING
## -----------------------------------------------------------

box_banner("STEP 10-13/15: Building annotations")

# Row annotations
annot_df <- data.frame(
    Sex      = sex_vec[colnames(comp_bin)],
    Age      = age_vec[colnames(comp_bin)],
    Celltype = ctype_vec[colnames(comp_bin)]
)
row_ha <- rowAnnotation(
    Age      = factor(annot_df$Age,      levels = AGE_LEVELS),
    Sex      = factor(annot_df$Sex,      levels = SEX_LEVELS),
    Celltype = factor(annot_df$Celltype, levels = names(CELLTYPE_COLS)),
    col = list(Age = AGE_COLS, Sex = SEX_COLS, Celltype = CELLTYPE_COLS)
)

# Column annotations function
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
    F <- matrix(0, nrow = nrow(A), ncol = ncol(A), dimnames = list(rownames(A), colnames(A)))
    nz <- rs > 0
    if (any(nz)) F[nz, ] <- sweep(A[nz, , drop = FALSE], 1, rs[nz], "/")
    F <- pmin(pmax(F, 0), 1)
    F <- F[ord_bins_ids, , drop = FALSE]
    t(F)
}

age_track      <- per_bin_active_strict(comp_bin, age_vec,   AGE_LEVELS,           ord_bins_ids)
sex_track      <- per_bin_active_strict(comp_bin, sex_vec,   SEX_LEVELS,           ord_bins_ids)
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
    Age = anno_barplot(age_mat_for_anno, which = "column", stack = TRUE,
                       gp = gpar(fill = unname(AGE_COLS[colnames(age_mat_for_anno)]), col = NA),
                       border = FALSE, axis = FALSE, ylim = c(0,1), bar_width = 2,
                       height = unit(22, "mm")),
    which = "column"
)
ha_sex <- HeatmapAnnotation(
    Sex = anno_barplot(sex_mat_for_anno, which = "column", stack = TRUE,
                       gp = gpar(fill = unname(SEX_COLS[colnames(sex_mat_for_anno)]), col = NA),
                       border = FALSE, axis = FALSE, ylim = c(0,1), bar_width = 2,
                       height = unit(16, "mm")),
    which = "column"
)
ha_cell <- HeatmapAnnotation(
    Celltype = anno_barplot(celltype_mat_for_anno, which = "column", stack = TRUE,
                            gp = gpar(fill = unname(CELLTYPE_COLS[colnames(celltype_mat_for_anno)]), col = NA),
                            border = FALSE, axis = FALSE, ylim = c(0,1), bar_width = 2,
                            height = unit(30, "mm")),
    which = "column"
)
top_anno_all <- c(ecg_blocks, ha_age, ha_sex, ha_cell)

# Optional tracks
chr_vec   <- bins_df$chr[match(ord_bins_ids, bins_df$bin_id)]
chr_mat   <- matrix(chr_vec, nrow = 1, dimnames = list("Chromosome", ord_bins_ids))
ht_chr <- Heatmap(chr_mat, name = "Chromosome", col = CHR_COLS,
                  cluster_rows = FALSE, cluster_columns = FALSE,
                  show_row_names = TRUE, show_column_names = FALSE,
                  height = unit(0.6,"cm"))

cpg_vec <- gc_dt[ord_bins_ids, "CpG_density", drop = TRUE]
at_vec  <- gc_dt[ord_bins_ids, "AT_content",  drop = TRUE]
cpg_brk <- quantile(cpg_vec, c(0.05, 0.50, 0.95), na.rm = TRUE)
at_brk  <- quantile(at_vec,  c(0.05, 0.50, 0.95), na.rm = TRUE)
cpg_mat <- matrix(cpg_vec, nrow = 1, dimnames = list("CpG_density", ord_bins_ids))
at_mat  <- matrix(at_vec,  nrow = 1, dimnames = list("AT_content",  ord_bins_ids))
ht_cpg <- Heatmap(cpg_mat, name = "CpG density",
                  col = circlize::colorRamp2(cpg_brk, c("white","lightpink","red")),
                  cluster_rows = FALSE, cluster_columns = FALSE,
                  show_row_names = TRUE, na_col = "white")
ht_at  <- Heatmap(at_mat, name = "AT content",
                  col = circlize::colorRamp2(at_brk, c("white","lightblue","purple")),
                  cluster_rows = FALSE, cluster_columns = FALSE,
                  show_row_names = TRUE, na_col = "white")

# ChromHMM panel
missing_cols <- setdiff(ord_bins_ids, rownames(log2fe_bin))
if (length(missing_cols)) {
    pad <- matrix(NA_real_, nrow = length(missing_cols), ncol = ncol(log2fe_bin),
                  dimnames = list(missing_cols, colnames(log2fe_bin)))
    log2fe_bin <- rbind(log2fe_bin, pad)
}
log2fe_bin <- log2fe_bin[ord_bins_ids, , drop = FALSE]
state_mat  <- t(log2fe_bin)
state_mat_z <- t(scale(t(state_mat)))
state_mat_z[!is.finite(state_mat_z)] <- NA_real_

state_order <- c(
    "mTSS","mTx","mTxEx","mTxWk","mTxEnh",
    "mEnhA","mEnhWk","mBivProm",
    "mPromF","mOpenC",
    "mReprPC","mReprPC_openC",
    "mQuies","mHET","mZNF"
)
state_order <- intersect(state_order, rownames(state_mat_z))
state_mat_z <- state_mat_z[state_order, , drop = FALSE]

cap_z <- 2.0
state_mat_z_clip <- pmax(pmin(state_mat_z, cap_z), -cap_z)
col_fun_chromhmm_z <- circlize::colorRamp2(c(-cap_z, 0, cap_z),
                                           c("#2166ac", "white", "#b2182b"))

ht_states <- Heatmap(
    state_mat_z_clip,
    name = "ChromHMM z",
    col  = col_fun_chromhmm_z,
    cluster_rows = FALSE, cluster_columns = FALSE,
    show_row_names = TRUE, show_column_names = FALSE,
    row_names_gp = gpar(fontsize = 8)
)

# Row ordering
comp_ord_numeric <- comp_bin[ord_bins_ids, , drop = FALSE]
mat_rowsamples   <- t(comp_ord_numeric)
col_w <- rep(1, ncol(mat_rowsamples))
mat_rowsamples_w <- sweep(mat_rowsamples, 2, col_w, "*")

row_age <- factor(as.character(age_vec[rownames(mat_chr_rot)]), levels = AGE_LEVELS)

k_max_per_age <- 6
order_rows_in_slice <- function(rids) {
    if (length(rids) <= 2) return(rids)
    X <- mat_rowsamples_w[rids, , drop = FALSE]
    k <- min(k_max_per_age, max(2, floor(sqrt(nrow(X)))))
    km <- kmeans(X, centers = k, nstart = 25)
    cl <- km$cluster; centers <- km$centers
    sw <- col_w / sum(col_w)
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

row_order_ids <- unlist(lapply(levels(row_age), function(ag){
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

message("  ✅ All annotations built")


## -----------------------------------------------------------
## STEP 14-15: CREATE AND SAVE HEATMAP
## -----------------------------------------------------------

box_banner("STEP 14-15/15: Creating heatmap")

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
        at = c("0","1"),
        labels = c("Repressive","Active"),
        labels_gp = gpar(fontsize = 14)
    ),
    column_split     = ecg_split,
    column_gap       = unit(0, "mm"),
    gap              = unit(0.5, "mm"),
    show_column_dend = FALSE,
    show_column_names= FALSE,
    top_annotation   = top_anno_all
)

ht_stack <- ht_main %v% ht_chr %v% ht_cpg %v% ht_at %v% ht_states

message("  Saving PDF...")
pdf_file <- file.path(CONFIG$outdir, "compartment_heatmap_full.pdf")
pdf(pdf_file, width = 24, height = 16)
draw(ht_stack, 
     heatmap_legend_side = "right",
     annotation_legend_side = "right",
     merge_legend = TRUE)
dev.off()

message("  Saving PNG...")
png_file <- file.path(CONFIG$outdir, "compartment_heatmap_full.png")
png(png_file, width = 4800, height = 3200, res = 200)
draw(ht_stack, 
     heatmap_legend_side = "right",
     annotation_legend_side = "right",
     merge_legend = TRUE)
dev.off()

message(sprintf("  ✅ Heatmap saved: %s", pdf_file))


## =========================================================
## STABILITY ANALYSIS (STEP 2)
## =========================================================

box_banner("STABILITY ANALYSIS: 5-Class Model")

# Stability functions
check_monotonicity <- function(states) {
    states <- states[!is.na(states)]
    if (length(states) < 2) {
        return(list(is_monotonic = NA, direction = "Insufficient", n_transitions = NA))
    }
    transitions <- states[-1] != states[-length(states)]
    n_trans <- sum(transitions)
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
    states[fracs >= t_on] <- "Active"
    states[fracs <= t_off] <- "Repressive"
    states[is.na(fracs)] <- NA
    valid_states <- states[states %in% c("Active", "Repressive")]
    if (length(valid_states) < 3) return("Insufficient_Data")
    n_active <- sum(valid_states == "Active")
    n_repressive <- sum(valid_states == "Repressive")
    n_total <- length(valid_states)
    if (n_active / n_total >= 0.8) return("Stable_Active")
    if (n_repressive / n_total >= 0.8) return("Stable_Repressive")
    mono_check <- check_monotonicity(valid_states)
    if (is.na(mono_check$is_monotonic)) return("Insufficient_Data")
    if (mono_check$is_monotonic) {
        if (mono_check$direction == "A_to_R") return("Monotonic_A_to_R")
        if (mono_check$direction == "R_to_A") return("Monotonic_R_to_A")
        return("Stable_Mixed")
    } else {
        return("Non_Monotonic")
    }
}

compute_stability_extended <- function(comp_bin, age_vec, sex_vec, ctype_vec,
                                       age_levels, t_on = 0.6, t_off = 0.4) {
    res <- list()
    sexes <- intersect(c("male", "female"), unique(as.character(sex_vec)))
    celltypes <- sort(unique(as.character(ctype_vec)))
    for (sx in sexes) {
        use_sx <- names(sex_vec)[as.character(sex_vec) == sx]
        for (ct in celltypes) {
            use <- intersect(use_sx, names(ctype_vec)[as.character(ctype_vec) == ct])
            if (length(use) < 2) next
            message(sprintf("  [STABILITY] %s - %s (%d samples)", sx, ct, length(use)))
            M <- comp_bin[, use, drop = FALSE]
            F <- frac_active_by(M, age_vec[use], age_levels)
            if (!is.matrix(F) || nrow(F) == 0 || ncol(F) < 2) next
            stability <- apply(F, 1, classify_stability_extended, t_on = t_on, t_off = t_off)
            res[[paste(sx, ct, sep = "|")]] <- data.table(
                bin_id = rownames(M),
                Sex = sx,
                Celltype = ct,
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
       file.path(CONFIG$outdir, "compartment_stability_5class.tsv"),
       sep = "\t")

message("\n[INFO] Stability distribution:")
print(stability_ext_dt[, .N, by = Stability][order(-N)])


## =========================================================
## STATE-PER-AGE DATA
## =========================================================

message("\n[INFO] Preparing state data per age...")

state_per_age_list <- list()
for (sample_name in colnames(comp_bin)) {
    sample_dt <- data.table(
        bin_id   = rownames(comp_bin),
        State    = ifelse(comp_bin[, sample_name] == 1L, "Active", "Repressive"),
        Age      = age_vec[sample_name],
        Sex      = sex_vec[sample_name],
        Celltype = ctype_vec[sample_name]
    )
    state_per_age_list[[sample_name]] <- sample_dt
}
state_per_age <- rbindlist(state_per_age_list)
state_per_age[, Age := factor(Age, levels = AGE_LEVELS)]

state_with_stability <- merge(
    state_per_age,
    stability_ext_dt[, .(bin_id, Sex, Celltype, Stability)],
    by = c("bin_id", "Sex", "Celltype"),
    all.x = TRUE
)

message(sprintf("  ✅ State data: %d rows", nrow(state_with_stability)))


## =========================================================
## STABILITY BARPLOTS
## =========================================================

message("\n[INFO] Creating stability barplots...")

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
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(CONFIG$outdir, "Barplot_Stability_5class.pdf"),
       p_5class, width = 12, height = 8)

message("  ✅ Stability barplots saved")


## =========================================================
## ALLUVIAL PLOTS - ALL VARIANTS
## =========================================================

box_banner("ALLUVIAL PLOTS: All Variants")

plot_data <- state_with_stability[!Stability %in% c("Missing", "Insufficient_Data")]
switching_data <- plot_data[Stability %in% SWITCHING_CLASSES]

message(sprintf("  Plotting data: %d rows", nrow(plot_data)))
message(sprintf("  Switching data: %d rows", nrow(switching_data)))


## -----------------------------------------------------------
## VARIANT 1: Per stability class - faceted by celltype
## -----------------------------------------------------------

banner("ALLUVIAL 1: Per stability class (faceted by celltype)")

for (stab in STABILITY_LEVELS) {
    stab_clean <- clean_name(stab)
    
    for (sx in c("male", "female")) {
        dat <- plot_data[Sex == sx & Stability == stab]
        if (nrow(dat) == 0) next
        
        dat_wide <- dcast(dat, bin_id + Celltype ~ Age, value.var = "State", 
                          fun.aggregate = function(x) x[1])
        dat_wide <- dat_wide[complete.cases(dat_wide)]
        if (nrow(dat_wide) == 0) next
        
        agg <- dat_wide[, .N, by = c("Celltype", AGE_LEVELS)]
        agg[, trajectory_id := do.call(paste, c(.SD, sep = "_")), .SDcols = c("Celltype", AGE_LEVELS)]
        
        agg_long <- melt(agg,
                         id.vars = c("N", "Celltype", "trajectory_id"),
                         measure.vars = AGE_LEVELS,
                         variable.name = "Age",
                         value.name = "State")
        agg_long[, Age := factor(Age, levels = AGE_LEVELS)]
        agg_long <- agg_long[!is.na(State)]
        
        ct_order <- agg_long[, .(total = sum(N)), by = Celltype][order(-total), Celltype]
        agg_long[, Celltype := factor(Celltype, levels = ct_order)]
        
        total_bins <- sum(agg$N)
        
        p <- ggplot(agg_long,
                    aes(x = Age, stratum = State, alluvium = trajectory_id,
                        y = N, fill = State)) +
            stat_alluvium(geom = "flow", alpha = 0.7, width = 0.4,
                          curve_type = "linear", na.rm = TRUE) +
            stat_stratum(width = 0.4, alpha = 0.9, color = "black",
                         linewidth = 0.3, na.rm = TRUE) +
            geom_text(stat = "stratum", aes(label = State),
                      size = 2.5, fontface = "bold", na.rm = TRUE) +
            facet_wrap(~ Celltype, ncol = 3, scales = "free_y") +
            scale_fill_manual(values = STATE_COLS, na.translate = FALSE) +
            labs(
                title = sprintf("%s - %s", gsub("_", " ", stab), tools::toTitleCase(sx)),
                subtitle = sprintf("Total: %s bins", format(total_bins, big.mark = ",")),
                x = "Age", y = "Number of 80kb Bins", fill = "State"
            ) +
            theme_bw(base_size = 11) +
            theme(
                axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 9),
                strip.text = element_text(size = 10, face = "bold"),
                strip.background = element_rect(fill = "gray95"),
                plot.title = element_text(face = "bold", size = 14),
                legend.position = "bottom"
            )
        
        fname <- sprintf("Alluvial_%s_%s_ByCelltype.pdf", stab_clean, tools::toTitleCase(sx))
        ggsave(file.path(CONFIG$outdir, fname), p, width = 14, height = 10)
        message(sprintf("  [SAVED] %s", fname))
    }
}


## -----------------------------------------------------------
## VARIANT 2: Celltype-colored (merged, no faceting)
## -----------------------------------------------------------

banner("ALLUVIAL 2: Celltype-colored merged plots")

for (stab in STABILITY_LEVELS) {
    stab_clean <- clean_name(stab)
    
    for (sx in c("male", "female")) {
        dat <- plot_data[Sex == sx & Stability == stab]
        if (nrow(dat) == 0) next
        
        dat_wide <- dcast(dat, bin_id + Celltype ~ Age, value.var = "State",
                          fun.aggregate = function(x) x[1])
        dat_wide <- dat_wide[complete.cases(dat_wide)]
        if (nrow(dat_wide) == 0) next
        
        agg <- dat_wide[, .N, by = c("Celltype", AGE_LEVELS)]
        agg[, trajectory_id := do.call(paste, c(.SD, sep = "_")), .SDcols = c("Celltype", AGE_LEVELS)]
        
        agg_long <- melt(agg,
                         id.vars = c("N", "Celltype", "trajectory_id"),
                         measure.vars = AGE_LEVELS,
                         variable.name = "Age",
                         value.name = "State")
        agg_long[, Age := factor(Age, levels = AGE_LEVELS)]
        agg_long <- agg_long[!is.na(State)]
        
        total_bins <- sum(agg$N)
        
        p <- ggplot(agg_long,
                    aes(x = Age, stratum = State, alluvium = trajectory_id,
                        y = N, fill = Celltype)) +
            stat_alluvium(geom = "flow", alpha = 0.6, width = 0.4,
                          curve_type = "linear", na.rm = TRUE) +
            stat_stratum(aes(fill = Celltype), width = 0.4, alpha = 0.9, 
                         color = "black", linewidth = 0.3, na.rm = TRUE) +
            geom_text(stat = "stratum", aes(label = State),
                      size = 3, fontface = "bold", na.rm = TRUE) +
            scale_fill_manual(values = CELLTYPE_COLS, na.translate = FALSE) +
            labs(
                title = sprintf("%s - %s (Colored by Celltype)",
                                gsub("_", " ", stab), tools::toTitleCase(sx)),
                subtitle = sprintf("Total: %s bins", format(total_bins, big.mark = ",")),
                x = "Age", y = "Number of 80kb Bins", fill = "Cell Type"
            ) +
            theme_bw(base_size = 12) +
            theme(
                axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
                plot.title = element_text(face = "bold", size = 14),
                legend.position = "right"
            )
        
        fname <- sprintf("Alluvial_%s_%s_CelltypeColored.pdf", stab_clean, tools::toTitleCase(sx))
        ggsave(file.path(CONFIG$outdir, fname), p, width = 12, height = 8)
        message(sprintf("  [SAVED] %s", fname))
    }
}


## -----------------------------------------------------------
## VARIANT 3: Merged alluvials (Male vs Female side-by-side)
## -----------------------------------------------------------

banner("ALLUVIAL 3: Merged male/female comparisons")

for (stab in STABILITY_LEVELS) {
    stab_clean <- clean_name(stab)
    
    dat <- plot_data[Stability == stab]
    if (nrow(dat) == 0) next
    
    dat_wide <- dcast(dat, bin_id + Sex ~ Age, value.var = "State", fun.aggregate = function(x) x[1])
    dat_wide <- dat_wide[complete.cases(dat_wide)]
    if (nrow(dat_wide) == 0) next
    
    agg <- dat_wide[, .N, by = c("Sex", AGE_LEVELS)]
    agg[, trajectory_id := do.call(paste, c(.SD, sep = "_")), .SDcols = c("Sex", AGE_LEVELS)]
    
    agg_long <- melt(agg,
                     id.vars = c("N", "Sex", "trajectory_id"),
                     measure.vars = AGE_LEVELS,
                     variable.name = "Age",
                     value.name = "State")
    agg_long[, Age := factor(Age, levels = AGE_LEVELS)]
    agg_long[, Sex := factor(Sex, levels = c("male", "female"), labels = c("Male", "Female"))]
    agg_long <- agg_long[!is.na(State)]
    
    p <- ggplot(agg_long,
                aes(x = Age, stratum = State, alluvium = trajectory_id,
                    y = N, fill = State)) +
        stat_alluvium(geom = "flow", alpha = 0.7, width = 0.4,
                      curve_type = "linear", na.rm = TRUE) +
        stat_stratum(width = 0.4, alpha = 0.9, color = "black",
                     linewidth = 0.3, na.rm = TRUE) +
        geom_text(stat = "stratum", aes(label = State),
                  size = 3, fontface = "bold", na.rm = TRUE) +
        facet_wrap(~ Sex, ncol = 2) +
        scale_fill_manual(values = STATE_COLS, na.translate = FALSE) +
        labs(
            title = gsub("_", " ", stab),
            subtitle = "All cell types merged | Male vs Female comparison",
            x = "Age", y = "Number of 80kb Bins", fill = "State"
        ) +
        theme_bw(base_size = 12) +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
            strip.text = element_text(size = 14, face = "bold"),
            strip.background = element_rect(fill = "gray90"),
            plot.title = element_text(face = "bold", size = 16)
        )
    
    fname <- sprintf("Alluvial_Merged_%s.pdf", stab_clean)
    ggsave(file.path(CONFIG$outdir, fname), p, width = 16, height = 5)
    message(sprintf("  [SAVED] %s", fname))
}


## -----------------------------------------------------------
## VARIANT 4: All stability classes combined (faceted by celltype)
## -----------------------------------------------------------

banner("ALLUVIAL 4: All classes combined by celltype")

for (sx in c("male", "female")) {
    
    dat <- plot_data[Sex == sx]
    
    dat_wide <- dcast(dat, bin_id + Celltype + Stability ~ Age, value.var = "State",
                      fun.aggregate = function(x) x[1])
    dat_wide <- dat_wide[complete.cases(dat_wide)]
    
    if (nrow(dat_wide) == 0) next
    
    agg <- dat_wide[, .N, by = c("Celltype", "Stability", AGE_LEVELS)]
    agg[, trajectory_id := do.call(paste, c(.SD, sep = "_")),
        .SDcols = c("Celltype", "Stability", AGE_LEVELS)]
    
    agg_long <- melt(agg,
                     id.vars = c("N", "Celltype", "Stability", "trajectory_id"),
                     measure.vars = AGE_LEVELS,
                     variable.name = "Age",
                     value.name = "State")
    agg_long[, Age := factor(Age, levels = AGE_LEVELS)]
    agg_long <- agg_long[!is.na(State)]
    
    ct_order <- names(CELLTYPE_COLS)[names(CELLTYPE_COLS) %in% unique(agg_long$Celltype)]
    agg_long[, Celltype := factor(Celltype, levels = ct_order)]
    
    p <- ggplot(agg_long,
                aes(x = Age, stratum = State, alluvium = trajectory_id,
                    y = N, fill = Stability)) +
        stat_alluvium(geom = "flow", alpha = 0.7, width = 0.4,
                      curve_type = "linear", na.rm = TRUE) +
        stat_stratum(width = 0.4, alpha = 0.9, color = "black",
                     linewidth = 0.3, na.rm = TRUE) +
        geom_text(stat = "stratum", aes(label = State),
                  size = 2, fontface = "bold", na.rm = TRUE) +
        facet_wrap(~ Celltype, ncol = 5, scales = "free_y") +
        scale_fill_manual(values = STABILITY_COLS, na.translate = FALSE) +
        labs(
            title = sprintf("All Stability Classes - %s", tools::toTitleCase(sx)),
            subtitle = sprintf("Total: %s bins", format(sum(agg$N), big.mark = ",")),
            x = "Age", y = "Number of 80kb Bins", fill = "Stability Class"
        ) +
        theme_bw(base_size = 10) +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 8),
            strip.text = element_text(size = 9, face = "bold"),
            strip.background = element_rect(fill = "gray95"),
            plot.title = element_text(face = "bold", size = 14),
            legend.position = "bottom"
        )
    
    fname <- sprintf("Alluvial_AllClasses_%s_ByCelltype.pdf", tools::toTitleCase(sx))
    ggsave(file.path(CONFIG$outdir, fname), p, width = 18, height = 14)
    message(sprintf("  [SAVED] %s", fname))
}


## -----------------------------------------------------------
## VARIANT 5: Switching bins only (by celltype)
## -----------------------------------------------------------

banner("ALLUVIAL 5: Switching bins only (by celltype)")

for (sx in c("male", "female")) {
    
    dat <- switching_data[Sex == sx]
    dat_wide <- dcast(dat, bin_id + Celltype + Stability ~ Age,
                      value.var = "State", fun.aggregate = function(x) x[1])
    dat_wide <- dat_wide[complete.cases(dat_wide)]
    
    n_bins <- nrow(dat_wide)
    message(sprintf("  %s: %d switching bins", sx, n_bins))
    
    if (n_bins == 0) next
    
    agg <- dat_wide[, .N, by = c("Celltype", "Stability", AGE_LEVELS)]
    agg[, trajectory_id := .I]
    
    agg_long <- melt(agg,
                     id.vars = c("N", "Celltype", "Stability", "trajectory_id"),
                     measure.vars = AGE_LEVELS,
                     variable.name = "Age",
                     value.name = "State")
    agg_long[, Age := factor(Age, levels = AGE_LEVELS)]
    agg_long <- agg_long[!is.na(State)]
    
    # Colored by CELLTYPE
    p_ct <- ggplot(agg_long,
                   aes(x = Age, stratum = State, alluvium = trajectory_id,
                       y = N, fill = Celltype)) +
        stat_alluvium(geom = "flow", alpha = 0.7, width = 0.4,
                      curve_type = "linear", na.rm = TRUE) +
        stat_stratum(width = 0.4, alpha = 0.9, color = "black",
                     linewidth = 0.3, na.rm = TRUE) +
        geom_text(stat = "stratum", aes(label = State),
                  size = 5, fontface = "bold", na.rm = TRUE) +
        scale_fill_manual(values = CELLTYPE_COLS, na.translate = FALSE) +
        scale_y_continuous(labels = scales::comma, expand = expansion(mult = c(0, 0.02))) +
        labs(
            title = sprintf("Switching Compartments - %s", tools::toTitleCase(sx)),
            subtitle = sprintf("Monotonic + Non-Monotonic | Total: %s bins", format(n_bins, big.mark = ",")),
            x = "Age", y = "Number of 80kb Bins", fill = "Cell Type"
        ) +
        theme_bw(base_size = 16) +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 14),
            plot.title = element_text(face = "bold", size = 20),
            legend.position = "right"
        )
    
    fname <- sprintf("Alluvial_Switching_ByCelltype_%s.pdf", tools::toTitleCase(sx))
    ggsave(file.path(CONFIG$outdir, fname), p_ct, width = 14, height = 10)
    message(sprintf("  [SAVED] %s", fname))
    
    # Colored by SWITCH TYPE
    p_sw <- ggplot(agg_long,
                   aes(x = Age, stratum = State, alluvium = trajectory_id,
                       y = N, fill = Stability)) +
        stat_alluvium(geom = "flow", alpha = 0.7, width = 0.4,
                      curve_type = "linear", na.rm = TRUE) +
        stat_stratum(width = 0.4, alpha = 0.9, color = "black",
                     linewidth = 0.3, na.rm = TRUE) +
        geom_text(stat = "stratum", aes(label = State),
                  size = 5, fontface = "bold", na.rm = TRUE) +
        scale_fill_manual(values = SWITCHING_COLS,
                          labels = c("Monotonic_A_to_R" = "A → R (Closing)",
                                     "Monotonic_R_to_A" = "R → A (Opening)",
                                     "Non_Monotonic" = "Non-Monotonic"),
                          na.translate = FALSE) +
        scale_y_continuous(labels = scales::comma, expand = expansion(mult = c(0, 0.02))) +
        labs(
            title = sprintf("Switching Compartments - %s", tools::toTitleCase(sx)),
            subtitle = sprintf("Monotonic + Non-Monotonic | Total: %s bins", format(n_bins, big.mark = ",")),
            x = "Age", y = "Number of 80kb Bins", fill = "Switch Type"
        ) +
        theme_bw(base_size = 16) +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 14),
            plot.title = element_text(face = "bold", size = 20),
            legend.position = "right"
        )
    
    fname <- sprintf("Alluvial_Switching_BySwitchType_%s.pdf", tools::toTitleCase(sx))
    ggsave(file.path(CONFIG$outdir, fname), p_sw, width = 14, height = 10)
    message(sprintf("  [SAVED] %s", fname))
}


## -----------------------------------------------------------
## VARIANT 6: Switching side-by-side (Male vs Female)
## -----------------------------------------------------------

banner("ALLUVIAL 6: Switching side-by-side comparison")

switching_both <- switching_data[Sex %in% c("male", "female")]

dat_wide_sw <- dcast(switching_both, bin_id + Sex + Celltype + Stability ~ Age,
                     value.var = "State", fun.aggregate = function(x) x[1])
dat_wide_sw <- dat_wide_sw[complete.cases(dat_wide_sw)]

n_male <- sum(dat_wide_sw$Sex == "male")
n_female <- sum(dat_wide_sw$Sex == "female")
message(sprintf("  Male: %d bins | Female: %d bins", n_male, n_female))

if (nrow(dat_wide_sw) > 0) {
    
    agg_sw <- dat_wide_sw[, .N, by = c("Sex", "Celltype", "Stability", AGE_LEVELS)]
    agg_sw[, trajectory_id := .I]
    
    agg_long_sw <- melt(agg_sw,
                        id.vars = c("N", "Sex", "Celltype", "Stability", "trajectory_id"),
                        measure.vars = AGE_LEVELS,
                        variable.name = "Age",
                        value.name = "State")
    agg_long_sw[, Age := factor(Age, levels = AGE_LEVELS)]
    agg_long_sw[, Sex := factor(Sex, levels = c("male", "female"), labels = c("Male", "Female"))]
    agg_long_sw <- agg_long_sw[!is.na(State)]
    
    # By Celltype
    p_ct <- ggplot(agg_long_sw,
                   aes(x = Age, stratum = State, alluvium = trajectory_id,
                       y = N, fill = Celltype)) +
        stat_alluvium(geom = "flow", alpha = 0.7, width = 0.4,
                      curve_type = "linear", na.rm = TRUE) +
        stat_stratum(width = 0.4, alpha = 0.9, color = "black",
                     linewidth = 0.3, na.rm = TRUE) +
        geom_text(stat = "stratum", aes(label = State),
                  size = 4, fontface = "bold", na.rm = TRUE) +
        facet_wrap(~ Sex, ncol = 2) +
        scale_fill_manual(values = CELLTYPE_COLS, na.translate = FALSE) +
        scale_y_continuous(labels = scales::comma, expand = expansion(mult = c(0, 0.02))) +
        labs(
            title = "Switching Compartments - Male vs Female",
            subtitle = sprintf("Monotonic + Non-Monotonic | Male: %s | Female: %s bins",
                               format(n_male, big.mark = ","), format(n_female, big.mark = ",")),
            x = "Age", y = "Number of 80kb Bins", fill = "Cell Type"
        ) +
        theme_bw(base_size = 14) +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 12),
            strip.text = element_text(size = 18, face = "bold"),
            strip.background = element_rect(fill = "gray95"),
            plot.title = element_text(face = "bold", size = 20),
            legend.position = "right"
        )
    
    ggsave(file.path(CONFIG$outdir, "Alluvial_Switching_SideBySide_ByCelltype.pdf"),
           p_ct, width = 23, height = 8)
    message("  [SAVED] Alluvial_Switching_SideBySide_ByCelltype.pdf")
    
    # By Switch Type
    p_sw <- ggplot(agg_long_sw,
                   aes(x = Age, stratum = State, alluvium = trajectory_id,
                       y = N, fill = Stability)) +
        stat_alluvium(geom = "flow", alpha = 0.7, width = 0.4,
                      curve_type = "linear", na.rm = TRUE) +
        stat_stratum(width = 0.4, alpha = 0.9, color = "black",
                     linewidth = 0.3, na.rm = TRUE) +
        geom_text(stat = "stratum", aes(label = State),
                  size = 4, fontface = "bold", na.rm = TRUE) +
        facet_wrap(~ Sex, ncol = 2) +
        scale_fill_manual(values = SWITCHING_COLS,
                          labels = c("Monotonic_A_to_R" = "A → R (Closing)",
                                     "Monotonic_R_to_A" = "R → A (Opening)",
                                     "Non_Monotonic" = "Non-Monotonic"),
                          na.translate = FALSE) +
        scale_y_continuous(labels = scales::comma, expand = expansion(mult = c(0, 0.02))) +
        labs(
            title = "Switching Compartments - Male vs Female",
            subtitle = sprintf("Monotonic + Non-Monotonic | Male: %s | Female: %s bins",
                               format(n_male, big.mark = ","), format(n_female, big.mark = ",")),
            x = "Age", y = "Number of 80kb Bins", fill = "Switch Type"
        ) +
        theme_bw(base_size = 14) +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 12),
            strip.text = element_text(size = 18, face = "bold"),
            strip.background = element_rect(fill = "gray95"),
            plot.title = element_text(face = "bold", size = 20),
            legend.position = "right"
        )
    
    ggsave(file.path(CONFIG$outdir, "Alluvial_Switching_SideBySide_BySwitchType.pdf"),
           p_sw, width = 23, height = 8)
    message("  [SAVED] Alluvial_Switching_SideBySide_BySwitchType.pdf")
}


## -----------------------------------------------------------
## VARIANT 7: Combined all (no faceting, colored by celltype/stability)
## -----------------------------------------------------------

banner("ALLUVIAL 7: Combined all (no faceting)")

for (sx in c("male", "female")) {
    
    dat <- plot_data[Sex == sx]
    
    dat_wide <- dcast(dat, bin_id + Celltype + Stability ~ Age, 
                      value.var = "State", fun.aggregate = function(x) x[1])
    dat_wide <- dat_wide[complete.cases(dat_wide)]
    
    n_bins <- nrow(dat_wide)
    message(sprintf("  %s: %d complete bins", sx, n_bins))
    
    if (n_bins == 0) next
    
    agg <- dat_wide[, .N, by = c("Celltype", "Stability", AGE_LEVELS)]
    agg[, trajectory_id := .I]
    
    agg_long <- melt(agg,
                     id.vars = c("N", "Celltype", "Stability", "trajectory_id"),
                     measure.vars = AGE_LEVELS,
                     variable.name = "Age",
                     value.name = "State")
    agg_long[, Age := factor(Age, levels = AGE_LEVELS)]
    agg_long <- agg_long[!is.na(State)]
    
    # By Celltype
    p_ct <- ggplot(agg_long,
                   aes(x = Age, stratum = State, alluvium = trajectory_id,
                       y = N, fill = Celltype)) +
        stat_alluvium(geom = "flow", alpha = 0.6, width = 0.4,
                      curve_type = "linear", na.rm = TRUE) +
        stat_stratum(width = 0.4, alpha = 0.9, color = "black",
                     linewidth = 0.2, na.rm = TRUE) +
        geom_text(stat = "stratum", aes(label = State),
                  size = 4, fontface = "bold", na.rm = TRUE) +
        scale_fill_manual(values = CELLTYPE_COLS, na.translate = FALSE) +
        labs(
            title = sprintf("Compartment Dynamics - %s", tools::toTitleCase(sx)),
            subtitle = sprintf("All stability classes + all celltypes | %s bins", 
                               format(n_bins, big.mark = ",")),
            x = "Age", y = "Number of 80kb Bins", fill = "Cell Type"
        ) +
        theme_bw(base_size = 14) +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 12),
            plot.title = element_text(face = "bold", size = 18),
            legend.position = "right"
        )
    
    fname <- sprintf("Alluvial_Combined_AllStability_%s_ByCelltype.pdf", tools::toTitleCase(sx))
    ggsave(file.path(CONFIG$outdir, fname), p_ct, width = 14, height = 10)
    message(sprintf("  [SAVED] %s", fname))
    
    # By Stability
    p_stab <- ggplot(agg_long,
                     aes(x = Age, stratum = State, alluvium = trajectory_id,
                         y = N, fill = Stability)) +
        stat_alluvium(geom = "flow", alpha = 0.6, width = 0.4,
                      curve_type = "linear", na.rm = TRUE) +
        stat_stratum(width = 0.4, alpha = 0.9, color = "black",
                     linewidth = 0.2, na.rm = TRUE) +
        geom_text(stat = "stratum", aes(label = State),
                  size = 4, fontface = "bold", na.rm = TRUE) +
        scale_fill_manual(values = STABILITY_COLS, na.translate = FALSE) +
        labs(
            title = sprintf("Compartment Dynamics - %s", tools::toTitleCase(sx)),
            subtitle = sprintf("All stability classes + all celltypes | %s bins",
                               format(n_bins, big.mark = ",")),
            x = "Age", y = "Number of 80kb Bins", fill = "Stability Class"
        ) +
        theme_bw(base_size = 14) +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 12),
            plot.title = element_text(face = "bold", size = 18),
            legend.position = "right"
        )
    
    fname <- sprintf("Alluvial_Combined_AllStability_%s_ByStability.pdf", tools::toTitleCase(sx))
    ggsave(file.path(CONFIG$outdir, fname), p_stab, width = 14, height = 10)
    message(sprintf("  [SAVED] %s", fname))
}

message("\n  ✅ ALL ALLUVIAL PLOTS SAVED")


## =========================================================
## FINAL SUMMARY
## =========================================================

end_time <- Sys.time()
duration <- difftime(end_time, start_time, units = "mins")

cat("\n")
cat("╔══════════════════════════════════════════════════════════════════╗\n")
cat("║                    PIPELINE COMPLETE                             ║\n")
cat("╠══════════════════════════════════════════════════════════════════╣\n")
cat(sprintf("║  Duration: %.1f minutes                                         ║\n", as.numeric(duration)))
cat(sprintf("║  Output: %-55s║\n", substr(CONFIG$outdir, 1, 55)))
cat("║                                                                  ║\n")
cat("║  Key outputs:                                                    ║\n")
cat("║    • compartment_heatmap_full.pdf/png                            ║\n")
cat("║    • compartments_binary_A1_R0.rds                               ║\n")
cat("║    • compartment_stability_5class.tsv                            ║\n")
cat("║    • Barplot_Stability_5class.pdf                                ║\n")
cat("║                                                                  ║\n")
cat("║  Alluvial plots (all variants):                                  ║\n")
cat("║    1. Alluvial_*_Male/Female_ByCelltype.pdf (faceted)            ║\n")
cat("║    2. Alluvial_*_Male/Female_CelltypeColored.pdf (merged)        ║\n")
cat("║    3. Alluvial_Merged_*.pdf (male vs female)                     ║\n")
cat("║    4. Alluvial_AllClasses_*_ByCelltype.pdf (all combined)        ║\n")
cat("║    5. Alluvial_Switching_ByCelltype/BySwitchType_*.pdf           ║\n")
cat("║    6. Alluvial_Switching_SideBySide_*.pdf                        ║\n")
cat("║    7. Alluvial_Combined_AllStability_*_By*.pdf                   ║\n")
cat("╚══════════════════════════════════════════════════════════════════╝\n")

message("\n[DONE] All analyses complete!")
