#!/usr/bin/env Rscript
## =========================================================
## COMPARTMENT ANALYSIS - MAIN PIPELINE
## 
## Steps:
##   1.  Load manifest
##   2.  Create genomic bins + blacklist + gap filter
##   3.  Summarize bigWigs → Z-score matrix
##   4.  HMM segmentation (A/B compartments)
##   5.  ChromHMM overlap → per-bin log2FE
##   6.  CpG density & AT content
##   7.  K-means clustering (ECG groups)
##   8.  Bin ordering by ECG and activity
##   9.  Relabel ECG clusters consecutively
##   10. Build annotation matrices
##   11. Build optional tracks
##   12. ChromHMM panel
##   13. Row ordering
##   14. Final validation
##   15. Create and save heatmap
##
## Prerequisites:
##   - Source compartment_00_config.R first
##   - BigWig files from ArchR (Step 9)
##
## Output:
##   - comp_bin: Compartment binary matrix (bins × samples)
##   - mat_z: Z-score signal matrix
##   - bins_df: Bin coordinates table
##   - Heatmap PDF/PNG
##
## =========================================================

box_banner("COMPARTMENT ANALYSIS PIPELINE - MAIN")

## =========================================================
## LOAD DEPENDENCIES
## =========================================================

# Source configuration (if not already loaded)
if (!exists("CONFIG")) {
    source("compartment_00_config.R")
}

# Load required libraries
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
})

addArchRThreads(threads = CONFIG$threads)


## =========================================================
## DIAGNOSTIC FUNCTIONS
## =========================================================

sanity_check_manifest <- function(man) {
    message("\n[SANITY CHECK 1] Manifest Loading")
    stopifnot("bigwig column missing" = "bigwig" %in% names(man))
    stopifnot("age column missing" = "age" %in% names(man))
    stopifnot("sex column missing" = "sex" %in% names(man))
    stopifnot("celltype column missing" = "celltype" %in% names(man))
    ages_found <- unique(as.character(man$age))
    missing_ages <- setdiff(ages_found, AGE_LEVELS)
    if (length(missing_ages) > 0) {
        warning(sprintf("Unexpected age values: %s", paste(missing_ages, collapse = ", ")))
    }
    sexes_found <- unique(as.character(man$sex))
    message(sprintf("  ✅ %d samples loaded", nrow(man)))
    message(sprintf("  ✅ Ages: %s", paste(ages_found, collapse = ", ")))
    message(sprintf("  ✅ Sexes: %s", paste(unique(sexes_found), collapse = ", ")))
    message(sprintf("  ✅ Celltypes: %d unique", length(unique(man$celltype))))
}

sanity_check_bins <- function(bins_use, bins_df) {
    message("\n[SANITY CHECK 2] Bin Creation")
    stopifnot("bins_use is empty" = length(bins_use) > 0)
    stopifnot("bins_df is empty" = nrow(bins_df) > 0)
    stopifnot("bin_id mismatch" = identical(bins_use$bin_id, bins_df$bin_id))
    message(sprintf("  ✅ %d bins created", length(bins_use)))
    message(sprintf("  ✅ Chromosomes: %d", length(unique(bins_df$chr))))
}

sanity_check_signal_matrix <- function(mat_z, bins_df, man) {
    message("\n[SANITY CHECK 3] Signal Matrix")
    stopifnot("mat_z rows != bins" = nrow(mat_z) == nrow(bins_df))
    stopifnot("mat_z cols != samples" = ncol(mat_z) == nrow(man))
    stopifnot("mat_z rownames missing" = !is.null(rownames(mat_z)))
    stopifnot("mat_z colnames missing" = !is.null(colnames(mat_z)))
    na_frac <- mean(is.na(mat_z))
    message(sprintf("  ✅ Matrix: %d bins × %d samples", nrow(mat_z), ncol(mat_z)))
    message(sprintf("  ✅ NA fraction: %.2f%%", 100 * na_frac))
    message(sprintf("  ✅ Z-score range: [%.2f, %.2f]", 
                    min(mat_z, na.rm = TRUE), max(mat_z, na.rm = TRUE)))
}

sanity_check_compartments <- function(comp_bin, mat_z) {
    message("\n[SANITY CHECK 4] Compartment Calling (HMM)")
    stopifnot("dimensions mismatch" = identical(dim(comp_bin), dim(mat_z)))
    stopifnot("rownames mismatch" = identical(rownames(comp_bin), rownames(mat_z)))
    stopifnot("colnames mismatch" = identical(colnames(comp_bin), colnames(mat_z)))
    unique_vals <- unique(as.vector(comp_bin))
    stopifnot("Invalid values" = all(unique_vals %in% c(0L, 1L, NA)))
    frac_active <- mean(comp_bin == 1L, na.rm = TRUE)
    message(sprintf("  ✅ Dimensions: %d × %d", nrow(comp_bin), ncol(comp_bin)))
    message(sprintf("  ✅ Active fraction: %.1f%%", 100 * frac_active))
    if (frac_active < 0.2 || frac_active > 0.8) {
        warning(sprintf("Unusual active fraction: %.1f%%", 100 * frac_active))
    }
}

sanity_check_kmeans <- function(groups_kept, comp_bin, k_bins) {
    message("\n[SANITY CHECK 5] K-means Clustering (ECG)")
    stopifnot("assignment count" = length(groups_kept) == nrow(comp_bin))
    stopifnot("names match" = identical(names(groups_kept), rownames(comp_bin)))
    cluster_sizes <- table(groups_kept)
    message(sprintf("  ✅ %d clusters created", length(cluster_sizes)))
    message(sprintf("  ✅ Size range: [%d, %d]", min(cluster_sizes), max(cluster_sizes)))
    message(sprintf("  ✅ Median size: %d", median(cluster_sizes)))
}

sanity_check_bin_ordering <- function(ord_bins_ids, ecg_split, comp_bin) {
    message("\n[SANITY CHECK 6] Bin Ordering")
    stopifnot("ord_bins_ids length" = length(ord_bins_ids) == nrow(comp_bin))
    stopifnot("ecg_split length" = length(ecg_split) == length(ord_bins_ids))
    stopifnot("all bins present" = all(ord_bins_ids %in% rownames(comp_bin)))
    stopifnot("no duplicates" = !any(duplicated(ord_bins_ids)))
    stopifnot("complete set" = setequal(ord_bins_ids, rownames(comp_bin)))
    message(sprintf("  ✅ %d bins ordered", length(ord_bins_ids)))
    message(sprintf("  ✅ %d ECG clusters", nlevels(ecg_split)))
}

sanity_check_rotation <- function(mat_chr_rot, comp_bin, ord_bins_ids) {
    message("\n[SANITY CHECK 7] Matrix Rotation")
    stopifnot("rows = samples" = nrow(mat_chr_rot) == ncol(comp_bin))
    stopifnot("cols = bins" = ncol(mat_chr_rot) == nrow(comp_bin))
    stopifnot("rownames are samples" = setequal(rownames(mat_chr_rot), colnames(comp_bin)))
    stopifnot("colnames = ord_bins_ids" = identical(colnames(mat_chr_rot), ord_bins_ids))
    message(sprintf("  ✅ Rotated: %d samples × %d bins", nrow(mat_chr_rot), ncol(mat_chr_rot)))
    message(sprintf("  ✅ Column alignment verified"))
}

sanity_check_annotation_matrices <- function(age_mat, sex_mat, cell_mat, mat_chr_rot) {
    message("\n[SANITY CHECK 8] Annotation Matrices")
    n_bins <- ncol(mat_chr_rot)
    stopifnot("age_mat rows" = nrow(age_mat) == n_bins)
    stopifnot("sex_mat rows" = nrow(sex_mat) == n_bins)
    stopifnot("cell_mat rows" = nrow(cell_mat) == n_bins)
    message(sprintf("  ✅ Age: %d × %d", nrow(age_mat), ncol(age_mat)))
    message(sprintf("  ✅ Sex: %d × %d", nrow(sex_mat), ncol(sex_mat)))
    message(sprintf("  ✅ Celltype: %d × %d", nrow(cell_mat), ncol(cell_mat)))
}

sanity_check_row_ordering <- function(row_order_idx, row_split_age, mat_chr_rot) {
    message("\n[SANITY CHECK 9] Row Ordering")
    n_samples <- nrow(mat_chr_rot)
    stopifnot("row_order length" = length(row_order_idx) == n_samples)
    stopifnot("row_split length" = length(row_split_age) == n_samples)
    stopifnot("indices in range" = all(row_order_idx >= 1 & row_order_idx <= n_samples))
    stopifnot("indices unique" = !any(duplicated(row_order_idx)))
    stopifnot("all rows present" = setequal(row_order_idx, seq_len(n_samples)))
    message(sprintf("  ✅ %d samples ordered", n_samples))
    message(sprintf("  ✅ %d age groups", nlevels(row_split_age)))
}

sanity_check_heatmap_final <- function(mat_chr_rot, row_split_age, row_order_idx, ecg_split) {
    message("\n[SANITY CHECK 10] Final Heatmap Inputs")
    stopifnot("row_split vs matrix" = length(row_split_age) == nrow(mat_chr_rot))
    stopifnot("row_order vs matrix" = length(row_order_idx) == nrow(mat_chr_rot))
    stopifnot("ecg_split vs matrix" = length(ecg_split) == ncol(mat_chr_rot))
    message(sprintf("  ✅ Matrix: %d × %d", nrow(mat_chr_rot), ncol(mat_chr_rot)))
    message(sprintf("  ✅ Row split: %d elements", length(row_split_age)))
    message(sprintf("  ✅ Row order: %d indices", length(row_order_idx)))
    message(sprintf("  ✅ Col split: %d elements", length(ecg_split)))
    message("\n  ════════════════════════════════════════")
    message("  ✅ ALL SANITY CHECKS PASSED")
    message("  ════════════════════════════════════════\n")
}

run_full_diagnostics <- function(comp_bin, mat_chr_rot, ord_bins_ids, ecg_split,
                                  row_split_age, row_order_idx, age_vec, sex_vec,
                                  ctype_vec, bins_df, outdir) {
    log_file <- file.path(outdir, "compartment_diagnostics_full.log")
    sink(log_file, split = TRUE)
    
    cat("\n╔══════════════════════════════════════════════════════════════════╗\n")
    cat("║        FULL DIAGNOSTIC REPORT                                    ║\n")
    cat("║        Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "                    ║\n")
    cat("╚══════════════════════════════════════════════════════════════════╝\n")
    
    all_pass <- TRUE
    
    cat("\n═══ 1. DIMENSION SUMMARY ═══\n")
    cat(sprintf("  comp_bin:      %d bins × %d samples\n", nrow(comp_bin), ncol(comp_bin)))
    cat(sprintf("  mat_chr_rot:   %d samples × %d bins\n", nrow(mat_chr_rot), ncol(mat_chr_rot)))
    cat(sprintf("  ord_bins_ids:  %d elements\n", length(ord_bins_ids)))
    cat(sprintf("  ecg_split:     %d elements, %d levels\n", length(ecg_split), nlevels(ecg_split)))
    cat(sprintf("  row_split_age: %d elements, %d levels\n", length(row_split_age), nlevels(row_split_age)))
    cat(sprintf("  row_order_idx: %d elements\n", length(row_order_idx)))
    
    cat("\n═══ 2. CRITICAL ALIGNMENT CHECKS ═══\n")
    
    p1 <- nrow(mat_chr_rot) == ncol(comp_bin)
    p2 <- ncol(mat_chr_rot) == nrow(comp_bin)
    p3 <- identical(colnames(mat_chr_rot), ord_bins_ids)
    p4 <- length(ecg_split) == ncol(mat_chr_rot)
    p5 <- length(row_split_age) == nrow(mat_chr_rot)
    p6 <- length(row_order_idx) == nrow(mat_chr_rot)
    p7 <- all(row_order_idx >= 1 & row_order_idx <= nrow(mat_chr_rot))
    p8 <- !any(duplicated(row_order_idx))
    
    cat(sprintf("  Transpose correct (rows=samples): %s\n", if(p1) "✅" else "❌"))
    cat(sprintf("  Transpose correct (cols=bins):    %s\n", if(p2) "✅" else "❌"))
    cat(sprintf("  colnames = ord_bins_ids:          %s\n", if(p3) "✅" else "❌ CRITICAL!"))
    cat(sprintf("  ecg_split aligned:                %s\n", if(p4) "✅" else "❌"))
    cat(sprintf("  row_split_age aligned:            %s\n", if(p5) "✅" else "❌"))
    cat(sprintf("  row_order_idx length:             %s\n", if(p6) "✅" else "❌"))
    cat(sprintf("  row_order_idx valid range:        %s\n", if(p7) "✅" else "❌"))
    cat(sprintf("  row_order_idx unique:             %s\n", if(p8) "✅" else "❌"))
    
    all_pass <- all(p1, p2, p3, p4, p5, p6, p7, p8)
    
    cat("\n═══ 3. DATA INTEGRITY ═══\n")
    
    frac_active <- mean(comp_bin == 1L, na.rm = TRUE)
    cat(sprintf("  Active compartment fraction: %.1f%%\n", 100 * frac_active))
    
    frac_active_all <- rowMeans(comp_bin == 1L, na.rm = TRUE)
    frac_active_ordered <- frac_active_all[ord_bins_ids]
    cluster_activity <- tapply(frac_active_ordered, ecg_split, mean, na.rm = TRUE)
    
    cat("\n  ECG Cluster Activity (should decrease):\n")
    for (i in seq_along(cluster_activity)) {
        bar <- paste(rep("█", round(cluster_activity[i] * 20)), collapse = "")
        cat(sprintf("    %2s: %s %.3f\n", names(cluster_activity)[i], bar, cluster_activity[i]))
    }
    
    cat("\n═══ 4. METADATA COVERAGE ═══\n")
    
    cat(sprintf("  Age groups: %s\n", paste(names(table(age_vec)), collapse = ", ")))
    cat(sprintf("  Sex groups: %s\n", paste(names(table(sex_vec)), collapse = ", ")))
    cat(sprintf("  Cell types: %d unique\n", length(unique(ctype_vec))))
    
    cat("\n═══ SUMMARY ═══\n")
    if (all_pass) {
        cat("\n  ╔════════════════════════════════════════╗\n")
        cat("  ║  ✅ ALL CHECKS PASSED                  ║\n")
        cat("  ╚════════════════════════════════════════╝\n")
    } else {
        cat("\n  ╔════════════════════════════════════════╗\n")
        cat("  ║  ❌ SOME CHECKS FAILED                 ║\n")
        cat("  ╚════════════════════════════════════════╝\n")
    }
    
    cat(sprintf("\n  Log saved: %s\n\n", log_file))
    sink()
    invisible(all_pass)
}


## =========================================================
## STEP 1: LOAD MANIFEST
## =========================================================

banner("STEP 1: Loading manifest")

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

sanity_check_manifest(man)


## =========================================================
## STEP 2: CREATE GENOMIC BINS + BLACKLIST + GAP FILTER
## =========================================================

banner("STEP 2: Creating genomic bins and filtering")

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

sanity_check_bins(bins_use, bins_df)


## =========================================================
## STEP 3: SUMMARIZE BIGWIGS → Z-SCORE MATRIX
## =========================================================

banner("STEP 3: Extracting signal from bigWigs")

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
colnames(mat) <- names(bw_vec); rownames(mat) <- bins_df$bin_id
mat_z <- scale(mat); mat_z[is.na(mat_z)] <- 0
saveRDS(mat_z, file.path(CONFIG$outdir, "matrix_signal_z.rds"))

sanity_check_signal_matrix(mat_z, bins_df, man)


## =========================================================
## STEP 4: HMM SEGMENTATION (A/B COMPARTMENTS)
## =========================================================

banner("STEP 4: Running HMM segmentation per chromosome")

joint_hmm_one_chr <- function(x_mat_chr) {
    n_bins   <- nrow(x_mat_chr); n_groups <- ncol(x_mat_chr)
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
colnames(comp_bin) <- colnames(mat_z); rownames(comp_bin) <- rownames(mat_z)
saveRDS(comp_bin, file.path(CONFIG$outdir, "compartments_binary_A1_R0.rds"))

sanity_check_compartments(comp_bin, mat_z)


## =========================================================
## STEP 5: CHROMHMM OVERLAP → PER-BIN LOG2FE
## =========================================================

banner("STEP 5: Computing ChromHMM state enrichments")

ol2 <- findOverlaps(bins_use, fs_gr, ignore.strand = TRUE)
q2  <- queryHits(ol2); s2 <- subjectHits(ol2)
ov2 <- width(pintersect(bins_use[q2], fs_gr[s2]))
cov_dt_mn <- data.table(bin = q2, mn = as.character(mcols(fs_gr)$mnemonic[s2]), w = ov2)
cov_wide_mn <- dcast(cov_dt_mn, bin ~ mn, value.var = "w", fill = 0L, fun.aggregate = sum)
bw_use <- width(bins_use); cov_wide_mn[, total := bw_use[bin]]
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
saveRDS(log2fe_bin, file.path(CONFIG$outdir, sprintf("log2fe_bin.nonGap%02d.rds", CONFIG$non_gap_thresh*100)))

message("  ✅ ChromHMM enrichments computed")


## =========================================================
## STEP 6: CPG DENSITY & AT CONTENT
## =========================================================

banner("STEP 6: Computing CpG density and AT content")

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


## =========================================================
## STEP 7: K-MEANS CLUSTERING (ECG GROUPS)
## =========================================================

banner("STEP 7: Running k-means clustering")

km <- kmeans(comp_bin, centers = CONFIG$k_bins, nstart = 25)
groups_kept <- km$cluster
names(groups_kept) <- rownames(comp_bin)

cluster_assignments <- data.frame(
    bin_id = rownames(comp_bin),
    ecg_cluster = groups_kept,
    stringsAsFactors = FALSE
)
fwrite(cluster_assignments, file.path(CONFIG$outdir, "kmeans_cluster_assignments.tsv"), sep = "\t")

sanity_check_kmeans(groups_kept, comp_bin, CONFIG$k_bins)


## =========================================================
## STEP 8: BIN ORDERING BY ECG AND ACTIVITY
## =========================================================

banner("STEP 8: Ordering bins by ECG cluster and activity")

# 8A: Compute per-bin active fraction
frac_active_all <- rowMeans(comp_bin == 1L, na.rm = TRUE)

# 8B: Initial ordering
ord_idx <- order(groups_kept, -frac_active_all, na.last = TRUE)
ord_bins_ids <- rownames(comp_bin)[ord_idx]
ecg_split <- factor(groups_kept[ord_bins_ids])
ecg_split <- droplevels(ecg_split)

# 8C: Re-order clusters by mean activity
cluster_activity <- tapply(frac_active_all[ord_bins_ids], ecg_split, mean)
cluster_order <- order(cluster_activity, decreasing = TRUE)
ecg_levels_ordered <- levels(ecg_split)[cluster_order]
ecg_split <- factor(ecg_split, levels = ecg_levels_ordered)

ord_idx2 <- order(as.integer(ecg_split), -frac_active_all[ord_bins_ids])
ord_bins_ids <- ord_bins_ids[ord_idx2]
ecg_split <- ecg_split[ord_idx2]

# 8D: Build rotated matrix
comp_ord <- comp_bin[ord_bins_ids, , drop = FALSE]
mat_chr_rot <- t(apply(comp_ord, 2, function(x) ifelse(x == 1L, "1", "0")))
colnames(mat_chr_rot) <- ord_bins_ids

sanity_check_bin_ordering(ord_bins_ids, ecg_split, comp_bin)
sanity_check_rotation(mat_chr_rot, comp_bin, ord_bins_ids)


## =========================================================
## STEP 9: RELABEL ECG CLUSTERS CONSECUTIVELY
## =========================================================

banner("STEP 9: Relabeling ECG clusters")

present_levels <- levels(ecg_split)
n_clusters <- length(present_levels)
new_labels <- as.character(1:n_clusters)
ecg_remap <- setNames(new_labels, present_levels)
ecg_split_consecutive <- factor(ecg_remap[as.character(ecg_split)], levels = new_labels)

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

ecg_split <- ecg_split_consecutive

message(sprintf("  ✅ ECG clusters relabeled: 1 → %d", n_clusters))

# Save diagnostic table
frac_dt <- data.frame(
    bin_id = ord_bins_ids,
    chr = bins_df$chr[match(ord_bins_ids, bins_df$bin_id)],
    ecg_cluster = as.character(ecg_split),
    frac_active = frac_active_all[ord_bins_ids],
    stringsAsFactors = FALSE
)
fwrite(frac_dt, file.path(CONFIG$outdir, "bin_ordering_by_ECG_and_activity.tsv"), sep = "\t")


## =========================================================
## STEP 10: BUILD ANNOTATION MATRICES
## =========================================================

banner("STEP 10: Building annotation matrices")

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

# Column annotations (stacked barplots)
age_track      <- per_bin_active_strict(comp_bin, age_vec,   AGE_LEVELS,           ord_bins_ids)
sex_track      <- per_bin_active_strict(comp_bin, sex_vec,   SEX_LEVELS,           ord_bins_ids)
celltype_track <- per_bin_active_strict(comp_bin, ctype_vec, names(CELLTYPE_COLS), ord_bins_ids)

age_mat_for_anno      <- t(age_track)
sex_mat_for_anno      <- t(sex_track)
celltype_mat_for_anno <- t(celltype_track)

sanity_check_annotation_matrices(age_mat_for_anno, sex_mat_for_anno, 
                                  celltype_mat_for_anno, mat_chr_rot)

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

message("  ✅ Annotations built")


## =========================================================
## STEP 11: BUILD OPTIONAL TRACKS
## =========================================================

banner("STEP 11: Building optional tracks")

# Chromosome track
chr_vec   <- bins_df$chr[match(ord_bins_ids, bins_df$bin_id)]
chr_mat   <- matrix(chr_vec, nrow = 1, dimnames = list("Chromosome", ord_bins_ids))
ht_chr <- Heatmap(chr_mat, name = "Chromosome", col = CHR_COLS,
                  cluster_rows = FALSE, cluster_columns = FALSE,
                  show_row_names = TRUE, show_column_names = FALSE,
                  height = unit(0.6,"cm"))

# CpG / AT tracks
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

message("  ✅ Optional tracks built")


## =========================================================
## STEP 12: CHROMHMM PANEL
## =========================================================

banner("STEP 12: Building ChromHMM panel")

log2fe_bin <- readRDS(file.path(CONFIG$outdir, sprintf("log2fe_bin.nonGap%02d.rds", CONFIG$non_gap_thresh * 100)))
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

message("  ✅ ChromHMM panel built")


## =========================================================
## STEP 13: ROW ORDERING
## =========================================================

banner("STEP 13: Computing row ordering")

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

sanity_check_row_ordering(row_order_idx, row_split_age, mat_chr_rot)


## =========================================================
## STEP 14: FINAL VALIDATION
## =========================================================

sanity_check_heatmap_final(mat_chr_rot, row_split_age, row_order_idx, ecg_split)

message("[DIAGNOSTICS] Running full diagnostic report...")
run_full_diagnostics(comp_bin, mat_chr_rot, ord_bins_ids, ecg_split,
                     row_split_age, row_order_idx, age_vec, sex_vec,
                     ctype_vec, bins_df, CONFIG$outdir)


## =========================================================
## STEP 15: CREATE AND SAVE HEATMAP
## =========================================================

banner("STEP 15: Assembling and saving heatmap")

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

# Stack all tracks
ht_stack <- ht_main %v% ht_chr %v% ht_cpg %v% ht_at %v% ht_states

# Save outputs
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


## =========================================================
## FINAL SUMMARY
## =========================================================

cat("\n")
cat("╔══════════════════════════════════════════════════════════════════╗\n")
cat("║                    PIPELINE COMPLETE                             ║\n")
cat("╠══════════════════════════════════════════════════════════════════╣\n")
cat(sprintf("║  Output directory: %-45s║\n", substr(CONFIG$outdir, 1, 45)))
cat("║                                                                  ║\n")
cat("║  Files generated:                                                ║\n")
cat("║    • compartment_heatmap_full.pdf                                ║\n")
cat("║    • compartment_heatmap_full.png                                ║\n")
cat("║    • compartment_diagnostics_full.log                            ║\n")
cat("║    • compartments_binary_A1_R0.rds                               ║\n")
cat("║    • matrix_signal_z.rds                                         ║\n")
cat("║    • kmeans_cluster_assignments.tsv                              ║\n")
cat("║    • bin_ordering_by_ECG_and_activity.tsv                        ║\n")
cat("╚══════════════════════════════════════════════════════════════════╝\n")

message("\n[INFO] Next step: Run compartment_02_stability_analysis.R")
