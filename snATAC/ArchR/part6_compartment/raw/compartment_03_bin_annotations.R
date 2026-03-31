#!/usr/bin/env Rscript
## =========================================================
## COMPARTMENT ANALYSIS - BIN-WISE ANNOTATION EXTRACTION
## 
## Consolidates ALL bin-level information into a single table
## for downstream annotation with ChIPseeker, GREAT, etc.
##
## Prerequisites:
##   - Source compartment_00_config.R
##   - Run compartment_01_main_analysis.R
##   - Run compartment_02_stability_analysis.R
##
## Output:
##   - Master table with columns for each analysis layer
##   - BED files for different subsets
##   - GRanges object for ChIPseeker
##
## =========================================================

box_banner("BIN-WISE ANNOTATION EXTRACTION")

## =========================================================
## LOAD DEPENDENCIES
## =========================================================

if (!exists("CONFIG")) {
    source("compartment_00_config.R")
}

suppressPackageStartupMessages({
    library(data.table)
    library(GenomicRanges)
})

# Check required objects
check_required_objects(c("comp_bin", "bins_df", "age_vec", "sex_vec", "ctype_vec",
                         "stability_ext_dt", "ord_bins_ids", "ecg_split", "groups_kept"))


## =========================================================
## SECTION 1: INITIALIZE MASTER BIN TABLE
## =========================================================

banner("STEP 1: Initializing master bin table")

master_bins <- data.table(
    bin_id = bins_df$bin_id,
    chr    = bins_df$chr,
    start  = bins_df$start,
    end    = bins_df$end
)

master_bins[, width := end - start + 1]
master_bins[, midpoint := as.integer((start + end) / 2)]

message(sprintf("  ✅ Initialized with %d bins", nrow(master_bins)))


## =========================================================
## SECTION 2: ADD ECG CLUSTER INFORMATION
## =========================================================

banner("STEP 2: Adding ECG cluster assignments")

ecg_dt <- data.table(
    bin_id = ord_bins_ids,
    ecg_cluster = as.character(ecg_split)
)
ecg_dt[, ecg_cluster_original := groups_kept[bin_id]]

master_bins <- merge(master_bins, ecg_dt, by = "bin_id", all.x = TRUE)
message(sprintf("  ✅ Added ECG clusters (1-%d)", nlevels(ecg_split)))


## =========================================================
## SECTION 3: ADD COMPARTMENT ACTIVITY METRICS
## =========================================================

banner("STEP 3: Adding compartment activity metrics")

frac_active_all <- rowMeans(comp_bin == 1L, na.rm = TRUE)
master_bins[, frac_active_all := frac_active_all[bin_id]]

# Active fraction by sex
frac_by_sex <- function(comp_bin, sex_vec, sex_val) {
    idx <- which(sex_vec == sex_val)
    if (length(idx) == 0) return(rep(NA_real_, nrow(comp_bin)))
    rowMeans(comp_bin[, idx, drop = FALSE] == 1L, na.rm = TRUE)
}

master_bins[, frac_active_male := frac_by_sex(comp_bin, sex_vec, "male")[bin_id]]
master_bins[, frac_active_female := frac_by_sex(comp_bin, sex_vec, "female")[bin_id]]
master_bins[, sex_diff_active := frac_active_male - frac_active_female]

message("  ✅ Added activity fractions (all, male, female, sex_diff)")


## =========================================================
## SECTION 4: ADD AGE-SPECIFIC ACTIVITY
## =========================================================

banner("STEP 4: Adding age-specific activity fractions")

frac_by_age <- frac_active_by(comp_bin, age_vec, AGE_LEVELS)

for (age in AGE_LEVELS) {
    col_name <- paste0("frac_active_", age)
    master_bins[, (col_name) := frac_by_age[bin_id, age]]
}

# Age trajectory slope
master_bins[, age_slope := {
    y <- c(frac_active_young, frac_active_mid_age, frac_active_old, 
           frac_active_pre_geriatric, frac_active_geriatric)
    x <- 1:5
    if (any(is.na(y))) NA_real_ else coef(lm(y ~ x))[2]
}, by = bin_id]

master_bins[, age_trend := fcase(
    age_slope > 0.05, "Increasing",
    age_slope < -0.05, "Decreasing",
    default = "Stable"
)]

message(sprintf("  ✅ Added %d age columns + slope + trend", length(AGE_LEVELS)))


## =========================================================
## SECTION 5: ADD STABILITY CLASSIFICATIONS
## =========================================================

banner("STEP 5: Adding stability classifications")

# Pivot stability data to wide format
stability_wide <- dcast(
    stability_ext_dt,
    bin_id ~ Sex,
    value.var = "Stability",
    fun.aggregate = function(x) x[1]
)
setnames(stability_wide, 
         old = c("male", "female"),
         new = c("stability_male", "stability_female"),
         skip_absent = TRUE)

master_bins <- merge(master_bins, stability_wide, by = "bin_id", all.x = TRUE)

# Consensus stability
master_bins[, stability_consensus := fcase(
    stability_male == stability_female, stability_male,
    is.na(stability_male) & !is.na(stability_female), stability_female,
    !is.na(stability_male) & is.na(stability_female), stability_male,
    default = "Sex_Discordant"
)]

# Simplified 3-class
master_bins[, stability_3class := fcase(
    stability_consensus %in% c("Stable_Active"), "Stable_Active",
    stability_consensus %in% c("Stable_Repressive"), "Stable_Repressive",
    stability_consensus %in% c("Monotonic_A_to_R", "Monotonic_R_to_A", "Non_Monotonic"), "Switching",
    default = stability_consensus
)]

message("  ✅ Added stability (male, female, consensus, 3-class)")


## =========================================================
## SECTION 6: ADD NON-MONOTONIC PATTERN CLASSIFICATION
## =========================================================

banner("STEP 6: Adding non-monotonic pattern classifications")

if (exists("state_with_stability") && exists("nonmono_wide") && nrow(nonmono_wide) > 0) {
    
    # Aggregate pattern by bin/sex
    pattern_by_bin_sex <- nonmono_wide[, .(
        nonmono_pattern = names(sort(table(pattern), decreasing = TRUE))[1],
        n_transitions_mean = mean(n_transitions),
        trajectory_example = trajectory[1]
    ), by = .(bin_id, Sex)]
    
    # Pivot to wide
    pattern_male <- pattern_by_bin_sex[Sex == "male", .(
        bin_id, 
        nonmono_pattern_male = nonmono_pattern,
        n_transitions_male = n_transitions_mean,
        trajectory_male = trajectory_example
    )]
    
    pattern_female <- pattern_by_bin_sex[Sex == "female", .(
        bin_id,
        nonmono_pattern_female = nonmono_pattern,
        n_transitions_female = n_transitions_mean,
        trajectory_female = trajectory_example
    )]
    
    master_bins <- merge(master_bins, pattern_male, by = "bin_id", all.x = TRUE)
    master_bins <- merge(master_bins, pattern_female, by = "bin_id", all.x = TRUE)
    
    # Consensus pattern
    master_bins[, nonmono_pattern_consensus := fcase(
        is.na(nonmono_pattern_male) & is.na(nonmono_pattern_female), NA_character_,
        is.na(nonmono_pattern_male), nonmono_pattern_female,
        is.na(nonmono_pattern_female), nonmono_pattern_male,
        nonmono_pattern_male == nonmono_pattern_female, nonmono_pattern_male,
        default = "Sex_Discordant"
    )]
    
    message(sprintf("  ✅ Added non-monotonic patterns for %d bins", 
                    sum(!is.na(master_bins$nonmono_pattern_consensus))))
} else {
    message("  [INFO] No non-monotonic data available")
}


## =========================================================
## SECTION 7: ADD SEQUENCE FEATURES
## =========================================================

banner("STEP 7: Adding sequence features")

if (exists("gc_dt")) {
    gc_data <- data.table(
        bin_id = rownames(gc_dt),
        CpG_density = gc_dt$CpG_density,
        AT_content = gc_dt$AT_content
    )
    gc_data[, GC_content := 1 - AT_content]
    
    cpg_quantiles <- quantile(gc_data$CpG_density, c(0.25, 0.75), na.rm = TRUE)
    gc_data[, CpG_category := fcase(
        CpG_density <= cpg_quantiles[1], "Low_CpG",
        CpG_density >= cpg_quantiles[2], "High_CpG",
        default = "Medium_CpG"
    )]
    
    master_bins <- merge(master_bins, gc_data, by = "bin_id", all.x = TRUE)
    message("  ✅ Added CpG_density, AT_content, GC_content, CpG_category")
} else {
    message("  [INFO] gc_dt not found - skipping sequence features")
}


## =========================================================
## SECTION 8: ADD CHROMHMM STATE ENRICHMENTS
## =========================================================

banner("STEP 8: Adding ChromHMM state enrichments")

if (exists("log2fe_bin")) {
    chromhmm_dt <- as.data.table(log2fe_bin, keep.rownames = "bin_id")
    old_names <- setdiff(names(chromhmm_dt), "bin_id")
    new_names <- paste0("chromHMM_", old_names)
    setnames(chromhmm_dt, old_names, new_names)
    
    master_bins <- merge(master_bins, chromhmm_dt, by = "bin_id", all.x = TRUE)
    
    # Dominant ChromHMM state
    chromhmm_cols <- grep("^chromHMM_", names(master_bins), value = TRUE)
    master_bins[, chromHMM_dominant := {
        vals <- unlist(.SD)
        if (all(is.na(vals))) NA_character_ else gsub("^chromHMM_", "", chromhmm_cols[which.max(vals)])
    }, .SDcols = chromhmm_cols, by = bin_id]
    
    message(sprintf("  ✅ Added %d ChromHMM state enrichments + dominant state", length(old_names)))
} else {
    message("  [INFO] log2fe_bin not found - skipping ChromHMM")
}


## =========================================================
## SECTION 9: ADD SIGNAL STATISTICS
## =========================================================

banner("STEP 9: Adding signal statistics")

if (exists("mat_z")) {
    signal_stats <- data.table(
        bin_id = rownames(mat_z),
        signal_mean = rowMeans(mat_z, na.rm = TRUE),
        signal_sd = apply(mat_z, 1, sd, na.rm = TRUE),
        signal_median = apply(mat_z, 1, median, na.rm = TRUE),
        signal_min = apply(mat_z, 1, min, na.rm = TRUE),
        signal_max = apply(mat_z, 1, max, na.rm = TRUE)
    )
    signal_stats[, signal_cv := signal_sd / abs(signal_mean)]
    signal_stats[!is.finite(signal_cv), signal_cv := NA_real_]
    
    master_bins <- merge(master_bins, signal_stats, by = "bin_id", all.x = TRUE)
    message("  ✅ Added signal statistics (mean, sd, median, min, max, cv)")
} else {
    message("  [INFO] mat_z not found - skipping signal stats")
}


## =========================================================
## SECTION 10: ADD COMPARTMENT SWITCH METRICS
## =========================================================

banner("STEP 10: Adding compartment switch metrics")

master_bins[, n_samples_active := rowSums(comp_bin[bin_id, ] == 1L, na.rm = TRUE)]
master_bins[, n_samples_repressive := rowSums(comp_bin[bin_id, ] == 0L, na.rm = TRUE)]
master_bins[, n_samples_total := ncol(comp_bin)]

# Compartment entropy
master_bins[, compartment_entropy := {
    p_active <- frac_active_all
    p_repressive <- 1 - p_active
    ifelse(p_active == 0 | p_active == 1, 0,
           -p_active * log2(p_active) - p_repressive * log2(p_repressive))
}]

master_bins[, compartment_variability := fcase(
    compartment_entropy < 0.5, "Low",
    compartment_entropy < 0.9, "Medium",
    default = "High"
)]

message("  ✅ Added sample counts and entropy metrics")


## =========================================================
## SECTION 11: CREATE GENOMICRANGES
## =========================================================

banner("STEP 11: Creating GenomicRanges object")

master_gr <- GRanges(
    seqnames = master_bins$chr,
    ranges = IRanges(start = master_bins$start, end = master_bins$end),
    strand = "*"
)

meta_cols <- setdiff(names(master_bins), c("chr", "start", "end"))
mcols(master_gr) <- master_bins[, ..meta_cols]

message(sprintf("  ✅ Created GRanges with %d bins and %d metadata columns",
                length(master_gr), ncol(mcols(master_gr))))


## =========================================================
## SECTION 12: CREATE BED FILES
## =========================================================

banner("STEP 12: Creating BED files for annotation tools")

bed_dir <- file.path(CONFIG$outdir, "bed_files_for_annotation")
ensure_dir(bed_dir)

# All bins
write_bed(master_bins, file.path(bed_dir, "all_bins.bed"))

# By stability class
if ("stability_consensus" %in% names(master_bins)) {
    for (stab in unique(na.omit(master_bins$stability_consensus))) {
        stab_clean <- clean_name(stab)
        subset_dt <- master_bins[stability_consensus == stab]
        if (nrow(subset_dt) > 0) {
            write_bed(subset_dt, file.path(bed_dir, sprintf("bins_%s.bed", stab_clean)))
        }
    }
}

# By ECG cluster
for (ecg in unique(na.omit(master_bins$ecg_cluster))) {
    subset_dt <- master_bins[ecg_cluster == ecg]
    write_bed(subset_dt, file.path(bed_dir, sprintf("bins_ECG_%s.bed", ecg)))
}

# High variability bins
if ("compartment_variability" %in% names(master_bins)) {
    high_var <- master_bins[compartment_variability == "High"]
    if (nrow(high_var) > 0) {
        write_bed(high_var, file.path(bed_dir, "bins_high_variability.bed"))
    }
}

# Sex-biased bins
if ("sex_diff_active" %in% names(master_bins)) {
    male_biased <- master_bins[sex_diff_active > 0.2]
    female_biased <- master_bins[sex_diff_active < -0.2]
    if (nrow(male_biased) > 0) write_bed(male_biased, file.path(bed_dir, "bins_male_biased.bed"))
    if (nrow(female_biased) > 0) write_bed(female_biased, file.path(bed_dir, "bins_female_biased.bed"))
}

# Age-trending bins
if ("age_trend" %in% names(master_bins)) {
    increasing <- master_bins[age_trend == "Increasing"]
    decreasing <- master_bins[age_trend == "Decreasing"]
    if (nrow(increasing) > 0) write_bed(increasing, file.path(bed_dir, "bins_age_increasing.bed"))
    if (nrow(decreasing) > 0) write_bed(decreasing, file.path(bed_dir, "bins_age_decreasing.bed"))
}

# Non-monotonic patterns
if ("nonmono_pattern_consensus" %in% names(master_bins)) {
    for (pattern in NONMONO_PATTERNS) {
        subset_dt <- master_bins[nonmono_pattern_consensus == pattern]
        if (nrow(subset_dt) > 0) {
            write_bed(subset_dt, file.path(bed_dir, sprintf("bins_nonmono_%s.bed", pattern)))
        }
    }
}

message(sprintf("  ✅ BED files saved to: %s", bed_dir))


## =========================================================
## SECTION 13: SAVE MASTER TABLE
## =========================================================

banner("STEP 13: Saving master annotation table")

master_file <- file.path(CONFIG$outdir, "master_bin_annotations.tsv")
fwrite(master_bins, master_file, sep = "\t")
message(sprintf("  ✅ TSV: %s", master_file))

master_rds <- file.path(CONFIG$outdir, "master_bin_annotations.rds")
saveRDS(master_bins, master_rds)
message(sprintf("  ✅ RDS: %s", master_rds))

gr_file <- file.path(CONFIG$outdir, "master_bin_annotations.GRanges.rds")
saveRDS(master_gr, gr_file)
message(sprintf("  ✅ GRanges: %s", gr_file))


## =========================================================
## SECTION 14: COLUMN DICTIONARY
## =========================================================

banner("STEP 14: Creating column dictionary")

col_dict <- data.table(
    column = names(master_bins),
    type = sapply(master_bins, class),
    n_unique = sapply(master_bins, function(x) length(unique(na.omit(x)))),
    n_na = sapply(master_bins, function(x) sum(is.na(x))),
    description = ""
)

# Add descriptions
col_dict[column == "bin_id", description := "Unique bin identifier (chr:start-end)"]
col_dict[column == "chr", description := "Chromosome"]
col_dict[column == "start", description := "Start position (1-based)"]
col_dict[column == "end", description := "End position"]
col_dict[column == "width", description := "Bin width in bp"]
col_dict[column == "ecg_cluster", description := "ECG cluster assignment (1-25, ordered by activity)"]
col_dict[column == "frac_active_all", description := "Fraction of all samples where bin is active"]
col_dict[column == "frac_active_male", description := "Fraction active in male samples"]
col_dict[column == "frac_active_female", description := "Fraction active in female samples"]
col_dict[column == "sex_diff_active", description := "Male - Female activity difference"]
col_dict[column == "age_slope", description := "Linear slope of activity across ages"]
col_dict[column == "age_trend", description := "Categorical age trend (Increasing/Decreasing/Stable)"]
col_dict[column == "stability_male", description := "5-class stability in males"]
col_dict[column == "stability_female", description := "5-class stability in females"]
col_dict[column == "stability_consensus", description := "Consensus stability"]
col_dict[column == "stability_3class", description := "Simplified 3-class stability"]
col_dict[column == "compartment_entropy", description := "Binary entropy of compartment state (0-1)"]
col_dict[column == "compartment_variability", description := "Variability category (Low/Medium/High)"]

fwrite(col_dict, file.path(CONFIG$outdir, "master_bin_column_dictionary.tsv"), sep = "\t")
message(sprintf("  ✅ Column dictionary: %d columns documented", nrow(col_dict)))


## =========================================================
## FINAL SUMMARY
## =========================================================

cat("\n")
cat("╔══════════════════════════════════════════════════════════════════╗\n")
cat("║              BIN-WISE EXTRACTION COMPLETE                        ║\n")
cat("╠══════════════════════════════════════════════════════════════════╣\n")
cat(sprintf("║  Total bins: %s                                           ║\n", 
            format(nrow(master_bins), big.mark = ",")))
cat(sprintf("║  Columns: %d                                                     ║\n",
            ncol(master_bins)))
cat("║                                                                  ║\n")
cat("║  Files created:                                                  ║\n")
cat("║    • master_bin_annotations.tsv                                  ║\n")
cat("║    • master_bin_annotations.rds                                  ║\n")
cat("║    • master_bin_annotations.GRanges.rds                          ║\n")
cat("║    • master_bin_column_dictionary.tsv                            ║\n")
cat("║    • bed_files_for_annotation/                                   ║\n")
cat("╚══════════════════════════════════════════════════════════════════╝\n")

message("\n[INFO] Next step: Run compartment_04_alluvial_plots.R")
