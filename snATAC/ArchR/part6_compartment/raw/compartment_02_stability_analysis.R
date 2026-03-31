#!/usr/bin/env Rscript
## =========================================================
## COMPARTMENT ANALYSIS - STABILITY CLASSIFICATION
## 
## Extended 5-Class Model:
##   1. Stable_Active       - Consistently active across aging
##   2. Stable_Repressive   - Consistently repressive across aging
##   3. Monotonic_A_to_R    - Clean closing (no reversals)
##   4. Monotonic_R_to_A    - Clean opening (no reversals)
##   5. Non_Monotonic       - Oscillating/complex patterns
##
## Prerequisites:
##   - Source compartment_00_config.R
##   - Run compartment_01_main_analysis.R (need comp_bin, age_vec, sex_vec, ctype_vec)
##
## Output:
##   - stability_ext_dt: Stability classifications per bin/sex/celltype
##   - state_per_age: State data per sample
##   - state_with_stability: Merged data for plotting
##   - Barplots and summary tables
##
## =========================================================

box_banner("COMPARTMENT STABILITY - 5-CLASS MODEL")

## =========================================================
## LOAD DEPENDENCIES
## =========================================================

# Source configuration (if not already loaded)
if (!exists("CONFIG")) {
    source("compartment_00_config.R")
}

# Load required libraries
suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(ggalluvial)
})

# Check required objects from main analysis
check_required_objects(c("comp_bin", "age_vec", "sex_vec", "ctype_vec", "CONFIG"))


## =========================================================
## SECTION 1: STABILITY CLASSIFICATION FUNCTIONS
## =========================================================

#' Check if trajectory is monotonic (no reversals)
#' @param states character vector of "Active"/"Repressive" across ages
#' @return list with: is_monotonic (logical), direction (character), n_transitions (integer)
check_monotonicity <- function(states) {
    
    # Remove NAs
    states <- states[!is.na(states)]
    if (length(states) < 2) {
        return(list(is_monotonic = NA, direction = "Insufficient", n_transitions = NA))
    }
    
    # Count transitions
    transitions <- states[-1] != states[-length(states)]
    n_trans <- sum(transitions)
    
    if (n_trans == 0) {
        return(list(is_monotonic = TRUE, direction = "Stable", n_transitions = 0))
    }
    
    # Convert to numeric: Active = 1, Repressive = 0
    numeric_states <- ifelse(states == "Active", 1L, 0L)
    
    # Check if monotonically non-increasing (A→R direction)
    is_monotonic_AR <- all(diff(numeric_states) <= 0) & any(diff(numeric_states) < 0)
    
    # Check if monotonically non-decreasing (R→A direction)
    is_monotonic_RA <- all(diff(numeric_states) >= 0) & any(diff(numeric_states) > 0)
    
    if (is_monotonic_AR) {
        return(list(is_monotonic = TRUE, direction = "A_to_R", n_transitions = n_trans))
    } else if (is_monotonic_RA) {
        return(list(is_monotonic = TRUE, direction = "R_to_A", n_transitions = n_trans))
    } else {
        return(list(is_monotonic = FALSE, direction = "Non_Monotonic", n_transitions = n_trans))
    }
}


#' Extended 5-class stability classifier
#' @param fracs numeric vector of fraction active per age
#' @param t_on threshold for "Active" (default 0.6)
#' @param t_off threshold for "Repressive" (default 0.4)
#' @return character: one of 5 classes
classify_stability_extended <- function(fracs, t_on = 0.6, t_off = 0.4) {
    
    if (all(is.na(fracs))) return("Missing")
    
    n_valid <- sum(!is.na(fracs))
    if (n_valid < 3) return("Insufficient_Data")
    
    # Convert fractions to states
    states <- rep("Intermediate", length(fracs))
    states[fracs >= t_on] <- "Active"
    states[fracs <= t_off] <- "Repressive"
    states[is.na(fracs)] <- NA
    
    # Remove intermediate for stability check
    valid_states <- states[states %in% c("Active", "Repressive")]
    
    if (length(valid_states) < 3) return("Insufficient_Data")
    
    # Check stability first (≥80% same state)
    n_active <- sum(valid_states == "Active")
    n_repressive <- sum(valid_states == "Repressive")
    n_total <- length(valid_states)
    
    if (n_active / n_total >= 0.8) return("Stable_Active")
    if (n_repressive / n_total >= 0.8) return("Stable_Repressive")
    
    # Not stable, so it's switching - check monotonicity
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


#' Compute extended stability for all bins
#' @param comp_bin Compartment binary matrix
#' @param age_vec Named vector of ages
#' @param sex_vec Named vector of sexes  
#' @param ctype_vec Named vector of celltypes
#' @param age_levels Age level order
#' @param t_on Active threshold
#' @param t_off Repressive threshold
#' @return data.table with stability per bin/sex/celltype
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
            
            message(sprintf("[STABILITY] %s - %s (%d samples)", sx, ct, length(use)))
            
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


## =========================================================
## SECTION 2: COMPUTE STABILITY CLASSIFICATIONS
## =========================================================

banner("STEP 1: Computing 5-class stability classifications")

stability_ext_dt <- compute_stability_extended(
    comp_bin, age_vec, sex_vec, ctype_vec, AGE_LEVELS,
    t_on = CONFIG$t_on, t_off = CONFIG$t_off
)

fwrite(stability_ext_dt,
       file.path(CONFIG$outdir, "compartment_stability_5class.tsv"),
       sep = "\t")

message("\n[INFO] Extended stability distribution:")
print(stability_ext_dt[, .N, by = Stability][order(-N)])


## =========================================================
## SECTION 3: CREATE 3-CLASS MAPPING
## =========================================================

banner("STEP 2: Creating 3-class mapping")

stability_ext_dt[, Stability_3class := fcase(
    Stability == "Stable_Active", "Stable_Active",
    Stability == "Stable_Repressive", "Stable_Repressive",
    Stability %in% c("Monotonic_A_to_R", "Monotonic_R_to_A", "Non_Monotonic"), "Switching",
    default = Stability
)]

# Comparison table
comparison_dt <- stability_ext_dt[, .N, by = .(Stability_3class, Stability)]
comparison_dt <- comparison_dt[order(Stability_3class, -N)]

message("\n[INFO] 3-class to 5-class breakdown:")
print(comparison_dt)

fwrite(comparison_dt,
       file.path(CONFIG$outdir, "stability_3class_vs_5class.tsv"),
       sep = "\t")


## =========================================================
## SECTION 4: PREPARE STATE-PER-AGE DATA
## =========================================================

banner("STEP 3: Preparing state data per age")

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

message(sprintf("  ✅ State data: %d rows", nrow(state_per_age)))


## =========================================================
## SECTION 5: MERGE STATE WITH STABILITY
## =========================================================

banner("STEP 4: Merging state data with stability")

state_with_stability <- merge(
    state_per_age,
    stability_ext_dt[, .(bin_id, Sex, Celltype, Stability)],
    by = c("bin_id", "Sex", "Celltype"),
    all.x = TRUE
)

message(sprintf("  ✅ Merged data: %d rows", nrow(state_with_stability)))


## =========================================================
## SECTION 6: NON-MONOTONIC TRAJECTORY ANALYSIS
## =========================================================

banner("STEP 5: Analyzing non-monotonic trajectories")

nonmono_dt <- state_with_stability[Stability == "Non_Monotonic"]

if (nrow(nonmono_dt) > 0) {
    
    nonmono_wide <- dcast(nonmono_dt, 
                          bin_id + Sex + Celltype ~ Age, 
                          value.var = "State")
    nonmono_wide <- nonmono_wide[complete.cases(nonmono_wide)]
    
    # Create trajectory string
    nonmono_wide[, trajectory := paste(young, mid_age, old, pre_geriatric, geriatric, sep = "→")]
    
    # Count transitions
    nonmono_wide[, n_transitions := {
        states <- c(young, mid_age, old, pre_geriatric, geriatric)
        sum(states[-1] != states[-length(states)])
    }, by = 1:nrow(nonmono_wide)]
    
    # Classify pattern
    nonmono_wide[, pattern := fcase(
        n_transitions == 2 & young == geriatric, "Transient_Return",
        n_transitions == 2, "Two_Switch",
        n_transitions == 3, "Triple_Switch",
        n_transitions >= 4, "Highly_Dynamic",
        default = "Other"
    )]
    
    # Summary
    message("\n[INFO] Non-monotonic pattern types:")
    print(nonmono_wide[, .N, by = pattern][order(-N)])
    
    # Top trajectories
    traj_counts <- nonmono_wide[, .N, by = .(Sex, Celltype, trajectory, pattern)]
    traj_counts <- traj_counts[order(Sex, Celltype, -N)]
    
    message("\n[INFO] Top 15 non-monotonic trajectories:")
    print(traj_counts[order(-N)][1:15])
    
    # Export
    fwrite(nonmono_wide, file.path(CONFIG$outdir, "nonmonotonic_trajectories_detail.tsv"), sep = "\t")
    fwrite(traj_counts, file.path(CONFIG$outdir, "nonmonotonic_trajectory_counts.tsv"), sep = "\t")
    
} else {
    message("  [INFO] No non-monotonic bins found")
    nonmono_wide <- data.table()
}


## =========================================================
## SECTION 7: STABILITY BARPLOTS
## =========================================================

banner("STEP 6: Creating stability barplots")

# Filter for plotting
plot_dt <- stability_ext_dt[!Stability %in% c("Missing", "Insufficient_Data")]

## --- 5-class barplot (counts) ---
p_5class_count <- ggplot(
    plot_dt[, .N, by = .(Sex, Celltype, Stability)],
    aes(x = Celltype, y = N, fill = Stability)
) +
    geom_col(position = "stack") +
    facet_wrap(~ Sex, ncol = 1) +
    scale_fill_manual(values = STABILITY_COLS) +
    labs(
        title = "Compartment Stability - Extended 5-Class Model",
        subtitle = "Monotonic switching distinguished from non-monotonic",
        y = "Number of 80 kb bins",
        x = "Cell Type",
        fill = "Stability Class"
    ) +
    theme_bw(base_size = 12) +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = 14, face = "bold")
    )

ggsave(file.path(CONFIG$outdir, "Barplot_Stability_5class_Counts.pdf"),
       p_5class_count, width = 12, height = 8)
message("[SAVED] Barplot_Stability_5class_Counts.pdf")

## --- 5-class barplot (proportions) ---
p_5class_prop <- ggplot(
    plot_dt[, .N, by = .(Sex, Celltype, Stability)],
    aes(x = Celltype, y = N, fill = Stability)
) +
    geom_col(position = "fill") +
    facet_wrap(~ Sex, ncol = 1) +
    scale_fill_manual(values = STABILITY_COLS) +
    labs(
        title = "Compartment Stability Proportions - Extended 5-Class Model",
        y = "Fraction of 80 kb bins",
        x = "Cell Type",
        fill = "Stability Class"
    ) +
    theme_bw(base_size = 12) +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = 14, face = "bold")
    )

ggsave(file.path(CONFIG$outdir, "Barplot_Stability_5class_Proportions.pdf"),
       p_5class_prop, width = 12, height = 8)
message("[SAVED] Barplot_Stability_5class_Proportions.pdf")

## --- Switching-only breakdown ---
switching_only <- plot_dt[Stability %in% SWITCHING_CLASSES]

p_switching_breakdown <- ggplot(
    switching_only[, .N, by = .(Sex, Celltype, Stability)],
    aes(x = Celltype, y = N, fill = Stability)
) +
    geom_col(position = "dodge") +
    facet_wrap(~ Sex, ncol = 1) +
    scale_fill_manual(values = STABILITY_COLS) +
    labs(
        title = "Switching Bins Breakdown: Monotonic vs Non-Monotonic",
        subtitle = "A-to-R (closing) | R-to-A (opening) | Non-monotonic (oscillating)",
        y = "Number of 80 kb bins",
        x = "Cell Type",
        fill = "Switch Type"
    ) +
    theme_bw(base_size = 12) +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = 14, face = "bold")
    )

ggsave(file.path(CONFIG$outdir, "Barplot_Switching_Breakdown.pdf"),
       p_switching_breakdown, width = 12, height = 8)
message("[SAVED] Barplot_Switching_Breakdown.pdf")

## --- Non-monotonic patterns barplot ---
if (nrow(nonmono_wide) > 0) {
    p_pattern <- ggplot(
        nonmono_wide[, .N, by = .(Sex, Celltype, pattern)],
        aes(x = Celltype, y = N, fill = pattern)
    ) +
        geom_col(position = "stack") +
        facet_wrap(~ Sex, ncol = 1) +
        scale_fill_manual(values = NONMONO_PATTERN_COLS) +
        labs(
            title = "Non-Monotonic Switching Patterns",
            subtitle = "Transient=returns to start | 2/3/4 switches = complexity levels",
            x = "Cell Type",
            y = "Number of 80kb Bins",
            fill = "Pattern"
        ) +
        theme_bw(base_size = 12) +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1),
            strip.text = element_text(size = 14, face = "bold")
        )
    
    ggsave(file.path(CONFIG$outdir, "Barplot_NonMonotonic_Patterns.pdf"),
           p_pattern, width = 12, height = 8)
    message("[SAVED] Barplot_NonMonotonic_Patterns.pdf")
}


## =========================================================
## SECTION 8: MERGED BARPLOTS (MALE VS FEMALE)
## =========================================================

banner("STEP 7: Creating merged comparison plots")

# Aggregate across all celltypes
merged_summary <- plot_dt[, .N, by = .(Sex, Stability)]
merged_summary[, Sex := factor(Sex, levels = c("male", "female"))]
merged_summary[, Stability := factor(Stability, levels = STABILITY_LEVELS)]

## --- Stacked counts ---
p_merged_counts <- ggplot(
    merged_summary,
    aes(x = Sex, y = N, fill = Stability)
) +
    geom_col(position = "stack", width = 0.7, color = "black", linewidth = 0.3) +
    scale_fill_manual(values = STABILITY_COLS) +
    scale_x_discrete(labels = c("male" = "Male", "female" = "Female")) +
    scale_y_continuous(labels = scales::comma, expand = expansion(mult = c(0, 0.05))) +
    labs(
        title = "Compartment Stability Distribution by Sex",
        subtitle = "All cell types merged",
        x = NULL,
        y = "Number of 80 kb bins",
        fill = "Stability Class"
    ) +
    theme_bw(base_size = 14) +
    theme(
        axis.text.x = element_text(size = 14, face = "bold"),
        legend.position = "right",
        plot.title = element_text(face = "bold", size = 16)
    )

ggsave(file.path(CONFIG$outdir, "Barplot_Merged_BySex_Counts.pdf"),
       p_merged_counts, width = 8, height = 6)
message("[SAVED] Barplot_Merged_BySex_Counts.pdf")

## --- Stacked proportions ---
p_merged_prop <- ggplot(
    merged_summary,
    aes(x = Sex, y = N, fill = Stability)
) +
    geom_col(position = "fill", width = 0.7, color = "black", linewidth = 0.3) +
    scale_fill_manual(values = STABILITY_COLS) +
    scale_x_discrete(labels = c("male" = "Male", "female" = "Female")) +
    scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.02))) +
    labs(
        title = "Compartment Stability Proportions by Sex",
        subtitle = "All cell types merged",
        x = NULL,
        y = "Fraction of 80 kb bins",
        fill = "Stability Class"
    ) +
    theme_bw(base_size = 14) +
    theme(
        axis.text.x = element_text(size = 14, face = "bold"),
        legend.position = "right",
        plot.title = element_text(face = "bold", size = 16)
    )

ggsave(file.path(CONFIG$outdir, "Barplot_Merged_BySex_Proportions.pdf"),
       p_merged_prop, width = 8, height = 6)
message("[SAVED] Barplot_Merged_BySex_Proportions.pdf")

## --- Switching-only merged ---
switching_merged <- switching_only[, .N, by = .(Sex, Stability)]
switching_merged[, Sex := factor(Sex, levels = c("male", "female"))]
switching_merged[, Stability := factor(Stability, levels = SWITCHING_CLASSES)]

p_switching_merged_counts <- ggplot(
    switching_merged,
    aes(x = Sex, y = N, fill = Stability)
) +
    geom_col(position = "stack", width = 0.7, color = "black", linewidth = 0.3) +
    scale_fill_manual(values = STABILITY_COLS) +
    scale_x_discrete(labels = c("male" = "Male", "female" = "Female")) +
    scale_y_continuous(labels = scales::comma, expand = expansion(mult = c(0, 0.05))) +
    labs(
        title = "Switching Bins by Sex",
        subtitle = "Monotonic vs Non-Monotonic (all cell types merged)",
        x = NULL,
        y = "Number of switching bins",
        fill = "Switch Type"
    ) +
    theme_bw(base_size = 14) +
    theme(
        axis.text.x = element_text(size = 14, face = "bold"),
        legend.position = "right",
        plot.title = element_text(face = "bold", size = 16)
    )

ggsave(file.path(CONFIG$outdir, "Barplot_Switching_Merged_BySex_Counts.pdf"),
       p_switching_merged_counts, width = 8, height = 6)
message("[SAVED] Barplot_Switching_Merged_BySex_Counts.pdf")


## =========================================================
## SECTION 9: SUMMARY STATISTICS
## =========================================================

banner("STEP 8: Exporting summary statistics")

# Summary by stability class
summary_5class <- stability_ext_dt[, .N, by = .(Sex, Celltype, Stability)]
summary_5class <- summary_5class[order(Sex, Celltype, -N)]

fwrite(summary_5class, file.path(CONFIG$outdir, "stability_5class_summary.tsv"), sep = "\t")

# Key findings
stable_n <- stability_ext_dt[Stability %in% c("Stable_Active", "Stable_Repressive"), .N]
mono_n <- stability_ext_dt[Stability %in% c("Monotonic_A_to_R", "Monotonic_R_to_A"), .N]
nonmono_n <- stability_ext_dt[Stability == "Non_Monotonic", .N]
total_n <- stability_ext_dt[!Stability %in% c("Missing", "Insufficient_Data"), .N]

message("\n[KEY FINDINGS]")
message(sprintf("  Stable bins:           %d (%.1f%%)", stable_n, 100*stable_n/total_n))
message(sprintf("  Monotonic switching:   %d (%.1f%%)", mono_n, 100*mono_n/total_n))
message(sprintf("  Non-monotonic:         %d (%.1f%%)", nonmono_n, 100*nonmono_n/total_n))


## =========================================================
## FINAL SUMMARY
## =========================================================

cat("\n")
cat("╔══════════════════════════════════════════════════════════════════╗\n")
cat("║           5-CLASS STABILITY ANALYSIS COMPLETE                    ║\n")
cat("╠══════════════════════════════════════════════════════════════════╣\n")
cat("║  Files generated:                                                ║\n")
cat("║    • compartment_stability_5class.tsv                            ║\n")
cat("║    • stability_3class_vs_5class.tsv                              ║\n")
cat("║    • stability_5class_summary.tsv                                ║\n")
cat("║    • nonmonotonic_trajectories_detail.tsv                        ║\n")
cat("║    • nonmonotonic_trajectory_counts.tsv                          ║\n")
cat("║    • Barplot_Stability_5class_Counts.pdf                         ║\n")
cat("║    • Barplot_Stability_5class_Proportions.pdf                    ║\n")
cat("║    • Barplot_Switching_Breakdown.pdf                             ║\n")
cat("║    • Barplot_NonMonotonic_Patterns.pdf                           ║\n")
cat("║    • Barplot_Merged_BySex_Counts.pdf                             ║\n")
cat("║    • Barplot_Merged_BySex_Proportions.pdf                        ║\n")
cat("║    • Barplot_Switching_Merged_BySex_Counts.pdf                   ║\n")
cat("╚══════════════════════════════════════════════════════════════════╝\n")

message("\n[INFO] Next step: Run compartment_03_bin_annotations.R")
