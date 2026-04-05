#!/usr/bin/env Rscript
# ==============================================================================
# Stability class distribution per ECG compartment
# ==============================================================================
#
# Description:
#   Computes and visualizes the distribution of 5-class stability labels
#   (Stable_Active, Stable_Repressive, Monotonic_A_to_R, Monotonic_R_to_A,
#   Non_Monotonic) within each ECG compartment. ECG clusters are reordered
#   by mean activity (ECG 1 = most active) to match the main heatmap.
#
# Input:
#   - compartments_binary_A1_R0.rds (from compartment_01_main_analysis.R)
#   - kmeans_cluster_assignments.tsv
#   - compartment_stability_5class.tsv
#
# Output:
#   - Barplot_Stability_per_ECG_Overall_Stacked.pdf/png
#   - stability_per_ECG_ordered.tsv (long format)
#   - stability_per_ECG_summary.tsv (wide format)
#
# ==============================================================================
# SETUP
# ==============================================================================

suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(scales)
})

# ------------------------------------------------------------------------------
# Configuration - Update this path
# ------------------------------------------------------------------------------

input_dir  <- "path/to/output/compartment_main_heatmap"
output_dir <- file.path(input_dir, "stability_per_ECG")

if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
}

# ------------------------------------------------------------------------------
# Constants
# ------------------------------------------------------------------------------

STABILITY_LEVELS <- c(
    "Stable_Active",
    "Stable_Repressive",
    "Monotonic_A_to_R",
    "Monotonic_R_to_A",
    "Non_Monotonic"
)

STABILITY_COLS <- c(
    "Stable_Active"     = "#FF6B6B",
    "Stable_Repressive" = "#4E79A7",
    "Monotonic_A_to_R"  = "#8B0000",
    "Monotonic_R_to_A"  = "#00008B",
    "Non_Monotonic"     = "#708238"
)

STABILITY_LABELS <- c(
    "Stable_Active"     = "Stable Active",
    "Stable_Repressive" = "Stable Repressive",
    "Monotonic_A_to_R"  = "A to R (Closing)",
    "Monotonic_R_to_A"  = "R to A (Opening)",
    "Non_Monotonic"     = "Non-Monotonic"
)


# ==============================================================================
# LOAD DATA
# ==============================================================================

message("[LOAD] Loading data...")

comp_bin <- readRDS(file.path(input_dir, "compartments_binary_A1_R0.rds"))
message(sprintf("  [OK] Compartments: %d bins x %d samples", nrow(comp_bin), ncol(comp_bin)))

cluster_dt <- fread(file.path(input_dir, "kmeans_cluster_assignments.tsv"))
message(sprintf("  [OK] Original clusters: %d bins", nrow(cluster_dt)))

stability_dt <- fread(file.path(input_dir, "compartment_stability_5class.tsv"))
message(sprintf("  [OK] Stability data: %d rows", nrow(stability_dt)))


# ==============================================================================
# RECOMPUTE ECG ORDERING (MATCHES MAIN PIPELINE)
# ==============================================================================

message("\n[COMPUTE] Reordering ECG clusters by activity...")

original_clusters <- setNames(cluster_dt$ecg_cluster, cluster_dt$bin_id)

# Mean fraction active per bin across all samples
frac_active_all <- rowMeans(comp_bin == 1L, na.rm = TRUE)
names(frac_active_all) <- rownames(comp_bin)

# Mean activity per original cluster
cluster_activity <- tapply(
    frac_active_all,
    original_clusters[names(frac_active_all)],
    mean, na.rm = TRUE
)

# Order: most active = ECG 1
cluster_order  <- order(cluster_activity, decreasing = TRUE)
original_levels <- as.character(sort(unique(original_clusters)))
ordered_levels  <- original_levels[cluster_order]

n_clusters <- length(ordered_levels)
new_labels <- as.character(1:n_clusters)
ecg_remap  <- setNames(new_labels, ordered_levels)

cluster_dt[, ECG_ordered := ecg_remap[as.character(ecg_cluster)]]
cluster_dt[, ECG_ordered := factor(ECG_ordered, levels = new_labels)]

message(sprintf("  [OK] Remapped %d clusters (ordered by activity)", n_clusters))

message("\n  Cluster activity ranking (original -> new label):")
for (i in 1:min(10, n_clusters)) {
    old <- ordered_levels[i]
    new <- new_labels[i]
    act <- cluster_activity[old]
    message(sprintf("    Original %s -> ECG %s (activity: %.2f)", old, new, act))
}
if (n_clusters > 10) message("    ...")


# ==============================================================================
# MERGE WITH STABILITY DATA
# ==============================================================================

message("\n[MERGE] Combining with stability data...")

stability_clean <- stability_dt[!Stability %in% c("Missing", "Insufficient_Data")]

merged_dt <- merge(
    stability_clean,
    cluster_dt[, .(bin_id, ECG_ordered)],
    by = "bin_id",
    all.x = TRUE
)
merged_dt <- merged_dt[!is.na(ECG_ordered)]
merged_dt[, Stability := factor(Stability, levels = STABILITY_LEVELS)]

message(sprintf("  [OK] Merged data: %d rows", nrow(merged_dt)))
message(sprintf("  [OK] ECG range: 1 to %d", n_clusters))


# ==============================================================================
# COMPUTE DISTRIBUTION
# ==============================================================================

message("\n[COMPUTE] Computing stability distribution per ECG...")

overall_dist <- merged_dt[, .(N = .N), by = .(ECG_ordered, Stability)]
overall_dist[, Total := sum(N), by = ECG_ordered]
overall_dist[, Percentage := N / Total * 100]

# Fill missing combinations with 0
all_combos <- CJ(
    ECG_ordered = factor(1:n_clusters, levels = 1:n_clusters),
    Stability   = factor(STABILITY_LEVELS, levels = STABILITY_LEVELS)
)
overall_dist <- merge(all_combos, overall_dist,
                      by = c("ECG_ordered", "Stability"), all.x = TRUE)
overall_dist[is.na(N), N := 0]
overall_dist[is.na(Percentage), Percentage := 0]
overall_dist[, Total := sum(N), by = ECG_ordered]

fwrite(overall_dist, file.path(output_dir, "stability_per_ECG_ordered.tsv"), sep = "\t")


# ==============================================================================
# STACKED BAR CHART
# ==============================================================================

message("\n[PLOT] Creating stacked bar chart...")

p <- ggplot(overall_dist, aes(x = ECG_ordered, y = Percentage, fill = Stability)) +
    geom_col(position = "stack", color = "black", linewidth = 0.2, width = 0.85) +
    scale_fill_manual(
        values = STABILITY_COLS,
        labels = STABILITY_LABELS,
        name = "Stability Class"
    ) +
    scale_x_discrete(drop = FALSE) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.02)),
                       breaks = seq(0, 100, 20)) +
    labs(
        title = "Stability Class Distribution Within Each ECG Compartment",
        subtitle = sprintf(
            "ECG 1 = Most Active | ECG %d = Least Active | All samples combined",
            n_clusters
        ),
        x = "ECG Compartment",
        y = "Percentage (%%)"
    ) +
    theme_bw(base_size = 14) +
    theme(
        text = element_text(family = "Arial"),
        axis.text.x  = element_text(face = "bold", size = 11),
        axis.text.y  = element_text(face = "bold", size = 11),
        axis.title   = element_text(face = "bold", size = 13),
        plot.title   = element_text(face = "bold", size = 16),
        plot.subtitle = element_text(size = 11, color = "gray30"),
        legend.position = "right",
        legend.title = element_text(face = "bold"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor   = element_blank()
    )

ggsave(file.path(output_dir, "Barplot_Stability_per_ECG_Overall_Stacked.pdf"),
       p, width = 16, height = 8, useDingbats = FALSE)
ggsave(file.path(output_dir, "Barplot_Stability_per_ECG_Overall_Stacked.png"),
       p, width = 16, height = 8, dpi = 300)

message("  [OK] Saved: Barplot_Stability_per_ECG_Overall_Stacked.pdf/png")


# ==============================================================================
# SUMMARY TABLE
# ==============================================================================

message("\n[SUMMARY] Stability distribution per ECG compartment:")

summary_wide <- dcast(overall_dist, ECG_ordered ~ Stability, value.var = "Percentage")
summary_wide[, Total_Bins := overall_dist[ECG_ordered == .BY[[1]], sum(N)], by = ECG_ordered]

cat("\nECG\tTotal\tStable_A\tStable_R\tA->R\tR->A\tNon-Mono\n")
cat("---\t-----\t--------\t--------\t----\t----\t--------\n")
for (i in 1:nrow(summary_wide)) {
    row <- summary_wide[i]
    cat(sprintf("%s\t%d\t%.1f%%\t\t%.1f%%\t\t%.1f%%\t%.1f%%\t%.1f%%\n",
                as.character(row$ECG_ordered),
                row$Total_Bins,
                row$Stable_Active,
                row$Stable_Repressive,
                row$Monotonic_A_to_R,
                row$Monotonic_R_to_A,
                row$Non_Monotonic))
}

fwrite(summary_wide, file.path(output_dir, "stability_per_ECG_summary.tsv"), sep = "\t")


# ==============================================================================
# SESSION INFO
# ==============================================================================

writeLines(
    capture.output(sessionInfo()),
    file.path(output_dir, "sessionInfo_stability_ECG.txt")
)


# ==============================================================================
# SUMMARY
# ==============================================================================

cat("\n")
cat("========================================================================\n")
cat("  COMPLETE\n")
cat("========================================================================\n")
cat(sprintf("  Output: %s\n", output_dir))
cat("\n")
cat("  Files:\n")
cat("    Barplot_Stability_per_ECG_Overall_Stacked.pdf/png\n")
cat("    stability_per_ECG_ordered.tsv\n")
cat("    stability_per_ECG_summary.tsv\n")
cat("    sessionInfo_stability_ECG.txt\n")
cat("========================================================================\n")
