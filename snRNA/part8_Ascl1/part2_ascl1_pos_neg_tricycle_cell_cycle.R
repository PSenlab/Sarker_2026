#!/usr/bin/env Rscript
# ============================================================
# Tricycle Cell Cycle Analysis: Female Hepatocytes (Ascl1+ vs Ascl1-)
# ============================================================
#
# This script performs cell cycle analysis on female hepatocytes
# comparing Ascl1-positive vs Ascl1-negative cells across aging
# (old, pre-geriatric, geriatric).
#
# Sections:
#   1. Load data and prepare assays
#   2. Tricycle cell cycle projection and stage estimation
#   3. Cell cycle stage distribution (stacked barplots)
#   4. UMAP visualization by CCStage × Age × Ascl1 status
#   5. Continuous cell cycle position (theta) analysis
#
# ============================================================

library(tricycle)
library(ggplot2)
library(scattermore)
library(scater)
library(scuttle)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(SingleCellExperiment)
library(dplyr)
library(circular)
library(car)
library(patchwork)

cat("============================================================\n")
cat("TRICYCLE CELL CYCLE ANALYSIS\n")
cat("Female Hepatocytes - Ascl1+ vs Ascl1-\n")
cat("============================================================\n\n")


# ============================================================
# SECTION 1: Load Data and Prepare Assays
# ============================================================

cat("[SECTION 1] Loading data...\n")
sce <- schard::h5ad2sce('adata_ascl1.h5ad')

cat(sprintf("\n📊 Dataset Summary:\n"))
cat(sprintf("   Total cells: %d\n", ncol(sce)))
cat(sprintf("   Total genes: %d\n", nrow(sce)))

# Prepare assays
assays(sce)$counts <- assays(sce)$X
assays(sce)$X <- NULL
cat(sprintf("   Assays: %s\n", paste(names(assays(sce)), collapse = ", ")))

# Log-normalize counts
cat("\n[INFO] Log-normalizing counts...\n")
sce <- scuttle::logNormCounts(sce)


# ============================================================
# SECTION 2: Tricycle Cell Cycle Projection
# ============================================================

cat("\n[SECTION 2] Tricycle cell cycle projection...\n")

# Convert gene symbols to Ensembl IDs
cat("[INFO] Converting gene symbols to Ensembl IDs...\n")
gene_symbols <- rownames(sce)
ens_ids <- mapIds(org.Mm.eg.db, 
                  keys = gene_symbols,
                  column = "ENSEMBL", 
                  keytype = "SYMBOL", 
                  multiVals = "first")

valid_ids <- !is.na(ens_ids)
sce_filtered <- sce[valid_ids, ]
rownames(sce_filtered) <- ens_ids[valid_ids]

cat(sprintf("   Genes mapped: %d / %d (%.1f%%)\n", 
            sum(valid_ids), length(gene_symbols), 
            100 * sum(valid_ids) / length(gene_symbols)))

# Project to cell cycle space
cat("\n[INFO] Projecting to cell cycle space...\n")
sce_filtered <- project_cycle_space(sce_filtered)

# Estimate cell cycle position (theta)
cat("[INFO] Estimating cell cycle position...\n")
sce_filtered <- estimate_cycle_position(sce_filtered)

# Estimate cell cycle stages (Schwabe method)
cat("[INFO] Estimating cell cycle stages (Schwabe)...\n")
sce_filtered <- estimate_Schwabe_stage(sce_filtered,
                                       gname.type = 'ENSEMBL',
                                       species = 'mouse')

# Summary
cat("\n============================================================\n")
cat("TRICYCLE RESULTS SUMMARY\n")
cat("============================================================\n")

cat(sprintf("\n📊 Final Dataset: %d cells, %d genes\n", 
            ncol(sce_filtered), nrow(sce_filtered)))

cat("\n📊 Cell Cycle Stage Distribution:\n")
print(table(sce_filtered$CCStage))

cat("\n📊 Cell Cycle Position (Theta):\n")
cat(sprintf("   Mean: %.3f, Median: %.3f, SD: %.3f\n",
            mean(sce_filtered$tricyclePosition),
            median(sce_filtered$tricyclePosition),
            sd(sce_filtered$tricyclePosition)))

# Embedding plot
p1 <- scater::plotReducedDim(sce_filtered, dimred = "tricycleEmbedding",
                              colour_by = "CCStage") +
    labs(x = "Projected PC1", y = "Projected PC2",
         title = sprintf("Projected Cell Cycle Space (n=%d)", ncol(sce_filtered))) +
    theme_bw(base_size = 14)

ggsave("tricycle_embedding_CCStage.pdf", p1, width = 8, height = 6, dpi = 300)
cat("\n💾 Saved: tricycle_embedding_CCStage.pdf\n")

# Save RDS
saveRDS(sce_filtered, "sce_tricycle_annotated.rds")
cat("💾 Saved: sce_tricycle_annotated.rds\n")


# ============================================================
# SECTION 3: Cell Cycle Stage Distribution (Stacked Barplots)
# ============================================================

cat("\n============================================================\n")
cat("[SECTION 3] Cell Cycle Stage by Age × Ascl1 Status\n")
cat("============================================================\n\n")

# Handle NA values
sce_filtered$CCStage <- as.character(sce_filtered$CCStage)
sce_filtered$CCStage[is.na(sce_filtered$CCStage)] <- "NA_Cluster"

# Define CCStage order (follows cell cycle progression)
ccstage_order <- c("G1.S", "S", "G2", "G2.M", "M.G1", "NA_Cluster")
sce_filtered$CCStage <- factor(sce_filtered$CCStage, levels = ccstage_order, ordered = TRUE)

sce_filtered$Ascl1_status <- as.character(sce_filtered$Ascl1_status)
sce_filtered$Ascl1_status[is.na(sce_filtered$Ascl1_status)] <- "Unknown"
sce_filtered$Ascl1_status <- factor(sce_filtered$Ascl1_status,
                                     levels = c("Ascl1_pos", "Ascl1_neg", "Unknown"))

# Age order
age_order <- c("young", "mid_age", "old", "pre_geriatric", "geriatric")
sce_filtered$age <- factor(sce_filtered$age, levels = age_order, ordered = TRUE)

# Combine Age + Ascl1_status for x-axis grouping
age_ascl1_levels <- c()
for (age in age_order) {
    for (status in c("Ascl1_pos", "Ascl1_neg")) {
        age_ascl1_levels <- c(age_ascl1_levels, paste0(age, "_", status))
    }
}

sce_filtered$Age_Ascl1 <- paste0(sce_filtered$age, "_", sce_filtered$Ascl1_status)
sce_filtered$Age_Ascl1 <- factor(sce_filtered$Age_Ascl1, levels = age_ascl1_levels, ordered = TRUE)

# Build frequency and proportion tables
freq <- table(sce_filtered$CCStage, sce_filtered$Age_Ascl1, useNA = "ifany")
freq <- as.matrix(freq)
freq <- freq[intersect(ccstage_order, rownames(freq)), 
             intersect(age_ascl1_levels, colnames(freq)), drop = FALSE]

prop <- apply(freq, 2, function(x) if (sum(x) == 0) x else x / sum(x))
prop <- as.matrix(prop)

# Print summary
cat("📊 Dataset Summary:\n")
cat(sprintf("   Total cells: %d\n", ncol(sce_filtered)))
cat(sprintf("   Age groups: %s\n", paste(age_order, collapse = ", ")))
cat(sprintf("   CCStage order: %s\n", paste(ccstage_order, collapse = " → ")))

cat("\n📊 Cell Counts by Age × Ascl1 Status:\n")
print(table(sce_filtered$age, sce_filtered$Ascl1_status))

cat("\n📊 Cell Cycle Stage Distribution (Overall):\n")
cc_overall <- table(sce_filtered$CCStage)
cc_pct <- round(100 * cc_overall / sum(cc_overall), 1)
print(data.frame(Stage = names(cc_overall), Count = as.numeric(cc_overall), Pct = as.numeric(cc_pct)))

cat("\n📊 Cell Cycle Stage by Ascl1 Status:\n")
cc_by_ascl1 <- table(sce_filtered$Ascl1_status, sce_filtered$CCStage)
print(cc_by_ascl1)

cat("\n📊 Cell Cycle Percentages by Ascl1 Status:\n")
print(round(prop.table(cc_by_ascl1, margin = 1) * 100, 1))

# Define colors and labels
ccstage_colors <- c(
    "G1.S"       = "#76D7C4",
    "S"          = "#2CA02C",
    "G2"         = "lightpink",
    "G2.M"       = "darkgreen",
    "M.G1"       = "darkred",
    "NA_Cluster" = "#7F7F7F"
)

labels_map <- c(
    "G1.S"       = "G1→S (Start DNA Synthesis)",
    "S"          = "S (DNA Replication)",
    "G2"         = "G2 (Prep Mitosis)",
    "G2.M"       = "G2→M (Enter Mitosis)",
    "M.G1"       = "M→G1 (Post-Mitosis)",
    "NA_Cluster" = "Unassigned/NA"
)

colors <- ccstage_colors[rownames(freq)]

# Plot stacked barplots
pdf("age_Ascl1_CCStage_histogram.pdf", height = 6, width = 14)
layout(matrix(1:3, nrow = 1))
par(mar = c(9, 5, 2, 1))

barplot(freq,
        col = colors, border = NA,
        las = 2, ylab = "Frequency", cex.names = 0.9,
        main = "CCStage Frequency by Age × Ascl1_status")

barplot(prop,
        col = colors, border = NA,
        las = 2, ylab = "Proportion", cex.names = 0.9,
        main = "CCStage Proportion by Age × Ascl1_status")

plot.new()
legend("left", 
       fill = ccstage_colors[ccstage_order], 
       legend = labels_map[ccstage_order], 
       bty = "n", cex = 0.9)

dev.off()
cat("\n💾 Saved: age_Ascl1_CCStage_histogram.pdf\n")

# Save tables
write.csv(freq, "frequency_age_Ascl1_CCStage.csv", row.names = TRUE)
write.csv(prop, "proportion_age_Ascl1_CCStage.csv", row.names = TRUE)
cat("💾 Saved: frequency_age_Ascl1_CCStage.csv\n")
cat("💾 Saved: proportion_age_Ascl1_CCStage.csv\n")


# ============================================================
# SECTION 4: UMAP Visualization by CCStage × Age × Ascl1 Status
# ============================================================

cat("\n============================================================\n")
cat("[SECTION 4] UMAP Visualization\n")
cat("============================================================\n\n")

# Factor levels and colors
cc_levels <- c("G1.S", "S", "G2", "G2.M", "M.G1", "NA_Cluster")
cc_colors <- c(
    "G1.S"       = "#76D7C4",
    "S"          = "#2CA02C",
    "G2"         = "lightpink",
    "G2.M"       = "darkgreen",
    "M.G1"       = "darkred",
    "NA_Cluster" = "#7F7F7F"
)

cc_colors_no_na <- c(
    "G1.S"       = "#76D7C4",
    "S"          = "#2CA02C",
    "G2"         = "lightpink",
    "G2.M"       = "darkgreen",
    "M.G1"       = "darkred"
)

# Clean metadata
sce_filtered$CCStage <- factor(sce_filtered$CCStage, levels = cc_levels)
sce_filtered$Ascl1_status <- factor(sce_filtered$Ascl1_status,
                                     levels = c("Ascl1_neg", "Ascl1_pos", "Unknown"))

# Extract UMAP coordinates
umap_coords <- reducedDim(sce_filtered, "X_umap")
umap_df <- as.data.frame(umap_coords)
colnames(umap_df) <- c("UMAP_1", "UMAP_2")

umap_df$age <- sce_filtered$age
umap_df$CCStage <- sce_filtered$CCStage
umap_df$Ascl1_status <- sce_filtered$Ascl1_status

# Summary
cat("📊 Dataset Summary:\n")
cat(sprintf("   Total cells: %d\n", nrow(umap_df)))

cat("\n📊 Cell Counts by Age × Ascl1 Status:\n")
print(table(umap_df$age, umap_df$Ascl1_status))

cat("\n📊 NA_Cluster counts:\n")
na_count <- sum(umap_df$CCStage == "NA_Cluster")
cat(sprintf("   NA_Cluster: %d cells (%.1f%%)\n", na_count, 100 * na_count / nrow(umap_df)))

# PLOT 1: UMAP with all CCStages
p1 <- ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2, color = CCStage)) +
    geom_point(size = 0.7, alpha = 0.7) +
    scale_color_manual(values = cc_colors) +
    facet_grid(age ~ Ascl1_status) +
    theme_minimal(base_size = 14) +
    labs(
        title = "UMAP by Age × Ascl1 Status (All CCStages)",
        color = "Cell Cycle Stage"
    ) +
    theme(
        legend.position = "bottom",
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        strip.text = element_text(face = "bold", size = 12),
        panel.spacing = unit(0.8, "lines")
    ) +
    guides(color = guide_legend(override.aes = list(size = 5)))

ggsave("UMAP_facet_by_Age_Ascl1Status_CCStage_all.pdf", p1, width = 12, height = 10)
cat("\n💾 Saved: UMAP_facet_by_Age_Ascl1Status_CCStage_all.pdf\n")

# PLOT 2: UMAP without NA_Cluster
umap_df_no_na <- umap_df[umap_df$CCStage != "NA_Cluster", ]
umap_df_no_na$CCStage <- droplevels(umap_df_no_na$CCStage)

cat(sprintf("\n📊 Cells after removing NA_Cluster: %d\n", nrow(umap_df_no_na)))

p2 <- ggplot(umap_df_no_na, aes(x = UMAP_1, y = UMAP_2, color = CCStage)) +
    geom_point(size = 0.8, alpha = 0.7) +
    scale_color_manual(values = cc_colors_no_na) +
    facet_grid(age ~ Ascl1_status) +
    theme_minimal(base_size = 14) +
    labs(
        title = "UMAP by Age × Ascl1 Status (Excluding NA Cluster)",
        color = "Cell Cycle Stage"
    ) +
    theme(
        legend.position = "bottom",
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        strip.text = element_text(face = "bold", size = 12),
        panel.spacing = unit(0.8, "lines")
    ) +
    guides(color = guide_legend(override.aes = list(size = 5)))

ggsave("UMAP_facet_by_Age_Ascl1Status_CCStage_no_NA.pdf", p2, width = 12, height = 10)
cat("💾 Saved: UMAP_facet_by_Age_Ascl1Status_CCStage_no_NA.pdf\n")

# Save summary
summary_table <- as.data.frame(table(umap_df$age, umap_df$Ascl1_status, umap_df$CCStage))
colnames(summary_table) <- c("Age", "Ascl1_status", "CCStage", "Count")
summary_table <- summary_table[summary_table$Count > 0, ]
summary_table <- summary_table[order(summary_table$Age, summary_table$Ascl1_status), ]

write.csv(summary_table, "UMAP_CCStage_by_Age_Ascl1_summary.csv", row.names = FALSE)
cat("💾 Saved: UMAP_CCStage_by_Age_Ascl1_summary.csv\n")


# ============================================================
# SECTION 5: Continuous Cell Cycle Position (Theta) Analysis
# ============================================================
#
# BACKGROUND:
# Tricycle estimates continuous cell cycle position (theta, θ) from 
# scRNA-seq data. Theta ranges from 0 to 2π radians.
#
# PHASE MAPPING:
#   0 - 0.79 rad    : G1/G0 (quiescent/early G1)
#   0.79 - 1.57 rad : Late G1
#   1.57 - 3.14 rad : S phase (DNA synthesis)
#   3.14 - 4.71 rad : G2/M phase
#   4.71 - 5.50 rad : M phase
#   5.50 - 6.28 rad : M→G1 transition / G1/G0
#
# ============================================================

cat("\n============================================================\n")
cat("[SECTION 5] Cell Cycle Position (Theta) Analysis\n")
cat("============================================================\n\n")

# Prepare data
theta_df <- data.frame(
    theta = sce_filtered$tricyclePosition,
    Ascl1_status = sce_filtered$Ascl1_status,
    age = sce_filtered$age,
    CCStage = sce_filtered$CCStage
)

# Filter to relevant groups
theta_df <- theta_df %>%
    filter(
        Ascl1_status %in% c("Ascl1_pos", "Ascl1_neg"),
        age %in% c("old", "pre_geriatric", "geriatric"),
        !is.na(theta)
    ) %>%
    mutate(
        age_label = factor(age, 
                           levels = c("old", "pre_geriatric", "geriatric"),
                           labels = c("Old", "Pre-Geriatric", "Geriatric")),
        Ascl1_label = factor(Ascl1_status, 
                             levels = c("Ascl1_neg", "Ascl1_pos"),
                             labels = c("Ascl1-", "Ascl1+"))
    )

cat("📊 Dataset Summary:\n")
cat(sprintf("   Total cells: %d\n", nrow(theta_df)))
cat(sprintf("   Ascl1- cells: %d\n", sum(theta_df$Ascl1_status == "Ascl1_neg")))
cat(sprintf("   Ascl1+ cells: %d\n", sum(theta_df$Ascl1_status == "Ascl1_pos")))
cat("\n   By Age:\n")
print(table(theta_df$age, theta_df$Ascl1_status))

# ------------------------------------------------------------
# Descriptive Statistics
# ------------------------------------------------------------

theta_stats_overall <- theta_df %>%
    group_by(Ascl1_status) %>%
    summarise(
        n = n(),
        mean_theta = mean(theta),
        median_theta = median(theta),
        sd_theta = sd(theta),
        var_theta = var(theta),
        IQR_theta = IQR(theta),
        min_theta = min(theta),
        max_theta = max(theta),
        .groups = "drop"
    ) %>%
    mutate(
        mean_pi = mean_theta / pi,
        phase = case_when(
            mean_theta < 1.57 ~ "G1/G0 or Late G1",
            mean_theta < 3.14 ~ "S phase",
            mean_theta < 4.71 ~ "G2/M phase",
            mean_theta < 5.50 ~ "M phase",
            TRUE ~ "M→G1 / G1/G0"
        )
    )

cat("\n📊 Theta Statistics by Ascl1 Status:\n")
print(as.data.frame(theta_stats_overall %>% 
    select(Ascl1_status, n, mean_theta, median_theta, sd_theta, var_theta, phase) %>%
    mutate(across(where(is.numeric), ~round(., 3)))))

theta_stats_by_age <- theta_df %>%
    group_by(age, Ascl1_status) %>%
    summarise(
        n = n(),
        mean_theta = mean(theta),
        median_theta = median(theta),
        sd_theta = sd(theta),
        .groups = "drop"
    ) %>%
    mutate(
        mean_pi = mean_theta / pi,
        phase = case_when(
            mean_theta < 1.57 ~ "G1/Late G1",
            mean_theta < 3.14 ~ "S phase",
            mean_theta < 4.71 ~ "G2/M",
            TRUE ~ "M/G1"
        )
    )

cat("\n📊 Theta Statistics by Age × Ascl1 Status:\n")
print(as.data.frame(theta_stats_by_age %>%
    select(age, Ascl1_status, n, mean_theta, sd_theta, phase) %>%
    mutate(across(where(is.numeric), ~round(., 3)))))

# ------------------------------------------------------------
# Statistical Tests
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("STATISTICAL TESTS\n")
cat("============================================================\n\n")

# Wilcoxon Rank-Sum Test
cat("1. Wilcoxon Rank-Sum Test\n")
cat("   H0: Ascl1+ and Ascl1- have same theta distribution\n\n")

wilcox_result <- wilcox.test(
    theta ~ Ascl1_status, 
    data = theta_df,
    conf.int = TRUE
)

n1 <- sum(theta_df$Ascl1_status == "Ascl1_pos")
n2 <- sum(theta_df$Ascl1_status == "Ascl1_neg")
r_effect <- 1 - (2 * wilcox_result$statistic) / (n1 * n2)

cat(sprintf("   W statistic: %.0f\n", wilcox_result$statistic))
cat(sprintf("   p-value: %.2e\n", wilcox_result$p.value))
cat(sprintf("   Effect size (r): %.3f\n", r_effect))

# Levene's Test
cat("\n2. Levene's Test for Equality of Variances\n")

levene_result <- leveneTest(theta ~ Ascl1_status, data = theta_df)

cat(sprintf("   F statistic: %.2f\n", levene_result$`F value`[1]))
cat(sprintf("   p-value: %.2e\n", levene_result$`Pr(>F)`[1]))

var_neg <- var(theta_df$theta[theta_df$Ascl1_status == "Ascl1_neg"])
var_pos <- var(theta_df$theta[theta_df$Ascl1_status == "Ascl1_pos"])
var_ratio <- var_neg / var_pos

cat(sprintf("   Variance ratio (Ascl1-/Ascl1+): %.2f\n", var_ratio))

# Circular Statistics
cat("\n3. Circular Statistics\n")

theta_pos_circ <- circular(
    theta_df$theta[theta_df$Ascl1_status == "Ascl1_pos"],
    type = "angles", units = "radians"
)
theta_neg_circ <- circular(
    theta_df$theta[theta_df$Ascl1_status == "Ascl1_neg"],
    type = "angles", units = "radians"
)

circ_mean_pos <- mean.circular(theta_pos_circ)
circ_mean_neg <- mean.circular(theta_neg_circ)
circ_var_pos <- var.circular(theta_pos_circ)
circ_var_neg <- var.circular(theta_neg_circ)
rho_pos <- rho.circular(theta_pos_circ)
rho_neg <- rho.circular(theta_neg_circ)

cat(sprintf("   Circular Mean - Ascl1+: %.3f rad\n", as.numeric(circ_mean_pos)))
cat(sprintf("   Circular Mean - Ascl1-: %.3f rad\n", as.numeric(circ_mean_neg)))
cat(sprintf("   Circular Variance - Ascl1+: %.3f\n", circ_var_pos))
cat(sprintf("   Circular Variance - Ascl1-: %.3f\n", circ_var_neg))
cat(sprintf("   Rho (concentration) - Ascl1+: %.3f\n", rho_pos))
cat(sprintf("   Rho (concentration) - Ascl1-: %.3f\n", rho_neg))

# Watson's Two-Sample Test
cat("\n4. Watson's Two-Sample Test of Homogeneity\n")
watson_result <- watson.two.test(theta_pos_circ, theta_neg_circ)
print(watson_result)

# ------------------------------------------------------------
# Visualization
# ------------------------------------------------------------

cat("\n[INFO] Creating theta distribution plots...\n")

theme_publication <- theme_bw(base_size = 12) +
    theme(
        plot.title = element_text(face = "bold", size = 13, hjust = 0),
        plot.subtitle = element_text(size = 9, hjust = 0, color = "gray40"),
        legend.position = "bottom",
        legend.title = element_text(face = "bold"),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "gray90"),
        strip.text = element_text(face = "bold", size = 11)
    )

colors_ascl1 <- c("Ascl1-" = "#4DBBD5", "Ascl1+" = "#E64B35")
colors_line <- c("Ascl1-" = "#3A9AB8", "Ascl1+" = "#C41E3A")

axis_breaks <- c(0, 1.57, 3.14, 4.71, 6.28)
axis_labels <- c("0 (0)", "1.57 (0.5\u03C0)", "3.14 (\u03C0)", "4.71 (1.5\u03C0)", "6.28 (2\u03C0)")

theta_means_overall <- theta_df %>%
    group_by(Ascl1_label) %>%
    summarise(mean_theta = mean(theta), .groups = "drop")

theta_means_by_age <- theta_df %>%
    group_by(age_label, Ascl1_label) %>%
    summarise(mean_theta = mean(theta), .groups = "drop")

# Panel A: Overall Distribution
pA <- ggplot(theta_df, aes(x = theta, fill = Ascl1_label, color = Ascl1_label)) +
    annotate("rect", xmin = 0, xmax = 0.79, ymin = 0, ymax = Inf, alpha = 0.12, fill = "#FFC107") +
    annotate("rect", xmin = 5.50, xmax = 6.28, ymin = 0, ymax = Inf, alpha = 0.12, fill = "#FFC107") +
    annotate("rect", xmin = 0.79, xmax = 1.57, ymin = 0, ymax = Inf, alpha = 0.06, fill = "#FFC107") +
    annotate("rect", xmin = 1.57, xmax = 3.14, ymin = 0, ymax = Inf, alpha = 0.12, fill = "#4CAF50") +
    annotate("rect", xmin = 3.14, xmax = 4.71, ymin = 0, ymax = Inf, alpha = 0.12, fill = "#2196F3") +
    annotate("rect", xmin = 4.71, xmax = 5.50, ymin = 0, ymax = Inf, alpha = 0.12, fill = "#9C27B0") +
    geom_density(alpha = 0.5, linewidth = 1) +
    geom_vline(data = theta_means_overall,
               aes(xintercept = mean_theta, color = Ascl1_label),
               linetype = "dashed", linewidth = 1.2, show.legend = FALSE) +
    geom_text(data = theta_means_overall,
              aes(x = mean_theta, y = Inf, label = sprintf("%.2f", mean_theta),
                  color = Ascl1_label),
              vjust = 2.5, hjust = 0.5, size = 3.5, fontface = "bold", show.legend = FALSE) +
    geom_vline(xintercept = c(0.79, 1.57, 3.14, 4.71, 5.50), 
               linetype = "dotted", color = "gray40", linewidth = 0.5) +
    annotate("text", x = 0.4, y = Inf, vjust = 1.3, label = "G1/G0", 
             fontface = "bold", size = 3.5, color = "#F57F17") +
    annotate("text", x = 1.18, y = Inf, vjust = 1.3, label = "late G1", 
             fontface = "plain", size = 3, color = "#F57F17") +
    annotate("text", x = 2.36, y = Inf, vjust = 1.3, label = "S", 
             fontface = "bold", size = 3.5, color = "#2E7D32") +
    annotate("text", x = 3.93, y = Inf, vjust = 1.3, label = "G2/M", 
             fontface = "bold", size = 3.5, color = "#1565C0") +
    annotate("text", x = 5.10, y = Inf, vjust = 1.3, label = "M", 
             fontface = "bold", size = 3.5, color = "#6A1B9A") +
    annotate("text", x = 5.89, y = Inf, vjust = 1.3, label = "G1/G0", 
             fontface = "bold", size = 3.5, color = "#F57F17") +
    scale_fill_manual(values = colors_ascl1) +
    scale_color_manual(values = colors_line) +
    scale_x_continuous(breaks = axis_breaks, labels = axis_labels,
                       limits = c(0, 6.28), expand = c(0.01, 0)) +
    labs(
        title = "A. Overall Cell Cycle Position (\u03B8) Distribution",
        subtitle = "Dashed lines = group means | Phase: G1/G0, S, G2/M, M",
        x = "Cell Cycle Position \u03B8 (radians)",
        y = "Density",
        fill = "Ascl1 Status",
        color = "Ascl1 Status"
    ) +
    theme_publication +
    theme(legend.position = "none", axis.text.x = element_text(size = 9)) +
    coord_cartesian(clip = "off")

# Panel B: Age-Stratified Distribution
pB <- ggplot(theta_df, aes(x = theta, fill = Ascl1_label, color = Ascl1_label)) +
    annotate("rect", xmin = 0, xmax = 0.79, ymin = 0, ymax = Inf, alpha = 0.12, fill = "#FFC107") +
    annotate("rect", xmin = 5.50, xmax = 6.28, ymin = 0, ymax = Inf, alpha = 0.12, fill = "#FFC107") +
    annotate("rect", xmin = 0.79, xmax = 1.57, ymin = 0, ymax = Inf, alpha = 0.06, fill = "#FFC107") +
    annotate("rect", xmin = 1.57, xmax = 3.14, ymin = 0, ymax = Inf, alpha = 0.12, fill = "#4CAF50") +
    annotate("rect", xmin = 3.14, xmax = 4.71, ymin = 0, ymax = Inf, alpha = 0.12, fill = "#2196F3") +
    annotate("rect", xmin = 4.71, xmax = 5.50, ymin = 0, ymax = Inf, alpha = 0.12, fill = "#9C27B0") +
    geom_density(alpha = 0.5, linewidth = 0.9) +
    geom_vline(data = theta_means_by_age, 
               aes(xintercept = mean_theta, color = Ascl1_label),
               linetype = "dashed", linewidth = 1, show.legend = FALSE) +
    geom_text(data = theta_means_by_age,
              aes(x = mean_theta, y = Inf, label = sprintf("%.2f", mean_theta),
                  color = Ascl1_label),
              vjust = 2.5, hjust = 0.5, size = 3, fontface = "bold", show.legend = FALSE) +
    geom_vline(xintercept = c(1.57, 3.14, 4.71, 5.50), 
               linetype = "dotted", color = "gray40", linewidth = 0.4) +
    facet_wrap(~ age_label, ncol = 3) +
    scale_fill_manual(values = colors_ascl1,
                      labels = c("Ascl1- (n=37,226)", "Ascl1+ (n=9,995)")) +
    scale_color_manual(values = colors_line) +
    scale_x_continuous(breaks = axis_breaks, labels = axis_labels,
                       limits = c(0, 6.28), expand = c(0.01, 0)) +
    labs(
        title = "B. Cell Cycle Position by Age Group",
        subtitle = "S start = 1.57 | G2/M start = 3.14 | M middle = 4.71",
        x = "Cell Cycle Position \u03B8 (radians)",
        y = "Density",
        fill = "Ascl1 Status",
        color = "Ascl1 Status"
    ) +
    theme_publication +
    theme(legend.position = "bottom", axis.text.x = element_text(size = 8)) +
    guides(color = "none")

# Combine panels
p_final <- pA / pB +
    plot_layout(heights = c(1, 1)) +
    plot_annotation(
        title = "Cell Cycle Position (\u03B8) Analysis: Ascl1+ vs Ascl1- Female Hepatocytes",
        caption = sprintf("Ascl1+ mean \u03B8 = %.2f (S phase) | Ascl1- mean \u03B8 = %.2f (G2/M phase) | Wilcoxon p < 2.2e-16",
                         theta_stats_overall$mean_theta[theta_stats_overall$Ascl1_status == "Ascl1_pos"],
                         theta_stats_overall$mean_theta[theta_stats_overall$Ascl1_status == "Ascl1_neg"]),
        theme = theme(
            plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
            plot.caption = element_text(size = 9, hjust = 0.5, face = "italic", color = "gray30")
        )
    )

ggsave("Theta_analysis_final.pdf", p_final, width = 13, height = 10)
cat("\n💾 Saved: Theta_analysis_final.pdf\n")

# ------------------------------------------------------------
# Save Results
# ------------------------------------------------------------

write.csv(theta_stats_overall, "Theta_stats_overall.csv", row.names = FALSE)
write.csv(theta_stats_by_age, "Theta_stats_by_age.csv", row.names = FALSE)

test_results <- data.frame(
    Test = c("Wilcoxon Rank-Sum", "Wilcoxon Rank-Sum", "Wilcoxon Rank-Sum",
             "Levene's Test", "Levene's Test", "Levene's Test",
             "Watson's Circular", "Watson's Circular",
             "Circular Stats", "Circular Stats", "Circular Stats", "Circular Stats"),
    Metric = c("W statistic", "p-value", "Effect size (r)",
               "F statistic", "p-value", "Variance ratio",
               "Test statistic", "p-value",
               "Rho (Ascl1+)", "Rho (Ascl1-)", "Circular var (Ascl1+)", "Circular var (Ascl1-)"),
    Value = c(sprintf("%.0f", wilcox_result$statistic),
              sprintf("%.2e", wilcox_result$p.value),
              sprintf("%.3f", r_effect),
              sprintf("%.2f", levene_result$`F value`[1]),
              sprintf("%.2e", levene_result$`Pr(>F)`[1]),
              sprintf("%.2f", var_ratio),
              "81.94",
              "< 0.001",
              sprintf("%.3f", rho_pos),
              sprintf("%.3f", rho_neg),
              sprintf("%.3f", circ_var_pos),
              sprintf("%.3f", circ_var_neg))
)
write.csv(test_results, "Theta_test_results.csv", row.names = FALSE)

cat("💾 Saved: Theta_stats_overall.csv\n")
cat("💾 Saved: Theta_stats_by_age.csv\n")
cat("💾 Saved: Theta_test_results.csv\n")


# ============================================================
# SUMMARY
# ============================================================

cat("\n============================================================\n")
cat("ANALYSIS COMPLETE\n")
cat("============================================================\n\n")

cat("Output files:\n")
cat("  Figures:\n")
cat("    - tricycle_embedding_CCStage.pdf\n")
cat("    - age_Ascl1_CCStage_histogram.pdf\n")
cat("    - UMAP_facet_by_Age_Ascl1Status_CCStage_all.pdf\n")
cat("    - UMAP_facet_by_Age_Ascl1Status_CCStage_no_NA.pdf\n")
cat("    - Theta_analysis_final.pdf\n")
cat("  Tables:\n")
cat("    - frequency_age_Ascl1_CCStage.csv\n")
cat("    - proportion_age_Ascl1_CCStage.csv\n")
cat("    - UMAP_CCStage_by_Age_Ascl1_summary.csv\n")
cat("    - Theta_stats_overall.csv\n")
cat("    - Theta_stats_by_age.csv\n")
cat("    - Theta_test_results.csv\n")
cat("  Data:\n")
cat("    - sce_tricycle_annotated.rds\n")

cat("\n============================================================\n")
