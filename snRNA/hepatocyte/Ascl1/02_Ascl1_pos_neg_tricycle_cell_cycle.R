#!/usr/bin/env Rscript
# ==============================================================================
# Tricycle Cell Cycle Analysis: Ascl1+ vs Ascl1- Female Hepatocytes
# ==============================================================================
#
# Description:
#   Uses the Tricycle package to project female hepatocytes onto the cell
#   cycle manifold, estimate discrete Schwabe stages (G1.S, S, G2, G2.M,
#   M.G1) and continuous cell-cycle position (theta), then compares the
#   distributions between Ascl1+ and Ascl1- cells across older age groups
#   (old, pre_geriatric, geriatric).
#
# Pipeline:
#   STEP 0: Load h5ad, subset to female hepatocytes, annotate Ascl1 status
#   STEP 1: Tricycle projection, theta estimation, Schwabe stage assignment
#   STEP 2: Stacked barplots of stage distribution by age x Ascl1 status
#   STEP 3: UMAP faceted by age x Ascl1 status colored by CCStage
#   STEP 4: Continuous theta analysis + circular statistics (Wilcoxon,
#           Levene, Watson's two-sample, circular mean/var/rho)
#
# Input:
#   - integrated_scvi.h5ad  (canonical annotated AnnData)
#
# Output:
#   Figures:
#     - tricycle_embedding_CCStage.pdf
#     - age_Ascl1_CCStage_histogram.pdf
#     - UMAP_facet_by_Age_Ascl1Status_CCStage_all.pdf
#     - UMAP_facet_by_Age_Ascl1Status_CCStage_no_NA.pdf
#     - Theta_analysis_final.pdf
#
#   Tables:
#     - frequency_age_Ascl1_CCStage.csv
#     - proportion_age_Ascl1_CCStage.csv
#     - UMAP_CCStage_by_Age_Ascl1_summary.csv
#     - Theta_stats_overall.csv
#     - Theta_stats_by_age.csv
#     - Theta_test_results.csv
#
#   Data:
#     - sce_tricycle_annotated.rds
#
#
# ==============================================================================

suppressPackageStartupMessages({
  library(tricycle)
  library(ggplot2)
  library(scattermore)
  library(scater)
  library(scuttle)
  library(org.Mm.eg.db)
  library(AnnotationDbi)
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(dplyr)
  library(circular)
  library(car)
  library(patchwork)
  library(schard)
})


# ==============================================================================
# CONFIGURATION - UPDATE THIS PATH
# ==============================================================================

H5AD_PATH <- "integrated_scvi.h5ad"

AGE_ORDER        <- c("young", "mid_age", "old", "pre_geriatric", "geriatric")
THETA_AGE_GROUPS <- c("old", "pre_geriatric", "geriatric")
CCSTAGE_ORDER    <- c("G1.S", "S", "G2", "G2.M", "M.G1", "NA_Cluster")

CCSTAGE_COLORS <- c(
  "G1.S"       = "#76D7C4",
  "S"          = "#2CA02C",
  "G2"         = "lightpink",
  "G2.M"       = "darkgreen",
  "M.G1"       = "darkred",
  "NA_Cluster" = "#7F7F7F"
)

CCSTAGE_LABELS <- c(
  "G1.S"       = "G1 to S (start DNA synthesis)",
  "S"          = "S (DNA replication)",
  "G2"         = "G2 (prep mitosis)",
  "G2.M"       = "G2 to M (enter mitosis)",
  "M.G1"       = "M to G1 (post-mitosis)",
  "NA_Cluster" = "Unassigned / NA"
)

ASCL1_FILL_COLORS <- c("Ascl1-" = "#4DBBD5", "Ascl1+" = "#E64B35")
ASCL1_LINE_COLORS <- c("Ascl1-" = "#3A9AB8", "Ascl1+" = "#C41E3A")


# ==============================================================================
# HELPERS
# ==============================================================================
banner <- function(text) {
  line <- paste(rep("=", 70), collapse = "")
  message("\n", line)
  message(text)
  message(line)
}


# ==============================================================================
# STEP 0: LOAD DATA AND SUBSET TO FEMALE HEPATOCYTES (Ascl1 annotated)
# ==============================================================================
banner("STEP 0: Load data and subset to female hepatocytes")

sce <- schard::h5ad2sce(H5AD_PATH)
message(sprintf("  [OK] Full dataset: %d cells x %d genes", ncol(sce), nrow(sce)))

# Subset to female hepatocytes
keep <- sce$sex == "female" & sce$celltype == "Hepatocyte"
sce <- sce[, keep]
message(sprintf("  [OK] Female hepatocytes: %d cells", ncol(sce)))

# Counts matrix
if ("counts" %in% assayNames(sce)) {
  # already present
} else if ("X" %in% assayNames(sce)) {
  assays(sce)$counts <- assays(sce)$X
  assays(sce)$X <- NULL
} else {
  assays(sce)$counts <- assay(sce, 1)
}
message(sprintf("  Assays: %s", paste(assayNames(sce), collapse = ", ")))

# Annotate Ascl1 status (pos / neg) from raw counts
if (!"Ascl1" %in% rownames(sce)) stop("Ascl1 gene not found in rownames(sce)")
ascl1_counts <- as.numeric(counts(sce)["Ascl1", ])
sce$Ascl1_status <- ifelse(ascl1_counts > 0, "Ascl1_pos", "Ascl1_neg")
message(sprintf("  [OK] Ascl1+: %d  |  Ascl1-: %d",
                sum(sce$Ascl1_status == "Ascl1_pos"),
                sum(sce$Ascl1_status == "Ascl1_neg")))

# Log-normalize
message("\n  Log-normalizing counts...")
sce <- scuttle::logNormCounts(sce)


# ==============================================================================
# STEP 1: TRICYCLE PROJECTION AND STAGE ESTIMATION
# ==============================================================================
banner("STEP 1: Tricycle projection, theta, and Schwabe stage")

# Convert gene symbols to Ensembl IDs (Tricycle reference uses Ensembl)
message("  Converting gene symbols to Ensembl IDs...")
gene_symbols <- rownames(sce)
ens_ids <- mapIds(
  org.Mm.eg.db,
  keys      = gene_symbols,
  column    = "ENSEMBL",
  keytype   = "SYMBOL",
  multiVals = "first"
)

valid_ids <- !is.na(ens_ids)
sce_filtered <- sce[valid_ids, ]
rownames(sce_filtered) <- ens_ids[valid_ids]
message(sprintf("  [OK] Genes mapped: %d / %d (%.1f%%)",
                sum(valid_ids), length(gene_symbols),
                100 * sum(valid_ids) / length(gene_symbols)))

# Project to cell cycle space, estimate position + Schwabe stage
message("\n  Running Tricycle project_cycle_space...")
sce_filtered <- project_cycle_space(sce_filtered)

message("  Running Tricycle estimate_cycle_position...")
sce_filtered <- estimate_cycle_position(sce_filtered)

message("  Running Tricycle estimate_Schwabe_stage...")
sce_filtered <- estimate_Schwabe_stage(
  sce_filtered,
  gname.type = "ENSEMBL",
  species    = "mouse"
)

message(sprintf("\n  [OK] Final SCE: %d cells x %d genes",
                ncol(sce_filtered), nrow(sce_filtered)))
message("\n  CCStage distribution:")
print(table(sce_filtered$CCStage, useNA = "ifany"))
message(sprintf("\n  Theta: mean=%.3f  median=%.3f  sd=%.3f",
                mean(sce_filtered$tricyclePosition, na.rm = TRUE),
                median(sce_filtered$tricyclePosition, na.rm = TRUE),
                sd(sce_filtered$tricyclePosition, na.rm = TRUE)))

# Embedding plot
p_embed <- scater::plotReducedDim(
  sce_filtered, dimred = "tricycleEmbedding", colour_by = "CCStage"
) +
  labs(
    x     = "Projected PC1",
    y     = "Projected PC2",
    title = sprintf("Projected cell cycle space (n = %d)", ncol(sce_filtered))
  ) +
  theme_bw(base_size = 14) +
  theme(text = element_text(family = "Arial"))

ggsave("tricycle_embedding_CCStage.pdf",
       p_embed, width = 8, height = 6, dpi = 300, useDingbats = FALSE)
message("  [OK] tricycle_embedding_CCStage.pdf")

saveRDS(sce_filtered, "sce_tricycle_annotated.rds")
message("  [OK] sce_tricycle_annotated.rds")


# ==============================================================================
# STEP 2: CCSTAGE DISTRIBUTION BY AGE x ASCL1 STATUS
# ==============================================================================
banner("STEP 2: CCStage distribution - stacked barplots")

# Clean and order metadata
sce_filtered$CCStage <- as.character(sce_filtered$CCStage)
sce_filtered$CCStage[is.na(sce_filtered$CCStage)] <- "NA_Cluster"
sce_filtered$CCStage <- factor(sce_filtered$CCStage, levels = CCSTAGE_ORDER, ordered = TRUE)

sce_filtered$Ascl1_status <- as.character(sce_filtered$Ascl1_status)
sce_filtered$Ascl1_status[is.na(sce_filtered$Ascl1_status)] <- "Unknown"
sce_filtered$Ascl1_status <- factor(
  sce_filtered$Ascl1_status,
  levels = c("Ascl1_pos", "Ascl1_neg", "Unknown")
)

sce_filtered$age <- factor(sce_filtered$age, levels = AGE_ORDER, ordered = TRUE)

# Combine age + Ascl1 status for x-axis grouping
age_ascl1_levels <- unlist(lapply(AGE_ORDER, function(a) {
  paste0(a, "_", c("Ascl1_pos", "Ascl1_neg"))
}))

sce_filtered$Age_Ascl1 <- paste0(sce_filtered$age, "_", sce_filtered$Ascl1_status)
sce_filtered$Age_Ascl1 <- factor(
  sce_filtered$Age_Ascl1,
  levels = age_ascl1_levels, ordered = TRUE
)

# Frequency and proportion tables
freq <- table(sce_filtered$CCStage, sce_filtered$Age_Ascl1, useNA = "ifany")
freq <- as.matrix(freq)
freq <- freq[intersect(CCSTAGE_ORDER, rownames(freq)),
             intersect(age_ascl1_levels, colnames(freq)), drop = FALSE]

prop <- apply(freq, 2, function(x) if (sum(x) == 0) x else x / sum(x))
prop <- as.matrix(prop)

message("  Cells by Age x Ascl1 status:")
print(table(sce_filtered$age, sce_filtered$Ascl1_status))

message("\n  CCStage by Ascl1 status (percentages):")
cc_by_ascl1 <- table(sce_filtered$Ascl1_status, sce_filtered$CCStage)
print(round(prop.table(cc_by_ascl1, margin = 1) * 100, 1))

bar_colors <- CCSTAGE_COLORS[rownames(freq)]

pdf("age_Ascl1_CCStage_histogram.pdf",
    height = 6, width = 14, useDingbats = FALSE)
layout(matrix(1:3, nrow = 1))
par(mar = c(9, 5, 2, 1), family = "Arial")

barplot(
  freq,
  col = bar_colors, border = NA,
  las = 2, ylab = "Frequency", cex.names = 0.9,
  main = "CCStage frequency by age x Ascl1 status"
)

barplot(
  prop,
  col = bar_colors, border = NA,
  las = 2, ylab = "Proportion", cex.names = 0.9,
  main = "CCStage proportion by age x Ascl1 status"
)

plot.new()
legend(
  "left",
  fill   = CCSTAGE_COLORS[CCSTAGE_ORDER],
  legend = CCSTAGE_LABELS[CCSTAGE_ORDER],
  bty = "n", cex = 0.9
)
dev.off()
message("\n  [OK] age_Ascl1_CCStage_histogram.pdf")

write.csv(freq, "frequency_age_Ascl1_CCStage.csv", row.names = TRUE)
write.csv(prop, "proportion_age_Ascl1_CCStage.csv", row.names = TRUE)
message("  [OK] frequency_age_Ascl1_CCStage.csv")
message("  [OK] proportion_age_Ascl1_CCStage.csv")


# ==============================================================================
# STEP 3: UMAP FACETED BY AGE x ASCL1 STATUS
# ==============================================================================
banner("STEP 3: UMAP faceted by age x Ascl1 status")

sce_filtered$CCStage <- factor(sce_filtered$CCStage, levels = CCSTAGE_ORDER)
sce_filtered$Ascl1_status <- factor(
  sce_filtered$Ascl1_status,
  levels = c("Ascl1_neg", "Ascl1_pos", "Unknown")
)

# Extract UMAP coordinates
if (!"X_umap" %in% reducedDimNames(sce_filtered)) {
  stop("reducedDim 'X_umap' not found. Make sure the h5ad has .obsm['X_umap'].")
}
umap_coords <- reducedDim(sce_filtered, "X_umap")
umap_df <- as.data.frame(umap_coords)
colnames(umap_df) <- c("UMAP_1", "UMAP_2")
umap_df$age          <- sce_filtered$age
umap_df$CCStage      <- sce_filtered$CCStage
umap_df$Ascl1_status <- sce_filtered$Ascl1_status

message(sprintf("  Cells in UMAP: %d", nrow(umap_df)))
na_count <- sum(umap_df$CCStage == "NA_Cluster")
message(sprintf("  NA_Cluster: %d cells (%.1f%%)",
                na_count, 100 * na_count / nrow(umap_df)))

umap_theme <- theme_minimal(base_size = 14) +
  theme(
    text            = element_text(family = "Arial"),
    legend.position = "bottom",
    legend.text     = element_text(size = 14),
    legend.title    = element_text(size = 16),
    strip.text      = element_text(face = "bold", size = 12),
    panel.spacing   = unit(0.8, "lines")
  )

# Panel 1: all stages
p_umap_all <- ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2, color = CCStage)) +
  geom_point(size = 0.7, alpha = 0.7) +
  scale_color_manual(values = CCSTAGE_COLORS) +
  facet_grid(age ~ Ascl1_status) +
  labs(
    title = "UMAP by age x Ascl1 status (all CCStages)",
    color = "Cell cycle stage"
  ) +
  umap_theme +
  guides(color = guide_legend(override.aes = list(size = 5)))

ggsave("UMAP_facet_by_Age_Ascl1Status_CCStage_all.pdf",
       p_umap_all, width = 12, height = 10, useDingbats = FALSE)
message("  [OK] UMAP_facet_by_Age_Ascl1Status_CCStage_all.pdf")

# Panel 2: drop NA_Cluster
umap_df_no_na <- umap_df[umap_df$CCStage != "NA_Cluster", ]
umap_df_no_na$CCStage <- droplevels(umap_df_no_na$CCStage)

p_umap_no_na <- ggplot(umap_df_no_na, aes(x = UMAP_1, y = UMAP_2, color = CCStage)) +
  geom_point(size = 0.8, alpha = 0.7) +
  scale_color_manual(values = CCSTAGE_COLORS[names(CCSTAGE_COLORS) != "NA_Cluster"]) +
  facet_grid(age ~ Ascl1_status) +
  labs(
    title = "UMAP by age x Ascl1 status (excluding NA cluster)",
    color = "Cell cycle stage"
  ) +
  umap_theme +
  guides(color = guide_legend(override.aes = list(size = 5)))

ggsave("UMAP_facet_by_Age_Ascl1Status_CCStage_no_NA.pdf",
       p_umap_no_na, width = 12, height = 10, useDingbats = FALSE)
message("  [OK] UMAP_facet_by_Age_Ascl1Status_CCStage_no_NA.pdf")

# Per-group summary table
summary_table <- as.data.frame(
  table(umap_df$age, umap_df$Ascl1_status, umap_df$CCStage)
)
colnames(summary_table) <- c("Age", "Ascl1_status", "CCStage", "Count")
summary_table <- summary_table[summary_table$Count > 0, ]
summary_table <- summary_table[order(summary_table$Age, summary_table$Ascl1_status), ]
write.csv(summary_table, "UMAP_CCStage_by_Age_Ascl1_summary.csv", row.names = FALSE)
message("  [OK] UMAP_CCStage_by_Age_Ascl1_summary.csv")


# ==============================================================================
# STEP 4: CONTINUOUS THETA ANALYSIS + CIRCULAR STATISTICS
# ==============================================================================
# Theta (tricyclePosition) ranges from 0 to 2*pi and encodes the continuous
# cell cycle position. Approximate phase mapping:
#   0    - 0.79  : G1/G0 (quiescent / early G1)
#   0.79 - 1.57  : Late G1
#   1.57 - 3.14  : S phase (DNA synthesis)
#   3.14 - 4.71  : G2/M phase
#   4.71 - 5.50  : M phase
#   5.50 - 6.28  : M -> G1 transition / G1/G0
# ==============================================================================
banner("STEP 4: Continuous theta analysis + circular statistics")

theta_df <- data.frame(
  theta        = sce_filtered$tricyclePosition,
  Ascl1_status = sce_filtered$Ascl1_status,
  age          = sce_filtered$age,
  CCStage      = sce_filtered$CCStage
)

theta_df <- theta_df %>%
  filter(
    Ascl1_status %in% c("Ascl1_pos", "Ascl1_neg"),
    age %in% THETA_AGE_GROUPS,
    !is.na(theta)
  ) %>%
  mutate(
    age_label = factor(
      age,
      levels = THETA_AGE_GROUPS,
      labels = c("Old", "Pre-geriatric", "Geriatric")
    ),
    Ascl1_label = factor(
      Ascl1_status,
      levels = c("Ascl1_neg", "Ascl1_pos"),
      labels = c("Ascl1-", "Ascl1+")
    )
  )

message(sprintf("  Cells: %d", nrow(theta_df)))
message(sprintf("    Ascl1-: %d", sum(theta_df$Ascl1_status == "Ascl1_neg")))
message(sprintf("    Ascl1+: %d", sum(theta_df$Ascl1_status == "Ascl1_pos")))
message("\n  By age:")
print(table(theta_df$age, theta_df$Ascl1_status))

# --- Descriptive stats ------------------------------------------------------
theta_stats_overall <- theta_df %>%
  group_by(Ascl1_status) %>%
  summarise(
    n            = n(),
    mean_theta   = mean(theta),
    median_theta = median(theta),
    sd_theta     = sd(theta),
    var_theta    = var(theta),
    IQR_theta    = IQR(theta),
    min_theta    = min(theta),
    max_theta    = max(theta),
    .groups = "drop"
  ) %>%
  mutate(
    mean_pi = mean_theta / pi,
    phase = case_when(
      mean_theta < 1.57 ~ "G1/G0 or late G1",
      mean_theta < 3.14 ~ "S phase",
      mean_theta < 4.71 ~ "G2/M phase",
      mean_theta < 5.50 ~ "M phase",
      TRUE              ~ "M to G1 / G1/G0"
    )
  )

message("\n  Theta statistics by Ascl1 status:")
print(as.data.frame(
  theta_stats_overall %>%
    select(Ascl1_status, n, mean_theta, median_theta, sd_theta, var_theta, phase) %>%
    mutate(across(where(is.numeric), ~round(., 3)))
))

theta_stats_by_age <- theta_df %>%
  group_by(age, Ascl1_status) %>%
  summarise(
    n            = n(),
    mean_theta   = mean(theta),
    median_theta = median(theta),
    sd_theta     = sd(theta),
    .groups = "drop"
  ) %>%
  mutate(
    mean_pi = mean_theta / pi,
    phase = case_when(
      mean_theta < 1.57 ~ "G1 / late G1",
      mean_theta < 3.14 ~ "S phase",
      mean_theta < 4.71 ~ "G2/M",
      TRUE              ~ "M / G1"
    )
  )

message("\n  Theta statistics by age x Ascl1 status:")
print(as.data.frame(
  theta_stats_by_age %>%
    select(age, Ascl1_status, n, mean_theta, sd_theta, phase) %>%
    mutate(across(where(is.numeric), ~round(., 3)))
))

# --- Statistical tests ------------------------------------------------------
banner("STATISTICAL TESTS")

# 1. Wilcoxon rank-sum
message("1. Wilcoxon rank-sum test")
wilcox_result <- wilcox.test(theta ~ Ascl1_status, data = theta_df, conf.int = TRUE)
n1 <- sum(theta_df$Ascl1_status == "Ascl1_pos")
n2 <- sum(theta_df$Ascl1_status == "Ascl1_neg")
r_effect <- 1 - (2 * wilcox_result$statistic) / (n1 * n2)
message(sprintf("   W statistic:     %.0f",  wilcox_result$statistic))
message(sprintf("   p-value:         %.2e",  wilcox_result$p.value))
message(sprintf("   Effect size (r): %.3f",  r_effect))

# 2. Levene's test for equality of variances
message("\n2. Levene's test (equality of variances)")
levene_result <- leveneTest(theta ~ Ascl1_status, data = theta_df)
var_neg <- var(theta_df$theta[theta_df$Ascl1_status == "Ascl1_neg"])
var_pos <- var(theta_df$theta[theta_df$Ascl1_status == "Ascl1_pos"])
var_ratio <- var_neg / var_pos
message(sprintf("   F statistic:     %.2f",  levene_result$`F value`[1]))
message(sprintf("   p-value:         %.2e",  levene_result$`Pr(>F)`[1]))
message(sprintf("   Var ratio (Ascl1-/Ascl1+): %.2f", var_ratio))

# 3. Circular statistics
message("\n3. Circular statistics")
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
circ_var_pos  <- var.circular(theta_pos_circ)
circ_var_neg  <- var.circular(theta_neg_circ)
rho_pos       <- rho.circular(theta_pos_circ)
rho_neg       <- rho.circular(theta_neg_circ)

message(sprintf("   Circular mean - Ascl1+: %.3f rad", as.numeric(circ_mean_pos)))
message(sprintf("   Circular mean - Ascl1-: %.3f rad", as.numeric(circ_mean_neg)))
message(sprintf("   Circular var  - Ascl1+: %.3f",     circ_var_pos))
message(sprintf("   Circular var  - Ascl1-: %.3f",     circ_var_neg))
message(sprintf("   Rho (conc.)   - Ascl1+: %.3f",     rho_pos))
message(sprintf("   Rho (conc.)   - Ascl1-: %.3f",     rho_neg))

# 4. Watson's two-sample test
message("\n4. Watson's two-sample test of homogeneity")
watson_result <- watson.two.test(theta_pos_circ, theta_neg_circ)
print(watson_result)
watson_stat <- tryCatch(
  as.numeric(watson_result$statistic),
  error = function(e) NA_real_
)


# --- Visualization ----------------------------------------------------------
banner("Generating theta distribution plots")

theme_publication <- theme_bw(base_size = 12) +
  theme(
    text             = element_text(family = "Arial"),
    plot.title       = element_text(face = "bold", size = 13, hjust = 0),
    plot.subtitle    = element_text(size = 9, hjust = 0, color = "gray40"),
    legend.position  = "bottom",
    legend.title     = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "gray90"),
    strip.text       = element_text(face = "bold", size = 11)
  )

axis_breaks <- c(0, 1.57, 3.14, 4.71, 6.28)
axis_labels <- c("0 (0)", "1.57 (0.5\u03C0)", "3.14 (\u03C0)",
                 "4.71 (1.5\u03C0)", "6.28 (2\u03C0)")

theta_means_overall <- theta_df %>%
  group_by(Ascl1_label) %>%
  summarise(mean_theta = mean(theta), .groups = "drop")

theta_means_by_age <- theta_df %>%
  group_by(age_label, Ascl1_label) %>%
  summarise(mean_theta = mean(theta), .groups = "drop")

# Shared phase shading
add_phase_shading <- function(p) {
  p +
    annotate("rect", xmin = 0,    xmax = 0.79, ymin = 0, ymax = Inf, alpha = 0.12, fill = "#FFC107") +
    annotate("rect", xmin = 5.50, xmax = 6.28, ymin = 0, ymax = Inf, alpha = 0.12, fill = "#FFC107") +
    annotate("rect", xmin = 0.79, xmax = 1.57, ymin = 0, ymax = Inf, alpha = 0.06, fill = "#FFC107") +
    annotate("rect", xmin = 1.57, xmax = 3.14, ymin = 0, ymax = Inf, alpha = 0.12, fill = "#4CAF50") +
    annotate("rect", xmin = 3.14, xmax = 4.71, ymin = 0, ymax = Inf, alpha = 0.12, fill = "#2196F3") +
    annotate("rect", xmin = 4.71, xmax = 5.50, ymin = 0, ymax = Inf, alpha = 0.12, fill = "#9C27B0")
}

# Panel A: overall distribution
pA <- ggplot(theta_df, aes(x = theta, fill = Ascl1_label, color = Ascl1_label))
pA <- add_phase_shading(pA) +
  geom_density(alpha = 0.5, linewidth = 1) +
  geom_vline(
    data = theta_means_overall,
    aes(xintercept = mean_theta, color = Ascl1_label),
    linetype = "dashed", linewidth = 1.2, show.legend = FALSE
  ) +
  geom_text(
    data = theta_means_overall,
    aes(x = mean_theta, y = Inf,
        label = sprintf("%.2f", mean_theta), color = Ascl1_label),
    vjust = 2.5, hjust = 0.5, size = 3.5, fontface = "bold", show.legend = FALSE
  ) +
  geom_vline(xintercept = c(0.79, 1.57, 3.14, 4.71, 5.50),
             linetype = "dotted", color = "gray40", linewidth = 0.5) +
  annotate("text", x = 0.40, y = Inf, vjust = 1.3, label = "G1/G0",
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
  scale_fill_manual(values = ASCL1_FILL_COLORS) +
  scale_color_manual(values = ASCL1_LINE_COLORS) +
  scale_x_continuous(
    breaks = axis_breaks, labels = axis_labels,
    limits = c(0, 6.28), expand = c(0.01, 0)
  ) +
  labs(
    title    = "A. Overall cell cycle position (theta) distribution",
    subtitle = "Dashed lines = group means; phase shading: G1/G0, S, G2/M, M",
    x        = "Cell cycle position theta (radians)",
    y        = "Density",
    fill     = "Ascl1 status",
    color    = "Ascl1 status"
  ) +
  theme_publication +
  theme(legend.position = "none", axis.text.x = element_text(size = 9)) +
  coord_cartesian(clip = "off")

# Panel B: age-stratified distribution
n_neg <- sum(theta_df$Ascl1_status == "Ascl1_neg")
n_pos <- sum(theta_df$Ascl1_status == "Ascl1_pos")

pB <- ggplot(theta_df, aes(x = theta, fill = Ascl1_label, color = Ascl1_label))
pB <- add_phase_shading(pB) +
  geom_density(alpha = 0.5, linewidth = 0.9) +
  geom_vline(
    data = theta_means_by_age,
    aes(xintercept = mean_theta, color = Ascl1_label),
    linetype = "dashed", linewidth = 1, show.legend = FALSE
  ) +
  geom_text(
    data = theta_means_by_age,
    aes(x = mean_theta, y = Inf,
        label = sprintf("%.2f", mean_theta), color = Ascl1_label),
    vjust = 2.5, hjust = 0.5, size = 3, fontface = "bold", show.legend = FALSE
  ) +
  geom_vline(xintercept = c(1.57, 3.14, 4.71, 5.50),
             linetype = "dotted", color = "gray40", linewidth = 0.4) +
  facet_wrap(~ age_label, ncol = 3) +
  scale_fill_manual(
    values = ASCL1_FILL_COLORS,
    labels = c(
      sprintf("Ascl1- (n=%s)", format(n_neg, big.mark = ",")),
      sprintf("Ascl1+ (n=%s)", format(n_pos, big.mark = ","))
    )
  ) +
  scale_color_manual(values = ASCL1_LINE_COLORS) +
  scale_x_continuous(
    breaks = axis_breaks, labels = axis_labels,
    limits = c(0, 6.28), expand = c(0.01, 0)
  ) +
  labs(
    title    = "B. Cell cycle position by age group",
    subtitle = "S start = 1.57  |  G2/M start = 3.14  |  M middle = 4.71",
    x        = "Cell cycle position theta (radians)",
    y        = "Density",
    fill     = "Ascl1 status",
    color    = "Ascl1 status"
  ) +
  theme_publication +
  theme(legend.position = "bottom", axis.text.x = element_text(size = 8)) +
  guides(color = "none")

# Combine panels
mean_pos_val <- theta_stats_overall$mean_theta[theta_stats_overall$Ascl1_status == "Ascl1_pos"]
mean_neg_val <- theta_stats_overall$mean_theta[theta_stats_overall$Ascl1_status == "Ascl1_neg"]

p_final <- pA / pB +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(
    title = "Cell cycle position (theta): Ascl1+ vs Ascl1- female hepatocytes",
    caption = sprintf(
      "Ascl1+ mean theta = %.2f | Ascl1- mean theta = %.2f | Wilcoxon p = %.2e",
      mean_pos_val, mean_neg_val, wilcox_result$p.value
    ),
    theme = theme(
      text         = element_text(family = "Arial"),
      plot.title   = element_text(face = "bold", size = 15, hjust = 0.5),
      plot.caption = element_text(size = 9, hjust = 0.5, face = "italic", color = "gray30")
    )
  )

ggsave("Theta_analysis_final.pdf",
       p_final, width = 13, height = 10, useDingbats = FALSE)
message("  [OK] Theta_analysis_final.pdf")


# --- Save results tables ----------------------------------------------------
write.csv(theta_stats_overall, "Theta_stats_overall.csv", row.names = FALSE)
write.csv(theta_stats_by_age,  "Theta_stats_by_age.csv",  row.names = FALSE)

test_results <- data.frame(
  Test = c(
    "Wilcoxon Rank-Sum", "Wilcoxon Rank-Sum", "Wilcoxon Rank-Sum",
    "Levene's Test", "Levene's Test", "Levene's Test",
    "Watson's Circular", "Circular Stats", "Circular Stats",
    "Circular Stats", "Circular Stats"
  ),
  Metric = c(
    "W statistic", "p-value", "Effect size (r)",
    "F statistic", "p-value", "Variance ratio",
    "Test statistic",
    "Rho (Ascl1+)", "Rho (Ascl1-)",
    "Circular var (Ascl1+)", "Circular var (Ascl1-)"
  ),
  Value = c(
    sprintf("%.0f", wilcox_result$statistic),
    sprintf("%.2e", wilcox_result$p.value),
    sprintf("%.3f", r_effect),
    sprintf("%.2f", levene_result$`F value`[1]),
    sprintf("%.2e", levene_result$`Pr(>F)`[1]),
    sprintf("%.2f", var_ratio),
    ifelse(is.na(watson_stat), "see Watson test output", sprintf("%.3f", watson_stat)),
    sprintf("%.3f", rho_pos),
    sprintf("%.3f", rho_neg),
    sprintf("%.3f", circ_var_pos),
    sprintf("%.3f", circ_var_neg)
  )
)
write.csv(test_results, "Theta_test_results.csv", row.names = FALSE)
message("  [OK] Theta_stats_overall.csv")
message("  [OK] Theta_stats_by_age.csv")
message("  [OK] Theta_test_results.csv")


# ==============================================================================
# DONE
# ==============================================================================
banner("ANALYSIS COMPLETE")
message("Figures:")
message("  tricycle_embedding_CCStage.pdf")
message("  age_Ascl1_CCStage_histogram.pdf")
message("  UMAP_facet_by_Age_Ascl1Status_CCStage_all.pdf")
message("  UMAP_facet_by_Age_Ascl1Status_CCStage_no_NA.pdf")
message("  Theta_analysis_final.pdf")
message("\nTables:")
message("  frequency_age_Ascl1_CCStage.csv")
message("  proportion_age_Ascl1_CCStage.csv")
message("  UMAP_CCStage_by_Age_Ascl1_summary.csv")
message("  Theta_stats_overall.csv")
message("  Theta_stats_by_age.csv")
message("  Theta_test_results.csv")
message("\nData:")
message("  sce_tricycle_annotated.rds")

sessionInfo()
