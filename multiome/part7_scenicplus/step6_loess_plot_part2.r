#!/usr/bin/env Rscript
# =============================================================================
# TF RSS Trajectory Clustering Analysis - MALE
# Input: age_sex_rss_matrix_gene_male_clean.csv
#        age_sex_rss_matrix_region_male_clean.csv
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(reshape2)
  library(ggrepel)
  library(patchwork)
})

set.seed(123)

# =============================================================================
# Constants
# =============================================================================
AGE_LEVELS <- c("young", "mid_age", "old", "pre_geriatric", "geriatric")
N_CLUSTERS <- 6
COLORS <- list(gene = "#1B9E77", region = "#D95F02")

# =============================================================================
# SECTION 1: Load Data
# =============================================================================
cat("=", rep("=", 59), "\n", sep = "")
cat("Loading cleaned RSS matrices (Male)\n")
cat("=", rep("=", 59), "\n\n", sep = "")

gene_mat <- read.csv("age_sex_rss_matrix_gene_male_clean.csv", 
                     row.names = 1, check.names = FALSE)
region_mat <- read.csv("age_sex_rss_matrix_region_male_clean.csv", 
                       row.names = 1, check.names = FALSE)

cat("Gene matrix:", dim(gene_mat)[1], "ages ×", dim(gene_mat)[2], "TFs\n")
cat("Region matrix:", dim(region_mat)[1], "ages ×", dim(region_mat)[2], "TFs\n")

# Match columns
common_tfs <- intersect(colnames(gene_mat), colnames(region_mat))
gene_mat <- gene_mat[, common_tfs, drop = FALSE]
region_mat <- region_mat[, common_tfs, drop = FALSE]
cat("Common TFs:", length(common_tfs), "\n\n")

# =============================================================================
# SECTION 2: Z-score Normalization
# =============================================================================
cat("Z-score normalization...\n")

z_gene <- t(scale(t(as.matrix(gene_mat))))
z_region <- t(scale(t(as.matrix(region_mat))))

# Transpose: TFs as rows × Ages as columns
z_gene <- t(z_gene)
z_region <- t(z_region)

cat("✅ Z-scoring complete. Shape:", dim(z_gene)[1], "TFs ×", dim(z_gene)[2], "ages\n\n")

# =============================================================================
# SECTION 3: K-means Clustering (Gene-based)
# =============================================================================
cat("K-means clustering (k =", N_CLUSTERS, ")...\n")

km <- kmeans(z_gene, centers = N_CLUSTERS)
gene_clusters <- data.frame(TF = rownames(z_gene), Cluster = factor(km$cluster))

cat("✅ Clustering complete\n\n")

# =============================================================================
# SECTION 4: Convert to Long Format
# =============================================================================
cat("Converting to long format...\n")

# Gene
gene_df <- as.data.frame(z_gene)
gene_df$TF <- rownames(gene_df)
gene_long <- melt(gene_df, id.vars = "TF", variable.name = "Age", value.name = "z_RSS")
gene_long$Source <- "Gene"

# Region
region_df <- as.data.frame(z_region)
region_df$TF <- rownames(region_df)
region_long <- melt(region_df, id.vars = "TF", variable.name = "Age", value.name = "z_RSS")
region_long$Source <- "Region"

# Combine
rss_long <- rbind(gene_long, region_long)
rss_long$Age <- factor(rss_long$Age, levels = AGE_LEVELS)
rss_long <- left_join(rss_long, gene_clusters, by = "TF")

cat("✅ Long format:", nrow(rss_long), "rows,", 
    length(unique(rss_long$TF)), "TFs,",
    length(unique(rss_long$Cluster)), "clusters\n\n")

# =============================================================================
# SECTION 5: Plot 1 - Individual TF Trajectories
# =============================================================================
cat("=", rep("=", 59), "\n", sep = "")
cat("Generating plots...\n")
cat("=", rep("=", 59), "\n\n")

last_age <- "geriatric"

# Gene-based trajectories
p1 <- ggplot(subset(rss_long, Source == "Gene"), 
             aes(x = Age, y = z_RSS, group = TF, color = TF)) +
  geom_line(linewidth = 0.8, alpha = 0.85) +
  geom_point(size = 1.6) +
  geom_text_repel(
    data = subset(rss_long, Source == "Gene" & Age == last_age),
    aes(label = TF),
    size = 3,
    segment.color = "gray60",
    max.overlaps = 100
  ) +
  facet_wrap(~ Cluster, scales = "free_y") +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold", size = 16)
  ) +
  labs(
    title = "TF Trajectories by Cluster (Gene-based, Male)",
    x = "Age", y = "z-scored RSS"
  )

ggsave("male_tf_trajectories_gene.pdf", plot = p1, width = 30, height = 7)
cat("✅ Saved: male_tf_trajectories_gene.pdf\n")

# Region-based trajectories
p2 <- ggplot(subset(rss_long, Source == "Region"), 
             aes(x = Age, y = z_RSS, group = TF, color = TF)) +
  geom_line(linewidth = 0.8, alpha = 0.85) +
  geom_point(size = 1.6) +
  geom_text_repel(
    data = subset(rss_long, Source == "Region" & Age == last_age),
    aes(label = TF),
    size = 3,
    segment.color = "gray60",
    max.overlaps = 100
  ) +
  facet_wrap(~ Cluster, scales = "free_y") +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold", size = 16)
  ) +
  labs(
    title = "TF Trajectories by Cluster (Region-based, Male)",
    x = "Age", y = "z-scored RSS"
  )

ggsave("male_tf_trajectories_region.pdf", plot = p2, width = 30, height = 7)
cat("✅ Saved: male_tf_trajectories_region.pdf\n")

# =============================================================================
# SECTION 6: Plot 2 - Gene vs Region Overlay
# =============================================================================
p3 <- ggplot(rss_long, aes(x = Age, y = z_RSS, group = TF, color = TF)) +
  geom_line(data = subset(rss_long, Source == "Gene"), linewidth = 1.1) +
  geom_line(data = subset(rss_long, Source == "Region"), 
            linewidth = 1, linetype = "dashed") +
  facet_wrap(~ Cluster, scales = "free_y") +
  theme_minimal(base_size = 16) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold", size = 18)
  ) +
  labs(
    title = "Gene (solid) vs Region (dashed) - Male",
    x = "Age group", y = "z-scored RSS"
  )

ggsave("male_gene_region_overlay.pdf", plot = p3, width = 30, height = 7)
cat("✅ Saved: male_gene_region_overlay.pdf\n")

# =============================================================================
# SECTION 7: Plot 3 - Mean Cluster Trajectories
# =============================================================================
summary_df <- rss_long %>%
  group_by(Cluster, Age, Source) %>%
  summarise(
    mean_zRSS = mean(z_RSS, na.rm = TRUE),
    sd_zRSS   = sd(z_RSS, na.rm = TRUE),
    n         = n(),
    .groups   = "drop"
  ) %>%
  mutate(
    se    = sd_zRSS / sqrt(n),
    lower = mean_zRSS - 1.96 * se,
    upper = mean_zRSS + 1.96 * se
  )

p4 <- ggplot(summary_df, aes(x = Age, y = mean_zRSS, 
                              group = Source, color = Source, linetype = Source)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = Source),
              alpha = 0.15, color = NA) +
  geom_smooth(method = "loess", se = FALSE, linewidth = 1.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  facet_wrap(~ Cluster, scales = "free_y") +
  scale_linetype_manual(values = c(Gene = "solid", Region = "dashed")) +
  scale_color_manual(values = c(Gene = COLORS$gene, Region = COLORS$region)) +
  scale_fill_manual(values = c(Gene = COLORS$gene, Region = COLORS$region)) +
  labs(
    title = "Mean RSS per Cluster (Gene vs Region) - Male",
    x = "Age group", y = "Mean z-scored RSS"
  ) +
  theme_minimal(base_size = 15) +
  theme(
    strip.text = element_text(face = "bold", size = 16),
    legend.position = "bottom",
    legend.title = element_blank()
  )

ggsave("male_cluster_means.pdf", plot = p4, width = 14, height = 7)
cat("✅ Saved: male_cluster_means.pdf\n")

# =============================================================================
# SECTION 8: Plot 4 - Clusters with TF Labels (Symmetric Scales)
# =============================================================================

# Compute symmetric limits per cluster
range_df <- summary_df %>%
  group_by(Cluster) %>%
  summarise(
    max_abs = max(abs(c(lower, upper)), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    max_abs = ifelse(is.finite(max_abs), max_abs, 0.5),
    ymax = ceiling(max_abs * 1.1),
    ymin = -ymax
  )

# Create TF labels
label_df <- rss_long %>%
  distinct(Cluster, TF) %>%
  group_by(Cluster) %>%
  summarise(
    TF_list = paste(sort(unique(TF)), collapse = "\n"),
    n_TFs = n(),
    .groups = "drop"
  ) %>%
  mutate(Age = "geriatric", mean_zRSS = 0)

# Base plot template
base_plot <- ggplot(summary_df, aes(x = Age, y = mean_zRSS,
                                     group = Source, color = Source, linetype = Source)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = Source),
              alpha = 0.15, color = NA) +
  geom_smooth(method = "loess", se = FALSE, linewidth = 1.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  scale_linetype_manual(values = c(Gene = "solid", Region = "dashed")) +
  scale_color_manual(values = c(Gene = COLORS$gene, Region = COLORS$region)) +
  scale_fill_manual(values = c(Gene = COLORS$gene, Region = COLORS$region)) +
  labs(x = "Age group", y = "Mean z-scored RSS") +
  theme_minimal(base_size = 15) +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    strip.text = element_text(face = "bold", size = 20)
  )

# Build each cluster facet
p_list <- list()
clusters <- sort(unique(summary_df$Cluster))

for (cl in clusters) {
  lims <- range_df[range_df$Cluster == cl, c("ymin", "ymax")]
  
  p_list[[as.character(cl)]] <-
    base_plot %+%
    filter(summary_df, Cluster == cl) +
    geom_text(
      data = filter(label_df, Cluster == cl),
      aes(x = Age, y = mean_zRSS, label = TF_list),
      inherit.aes = FALSE,
      hjust = 0, vjust = 0.5,
      nudge_x = 0.25, size = 3.2, lineheight = 0.9, color = "black"
    ) +
    coord_cartesian(ylim = c(lims$ymin, lims$ymax), clip = "off") +
    ggtitle(paste("Cluster", cl))
}

# Combine all clusters
p5 <- wrap_plots(p_list, ncol = N_CLUSTERS) +
  plot_annotation(
    title = "Mean z-scored RSS per Cluster (Male)",
    theme = theme(plot.title = element_text(face = "bold", size = 18, hjust = 0.5))
  )

ggsave("male_cluster_means_labeled.pdf", plot = p5, 
       width = N_CLUSTERS * 6.5, height = 8, device = cairo_pdf)
cat("✅ Saved: male_cluster_means_labeled.pdf\n")

# =============================================================================
# SECTION 9: Save Cluster Summary
# =============================================================================
cluster_summary <- gene_clusters %>%
  group_by(Cluster) %>%
  summarise(
    n_TFs = n(),
    TFs = paste(sort(TF), collapse = ", "),
    .groups = "drop"
  )

write.csv(cluster_summary, "male_cluster_summary.csv", row.names = FALSE)
cat("✅ Saved: male_cluster_summary.csv\n")

# =============================================================================
# SECTION 10: Summary
# =============================================================================
cat("\n")
cat("=", rep("=", 59), "\n", sep = "")
cat("ANALYSIS COMPLETE - MALE\n")
cat("=", rep("=", 59), "\n\n", sep = "")

cat("Output files:\n")
cat("  • male_tf_trajectories_gene.pdf\n")
cat("  • male_tf_trajectories_region.pdf\n")
cat("  • male_gene_region_overlay.pdf\n")
cat("  • male_cluster_means.pdf\n")
cat("  • male_cluster_means_labeled.pdf\n")
cat("  • male_cluster_summary.csv\n\n")

cat("Cluster distribution:\n")
print(table(gene_clusters$Cluster))

cat("\nDone! 🎉\n")

#!/usr/bin/env Rscript
# =============================================================================
# TF RSS Trajectory Clustering Analysis - FEMALE
# Input: age_sex_rss_matrix_gene_female_clean.csv
#        age_sex_rss_matrix_region_female_clean.csv
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(reshape2)
  library(ggrepel)
  library(patchwork)
})

set.seed(123)

# =============================================================================
# Constants
# =============================================================================
AGE_LEVELS <- c("young", "mid_age", "old", "pre_geriatric", "geriatric")
N_CLUSTERS <- 6
COLORS <- list(gene = "#1B9E77", region = "#D95F02")

# =============================================================================
# SECTION 1: Load Data
# =============================================================================
cat("=", rep("=", 59), "\n", sep = "")
cat("Loading cleaned RSS matrices (Female)\n")
cat("=", rep("=", 59), "\n\n", sep = "")

gene_mat <- read.csv("age_sex_rss_matrix_gene_female_clean.csv", 
                     row.names = 1, check.names = FALSE)
region_mat <- read.csv("age_sex_rss_matrix_region_female_clean.csv", 
                       row.names = 1, check.names = FALSE)

cat("Gene matrix:", dim(gene_mat)[1], "ages ×", dim(gene_mat)[2], "TFs\n")
cat("Region matrix:", dim(region_mat)[1], "ages ×", dim(region_mat)[2], "TFs\n")

# Match columns
common_tfs <- intersect(colnames(gene_mat), colnames(region_mat))
gene_mat <- gene_mat[, common_tfs, drop = FALSE]
region_mat <- region_mat[, common_tfs, drop = FALSE]
cat("Common TFs:", length(common_tfs), "\n\n")

# =============================================================================
# SECTION 2: Z-score Normalization
# =============================================================================
cat("Z-score normalization...\n")

z_gene <- t(scale(t(as.matrix(gene_mat))))
z_region <- t(scale(t(as.matrix(region_mat))))

# Transpose: TFs as rows × Ages as columns
z_gene <- t(z_gene)
z_region <- t(z_region)

cat("✅ Z-scoring complete. Shape:", dim(z_gene)[1], "TFs ×", dim(z_gene)[2], "ages\n\n")

# =============================================================================
# SECTION 3: K-means Clustering (Gene-based)
# =============================================================================
cat("K-means clustering (k =", N_CLUSTERS, ")...\n")

km <- kmeans(z_gene, centers = N_CLUSTERS)
gene_clusters <- data.frame(TF = rownames(z_gene), Cluster = factor(km$cluster))

cat("✅ Clustering complete\n\n")

# =============================================================================
# SECTION 4: Convert to Long Format
# =============================================================================
cat("Converting to long format...\n")

# Gene
gene_df <- as.data.frame(z_gene)
gene_df$TF <- rownames(gene_df)
gene_long <- melt(gene_df, id.vars = "TF", variable.name = "Age", value.name = "z_RSS")
gene_long$Source <- "Gene"

# Region
region_df <- as.data.frame(z_region)
region_df$TF <- rownames(region_df)
region_long <- melt(region_df, id.vars = "TF", variable.name = "Age", value.name = "z_RSS")
region_long$Source <- "Region"

# Combine
rss_long <- rbind(gene_long, region_long)
rss_long$Age <- factor(rss_long$Age, levels = AGE_LEVELS)
rss_long <- left_join(rss_long, gene_clusters, by = "TF")

cat("✅ Long format:", nrow(rss_long), "rows,", 
    length(unique(rss_long$TF)), "TFs,",
    length(unique(rss_long$Cluster)), "clusters\n\n")

# =============================================================================
# SECTION 5: Plot 1 - Individual TF Trajectories
# =============================================================================
cat("=", rep("=", 59), "\n", sep = "")
cat("Generating plots...\n")
cat("=", rep("=", 59), "\n\n", sep = "")

last_age <- "geriatric"

# Gene-based trajectories
p1 <- ggplot(subset(rss_long, Source == "Gene"), 
             aes(x = Age, y = z_RSS, group = TF, color = TF)) +
  geom_line(linewidth = 0.8, alpha = 0.85) +
  geom_point(size = 1.6) +
  geom_text_repel(
    data = subset(rss_long, Source == "Gene" & Age == last_age),
    aes(label = TF),
    size = 3,
    segment.color = "gray60",
    max.overlaps = 100
  ) +
  facet_wrap(~ Cluster, scales = "free_y") +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold", size = 16)
  ) +
  labs(
    title = "TF Trajectories by Cluster (Gene-based, Female)",
    x = "Age", y = "z-scored RSS"
  )

ggsave("female_tf_trajectories_gene.pdf", plot = p1, width = 30, height = 7)
cat("✅ Saved: female_tf_trajectories_gene.pdf\n")

# Region-based trajectories
p2 <- ggplot(subset(rss_long, Source == "Region"), 
             aes(x = Age, y = z_RSS, group = TF, color = TF)) +
  geom_line(linewidth = 0.8, alpha = 0.85) +
  geom_point(size = 1.6) +
  geom_text_repel(
    data = subset(rss_long, Source == "Region" & Age == last_age),
    aes(label = TF),
    size = 3,
    segment.color = "gray60",
    max.overlaps = 100
  ) +
  facet_wrap(~ Cluster, scales = "free_y") +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold", size = 16)
  ) +
  labs(
    title = "TF Trajectories by Cluster (Region-based, Female)",
    x = "Age", y = "z-scored RSS"
  )

ggsave("female_tf_trajectories_region.pdf", plot = p2, width = 30, height = 7)
cat("✅ Saved: female_tf_trajectories_region.pdf\n")

# =============================================================================
# SECTION 6: Plot 2 - Gene vs Region Overlay
# =============================================================================
p3 <- ggplot(rss_long, aes(x = Age, y = z_RSS, group = TF, color = TF)) +
  geom_line(data = subset(rss_long, Source == "Gene"), linewidth = 1.1) +
  geom_line(data = subset(rss_long, Source == "Region"), 
            linewidth = 1, linetype = "dashed") +
  facet_wrap(~ Cluster, scales = "free_y") +
  theme_minimal(base_size = 16) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold", size = 18)
  ) +
  labs(
    title = "Gene (solid) vs Region (dashed) - Female",
    x = "Age group", y = "z-scored RSS"
  )

ggsave("female_gene_region_overlay.pdf", plot = p3, width = 30, height = 7)
cat("✅ Saved: female_gene_region_overlay.pdf\n")

# =============================================================================
# SECTION 7: Plot 3 - Mean Cluster Trajectories
# =============================================================================
summary_df <- rss_long %>%
  group_by(Cluster, Age, Source) %>%
  summarise(
    mean_zRSS = mean(z_RSS, na.rm = TRUE),
    sd_zRSS   = sd(z_RSS, na.rm = TRUE),
    n         = n(),
    .groups   = "drop"
  ) %>%
  mutate(
    se    = sd_zRSS / sqrt(n),
    lower = mean_zRSS - 1.96 * se,
    upper = mean_zRSS + 1.96 * se
  )

p4 <- ggplot(summary_df, aes(x = Age, y = mean_zRSS, 
                              group = Source, color = Source, linetype = Source)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = Source),
              alpha = 0.15, color = NA) +
  geom_smooth(method = "loess", se = FALSE, linewidth = 1.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  facet_wrap(~ Cluster, scales = "free_y") +
  scale_linetype_manual(values = c(Gene = "solid", Region = "dashed")) +
  scale_color_manual(values = c(Gene = COLORS$gene, Region = COLORS$region)) +
  scale_fill_manual(values = c(Gene = COLORS$gene, Region = COLORS$region)) +
  labs(
    title = "Mean RSS per Cluster (Gene vs Region) - Female",
    x = "Age group", y = "Mean z-scored RSS"
  ) +
  theme_minimal(base_size = 15) +
  theme(
    strip.text = element_text(face = "bold", size = 16),
    legend.position = "bottom",
    legend.title = element_blank()
  )

ggsave("female_cluster_means.pdf", plot = p4, width = 14, height = 7)
cat("✅ Saved: female_cluster_means.pdf\n")

# =============================================================================
# SECTION 8: Plot 4 - Clusters with TF Labels (Symmetric Scales)
# =============================================================================

# Compute symmetric limits per cluster
range_df <- summary_df %>%
  group_by(Cluster) %>%
  summarise(
    max_abs = max(abs(c(lower, upper)), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    max_abs = ifelse(is.finite(max_abs), max_abs, 0.5),
    ymax = ceiling(max_abs * 1.1),
    ymin = -ymax
  )

# Create TF labels
label_df <- rss_long %>%
  distinct(Cluster, TF) %>%
  group_by(Cluster) %>%
  summarise(
    TF_list = paste(sort(unique(TF)), collapse = "\n"),
    n_TFs = n(),
    .groups = "drop"
  ) %>%
  mutate(Age = "geriatric", mean_zRSS = 0)

# Base plot template
base_plot <- ggplot(summary_df, aes(x = Age, y = mean_zRSS,
                                     group = Source, color = Source, linetype = Source)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = Source),
              alpha = 0.15, color = NA) +
  geom_smooth(method = "loess", se = FALSE, linewidth = 1.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  scale_linetype_manual(values = c(Gene = "solid", Region = "dashed")) +
  scale_color_manual(values = c(Gene = COLORS$gene, Region = COLORS$region)) +
  scale_fill_manual(values = c(Gene = COLORS$gene, Region = COLORS$region)) +
  labs(x = "Age group", y = "Mean z-scored RSS") +
  theme_minimal(base_size = 15) +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    strip.text = element_text(face = "bold", size = 20)
  )

# Build each cluster facet
p_list <- list()
clusters <- sort(unique(summary_df$Cluster))

for (cl in clusters) {
  lims <- range_df[range_df$Cluster == cl, c("ymin", "ymax")]
  
  p_list[[as.character(cl)]] <-
    base_plot %+%
    filter(summary_df, Cluster == cl) +
    geom_text(
      data = filter(label_df, Cluster == cl),
      aes(x = Age, y = mean_zRSS, label = TF_list),
      inherit.aes = FALSE,
      hjust = 0, vjust = 0.5,
      nudge_x = 0.25, size = 3.2, lineheight = 0.9, color = "black"
    ) +
    coord_cartesian(ylim = c(lims$ymin, lims$ymax), clip = "off") +
    ggtitle(paste("Cluster", cl))
}

# Combine all clusters
p5 <- wrap_plots(p_list, ncol = N_CLUSTERS) +
  plot_annotation(
    title = "Mean z-scored RSS per Cluster (Female)",
    theme = theme(plot.title = element_text(face = "bold", size = 18, hjust = 0.5))
  )

ggsave("female_cluster_means_labeled.pdf", plot = p5, 
       width = N_CLUSTERS * 6.5, height = 8, device = cairo_pdf)
cat("✅ Saved: female_cluster_means_labeled.pdf\n")

# =============================================================================
# SECTION 9: Save Cluster Summary
# =============================================================================
cluster_summary <- gene_clusters %>%
  group_by(Cluster) %>%
  summarise(
    n_TFs = n(),
    TFs = paste(sort(TF), collapse = ", "),
    .groups = "drop"
  )

write.csv(cluster_summary, "female_cluster_summary.csv", row.names = FALSE)
cat("✅ Saved: female_cluster_summary.csv\n")

# =============================================================================
# SECTION 10: Summary
# =============================================================================
cat("\n")
cat("=", rep("=", 59), "\n", sep = "")
cat("ANALYSIS COMPLETE - FEMALE\n")
cat("=", rep("=", 59), "\n\n", sep = "")

cat("Output files:\n")
cat("  • female_tf_trajectories_gene.pdf\n")
cat("  • female_tf_trajectories_region.pdf\n")
cat("  • female_gene_region_overlay.pdf\n")
cat("  • female_cluster_means.pdf\n")
cat("  • female_cluster_means_labeled.pdf\n")
cat("  • female_cluster_summary.csv\n\n")

cat("Cluster distribution:\n")
print(table(gene_clusters$Cluster))

cat("\nDone! 🎉\n")
