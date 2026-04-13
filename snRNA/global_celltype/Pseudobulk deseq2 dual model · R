#!/usr/bin/env Rscript
# ==============================================================================
# Pseudobulk DESeq2 Dual-Model LRT Comparison (Hepatocyte)
# ==============================================================================
#
# Description:
#   Runs two DESeq2 Likelihood Ratio Tests side-by-side on sample-level
#   pseudobulk counts from one cell type (default: Hepatocyte) to compare:
#
#     Model 1 (with sex): full = ~ sex + age   reduced = ~ sex
#                         -> tests age effect controlling for sex
#     Model 2 (age only): full = ~ age          reduced = ~ 1
#                         -> tests age effect without sex adjustment
#
#   For each model it produces: LRT results, pairwise age contrasts, a VST
#   PCA, and a ComplexHeatmap of age-significant genes (k-means clustered)
#   with a per-cluster line plot of mean expression across age split by sex.
#   Finally it compares the two significant-gene sets (overlap + Venn diagram).
#
# Input:
#   - Annotated Seurat .rds OR AnnData .h5ad with celltype, sex, age, sample
#
# Output:
#   deseq2_<celltype>/
#     model_with_sex/
#       deseq2_LRT_age_results.csv
#       deseq2_LRT_age_significant.csv
#       deseq2_age_<Age2>_vs_<Age1>.csv   (all pairwise contrasts)
#       deseq2_all_pairwise_contrasts.rds
#       PCA_age_sex.pdf
#       heatmap_zscore.pdf
#       genes_cluster_assignments.csv
#       genes_cluster_<k>.csv
#       genes_clusters_with_stats.csv
#       lineplot_clusters_by_sex.pdf
#     model_without_sex/
#       (same set of outputs)
#     model_comparison_summary.csv
#     sig_gene_overlap.csv
#     venn_sig_overlap.png
#
# Reference:
#   Sarker et al. (2026) Cell Metabolism
#
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(DESeq2)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(ComplexHeatmap)
  library(circlize)
  library(RColorBrewer)
  library(VennDiagram)
  library(schard)
})


# ==============================================================================
# CONFIGURATION - UPDATE THIS PATH
# ==============================================================================
INPUT_PATH      <- "integrated_scvi.h5ad"   # .h5ad or .rds accepted
TARGET_CELLTYPE <- "Hepatocyte"

BASE_OUT_DIR    <- paste0("deseq2_", tolower(TARGET_CELLTYPE))
if (!dir.exists(BASE_OUT_DIR)) dir.create(BASE_OUT_DIR, recursive = TRUE)

AGE_LEVELS <- c("young", "mid_age", "old", "pre_geriatric", "geriatric")

AGE_COLORS <- c(
  "young"         = "#1ABC9C",
  "mid_age"       = "#F1C40F",
  "old"           = "#C39BD3",
  "pre_geriatric" = "#2980B9",
  "geriatric"     = "#E84393"
)
SEX_COLORS <- c("male" = "navy", "female" = "deeppink")

K_CLUSTERS   <- 4
MIN_COUNT    <- 10     # row filter: at least MIN_COUNT in >= MIN_SAMPLES
MIN_SAMPLES  <- 2
SEED         <- 123    # k-means reproducibility


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
# STEP 0: LOAD DATA AND SUBSET TO TARGET CELL TYPE
# ==============================================================================
banner(sprintf("STEP 0: Load data and subset to %s", TARGET_CELLTYPE))

if (grepl("\\.h5ad$", INPUT_PATH, ignore.case = TRUE)) {
  message("  Loading via schard::h5ad2seurat: ", INPUT_PATH)
  seuratObj <- schard::h5ad2seurat(INPUT_PATH)
} else if (grepl("\\.rds$", INPUT_PATH, ignore.case = TRUE)) {
  message("  Loading RDS: ", INPUT_PATH)
  seuratObj <- readRDS(INPUT_PATH)
} else {
  stop("INPUT_PATH must end in .h5ad or .rds")
}

data <- subset(seuratObj, subset = celltype == TARGET_CELLTYPE)
DefaultAssay(data) <- "RNA"
message(sprintf("  [OK] %s cells: %d", TARGET_CELLTYPE, ncol(data)))


# ==============================================================================
# STEP 1: PREPARE SAMPLE-LEVEL METADATA
# ==============================================================================
banner("STEP 1: Prepare metadata")

meta <- data@meta.data
required <- c("sample", "age", "sex")
missing  <- setdiff(required, colnames(meta))
if (length(missing) > 0) {
  stop("Missing metadata columns: ", paste(missing, collapse = ", "))
}
meta$age <- factor(meta$age, levels = AGE_LEVELS)
meta$group_id <- meta$sample
message(sprintf("  [OK] %d cells, %d samples", nrow(meta),
                length(unique(meta$sample))))


# ==============================================================================
# STEP 2: PSEUDOBULK AGGREGATION
# ==============================================================================
banner("STEP 2: Pseudobulk aggregation (sum counts per sample)")

rna_counts <- GetAssayData(data, assay = "RNA", layer = "counts")[, rownames(meta)]

groups <- split(rownames(meta), meta$group_id)
bulk_counts <- sapply(groups, function(cells) {
  Matrix::rowSums(rna_counts[, cells, drop = FALSE])
})
bulk_counts <- as.matrix(bulk_counts)

# Sample-level metadata
sample_info <- data.frame(
  sample = colnames(bulk_counts), stringsAsFactors = FALSE
)
rownames(sample_info) <- sample_info$sample
sample_info$age <- sapply(sample_info$sample, function(x) {
  unique(meta$age[meta$sample == x])[1]
})
sample_info$sex <- sapply(sample_info$sample, function(x) {
  unique(meta$sex[meta$sample == x])[1]
})
sample_info$age <- factor(sample_info$age, levels = AGE_LEVELS)
sample_info$sex <- factor(sample_info$sex)

message(sprintf("  [OK] Pseudobulk: %d genes x %d samples",
                nrow(bulk_counts), ncol(bulk_counts)))


# ==============================================================================
# HELPER: Run a single DESeq2 LRT model + all downstream analyses
# ==============================================================================
run_deseq_model <- function(model_name, full_design, reduced_design,
                            use_sex_in_valid = TRUE) {

  out_dir <- file.path(BASE_OUT_DIR, model_name)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  banner(sprintf("MODEL: %s", model_name))
  message(sprintf("  Full design:    %s", deparse(full_design)))
  message(sprintf("  Reduced design: %s", deparse(reduced_design)))

  # Valid samples
  valid <- !is.na(sample_info$age)
  if (use_sex_in_valid) valid <- valid & !is.na(sample_info$sex)

  if (sum(valid) < 2) {
    message("  [SKIP] Fewer than 2 valid samples")
    return(NULL)
  }

  dds <- DESeqDataSetFromMatrix(
    countData = round(bulk_counts[, valid]),
    colData   = sample_info[valid, , drop = FALSE],
    design    = full_design
  )

  # Row filter
  keep <- rowSums(counts(dds) >= MIN_COUNT) >= MIN_SAMPLES
  dds <- dds[keep, ]
  message(sprintf("  Genes after filtering (>= %d in >= %d samples): %d",
                  MIN_COUNT, MIN_SAMPLES, nrow(dds)))

  # LRT
  dds <- DESeq(dds, test = "LRT", reduced = reduced_design)

  # ---------- LRT results ----------
  lrt_res <- results(dds)
  lrt_df  <- as.data.frame(lrt_res)
  write.csv(lrt_df, file.path(out_dir, "deseq2_LRT_age_results.csv"))

  lrt_sig <- subset(lrt_df, padj < 0.05)
  lrt_sig <- lrt_sig[order(lrt_sig$padj), ]
  write.csv(lrt_sig, file.path(out_dir, "deseq2_LRT_age_significant.csv"))
  message(sprintf("  [OK] Significant genes (padj < 0.05): %d", nrow(lrt_sig)))

  # ---------- Pairwise age contrasts ----------
  message("  Running pairwise age contrasts...")
  contrast_results <- list()
  for (i in 1:(length(AGE_LEVELS) - 1)) {
    for (j in (i + 1):length(AGE_LEVELS)) {
      level1 <- AGE_LEVELS[j]
      level2 <- AGE_LEVELS[i]
      res_pair <- tryCatch(
        results(dds, contrast = c("age", level1, level2)),
        error = function(e) NULL
      )
      if (is.null(res_pair)) next
      cname <- paste0("age_", level1, "_vs_", level2)
      contrast_results[[cname]] <- as.data.frame(res_pair)
      write.csv(as.data.frame(res_pair),
                file.path(out_dir, paste0("deseq2_", cname, ".csv")))
    }
  }
  saveRDS(contrast_results,
          file.path(out_dir, "deseq2_all_pairwise_contrasts.rds"))

  # ---------- VST + PCA ----------
  vsd <- vst(dds, blind = FALSE)

  pca_df <- plotPCA(vsd, intgroup = c("age", "sex"), returnData = TRUE)
  percentVar <- round(100 * attr(pca_df, "percentVar"))
  p_pca <- ggplot(pca_df, aes(PC1, PC2, color = age, shape = sex)) +
    geom_point(size = 4) +
    scale_color_manual(values = AGE_COLORS) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    theme_classic(base_size = 14) +
    theme(text = element_text(family = "Arial")) +
    labs(title = sprintf("PCA - %s (%s)", TARGET_CELLTYPE, model_name))

  ggsave(file.path(out_dir, "PCA_age_sex.pdf"),
         p_pca, width = 7, height = 5, useDingbats = FALSE)
  message("  [OK] PCA_age_sex.pdf")

  # ---------- Heatmap + k-means clustering ----------
  if (nrow(lrt_sig) == 0) {
    message("  [SKIP] No significant genes -> no heatmap")
    return(list(dds = dds, lrt_df = lrt_df, lrt_sig = lrt_sig))
  }

  all_sig_genes <- rownames(lrt_sig)
  mat_all <- assay(vsd)[all_sig_genes, , drop = FALSE]
  mat_all <- t(scale(t(mat_all)))

  # Column annotation + ordering: male first (by age), then female (by age)
  anno_col <- sample_info[colnames(mat_all), c("age", "sex")]
  anno_col$age <- factor(anno_col$age, levels = AGE_LEVELS)
  anno_col$sex <- factor(anno_col$sex, levels = c("male", "female"))
  annotation_colors <- list(age = AGE_COLORS, sex = SEX_COLORS)

  male_samples   <- rownames(anno_col)[anno_col$sex == "male"]
  female_samples <- rownames(anno_col)[anno_col$sex == "female"]
  male_samples   <- male_samples[order(anno_col[male_samples, "age"])]
  female_samples <- female_samples[order(anno_col[female_samples, "age"])]

  expr_combined <- cbind(
    mat_all[, male_samples, drop = FALSE],
    mat_all[, female_samples, drop = FALSE]
  )
  anno_combined <- rbind(
    anno_col[male_samples, , drop = FALSE],
    anno_col[female_samples, , drop = FALSE]
  )

  # K-means (seeded) on genes
  set.seed(SEED)
  km_res <- kmeans(expr_combined, centers = K_CLUSTERS, nstart = 25)
  gene_clusters <- km_res$cluster
  message(sprintf("  K-means cluster sizes: %s",
                  paste(table(gene_clusters), collapse = ", ")))

  # Export gene-cluster assignments
  gene_cluster_table <- data.frame(
    Gene = names(gene_clusters), Cluster = gene_clusters,
    stringsAsFactors = FALSE
  )
  write.csv(gene_cluster_table,
            file.path(out_dir, "genes_cluster_assignments.csv"),
            row.names = FALSE)

  cluster_gene_list <- split(names(gene_clusters), gene_clusters)
  for (cl in names(cluster_gene_list)) {
    write.csv(
      data.frame(Gene = cluster_gene_list[[cl]]),
      file.path(out_dir, paste0("genes_cluster_", cl, ".csv")),
      row.names = FALSE, quote = FALSE
    )
    message(sprintf("    Cluster %s: %d genes",
                    cl, length(cluster_gene_list[[cl]])))
  }

  # Gene stats + cluster assignment combined
  gene_cluster_full <- data.frame(
    Gene = names(gene_clusters), Cluster = gene_clusters,
    stringsAsFactors = FALSE
  ) %>%
    left_join(
      lrt_sig %>% rownames_to_column("Gene") %>%
        select(Gene, baseMean, log2FoldChange, pvalue, padj),
      by = "Gene"
    )
  write.csv(gene_cluster_full,
            file.path(out_dir, "genes_clusters_with_stats.csv"),
            row.names = FALSE)

  # Z-score color ramp (1st-99th percentile)
  zlim <- quantile(expr_combined, c(0.01, 0.99), na.rm = TRUE)
  col_fun_z <- colorRamp2(
    seq(zlim[1], zlim[2], length.out = 11),
    rev(brewer.pal(11, "RdYlBu"))
  )

  ht <- Heatmap(
    expr_combined,
    name = "z-score",
    col  = col_fun_z,
    show_row_names    = FALSE,
    show_column_names = FALSE,
    top_annotation    = HeatmapAnnotation(df = anno_combined,
                                          col = annotation_colors),
    cluster_rows      = TRUE,
    cluster_columns   = FALSE,
    row_split         = gene_clusters
  )
  pdf(file.path(out_dir, "heatmap_zscore.pdf"),
      width = 5, height = 10, useDingbats = FALSE)
  draw(ht, heatmap_legend_side = "right",
       annotation_legend_side = "right")
  dev.off()
  message("  [OK] heatmap_zscore.pdf")

  # Line plot of cluster mean expression by sex x age
  mat_long <- as.data.frame(expr_combined)
  mat_long$Cluster <- gene_clusters
  mat_long$Gene    <- rownames(expr_combined)
  mat_long_tidy <- mat_long %>%
    pivot_longer(-c(Cluster, Gene), names_to = "Sample", values_to = "Expression") %>%
    left_join(anno_combined %>% rownames_to_column("Sample"), by = "Sample")

  cluster_pattern_sex <- mat_long_tidy %>%
    group_by(Cluster, age, sex) %>%
    summarise(
      mean_expr = mean(Expression, na.rm = TRUE),
      sd_expr   = sd(Expression, na.rm = TRUE),
      n         = n(),
      se        = sd_expr / sqrt(n),
      .groups = "drop"
    ) %>%
    mutate(age = factor(age, levels = AGE_LEVELS))

  p_sex <- ggplot(
    cluster_pattern_sex,
    aes(x = age, y = mean_expr, color = sex, group = sex)
  ) +
    geom_ribbon(aes(ymin = mean_expr - se, ymax = mean_expr + se, fill = sex),
                alpha = 0.2, color = NA) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 2.5) +
    scale_color_manual(values = SEX_COLORS) +
    scale_fill_manual(values = SEX_COLORS) +
    facet_wrap(~ Cluster, scales = "free_y", ncol = 2) +
    theme_minimal(base_size = 14) +
    theme(
      text            = element_text(family = "Arial"),
      axis.text.x     = element_text(angle = 45, hjust = 1),
      strip.text      = element_text(face = "bold", size = 12),
      legend.position = "bottom"
    ) +
    labs(
      title = sprintf("Gene cluster patterns - %s (%s)",
                      TARGET_CELLTYPE, model_name),
      x     = "Age group", y = "Mean z-score",
      color = "Sex", fill = "Sex"
    )
  ggsave(file.path(out_dir, "lineplot_clusters_by_sex.pdf"),
         p_sex, width = 10, height = 8, useDingbats = FALSE)
  message("  [OK] lineplot_clusters_by_sex.pdf")

  return(list(dds = dds, lrt_df = lrt_df, lrt_sig = lrt_sig))
}


# ==============================================================================
# STEP 3: RUN BOTH MODELS
# ==============================================================================
res_with_sex <- run_deseq_model(
  model_name       = "model_with_sex",
  full_design      = ~ sex + age,
  reduced_design   = ~ sex,
  use_sex_in_valid = TRUE
)

res_without_sex <- run_deseq_model(
  model_name       = "model_without_sex",
  full_design      = ~ age,
  reduced_design   = ~ 1,
  use_sex_in_valid = FALSE
)


# ==============================================================================
# STEP 4: COMPARE MODELS
# ==============================================================================
banner("STEP 4: Compare the two models")

sig_with    <- if (!is.null(res_with_sex))    rownames(res_with_sex$lrt_sig)    else character(0)
sig_without <- if (!is.null(res_without_sex)) rownames(res_without_sex$lrt_sig) else character(0)

comparison <- data.frame(
  Model = c("with_sex (~ sex + age)", "without_sex (~ age)"),
  N_sig_padj_0.05 = c(length(sig_with), length(sig_without)),
  Shared = length(intersect(sig_with, sig_without)),
  Unique = c(
    length(setdiff(sig_with, sig_without)),
    length(setdiff(sig_without, sig_with))
  )
)
print(comparison)
write.csv(
  comparison,
  file.path(BASE_OUT_DIR, "model_comparison_summary.csv"),
  row.names = FALSE
)

# Gene-level overlap
all_genes <- union(sig_with, sig_without)
overlap_df <- data.frame(
  Gene            = all_genes,
  Sig_with_sex    = all_genes %in% sig_with,
  Sig_without_sex = all_genes %in% sig_without
)
write.csv(overlap_df,
          file.path(BASE_OUT_DIR, "sig_gene_overlap.csv"),
          row.names = FALSE)

message(sprintf("\n  Shared sig genes:          %d",
                length(intersect(sig_with, sig_without))))
message(sprintf("  Only in with_sex model:    %d",
                length(setdiff(sig_with, sig_without))))
message(sprintf("  Only in without_sex model: %d",
                length(setdiff(sig_without, sig_with))))

# Venn of the two significant-gene sets
venn_out <- file.path(BASE_OUT_DIR, "venn_sig_overlap.png")
if (length(sig_with) > 0 && length(sig_without) > 0) {
  tryCatch({
    venn.diagram(
      x = list(with_sex = sig_with, without_sex = sig_without),
      filename  = venn_out,
      fill      = c("#E74C3C", "#3498DB"),
      alpha     = 0.5,
      cex       = 1.5,
      cat.cex   = 1.2,
      cat.pos   = c(-60, 60),
      cat.dist  = c(0.08, 0.08),
      margin    = 0.15,
      scaled    = FALSE,
      height    = 2000,
      width     = 2000,
      resolution = 400
    )
    message(sprintf("  [OK] %s", venn_out))
  }, error = function(e) {
    message(sprintf("  [WARN] VennDiagram failed: %s", e$message))
  })
} else {
  message("  [SKIP] Venn diagram (one or both sets are empty)")
}


# ==============================================================================
# DONE
# ==============================================================================
banner("PIPELINE COMPLETE")
message(sprintf("  Base output directory: %s/", BASE_OUT_DIR))
message("    model_with_sex/       (all outputs for ~ sex + age)")
message("    model_without_sex/    (all outputs for ~ age)")
message("    model_comparison_summary.csv")
message("    sig_gene_overlap.csv")
message("    venn_sig_overlap.png")

sessionInfo()
