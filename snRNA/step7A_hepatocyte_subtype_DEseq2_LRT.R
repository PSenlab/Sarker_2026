#!/usr/bin/env Rscript
# ==============================================================
# hepatocyte_complete_R_pipeline.R
# ==============================================================
# Complete R pipeline for hepatocyte sub-cluster analysis:
#
#   STEP 0: Load h5ad → Seurat, subset to Hepatocytes
#   STEP 1: Pseudobulk aggregation + DESeq2 LRT per celltype2
#   STEP 2: Extract significant DEGs (FDR < 0.05)
#   STEP 3: Z-scored heatmaps (K-means clustered, log2FC annotation)
#   STEP 4: Export cluster gene lists for Reactome enrichment
#   STEP 5: Merged Excel — all FDR < 0.05 DEGs (one tab per celltype)
#   STEP 6: Merged Excel — heatmap genes stat > 40 with cluster column
#
# Outputs:
#   - bulk_counts_<ct>.rds / sample_info_<ct>.rds
#   - deseq2_LRT_age_<ct>.csv (full results)
#   - deseq2_LRT_age_<ct>_sig_FDR_lt_0.05.csv
#   - plotMA_<ct>.pdf
#   - <ct>_heatmap_stat<threshold>.pdf / .png
#   - <ct>_cluster_<k>_genes_cluster.csv
#   - <ct>_all_cluster_assignments.csv
#   - DESeq2_LRT_all_celltypes_FDR005.xlsx
#   - Heatmap_genes_stat40_with_clusters.xlsx
# ==============================================================

# ============================================================
# Libraries
# ============================================================
library(Seurat)
library(Matrix)
library(DESeq2)
library(dplyr)
library(pheatmap)
library(circlize)
library(ComplexHeatmap)
library(RColorBrewer)
library(gridExtra)
library(grid)
library(ggplot2)
library(openxlsx)

# ============================================================
# Configuration
# ============================================================
H5AD_PATH       <- "adata_with_resolutions2.h5ad"
CELLTYPES       <- c("Hep-01", "Hep-02", "Hep-03", "Hep-04",
                      "Hep-05", "Hep-06", "Hep-07")
AGE_LEVELS      <- c("young", "mid_age", "old", "pre_geriatric", "geriatric")
STAT_THRESHOLD  <- 40
KMEANS_K        <- 2
LFC_CAP         <- 7
MIN_GENES       <- 5
FDR_CUTOFF      <- 0.05

# Color palettes
AGE_COLORS <- c(
  "young"          = "lightblue",
  "mid_age"        = "mediumseagreen",
  "old"            = "goldenrod",
  "pre_geriatric"  = "lightpink",
  "geriatric"      = "red"
)
SEX_COLORS <- c("male" = "navy", "female" = "deeppink")

CLUSTER_COLORS <- c(
  "1" = "forestgreen", "2" = "firebrick",
  "3" = "steelblue",   "4" = "darkorange",
  "5" = "purple",       "6" = "gold"
)


# ============================================================
# STEP 0: Load h5ad → Seurat & subset to Hepatocytes
# ============================================================
message("\n", paste(rep("=", 60), collapse = ""))
message("STEP 0: Loading data & subsetting to Hepatocytes")
message(paste(rep("=", 60), collapse = ""))

seuratObj <- schard::h5ad2seurat(H5AD_PATH)
message("  Full dataset: ", ncol(seuratObj), " cells x ", nrow(seuratObj), " genes")

data <- subset(seuratObj, subset = celltype == "Hepatocyte")
message("  Hepatocytes:  ", ncol(data), " cells")
message("  Sub-clusters: ", paste(sort(unique(data$celltype2)), collapse = ", "))
message("  Sex:          ", paste(names(table(data$sex)), table(data$sex),
                                   sep = "=", collapse = ", "))


# ============================================================
# STEP 1: Pseudobulk aggregation + DESeq2 LRT (per celltype2)
# ============================================================
message("\n", paste(rep("=", 60), collapse = ""))
message("STEP 1: Pseudobulk aggregation + DESeq2 LRT")
message(paste(rep("=", 60), collapse = ""))

for (ct in CELLTYPES) {

  message("\n--- Processing celltype2: ", ct, " ---")

  # 1a. Subset
  data_ct <- tryCatch(
    subset(data, subset = celltype2 == ct),
    error = function(e) { message("  Skipping ", ct, ": subset failed"); NULL }
  )
  if (is.null(data_ct) || ncol(data_ct) == 0) {
    message("  Skipping ", ct, ": no cells found")
    next
  }

  meta_ct <- data_ct@meta.data
  if (length(unique(meta_ct$sample)) < 2) {
    message("  Skipping ", ct, ": not enough samples (need >= 2)")
    next
  }

  meta_ct$age <- factor(meta_ct$age, levels = AGE_LEVELS)

  # 1b. Aggregate raw counts per sample (pseudobulk)
  rna_counts <- GetAssayData(data_ct, assay = "RNA", slot = "counts")
  groups <- split(rownames(meta_ct), meta_ct$sample)
  bulk_counts <- sapply(groups, function(cells) {
    Matrix::rowSums(rna_counts[, cells, drop = FALSE])
  })
  bulk_counts <- as.matrix(bulk_counts)

  # 1c. Build sample-level metadata
  sample_info <- data.frame(sample = colnames(bulk_counts), stringsAsFactors = FALSE)
  rownames(sample_info) <- sample_info$sample
  sample_info$age <- sapply(sample_info$sample, function(x) {
    unique(meta_ct$age[meta_ct$sample == x])
  })
  sample_info$sex <- sapply(sample_info$sample, function(x) {
    unique(meta_ct$sex[meta_ct$sample == x])
  })
  sample_info$age <- factor(sample_info$age, levels = AGE_LEVELS)
  sample_info$sex <- factor(sample_info$sex)

  valid <- !is.na(sample_info$age) & !is.na(sample_info$sex)
  if (sum(valid) < 2) {
    message("  Skipping ", ct, ": not enough valid samples after filtering")
    next
  }

  # 1d. Save pseudobulk intermediates
  saveRDS(bulk_counts, paste0("bulk_counts_", ct, ".rds"))
  saveRDS(sample_info, paste0("sample_info_", ct, ".rds"))
  message("  Saved: bulk_counts_", ct, ".rds + sample_info_", ct, ".rds")

  # 1e. DESeq2 LRT: full ~sex+age vs reduced ~sex
  dds <- DESeqDataSetFromMatrix(
    countData = round(bulk_counts[, valid]),
    colData   = sample_info[valid, ],
    design    = ~ sex + age
  )
  dds <- DESeq(dds, test = "LRT", reduced = ~ sex)
  res <- results(dds)

  out_file <- paste0("deseq2_LRT_age_", ct, ".csv")
  write.csv(as.data.frame(res), file = out_file)
  message("  Saved: ", out_file)

  # MA plot
  pdf(paste0("plotMA_", ct, ".pdf"))
  plotMA(res, ylim = c(-5, 5), main = ct)
  dev.off()
  message("  Saved: plotMA_", ct, ".pdf")
}


# ============================================================
# STEP 2: Extract significant DEGs (FDR < 0.05)
# ============================================================
message("\n", paste(rep("=", 60), collapse = ""))
message("STEP 2: Extracting significant DEGs (FDR < ", FDR_CUTOFF, ")")
message(paste(rep("=", 60), collapse = ""))

result_files <- list.files(pattern = "^deseq2_LRT_age_Hep-0[1-7]\\.csv$")

for (file in result_files) {
  message("\n  Processing: ", file)
  res_df <- read.csv(file, row.names = 1)

  if (!"padj" %in% colnames(res_df)) {
    warning("  Skipping ", file, ": 'padj' column missing")
    next
  }

  res_sig <- res_df[which(res_df$padj < FDR_CUTOFF), ]
  res_sig <- res_sig[order(res_sig$padj), ]

  if (nrow(res_sig) == 0) {
    message("  No significant genes in ", file)
    next
  }

  output_file <- sub("\\.csv$", "_sig_FDR_lt_0.05.csv", file)
  write.csv(res_sig, output_file, quote = FALSE)
  message("  Saved: ", output_file, " (", nrow(res_sig), " genes)")
}


# ============================================================
# STEP 3: Heatmaps — z-scored log1p pseudobulk, K-means
# ============================================================
message("\n", paste(rep("=", 60), collapse = ""))
message("STEP 3: Generating heatmaps (z-scored log1p expression)")
message(paste(rep("=", 60), collapse = ""))

for (ct in CELLTYPES) {

  message("\n--- Heatmap for: ", ct, " ---")

  # 3a. Load significant DE results
  res_file <- paste0("deseq2_LRT_age_", ct, "_sig_FDR_lt_0.05.csv")
  if (!file.exists(res_file)) {
    message("  Skipping: ", res_file, " not found")
    next
  }
  res_sig <- read.csv(res_file, row.names = 1)

  strong_sig_genes <- rownames(res_sig[res_sig$stat > STAT_THRESHOLD, ])
  if (length(strong_sig_genes) < MIN_GENES) {
    message("  Skipping ", ct, ": only ", length(strong_sig_genes),
            " genes pass stat > ", STAT_THRESHOLD)
    next
  }

  # 3b. Load pseudobulk
  bulk_file <- paste0("bulk_counts_", ct, ".rds")
  info_file <- paste0("sample_info_", ct, ".rds")
  if (!file.exists(bulk_file) || !file.exists(info_file)) {
    message("  Skipping: missing .rds files for ", ct)
    next
  }
  bulk_counts <- readRDS(bulk_file)
  sample_info <- readRDS(info_file)

  # 3c. Z-scored expression matrix
  strong_sig_genes <- intersect(strong_sig_genes, rownames(bulk_counts))
  if (length(strong_sig_genes) < MIN_GENES) {
    message("  Skipping: not enough genes after intersect")
    next
  }

  expr_mat <- bulk_counts[strong_sig_genes, ]
  expr_scaled <- t(scale(t(log1p(expr_mat))))

  message("  Matrix: ", nrow(expr_scaled), " genes x ", ncol(expr_scaled), " samples")

  # 3d. Column annotation
  anno_col <- sample_info[, c("age", "sex")]
  anno_col$age <- factor(anno_col$age, levels = AGE_LEVELS)
  anno_col$sex <- factor(anno_col$sex, levels = c("male", "female"))

  annotation_colors <- list(age = AGE_COLORS, sex = SEX_COLORS)

  # 3e. Order columns: male (by age) | female (by age)
  male_samples   <- rownames(anno_col)[anno_col$sex == "male"]
  female_samples <- rownames(anno_col)[anno_col$sex == "female"]
  male_samples   <- male_samples[order(anno_col[male_samples, "age"])]
  female_samples <- female_samples[order(anno_col[female_samples, "age"])]

  expr_combined <- cbind(expr_scaled[, male_samples, drop = FALSE],
                         expr_scaled[, female_samples, drop = FALSE])
  anno_combined <- rbind(anno_col[male_samples, , drop = FALSE],
                         anno_col[female_samples, , drop = FALSE])

  # 3f. K-means row clustering
  heatmap_km <- pheatmap::pheatmap(
    expr_combined,
    annotation_col           = anno_combined,
    annotation_colors        = annotation_colors,
    show_colnames            = FALSE,
    show_rownames            = TRUE,
    cluster_rows             = TRUE,
    cluster_cols             = FALSE,
    kmeans_k                 = KMEANS_K,
    clustering_distance_rows = "correlation",
    gaps_col                 = length(male_samples),
    fontsize_row             = 5,
    main                     = paste0(ct, ": K-means pre-clustering"),
    silent                   = TRUE
  )

  gene_clusters  <- heatmap_km$kmeans$cluster
  ordered_genes  <- names(sort(gene_clusters))
  expr_reordered <- expr_combined[ordered_genes, ]

  # 3g. Row annotation: cluster + capped log2FC
  lfc_raw    <- res_sig[ordered_genes, "log2FoldChange"]
  lfc_capped <- pmax(pmin(lfc_raw, LFC_CAP), -LFC_CAP)

  anno_row <- data.frame(
    cluster = factor(gene_clusters[ordered_genes]),
    log2FC  = lfc_capped
  )
  rownames(anno_row) <- ordered_genes

  cluster_cols <- CLUSTER_COLORS[as.character(seq_len(KMEANS_K))]
  names(cluster_cols) <- as.character(seq_len(KMEANS_K))

  col_fun_lfc      <- colorRamp2(c(-LFC_CAP, 0, LFC_CAP),
                                  c("darkblue", "white", "darkred"))
  lfc_strip_colors <- colorRampPalette(c("darkblue", "white", "darkred"))(100)

  annotation_colors_all <- c(
    annotation_colors,
    list(cluster = cluster_cols),
    list(log2FC  = lfc_strip_colors)
  )

  # 3h. Final heatmap
  heatmap_final <- pheatmap::pheatmap(
    expr_reordered,
    annotation_col    = anno_combined,
    annotation_row    = anno_row[, c("cluster", "log2FC")],
    annotation_colors = annotation_colors_all,
    show_colnames     = FALSE,
    show_rownames     = TRUE,
    cluster_rows      = FALSE,
    cluster_cols      = FALSE,
    gaps_col          = length(male_samples),
    fontsize_row      = 5,
    fontsize_col      = 10,
    main              = paste0(ct, ": Male vs Female across Age\n",
                               "(K-means=", KMEANS_K,
                               ", LRT stat > ", STAT_THRESHOLD,
                               ", FDR < ", FDR_CUTOFF, ")"),
    silent            = TRUE
  )

  # 3i. External log2FC legend
  log2fc_legend <- Legend(
    title          = "log2FC",
    col_fun        = col_fun_lfc,
    at             = c(-LFC_CAP, -LFC_CAP / 2, 0, LFC_CAP / 2, LFC_CAP),
    direction      = "vertical",
    title_position = "topcenter",
    legend_height  = unit(4, "cm")
  )

  legend_grob <- grid.grabExpr(draw(log2fc_legend))
  final_plot  <- grid.arrange(heatmap_final$gtable, legend_grob,
                              ncol = 2, widths = c(15, 2))

  # 3j. Save
  safe_name <- tolower(gsub("-", "_", ct))
  pdf_out <- paste0(safe_name, "_heatmap_stat", STAT_THRESHOLD, ".pdf")
  png_out <- paste0(safe_name, "_heatmap_stat", STAT_THRESHOLD, ".png")

  ggsave(pdf_out, final_plot, width = 10, height = 15, dpi = 300)
  ggsave(png_out, final_plot, width = 10, height = 15, dpi = 300)
  message("  Saved: ", pdf_out, " + ", png_out)
}


# ============================================================
# STEP 4: Export cluster gene lists (CSVs) for Reactome
# ============================================================
message("\n", paste(rep("=", 60), collapse = ""))
message("STEP 4: Exporting cluster gene CSVs for Reactome")
message(paste(rep("=", 60), collapse = ""))

for (ct in CELLTYPES) {

  message("\n--- Exporting clusters for: ", ct, " ---")

  res_file  <- paste0("deseq2_LRT_age_", ct, "_sig_FDR_lt_0.05.csv")
  bulk_file <- paste0("bulk_counts_", ct, ".rds")
  info_file <- paste0("sample_info_", ct, ".rds")

  if (!file.exists(res_file) || !file.exists(bulk_file) || !file.exists(info_file)) {
    message("  Skipping: required files not found")
    next
  }

  res_sig     <- read.csv(res_file, row.names = 1)
  bulk_counts <- readRDS(bulk_file)
  sample_info <- readRDS(info_file)

  strong_genes <- rownames(res_sig[res_sig$stat > STAT_THRESHOLD, ])
  strong_genes <- intersect(strong_genes, rownames(bulk_counts))

  if (length(strong_genes) < MIN_GENES) {
    message("  Skipping: only ", length(strong_genes), " genes pass threshold")
    next
  }

  expr_mat    <- bulk_counts[strong_genes, ]
  expr_scaled <- t(scale(t(log1p(expr_mat))))

  anno_col <- sample_info[, c("age", "sex")]
  anno_col$age <- factor(anno_col$age, levels = AGE_LEVELS)
  anno_col$sex <- factor(anno_col$sex, levels = c("male", "female"))

  male_samples   <- rownames(anno_col)[anno_col$sex == "male"]
  female_samples <- rownames(anno_col)[anno_col$sex == "female"]
  male_samples   <- male_samples[order(anno_col[male_samples, "age"])]
  female_samples <- female_samples[order(anno_col[female_samples, "age"])]
  expr_combined  <- cbind(expr_scaled[, male_samples, drop = FALSE],
                          expr_scaled[, female_samples, drop = FALSE])

  hm <- pheatmap::pheatmap(
    expr_combined,
    cluster_rows = TRUE, cluster_cols = FALSE,
    kmeans_k = KMEANS_K, clustering_distance_rows = "correlation",
    silent = TRUE
  )
  gene_clusters <- hm$kmeans$cluster
  ct_safe <- gsub("-", "_", ct)

  for (k in seq_len(KMEANS_K)) {
    genes_k <- names(gene_clusters[gene_clusters == k])
    if (length(genes_k) == 0) next
    gene_info <- data.frame(
      gene = genes_k, log2FC = res_sig[genes_k, "log2FoldChange"],
      stat = res_sig[genes_k, "stat"], padj = res_sig[genes_k, "padj"],
      cluster = k, stringsAsFactors = FALSE
    )
    gene_info <- gene_info[order(-abs(gene_info$log2FC)), ]
    out_file <- paste0(ct_safe, "_cluster_", k, "_genes_cluster.csv")
    write.csv(gene_info, out_file, row.names = FALSE)
    message("  Cluster ", k, ": ", length(genes_k), " genes -> ", out_file)
  }

  all_clusters <- data.frame(
    gene = names(gene_clusters), cluster = gene_clusters,
    log2FC = res_sig[names(gene_clusters), "log2FoldChange"],
    stat = res_sig[names(gene_clusters), "stat"],
    padj = res_sig[names(gene_clusters), "padj"],
    stringsAsFactors = FALSE
  )
  combined_out <- paste0(ct_safe, "_all_cluster_assignments.csv")
  write.csv(all_clusters, combined_out, row.names = FALSE)
  message("  Combined: ", combined_out)
}


# ============================================================
# STEP 5: Merged Excel — all FDR < 0.05 DEGs per celltype
# ============================================================
message("\n", paste(rep("=", 60), collapse = ""))
message("STEP 5: Merged Excel — all FDR < 0.05 DEGs")
message(paste(rep("=", 60), collapse = ""))

wb1 <- createWorkbook()

for (ct in CELLTYPES) {

  res_file <- paste0("deseq2_LRT_age_", ct, "_sig_FDR_lt_0.05.csv")
  if (!file.exists(res_file)) {
    message("  Skipping ", ct, ": ", res_file, " not found")
    next
  }

  res_sig <- read.csv(res_file, row.names = 1)

  # Gene column from rownames, reorder gene first
  res_sig$gene <- rownames(res_sig)
  col_order <- c("gene", setdiff(colnames(res_sig), "gene"))
  res_sig <- res_sig[, col_order]
  res_sig <- res_sig[order(res_sig$padj), ]

  # Write tab
  addWorksheet(wb1, ct)
  writeData(wb1, ct, res_sig, rowNames = FALSE)

  # Format: auto-width + bold blue header
  setColWidths(wb1, ct, cols = seq_along(col_order), widths = "auto")
  headerStyle <- createStyle(textDecoration = "Bold", border = "Bottom",
                              fgFill = "#D9E1F2")
  addStyle(wb1, ct, headerStyle, rows = 1, cols = seq_along(col_order))

  # Highlight padj < 0.001 rows in yellow
  if (any(res_sig$padj < 0.001, na.rm = TRUE)) {
    highlight_rows <- which(res_sig$padj < 0.001) + 1
    highlightStyle <- createStyle(fgFill = "#FFF2CC")
    addStyle(wb1, ct, highlightStyle, rows = highlight_rows,
             cols = seq_along(col_order), gridExpand = TRUE, stack = TRUE)
  }

  message("  Tab: ", ct, " — ", nrow(res_sig), " genes")
}

out1 <- "DESeq2_LRT_all_celltypes_FDR005.xlsx"
saveWorkbook(wb1, out1, overwrite = TRUE)
message("  Saved: ", out1)


# ============================================================
# STEP 6: Merged Excel — heatmap genes (stat > 40) + cluster
# ============================================================
message("\n", paste(rep("=", 60), collapse = ""))
message("STEP 6: Merged Excel — heatmap genes (stat > ", STAT_THRESHOLD, ") + clusters")
message(paste(rep("=", 60), collapse = ""))

wb2 <- createWorkbook()

# Cluster row fill colors
cluster_fill <- c("1" = "#E8F5E9", "2" = "#FFEBEE", "3" = "#E3F2FD",
                   "4" = "#FFF3E0", "5" = "#F3E5F5", "6" = "#FFFDE7")

for (ct in CELLTYPES) {

  message("\n--- ", ct, " ---")

  res_file  <- paste0("deseq2_LRT_age_", ct, "_sig_FDR_lt_0.05.csv")
  bulk_file <- paste0("bulk_counts_", ct, ".rds")
  info_file <- paste0("sample_info_", ct, ".rds")

  if (!file.exists(res_file) || !file.exists(bulk_file) || !file.exists(info_file)) {
    message("  Skipping: required files not found")
    next
  }

  res_sig     <- read.csv(res_file, row.names = 1)
  bulk_counts <- readRDS(bulk_file)
  sample_info <- readRDS(info_file)

  # Filter stat > threshold
  strong_genes <- rownames(res_sig[res_sig$stat > STAT_THRESHOLD, ])
  strong_genes <- intersect(strong_genes, rownames(bulk_counts))

  if (length(strong_genes) < MIN_GENES) {
    message("  Skipping: only ", length(strong_genes), " genes pass stat > ", STAT_THRESHOLD)
    next
  }

  # Z-score (same as heatmap)
  expr_mat    <- bulk_counts[strong_genes, ]
  expr_scaled <- t(scale(t(log1p(expr_mat))))

  # Order samples
  anno_col <- sample_info[, c("age", "sex")]
  anno_col$age <- factor(anno_col$age, levels = AGE_LEVELS)
  anno_col$sex <- factor(anno_col$sex, levels = c("male", "female"))

  male_samples   <- rownames(anno_col)[anno_col$sex == "male"]
  female_samples <- rownames(anno_col)[anno_col$sex == "female"]
  male_samples   <- male_samples[order(anno_col[male_samples, "age"])]
  female_samples <- female_samples[order(anno_col[female_samples, "age"])]
  expr_combined  <- cbind(expr_scaled[, male_samples, drop = FALSE],
                          expr_scaled[, female_samples, drop = FALSE])

  # K-means (same as heatmap)
  hm <- pheatmap::pheatmap(
    expr_combined,
    cluster_rows = TRUE, cluster_cols = FALSE,
    kmeans_k = KMEANS_K, clustering_distance_rows = "correlation",
    silent = TRUE
  )
  gene_clusters <- hm$kmeans$cluster

  # Build output
  out_df <- data.frame(
    gene           = names(gene_clusters),
    baseMean       = res_sig[names(gene_clusters), "baseMean"],
    log2FoldChange = res_sig[names(gene_clusters), "log2FoldChange"],
    lfcSE          = res_sig[names(gene_clusters), "lfcSE"],
    stat           = res_sig[names(gene_clusters), "stat"],
    pvalue         = res_sig[names(gene_clusters), "pvalue"],
    padj           = res_sig[names(gene_clusters), "padj"],
    cluster        = gene_clusters,
    stringsAsFactors = FALSE
  )
  out_df <- out_df[order(out_df$cluster, -out_df$stat), ]

  # Write tab
  addWorksheet(wb2, ct)
  writeData(wb2, ct, out_df, rowNames = FALSE)

  # Format: auto-width + bold blue header
  setColWidths(wb2, ct, cols = seq_along(out_df), widths = "auto")
  headerStyle <- createStyle(textDecoration = "Bold", border = "Bottom",
                              fgFill = "#D9E1F2")
  addStyle(wb2, ct, headerStyle, rows = 1, cols = seq_along(out_df))

  # Color-code entire row by cluster
  for (k in seq_len(KMEANS_K)) {
    cluster_rows <- which(out_df$cluster == k) + 1
    if (length(cluster_rows) == 0) next
    fill <- cluster_fill[as.character(k)]
    if (is.na(fill)) fill <- "#F5F5F5"
    rowStyle <- createStyle(fgFill = fill)
    addStyle(wb2, ct, rowStyle, rows = cluster_rows,
             cols = seq_along(out_df), gridExpand = TRUE, stack = TRUE)
  }

  # Summary
  cluster_summary <- table(out_df$cluster)
  summary_str <- paste(paste0("C", names(cluster_summary), "=",
                               cluster_summary), collapse = ", ")
  message("  Tab: ", ct, " — ", nrow(out_df), " genes (", summary_str, ")")
}

out2 <- "Heatmap_genes_stat40_with_clusters.xlsx"
saveWorkbook(wb2, out2, overwrite = TRUE)
message("  Saved: ", out2)


# ============================================================
# Done
# ============================================================
message("\n", paste(rep("=", 60), collapse = ""))
message("R PIPELINE COMPLETE — ALL 7 STEPS")
message(paste(rep("-", 60), collapse = ""))
message("  STEP 0: h5ad -> Seurat -> Hepatocyte subset")
message("  STEP 1: Pseudobulk + DESeq2 LRT        -> deseq2_LRT_age_*.csv")
message("  STEP 2: FDR < 0.05 filtering            -> *_sig_FDR_lt_0.05.csv")
message("  STEP 3: Z-scored heatmaps               -> *_heatmap_stat", STAT_THRESHOLD, ".pdf/.png")
message("  STEP 4: Cluster gene CSVs               -> *_cluster_*_genes_cluster.csv")
message("  STEP 5: Merged Excel (FDR < 0.05)       -> ", out1)
message("  STEP 6: Merged Excel (stat > ", STAT_THRESHOLD, " + clusters) -> ", out2)
message(paste(rep("-", 60), collapse = ""))
message("  Next: Run hepatocyte_complete_python_pipeline.py")
message(paste(rep("=", 60), collapse = ""))
