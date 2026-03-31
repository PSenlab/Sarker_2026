# ======================================================
# ENTROPY ANALYSIS - COMPLETE PIPELINE
# Part 1: Entropy Calculation
# Part 2: Ridge Plots + ANOVA Statistics
# ======================================================

# ===============================
# 0. Load Libraries
# ===============================
suppressPackageStartupMessages({
  library(Entropy)
  library(Seurat)
  library(ggplot2)
  library(ggridges)
  library(dplyr)
  library(stringr)
})

# ===============================
# 1. Configuration
# ===============================
# Input/Output paths
input_rds <- "/data/sarkern2/multiome_liver/Seurat/analyzed_data/Figure3/combined_rna_atac_seurat_ccans.rds"
output_rds <- "seurat_with_entropy_all.rds"
output_csv <- "seurat_metadata_with_entropy_all.csv"

# Output directory for plots
out_dir <- "entropy_analysis_results"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Age levels and colors
age_levels <- c("young", "mid_age", "old", "pre_geriatric", "geriatric")

age_colors <- c(
  "young"         = "#1ABC9C",
  "mid_age"       = "#F1C40F",
  "old"           = "#C39BD3",

  "pre_geriatric" = "#2980B9",
  "geriatric"     = "#E84393"
)

# Celltype order
celltype_levels <- c(
  "Hepatocyte", "Stellate",
  "Endothelial-01", "Endothelial-02",
  "Cholangiocyte-01", "Cholangiocyte-02",
  "lymp_B", "lymp_T", "Kupffer", "MoMFs"
)

min_cells <- 20

# ======================================================
# PART 1: ENTROPY CALCULATION
# ======================================================
message("\n========================================")
message("PART 1: ENTROPY CALCULATION")
message("========================================\n")

# ===============================
# 1.1 Load Seurat Object
# ===============================
message("Loading Seurat object...")
seuratObj <- readRDS(input_rds)

# Set RNA as the default assay
DefaultAssay(seuratObj) <- "RNA"
message(sprintf("Loaded: %d cells, %d features", ncol(seuratObj), nrow(seuratObj)))

# ===============================
# 1.2 Run PCA and Neighbors (REQUIRED for entropy)
# ===============================
message("Running PCA and FindNeighbors...")
seuratObj <- NormalizeData(seuratObj)
seuratObj <- FindVariableFeatures(seuratObj)
seuratObj <- ScaleData(seuratObj)
seuratObj <- RunPCA(seuratObj, npcs = 30, verbose = FALSE)
seuratObj <- FindNeighbors(seuratObj, reduction = "pca", dims = 1:30, verbose = FALSE)
message("PCA and neighbors complete")

# ===============================
# 1.3 Run Entropy Estimation
# ===============================
message("Running entropy estimation...")
seuratObj <- run_entropy(
  seu = seuratObj,
  assay = "RNA",
  nn_list = NULL,
  output_path = NULL,
  add_assay = TRUE
)
message("Entropy estimation complete")

# ===============================
# 1.4 Compute Mean Entropy Per Cell
# ===============================
message("Computing mean entropy per cell...")
entropy_means <- colMeans(GetAssayData(seuratObj[["RNA_entropy"]]))
seuratObj$entropy_score <- entropy_means

# ===============================
# 1.5 Save Outputs
# ===============================
message("Saving outputs...")
metadata_df <- seuratObj@meta.data
write.csv(metadata_df, file = output_csv, row.names = TRUE)
saveRDS(seuratObj, output_rds)
message(sprintf("Saved: %s", output_csv))
message(sprintf("Saved: %s", output_rds))

# ======================================================
# PART 2: RIDGE PLOTS + ANOVA
# ======================================================
message("\n========================================")
message("PART 2: RIDGE PLOTS + ANOVA")
message("========================================\n")

# ===============================
# 2.1 Clean Metadata
# ===============================
message("Cleaning metadata...")
df <- seuratObj@meta.data

# Standardize age
df$age <- tolower(str_trim(df$age))
df$age <- dplyr::recode(df$age,
  "midage"       = "mid_age",
  "pregeriatric" = "pre_geriatric"
)
df$age <- factor(df$age, levels = age_levels, ordered = TRUE)

# Standardize sex and celltype
df$sex <- factor(tolower(str_trim(as.character(df$sex))))
df$celltype <- str_trim(as.character(df$celltype))

# Apply celltype order
present_ct <- intersect(celltype_levels, unique(df$celltype))
extra_ct <- setdiff(unique(df$celltype), celltype_levels)
if (length(extra_ct) > 0) {
  message("Additional celltypes found: ", paste(extra_ct, collapse = ", "))
  present_ct <- c(present_ct, extra_ct)
}
df$celltype <- factor(df$celltype, levels = present_ct)

message(sprintf("Data: %d cells, %d celltypes, %d samples",
                nrow(df), length(unique(df$celltype)), length(unique(df$sample))))

# ===============================
# 2.2 Pseudobulk Aggregation for ANOVA
# ===============================
message("Aggregating to sample level for ANOVA...")
df_sample <- df %>%
  filter(!is.na(entropy_score), !is.na(age), !is.na(sex), !is.na(celltype)) %>%
  group_by(sample, age, sex, celltype) %>%
  summarise(entropy = mean(entropy_score, na.rm = TRUE), .groups = "drop")

# ===============================
# 2.3 Ridge Plot Function
# ===============================
draw_ridge_plot <- function(df, target_sex, ncol_facets = 5) {
  
  df_sex <- df %>%
    filter(sex == target_sex, !is.na(age), !is.na(entropy_score))
  
  if (nrow(df_sex) < min_cells) {
    message("Skipping ", target_sex, ": too few cells")
    return(NULL)
  }
  
  # Keep celltypes with enough cells
  keep_ct <- df_sex %>%
    count(celltype) %>%
    filter(n >= min_cells) %>%
    pull(celltype)
  
  df_sex <- df_sex %>% filter(celltype %in% keep_ct)
  df_sex$celltype <- factor(df_sex$celltype, levels = present_ct[present_ct %in% keep_ct])
  
  p <- ggplot(df_sex, aes(x = entropy_score, y = age, fill = age)) +
    geom_density_ridges(
      scale = 1.5,
      rel_min_height = 0.01,
      alpha = 0.85,
      color = "grey20",
      linewidth = 0.2
    ) +
    scale_fill_manual(values = age_colors, drop = FALSE) +
    scale_y_discrete(limits = rev(age_levels), drop = FALSE) +
    facet_wrap(~ celltype, ncol = ncol_facets) +
    theme_ridges() +
    theme(
      legend.position = "none",
      strip.text = element_text(size = 10, face = "bold"),
      axis.text.y = element_text(size = 9),
      axis.text.x = element_text(size = 8),
      plot.title = element_text(size = 14, face = "bold")
    ) +
    labs(
      title = paste0("Entropy Distribution by Age — ", tools::toTitleCase(target_sex)),
      x = "Entropy Score",
      y = "Age"
    )
  
  # Dynamic sizing
  n_ct <- length(unique(df_sex$celltype))
  n_row <- ceiling(n_ct / ncol_facets)
  
  outfile <- file.path(out_dir, paste0("Ridge_Entropy_", target_sex, ".pdf"))
  ggsave(outfile, p, width = ncol_facets * 3, height = max(4, n_row * 3.5))
  message("Saved: ", outfile)
  
  return(p)
}

# ===============================
# 2.4 Generate Ridge Plots
# ===============================
message("\nCreating ridge plots...")
p_ridge_female <- draw_ridge_plot(df, "female")
p_ridge_male <- draw_ridge_plot(df, "male")

# ===============================
# 2.5 Run ANOVA
# ===============================
message("\nRunning ANOVA...")

anova_results <- df_sample %>%
  group_by(celltype, sex) %>%
  summarise(
    n_samples = n(),
    n_ages = n_distinct(age),
    p_ANOVA = tryCatch(
      summary(aov(entropy ~ age))[[1]][["Pr(>F)"]][1],
      error = function(e) NA_real_
    ),
    .groups = "drop"
  ) %>%
  mutate(
    padj_ANOVA = p.adjust(p_ANOVA, method = "BH"),
    significance = case_when(
      is.na(padj_ANOVA) ~ "",
      padj_ANOVA < 0.001 ~ "***",
      padj_ANOVA < 0.01  ~ "**",
      padj_ANOVA < 0.05  ~ "*",
      TRUE ~ "ns"
    )
  )

# Save ANOVA results
anova_file <- file.path(out_dir, "Entropy_ANOVA_results.csv")
write.csv(anova_results, anova_file, row.names = FALSE)
message("Saved: ", anova_file)

# ======================================================
# SUMMARY
# ======================================================
message("\n========================================")
message("PIPELINE COMPLETE")
message("========================================")
message("\nOutputs:")
message(sprintf("  • %s", output_rds))
message(sprintf("  • %s", output_csv))
message(sprintf("  • %s/Ridge_Entropy_female.pdf", out_dir))
message(sprintf("  • %s/Ridge_Entropy_male.pdf", out_dir))
message(sprintf("  • %s/Entropy_ANOVA_results.csv", out_dir))

# Print ANOVA summary
message("\nANOVA Summary (significant results):")
sig_results <- anova_results %>% filter(significance != "ns" & significance != "")
if (nrow(sig_results) > 0) {
  print(sig_results)
} else {
  message("No significant results found")
}