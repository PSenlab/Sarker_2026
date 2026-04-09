#!/usr/bin/env Rscript
# ==============================================================================
# Transcriptional Entropy Analysis Across Aging
# ==============================================================================
#
# Description:
#   Computes per-cell transcriptional entropy from the RNA modality of the
#   snMultiome object using the Entropy package, then visualizes age-related
#   entropy distributions per cell type per sex (ridge plots) and tests
#   significance via per-sample pseudobulked one-way ANOVA + BH FDR.
#
# Pipeline:
#   PART 1: Entropy calculation
#     1.1  Load Seurat object (or h5ad -> Seurat)
#     1.2  PCA + neighbors (required by run_entropy)
#     1.3  Run entropy estimation
#     1.4  Compute mean entropy per cell -> entropy_score
#     1.5  Save annotated object + metadata
#
#   PART 2: Ridge plots + ANOVA
#     2.1  Clean and harmonize metadata
#     2.2  Pseudobulk aggregation per sample x celltype
#     2.3  Ridge plot per sex (one panel per cell type)
#     2.4  Per-celltype, per-sex one-way ANOVA + BH FDR
#
# Input:
#   - Annotated Seurat .rds OR AnnData .h5ad with celltype, sex, age, sample
#
# Output:
#   - seurat_with_entropy_all.rds
#   - seurat_metadata_with_entropy_all.csv
#   - Ridge_Entropy_female.pdf
#   - Ridge_Entropy_male.pdf
#   - Entropy_ANOVA_results.csv
#
# Reference:
#   Sarker et al. (2026) Cell Metabolism
#
# ==============================================================================

suppressPackageStartupMessages({
  library(Entropy)
  library(Seurat)
  library(schard)
  library(ggplot2)
  library(ggridges)
  library(dplyr)
  library(stringr)
})


# ==============================================================================
# CONFIGURATION - UPDATE THIS PATH
# ==============================================================================

# Input: either an h5ad (will be loaded via schard) or a Seurat .rds
INPUT_PATH  <- "integrated_scvi.h5ad"

# Output object + metadata
OUTPUT_RDS  <- "seurat_with_entropy_all.rds"
OUTPUT_CSV  <- "seurat_metadata_with_entropy_all.csv"

# Output directory for plots and tables
OUT_DIR     <- "entropy_analysis_results"
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

# Age levels and canonical color palette
AGE_LEVELS <- c("young", "mid_age", "old", "pre_geriatric", "geriatric")
AGE_COLORS <- c(
  "young"         = "#1ABC9C",
  "mid_age"       = "#F1C40F",
  "old"           = "#C39BD3",
  "pre_geriatric" = "#2980B9",
  "geriatric"     = "#E84393"
)

# Cell type display order
CELLTYPE_LEVELS <- c(
  "Hepatocyte", "Stellate",
  "Endothelial-01", "Endothelial-02",
  "Cholangiocyte-01", "Cholangiocyte-02",
  "lymp_B", "lymp_T", "Kupffer", "MoMFs"
)

# Age aliases for harmonization
AGE_ALIASES <- c("midage" = "mid_age", "pregeriatric" = "pre_geriatric")

MIN_CELLS <- 20


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
# PART 1: ENTROPY CALCULATION
# ==============================================================================
banner("PART 1: Entropy calculation")

# 1.1 Load object (h5ad or rds)
message("Loading input object: ", INPUT_PATH)
if (grepl("\\.h5ad$", INPUT_PATH, ignore.case = TRUE)) {
  seuratObj <- schard::h5ad2seurat(INPUT_PATH)
} else if (grepl("\\.rds$", INPUT_PATH, ignore.case = TRUE)) {
  seuratObj <- readRDS(INPUT_PATH)
} else {
  stop("INPUT_PATH must end in .h5ad or .rds")
}
DefaultAssay(seuratObj) <- "RNA"
message(sprintf("  [OK] Loaded: %d cells, %d features",
                ncol(seuratObj), nrow(seuratObj)))

# 1.2 PCA + neighbors (required by run_entropy)
message("\nRunning Normalize / VariableFeatures / Scale / PCA / FindNeighbors...")
seuratObj <- NormalizeData(seuratObj, verbose = FALSE)
seuratObj <- FindVariableFeatures(seuratObj, verbose = FALSE)
seuratObj <- ScaleData(seuratObj, verbose = FALSE)
seuratObj <- RunPCA(seuratObj, npcs = 30, verbose = FALSE)
seuratObj <- FindNeighbors(seuratObj, reduction = "pca", dims = 1:30, verbose = FALSE)
message("  [OK] PCA + neighbors complete")

# 1.3 Run Entropy estimation
message("\nRunning entropy estimation (Entropy::run_entropy)...")
seuratObj <- run_entropy(
  seu         = seuratObj,
  assay       = "RNA",
  nn_list     = NULL,
  output_path = NULL,
  add_assay   = TRUE
)
message("  [OK] Entropy estimation complete")

# 1.4 Compute mean entropy per cell
message("\nComputing mean entropy per cell...")
entropy_means <- colMeans(GetAssayData(seuratObj[["RNA_entropy"]]))
seuratObj$entropy_score <- entropy_means

# 1.5 Save outputs
message("\nSaving entropy-annotated object and metadata...")
write.csv(seuratObj@meta.data, file = OUTPUT_CSV, row.names = TRUE)
saveRDS(seuratObj, OUTPUT_RDS)
message("  [OK] ", OUTPUT_CSV)
message("  [OK] ", OUTPUT_RDS)


# ==============================================================================
# PART 2: RIDGE PLOTS + ANOVA
# ==============================================================================
banner("PART 2: Ridge plots + ANOVA")

# 2.1 Clean metadata
message("Cleaning metadata...")
df <- seuratObj@meta.data

# Harmonize age labels
df$age <- tolower(str_trim(df$age))
df$age <- dplyr::recode(df$age, !!!AGE_ALIASES)
df$age <- factor(df$age, levels = AGE_LEVELS, ordered = TRUE)

# Sex and celltype
df$sex      <- factor(tolower(str_trim(as.character(df$sex))))
df$celltype <- str_trim(as.character(df$celltype))

# Apply celltype order, append any extras at the end
present_ct <- intersect(CELLTYPE_LEVELS, unique(df$celltype))
extra_ct   <- setdiff(unique(df$celltype), CELLTYPE_LEVELS)
if (length(extra_ct) > 0) {
  message("  [INFO] Additional celltypes found: ", paste(extra_ct, collapse = ", "))
  present_ct <- c(present_ct, extra_ct)
}
df$celltype <- factor(df$celltype, levels = present_ct)

message(sprintf("  Data: %d cells, %d celltypes, %d samples",
                nrow(df), length(unique(df$celltype)),
                length(unique(df$sample))))

# 2.2 Pseudobulk aggregation for ANOVA
message("\nPseudobulk aggregation per sample x celltype...")
df_sample <- df %>%
  filter(!is.na(entropy_score), !is.na(age), !is.na(sex), !is.na(celltype)) %>%
  group_by(sample, age, sex, celltype) %>%
  summarise(entropy = mean(entropy_score, na.rm = TRUE), .groups = "drop")
message(sprintf("  [OK] %d sample-celltype rows", nrow(df_sample)))


# 2.3 Ridge plot helper
draw_ridge_plot <- function(df, target_sex, ncol_facets = 5) {

  df_sex <- df %>%
    filter(sex == target_sex, !is.na(age), !is.na(entropy_score))

  if (nrow(df_sex) < MIN_CELLS) {
    message("  [SKIP] ", target_sex, ": too few cells")
    return(NULL)
  }

  # Keep celltypes with enough cells
  keep_ct <- df_sex %>%
    count(celltype) %>%
    filter(n >= MIN_CELLS) %>%
    pull(celltype)

  df_sex <- df_sex %>% filter(celltype %in% keep_ct)
  df_sex$celltype <- factor(
    df_sex$celltype,
    levels = present_ct[present_ct %in% keep_ct]
  )

  p <- ggplot(df_sex, aes(x = entropy_score, y = age, fill = age)) +
    geom_density_ridges(
      scale          = 1.5,
      rel_min_height = 0.01,
      alpha          = 0.85,
      color          = "grey20",
      linewidth      = 0.2
    ) +
    scale_fill_manual(values = AGE_COLORS, drop = FALSE) +
    scale_y_discrete(limits = rev(AGE_LEVELS), drop = FALSE) +
    facet_wrap(~ celltype, ncol = ncol_facets) +
    theme_ridges() +
    theme(
      legend.position = "none",
      strip.text      = element_text(size = 10, face = "bold", family = "Arial"),
      axis.text.y     = element_text(size = 9,  family = "Arial"),
      axis.text.x     = element_text(size = 8,  family = "Arial"),
      plot.title      = element_text(size = 14, face = "bold", family = "Arial"),
      text            = element_text(family = "Arial")
    ) +
    labs(
      title = paste0("Entropy distribution by age - ", tools::toTitleCase(target_sex)),
      x     = "Entropy score",
      y     = "Age"
    )

  # Dynamic sizing
  n_ct  <- length(unique(df_sex$celltype))
  n_row <- ceiling(n_ct / ncol_facets)

  outfile <- file.path(OUT_DIR, paste0("Ridge_Entropy_", target_sex, ".pdf"))
  ggsave(
    outfile, p,
    width = ncol_facets * 3, height = max(4, n_row * 3.5),
    useDingbats = FALSE
  )
  message("  [OK] ", outfile)

  return(p)
}


# 2.4 Generate ridge plots
message("\nCreating ridge plots...")
p_ridge_female <- draw_ridge_plot(df, "female")
p_ridge_male   <- draw_ridge_plot(df, "male")


# 2.5 Per-celltype, per-sex one-way ANOVA on per-sample pseudobulked entropy
message("\nRunning per-celltype x per-sex ANOVA (BH FDR)...")
anova_results <- df_sample %>%
  group_by(celltype, sex) %>%
  summarise(
    n_samples = n(),
    n_ages    = n_distinct(age),
    p_ANOVA   = tryCatch(
      summary(aov(entropy ~ age))[[1]][["Pr(>F)"]][1],
      error = function(e) NA_real_
    ),
    .groups = "drop"
  ) %>%
  mutate(
    padj_ANOVA = p.adjust(p_ANOVA, method = "BH"),
    significance = case_when(
      is.na(padj_ANOVA)  ~ "",
      padj_ANOVA < 0.001 ~ "***",
      padj_ANOVA < 0.01  ~ "**",
      padj_ANOVA < 0.05  ~ "*",
      TRUE               ~ "n.s."
    )
  )

anova_file <- file.path(OUT_DIR, "Entropy_ANOVA_results.csv")
write.csv(anova_results, anova_file, row.names = FALSE)
message("  [OK] ", anova_file)


# ==============================================================================
# SUMMARY
# ==============================================================================
banner("PIPELINE COMPLETE")
message("Outputs:")
message("  ", OUTPUT_RDS)
message("  ", OUTPUT_CSV)
message("  ", file.path(OUT_DIR, "Ridge_Entropy_female.pdf"))
message("  ", file.path(OUT_DIR, "Ridge_Entropy_male.pdf"))
message("  ", anova_file)

message("\nANOVA significant results:")
sig_results <- anova_results %>%
  filter(significance %in% c("*", "**", "***"))
if (nrow(sig_results) > 0) {
  print(as.data.frame(sig_results))
} else {
  message("  None")
}

# Reproducibility
sessionInfo()
