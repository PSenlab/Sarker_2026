#!/usr/bin/env Rscript
# ==============================================================================
# Augur Cell-Type Prioritization Across Aging (Young vs. Each Older Group)
# ==============================================================================
#
# Description:
#   For each sex (male, female) and each older age group (mid_age, old,
#   pre_geriatric, geriatric), runs Augur to compute per-cell-type AUC
#   scores for distinguishing young vs. that older group from the RNA
#   modality of the snMultiome object.
#
# Input:
#   - Annotated AnnData (.h5ad) with celltype, sex, age labels
#
# Output:
#   - augur_auc_<sex>_young_vs_<age>.csv  (per-cell-type AUC values)
#
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Augur)
  library(schard)
})


# ==============================================================================
# CONFIGURATION - UPDATE THIS PATH
# ==============================================================================
H5AD_PATH <- "integrated_scvi.h5ad"

SEXES           <- c("male", "female")
AGES_TO_COMPARE <- c("mid_age", "old", "pre_geriatric", "geriatric")
N_THREADS       <- 15


# ==============================================================================
# STEP 1: LOAD H5AD VIA SCHARD
# ==============================================================================
message("\n", paste(rep("=", 70), collapse = ""))
message("STEP 1: Load h5ad and convert to Seurat")
message(paste(rep("=", 70), collapse = ""))

seuratObj <- schard::h5ad2seurat(H5AD_PATH)
DefaultAssay(seuratObj) <- "RNA"
message("  [OK] Loaded: ", ncol(seuratObj), " cells x ", nrow(seuratObj), " genes")


# ==============================================================================
# STEP 2: AUGUR PER SEX x AGE COMPARISON
# ==============================================================================
message("\n", paste(rep("=", 70), collapse = ""))
message("STEP 2: Augur AUC per (sex, young vs older group)")
message(paste(rep("=", 70), collapse = ""))

for (sx in SEXES) {
  for (ag in AGES_TO_COMPARE) {

    message("\n--- ", sx, ": young vs ", ag, " ---")

    # Subset to one sex and two age groups
    sub_obj <- subset(seuratObj,
                      subset = sex == sx & age %in% c("young", ag))
    message("  Cells: ", ncol(sub_obj))

    if (ncol(sub_obj) == 0) {
      message("  [SKIP] No cells for this combination")
      next
    }

    # Normalized expression matrix
    expr <- as.matrix(GetAssayData(sub_obj, assay = "RNA", layer = "data"))

    # Metadata: cell type and label
    meta <- sub_obj@meta.data[, c("celltype", "age")]
    meta$age <- factor(meta$age, levels = c("young", ag))

    # Run Augur
    augur_res <- calculate_auc(
      expr,
      meta,
      cell_type_col = "celltype",
      label_col     = "age",
      n_threads     = N_THREADS
    )

    # Save AUC table
    out_file <- paste0("augur_auc_", sx, "_young_vs_", ag, ".csv")
    write.csv(augur_res$AUC, out_file, row.names = FALSE)
    message("  [OK] ", out_file)
  }
}


# ==============================================================================
# DONE
# ==============================================================================
message("\n", paste(rep("=", 70), collapse = ""))
message("AUGUR PIPELINE COMPLETE")
message(paste(rep("=", 70), collapse = ""))

sessionInfo()
