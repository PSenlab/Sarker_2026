library(Seurat)
library(Augur)
library(schard)

# ===============================
# 1. Load via h5ad Conversion
# ===============================
h5ad_path <- "/data/sarkern2/multiome_liver/final_object/final_rna_wnn.h5ad"

seuratObj <- schard::h5ad2seurat(h5ad_path)

DefaultAssay(seuratObj) <- "RNA"

# ===============================
# 2. Define Sexes & Age Comparisons
# ===============================
sexes          <- c("male", "female")
ages_to_compare <- c("mid_age", "old", "pre_geriatric", "geriatric")

# ===============================
# 3. Loop Through Each Combination
# ===============================
for (sx in sexes) {
  for (ag in ages_to_compare) {

    message("Running Augur for: ", sx, " (young vs ", ag, ")")

    # Subset to one sex and two age groups
    sub_obj <- subset(seuratObj,
                      subset = sex == sx & age %in% c("young", ag))

    # Extract normalized expression matrix
    expr <- as.matrix(GetAssayData(sub_obj, assay = "RNA", layer = "data"))

    # Metadata: must contain cell type and label columns
    meta <- sub_obj@meta.data[, c("celltype", "age")]

    # Confirm label levels
    meta$age <- factor(meta$age, levels = c("young", ag))

    # Run Augur
    augur_res <- calculate_auc(
      expr,
      meta,
      cell_type_col = "celltype",
      label_col     = "age",
      n_threads     = 15
    )

    # Save AUC results
    out_file <- paste0("augur_auc_", sx, "_young_vs_", ag, ".csv")
    write.csv(augur_res$AUC, out_file, row.names = FALSE)

    message("Saved: ", out_file)
  }
}