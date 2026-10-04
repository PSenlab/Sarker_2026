#!/usr/bin/env Rscript
# ==============================================================================
# Module Score UMAPs (WNN embedding)
# ==============================================================================
#
# Description:
#   Scores gene sets with Seurat::AddModuleScore ("combined module score") and
#   draws them on the WNN UMAP with scCustomize. A set is drawn either as one
#   UMAP, or split by age with a male row and a female row.
#
#   Gene sets:
#     1. Cell type markers, all cells (one UMAP per gene):
#          Abcc2, Spp1, Dcn, Ptprb, Cd5l, Mctp1, Ms4a4b, Ebf1
#     2. Hepatocyte P2G modules C1-C4 (one UMAP per module), gene lists from
#        03_atac/09_archr_p2g_analysis.R (hepatocyte run)
#     3. Hepatocyte zonation markers (one UMAP per gene):
#          Cyp2f2, Cdh1, Hal, Hamp2, Cyp2e1, Glul, Cyp7a1
#     4. Hepatocyte hep_05 cluster genes by age and sex (one grid per gene):
#          Cd14, Tox2, Anxa2
#     5. Ascl1 figure genes by age and sex (one grid per gene):
#          Ascl1, Mki67, Vwf
#
# Input:
#   - rna_wnn.h5ad (from 04_multiome_integration/02_wnn_integration.py)
#     with cell_type, sex, age and the X_wnn embedding
#   - P2G_C1_genes.csv ... P2G_C4_genes.csv
#     (from 03_atac/09_archr_p2g_analysis.R, hepatocyte run)
#
# Output (in OUTPUT_DIR):
#   - markers_<gene>.pdf
#   - hep_P2G_C<k>.pdf
#   - hep_zonation_<gene>.pdf
#   - hep_<gene>_by_age_sex.pdf
#
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(scCustomize)
  library(schard)
  library(ggplot2)
  library(patchwork)
})

set.seed(1)


# ==============================================================================
# CONFIGURATION - UPDATE THESE PATHS
# ==============================================================================
H5AD_PATH  <- "rna_wnn.h5ad"
P2G_DIR    <- "Step12_P2G_Analysis"   # folder with P2G_C*_genes.csv (hepatocyte run)
OUTPUT_DIR <- "module_score_umaps"
REDUCTION  <- "Xwnn_"                 # schard name for obsm["X_wnn"]

AGE_LEVELS <- c("young", "mid_age", "old", "pre_geriatric", "geriatric")
SEX_LEVELS <- c("male", "female")

MARKER_GENES   <- c("Abcc2", "Spp1", "Dcn", "Ptprb", "Cd5l", "Mctp1", "Ms4a4b", "Ebf1")
ZONATION_GENES <- c("Cyp2f2", "Cdh1", "Hal", "Hamp2", "Cyp2e1", "Glul", "Cyp7a1")
HEP05_GENES    <- c("Cd14", "Tox2", "Anxa2")    # hep_05 cluster genes
ASCL1_GENES    <- c("Ascl1", "Mki67", "Vwf")
P2G_CLUSTERS   <- 1:4

dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)


# ==============================================================================
# HELPERS
# ==============================================================================
banner <- function(text) {
  line <- paste(rep("=", 70), collapse = "")
  message("\n", line)
  message(text)
  message(line)
}

# Add a module score for one gene set; returns the object with the new
# column "<name>1", or NULL if none of the genes are present
add_score <- function(obj, genes, name) {
  present <- genes[genes %in% rownames(obj)]
  missing <- setdiff(genes, present)
  if (length(missing) > 0) {
    message("    not found: ", paste(missing, collapse = ", "))
  }
  if (length(present) == 0) return(NULL)
  AddModuleScore(object = obj, features = list(present), name = name)
}

# Single UMAP
plot_one <- function(obj, feature, title, out_pdf) {
  p <- FeaturePlot_scCustom(
    obj,
    reduction = REDUCTION,
    features  = feature,
    order     = TRUE
  ) + ggtitle(title)
  ggsave(out_pdf, p, width = 5.5, height = 5)
  message("  [OK] ", out_pdf)
}

# One row per sex (male top, female bottom), split by age
plot_age_sex <- function(obj, feature, title, out_pdf) {
  p_male <- FeaturePlot_scCustom(
    subset(obj, subset = sex == "male"),
    reduction = REDUCTION,
    features  = feature,
    split.by  = "age",
    order     = TRUE
  )
  p_female <- FeaturePlot_scCustom(
    subset(obj, subset = sex == "female"),
    reduction = REDUCTION,
    features  = feature,
    split.by  = "age",
    order     = TRUE
  )
  p <- wrap_plots(p_male, p_female, ncol = 1) + plot_annotation(title = title)
  ggsave(out_pdf, p, width = 20, height = 8)
  message("  [OK] ", out_pdf)
}

safe_name <- function(x) gsub("[^A-Za-z0-9]", "_", x)


# ==============================================================================
# STEP 1: LOAD DATA
# ==============================================================================
banner("STEP 1: Load h5ad and convert to Seurat")

seurat_obj <- schard::h5ad2seurat(H5AD_PATH)
DefaultAssay(seurat_obj) <- "RNA"

seurat_obj$age <- factor(seurat_obj$age, levels = AGE_LEVELS)
seurat_obj$sex <- factor(seurat_obj$sex, levels = SEX_LEVELS)

if (!REDUCTION %in% names(seurat_obj@reductions)) {
  warning(sprintf("'%s' not found; using '%s'",
                  REDUCTION, names(seurat_obj@reductions)[1]))
  REDUCTION <- names(seurat_obj@reductions)[1]
}

hep <- subset(seurat_obj, subset = cell_type == "hepatocyte")

message("  [OK] All cells:   ", ncol(seurat_obj))
message("  [OK] Hepatocytes: ", ncol(hep))
message("  [OK] Reduction:   ", REDUCTION)


# ==============================================================================
# STEP 2: CELL TYPE MARKERS (ALL CELLS)
# ==============================================================================
banner("STEP 2: Cell type markers (all cells)")

for (g in MARKER_GENES) {
  message("  ", g)
  nm  <- paste0("marker_", safe_name(g), "_")
  res <- add_score(seurat_obj, g, nm)
  if (is.null(res)) next
  seurat_obj <- res
  plot_one(seurat_obj, paste0(nm, "1"), g,
           file.path(OUTPUT_DIR, paste0("markers_", safe_name(g), ".pdf")))
}


# ==============================================================================
# STEP 3: HEPATOCYTE P2G MODULES C1-C4
# ==============================================================================
banner("STEP 3: Hepatocyte P2G modules")

for (k in P2G_CLUSTERS) {
  f <- file.path(P2G_DIR, sprintf("P2G_C%d_genes.csv", k))
  if (!file.exists(f)) {
    message("  [SKIP] not found: ", f)
    next
  }
  genes <- unique(na.omit(read.csv(f, stringsAsFactors = FALSE)$gene))
  message(sprintf("  P2G C%d: %d genes", k, length(genes)))
  nm  <- sprintf("P2G_C%d_", k)
  res <- add_score(hep, genes, nm)
  if (is.null(res)) next
  hep <- res
  plot_one(hep, paste0(nm, "1"), sprintf("P2G C%d", k),
           file.path(OUTPUT_DIR, sprintf("hep_P2G_C%d.pdf", k)))
}


# ==============================================================================
# STEP 4: HEPATOCYTE ZONATION MARKERS
# ==============================================================================
banner("STEP 4: Hepatocyte zonation markers")

for (g in ZONATION_GENES) {
  message("  ", g)
  nm  <- paste0("zonation_", safe_name(g), "_")
  res <- add_score(hep, g, nm)
  if (is.null(res)) next
  hep <- res
  plot_one(hep, paste0(nm, "1"), g,
           file.path(OUTPUT_DIR, paste0("hep_zonation_", safe_name(g), ".pdf")))
}


# ==============================================================================
# STEP 5: HEPATOCYTE HEP_05 AND ASCL1 FIGURE GENES BY AGE AND SEX
# ==============================================================================
banner("STEP 5: Hepatocyte hep_05 and Ascl1 figure genes by age and sex")

for (g in c(HEP05_GENES, ASCL1_GENES)) {
  message("  ", g)
  nm  <- paste0("hep_", safe_name(g), "_")
  res <- add_score(hep, g, nm)
  if (is.null(res)) next
  hep <- res
  plot_age_sex(hep, paste0(nm, "1"), g,
               file.path(OUTPUT_DIR, paste0("hep_", safe_name(g), "_by_age_sex.pdf")))
}


# ==============================================================================
# DONE
# ==============================================================================
banner("MODULE SCORE UMAPS COMPLETE")
message("  Output: ", OUTPUT_DIR, "/")

sessionInfo()
