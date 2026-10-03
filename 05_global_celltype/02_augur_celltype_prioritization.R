#!/usr/bin/env Rscript
# ==============================================================================
# Augur Cell-Type Prioritization Across Aging (Young vs. Each Older Group)
# ==============================================================================
#
# Description:
#   For each sex (male, female) and each older age group (mid_age, old,
#   pre_geriatric, geriatric), runs Augur to compute per-cell-type AUC
#   scores for distinguishing young vs. that older group from the RNA
#   modality of the snMultiome object. The AUC tables are then combined
#   and drawn as a dot plot (dot size = AUC, color = comparison, by sex).
#
# Input:
#   - rna_wnn.h5ad (from 04_multiome_integration/02_wnn_integration.py)
#     with cell_type, sex, age labels
#
# Output:
#   - augur_auc_<sex>_young_vs_<age>.csv  (per-cell-type AUC values)
#   - augur_auc_all.csv                   (all comparisons in one table)
#   - augur_auc_dotplot.pdf
#
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Augur)
  library(schard)
  library(dplyr)
  library(ggplot2)
})


# ==============================================================================
# CONFIGURATION - UPDATE THIS PATH
# ==============================================================================
H5AD_PATH <- "rna_wnn.h5ad"

SEXES           <- c("male", "female")
AGES_TO_COMPARE <- c("mid_age", "old", "pre_geriatric", "geriatric")
N_THREADS       <- 15

# Dot plot: comparison labels and colors
COMP_LABELS <- c(
  mid_age       = "mid-age vs young",
  old           = "old vs young",
  pre_geriatric = "pre-geriatric vs young",
  geriatric     = "geriatric vs young"
)
COMP_COLORS <- c(
  "mid-age vs young"       = "#5B9BD5",
  "old vs young"           = "#4DAF4A",
  "pre-geriatric vs young" = "#F0A73A",
  "geriatric vs young"     = "#E8483F"
)

# Dot plot: cell type order (top to bottom)
CELLTYPE_ORDER <- c(
  "mesenchymal",
  "non-resident myeloid",
  "T/ILC cells",
  "B cells",
  "Kupffer 01",
  "hepatocyte",
  "Kupffer 02",
  "endothelial",
  "cholangiocyte 01"
)


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
    meta <- sub_obj@meta.data[, c("cell_type", "age")]
    meta$age <- factor(meta$age, levels = c("young", ag))

    # Run Augur
    augur_res <- calculate_auc(
      expr,
      meta,
      cell_type_col = "cell_type",
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
# STEP 3: COMBINE AUGUR AUC TABLES
# ==============================================================================
message("\n", paste(rep("=", 70), collapse = ""))
message("STEP 3: Combine Augur AUC tables")
message(paste(rep("=", 70), collapse = ""))

auc_list <- list()
for (sx in SEXES) {
  for (ag in AGES_TO_COMPARE) {
    f <- paste0("augur_auc_", sx, "_young_vs_", ag, ".csv")
    if (!file.exists(f)) {
      message("  [SKIP] Not found: ", f)
      next
    }
    df <- read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
    stopifnot(all(c("cell_type", "auc") %in% colnames(df)))
    df$sex        <- sx
    df$comparison <- COMP_LABELS[[ag]]
    auc_list[[paste(sx, ag)]] <- df[, c("sex", "comparison", "cell_type", "auc")]
  }
}

if (length(auc_list) == 0) stop("No Augur AUC tables found")

auc_all <- bind_rows(auc_list)
write.csv(auc_all, "augur_auc_all.csv", row.names = FALSE)
message("  [OK] augur_auc_all.csv (", nrow(auc_all), " rows)")

not_plotted <- setdiff(unique(auc_all$cell_type), CELLTYPE_ORDER)
if (length(not_plotted) > 0) {
  message("  Cell types in the tables but not in CELLTYPE_ORDER (not plotted): ",
          paste(not_plotted, collapse = ", "))
}


# ==============================================================================
# STEP 4: AUGUR AUC DOT PLOT
# ==============================================================================
message("\n", paste(rep("=", 70), collapse = ""))
message("STEP 4: Augur AUC dot plot")
message(paste(rep("=", 70), collapse = ""))

plot_df <- auc_all %>%
  filter(cell_type %in% CELLTYPE_ORDER) %>%
  mutate(
    cell_type  = factor(cell_type, levels = rev(CELLTYPE_ORDER)),
    comparison = factor(comparison, levels = unname(COMP_LABELS)),
    sex        = factor(sex, levels = SEXES),
    auc_plot   = pmin(pmax(auc, 0.5), 1)   # keep AUCs below 0.5 on the scale
  )

p <- ggplot(plot_df, aes(x = comparison, y = cell_type, size = auc_plot, fill = comparison)) +
  geom_point(shape = 21, color = "black", stroke = 0.4) +
  facet_wrap(~ sex, nrow = 1) +
  scale_fill_manual(values = COMP_COLORS, name = NULL) +
  scale_radius(
    name   = "Augur AUC score",
    limits = c(0.5, 1),
    breaks = c(0.5, 0.75, 1),
    labels = c("0.5", "0.75", "1"),
    range  = c(0.5, 11)
  ) +
  guides(
    fill = guide_legend(override.aes = list(size = 4), ncol = 1, order = 1),
    size = guide_legend(ncol = 1, order = 2, title.position = "top")
  ) +
  labs(x = NULL, y = NULL) +
  theme_bw(base_size = 12, base_family = "Arial") +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    axis.text.x        = element_blank(),
    axis.ticks.x       = element_blank(),
    strip.background   = element_rect(fill = "grey90", color = "black"),
    legend.position    = "bottom",
    legend.box         = "horizontal",
    legend.key         = element_blank()
  )

ggsave("augur_auc_dotplot.pdf", p, width = 7, height = 6, device = cairo_pdf)
message("  [OK] augur_auc_dotplot.pdf")


# ==============================================================================
# DONE
# ==============================================================================
message("\n", paste(rep("=", 70), collapse = ""))
message("AUGUR PIPELINE COMPLETE")
message(paste(rep("=", 70), collapse = ""))

sessionInfo()
