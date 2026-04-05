#!/usr/bin/env Rscript
# ==============================================================================
# Step 12b: Cross-Celltype P2G Annotation Visualization
# ==============================================================================
#
# Description:
#   This script combines P2G annotation results from cell types 
#
# Prerequisites:
#   - Run Step 12 for EACH cell type to generate per-celltype CSVs:
#     - P2G_ChIPseeker_annotation.csv
#     - P2G_ENCODE_cCRE_annotation.csv
#
# Input:
#   - Per-celltype annotation CSVs (renamed with celltype suffix)
#     e.g., P2G_ChIPseeker_Hepatocyte.csv, P2G_ChIPseeker_Endothelial_01.csv
#
# Output:
#   - Combined stacked barplots (PDF)
#   - Combined annotation tables (CSV)
#
# ==============================================================================

suppressPackageStartupMessages({
    library(tidyverse)
    library(patchwork)
})

# ==============================================================================
# USER CONFIGURATION
# ==============================================================================

# Cell type order (left → right on x-axis)
CELL_ORDER <- c(
    "Hepatocyte",
    "Endothelial_01",
    "Endothelial_02",
    "Stellate",
    "Cholangiocyte_01",
    "Cholangiocyte_02",
    "Kupffer",
    "MoMFs",
    "Tcells",
    "Bcells"
)

# ==============================================================================
# HELPER FUNCTION
# ==============================================================================

banner <- function(text) {
    line <- paste(rep("=", 70), collapse = "")
    message("\n", line)
    message(text)
    message(line, "\n")
}

# ==============================================================================
# SECTION 1: CHIPSEEKER CROSS-CELLTYPE VISUALIZATION
# ==============================================================================

banner("SECTION 1: ChIPseeker Cross-Celltype Visualization")

# ------------------------------------------------------------
# 1.1 Read all per-celltype ChIPseeker CSVs
# ------------------------------------------------------------
chipseeker_files <- list.files(
    pattern = "^P2G_ChIPseeker_.*\\.csv$",
    full.names = TRUE
)

if (length(chipseeker_files) == 0) {
    message("[WARNING] No ChIPseeker files found. Skipping.")
} else {
    
    chipseeker_df <- purrr::map_dfr(chipseeker_files, function(f) {
        df <- read.csv(f, stringsAsFactors = FALSE)
        celltype <- gsub("^P2G_ChIPseeker_|\\.csv$", "", basename(f))
        df %>% mutate(celltype = celltype)
    })
    
    message(sprintf("[1.1] Loaded %d ChIPseeker files", length(chipseeker_files)))
    message(sprintf("[1.1] Cell types: %s", paste(unique(chipseeker_df$celltype), collapse = ", ")))
    
    # ------------------------------------------------------------
    # 1.2 Lock annotation order (bottom → top)
    # ------------------------------------------------------------
    annotation_order <- c(
        "Distal Intergenic",
        "Intron",
        "Exon",
        "5' UTR",
        "3' UTR",
        "Downstream",
        "Promoter (2-3kb)",
        "Promoter (1-2kb)",
        "Promoter (<=1kb)"
    )
    
    # ------------------------------------------------------------
    # 1.3 Standardize labels + complete missing categories
    # ------------------------------------------------------------
    chipseeker_df <- chipseeker_df %>%
        mutate(
            chipseeker_annotation = gsub("_", " ", chipseeker_annotation),
            chipseeker_annotation = trimws(chipseeker_annotation)
        ) %>%
        tidyr::complete(
            celltype,
            chipseeker_annotation = annotation_order,
            fill = list(n_peaks = 0, pct_of_all_p2g = 0)
        )
    
    # ------------------------------------------------------------
    # 1.4 Renormalize percentages per celltype
    # ------------------------------------------------------------
    chipseeker_df <- chipseeker_df %>%
        group_by(celltype) %>%
        mutate(
            pct_of_all_p2g = pct_of_all_p2g / sum(pct_of_all_p2g) * 100
        ) %>%
        ungroup()
    
    # ------------------------------------------------------------
    # 1.5 Set factor levels
    # ------------------------------------------------------------
    chipseeker_df$chipseeker_annotation <- factor(
        chipseeker_df$chipseeker_annotation,
        levels = annotation_order
    )
    
    chipseeker_df$celltype <- factor(
        chipseeker_df$celltype,
        levels = intersect(CELL_ORDER, unique(chipseeker_df$celltype))
    )
    
    # ------------------------------------------------------------
    # 1.6 Color palette
    # ------------------------------------------------------------
    annotation_colors <- c(
        "Distal Intergenic" = "#1B9E77",
        "Intron"            = "#66C2A5",
        "Exon"              = "#1F78B4",
        "5' UTR"            = "#8C8C8C",
        "3' UTR"            = "#B2ABD2",
        "Downstream"        = "#8E6BBE",
        "Promoter (2-3kb)"  = "#FBB4B9",
        "Promoter (1-2kb)"  = "#F768A1",
        "Promoter (<=1kb)"  = "#C51B7D"
    )
    
    # ------------------------------------------------------------
    # 1.7 Compute total peaks per celltype for x-axis labels
    # ------------------------------------------------------------
    celltype_totals <- chipseeker_df %>%
        group_by(celltype) %>%
        summarise(total_peaks = sum(n_peaks), .groups = "drop")
    
    celltype_labels <- setNames(
        paste0(
            celltype_totals$celltype,
            "\n(n=",
            format(celltype_totals$total_peaks, big.mark = ","),
            ")"
        ),
        celltype_totals$celltype
    )
    
    # ------------------------------------------------------------
    # 1.8 Build stacked barplot
    # ------------------------------------------------------------
    p_chipseeker <- ggplot(
        chipseeker_df,
        aes(x = celltype, y = pct_of_all_p2g, fill = chipseeker_annotation)
    ) +
        geom_bar(stat = "identity", width = 0.8, color = "black", linewidth = 0.25) +
        scale_fill_manual(values = annotation_colors, drop = FALSE) +
        scale_x_discrete(labels = celltype_labels) +
        scale_y_continuous(
            breaks = seq(0, 100, 25),
            labels = function(x) paste0(x, "%"),
            expand = c(0, 0)
        ) +
        labs(
            x = NULL,
            y = "% of P2G-linked peaks",
            fill = "ChIPseeker annotation"
        ) +
        theme_classic(base_size = 12) +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, lineheight = 0.9),
            axis.line = element_line(color = "black"),
            legend.position = "right",
            legend.title = element_text(size = 11),
            legend.text = element_text(size = 10)
        )
    
    print(p_chipseeker)
    
    ggsave(
        "P2G_ChIPseeker_all_celltypes_stacked_barplot.pdf",
        p_chipseeker,
        width = 9, height = 4.8,
        useDingbats = FALSE
    )
    message("[1.8] Saved: P2G_ChIPseeker_all_celltypes_stacked_barplot.pdf")
    
    # Save combined table
    write.csv(chipseeker_df, "P2G_ChIPseeker_all_celltypes_combined.csv", row.names = FALSE)
    message("[1.8] Saved: P2G_ChIPseeker_all_celltypes_combined.csv")
}


# ==============================================================================
# SECTION 2: ENCODE cCRE CROSS-CELLTYPE VISUALIZATION
# ==============================================================================

banner("SECTION 2: ENCODE cCRE Cross-Celltype Visualization")

# ------------------------------------------------------------
# 2.1 Read all per-celltype ENCODE cCRE CSVs
# ------------------------------------------------------------
ccre_files <- list.files(
    pattern = "^P2G_ENCODE_cCRE_.*\\.csv$",
    full.names = TRUE
)

if (length(ccre_files) == 0) {
    message("[WARNING] No ENCODE cCRE files found. Skipping.")
} else {
    
    ccre_df <- purrr::map_dfr(ccre_files, function(f) {
        df <- read.csv(f, stringsAsFactors = FALSE)
        celltype <- gsub("^P2G_ENCODE_cCRE_|\\.csv$", "", basename(f))
        df %>% mutate(celltype = celltype)
    })
    
    message(sprintf("[2.1] Loaded %d ENCODE cCRE files", length(ccre_files)))
    message(sprintf("[2.1] Cell types: %s", paste(unique(ccre_df$celltype), collapse = ", ")))
    
    # ------------------------------------------------------------
    # 2.2 Lock cCRE class order (bottom → top)
    # ------------------------------------------------------------
    ccre_class_order <- c(
        "not_annotated",
        "dELS",
        "dELS,CTCF-bound",
        "pELS",
        "pELS,CTCF-bound",
        "PLS",
        "PLS,CTCF-bound",
        "DNase-H3K4me3",
        "DNase-H3K4me3,CTCF-bound",
        "CTCF-only,CTCF-bound"
    )
    
    # ------------------------------------------------------------
    # 2.3 Set factor levels
    # ------------------------------------------------------------
    ccre_df$ccre_class <- factor(ccre_df$ccre_class, levels = ccre_class_order)
    
    ccre_df$celltype <- factor(
        ccre_df$celltype,
        levels = intersect(CELL_ORDER, unique(ccre_df$celltype))
    )
    
    # ------------------------------------------------------------
    # 2.4 Color palette
    # ------------------------------------------------------------
    ccre_colors <- c(
        "not_annotated"           = "#C7C7C7",
        "dELS"                    = "#A50F15",
        "dELS,CTCF-bound"         = "#FC9272",
        "pELS"                    = "lightpink",
        "pELS,CTCF-bound"         = "#A6D854",
        "PLS"                     = "#E6AB02",
        "PLS,CTCF-bound"          = "#FFD92F",
        "DNase-H3K4me3"           = "#E5C494",
        "DNase-H3K4me3,CTCF-bound"= "#B8A136",
        "CTCF-only,CTCF-bound"    = "#6A3D9A"
    )
    
    # ------------------------------------------------------------
    # 2.5 Compute total peaks per celltype for x-axis labels
    # ------------------------------------------------------------
    celltype_totals <- ccre_df %>%
        group_by(celltype) %>%
        summarise(total_peaks = sum(n_peaks), .groups = "drop")
    
    celltype_labels <- setNames(
        paste0(
            celltype_totals$celltype,
            "\n(n=",
            format(celltype_totals$total_peaks, big.mark = ","),
            ")"
        ),
        celltype_totals$celltype
    )
    
    # ------------------------------------------------------------
    # 2.6 Build stacked barplot
    # ------------------------------------------------------------
    p_ccre <- ggplot(
        ccre_df,
        aes(x = celltype, y = pct_of_all_p2g, fill = ccre_class)
    ) +
        geom_bar(stat = "identity", width = 0.8, color = "black", linewidth = 0.25) +
        scale_fill_manual(values = ccre_colors, drop = FALSE) +
        scale_x_discrete(labels = celltype_labels) +
        scale_y_continuous(
            breaks = seq(0, 100, 25),
            labels = function(x) paste0(x, "%"),
            expand = c(0, 0)
        ) +
        labs(
            x = NULL,
            y = "% of P2G-linked peaks",
            fill = "ENCODE cCRE class"
        ) +
        theme_classic(base_size = 12) +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, lineheight = 0.9),
            axis.line = element_line(color = "black"),
            legend.position = "right",
            legend.title = element_text(size = 11),
            legend.text = element_text(size = 10)
        )
    
    print(p_ccre)
    
    ggsave(
        "P2G_ENCODE_cCRE_all_celltypes_stacked_barplot.pdf",
        p_ccre,
        width = 9, height = 5,
        dpi = 300,
        useDingbats = FALSE
    )
    message("[2.6] Saved: P2G_ENCODE_cCRE_all_celltypes_stacked_barplot.pdf")
    
    # Save combined table
    write.csv(ccre_df, "P2G_ENCODE_cCRE_all_celltypes_combined.csv", row.names = FALSE)
    message("[2.6] Saved: P2G_ENCODE_cCRE_all_celltypes_combined.csv")
}


# ==============================================================================
# SUMMARY
# ==============================================================================

banner("STEP 12b COMPLETE")

message("Output files:")
message("  ✓ P2G_ChIPseeker_all_celltypes_stacked_barplot.pdf")
message("  ✓ P2G_ChIPseeker_all_celltypes_combined.csv")
message("  ✓ P2G_ENCODE_cCRE_all_celltypes_stacked_barplot.pdf")
message("  ✓ P2G_ENCODE_cCRE_all_celltypes_combined.csv")

message("\n============================================================")
