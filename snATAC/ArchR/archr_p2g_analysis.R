#!/usr/bin/env Rscript
# ==============================================================================
# Step 12: P2G Extraction, Clustering, Annotation & Cross-Celltype Visualization
# ==============================================================================
#
# Description:
#   Extracts Peak-to-Gene (P2G) linkages from ArchR projects, performs k-means
#   clustering on P2G heatmaps, annotates peaks using ENCODE cCREs and
#   ChIPseeker, and generates cross-celltype stacked barplots.
#
# Prerequisites:
#   - Run Step 11 (archr_peak_calling_p2g.R) to generate P2G linkages
#   - Download ENCODE cCREs: mm10-cCREs.bed from ENCODE SCREEN
#   - Run this script for EACH cell type, then run Section 13-14 once
#     after all cell types are processed
#
# Input:
#   - ArchR project with P2G linkages (from Step 11b)
#   - mm10-cCREs.bed (ENCODE candidate cis-Regulatory Elements)
#   - Per-celltype annotation CSVs (for cross-celltype visualization)
#
# Output:
#   - P2G cluster assignments (k-means)
#   - Per-cluster genes, peaks, links (CSV + Excel)
#   - BED files for external analysis
#   - ENCODE cCRE annotation summary
#   - ChIPseeker genomic annotation summary
#   - Cross-celltype stacked barplots (ChIPseeker + ENCODE cCRE)
#
# Pipeline Overview:
#   12.1:  Load ArchR Project
#   12.2:  Extract P2G Links
#   12.3:  P2G Heatmap with K-means Clustering
#   12.4:  Export Per-Cluster Data (CSV + Excel)
#   12.5:  Export BED Files
#   12.6:  ENCODE cCRE Annotation
#   12.7:  ChIPseeker Annotation
#   12.8:  Save RDS and Final Summary
#   12.9:  ChIPseeker Stacked Barplot (Single Cell Type)
#   12.10: ENCODE cCRE Stacked Barplot (Single Cell Type)
#   12.11: ChIPseeker Cross-Celltype Stacked Barplot
#   12.12: ENCODE cCRE Cross-Celltype Stacked Barplot
#
# ==============================================================================

# ==============================================================================
# SETUP AND CONFIGURATION
# ==============================================================================

suppressPackageStartupMessages({
    library(ArchR)
    library(GenomicRanges)
    library(dplyr)
    library(tidyr)
    library(stringr)
    library(openxlsx)
    library(data.table)
    library(ggplot2)
    library(patchwork)
    library(ChIPseeker)
    library(TxDb.Mmusculus.UCSC.mm10.knownGene)
    library(org.Mm.eg.db)
    library(BSgenome.Mmusculus.UCSC.mm10)
})

# Set Global Parameters
addArchRThreads(threads = 60)
addArchRGenome("mm10")

# ==============================================================================
# USER CONFIGURATION - MODIFY THESE PATHS
# ==============================================================================

# Path to ArchR project from Step 11b (with P2G linkages)
ARCHR_PROJECT_PATH <- "path/to/Step11b_CellType_P2G"

# Path to ENCODE cCREs file
ENCODE_CCRE_PATH <- "mm10-cCREs.bed"

# Output directory
OUTPUT_DIR <- "Step12_P2G_Analysis"

# Analysis parameters
K_CLUSTERS <- 3                  # Number of k-means clusters for P2G heatmap
P2G_COR_CUTOFF <- 0.45           # Correlation cutoff for P2G extraction
COA_COR_CUTOFF <- 0.5            # Correlation cutoff for co-accessibility

# Cell type order for cross-celltype visualization (left → right on x-axis)
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
# HELPER FUNCTIONS
# ==============================================================================

banner <- function(text) {
    line <- paste(rep("=", 70), collapse = "")
    message("\n", line)
    message(text)
    message(line, "\n")
}

ensure_dir <- function(dir_path) {
    if (!dir.exists(dir_path)) {
        dir.create(dir_path, recursive = TRUE)
    }
    invisible(dir_path)
}

# ==============================================================================
# 12.1: LOAD ARCHR PROJECT
# ==============================================================================

banner("12.1: Load ArchR Project")

if (!dir.exists(ARCHR_PROJECT_PATH)) {
    stop(sprintf("[FATAL] Project path does not exist: %s", ARCHR_PROJECT_PATH))
}

proj <- loadArchRProject(path = ARCHR_PROJECT_PATH)
message(sprintf("[12.1] Loaded project: %s", ARCHR_PROJECT_PATH))
message(sprintf("[12.1] Total cells: %d", nCells(proj)))

# Create output directory
ensure_dir(OUTPUT_DIR)
setwd(OUTPUT_DIR)

# ==============================================================================
# 12.2: EXTRACT P2G LINKS
# ==============================================================================

banner("12.2: Extract P2G Links")

# Get co-accessibility
cA <- getCoAccessibility(
    ArchRProj = proj,
    corCutOff = COA_COR_CUTOFF,
    resolution = 1,
    returnLoops = FALSE
)
message(sprintf("[12.2] Co-accessibility links: %d", length(cA)))

# Get P2G links
p2g <- getPeak2GeneLinks(
    ArchRProj = proj,
    corCutOff = P2G_COR_CUTOFF,
    resolution = 1,
    returnLoops = FALSE
)

# Add gene and peak names
p2g$geneName <- mcols(metadata(p2g)$geneSet)$name[p2g$idxRNA]
p2g$peakName <- (metadata(p2g)$peakSet %>% {
    paste0(seqnames(.), ":", start(.), "-", end(.))
})[p2g$idxATAC]

message(sprintf("[12.2] P2G links extracted: %d", length(p2g)))

# ==============================================================================
# 12.3: P2G HEATMAP WITH K-MEANS CLUSTERING
# ==============================================================================

banner("12.3: P2G Heatmap with K-means Clustering")

# Create sex_age grouping if not exists
if (!"sex_age" %in% colnames(getCellColData(proj))) {
    cellMeta <- getCellColData(proj)
    proj$sex_age <- paste0(cellMeta$sex, "_", cellMeta$age)
}

# Generate P2G heatmap
p2g_heatmap <- plotPeak2GeneHeatmap(
    ArchRProj = proj,
    groupBy = "sex_age",
    k = K_CLUSTERS
)

# Save heatmap
pdf(sprintf("P2G_heatmap_k%d.pdf", K_CLUSTERS), width = 8, height = 4)
print(p2g_heatmap)
dev.off()
message(sprintf("[12.3] Saved: P2G_heatmap_k%d.pdf", K_CLUSTERS))

# Get matrices with cluster assignments
p2g_data <- plotPeak2GeneHeatmap(
    ArchRProj = proj,
    groupBy = "sex_age",
    k = K_CLUSTERS,
    returnMatrices = TRUE
)

# Extract cluster assignments
clusters <- p2g_data$ATAC$kmeansId
message(sprintf("[12.3] Cluster distribution:"))
print(table(clusters))

# ==============================================================================
# 12.4: EXPORT PER-CLUSTER DATA (CSV + EXCEL)
# ==============================================================================

banner("12.4: Export Per-Cluster Data")

# Get P2G links as dataframe
p2g_links <- as.data.frame(p2g_data$Peak2GeneLinks)
p2g_links$cluster <- clusters

message(sprintf("[12.4] Total P2G links: %d", nrow(p2g_links)))

# Parse peak coordinates
p2g_links$peak_chr <- gsub(":.*", "", p2g_links$peak)
p2g_links$peak_start <- as.integer(gsub(".*:(\\d+)-.*", "\\1", p2g_links$peak))
p2g_links$peak_end <- as.integer(gsub(".*-", "", p2g_links$peak))

# Get peak annotations from ArchR
peakSet <- getPeakSet(proj)
peak_annot <- data.frame(
    peak = paste0(seqnames(peakSet), ":", start(peakSet), "-", end(peakSet)),
    peakType = peakSet$peakType,
    nearestGene = peakSet$nearestGene,
    distToTSS = peakSet$distToTSS
)

# Add annotations to p2g_links
p2g_links <- p2g_links %>%
    left_join(peak_annot, by = "peak")

message("[12.4] Peak type distribution:")
print(table(p2g_links$peakType))

# Initialize storage
cluster_summary <- data.frame()
all_genes_list <- list()
all_peaks_list <- list()
all_links_list <- list()

# Process each cluster
for (c in sort(unique(clusters))) {
    
    p2g_c <- p2g_links[p2g_links$cluster == c, ]
    
    # Unique genes and peaks
    genes_c <- unique(p2g_c$gene)
    peaks_c <- p2g_c %>%
        select(peak, peak_chr, peak_start, peak_end, peakType, nearestGene, distToTSS) %>%
        distinct()
    
    # Peak type counts
    peak_types <- table(peaks_c$peakType)
    
    # Store for Excel
    all_genes_list[[paste0("Cluster_", c)]] <- data.frame(gene = genes_c, cluster = c)
    all_peaks_list[[paste0("Cluster_", c)]] <- data.frame(peak = unique(p2g_c$peak), cluster = c)
    all_links_list[[paste0("Cluster_", c)]] <- p2g_c[, c("peak", "gene", "Correlation", "FDR", "peakType", "distToTSS", "cluster")]
    
    # Save individual CSVs
    write.csv(data.frame(gene = genes_c),
              sprintf("P2G_C%d_genes.csv", c), row.names = FALSE)
    write.csv(peaks_c,
              sprintf("P2G_C%d_peaks_annotated.csv", c), row.names = FALSE)
    write.csv(p2g_c[, c("peak", "gene", "Correlation", "FDR", "peakType", "distToTSS")],
              sprintf("P2G_C%d_links.csv", c), row.names = FALSE)
    
    # Summary
    cluster_summary <- rbind(cluster_summary, data.frame(
        cluster = c,
        n_links = nrow(p2g_c),
        n_peaks = nrow(peaks_c),
        n_genes = length(genes_c),
        n_promoter = ifelse("Promoter" %in% names(peak_types), peak_types["Promoter"], 0),
        n_distal = ifelse("Distal" %in% names(peak_types), peak_types["Distal"], 0),
        n_intronic = ifelse("Intronic" %in% names(peak_types), peak_types["Intronic"], 0),
        n_exonic = ifelse("Exonic" %in% names(peak_types), peak_types["Exonic"], 0),
        mean_cor = round(mean(p2g_c$Correlation), 3)
    ))
    
    message(sprintf("  C%d: %d links | %d peaks | %d genes", c, nrow(p2g_c), nrow(peaks_c), length(genes_c)))
}

# Add percentages
cluster_summary$pct_promoter <- round(cluster_summary$n_promoter / cluster_summary$n_peaks * 100, 1)
cluster_summary$pct_distal <- round(cluster_summary$n_distal / cluster_summary$n_peaks * 100, 1)

# Save cluster summary
write.csv(cluster_summary, "P2G_cluster_summary.csv", row.names = FALSE)

# ------------------------------------------------------------
# Create Excel workbooks
# ------------------------------------------------------------

# GENES Excel
wb_genes <- createWorkbook()
all_genes_combined <- do.call(rbind, all_genes_list) %>% arrange(cluster, gene)
addWorksheet(wb_genes, "All_Genes")
writeData(wb_genes, "All_Genes", all_genes_combined)

for (c in sort(unique(clusters))) {
    addWorksheet(wb_genes, paste0("C", c, "_Genes"))
    writeData(wb_genes, paste0("C", c, "_Genes"), all_genes_list[[paste0("Cluster_", c)]])
}
addWorksheet(wb_genes, "Summary")
writeData(wb_genes, "Summary", cluster_summary)
saveWorkbook(wb_genes, "P2G_AllClusters_GENES.xlsx", overwrite = TRUE)
message("[12.4] Saved: P2G_AllClusters_GENES.xlsx")

# PEAKS Excel
wb_peaks <- createWorkbook()
all_peaks_combined <- do.call(rbind, all_peaks_list) %>% arrange(cluster, peak)
addWorksheet(wb_peaks, "All_Peaks")
writeData(wb_peaks, "All_Peaks", all_peaks_combined)

for (c in sort(unique(clusters))) {
    addWorksheet(wb_peaks, paste0("C", c, "_Peaks"))
    writeData(wb_peaks, paste0("C", c, "_Peaks"), all_peaks_list[[paste0("Cluster_", c)]])
}
addWorksheet(wb_peaks, "Summary")
writeData(wb_peaks, "Summary", cluster_summary)
saveWorkbook(wb_peaks, "P2G_AllClusters_PEAKS.xlsx", overwrite = TRUE)
message("[12.4] Saved: P2G_AllClusters_PEAKS.xlsx")

# LINKS Excel
wb_links <- createWorkbook()
all_links_combined <- do.call(rbind, all_links_list) %>% arrange(cluster, desc(Correlation))
addWorksheet(wb_links, "All_Links")
writeData(wb_links, "All_Links", all_links_combined)

for (c in sort(unique(clusters))) {
    addWorksheet(wb_links, paste0("C", c, "_Links"))
    writeData(wb_links, paste0("C", c, "_Links"), 
              all_links_list[[paste0("Cluster_", c)]] %>% arrange(desc(Correlation)))
}
addWorksheet(wb_links, "Summary")
writeData(wb_links, "Summary", cluster_summary)
saveWorkbook(wb_links, "P2G_AllClusters_LINKS.xlsx", overwrite = TRUE)
message("[12.4] Saved: P2G_AllClusters_LINKS.xlsx")

# ==============================================================================
# 12.5: EXPORT BED FILES
# ==============================================================================

banner("12.5: Export BED Files")

# Per-cluster BED files
for (c in sort(unique(clusters))) {
    p2g_c <- p2g_links[p2g_links$cluster == c, ]
    bed_c <- p2g_c %>%
        select(peak_chr, peak_start, peak_end, peak, peakType) %>%
        distinct()
    
    write.table(bed_c,
                sprintf("P2G_C%d_peaks.bed", c),
                sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
}

# Combined BED file (all peaks)
all_peaks <- p2g_links %>%
    select(peak_chr, peak_start, peak_end, peak, cluster) %>%
    distinct()

write.table(all_peaks, "P2G_all_peaks.bed",
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

message(sprintf("[12.5] Saved BED files: %d clusters + combined", length(unique(clusters))))

# ==============================================================================
# 12.6: ENCODE cCRE ANNOTATION
# ==============================================================================

banner("12.6: ENCODE cCRE Annotation")

if (!file.exists(ENCODE_CCRE_PATH)) {
    warning(sprintf("[12.6] ENCODE cCRE file not found: %s", ENCODE_CCRE_PATH))
    warning("[12.6] Skipping ENCODE annotation. Download from ENCODE SCREEN.")
} else {
    
    # Load ENCODE cCREs
    ccre_df <- fread(ENCODE_CCRE_PATH, sep = "\t", header = FALSE)
    colnames(ccre_df) <- c("chr", "start", "end", "dnase_id", "enhancer_id", "ccre_class")
    
    encode_ccres <- GRanges(
        seqnames = ccre_df$chr,
        ranges = IRanges(start = ccre_df$start + 1, end = ccre_df$end),
        ccre_class = ccre_df$ccre_class
    )
    message(sprintf("[12.6] ENCODE cCREs loaded: %d", length(encode_ccres)))
    
    # Create GRanges for P2G peaks
    p2g_peaks <- GRanges(
        seqnames = p2g_links$peak_chr,
        ranges = IRanges(start = p2g_links$peak_start, end = p2g_links$peak_end)
    )
    p2g_peaks <- unique(p2g_peaks)
    
    total_p2g_peaks <- length(p2g_peaks)
    message(sprintf("[12.6] Unique P2G peaks: %d", total_p2g_peaks))
    
    # Harmonize chromosome style
    seqlevelsStyle(encode_ccres) <- "UCSC"
    seqlevelsStyle(p2g_peaks) <- "UCSC"
    
    # Find overlaps
    ov <- findOverlaps(p2g_peaks, encode_ccres, ignore.strand = TRUE)
    
    overlap_df <- data.frame(
        peak_idx = queryHits(ov),
        ccre_class = mcols(encode_ccres)$ccre_class[subjectHits(ov)],
        stringsAsFactors = FALSE
    )
    
    # Define priority (one class per peak)
    ccre_priority <- c(
        "PLS", "PLS,CTCF-bound",
        "DNase-H3K4me3", "DNase-H3K4me3,CTCF-bound",
        "pELS", "pELS,CTCF-bound",
        "dELS", "dELS,CTCF-bound",
        "CTCF-only,CTCF-bound"
    )
    
    overlap_df$priority <- match(overlap_df$ccre_class, ccre_priority)
    
    overlap_df_unique <- overlap_df %>%
        filter(!is.na(priority)) %>%
        arrange(peak_idx, priority) %>%
        group_by(peak_idx) %>%
        slice(1) %>%
        ungroup()
    
    annotated_peaks <- nrow(overlap_df_unique)
    non_annotated_peaks <- total_p2g_peaks - annotated_peaks
    
    # Build summary
    ccre_summary <- overlap_df_unique %>%
        count(ccre_class, name = "n_peaks") %>%
        mutate(pct_of_all_p2g = round(100 * n_peaks / total_p2g_peaks, 2))
    
    ccre_summary <- bind_rows(
        ccre_summary,
        tibble(
            ccre_class = "not_annotated",
            n_peaks = non_annotated_peaks,
            pct_of_all_p2g = round(100 * non_annotated_peaks / total_p2g_peaks, 2)
        )
    ) %>%
        arrange(desc(n_peaks))
    
    # Print and save
    message("\n[12.6] ENCODE cCRE Summary:")
    message(sprintf("  Annotated: %d (%.1f%%)", annotated_peaks, 100 * annotated_peaks / total_p2g_peaks))
    message(sprintf("  Not annotated: %d (%.1f%%)", non_annotated_peaks, 100 * non_annotated_peaks / total_p2g_peaks))
    print(as.data.frame(ccre_summary))
    
    write.csv(ccre_summary, "P2G_ENCODE_cCRE_annotation.csv", row.names = FALSE)
    message("[12.6] Saved: P2G_ENCODE_cCRE_annotation.csv")
}

# ==============================================================================
# 12.7: CHIPSEEKER ANNOTATION
# ==============================================================================

banner("12.7: ChIPseeker Annotation")

# Load TxDb
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

# Create GRanges for annotation
p2g_peaks_gr <- GRanges(
    seqnames = p2g_links$peak_chr,
    ranges = IRanges(start = p2g_links$peak_start, end = p2g_links$peak_end)
)
p2g_peaks_gr <- unique(p2g_peaks_gr)

# Annotate
peak_anno <- annotatePeak(
    peak = p2g_peaks_gr,
    TxDb = txdb,
    annoDb = "org.Mm.eg.db",
    tssRegion = c(-3000, 3000),
    verbose = FALSE
)

peak_anno_df <- as.data.frame(peak_anno)

# Hybrid annotation (collapse Intron/Exon/Downstream numbers)
peak_anno_df$chipseeker_annotation <- vapply(
    peak_anno_df$annotation,
    function(x) {
        x <- as.character(x)[1]
        if (grepl("^Intron", x)) {
            "Intron"
        } else if (grepl("^Exon", x)) {
            "Exon"
        } else if (grepl("^Downstream", x)) {
            "Downstream"
        } else {
            x
        }
    },
    character(1)
)

# Summary
chipseeker_summary <- as.data.frame(
    sort(table(peak_anno_df$chipseeker_annotation), decreasing = TRUE),
    stringsAsFactors = FALSE
)
colnames(chipseeker_summary) <- c("chipseeker_annotation", "n_peaks")

total_peaks <- sum(chipseeker_summary$n_peaks)
chipseeker_summary$pct_of_all_p2g <- round(100 * chipseeker_summary$n_peaks / total_peaks, 2)

# Print and save
message("\n[12.7] ChIPseeker Summary:")
print(chipseeker_summary)

write.csv(chipseeker_summary, "P2G_ChIPseeker_annotation.csv", row.names = FALSE)
message("[12.7] Saved: P2G_ChIPseeker_annotation.csv")

# ==============================================================================
# 12.8: SAVE RDS AND FINAL SUMMARY
# ==============================================================================

banner("12.8: Save RDS and Summary")

# Save full annotated data
saveRDS(p2g_links, "P2G_full_annotated.rds")
message("[12.8] Saved: P2G_full_annotated.rds")

# Peak type summary by cluster
peak_type_summary <- p2g_links %>%
    group_by(cluster, peakType) %>%
    summarise(n = n_distinct(peak), .groups = "drop") %>%
    pivot_wider(names_from = peakType, values_from = n, values_fill = 0)

write.csv(peak_type_summary, "P2G_cluster_peakType_summary.csv", row.names = FALSE)

# Top regulated genes
top_genes <- p2g_links %>%
    group_by(gene) %>%
    summarise(
        n_links = n(),
        mean_cor = mean(Correlation),
        .groups = "drop"
    ) %>%
    arrange(desc(n_links)) %>%
    head(20)

write.csv(top_genes, "P2G_top_regulated_genes.csv", row.names = FALSE)

# ==============================================================================
# 12.9: CHIPSEEKER STACKED BARPLOT (SINGLE CELL TYPE)
# ==============================================================================

banner("12.9: ChIPseeker Stacked Barplot")

# Lock annotation order (bottom → top)
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

# Publication-safe color palette
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

# Prepare data
chipseeker_plot_df <- chipseeker_summary %>%
    mutate(
        chipseeker_annotation = factor(chipseeker_annotation, levels = annotation_order)
    )

total_peaks_label <- sum(chipseeker_plot_df$n_peaks)

# Build barplot
p_chipseeker <- ggplot(
    chipseeker_plot_df,
    aes(x = "P2G Peaks", y = pct_of_all_p2g, fill = chipseeker_annotation)
) +
    geom_bar(stat = "identity", width = 0.6, color = "black", linewidth = 0.25) +
    scale_fill_manual(values = annotation_colors, drop = FALSE) +
    scale_y_continuous(
        breaks = seq(0, 100, 25),
        labels = function(x) paste0(x, "%"),
        expand = c(0, 0)
    ) +
    labs(
        x = sprintf("n = %s", format(total_peaks_label, big.mark = ",")),
        y = "% of P2G-linked peaks",
        fill = "ChIPseeker annotation"
    ) +
    theme_classic(base_size = 12) +
    theme(
        axis.text.x = element_text(size = 11),
        axis.line = element_line(color = "black"),
        legend.position = "right",
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 10)
    )

ggsave("P2G_ChIPseeker_stacked_barplot.pdf", p_chipseeker,
       width = 6, height = 5, useDingbats = FALSE)
message("[12.9] Saved: P2G_ChIPseeker_stacked_barplot.pdf")

# ==============================================================================
# 12.10: ENCODE cCRE STACKED BARPLOT (SINGLE CELL TYPE)
# ==============================================================================

banner("12.10: ENCODE cCRE Stacked Barplot")

if (exists("ccre_summary")) {
    
    # Lock cCRE class order (bottom → top)
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
    
    # Prepare data
    ccre_plot_df <- ccre_summary %>%
        mutate(ccre_class = factor(ccre_class, levels = ccre_class_order))
    
    total_ccre_peaks <- sum(ccre_plot_df$n_peaks)
    
    # Build barplot
    p_ccre <- ggplot(
        ccre_plot_df,
        aes(x = "P2G Peaks", y = pct_of_all_p2g, fill = ccre_class)
    ) +
        geom_bar(stat = "identity", width = 0.6, color = "black", linewidth = 0.25) +
        scale_fill_manual(values = ccre_colors, drop = FALSE) +
        scale_y_continuous(
            breaks = seq(0, 100, 25),
            labels = function(x) paste0(x, "%"),
            expand = c(0, 0)
        ) +
        labs(
            x = sprintf("n = %s", format(total_ccre_peaks, big.mark = ",")),
            y = "% of P2G-linked peaks",
            fill = "ENCODE cCRE class"
        ) +
        theme_classic(base_size = 12) +
        theme(
            axis.text.x = element_text(size = 11),
            axis.line = element_line(color = "black"),
            legend.position = "right",
            legend.title = element_text(size = 11),
            legend.text = element_text(size = 10)
        )
    
    ggsave("P2G_ENCODE_cCRE_stacked_barplot.pdf", p_ccre,
           width = 6, height = 5, useDingbats = FALSE)
    message("[12.10] Saved: P2G_ENCODE_cCRE_stacked_barplot.pdf")
    
} else {
    message("[12.10] Skipping ENCODE cCRE plot (annotation not performed)")
}

# ==============================================================================
# 12.11: CHIPSEEKER CROSS-CELLTYPE STACKED BARPLOT
# ==============================================================================

banner("12.11: ChIPseeker Cross-Celltype Stacked Barplot")

chipseeker_files <- list.files(
    pattern = "^P2G_ChIPseeker_.*\\.csv$",
    full.names = TRUE
)

if (length(chipseeker_files) < 2) {
    message("[12.11] Fewer than 2 ChIPseeker files found. Skipping cross-celltype plot.")
    message("[12.11] Run this script for each cell type first, then re-run this section.")
} else {
    
    chipseeker_cross_df <- purrr::map_dfr(chipseeker_files, function(f) {
        df <- read.csv(f, stringsAsFactors = FALSE)
        celltype <- gsub("^P2G_ChIPseeker_|\\.csv$", "", basename(f))
        df %>% mutate(celltype = celltype)
    })
    
    message(sprintf("[12.11] Loaded %d ChIPseeker files", length(chipseeker_files)))
    message(sprintf("[12.11] Cell types: %s", paste(unique(chipseeker_cross_df$celltype), collapse = ", ")))
    
    # Standardize labels + complete missing categories
    chipseeker_cross_df <- chipseeker_cross_df %>%
        mutate(
            chipseeker_annotation = gsub("_", " ", chipseeker_annotation),
            chipseeker_annotation = trimws(chipseeker_annotation)
        ) %>%
        tidyr::complete(
            celltype,
            chipseeker_annotation = annotation_order,
            fill = list(n_peaks = 0, pct_of_all_p2g = 0)
        )
    
    # Renormalize percentages per celltype
    chipseeker_cross_df <- chipseeker_cross_df %>%
        group_by(celltype) %>%
        mutate(
            pct_of_all_p2g = pct_of_all_p2g / sum(pct_of_all_p2g) * 100
        ) %>%
        ungroup()
    
    # Set factor levels
    chipseeker_cross_df$chipseeker_annotation <- factor(
        chipseeker_cross_df$chipseeker_annotation,
        levels = annotation_order
    )
    
    chipseeker_cross_df$celltype <- factor(
        chipseeker_cross_df$celltype,
        levels = intersect(CELL_ORDER, unique(chipseeker_cross_df$celltype))
    )
    
    # Compute total peaks per celltype for x-axis labels
    celltype_totals <- chipseeker_cross_df %>%
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
    
    # Build stacked barplot
    p_chipseeker_cross <- ggplot(
        chipseeker_cross_df,
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
    
    ggsave(
        "P2G_ChIPseeker_all_celltypes_stacked_barplot.pdf",
        p_chipseeker_cross,
        width = 9, height = 4.8,
        useDingbats = FALSE
    )
    message("[12.11] Saved: P2G_ChIPseeker_all_celltypes_stacked_barplot.pdf")
    
    write.csv(chipseeker_cross_df, "P2G_ChIPseeker_all_celltypes_combined.csv", row.names = FALSE)
    message("[12.11] Saved: P2G_ChIPseeker_all_celltypes_combined.csv")
}

# ==============================================================================
# 12.12: ENCODE cCRE CROSS-CELLTYPE STACKED BARPLOT
# ==============================================================================

banner("12.12: ENCODE cCRE Cross-Celltype Stacked Barplot")

ccre_files <- list.files(
    pattern = "^P2G_ENCODE_cCRE_.*\\.csv$",
    full.names = TRUE
)

if (length(ccre_files) < 2) {
    message("[12.12] Fewer than 2 ENCODE cCRE files found. Skipping cross-celltype plot.")
    message("[12.12] Run this script for each cell type first, then re-run this section.")
} else {
    
    ccre_cross_df <- purrr::map_dfr(ccre_files, function(f) {
        df <- read.csv(f, stringsAsFactors = FALSE)
        celltype <- gsub("^P2G_ENCODE_cCRE_|\\.csv$", "", basename(f))
        df %>% mutate(celltype = celltype)
    })
    
    message(sprintf("[12.12] Loaded %d ENCODE cCRE files", length(ccre_files)))
    message(sprintf("[12.12] Cell types: %s", paste(unique(ccre_cross_df$celltype), collapse = ", ")))
    
    # cCRE class order and colors (reuse from 12.10)
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
    
    # Set factor levels
    ccre_cross_df$ccre_class <- factor(ccre_cross_df$ccre_class, levels = ccre_class_order)
    
    ccre_cross_df$celltype <- factor(
        ccre_cross_df$celltype,
        levels = intersect(CELL_ORDER, unique(ccre_cross_df$celltype))
    )
    
    # Compute total peaks per celltype for x-axis labels
    celltype_totals <- ccre_cross_df %>%
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
    
    # Build stacked barplot
    p_ccre_cross <- ggplot(
        ccre_cross_df,
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
    
    ggsave(
        "P2G_ENCODE_cCRE_all_celltypes_stacked_barplot.pdf",
        p_ccre_cross,
        width = 9, height = 5,
        dpi = 300,
        useDingbats = FALSE
    )
    message("[12.12] Saved: P2G_ENCODE_cCRE_all_celltypes_stacked_barplot.pdf")
    
    write.csv(ccre_cross_df, "P2G_ENCODE_cCRE_all_celltypes_combined.csv", row.names = FALSE)
    message("[12.12] Saved: P2G_ENCODE_cCRE_all_celltypes_combined.csv")
}

# ==============================================================================
# FINAL SUMMARY
# ==============================================================================

banner("STEP 12 COMPLETE")

message(sprintf("Total P2G links: %d", nrow(p2g_links)))
message(sprintf("Unique peaks: %d", n_distinct(p2g_links$peak)))
message(sprintf("Unique genes: %d", n_distinct(p2g_links$gene)))
message(sprintf("K-means clusters: %d", length(unique(clusters))))

message("\nCluster Summary:")
print(cluster_summary)

message("\n=== OUTPUT FILES ===")
message("Figures:")
message(sprintf("  P2G_heatmap_k%d.pdf", K_CLUSTERS))
message("  P2G_ChIPseeker_stacked_barplot.pdf")
message("  P2G_ENCODE_cCRE_stacked_barplot.pdf")
message("  P2G_ChIPseeker_all_celltypes_stacked_barplot.pdf")
message("  P2G_ENCODE_cCRE_all_celltypes_stacked_barplot.pdf")
message("\nPer-cluster CSVs:")
message("  P2G_C[1-k]_genes.csv")
message("  P2G_C[1-k]_peaks_annotated.csv")
message("  P2G_C[1-k]_links.csv")
message("  P2G_C[1-k]_peaks.bed")
message("\nCombined Excel files:")
message("  P2G_AllClusters_GENES.xlsx")
message("  P2G_AllClusters_PEAKS.xlsx")
message("  P2G_AllClusters_LINKS.xlsx")
message("\nAnnotation summaries:")
message("  P2G_ENCODE_cCRE_annotation.csv")
message("  P2G_ChIPseeker_annotation.csv")
message("  P2G_ChIPseeker_all_celltypes_combined.csv")
message("  P2G_ENCODE_cCRE_all_celltypes_combined.csv")
message("\nSummary files:")
message("  P2G_cluster_summary.csv")
message("  P2G_cluster_peakType_summary.csv")
message("  P2G_top_regulated_genes.csv")
message("  P2G_all_peaks.bed")
message("  P2G_full_annotated.rds")

message("\n====================================================================")
