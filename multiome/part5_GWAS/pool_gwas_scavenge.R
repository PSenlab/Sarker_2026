#!/usr/bin/env Rscript
#===============================================================================
# POOL GWAS-SCAVENGE RESULTS ACROSS ALL CELL TYPES (Excluding Stellate)
#===============================================================================

library(data.table)
library(dplyr)
library(tidyr)

#===============================================================================
# CONFIGURATION
#===============================================================================
base_dir <- "/data/sarkern2/multiome_liver/Seurat/epigenome/gwas_peak_gene"
output_dir <- file.path(base_dir, "pooled_results")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Cell types (excluding Stellate - insufficient P2G data: only 169 peaks)
cell_types <- list(
  "Hepatocyte" = "hepatocyte",
  "Endothelial_01" = "endothelial_01",
  "Endothelial_02" = "endothelial_02",
  "Kupffer" = "kupffer",
  "MoMFs" = "momfs",
  "Cholangiocyte_01" = "cholangiocyte_01",
  "Cholangiocyte_02" = "cholangiocyte_02",
  "Lymp_B" = "lymp_b",
  "Lymp_T" = "lymp_t"
)

cat("===============================================================================\n")
cat("POOLING GWAS-SCAVENGE RESULTS (9 Cell Types, Excluding Stellate)\n")
cat("===============================================================================\n")
cat("Cell types:", length(cell_types), "\n")
cat("Output:", output_dir, "\n\n")

#===============================================================================
# 1. POOL FILTERED RESULTS (FDR < 0.05 cells)
#===============================================================================
cat("====== 1. Pooling Filtered Results ======\n")

filtered_list <- list()

for (cell_type in names(cell_types)) {
  dir_name <- cell_types[[cell_type]]
  results_dir <- file.path(base_dir, dir_name, "gwas_p2g/Liver/results")
  
  # Try different possible file naming patterns
  possible_files <- c(
    file.path(results_dir, paste0("combined_ALL_filtered_", cell_type, ".csv")),
    file.path(results_dir, paste0("combined_ALL_filtered_", dir_name, ".csv")),
    file.path(results_dir, "combined_ALL_filtered.csv")
  )
  
  filtered_file <- NULL
  for (f in possible_files) {
    if (file.exists(f)) {
      filtered_file <- f
      break
    }
  }
  
  if (!is.null(filtered_file)) {
    df <- fread(filtered_file)
    df$cell_type <- cell_type
    filtered_list[[cell_type]] <- df
    cat(sprintf("  %-20s: %7d cells (FDR < 0.05) | %3d traits\n", 
                cell_type, nrow(df), length(unique(df$trait))))
  } else {
    cat(sprintf("  %-20s: NOT FOUND\n", cell_type))
  }
}

# Combine all filtered results
if (length(filtered_list) > 0) {
  pooled_filtered <- rbindlist(filtered_list, fill = TRUE)
  fwrite(pooled_filtered, file.path(output_dir, "pooled_ALL_filtered.csv"))
  
  cat(sprintf("\n  TOTAL: %d cells across %d cell types, %d unique traits\n", 
              nrow(pooled_filtered), 
              length(unique(pooled_filtered$cell_type)),
              length(unique(pooled_filtered$trait))))
}

#===============================================================================
# 2. POOL SUMMARY FILES
#===============================================================================
cat("\n====== 2. Pooling Summary Files ======\n")

summary_list <- list()

for (cell_type in names(cell_types)) {
  dir_name <- cell_types[[cell_type]]
  results_dir <- file.path(base_dir, dir_name, "gwas_p2g/Liver/results")
  
  possible_files <- c(
    file.path(results_dir, paste0("summary_Liver_P2G_", cell_type, ".csv")),
    file.path(results_dir, paste0("summary_Liver_P2G_", dir_name, ".csv")),
    file.path(results_dir, "summary_Liver_P2G.csv")
  )
  
  summary_file <- NULL
  for (f in possible_files) {
    if (file.exists(f)) {
      summary_file <- f
      break
    }
  }
  
  if (!is.null(summary_file)) {
    df <- fread(summary_file)
    df$cell_type <- cell_type
    summary_list[[cell_type]] <- df
    
    # Handle Inf values in median_TRS
    med_trs <- median(df$median_TRS[is.finite(df$median_TRS)], na.rm = TRUE)
    
    cat(sprintf("  %-20s: %3d traits | median FDR<0.05: %5.1f%% | median TRS: %.4f\n", 
                cell_type, nrow(df), 
                median(df$pct_fdr_sig, na.rm = TRUE),
                med_trs))
  } else {
    cat(sprintf("  %-20s: NOT FOUND\n", cell_type))
  }
}

# Combine all summaries
if (length(summary_list) > 0) {
  pooled_summary <- rbindlist(summary_list, fill = TRUE)
  
  # Clean up Inf values
  pooled_summary[!is.finite(median_TRS), median_TRS := NA]
  pooled_summary[!is.finite(mean_TRS), mean_TRS := NA]
  
  fwrite(pooled_summary, file.path(output_dir, "pooled_summary.csv"))
  cat(sprintf("\n  TOTAL: %d trait-celltype combinations\n", nrow(pooled_summary)))
}

#===============================================================================
# 3. CREATE CROSS-CELL-TYPE MATRICES
#===============================================================================
cat("\n====== 3. Creating Cross-Cell-Type Matrices ======\n")

if (exists("pooled_summary") && nrow(pooled_summary) > 0) {
  
  # 3a. % FDR-significant cells matrix
  pct_fdr_matrix <- pooled_summary %>%
    select(trait, cell_type, pct_fdr_sig) %>%
    pivot_wider(names_from = cell_type, values_from = pct_fdr_sig, values_fill = 0) %>%
    as.data.frame()
  fwrite(pct_fdr_matrix, file.path(output_dir, "matrix_pct_fdr_sig.csv"))
  cat(sprintf("  pct_fdr_sig matrix: %d traits x %d cell types\n", nrow(pct_fdr_matrix), ncol(pct_fdr_matrix) - 1))
  
  # 3b. Median TRS matrix
  trs_matrix <- pooled_summary %>%
    select(trait, cell_type, median_TRS) %>%
    pivot_wider(names_from = cell_type, values_from = median_TRS, values_fill = NA) %>%
    as.data.frame()
  fwrite(trs_matrix, file.path(output_dir, "matrix_median_TRS.csv"))
  cat(sprintf("  median_TRS matrix: %d traits x %d cell types\n", nrow(trs_matrix), ncol(trs_matrix) - 1))
  
  # 3c. N FDR-significant cells matrix
  n_fdr_matrix <- pooled_summary %>%
    select(trait, cell_type, n_fdr_sig) %>%
    pivot_wider(names_from = cell_type, values_from = n_fdr_sig, values_fill = 0) %>%
    as.data.frame()
  fwrite(n_fdr_matrix, file.path(output_dir, "matrix_n_fdr_sig.csv"))
  cat(sprintf("  n_fdr_sig matrix: %d traits x %d cell types\n", nrow(n_fdr_matrix), ncol(n_fdr_matrix) - 1))
}

#===============================================================================
# 4. SUMMARY BY CELL TYPE
#===============================================================================
cat("\n====== 4. Summary by Cell Type ======\n")

if (exists("pooled_summary") && nrow(pooled_summary) > 0) {
  
  celltype_stats <- pooled_summary %>%
    group_by(cell_type) %>%
    summarise(
      n_traits = n(),
      total_sig_cells = sum(n_fdr_sig, na.rm = TRUE),
      mean_pct_fdr = mean(pct_fdr_sig, na.rm = TRUE),
      median_pct_fdr = median(pct_fdr_sig, na.rm = TRUE),
      mean_TRS = mean(median_TRS, na.rm = TRUE),
      median_TRS = median(median_TRS, na.rm = TRUE),
      traits_gt5pct = sum(pct_fdr_sig > 5, na.rm = TRUE),
      traits_gt10pct = sum(pct_fdr_sig > 10, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(total_sig_cells))
  
  fwrite(celltype_stats, file.path(output_dir, "summary_by_celltype.csv"))
  
  cat("\n")
  print(as.data.frame(celltype_stats), row.names = FALSE)
}

#===============================================================================
# 5. SUMMARY BY TRAIT
#===============================================================================
cat("\n====== 5. Summary by Trait ======\n")

if (exists("pooled_summary") && nrow(pooled_summary) > 0) {
  
  trait_stats <- pooled_summary %>%
    filter(!is.na(trait)) %>%
    group_by(trait) %>%
    reframe(
      n_celltypes = n(),
      total_sig_cells = sum(n_fdr_sig, na.rm = TRUE),
      mean_pct_fdr = mean(pct_fdr_sig, na.rm = TRUE),
      max_pct_fdr = max(pct_fdr_sig, na.rm = TRUE),
      mean_TRS = mean(median_TRS, na.rm = TRUE),
      max_TRS = max(median_TRS, na.rm = TRUE),
      best_celltype_pct = cell_type[which.max(pct_fdr_sig)][1],
      best_celltype_TRS = cell_type[which.max(median_TRS)][1]
    ) %>%
    arrange(desc(total_sig_cells))
  
  fwrite(trait_stats, file.path(output_dir, "summary_by_trait.csv"))
  
  cat("\nTop 20 Traits by Total Significant Cells:\n")
  top20 <- head(trait_stats, 20) %>%
    select(trait, n_celltypes, total_sig_cells, max_pct_fdr, best_celltype_pct)
  print(as.data.frame(top20), row.names = FALSE)
}

#===============================================================================
# 6. CELL-TYPE SPECIFICITY
#===============================================================================
cat("\n====== 6. Cell-Type Specificity ======\n")

if (exists("pooled_summary") && nrow(pooled_summary) > 0) {
  
  specificity <- pooled_summary %>%
    filter(!is.na(trait)) %>%
    group_by(trait) %>%
    mutate(
      mean_pct_across = mean(pct_fdr_sig, na.rm = TRUE),
      specificity_score = ifelse(mean_pct_across > 0, pct_fdr_sig / mean_pct_across, 0)
    ) %>%
    ungroup() %>%
    select(trait, cell_type, pct_fdr_sig, mean_pct_across, specificity_score, median_TRS) %>%
    arrange(desc(specificity_score))
  
  fwrite(specificity, file.path(output_dir, "trait_celltype_specificity.csv"))
  
  # Top specific per cell type
  cat("\nTop 3 Most Specific Traits per Cell Type (min 5% FDR-sig):\n\n")
  
  top_per_ct <- specificity %>%
    filter(pct_fdr_sig >= 5) %>%
    group_by(cell_type) %>%
    slice_max(specificity_score, n = 3) %>%
    select(cell_type, trait, pct_fdr_sig, specificity_score) %>%
    arrange(cell_type, desc(specificity_score))
  
  fwrite(top_per_ct, file.path(output_dir, "top_specific_traits_per_celltype.csv"))
  print(as.data.frame(top_per_ct), row.names = FALSE)
}

#===============================================================================
# 7. QUICK STATS
#===============================================================================
cat("\n===============================================================================\n")
cat("POOLING SUMMARY\n")
cat("===============================================================================\n")

if (exists("pooled_filtered")) {
  cat(sprintf("Total FDR<0.05 cells: %s\n", format(nrow(pooled_filtered), big.mark = ",")))
}
if (exists("pooled_summary")) {
  cat(sprintf("Total trait-celltype combinations: %d\n", nrow(pooled_summary)))
  cat(sprintf("Unique traits: %d\n", length(unique(pooled_summary$trait[!is.na(pooled_summary$trait)]))))
}
cat(sprintf("Cell types: %d\n", length(cell_types)))

cat("\n===============================================================================\n")
cat("OUTPUT FILES\n")
cat("===============================================================================\n")
cat(sprintf("Directory: %s\n\n", output_dir))
cat("  pooled_ALL_filtered.csv           - All FDR<0.05 cells\n")
cat("  pooled_summary.csv                - All trait summaries\n")
cat("  matrix_pct_fdr_sig.csv            - Traits x CellTypes (% FDR-sig)\n")
cat("  matrix_median_TRS.csv             - Traits x CellTypes (median TRS)\n")
cat("  matrix_n_fdr_sig.csv              - Traits x CellTypes (N cells)\n")
cat("  summary_by_celltype.csv           - Stats per cell type\n")
cat("  summary_by_trait.csv              - Stats per trait\n")
cat("  trait_celltype_specificity.csv    - Specificity scores\n")
cat("  top_specific_traits_per_celltype.csv - Top specific traits\n")

cat("\n===============================================================================\n")
cat("COMPLETE\n")
cat("===============================================================================\n")

#!/usr/bin/env Rscript
#===============================================================================
# DOTPLOT: Median TRS (color) and % FDR-sig (size) by Cell Type
#===============================================================================

library(data.table)
library(ggplot2)
library(dplyr)
library(tidyr)

#===============================================================================
# LOAD DATA
#===============================================================================
output_dir <- "/data/sarkern2/multiome_liver/Seurat/epigenome/gwas_peak_gene/pooled_results"

# Load pooled summary (has both median_TRS and pct_fdr_sig)
summary_df <- fread(file.path(output_dir, "pooled_summary.csv"))

cat("Loaded:", nrow(summary_df), "trait-celltype combinations\n")
cat("Unique traits:", length(unique(summary_df$trait)), "\n")
cat("Cell types:", length(unique(summary_df$cell_type)), "\n")

#===============================================================================
# DEFINE TRAIT CATEGORIES
#===============================================================================
category_patterns <- list(
  "Liver Enzymes" = c(
    "\\bALT\\b", "\\bAST\\b", "\\bGGT\\b", "\\bALP\\b",
    "alanine aminotransferase", "alanine transaminase",
    "aspartate aminotransferase", "aspartate transaminase",
    "gamma.glutamyl", "gamma glutamyl", "gamma-glutamyl",
    "alkaline phosphatase"
  ),
  "Bilirubin" = c("bilirubin"),
  "Albumin/Globulins" = c("albumin", "globulin", "non-albumin protein"),
  "Complement/Coagulation" = c(
    "complement C", "complement component",
    "fibrinogen", "prothrombin", "coagulation factor",
    "factor VII", "factor VIII", "factor IX", "factor X"
  ),
  "Apolipoproteins" = c("apolipoprotein", "\\bAPO[ABCDEFHLM]\\b", "\\bAPOM\\b"),
  "VLDL/Lipoproteins" = c("\\bVLDL\\b", "lipoprotein"),
  "Iron Metabolism" = c("ferritin", "transferrin", "hemochromatosis", "hepcidin", "iron\\b", "\\biron"),
  "Other Liver Proteins" = c(
    "ceruloplasmin", "haptoglobin", "alpha.1.antitrypsin", "alpha-1-antitrypsin",
    "angiotensinogen", "protein quantitative trait"
  ),
  "NAFLD/Steatosis" = c("NAFLD", "NASH", "steatosis", "steatohepatitis", "fatty liver", "hepatic fat", "liver fat"),
  "Fibrosis/Cirrhosis" = c("fibrosis", "cirrhosis", "liver condition"),
  "Liver Cancer" = c("hepatocellular carcinoma", "hepatic cancer", "liver cancer", "cholangiocarcinoma"),
  "Viral Hepatitis" = c("hepatitis B", "hepatitis C", "hepatitis E", "\\bHBV\\b", "\\bHCV\\b", "anti-hepatitis"),
  "Autoimmune Liver" = c("autoimmune hepatitis", "primary biliary", "primary sclerosing", "\\bPBC\\b", "\\bPSC\\b"),
  "Biliary/Cholestasis" = c("biliary", "bile ", "cholestasis", "cholangitis", "gallbladder"),
  "Drug-Induced" = c("drug.induced", "hepatotoxicity", "alcohol"),
  "Pregnancy-Related" = c("pregnancy", "childbirth", "puerperium", "neonatal"),
  "Metabolic Syndrome" = c("metabolic syndrome", "obesity", "BMI", "diabetes")
)

categorize_trait <- function(trait, patterns_list) {
  if (is.na(trait) || trait == "") return("Other")
  trait_lower <- tolower(trait)
  for (category in names(patterns_list)) {
    for (pattern in patterns_list[[category]]) {
      if (grepl(tolower(pattern), trait_lower, ignore.case = TRUE)) return(category)
    }
  }
  return("Other")
}

#===============================================================================
# PREPARE DATA
#===============================================================================
# Remove NA traits
df <- summary_df[!is.na(trait) & trait != ""]

# Add categories
df$category <- sapply(df$trait, function(t) categorize_trait(t, category_patterns))

# Handle Inf/NA in median_TRS
df$median_TRS[!is.finite(df$median_TRS)] <- NA

# Category order
category_order <- c("Liver Enzymes", "Bilirubin", "Albumin/Globulins",
                    "Complement/Coagulation", "Apolipoproteins", "VLDL/Lipoproteins",
                    "Iron Metabolism", "Other Liver Proteins", "NAFLD/Steatosis",
                    "Fibrosis/Cirrhosis", "Liver Cancer", "Viral Hepatitis",
                    "Autoimmune Liver", "Biliary/Cholestasis", "Drug-Induced",
                    "Pregnancy-Related", "Metabolic Syndrome", "Other")

df$category <- factor(df$category, levels = category_order)

# Cell type order
celltype_order <- c("Hepatocyte", "Cholangiocyte_01", "Cholangiocyte_02",
                    "Endothelial_01", "Endothelial_02", "Kupffer", "MoMFs",
                    "Lymp_B", "Lymp_T")
df$cell_type <- factor(df$cell_type, levels = celltype_order)

# Order traits by category, then by max TRS
trait_order <- df %>%
  group_by(trait, category) %>%
  summarise(max_trs = max(median_TRS, na.rm = TRUE), .groups = "drop") %>%
  arrange(category, desc(max_trs)) %>%
  pull(trait)

df$trait <- factor(df$trait, levels = rev(trait_order))  # Reverse for y-axis

cat("Traits:", length(unique(df$trait)), "\n")
cat("Categories:\n")
print(table(df$category))#===============================================================================
# DOTPLOT - ALL TRAITS
#===============================================================================
cat("\nGenerating dotplot...\n")

n_traits <- length(unique(df$trait))
plot_height <- max(15, n_traits * 0.15)

p <- ggplot(df, aes(x = cell_type, y = trait)) +
  geom_point(aes(size = pct_fdr_sig, color = median_TRS)) +
  
  # Color scale - red gradient (low to high)
  scale_color_gradient(
    low = "white", 
    high = "#CC0000",
    na.value = "grey80",
    name = "Median TRS",
    limits = c(0, max(df$median_TRS, na.rm = TRUE))
  ) +
  
  # Size scale
  scale_size_continuous(
    name = "% FDR < 0.05",
    range = c(0.5, 6),
    breaks = c(0, 2, 5, 10, 15),
    limits = c(0, max(df$pct_fdr_sig, na.rm = TRUE))
  ) +
  
  # Facet by category
  facet_grid(category ~ ., scales = "free_y", space = "free_y") +
  
  # Theme
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 5),
    axis.title = element_blank(),
    strip.text.y = element_text(angle = 0, size = 7, face = "bold"),
    strip.background = element_rect(fill = "grey90"),
    panel.grid.major = element_line(color = "grey90", size = 0.3),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 8),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
  ) +
  
  labs(title = "GWAS-SCAVENGE: Liver Trait Enrichment by Cell Type")

# Save
ggsave(file.path(output_dir, "dotplot_TRS_pctFDR_all_traits.pdf"),
       p, width = 8, height = plot_height, limitsize = FALSE)


cat("Saved: dotplot_TRS_pctFDR_all_traits.pdf/png\n")
cat("Dimensions: 12 x", plot_height, "inches\n")



#===============================================================================
# SUMMARY
#===============================================================================
cat("\n===============================================================================\n")
cat("DOTPLOTS CREATED\n")
cat("===============================================================================\n")
cat("1. dotplot_TRS_pctFDR_all_traits.pdf/png       - All traits\n")
cat("2. dotplot_TRS_pctFDR_top5_per_category.pdf/png - Top 5 per category\n")
cat("\nDot color = Median TRS (white → red)\n")
cat("Dot size  = % cells with FDR < 0.05\n")
cat("\nOutput:", output_dir, "\n")
cat("===============================================================================\n")


