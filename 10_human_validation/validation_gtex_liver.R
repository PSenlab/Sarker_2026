#!/usr/bin/env Rscript
# ==============================================================================
# GTEx Liver Validation - ASCL1 Male vs Female
# ==============================================================================
#
# Description:
#   Validates Ascl1 sex dimorphism in the GTEx v10 liver cohort using
#   raw gene-level counts and DESeq2 median-of-ratios normalization.
#   Sex (1 = male, 2 = female) and age group are annotated in the GTEx
#   subject phenotype file. Subjects are restricted to death-circumstance
#   Hardy scales 0-3 (sudden or natural deaths) which are standard GTEx
#   quality-filter strata for RNA-seq expression analyses.
#
# Dataset-specific quirks:
#   - GCT counts file (GTEx v10 liver)
#   - Sample ID rewrite: dotted "GTEX.XXXXX.YYYY..." -> "GTEX-XXXXX"
#     (first two hyphen-delimited tokens only, since the phenotype table
#     keys by SUBJID not SAMPID)
#   - DESeq2 normalized counts (size-factor scaled) -> log2(x + 1)
#   - Hardy scale filter: DTHHRDY in {0, 1, 2, 3}
#   - Low-count filter: rowSums > 10
#   - ASCL1 looked up by versioned Ensembl ID (ENSG00000139352.4)
#
# Inputs:
#   GTEx_Analysis_v10_Annotations_SubjectPhenotypesDS.txt
#   gene_reads_v10_liver.gct
#
# Output:
#   verify_GTEx_liver_ASCL1_barplot.pdf
#
# ==============================================================================

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(dplyr)
})


# ==============================================================================
# CONFIGURATION - UPDATE THESE PATHS
# ==============================================================================

PHENOTYPE_FILE <- "GTEx_Analysis_v10_Annotations_SubjectPhenotypesDS.txt"
COUNTS_FILE    <- "gene_reads_v10_liver.gct"

ASCL1_ENSG     <- "ENSG00000139352.4"
HARDY_KEEP     <- c(0, 1, 2, 3)     # death-circumstance filter
MIN_ROW_COUNT  <- 10                # filter: rowSums > MIN_ROW_COUNT

# Plot palette - matches the rest of the validation folder
SEX_COLORS     <- c("Male" = "#3498db", "Female" = "#e74c3c")

OUTPUT_PDF     <- "verify_GTEx_liver_ASCL1_barplot.pdf"


# ==============================================================================
# HELPERS
# ==============================================================================

banner <- function(text) {
  line <- paste(rep("=", 70), collapse = "")
  message("\n", line)
  message(text)
  message(line)
}

sig_label <- function(p) {
  if (p < 0.001) return("***")
  if (p < 0.01)  return("**")
  if (p < 0.05)  return("*")
  return("n.s.")
}

# Strip trailing sample-run tokens so that column IDs match GTEx SUBJID
clean_sample_ids <- function(x) {
  parts <- strsplit(gsub("\\.", "-", x), "-", fixed = FALSE)
  sapply(parts, function(p) paste(p[1:2], collapse = "-"))
}


# ==============================================================================
# STEP 1: Load phenotypes and counts
# ==============================================================================
banner("STEP 1: Load phenotypes + counts (GCT)")

phenotypes <- read.table(PHENOTYPE_FILE, sep = "\t", header = TRUE)
phenotypes_sub <- subset(phenotypes, DTHHRDY %in% HARDY_KEEP)
message(sprintf("  Phenotype table: %d subjects (Hardy 0-3 subset: %d)",
                nrow(phenotypes), nrow(phenotypes_sub)))

counts_raw <- read.table(COUNTS_FILE, skip = 2, header = TRUE)
rownames(counts_raw) <- counts_raw$Name
counts_raw$Name        <- NULL
counts_raw$Description <- NULL
message(sprintf("  Counts matrix:   %d genes x %d samples",
                nrow(counts_raw), ncol(counts_raw)))


# ==============================================================================
# STEP 2: Clean sample IDs and intersect with phenotype SUBJIDs
# ==============================================================================
banner("STEP 2: Sample ID harmonization")

colnames(counts_raw) <- clean_sample_ids(colnames(counts_raw))
common <- intersect(colnames(counts_raw), phenotypes_sub$SUBJID)
message(sprintf("  Common samples (counts AND Hardy 0-3 phenotype): %d",
                length(common)))

counts_matched <- counts_raw[, common, drop = FALSE]

# Drop low-count genes
counts_matched <- counts_matched[rowSums(counts_matched) > MIN_ROW_COUNT, , drop = FALSE]
message(sprintf("  Genes after filter (rowSums > %d): %d",
                MIN_ROW_COUNT, nrow(counts_matched)))

# Align phenotypes to counts column order
pheno_matched <- phenotypes_sub[phenotypes_sub$SUBJID %in% colnames(counts_matched), ]
pheno_matched <- pheno_matched[match(colnames(counts_matched), pheno_matched$SUBJID), ]
stopifnot(all(pheno_matched$SUBJID == colnames(counts_matched)))
pheno_matched$AGE <- as.factor(as.character(pheno_matched$AGE))


# ==============================================================================
# STEP 3: DESeq2 normalization
# ==============================================================================
banner("STEP 3: DESeq2 median-of-ratios normalization")

dds <- DESeqDataSetFromMatrix(
  countData = counts_matched,
  colData   = pheno_matched,
  design    = ~ AGE
)
dds <- DESeq(dds)
norm_counts <- counts(dds, normalized = TRUE)
message(sprintf("  [OK] Normalized counts: %d genes x %d samples",
                nrow(norm_counts), ncol(norm_counts)))


# ==============================================================================
# STEP 4: ASCL1 extraction
# ==============================================================================
banner("STEP 4: ASCL1 extraction")

if (!(ASCL1_ENSG %in% rownames(norm_counts))) {
  stop(sprintf("ASCL1 (%s) not found in normalized count matrix", ASCL1_ENSG))
}

ascl1_expr <- as.numeric(norm_counts[ASCL1_ENSG, ])

plot_df <- data.frame(
  SUBJID     = pheno_matched$SUBJID,
  AGE        = pheno_matched$AGE,
  SEX        = factor(pheno_matched$SEX, levels = c(1, 2),
                       labels = c("Male", "Female")),
  DTHHRDY    = pheno_matched$DTHHRDY,
  Expression = ascl1_expr,
  log2_expr  = log2(ascl1_expr + 1)
)

message(sprintf("  log2(norm+1) range: %.3f - %.3f",
                min(plot_df$log2_expr), max(plot_df$log2_expr)))
message(sprintf("  Non-zero ASCL1: %d / %d",
                sum(plot_df$Expression > 0), nrow(plot_df)))


# ==============================================================================
# STEP 5: Wilcoxon rank-sum (Male vs Female)
# ==============================================================================
banner("STEP 5: Wilcoxon rank-sum test (Male vs Female)")

sex_counts <- table(plot_df$SEX)
n_male   <- sex_counts["Male"]
n_female <- sex_counts["Female"]

wilcox_res <- wilcox.test(log2_expr ~ SEX, data = plot_df)
p_val <- wilcox_res$p.value
sig   <- sig_label(p_val)

mean_m <- mean(plot_df$log2_expr[plot_df$SEX == "Male"])
mean_f <- mean(plot_df$log2_expr[plot_df$SEX == "Female"])
sd_m   <- sd(plot_df$log2_expr[plot_df$SEX == "Male"])
sd_f   <- sd(plot_df$log2_expr[plot_df$SEX == "Female"])

message(sprintf("  Male:   n=%d, mean=%.4f +/- %.4f", n_male,   mean_m, sd_m))
message(sprintf("  Female: n=%d, mean=%.4f +/- %.4f", n_female, mean_f, sd_f))
message(sprintf("  W = %.1f, p = %.4e %s",
                wilcox_res$statistic, p_val, sig))


# ==============================================================================
# STEP 6: Bar + jitter + significance bracket
# ==============================================================================
banner("STEP 6: Generating barplot")

summary_df <- plot_df %>%
  group_by(SEX) %>%
  summarise(
    mean_expr = mean(log2_expr),
    sem       = sd(log2_expr) / sqrt(n()),
    .groups   = "drop"
  )

# Dynamic bracket heights
y_data_max <- max(plot_df$log2_expr)
y_sig_line <- y_data_max * 1.08
y_sig_text <- y_data_max * 1.12
y_axis_max <- y_data_max * 1.18

p <- ggplot(summary_df, aes(x = SEX, y = mean_expr, fill = SEX)) +
  geom_bar(stat = "identity", width = 0.55,
           color = "black", linewidth = 0.6) +
  geom_errorbar(aes(ymin = mean_expr - sem, ymax = mean_expr + sem),
                width = 0.18, linewidth = 0.6) +
  geom_jitter(
    data = plot_df,
    aes(x = SEX, y = log2_expr),
    width = 0.12, alpha = 0.35, size = 1.2,
    color = "black", inherit.aes = FALSE
  ) +
  annotate("segment",
           x = 1, xend = 2,
           y = y_sig_line, yend = y_sig_line,
           linewidth = 0.8) +
  annotate("text",
           x = 1.5, y = y_sig_text,
           label = sig, size = 5.5, fontface = "bold") +
  scale_fill_manual(values = SEX_COLORS) +
  scale_y_continuous(limits = c(0, y_axis_max), expand = c(0, 0)) +
  scale_x_discrete(
    labels = c(
      sprintf("Male\n(n=%d)",   n_male),
      sprintf("Female\n(n=%d)", n_female)
    )
  ) +
  labs(
    title    = "GTEx liver - ASCL1: Male vs Female",
    subtitle = sprintf("Wilcoxon p = %.2e %s", p_val, sig),
    x = NULL,
    y = "ASCL1 expression\n(log2(DESeq2-normalized counts))"
  ) +
  theme_classic(base_size = 13) +
  theme(
    text            = element_text(family = "Arial"),
    plot.title      = element_text(face = "bold", hjust = 0.5),
    plot.subtitle   = element_text(hjust = 0.5, color = "gray30", size = 10),
    legend.position = "none"
  )

pdf(OUTPUT_PDF, width = 5, height = 5.5, useDingbats = FALSE)
print(p)
dev.off()
message(sprintf("  [OK] %s", OUTPUT_PDF))


# ==============================================================================
# SUMMARY
# ==============================================================================
banner("PIPELINE SUMMARY")
message(sprintf("  Platform:       GTEx v10 bulk RNA-seq (liver)"))
message(sprintf("  Gene IDs:       Versioned Ensembl (ENSG...) from GCT"))
message(sprintf("  Normalization:  DESeq2 median-of-ratios -> log2(x + 1)"))
message(sprintf("  Sex:            ANNOTATED (GTEx subject phenotype, 1=Male/2=Female)"))
message(sprintf("  Hardy filter:   DTHHRDY in {%s}",
                paste(HARDY_KEEP, collapse = ", ")))
message(sprintf("  Samples:        %d total (%dM / %dF)",
                nrow(plot_df), n_male, n_female))
message(sprintf("  ASCL1 result:   Wilcoxon, p = %.4e %s", p_val, sig))

sessionInfo()
