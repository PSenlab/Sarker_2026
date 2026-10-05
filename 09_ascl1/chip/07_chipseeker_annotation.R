#!/usr/bin/env Rscript
# ------------------------------------------------------------------
# ChIPseeker annotation: ASCL1 peaks, male vs female
# Figure 7, panels a (annotation pies) and d (genes by sex
# specificity)
#
# Input : Ascl1_M_narrowpeaks.bed, Ascl1_F_narrowpeaks.bed
#         BED3, tab-delimited, WITH a header row (chr / start / end)
# Output (chipseeker_out/):
#   annotation_distribution.pdf   bar, TSS distance, pies (panel a)
#   sex_specificity.pdf           genes by sex specificity (panel d)
#   annotated_peaks.xlsx          per-sex, merged, promoter,
#                                 sex-specificity, summary
#   all_peaks.csv                 merged annotation, input for
#                                 08_auc_by_annotation.py and
#                                 09_genelist_peak_overlap.py
#   sex_specificity.csv
#   peakAnnoList.rds, specificity.rds
# ------------------------------------------------------------------

suppressPackageStartupMessages({
  library(ChIPseeker)
  library(GenomicRanges)
  library(GenomeInfoDb)
  library(openxlsx)
  library(ggplot2)
  library(systemfonts)
})

# ----------------------- CONFIG -----------------------------------
GENOME <- "mm10"   # "mm10" or "hg38"

if (GENOME == "mm10") {
  library(TxDb.Mmusculus.UCSC.mm10.knownGene)
  library(org.Mm.eg.db)
  txdb   <- TxDb.Mmusculus.UCSC.mm10.knownGene
  annoDb <- "org.Mm.eg.db"
  n_auto <- 19
} else {
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(org.Hs.eg.db)
  txdb   <- TxDb.Hsapiens.UCSC.hg38.knownGene
  annoDb <- "org.Hs.eg.db"
  n_auto <- 22
}

std_chroms <- paste0("chr", c(seq_len(n_auto), "X", "Y"))

# List order defines plot order and worksheet order downstream
peakfiles <- list(
  male   = "Ascl1_M_narrowpeaks.bed",
  female = "Ascl1_F_narrowpeaks.bed"
)
sex_levels <- names(peakfiles)

TSS_WINDOW <- c(-3000, 3000)   # promoter definition
BED_0BASED <- TRUE             # TRUE if written by a peak caller (BED convention)
OUTDIR     <- "chipseeker_out"
dir.create(OUTDIR, showWarnings = FALSE)

# ----------------------- FONT RESOLUTION ---------------------------
# pdf() uses the PostScript font database, which has no Arial. cairo_pdf()
# goes through the system font stack. Even so, Arial is often absent on
# Linux nodes -- fall back to a metric-compatible substitute rather than
# letting Cairo silently pick something wider.
resolve_font <- function() {
  cands <- c("Arial", "Liberation Sans", "Nimbus Sans", "Helvetica")
  for (f in cands) {
    p <- tryCatch(systemfonts::match_font(f)$path, error = function(e) "")
    if (nzchar(p) &&
        grepl(gsub(" ", "", f), gsub(" ", "", basename(p)), ignore.case = TRUE)) {
      return(f)
    }
  }
  "sans"
}
FONT <- resolve_font()
message("Using font family: ", FONT)
if (FONT != "Arial") {
  message("  NOTE: Arial not found. '", FONT, "' substituted. ",
          "Liberation Sans is metrically identical to Arial; others are not.")
}

theme_pub <- theme(text = element_text(family = FONT))

# ----------------------- 1. READ ----------------------------------
# readPeakFile() cannot handle a header row -- it evaluates peak.df[,2] + 1
# on the literal string "start". Read manually.
read_bed3 <- function(f) {
  if (!file.exists(f)) stop("File not found: ", f)
  df <- read.table(f, header = TRUE, sep = "\t",
                   stringsAsFactors = FALSE,
                   col.names = c("chr", "start", "end"))
  if (!is.numeric(df$start) || !is.numeric(df$end)) {
    stop("start/end not numeric in ", f,
         " -- check for extra header or track lines")
  }
  makeGRangesFromDataFrame(
    df,
    seqnames.field          = "chr",
    start.field             = "start",
    end.field               = "end",
    starts.in.df.are.0based = BED_0BASED
  )
}

peaks <- lapply(peakfiles, read_bed3)

# Harmonize chromosome naming, drop scaffolds and chrM
peaks <- lapply(peaks, function(gr) {
  if (!any(grepl("^chr", as.character(seqnames(gr))))) {
    seqlevelsStyle(gr) <- "UCSC"
  }
  gr <- gr[as.character(seqnames(gr)) %in% std_chroms]
  seqlevels(gr) <- seqlevelsInUse(gr)
  sort(gr)
})

message("Peak counts after filtering:")
for (nm in sex_levels) {
  message("  ", nm, ": ", length(peaks[[nm]]),
          "  (median width ", median(width(peaks[[nm]])), " bp)")
}

# ----------------------- 2. ANNOTATE -------------------------------
peakAnnoList <- lapply(peaks, annotatePeak,
                       TxDb      = txdb,
                       tssRegion = TSS_WINDOW,
                       annoDb    = annoDb,
                       verbose   = FALSE)

anno_dfs <- lapply(peakAnnoList, as.data.frame)

# ----------------------- 3. DISTRIBUTION PLOTS ---------------------
# cairo_pdf, not pdf(). onefile=TRUE or you get one file per page.
cairo_pdf(file.path(OUTDIR, "annotation_distribution.pdf"),
          width = 9, height = 5, onefile = TRUE)

# ggplot-based: theme() works
print(plotAnnoBar(peakAnnoList,
                  title = "Ascl1 peak feature distribution") + theme_pub)

print(plotDistToTSS(peakAnnoList,
                    title = "Ascl1 peak distribution relative to TSS") + theme_pub)

# plotAnnoPie is base graphics -- theme() is ignored, set via par()
old_par <- par(no.readonly = TRUE)
for (nm in sex_levels) {
  par(family = FONT)
  plotAnnoPie(peakAnnoList[[nm]], main = nm)
}
par(old_par)

dev.off()

# ----------------------- 4. SEX-SPECIFICITY + VENN -----------------
male_g   <- unique(na.omit(anno_dfs$male$SYMBOL))
female_g <- unique(na.omit(anno_dfs$female$SYMBOL))
m_only <- setdiff(male_g, female_g)
f_only <- setdiff(female_g, male_g)
shared <- intersect(male_g, female_g)

specificity <- data.frame(
  SYMBOL = c(m_only, f_only, shared),
  class  = c(rep("male_specific",   length(m_only)),
             rep("female_specific", length(f_only)),
             rep("shared",          length(shared))),
  stringsAsFactors = FALSE
)

# ChIPseeker's vennplot() is grid-based and ignores both theme() and par(),
# so its labels would render in the default font. Draw it as a ggplot instead.
venn_df <- data.frame(
  class = factor(c("male_specific", "shared", "female_specific"),
                 levels = c("male_specific", "shared", "female_specific")),
  n     = c(length(m_only), length(shared), length(f_only))
)

p_venn <- ggplot(venn_df, aes(x = class, y = n, fill = class)) +
  geom_col(width = 0.65) +
  geom_text(aes(label = n), vjust = -0.4, family = FONT) +
  scale_fill_manual(values = c(male_specific   = "#89d9e1",
                               shared          = "grey85",
                               female_specific = "#fb9a99")) +
  labs(title = "Annotated genes by sex specificity",
       x = NULL, y = "Number of genes") +
  theme_classic() +
  theme_pub +
  theme(legend.position = "none")

cairo_pdf(file.path(OUTDIR, "sex_specificity.pdf"),
          width = 5, height = 4.5, onefile = TRUE)
print(p_venn)
dev.off()

# ----------------------- 5. EXCEL ----------------------------------
wb <- createWorkbook()

for (nm in sex_levels) {
  addWorksheet(wb, nm)
  writeData(wb, nm, anno_dfs[[nm]])
  freezePane(wb, nm, firstRow = TRUE)
}

merged <- do.call(rbind, lapply(sex_levels, function(nm) {
  df <- anno_dfs[[nm]]
  df$sex <- nm
  df
}))
merged$sex <- factor(merged$sex, levels = sex_levels)

addWorksheet(wb, "All_peaks")
writeData(wb, "All_peaks", merged)
freezePane(wb, "All_peaks", firstRow = TRUE)

prom <- merged[grepl("^Promoter", merged$annotation), ]
addWorksheet(wb, "Promoter_only")
writeData(wb, "Promoter_only", prom)
freezePane(wb, "Promoter_only", firstRow = TRUE)

addWorksheet(wb, "Sex_specificity")
writeData(wb, "Sex_specificity", specificity)
freezePane(wb, "Sex_specificity", firstRow = TRUE)

summary_df <- data.frame(
  sex          = sex_levels,
  n_peaks      = vapply(peaks[sex_levels], length, integer(1)),
  median_width = vapply(peaks[sex_levels], function(g) median(width(g)), numeric(1)),
  n_genes      = vapply(anno_dfs[sex_levels],
                        function(d) length(unique(na.omit(d$geneId))), integer(1)),
  pct_promoter = vapply(anno_dfs[sex_levels], function(d)
                        round(100 * mean(grepl("^Promoter", d$annotation)), 1),
                        numeric(1)),
  row.names = NULL
)
addWorksheet(wb, "Summary")
writeData(wb, "Summary", summary_df)

saveWorkbook(wb, file.path(OUTDIR, "annotated_peaks.xlsx"), overwrite = TRUE)

# Saved so the annotation does not need to be repeated
saveRDS(peakAnnoList, file.path(OUTDIR, "peakAnnoList.rds"))
saveRDS(specificity,  file.path(OUTDIR, "specificity.rds"))

# CSV copies for the downstream Python scripts
write.csv(merged, file.path(OUTDIR, "all_peaks.csv"), row.names = FALSE)
write.csv(specificity, file.path(OUTDIR, "sex_specificity.csv"), row.names = FALSE)

print(summary_df)
message("Done. Outputs in: ", normalizePath(OUTDIR))
