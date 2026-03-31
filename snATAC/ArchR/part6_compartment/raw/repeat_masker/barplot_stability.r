#!/usr/bin/env Rscript
## =========================================================
## PLOT: Repeat Class Percentage by Stability Class and Sex
## 
## Single stacked bar plot showing % of each repeat class
## within each stability class, faceted by sex
## =========================================================

cat("\n")
cat("╔══════════════════════════════════════════════════════════════════╗\n")
cat("║   REPEAT CLASS PERCENTAGE BY STABILITY CLASS (BY SEX)            ║\n")
cat("╚══════════════════════════════════════════════════════════════════╝\n\n")

## =========================================================
## CONFIGURATION
## =========================================================

rmsk_file <- "/data/sarkern2/multiome_liver/Seurat/epigenome/repeatmasker/mm10_repeatmasker_processed.rds"
base_outdir <- "/data/sarkern2/multiome_liver/Seurat/epigenome/downstream_stability"

# Choose version: "default" or "relaxed"
VERSION <- "default"

p2g_dir <- file.path(base_outdir, VERSION, "peak2gene_stability_analysis")
output_dir <- "/data/sarkern2/multiome_liver/Seurat/epigenome/repeatmasker_enrichment"

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)


## =========================================================
## LOAD PACKAGES
## =========================================================

suppressPackageStartupMessages({
    library(GenomicRanges)
    library(data.table)
    library(ggplot2)
})


## =========================================================
## CONSTANTS - BIOLOGICAL ORDER
## =========================================================

# Stability classes in biological order:
# Active → Dynamic/Switching → Repressive
STABILITY_CLASSES <- c(
    "Stable_Active",
    "Non_Monotonic",
    "Monotonic_R_to_A",
    "Monotonic_A_to_R",
    "Stable_Repressive"
)

# Repeat classes in biological order:
# Euchromatin-associated → Heterochromatin-associated
MAJOR_REPEAT_CLASSES <- c(
    "SINE",            # Euchromatin (A compartment)
    "Simple_repeat",   # Neutral
    "Low_complexity",  # Neutral
    "RNA",             # Functional RNA
    "Rolling_Circle",  # DNA transposon-like
    "DNA_Transposon",  # DNA transposons
    "LTR",             # ERVs (can be regulatory)
    "LINE",            # Heterochromatin (B compartment)
    "Satellite",       # Constitutive heterochromatin
    "Other"            # Uncertain/Unknown
)

# Paper-grade, colorblind-safe, muted palette
REPEAT_COLS <- c(
    "SINE"           = "#76B7B2",
    "Simple_repeat"  = "#B07AA1",
    "LINE"           = "#E15759",
    "LTR"            = "#EDC948",
    "Low_complexity" = "#FF9DA7",
    "DNA_Transposon" = "#59A14F",
    "Satellite"      = "#F28E2B",
    "Other"          = "#9C9C9C",
    "RNA"            = "#4E79A7",
    "Rolling_Circle" = "#1F9ED4"
)

# Stability labels
STABILITY_LABELS <- c(
    "Stable_Active"     = "Stable\nActive",
    "Non_Monotonic"     = "Non-\nMonotonic",
    "Monotonic_R_to_A"  = "R → A\n(Opening)",
    "Monotonic_A_to_R"  = "A → R\n(Closing)",
    "Stable_Repressive" = "Stable\nRepressive"
)


## =========================================================
## LOAD DATA
## =========================================================

message("[1] Loading data...")

# RepeatMasker
rmsk_gr <- readRDS(rmsk_file)
message(sprintf("    RepeatMasker: %s elements", format(length(rmsk_gr), big.mark = ",")))

# P2G peaks
peaks_gr <- readRDS(file.path(p2g_dir, "hepatocyte_p2g_peaks_gr.rds"))
message(sprintf("    P2G peaks: %d", length(peaks_gr)))

# Stability
stability_gr_list <- readRDS(file.path(p2g_dir, "stability_gr_list.rds"))
message(sprintf("    Stability: %s", paste(names(stability_gr_list), collapse = ", ")))


## =========================================================
## OVERLAP PEAKS WITH REPEATMASKER
## =========================================================

message("\n[2] Overlapping peaks with RepeatMasker...")

olaps <- findOverlaps(peaks_gr, rmsk_gr, ignore.strand = TRUE)

olap_dt <- data.table(
    peak_idx = queryHits(olaps),
    peak_id = mcols(peaks_gr)$peak_id[queryHits(olaps)],
    Major_Class = mcols(rmsk_gr)$Major_Class[subjectHits(olaps)]
)

olap_dt <- olap_dt[Major_Class %in% MAJOR_REPEAT_CLASSES]
message(sprintf("    Overlaps: %s", format(nrow(olap_dt), big.mark = ",")))


## =========================================================
## CALCULATE PERCENTAGES FOR EACH SEX × STABILITY
## =========================================================

message("\n[3] Calculating percentages...")

results_list <- list()

for (sex in c("male", "female")) {
    
    if (!sex %in% names(stability_gr_list)) next
    
    for (stab_class in STABILITY_CLASSES) {
        
        stab_bins <- stability_gr_list[[sex]][[stab_class]]
        if (is.null(stab_bins) || length(stab_bins) == 0) next
        
        # Get peaks in this stability class
        hits <- findOverlaps(peaks_gr, stab_bins, ignore.strand = TRUE)
        peak_ids <- unique(mcols(peaks_gr)$peak_id[queryHits(hits)])
        
        if (length(peak_ids) == 0) next
        
        # Count repeats
        stab_olaps <- olap_dt[peak_id %in% peak_ids]
        
        if (nrow(stab_olaps) == 0) next
        
        # Count by repeat class
        counts <- stab_olaps[, .(N = .N), by = Major_Class]
        total <- sum(counts$N)
        counts[, Percentage := N / total * 100]
        counts[, Sex := sex]
        counts[, Stability := stab_class]
        counts[, N_peaks := length(peak_ids)]
        counts[, Total_overlaps := total]
        
        results_list[[paste(sex, stab_class, sep = "_")]] <- counts
    }
}

plot_dt <- rbindlist(results_list, fill = TRUE)

# Ensure all combinations exist
all_combos <- CJ(
    Sex = c("male", "female"),
    Stability = STABILITY_CLASSES,
    Major_Class = MAJOR_REPEAT_CLASSES
)
plot_dt <- merge(all_combos, plot_dt, by = c("Sex", "Stability", "Major_Class"), all.x = TRUE)
plot_dt[is.na(Percentage), Percentage := 0]

# Factor ordering (biological order)
plot_dt[, Sex := factor(Sex, levels = c("male", "female"), labels = c("Male", "Female"))]
plot_dt[, Stability := factor(Stability, levels = STABILITY_CLASSES)]
plot_dt[, Major_Class := factor(Major_Class, levels = MAJOR_REPEAT_CLASSES)]

message("    Done!")


## =========================================================
## CREATE PLOT
## =========================================================

message("\n[4] Creating plot...")

p <- ggplot(plot_dt, aes(x = Stability, y = Percentage, fill = Major_Class)) +
    geom_col(position = "stack", color = "black", linewidth = 0.7, width = 1.0) +
    facet_wrap(~ Sex, ncol = 2) +
    scale_fill_manual(values = REPEAT_COLS, name = "Repeat Class") +
    scale_x_discrete(labels = STABILITY_LABELS) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.02)), 
                       breaks = seq(0, 100, 20)) +
    labs(
        title = "Repeat Element Composition by Stability Class",
        subtitle = sprintf("Hepatocyte P2G Peaks (%s) — Percentage of repeat overlaps", toupper(VERSION)),
        x = NULL,
        y = "Percentage (%)"
    ) +
    theme_minimal(base_size = 14) +
    theme(
        plot.title = element_text(face = "bold", size = 18, hjust = 0.5, 
                                  margin = margin(b = 5)),
        plot.subtitle = element_text(size = 12, hjust = 0.5, color = "gray40",
                                     margin = margin(b = 15)),
        strip.text = element_text(face = "bold", size = 16, color = "gray20"),
        strip.background = element_rect(fill = "gray95", color = NA),
        axis.text.x = element_text(face = "bold", size = 10, color = "gray20",
                                   lineheight = 0.9),
        axis.text.y = element_text(size = 12, color = "gray30"),
        axis.title.y = element_text(face = "bold", size = 13, color = "gray20",
                                    margin = margin(r = 10)),
        legend.position = "right",
        legend.title = element_text(face = "bold", size = 12),
        legend.text = element_text(size = 11),
        legend.key.size = unit(0.5, "cm"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_line(color = "gray90", linewidth = 0.5),
        panel.spacing = unit(1, "cm"),
        plot.margin = margin(20, 20, 20, 20),
        plot.background = element_rect(fill = "white", color = NA)
    )


## =========================================================
## SAVE PLOT
## =========================================================

outfile_pdf <- file.path(output_dir, sprintf("RepeatClass_Percentage_by_Stability_Sex_%s.pdf", VERSION))
outfile_png <- file.path(output_dir, sprintf("RepeatClass_Percentage_by_Stability_Sex_%s.png", VERSION))

ggsave(outfile_pdf, p, width = 10, height = 6)


message(sprintf("\n✅ Saved: %s", outfile_pdf))
message(sprintf("✅ Saved: %s", outfile_png))


## =========================================================
## PRINT SUMMARY TABLE
## =========================================================

message("\n[5] Summary table:\n")

summary_wide <- dcast(plot_dt, Sex + Stability ~ Major_Class, value.var = "Percentage")
print(summary_wide)

# Save table
fwrite(summary_wide, file.path(output_dir, sprintf("RepeatClass_Percentage_Summary_%s.tsv", VERSION)), sep = "\t")

message("\n✅ Done!")



