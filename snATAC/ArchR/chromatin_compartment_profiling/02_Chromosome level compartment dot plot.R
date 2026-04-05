#!/usr/bin/env Rscript
# ==============================================================================
# Chromosome level compartment dot plots
# ==============================================================================
#
# Description:
#   Generates chromosome-level visualizations of A/B compartment activity
#   across aging, sex, and cell type. Exports underlying data as a single
#   Excel workbook with one tab per panel.
#
# Input:
#   - compartments_binary_A1_R0.rds (from compartment_01_main_analysis.R)
#   - manifest.normalized.tsv
#   - bins_table.tsv
#   - compartment_stability_5class.tsv
#
# Output:
#   Figures (PDF + PNG):
#     Panel A - Age-domain vs non-age-domain proportion per chromosome
#     Panel B - Dot plot: Chr x Celltype, fraction active (by sex)
#     Panel C - Dot plot: Chr x Age, fraction active (by sex) - all cell types
#     Panel D - Dot plot: Chr x Age, fraction active (by sex) - Hepatocyte
#
#   Data:
#     compartment_dotplot_data.xlsx (4 tabs matching panels A-D)
#
# Note:
#   Color scale uses data-driven midpoint (center of observed range) rather
#   than fixed 0.5, so the full blue-white-red gradient spans actual values.
#
#
# ==============================================================================


# ==============================================================================
# SETUP
# ==============================================================================

suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(scales)
    library(openxlsx)
})

# ------------------------------------------------------------------------------
# Configuration - Update this path
# ------------------------------------------------------------------------------

outdir <- "path/to/output/compartment_main_heatmap"

# ------------------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------------------

message("[SETUP] Loading data...")

comp_bin     <- readRDS(file.path(outdir, "compartments_binary_A1_R0.rds"))
manifest     <- fread(file.path(outdir, "manifest.normalized.tsv"))
bins_df      <- fread(file.path(outdir, "bins_table.tsv"))
stability_dt <- fread(file.path(outdir, "compartment_stability_5class.tsv"))

age_vec   <- setNames(manifest$age, manifest$group)
sex_vec   <- setNames(manifest$sex, manifest$group)
ctype_vec <- setNames(manifest$celltype, manifest$group)

# Constants
AGE_LEVELS        <- c("young", "mid_age", "old", "pre_geriatric", "geriatric")
CHR_ORDER         <- paste0("chr", 1:19)
SWITCHING_CLASSES <- c("Monotonic_A_to_R", "Monotonic_R_to_A", "Non_Monotonic")

message("  [OK] Data loaded")


# ==============================================================================
# HELPER FUNCTION
# ==============================================================================

# Compute data-driven color scale limits and midpoint.
# Using midpoint = center of data range ensures the full blue-white-red
# gradient spans the observed values, avoiding washed-out colors when
# the range is narrow (e.g. 55%-65% would look uniformly pale with
# a fixed midpoint of 0.5).
make_color_scale <- function(frac_vec) {
    lo  <- floor(min(frac_vec, na.rm = TRUE) * 100) / 100
    hi  <- ceiling(max(frac_vec, na.rm = TRUE) * 100) / 100
    mid <- (lo + hi) / 2
    list(lo = lo, hi = hi, mid = mid)
}

# Shared ggplot theme for dot plots
theme_dotplot <- function(base_size = 11) {
    theme_bw(base_size = base_size) %+replace%
        theme(
            text = element_text(family = "Arial"),
            axis.text.x  = element_text(angle = 45, hjust = 1, face = "bold", size = 9),
            axis.text.y  = element_text(face = "bold", size = 9),
            strip.text   = element_text(face = "bold", size = 12),
            strip.background = element_rect(fill = "gray95"),
            panel.grid.minor = element_blank(),
            panel.grid.major = element_line(linewidth = 0.3, color = "gray90"),
            legend.position  = "right",
            plot.title = element_text(face = "bold", hjust = 0.5, size = 13)
        )
}


# ==============================================================================
# BUILD LONG-FORMAT DATA
# ==============================================================================

message("[DATA] Building long-format compartment data...")

dt_list <- list()
for (sn in colnames(comp_bin)) {
    dt_list[[sn]] <- data.table(
        bin_id   = rownames(comp_bin),
        state    = comp_bin[, sn],
        age      = age_vec[sn],
        sex      = sex_vec[sn],
        celltype = ctype_vec[sn]
    )
}
all_data <- rbindlist(dt_list)
all_data <- merge(all_data, bins_df[, .(bin_id, chr)], by = "bin_id")
all_data <- all_data[chr %in% CHR_ORDER]

message(sprintf("  [OK] %d rows", nrow(all_data)))


# ==============================================================================
# PANEL A: AGE-DOMAIN VS NON-AGE-DOMAIN PROPORTION PER CHROMOSOME
# ==============================================================================

banner_a <- "\n======================================================================
PANEL A: Age-domain proportion per chromosome
======================================================================"
message(banner_a)

bin_stability <- stability_dt[, .(
    is_age_domain = any(Stability %in% SWITCHING_CLASSES)
), by = bin_id]

bin_stability <- merge(bin_stability, bins_df[, .(bin_id, chr)], by = "bin_id")
bin_stability <- bin_stability[chr %in% CHR_ORDER]

chr_proportion <- bin_stability[, .(
    n_total          = .N,
    n_age_domain     = sum(is_age_domain),
    n_non_age_domain = sum(!is_age_domain)
), by = chr]

chr_proportion[, frac_age_domain     := n_age_domain / n_total]
chr_proportion[, frac_non_age_domain := n_non_age_domain / n_total]

# Order chromosomes by age-domain fraction (descending)
chr_order_by_prop <- chr_proportion[order(-frac_age_domain), as.character(chr)]
chr_proportion[, chr := factor(chr, levels = chr_order_by_prop)]

chr_prop_long <- melt(
    chr_proportion,
    id.vars = "chr",
    measure.vars = c("frac_age_domain", "frac_non_age_domain"),
    variable.name = "Domain",
    value.name = "Proportion"
)
chr_prop_long[, Domain := fifelse(
    Domain == "frac_age_domain", "age-domain", "non-age-domain"
)]
chr_prop_long[, Domain := factor(Domain, levels = c("non-age-domain", "age-domain"))]

pA <- ggplot(chr_prop_long, aes(x = Proportion, y = chr, fill = Domain)) +
    geom_col(position = "stack", color = NA, width = 0.8) +
    scale_fill_manual(
        values = c("non-age-domain" = "#BDBDBD", "age-domain" = "#E15759"),
        name = NULL
    ) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.02)),
                       breaks = seq(0, 1, 0.25)) +
    labs(x = "proportion", y = "chromosome") +
    theme_minimal(base_size = 12) +
    theme(
        text = element_text(family = "Arial"),
        axis.text.y  = element_text(face = "bold", size = 10),
        axis.text.x  = element_text(size = 10),
        axis.title   = element_text(face = "bold"),
        legend.position = "top",
        legend.text  = element_text(size = 10),
        panel.grid.major.y = element_blank(),
        panel.grid.minor   = element_blank()
    )

ggsave(file.path(outdir, "PanelA_AgeDomain_Proportion_Chr.pdf"),
       pA, width = 5, height = 7, useDingbats = FALSE)
ggsave(file.path(outdir, "PanelA_AgeDomain_Proportion_Chr.png"),
       pA, width = 5, height = 7, dpi = 300)

message("  [OK] Panel A saved")


# ==============================================================================
# PANEL B: DOT PLOT - CHR x CELLTYPE (BY SEX)
# ==============================================================================

banner_b <- "\n======================================================================
PANEL B: Fraction active per Chr x Celltype
======================================================================"
message(banner_b)

chr_ct_summary <- all_data[, .(
    frac_active = mean(state == 1L, na.rm = TRUE),
    n_bins      = .N
), by = .(sex, chr, celltype)]

chr_ct_summary[, chr := factor(chr, levels = rev(CHR_ORDER))]
chr_ct_summary[, sex := factor(sex, levels = c("male", "female"))]

cs_B <- make_color_scale(chr_ct_summary$frac_active)

pB <- ggplot(chr_ct_summary,
             aes(x = celltype, y = chr,
                 size = frac_active, color = frac_active)) +
    geom_point() +
    facet_wrap(~ sex, ncol = 2) +
    scale_size_continuous(
        range = c(2, 10), name = "fraction active",
        labels = percent, limits = c(cs_B$lo, cs_B$hi)
    ) +
    scale_color_gradient2(
        low = "#2166AC", mid = "white", high = "#B2182B",
        midpoint = cs_B$mid, limits = c(cs_B$lo, cs_B$hi),
        name = "fraction active", labels = percent
    ) +
    labs(title = "all cell types", x = NULL, y = "chromosome") +
    theme_dotplot()

ggsave(file.path(outdir, "PanelB_DotPlot_Chr_Celltype_BySex.pdf"),
       pB, width = 14, height = 8, useDingbats = FALSE)
ggsave(file.path(outdir, "PanelB_DotPlot_Chr_Celltype_BySex.png"),
       pB, width = 14, height = 8, dpi = 300)

message(sprintf("  [OK] Panel B saved [range: %d%%-%d%%, midpoint: %d%%]",
                round(cs_B$lo * 100), round(cs_B$hi * 100), round(cs_B$mid * 100)))


# ==============================================================================
# PANEL C: DOT PLOT - CHR x AGE (BY SEX) - ALL CELL TYPES
# ==============================================================================

banner_c <- "\n======================================================================
PANEL C: Fraction active per Chr x Age (all cell types)
======================================================================"
message(banner_c)

chr_age_summary <- all_data[, .(
    frac_active = mean(state == 1L, na.rm = TRUE),
    n_bins      = .N
), by = .(sex, chr, age)]

chr_age_summary[, chr := factor(chr, levels = rev(CHR_ORDER))]
chr_age_summary[, age := factor(age, levels = AGE_LEVELS)]
chr_age_summary[, sex := factor(sex, levels = c("male", "female"))]

cs_C <- make_color_scale(chr_age_summary$frac_active)

pC <- ggplot(chr_age_summary,
             aes(x = age, y = chr,
                 size = frac_active, color = frac_active)) +
    geom_point() +
    facet_wrap(~ sex, ncol = 2) +
    scale_size_continuous(
        range = c(3, 12), name = "fraction active",
        labels = percent, limits = c(cs_C$lo, cs_C$hi)
    ) +
    scale_color_gradient2(
        low = "#2166AC", mid = "white", high = "#B2182B",
        midpoint = cs_C$mid, limits = c(cs_C$lo, cs_C$hi),
        name = "fraction active", labels = percent
    ) +
    labs(title = "all cell types", x = NULL, y = "chromosome") +
    theme_dotplot() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 11))

ggsave(file.path(outdir, "PanelC_DotPlot_Chr_Age_AllCelltypes_BySex.pdf"),
       pC, width = 12, height = 8, useDingbats = FALSE)
ggsave(file.path(outdir, "PanelC_DotPlot_Chr_Age_AllCelltypes_BySex.png"),
       pC, width = 12, height = 8, dpi = 300)

message(sprintf("  [OK] Panel C saved [range: %d%%-%d%%, midpoint: %d%%]",
                round(cs_C$lo * 100), round(cs_C$hi * 100), round(cs_C$mid * 100)))


# ==============================================================================
# PANEL D: DOT PLOT - CHR x AGE (BY SEX) - HEPATOCYTE
# ==============================================================================

banner_d <- "\n======================================================================
PANEL D: Fraction active per Chr x Age (Hepatocyte)
======================================================================"
message(banner_d)

hep_data <- all_data[celltype == "Hepatocyte"]

hep_chr_age <- hep_data[, .(
    frac_active = mean(state == 1L, na.rm = TRUE),
    n_bins      = .N
), by = .(sex, chr, age)]

hep_chr_age[, chr := factor(chr, levels = rev(CHR_ORDER))]
hep_chr_age[, age := factor(age, levels = AGE_LEVELS)]
hep_chr_age[, sex := factor(sex, levels = c("male", "female"))]

cs_D <- make_color_scale(hep_chr_age$frac_active)

pD <- ggplot(hep_chr_age,
             aes(x = age, y = chr,
                 size = frac_active, color = frac_active)) +
    geom_point() +
    facet_wrap(~ sex, ncol = 2) +
    scale_size_continuous(
        range = c(3, 12), name = "fraction active",
        labels = percent, limits = c(cs_D$lo, cs_D$hi)
    ) +
    scale_color_gradient2(
        low = "#2166AC", mid = "white", high = "#B2182B",
        midpoint = cs_D$mid, limits = c(cs_D$lo, cs_D$hi),
        name = "fraction active", labels = percent
    ) +
    labs(title = "Hepatocyte", x = NULL, y = "chromosome") +
    theme_dotplot() +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 11),
        plot.title  = element_text(face = "bold", hjust = 0.5, size = 14)
    )

ggsave(file.path(outdir, "PanelD_DotPlot_Chr_Age_Hepatocyte_BySex.pdf"),
       pD, width = 12, height = 8, useDingbats = FALSE)
ggsave(file.path(outdir, "PanelD_DotPlot_Chr_Age_Hepatocyte_BySex.png"),
       pD, width = 12, height = 8, dpi = 300)

message(sprintf("  [OK] Panel D saved [range: %d%%-%d%%, midpoint: %d%%]",
                round(cs_D$lo * 100), round(cs_D$hi * 100), round(cs_D$mid * 100)))


# ==============================================================================
# EXCEL EXPORT (4 TABS)
# ==============================================================================

banner_e <- "\n======================================================================
EXCEL EXPORT: compartment_dotplot_data.xlsx
======================================================================"
message(banner_e)

# Tab 1: Panel A
tab1 <- chr_proportion[, .(
    chr, n_total, n_age_domain, n_non_age_domain,
    pct_age_domain     = round(frac_age_domain * 100, 2),
    pct_non_age_domain = round(frac_non_age_domain * 100, 2)
)]
tab1[, chr := factor(chr, levels = CHR_ORDER)]
tab1 <- tab1[order(-pct_age_domain)]

# Tab 2: Panel B
tab2 <- chr_ct_summary[, .(
    chr, celltype, sex,
    pct_active = round(frac_active * 100, 2),
    n_bins
)]
tab2[, chr := factor(chr, levels = CHR_ORDER)]
tab2 <- tab2[order(chr, celltype, sex)]

# Tab 3: Panel C
tab3 <- chr_age_summary[, .(
    chr, age, sex,
    pct_active = round(frac_active * 100, 2),
    n_bins
)]
tab3[, chr := factor(chr, levels = CHR_ORDER)]
tab3[, age := factor(age, levels = AGE_LEVELS)]
tab3 <- tab3[order(chr, age, sex)]

# Tab 4: Panel D
tab4 <- hep_chr_age[, .(
    chr, age, sex,
    pct_active = round(frac_active * 100, 2),
    n_bins
)]
tab4[, chr := factor(chr, levels = CHR_ORDER)]
tab4[, age := factor(age, levels = AGE_LEVELS)]
tab4 <- tab4[order(chr, age, sex)]

# Write workbook
wb <- createWorkbook()

addWorksheet(wb, "PanelA_AgeDomain")
addWorksheet(wb, "PanelB_Chr_Celltype")
addWorksheet(wb, "PanelC_Chr_Age_All")
addWorksheet(wb, "PanelD_Chr_Age_Hep")

hs <- createStyle(textDecoration = "bold", fontName = "Arial", fontSize = 11)

writeData(wb, "PanelA_AgeDomain",    tab1, headerStyle = hs)
writeData(wb, "PanelB_Chr_Celltype", tab2, headerStyle = hs)
writeData(wb, "PanelC_Chr_Age_All",  tab3, headerStyle = hs)
writeData(wb, "PanelD_Chr_Age_Hep",  tab4, headerStyle = hs)

setColWidths(wb, "PanelA_AgeDomain",    cols = 1:6, widths = "auto")
setColWidths(wb, "PanelB_Chr_Celltype", cols = 1:5, widths = "auto")
setColWidths(wb, "PanelC_Chr_Age_All",  cols = 1:5, widths = "auto")
setColWidths(wb, "PanelD_Chr_Age_Hep",  cols = 1:5, widths = "auto")

outfile <- file.path(outdir, "compartment_dotplot_data.xlsx")
saveWorkbook(wb, outfile, overwrite = TRUE)

message(sprintf("  [OK] Saved: %s", outfile))


# ==============================================================================
# SESSION INFO
# ==============================================================================

writeLines(
    capture.output(sessionInfo()),
    file.path(outdir, "sessionInfo_dotplot.txt")
)


# ==============================================================================
# SUMMARY
# ==============================================================================

cat("\n")
cat("========================================================================\n")
cat("  COMPLETE\n")
cat("========================================================================\n")
cat("\n")
cat("  Figures (PDF + PNG):\n")
cat("    PanelA_AgeDomain_Proportion_Chr\n")
cat("    PanelB_DotPlot_Chr_Celltype_BySex\n")
cat("    PanelC_DotPlot_Chr_Age_AllCelltypes_BySex\n")
cat("    PanelD_DotPlot_Chr_Age_Hepatocyte_BySex\n")
cat("\n")
cat("  Data:\n")
cat("    compartment_dotplot_data.xlsx (4 tabs)\n")
cat("\n")
cat("  Color scales (data-driven midpoint):\n")
cat(sprintf("    Panel B: %d%%-%d%% (mid=%d%%)\n",
            round(cs_B$lo * 100), round(cs_B$hi * 100), round(cs_B$mid * 100)))
cat(sprintf("    Panel C: %d%%-%d%% (mid=%d%%)\n",
            round(cs_C$lo * 100), round(cs_C$hi * 100), round(cs_C$mid * 100)))
cat(sprintf("    Panel D: %d%%-%d%% (mid=%d%%)\n",
            round(cs_D$lo * 100), round(cs_D$hi * 100), round(cs_D$mid * 100)))
cat("\n")
cat("  Reproducibility:\n")
cat("    sessionInfo_dotplot.txt\n")
cat("========================================================================\n")
