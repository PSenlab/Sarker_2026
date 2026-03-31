#!/usr/bin/env Rscript
## =========================================================
## COMPARTMENT FIGURE: 4 Panels for GitHub
##
## 1. Panel A: Age-domain vs Non-age-domain proportion per chr
## 2. Panel B: Dot plot - Chr × Celltype, fraction active (by sex)
## 3. Panel C: Dot plot - Chr × Age, fraction active (by sex) - All celltypes
## 4. Panel D: Dot plot - Chr × Age, fraction active (by sex) - Hepatocyte
##
## Color fix: midpoint = center of DATA range (not 0.5)
##   so full blue-white-red gradient spans your actual values
## =========================================================

suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(scales)
})

## =========================================================
## CONFIGURATION
## =========================================================

outdir <- "/data/sarkern2/multiome_liver/Seurat/epigenome/compartment_main_heatmap"

comp_bin      <- readRDS(file.path(outdir, "compartments_binary_A1_R0.rds"))
manifest      <- fread(file.path(outdir, "manifest.normalized.tsv"))
bins_df       <- fread(file.path(outdir, "bins_table.tsv"))
stability_dt  <- fread(file.path(outdir, "compartment_stability_5class.tsv"))

age_vec   <- setNames(manifest$age, manifest$group)
sex_vec   <- setNames(manifest$sex, manifest$group)
ctype_vec <- setNames(manifest$celltype, manifest$group)

age_levels <- c("young", "mid_age", "old", "pre_geriatric", "geriatric")
chr_order  <- paste0("chr", 1:19)

switching_classes <- c("Monotonic_A_to_R", "Monotonic_R_to_A", "Non_Monotonic")

## =========================================================
## HELPER: compute data-driven limits + midpoint
##
## KEY INSIGHT: gradient2 with midpoint=0.5 makes 60% look
## salmon because it's only 20% of the way from white→red.
## Setting midpoint = center of data range means the FULL
## gradient spans your actual values → deep saturated colors.
## =========================================================

make_color_scale <- function(frac_vec) {
    lo  <- floor(min(frac_vec, na.rm = TRUE) * 100) / 100
    hi  <- ceiling(max(frac_vec, na.rm = TRUE) * 100) / 100
    mid <- (lo + hi) / 2
    list(lo = lo, hi = hi, mid = mid)
}

## =========================================================
## BUILD LONG-FORMAT DATA
## =========================================================

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
all_data <- all_data[chr %in% chr_order]

message("    ✓ ", nrow(all_data), " rows")


## =========================================================
## PANEL A: Age-domain vs Non-age-domain proportion per chr
## =========================================================

message("[Panel A] Age-domain proportion per chromosome...")

bin_stability <- stability_dt[, .(
    is_age_domain = any(Stability %in% switching_classes)
), by = bin_id]

bin_stability <- merge(bin_stability, bins_df[, .(bin_id, chr)], by = "bin_id")
bin_stability <- bin_stability[chr %in% chr_order]

chr_proportion <- bin_stability[, .(
    n_total          = .N,
    n_age_domain     = sum(is_age_domain),
    n_non_age_domain = sum(!is_age_domain)
), by = chr]

chr_proportion[, frac_age_domain     := n_age_domain / n_total]
chr_proportion[, frac_non_age_domain := n_non_age_domain / n_total]

chr_order_by_prop <- chr_proportion[order(-frac_age_domain), as.character(chr)]
chr_proportion[, chr := factor(chr, levels = chr_order_by_prop)]

chr_prop_long <- melt(chr_proportion,
                      id.vars = "chr",
                      measure.vars = c("frac_age_domain", "frac_non_age_domain"),
                      variable.name = "Domain",
                      value.name = "Proportion")

chr_prop_long[, Domain := fifelse(Domain == "frac_age_domain", "age-domain", "non-age-domain")]
chr_prop_long[, Domain := factor(Domain, levels = c("non-age-domain", "age-domain"))]

pA <- ggplot(chr_prop_long, aes(x = Proportion, y = chr, fill = Domain)) +
    geom_col(position = "stack", color = NA, width = 0.8) +
    scale_fill_manual(values = c("non-age-domain" = "#BDBDBD",
                                  "age-domain"     = "#E15759"),
                      name = NULL) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.02)),
                       breaks = seq(0, 1, 0.25)) +
    labs(x = "proportion", y = "chromosome") +
    theme_minimal(base_size = 12) +
    theme(
        axis.text.y  = element_text(face = "bold", size = 10),
        axis.text.x  = element_text(size = 10),
        axis.title   = element_text(face = "bold"),
        legend.position = "top",
        legend.text  = element_text(size = 10),
        panel.grid.major.y = element_blank(),
        panel.grid.minor   = element_blank()
    )

ggsave(file.path(outdir, "PanelA_AgeDomain_Proportion_Chr.pdf"),
       pA, width = 5, height = 7)
ggsave(file.path(outdir, "PanelA_AgeDomain_Proportion_Chr.png"),
       pA, width = 5, height = 7, dpi = 300)

message("    ✅ Panel A saved")


## =========================================================
## PANEL B: Dot plot - Chr × Celltype (by sex)
## =========================================================

message("[Panel B] Fraction active per Chr × Celltype...")

chr_ct_summary <- all_data[, .(
    frac_active = mean(state == 1L, na.rm = TRUE),
    n_bins      = .N
), by = .(sex, chr, celltype)]

chr_ct_summary[, chr := factor(chr, levels = rev(chr_order))]
chr_ct_summary[, sex := factor(sex, levels = c("male", "female"))]

cs_B <- make_color_scale(chr_ct_summary$frac_active)

pB <- ggplot(chr_ct_summary,
             aes(x = celltype, y = chr,
                 size = frac_active, color = frac_active)) +
    geom_point() +
    facet_wrap(~ sex, ncol = 2) +
    scale_size_continuous(range  = c(2, 10),
                          name   = "fraction active",
                          labels = percent,
                          limits = c(cs_B$lo, cs_B$hi)) +
    scale_color_gradient2(low      = "#2166AC",
                          mid      = "white",
                          high     = "#B2182B",
                          midpoint = cs_B$mid,
                          limits   = c(cs_B$lo, cs_B$hi),
                          name     = "fraction active",
                          labels   = percent) +
    labs(title = "all cell types",
         x = NULL, y = "chromosome") +
    theme_bw(base_size = 11) +
    theme(
        axis.text.x  = element_text(angle = 45, hjust = 1, face = "bold", size = 9),
        axis.text.y  = element_text(face = "bold", size = 9),
        strip.text   = element_text(face = "bold", size = 12),
        strip.background = element_rect(fill = "gray95"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linewidth = 0.3, color = "gray90"),
        legend.position  = "right",
        plot.title = element_text(face = "bold", hjust = 0.5, size = 13)
    )

ggsave(file.path(outdir, "PanelB_DotPlot_Chr_Celltype_BySex.pdf"),
       pB, width = 14, height = 8)
ggsave(file.path(outdir, "PanelB_DotPlot_Chr_Celltype_BySex.png"),
       pB, width = 14, height = 8, dpi = 300)

message("    ✅ Panel B saved  [range: ", round(cs_B$lo*100),
        "% – ", round(cs_B$hi*100), "%, midpoint: ", round(cs_B$mid*100), "%]")


## =========================================================
## PANEL C: Dot plot - Chr × Age (by sex) - ALL CELLTYPES
## =========================================================

message("[Panel C] Fraction active per Chr × Age (all celltypes)...")

chr_age_summary <- all_data[, .(
    frac_active = mean(state == 1L, na.rm = TRUE),
    n_bins      = .N
), by = .(sex, chr, age)]

chr_age_summary[, chr := factor(chr, levels = rev(chr_order))]
chr_age_summary[, age := factor(age, levels = age_levels)]
chr_age_summary[, sex := factor(sex, levels = c("male", "female"))]

cs_C <- make_color_scale(chr_age_summary$frac_active)

pC <- ggplot(chr_age_summary,
             aes(x = age, y = chr,
                 size = frac_active, color = frac_active)) +
    geom_point() +
    facet_wrap(~ sex, ncol = 2) +
    scale_size_continuous(range  = c(3, 12),
                          name   = "fraction active",
                          labels = percent,
                          limits = c(cs_C$lo, cs_C$hi)) +
    scale_color_gradient2(low      = "#2166AC",
                          mid      = "white",
                          high     = "#B2182B",
                          midpoint = cs_C$mid,
                          limits   = c(cs_C$lo, cs_C$hi),
                          name     = "fraction active",
                          labels   = percent) +
    labs(title = "all cell types",
         x = NULL, y = "chromosome") +
    theme_bw(base_size = 11) +
    theme(
        axis.text.x  = element_text(angle = 45, hjust = 1, face = "bold", size = 11),
        axis.text.y  = element_text(face = "bold", size = 9),
        strip.text   = element_text(face = "bold", size = 12),
        strip.background = element_rect(fill = "gray95"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linewidth = 0.3, color = "gray90"),
        legend.position  = "right",
        plot.title = element_text(face = "bold", hjust = 0.5, size = 13)
    )

ggsave(file.path(outdir, "PanelC_DotPlot_Chr_Age_AllCelltypes_BySex.pdf"),
       pC, width = 12, height = 8)
ggsave(file.path(outdir, "PanelC_DotPlot_Chr_Age_AllCelltypes_BySex.png"),
       pC, width = 12, height = 8, dpi = 300)

message("    ✅ Panel C saved  [range: ", round(cs_C$lo*100),
        "% – ", round(cs_C$hi*100), "%, midpoint: ", round(cs_C$mid*100), "%]")


## =========================================================
## PANEL D: Dot plot - Chr × Age (by sex) - HEPATOCYTE
## =========================================================

message("[Panel D] Fraction active per Chr × Age (Hepatocyte)...")

hep_data <- all_data[celltype == "Hepatocyte"]

hep_chr_age <- hep_data[, .(
    frac_active = mean(state == 1L, na.rm = TRUE),
    n_bins      = .N
), by = .(sex, chr, age)]

hep_chr_age[, chr := factor(chr, levels = rev(chr_order))]
hep_chr_age[, age := factor(age, levels = age_levels)]
hep_chr_age[, sex := factor(sex, levels = c("male", "female"))]

cs_D <- make_color_scale(hep_chr_age$frac_active)

pD <- ggplot(hep_chr_age,
             aes(x = age, y = chr,
                 size = frac_active, color = frac_active)) +
    geom_point() +
    facet_wrap(~ sex, ncol = 2) +
    scale_size_continuous(range  = c(3, 12),
                          name   = "fraction active",
                          labels = percent,
                          limits = c(cs_D$lo, cs_D$hi)) +
    scale_color_gradient2(low      = "#2166AC",
                          mid      = "white",
                          high     = "#B2182B",
                          midpoint = cs_D$mid,
                          limits   = c(cs_D$lo, cs_D$hi),
                          name     = "fraction active",
                          labels   = percent) +
    labs(title = "Hepatocyte",
         x = NULL, y = "chromosome") +
    theme_bw(base_size = 11) +
    theme(
        axis.text.x  = element_text(angle = 45, hjust = 1, face = "bold", size = 11),
        axis.text.y  = element_text(face = "bold", size = 9),
        strip.text   = element_text(face = "bold", size = 12),
        strip.background = element_rect(fill = "gray95"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linewidth = 0.3, color = "gray90"),
        legend.position  = "right",
        plot.title = element_text(face = "bold", hjust = 0.5, size = 14)
    )

ggsave(file.path(outdir, "PanelD_DotPlot_Chr_Age_Hepatocyte_BySex.pdf"),
       pD, width = 12, height = 8)
ggsave(file.path(outdir, "PanelD_DotPlot_Chr_Age_Hepatocyte_BySex.png"),
       pD, width = 12, height = 8, dpi = 300)

message("    ✅ Panel D saved  [range: ", round(cs_D$lo*100),
        "% – ", round(cs_D$hi*100), "%, midpoint: ", round(cs_D$mid*100), "%]")


## =========================================================
## DONE
## =========================================================

cat("\n========================================\n")
cat("✅ 4 plots saved (PDF + PNG):\n")
cat("  1. PanelA_AgeDomain_Proportion_Chr\n")
cat("  2. PanelB_DotPlot_Chr_Celltype_BySex\n")
cat("  3. PanelC_DotPlot_Chr_Age_AllCelltypes_BySex\n")
cat("  4. PanelD_DotPlot_Chr_Age_Hepatocyte_BySex\n")
cat("\nColor scales (data-driven, midpoint = center of data range):\n")
cat(sprintf("  Panel B: %d%% – %d%% (mid=%d%%)\n", round(cs_B$lo*100), round(cs_B$hi*100), round(cs_B$mid*100)))
cat(sprintf("  Panel C: %d%% – %d%% (mid=%d%%)\n", round(cs_C$lo*100), round(cs_C$hi*100), round(cs_C$mid*100)))
cat(sprintf("  Panel D: %d%% – %d%% (mid=%d%%)\n", round(cs_D$lo*100), round(cs_D$hi*100), round(cs_D$mid*100)))
cat("========================================\n")

#!/usr/bin/env Rscript
## =========================================================
## EXPORT: One Excel file with 4 tabs (chr1-19, no chrX)
##
## Tab 1: PanelA_AgeDomain       — age-domain proportion per chr
## Tab 2: PanelB_Chr_Celltype    — fraction active per chr × celltype × sex
## Tab 3: PanelC_Chr_Age_All     — fraction active per chr × age × sex (all celltypes)
## Tab 4: PanelD_Chr_Age_Hep     — fraction active per chr × age × sex (Hepatocyte)
## =========================================================

suppressPackageStartupMessages({
    library(data.table)
    library(openxlsx)
})

outdir <- "/data/sarkern2/multiome_liver/Seurat/epigenome/compartment_main_heatmap"

comp_bin      <- readRDS(file.path(outdir, "compartments_binary_A1_R0.rds"))
manifest      <- fread(file.path(outdir, "manifest.normalized.tsv"))
bins_df       <- fread(file.path(outdir, "bins_table.tsv"))
stability_dt  <- fread(file.path(outdir, "compartment_stability_5class.tsv"))

age_vec   <- setNames(manifest$age, manifest$group)
sex_vec   <- setNames(manifest$sex, manifest$group)
ctype_vec <- setNames(manifest$celltype, manifest$group)

age_levels <- c("young", "mid_age", "old", "pre_geriatric", "geriatric")
chr_order  <- paste0("chr", 1:19)

switching_classes <- c("Monotonic_A_to_R", "Monotonic_R_to_A", "Non_Monotonic")

## =========================================================
## BUILD LONG-FORMAT DATA
## =========================================================

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
all_data <- all_data[chr %in% chr_order]


## =========================================================
## TAB 1: Panel A — Age-domain proportion per chromosome
## =========================================================

message("[Tab 1] Panel A: Age-domain proportion...")

bin_stability <- stability_dt[, .(
    is_age_domain = any(Stability %in% switching_classes)
), by = bin_id]
bin_stability <- merge(bin_stability, bins_df[, .(bin_id, chr)], by = "bin_id")
bin_stability <- bin_stability[chr %in% chr_order]

tab1 <- bin_stability[, .(
    n_total            = .N,
    n_age_domain       = sum(is_age_domain),
    n_non_age_domain   = sum(!is_age_domain),
    pct_age_domain     = round(sum(is_age_domain) / .N * 100, 2),
    pct_non_age_domain = round(sum(!is_age_domain) / .N * 100, 2)
), by = chr]
tab1[, chr := factor(chr, levels = chr_order)]
tab1 <- tab1[order(-pct_age_domain)]


## =========================================================
## TAB 2: Panel B — Fraction active per chr × celltype × sex
## =========================================================

message("[Tab 2] Panel B: Chr × Celltype × Sex...")

tab2 <- all_data[, .(
    pct_active = round(mean(state == 1L, na.rm = TRUE) * 100, 2),
    n_bins     = .N
), by = .(chr, celltype, sex)]
tab2[, chr := factor(chr, levels = chr_order)]
tab2 <- tab2[order(chr, celltype, sex)]


## =========================================================
## TAB 3: Panel C — Fraction active per chr × age × sex (all celltypes)
## =========================================================

message("[Tab 3] Panel C: Chr × Age × Sex (all celltypes)...")

tab3 <- all_data[, .(
    pct_active = round(mean(state == 1L, na.rm = TRUE) * 100, 2),
    n_bins     = .N
), by = .(chr, age, sex)]
tab3[, chr := factor(chr, levels = chr_order)]
tab3[, age := factor(age, levels = age_levels)]
tab3 <- tab3[order(chr, age, sex)]


## =========================================================
## TAB 4: Panel D — Fraction active per chr × age × sex (Hepatocyte)
## =========================================================

message("[Tab 4] Panel D: Chr × Age × Sex (Hepatocyte)...")

tab4 <- all_data[celltype == "Hepatocyte", .(
    pct_active = round(mean(state == 1L, na.rm = TRUE) * 100, 2),
    n_bins     = .N
), by = .(chr, age, sex)]
tab4[, chr := factor(chr, levels = chr_order)]
tab4[, age := factor(age, levels = age_levels)]
tab4 <- tab4[order(chr, age, sex)]


## =========================================================
## WRITE SINGLE EXCEL FILE
## =========================================================

message("[EXPORT] Writing Excel file...")

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
setColWidths(wb, "PanelB_Chr_Celltype", cols = 1:4, widths = "auto")
setColWidths(wb, "PanelC_Chr_Age_All",  cols = 1:4, widths = "auto")
setColWidths(wb, "PanelD_Chr_Age_Hep",  cols = 1:4, widths = "auto")

outfile <- file.path(outdir, "compartment_dotplot_data.xlsx")
saveWorkbook(wb, outfile, overwrite = TRUE)

message("\n✅ Saved: ", outfile)
message("\nTabs:")
message("  1. PanelA_AgeDomain      — chr, n_total, n_age_domain, pct_age_domain")
message("  2. PanelB_Chr_Celltype   — chr, celltype, sex, pct_active, n_bins")
message("  3. PanelC_Chr_Age_All    — chr, age, sex, pct_active, n_bins")
message("  4. PanelD_Chr_Age_Hep    — chr, age, sex, pct_active, n_bins")
