#!/usr/bin/env Rscript
#===============================================================================
# DOTPLOT: liver-trait TRS across age, gaining vs losing seeding
#   Reconstruction of the two-panel figure. Reads the SCAVENGE combined_ALL_raw
#   outputs from both direction folders, aggregates per trait x age x direction,
#   and draws traits (y) x age (x), faceted by category (rows) and direction
#   (columns). Dot size = % FDR-significant cells; fill = median TRS.
#
#   NOTE: this figure is a visualization of the cell-level SCAVENGE output only.
#   It is not a test of disease specificity (see analysis notes).
#===============================================================================

suppressPackageStartupMessages({ library(data.table); library(ggplot2) })

base       <- "/data/sarkern2/multiome_liver/Seurat/epigenome/gwas_peak_gene"
age_levels <- c("young","mid_age","old","pre_geriatric","geriatric")

#-------------------------------------------------------------------------------
# category scheme (matches the paper dotplot)
#-------------------------------------------------------------------------------
category_patterns <- list(
  "Liver Enzymes" = c("\\bALT\\b","\\bAST\\b","\\bGGT\\b","\\bALP\\b","alanine aminotransferase",
    "aspartate aminotransferase","gamma.glutamyl","gamma glutamyl","gamma-glutamyl","alkaline phosphatase"),
  "Bilirubin" = c("bilirubin"),
  "Albumin/Globulins" = c("albumin","globulin","non-albumin protein"),
  "Complement/Coagulation" = c("complement C","complement component","fibrinogen","prothrombin",
    "coagulation factor","factor VII","factor VIII","factor IX","factor X"),
  "Apolipoproteins" = c("apolipoprotein","\\bAPO[ABCDEFHLM]\\b","\\bAPOM\\b"),
  "VLDL/Lipoproteins" = c("\\bVLDL\\b","lipoprotein"),
  "Iron Metabolism" = c("ferritin","transferrin","hemochromatosis","hepcidin","iron\\b","\\biron"),
  "Other Liver Proteins" = c("ceruloplasmin","haptoglobin","alpha.1.antitrypsin","alpha-1-antitrypsin",
    "angiotensinogen","protein quantitative trait"),
  "NAFLD/Steatosis" = c("NAFLD","NASH","steatosis","steatohepatitis","fatty liver","hepatic fat","liver fat"),
  "Fibrosis/Cirrhosis" = c("fibrosis","cirrhosis","liver condition"),
  "Liver Cancer" = c("hepatocellular carcinoma","hepatic cancer","liver cancer","cholangiocarcinoma"),
  "Viral Hepatitis" = c("hepatitis B","hepatitis C","hepatitis E","\\bHBV\\b","\\bHCV\\b","anti-hepatitis"),
  "Autoimmune Liver" = c("autoimmune hepatitis","primary biliary","primary sclerosing","\\bPBC\\b","\\bPSC\\b"),
  "Biliary/Cholestasis" = c("biliary","bile ","cholestasis","cholangitis","gallbladder"),
  "Drug-Induced" = c("drug.induced","hepatotoxicity","alcohol"),
  "Pregnancy-Related" = c("pregnancy","childbirth","puerperium","neonatal"),
  "Metabolic Syndrome" = c("metabolic syndrome","obesity","BMI","diabetes")
)
categorize_trait <- function(trait, patterns_list) {
  if (is.na(trait) || trait == "") return("Other")
  tl <- tolower(trait)
  for (cat in names(patterns_list))
    for (p in patterns_list[[cat]])
      if (grepl(tolower(p), tl, ignore.case = TRUE)) return(cat)
  "Other"
}

#-------------------------------------------------------------------------------
# load both directions, aggregate per trait x age
#-------------------------------------------------------------------------------
load_dir <- function(dir, lab) {
  r <- fread(file.path(base, dir, "gwas_p2g/Liver/results/combined_ALL_raw.csv"))
  a <- r[, .(median_TRS  = median(TRS, na.rm = TRUE),
             pct_fdr_sig = 100 * mean(fdr < 0.05, na.rm = TRUE),
             n_cells     = .N),
         by = .(trait, age)]
  a$direction <- lab; a
}
agg <- rbind(load_dir("hepatocyte_gaining","gaining"),
             load_dir("hepatocyte_losing","losing"))
agg <- agg[n_cells >= 20 & is.finite(median_TRS)]
agg$category <- sapply(agg$trait, function(t) categorize_trait(t, category_patterns))
agg <- agg[!is.na(category) & category != "Other"]

cat_order <- c("Liver Enzymes","Bilirubin","Albumin/Globulins","Complement/Coagulation",
               "Apolipoproteins","VLDL/Lipoproteins","Iron Metabolism","Other Liver Proteins",
               "NAFLD/Steatosis","Fibrosis/Cirrhosis","Liver Cancer","Viral Hepatitis",
               "Autoimmune Liver","Biliary/Cholestasis","Drug-Induced","Pregnancy-Related",
               "Metabolic Syndrome")
agg$category  <- factor(agg$category,  levels = cat_order)
agg$age       <- factor(agg$age,       levels = age_levels)
agg$direction <- factor(agg$direction, levels = c("gaining","losing"))

#-------------------------------------------------------------------------------
# plot
#-------------------------------------------------------------------------------
p <- ggplot(agg, aes(x = age, y = trait)) +
  geom_point(aes(size = pct_fdr_sig, fill = median_TRS),
             shape = 21, color = "grey40", stroke = 0.2) +
  scale_fill_gradientn(colours = c("#fff5f0","#fcbba1","#fb9a99","#cb181d"),
                       name = "median TRS") +
  scale_size_continuous(range = c(1, 6), name = "% cells FDR < 0.05") +
  facet_grid(category ~ direction, scales = "free_y", space = "free_y") +
  labs(x = NULL, y = NULL, title = "Liver-trait TRS across age: gaining vs losing seeding") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        axis.text.y = element_text(size = 5),
        strip.text.y = element_text(angle = 0, size = 7, face = "bold"),
        strip.text.x = element_text(size = 9, face = "bold"),
        strip.background = element_rect(fill = "grey90"),
        panel.spacing = unit(0.15, "lines"))

ggsave(file.path(base, "dotplot_TRS_by_age_gain_vs_lose.pdf"), p,
       width = 11, height = max(15, length(unique(agg$trait)) * 0.16), limitsize = FALSE)
cat("saved dotplot;", length(unique(agg$trait)), "traits,",
    nrow(agg), "trait x age x direction points\n")
