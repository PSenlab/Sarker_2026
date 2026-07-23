# =============================================================================
# Hepatocyte transcriptional entropy x senescence module score — alluvial
#
# Module score : 42-gene panel = RNA C2 DE cluster ∩ SHGS
# Bins         : equal-width thirds of [-1,1]-rescaled scores
#                cut points computed on ALL hepatocytes pooled (global ruler)
#
# Sarker et al. — liver aging multiome
# =============================================================================

library(Seurat)
library(dplyr)
library(ggalluvial)
library(ggplot2)
library(patchwork)

# ---- config -----------------------------------------------------------------
IN_RDS     <- "seurat_with_entropy_merged_rep.rds"
AGE_LEVELS <- c("young", "mid_age", "old", "pre_geriatric", "geriatric")

SEX_PAL <- c(female = "#fb9a99", male = "#89d9e1")
ENT_PAL <- c("E-Low" = "#d9e8de", "E-Mid" = "#6ba585", "E-High" = "#1d5138")

CAP_SHARED <- paste0(
  "Equal-width thirds of [\u22121,1]-rescaled scores, ",
  "cut points computed on all hepatocytes pooled (identical across sex and age) \u00b7 ",
  "each panel normalised to 100%")

set.seed(1)

# ---- 1. load and subset to hepatocytes --------------------------------------
seu <- readRDS(IN_RDS)
hep <- subset(seu, subset = celltype == "Hepatocyte")
DefaultAssay(hep) <- "RNA"

table(hep$sex, hep$age)

# ---- 2. module score --------------------------------------------------------
# 42 genes: intersection of the RNA C2 up-regulated DE cluster with the
# senescent hepatocyte gene signature (SHGS).
gene_list_combined <- c(
  "Abhd2", "Arl4c", "Bmp8b", "Capn2", "Cd151", "Cd44", "Cdkn2b",
  "Cxcl9", "Ddit4l", "Dnajc10", "Elovl7", "Enc1", "Flna", "Fosl1",
  "Gbp2", "Gls", "Ier3", "Ifngr1", "Igdcc4", "Inpp4b", "Iqgap1",
  "Itgav", "Lgals3", "Ltb", "Myof", "Panx1", "Ppfibp1", "Rab31",
  "Rhbdf1", "S100a6", "Serinc2", "Sh3bgrl3", "Slc48a1", "Slpi",
  "Srxn1", "Sytl5", "Tagln2", "Tanc1", "Tgfbr2", "Tnfrsf12a",
  "Trim35", "Vim")

genes_present <- intersect(gene_list_combined, rownames(hep))
message(length(genes_present), "/", length(gene_list_combined), " genes found")
if (length(genes_present) < length(gene_list_combined))
  message("missing: ", paste(setdiff(gene_list_combined, genes_present), collapse = ", "))

hep <- AddModuleScore(hep, features = list(genes_present),
                      name = "Combined_ModuleScore")

# ---- 3. global [-1,1] rescale + equal-width tertile bins ---------------------
rescale11 <- function(x) {
  rng <- range(x, na.rm = TRUE)
  2 * (x - rng[1]) / (rng[2] - rng[1]) - 1
}
eqbin <- function(x) cut(x, breaks = c(-1, -1/3, 1/3, 1),
                         labels = c("Low", "Mid", "High"), include.lowest = TRUE)

af <- hep@meta.data %>%
  transmute(sample, sex, age,
            entropy = rescale11(entropy_score),
            module  = rescale11(Combined_ModuleScore1)) %>%
  mutate(entropy_bin = factor(paste0("E-", eqbin(entropy)),
                              levels = c("E-Low", "E-Mid", "E-High")),
         module_bin  = factor(paste0("M-", eqbin(module)),
                              levels = c("M-Low", "M-Mid", "M-High")),
         age = factor(age, levels = AGE_LEVELS))

flow <- af %>% count(age, entropy_bin, module_bin, sex, name = "n")

# ---- 4. plotting helper -----------------------------------------------------
make_alluvial <- function(dat, fill_var, palette, fill_name, cap,
                          title = NULL, show_n = FALSE) {

  dat <- dat %>%
    group_by(age) %>% mutate(frac = n / sum(n)) %>% ungroup() %>%
    mutate(age = factor(age, levels = AGE_LEVELS))

  p <- ggplot(dat, aes(axis1 = entropy_bin, axis2 = module_bin, y = frac)) +
    geom_alluvium(aes(fill = .data[[fill_var]]), alpha = 0.65, width = 0.2,
                  color = "white", linewidth = 0.12, curve_type = "sigmoid") +
    geom_stratum(width = 0.2, fill = "grey92", color = "grey45",
                 linewidth = 0.25) +
    geom_text(stat = "stratum",
              aes(label = sub("^._", "", after_stat(stratum))),
              family = "Arial", size = 2.6) +
    facet_wrap(~ age, nrow = 1) +
    scale_x_discrete(limits = c("Entropy", "Module"), expand = c(0.15, 0.15),
                     position = "top") +
    scale_fill_manual(values = palette, name = fill_name) +
    scale_y_continuous(labels = scales::percent,
                       expand = expansion(mult = c(0, 0.08))) +
    coord_cartesian(clip = "off") +
    labs(title = title, y = "Proportion of hepatocytes", caption = cap) +
    theme_minimal(base_family = "Arial") +
    theme(panel.grid    = element_blank(),
          axis.title.x  = element_blank(),
          axis.text.x   = element_text(size = 9, face = "bold"),
          axis.text.y   = element_text(size = 7, colour = "grey45"),
          strip.text    = element_text(size = 10, face = "bold"),
          plot.title    = element_text(size = 11, face = "bold"),
          panel.spacing = unit(1, "lines"),
          plot.caption  = element_text(size = 7.5, colour = "grey50", hjust = 0))

  if (show_n) {
    pn <- dat %>% count(age, wt = n, name = "N") %>%
      mutate(age = factor(age, levels = AGE_LEVELS),
             lab = paste0("n = ", scales::comma(N)))
    p <- p + geom_text(data = pn, inherit.aes = FALSE,
                       aes(x = 1.5, y = 1.045, label = lab),
                       family = "Arial", size = 2.4, colour = "grey45")
  }
  p
}

# ---- 5. figures -------------------------------------------------------------
# combined, both sexes, sex as fill
p_prop <- make_alluvial(flow, "sex", SEX_PAL, "Sex", CAP_SHARED)

# Variant A — per sex, sex palette (consistent with combined)
pA_female <- make_alluvial(filter(flow, sex == "female"), "sex", SEX_PAL,
                           "Sex", CAP_SHARED, "Female")
pA_male   <- make_alluvial(filter(flow, sex == "male"),   "sex", SEX_PAL,
                           "Sex", CAP_SHARED, "Male")
pA_stack  <- (pA_female + labs(caption = NULL)) / pA_male +
  plot_layout(guides = "collect") & theme(legend.position = "right")

# Variant B — per sex, entropy-tertile palette (flows traceable by origin)
pB_female <- make_alluvial(filter(flow, sex == "female"), "entropy_bin",
                           ENT_PAL, "Entropy tertile", CAP_SHARED, "Female")
pB_male   <- make_alluvial(filter(flow, sex == "male"),   "entropy_bin",
                           ENT_PAL, "Entropy tertile", CAP_SHARED, "Male")
pB_stack  <- (pB_female + labs(caption = NULL)) / pB_male +
  plot_layout(guides = "collect") & theme(legend.position = "right")

# ---- 6. export --------------------------------------------------------------
ggsave("alluvial_pm1_byage_prop_combined.pdf", p_prop,    width = 14, height = 5, device = cairo_pdf)
ggsave("alluvial_sexpal_female.pdf",           pA_female, width = 14, height = 5, device = cairo_pdf)
ggsave("alluvial_sexpal_male.pdf",             pA_male,   width = 14, height = 5, device = cairo_pdf)
ggsave("alluvial_sexpal_stack.pdf",            pA_stack,  width = 14, height = 9, device = cairo_pdf)
ggsave("alluvial_entpal_female.pdf",           pB_female, width = 14, height = 5, device = cairo_pdf)
ggsave("alluvial_entpal_male.pdf",             pB_male,   width = 14, height = 5, device = cairo_pdf)
ggsave("alluvial_entpal_stack.pdf",            pB_stack,  width = 14, height = 9, device = cairo_pdf)

# ---- 7. source data ---------------------------------------------------------
# marginals actually underlying the figures (from `af`, not the old `df`)
af %>% count(age, sex, entropy_bin) %>%
  group_by(age, sex) %>% mutate(pct = round(100 * n / sum(n), 1)) %>%
  ungroup() %>% arrange(entropy_bin, age, sex) %>%
  readr::write_csv("source_data_entropy_bin_marginals.csv")

flow %>% readr::write_csv("source_data_alluvial_flow.csv")

sessionInfo()
