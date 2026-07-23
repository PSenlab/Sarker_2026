# =============================================================================
# Hepatocyte transcriptional entropy x nuclear ploidy — alluvial figures
#
# Entropy : whole-transcriptome, RNA_entropy assay
# Ploidy  : scploidy ATAC-based per-cell calls (2n / 4n / 8n)
# Bins    : global tertiles of the [-1,1]-rescaled entropy score
#           (ploidy is natively categorical and is not binned)
#
# Sarker et al. — liver aging multiome
# =============================================================================

library(Seurat)
library(dplyr)
library(stringr)
library(ggalluvial)
library(ggplot2)
library(lmerTest)

# ---- config -----------------------------------------------------------------
IN_RDS     <- "seurat_with_entropy_merged_rep.rds"
WORK       <- "/data/sarkern2/multiome_liver/scploidy"
PLOIDY_CSV <- file.path(WORK, "aggregate", "ploidy_annotated_percell_CLEAN.csv")

AGE_LEVELS <- c("young", "mid_age", "old", "pre_geriatric", "geriatric")


SEX_PAL <- c(female = "#fb9a99", male = "#89d9e1")
ENT_PAL <- c("E-Low" = "#9a7fc4", "E-Mid" = "#6ba585", "E-High" = "#1d5138")

# ---- 1. load and subset to hepatocytes --------------------------------------
seu <- readRDS(IN_RDS)
hep <- subset(seu, subset = celltype == "Hepatocyte")
DefaultAssay(hep) <- "RNA"

# ---- 2. attach scploidy calls -----------------------------------------------
# scploidy cellkey = "<sample>#<barcode>"; Seurat = "<barcode>-<sample>",
# with leading zeros stripped from the sample index (geriatric_05 -> geriatric_5)
ploidy <- readr::read_csv(PLOIDY_CSV, show_col_types = FALSE) %>%
  mutate(samp_raw   = str_extract(cellkey, "^[^#]+"),
         bc         = str_extract(cellkey, "(?<=#)[ACGT]+"),
         samp       = str_replace(samp_raw, "_0*(\\d+)$", "_\\1"),
         seurat_key = paste0(bc, "-", samp))

message(sum(ploidy$seurat_key %in% colnames(hep)), " / ", ncol(hep),
        " hepatocytes matched to a ploidy call")

idx <- match(colnames(hep), ploidy$seurat_key)
hep$ploidy_call <- ploidy$ploidy_call[idx]
hep$nFrags      <- ploidy$nFrags[idx]
hep$has_ploidy  <- !is.na(hep$ploidy_call)

# ---- 3. dropout diagnostics -------------------------------------------------
# is the ~10% of hepatocytes without a ploidy call uniform across groups?
hep@meta.data %>%
  group_by(age, sex) %>%
  summarise(pct_matched = round(100 * mean(has_ploidy), 1), n = n(),
            .groups = "drop") %>%
  arrange(factor(age, levels = AGE_LEVELS), sex) %>%
  print(n = 20)

# do unmatched cells differ technically?
hep@meta.data %>%
  group_by(has_ploidy) %>%
  summarise(nFeature    = mean(nFeature_RNA),
            nCount_ATAC = mean(nCount_ATAC, na.rm = TRUE),
            entropy     = mean(entropy_score),
            n = n())

# entropy / depth by ploidy class
hep@meta.data %>%
  filter(has_ploidy) %>%
  group_by(ploidy_call) %>%
  summarise(mean_nFeature = mean(nFeature_RNA),
            mean_nCount   = mean(nCount_RNA),
            mean_entropy  = mean(entropy_score),
            n = n())

# ---- 4. does the sex effect on entropy survive adjustment for ploidy? --------
sex_by_age <- purrr::map_dfr(AGE_LEVELS, function(a) {
  d <- filter(hep@meta.data, age == a, has_ploidy)
  m <- lmerTest::lmer(entropy_score ~ sex + factor(ploidy_call) +
                        scale(nFeature_RNA) + scale(nCount_RNA) + (1 | sample),
                      data = d)
  s <- summary(m)$coefficients["sexmale", ]
  tibble::tibble(age = a, beta = s[1], se = s[2], df = s[3], p = s[5])
})
print(sex_by_age)

# young only: is a sex x ploidy interaction warranted?
yng    <- filter(hep@meta.data, age == "young", has_ploidy)
m_full <- lmer(entropy_score ~ sex + factor(ploidy_call) + scale(nFeature_RNA) +
                 scale(nCount_RNA) + (1 | sample), data = yng)
m_int  <- lmer(entropy_score ~ sex * factor(ploidy_call) + scale(nFeature_RNA) +
                 scale(nCount_RNA) + (1 | sample), data = yng)
print(anova(m_int, m_full))

# ---- 5. binning -------------------------------------------------------------
rescale11 <- function(x) {
  rng <- range(x, na.rm = TRUE); 2 * (x - rng[1]) / (rng[2] - rng[1]) - 1
}

pdat <- hep@meta.data %>%
  filter(has_ploidy) %>%
  transmute(sex     = factor(sex, levels = c("female", "male")),
            age     = factor(age, levels = AGE_LEVELS),
            entropy = rescale11(entropy_score),
            ploidy  = factor(paste0(ploidy_call, "n"),
                             levels = c("2n", "4n", "8n")))

pdat <- pdat %>%
  mutate(entropy_bin = cut(entropy, quantile(entropy, 0:3/3, na.rm = TRUE),
                           labels = c("E-Low", "E-Mid", "E-High"),
                           include.lowest = TRUE))

CAP_PL <- paste0("Entropy: global tertiles of the [\u22121,1]-rescaled score, ",
                 "cut points computed on all ploidy-called hepatocytes pooled \u00b7 ",
                 "ploidy: scploidy ATAC-based calls \u00b7 ",
                 "each panel normalised to 100% \u00b7 ",
                 scales::comma(nrow(pdat)), " cells with ploidy calls")

# ---- 6. plotting helper -----------------------------------------------------
make_ploidy_alluvial <- function(dat, fill_var, palette, fill_name, cap,
                                 title = NULL, by_sex = FALSE) {

  dat <- if (by_sex) dat %>% group_by(age, sex) else dat %>% group_by(age)
  dat <- dat %>% mutate(frac = n / sum(n)) %>% ungroup()

  p <- ggplot(dat, aes(axis1 = entropy_bin, axis2 = ploidy, y = frac)) +
    geom_alluvium(aes(fill = .data[[fill_var]]), alpha = 0.65, width = 0.2,
                  color = "white", linewidth = 0.12, curve_type = "sigmoid") +
    geom_stratum(width = 0.2, fill = "grey92", color = "grey45",
                 linewidth = 0.25) +
    geom_text(stat = "stratum", aes(label = after_stat(stratum)),
              family = "Arial", size = 2.5) +
    scale_x_discrete(limits = c("Entropy", "Ploidy"), expand = c(0.15, 0.15),
                     position = "top") +
    scale_fill_manual(values = palette, name = fill_name) +
    scale_y_continuous(labels = scales::percent,
                       expand = expansion(mult = c(0, 0.08))) +
    coord_cartesian(clip = "off") +
    labs(title = title, y = "Proportion of hepatocytes", caption = cap) +
    theme_minimal(base_family = "Arial") +
    theme(panel.grid        = element_blank(),
          axis.title.x      = element_blank(),
          axis.text.x       = element_text(size = 9, face = "bold"),
          axis.text.y       = element_text(size = 7, colour = "grey45"),
          strip.text.x      = element_text(size = 10, face = "bold"),
          strip.text.y.left = element_text(size = 11, face = "bold", angle = 0),
          strip.placement   = "outside",
          plot.title        = element_text(size = 11, face = "bold"),
          panel.spacing     = unit(0.9, "lines"),
          plot.caption      = element_text(size = 7.5, colour = "grey50",
                                           hjust = 0))

  if (by_sex) p + facet_grid(sex ~ age, switch = "y")
  else        p + facet_wrap(~ age, nrow = 1)
}

# ---- 7. figures -------------------------------------------------------------
flow_pl     <- pdat %>% count(age, entropy_bin, ploidy, sex, name = "n")
flow_pl_sex <- pdat %>% count(age, sex, entropy_bin, ploidy, name = "n")

# combined: both sexes pooled per age panel, sex as fill
p_pl <- make_ploidy_alluvial(flow_pl, "sex", SEX_PAL, "Sex", CAP_PL)

# Variant A: sex on facet rows, sex palette (consistent with combined)
pA_pl <- make_ploidy_alluvial(
  flow_pl_sex, "sex", SEX_PAL, "Sex", CAP_PL,
  "Hepatocyte entropy and ploidy across age, by sex", by_sex = TRUE)

# Variant B: sex on facet rows, entropy palette (flows traceable by origin)
pB_pl <- make_ploidy_alluvial(
  flow_pl_sex, "entropy_bin", ENT_PAL, "Entropy stratum", CAP_PL,
  "Hepatocyte entropy and ploidy across age, by sex", by_sex = TRUE)

# ---- 8. export --------------------------------------------------------------
ggsave("alluvial_entropy_ploidy_byage.pdf",       p_pl,  width = 14, height = 5, device = cairo_pdf)
ggsave("alluvial_entropy_ploidy_sexpal_rows.pdf", pA_pl, width = 14, height = 8, device = cairo_pdf)
ggsave("alluvial_entropy_ploidy_entpal_rows.pdf", pB_pl, width = 14, height = 8, device = cairo_pdf)

# ---- 9. source data ---------------------------------------------------------
readr::write_csv(flow_pl_sex, "source_data_entropy_ploidy_flow.csv")
readr::write_csv(sex_by_age,  "source_data_sex_effect_by_age_ploidy_adjusted.csv")

sessionInfo()

print(p_pl); print(pA_pl); print(pB_pl)
