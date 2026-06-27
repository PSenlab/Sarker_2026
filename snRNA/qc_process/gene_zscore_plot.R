library(ggplot2); library(dplyr)

# ── endo orders ──
celltype_order <- c("LSEC","LSEC_cycling","Macrovascular_portal",
                    "Macrovascular_central","endo_02")

gene_order <- c("Stab2","Clec4g","Gata4","Oit3","Dnase1l3","Mrc1",   # Pan-LSEC
                "Mki67","Top2a","Cenpf",                              # Cycling
                "Vwf","Sox17","Efnb2",                                # Portal
                "Rspo3","Wnt9b","Wnt2","Thbd",                        # Central
                "Clec4f","Vsig4","Timd4","Cd5l","C1qa","Marco","Csf1r") # KC2/myeloid

# ── helpers ──
clamp20_100 <- function(x) pmin(pmax(x, 20), 100)
size_scale_20_100 <- scale_size_continuous(
  range  = c(2, 7),
  limits = c(20, 100),
  breaks = seq(20, 100, by = 20),
  name   = "fraction of cells in group (%)"
)

# ── load endo DE (LogFC) + expr stats (PctPos) exported from python ──
de_path <- "/data/sarkern2/multiome_liver/endo_DE_csvs"
de <- bind_rows(lapply(celltype_order, function(ct){
  d <- read.csv(file.path(de_path, paste0("DE_Cluster_", ct, ".csv")))
  d$celltype <- ct; d
}))
expr_stats <- read.delim(file.path(de_path, "endo_expr_stats.tsv"))

# build plot_df: Zscore per gene across groups + PctPos
plot_df <- de %>%
  select(gene, celltype, LogFC = avg_log2FC) %>%
  left_join(expr_stats, by = c("gene","celltype")) %>%
  filter(gene %in% gene_order) %>%
  group_by(gene) %>%
  mutate(Zscore = (LogFC - mean(LogFC)) / (sd(LogFC) + 1e-9)) %>%   # z across groups
  ungroup() %>%
  rename(PctPos = pct_expressed)

# ── ATAC/RNA panel ──
plot_df_endo <- plot_df %>%
  mutate(
    celltype    = factor(celltype, levels = celltype_order),
    gene        = factor(gene,     levels = gene_order),
    PctPos_plot = clamp20_100(PctPos)
  )

p_endo <- ggplot(plot_df_endo, aes(x = gene, y = celltype)) +
  geom_point(aes(size = PctPos_plot, color = Zscore)) +
  size_scale_20_100 +
  scale_color_gradient2(low = "blue", mid = "white", high = "red",
                        midpoint = 0, name = "Z-score") +
  scale_y_discrete(limits = rev(celltype_order)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
        axis.text.y = element_text(size = 9),
        axis.title  = element_blank()) +
  ggtitle("Endo subclusters — size: % cells (20–100), color: Z-score")

ggsave("/data/sarkern2/multiome_liver/endo_DE_csvs/PanelC_endo_bubble_size20_100.pdf",
       p_endo, width = 8, height = 3, limitsize = FALSE)
