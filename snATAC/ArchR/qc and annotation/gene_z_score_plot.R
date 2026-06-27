library(ArchR); library(dplyr); library(ggplot2)

# ── orders ──
celltype_order <- c("LSEC","LSEC_cycling","Macrovascular_portal",
                    "Macrovascular_central","endo_02")
gene_order <- c("Stab2","Clec4g","Gata4","Oit3","Dnase1l3","Mrc1",
                "Mki67","Top2a","Cenpf","Vwf","Sox17","Efnb2",
                "Rspo3","Wnt9b","Wnt2","Thbd",
                "Clec4f","Vsig4","Timd4","Cd5l","C1qa","Marco","Csf1r")

# ── load the plot_df we already built from getMarkerFeatures ──
# (this has gene, group, Log2FC, Mean, PctPos — the bias-corrected version)
plot_df <- read.csv("/data/sarkern2/multiome_liver/atac_endo_figs/plot_df_endo_sub.csv",
                    stringsAsFactors = FALSE)
colnames(plot_df)[colnames(plot_df) == "group"] <- "celltype"

# ── color variable = Log2FC (bias-corrected), z-scored per gene across celltypes ──
plot_df <- plot_df %>%
  filter(gene %in% gene_order) %>%
  group_by(gene) %>%
  mutate(Zscore = as.numeric(scale(Log2FC))) %>%   # <- z-score of Log2FC, not raw Mean
  ungroup() %>%
  mutate(celltype = factor(celltype, levels = celltype_order),
         gene     = factor(gene,     levels = gene_order))

# ── PanelC style ──
clamp20_100 <- function(x) pmin(pmax(x, 20), 100)
size_scale_20_100 <- scale_size_continuous(
  range = c(2, 7), limits = c(20, 100), breaks = seq(20, 100, by = 20),
  name = "fraction of cells in group (%)")
plot_df$PctPos_plot <- clamp20_100(plot_df$PctPos)

p_atac <- ggplot(plot_df, aes(x = gene, y = celltype)) +
  geom_point(aes(size = PctPos_plot, color = Zscore)) +
  size_scale_20_100 +
  scale_color_gradient2(low = "blue", mid = "white", high = "red",
                        midpoint = 0, name = "Z-score") +
  scale_y_discrete(limits = rev(celltype_order)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9, face = "italic"),
        axis.text.y = element_text(size = 9),
        axis.title  = element_blank()) +
  ggtitle("Endo ATAC GeneScore markers — size: % cells (20–100), color: Z-score")

ggsave("/data/sarkern2/multiome_liver/endo_DE_csvs/PanelC_endo_ATAC_bubble.pdf",
       p_atac, width = 11, height = 4, limitsize = FALSE)
