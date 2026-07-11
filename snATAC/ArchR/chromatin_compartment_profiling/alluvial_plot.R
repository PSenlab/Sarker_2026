#!/usr/bin/env Rscript
# ==============================================================================
# Panel g : compartment switching alluvial  (TRAJECTORY-BASED, matches manuscript)
# ==============================================================================
# Reconstructed from the trajectory-based method used for the paper (not the
# broken interaction(Stability, State) stub in the master's Step 17).
#
# Reads ONLY files the master already wrote to reviewer/:
#   compartments_binary_A1_R0.rds, compartment_stability_5class.tsv,
#   manifest.normalized.tsv
# ==============================================================================

suppressPackageStartupMessages({
    library(data.table); library(ggplot2); library(ggalluvial)
})

base   <- "/data/sarkern2/multiome_liver/Seurat/epigenome/reviewer"
outdir <- base

AGE_LEVELS        <- c("young", "mid_age", "old", "pre_geriatric", "geriatric")
SWITCHING_CLASSES <- c("Monotonic_A_to_R", "Monotonic_R_to_A", "Non_Monotonic")

# switch-type colours + legend labels (match panel g)
STAB_COLS <- c("Monotonic_A_to_R" = "darkred",
               "Monotonic_R_to_A" = "darkblue",
               "Non_Monotonic"    = "#708238")
STAB_LABS <- c("Monotonic_A_to_R" = "monotonic AR",
               "Monotonic_R_to_A" = "monotonic RA",
               "Non_Monotonic"    = "non-monotonic")

plot_font <- tryCatch(
    if (requireNamespace("systemfonts", quietly = TRUE) &&
        "Arial" %in% systemfonts::system_fonts()$family) "Arial" else "sans",
    error = function(e) "sans")

# ---- load ----
comp_bin     <- readRDS(file.path(base, "compartments_binary_A1_R0.rds"))
manifest     <- fread(file.path(base, "manifest.normalized.tsv"))
stability_dt <- fread(file.path(base, "compartment_stability_5class.tsv"))

age_vec   <- setNames(as.character(manifest$age),      manifest$group)
sex_vec   <- setNames(as.character(manifest$sex),      manifest$group)
ctype_vec <- setNames(as.character(manifest$celltype), manifest$group)

# ---- state per age (per bin/sample) + attach stability ----
state_per_age <- rbindlist(lapply(colnames(comp_bin), function(g) data.table(
    bin_id   = rownames(comp_bin),
    State    = ifelse(comp_bin[, g] == 1L, "Active", "Repressive"),
    Age      = age_vec[g], Sex = sex_vec[g], Celltype = ctype_vec[g]
)))
state_per_age[, Age := factor(Age, levels = AGE_LEVELS)]

state_with_stability <- merge(
    state_per_age,
    stability_dt[, .(bin_id, Sex, Celltype, Stability)],
    by = c("bin_id", "Sex", "Celltype"), all.x = TRUE
)

# ==============================================================================
# TRAJECTORY-BASED PREP (the correct approach)
# ==============================================================================
switching_data <- state_with_stability[Stability %in% SWITCHING_CLASSES]

# wide: one column per age; majority vote if multiple samples per age
state_wide <- dcast(
    switching_data,
    bin_id + Sex + Celltype + Stability ~ Age,
    value.var     = "State",
    fun.aggregate = function(x) { tbl <- table(x); names(tbl)[which.max(tbl)] }
)
state_wide <- state_wide[complete.cases(state_wide)]

# trajectory id = full 5-age path (+ class/celltype/sex); N counts bins on it
agg <- state_wide[, .N, by = .(Sex, Celltype, Stability,
                               young, mid_age, old, pre_geriatric, geriatric)]
agg[, trajectory_id := paste(Sex, Celltype, Stability,
                             young, mid_age, old, pre_geriatric, geriatric,
                             sep = "_")]

# back to long for ggalluvial
agg_long <- melt(agg,
    id.vars     = c("N", "Sex", "Celltype", "Stability", "trajectory_id"),
    measure.vars = AGE_LEVELS, variable.name = "Age", value.name = "State")
agg_long[, Age := factor(Age, levels = AGE_LEVELS)]
agg_long <- agg_long[!is.na(State)]
agg_long[, Stability := factor(Stability, levels = SWITCHING_CLASSES)]
agg_long[, Sex := factor(Sex, levels = c("male", "female"))]

# ==============================================================================
# PLOT panel g
# ==============================================================================
p_g <- ggplot(agg_long,
       aes(x = Age, stratum = State, alluvium = trajectory_id,
           y = N, fill = Stability)) +
    stat_alluvium(geom = "flow", alpha = 0.7, width = 0.4,
                  curve_type = "linear", na.rm = TRUE) +
    stat_stratum(width = 0.4, alpha = 0.9, colour = "black",
                 linewidth = 0.2, fill = "white", na.rm = TRUE) +
    geom_text(stat = "stratum", aes(label = after_stat(stratum)),
              angle = 90, size = 3, na.rm = TRUE) +
    scale_fill_manual(values = STAB_COLS, labels = STAB_LABS,
                      breaks = SWITCHING_CLASSES, name = "switch type") +
    facet_wrap(~ Sex, ncol = 1) +
    labs(x = NULL, y = "number of 80kb bins") +
    theme_bw(base_size = 12) +
    theme(text = element_text(family = plot_font),
          axis.text.x = element_text(angle = 45, hjust = 1),
          strip.text = element_text(hjust = 0),
          legend.position = "bottom")

ggsave(file.path(outdir, "Alluvial_panelG.pdf"), p_g, width = 7, height = 10, device = cairo_pdf)
ggsave(file.path(outdir, "Alluvial_panelG.png"), p_g, width = 7, height = 10, dpi = 300)

cat(sprintf("Trajectories: %d | long rows: %d\n", nrow(agg), nrow(agg_long)))
cat("Saved Alluvial_panelG.pdf/png\n")
