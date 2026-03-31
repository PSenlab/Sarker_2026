#!/usr/bin/env Rscript
## =========================================================
## COMPARTMENT ANALYSIS - ALLUVIAL PLOTS
## 
## All alluvial plot variants:
##   1. Per stability class - faceted by celltype
##   2. Per stability class - celltype-colored (merged)
##   3. Male/Female side-by-side comparisons
##   4. All stability classes combined
##   5. Switching bins only
##   6. Merged alluvials (all celltypes combined)
##
## Prerequisites:
##   - Source compartment_00_config.R
##   - Run compartment_01_main_analysis.R
##   - Run compartment_02_stability_analysis.R
##
## =========================================================

box_banner("ALLUVIAL PLOTS FOR COMPARTMENT DYNAMICS")

## =========================================================
## LOAD DEPENDENCIES
## =========================================================

if (!exists("CONFIG")) {
    source("compartment_00_config.R")
}

suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(ggalluvial)
    library(cowplot)
})

# Check required objects
check_required_objects(c("state_with_stability", "stability_ext_dt", "CONFIG"))


## =========================================================
## SECTION 1: VALIDATION FUNCTIONS
## =========================================================

#' Validate bin counts for alluvial plot
validate_alluvial_bins <- function(agg_data, stability_class, sex, expected_total = NULL) {
    actual_total <- sum(agg_data$N, na.rm = TRUE)
    message(sprintf("  [CHECK] %s %s: %d bins in plot", sex, stability_class, actual_total))
    
    if (!is.null(expected_total) && actual_total != expected_total) {
        warning(sprintf("    Expected %d bins, got %d", expected_total, actual_total))
        return(FALSE)
    }
    return(TRUE)
}

#' Comprehensive validation report
validate_bin_data <- function(state_data, stability_data) {
    message("\n[VALIDATION] Running comprehensive bin data checks...")
    
    # Check state values
    unique_states <- unique(state_data$State)
    valid_states <- c("Active", "Repressive")
    invalid_states <- setdiff(unique_states, valid_states)
    
    if (length(invalid_states) > 0) {
        warning(sprintf("Invalid states found: %s", paste(invalid_states, collapse = ", ")))
    }
    
    # Check age coverage
    ages_found <- unique(state_data$Age)
    missing_ages <- setdiff(AGE_LEVELS, ages_found)
    
    if (length(missing_ages) > 0) {
        warning(sprintf("Missing ages: %s", paste(missing_ages, collapse = ", ")))
    }
    
    # Summary
    message(sprintf("  Unique bins in state data: %d", length(unique(state_data$bin_id))))
    message(sprintf("  Celltypes: %d", length(unique(state_data$Celltype))))
    message("  ✅ Validation complete")
}


## =========================================================
## SECTION 2: ALLUVIAL PLOT FUNCTIONS
## =========================================================

#' Create alluvial plot for a specific stability class (faceted by celltype)
#' @param data Merged state+stability data
#' @param stability_class Which stability class to plot
#' @param sex_select Which sex ("male" or "female")
#' @param title_suffix Additional title text
#' @return ggplot object
create_stability_alluvial_by_celltype <- function(data, stability_class, sex_select,
                                                   title_suffix = "") {
    
    dat <- data[Sex == sex_select & Stability == stability_class]
    
    if (nrow(dat) == 0) {
        message(sprintf("  [SKIP] No data for %s %s", sex_select, stability_class))
        return(NULL)
    }
    
    # Wide format
    dat_wide <- dcast(dat, bin_id + Celltype ~ Age, value.var = "State", 
                      fun.aggregate = function(x) x[1])
    dat_wide <- dat_wide[complete.cases(dat_wide)]
    
    if (nrow(dat_wide) == 0) {
        message(sprintf("  [SKIP] No complete trajectories for %s %s", sex_select, stability_class))
        return(NULL)
    }
    
    # Aggregate by trajectory
    agg <- dat_wide[, .N, by = c("Celltype", AGE_LEVELS)]
    agg[, trajectory_id := do.call(paste, c(.SD, sep = "_")), .SDcols = c("Celltype", AGE_LEVELS)]
    
    # Validate
    validate_alluvial_bins(agg, stability_class, sex_select)
    
    # Long format for ggalluvial
    agg_long <- melt(agg,
                     id.vars = c("N", "Celltype", "trajectory_id"),
                     measure.vars = AGE_LEVELS,
                     variable.name = "Age",
                     value.name = "State")
    agg_long[, Age := factor(Age, levels = AGE_LEVELS)]
    agg_long <- agg_long[!is.na(State)]
    
    # Order celltypes by total bins
    ct_order <- agg_long[, .(total = sum(N)), by = Celltype][order(-total), Celltype]
    agg_long[, Celltype := factor(Celltype, levels = ct_order)]
    
    total_bins <- sum(agg$N)
    
    p <- ggplot(agg_long,
                aes(x = Age, stratum = State, alluvium = trajectory_id,
                    y = N, fill = State)) +
        stat_alluvium(geom = "flow", alpha = 0.7, width = 0.4,
                      curve_type = "linear", na.rm = TRUE) +
        stat_stratum(width = 0.4, alpha = 0.9, color = "black",
                     linewidth = 0.3, na.rm = TRUE) +
        geom_text(stat = "stratum", aes(label = State),
                  size = 2.5, fontface = "bold", na.rm = TRUE) +
        facet_wrap(~ Celltype, ncol = 3, scales = "free_y") +
        scale_fill_manual(values = STATE_COLS, na.translate = FALSE) +
        labs(
            title = sprintf("%s - %s%s", 
                            gsub("_", " ", stability_class),
                            tools::toTitleCase(sex_select),
                            title_suffix),
            subtitle = sprintf("Total: %s bins across %d celltypes", 
                               format(total_bins, big.mark = ","),
                               length(unique(agg_long$Celltype))),
            x = "Age",
            y = "Number of 80kb Bins",
            fill = "Compartment State"
        ) +
        theme_bw(base_size = 11) +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 9),
            panel.grid.minor = element_blank(),
            strip.text = element_text(size = 10, face = "bold"),
            strip.background = element_rect(fill = "gray95"),
            plot.title = element_text(face = "bold", size = 14),
            plot.subtitle = element_text(size = 10, color = "gray40"),
            legend.position = "bottom"
        )
    
    return(p)
}


#' Create alluvial colored by celltype (all celltypes merged)
#' @param data Merged state+stability data
#' @param stability_class Which stability class to plot
#' @param sex_select Which sex
#' @return ggplot object
create_celltype_colored_alluvial <- function(data, stability_class, sex_select) {
    
    dat <- data[Sex == sex_select & Stability == stability_class]
    
    if (nrow(dat) == 0) return(NULL)
    
    dat_wide <- dcast(dat, bin_id + Celltype ~ Age, value.var = "State",
                      fun.aggregate = function(x) x[1])
    dat_wide <- dat_wide[complete.cases(dat_wide)]
    
    if (nrow(dat_wide) == 0) return(NULL)
    
    agg <- dat_wide[, .N, by = c("Celltype", AGE_LEVELS)]
    agg[, trajectory_id := do.call(paste, c(.SD, sep = "_")), .SDcols = c("Celltype", AGE_LEVELS)]
    
    validate_alluvial_bins(agg, stability_class, sex_select)
    
    agg_long <- melt(agg,
                     id.vars = c("N", "Celltype", "trajectory_id"),
                     measure.vars = AGE_LEVELS,
                     variable.name = "Age",
                     value.name = "State")
    agg_long[, Age := factor(Age, levels = AGE_LEVELS)]
    agg_long <- agg_long[!is.na(State)]
    
    total_bins <- sum(agg$N)
    
    p <- ggplot(agg_long,
                aes(x = Age, stratum = State, alluvium = trajectory_id,
                    y = N, fill = Celltype)) +
        stat_alluvium(geom = "flow", alpha = 0.6, width = 0.4,
                      curve_type = "linear", na.rm = TRUE) +
        stat_stratum(aes(fill = Celltype), width = 0.4, alpha = 0.9, 
                     color = "black", linewidth = 0.3, na.rm = TRUE) +
        geom_text(stat = "stratum", aes(label = State),
                  size = 3, fontface = "bold", na.rm = TRUE) +
        scale_fill_manual(values = CELLTYPE_COLS, na.translate = FALSE) +
        labs(
            title = sprintf("%s - %s (Colored by Celltype)",
                            gsub("_", " ", stability_class),
                            tools::toTitleCase(sex_select)),
            subtitle = sprintf("Total: %s bins", format(total_bins, big.mark = ",")),
            x = "Age",
            y = "Number of 80kb Bins",
            fill = "Cell Type"
        ) +
        theme_bw(base_size = 12) +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
            panel.grid.minor = element_blank(),
            plot.title = element_text(face = "bold", size = 14),
            legend.position = "right"
        )
    
    return(p)
}


#' Create merged alluvial for a specific stability class (all celltypes merged)
#' @param data Merged state+stability data
#' @param stability_class Which stability class to plot
#' @param title_label Title for the plot
#' @return ggplot object
create_merged_alluvial <- function(data, stability_class, title_label) {
    
    dat <- data[Stability == stability_class]
    if (nrow(dat) == 0) return(NULL)
    
    dat_wide <- dcast(dat, bin_id + Sex ~ Age, value.var = "State", fun.aggregate = function(x) x[1])
    dat_wide <- dat_wide[complete.cases(dat_wide)]
    if (nrow(dat_wide) == 0) return(NULL)
    
    message(sprintf("[INFO] %s merged: %d bins", stability_class, nrow(dat_wide)))
    
    agg <- dat_wide[, .N, by = c("Sex", AGE_LEVELS)]
    agg[, trajectory_id := do.call(paste, c(.SD, sep = "_")), .SDcols = c("Sex", AGE_LEVELS)]
    
    agg_long <- melt(agg,
                     id.vars = c("N", "Sex", "trajectory_id"),
                     measure.vars = AGE_LEVELS,
                     variable.name = "Age",
                     value.name = "State")
    agg_long[, Age := factor(Age, levels = AGE_LEVELS)]
    agg_long[, Sex := factor(Sex, levels = c("male", "female"), labels = c("Male", "Female"))]
    agg_long <- agg_long[!is.na(State)]
    
    p <- ggplot(agg_long,
                aes(x = Age, stratum = State, alluvium = trajectory_id,
                    y = N, fill = State)) +
        stat_alluvium(geom = "flow", alpha = 0.7, width = 0.4,
                      curve_type = "linear", na.rm = TRUE) +
        stat_stratum(width = 0.4, alpha = 0.9, color = "black",
                     linewidth = 0.3, na.rm = TRUE) +
        geom_text(stat = "stratum", aes(label = State),
                  size = 3, fontface = "bold", na.rm = TRUE) +
        facet_wrap(~ Sex, ncol = 2) +
        scale_fill_manual(values = STATE_COLS, na.translate = FALSE) +
        labs(title = title_label,
             subtitle = "All cell types merged | Male vs Female comparison",
             x = "Age", y = "Number of 80kb Bins", fill = "State") +
        theme_bw(base_size = 12) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
              panel.grid.minor = element_blank(),
              strip.text = element_text(size = 14, face = "bold"),
              strip.background = element_rect(fill = "gray90"),
              plot.title = element_text(face = "bold", size = 16))
    
    return(p)
}


#' Create side-by-side male/female comparison
#' @param data Merged data
#' @param stability_class Stability class
#' @param plot_type "faceted" or "colored"
#' @return Combined ggplot
create_sex_comparison <- function(data, stability_class, plot_type = "faceted") {
    
    if (plot_type == "faceted") {
        p_male <- create_stability_alluvial_by_celltype(data, stability_class, "male")
        p_female <- create_stability_alluvial_by_celltype(data, stability_class, "female")
    } else {
        p_male <- create_celltype_colored_alluvial(data, stability_class, "male")
        p_female <- create_celltype_colored_alluvial(data, stability_class, "female")
    }
    
    if (is.null(p_male) && is.null(p_female)) return(NULL)
    if (is.null(p_male)) return(p_female + labs(title = paste(stability_class, "- Female Only")))
    if (is.null(p_female)) return(p_male + labs(title = paste(stability_class, "- Male Only")))
    
    combined <- cowplot::plot_grid(
        p_male + theme(legend.position = "none"),
        p_female + theme(legend.position = "none"),
        ncol = 2, labels = c("Male", "Female"), label_size = 14
    )
    
    legend <- cowplot::get_legend(p_male + theme(legend.position = "bottom"))
    
    final <- cowplot::plot_grid(combined, legend, ncol = 1, rel_heights = c(1, 0.1))
    
    return(final)
}


## =========================================================
## SECTION 3: RUN VALIDATION
## =========================================================

banner("STEP 1: Validating input data")

validate_bin_data(state_with_stability, stability_ext_dt)

# Filter valid data
plot_data <- state_with_stability[!Stability %in% c("Missing", "Insufficient_Data")]

message(sprintf("  Plotting %d rows with valid stability", nrow(plot_data)))


## =========================================================
## SECTION 4: PER-STABILITY CLASS PLOTS (FACETED BY CELLTYPE)
## =========================================================

banner("STEP 2: Creating faceted plots by celltype")

for (stab in STABILITY_LEVELS) {
    stab_clean <- clean_name(stab)
    
    for (sx in c("male", "female")) {
        p <- create_stability_alluvial_by_celltype(plot_data, stab, sx)
        
        if (!is.null(p)) {
            fname <- sprintf("Alluvial_%s_%s_ByCelltype.pdf", stab_clean, tools::toTitleCase(sx))
            ggsave(file.path(CONFIG$outdir, fname), p, width = 14, height = 10)
            message(sprintf("  [SAVED] %s", fname))
        }
    }
}


## =========================================================
## SECTION 5: CELLTYPE-COLORED PLOTS (MERGED)
## =========================================================

banner("STEP 3: Creating celltype-colored merged plots")

for (stab in STABILITY_LEVELS) {
    stab_clean <- clean_name(stab)
    
    for (sx in c("male", "female")) {
        p <- create_celltype_colored_alluvial(plot_data, stab, sx)
        
        if (!is.null(p)) {
            fname <- sprintf("Alluvial_%s_%s_CelltypeColored.pdf", stab_clean, tools::toTitleCase(sx))
            ggsave(file.path(CONFIG$outdir, fname), p, width = 12, height = 8)
            message(sprintf("  [SAVED] %s", fname))
        }
    }
}


## =========================================================
## SECTION 6: MALE/FEMALE SIDE-BY-SIDE COMPARISONS
## =========================================================

banner("STEP 4: Creating male/female comparison plots")

for (stab in STABILITY_LEVELS) {
    stab_clean <- clean_name(stab)
    
    # Faceted version
    p_faceted <- create_sex_comparison(plot_data, stab, "faceted")
    if (!is.null(p_faceted)) {
        fname <- sprintf("Alluvial_%s_MaleFemale_Faceted.pdf", stab_clean)
        ggsave(file.path(CONFIG$outdir, fname), p_faceted, width = 24, height = 12)
        message(sprintf("  [SAVED] %s", fname))
    }
    
    # Celltype-colored version
    p_colored <- create_sex_comparison(plot_data, stab, "colored")
    if (!is.null(p_colored)) {
        fname <- sprintf("Alluvial_%s_MaleFemale_CelltypeColored.pdf", stab_clean)
        ggsave(file.path(CONFIG$outdir, fname), p_colored, width = 20, height = 8)
        message(sprintf("  [SAVED] %s", fname))
    }
}


## =========================================================
## SECTION 7: MERGED ALLUVIAL PLOTS
## =========================================================

banner("STEP 5: Creating merged alluvial plots")

# Per stability class
for (stab in STABILITY_LEVELS) {
    stab_clean <- clean_name(stab)
    title_label <- gsub("_", " ", stab)
    
    p <- create_merged_alluvial(plot_data, stab, title_label)
    if (!is.null(p)) {
        fname <- sprintf("Alluvial_Merged_%s.pdf", stab_clean)
        ggsave(file.path(CONFIG$outdir, fname), p, width = 16, height = 5)
        message(sprintf("  [SAVED] %s", fname))
    }
}


## =========================================================
## SECTION 8: ALL CLASSES COMBINED
## =========================================================

banner("STEP 6: Creating combined all-classes plots")

for (sx in c("male", "female")) {
    
    dat <- plot_data[Sex == sx]
    
    dat_wide <- dcast(dat, bin_id + Celltype + Stability ~ Age, value.var = "State",
                      fun.aggregate = function(x) x[1])
    dat_wide <- dat_wide[complete.cases(dat_wide)]
    
    if (nrow(dat_wide) == 0) next
    
    agg <- dat_wide[, .N, by = c("Celltype", "Stability", AGE_LEVELS)]
    agg[, trajectory_id := do.call(paste, c(.SD, sep = "_")),
        .SDcols = c("Celltype", "Stability", AGE_LEVELS)]
    
    agg_long <- melt(agg,
                     id.vars = c("N", "Celltype", "Stability", "trajectory_id"),
                     measure.vars = AGE_LEVELS,
                     variable.name = "Age",
                     value.name = "State")
    agg_long[, Age := factor(Age, levels = AGE_LEVELS)]
    agg_long <- agg_long[!is.na(State)]
    
    ct_order <- names(CELLTYPE_COLS)[names(CELLTYPE_COLS) %in% unique(agg_long$Celltype)]
    agg_long[, Celltype := factor(Celltype, levels = ct_order)]
    
    p <- ggplot(agg_long,
                aes(x = Age, stratum = State, alluvium = trajectory_id,
                    y = N, fill = Stability)) +
        stat_alluvium(geom = "flow", alpha = 0.7, width = 0.4,
                      curve_type = "linear", na.rm = TRUE) +
        stat_stratum(width = 0.4, alpha = 0.9, color = "black",
                     linewidth = 0.3, na.rm = TRUE) +
        geom_text(stat = "stratum", aes(label = State),
                  size = 2, fontface = "bold", na.rm = TRUE) +
        facet_wrap(~ Celltype, ncol = 5, scales = "free_y") +
        scale_fill_manual(values = STABILITY_COLS, na.translate = FALSE) +
        labs(
            title = sprintf("All Stability Classes - %s (5-Class Model)", tools::toTitleCase(sx)),
            subtitle = sprintf("Total: %s bins", format(sum(agg$N), big.mark = ",")),
            x = "Age", y = "Number of 80kb Bins", fill = "Stability Class"
        ) +
        theme_bw(base_size = 10) +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 8),
            panel.grid.minor = element_blank(),
            strip.text = element_text(size = 9, face = "bold"),
            strip.background = element_rect(fill = "gray95"),
            plot.title = element_text(face = "bold", size = 14),
            legend.position = "bottom"
        )
    
    fname <- sprintf("Alluvial_AllClasses_%s_ByCelltype.pdf", tools::toTitleCase(sx))
    ggsave(file.path(CONFIG$outdir, fname), p, width = 18, height = 14)
    message(sprintf("  [SAVED] %s", fname))
}


## =========================================================
## SECTION 9: SWITCHING BINS ONLY
## =========================================================

banner("STEP 7: Creating switching-only plots")

switching_data <- plot_data[Stability %in% SWITCHING_CLASSES]

message("[CHECK] Switching bins by sex and stability:")
print(switching_data[, .(n_bins = length(unique(bin_id))), by = .(Sex, Stability)][order(Sex, Stability)])

for (sx in c("male", "female")) {
    
    dat <- switching_data[Sex == sx]
    
    dat_wide <- dcast(dat, bin_id + Celltype + Stability ~ Age,
                      value.var = "State",
                      fun.aggregate = function(x) x[1])
    dat_wide <- dat_wide[complete.cases(dat_wide)]
    
    n_bins <- nrow(dat_wide)
    message(sprintf("\n[%s] %d switching bins", toupper(sx), n_bins))
    
    if (n_bins == 0) next
    
    agg <- dat_wide[, .N, by = c("Celltype", "Stability", AGE_LEVELS)]
    agg[, trajectory_id := .I]
    
    agg_long <- melt(agg,
                     id.vars = c("N", "Celltype", "Stability", "trajectory_id"),
                     measure.vars = AGE_LEVELS,
                     variable.name = "Age",
                     value.name = "State")
    agg_long[, Age := factor(Age, levels = AGE_LEVELS)]
    agg_long <- agg_long[!is.na(State)]
    
    # Colored by CELLTYPE
    p_ct <- ggplot(agg_long,
                   aes(x = Age, stratum = State, alluvium = trajectory_id,
                       y = N, fill = Celltype)) +
        stat_alluvium(geom = "flow", alpha = 0.7, width = 0.4,
                      curve_type = "linear", na.rm = TRUE) +
        stat_stratum(width = 0.4, alpha = 0.9, color = "black",
                     linewidth = 0.3, na.rm = TRUE) +
        geom_text(stat = "stratum", aes(label = State),
                  size = 5, fontface = "bold", na.rm = TRUE) +
        scale_fill_manual(values = CELLTYPE_COLS, na.translate = FALSE) +
        scale_y_continuous(labels = scales::comma, expand = expansion(mult = c(0, 0.02))) +
        labs(
            title = sprintf("Switching Compartments - %s", tools::toTitleCase(sx)),
            subtitle = sprintf("Monotonic + Non-Monotonic | Total: %s bins", format(n_bins, big.mark = ",")),
            x = "Age", y = "Number of 80kb Bins", fill = "Cell Type"
        ) +
        theme_bw(base_size = 16) +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 14),
            panel.grid.minor = element_blank(),
            plot.title = element_text(face = "bold", size = 20),
            legend.position = "right"
        )
    
    fname <- sprintf("Alluvial_Switching_ByCelltype_%s.pdf", tools::toTitleCase(sx))
    ggsave(file.path(CONFIG$outdir, fname), p_ct, width = 14, height = 10)
    message(sprintf("  [SAVED] %s", fname))
    
    # Colored by SWITCH TYPE
    p_sw <- ggplot(agg_long,
                   aes(x = Age, stratum = State, alluvium = trajectory_id,
                       y = N, fill = Stability)) +
        stat_alluvium(geom = "flow", alpha = 0.7, width = 0.4,
                      curve_type = "linear", na.rm = TRUE) +
        stat_stratum(width = 0.4, alpha = 0.9, color = "black",
                     linewidth = 0.3, na.rm = TRUE) +
        geom_text(stat = "stratum", aes(label = State),
                  size = 5, fontface = "bold", na.rm = TRUE) +
        scale_fill_manual(values = SWITCHING_COLS,
                          labels = c("Monotonic_A_to_R" = "A → R (Closing)",
                                     "Monotonic_R_to_A" = "R → A (Opening)",
                                     "Non_Monotonic" = "Non-Monotonic"),
                          na.translate = FALSE) +
        scale_y_continuous(labels = scales::comma, expand = expansion(mult = c(0, 0.02))) +
        labs(
            title = sprintf("Switching Compartments - %s", tools::toTitleCase(sx)),
            subtitle = sprintf("Monotonic + Non-Monotonic | Total: %s bins", format(n_bins, big.mark = ",")),
            x = "Age", y = "Number of 80kb Bins", fill = "Switch Type"
        ) +
        theme_bw(base_size = 16) +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 14),
            panel.grid.minor = element_blank(),
            plot.title = element_text(face = "bold", size = 20),
            legend.position = "right"
        )
    
    fname <- sprintf("Alluvial_Switching_BySwitchType_%s.pdf", tools::toTitleCase(sx))
    ggsave(file.path(CONFIG$outdir, fname), p_sw, width = 14, height = 10)
    message(sprintf("  [SAVED] %s", fname))
}


## =========================================================
## SECTION 10: SWITCHING SIDE-BY-SIDE
## =========================================================

banner("STEP 8: Creating switching side-by-side plots")

switching_both <- switching_data[Sex %in% c("male", "female")]

dat_wide_sw <- dcast(switching_both, bin_id + Sex + Celltype + Stability ~ Age,
                     value.var = "State", fun.aggregate = function(x) x[1])
dat_wide_sw <- dat_wide_sw[complete.cases(dat_wide_sw)]

n_male <- sum(dat_wide_sw$Sex == "male")
n_female <- sum(dat_wide_sw$Sex == "female")
message(sprintf("\n[DATA] Male: %d bins | Female: %d bins", n_male, n_female))

if (nrow(dat_wide_sw) > 0) {
    
    agg_sw <- dat_wide_sw[, .N, by = c("Sex", "Celltype", "Stability", AGE_LEVELS)]
    agg_sw[, trajectory_id := .I]
    
    agg_long_sw <- melt(agg_sw,
                        id.vars = c("N", "Sex", "Celltype", "Stability", "trajectory_id"),
                        measure.vars = AGE_LEVELS,
                        variable.name = "Age",
                        value.name = "State")
    agg_long_sw[, Age := factor(Age, levels = AGE_LEVELS)]
    agg_long_sw[, Sex := factor(Sex, levels = c("male", "female"), labels = c("Male", "Female"))]
    agg_long_sw <- agg_long_sw[!is.na(State)]
    
    # By Celltype
    p_ct <- ggplot(agg_long_sw,
                   aes(x = Age, stratum = State, alluvium = trajectory_id,
                       y = N, fill = Celltype)) +
        stat_alluvium(geom = "flow", alpha = 0.7, width = 0.4,
                      curve_type = "linear", na.rm = TRUE) +
        stat_stratum(width = 0.4, alpha = 0.9, color = "black",
                     linewidth = 0.3, na.rm = TRUE) +
        geom_text(stat = "stratum", aes(label = State),
                  size = 4, fontface = "bold", na.rm = TRUE) +
        facet_wrap(~ Sex, ncol = 2) +
        scale_fill_manual(values = CELLTYPE_COLS, na.translate = FALSE) +
        scale_y_continuous(labels = scales::comma, expand = expansion(mult = c(0, 0.02))) +
        labs(
            title = "Switching Compartments - Male vs Female",
            subtitle = sprintf("Monotonic + Non-Monotonic | Male: %s | Female: %s bins",
                               format(n_male, big.mark = ","), format(n_female, big.mark = ",")),
            x = "Age", y = "Number of 80kb Bins", fill = "Cell Type"
        ) +
        theme_bw(base_size = 14) +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 12),
            strip.text = element_text(size = 18, face = "bold"),
            strip.background = element_rect(fill = "gray95"),
            panel.grid.minor = element_blank(),
            plot.title = element_text(face = "bold", size = 20),
            legend.position = "right"
        )
    
    ggsave(file.path(CONFIG$outdir, "Alluvial_Switching_SideBySide_ByCelltype.pdf"),
           p_ct, width = 23, height = 8)
    message("  [SAVED] Alluvial_Switching_SideBySide_ByCelltype.pdf")
    
    # By Switch Type
    p_sw <- ggplot(agg_long_sw,
                   aes(x = Age, stratum = State, alluvium = trajectory_id,
                       y = N, fill = Stability)) +
        stat_alluvium(geom = "flow", alpha = 0.7, width = 0.4,
                      curve_type = "linear", na.rm = TRUE) +
        stat_stratum(width = 0.4, alpha = 0.9, color = "black",
                     linewidth = 0.3, na.rm = TRUE) +
        geom_text(stat = "stratum", aes(label = State),
                  size = 4, fontface = "bold", na.rm = TRUE) +
        facet_wrap(~ Sex, ncol = 2) +
        scale_fill_manual(values = SWITCHING_COLS,
                          labels = c("Monotonic_A_to_R" = "A → R (Closing)",
                                     "Monotonic_R_to_A" = "R → A (Opening)",
                                     "Non_Monotonic" = "Non-Monotonic"),
                          na.translate = FALSE) +
        scale_y_continuous(labels = scales::comma, expand = expansion(mult = c(0, 0.02))) +
        labs(
            title = "Switching Compartments - Male vs Female",
            subtitle = sprintf("Monotonic + Non-Monotonic | Male: %s | Female: %s bins",
                               format(n_male, big.mark = ","), format(n_female, big.mark = ",")),
            x = "Age", y = "Number of 80kb Bins", fill = "Switch Type"
        ) +
        theme_bw(base_size = 14) +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 12),
            strip.text = element_text(size = 18, face = "bold"),
            strip.background = element_rect(fill = "gray95"),
            panel.grid.minor = element_blank(),
            plot.title = element_text(face = "bold", size = 20),
            legend.position = "right"
        )
    
    ggsave(file.path(CONFIG$outdir, "Alluvial_Switching_SideBySide_BySwitchType.pdf"),
           p_sw, width = 23, height = 8)
    message("  [SAVED] Alluvial_Switching_SideBySide_BySwitchType.pdf")
}


## =========================================================
## SECTION 11: COMBINED ALL (NO FACETING)
## =========================================================

banner("STEP 9: Creating combined all plots (no faceting)")

for (sx in c("male", "female")) {
    
    dat <- plot_data[Sex == sx]
    
    dat_wide <- dcast(dat, bin_id + Celltype + Stability ~ Age, 
                      value.var = "State", fun.aggregate = function(x) x[1])
    dat_wide <- dat_wide[complete.cases(dat_wide)]
    
    n_bins <- nrow(dat_wide)
    message(sprintf("  %s: %d complete bins", sx, n_bins))
    
    if (n_bins == 0) next
    
    agg <- dat_wide[, .N, by = c("Celltype", "Stability", AGE_LEVELS)]
    agg[, trajectory_id := .I]
    
    agg_long <- melt(agg,
                     id.vars = c("N", "Celltype", "Stability", "trajectory_id"),
                     measure.vars = AGE_LEVELS,
                     variable.name = "Age",
                     value.name = "State")
    agg_long[, Age := factor(Age, levels = AGE_LEVELS)]
    agg_long <- agg_long[!is.na(State)]
    
    # By Celltype
    p_ct <- ggplot(agg_long,
                   aes(x = Age, stratum = State, alluvium = trajectory_id,
                       y = N, fill = Celltype)) +
        stat_alluvium(geom = "flow", alpha = 0.6, width = 0.4,
                      curve_type = "linear", na.rm = TRUE) +
        stat_stratum(width = 0.4, alpha = 0.9, color = "black",
                     linewidth = 0.2, na.rm = TRUE) +
        geom_text(stat = "stratum", aes(label = State),
                  size = 4, fontface = "bold", na.rm = TRUE) +
        scale_fill_manual(values = CELLTYPE_COLS, na.translate = FALSE) +
        labs(
            title = sprintf("Compartment Dynamics - %s", tools::toTitleCase(sx)),
            subtitle = sprintf("All stability classes + all celltypes | %s bins", 
                               format(n_bins, big.mark = ",")),
            x = "Age", y = "Number of 80kb Bins", fill = "Cell Type"
        ) +
        theme_bw(base_size = 14) +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 12),
            panel.grid.minor = element_blank(),
            plot.title = element_text(face = "bold", size = 18),
            legend.position = "right"
        )
    
    fname <- sprintf("Alluvial_Combined_AllStability_AllCelltypes_%s_ByCelltype.pdf", 
                     tools::toTitleCase(sx))
    ggsave(file.path(CONFIG$outdir, fname), p_ct, width = 14, height = 10)
    message(sprintf("  [SAVED] %s", fname))
    
    # By Stability
    p_stab <- ggplot(agg_long,
                     aes(x = Age, stratum = State, alluvium = trajectory_id,
                         y = N, fill = Stability)) +
        stat_alluvium(geom = "flow", alpha = 0.6, width = 0.4,
                      curve_type = "linear", na.rm = TRUE) +
        stat_stratum(width = 0.4, alpha = 0.9, color = "black",
                     linewidth = 0.2, na.rm = TRUE) +
        geom_text(stat = "stratum", aes(label = State),
                  size = 4, fontface = "bold", na.rm = TRUE) +
        scale_fill_manual(values = STABILITY_COLS, na.translate = FALSE) +
        labs(
            title = sprintf("Compartment Dynamics - %s", tools::toTitleCase(sx)),
            subtitle = sprintf("All stability classes + all celltypes | %s bins",
                               format(n_bins, big.mark = ",")),
            x = "Age", y = "Number of 80kb Bins", fill = "Stability Class"
        ) +
        theme_bw(base_size = 14) +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 12),
            panel.grid.minor = element_blank(),
            plot.title = element_text(face = "bold", size = 18),
            legend.position = "right"
        )
    
    fname <- sprintf("Alluvial_Combined_AllStability_AllCelltypes_%s_ByStability.pdf",
                     tools::toTitleCase(sx))
    ggsave(file.path(CONFIG$outdir, fname), p_stab, width = 14, height = 10)
    message(sprintf("  [SAVED] %s", fname))
}


## =========================================================
## FINAL SUMMARY
## =========================================================

cat("\n")
cat("╔══════════════════════════════════════════════════════════════════╗\n")
cat("║                  ALLUVIAL PLOTS COMPLETE                         ║\n")
cat("╠══════════════════════════════════════════════════════════════════╣\n")
cat("║  Per stability class (faceted by celltype):                      ║\n")
cat("║    • Alluvial_*_Male/Female_ByCelltype.pdf                       ║\n")
cat("║                                                                  ║\n")
cat("║  Celltype-colored (merged):                                      ║\n")
cat("║    • Alluvial_*_Male/Female_CelltypeColored.pdf                  ║\n")
cat("║                                                                  ║\n")
cat("║  Male/Female side-by-side:                                       ║\n")
cat("║    • Alluvial_*_MaleFemale_Faceted.pdf                           ║\n")
cat("║    • Alluvial_*_MaleFemale_CelltypeColored.pdf                   ║\n")
cat("║                                                                  ║\n")
cat("║  Merged alluvials:                                               ║\n")
cat("║    • Alluvial_Merged_*.pdf                                       ║\n")
cat("║                                                                  ║\n")
cat("║  All classes combined:                                           ║\n")
cat("║    • Alluvial_AllClasses_Male/Female_ByCelltype.pdf              ║\n")
cat("║                                                                  ║\n")
cat("║  Switching only:                                                 ║\n")
cat("║    • Alluvial_Switching_ByCelltype/BySwitchType_Male/Female.pdf  ║\n")
cat("║    • Alluvial_Switching_SideBySide_*.pdf                         ║\n")
cat("║                                                                  ║\n")
cat("║  Combined all (no faceting):                                     ║\n")
cat("║    • Alluvial_Combined_AllStability_AllCelltypes_*.pdf           ║\n")
cat("╚══════════════════════════════════════════════════════════════════╝\n")
