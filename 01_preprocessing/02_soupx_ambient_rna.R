#!/usr/bin/env Rscript
#===============================================================================
# Step 02: SoupX ambient RNA removal (snRNA-seq modality of 10x Multiome)
#===============================================================================
# Description: Estimates and removes ambient RNA contamination from the
#              snRNA-seq counts of each multiome sample using SoupX with
#              automated contamination estimation (autoEstCont).
#
# Input:       Step 01 CellRanger ARC output, one folder per sample:
#                <input_dir>/<sample_id>/outs/
#                  raw_feature_bc_matrix/       all barcodes (ambient profile)
#                  filtered_feature_bc_matrix/  called cells
#                  analysis/clustering/gex/graphclust/clusters.csv
#              SoupX::load10X reads these directly and keeps only
#              Gene Expression features (ATAC peaks are not used here).
#
# Output:      <output_dir>/<sample_id>/                  corrected counts (10x v3)
#              <output_dir>/qc/<sample_id>_autoEstCont.pdf contamination fit plot
#              <output_dir>/soupx_contamination_summary.csv
#                  estimated contamination fraction (rho) and cell count per sample
#              <output_dir>/sessionInfo.txt               R and package versions
#
# Usage:       Rscript 02_soupx_ambient_rna.R
#
# Reference:   Young MD, Behjati S (2020). SoupX removes ambient RNA
#              contamination from droplet-based single-cell RNA sequencing data.
#              GigaScience, 9(12):giaa151
#===============================================================================

suppressPackageStartupMessages({
  library(SoupX)
  library(DropletUtils)
})

#-------------------------------------------------------------------------------
# Configuration: update paths for your environment
#-------------------------------------------------------------------------------
input_dir  <- "path/to/cellranger_arc/outputs"   # Step 01 output folder
output_dir <- "path/to/soupx/corrected"

#-------------------------------------------------------------------------------
# Sample sheet: same naming as Step 01 (5 age groups x 8 animals = 40 samples)
#-------------------------------------------------------------------------------
age_groups <- c(Y  = "Young",
                MA = "Mid-age",
                O  = "Old",
                PG = "Pre-geriatric",
                G  = "Geriatric")
n_per_group <- 8

samples <- data.frame(
  sample_id = paste0(rep(names(age_groups), each = n_per_group),
                     rep(seq_len(n_per_group), times = length(age_groups))),
  age_group = rep(unname(age_groups), each = n_per_group),
  stringsAsFactors = FALSE
)

#-------------------------------------------------------------------------------
# SoupX processing function
#-------------------------------------------------------------------------------
run_soupx <- function(sample_id, input_dir, output_dir, qc_dir) {

  message("Processing: ", sample_id, " - ", Sys.time())

  sample_path <- file.path(input_dir, sample_id, "outs")
  if (!dir.exists(sample_path)) {
    stop("CellRanger ARC output not found: ", sample_path)
  }

  # Load raw and filtered matrices plus the ARC GEX graph-based clusters
  sc <- load10X(sample_path)

  # autoEstCont needs clusters; stop clearly if load10X could not find them
  if (!"clusters" %in% colnames(sc$metaData)) {
    stop("No clusters loaded for ", sample_id,
         "; expected analysis/clustering/gex/graphclust/clusters.csv")
  }

  # Estimate contamination fraction automatically and save the diagnostic plot
  pdf(file.path(qc_dir, paste0(sample_id, "_autoEstCont.pdf")))
  sc <- tryCatch(autoEstCont(sc), finally = dev.off())
  rho <- sc$metaData$rho[1]

  # Adjust counts using the subtraction method, rounded to integers
  sc_corrected <- adjustCounts(sc, method = "subtraction", roundToInt = TRUE)

  # Export corrected counts (overwrite = TRUE so the step can be rerun)
  write10xCounts(file.path(output_dir, sample_id), sc_corrected,
                 version = "3", overwrite = TRUE)

  message("Completed: ", sample_id, " (rho = ", signif(rho, 3), ")")

  data.frame(sample_id = sample_id,
             rho       = rho,
             n_cells   = ncol(sc_corrected),
             stringsAsFactors = FALSE)
}

#-------------------------------------------------------------------------------
# Process all samples
#-------------------------------------------------------------------------------
qc_dir <- file.path(output_dir, "qc")
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)

results <- list()
failed  <- character(0)

for (i in seq_len(nrow(samples))) {
  sample_id <- samples$sample_id[i]

  # Announce each age group once
  if (i == 1 || samples$age_group[i] != samples$age_group[i - 1]) {
    message("\n=== Processing ", samples$age_group[i], " samples ===\n")
  }

  res <- tryCatch(
    run_soupx(sample_id, input_dir, output_dir, qc_dir),
    error = function(e) {
      message("Error processing ", sample_id, ": ", conditionMessage(e))
      NULL
    }
  )

  if (is.null(res)) {
    failed <- c(failed, sample_id)
  } else {
    results[[sample_id]] <- res
  }
}

#-------------------------------------------------------------------------------
# Summary, versions and exit status
#-------------------------------------------------------------------------------
if (length(results) > 0) {
  summary_df <- do.call(rbind, results)
  summary_df$age_group <- samples$age_group[match(summary_df$sample_id,
                                                  samples$sample_id)]
  summary_df <- summary_df[, c("sample_id", "age_group", "rho", "n_cells")]
  write.csv(summary_df,
            file.path(output_dir, "soupx_contamination_summary.csv"),
            row.names = FALSE)
}

writeLines(capture.output(sessionInfo()),
           file.path(output_dir, "sessionInfo.txt"))

# Fail loudly so missing samples are never passed silently to the next step
if (length(failed) > 0) {
  stop(length(failed), " of ", nrow(samples), " samples failed: ",
       paste(failed, collapse = ", "))
}

message("\nPipeline complete: ", Sys.time())
