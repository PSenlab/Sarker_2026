#!/usr/bin/env Rscript
#===============================================================================
# SoupX ambient RNA removal for single-nucleus multi-ome data
#===============================================================================
# Description: Removes ambient RNA contamination from snRNA-seq data using 
#              SoupX with automated contamination estimation
#
# Input:       Cell Ranger ARC output directories (outs folder with raw and 
#              filtered matrices)
# Output:      Corrected count matrices in 10x-compatible format
#
# Reference:   Young MD, Behjati S (2020). SoupX removes ambient RNA 
#              contamination from droplet-based single-cell RNA sequencing data.
#              GigaScience, 9(12):giaa151
#===============================================================================

library(SoupX)
library(DropletUtils)

#-------------------------------------------------------------------------------
# Configuration - Update paths for your environment
#-------------------------------------------------------------------------------
input_dir <- "path/to/cellranger_arc/outputs"
output_dir <- "path/to/soupx/corrected"

# Sample prefixes by age group
age_groups <- list(
  young       = paste0("Y", 1:8),
  mid_age     = paste0("MA", 1:8),
  old         = paste0("O", 1:8),
  pre_ger     = paste0("PG", 1:8),
  geriatric   = paste0("G", 1:8)
)

#-------------------------------------------------------------------------------
# SoupX Processing Function
#-------------------------------------------------------------------------------
run_soupx <- function(sample_id, input_dir, output_dir) {
  
  message(paste0("Processing: ", sample_id, " - ", Sys.time()))
  
  # Load 10X data (requires raw and filtered matrices)
  sample_path <- file.path(input_dir, sample_id, "outs")
  sc <- load10X(sample_path)
  
  # Estimate contamination fraction automatically
  sc <- autoEstCont(sc)
  
  # Adjust counts using subtraction method
  sc_corrected <- adjustCounts(sc, method = "subtraction", roundToInt = TRUE)
  
  # Export corrected counts
  output_path <- file.path(output_dir, sample_id)
  write10xCounts(output_path, sc_corrected, version = "3")
  
  message(paste0("Completed: ", sample_id))
}

#-------------------------------------------------------------------------------
# Process All Samples
#-------------------------------------------------------------------------------
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

for (group in names(age_groups)) {
  message(paste0("\n=== Processing ", group, " samples ===\n"))
  
  for (sample_id in age_groups[[group]]) {
    tryCatch({
      run_soupx(sample_id, input_dir, output_dir)
    }, error = function(e) {
      message(paste0("Error processing ", sample_id, ": ", e$message))
    })
  }
}

message(paste0("\nPipeline complete: ", Sys.time()))
