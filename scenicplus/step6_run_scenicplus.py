#!/usr/bin/env python3
# ==============================================================================
# SCENIC+ Step 5: Launch the Snakemake Pipeline
# ==============================================================================
#
# Description:
#   Launches the SCENIC+ Snakemake workflow. Assumes the Snakemake folder
#   has been initialized via:
#       scenicplus init_snakemake --out_dir scplus_pipeline
#   and that step5_config.yaml has been copied in as:
#       scplus_pipeline/Snakemake/config/config.yaml
#
# Usage:
#   cd scplus_pipeline/Snakemake
#   python /path/to/step5_run_scenicplus.py
#
#
# ==============================================================================

import os

# Number of cores to allocate to the Snakemake DAG
N_CORES = 50

# Unlock any stale lock from a prior interrupted run
os.system("snakemake --unlock")

# Run the pipeline
os.system(f"snakemake --cores {N_CORES}")
