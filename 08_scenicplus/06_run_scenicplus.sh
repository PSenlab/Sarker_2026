#!/bin/bash
#SBATCH --job-name=scenicplus
#SBATCH --output=scenicplus_%j.log
#SBATCH --error=scenicplus_%j.err
#SBATCH --time=148:00:00
#SBATCH --partition=largemem
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=50
#SBATCH --mem=1600G
# ==============================================================================
# Launch the SCENIC+ Snakemake Pipeline (Step 6)
# ==============================================================================
#
# Description:
#   Submits the full SCENIC+ Snakemake workflow to SLURM. Scale the
#   --time, --cpus-per-task, and --mem directives to your allocation.
#
# Prerequisites:
#   - Snakemake pipeline already initialized:
#         scenicplus init_snakemake --out_dir scplus_pipeline
#   - config.yaml written by 05_make_scenicplus_configs.py copied in as:
#         scplus_pipeline/Snakemake/config/config.yaml
#
# Submit:
#   sbatch 06_run_scenicplus.sh
# ==============================================================================

set -euo pipefail

ulimit -u 4096
ulimit -n 8192

# Activate conda environment containing the scenicplus package
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate scenicplus

# Run from inside the initialized Snakemake folder
cd scplus_pipeline/Snakemake/

python ../../06_run_scenicplus.py
