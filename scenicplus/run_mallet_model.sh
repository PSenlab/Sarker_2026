#!/bin/bash
# ==============================================================================
# Launch MALLET LDA Model Fitting (Step 2)
# ==============================================================================
#
# Description:
#   Sets up the MALLET / Java environment and runs step2_mallet_model.py.
#   Adjust memory, conda environment, and Java paths to your system.
#
# Example submission (SLURM):
#   sbatch --partition=largemem --cpus-per-task=30 --mem=560g --time=12:00:00 \
#          run_mallet_model.sh
# ==============================================================================

set -euo pipefail

# Increase process / file descriptor limits
ulimit -u 4096
ulimit -n 8192

# Activate conda environment containing pycisTopic + MALLET
# (adjust to your install)
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate scenicplus

# Java 8 required for MALLET; set JAVA_HOME to the absolute path of your JDK
export JAVA_HOME="path/to/java-1.8.0-openjdk"
export PATH="${JAVA_HOME}/bin:${PATH}"

# MALLET installation - point to the unpacked Mallet-202108 directory
export MALLET_HOME="path/to/Mallet-202108"
export PATH="${MALLET_HOME}/bin:${PATH}"

# MALLET / JVM memory - tune to your allocation
export MALLET_MEMORY=500G
export JAVA_TOOL_OPTIONS="-Xms500g -Xmx500g -XX:+UseG1GC"

echo "JAVA_HOME:         ${JAVA_HOME}"
echo "MALLET_HOME:       ${MALLET_HOME}"
echo "JAVA_TOOL_OPTIONS: ${JAVA_TOOL_OPTIONS}"

python step2_mallet_model.py
