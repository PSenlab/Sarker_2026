#!/bin/bash
#===============================================================================
# CellRanger Arc processing pipeline for single-nucleus multi-ome samples
#===============================================================================
# Description: Processes snRNA-seq + snATAC-seq multiome libraries using 
#              CellRanger Arc for mouse liver aging study
#
# Samples:     40 total across 5 age groups (young, mid-age, old, pre-geriatric, 
#              geriatric); n=8 per group
#
# Reference:   mm10 (refdata-cellranger-arc-mm10-2020-A)
#===============================================================================
#SBATCH --job-name=cellranger_arc
#SBATCH --partition=largemem
#SBATCH --cpus-per-task=50
#SBATCH --time=168:00:00
#SBATCH --mem=600g
set -euo pipefail
ulimit -u 4096
module load cellranger-arc
#-------------------------------------------------------------------------------
# Configuration - Update paths for your environment
#-------------------------------------------------------------------------------
REFERENCE="${CELLRANGER_ARC_REF}/refdata-cellranger-arc-mm10-2020-A"
LOCALCORES=50
LOCALMEM=600
# FASTQ directories (update these paths)
FASTQ_GEX="path/to/gene_expression/fastqs"
FASTQ_ATAC="path/to/chromatin_accessibility/fastqs"
#-------------------------------------------------------------------------------
# Processing Function
#-------------------------------------------------------------------------------
run_cellranger_arc() {
    local sample_id=$1
    local gex_sample=$2
    local atac_sample=$3
    
    echo "Processing: ${sample_id} - $(date)"
    
    cat > libraries_${sample_id}.csv <<EOF
fastqs,sample,library_type
${FASTQ_GEX},${gex_sample},Gene Expression
${FASTQ_ATAC},${atac_sample},Chromatin Accessibility
EOF
    cellranger-arc count \
        --id="${sample_id}" \
        --reference="${REFERENCE}" \
        --libraries="libraries_${sample_id}.csv" \
        --localcores="${LOCALCORES}" \
        --localmem="${LOCALMEM}"
    
    rm -f "libraries_${sample_id}.csv"
}
#-------------------------------------------------------------------------------
# Sample Processing
# Format: run_cellranger_arc <output_id> <gex_sample_name> <atac_sample_name>
#-------------------------------------------------------------------------------
# Young (n=8)
for i in {1..8}; do
    run_cellranger_arc "Y${i}" "snRNA_Y${i}" "snATAC_Y${i}"
done
# Mid-age (n=8)
for i in {1..8}; do
    run_cellranger_arc "MA${i}" "snRNA_MA${i}" "snATAC_MA${i}"
done
# Old (n=8)
for i in {1..8}; do
    run_cellranger_arc "O${i}" "snRNA_O${i}" "snATAC_O${i}"
done
# Pre-geriatric (n=8)
for i in {1..8}; do
    run_cellranger_arc "PG${i}" "snRNA_PG${i}" "snATAC_PG${i}"
done
# Geriatric (n=8)
for i in {1..8}; do
    run_cellranger_arc "G${i}" "snRNA_G${i}" "snATAC_G${i}"
done
echo "Pipeline complete: $(date)"
