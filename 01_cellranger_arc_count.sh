#!/bin/bash
#===============================================================================
# Step 01: CellRanger ARC alignment and feature counting
#          (snRNA-seq + snATAC-seq, 10x Multiome)
#===============================================================================
# Description: Aligns and counts paired snRNA-seq + snATAC-seq multiome
#              libraries with CellRanger ARC for the mouse liver aging study.
#
# Samples:     40 total across 5 age groups, n=8 per group
#                Young (Y)            <age> months
#                Mid-age (MA)         <age> months
#                Old (O)              <age> months
#                Pre-geriatric (PG)   <age> months
#                Geriatric (G)        <age> months
#
# Software:    CellRanger ARC (version set in CELLRANGER_ARC_VERSION below)
# Reference:   mm10 (refdata-cellranger-arc-mm10-2020-A)
#
# Usage:
#   All 40 samples sequentially in one job:
#       sbatch 01_cellranger_arc_count.sh
#   One sample per array task (runs samples in parallel):
#       sbatch --array=1-40 01_cellranger_arc_count.sh
#
# Outputs:     <sample_id>/outs/                    per-sample CellRanger ARC output
#              libraries/libraries_<sample_id>.csv  FASTQ pairing used for each sample
#===============================================================================
#SBATCH --job-name=cellranger_arc
#SBATCH --cpus-per-task=50
#SBATCH --time=168:00:00
#SBATCH --mem=600g
## If your cluster requires a partition for high-memory jobs, remove one "#"
## from the line below and set the partition name:
##SBATCH --partition=<high_memory_partition>

set -euo pipefail
ulimit -u 4096   # CellRanger ARC launches many processes; raise the per-user limit

#-------------------------------------------------------------------------------
# Software
#-------------------------------------------------------------------------------
CELLRANGER_ARC_VERSION="X.Y.Z"   # set to the exact version used for the paper
# Load via environment modules if available; otherwise cellranger-arc must be on PATH
if command -v module >/dev/null 2>&1; then
    module load "cellranger-arc/${CELLRANGER_ARC_VERSION}"
fi

#-------------------------------------------------------------------------------
# Configuration: update paths for your environment
#-------------------------------------------------------------------------------
# REF_DIR: folder containing refdata-cellranger-arc-mm10-2020-A
#          (download from 10x Genomics), e.g. export REF_DIR=/path/to/references
REFERENCE="${REF_DIR:?REF_DIR is not set; export REF_DIR=/path/to/references}/refdata-cellranger-arc-mm10-2020-A"
LOCALCORES="${SLURM_CPUS_PER_TASK:-50}"
LOCALMEM=560   # GB; kept below the 600 GB SLURM allocation to leave headroom

# FASTQ directories (update these paths)
FASTQ_GEX="path/to/gene_expression/fastqs"
FASTQ_ATAC="path/to/chromatin_accessibility/fastqs"

# Library CSVs are kept as a record of which FASTQs were paired per sample
LIB_DIR="libraries"
mkdir -p "${LIB_DIR}"

#-------------------------------------------------------------------------------
# Sample sheet: 5 age groups x 8 animals = 40 samples
# Naming: output ID <group><n>, GEX sample snRNA_<group><n>,
#         ATAC sample snATAC_<group><n>
#-------------------------------------------------------------------------------
AGE_GROUPS=(Y MA O PG G)   # Young, Mid-age, Old, Pre-geriatric, Geriatric
N_PER_GROUP=8

SAMPLES=()
for group in "${AGE_GROUPS[@]}"; do
    for i in $(seq 1 "${N_PER_GROUP}"); do
        SAMPLES+=("${group}${i}")
    done
done

#-------------------------------------------------------------------------------
# Processing function
#-------------------------------------------------------------------------------
run_cellranger_arc() {
    local sample_id=$1
    local gex_sample="snRNA_${sample_id}"
    local atac_sample="snATAC_${sample_id}"
    local lib_csv="${LIB_DIR}/libraries_${sample_id}.csv"

    # Skip samples that already finished (safe to resubmit after a timeout)
    if [[ -f "${sample_id}/outs/summary.csv" ]]; then
        echo "Skipping ${sample_id}: output already complete"
        return 0
    fi

    echo "Processing: ${sample_id} - $(date)"

    cat > "${lib_csv}" <<EOF
fastqs,sample,library_type
${FASTQ_GEX},${gex_sample},Gene Expression
${FASTQ_ATAC},${atac_sample},Chromatin Accessibility
EOF

    cellranger-arc count \
        --id="${sample_id}" \
        --reference="${REFERENCE}" \
        --libraries="${lib_csv}" \
        --localcores="${LOCALCORES}" \
        --localmem="${LOCALMEM}"

    echo "Finished: ${sample_id} - $(date)"
}

#-------------------------------------------------------------------------------
# Run: one sample per array task, or all samples sequentially
#-------------------------------------------------------------------------------
if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    idx=$((SLURM_ARRAY_TASK_ID - 1))
    if (( idx < 0 || idx >= ${#SAMPLES[@]} )); then
        echo "Array task ${SLURM_ARRAY_TASK_ID} is out of range (1-${#SAMPLES[@]})" >&2
        exit 1
    fi
    run_cellranger_arc "${SAMPLES[$idx]}"
else
    for sample_id in "${SAMPLES[@]}"; do
        run_cellranger_arc "${sample_id}"
    done
fi

echo "Job complete: $(date)"
