#!/bin/bash
#===============================================================================
# Run GWAS-SCAVENGE for All Liver Cell Types
#
# Launches run_gwas_scavenge.R for each cell type sequentially.
# Modify THREADS/CORES based on available compute resources.
#
# Usage:
#   bash run_all_celltypes.sh
#   nohup bash run_all_celltypes.sh > run_all.log 2>&1 &
#===============================================================================

set -euo pipefail

#--- PATHS (edit these) --------------------------------------------------------
BASE_DIR="/data/sarkern2/multiome_liver/Seurat/epigenome/gwas_peak_gene"
ARCHR_BASE="/data/sarkern2/multiome_liver/Seurat/archR/ArchR_Projects"
GWAS_CATALOG="${BASE_DIR}/gwas_catalog_v1.0.2-associations.tsv"
CHAIN_FILE="${BASE_DIR}/hg38ToMm10.over.chain"
SCRIPT="run_gwas_scavenge.R"

#--- RESOURCES -----------------------------------------------------------------
THREADS=60
SCAVENGE_CORES=55

#--- CELL TYPES ----------------------------------------------------------------
# Format: "CellTypeName|directory_name|ArchR_project_folder"
# Stellate excluded: only 169 P2G peaks → insufficient for SCAVENGE

CELL_TYPES=(
  "Hepatocyte|hepatocyte|Step8_Hepatocyte_CCAN_P2G"
  "Endothelial_01|endothelial_01|Step8_Endothelial_01_CCAN_P2G"
  "Endothelial_02|endothelial_02|Step8_Endothelial_02_CCAN_P2G"
  "Kupffer|kupffer|Step8_Kupffer_CCAN_P2G"
  "MoMFs|momfs|Step8_MoMFs_CCAN_P2G"
  "Cholangiocyte_01|cholangiocyte_01|Step8_Cholangiocyte_01_CCAN_P2G"
  "Cholangiocyte_02|cholangiocyte_02|Step8_Cholangiocyte_02_CCAN_P2G"
  "Lymp_B|lymp_b|Step8_Lymp_B_CCAN_P2G"
  "Lymp_T|lymp_t|Step8_Lymp_T_CCAN_P2G"
)

#--- RUN -----------------------------------------------------------------------
echo "================================================================"
echo "GWAS-SCAVENGE: Running ${#CELL_TYPES[@]} cell types"
echo "================================================================"
echo ""

for entry in "${CELL_TYPES[@]}"; do
  IFS='|' read -r CT_NAME DIR_NAME ARCHR_PROJ <<< "$entry"
  
  OUTPUT_DIR="${BASE_DIR}/${DIR_NAME}"
  ARCHR_PATH="${ARCHR_BASE}/${ARCHR_PROJ}"
  
  echo "--------------------------------------------------------------"
  echo "[$(date)] Starting: ${CT_NAME}"
  echo "  ArchR:  ${ARCHR_PATH}"
  echo "  Output: ${OUTPUT_DIR}"
  echo "--------------------------------------------------------------"
  
  Rscript "${SCRIPT}" \
    --cell_type       "${CT_NAME}" \
    --archr_project   "${ARCHR_PATH}" \
    --output_dir      "${OUTPUT_DIR}" \
    --gwas_catalog    "${GWAS_CATALOG}" \
    --chain_file      "${CHAIN_FILE}" \
    --threads         "${THREADS}" \
    --scavenge_cores  "${SCAVENGE_CORES}"
  
  echo "[$(date)] Finished: ${CT_NAME}"
  echo ""
done

echo "================================================================"
echo "[$(date)] All cell types complete"
echo "================================================================"

# Pool results
echo "Pooling results..."
Rscript pool_gwas_scavenge.R
echo "Done."
