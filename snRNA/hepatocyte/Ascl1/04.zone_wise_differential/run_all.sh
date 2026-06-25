#!/usr/bin/env bash
# ============================================================
# run_all.sh  -  export once, then MAST for global + all zones
# ============================================================
set -euo pipefail

# 1) export (needs your Python env with scanpy + adata path set in the script)
python 01_export_for_mast.py

# 2) MAST per subset (serial; run on a node with RAM headroom)
for subset in allfemhep periportal midlobular pericentral; do
    echo ">>> MAST: ${subset}"
    Rscript 02_run_mast.R "${subset}"
done

echo "All subsets done. Results in ascl1_de_by_zone_MAST/"
