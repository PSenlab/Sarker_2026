# SCENIC+ GRN Inference Pipeline (Mouse Liver Hepatocytes)

End-to-end SCENIC+ workflow that infers enhancer-driven gene regulatory networks (eGRNs) from the single-nucleus multi-ome liver aging dataset. The pipeline runs from raw fragments through pycisTopic topic modeling to SCENIC+ eRegulon inference and downstream trajectory analysis.

Reference: Sarker et al. (2026) *Cell Metabolism*.

## Pipeline overview

```
   fragments.tsv.gz + RNA cellranger matrix
                  |
     step1_pycistopic_preprocessing.py
        (per-sample QC, doublets, consensus
         peak calling via MACS2, cisTopic obj)
                  |
     step2_mallet_model.py
        (MALLET LDA fits at multiple topic counts)
                  |
     step3_pycistopic_lda_analysis.py
        (model selection, topic binarization,
         DARs per age group, region-set export)
                  |
     step4_prepare_scenicplus_input.py
        (harmonize RNA h5ad + cisTopic obj into
         SCENIC+-compatible inputs)
                  |
     step5_config.yaml + step5_run_scenicplus.py
        (Snakemake-driven motif enrichment,
         TF-to-gene + region-to-gene, eRegulon
         inference, AUCell)
                  |
     step6_downstream_analysis.py
     step6_loess_trajectory_plots.R
        (eRegulon dot plots, TF perturbation,
         age-trajectory LOESS, hub-TF networks)
```

## Scripts

| Step | Script | Description |
|---|---|---|
| 1 | `step1_pycistopic_preprocessing.py` | Per-sample QC, Scrublet doublet detection, cell annotation import, MACS2 consensus peak calling, cisTopic object construction |
| 2 | `step2_mallet_model.py` | MALLET LDA topic modeling (fits n_topics = 300, 350, 400 by default) |
|   | `run_mallet_model.sh` | SLURM launcher with Java/MALLET environment setup for step 2 |
| 3 | `step3_pycistopic_lda_analysis.py` | Model selection, topic binarization (Otsu + top-3k), age DAR calling, region-set BED export |
| 4 | `step4_prepare_scenicplus_input.py` | Harmonizes RNA AnnData + cisTopic object into SCENIC+ input files |
| 5 | `step5_config.yaml` | SCENIC+ Snakemake configuration (paths, motif parameters, inference thresholds) |
|   | `step5_run_scenicplus.py` | Snakemake launcher (wraps `snakemake --cores N`) |
|   | `run_scenicplus.sh` | SLURM launcher for the Snakemake pipeline |
| 6 | `step6_downstream_analysis.py` | eRegulon dot plots, motif hit reports, TF perturbation simulations, hub-TF module networks |
|   | `step6_loess_trajectory_plots.R` | LOESS-smoothed TF/target trajectories across age groups |

## Quick-start

```bash
# --- Step 1: pycisTopic preprocessing ---
# Edit paths at the top of the script, then:
python step1_pycistopic_preprocessing.py

# --- Step 2: MALLET LDA ---
# Update MALLET_BIN / MALLET_HOME / JAVA_HOME in run_mallet_model.sh
sbatch --partition=largemem --cpus-per-task=30 --mem=560g --time=12:00:00 \
       run_mallet_model.sh

# --- Step 3: LDA model selection + DARs + region sets ---
python step3_pycistopic_lda_analysis.py

# --- Step 4: Prepare SCENIC+ input ---
python step4_prepare_scenicplus_input.py \
    --rna_anndata integrated_scvi.h5ad \
    --cistopic_obj cisTopicObject_lda_hep_all_cluster.pkl \
    --output_dir scenicplus_input/

# --- Step 5: SCENIC+ Snakemake pipeline ---
# Initialize and configure:
scenicplus init_snakemake --out_dir scplus_pipeline
cp step5_config.yaml scplus_pipeline/Snakemake/config/config.yaml
# Edit input_data: paths inside the copied config, then:
sbatch run_scenicplus.sh

# --- Step 6: Downstream analysis ---
python step6_downstream_analysis.py
Rscript step6_loess_trajectory_plots.R
```

## Required external resources

| Resource | Where to get it |
|---|---|
| MALLET 202108 (or compatible) | http://mallet.cs.umass.edu/download.php |
| Java 8 JDK (MALLET dependency) | System package or OpenJDK 1.8 |
| mm10 blacklist | https://github.com/Boyle-Lab/Blacklist/blob/master/lists/mm10-blacklist.v2.bed.gz |
| mm10 SCREEN v10 motif database | https://resources.aertslab.org/cistarget/databases/ |
| Motif annotations (`motifs-v10nr_clust-nr.mgi`) | https://resources.aertslab.org/cistarget/motif2tf/ |

All paths to these resources are parameterized at the top of each script / config file as `path/to/...` placeholders - update them to your local installation before running.

## Input data

All scripts assume the canonical multi-ome object for this paper:

- `integrated_scvi.h5ad` - scVI-integrated AnnData with `celltype`, `celltype2`, `sex`, `age`, `sample` columns in `.obs`

Additionally, step 1 expects per-sample fragment files (`fragments.tsv.gz`) and a cell annotation TSV mapping barcodes to broad cell types.

## Compute requirements

| Step | Memory | CPUs | Time |
|---|---|---|---|
| Step 1 | ~60 GB | 20 | ~4 h |
| Step 2 (MALLET) | ~500 GB | 30+ | ~8-12 h |
| Step 3 | ~60 GB | 20 | ~2 h |
| Step 4 | ~30 GB | 4 | ~30 min |
| Step 5 (SCENIC+) | ~1.5 TB | 50+ | ~24-48 h |
| Step 6 | ~60 GB | 8 | ~1 h |

MALLET and the SCENIC+ inference step are the expensive ones. MALLET is memory-bound (step 2), SCENIC+ is both memory- and compute-bound (step 5). All other steps are feasible on a standard node.

## Output structure

```
outs/
├── qc/                                      # Step 1 per-sample QC plots
├── individual_cistopic_objects/             # Step 1 per-sample cisTopic obj
├── consensus_peak_calling/                  # Step 1 MACS2 consensus peaks
├── cisTopicObject_filtered_annotated.pkl    # Step 1 final obj
├── mal_result/mallet_model_<N>.pkl          # Step 2 LDA models
├── plots/                                   # Step 3 model-selection + DAR plots
├── region_sets/                             # Step 3 binarized topics / DARs (BED)
├── cisTopicObject_lda_hep_all_cluster.pkl   # Step 3 final annotated obj
scenicplus_input/
├── Hep_all_GEX_anndata_all.h5ad             # Step 4
├── Hep_all_cisTopic_obj.pkl                 # Step 4
└── region_sets_heps_all/                    # Step 4 (copied from Step 3)
scplus_pipeline/Snakemake/outs/
├── scplusmdata.h5mu                         # Step 5 final SCENIC+ mudata
├── eRegulon_direct.tsv                      # Step 5 direct eRegulons
├── eRegulons_extended.tsv                   # Step 5 extended eRegulons
├── AUCell_direct.h5mu                       # Step 5 AUCell scores
└── ...
figures/                                     # Step 6 dot plots, networks, trajectories
results/                                     # Step 6 tables
```

## Dependencies

### Python (step 1-6 Python)

```bash
# Core SCENIC+ env (contains pycisTopic, scenicplus, chromVAR, etc.)
mamba create -n scenicplus -c conda-forge -c bioconda python=3.11
mamba activate scenicplus
pip install scenicplus pycisTopic
pip install scanpy anndata polars pyranges scrublet ray joblib \
            gseapy mygene seaborn matplotlib
```

### R (step 6 LOESS plots)

```r
install.packages(c("ggplot2", "dplyr", "tidyr", "readr", "purrr"))
```

## Notes

- **Reproducibility**: All stochastic steps use fixed seeds (MALLET `random_state=555`, SCENIC+ `seed=666`).
- **Path placeholders**: Every `path/to/...` string in the scripts and config must be updated to match your local installation before running.
- **The SCENIC+ Snakemake folder is not included** in this repo - it is generated fresh on each run via `scenicplus init_snakemake`. Only the `step5_config.yaml` template is tracked.

## Reference

Bravo Gonzalez-Blas C, et al. (2023). SCENIC+: single-cell multiomic inference of enhancers and gene regulatory networks. *Nature Methods* 20, 1355-1367.

Bravo Gonzalez-Blas C, et al. (2019). cisTopic: cis-regulatory topic modeling on single-cell ATAC-seq data. *Nature Methods* 16, 397-400.
