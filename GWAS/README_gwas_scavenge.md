# GWAS-SCAVENGE: Cell type specific GWAS trait enrichment in mouse liver

Integration of human GWAS Catalog variants with single-nucleus multi-ome (RNA + ATAC-seq) data from the aging mouse liver to identify cell-type-specific enrichment of liver-disease-associated genetic variants using the [SCAVENGE](https://github.com/sankaranlab/SCAVENGE) framework.

## Overview

This pipeline maps human GWAS SNPs to the mouse genome via UCSC liftOver, overlaps them with cell-type-specific Peak-to-Gene (P2G) linked regulatory elements derived from ArchR, and uses SCAVENGE's network propagation approach to identify cells enriched for trait-associated chromatin accessibility.

### Pipeline overview

```
GWAS Catalog (hg38) -> liftOver -> mm10 SNPs -> LD expansion (+/- 25 kb)
                                                    |
ArchR P2G peaks (per cell type) <- overlap -> Per-trait BED files
                                                    |
                                 SCAVENGE (per trait):
                                 chromVAR z-scores -> seed cells (top 5%)
                                 TF-IDF -> LSI (30 dims) -> mutual KNN (k=30)
                                 Random walk with restart (gamma=0.05)
                                 TRS = 95th-pct capped, min-max scaled,
                                       z-score-derived scale factor
                                 1,000 degree-matched permutations -> BH FDR
```

## Scripts

| Script | Description |
|---|---|
| `run_gwas_scavenge.R` | Main pipeline - runs SCAVENGE for a single cell type (parameterized) |
| `run_all_celltypes.sh` | Batch runner - launches the pipeline for all 9 cell types |
| `pool_gwas_scavenge.R` | Pools results across cell types into publication tables |
| `GWAS_download.R` | Helper to download the NHGRI-EBI GWAS Catalog and liftOver chain |

## Usage

### Single cell type

```bash
Rscript run_gwas_scavenge.R \
  --cell_type Hepatocyte \
  --archr_project path/to/Step8_Hepatocyte_CCAN_P2G \
  --output_dir path/to/output/hepatocyte \
  --gwas_catalog path/to/gwas_catalog_v1.0.2-associations.tsv \
  --chain_file path/to/hg38ToMm10.over.chain \
  --threads 60 \
  --scavenge_cores 55
```

### All cell types

```bash
# Update the path placeholders inside run_all_celltypes.sh, then:
nohup bash run_all_celltypes.sh > run_all.log 2>&1 &
```

### Pool results

```bash
# After all cell types finish:
Rscript pool_gwas_scavenge.R
```

## Parameters

| Parameter | Default | Description |
|---|---|---|
| `--cell_type` | *required* | Cell type label (e.g., Hepatocyte, Kupffer) |
| `--archr_project` | *required* | Path to ArchR project with P2G links |
| `--output_dir` | *required* | Output directory for this cell type |
| `--gwas_catalog` | *required* | NHGRI-EBI GWAS Catalog TSV |
| `--chain_file` | *required* | hg38 -> mm10 liftOver chain file |
| `--threads` | 60 | ArchR threads |
| `--scavenge_cores` | 55 | Parallel cores for permutation test |
| `--peak_extension` | 500 | Peak extension (bp) for GWAS-peak overlap |
| `--ld_window` | 25000 | LD proxy expansion window (bp) |
| `--min_regions` | 2 | Minimum GWAS regions per trait |
| `--min_seed_cells` | 5 | Minimum seed cells for SCAVENGE |
| `--min_valid_cells` | 10 | Minimum valid cells after random walk |
| `--fdr_threshold` | 0.05 | BH FDR significance threshold |
| `--permutations` | 1000 | Permutation test iterations |
| `--seed` | 42 | Random seed |

## Output structure

```
<cell_type>/
└── gwas_p2g/Liver/
    ├── checkpoints/              # Resumable pipeline checkpoints
    │   ├── step2_SE_Data.rds     # SummarizedExperiment (P2G peaks x cells)
    │   ├── step2_peaks_p2g.rds
    │   ├── step3_gwas_filtered.rds
    │   ├── step3_liver_traits.rds
    │   ├── step4_gwas_mm10.rds   # liftOver results
    │   ├── step5_gwas_mm10_ld.rds
    │   ├── step6_trait_counts.rds
    │   └── step7_bed_files.rds
    ├── trait_beds/               # One BED file per trait
    └── results/
        ├── results_Liver_P2G_<CellType>_ALL.rds
        ├── summary_Liver_P2G_<CellType>.csv
        ├── skipped_Liver_P2G_<CellType>.csv
        ├── combined_ALL_raw_<CellType>.csv
        ├── combined_ALL_filtered_<CellType>.csv
        ├── raw/                  # Per-trait all-cell results
        └── filtered/             # Per-trait FDR < 0.05 results
```

### Pooled cross-celltype output

```
pooled_results/
├── pooled_ALL_filtered.csv                 # All FDR < 0.05 cells across celltypes
├── pooled_summary.csv                      # All trait-celltype summaries
├── matrix_pct_fdr_sig.csv                  # Traits x celltypes (% FDR-sig)
├── matrix_median_TRS.csv                   # Traits x celltypes (median TRS)
├── matrix_n_fdr_sig.csv                    # Traits x celltypes (n cells)
├── summary_by_celltype.csv                 # Stats per cell type
├── summary_by_trait.csv                    # Stats per trait
├── trait_celltype_specificity.csv          # Specificity scores
└── top_specific_traits_per_celltype.csv    # Top specific traits
```

## Cell types analyzed

| Cell type | ArchR project | Status |
|---|---|---|
| Hepatocyte | Step8_Hepatocyte_CCAN_P2G | included |
| Endothelial_01 | Step8_Endothelial_01_CCAN_P2G | included |
| Endothelial_02 | Step8_Endothelial_02_CCAN_P2G | included |
| Kupffer | Step8_Kupffer_CCAN_P2G | included |
| MoMFs | Step8_MoMFs_CCAN_P2G | included |
| Cholangiocyte_01 | Step8_Cholangiocyte_01_CCAN_P2G | included |
| Cholangiocyte_02 | Step8_Cholangiocyte_02_CCAN_P2G | included |
| Lymp_B | Step8_Lymp_B_CCAN_P2G | included |
| Lymp_T | Step8_Lymp_T_CCAN_P2G | included |
| Stellate | - | excluded (~169 P2G peaks, insufficient for SCAVENGE) |

## GWAS trait curation

Liver-relevant traits are identified from the GWAS Catalog using a two-stage filter (regex patterns defined in `run_gwas_scavenge.R`):

1. **Inclusion** (~66 patterns): core liver terms (`hepat`, `liver`), liver diseases (NAFLD, NASH, steatosis, cirrhosis), liver enzymes (ALT, AST, GGT, ALP, bilirubin), liver-synthesized proteins (albumin, fibrinogen, complement C, haptoglobin, transferrin, ceruloplasmin), lipoproteins (apolipoprotein, VLDL), iron metabolism (ferritin, hemochromatosis, hepcidin), viral hepatitis, autoimmune liver disease, drug-induced liver injury, biliary pathology.

2. **Exclusion** (~90 patterns): brain iron MRI (putamen, pallidum, etc.), gamma-glutamyl amino acid metabolites (not the GGT enzyme), VLDL subfractions, APOE / Alzheimer's, urinary albumin, CKD / kidney disease, cystic and pulmonary fibrosis, non-liver autoimmune diseases, and several miscellaneous non-liver traits.

## Features

- **Checkpoint / resume**: every major step saves RDS checkpoints; interrupted runs resume automatically.
- **Parameterized**: all paths and thresholds configurable via command-line flags.
- **Per-trait outputs**: both raw (all cells) and filtered (FDR < 0.05) CSVs per trait.
- **Cross-cell-type pooling**: matrices ready for heatmap or dot-plot visualization.

## Dependencies

```r
# Bioconductor
BiocManager::install(c(
  "ArchR", "GenomicRanges", "SummarizedExperiment",
  "BSgenome.Mmusculus.UCSC.mm10", "chromVAR", "rtracklayer"
))

# GitHub
devtools::install_github("sankaranlab/SCAVENGE")
devtools::install_github("caleblareau/gchromVAR")

# CRAN
install.packages(c("data.table", "dplyr", "tidyr", "Matrix"))
```

## External data

- NHGRI-EBI GWAS Catalog associations TSV: https://www.ebi.ac.uk/gwas/docs/file-downloads
- UCSC hg38 to mm10 liftOver chain: https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToMm10.over.chain.gz

Both can be downloaded via the companion `GWAS_download.R` helper.

## Reference

Yu W, et al. (2022). SCAVENGE: Single Cell Analysis of Variant Enrichment through Network propagation of GEnetic associations. *Nature Biotechnology* 40, 1443-1450.
