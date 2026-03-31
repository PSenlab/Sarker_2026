# GWAS-SCAVENGE: Cell-Type-Specific GWAS Trait Enrichment in Mouse Liver

Integration of human GWAS Catalog variants with single-cell multiome (RNA + ATAC-seq) data to identify cell-type-specific enrichment of liver disease-associated genetic variants using the [SCAVENGE](https://github.com/sankaranlab/SCAVENGE) framework.

## Overview

This pipeline maps human GWAS SNPs to the mouse genome via liftOver, overlaps them with cell-type-specific Peak-to-Gene (P2G) linked regulatory elements from ArchR, and uses SCAVENGE's network propagation approach to identify cells enriched for trait-associated chromatin accessibility.

### Pipeline Steps

```
GWAS Catalog (hg38) → LiftOver → mm10 SNPs → LD expansion (±25kb)
                                                    ↓
ArchR P2G peaks (cell-type) ← Overlap → Trait BED files
                                                    ↓
                              SCAVENGE (per trait):
                                chromVAR z-scores → Seed cells
                                TF-IDF → LSI → KNN graph
                                Random walk → TRS scores
                                Permutation test → FDR
```

## Scripts

| Script | Description |
|---|---|
| `run_gwas_scavenge.R` | Main pipeline — runs SCAVENGE for one cell type (parameterized) |
| `run_all_celltypes.sh` | Batch runner — launches pipeline for all 9 cell types |
| `pool_gwas_scavenge.R` | Pools results across cell types into publication tables |

## Usage

### Single cell type

```bash
Rscript run_gwas_scavenge.R \
  --cell_type Hepatocyte \
  --archr_project /path/to/Step8_Hepatocyte_CCAN_P2G \
  --output_dir /path/to/output/hepatocyte \
  --gwas_catalog /path/to/gwas_catalog_v1.0.2-associations.tsv \
  --chain_file /path/to/hg38ToMm10.over.chain \
  --threads 60 \
  --scavenge_cores 55
```

### All cell types

```bash
# Edit paths in run_all_celltypes.sh, then:
nohup bash run_all_celltypes.sh > run_all.log 2>&1 &
```

### Pool results

```bash
# After all cell types complete:
Rscript pool_gwas_scavenge.R
```

## Parameters

| Parameter | Default | Description |
|---|---|---|
| `--cell_type` | *required* | Cell type label (e.g., Hepatocyte, Kupffer) |
| `--archr_project` | *required* | Path to ArchR project with P2G links |
| `--output_dir` | *required* | Output directory for this cell type |
| `--gwas_catalog` | *required* | NHGRI-EBI GWAS Catalog TSV |
| `--chain_file` | *required* | hg38→mm10 liftOver chain file |
| `--threads` | 60 | ArchR threads |
| `--scavenge_cores` | 55 | Parallel cores for permutation test |
| `--peak_extension` | 500 | Peak extension for overlap (bp) |
| `--ld_window` | 25000 | LD proxy expansion window (bp) |
| `--min_regions` | 2 | Minimum GWAS regions per trait |
| `--fdr_threshold` | 0.05 | FDR significance threshold |
| `--permutations` | 1000 | Permutation test iterations |

## Output Structure

```
<cell_type>/
└── gwas_p2g/Liver/
    ├── checkpoints/              # Resumable pipeline checkpoints
    │   ├── step2_SE_Data.rds     # SummarizedExperiment (P2G peaks × cells)
    │   ├── step2_peaks_p2g.rds   # P2G peak GRanges
    │   ├── step3_gwas_filtered.rds
    │   ├── step3_liver_traits.rds
    │   ├── step4_gwas_mm10.rds   # LiftOver results
    │   ├── step5_gwas_mm10_ld.rds # LD-expanded regions
    │   ├── step6_trait_counts.rds # Per-trait overlap counts
    │   └── step7_bed_files.rds   # BED file paths
    ├── trait_beds/                # One BED per trait
    └── results/
        ├── results_Liver_P2G_<CellType>_ALL.rds   # Full results (list of DFs)
        ├── summary_Liver_P2G_<CellType>.csv        # Per-trait summary
        ├── skipped_Liver_P2G_<CellType>.csv        # Failed traits
        ├── combined_ALL_raw_<CellType>.csv         # All cells combined
        ├── combined_ALL_filtered_<CellType>.csv    # FDR < 0.05 cells
        ├── raw/                   # Per-trait all-cell results
        └── filtered/              # Per-trait FDR < 0.05 results
```

### Pooled output (from `pool_gwas_scavenge.R`)

```
pooled_results/
├── Table_S1_celltype_summary.csv         # Per-cell-type overview
├── Table_S2_trait_celltype_detail.csv     # Trait × cell type results
├── matrix_pct_fdr_sig.csv                # Heatmap: % FDR < 0.05
├── matrix_median_TRS.csv                 # Heatmap: median TRS
├── matrix_n_fdr_sig.csv                  # Heatmap: N significant cells
├── matrix_median_zscore.csv              # Heatmap: median z-score
├── summary_by_trait.csv                  # Per-trait cross-celltype summary
├── trait_celltype_specificity.csv        # Specificity scores
└── top_specific_traits_per_celltype.csv  # Top cell-type-specific traits
```

## Cell Types Analyzed

| Cell Type | ArchR Project | Status |
|---|---|---|
| Hepatocyte | Step8_Hepatocyte_CCAN_P2G | ✓ |
| Endothelial_01 | Step8_Endothelial_01_CCAN_P2G | ✓ |
| Endothelial_02 | Step8_Endothelial_02_CCAN_P2G | ✓ |
| Kupffer | Step8_Kupffer_CCAN_P2G | ✓ |
| MoMFs | Step8_MoMFs_CCAN_P2G | ✓ |
| Cholangiocyte_01 | Step8_Cholangiocyte_01_CCAN_P2G | ✓ |
| Cholangiocyte_02 | Step8_Cholangiocyte_02_CCAN_P2G | ✓ |
| Lymp_B | Step8_Lymp_B_CCAN_P2G | ✓ |
| Lymp_T | Step8_Lymp_T_CCAN_P2G | ✓ |
| Stellate | — | Excluded (169 P2G peaks) |

## GWAS Trait Curation

Liver-relevant traits are identified from the GWAS Catalog using a two-stage filter:

1. **Inclusion** (66 regex patterns): Core liver terms, liver diseases, enzymes (ALT, AST, GGT, ALP), liver-synthesized proteins, lipoproteins, iron metabolism, viral hepatitis, autoimmune liver disease, biliary pathology
2. **Exclusion** (90+ regex patterns): Brain iron MRI, gamma-glutamyl amino acid metabolites, VLDL subfractions, APOE/Alzheimer's, urinary albumin, CKD/kidney, cystic fibrosis, non-liver autoimmune diseases

## Features

- **Checkpoint/resume**: Every major step saves RDS checkpoints; interrupted runs resume automatically
- **Parameterized**: All paths and thresholds configurable via command line
- **Per-trait outputs**: Both raw (all cells) and filtered (FDR < 0.05) CSVs per trait
- **Cross-cell-type pooling**: Matrices ready for heatmap visualization

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

## Reference

Yu et al. (2022) SCAVENGE: Single Cell Analysis of Variant Enrichment through Network propagation of GEnetic associations. *Nature Biotechnology*.
