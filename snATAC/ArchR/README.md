# ArchR Preprocessing Pipeline for Single-Cell Multiome ATAC-seq

[![R](https://img.shields.io/badge/R-%3E%3D4.0-blue)](https://www.r-project.org/)
[![ArchR](https://img.shields.io/badge/ArchR-%3E%3D1.0.2-green)](https://www.archrproject.com/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Overview

This repository contains the preprocessing pipeline for single-cell ATAC-seq data from a 10x Genomics Multiome experiment, designed to study chromatin accessibility changes during liver aging in mice.

## Study Design

| Parameter | Value |
|-----------|-------|
| **Species** | Mouse (*Mus musculus*) |
| **Genome** | mm10 |
| **Tissue** | Liver |
| **Age Groups** | Young, Middle-age, Old, Pre-geriatric, Geriatric |
| **Replicates** | 8 biological replicates per group |
| **Total Samples** | 40 |

## Pipeline Steps

```
┌─────────────────────────────────────────────────────────────────────────┐
│                        PREPROCESSING PIPELINE                           │
├─────────────────────────────────────────────────────────────────────────┤
│                                                                         │
│  Step 1: Project Creation                                               │
│     └── Create ArchR project from Arrow files                           │
│     └── Add age group metadata                                          │
│                                                                         │
│  Step 2: Doublet Filtering                                              │
│     └── k = 10, knnMethod = "UMAP", LSIMethod = 1                       │
│     └── Filter detected doublets                                        │
│                                                                         │
│  Step 3: Dimensionality Reduction (Iterative LSI)                       │
│     └── TileMatrix, 2 iterations                                        │
│     └── 25,000 variable features                                        │
│     └── Dimensions 1-30                                                 │
│                                                                         │
│  Step 4: Batch Correction (Harmony)                                     │
│     └── Correct for age group batch effects                             │
│                                                                         │
│  Step 5: Seurat Metadata Integration                                    │
│     └── Transfer cell type annotations                                  │
│     └── Transfer sample, age, sex metadata                              │
│                                                                         │
│  Step 6: WNN UMAP Transfer                                              │
│     └── Import multiome-derived UMAP embedding                          │
│                                                                         │
└─────────────────────────────────────────────────────────────────────────┘
```

## Requirements

### Software

- R ≥ 4.0
- ArchR ≥ 1.0.2
- Seurat ≥ 4.0
- harmony
- dplyr
- stringr
- EnsDb.Mmusculus.v79

### Installation

```r
# Install ArchR
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
devtools::install_github("GreenleafLab/ArchR")

# Install other dependencies
install.packages(c("Seurat", "dplyr", "stringr"))
BiocManager::install("EnsDb.Mmusculus.v79")
```

## Usage

### 1. Prepare Input Files

Ensure the following are available:
- Pre-generated Arrow files (one per sample)
- Integrated Seurat object with cell type annotations and WNN UMAP

### 2. Configure Paths

Edit the following variables in the script:

```r
WORK_DIR <- "/path/to/arrow/files"
OUTPUT_BASE <- "/path/to/output"
SEURAT_OBJ_PATH <- "/path/to/seurat_object.rds"
```

### 3. Run Pipeline

```bash
Rscript ArchR_Preprocessing_Pipeline.R
```

## Parameters

### Doublet Detection

| Parameter | Value | Description |
|-----------|-------|-------------|
| k | 10 | Number of nearest neighbors |
| knnMethod | "UMAP" | Method for kNN calculation |
| LSIMethod | 1 | LSI method for dimensionality reduction |

### Iterative LSI

| Parameter | Value | Description |
|-----------|-------|-------------|
| useMatrix | "TileMatrix" | Matrix type for LSI |
| iterations | 2 | Number of LSI iterations |
| varFeatures | 25,000 | Number of variable features |
| dimsToUse | 1:30 | Dimensions to retain |
| resolution | 0.1 | Clustering resolution |
| sampleCells | 10,000 | Cells sampled per iteration |

### Harmony Batch Correction

| Parameter | Value | Description |
|-----------|-------|-------------|
| reducedDims | "IterativeLSI" | Input dimensionality reduction |
| groupBy | "DatasetGroup" | Variable for batch correction |

## Output Structure

```
ArchR_Projects/
├── Step1_Project/
│   └── ArchRProject.rds
├── Step2_DoubletsFiltered/
│   └── ArchRProject.rds
├── Step3_LSI/
│   └── ArchRProject.rds
├── Step4_Harmony/
│   └── ArchRProject.rds
├── Step5_MetadataAdded/
│   └── ArchRProject.rds
└── Step6_Xwnn_UMAP/
    └── ArchRProject.rds          # Final preprocessed project
```

## Metadata Columns

The final ArchR project contains the following cell-level metadata:

| Column | Description |
|--------|-------------|
| Sample | Sample identifier |
| DatasetGroup | Age group (Young, Middle_age, Old, Pre_Geriatric, Geriatric) |
| celltype | Cell type annotation (from Seurat) |
| celltype2 | Alternative cell type annotation |
| age | Age information |
| sex | Sex (M/F) |
| Xwnn_UMAP_1 | WNN UMAP dimension 1 |
| Xwnn_UMAP_2 | WNN UMAP dimension 2 |

## Next Steps

After preprocessing, the ArchR project is ready for:

1. **Peak Calling** - `addGroupCoverages()` + `addReproduciblePeakSet()`
2. **Differential Accessibility** - `getMarkerFeatures()`
3. **Motif Enrichment** - `addMotifAnnotations()` + `enrichMotifs()`
4. **Peak-to-Gene Linkage** - `addPeak2GeneLinks()`
5. **Trajectory Analysis** - `addTrajectory()`

## Citation

If you use this pipeline, please cite:

- **ArchR**: Granja JM et al. ArchR is a scalable software package for integrative single-cell chromatin accessibility analysis. *Nature Genetics* (2021). https://doi.org/10.1038/s41588-021-00790-6

- **Harmony**: Korsunsky I et al. Fast, sensitive and accurate integration of single-cell data with Harmony. *Nature Methods* (2019). https://doi.org/10.1038/s41592-019-0619-0

- **Seurat**: Hao Y et al. Integrated analysis of multimodal single-cell data. *Cell* (2021). https://doi.org/10.1016/j.cell.2021.04.048

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contact

[Your Name]  
[Your Email]  
[Your Institution]

---

*This pipeline was developed as part of a study investigating chromatin accessibility dynamics during liver aging.*
