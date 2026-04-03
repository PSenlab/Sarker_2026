# Hepatomap: A Multimodal Atlas of the Aging Mouse Liver
## Overview

This repository contains the complete analysis code for the manuscript:

**A lifespan-resolved multi-omic study reveals dynamic transcription factor networks and the regulatory logic of liver aging
> Nishat Sarker, Na Yang, Michel Bernier, Nirad Banskota, Minjung Kwon, Yaohui Chen, Amit Singh, Robert W. Maul, Sadia Afrin, Lin Wang, Miguel Aon, Nathan L. Price, Chang-Yi Cui, Jinshui Fan, Supriyo De, Mary Kaileh, Jeffrey Albrecht, Ranjan Sen, Myriam Gorospe, Rafael de Cabo, Payel Sen1.
> 2026.

We generated paired single-nucleus RNA-seq and ATAC-seq (10X Multiome) from **40 mouse liver samples** spanning **five age groups** (young, mid-age, old, pre-geriatric, geriatric) and **both sexes**, yielding ~300,000 high-quality nuclei with paired modalities.

## Repository Structure

```
liver-aging-multiome/
├── README.md
├── LICENSE
├── .gitignore
├── environment.yml                     # Conda environment (Python)
├── renv.lock                           # R package versions (renv)
├── requirements.txt                    # Python pip dependencies
│
├── config/
│   ├── config.py                       # Python: paths, palettes, constants
│   └── config.R                        # R: paths, palettes, constants
│
├── scripts/
│   ├── 01_preprocessing/
│   │   ├── 01_cellranger_arc.sh        # Cell Ranger Arc alignment
│   │   ├── 02_soupx_decontamination.R  # Ambient RNA removal
│   │   └── 03_qc_filtering.py          # QC and nuclei filtering
│   │
│   ├── 02_snrna_processing/
│   │   ├── 01_scvi_integration.py      # scVI batch integration
│   │   ├── 02_clustering_annotation.py # Leiden clustering & cell-type annotation
│   │   └── 03_pseudobulk_deseq2.R     # DESeq2 LRT across age (sex covariate)
│   │
│   ├── 03_snatac_processing/
│   │   ├── 01_archr_preprocessing.R    # ArchR 6-step pipeline
│   │   ├── 02_peak_calling.R           # MACS2 peak calling per group
│   │   └── 03_motif_enrichment.R       # chromVAR motif deviations
│   │
│   ├── 04_wnn_integration/
│   │   ├── 01_muon_wnn.py             # Weighted nearest neighbors (Muon)
│   │   └── 02_schard_conversion.py    # h5ad ↔ ArchR conversion (schard)
│   │
│   ├── 05_compartment_switching/
│   │   ├── 01_ab_compartment_calling.R # A/B compartment eigen decomposition
│   │   ├── 02_hmm_stability_model.R   # 5-class HMM stability classification
│   │   └── 03_compartment_figures.R    # Heatmaps, alluvial plots, summaries
│   │
│   ├── 06_scenic_plus/
│   │   ├── 01_scenic_plus_grn.py      # SCENIC+ gene regulatory networks
│   │   └── 02_eregulon_analysis.py    # eRegulon scoring & visualization
│   │
│   ├── 07_gwas_scavenge/
│   │   ├── 01_gwas_scavenge.R         # SCAVENGE enrichment (liver/metabolic/aging/immune)
│   │   └── 02_ewas_integration.py     # EWAS CpG overlap analysis
│   │
│   ├── 08_peak_gene_linkage/
│   │   ├── 01_peak_to_gene.R          # ArchR peak-to-gene correlation
│   │   └── 02_ccans.R                 # Cis-coaccessibility networks
│   │
│   ├── 09_celltype_composition/
│   │   ├── 01_celltype_proportions.py # Sample-level proportions & ANOVA
│   │   └── 02_composition_figures.py  # Boxplots, stacked bars, UMAPs
│   │
│   ├── 10_hepatocyte_zonation/
│   │   ├── 01_zonation_assignment.py  # Periportal/Midlobular/Pericentral scoring
│   │   ├── 02_zonation_figures.py     # Zonation UMAPs, boxplots, stacked bars
│   │   └── 03_gene_module_scoring.py  # Zonation gene-set heatmaps
│   │
│   ├── 11_transcriptional_entropy/
│   │   └── 01_augur_entropy.py        # Augur cell-type prioritization & entropy
│   │
│   └── 12_ascl1_bcl6_analysis/
│       ├── 01_ascl1_female_eregulon.py # Ascl1 female-specific eRegulon
│       ├── 02_bcl6_male_eregulon.py    # Bcl6 male hepatocyte aging signal
│       └── 03_human_validation.py      # MASLD cohort & TCGA-LIHC validation
│
├── data/                               # Input data (not tracked by git)
│   └── .gitkeep
├── figures/                            # Generated figures (not tracked)
│   └── .gitkeep
└── results/                            # Generated tables/CSVs (not tracked)
    └── .gitkeep
```

## Setup

### Python environment

```bash
git clone https://github.com/SenLab-NIA/liver-aging-multiome.git
cd liver-aging-multiome

# Option A: Conda (recommended)
conda env create -f environment.yml
conda activate liver_multiome

# Option B: pip
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
```

### R environment

```r
# Install renv if needed
install.packages("renv")
renv::restore()
```

## Data Availability

- **Raw and processed data**: GEO accession [GSExxxxxx]
- **Interactive browser**: [Hepatomap Shiny App URL]
- **Processed AnnData**: Zenodo DOI [10.5281/zenodo.xxxxxxx]

### Required input files

| File | Description | Location |
|------|-------------|----------|
| `adata_wnn.h5ad` | WNN-integrated AnnData (RNA + ATAC) | `data/` |
| `adata_rna.h5ad` | snRNA-seq processed AnnData | `data/` |
| ArchR project | ArchR SavedProject directory | `data/archr_project/` |
| `pseudobulk_counts/` | Pseudobulk count matrices per cell type | `data/` |
| GWAS summary statistics | Formatted GWAS files | `data/gwas/` |

## Running the Pipeline

Scripts are numbered and can be run sequentially:

```bash
# Example: run all Python scripts in order
for dir in scripts/*/; do
    for script in "$dir"*.py; do
        [ -f "$script" ] && python "$script"
    done
done

# Example: run all R scripts
for dir in scripts/*/; do
    for script in "$dir"*.R; do
        [ -f "$script" ] && Rscript "$script"
    done
done
```

Or run individual analyses:

```bash
python scripts/09_celltype_composition/01_celltype_proportions.py
Rscript scripts/05_compartment_switching/01_ab_compartment_calling.R
```

## Key Parameters

| Parameter | Value | Description |
|-----------|-------|-------------|
| Tile size | 80 kb | Compartment calling resolution |
| scVI latent dims | 30 | Integration dimensionality |
| Leiden resolution | 0.6 | Primary clustering resolution |
| Peak calling | MACS2 q < 0.01 | Reproducible peak set |
| DESeq2 design | ~ sex + age (LRT) | Pseudobulk differential expression |
| HMM states | 5 | Stable_active, Stable_repressive, AR, RA, Non_nonotonic |

## Hepatocyte Zonation Mapping

| Cluster | Zonation | Key Markers |
|---------|----------|-------------|
| Hep-01 | Periportal | *Sds*, *Hal*, *Cyp2f2* |
| Hep-02 | Midlobular | *Hamp*, *Igfbp2* |
| Hep-03 | Pericentral | *Glul*, *Oat*, *Cyp2e1* |
| Hep-04 | Pericentral | *Cyp1a2*, *Aldh1a1* |
| Hep-05 | Periportal | *Arg1*, *Ass1* |
| Hep-06 | Midlobular | *Mup20*, *Serpina1e* |
| Hep-07 | Midlobular | *Apoa4*, *Ttr* |

## Contact

- **Nishat Sarker** — [nishat.sarker@nih.gov](mailto:nishat.sarker@nih.gov)
- **Payel Sen** (PI) — [payel.sen@nih.gov](mailto:payel.sen@nih.gov)
- Functional Epigenomics Unit, Laboratory of Genetics and Genomics, NIA/NIH
