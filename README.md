# Hepatomap: A Multimodal Atlas of the Aging Mouse Liver
## Overview

This repository contains the complete analysis code for the manuscript:

A lifespan-resolved multi-omic study reveals dynamic transcription factor networks and the regulatory logic of liver aging
> Nishat Sarker, Na Yang, Michel Bernier, Nirad Banskota, Minjung Kwon, Yaohui Chen, Amit Singh, Robert W. Maul, Sadia Afrin, Lin Wang, Miguel Aon, Nathan L. Price, Chang-Yi Cui, Jinshui Fan, Supriyo De, Mary Kaileh, Jeffrey Albrecht, Ranjan Sen, Myriam Gorospe, Rafael de Cabo, Payel Sen.


We generated paired single-nucleus RNA-seq and ATAC-seq (10X Multiome) from **40 mouse liver samples** spanning **five age groups** (young, mid-age, old, pre-geriatric, geriatric) and **both sexes**, yielding ~300,000 high-quality nuclei with paired modalities.

## Repository Structure



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



## Contact

- **Nishat Sarker** — [nishat.sarker@nih.gov](mailto:nishat.sarker@nih.gov)
- **Payel Sen** (PI) — [payel.sen@nih.gov](mailto:payel.sen@nih.gov)
- Functional Epigenomics Unit, Laboratory of Genetics and Genomics, NIA/NIH
