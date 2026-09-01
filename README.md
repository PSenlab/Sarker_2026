# A lifespan-resolved multi-omic study reveals dynamic transcription factor networks and the regulatory logic of liver aging
## Overview

This repository contains the complete analysis code for the manuscript:

A lifespan-resolved multi-omic study reveals dynamic transcription factor networks and the regulatory logic of liver aging
> Nishat Sarker, Na Yang, Sadia Afrin, Michel Bernier, Nirad Banskota, Minjung Kwon, Yaohui Chen, Amit Singh, Nathan L. Price, Robert W. Maul, Lin Wang, Miguel Aon, Chang-Yi Cui, Jinshui Fan, Supriyo De, Mary Kaileh, Jeffrey Albrecht, Ranjan Sen, Myriam Gorospe, Rafael de Cabo, Payel Sen.


We generated paired single-nucleus RNA-seq and ATAC-seq (10X Multiome) from **40 mouse liver samples** spanning **five age groups** (young, mid-age, old, pre-geriatric, geriatric) and **both sexes**, yielding ~300,000 high-quality nuclei with paired modalities. Additionally, ASCL1 ChIP-seq data was generated from **three biological replicates** in pre-geriatric mice of **both sexes**.

## Key analyses:

- Integrated RNA + ATAC atlas 
- Age-associated differential expression and chromatin accessibility (pseudobulk DESeq2 LRT, k-means clustering)
- Chromatin compartment dynamics (HMM-based stability classification)
- Peak-to-gene linkage via ArchR CCAN co-accessibility
- SCENIC+ eGRN inference with hub TF identification 
- Cell-type-resolved GWAS trait enrichment (SCAVENGE)
- Transcriptional entropy, cell-cycle (Tricycle), and senescence scoring
- Human MASLD validation across four cohorts (Duke, German, Japanese, European)
- ChIP-seq for ASCL1 in male and female pre-geriatric livers


## Data Availability

- **Raw and processed data**: GEO accession [GSE341533]
- **Interactive browser**: [https://sciapps.nia.nih.gov/hepatomap/]
- **Processed AnnData**: Zenodo DOI [10.5281/zenodo.22218673]


## Contact

- **Nishat Sarker** — [nishat.sarker@nih.gov](mailto:nishat.sarker@nih.gov)
- **Payel Sen** (PI) — [payel.sen@nih.gov](mailto:payel.sen@nih.gov)
- Functional Epigenomics Unit, Laboratory of Genetics and Genomics, NIA/NIH
