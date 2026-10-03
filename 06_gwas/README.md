# GWAS-SCAVENGE: Cell type specific GWAS trait enrichment in mouse liver

Integration of human GWAS Catalog variants with single-nucleus multi-ome (RNA + ATAC-seq) data from the aging mouse liver to identify cell-type-specific enrichment of liver-disease-associated genetic variants using the [SCAVENGE](https://github.com/sankaranlab/SCAVENGE) framework.

## Overview

This pipeline maps human GWAS SNPs to the mouse genome via UCSC liftOver, overlaps them with cell-type-specific Peak-to-Gene (P2G) linked regulatory elements derived from ArchR, and uses SCAVENGE's network propagation approach to identify cells enriched for trait-associated chromatin accessibility.

## Usage

### Single cell type

```bash
Rscript 01_run_gwas_scavenge.R \
  --cell_type hepatocyte \
  --archr_project path/to/Step11b_hepatocyte_P2G \
  --output_dir path/to/output/hepatocyte \
  --gwas_catalog path/to/gwas_catalog_v1.0.2-associations.tsv \
  --chain_file path/to/hg38ToMm10.over.chain \
  --threads 60 \
  --scavenge_cores 55
```

## External data

- NHGRI-EBI GWAS Catalog associations TSV: https://www.ebi.ac.uk/gwas/docs/file-downloads
- UCSC hg38 to mm10 liftOver chain: https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToMm10.over.chain.gz

## Reference

Yu W, et al. (2022). SCAVENGE: Single Cell Analysis of Variant Enrichment through Network propagation of GEnetic associations. *Nature Biotechnology* 40, 1443-1450.
