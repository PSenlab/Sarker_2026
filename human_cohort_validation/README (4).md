# Human Cohort Validation (MASLD / MASH / HCC)

Independent validation of the mouse Ascl1 sex-dimorphic aging signal across 8 publicly available human liver datasets spanning healthy reference tissue (GTEx), MASLD, MASH, bariatric-surgery intervention, statin treatment, and HCC. Each script runs the same core comparison - ASCL1 expression in male vs female livers via Mann-Whitney / Wilcoxon rank-sum test, with a bar + jitter + significance-bracket plot - on a different dataset, with dataset-specific preprocessing handled per script.

Reference: Sarker et al. (2026) *Cell Metabolism*.

## Scripts

| Script | Cohort | Platform | Sex source | n |
|---|---|---|---|---|
| `validation_gtex_liver.R` | GTEx v10 (healthy liver reference) | Bulk RNA-seq (GCT) | Annotated (GTEx phenotype) | ~200 |
| `validation_german_GSE33814.py` | German MASLD | Illumina HumanWG-6 v3 | Inferred (XIST + RPS4Y1) | ~60 |
| `validation_european_GSE135251.py` | European MASLD | RNA-seq (per-sample counts) | Inferred (XIST + 7 Y-genes) | ~200 |
| `validation_japanese_GSE167523.py` | Japanese MASLD | RNA-seq (counts matrix) | Annotated (M/F) | ~100 |
| `validation_japanese_hcc_GSE193066.py` | Japanese MASLD-HCC | RNA-seq (GCT) | Annotated (male/female) | 164 |
| `validation_statin_GSE130991.py` | Statin treatment | Affy HuGene 2.0 ST | Annotated (gender) | ~200 |
| `validation_bariatric_GSE83452.py` | Bariatric surgery | Affy HuGene 2.0 ST | Annotated (gender) | ~100 |
| `validation_tcga_lihc.py` | TCGA-LIHC (HCC) | RNA-seq (RSEM) | Annotated (TCGA-CDR) | ~400 |
| `validation_utils.py` | *(shared Python helpers)* | - | - | - |

## Shared pipeline

All scripts follow the same logical flow (implemented in `validation_utils.py`):

```
series matrix / counts file
    |
  dataset-specific preprocessing
  (probe -> gene, version strip, GCT parse, etc.)
    |
  gene-level matrix (ASCL1 accessible)
    |
  sex annotation (look up OR KMeans-infer from XIST + Y-genes)
    |
  ASCL1 extraction (dataset-specific transform)
    |
  Mann-Whitney U test (Male vs Female)
    |
  bar + jitter + significance bracket  ->  verify_<GSEID>_ASCL1_barplot.{png,pdf}
```

## Dataset-specific preprocessing at a glance

| Cohort | Key quirk |
|---|---|
| GTEx v10 | Bulk RNA-seq GCT; Hardy scale filter (DTHHRDY in {0,1,2,3}); sample IDs trimmed to SUBJID; **DESeq2 size-factor normalization + log2(x+1)**; R/DESeq2 pipeline (not Python) |
| GSE33814 | Illumina BeadArray; multi-probe collapse by **MAX**; sex inferred with single Y-gene |
| GSE135251 | Per-sample Ensembl-ID TSVs; version numbers stripped; CPM + log2 |
| GSE167523 | Single counts matrix with gene symbols; sex pre-annotated M/F |
| GSE193066 | **GCT format**; sample matching via **Title** (not GSM); pseudocount **0.1** |
| GSE130991 | GPL20265 SOFT file; gene symbol from `gene_assignment` field (// delimited); RMA direct |
| GSE83452 | Pre-computed probe->gene CSV; longitudinal (baseline + follow-up) |
| TCGA-LIHC | `GENE|ENTREZ` index; log2(RSEM+1); **barcode -01 = tumor / -11 = normal**; TCGA-CDR clinical |

## Quick-start

```bash
# All scripts look for data in <CohortID>_raw/ by default.
# Each script has a CONFIGURATION block at the top - update paths there.

Rscript validation_gtex_liver.R
python validation_german_GSE33814.py
python validation_european_GSE135251.py
python validation_japanese_GSE167523.py
python validation_japanese_hcc_GSE193066.py
python validation_statin_GSE130991.py
python validation_bariatric_GSE83452.py
python validation_tcga_lihc.py
```

Each script prints a step-by-step log and writes:
- `verify_<cohort>_ASCL1_barplot.png`
- `verify_<cohort>_ASCL1_barplot.pdf`

## Where to get the raw data

| Cohort | Accession | URL |
|---|---|---|
| GTEx v10 | GTEx v10 liver | https://gtexportal.org/home/downloads/adult-gtex/bulk_tissue_expression |
| GSE33814 | GSE33814 | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE33814 |
| GSE135251 | GSE135251 | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE135251 |
| GSE167523 | GSE167523 | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE167523 |
| GSE193066 | GSE193066 | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE193066 |
| GSE130991 | GSE130991 | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE130991 |
| GSE83452 | GSE83452 | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE83452 |
| TCGA-LIHC | TCGA Pan-Cancer | https://gdc.cancer.gov/about-data/publications/pancanatlas |

## Dependencies

### Python (all validation_*.py scripts)

```bash
pip install numpy pandas scipy scikit-learn matplotlib openpyxl
```

### R (validation_gtex_liver.R only)

```r
BiocManager::install(c("DESeq2"))
install.packages(c("ggplot2", "dplyr"))
```
