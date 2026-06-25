# Ascl1 zone-stratified MAST DE (female hepatocytes)

CDR-controlled differential expression of Ascl1-positive vs Ascl1-negative
female hepatocytes, run globally and per hepatic zone. MAST includes the
cellular detection rate (`cngeneson`) as a covariate to control for the
higher sequencing depth of Ascl1-detected cells.

## Pipeline

1. `01_export_for_mast.py` — load adata, subset to female hepatocytes, map
   sub-clusters (Hep-01..07) to zones, binarize Ascl1 (detection > 0), and
   export log-normalized expression + metadata per subset to `mast_export/`.
2. `02_run_mast.R <subset>` — CDR-controlled MAST hurdle model for one subset.
   `subset` is one of `allfemhep`, `periportal`, `midlobular`, `pericentral`.
3. `run_all.sh` — export once, then loop MAST over all four subsets.

## Run

```bash
# set ADATA_PATH at the top of 01_export_for_mast.py first
python 01_export_for_mast.py

# serial MAST; global (~90k cells) needs the most memory
Rscript 02_run_mast.R allfemhep      # sinteractive --mem=128g
Rscript 02_run_mast.R periportal     # sinteractive --mem=96g
Rscript 02_run_mast.R midlobular
Rscript 02_run_mast.R pericentral
```

## Outputs (`ascl1_de_by_zone_MAST/`)

Per subset: `_all.csv`, `_significant.csv`, `_up.csv`, `_down.csv`.
Columns: `primerid, p_value, logFC, ci_lo, ci_hi, fdr`.
`logFC` is on the natural-log scale (log1p input); divide by `log(2)` for log2.
Global subset is written as `..._female_hepatocytes_*`.

## Notes

- `.X` must be the log-normalized layer (matches `use_raw=False`).
- Serial execution (`SerialParam`) avoids MAST fork-death / memory blowup on
  large cell counts; bump `MIN_FRAC` in the export to 0.05 if memory is tight.
- Downstream enrichment uses a strict gene filter (`fdr < 0.05 & |logFC| > 0.21`)
  before Enrichr/GSEA to remove residual detection-depth genes.
