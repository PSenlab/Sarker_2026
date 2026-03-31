#!/usr/bin/env Rscript
## =========================================================
## GWAS CATALOG TO MOUSE ORTHOLOG MAPPING
## 
## This script:
## 1. Loads GWAS catalog (via Bioconductor gwascat)
## 2. Extracts genes associated with traits
## 3. Maps human genes to mouse orthologs
## 4. Creates annotation file for compartment bins
##
## Output: gwas_mouse_orthologs.tsv (for database_annotations.R)
## =========================================================

suppressPackageStartupMessages({
  library(data.table)
  library(biomaRt)
  library(GenomicRanges)
  library(gwascat)
})

message("\n")
message("╔══════════════════════════════════════════════════════════════════╗")
message("║     GWAS CATALOG → MOUSE ORTHOLOG MAPPING                        ║")
message("╚══════════════════════════════════════════════════════════════════╝\n")

## =========================================================
## CONFIGURATION
## =========================================================

outdir <- "/data/sarkern2/multiome_liver/Seurat/epigenome/external_annotations"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# Traits of interest for liver/aging research (optional filtering)
# Set to NULL to keep all traits
liver_aging_keywords <- c(
  # Liver-related
  "liver", "hepat", "biliary", "bile", "cholesterol", "triglyceride",
  "lipid", "fatty", "NAFLD", "NASH", "cirrhosis", "fibrosis",
  "ALT", "AST", "GGT", "albumin", "bilirubin",
  # Metabolic
  "metabol", "insulin", "glucose", "diabetes", "obesity", "BMI",
  "adipos", "weight", "glycem",
  # Aging-related
  "aging", "ageing", "longevity", "lifespan", "senescen",
  "telomere", "frailty", "sarcopenia",
  # Inflammation
  "inflam", "CRP", "IL-6", "TNF", "cytokine",
  # General
  "blood", "serum", "plasma"
)


## =========================================================
## STEP 1: LOAD GWAS CATALOG
## =========================================================

message("[STEP 1] Loading GWAS catalog...")

gwas_file <- file.path(outdir, "gwas_catalog_v1.0.2-associations.tsv")

if (!file.exists(gwas_file) || file.size(gwas_file) < 1000000) {

  # Try live download first, then fall back to bundled data
  gwas_cat <- tryCatch({
    message("  Trying makeCurrentGwascat()...")
    makeCurrentGwascat()
  }, error = function(e) {
    message(sprintf("  Live download failed: %s", e$message))
    message("  Checking bundled datasets...")
    for (ds in c("ebicat38", "ebicat_2020_04_15", "ebicat37")) {
      result <- tryCatch({
        data(list = ds, package = "gwascat", envir = environment())
        get(ds, envir = environment())
      }, warning = function(w) NULL, error = function(e) NULL)
      if (!is.null(result)) {
        message(sprintf("  Found bundled dataset: %s", ds))
        return(result)
      }
    }
    stop("No GWAS catalog data available. Run: data(package='gwascat') to check available datasets")
  })

  gwas_df <- as.data.frame(mcols(gwas_cat))
  gwas_df$CHR_ID <- as.character(seqnames(gwas_cat))
  gwas_df$CHR_POS <- start(gwas_cat)

  fwrite(as.data.table(gwas_df), gwas_file, sep = "\t")
  message(sprintf("  ✅ Saved GWAS catalog (%d associations)", nrow(gwas_df)))
  message("  ⚠️  For publication, replace with full catalog from:")
  message("     https://www.ebi.ac.uk/gwas/docs/file-downloads")

} else {
  file_size <- file.size(gwas_file) / 1e6
  message(sprintf("  ✅ GWAS catalog already exists (%.1f MB)", file_size))
}


## =========================================================
## STEP 2: LOAD AND PARSE GWAS DATA
## =========================================================

message("\n[STEP 2] Loading GWAS catalog...")

# Read header separately to get correct column names
header_line <- readLines(gwas_file, n = 1)
header_cols <- strsplit(header_line, "\t")[[1]]

message(sprintf("  Header has %d columns", length(header_cols)))

# Read the data skipping header, then assign column names
gwas <- fread(
  gwas_file, 
  header = FALSE,
  skip = 1,
  sep = "\t",
  quote = "",
  fill = TRUE
)

# Assign proper column names (handle if data has more/fewer cols than header)
n_cols <- min(ncol(gwas), length(header_cols))
setnames(gwas, 1:n_cols, header_cols[1:n_cols])

message(sprintf("  ✅ Loaded %d associations", nrow(gwas)))

# Now the columns should be correctly named
gwas_cols <- names(gwas)
message("  Key columns found:")
message(sprintf("     DISEASE/TRAIT: %s", "DISEASE/TRAIT" %in% gwas_cols))
message(sprintf("     MAPPED_GENE: %s", "MAPPED_GENE" %in% gwas_cols))
message(sprintf("     P-VALUE: %s", "P-VALUE" %in% gwas_cols))
message(sprintf("     SNPS: %s", "SNPS" %in% gwas_cols))
message(sprintf("     PUBMEDID: %s", "PUBMEDID" %in% gwas_cols))

# Extract relevant columns using exact names from the file
gwas_clean <- gwas[, .(
  trait = `DISEASE/TRAIT`,
  mapped_gene = MAPPED_GENE,
  reported_gene = `REPORTED GENE(S)`,
  pvalue = as.numeric(`P-VALUE`),
  snp = SNPS,
  pubmed = as.character(PUBMEDID),
  study = STUDY,
  chr = CHR_ID,
  pos = CHR_POS
)]

# Remove rows without gene information
gwas_clean <- gwas_clean[!is.na(mapped_gene) & mapped_gene != "" & mapped_gene != "NR" & mapped_gene != "-"]

message(sprintf("\n  ✅ %d associations with mapped genes", nrow(gwas_clean)))


## =========================================================
## STEP 3: FILTER FOR RELEVANT TRAITS (OPTIONAL)
## =========================================================

message("\n[STEP 3] Filtering for relevant traits...")

if (!is.null(liver_aging_keywords) && length(liver_aging_keywords) > 0) {
  
  # Create regex pattern
  pattern <- paste(liver_aging_keywords, collapse = "|")
  
  # Filter traits
  gwas_filtered <- gwas_clean[grepl(pattern, trait, ignore.case = TRUE)]
  
  message(sprintf("  ✅ Filtered to %d associations matching keywords", nrow(gwas_filtered)))
  message(sprintf("     (from %d total)", nrow(gwas_clean)))
  
  # Show top traits
  top_traits <- gwas_filtered[, .N, by = trait][order(-N)][1:20]
  message("\n  Top 20 traits:")
  for (i in 1:nrow(top_traits)) {
    message(sprintf("     %d. %s (n=%d)", i, top_traits$trait[i], top_traits$N[i]))
  }
  
} else {
  gwas_filtered <- gwas_clean
  message("  ⚠️ No filtering applied - using all traits")
}
