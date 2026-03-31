#!/usr/bin/env Rscript
#===============================================================================
# GWAS-SCAVENGE Integration via Peak-to-Gene Links
#
# Identifies cell-type-specific enrichment of human GWAS trait-associated
# variants within single-cell ATAC-seq data using the SCAVENGE framework.
#
# Pipeline:
#   1. Load cell-type-specific ArchR project with Peak-to-Gene (P2G) links
#   2. Build SummarizedExperiment from P2G-linked peaks (binarized)
#   3. Filter GWAS Catalog for liver-relevant traits (inclusion/exclusion)
#   4. LiftOver GWAS SNPs from hg38 → mm10
#   5. LD-proxy expansion (±25 kb) and merge overlapping regions per trait
#   6. Overlap expanded GWAS regions with P2G peaks
#   7. Generate weighted BED files per passing trait
#   8. Run SCAVENGE per trait:
#        - GC bias correction → background peaks → weighted deviations
#        - Seed cell identification → KNN graph → random walk → TRS
#        - Permutation test (1000×) → empirical p-values → BH FDR
#   9. Save per-trait raw + filtered (FDR < 0.05) results
#  10. Combine into summary tables
#
# Usage:
#   Rscript run_gwas_scavenge.R \
#     --cell_type Hepatocyte \
#     --archr_project /path/to/ArchR_Project \
#     --output_dir /path/to/output \
#     --gwas_catalog /path/to/gwas_catalog_v1.0.2-associations.tsv \
#     --chain_file /path/to/hg38ToMm10.over.chain \
#     --threads 60 \
#     --scavenge_cores 55
#
# Required packages:
#   ArchR, SCAVENGE, chromVAR, gchromVAR, GenomicRanges, rtracklayer,
#   SummarizedExperiment, BSgenome.Mmusculus.UCSC.mm10, Matrix,
#   data.table, dplyr, tidyr, parallel
#
# Reference:
#   Yu et al. (2022) SCAVENGE: Single Cell Analysis of Variant Enrichment
#   through Network propagation of GEnetic associations. Nat Biotechnol.
#===============================================================================

suppressPackageStartupMessages({
  library(ArchR)
  library(GenomicRanges)
  library(dplyr)
  library(tidyr)
  library(Matrix)
  library(SummarizedExperiment)
  library(BSgenome.Mmusculus.UCSC.mm10)
  library(SCAVENGE)
  library(chromVAR)
  library(gchromVAR)
  library(parallel)
  library(data.table)
  library(rtracklayer)
})

#===============================================================================
# COMMAND-LINE ARGUMENTS
#===============================================================================
parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  
  # Defaults
  params <- list(
    cell_type       = NULL,
    archr_project   = NULL,
    output_dir      = NULL,
    gwas_catalog    = NULL,
    chain_file      = NULL,
    threads         = 60,
    scavenge_cores  = 55,
    peak_extension  = 500,
    ld_window       = 25000,
    min_regions     = 2,
    min_seed_cells  = 5,
    min_valid_cells = 10,
    fdr_threshold   = 0.05,
    permutations    = 1000,
    seed            = 42
  )
  
  i <- 1
  while (i <= length(args)) {
    flag <- args[i]
    val  <- if (i + 1 <= length(args)) args[i + 1] else NULL
    switch(flag,
      "--cell_type"       = { params$cell_type       <- val;           i <- i + 2 },
      "--archr_project"   = { params$archr_project   <- val;           i <- i + 2 },
      "--output_dir"      = { params$output_dir      <- val;           i <- i + 2 },
      "--gwas_catalog"    = { params$gwas_catalog    <- val;           i <- i + 2 },
      "--chain_file"      = { params$chain_file      <- val;           i <- i + 2 },
      "--threads"         = { params$threads         <- as.integer(val); i <- i + 2 },
      "--scavenge_cores"  = { params$scavenge_cores  <- as.integer(val); i <- i + 2 },
      "--peak_extension"  = { params$peak_extension  <- as.integer(val); i <- i + 2 },
      "--ld_window"       = { params$ld_window       <- as.integer(val); i <- i + 2 },
      "--min_regions"     = { params$min_regions     <- as.integer(val); i <- i + 2 },
      "--min_seed_cells"  = { params$min_seed_cells  <- as.integer(val); i <- i + 2 },
      "--min_valid_cells" = { params$min_valid_cells <- as.integer(val); i <- i + 2 },
      "--fdr_threshold"   = { params$fdr_threshold   <- as.numeric(val); i <- i + 2 },
      "--permutations"    = { params$permutations    <- as.integer(val); i <- i + 2 },
      "--seed"            = { params$seed            <- as.integer(val); i <- i + 2 },
      "--help"            = { print_usage(); q(save = "no") },
      { cat("Unknown argument:", flag, "\n"); print_usage(); q(save = "no", status = 1) }
    )
  }
  
  # Validate required arguments
  required <- c("cell_type", "archr_project", "output_dir", "gwas_catalog", "chain_file")
  missing <- required[sapply(required, function(x) is.null(params[[x]]))]
  if (length(missing) > 0) {
    cat("ERROR: Missing required arguments:", paste0("--", missing, collapse = ", "), "\n\n")
    print_usage()
    q(save = "no", status = 1)
  }
  
  return(params)
}

print_usage <- function() {
  cat("
Usage:
  Rscript run_gwas_scavenge.R \\
    --cell_type       <name>   Cell type label (e.g., Hepatocyte, Kupffer)
    --archr_project   <path>   Path to ArchR project with P2G links
    --output_dir      <path>   Output directory for this cell type
    --gwas_catalog    <path>   GWAS Catalog associations TSV file
    --chain_file      <path>   hg38ToMm10 liftOver chain file

  Optional:
    --threads         <int>    ArchR threads (default: 60)
    --scavenge_cores  <int>    SCAVENGE permutation cores (default: 55)
    --peak_extension  <int>    Peak extension for overlap, bp (default: 500)
    --ld_window       <int>    LD proxy window, bp (default: 25000)
    --min_regions     <int>    Min GWAS regions per trait (default: 2)
    --min_seed_cells  <int>    Min seed cells for SCAVENGE (default: 5)
    --min_valid_cells <int>    Min valid cells after random walk (default: 10)
    --fdr_threshold   <num>    FDR significance threshold (default: 0.05)
    --permutations    <int>    Permutation test iterations (default: 1000)
    --seed            <int>    Random seed (default: 42)
    --help                     Show this message
\n")
}

params <- parse_args()

#===============================================================================
# SETUP
#===============================================================================
CELL_TYPE      <- params$cell_type
ARCHR_PROJECT  <- params$archr_project
OUTPUT_DIR     <- params$output_dir
GWAS_FILE      <- params$gwas_catalog
CHAIN_FILE     <- params$chain_file

addArchRThreads(params$threads)
addArchRGenome("mm10")
set.seed(params$seed)

# Directory structure
gwas_dir       <- file.path(OUTPUT_DIR, "gwas_p2g/Liver")
checkpoint_dir <- file.path(gwas_dir, "checkpoints")
bed_dir        <- file.path(gwas_dir, "trait_beds")
results_dir    <- file.path(gwas_dir, "results")
raw_dir        <- file.path(results_dir, "raw")
filtered_dir   <- file.path(results_dir, "filtered")

for (d in c(gwas_dir, checkpoint_dir, bed_dir, results_dir, raw_dir, filtered_dir)) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

log_file <- file.path(results_dir,
                      paste0("log_", CELL_TYPE, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt"))

#===============================================================================
# HELPER FUNCTIONS
#===============================================================================
write_log <- function(msg) {
  cat(msg, "\n")
  cat(paste0("[", Sys.time(), "] ", msg, "\n"), file = log_file, append = TRUE)
}

save_checkpoint <- function(obj, name) {
  saveRDS(obj, file.path(checkpoint_dir, paste0(name, ".rds")))
  write_log(paste("  Checkpoint saved:", name))
}

load_checkpoint <- function(name) {
  fp <- file.path(checkpoint_dir, paste0(name, ".rds"))
  if (file.exists(fp)) { write_log(paste("  Checkpoint loaded:", name)); return(readRDS(fp)) }
  return(NULL)
}

checkpoint_exists <- function(name) {
  file.exists(file.path(checkpoint_dir, paste0(name, ".rds")))
}

clean_trait_name <- function(trait) {
  clean <- gsub("[^a-zA-Z0-9]", "_", trait)
  clean <- gsub("_+", "_", clean)
  clean <- gsub("^_|_$", "", clean)
  substr(clean, 1, 100)
}

#===============================================================================
# LIVER TRAIT INCLUSION PATTERNS
#
# Curated regex patterns to identify liver-relevant GWAS traits from the
# NHGRI-EBI GWAS Catalog DISEASE/TRAIT column. Covers: core liver terms,
# liver diseases, liver enzymes, liver-synthesized proteins, lipoproteins,
# iron metabolism, viral hepatitis, autoimmune liver disease, drug-induced
# liver injury, and biliary pathology.
#===============================================================================
LIVER_PATTERNS <- c(
  # Core liver
  "\\bliver\\b", "hepat",
  # Liver diseases
  "NAFLD", "NASH", "steatosis", "cirrhosis",
  # Liver enzymes
  "\\bALT\\b", "\\bAST\\b", "\\bGGT\\b", "\\bALP\\b",
  "alanine aminotransferase", "aspartate aminotransferase",
  "gamma glutamyl transferase", "gamma glutamyltransferase",
  "gamma-glutamyl transpeptidase", "alkaline phosphatase", "bilirubin",
  # Liver-synthesized proteins
  "albumin", "fibrinogen", "prothrombin", "coagulation factor",
  "complement C", "haptoglobin", "transferrin", "ceruloplasmin",
  "alpha.1.antitrypsin", "angiotensinogen",
  # Lipoproteins (liver-produced)
  "apolipoprotein", "\\bAPO[ABCDEFHLM]\\b", "\\bVLDL\\b",
  # Iron metabolism
  "ferritin", "hemochromatosis", "hepcidin", "\\bliver iron\\b",
  # Viral hepatitis
  "hepatitis B", "hepatitis C", "\\bHBV\\b", "\\bHCV\\b",
  # Autoimmune liver
  "autoimmune hepatitis", "primary biliary", "primary sclerosing",
  "\\bPBC\\b", "\\bPSC\\b",
  # Drug-induced liver injury
  "drug.induced liver", "hepatotoxicity",
  # Biliary
  "cholangiocarcinoma", "biliary", "bile duct", "gallbladder"
)

#===============================================================================
# EXCLUSION PATTERNS
#
# Remove non-liver traits that match inclusion patterns due to shared
# keywords (e.g., brain iron MRI studies matching "iron", urinary albumin
# matching "albumin", gamma-glutamyl amino acid metabolites matching "GGT").
#===============================================================================
EXCLUSION_PATTERNS <- c(
  # Brain iron MRI
  "putamen iron", "pallidum iron", "caudate iron", "substantia nigra iron",
  "thalamus iron", "hippocampus iron", "accumbens iron", "amygdala iron",
  "\\bR2\\* MRI", "quantitative susceptibility mapping",
  # Environmental
  "environmental pollutants", "environmental tobacco", "environmental stress",
  "environmental factor", "environmentalism", "environmental exposures",
  "occupational environmental",
  # Gamma-glutamyl amino acids (metabolites, not GGT enzyme)
  "gamma-glutamylglutamate", "gamma-glutamylglutamine", "gamma-glutamylleucine",
  "gamma-glutamylvaline", "gamma-glutamyltyrosine", "gamma-glutamylphenylalanine",
  "gamma-glutamylthreonine", "gamma-glutamylisoleucine", "gamma-glutamylmethionine",
  "gamma-glutamylhistidine", "gamma-glutamylalanine", "gamma-glutamyltryptophan",
  "gamma-glutamylglycine", "gamma-glutamyl-alpha-lysine", "gamma-glutamyl-epsilon-lysine",
  "gamma-glutamyl-2-aminobutyrate", "gamma-glutamylcitrulline",
  "gamma-glutamylaminecyclotransferase", "gamma-glutamylcyclotransferase",
  "gamma-glutamyl hydrolase", "Glutathione-specific gamma-glutamylcyclotransferase",
  "Protein-glutamine gamma-glutamyltransferase", "Inactive gamma-glutamyltranspeptidase",
  "gamma-glutamyltransferase 5",
  # Redundant VLDL subfractions
  "very large VLDL", "very small VLDL", "large VLDL", "medium VLDL", "small VLDL",
  "very very large VLDL", "chylomicrons and extremely large VLDL",
  "chylomicron and extremely large VLDL", "Fasting.*VLDL", "Postprandial.*VLDL",
  "VLDL particle concentration", "VLDL particles", "diameter.*VLDL",
  "to Total Lipids.*VLDL percentage", "to total lipids ratio.*VLDL", "UKB data field 23",
  # APOE / Alzheimer's
  "Alzheimer", "APOE e4", "APOEe4", "APOE E4", "dementia", "amyloid",
  "cerebral", "cerebrospinal",
  # Apoptosis proteins
  "apoptosis", "apoptotic", "APOBEC",
  # Cystic / pulmonary fibrosis
  "cystic fibrosis", "pulmonary fibrosis", "idiopathic pulmonary",
  "lung disease", "lung function", "FEV1", "FVC", "meconium ileus", "CFTR",
  # Intestinal / placental ALP
  "intestinal.type alkaline phosphatase", "placental.type alkaline phosphatase",
  "placental alkaline phosphatase", "alkaline phosphatase, placental",
  "alkaline phosphatase, tissue-nonspecific",
  # Urinary albumin
  "urinary albumin", "microalbuminuria", "macroalbuminuria", "albuminuria",
  "albumin-to-creatinine ratio", "albumin excretion", "microalbumin",
  # CKD / kidney
  "Chronic kidney disease", "chronic kidney disease", "renal", "kidney disease",
  "end.stage.*kidney", "end stage.*kidney",
  # Other non-liver
  "interstitial cystitis", "multiple sclerosis", "rheumatoid arthritis",
  "inflammatory bowel disease", "Crohn's disease", "ulcerative colitis",
  "ankylosing spondylitis", "psoriasis", "colorectal cancer", "colorectal carcinoma",
  "acute lymphoblastic leukemia", "B-cell lymphoblastic leukemia",
  # Misc
  "elite athletes", "mercapturic acid", "vaginal microbiome", "urea cycle",
  "Biological.*Liver Condition", "Liver Liking", "pancreas iron", "spleen iron"
)

#===============================================================================
# START
#===============================================================================
write_log("================================================================")
write_log(paste("GWAS-SCAVENGE:", CELL_TYPE))
write_log("================================================================")
write_log(paste("ArchR project: ", ARCHR_PROJECT))
write_log(paste("Output:        ", OUTPUT_DIR))
write_log(paste("GWAS catalog:  ", GWAS_FILE))
write_log(paste("Chain file:    ", CHAIN_FILE))
write_log(paste("Threads:       ", params$threads))
write_log(paste("SCAVENGE cores:", params$scavenge_cores))
write_log(paste("LD window:     ", params$ld_window, "bp"))
write_log(paste("Peak extension:", params$peak_extension, "bp"))
write_log(paste("FDR threshold: ", params$fdr_threshold))
write_log(paste("Permutations:  ", params$permutations))

#===============================================================================
# STEP 1: Load ArchR Project
#===============================================================================
write_log("\n====== Step 1: Load ArchR Project ======")

proj <- tryCatch(
  loadArchRProject(ARCHR_PROJECT),
  error = function(e) { write_log(paste("FATAL:", e$message)); stop(e) }
)
write_log(paste("Cells:", length(getCellNames(proj))))

#===============================================================================
# STEP 2: Build SummarizedExperiment from P2G Peaks
#===============================================================================
write_log("\n====== Step 2: Build SE_Data from P2G Peaks ======")

if (checkpoint_exists("step2_SE_Data") && checkpoint_exists("step2_peaks_p2g")) {
  SE_Data  <- load_checkpoint("step2_SE_Data")
  peaks_p2g <- load_checkpoint("step2_peaks_p2g")
  write_log(paste("  Dimensions:", nrow(SE_Data), "peaks ×", ncol(SE_Data), "cells"))
} else {
  
  # Extract P2G links
  p2g <- getPeak2GeneLinks(ArchRProj = proj, returnLoops = FALSE)
  p2g_peak_idx <- unique(p2g$idxATAC)
  p2g_gene_idx <- unique(p2g$idxRNA)
  
  write_log(paste("  P2G:", nrow(p2g), "links |", length(p2g_peak_idx), "peaks |",
                  length(p2g_gene_idx), "genes"))
  
  # Build binarized peak matrix subset
  peakMat <- getMatrixFromProject(ArchRProj = proj, useMatrix = "PeakMatrix", binarize = TRUE)
  write_log(paste("  Full PeakMatrix:", nrow(peakMat), "×", ncol(peakMat)))
  
  SE_Data <- peakMat[p2g_peak_idx, ]
  if (!"counts" %in% assayNames(SE_Data)) assayNames(SE_Data)[1] <- "counts"
  peaks_p2g <- rowRanges(SE_Data)
  
  save_checkpoint(SE_Data, "step2_SE_Data")
  save_checkpoint(peaks_p2g, "step2_peaks_p2g")
  
  write_log(paste("  SE_Data:", nrow(SE_Data), "peaks ×", ncol(SE_Data), "cells"))
  write_log(paste("  P2G peak fraction:", round(nrow(SE_Data) / nrow(peakMat) * 100, 1), "%"))
  rm(peakMat); gc(verbose = FALSE)
}

#===============================================================================
# STEP 3: Load GWAS Catalog & Filter Liver Traits
#===============================================================================
write_log("\n====== Step 3: Filter GWAS Catalog for Liver Traits ======")

if (checkpoint_exists("step3_gwas_filtered") && checkpoint_exists("step3_liver_traits")) {
  gwas_filtered <- load_checkpoint("step3_gwas_filtered")
  liver_traits  <- load_checkpoint("step3_liver_traits")
} else {
  
  if (!file.exists(GWAS_FILE)) stop("GWAS Catalog not found: ", GWAS_FILE)
  
  gwas <- fread(GWAS_FILE, quote = "")
  write_log(paste("  GWAS Catalog:", nrow(gwas), "associations |",
                  length(unique(gwas$`DISEASE/TRAIT`)), "unique traits"))
  
  unique_traits <- unique(gwas$`DISEASE/TRAIT`)
  
  # Inclusion
  liver_traits_raw <- unique(unlist(lapply(LIVER_PATTERNS, function(p) {
    grep(p, unique_traits, value = TRUE, ignore.case = TRUE)
  })))
  
  # Exclusion
  exclusion_regex <- paste(EXCLUSION_PATTERNS, collapse = "|")
  liver_traits <- liver_traits_raw[!grepl(exclusion_regex, liver_traits_raw, ignore.case = TRUE)]
  
  write_log(paste("  Matched:", length(liver_traits_raw), "→ After exclusion:", length(liver_traits)))
  
  # Log trait list
  for (t in sort(liver_traits)) {
    write_log(sprintf("    [%4d] %s", nrow(gwas[`DISEASE/TRAIT` == t]), t))
  }
  
  gwas_filtered <- gwas[`DISEASE/TRAIT` %in% liver_traits]
  
  save_checkpoint(gwas_filtered, "step3_gwas_filtered")
  save_checkpoint(liver_traits, "step3_liver_traits")
  rm(gwas); gc(verbose = FALSE)
}

write_log(paste("  GWAS entries for liver traits:", nrow(gwas_filtered)))

#===============================================================================
# STEP 4: LiftOver hg38 → mm10
#===============================================================================
write_log("\n====== Step 4: LiftOver hg38 → mm10 ======")

if (checkpoint_exists("step4_gwas_mm10")) {
  gwas_mm10 <- load_checkpoint("step4_gwas_mm10")
} else {
  
  gwas_clean <- gwas_filtered[!is.na(CHR_ID) & CHR_ID != "" & !is.na(CHR_POS) & CHR_POS != ""]
  gwas_clean <- gwas_clean[!grepl(";|x", CHR_ID, ignore.case = TRUE)]
  gwas_clean[, CHR_POS := as.numeric(CHR_POS)]
  gwas_clean <- gwas_clean[!is.na(CHR_POS)]
  
  gwas_gr <- GRanges(
    seqnames = paste0("chr", gwas_clean$CHR_ID),
    ranges   = IRanges(start = gwas_clean$CHR_POS, width = 1),
    snp      = gwas_clean$SNPS,
    trait    = gwas_clean$`DISEASE/TRAIT`,
    pvalue   = gwas_clean$`P-VALUE`,
    gene     = gwas_clean$MAPPED_GENE
  )
  gwas_gr <- keepStandardChromosomes(gwas_gr, pruning.mode = "coarse")
  
  # Handle chain file (support .chain or .chain.gz)
  chain_path <- CHAIN_FILE
  if (grepl("\\.gz$", chain_path)) {
    unzipped <- sub("\\.gz$", "", chain_path)
    if (!file.exists(unzipped)) system(paste("gunzip -k", chain_path))
    chain_path <- unzipped
  }
  
  chain <- import.chain(chain_path)
  gwas_mm10 <- unlist(liftOver(gwas_gr, chain))
  
  write_log(paste("  Mapped:", length(gwas_mm10), "/", length(gwas_gr),
                  "(", round(length(gwas_mm10) / length(gwas_gr) * 100, 1), "%)"))
  save_checkpoint(gwas_mm10, "step4_gwas_mm10")
}

#===============================================================================
# STEP 5: LD-Proxy Expansion
#===============================================================================
write_log("\n====== Step 5: LD-Proxy Expansion ======")

if (checkpoint_exists("step5_gwas_mm10_ld")) {
  gwas_mm10_ld <- load_checkpoint("step5_gwas_mm10_ld")
} else {
  
  gwas_expanded <- gwas_mm10
  start(gwas_expanded) <- pmax(1, start(gwas_mm10) - params$ld_window)
  end(gwas_expanded)   <- end(gwas_mm10) + params$ld_window
  
  # Merge overlapping regions per trait
  traits <- unique(gwas_expanded$trait)
  merged <- vector("list", length(traits))
  
  for (i in seq_along(traits)) {
    tr_gr <- gwas_expanded[gwas_expanded$trait == traits[i]]
    tr_red <- reduce(tr_gr)
    if (length(tr_red) > 0) {
      mcols(tr_red)$trait <- traits[i]
      mcols(tr_red)$n_snps_merged <- countOverlaps(tr_red, tr_gr)
      merged[[i]] <- tr_red
    }
  }
  
  merged <- merged[!sapply(merged, is.null)]
  gwas_mm10_ld <- unlist(GRangesList(merged), use.names = FALSE)
  
  write_log(paste("  LD-expanded regions:", length(gwas_mm10_ld)))
  save_checkpoint(gwas_mm10_ld, "step5_gwas_mm10_ld")
}

#===============================================================================
# STEP 6: Overlap GWAS Regions with P2G Peaks
#===============================================================================
write_log("\n====== Step 6: GWAS–Peak Overlap ======")

if (checkpoint_exists("step6_trait_counts")) {
  trait_counts <- load_checkpoint("step6_trait_counts")
} else {
  
  peaks_ext <- peaks_p2g
  start(peaks_ext) <- pmax(1, start(peaks_p2g) - params$peak_extension)
  end(peaks_ext)   <- end(peaks_p2g) + params$peak_extension
  
  ov <- findOverlaps(peaks_ext, gwas_mm10_ld, ignore.strand = TRUE)
  
  if (length(ov) > 0) {
    ov_dt <- data.table(
      peak_idx = queryHits(ov),
      gwas_idx = subjectHits(ov),
      trait    = gwas_mm10_ld$trait[subjectHits(ov)]
    )
    trait_counts <- ov_dt[, .(n_regions = uniqueN(gwas_idx), n_peaks = uniqueN(peak_idx)), by = trait]
    trait_counts[, passes_filter := n_regions >= params$min_regions]
  } else {
    trait_counts <- data.table(trait = character(0), n_regions = integer(0),
                               n_peaks = integer(0), passes_filter = logical(0))
  }
  
  save_checkpoint(trait_counts, "step6_trait_counts")
}

passing_traits <- trait_counts[passes_filter == TRUE, trait]
write_log(paste("  Traits passing (≥", params$min_regions, "regions):", length(passing_traits)))

#===============================================================================
# STEP 7: Create BED Files
#===============================================================================
write_log("\n====== Step 7: Create BED Files ======")

if (checkpoint_exists("step7_bed_files")) {
  trait_bed_files <- load_checkpoint("step7_bed_files")
  write_log(paste("  BED files:", length(trait_bed_files)))
} else {
  
  create_bed <- function(trait_name, gwas_gr, out_dir) {
    idx <- which(gwas_gr$trait == trait_name)
    if (length(idx) == 0) return(NULL)
    
    gr <- gwas_gr[idx]
    raw_score <- if (!is.null(mcols(gr)$n_snps_merged)) mcols(gr)$n_snps_merged else rep(1, length(gr))
    score <- if (length(unique(raw_score)) > 1 && max(raw_score) > 1) log1p(raw_score) else raw_score
    
    bed <- data.frame(
      chr = as.character(seqnames(gr)), start = start(gr), end = end(gr),
      name = paste0(clean_trait_name(trait_name), "_", seq_along(gr)),
      score = score, strand = ".", thickStart = start(gr), thickEnd = end(gr)
    )
    
    fp <- file.path(out_dir, paste0(clean_trait_name(trait_name), ".bed"))
    write.table(bed, fp, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
    return(fp)
  }
  
  trait_bed_files <- list()
  for (i in seq_along(passing_traits)) {
    bf <- tryCatch(create_bed(passing_traits[i], gwas_mm10_ld, bed_dir), error = function(e) NULL)
    if (!is.null(bf)) trait_bed_files[[passing_traits[i]]] <- bf
    if (i %% 50 == 0) write_log(paste("  Created", i, "/", length(passing_traits)))
  }
  
  write_log(paste("  Total BED files:", length(trait_bed_files)))
  save_checkpoint(trait_bed_files, "step7_bed_files")
}

#===============================================================================
# STEP 8: SCAVENGE Function
#
# For each trait:
#   1. GC bias correction on P2G peak SummarizedExperiment
#   2. Compute background peaks (200 iterations)
#   3. Build weighted trait-import vector from BED-peak overlaps
#   4. chromVAR weighted deviations → per-cell z-scores
#   5. Seed cell identification (top 5% z-score)
#   6. TF-IDF → LSI (30 dims) → mutual KNN (k=30) graph
#   7. Random walk with restart (gamma=0.05)
#   8. Trait Relevance Score (TRS) = capped + min-max-scaled network scores
#   9. Permutation test (1000×) → empirical p-values → BH FDR
#===============================================================================

run_scavenge <- function(SE_Data, trait_file, mycores, min_seed, min_cells, n_perm) {
  
  result <- tryCatch({
    
    SE_sub <- addGCBias(SE_Data, genome = BSgenome.Mmusculus.UCSC.mm10)
    SE_bg  <- getBackgroundPeaks(SE_sub, niterations = 200)
    
    if (!is.matrix(SE_bg) || nrow(SE_bg) != nrow(SE_sub))
      return(list(success = FALSE, error = "Background peaks failed"))
    
    # Load BED and find overlaps
    bed_df <- read.table(trait_file, sep = "\t", stringsAsFactors = FALSE)
    if (nrow(bed_df) == 0) return(list(success = FALSE, error = "Empty BED"))
    colnames(bed_df) <- c("chr","start","end","name","score","strand","thickStart","thickEnd")
    
    bed_gr <- GRanges(seqnames = bed_df$chr,
                      ranges = IRanges(start = bed_df$start, end = bed_df$end),
                      score = bed_df$score)
    peaks <- rowRanges(SE_sub)
    ov <- findOverlaps(peaks, bed_gr)
    if (length(ov) == 0) return(list(success = FALSE, error = "No overlaps"))
    
    # Weighted trait-import vector
    trait_import <- rep(0, length(peaks))
    scores_by_peak <- tapply(bed_gr$score[subjectHits(ov)], queryHits(ov), sum)
    trait_import[as.integer(names(scores_by_peak))] <- scores_by_peak
    
    n_peaks_signal <- sum(trait_import > 0)
    if (n_peaks_signal < 1) return(list(success = FALSE, error = "No peaks with signal"))
    
    trait_import <- Matrix(trait_import, ncol = 1, sparse = TRUE)
    colnames(trait_import) <- "trait"
    
    # Weighted deviations → z-scores
    SE_DEV <- computeWeightedDeviations(SE_sub, trait_import, background_peaks = SE_bg)
    z <- assays(SE_DEV)[["z"]]
    
    if (sum(is.na(z)) / length(z) > 0.9)
      return(list(success = FALSE, error = paste("NA rate:", round(sum(is.na(z))/length(z)*100,1), "%")))
    
    z_mat <- data.frame(cell = colnames(z), z_score = as.numeric(z),
                        colData(SE_sub)[colnames(z), , drop = FALSE],
                        row.names = colnames(z))
    z_mat <- z_mat[complete.cases(z_mat), ]
    
    if (nrow(z_mat) < min_cells)
      return(list(success = FALSE, error = paste("Only", nrow(z_mat), "cells")))
    
    # Seed cells and scale factor
    seed_idx <- seedindex(z_mat$z_score, 0.05)
    if (sum(seed_idx) < min_seed)
      return(list(success = FALSE, error = paste("Only", sum(seed_idx), "seeds")))
    
    scale_factor <- cal_scalefactor(z_score = z_mat$z_score, 0.01)
    
    # KNN graph from P2G peak accessibility
    peak_cell_mat <- assay(SE_sub)[, rownames(z_mat)]
    tfidf_mat <- tfidf(bmat = peak_cell_mat, mat_binary = TRUE, TF = TRUE, log_TF = TRUE)
    lsi_mat   <- do_lsi(tfidf_mat, dims = 30)
    knn30     <- getmutualknn(lsi_mat, 30)
    
    seed_cells <- intersect(rownames(z_mat)[seed_idx], rownames(knn30))
    if (length(seed_cells) < min_seed)
      return(list(success = FALSE, error = paste("Only", length(seed_cells), "valid seeds")))
    
    # Random walk
    np_score <- randomWalk_sparse(intM = knn30, queryCells = seed_cells, gamma = 0.05)
    active <- np_score != 0
    if (sum(active) < min_cells)
      return(list(success = FALSE, error = paste("Only", sum(active), "after walk")))
    
    knn30    <- knn30[active, active]
    np_score <- np_score[active]
    TRS <- np_score %>% capOutlierQuantile(., 0.95) %>% max_min_scale
    TRS <- TRS * scale_factor
    
    med_TRS <- median(TRS, na.rm = TRUE)
    if (med_TRS < 0.001)
      return(list(success = FALSE, error = paste("Walk collapse, median TRS:", round(med_TRS, 6))))
    
    filt_cells <- names(np_score)
    mono_mat <- data.frame(z_mat[filt_cells, ],
                           seed_idx = filt_cells %in% seed_cells,
                           np_score = np_score,
                           TRS = TRS)
    
    # Permutation test
    perm_test <- function(knn, seed_idx, np_scores, n_perm, mycores, rw_gamma = 0.05) {
      cell_mat <- data.frame(cell = 1:nrow(knn), degree = colSums(knn))
      seed_mat <- data.frame(seed = which(seed_idx),
                             degree = colSums(knn[, which(seed_idx), drop = FALSE]))
      seed_tbl <- data.frame(table(seed_mat$degree))
      
      if (nrow(seed_tbl) == 0)
        return(data.frame(seed_idx = seed_idx, true_cell_top_idx = FALSE, empirical_pval = 1))
      
      degree_groups <- tapply(cell_mat[,1], cell_mat[,2], list)
      matched <- degree_groups[names(degree_groups) %in% seed_tbl$Var1]
      if (length(matched) == 0)
        return(data.frame(seed_idx = seed_idx, true_cell_top_idx = FALSE, empirical_pval = 1))
      
      perm_scores <- mclapply(1:n_perm, mc.cores = mycores, function(i) {
        tryCatch({
          sampled <- matched %>% mapply(sample, ., seed_tbl$Freq) %>% unlist() %>% sort()
          randomWalk_sparse(intM = knn, queryCells = rownames(knn)[as.numeric(sampled)], gamma = rw_gamma)
        }, error = function(e) rep(NA, nrow(knn)))
      })
      
      perm_df <- data.frame(sapply(perm_scores, c))
      valid <- colSums(is.na(perm_df)) == 0
      if (sum(valid) == 0)
        return(data.frame(seed_idx = seed_idx, true_cell_top_idx = FALSE, empirical_pval = 1))
      
      perm_df <- perm_df[, valid, drop = FALSE]
      gt_mat <- apply(perm_df, 2, function(x) x > np_scores)
      if (is.null(dim(gt_mat))) gt_mat <- matrix(gt_mat, ncol = 1)
      
      emp_pval <- rowSums(gt_mat) / ncol(gt_mat)
      data.frame(seed_idx, true_cell_top_idx = emp_pval <= 0.05, empirical_pval = emp_pval)
    }
    
    perm_res <- perm_test(knn30, mono_mat$seed_idx, mono_mat$np_score, n_perm, mycores)
    mono_final <- data.frame(mono_mat, perm_res)
    mono_final$fdr <- p.adjust(mono_final$empirical_pval, method = "BH")
    
    list(success = TRUE, result = mono_final,
         n_significant = sum(mono_final$true_cell_top_idx),
         n_cells = nrow(mono_final), n_peaks = nrow(SE_sub),
         n_peaks_with_signal = n_peaks_signal,
         n_seed = length(seed_cells), median_TRS = med_TRS)
    
  }, error = function(e) list(success = FALSE, error = e$message))
  
  return(result)
}

#===============================================================================
# STEP 9: Run SCAVENGE Per Trait
#===============================================================================
write_log("\n====== Step 9: Run SCAVENGE ======")

results_file <- file.path(results_dir, paste0("results_Liver_P2G_", CELL_TYPE, "_ALL.rds"))
summary_file <- file.path(results_dir, paste0("summary_Liver_P2G_", CELL_TYPE, ".csv"))
skipped_file <- file.path(results_dir, paste0("skipped_Liver_P2G_", CELL_TYPE, ".csv"))

# Resume from previous run
if (file.exists(results_file)) {
  scavenge_results <- readRDS(results_file)
  scavenge_summary <- fread(summary_file)
  skipped_runs <- if (file.exists(skipped_file)) fread(skipped_file) else data.table()
  write_log(paste("  Resuming:", nrow(scavenge_summary), "completed"))
} else {
  scavenge_results <- list()
  scavenge_summary <- data.table()
  skipped_runs     <- data.table()
}

completed <- scavenge_summary$trait
remaining <- setdiff(names(trait_bed_files), completed)

write_log(paste("  Completed:", length(completed), "| Remaining:", length(remaining)))

if (length(remaining) > 0) {
  
  t0 <- Sys.time()
  
  for (i in seq_along(remaining)) {
    
    trait <- remaining[i]
    write_log(paste0("\n  [", i, "/", length(remaining), "] ", trait))
    
    res <- tryCatch(
      run_scavenge(SE_Data, trait_bed_files[[trait]],
                   params$scavenge_cores, params$min_seed_cells,
                   params$min_valid_cells, params$permutations),
      error = function(e) list(success = FALSE, error = e$message)
    )
    
    if (res$success) {
      df <- res$result
      scavenge_results[[trait]] <- df
      
      n_fdr <- sum(df$fdr < params$fdr_threshold, na.rm = TRUE)
      
      # Get region count for this trait
      tc_row <- trait_counts[trait_counts$trait == trait]
      n_reg <- if (nrow(tc_row) > 0) tc_row$n_regions[1] else NA
      
      scavenge_summary <- rbind(scavenge_summary, data.table(
        trait               = trait,
        cell_type           = CELL_TYPE,
        n_regions           = n_reg,
        n_peaks             = res$n_peaks,
        n_peaks_with_signal = res$n_peaks_with_signal,
        n_cells             = res$n_cells,
        n_seed              = res$n_seed,
        mean_zscore         = round(mean(df$z_score, na.rm = TRUE), 4),
        n_fdr_sig           = n_fdr,
        pct_fdr_sig         = round(n_fdr / res$n_cells * 100, 2),
        n_true_cell         = sum(df$true_cell_top_idx, na.rm = TRUE),
        pct_true_cell       = round(sum(df$true_cell_top_idx, na.rm = TRUE) / res$n_cells * 100, 2),
        median_TRS          = round(res$median_TRS, 6),
        mean_TRS            = round(mean(df$TRS, na.rm = TRUE), 6)
      ), fill = TRUE)
      
      write_log(paste("    OK | FDR<0.05:", n_fdr, "| TRS:", round(res$median_TRS, 4)))
      
      # Save per-trait raw & filtered
      df$trait <- trait; df$cell_type <- CELL_TYPE
      fwrite(df, file.path(raw_dir, paste0(clean_trait_name(trait), "_raw.csv")))
      
      df_filt <- df[df$fdr < params$fdr_threshold, ]
      if (nrow(df_filt) > 0) {
        fwrite(df_filt, file.path(filtered_dir, paste0(clean_trait_name(trait), "_filtered.csv")))
      }
      
    } else {
      write_log(paste("    SKIP:", res$error))
      skipped_runs <- rbind(skipped_runs, data.table(
        trait = trait, cell_type = CELL_TYPE, reason = res$error
      ), fill = TRUE)
    }
    
    # Periodic checkpoint every 5 traits
    if (i %% 5 == 0 || i == length(remaining)) {
      saveRDS(scavenge_results, results_file)
      fwrite(scavenge_summary, summary_file)
      fwrite(skipped_runs, skipped_file)
      elapsed <- round(as.numeric(difftime(Sys.time(), t0, units = "hours")), 2)
      write_log(paste("    [checkpoint] elapsed:", elapsed, "h"))
    }
  }
}

# Final save
saveRDS(scavenge_results, results_file)
fwrite(scavenge_summary, summary_file)
fwrite(skipped_runs, skipped_file)

#===============================================================================
# STEP 10: Combine Per-Trait Results
#===============================================================================
write_log("\n====== Step 10: Combine Results ======")

raw_files <- list.files(raw_dir, "_raw.csv$", full.names = TRUE)
if (length(raw_files) > 0) {
  all_raw <- rbindlist(lapply(raw_files, fread), fill = TRUE)
  fwrite(all_raw, file.path(results_dir, paste0("combined_ALL_raw_", CELL_TYPE, ".csv")))
  write_log(paste("  Raw combined:", nrow(all_raw), "rows from", length(raw_files), "traits"))
}

filt_files <- list.files(filtered_dir, "_filtered.csv$", full.names = TRUE)
if (length(filt_files) > 0) {
  all_filt <- rbindlist(lapply(filt_files, fread), fill = TRUE)
  fwrite(all_filt, file.path(results_dir, paste0("combined_ALL_filtered_", CELL_TYPE, ".csv")))
  write_log(paste("  Filtered combined:", nrow(all_filt), "rows from", length(filt_files), "traits"))
}

#===============================================================================
# SUMMARY
#===============================================================================
write_log("\n================================================================")
write_log(paste("COMPLETE:", CELL_TYPE))
write_log("================================================================")
write_log(paste("  Successful:", nrow(scavenge_summary)))
write_log(paste("  Skipped:   ", nrow(skipped_runs)))
write_log(paste("  Results:   ", results_file))
write_log(paste("  Summary:   ", summary_file))

if (nrow(scavenge_summary) > 0) {
  write_log("\n  Top 10 traits by % FDR-significant:")
  top <- scavenge_summary[order(-pct_fdr_sig)][1:min(10, .N)]
  for (j in 1:nrow(top)) {
    write_log(sprintf("    %2d. %-45s  FDR<0.05: %5d (%5.1f%%)  TRS: %.4f",
                      j, substr(top$trait[j], 1, 45),
                      top$n_fdr_sig[j], top$pct_fdr_sig[j], top$median_TRS[j]))
  }
}

gc()
write_log(paste("\nFinished:", Sys.time()))
write_log("================================================================")
