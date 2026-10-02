# ==============================================================================
# export_soybean_standardized_to_anndata.R
# Converts Standardized Soybean Multiomics Datasets (Seurat RDS) to AnnData / Matrix Market
# Supports all 7 Soybean tissues:
#   Cotyledon_stage_seeds, Early_maturation_stage_seeds, Early_nodule,
#   Globular_stage_seeds, Heart_stage_seeds, Hypocotyl, Root
# Optimized for HP 280 G4 Workstation (R 4.6.1 + Python sc_env)
# ==============================================================================

suppressPackageStartupMessages({
  library(Matrix)
})

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  
  script_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", script_args, value = TRUE)
  this_dir <- if (length(file_arg) > 0) normalizePath(dirname(sub("^--file=", "", file_arg[1]))) else getwd()
  drive_prefix <- substr(this_dir, 1, 2)
  if (!grepl("^[A-Za-z]:", drive_prefix)) drive_prefix <- "G:"
  
  default_input <- sprintf("%s/PhD/sc_datasets/soybean/featured_datasets_standardized", drive_prefix)
  default_output <- sprintf("%s/PhD/sc_datasets/soybean/featured_datasets_anndata", drive_prefix)
  
  all_soytissues <- c(
    "Cotyledon_stage_seeds",
    "Early_maturation_stage_seeds",
    "Early_nodule",
    "Globular_stage_seeds",
    "Heart_stage_seeds",
    "Hypocotyl",
    "Root"
  )
  
  opt <- list(
    input_dir  = default_input,
    output_dir = default_output,
    tissues    = all_soytissues
  )
  
  i <- 1
  while (i <= length(args)) {
    arg <- args[i]
    if (arg %in% c("-i", "--input_dir") && i < length(args)) {
      opt$input_dir <- args[i + 1]
      i <- i + 2
    } else if (arg %in% c("-o", "--output_dir") && i < length(args)) {
      opt$output_dir <- args[i + 1]
      i <- i + 2
    } else if (arg %in% c("-t", "--tissue") && i < length(args)) {
      val <- args[i + 1]
      if (val == "all") {
        opt$tissues <- all_soytissues
      } else {
        opt$tissues <- strsplit(val, ",")[[1]]
      }
      i <- i + 2
    } else {
      i <- i + 1
    }
  }
  return(opt)
}

opt <- parse_args()

message("==================================================================")
message("Exporting Standardized Soybean Datasets for PAGA Analysis")
message("Input directory:  ", opt$input_dir)
message("Output directory: ", opt$output_dir)
message("Tissues:          ", paste(opt$tissues, collapse = ", "))
message("==================================================================")

dir.create(opt$output_dir, recursive = TRUE, showWarnings = FALSE)
has_zellkonverter <- requireNamespace("zellkonverter", quietly = TRUE)

export_single_tissue <- function(tissue, input_dir, output_dir) {
  message(sprintf("\n>>> Processing Soybean Tissue: %s <<<", tissue))
  tissue_dir <- file.path(input_dir, tissue)
  if (!dir.exists(tissue_dir)) {
    warning("Directory does not exist: ", tissue_dir)
    return(FALSE)
  }
  
  counts_rds   <- file.path(tissue_dir, sprintf("soy_processed_%s_rna_counts.rds", tissue))
  rna_rds      <- file.path(tissue_dir, sprintf("soy_processed_%s_rna.rds", tissue))
  metadata_csv <- file.path(tissue_dir, sprintf("soy_processed_%s_rna_metadata.csv", tissue))
  celltypes_csv<- file.path(tissue_dir, sprintf("soy_processed_%s_cell_types.csv", tissue))
  
  counts <- NULL
  meta   <- NULL
  
  if (file.exists(counts_rds)) {
    message("  Loading RNA counts: ", basename(counts_rds))
    counts <- readRDS(counts_rds)
  } else if (file.exists(rna_rds)) {
    message("  Loading Seurat object: ", basename(rna_rds))
    suppressPackageStartupMessages(library(Seurat))
    obj <- readRDS(rna_rds)
    counts <- GetAssayData(obj, assay = "RNA", slot = "counts")
    meta <- obj@meta.data
  } else {
    warning("  Counts not found for: ", tissue)
    return(FALSE)
  }
  
  if (is.null(meta) && file.exists(metadata_csv)) {
    message("  Loading metadata: ", basename(metadata_csv))
    meta <- read.csv(metadata_csv, row.names = 1, check.names = FALSE)
  }
  
  if (!"cell_type" %in% colnames(meta) && file.exists(celltypes_csv)) {
    ct_df <- read.csv(celltypes_csv)
    rownames(ct_df) <- ct_df$CellBarcode
    meta$cell_type <- ct_df[rownames(meta), "CellType"]
  }
  
  common_cells <- intersect(colnames(counts), rownames(meta))
  if (length(common_cells) == 0) {
    stop("No overlapping barcodes between counts and metadata for ", tissue)
  }
  
  counts <- counts[, common_cells, drop = FALSE]
  meta   <- meta[common_cells, , drop = FALSE]
  
  message(sprintf("  Validated dimensions: %d genes x %d cells across %d cell types",
                  nrow(counts), ncol(counts), length(unique(meta$cell_type))))
                  
  tissue_out <- file.path(output_dir, tissue)
  dir.create(tissue_out, recursive = TRUE, showWarnings = FALSE)
  
  exported_h5ad <- FALSE
  h5ad_file <- file.path(tissue_out, sprintf("soy_processed_%s.h5ad", tissue))
  if (has_zellkonverter && requireNamespace("SingleCellExperiment", quietly = TRUE)) {
    tryCatch({
      sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = counts),
        colData = meta
      )
      zellkonverter::writeH5AD(sce, file = h5ad_file)
      exported_h5ad <- TRUE
      message("  Wrote H5AD: ", basename(h5ad_file))
    }, error = function(e) {})
  }
  
  mm_dir <- file.path(tissue_out, "matrix_market")
  dir.create(mm_dir, recursive = TRUE, showWarnings = FALSE)
  Matrix::writeMM(counts, file = file.path(mm_dir, "matrix.mtx"))
  writeLines(rownames(counts), con = file.path(mm_dir, "genes.tsv"))
  writeLines(colnames(counts), con = file.path(mm_dir, "barcodes.tsv"))
  write.csv(meta, file = file.path(mm_dir, "metadata.csv"), quote = TRUE)
  
  manifest_file <- file.path(tissue_out, "dataset_manifest.json")
  json_str <- sprintf(
    '{\n  "species": "Glycine max",\n  "tissue": "%s",\n  "n_genes": %d,\n  "n_cells": %d,\n  "n_cell_types": %d,\n  "h5ad_exported": %s\n}',
    tissue, nrow(counts), ncol(counts), length(unique(meta$cell_type)), tolower(as.character(exported_h5ad))
  )
  writeLines(json_str, con = manifest_file)
  return(TRUE)
}

results <- sapply(opt$tissues, function(t) {
  export_single_tissue(t, opt$input_dir, opt$output_dir)
})

message("\n==================================================================")
message("Soybean Export Summary:")
for (i in seq_along(results)) {
  message(sprintf("  %-30s : %s", names(results)[i], if (results[i]) "SUCCESS" else "FAILED"))
}
message("Output directory: ", opt$output_dir)
message("==================================================================")
