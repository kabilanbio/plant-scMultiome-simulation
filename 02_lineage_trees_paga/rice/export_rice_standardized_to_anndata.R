# ==============================================================================
# export_rice_standardized_to_anndata.R
# Converts Standardized Rice Multiomics Datasets (Seurat RDS) to AnnData / Matrix Market
# Supports all 8 Rice tissues: Bud, Flag, Leaf, Root, SAM, Seed, SP, ST
# Ensures 100% cross-platform compatibility between R and Python / Scanpy / PAGA
# ==============================================================================

suppressPackageStartupMessages({
  library(Matrix)
})

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  
  # Default paths
  default_input <- "G:/PhD/sc_datasets/rice/featured_datasets_standardized"
  default_output <- "G:/PhD/sc_datasets/rice/featured_datasets_anndata"
  
  opt <- list(
    input_dir  = default_input,
    output_dir = default_output,
    tissues    = c("Bud", "Flag", "Leaf", "Root", "SAM", "Seed", "SP", "ST")
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
        opt$tissues <- c("Bud", "Flag", "Leaf", "Root", "SAM", "Seed", "SP", "ST")
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
message("Exporting Standardized Rice Datasets for PAGA Analysis")
message("Input directory:  ", opt$input_dir)
message("Output directory: ", opt$output_dir)
message("Tissues:          ", paste(opt$tissues, collapse = ", "))
message("==================================================================")

dir.create(opt$output_dir, recursive = TRUE, showWarnings = FALSE)

# Check for direct H5AD export support via zellkonverter or reticulate/anndata
has_zellkonverter <- requireNamespace("zellkonverter", quietly = TRUE)
has_reticulate    <- requireNamespace("reticulate", quietly = TRUE)

export_single_tissue <- function(tissue, input_dir, output_dir) {
  message(sprintf("\n>>> Processing Tissue: %s <<<", tissue))
  
  tissue_dir <- file.path(input_dir, tissue)
  if (!dir.exists(tissue_dir)) {
    warning(sprintf("Tissue directory does not exist: %s", tissue_dir))
    return(FALSE)
  }
  
  # Candidate files
  rna_rds      <- file.path(tissue_dir, sprintf("rice_processed_%s_rna.rds", tissue))
  counts_rds   <- file.path(tissue_dir, sprintf("rice_processed_%s_rna_counts.rds", tissue))
  metadata_csv <- file.path(tissue_dir, sprintf("rice_processed_%s_metadata.csv", tissue))
  celltype_csv <- file.path(tissue_dir, sprintf("rice_processed_%s_cell_types.csv", tissue))
  
  # Load count matrix
  counts <- NULL
  meta   <- NULL
  
  if (file.exists(counts_rds)) {
    message("  Loading RNA counts matrix: ", basename(counts_rds))
    counts <- readRDS(counts_rds)
  } else if (file.exists(rna_rds)) {
    message("  Loading Seurat object: ", basename(rna_rds))
    suppressPackageStartupMessages(library(Seurat))
    obj <- readRDS(rna_rds)
    counts <- GetAssayData(obj, assay = "RNA", slot = "counts")
    meta <- obj@meta.data
  } else {
    warning("  No RNA counts or Seurat object found for tissue: ", tissue)
    return(FALSE)
  }
  
  # Load metadata
  if (is.null(meta) && file.exists(metadata_csv)) {
    message("  Loading metadata: ", basename(metadata_csv))
    meta <- read.csv(metadata_csv, row.names = 1, check.names = FALSE)
  }
  
  # Ensure cell_type and batch columns exist
  if (!"cell_type" %in% colnames(meta) && file.exists(celltype_csv)) {
    ct_df <- read.csv(celltype_csv)
    rownames(ct_df) <- ct_df$CellBarcode
    meta$cell_type <- ct_df[rownames(meta), "CellType"]
  }
  
  # Harmonize barcodes
  common_cells <- intersect(colnames(counts), rownames(meta))
  if (length(common_cells) == 0) {
    stop("No matching cell barcodes between counts matrix and metadata!")
  }
  
  counts <- counts[, common_cells, drop = FALSE]
  meta   <- meta[common_cells, , drop = FALSE]
  
  message(sprintf("  Validated dimensions: %d genes x %d cells across %d cell types", 
                  nrow(counts), ncol(counts), length(unique(meta$cell_type))))
  
  # Create tissue output folder
  tissue_out <- file.path(output_dir, tissue)
  dir.create(tissue_out, recursive = TRUE, showWarnings = FALSE)
  
  # Strategy 1: Direct H5AD if zellkonverter is present
  exported_h5ad <- FALSE
  h5ad_file <- file.path(tissue_out, sprintf("rice_processed_%s.h5ad", tissue))
  
  if (has_zellkonverter && requireNamespace("SingleCellExperiment", quietly = TRUE)) {
    tryCatch({
      message("  Exporting directly to H5AD via zellkonverter...")
      sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = counts),
        colData = meta
      )
      zellkonverter::writeH5AD(sce, file = h5ad_file)
      message("  Successfully wrote H5AD: ", h5ad_file)
      exported_h5ad <- TRUE
    }, error = function(e) {
      message("  zellkonverter export failed: ", e$message)
    })
  }
  
  # Strategy 2: Export 10x Matrix Market format + Metadata CSV (100% universal)
  # Python Scanpy reads this instantly with sc.read_10x_mtx or scanpy.read_mtx
  mm_dir <- file.path(tissue_out, "matrix_market")
  dir.create(mm_dir, recursive = TRUE, showWarnings = FALSE)
  
  message("  Exporting Matrix Market archive for Python Scanpy/AnnData...")
  Matrix::writeMM(counts, file = file.path(mm_dir, "matrix.mtx"))
  writeLines(rownames(counts), con = file.path(mm_dir, "genes.tsv"))
  writeLines(colnames(counts), con = file.path(mm_dir, "barcodes.tsv"))
  write.csv(meta, file = file.path(mm_dir, "metadata.csv"), quote = TRUE)
  
  # Write a helper JSON manifest
  manifest <- list(
    tissue = tissue,
    n_genes = nrow(counts),
    n_cells = ncol(counts),
    cell_types = unique(as.character(meta$cell_type)),
    batches = if ("batch" %in% colnames(meta)) unique(as.character(meta$batch)) else character(0),
    h5ad_exported = exported_h5ad
  )
  manifest_file <- file.path(tissue_out, "dataset_manifest.json")
  # Simple JSON writer
  json_str <- sprintf(
    '{\n  "tissue": "%s",\n  "n_genes": %d,\n  "n_cells": %d,\n  "n_cell_types": %d,\n  "h5ad_exported": %s\n}',
    tissue, nrow(counts), ncol(counts), length(unique(meta$cell_type)), tolower(as.character(exported_h5ad))
  )
  writeLines(json_str, con = manifest_file)
  
  message(sprintf("  Successfully exported %s to: %s", tissue, tissue_out))
  return(TRUE)
}

# Run loop over requested tissues
results <- sapply(opt$tissues, function(t) {
  export_single_tissue(t, opt$input_dir, opt$output_dir)
})

message("\n==================================================================")
message("Export Summary:")
for (i in seq_along(results)) {
  message(sprintf("  %-6s : %s", names(results)[i], if (results[i]) "SUCCESS" else "FAILED"))
}
message("Output directory: ", opt$output_dir)
message("==================================================================")
