# ==============================================================================
# export_arab_standardized_to_anndata.R
# Converts Standardized Arabidopsis Root Dataset (Seurat RDS) to AnnData / Matrix Market
# Automatically applies curated biological cell type and tissue layer annotations:
#   - 21 Clusters (1-21) -> 6 Biological Cell Types -> 4 Tissue Layers
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
  
  default_input <- sprintf("%s/PhD/sc_datasets/arabidopsis_root/featured_datasets_standardized", drive_prefix)
  default_output <- sprintf("%s/PhD/sc_datasets/arabidopsis_root/featured_datasets_anndata", drive_prefix)
  
  opt <- list(
    input_dir  = default_input,
    output_dir = default_output
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
    } else {
      i <- i + 1
    }
  }
  return(opt)
}

opt <- parse_args()

message("==================================================================")
message("Exporting Standardized Arabidopsis Root with Informative Cell Types")
message("Input directory:  ", opt$input_dir)
message("Output directory: ", opt$output_dir)
message("==================================================================")

dir.create(opt$output_dir, recursive = TRUE, showWarnings = FALSE)

tissue <- "Root"
tissue_dir <- file.path(opt$input_dir, tissue)
if (!dir.exists(tissue_dir)) {
  if (file.exists(file.path(opt$input_dir, "arab_processed_Root_rna_counts.rds"))) {
    tissue_dir <- opt$input_dir
  } else {
    stop("Could not locate Arabidopsis Root directory at: ", tissue_dir)
  }
}

counts_rds   <- file.path(tissue_dir, "arab_processed_Root_rna_counts.rds")
rna_rds      <- file.path(tissue_dir, "arab_processed_Root_rna.rds")
metadata_csv <- file.path(tissue_dir, "arab_processed_Root_rna_metadata.csv")

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
  stop("Counts matrix not found at: ", counts_rds)
}

if (is.null(meta) && file.exists(metadata_csv)) {
  message("  Loading metadata: ", basename(metadata_csv))
  meta <- read.csv(metadata_csv, row.names = 1, check.names = FALSE)
}

common_cells <- intersect(colnames(counts), rownames(meta))
if (length(common_cells) == 0) {
  stop("No overlapping barcodes between counts and metadata!")
}

counts <- counts[, common_cells, drop = FALSE]
meta   <- meta[common_cells, , drop = FALSE]

# ------------------------------------------------------------------------------
# Biological Annotation Mapping (Curated from Farmer et al., 2021)
# ------------------------------------------------------------------------------
cluster_to_celltype_str <- c(
  'k0_11' = 'Cortical cell',        'o1_15' = 'Endodermal cell',  
  'e2_5'  = 'Atrichoblast',         'r3_18' = 'Stele cell',  
  'g4_7'  = 'Atrichoblast',         'b5_2'  = 'Trichoblast',  
  'h6_8'  = 'Meristematic cell',    'q7_17' = 'Stele cell',  
  'f8_6'  = 'Atrichoblast',         'i9_9'  = 'Meristematic cell',  
  'l10_12'= 'Cortical cell',        'm11_13'= 'Endodermal cell',  
  'd12_4' = 'Atrichoblast',         'c13_3' = 'Trichoblast',  
  'a14_1' = 'Trichoblast',          'u15_21'= 'Stele cell',  
  'p16_16'= 'Endodermal cell',      'n17_14'= 'Endodermal cell',  
  's18_19'= 'Stele cell',           't19_20'= 'Stele cell',  
  'j20_10'= 'Meristematic cell'
)

celltype_to_layer_str <- c(
  'Trichoblast'       = 'Epidermis',
  'Atrichoblast'      = 'Epidermis',
  'Meristematic cell' = 'Meristem',
  'Cortical cell'     = 'Ground tissue',
  'Endodermal cell'   = 'Ground tissue',
  'Stele cell'        = 'Vasculature'
)

# Extract cluster number (1-21)
cluster_key <- if ("seurat_clusters_renamed" %in% colnames(meta)) {
  as.character(meta$seurat_clusters_renamed)
} else {
  as.character(meta$cell_type)
}

# Add informative annotations
meta$seurat_id        <- cluster_key
meta$cell_type        <- cluster_to_celltype_str[cluster_key]
meta$tissue_layer     <- celltype_to_layer_str[meta$cell_type]
meta$cluster_num      <- as.integer(sub(".*_([0-9]+)$", "\\1", cluster_key))
meta$cluster_tip      <- paste0("cluster", meta$cluster_num)
meta$informative_label<- paste0(gsub(" ", "_", meta$cell_type), "_c", meta$cluster_num)

message(sprintf("  Successfully annotated %d cells across 21 developmental clusters:", nrow(meta)))
for (ct in unique(meta$cell_type)) {
  n_ct <- sum(meta$cell_type == ct)
  message(sprintf("    * %-18s (%s) : %d cells", ct, celltype_to_layer_str[ct], n_ct))
}

tissue_out <- file.path(opt$output_dir, tissue)
dir.create(tissue_out, recursive = TRUE, showWarnings = FALSE)

# Export H5AD if zellkonverter is present
has_zellkonverter <- requireNamespace("zellkonverter", quietly = TRUE)
exported_h5ad <- FALSE
h5ad_file <- file.path(tissue_out, "arab_processed_Root.h5ad")

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

# Export Matrix Market
mm_dir <- file.path(tissue_out, "matrix_market")
dir.create(mm_dir, recursive = TRUE, showWarnings = FALSE)
Matrix::writeMM(counts, file = file.path(mm_dir, "matrix.mtx"))
writeLines(rownames(counts), con = file.path(mm_dir, "genes.tsv"))
writeLines(colnames(counts), con = file.path(mm_dir, "barcodes.tsv"))
write.csv(meta, file = file.path(mm_dir, "metadata.csv"), quote = TRUE)

manifest_file <- file.path(tissue_out, "dataset_manifest.json")
json_str <- sprintf(
  '{\n  "species": "Arabidopsis thaliana",\n  "tissue": "Root",\n  "n_genes": %d,\n  "n_cells": %d,\n  "n_clusters": 21,\n  "n_cell_types": 6,\n  "h5ad_exported": %s\n}',
  nrow(counts), ncol(counts), tolower(as.character(exported_h5ad))
)
writeLines(json_str, con = manifest_file)

message("\n==================================================================")
message("SUCCESS! Exported Arabidopsis Root dataset to: ", tissue_out)
message("==================================================================")
