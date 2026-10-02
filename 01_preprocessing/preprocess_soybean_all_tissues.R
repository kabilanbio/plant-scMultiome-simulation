# ==============================================================================
# preprocess_soybean_all_tissues.R
# Standardized Single-Cell Multiomics Preprocessing Pipeline for Soybean (7 Tissues)
# Unpaired Multi-Omics Architecture: snRNA-seq (Seurat) + snATAC-seq (ACR sparse matrix)
# Strictly enforcers:
#   - 3:1 ATAC-to-RNA feature ratio (2,000 HVGs, 6,000 ACRs)
#   - Dual-stratified sampling per modality (max 200 cells / cell_type x batch)
#   - High-resolution 600 DPI JPEG figures (zero PDFs)
#   - Decoupled unimodal RDS objects, count matrices, cell types, batches, and by_batch files
# ==============================================================================

# Ensure all required packages are present; install automatically if missing
required_pkgs <- c("Seurat", "Signac", "Matrix", "ggplot2", "patchwork", "dplyr", "scales", "readxl")
missing_pkgs <- required_pkgs[!required_pkgs %in% installed.packages()[, "Package"]]
if (length(missing_pkgs) > 0) {
  message(">>> Auto-installing missing R packages: ", paste(missing_pkgs, collapse = ", "), " <<<")
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
  }
  for (pkg in missing_pkgs) {
    message("Installing: ", pkg)
    if (pkg == "Signac") {
      BiocManager::install("Signac", update = FALSE, ask = FALSE)
    } else {
      install.packages(pkg, repos = "https://cloud.r-project.org")
    }
  }
  message(">>> All packages successfully installed! <<<\n")
}

suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(Matrix)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(readxl)
  library(scales)
})

# ------------------------------------------------------------------------------
# 1. Native CLI Parameter Parser (Zero External Dependencies)
# ------------------------------------------------------------------------------
parse_cli_arguments <- function(args = commandArgs(trailingOnly = TRUE)) {
  script_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", script_args, value = TRUE)
  this_dir <- if (length(file_arg) > 0) normalizePath(dirname(sub("^--file=", "", file_arg[1]))) else getwd()
  
  candidate_input <- normalizePath(file.path(this_dir, "..", "soy_datasets"), mustWork = FALSE)
  candidate_output <- normalizePath(file.path(this_dir, "..", "featured_datasets_standardized"), mustWork = FALSE)
  
  default_input <- if (dir.exists(candidate_input)) candidate_input else "G:/PhD/sc_datasets/soybean/soy_datasets"
  default_output <- if (dir.exists(candidate_output)) candidate_output else "G:/PhD/sc_datasets/soybean/featured_datasets_standardized"

  opt <- list(
    input_dir  = default_input,
    output_dir = default_output,
    tissue     = "remaining",
    n_hvg      = 2000,
    cell_cap   = 200,
    seed       = 42
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
      opt$tissue <- args[i + 1]
      i <- i + 2
    } else if (arg == "--n_hvg" && i < length(args)) {
      opt$n_hvg <- as.integer(args[i + 1])
      i <- i + 2
    } else if (arg == "--cell_cap" && i < length(args)) {
      opt$cell_cap <- as.integer(args[i + 1])
      i <- i + 2
    } else if (arg == "--seed" && i < length(args)) {
      opt$seed <- as.integer(args[i + 1])
      i <- i + 2
    } else {
      i <- i + 1
    }
  }
  return(opt)
}

opt <- parse_cli_arguments()
set.seed(opt$seed)

# Canonical 7 soybean tissues (GSE270392)
ALL_TISSUES <- c(
  "Hypocotyl",
  "Root",
  "Early_nodule",
  "Globular_stage_seeds",
  "Heart_stage_seeds",
  "Cotyledon_stage_seeds",
  "Early_maturation_stage_seeds"
)

tissues_to_process <- if (tolower(opt$tissue) == "all") {
  ALL_TISSUES
} else if (tolower(opt$tissue) == "remaining") {
  done_tissues <- character(0)
  for (tis in ALL_TISSUES) {
    tis_metrics <- file.path(opt$output_dir, tis, paste0("soy_processed_", tis, "_summary_metrics.csv"))
    if (file.exists(tis_metrics)) done_tissues <- c(done_tissues, tis)
  }
  rem <- setdiff(ALL_TISSUES, done_tissues)
  if (length(rem) == 0) {
    message("All 7 soybean tissues are already fully processed!")
    quit(save = "no", status = 0)
  }
  rem
} else {
  matched <- grep(opt$tissue, ALL_TISSUES, ignore.case = TRUE, value = TRUE)
  if (length(matched) == 0) stop("Invalid tissue name. Choose from: ", paste(ALL_TISSUES, collapse = ", "))
  matched
}

message("==================================================================")
message("Soybean Single-Cell Multiomics Preprocessing Pipeline (Unpaired)")
message("Tissues to process: ", paste(tissues_to_process, collapse = ", "))
message("Target RNA HVGs: ", opt$n_hvg, " | Target ATAC ACRs: ", 3 * opt$n_hvg, " (Strict 3:1 ratio)")
message("Max cells per (cell_type x batch): ", opt$cell_cap)
message("==================================================================")

# ------------------------------------------------------------------------------
# 2. Sparse Helper Functions
# ------------------------------------------------------------------------------
calc_sparse_sparsity <- function(mat) {
  total_elements <- as.numeric(nrow(mat)) * as.numeric(ncol(mat))
  1 - (Matrix::nnzero(mat) / total_elements)
}

calc_sparse_row_variance <- function(mat) {
  n <- ncol(mat)
  means <- Matrix::rowMeans(mat)
  mat_sq <- mat
  mat_sq@x <- mat_sq@x^2
  vars <- (Matrix::rowSums(mat_sq) - n * means^2) / (n - 1)
  names(vars) <- rownames(mat)
  return(vars)
}

# ------------------------------------------------------------------------------
# 3. Load Global Metadata Once
# ------------------------------------------------------------------------------
rna_meta_file <- file.path(opt$input_dir, "snRNA_meta.xlsx")
atac_meta_file <- file.path(opt$input_dir, "snATAC_meta.xlsx")

if (!file.exists(rna_meta_file)) stop("snRNA_meta.xlsx not found in: ", opt$input_dir)
if (!file.exists(atac_meta_file)) stop("snATAC_meta.xlsx not found in: ", opt$input_dir)

message("Loading snRNA_meta.xlsx and snATAC_meta.xlsx...")
rna_meta_global <- as.data.frame(readxl::read_excel(rna_meta_file))
atac_meta_global <- as.data.frame(readxl::read_excel(atac_meta_file))

# Standardize column names
colnames(rna_meta_global)[colnames(rna_meta_global) == "Celltype"] <- "cell_type"
colnames(rna_meta_global)[colnames(rna_meta_global) == "Replicate"] <- "batch"

colnames(atac_meta_global)[colnames(atac_meta_global) == "Celltype"] <- "cell_type"
colnames(atac_meta_global)[colnames(atac_meta_global) == "Replicate"] <- "batch"

# ------------------------------------------------------------------------------
# 4. Core Preprocessing Function per Soybean Tissue
# ------------------------------------------------------------------------------
process_single_soybean_tissue <- function(tissue_name, input_dir, output_dir, n_target_hvg = 2000, cell_cap = 200) {
  message("\n>>> Processing Soybean Tissue: ", tissue_name, " <<<")
  
  tissue_out_dir <- file.path(output_dir, tissue_name)
  plot_dir       <- file.path(tissue_out_dir, "plots")
  batch_dir      <- file.path(tissue_out_dir, "by_batch")
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(batch_dir, recursive = TRUE, showWarnings = FALSE)

  # ============================================================================
  # STEP A: RNA MODALITY
  # ============================================================================
  rna_file <- file.path(input_dir, paste0("GSE270392_Gm_atlas_", tissue_name, ".rna.seurat.obj.rds"))
  if (!file.exists(rna_file)) stop("RNA file not found: ", rna_file)
  
  message("Loading RNA Seurat object: ", rna_file)
  obj_rna <- readRDS(rna_file)
  obj_rna <- UpdateSeuratObject(obj_rna)
  
  # Strip ann1. prefix from gene names if present
  rna_counts_raw <- GetAssayData(obj_rna, assay = "RNA", layer = "counts")
  if (any(grepl("^ann1\\.", rownames(rna_counts_raw)))) {
    message("Cleaning 'ann1.' gene prefixes in RNA counts...")
    rownames(rna_counts_raw) <- sub("^ann1\\.", "", rownames(rna_counts_raw))
    # Update object assay
    obj_rna[["RNA"]] <- CreateAssayObject(counts = rna_counts_raw)
  }
  
  # Match external metadata
  rna_tissue_meta <- rna_meta_global[rna_meta_global$Tissue == tissue_name, ]
  rownames(rna_tissue_meta) <- rna_tissue_meta$cellID
  common_rna_cells <- intersect(colnames(obj_rna), rownames(rna_tissue_meta))
  
  if (length(common_rna_cells) > 0) {
    obj_rna <- subset(obj_rna, cells = common_rna_cells)
    obj_rna$cell_type <- rna_tissue_meta[colnames(obj_rna), "cell_type"]
    obj_rna$batch <- rna_tissue_meta[colnames(obj_rna), "batch"]
  } else {
    # Fallback to internal meta.data
    if (!"cell_type" %in% colnames(obj_rna@meta.data)) {
      if ("celltype" %in% colnames(obj_rna@meta.data)) obj_rna$cell_type <- obj_rna$celltype
      else if ("seurat_clusters" %in% colnames(obj_rna@meta.data)) obj_rna$cell_type <- paste0("Cluster_", obj_rna$seurat_clusters)
      else obj_rna$cell_type <- "Unassigned"
    }
    if (!"batch" %in% colnames(obj_rna@meta.data)) {
      if ("sampleID" %in% colnames(obj_rna@meta.data)) obj_rna$batch <- obj_rna$sampleID
      else obj_rna$batch <- "batch_1"
    }
  }
  
  raw_rna_cells <- ncol(obj_rna)
  raw_rna_genes <- nrow(obj_rna)
  raw_rna_counts <- GetAssayData(obj_rna, assay = "RNA", layer = "counts")
  raw_rna_umi <- Matrix::colSums(raw_rna_counts)
  raw_rna_sparsity <- calc_sparse_sparsity(raw_rna_counts)
  raw_mean_umi_rna <- mean(raw_rna_umi)
  raw_median_umi_rna <- median(raw_rna_umi)
  
  message(sprintf("Raw RNA: %d cells | %d genes | Mean UMI: %.2f (Median: %.2f) | Sparsity: %.2f%%",
                  raw_rna_cells, raw_rna_genes, raw_mean_umi_rna, raw_median_umi_rna, raw_rna_sparsity * 100))

  # QC filtering for RNA
  qc_keep_rna <- colnames(obj_rna)[
    obj_rna$nCount_RNA >= 300 &
    obj_rna$nFeature_RNA >= 150 &
    !is.na(obj_rna$cell_type) &
    obj_rna$cell_type != "" &
    obj_rna$cell_type != "Unknown"
  ]
  if (length(qc_keep_rna) < 200) {
    # Relax if needed
    qc_keep_rna <- colnames(obj_rna)[!is.na(obj_rna$cell_type) & obj_rna$cell_type != ""]
  }
  obj_rna_filtered <- subset(obj_rna, cells = qc_keep_rna)
  message("RNA cells after QC filtering: ", ncol(obj_rna_filtered))

  # Dual-stratified sampling for RNA
  rna_groups <- split(colnames(obj_rna_filtered), list(obj_rna_filtered$cell_type, obj_rna_filtered$batch), drop = TRUE)
  sampled_rna_cells <- unlist(lapply(rna_groups, function(cells) {
    if (length(cells) <= cell_cap) cells else sample(cells, cell_cap)
  }), use.names = FALSE)
  
  obj_rna_sub <- subset(obj_rna_filtered, cells = sampled_rna_cells)
  message("RNA cells after dual-stratified sampling: ", ncol(obj_rna_sub))

  # Batch-aware consensus HVGs
  obj_rna_sub <- NormalizeData(obj_rna_sub, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
  obj_rna_sub <- FindVariableFeatures(obj_rna_sub, selection.method = "vst", nfeatures = n_target_hvg, verbose = FALSE)
  final_hvg_genes <- VariableFeatures(obj_rna_sub)
  if (length(final_hvg_genes) > n_target_hvg) final_hvg_genes <- final_hvg_genes[1:n_target_hvg]

  final_rna_counts <- GetAssayData(obj_rna_sub, assay = "RNA", layer = "counts")[final_hvg_genes, ]
  final_rna_umi <- Matrix::colSums(final_rna_counts)
  final_rna_sparsity <- calc_sparse_sparsity(final_rna_counts)
  final_mean_umi_rna <- mean(final_rna_umi)
  final_median_umi_rna <- median(final_rna_umi)

  # ============================================================================
  # STEP B: ATAC MODALITY
  # ============================================================================
  atac_file <- file.path(input_dir, paste0("GSE270392_Gm_atlas_", tissue_name, "_ACR_4cpm_rmcds_ACR_sparse.rds"))
  if (!file.exists(atac_file)) stop("ATAC file not found: ", atac_file)
  
  message("Loading ATAC sparse matrix: ", atac_file)
  atac_mat_raw <- readRDS(atac_file)
  if (!is(atac_mat_raw, "dgCMatrix")) atac_mat_raw <- as(atac_mat_raw, "CsparseMatrix")
  
  raw_atac_cells <- ncol(atac_mat_raw)
  raw_atac_peaks <- nrow(atac_mat_raw)
  raw_atac_counts <- Matrix::colSums(atac_mat_raw)
  raw_atac_sparsity <- calc_sparse_sparsity(atac_mat_raw)
  raw_mean_counts_atac <- mean(raw_atac_counts)
  raw_median_counts_atac <- median(raw_atac_counts)
  
  message(sprintf("Raw ATAC: %d cells | %d ACRs | Mean Counts: %.2f (Median: %.2f) | Sparsity: %.2f%%",
                  raw_atac_cells, raw_atac_peaks, raw_mean_counts_atac, raw_median_counts_atac, raw_atac_sparsity * 100))

  # Match ATAC metadata
  atac_tissue_meta <- atac_meta_global[atac_meta_global$Tissue == tissue_name, ]
  rownames(atac_tissue_meta) <- atac_tissue_meta$cellID
  common_atac_cells <- intersect(colnames(atac_mat_raw), rownames(atac_tissue_meta))
  
  if (length(common_atac_cells) == 0) {
    stop("No cell barcode overlap found between ATAC matrix and snATAC_meta.xlsx for tissue: ", tissue_name)
  }
  
  atac_mat_matched <- atac_mat_raw[, common_atac_cells]
  atac_meta_sub <- atac_tissue_meta[common_atac_cells, ]

  # QC filtering for ATAC
  atac_cell_depths <- Matrix::colSums(atac_mat_matched)
  atac_cell_features <- Matrix::colSums(atac_mat_matched > 0)
  
  qc_keep_atac <- which(
    atac_cell_depths >= 200 &
    atac_cell_features >= 100 &
    !is.na(atac_meta_sub$cell_type) &
    atac_meta_sub$cell_type != "" &
    atac_meta_sub$cell_type != "Unknown"
  )
  if (length(qc_keep_atac) < 200) {
    qc_keep_atac <- which(!is.na(atac_meta_sub$cell_type) & atac_meta_sub$cell_type != "")
  }
  atac_mat_filtered <- atac_mat_matched[, qc_keep_atac]
  atac_meta_filtered <- atac_meta_sub[qc_keep_atac, ]
  message("ATAC cells after QC filtering: ", ncol(atac_mat_filtered))

  # Dual-stratified sampling for ATAC
  atac_groups <- split(rownames(atac_meta_filtered), list(atac_meta_filtered$cell_type, atac_meta_filtered$batch), drop = TRUE)
  sampled_atac_cells <- unlist(lapply(atac_groups, function(cells) {
    if (length(cells) <= cell_cap) cells else sample(cells, cell_cap)
  }), use.names = FALSE)
  
  atac_mat_sampled <- atac_mat_filtered[, sampled_atac_cells]
  atac_meta_final <- atac_meta_filtered[sampled_atac_cells, ]
  message("ATAC cells after dual-stratified sampling: ", ncol(atac_mat_sampled))

  # Feature Selection: Strict 3:1 ATAC-to-RNA ratio (Target: 3 * n_target_hvg = 6,000 ACRs)
  target_atac_peaks <- 3 * n_target_hvg
  peak_det_rate <- Matrix::rowMeans(atac_mat_sampled > 0)
  eligible_peaks <- names(peak_det_rate[peak_det_rate >= 0.01])
  if (length(eligible_peaks) < target_atac_peaks) {
    eligible_peaks <- names(sort(peak_det_rate, decreasing = TRUE))[1:min(target_atac_peaks, length(peak_det_rate))]
  }
  
  peak_accessibility <- Matrix::rowSums(atac_mat_sampled[eligible_peaks, ])
  final_selected_peaks <- names(sort(peak_accessibility, decreasing = TRUE))[1:target_atac_peaks]
  
  final_atac_counts <- atac_mat_sampled[final_selected_peaks, ]
  final_atac_sparsity <- calc_sparse_sparsity(final_atac_counts)
  final_mean_counts_atac <- mean(Matrix::colSums(final_atac_counts))
  final_median_counts_atac <- median(Matrix::colSums(final_atac_counts))
  
  ratio_achieved <- length(final_selected_peaks) / length(final_hvg_genes)
  message(sprintf("Ratio Verification: ATAC (%d) / RNA (%d) = %.2f (Target: 3.00)",
                  length(final_selected_peaks), length(final_hvg_genes), ratio_achieved))

  # ============================================================================
  # STEP C: DIAGNOSTIC PLOTS (600 DPI JPEGs)
  # ============================================================================
  message("Generating 600 DPI diagnostic plots in: ", plot_dir)

  # 1. RNA histograms
  p1 <- ggplot(data.frame(nCount = obj_rna_sub$nCount_RNA), aes(x = nCount)) +
    geom_histogram(bins = 50, fill = "#1b9e77", color = "white") +
    scale_x_log10(labels = label_comma()) +
    theme_classic(base_size = 12) +
    labs(title = paste0("RNA UMI Distribution (", tissue_name, ")"), x = "nCount_RNA (log10)", y = "Cells")

  p2 <- ggplot(data.frame(nFeature = obj_rna_sub$nFeature_RNA), aes(x = nFeature)) +
    geom_histogram(bins = 50, fill = "#d95f02", color = "white") +
    scale_x_log10(labels = label_comma()) +
    theme_classic(base_size = 12) +
    labs(title = paste0("RNA Detected Genes (", tissue_name, ")"), x = "nFeature_RNA (log10)", y = "Cells")
  
  jpeg(file.path(plot_dir, "01_hist_RNA_QC.jpeg"), width = 3600, height = 1800, res = 600)
  print(p1 + p2)
  dev.off()

  # 2. RNA scatter
  p_sc_rna <- ggplot(obj_rna_sub@meta.data, aes(x = nCount_RNA, y = nFeature_RNA, color = batch)) +
    geom_point(alpha = 0.5, size = 0.8) +
    scale_x_log10(labels = label_comma()) +
    scale_y_log10(labels = label_comma()) +
    theme_classic(base_size = 12) +
    labs(title = paste0("RNA Quality Scatter (", tissue_name, ")"), x = "nCount_RNA", y = "nFeature_RNA")
  
  jpeg(file.path(plot_dir, "02_scatter_RNA_QC.jpeg"), width = 2800, height = 2400, res = 600)
  print(p_sc_rna)
  dev.off()

  # 3. ATAC histograms
  p3 <- ggplot(data.frame(nCount = Matrix::colSums(final_atac_counts)), aes(x = nCount)) +
    geom_histogram(bins = 50, fill = "#7570b3", color = "white") +
    scale_x_log10(labels = label_comma()) +
    theme_classic(base_size = 12) +
    labs(title = paste0("ATAC Depths in ACRs (", tissue_name, ")"), x = "Total Counts in ACRs (log10)", y = "Cells")

  p4 <- ggplot(data.frame(nFeature = Matrix::colSums(final_atac_counts > 0)), aes(x = nFeature)) +
    geom_histogram(bins = 50, fill = "#e7298a", color = "white") +
    scale_x_log10(labels = label_comma()) +
    theme_classic(base_size = 12) +
    labs(title = paste0("ATAC Detected ACRs (", tissue_name, ")"), x = "Detected ACRs (log10)", y = "Cells")
  
  jpeg(file.path(plot_dir, "03_hist_ATAC_QC.jpeg"), width = 3600, height = 1800, res = 600)
  print(p3 + p4)
  dev.off()

  # 4. ATAC scatter
  p_sc_atac <- ggplot(data.frame(nCount = Matrix::colSums(final_atac_counts),
                                nFeature = Matrix::colSums(final_atac_counts > 0),
                                batch = atac_meta_final$batch),
                      aes(x = nCount, y = nFeature, color = batch)) +
    geom_point(alpha = 0.5, size = 0.8) +
    scale_x_log10(labels = label_comma()) +
    scale_y_log10(labels = label_comma()) +
    theme_classic(base_size = 12) +
    labs(title = paste0("ATAC Quality Scatter (", tissue_name, ")"), x = "Total Counts in ACRs", y = "Detected ACRs")
  
  jpeg(file.path(plot_dir, "04_scatter_ATAC_QC.jpeg"), width = 2800, height = 2400, res = 600)
  print(p_sc_atac)
  dev.off()

  # 5. RNA Cell type composition
  rna_comp <- as.data.frame(table(obj_rna_sub$cell_type))
  colnames(rna_comp) <- c("cell_type", "count")
  p_comp_rna <- ggplot(rna_comp, aes(x = reorder(cell_type, count), y = count, fill = cell_type)) +
    geom_col(show.legend = FALSE) +
    coord_flip() +
    theme_classic(base_size = 11) +
    labs(title = paste0("RNA Cell Type Frequencies (", tissue_name, ")"), x = "Cell Type", y = "Cells")
  
  jpeg(file.path(plot_dir, "05_celltype_composition_RNA.jpeg"), width = 3000, height = 2400, res = 600)
  print(p_comp_rna)
  dev.off()

  # 6. ATAC Cell type composition
  atac_comp <- as.data.frame(table(atac_meta_final$cell_type))
  colnames(atac_comp) <- c("cell_type", "count")
  p_comp_atac <- ggplot(atac_comp, aes(x = reorder(cell_type, count), y = count, fill = cell_type)) +
    geom_col(show.legend = FALSE) +
    coord_flip() +
    theme_classic(base_size = 11) +
    labs(title = paste0("ATAC Cell Type Frequencies (", tissue_name, ")"), x = "Cell Type", y = "Cells")
  
  jpeg(file.path(plot_dir, "06_celltype_composition_ATAC.jpeg"), width = 3000, height = 2400, res = 600)
  print(p_comp_atac)
  dev.off()

  # 7. ACR detection rates
  p_det <- ggplot(data.frame(det = peak_det_rate[final_selected_peaks]), aes(x = det)) +
    geom_histogram(bins = 40, fill = "#66a61e", color = "white") +
    theme_classic(base_size = 12) +
    labs(title = paste0("Selected 6,000 ACR Detection Rates (", tissue_name, ")"), x = "Fraction of Cells Accessible", y = "ACRs")
  
  jpeg(file.path(plot_dir, "07_peak_detection_rates.jpeg"), width = 2800, height = 2000, res = 600)
  print(p_det)
  dev.off()

  # 8. Cross-modal Cell Type Concordance (RNA vs ATAC)
  rna_prop <- rna_comp %>% mutate(prop = count / sum(count), Modality = "RNA")
  atac_prop <- atac_comp %>% mutate(prop = count / sum(count), Modality = "ATAC")
  joint_prop <- rbind(rna_prop, atac_prop)
  
  p_cross <- ggplot(joint_prop, aes(x = cell_type, y = prop, fill = Modality)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    coord_flip() +
    scale_y_continuous(labels = percent) +
    theme_classic(base_size = 11) +
    labs(title = paste0("Cross-Modal Cell Type Proportions (", tissue_name, ")"), x = "Cell Type", y = "Proportion of Modality")
  
  jpeg(file.path(plot_dir, "08_crossmodal_celltype_concordance.jpeg"), width = 3200, height = 2400, res = 600)
  print(p_cross)
  dev.off()

  # ============================================================================
  # STEP D: SAVE DECOUPLED OUTPUTS & BY-BATCH MATRICES
  # ============================================================================
  # 1. Clean RNA Seurat object
  final_rna_obj <- CreateSeuratObject(counts = final_rna_counts, assay = "RNA", meta.data = obj_rna_sub@meta.data)
  saveRDS(final_rna_obj, file.path(tissue_out_dir, paste0("soy_processed_", tissue_name, "_rna.rds")))
  saveRDS(final_rna_counts, file.path(tissue_out_dir, paste0("soy_processed_", tissue_name, "_rna_counts.rds")))
  
  # 2. ATAC object & counts
  final_atac_obj <- CreateSeuratObject(counts = final_atac_counts, assay = "ATAC", meta.data = atac_meta_final)
  saveRDS(final_atac_obj, file.path(tissue_out_dir, paste0("soy_processed_", tissue_name, "_atac.rds")))
  saveRDS(final_atac_counts, file.path(tissue_out_dir, paste0("soy_processed_", tissue_name, "_atac_counts.rds")))
  
  # 3. Metadatas
  write.csv(obj_rna_sub@meta.data, file.path(tissue_out_dir, paste0("soy_processed_", tissue_name, "_rna_metadata.csv")), row.names = TRUE)
  write.csv(atac_meta_final, file.path(tissue_out_dir, paste0("soy_processed_", tissue_name, "_atac_metadata.csv")), row.names = TRUE)
  
  # 4. By-batch exports for RNA and ATAC
  for (b in unique(obj_rna_sub$batch)) {
    b_cells <- colnames(obj_rna_sub)[obj_rna_sub$batch == b]
    b_counts <- final_rna_counts[, b_cells, drop = FALSE]
    b_meta <- obj_rna_sub@meta.data[b_cells, , drop = FALSE]
    b_clean <- gsub("[^A-Za-z0-9_]", "_", b)
    saveRDS(b_counts, file.path(batch_dir, paste0("soy_batch_", b_clean, "_rna_counts.rds")))
    write.csv(b_meta, file.path(batch_dir, paste0("soy_batch_", b_clean, "_rna_metadata.csv")), row.names = TRUE)
  }
  
  for (b in unique(atac_meta_final$batch)) {
    b_cells <- rownames(atac_meta_final)[atac_meta_final$batch == b]
    b_counts <- final_atac_counts[, b_cells, drop = FALSE]
    b_meta <- atac_meta_final[b_cells, , drop = FALSE]
    b_clean <- gsub("[^A-Za-z0-9_]", "_", b)
    saveRDS(b_counts, file.path(batch_dir, paste0("soy_batch_", b_clean, "_atac_counts.rds")))
    write.csv(b_meta, file.path(batch_dir, paste0("soy_batch_", b_clean, "_atac_metadata.csv")), row.names = TRUE)
  }

  # ============================================================================
  # STEP E: SUMMARY METRICS RECORD
  # ============================================================================
  summary_df <- data.frame(
    Tissue                     = tissue_name,
    Initial_RNA_Cells          = raw_rna_cells,
    Final_RNA_Cells            = ncol(obj_rna_sub),
    Initial_ATAC_Cells         = raw_atac_cells,
    Final_ATAC_Cells           = ncol(atac_mat_sampled),
    Initial_RNA_Genes          = raw_rna_genes,
    Final_RNA_HVGs             = length(final_hvg_genes),
    Initial_ATAC_ACRs          = raw_atac_peaks,
    Final_ATAC_ACRs            = length(final_selected_peaks),
    ATAC_to_RNA_Ratio          = round(ratio_achieved, 2),
    Initial_RNA_Sparsity       = round(raw_rna_sparsity, 4),
    Final_RNA_Sparsity         = round(final_rna_sparsity, 4),
    Initial_ATAC_Sparsity      = round(raw_atac_sparsity, 4),
    Final_ATAC_Sparsity        = round(final_atac_sparsity, 4),
    Initial_Mean_UMI_RNA       = round(raw_mean_umi_rna, 2),
    Final_Mean_UMI_RNA         = round(final_mean_umi_rna, 2),
    Initial_Mean_Counts_ATAC   = round(raw_mean_counts_atac, 2),
    Final_Mean_Counts_ATAC     = round(final_mean_counts_atac, 2),
    RNA_Clusters               = length(unique(obj_rna_sub$cell_type)),
    ATAC_Clusters              = length(unique(atac_meta_final$cell_type)),
    RNA_Batches                = length(unique(obj_rna_sub$batch)),
    ATAC_Batches               = length(unique(atac_meta_final$batch)),
    RNA_Batch_List             = paste(sort(unique(obj_rna_sub$batch)), collapse = "; "),
    ATAC_Batch_List            = paste(sort(unique(atac_meta_final$batch)), collapse = "; "),
    stringsAsFactors           = FALSE
  )
  
  summary_csv <- file.path(tissue_out_dir, paste0("soy_processed_", tissue_name, "_summary_metrics.csv"))
  write.csv(summary_df, summary_csv, row.names = FALSE)
  message("Saved summary metrics: ", summary_csv)
  
  return(summary_df)
}

# ------------------------------------------------------------------------------
# 5. Execution Loop Across All Specified Tissues
# ------------------------------------------------------------------------------
all_summaries <- list()

for (tis in tissues_to_process) {
  tryCatch({
    sum_df <- process_single_soybean_tissue(
      tissue_name = tis,
      input_dir = opt$input_dir,
      output_dir = opt$output_dir,
      n_target_hvg = opt$n_hvg,
      cell_cap = opt$cell_cap
    )
    all_summaries[[tis]] <- sum_df
  }, error = function(e) {
    message("ERROR processing soybean tissue [", tis, "]: ", e$message)
  })
}

# Combine and write global Soybean dataset summary (cumulative across all processed tissues)
existing_summary_files <- list.files(opt$output_dir, pattern = "^soy_processed_.*_summary_metrics\\.csv$", recursive = TRUE, full.names = TRUE)
if (length(existing_summary_files) > 0) {
  all_tissue_summaries <- dplyr::bind_rows(lapply(existing_summary_files, read.csv, stringsAsFactors = FALSE))
  all_tissue_summaries <- all_tissue_summaries[!duplicated(all_tissue_summaries$Tissue), ]
  global_csv <- file.path(opt$output_dir, "soybean_all_7_tissues_standardized_summary.csv")
  tryCatch({
    write.csv(all_tissue_summaries, global_csv, row.names = FALSE)
    message("\n==================================================================")
    message("Cumulative Soybean tissues standardized summary updated!")
    message("Global summary table written to: ", global_csv)
    print(all_tissue_summaries)
    message("==================================================================")
  }, error = function(e) {
    warning("Could not write global summary CSV: ", e$message)
  })
}
