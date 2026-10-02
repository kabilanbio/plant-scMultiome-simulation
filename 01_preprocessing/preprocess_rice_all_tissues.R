# ==============================================================================
# preprocess_rice_all_tissues.R
# Standardized Single-Cell Multiomics Preprocessing Pipeline for Rice (8 Tissues)
# Supports: scMultiSim (strictly 3x ATAC:RNA ratio), scDesign3, MOSim, 
#           Matilda, scCross, and scMoMtF
# ==============================================================================

# Ensure all required packages are present; install automatically if missing
required_pkgs <- c("Seurat", "Signac", "Matrix", "ggplot2", "patchwork", "dplyr", "scales")
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
})

# ------------------------------------------------------------------------------
# 1. Native CLI Parameter Parser (Zero External Dependencies)
# ------------------------------------------------------------------------------
parse_cli_arguments <- function(args = commandArgs(trailingOnly = TRUE)) {
  # Auto-detect base directory if script is moved to external hard drive
  script_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", script_args, value = TRUE)
  this_dir <- if (length(file_arg) > 0) normalizePath(dirname(sub("^--file=", "", file_arg[1]))) else getwd()
  
  candidate_input <- normalizePath(file.path(this_dir, "..", "original_dataset", "rice_original_dataset", "rice_original_tissues"), mustWork = FALSE)
  candidate_output <- normalizePath(file.path(this_dir, "..", "featured_datasets_standardized"), mustWork = FALSE)
  
  default_input <- if (dir.exists(candidate_input)) candidate_input else "G:/PhD/sc_datasets/rice/original_dataset/rice_original_dataset/rice_original_tissues"
  default_output <- if (dir.exists(candidate_output)) candidate_output else "G:/PhD/sc_datasets/rice/featured_datasets_standardized"

  opt <- list(
    input_dir  = default_input,
    output_dir = default_output,
    tissue     = "remaining",
    n_hvg      = 2000,
    cell_cap   = 400,
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

# All 8 canonical rice tissues
ALL_TISSUES <- c("Bud", "Flag", "Leaf", "Root", "SAM", "Seed", "SP", "ST")

tissues_to_process <- if (tolower(opt$tissue) == "all") {
  ALL_TISSUES
} else if (tolower(opt$tissue) == "remaining") {
  done_tissues <- character(0)
  for (tis in ALL_TISSUES) {
    tis_metrics <- file.path(opt$output_dir, tis, paste0("rice_processed_", tis, "_summary_metrics.csv"))
    if (file.exists(tis_metrics)) done_tissues <- c(done_tissues, tis)
  }
  rem <- setdiff(ALL_TISSUES, done_tissues)
  if (length(rem) == 0) {
    message("All 8 tissues are already fully processed!")
    quit(save = "no", status = 0)
  }
  rem
} else {
  matched <- grep(opt$tissue, ALL_TISSUES, ignore.case = TRUE, value = TRUE)
  if (length(matched) == 0) stop("Invalid tissue name. Choose from: ", paste(ALL_TISSUES, collapse = ", "))
  matched
}

message("==================================================================")
message("Rice Single-Cell Multiomics Preprocessing Pipeline")
message("Tissues to process: ", paste(tissues_to_process, collapse = ", "))
message("Target HVGs: ", opt$n_hvg, " | Target ATAC peaks: ", 3 * opt$n_hvg, " (Strict 3:1 ratio)")
message("==================================================================")

# ------------------------------------------------------------------------------
# 2. Sparse Helper Functions
# ------------------------------------------------------------------------------
# Safe sparse row variance without memory duplication
calc_sparse_row_variance <- function(mat) {
  n <- ncol(mat)
  means <- Matrix::rowMeans(mat)
  mat_sq <- mat
  mat_sq@x <- mat_sq@x^2
  vars <- (Matrix::rowSums(mat_sq) - n * means^2) / (n - 1)
  names(vars) <- rownames(mat)
  return(vars)
}

# Safe sparse sparsity calculation
calc_sparse_sparsity <- function(mat) {
  total_elements <- as.numeric(nrow(mat)) * as.numeric(ncol(mat))
  1 - (Matrix::nnzero(mat) / total_elements)
}

# ------------------------------------------------------------------------------
# 3. Core Preprocessing Function per Tissue
# ------------------------------------------------------------------------------
process_single_rice_tissue <- function(tissue_name, input_dir, output_dir, n_target_hvg = 2000, cell_cap = 400) {
  message("\n>>> Processing Tissue: ", tissue_name, " <<<")
  
  # Setup directory tree
  tissue_out_dir <- file.path(output_dir, tissue_name)
  plot_dir       <- file.path(tissue_out_dir, "plots")
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Locate input file
  input_file <- file.path(input_dir, paste0("rice_original_", tissue_name, ".rds"))
  if (!file.exists(input_file)) {
    # Check lowercase fallback
    input_file_alt <- file.path(input_dir, paste0("rice_original_", tolower(tissue_name), ".rds"))
    if (file.exists(input_file_alt)) input_file <- input_file_alt
    else stop("File not found for tissue: ", input_file)
  }
  
  message("Loading input: ", input_file)
  obj <- UpdateSeuratObject(readRDS(input_file))
  
  # Verify required assays
  if (!"RNA" %in% names(obj@assays)) stop("RNA assay missing from Seurat object!")
  atac_assay_name <- if ("ATAC" %in% names(obj@assays)) "ATAC" else if ("peaks" %in% names(obj@assays)) "peaks" else stop("ATAC/peaks assay missing!")
  
  raw_cells <- ncol(obj)
  raw_genes <- nrow(obj[["RNA"]])
  raw_peaks <- nrow(obj[[atac_assay_name]])
  message(sprintf("Raw Dimensions: %d cells | %d genes | %d peaks", raw_cells, raw_genes, raw_peaks))
  
  meta <- obj@meta.data
  
  # Compute Initial Raw Library Depths and Sparsity
  raw_mean_umi_rna <- round(if ("nCount_RNA" %in% colnames(meta)) mean(meta$nCount_RNA) else mean(Matrix::colSums(GetAssayData(obj, assay = "RNA", layer = "counts"))), 2)
  raw_median_umi_rna <- round(if ("nCount_RNA" %in% colnames(meta)) median(meta$nCount_RNA) else median(Matrix::colSums(GetAssayData(obj, assay = "RNA", layer = "counts"))), 2)
  
  ncount_atac_col <- paste0("nCount_", atac_assay_name)
  raw_mean_counts_atac <- round(if (ncount_atac_col %in% colnames(meta)) mean(meta[[ncount_atac_col]]) else mean(Matrix::colSums(GetAssayData(obj, assay = atac_assay_name, layer = "counts"))), 2)
  raw_median_counts_atac <- round(if (ncount_atac_col %in% colnames(meta)) median(meta[[ncount_atac_col]]) else median(Matrix::colSums(GetAssayData(obj, assay = atac_assay_name, layer = "counts"))), 2)
  
  raw_rna_sparsity <- round(calc_sparse_sparsity(GetAssayData(obj, assay = "RNA", layer = "counts")), 4)
  raw_atac_sparsity <- round(calc_sparse_sparsity(GetAssayData(obj, assay = atac_assay_name, layer = "counts")), 4)
  
  message(sprintf("Raw Depths: Mean RNA UMI = %.2f (Median = %.2f) | Mean ATAC = %.2f (Median = %.2f)", 
                  raw_mean_umi_rna, raw_median_umi_rna, raw_mean_counts_atac, raw_median_counts_atac))
  message(sprintf("Raw Sparsity: RNA = %.2f%% | ATAC = %.2f%%", raw_rna_sparsity * 100, raw_atac_sparsity * 100))
  
  # ----------------------------------------------------------------------------
  # QC Plot Suite 1: Pre-filtering distributions
  # ----------------------------------------------------------------------------
  message("Generating initial QC distribution plots...")
  p_rna1 <- ggplot(meta, aes(x = nCount_RNA)) +
    geom_histogram(bins = 100, fill = "#1f77b4", alpha = 0.8) +
    scale_x_log10() +
    geom_vline(xintercept = c(800, 30000), color = "red", linetype = "dashed") +
    labs(title = paste0(tissue_name, " - nCount_RNA"), x = "UMI Count (log10)", y = "Cells") +
    theme_minimal(base_size = 12)

  p_rna2 <- ggplot(meta, aes(x = nFeature_RNA)) +
    geom_histogram(bins = 100, fill = "#ff7f0e", alpha = 0.8) +
    scale_x_log10() +
    geom_vline(xintercept = c(500, 5000), color = "red", linetype = "dashed") +
    labs(title = paste0(tissue_name, " - nFeature_RNA"), x = "Detected Genes (log10)", y = "Cells") +
    theme_minimal(base_size = 12)
  
  p_scatter_rna <- ggplot(meta, aes(x = nCount_RNA, y = nFeature_RNA)) +
    geom_point(size = 0.3, alpha = 0.25, color = "grey30") +
    scale_x_log10() + scale_y_log10() +
    labs(title = "RNA QC: UMI Count vs Detected Genes", x = "nCount_RNA", y = "nFeature_RNA") +
    theme_minimal(base_size = 12)
  
  ggsave(file.path(plot_dir, "01_hist_RNA_QC.jpeg"), plot = p_rna1 + p_rna2, width = 10, height = 4.5, dpi = 600)
  ggsave(file.path(plot_dir, "02_scatter_RNA_QC.jpeg"), plot = p_scatter_rna, width = 6, height = 5, dpi = 600)
  
  # ATAC QC plots (nFeature_ATAC / TSS enrichment if present)
  nfeat_atac_col <- paste0("nFeature_", atac_assay_name)
  if (nfeat_atac_col %in% colnames(meta)) {
    p_atac1 <- ggplot(meta, aes(x = .data[[nfeat_atac_col]])) +
      geom_histogram(bins = 100, fill = "#2ca02c", alpha = 0.8) +
      scale_x_log10() +
      geom_vline(xintercept = c(300, 10000), color = "red", linetype = "dashed") +
      labs(title = paste0(tissue_name, " - Peaks per cell"), x = "Peaks (log10)", y = "Cells") +
      theme_minimal(base_size = 12)
    
    if ("TSS.enrichment" %in% colnames(meta)) {
      p_atac2 <- ggplot(meta, aes(x = TSS.enrichment)) +
        geom_histogram(bins = 100, fill = "#9467bd", alpha = 0.8) +
        geom_vline(xintercept = c(1.5, 4.0), color = "red", linetype = "dashed") +
        labs(title = "TSS Enrichment", x = "Enrichment Score", y = "Cells") +
        theme_minimal(base_size = 12)
      
      p_scatter_atac <- ggplot(meta, aes(x = .data[[nfeat_atac_col]], y = TSS.enrichment)) +
        geom_point(size = 0.3, alpha = 0.25, color = "#2ca02c") +
        scale_x_log10() +
        labs(title = "ATAC QC: Peaks vs TSS Enrichment", x = "Peaks per cell", y = "TSS Enrichment") +
        theme_minimal(base_size = 12)
      
      ggsave(file.path(plot_dir, "03_hist_ATAC_QC.jpeg"), plot = p_atac1 + p_atac2, width = 10, height = 4.5, dpi = 600)
      ggsave(file.path(plot_dir, "04_scatter_ATAC_QC.jpeg"), plot = p_scatter_atac, width = 6, height = 5, dpi = 600)
    } else {
      ggsave(file.path(plot_dir, "03_hist_ATAC_QC.jpeg"), plot = p_atac1, width = 6, height = 4.5, dpi = 600)
    }
  }

  # ----------------------------------------------------------------------------
  # QC Step 1: Organelle (Mitochondrial & Chloroplast) Read Filtering
  # ----------------------------------------------------------------------------
  rna_features <- rownames(obj[["RNA"]])
  # Rice nuclear genes begin with "Os" (or LOC_Os). Organelle genes are non-Os or OsMT/OsCP
  organelle_genes <- grep("^OsMT|^OSMT|^Mt|^mt\\.|^rrn|^nad|^cox|^atp|^ccb|^OsCP|^OSCP|^Pt|^pt\\.|^psa|^psb|^rbc|^rpo", 
                          rna_features, value = TRUE, ignore.case = TRUE)
  if (length(organelle_genes) < 5) {
    # Fallback to all non-Os features if explicit MT/CP prefixes were stripped
    organelle_genes <- grep("^Os", rna_features, invert = TRUE, value = TRUE)
  }
  
  obj[["percent.organelle"]] <- PercentageFeatureSet(obj, features = organelle_genes, assay = "RNA")
  
  p_org <- ggplot(obj@meta.data, aes(x = percent.organelle)) +
    geom_histogram(bins = 80, fill = "#d62728", alpha = 0.8) +
    geom_vline(xintercept = 5, color = "black", linetype = "dashed") +
    labs(title = paste0(tissue_name, " - Organelle Gene Read %"), x = "Percent Organelle Reads (%)", y = "Cells") +
    theme_minimal(base_size = 12)
  ggsave(file.path(plot_dir, "05_organelle_percent.jpeg"), plot = p_org, width = 6, height = 4.5, dpi = 600)

  # ----------------------------------------------------------------------------
  # QC Step 2: Quality Filtering Thresholds
  # ----------------------------------------------------------------------------
  message("Applying strict QC filtering thresholds...")
  # RNA thresholds
  obj <- subset(obj, nCount_RNA >= 800 & nCount_RNA <= 30000)
  obj <- subset(obj, nFeature_RNA >= 500 & nFeature_RNA <= 5000)
  obj <- subset(obj, percent.organelle < 5)
  
  # ATAC thresholds (using cells = ... to prevent NSE lookup errors)
  if (nfeat_atac_col %in% colnames(obj@meta.data)) {
    valid_atac_cells <- colnames(obj)[obj@meta.data[[nfeat_atac_col]] >= 300 & obj@meta.data[[nfeat_atac_col]] <= 10000]
    obj <- subset(obj, cells = valid_atac_cells)
  }
  if ("TSS.enrichment" %in% colnames(obj@meta.data)) {
    valid_tss_cells <- colnames(obj)[obj$TSS.enrichment >= 1.5 & obj$TSS.enrichment <= 4.0]
    obj <- subset(obj, cells = valid_tss_cells)
  }
  
  # Standardize cell type and batch columns
  cluster_col <- intersect(c("cluster_names", "tissue_cluster_names", "cluster", "cell_type"), colnames(obj@meta.data))[1]
  batch_col   <- intersect(c("sample", "batch", "orig.ident"), colnames(obj@meta.data))[1]
  
  obj$cell_type <- if (!is.na(cluster_col)) as.character(obj@meta.data[[cluster_col]]) else "Cluster_1"
  obj$batch     <- if (!is.na(batch_col)) as.character(obj@meta.data[[batch_col]]) else "Batch_1"
  
  # Remove unannotated / unknown clusters if present
  unknown_mask <- grepl("unknown|unassigned", obj$cell_type, ignore.case = TRUE)
  if (any(unknown_mask)) {
    message("Removing ", sum(unknown_mask), " unknown/unassigned cells...")
    obj <- obj[, !unknown_mask]
  }
  message("Cells remaining after QC filtering: ", ncol(obj))

  # ----------------------------------------------------------------------------
  # QC Step 3: Dual-Stratified Subsampling (Balanced Cell Type & Batch Sampling)
  # ----------------------------------------------------------------------------
  n_batches <- length(unique(obj$batch))
  per_batch_cap <- max(10, floor(cell_cap / n_batches))
  message(sprintf("Dual-stratified sampling: Capping at %d cells per (CellType x Batch) [%d batches]", 
                  per_batch_cap, n_batches))
  
  # Pre-subsampling composition by cell type and batch
  df_pre <- as.data.frame(table(CellType = obj$cell_type, Batch = obj$batch))
  
  cell_batch_pairs <- paste(obj$cell_type, obj$batch, sep = "___")
  keep_cells <- unlist(lapply(unique(cell_batch_pairs), function(pair) {
    cells_in_pair <- colnames(obj)[cell_batch_pairs == pair]
    if (length(cells_in_pair) <= per_batch_cap) {
      return(cells_in_pair)
    } else {
      return(sample(cells_in_pair, size = per_batch_cap))
    }
  }))
  
  obj_sub <- subset(obj, cells = keep_cells)
  
  # Post-subsampling composition plots
  df_post <- as.data.frame(table(CellType = obj_sub$cell_type))
  df_pre_ct <- as.data.frame(table(CellType = obj$cell_type))
  df_comp <- merge(df_pre_ct, df_post, by = "CellType", suffixes = c("_Before", "_After"))
  df_comp_long <- reshape2::melt(df_comp, id.vars = "CellType", 
                                 variable.name = "Stage", value.name = "CellCount")
  
  p_comp <- ggplot(df_comp_long, aes(x = CellType, y = CellCount, fill = Stage)) +
    geom_bar(stat = "identity", position = "dodge", alpha = 0.85) +
    scale_fill_manual(values = c("Freq_Before" = "#7f7f7f", "Freq_After" = "#1f77b4")) +
    coord_flip() +
    labs(title = paste0(tissue_name, " - Cell Type Representation (Cap = ", cell_cap, ")"),
         x = "Cell Type / Cluster", y = "Number of Cells") +
    theme_minimal(base_size = 11)
  ggsave(file.path(plot_dir, "06_celltype_composition.jpeg"), plot = p_comp, width = 8, height = max(4, nrow(df_pre_ct) * 0.35), dpi = 600)
  
  # Batch representation plot per cell type
  df_batch <- as.data.frame(table(CellType = obj_sub$cell_type, Batch = obj_sub$batch))
  p_batch <- ggplot(df_batch, aes(x = CellType, y = Freq, fill = Batch)) +
    geom_bar(stat = "identity", position = "stack", alpha = 0.85) +
    scale_fill_brewer(palette = "Set2") +
    coord_flip() +
    labs(title = paste0(tissue_name, " - Batch Representation across Cell Types"),
         x = "Cell Type", y = "Cells", fill = "Batch") +
    theme_minimal(base_size = 11)
  ggsave(file.path(plot_dir, "08_batch_composition.jpeg"), plot = p_batch, width = 8, height = max(4, length(unique(obj_sub$cell_type)) * 0.35), dpi = 600)
  
  message("Cells after dual-stratified sampling: ", ncol(obj_sub))
  message("Batches detected: ", length(unique(obj_sub$batch)), " (", paste(sort(unique(obj_sub$batch)), collapse = ", "), ")")

  # ----------------------------------------------------------------------------
  # Feature Selection: Batch-Aware RNA Highly Variable Genes (HVGs)
  # ----------------------------------------------------------------------------
  message("Selecting top Highly Variable Genes (Batch-Aware Consensus HVGs)...")
  DefaultAssay(obj_sub) <- "RNA"
  obj_sub <- NormalizeData(obj_sub, verbose = FALSE)
  
  # Use split.by = 'batch' if multiple batches exist to isolate true biological markers from technical batch noise
  if (length(unique(obj_sub$batch)) > 1) {
    obj_sub <- FindVariableFeatures(obj_sub, selection.method = "vst", split.by = "batch", 
                                    nfeatures = n_target_hvg * 1.5, verbose = FALSE)
  } else {
    obj_sub <- FindVariableFeatures(obj_sub, selection.method = "vst", 
                                    nfeatures = n_target_hvg * 1.5, verbose = FALSE)
  }
  
  cand_hvgs <- VariableFeatures(obj_sub)
  
  # Remove all-zero genes in the subsampled cohort
  rna_counts <- GetAssayData(obj_sub, assay = "RNA", layer = "counts")[cand_hvgs, ]
  non_zero_hvgs <- cand_hvgs[Matrix::rowSums(rna_counts) > 0]
  
  # Strictly exclude organelle genes from simulated HVGs
  clean_hvgs <- non_zero_hvgs[!non_zero_hvgs %in% organelle_genes]
  
  # Finalize to exactly n_target_hvg (or maximum available clean HVGs)
  final_n_rna <- min(n_target_hvg, length(clean_hvgs))
  selected_hvgs <- clean_hvgs[1:final_n_rna]
  
  # Subset RNA assay
  obj_sub[["RNA"]] <- subset(obj_sub[["RNA"]], features = selected_hvgs)
  VariableFeatures(obj_sub, assay = "RNA") <- selected_hvgs
  message("Final Selected Batch-Aware RNA HVGs: ", length(selected_hvgs))

  # ----------------------------------------------------------------------------
  # Feature Selection: ATAC Peaks (STRICT 3x RNA RATIO FOR scMultiSim)
  # ----------------------------------------------------------------------------
  target_atac_peaks <- 3 * length(selected_hvgs)
  message(sprintf("Target ATAC Peaks for scMultiSim (3 * %d): %d peaks", length(selected_hvgs), target_atac_peaks))
  
  DefaultAssay(obj_sub) <- atac_assay_name
  atac_counts_all <- GetAssayData(obj_sub, assay = atac_assay_name, layer = "counts")
  
  # 1. Peak detection rate filter (accessible in >= 1% of cells)
  peak_detection_rate <- Matrix::rowMeans(atac_counts_all > 0)
  
  # Diagnostic plot of peak detection rates
  jpeg(file.path(plot_dir, "07_peak_detection_rates.jpeg"), width = 7, height = 4.5, units = "in", res = 600)
  hist(log10(peak_detection_rate + 1e-5), breaks = 100, col = "#2ca02c", border = "white",
       main = paste0(tissue_name, " - Peak Detection Rates"), xlab = "log10(Fraction of Cells)")
  abline(v = log10(0.01), col = "red", lty = 2, lwd = 1.5)
  dev.off()
  
  accessible_peaks <- names(peak_detection_rate[peak_detection_rate >= 0.01])
  message("Peaks passing 1% detection rate: ", length(accessible_peaks))
  
  if (length(accessible_peaks) < target_atac_peaks) {
    warning("Fewer than target peaks pass 1% filter (", length(accessible_peaks), 
            " < ", target_atac_peaks, "). Relaxing detection threshold.")
    accessible_peaks <- names(sort(peak_detection_rate, decreasing = TRUE)[1:min(target_atac_peaks * 2, length(peak_detection_rate))])
  }
  
  # Subset ATAC assay to accessible peaks
  obj_sub[[atac_assay_name]] <- subset(obj_sub[[atac_assay_name]], features = accessible_peaks)
  
  # 2. Run TF-IDF to find top variable peaks
  obj_sub <- RunTFIDF(obj_sub, assay = atac_assay_name, verbose = FALSE)
  tfidf_data <- GetAssayData(obj_sub, assay = atac_assay_name, layer = "data")
  
  # Compute sparse row variance of TF-IDF
  peak_variances <- calc_sparse_row_variance(tfidf_data)
  
  # Rank and take top 3x peaks
  ranked_peaks <- names(sort(peak_variances, decreasing = TRUE))
  selected_peaks <- ranked_peaks[1:target_atac_peaks]
  
  obj_sub[[atac_assay_name]] <- subset(obj_sub[[atac_assay_name]], features = selected_peaks)
  VariableFeatures(obj_sub, assay = atac_assay_name) <- selected_peaks
  message("Final Selected ATAC Peaks: ", length(selected_peaks))

  # ----------------------------------------------------------------------------
  # Verification & Sparsity Audit
  # ----------------------------------------------------------------------------
  final_rna_mat  <- GetAssayData(obj_sub, assay = "RNA", layer = "counts")
  final_atac_mat <- GetAssayData(obj_sub, assay = atac_assay_name, layer = "counts")
  
  rna_sparsity  <- calc_sparse_sparsity(final_rna_mat)
  atac_sparsity <- calc_sparse_sparsity(final_atac_mat)
  
  ratio_check <- nrow(final_atac_mat) / nrow(final_rna_mat)
  message(sprintf("Ratio Verification: ATAC (%d) / RNA (%d) = %.2f (Target: 3.00)", 
                  nrow(final_atac_mat), nrow(final_rna_mat), ratio_check))
  message(sprintf("RNA Sparsity: %.4f | ATAC Sparsity: %.4f", rna_sparsity, atac_sparsity))

  # ----------------------------------------------------------------------------
  # Save RNA and ATAC Separately
  # ----------------------------------------------------------------------------
  # 1. Create separate unimodal Seurat objects
  clean_meta <- obj_sub@meta.data
  
  obj_rna <- CreateSeuratObject(
    counts = final_rna_mat,
    assay = "RNA",
    meta.data = clean_meta
  )
  VariableFeatures(obj_rna) <- selected_hvgs
  
  # For ATAC, create ChromatinAssay or Seurat object
  atac_assay_clean <- CreateChromatinAssay(
    counts = final_atac_mat,
    sep = c("-", "-"),
    genome = NULL
  )
  obj_atac <- CreateSeuratObject(
    counts = atac_assay_clean,
    assay = "ATAC",
    meta.data = clean_meta
  )
  VariableFeatures(obj_atac) <- selected_peaks

  # 2. File paths for separate outputs
  out_rna_rds         <- file.path(tissue_out_dir, paste0("rice_processed_", tissue_name, "_rna.rds"))
  out_atac_rds        <- file.path(tissue_out_dir, paste0("rice_processed_", tissue_name, "_atac.rds"))
  out_rna_counts_rds  <- file.path(tissue_out_dir, paste0("rice_processed_", tissue_name, "_rna_counts.rds"))
  out_atac_counts_rds <- file.path(tissue_out_dir, paste0("rice_processed_", tissue_name, "_atac_counts.rds"))
  out_metadata_csv    <- file.path(tissue_out_dir, paste0("rice_processed_", tissue_name, "_metadata.csv"))

  # Save separate Seurat objects
  saveRDS(obj_rna,  out_rna_rds)
  saveRDS(obj_atac, out_atac_rds)
  message("Saved separate RNA object:  ", out_rna_rds)
  message("Saved separate ATAC object: ", out_atac_rds)

  # Save separate raw count matrices (dgCMatrix) for fast simulator ingestion
  saveRDS(final_rna_mat,  out_rna_counts_rds)
  saveRDS(final_atac_mat, out_atac_counts_rds)
  message("Saved RNA count matrix:     ", out_rna_counts_rds)
  message("Saved ATAC count matrix:    ", out_atac_counts_rds)

  # Save matched cell metadata
  write.csv(clean_meta, out_metadata_csv, row.names = TRUE)
  message("Saved matched cell metadata:", out_metadata_csv)

  # Save dedicated Cell Types and Batches annotation files
  out_celltypes_csv <- file.path(tissue_out_dir, paste0("rice_processed_", tissue_name, "_cell_types.csv"))
  out_batches_csv   <- file.path(tissue_out_dir, paste0("rice_processed_", tissue_name, "_batches.csv"))
  
  df_celltypes <- data.frame(CellBarcode = colnames(obj_sub), CellType = obj_sub$cell_type, stringsAsFactors = FALSE)
  df_batches   <- data.frame(CellBarcode = colnames(obj_sub), Batch = obj_sub$batch, stringsAsFactors = FALSE)
  
  write.csv(df_celltypes, out_celltypes_csv, row.names = FALSE)
  write.csv(df_batches,   out_batches_csv,   row.names = FALSE)
  message("Saved dedicated cell types: ", out_celltypes_csv)
  message("Saved dedicated batches:    ", out_batches_csv)

  # Save individual per-batch subsets (for batch integration tests, single-batch simulation)
  batch_out_dir <- file.path(tissue_out_dir, "by_batch")
  dir.create(batch_out_dir, recursive = TRUE, showWarnings = FALSE)
  for (b in unique(obj_sub$batch)) {
    b_cells <- colnames(obj_sub)[obj_sub$batch == b]
    b_safe  <- gsub("[^A-Za-z0-9_]", "_", b)
    
    saveRDS(final_rna_mat[, b_cells, drop = FALSE], 
            file.path(batch_out_dir, paste0(tissue_name, "_batch_", b_safe, "_rna_counts.rds")))
    saveRDS(final_atac_mat[, b_cells, drop = FALSE], 
            file.path(batch_out_dir, paste0(tissue_name, "_batch_", b_safe, "_atac_counts.rds")))
    write.csv(clean_meta[b_cells, , drop = FALSE], 
              file.path(batch_out_dir, paste0(tissue_name, "_batch_", b_safe, "_metadata.csv")), row.names = TRUE)
  }
  message("Saved individual per-batch files in: ", batch_out_dir)
  
  # Save Tissue Summary Metadata CSV
  summary_df <- data.frame(
    Tissue                   = tissue_name,
    Initial_Cells            = raw_cells,
    Final_Cells              = ncol(obj_sub),
    Initial_RNA_Genes        = raw_genes,
    Final_RNA_HVGs           = nrow(final_rna_mat),
    Initial_ATAC_Peaks       = raw_peaks,
    Final_ATAC_Peaks         = nrow(final_atac_mat),
    ATAC_to_RNA_Ratio        = round(ratio_check, 2),
    Initial_RNA_Sparsity     = raw_rna_sparsity,
    Final_RNA_Sparsity       = round(rna_sparsity, 4),
    Initial_ATAC_Sparsity    = raw_atac_sparsity,
    Final_ATAC_Sparsity      = round(atac_sparsity, 4),
    Initial_Mean_UMI_RNA     = raw_mean_umi_rna,
    Final_Mean_UMI_RNA       = round(mean(Matrix::colSums(final_rna_mat)), 2),
    Initial_Mean_Counts_ATAC = raw_mean_counts_atac,
    Final_Mean_Counts_ATAC   = round(mean(Matrix::colSums(final_atac_mat)), 2),
    Clusters                 = length(unique(obj_sub$cell_type)),
    Batches                  = length(unique(obj_sub$batch)),
    Batch_List               = paste(sort(unique(obj_sub$batch)), collapse = "; "),
    stringsAsFactors         = FALSE
  )
  
  summary_csv <- file.path(tissue_out_dir, paste0("rice_processed_", tissue_name, "_summary_metrics.csv"))
  tryCatch({
    write.csv(summary_df, summary_csv, row.names = FALSE)
    message("Saved summary metrics:      ", summary_csv)
  }, error = function(e) {
    warning("Could not overwrite ", summary_csv, " (file may be open in Excel): ", e$message)
  })
  
  return(summary_df)
}

# ------------------------------------------------------------------------------
# 4. Execution Loop Across All Specified Tissues
# ------------------------------------------------------------------------------
all_summaries <- list()

for (tis in tissues_to_process) {
  tryCatch({
    sum_df <- process_single_rice_tissue(
      tissue_name = tis,
      input_dir = opt$input_dir,
      output_dir = opt$output_dir,
      n_target_hvg = opt$n_hvg,
      cell_cap = opt$cell_cap
    )
    all_summaries[[tis]] <- sum_df
  }, error = function(e) {
    message("ERROR processing tissue [", tis, "]: ", e$message)
  })
}

# Combine and write global Rice dataset summary (cumulative across all processed tissues)
existing_summary_files <- list.files(opt$output_dir, pattern = "^rice_processed_.*_summary_metrics\\.csv$", recursive = TRUE, full.names = TRUE)
if (length(existing_summary_files) > 0) {
  all_tissue_summaries <- dplyr::bind_rows(lapply(existing_summary_files, read.csv, stringsAsFactors = FALSE))
  all_tissue_summaries <- all_tissue_summaries[!duplicated(all_tissue_summaries$Tissue), ]
  global_csv <- file.path(opt$output_dir, "rice_all_8_tissues_standardized_summary.csv")
  tryCatch({
    write.csv(all_tissue_summaries, global_csv, row.names = FALSE)
    message("\n==================================================================")
    message("Cumulative Rice tissues standardized summary updated!")
    message("Global summary table written to: ", global_csv)
    print(all_tissue_summaries)
    message("==================================================================")
  }, error = function(e) {
    warning("Could not write global summary CSV (file may be open in Excel): ", e$message)
  })
}
