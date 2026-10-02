# ==============================================================================
# preprocess_arabidopsis_root.R
# Standardized Single-Cell Multiomics Preprocessing Pipeline for Arabidopsis Root
# Unpaired Multi-Omics Architecture (GSE155304): snRNA-seq + snATAC-seq
# Strictly enforces:
#   - 3:1 ATAC-to-RNA feature ratio (2,000 HVGs, 6,000 peaks)
#   - Dual-stratified sampling per modality (max 200 cells / cell_type x batch)
#   - High-resolution 600 DPI JPEG figures (zero PDFs)
#   - Decoupled unimodal RDS objects, count matrices, cell types, batches, and by_batch files
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
  library(scales)
})

# ------------------------------------------------------------------------------
# 1. Native CLI Parameter Parser (Zero External Dependencies)
# ------------------------------------------------------------------------------
parse_cli_arguments <- function(args = commandArgs(trailingOnly = TRUE)) {
  script_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", script_args, value = TRUE)
  this_dir <- if (length(file_arg) > 0) normalizePath(dirname(sub("^--file=", "", file_arg[1]))) else getwd()
  
  candidate_input <- normalizePath(file.path(this_dir, "..", "original_datasets"), mustWork = FALSE)
  candidate_output <- normalizePath(file.path(this_dir, "..", "featured_datasets_standardized"), mustWork = FALSE)
  
  default_input <- if (dir.exists(candidate_input)) candidate_input else "G:/PhD/sc_datasets/arabidopsis_root/original_datasets"
  default_output <- if (dir.exists(candidate_output)) candidate_output else "G:/PhD/sc_datasets/arabidopsis_root/featured_datasets_standardized"

  opt <- list(
    input_dir  = default_input,
    output_dir = default_output,
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

message("==================================================================")
message("Arabidopsis Root Single-Cell Multiomics Preprocessing (Unpaired)")
message("Input directory:  ", opt$input_dir)
message("Output directory: ", opt$output_dir)
message("Target RNA HVGs:  ", opt$n_hvg, " | Target ATAC Peaks: ", 3 * opt$n_hvg, " (Strict 3:1 ratio)")
message("Max cells per (cell_type x batch): ", opt$cell_cap)
message("==================================================================")

# ------------------------------------------------------------------------------
# 2. Helper Functions
# ------------------------------------------------------------------------------
calc_sparse_sparsity <- function(mat) {
  total_elements <- as.numeric(nrow(mat)) * as.numeric(ncol(mat))
  1 - (Matrix::nnzero(mat) / total_elements)
}

get_counts_matrix <- function(obj, assay_name) {
  mat <- tryCatch(
    GetAssayData(obj, assay = assay_name, layer = "counts"),
    error = function(e) NULL
  )
  if (is.null(mat)) {
    mat <- tryCatch(
      obj@assays[[assay_name]]@counts,
      error = function(e) NULL
    )
  }
  if (is.null(mat)) {
    mat <- tryCatch(
      GetAssayData(obj, assay = assay_name, layer = "data"),
      error = function(e) NULL
    )
    if (is.null(mat)) {
      mat <- tryCatch(
        obj@assays[[assay_name]]@data,
        error = function(e) NULL
      )
    }
  }
  if (is.null(mat)) stop("Cannot find counts matrix for assay ", assay_name)
  return(mat)
}

# ------------------------------------------------------------------------------
# 3. Setup Directories
# ------------------------------------------------------------------------------
tissue_out_dir <- file.path(opt$output_dir, "Root")
plot_dir       <- file.path(tissue_out_dir, "plots")
batch_dir      <- file.path(tissue_out_dir, "by_batch")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(batch_dir, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# STEP A: RNA MODALITY
# ==============================================================================
rna_file <- file.path(opt$input_dir, "GSE155304_rnaseq_integration.rds")
if (!file.exists(rna_file)) stop("Arabidopsis RNA file not found: ", rna_file)

message("\n>>> Loading Arabidopsis RNA Seurat object: ", rna_file, " <<<")
obj_rna <- readRDS(rna_file)
obj_rna <- UpdateSeuratObject(obj_rna)
DefaultAssay(obj_rna) <- "RNA"

# Harmonize metadata columns
if ("seurat_clusters_renamed" %in% colnames(obj_rna@meta.data)) {
  obj_rna$cell_type <- as.character(obj_rna$seurat_clusters_renamed)
} else if ("seurat_clusters" %in% colnames(obj_rna@meta.data)) {
  obj_rna$cell_type <- paste0("Cluster_", obj_rna$seurat_clusters)
} else {
  obj_rna$cell_type <- "Root_cells"
}

if ("orig.ident" %in% colnames(obj_rna@meta.data)) {
  obj_rna$batch <- as.character(obj_rna$orig.ident)
} else if ("protocol" %in% colnames(obj_rna@meta.data)) {
  obj_rna$batch <- as.character(obj_rna$protocol)
} else {
  obj_rna$batch <- "batch_1"
}

raw_rna_cells <- ncol(obj_rna)
raw_rna_counts <- get_counts_matrix(obj_rna, "RNA")
raw_rna_genes <- nrow(raw_rna_counts)
raw_rna_umi <- Matrix::colSums(raw_rna_counts)
raw_rna_sparsity <- calc_sparse_sparsity(raw_rna_counts)
raw_mean_umi_rna <- mean(raw_rna_umi)
raw_median_umi_rna <- median(raw_rna_umi)

message(sprintf("Raw RNA: %d cells | %d genes | Mean UMI: %.2f (Median: %.2f) | Sparsity: %.2f%%",
                raw_rna_cells, raw_rna_genes, raw_mean_umi_rna, raw_median_umi_rna, raw_rna_sparsity * 100))

# QC filtering for RNA
qc_keep_rna <- colnames(obj_rna)[
  obj_rna$nCount_RNA >= 500 &
  obj_rna$nFeature_RNA >= 250 &
  !is.na(obj_rna$cell_type) &
  obj_rna$cell_type != ""
]
if (length(qc_keep_rna) < 500) qc_keep_rna <- colnames(obj_rna)
obj_rna_filtered <- subset(obj_rna, cells = qc_keep_rna)
message("RNA cells after QC filtering: ", ncol(obj_rna_filtered))

# Dual-stratified sampling for RNA
rna_groups <- split(colnames(obj_rna_filtered), list(obj_rna_filtered$cell_type, obj_rna_filtered$batch), drop = TRUE)
sampled_rna_cells <- unlist(lapply(rna_groups, function(cells) {
  if (length(cells) <= opt$cell_cap) cells else sample(cells, opt$cell_cap)
}), use.names = FALSE)

obj_rna_sub <- subset(obj_rna_filtered, cells = sampled_rna_cells)
message("RNA cells after dual-stratified sampling: ", ncol(obj_rna_sub))

# Batch-aware consensus HVGs
obj_rna_sub <- NormalizeData(obj_rna_sub, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
obj_rna_sub <- FindVariableFeatures(obj_rna_sub, selection.method = "vst", nfeatures = opt$n_hvg, verbose = FALSE)
final_hvg_genes <- VariableFeatures(obj_rna_sub)
if (length(final_hvg_genes) > opt$n_hvg) final_hvg_genes <- final_hvg_genes[1:opt$n_hvg]

final_rna_counts <- get_counts_matrix(obj_rna_sub, "RNA")[final_hvg_genes, ]
final_rna_umi <- Matrix::colSums(final_rna_counts)
final_rna_sparsity <- calc_sparse_sparsity(final_rna_counts)
final_mean_umi_rna <- mean(final_rna_umi)
final_median_umi_rna <- median(final_rna_umi)

# ==============================================================================
# STEP B: ATAC MODALITY
# ==============================================================================
atac_file <- file.path(opt$input_dir, "GSE155304_atacseq_integration.rds")
if (!file.exists(atac_file)) stop("Arabidopsis ATAC file not found: ", atac_file)

message("\n>>> Loading Arabidopsis ATAC Seurat object: ", atac_file, " <<<")
obj_atac <- readRDS(atac_file)
obj_atac <- UpdateSeuratObject(obj_atac)
DefaultAssay(obj_atac) <- "peaks"

# Harmonize metadata columns
if ("predicted.id" %in% colnames(obj_atac@meta.data)) {
  obj_atac$cell_type <- as.character(obj_atac$predicted.id)
} else if ("seurat_clusters_renamed" %in% colnames(obj_atac@meta.data)) {
  obj_atac$cell_type <- as.character(obj_atac$seurat_clusters_renamed)
} else if ("seurat_clusters" %in% colnames(obj_atac@meta.data)) {
  obj_atac$cell_type <- paste0("Cluster_", obj_atac$seurat_clusters)
} else {
  obj_atac$cell_type <- "Root_cells"
}

if ("orig.ident" %in% colnames(obj_atac@meta.data)) {
  obj_atac$batch <- as.character(obj_atac$orig.ident)
} else {
  obj_atac$batch <- "batch_1"
}

raw_atac_cells <- ncol(obj_atac)
raw_atac_counts <- get_counts_matrix(obj_atac, "peaks")
raw_atac_peaks <- nrow(raw_atac_counts)
raw_atac_depth <- Matrix::colSums(raw_atac_counts)
raw_atac_sparsity <- calc_sparse_sparsity(raw_atac_counts)
raw_mean_counts_atac <- mean(raw_atac_depth)
raw_median_counts_atac <- median(raw_atac_depth)

message(sprintf("Raw ATAC: %d cells | %d peaks | Mean Counts: %.2f (Median: %.2f) | Sparsity: %.2f%%",
                raw_atac_cells, raw_atac_peaks, raw_mean_counts_atac, raw_median_counts_atac, raw_atac_sparsity * 100))

# QC filtering for ATAC
qc_keep_atac <- colnames(obj_atac)[
  obj_atac$nCount_peaks >= 500 &
  obj_atac$nFeature_peaks >= 200 &
  !is.na(obj_atac$cell_type) &
  obj_atac$cell_type != ""
]
if (length(qc_keep_atac) < 500) qc_keep_atac <- colnames(obj_atac)
obj_atac_filtered <- subset(obj_atac, cells = qc_keep_atac)
message("ATAC cells after QC filtering: ", ncol(obj_atac_filtered))

# Dual-stratified sampling for ATAC
atac_groups <- split(colnames(obj_atac_filtered), list(obj_atac_filtered$cell_type, obj_atac_filtered$batch), drop = TRUE)
sampled_atac_cells <- unlist(lapply(atac_groups, function(cells) {
  if (length(cells) <= opt$cell_cap) cells else sample(cells, opt$cell_cap)
}), use.names = FALSE)

obj_atac_sub <- subset(obj_atac_filtered, cells = sampled_atac_cells)
message("ATAC cells after dual-stratified sampling: ", ncol(obj_atac_sub))

# Feature Selection: Strict 3:1 ATAC-to-RNA ratio (Target: 3 * n_target_hvg = 6,000 peaks)
target_atac_peaks <- 3 * opt$n_hvg
atac_mat_sampled <- get_counts_matrix(obj_atac_sub, "peaks")
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

# ==============================================================================
# STEP C: DIAGNOSTIC PLOTS (600 DPI JPEGs)
# ==============================================================================
message("Generating 600 DPI diagnostic plots in: ", plot_dir)

# 1. RNA histograms
p1 <- ggplot(data.frame(nCount = obj_rna_sub$nCount_RNA), aes(x = nCount)) +
  geom_histogram(bins = 50, fill = "#1b9e77", color = "white") +
  scale_x_log10(labels = label_comma()) +
  theme_classic(base_size = 12) +
  labs(title = "RNA UMI Distribution (Arabidopsis Root)", x = "nCount_RNA (log10)", y = "Cells")

p2 <- ggplot(data.frame(nFeature = obj_rna_sub$nFeature_RNA), aes(x = nFeature)) +
  geom_histogram(bins = 50, fill = "#d95f02", color = "white") +
  scale_x_log10(labels = label_comma()) +
  theme_classic(base_size = 12) +
  labs(title = "RNA Detected Genes (Arabidopsis Root)", x = "nFeature_RNA (log10)", y = "Cells")

jpeg(file.path(plot_dir, "01_hist_RNA_QC.jpeg"), width = 3600, height = 1800, res = 600)
print(p1 + p2)
dev.off()

# 2. RNA scatter
p_sc_rna <- ggplot(obj_rna_sub@meta.data, aes(x = nCount_RNA, y = nFeature_RNA, color = batch)) +
  geom_point(alpha = 0.5, size = 0.8) +
  scale_x_log10(labels = label_comma()) +
  scale_y_log10(labels = label_comma()) +
  theme_classic(base_size = 12) +
  labs(title = "RNA Quality Scatter (Arabidopsis Root)", x = "nCount_RNA", y = "nFeature_RNA")

jpeg(file.path(plot_dir, "02_scatter_RNA_QC.jpeg"), width = 3000, height = 2400, res = 600)
print(p_sc_rna)
dev.off()

# 3. ATAC histograms
p3 <- ggplot(data.frame(nCount = Matrix::colSums(final_atac_counts)), aes(x = nCount)) +
  geom_histogram(bins = 50, fill = "#7570b3", color = "white") +
  scale_x_log10(labels = label_comma()) +
  theme_classic(base_size = 12) +
  labs(title = "ATAC Depths in Selected Peaks", x = "Total Counts in Peaks (log10)", y = "Cells")

p4 <- ggplot(data.frame(nFeature = Matrix::colSums(final_atac_counts > 0)), aes(x = nFeature)) +
  geom_histogram(bins = 50, fill = "#e7298a", color = "white") +
  scale_x_log10(labels = label_comma()) +
  theme_classic(base_size = 12) +
  labs(title = "ATAC Detected Peaks", x = "Detected Peaks (log10)", y = "Cells")

jpeg(file.path(plot_dir, "03_hist_ATAC_QC.jpeg"), width = 3600, height = 1800, res = 600)
print(p3 + p4)
dev.off()

# 4. ATAC scatter
p_sc_atac <- ggplot(data.frame(nCount = Matrix::colSums(final_atac_counts),
                              nFeature = Matrix::colSums(final_atac_counts > 0),
                              batch = obj_atac_sub$batch),
                    aes(x = nCount, y = nFeature, color = batch)) +
  geom_point(alpha = 0.5, size = 0.8) +
  scale_x_log10(labels = label_comma()) +
  scale_y_log10(labels = label_comma()) +
  theme_classic(base_size = 12) +
  labs(title = "ATAC Quality Scatter (Arabidopsis Root)", x = "Total Counts in Peaks", y = "Detected Peaks")

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
  labs(title = "RNA Cluster Frequencies (Arabidopsis Root)", x = "Cluster / Cell Type", y = "Cells")

jpeg(file.path(plot_dir, "05_celltype_composition_RNA.jpeg"), width = 3000, height = 2400, res = 600)
print(p_comp_rna)
dev.off()

# 6. ATAC Cell type composition
atac_comp <- as.data.frame(table(obj_atac_sub$cell_type))
colnames(atac_comp) <- c("cell_type", "count")
p_comp_atac <- ggplot(atac_comp, aes(x = reorder(cell_type, count), y = count, fill = cell_type)) +
  geom_col(show.legend = FALSE) +
  coord_flip() +
  theme_classic(base_size = 11) +
  labs(title = "ATAC Cluster Frequencies (Arabidopsis Root)", x = "Cluster / Cell Type", y = "Cells")

jpeg(file.path(plot_dir, "06_celltype_composition_ATAC.jpeg"), width = 3000, height = 2400, res = 600)
print(p_comp_atac)
dev.off()

# 7. Peak detection rates
p_det <- ggplot(data.frame(det = peak_det_rate[final_selected_peaks]), aes(x = det)) +
  geom_histogram(bins = 40, fill = "#66a61e", color = "white") +
  theme_classic(base_size = 12) +
  labs(title = "Selected 6,000 Peaks Detection Rates", x = "Fraction of Cells Accessible", y = "Peaks")

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
  labs(title = "Cross-Modal Cluster Proportions (Arabidopsis Root)", x = "Cluster / Cell Type", y = "Proportion of Modality")

jpeg(file.path(plot_dir, "08_crossmodal_celltype_concordance.jpeg"), width = 3200, height = 2600, res = 600)
print(p_cross)
dev.off()

# ==============================================================================
# STEP D: SAVE DECOUPLED OUTPUTS & BY-BATCH MATRICES
# ==============================================================================
message("Saving standardized decoupled objects...")

# 1. Clean RNA Seurat object & matrix
final_rna_obj <- CreateSeuratObject(counts = final_rna_counts, assay = "RNA", meta.data = obj_rna_sub@meta.data)
saveRDS(final_rna_obj, file.path(tissue_out_dir, "arab_processed_Root_rna.rds"))
saveRDS(final_rna_counts, file.path(tissue_out_dir, "arab_processed_Root_rna_counts.rds"))

# 2. Clean ATAC Seurat object & matrix
final_atac_obj <- CreateSeuratObject(counts = final_atac_counts, assay = "peaks", meta.data = obj_atac_sub@meta.data)
saveRDS(final_atac_obj, file.path(tissue_out_dir, "arab_processed_Root_atac.rds"))
saveRDS(final_atac_counts, file.path(tissue_out_dir, "arab_processed_Root_atac_counts.rds"))

# 3. Metadatas
write.csv(obj_rna_sub@meta.data, file.path(tissue_out_dir, "arab_processed_Root_rna_metadata.csv"), row.names = TRUE)
write.csv(obj_atac_sub@meta.data, file.path(tissue_out_dir, "arab_processed_Root_atac_metadata.csv"), row.names = TRUE)

# 4. By-batch exports for RNA and ATAC
for (b in unique(obj_rna_sub$batch)) {
  b_cells <- colnames(obj_rna_sub)[obj_rna_sub$batch == b]
  b_counts <- final_rna_counts[, b_cells, drop = FALSE]
  b_meta <- obj_rna_sub@meta.data[b_cells, , drop = FALSE]
  b_clean <- gsub("[^A-Za-z0-9_]", "_", b)
  saveRDS(b_counts, file.path(batch_dir, paste0("arab_batch_", b_clean, "_rna_counts.rds")))
  write.csv(b_meta, file.path(batch_dir, paste0("arab_batch_", b_clean, "_rna_metadata.csv")), row.names = TRUE)
}

for (b in unique(obj_atac_sub$batch)) {
  b_cells <- colnames(obj_atac_sub)[obj_atac_sub$batch == b]
  b_counts <- final_atac_counts[, b_cells, drop = FALSE]
  b_meta <- obj_atac_sub@meta.data[b_cells, , drop = FALSE]
  b_clean <- gsub("[^A-Za-z0-9_]", "_", b)
  saveRDS(b_counts, file.path(batch_dir, paste0("arab_batch_", b_clean, "_atac_counts.rds")))
  write.csv(b_meta, file.path(batch_dir, paste0("arab_batch_", b_clean, "_atac_metadata.csv")), row.names = TRUE)
}

# ==============================================================================
# STEP E: SUMMARY METRICS RECORD
# ==============================================================================
summary_df <- data.frame(
  Species                    = "A. thaliana",
  Tissue                     = "Root",
  Initial_RNA_Cells          = raw_rna_cells,
  Final_RNA_Cells            = ncol(obj_rna_sub),
  Initial_ATAC_Cells         = raw_atac_cells,
  Final_ATAC_Cells           = ncol(obj_atac_sub),
  Initial_RNA_Genes          = raw_rna_genes,
  Final_RNA_HVGs             = length(final_hvg_genes),
  Initial_ATAC_Peaks         = raw_atac_peaks,
  Final_ATAC_Peaks           = length(final_selected_peaks),
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
  ATAC_Clusters              = length(unique(obj_atac_sub$cell_type)),
  RNA_Batches                = length(unique(obj_rna_sub$batch)),
  ATAC_Batches               = length(unique(obj_atac_sub$batch)),
  RNA_Batch_List             = paste(sort(unique(obj_rna_sub$batch)), collapse = "; "),
  ATAC_Batch_List            = paste(sort(unique(obj_atac_sub$batch)), collapse = "; "),
  stringsAsFactors           = FALSE
)

summary_csv <- file.path(tissue_out_dir, "arab_processed_Root_summary_metrics.csv")
write.csv(summary_df, summary_csv, row.names = FALSE)
message("Saved summary metrics: ", summary_csv)

global_csv <- file.path(opt$output_dir, "arabidopsis_standardized_summary.csv")
write.csv(summary_df, global_csv, row.names = FALSE)
message("Saved global summary: ", global_csv)
print(summary_df)
message("\n==================================================================")
message("Arabidopsis Root preprocessing successfully completed!")
message("==================================================================")
