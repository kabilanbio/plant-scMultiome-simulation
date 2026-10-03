# ==============================================================================
# rice_plantFigR_GRN.R
# Scientific Gene Regulatory Network (GRN) Inference for Rice Multi-Omics
# Using plantFigR (Adapted from FigR, Buenrostro Lab)
# Generates Trajectory-Conditioned GRNs Formatted Specifically for scMultiSim
#
# Processed Tissues:
#   Bud, Flag, Leaf, Root, SAM, Seed, SP, ST
# ==============================================================================

suppressPackageStartupMessages({
  library(Matrix)
  library(SummarizedExperiment)
  library(GenomicRanges)
  library(IRanges)
  library(GenomeInfoDb)
  library(ggplot2)
  library(dplyr)
  library(scales)
  library(FNN)
  library(doParallel)
  library(foreach)
})

# Ensure required Bioconductor and specialized packages are available
ensure_package <- function(pkg, is_bioc = FALSE) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(">>> Package '", pkg, "' is not installed. Installing now... <<<")
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager", repos = "https://cloud.r-project.org")
    }
    if (is_bioc) {
      BiocManager::install(pkg, update = FALSE, ask = FALSE)
    } else {
      install.packages(pkg, repos = "https://cloud.r-project.org")
    }
  }
}

ensure_package("chromVAR", is_bioc = TRUE)
ensure_package("BSgenome.Osativa.MSU.MSU7", is_bioc = TRUE)
ensure_package("plantFigR", is_bioc = FALSE)
ensure_package("FigR", is_bioc = FALSE)
ensure_package("inflection", is_bioc = FALSE)
ensure_package("ggrepel", is_bioc = FALSE)

suppressPackageStartupMessages({
  library(chromVAR)
  library(BSgenome.Osativa.MSU.MSU7)
  library(plantFigR)
  library(FigR)
  library(inflection)
  library(ggrepel)
})

# ------------------------------------------------------------------------------
# 1. Native CLI Parameter Parser
# ------------------------------------------------------------------------------
parse_cli_arguments <- function(args = commandArgs(trailingOnly = TRUE)) {
  script_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", script_args, value = TRUE)
  this_dir <- if (length(file_arg) > 0) normalizePath(dirname(sub("^--file=", "", file_arg[1]))) else getwd()
  drive_prefix <- substr(this_dir, 1, 2)
  
  params <- list(
    tissue       = "all",
    data_dir     = file.path(drive_prefix, "PhD", "sc_datasets", "rice", "featured_datasets_standardized"),
    ref_dir      = file.path(drive_prefix, "PhD", "sc_datasets", "rice"),
    full_rna     = file.path(drive_prefix, "PhD", "sc_datasets", "rice", "original_dataset", "rice_original_dataset", "GSE232863_scRNA_omics.Rds"),
    output_dir   = file.path(drive_prefix, "PhD", "data_simulation_comparison", "scMultiSim", "rice_datasets", "03_grn_plantFigR"),
    window       = 50000,
    n_bg         = 50,
    ncores       = 10,
    score_cutoff = 1.0,
    dpi          = 600
  )
  
  i <- 1
  while (i <= length(args)) {
    arg <- args[i]
    if (arg == "--tissue" && i < length(args)) {
      params$tissue <- args[i + 1]; i <- i + 2
    } else if (arg == "--data_dir" && i < length(args)) {
      params$data_dir <- args[i + 1]; i <- i + 2
    } else if (arg == "--ref_dir" && i < length(args)) {
      params$ref_dir <- args[i + 1]; i <- i + 2
    } else if (arg == "--full_rna" && i < length(args)) {
      params$full_rna <- args[i + 1]; i <- i + 2
    } else if (arg == "--output_dir" && i < length(args)) {
      params$output_dir <- args[i + 1]; i <- i + 2
    } else if (arg == "--window" && i < length(args)) {
      params$window <- as.numeric(args[i + 1]); i <- i + 2
    } else if (arg == "--n_bg" && i < length(args)) {
      params$n_bg <- as.numeric(args[i + 1]); i <- i + 2
    } else if (arg == "--ncores" && i < length(args)) {
      params$ncores <- as.integer(args[i + 1]); i <- i + 2
    } else if (arg == "--score_cutoff" && i < length(args)) {
      params$score_cutoff <- as.numeric(args[i + 1]); i <- i + 2
    } else if (arg == "--dpi" && i < length(args)) {
      params$dpi <- as.numeric(args[i + 1]); i <- i + 2
    } else {
      i <- i + 1
    }
  }
  
  return(params)
}

# ------------------------------------------------------------------------------
# 2. Peak Parser & SummarizedExperiment Builder
# ------------------------------------------------------------------------------
build_rice_atac_se <- function(atac_counts, genome = BSgenome.Osativa.MSU.MSU7) {
  peak_names <- rownames(atac_counts)
  
  # Support both 'chr1-start-end' and 'Chr1:start-end' formats
  if (all(grepl(":", peak_names) & grepl("-", peak_names))) {
    peak_chr_raw <- sub(":.*", "", peak_names)
    coords <- sub(".*:", "", peak_names)
    peak_start <- as.integer(sub("-.*", "", coords))
    peak_end   <- as.integer(sub(".*-", "", coords))
  } else {
    parts <- strsplit(peak_names, "-")
    peak_chr_raw <- sapply(parts, `[`, 1)
    peak_start   <- as.integer(sapply(parts, `[`, 2))
    peak_end     <- as.integer(sapply(parts, `[`, 3))
  }
  
  # Normalize to 'Chr1' to 'Chr12' format matching BSgenome
  peak_chr <- gsub("^chr", "Chr", peak_chr_raw)
  
  peak_gr <- GRanges(
    seqnames = peak_chr,
    ranges   = IRanges(start = peak_start, end = peak_end)
  )
  names(peak_gr) <- peak_names
  
  # Filter to canonical chromosomes Chr1 - Chr12
  main_chrs <- paste0("Chr", 1:12)
  keep_peaks <- as.logical(seqnames(peak_gr) %in% main_chrs)
  
  atac_counts_filt <- atac_counts[keep_peaks, ]
  peak_gr_filt     <- peak_gr[keep_peaks]
  peak_gr_filt     <- keepSeqlevels(peak_gr_filt, main_chrs, pruning.mode = "coarse")
  
  # Assign sequence lengths and genome
  sl_names <- intersect(seqlevels(peak_gr_filt), seqlevels(genome))
  seqlengths(peak_gr_filt)[sl_names] <- seqlengths(genome)[sl_names]
  genome(peak_gr_filt) <- unique(genome(genome))
  
  # Trim out-of-bound peaks
  oob <- end(peak_gr_filt) > seqlengths(peak_gr_filt)[as.character(seqnames(peak_gr_filt))]
  if (any(oob)) {
    atac_counts_filt <- atac_counts_filt[!oob, ]
    peak_gr_filt     <- peak_gr_filt[!oob]
  }
  
  # Standardize peak row names
  rownames(atac_counts_filt) <- paste0(seqnames(peak_gr_filt), ":", start(peak_gr_filt), "-", end(peak_gr_filt))
  names(peak_gr_filt)        <- rownames(atac_counts_filt)
  
  # Construct SummarizedExperiment
  ATAC.se <- SummarizedExperiment(
    assays    = list(counts = atac_counts_filt),
    rowRanges = peak_gr_filt
  )
  
  # Add GC bias for background matching
  message("  [ATAC] Computing peak GC bias with chromVAR...")
  ATAC.se <- chromVAR::addGCBias(ATAC.se, genome = genome)
  
  bias_vals <- rowData(ATAC.se)$bias
  if (any(is.na(bias_vals))) {
    rowData(ATAC.se)$bias[is.na(bias_vals)] <- median(bias_vals, na.rm = TRUE)
  }
  
  return(ATAC.se)
}

# ------------------------------------------------------------------------------
# 3. Core plantFigR Pipeline per Rice Tissue
# ------------------------------------------------------------------------------
process_rice_tissue_grn <- function(tissue, params, rice_TSSg, rice_motifs, full_rna_obj = NULL) {
  cat("\n==============================================================================\n")
  cat(sprintf("[plantFigR GRN Inference] Processing Tissue: %s\n", tissue))
  cat("==============================================================================\n")
  
  t0 <- Sys.time()
  tissue_in_dir  <- file.path(params$data_dir, tissue)
  tissue_out_dir <- file.path(params$output_dir, tissue)
  fig_dir        <- file.path(tissue_out_dir, "figures")
  dir.create(tissue_out_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(fig_dir,        recursive = TRUE, showWarnings = FALSE)
  
  # 1. Load Standardized Input Multi-Omics Datasets
  atac_rds_file <- file.path(tissue_in_dir, sprintf("rice_processed_%s_atac.rds", tissue))
  rna_rds_file  <- file.path(tissue_in_dir, sprintf("rice_processed_%s_rna.rds", tissue))
  
  if (!file.exists(atac_rds_file) || !file.exists(rna_rds_file)) {
    warning(sprintf("Standardized RDS files missing for tissue '%s' in %s", tissue, tissue_in_dir))
    return(NULL)
  }
  
  message("  Loading standardized ATAC and RNA objects...")
  atac_obj <- readRDS(atac_rds_file)
  rna_obj  <- readRDS(rna_rds_file)
  
  # Extract raw counts and log-normalized data
  if (inherits(atac_obj, "Seurat")) {
    atac_counts <- Seurat::GetAssayData(atac_obj, assay = "ATAC", layer = "counts")
  } else {
    atac_counts <- atac_obj
  }
  
  if (inherits(rna_obj, "Seurat")) {
    rna_data <- Seurat::GetAssayData(rna_obj, assay = "RNA", layer = "data")
  } else {
    rna_data <- rna_obj
  }
  
  # Ensure cells match 1:1 between modalities
  common_cells <- intersect(colnames(atac_counts), colnames(rna_data))
  cat(sprintf("  Matched multiome cells: %d\n", length(common_cells)))
  atac_counts <- atac_counts[, common_cells]
  rna_data    <- rna_data[, common_cells]
  
  # 2. Build ATAC SummarizedExperiment & GC Bias
  ATAC.se <- build_rice_atac_se(atac_counts, genome = BSgenome.Osativa.MSU.MSU7)
  cat(sprintf("  ATAC.se constructed: %d peaks x %d cells\n", nrow(ATAC.se), ncol(ATAC.se)))
  
  # 3. Intersect RNA with TSSg
  rna_genes    <- intersect(rownames(rna_data), names(rice_TSSg))
  rna_mat_filt <- rna_data[rna_genes, ]
  cat(sprintf("  RNA genes overlapping with TSSg: %d / %d\n", length(rna_genes), nrow(rna_data)))
  
  # 4. Step A: Gene-Peak Correlation (DORC Calling)
  message("  Running gene-peak correlation testing (runGenePeakcorr_plant)...")
  dorcTab <- plantFigR::runGenePeakcorr_plant(
    ATAC.se               = ATAC.se,
    RNAmat                = rna_mat_filt,
    TSSg                  = rice_TSSg,
    bsgenome              = BSgenome.Osativa.MSU.MSU7,
    windowPadSize         = params$window,
    normalizeATACmat      = TRUE,
    nCores                = params$ncores,
    keepPosCorOnly        = TRUE,
    keepMultiMappingPeaks = FALSE,
    n_bg                  = params$n_bg,
    p.cut                 = NULL
  )
  
  cat(sprintf("  Total gene-peak pairs tested: %d\n", nrow(dorcTab)))
  saveRDS(dorcTab, file.path(tissue_out_dir, sprintf("%s_dorcTab_all.rds", tissue)))
  
  # 5. Step B: Filter Significant Pairs & Knee Point Cutoff
  dorcTab_sig <- dorcTab[dorcTab$pvalZ < 0.05, ]
  dorcTab_sig <- dorcTab_sig[!is.na(dorcTab_sig$Gene), ]
  cat(sprintf("  Significant gene-peak pairs (pvalZ < 0.05): %d\n", nrow(dorcTab_sig)))
  saveRDS(dorcTab_sig, file.path(tissue_out_dir, sprintf("%s_dorcTab_sig.rds", tissue)))
  
  numDorcs <- dorcTab_sig %>%
    dplyr::group_by(Gene) %>%
    dplyr::tally() %>%
    dplyr::arrange(desc(n))
  
  knee <- tryCatch({
    inflection::uik(x = 1:nrow(numDorcs), y = numDorcs$n)
  }, error = function(e) {
    min(5, max(1, round(median(numDorcs$n))))
  })
  
  auto_cutoff <- max(1, numDorcs$n[knee])
  cat(sprintf("  DORC Knee-point cutoff: %d peaks per gene\n", auto_cutoff))
  
  # DORC J-Plot Figure
  p_jplot <- ggplot(numDorcs %>% dplyr::mutate(rank = 1:n()), aes(x = rank, y = n)) +
    geom_point(size = 0.8, color = "gray40") +
    geom_hline(yintercept = auto_cutoff, linetype = "dashed", color = "firebrick", linewidth = 0.8) +
    annotate("text", x = nrow(numDorcs)*0.6, y = auto_cutoff + 1, 
             label = sprintf("Cutoff = %d peaks (knee)", auto_cutoff), color = "firebrick", fontface = "bold") +
    theme_classic(base_size = 12) +
    labs(x = "Gene Rank", y = "Significant Peaks per Gene", 
         title = sprintf("Rice %s: DORC Identification (J-Plot)", tissue))
  
  ggsave(file.path(fig_dir, sprintf("Fig1_%s_DORC_Jplot.png", tissue)),
         plot = p_jplot, width = 7, height = 5, dpi = params$dpi, bg = "white")
  
  dorcGenes <- numDorcs$Gene[numDorcs$n >= auto_cutoff]
  cat(sprintf("  Identified DORC Genes: %d\n", length(dorcGenes)))
  saveRDS(dorcGenes, file.path(tissue_out_dir, sprintf("%s_dorcGenes.rds", tissue)))
  
  # Volcano Figure (Correlation vs Significance)
  p_volcano <- ggplot(dorcTab, aes(x = rObs, y = -log10(pvalZ + 1e-300))) +
    geom_point(size = 0.4, alpha = 0.3, color = "gray50") +
    geom_point(data = dorcTab_sig, size = 0.6, alpha = 0.5, color = "dodgerblue4") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "firebrick") +
    theme_classic(base_size = 12) +
    labs(x = "Observed Correlation (rObs)", y = "-log10(pvalZ)",
         title = sprintf("Rice %s: Gene-Peak Association Significance", tissue))
  
  ggsave(file.path(fig_dir, sprintf("Fig2_%s_gene_peak_volcano.png", tissue)),
         plot = p_volcano, width = 7, height = 5, dpi = params$dpi, bg = "white")
  
  # 6. Step C: Per-Cell DORC Accessibility Scores
  message("  Computing per-cell DORC scores (FigR::getDORCScores)...")
  dorcMat_raw <- FigR::getDORCScores(
    ATAC.se          = ATAC.se,
    dorcTab          = dorcTab_sig,
    normalizeATACmat = TRUE,
    nCores           = params$ncores
  )
  
  # 7. Step D: Cell KNN Graph & Smoothing
  message("  Constructing KNN graph & smoothing accessibility/expression scores...")
  # Compute PCA / SVD embedding for KNN
  if (inherits(atac_obj, "Seurat") && "lsi" %in% names(atac_obj@reductions)) {
    embed_mat <- Seurat::Embeddings(atac_obj, "lsi")[common_cells, 2:min(30, ncol(Seurat::Embeddings(atac_obj, "lsi")))]
  } else {
    svd_res   <- irlba::irlba(Matrix::t(atac_counts), nv = min(30, ncol(atac_counts) - 1))
    embed_mat <- svd_res$u %*% diag(svd_res$d)
    rownames(embed_mat) <- common_cells
  }
  
  knn_res <- FNN::get.knn(data = embed_mat, k = 30)
  NNmat   <- knn_res$nn.index
  rownames(NNmat) <- common_cells
  
  dorcMat_smooth <- FigR::smoothScoresNN(NNmat = NNmat, mat = dorcMat_raw, nCores = params$ncores)
  
  # 8. Step E: Transcription Factor Expression Extraction
  tf_names <- names(rice_motifs)
  
  # Check if full RNA matrix is available to retrieve lowly expressed TFs
  if (!is.null(full_rna_obj) && inherits(full_rna_obj, "Seurat")) {
    full_rna_data <- Seurat::GetAssayData(full_rna_obj, assay = "RNA", layer = "data")
    avail_cells   <- intersect(common_cells, colnames(full_rna_data))
    matched_tfs   <- intersect(rownames(full_rna_data), tf_names)
    rna_tf_raw    <- full_rna_data[matched_tfs, avail_cells]
  } else {
    matched_tfs   <- intersect(rownames(rna_data), tf_names)
    rna_tf_raw    <- rna_data[matched_tfs, common_cells]
  }
  
  cat(sprintf("  Transcription Factors with motifs available in expression matrix: %d\n", length(matched_tfs)))
  
  if (length(matched_tfs) < 5) {
    warning("  Very few TFs found in expression matrix! Expanding with all available regulators.")
  }
  
  rna_tf_smooth <- FigR::smoothScoresNN(NNmat = NNmat, mat = rna_tf_raw, nCores = params$ncores)
  rna_tf_dense  <- as.matrix(rna_tf_smooth)
  
  # 9. Step F: Run plantFigR GRN Inference
  message("  Inferring Gene Regulatory Network (runFigRGRN_plant)...")
  figR_res <- plantFigR::runFigRGRN_plant(
    ATAC.se   = ATAC.se,
    dorcK     = 30,
    dorcTab   = dorcTab_sig,
    n_bg      = params$n_bg,
    bsgenome  = BSgenome.Osativa.MSU.MSU7,
    pwm       = rice_motifs,
    dorcMat   = dorcMat_smooth,
    rnaMat    = rna_tf_dense,
    dorcGenes = dorcGenes,
    nCores    = params$ncores
  )
  
  saveRDS(figR_res, file.path(tissue_out_dir, sprintf("%s_figR_GRN.rds", tissue)))
  cat(sprintf("  Raw FigR GRN interactions: %d\n", nrow(figR_res)))
  
  # 10. Step G: Format specifically for scMultiSim
  # scMultiSim requires a data.frame with: TF, Target, Effect
  figR_clean <- figR_res %>%
    dplyr::filter(!is.na(Score) & Score != 0) %>%
    dplyr::filter(abs(Score) >= params$score_cutoff) %>%
    dplyr::select(
      TF     = Motif,
      Target = DORC,
      Effect = Score,
      Corr,
      Corr.P,
      Enrichment.P
    ) %>%
    dplyr::arrange(desc(abs(Effect)))
  
  cat(sprintf("  scMultiSim filtered GRN edges (|Score| >= %.1f): %d (TFs: %d, Targets: %d)\n",
              params$score_cutoff, nrow(figR_clean), length(unique(figR_clean$TF)), length(unique(figR_clean$Target))))
  
  # Save scMultiSim-ready tables
  scMultiSim_df <- figR_clean[, c("TF", "Target", "Effect")]
  write.csv(scMultiSim_df, file.path(tissue_out_dir, sprintf("%s_scMultiSim_GRN.csv", tissue)), row.names = FALSE)
  saveRDS(scMultiSim_df,   file.path(tissue_out_dir, sprintf("%s_scMultiSim_GRN.rds", tissue)))
  
  # 11. Summary Metrics CSV
  summary_row <- data.frame(
    Tissue            = tissue,
    Cells             = length(common_cells),
    ATAC_Peaks        = nrow(ATAC.se),
    RNA_Genes         = nrow(rna_mat_filt),
    Sig_Gene_Peak     = nrow(dorcTab_sig),
    DORC_Genes        = length(dorcGenes),
    Active_TFs        = length(unique(figR_clean$TF)),
    scMultiSim_Edges  = nrow(scMultiSim_df),
    Activation_Edges  = sum(scMultiSim_df$Effect > 0),
    Repression_Edges  = sum(scMultiSim_df$Effect < 0),
    Runtime_Mins      = round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2)
  )
  write.csv(summary_row, file.path(tissue_out_dir, sprintf("%s_grn_summary.csv", tissue)), row.names = FALSE)
  
  cat(sprintf(">>> Tissue '%s' complete in %.2f minutes! Outputs saved in: %s <<<\n",
              tissue, summary_row$Runtime_Mins, tissue_out_dir))
  
  return(summary_row)
}

# ------------------------------------------------------------------------------
# 4. Master Execution Loop
# ------------------------------------------------------------------------------
main <- function() {
  params <- parse_cli_arguments()
  
  cat("==============================================================================\n")
  cat("[Rice Multi-Omics plantFigR Gene Regulatory Network (GRN) Pipeline]\n")
  cat(sprintf("  Standardized Datasets Directory: %s\n", params$data_dir))
  cat(sprintf("  Reference Annotations Directory: %s\n", params$ref_dir))
  cat(sprintf("  Output Directory:               %s\n", params$output_dir))
  cat(sprintf("  CPU Cores:                      %d\n", params$ncores))
  cat(sprintf("  Target Tissue(s):               %s\n", params$tissue))
  cat("==============================================================================\n\n")
  
  # Load Master Reference Files
  tss_file   <- file.path(params$ref_dir, "rice_MSU7_TSSg.rds")
  motif_file <- file.path(params$ref_dir, "rice_motifs_CISBP_annotated_RAP.rds")
  
  if (!file.exists(tss_file)) stop("TSS reference file missing: ", tss_file)
  if (!file.exists(motif_file)) stop("Motif reference file missing: ", motif_file)
  
  message("Loading rice TSS genomic ranges (MSU7)...")
  rice_TSSg <- readRDS(tss_file)
  if (!"gene_name" %in% names(GenomicRanges::mcols(rice_TSSg))) {
    rice_TSSg$gene_name <- as.character(rice_TSSg$gene_id)
  }
  
  message("Loading rice CIS-BP motif library...")
  rice_motifs <- readRDS(motif_file)
  
  # Optional: Full RNA matrix for comprehensive TF expression
  full_rna_obj <- NULL
  if (file.exists(params$full_rna)) {
    message("Loading full RNA reference object for maximum TF coverage...")
    full_rna_obj <- tryCatch(readRDS(params$full_rna), error = function(e) NULL)
  }
  
  # Determine target tissues
  all_tissues <- c("Bud", "Flag", "Leaf", "Root", "SAM", "Seed", "SP", "ST")
  if (tolower(params$tissue) == "all") {
    target_tissues <- all_tissues
  } else {
    target_tissues <- intersect(params$tissue, all_tissues)
    if (length(target_tissues) == 0) {
      stop("Invalid tissue name: '", params$tissue, "'. Available: ", paste(all_tissues, collapse = ", "))
    }
  }
  
  all_summaries <- list()
  for (t in target_tissues) {
    res <- tryCatch({
      process_rice_tissue_grn(t, params, rice_TSSg, rice_motifs, full_rna_obj)
    }, error = function(e) {
      message(sprintf("\n[ERROR] Processing tissue '%s' failed: %s\n", t, e$message))
      NULL
    })
    if (!is.null(res)) all_summaries[[t]] <- res
  }
  
  if (length(all_summaries) > 0) {
    master_summary <- do.call(rbind, all_summaries)
    write.csv(master_summary, file.path(params$output_dir, "rice_all_tissues_grn_summary.csv"), row.names = FALSE)
    cat("\n==============================================================================\n")
    cat("Rice plantFigR GRN Pipeline Complete for all tissues!\n")
    print(master_summary)
    cat("==============================================================================\n")
  }
}

main()
