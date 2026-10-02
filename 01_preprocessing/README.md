# Preprocessing Pipelines

Standardized preprocessing scripts for single-cell multi-omics (paired and unpaired snRNA-seq + snATAC-seq) across three plant species:
* **Rice (*Oryza sativa*)**: 8 distinct tissues (`preprocess_rice_all_tissues.R`)
* **Soybean (*Glycine max*)**: 7 tissues and developmental stages (`preprocess_soybean_all_tissues.R`)
* **Arabidopsis (*Arabidopsis thaliana*)**: Root single-cell atlas (`preprocess_arabidopsis_root.R`)

### Key Parameters:
- **3:1 Multiomic Ratio**: 2,000 RNA HVGs and 6,000 accessible chromatin peaks.
- **Dual-Stratified Sampling**: Max 200 cells per `cell_type` $\times$ `batch`.
- **Outputs**: Decoupled count matrices (`counts_rna.rds`, `counts_atac.rds`), `cell_types.rds`, and `batches.rds`.
