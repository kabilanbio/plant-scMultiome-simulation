# Standardized Single-Cell Multi-Omics Preprocessing Pipeline

This directory contains reproducible, production-grade preprocessing scripts for single-cell multiomics (paired and unpaired snRNA-seq + snATAC-seq) across three plant species:
1. **Rice (*Oryza sativa*)**: 8 distinct tissues (Bud, Flag leaf, Leaf, Root, SAM, Seed, Spikelet/Panicle, Stem)
2. **Soybean (*Glycine max*)**: 7 tissues/developmental stages (Cotyledon seeds, Early maturation seeds, Early nodule, Globular seeds, Heart seeds, Hypocotyl, Root)
3. **Arabidopsis (*Arabidopsis thaliana*)**: Comprehensive Root cell atlas across 21 developmental clusters and 6 canonical cell types

---

## 1. Standardization Principles & Constraints

To ensure computational feasibility, uniform cross-species comparison, and exact compatibility with `scMultiSim` ([Zhang et al., *Nature Methods* 2023](https://doi.org/10.1038/s41592-023-02051-0)):

1. **3:1 Multiomic Feature Ratio**:
   - **RNA HVGs**: 2,000 highly variable genes selected via variance-stabilizing transformation (`vst`).
   - **ATAC Peaks**: 6,000 top accessible chromatin peaks selected by frequency and dispersion.
   - Strictly enforces the biological and mathematical 3:1 ratio required for coupled multimodal simulation.

2. **Dual-Stratified Subsampling**:
   - Caps cells at a maximum of 200 cells per `cell_type` $\times$ `batch` combination.
   - Prevents dominant cell types from causing out-of-memory errors on standard workstations (16 GB RAM).
   - Preserves rare transitional populations and statistical representation of every lineage.

3. **Decoupled Architecture**:
   - Each tissue yields decoupled count matrices (`counts_rna.rds`, `counts_atac.rds`), unified cell metadata (`cell_types.rds`, `batches.rds`), and subsetted Seurat objects.

---

## 2. Scripts Description

| Script | Species | Architecture | Features | Description |
| :--- | :--- | :--- | :--- | :--- |
| `preprocess_rice_all_tissues.R` | Rice (*O. sativa*) | Paired multiome | 2,000 HVGs, 6,000 peaks | Processes all 8 rice tissues with automated CLI flags |
| `preprocess_soybean_all_tissues.R` | Soybean (*G. max*) | Paired multiome | 2,000 HVGs, 6,000 peaks | Processes all 7 soybean tissues/developmental stages |
| `preprocess_arabidopsis_root.R` | Arabidopsis (*A. thaliana*) | Unpaired multiome (GSE155304) | 2,000 HVGs, 6,000 peaks | Harmonizes snRNA + snATAC into standardized 21-cluster root atlas |

---

## 3. Usage Instructions

Run locally in R (version $\ge$ 4.2):

```bash
# Process all rice tissues
Rscript preprocess_rice_all_tissues.R --tissue all

# Process all soybean tissues
Rscript preprocess_soybean_all_tissues.R --tissue all

# Process Arabidopsis root dataset
Rscript preprocess_arabidopsis_root.R
```
