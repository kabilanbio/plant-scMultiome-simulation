# plant-scMultiome-simulation

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23068394.svg)](https://doi.org/10.5281/zenodo.23068394)
[![Language: R](https://img.shields.io/badge/Language-R%20%E2%89%A5%204.2-blue.svg)](https://www.r-project.org/)
[![Language: Python](https://img.shields.io/badge/Language-Python%20%E2%89%A5%203.10-blue.svg)](https://www.python.org/)

This repository contains custom analysis scripts supporting the manuscript:
> **Benchmarking Single-Cell Multi-Omics Simulation and Integration Methods in Plant Systems**  
> *Authors: Kabilan Sakthivel, Dwijesh Chandra Mishra, et al.*  
> ICAR - Indian Agricultural Research Institute (IARI) & ICAR - Indian Agricultural Statistics Research Institute (IASRI)

This repository serves as a transparent code archive detailing the computational workflows, parameterizations, and data processing routines utilized in this study.

---

## Repository Structure & Script Catalog

```
plant-scMultiome-simulation/
├── 01_preprocessing/                          # Multi-omics standardization workflows
│   ├── preprocess_arabidopsis_root.R          # Arabidopsis Root atlas (21 developmental clusters)
│   ├── preprocess_rice_all_tissues.R          # Rice (8 organs/tissues)
│   └── preprocess_soybean_all_tissues.R       # Soybean (7 tissues & developmental stages)
└── 02_lineage_trees_paga/                     # Empirical PAGA lineage trajectory extraction
    ├── arabidopsis/
    │   ├── export_arab_standardized_to_anndata.R
    │   └── arab_paga_lineage_tree.py
    ├── rice/
    │   ├── export_rice_standardized_to_anndata.R
    │   └── rice_paga_lineage_tree.py
    └── soybean/
        ├── export_soybean_standardized_to_anndata.R
        └── soybean_paga_lineage_tree.py
```

---

## Script Descriptions & Methodological Roles

### 1. Data Preprocessing & Standardization (`01_preprocessing/`)
Implements uniform feature selection and dual-stratified subsampling across 16 plant tissues to establish consistent 3:1 multiomic feature ratios (2,000 RNA HVGs to 6,000 accessible chromatin peaks) while preserving cell-type representations:

* **`preprocess_arabidopsis_root.R`**: Integrates and standardizes unpaired single-nucleus RNA-seq and ATAC-seq data from the *Arabidopsis thaliana* root atlas (GEO: GSE155304), resolving 21 developmental clusters across 6 canonical cell types.
* **`preprocess_rice_all_tissues.R`**: Standardizes paired single-cell multiomics data across 8 *Oryza sativa* tissues: Bud, Flag leaf, Young leaf, Root, Shoot Apical Meristem (SAM), Developing seed, Spikelet/Panicle (SP), and Stem (ST).
* **`preprocess_soybean_all_tissues.R`**: Processes paired multiomics datasets across 7 *Glycine max* tissues and seed developmental stages (Cotyledon, Early maturation, Early nodule, Globular, Heart, Hypocotyl, and Root).

### 2. Empirical Lineage Trajectory Inference (`02_lineage_trees_paga/`)
Derives data-driven cell lineage trees from empirical high-dimensional transcriptomic geometry using Partition-based Graph Abstraction (PAGA) and Diffusion Pseudotime (DPT) to condition multiomics simulations in `scMultiSim`:

* **`export_*_standardized_to_anndata.R`**: Translates standardized Seurat objects into Scanpy-compatible AnnData matrix formats.
* **`*_paga_lineage_tree.py`**: 
  - Computes statistical connectivity between cell clusters and identifies biological progenitor roots (e.g., Root Apical Meristem, SAM, Procambium).
  - Calculates Diffusion Pseudotime (DPT) and Maximum Spanning Tree (MST) trajectory backbones.
  - Exports dual Newick formats: `all_tips.nwk` (guaranteeing exact 1-to-1 tip matching for `scMultiSim` population vectors) and `classic.nwk` (macro-topology).

---

## Software Dependencies & Environment

The analysis was conducted using the following computational environments:

* **R (v4.6.1 / ≥ 4.2.0)**: `Seurat` (≥ 5.0.0), `Signac` (≥ 1.12.0), `Matrix`, `ggplot2`, `dplyr`, `ape`
* **Python (v3.12 / ≥ 3.10)**: `scanpy` (1.12.4), `anndata` (0.13.4), `igraph` (1.0.0), `leidenalg` (0.12.0), `networkx` (3.4), `biopython` (1.85)

Dependency specifications are provided in `environment.yml` and `requirements.txt`.

---

## Data Availability

All standardized matrices, processed AnnData (`.h5ad`) objects, empirical Newick lineage trees (`.nwk`), and high-resolution trajectory figures are deposited in Zenodo:
* **Zenodo Record**: [https://zenodo.org/records/23068394](https://zenodo.org/records/23068394)
* **Persistent DOI**: [10.5281/zenodo.23068394](https://doi.org/10.5281/zenodo.23068394)

---

## License

The code in this repository is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
