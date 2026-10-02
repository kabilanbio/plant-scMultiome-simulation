# plant-scMultiome-simulation

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23068394.svg)](https://doi.org/10.5281/zenodo.23068394)
[![Python](https://img.shields.io/badge/Python-3.10%20%7C%203.11%20%7C%203.12-blue.svg)](https://www.python.org/)
[![R](https://img.shields.io/badge/R-%E2%89%A5%204.2.0-blue.svg)](https://www.r-project.org/)

Standardized single-cell multi-omics preprocessing and empirical PAGA lineage tree reconstruction pipelines across 16 plant tissues (*Oryza sativa*, *Glycine max*, and *Arabidopsis thaliana*).

---

## Repository Structure

```
plant-scMultiome-simulation/
├── 01_preprocessing/                          # Standardized 3:1 multiomic preprocessing (Seurat/Signac)
│   ├── preprocess_arabidopsis_root.R          # Arabidopsis Root atlas (21 clusters, GSE155304)
│   ├── preprocess_rice_all_tissues.R          # Rice 8 organs (Bud, Flag, Leaf, Root, SAM, Seed, SP, ST)
│   └── preprocess_soybean_all_tissues.R       # Soybean 7 tissues & developmental stages
└── 02_lineage_trees_paga/                     # Empirical PAGA lineage trajectory extraction (Scanpy)
    ├── arabidopsis/
    │   ├── export_arab_standardized_to_anndata.R
    │   ├── arab_paga_lineage_tree.py
    │   └── run_all_arab_paga.bat
    ├── rice/
    │   ├── export_rice_standardized_to_anndata.R
    │   ├── rice_paga_lineage_tree.py
    │   ├── run_all_rice_paga.bat
    │   └── run_rice_paga_cluster.pbs
    └── soybean/
        ├── export_soybean_standardized_to_anndata.R
        ├── soybean_paga_lineage_tree.py
        └── run_all_soybean_paga.bat
```

---

## Installation

Ensure R ($\ge$ 4.2.0) and Python ($\ge$ 3.10) are installed:

```bash
# Python dependencies
pip install -r requirements.txt

# Or via conda
conda env create -f environment.yml
```

---

## Usage

### 1. Preprocessing (Standardized 3:1 Feature Ratio)
Extracts 2,000 RNA HVGs $\times$ 6,000 accessible chromatin peaks with dual-stratified cell subsampling (max 200 cells per cell type $\times$ batch):

```bash
# Rice (all 8 tissues)
Rscript 01_preprocessing/preprocess_rice_all_tissues.R --tissue all

# Soybean (all 7 tissues)
Rscript 01_preprocessing/preprocess_soybean_all_tissues.R --tissue all

# Arabidopsis Root (21 clusters)
Rscript 01_preprocessing/preprocess_arabidopsis_root.R
```

### 2. Empirical Lineage Tree Reconstruction
Infers data-driven PAGA graphs, Diffusion Pseudotime (DPT), and Newick lineage trees (`all_tips.nwk` for `scMultiSim` and `classic.nwk`):

```bash
# Arabidopsis Root
cd 02_lineage_trees_paga/arabidopsis && python arab_paga_lineage_tree.py

# Rice
cd 02_lineage_trees_paga/rice && python rice_paga_lineage_tree.py --tissue all

# Soybean
cd 02_lineage_trees_paga/soybean && python soybean_paga_lineage_tree.py --tissue all
```

---

## Data Availability

Standardized count matrices, processed AnnData (`.h5ad`) objects, Newick lineage trees (`.nwk`), and 600 DPI publication figures are archived on Zenodo:
* **DOI**: [10.5281/zenodo.23068394](https://doi.org/10.5281/zenodo.23068394)

---

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
