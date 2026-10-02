# Master Zenodo Upload Repository & Artifact Archive

**Local Path**: `G:\PhD\zenodo_upload_files\`  
**Target Zenodo Record**: [https://zenodo.org/records/23068394](https://zenodo.org/records/23068394)  
**Persistent DOI**: [10.5281/zenodo.23068394](https://doi.org/10.5281/zenodo.23068394)  
**Associated GitHub**: [https://github.com/kabilanbio/plant-scMultiome-simulation](https://github.com/kabilanbio/plant-scMultiome-simulation)

---

## 1. Directory Structure

```
G:\PhD\zenodo_upload_files\
├── README.md                              # This master archive guide
├── ZENODO_DATASET_MANIFEST.md             # Scientific manifest with SHA-256 hashes & specs
├── upload_to_zenodo.py                    # 1-Command REST API uploader script
│
├── zenodo_ready_zip_packages/             # READY TO UPLOAD DIRECTLY TO ZENODO
│   ├── paga_lineage_trees_and_metadata.zip   # 33.4 KB: All 16 Newick trees & metadata
│   ├── paga_processed_anndata_h5ad.zip       # 53.39 MB: 16 Scanpy AnnData with PAGA graphs
│   └── paga_publication_figures_600dpi.zip   # 114.36 MB: Publication 600 DPI figures
│
├── 01_preprocessing_standardized_datasets/# Uncompressed preprocessing outputs
│   ├── arabidopsis_root/                     # Standardized count matrices (RNA & ATAC)
│   ├── rice/                                 # Rice 8-tissue standardized summary metrics
│   └── soybean/                              # Soybean 7-tissue standardized summary metrics
│
└── 02_paga_lineage_trees_and_metadata/    # Uncompressed PAGA trajectory trees & graphs
    ├── arabidopsis/                          # 21-cluster Newick trees, metadata, connectivity
    ├── rice/                                 # 8 tissue subfolders (.nwk, .csv) + master summary
    └── soybean/                              # 7 tissue subfolders (.nwk, .csv) + master summary
```

---

## 2. Summary of Scientific Deliverables

### A. Preprocessing Standardization (3:1 Feature Ratio)
- **RNA HVGs**: 2,000 highly variable genes
- **ATAC Peaks**: 6,000 top accessible chromatin peaks
- **Dual-Stratified Sampling**: Max 200 cells per `cell_type` $\times$ `batch`
- **Tissues**: 16 distinct plant systems across Rice (*Oryza sativa*), Soybean (*Glycine max*), and Arabidopsis (*Arabidopsis thaliana*).

### B. Empirical Lineage Inference (PAGA + Diffusion Pseudotime)
- Empirical differentiation backbones computed from statistical connectivity graphs.
- **Botanical Progenitor Rooting**: Pinned to stem cells, apical meristems (RAM/SAM), or embryonic axes.
- **Dual Newick Tree Architecture**:
  * `*_paga_tree_all_tips.nwk` (or `arab_root_paga_tree_21tips.nwk`): Guaranteed 100% terminal tip matching for `scMultiSim` population allocation vectors (`discrete.pop.size`).
  * `*_paga_tree_classic.nwk`: Traditional hierarchical branching topology.

---

## 3. How to Upload to Zenodo Record 23068394

### Method 1: Web Interface Drag & Drop (Fastest & Simplest)
1. Open [https://zenodo.org/records/23068394](https://zenodo.org/records/23068394).
2. Click **Edit** (or **New Version** if editing a published version).
3. Drag and drop the 3 zip files from `G:\PhD\zenodo_upload_files\zenodo_ready_zip_packages\` and `ZENODO_DATASET_MANIFEST.md`.
4. Click **Save** and **Publish**.

### Method 2: Automated via CLI
```powershell
python G:\PhD\zenodo_upload_files\upload_to_zenodo.py --token YOUR_ZENODO_PERSONAL_ACCESS_TOKEN
```
