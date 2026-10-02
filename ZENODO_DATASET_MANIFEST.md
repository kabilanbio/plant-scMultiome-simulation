# Zenodo Data Manifest & Description

**Zenodo Record**: [https://zenodo.org/records/23068394](https://zenodo.org/records/23068394)  
**Persistent DOI**: [10.5281/zenodo.23068394](https://doi.org/10.5281/zenodo.23068394)  
**Title**: *Single cell multiomics data simulation analysis results*  
**Authors**: SAKTHIVEL, KABILAN & Mishra, Dwijesh Chandra  
**Affiliations**: ICAR - Indian Agricultural Research Institute (IARI) & ICAR - Indian Agricultural Statistics Research Institute (IASRI)

---

## 1. Overview of Deposited Packages

This deposition contains empirical cell lineage topologies, standardized single-cell multiomics graphs, and high-resolution trajectory visualizations generated for cross-species *in silico* benchmarking in plant biology using **scMultiSim** ([Zhang et al., *Nature Methods* 2023](https://doi.org/10.1038/s41592-023-02051-0)) and **PAGA** ([Wolf et al., *Genome Biology* 2019](https://doi.org/10.1186/s13059-019-1663-x)).

### Package Files & Integrity Hashes

| Filename | File Size | SHA-256 Checksum | Description |
| :--- | :--- | :--- | :--- |
| `paga_lineage_trees_and_metadata.zip` | 33.4 KB | `d85e73d181fdb6f8da26017c2722f34b60894dc789f355260080c6b67b57173d` | Empirical Newick lineage trees (`.nwk`), connectivity matrices (`.csv`), and cell-type metadata across 16 plant tissues. |
| `paga_processed_anndata_h5ad.zip` | 53.39 MB | `2f910497226893db7e3b2dc769800b4109b59acd2bb0de98ec9f4c588b426160` | 16 standardized Scanpy AnnData (`.h5ad`) objects with PAGA graph abstraction, Diffusion Pseudotime (DPT), and UMAP embeddings. |
| `paga_publication_figures_600dpi.zip` | 114.36 MB | `e6645dcd391cce3de2075abde38686769e9536231063e8a4fa6850104a3846bf` | Publication-ready 600 DPI figures (PNG/PDF) showing PAGA graphs, single-cell UMAPs, pseudotime gradients, and MST layouts. |

---

## 2. Detailed Dataset Structure

### A. Arabidopsis (*Arabidopsis thaliana* - Root Atlas, GSE155304)
- **Clusters**: 21 developmental clusters (`a14_1` to `u15_21`).
- **Cell Types (6 canonical classes across 4 tissue layers)**:
  * Epidermis: `Trichoblast` (c1–c3), `Atrichoblast` (c4–c7)
  * Meristem: `Meristematic cell` (c8–c10, **Root Apical Meristem Progenitor**, pseudotime = 0.1122)
  * Ground tissue: `Cortical cell` (c11–c12), `Endodermal cell` (c13–c16)
  * Vasculature: `Stele cell` (c17–c21)
- **Dual Newick Architecture**:
  * `arab_root_paga_tree_21tips.nwk`: Exactly 21 terminal tips matching `scMultiSim` population allocation vectors.
  * `arab_root_paga_tree_informative.nwk`: Detailed lineage labels with cell-type and cluster names.

### B. Rice (*Oryza sativa* - 8 Tissues)
- **Tissues**: Bud (8 clusters), Flag leaf (5 clusters), Young leaf (6 clusters), Root (9 clusters), Shoot Apical Meristem (SAM; 7 clusters), Seed (4 clusters), Spikelet/Panicle (SP; 9 clusters), Stem (ST; 6 clusters).
- **Progenitor Rooting**: Empirically pinned to stem/meristematic populations (e.g., Root Apical Meristem, SAM, Plumule, Procambium).
- **Deliverables per tissue**: `<tissue>_paga_tree_all_tips.nwk`, `<tissue>_paga_tree_classic.nwk`, `<tissue>_paga_connectivity.csv`, `<tissue>_celltype_metadata.csv`, and `<tissue>_paga_processed.h5ad`.

### C. Soybean (*Glycine max* - 7 Tissues & Developmental Stages)
- **Tissues / Stages**: Cotyledon stage seeds, Early maturation stage seeds, Early nodule, Globular stage seeds, Heart stage seeds, Hypocotyl, Root.
- **Progenitor Rooting**: Zygotic/embryonic progenitor states and apical meristems.
- **Deliverables per tissue**: All tips Newick tree, classic tree, connectivity matrix, metadata, and AnnData.

---

## 3. How to Use in `scMultiSim`

In R, load the tree directly using `ape` and feed it into `scMultiSim_input`:

```r
library(scMultiSim)
library(ape)

# 1. Load the empirical Newick tree
tree <- read.tree("arab_root_paga_tree_21tips.nwk")

# 2. Confirm tip labels
print(tree$tip.label)

# 3. Define population sizes matching tips 1-to-1
meta <- read.csv("cluster_metadata_21tips.csv")
pop_size <- setNames(meta$n_cells, meta$cluster_id)[tree$tip.label]

# 4. Configure scMultiSim
sim_options <- list(
  rand.seed = 2026,
  GRN = grn_network,                # Inferred from plantFigR
  tree = tree,                      # Empirical PAGA lineage tree
  num.cells = sum(pop_size),
  discrete.pop.size = pop_size,
  diff.mode = "continuous"
)
```

---

## 4. Associated Source Code Repository

Source code and reproducible execution scripts are hosted on GitHub:
- **Repository**: [https://github.com/kabilanbio/plant-scMultiome-simulation](https://github.com/kabilanbio/plant-scMultiome-simulation)
