# Empirical PAGA Cell Lineage Tree Reconstruction for Plant Single-Cell Multi-Omics

## 1. Scientific Framework & Publication Grounding

In single-cell simulation with **scMultiSim** ([Zhang et al., *Nature Methods* 2023](https://doi.org/10.1038/s41592-023-02051-0)), cell differentiation trajectories are modeled as continuous branching processes governed by phylogenetic cell lineage trees.

Rather than generating arbitrary synthetic trees, our framework computes empirical, data-driven lineage trees directly from single-cell transcriptomic geometries using **PAGA (Partition-based Graph Abstraction)**:
* **Theory & Methodology**: [Wolf, F.A., Hamey, F.K., Plass, M. et al. PAGA: graph abstraction generates coarse-grained topologies of high-dimensional single-cell data. *Genome Biology* 20, 59 (2019)](https://doi.org/10.1186/s13059-019-1663-x).
* **Official Implementation**: [Theis Lab PAGA Repository](https://github.com/theislab/paga) and `scanpy.tl.paga`.

PAGA computes the statistical connectivity $c_{ij}$ between cell clusters:
$$c_{ij} = \frac{e_{ij} - \bar{e}_{ij}}{e_{ij}^{\max} - \bar{e}_{ij}}$$
where $e_{ij}$ is the observed inter-cluster edge count in the $k$-NN graph and $\bar{e}_{ij}$ is the expected random edge frequency. 

From the connectivity graph $G$, the **Maximum Spanning Tree (MST)** $T = \text{MST}(G, \text{weight}=c_{ij})$ provides the developmental differentiation backbone.

---

## 2. Botanical Progenitor Rooting Across 3 Plant Systems

Rooting the differentiation graph requires biologically grounded initial states:

### A. Rice (*Oryza sativa* - 8 Tissues)
- **Bud**: Stem cell `Bud_primordium` / `Meristem`
- **Flag Leaf**: `Vascular_cylinder` / `Mesophyll`
- **Leaf**: `Vascular_cylinder` / `Parenchyma_cell`
- **Root**: Root Apical Meristem (`Meristem`)
- **SAM**: Shoot Apical Meristem (`Meristem` / `Proliferating_cell`)
- **Seed**: Embryonic axis (`Plumule` / `Scutellum`)
- **SP (Spikelet)**: Inflorescence stem state (`Primordia_apex` / `Ovule`)
- **ST (Stem)**: Vascular cambium progenitor (`Procambium`)

### B. Soybean (*Glycine max* - 7 Tissues & Developmental Stages)
- **Cotyledon / Heart / Globular / Early maturation**: Zygotic / embryonic progenitor states
- **Early nodule**: Symbiotic nitrogen-fixing meristematic tissue
- **Hypocotyl & Root**: Procambium / Root apical meristematic lineages

### C. Arabidopsis (*Arabidopsis thaliana* - Root Atlas, GSE155304)
- **21 developmental clusters** mapped to **6 canonical cell types**:
  * `Trichoblast` (c1–c3, Root hair epidermis)
  * `Atrichoblast` (c4–c7, Non-hair epidermis)
  * `Meristematic cell` (c8–c10, **Botanical Progenitor Root**, Root Apical Meristem / QC)
  * `Cortical cell` (c11–c12, Ground tissue cortex)
  * `Endodermal cell` (c13–c16, Ground tissue Casparian strip / endodermis)
  * `Stele cell` (c17–c21, Vascular cylinder / xylem / phloem / pericycle)

---

## 3. Resolving the `scMultiSim` Tip Matching Bottleneck

### The Critical Bottleneck:
Classic phylogenetic Newick trees format intermediate ancestral clusters as internal nodes:
```text
((c1:0.5, c2:0.5)c_parent:0.2)root:0.0;
```
Standard phylogenetic parsers (`ape::read.tree()` in R) assign only terminal leaves `c1` and `c2` to `tree$tip.label`, placing `c_parent` into `tree$node.label`. In `scMultiSim`:
* Discrete population sizes (`discrete.pop.size`) requires an exact 1-to-1 match with `tree$tip.label`.
* Any intermediate cell type relegated to `node.label` cannot be assigned cell populations, crashing the simulation or omitting key cell types.

### The Solution: Dual Newick Architecture
Our pipeline exports two complementary trees per tissue:
1. **`*_paga_tree_all_tips.nwk` (Default for scMultiSim)**:
   Attaches a zero-length terminal tip to every developmental node. Every cluster is an explicit tip in `tree$tip.label`, enabling 100% population allocation in `scMultiSim`.
2. **`*_paga_tree_classic.nwk`**:
   Traditional phylogenetic hierarchical branching tree for macro-topology visualization.

---

## 4. Repository Structure

```
02_lineage_trees_paga/
├── README.md                                    # This scientific documentation
├── arabidopsis/
│   ├── export_arab_standardized_to_anndata.R    # Seurat RDS to Matrix Market / AnnData
│   ├── arab_paga_lineage_tree.py                # 21-cluster PAGA, DPT, and Newick tree generator
│   └── run_all_arab_paga.bat                    # 1-Click execution script
├── rice/
│   ├── export_rice_standardized_to_anndata.R    # 8-tissue Seurat to AnnData converter
│   ├── rice_paga_lineage_tree.py                # 8-tissue PAGA pipeline
│   ├── run_all_rice_paga.bat                    # 1-Click Windows execution script
│   └── run_rice_paga_cluster.pbs                # HPC PBS job submission script
└── soybean/
    ├── export_soybean_standardized_to_anndata.R # 7-tissue Seurat to AnnData converter
    ├── soybean_paga_lineage_tree.py             # 7-tissue PAGA pipeline
    └── run_all_soybean_paga.bat                 # 1-Click execution script
```

---

## 5. Execution Instructions

Ensure Python 3.10+ with `scanpy`, `anndata`, `igraph`, and `leidenalg` is installed (e.g., via `conda env create -f ../environment.yml` or `pip install -r ../requirements.txt`).

```bash
# Arabidopsis Root:
cd arabidopsis
python arab_paga_lineage_tree.py --dpi 600

# Rice (all 8 tissues):
cd ../rice
python rice_paga_lineage_tree.py --tissue all --dpi 600

# Soybean (all 7 tissues):
cd ../soybean
python soybean_paga_lineage_tree.py --tissue all --dpi 600
```
