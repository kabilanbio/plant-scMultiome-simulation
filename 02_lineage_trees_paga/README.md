# Empirical PAGA Lineage Tree Reconstruction

PAGA-based cell lineage tree inference pipelines for single-cell multiomics datasets:
* **`arabidopsis/`**: 21 developmental clusters mapped to 6 canonical cell types, rooted at Root Apical Meristematic cells.
* **`rice/`**: 8 tissues with empirical progenitor rooting (e.g., RAM, SAM, Plumule, Procambium).
* **`soybean/`**: 7 tissues and seed developmental stages.

### Tree Formats:
* `*_paga_tree_all_tips.nwk` (or `arab_root_paga_tree_21tips.nwk`): Guaranteed 100% terminal tip matching for `scMultiSim` population allocation vectors (`discrete.pop.size`).
* `*_paga_tree_classic.nwk`: Standard hierarchical phylogenetic tree.
