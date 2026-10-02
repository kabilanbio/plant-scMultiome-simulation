# plant-scMultiome-simulation

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23068394.svg)](https://doi.org/10.5281/zenodo.23068394)

Scripts for single-cell multi-omics data preprocessing and empirical PAGA lineage tree reconstruction across plant species (*Oryza sativa*, *Glycine max*, and *Arabidopsis thaliana*).

---

## Repository Structure

```
plant-scMultiome-simulation/
├── 01_preprocessing/                          # Multi-omics standardization scripts
│   ├── preprocess_arabidopsis_root.R          # Arabidopsis Root atlas (21 developmental clusters)
│   ├── preprocess_rice_all_tissues.R          # Rice (8 tissues)
│   └── preprocess_soybean_all_tissues.R       # Soybean (7 tissues & developmental stages)
└── 02_lineage_trees_paga/                     # Empirical PAGA lineage tree inference
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

## Data Availability

All datasets, processed AnnData objects, and Newick lineage trees are available on Zenodo:
* **DOI**: [10.5281/zenodo.23068394](https://doi.org/10.5281/zenodo.23068394)

---

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
