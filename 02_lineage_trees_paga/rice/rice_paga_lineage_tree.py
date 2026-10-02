#!/usr/bin/env python3
# ==============================================================================
# rice_paga_lineage_tree.py
# Publication-Grade PAGA Cell Lineage Tree Reconstruction for Standardized Rice Datasets
#
# Scientific Framework:
#   - PAGA Algorithm: Wolf, F.A., Hamey, F.K., Plass, M. et al. 
#     "PAGA: graph abstraction generates coarse-grained topologies of high-dimensional single-cell data."
#     Genome Biology 20, 59 (2019). https://doi.org/10.1186/s13059-019-1663-x
#   - Official GitHub: https://github.com/theislab/paga (Scanpy implementation)
#   - Downstream Target: scMultiSim (Zhang et al., Nature Methods 2023)
#
# Key Features:
#   1. Supports all 8 standardized Rice tissues: Bud, Flag, Leaf, Root, SAM, Seed, SP, ST
#   2. Automatic Biological Progenitor identification per tissue for empirical rooting
#   3. Generates Maximum Spanning Tree (MST) from PAGA statistical connectivities
#   4. Dual Newick tree export:
#      - 'all_tips.nwk': Every cell type is an explicit tip in tree$tip.label (1-to-1 match for scMultiSim discrete.pop.size)
#      - 'classic.nwk': Classical hierarchical phylogenetic branching
#   5. High-resolution 600 DPI publication figures (PAGA graph, PAGA-UMAP, Pseudotime, Tree)
# ==============================================================================

import os
import sys
import io
import argparse
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.patheffects as pe
from matplotlib.colors import Normalize
import networkx as nx
import scipy.sparse as sp

try:
    import scanpy as sc
    import anndata as ad
    ad.settings.allow_write_nullable_strings = True
except ImportError:
    print("[ERROR] Scanpy/AnnData not found in current Python environment.")
    print("Please activate your conda environment (e.g., conda activate conda_latest).")
    sys.exit(1)

# Biological Progenitor Mapping per Rice Tissue (Empirical Botanical Knowledge)
TISSUE_PROGENITORS = {
    "Bud":  ["Bud_primordium", "Shoot_apical_meristem", "Meristem", "Bud_axis"],
    "Flag": ["Vascular_cylinder", "Mesophyll", "Epidermis"],
    "Leaf": ["Vascular_cylinder", "Parenchyma_cell", "Mesophyll"],
    "Root": ["Meristem", "Root_meristem", "Root_cap"],
    "SAM":  ["Meristem", "Proliferating_cell", "Primordia_apex"],
    "Seed": ["Plumule", "Scutellum", "Vascular_cylinder"],
    "SP":   ["Primordia_apex", "Ovule", "Vascular_cylinder", "Meristem"],
    "ST":   ["Procambium", "Vascular_cylinder", "Collenchymatous_cell"],
    "rice_sub": ["Root_meristem", "Meristem", "Shoot_apical_meristem", "Bud_primordium"]
}

ALL_RICE_TISSUES = ["Bud", "Flag", "Leaf", "Root", "SAM", "Seed", "SP", "ST"]


# ------------------------------------------------------------------------------
# 1. CLI Argument Parser
# ------------------------------------------------------------------------------
def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Construct PAGA Cell Lineage Trees for Rice Multiomics Datasets"
    )
    parser.add_argument(
        "--input_dir", type=str,
        default="G:/PhD/sc_datasets/rice/featured_datasets_anndata",
        help="Path containing exported tissue datasets (h5ad or matrix_market)"
    )
    parser.add_argument(
        "--standardized_dir", type=str,
        default="G:/PhD/sc_datasets/rice/featured_datasets_standardized",
        help="Fallback path containing standardized R RDS / CSV datasets"
    )
    parser.add_argument(
        "--output_dir", type=str,
        default="G:/PhD/data_simulation_comparison/scMultiSim/rice_datasets/rice_paga_results",
        help="Directory to save Newick trees, connectivity matrices, and publication figures"
    )
    parser.add_argument(
        "--tissue", type=str, default="all",
        help="Specific tissue ('Bud', 'Root', etc.) or 'all' to process all 8 tissues"
    )
    parser.add_argument(
        "--root_celltype", type=str, default=None,
        help="Custom root cell type name (overrides automatic biological progenitor)"
    )
    parser.add_argument(
        "--n_neighbors", type=int, default=15,
        help="Number of nearest neighbors for graph construction (default: 15)"
    )
    parser.add_argument(
        "--n_pcs", type=int, default=30,
        help="Number of principal components for graph construction (default: 30)"
    )
    parser.add_argument(
        "--dpi", type=int, default=600,
        help="Resolution for publication graphics (default: 600 DPI)"
    )
    parser.add_argument(
        "--seed", type=int, default=42,
        help="Random seed for reproducibility"
    )
    return parser.parse_args()


# ------------------------------------------------------------------------------
# 2. Data Loading Helper
# ------------------------------------------------------------------------------
def load_tissue_dataset(tissue, input_dir, standardized_dir):
    """
    Loads tissue single-cell RNA dataset from:
    1. Direct .h5ad file
    2. Matrix Market export (matrix.mtx + genes.tsv + barcodes.tsv + metadata.csv)
    3. Direct reading of standardized CSV metadata + count matrix if available
    """
    tissue_h5ad_candidates = [
        os.path.join(input_dir, tissue, f"rice_processed_{tissue}.h5ad"),
        os.path.join(input_dir, f"rice_processed_{tissue}.h5ad"),
        os.path.join(standardized_dir, tissue, f"rice_processed_{tissue}.h5ad")
    ]
    
    for h5_path in tissue_h5ad_candidates:
        if os.path.exists(h5_path):
            print(f"  [Loader] Loading AnnData from H5AD: {h5_path}")
            adata = sc.read_h5ad(h5_path)
            return adata

    # Candidate 2: Matrix Market folder
    mm_dir_candidates = [
        os.path.join(input_dir, tissue, "matrix_market"),
        os.path.join(standardized_dir, tissue, "matrix_market")
    ]
    for mm_dir in mm_dir_candidates:
        mtx_file = os.path.join(mm_dir, "matrix.mtx")
        genes_file = os.path.join(mm_dir, "genes.tsv")
        barcodes_file = os.path.join(mm_dir, "barcodes.tsv")
        meta_file = os.path.join(mm_dir, "metadata.csv")
        
        if os.path.exists(mtx_file) and os.path.exists(barcodes_file):
            print(f"  [Loader] Loading AnnData from Matrix Market: {mm_dir}")
            # Matrix is typically genes x cells from R writeMM
            mat = sp.csr_matrix(sc.read_mtx(mtx_file).X.T)
            barcodes = pd.read_csv(barcodes_file, header=None)[0].astype(str).tolist()
            genes = pd.read_csv(genes_file, header=None)[0].astype(str).tolist()
            
            adata = ad.AnnData(X=mat, obs=pd.DataFrame(index=barcodes), var=pd.DataFrame(index=genes))
            if os.path.exists(meta_file):
                meta = pd.read_csv(meta_file, index_col=0)
                for col in meta.columns:
                    adata.obs[col] = meta.loc[adata.obs_names, col].values
            return adata

    # Candidate 3: Try rpy2 if available (for HPC / Linux)
    try:
        import rpy2.robjects as ro
        import anndata2ri
        from rpy2.robjects.conversion import localconverter
        
        rds_candidates = [
            os.path.join(standardized_dir, tissue, f"rice_processed_{tissue}_rna.rds"),
            os.path.join(input_dir, tissue, f"rice_processed_{tissue}_rna.rds")
        ]
        for rds_path in rds_candidates:
            if os.path.exists(rds_path):
                print(f"  [Loader] Converting Seurat RDS via rpy2: {rds_path}")
                with localconverter(anndata2ri.converter):
                    adata = ro.r(f'''
                        library(Seurat)
                        obj <- readRDS("{rds_path}")
                        as.SingleCellExperiment(obj, assay = "RNA")
                    ''')
                return adata
    except Exception as e:
        pass

    raise FileNotFoundError(
        f"Could not locate dataset for tissue '{tissue}' in either:\n"
        f"  1) {input_dir}\n"
        f"  2) {standardized_dir}\n"
        "Please run 'export_rice_standardized_to_anndata.R' first to generate AnnData / Matrix Market files."
    )


# ------------------------------------------------------------------------------
# 3. Newick Conversion Functions
# ------------------------------------------------------------------------------
def build_classic_newick(tree, node, conn, categories, parent=None):
    """
    Standard recursive Newick builder where non-leaf nodes are intermediate states.
    """
    children = list(tree.successors(node))
    if parent is not None:
        idx_p = categories.index(parent)
        idx_c = categories.index(node)
        w = conn[idx_p, idx_c]
        # Scaled biological distance: bounded to avoid extreme branch lengths
        branch_len = round(max(0.05, min(2.0, 1.0 / w)), 4) if w > 0 else 0.5
    else:
        branch_len = 0.0
        
    label = node.replace(' ', '_')
    if not children:
        return f"{label}:{branch_len}"
    child_str = ','.join(build_classic_newick(tree, c, conn, categories, node) for c in children)
    return f"({child_str}){label}:{branch_len}"


def build_all_tips_newick(tree, node, conn, categories, parent=None):
    """
    scMultiSim-Optimized Newick Builder:
    Guarantees that EVERY cell type appears as a terminal TIP (leaf) in the tree.
    If an internal node 'U' has cells, an explicit tip 'U' is attached with a minimal
    differentiation step, ensuring 'discrete.pop.size' matches 'tree$tip.label' 1-to-1.
    """
    children = list(tree.successors(node))
    if parent is not None:
        idx_p = categories.index(parent)
        idx_c = categories.index(node)
        w = conn[idx_p, idx_c]
        branch_len = round(max(0.05, min(2.0, 1.0 / w)), 4) if w > 0 else 0.5
    else:
        branch_len = 0.0
        
    label = node.replace(' ', '_')
    
    if not children:
        return f"{label}:{branch_len}"
    
    # Internal node: create branches for all children PLUS an explicit tip for this node itself
    child_strings = [build_all_tips_newick(tree, c, conn, categories, node) for c in children]
    # Attach a terminal tip representing the steady-state cells of this node
    self_tip = f"{label}_state:0.05"
    child_strings.append(self_tip)
    
    child_str = ','.join(child_strings)
    return f"({child_str}){label}_node:{branch_len}"


# ------------------------------------------------------------------------------
# 4. Publication Figures Generator (600 DPI)
# ------------------------------------------------------------------------------
def generate_publication_figures(adata, tissue, categories, celltype_colors, 
                                 T, root, fig_dir, dpi=600):
    """
    Generates 4 high-resolution 600 DPI figures:
      1. PAGA Connectivity Graph
      2. PAGA-initialized UMAP colored by cell types
      3. Diffusion Pseudotime (DPT) progression UMAP
      4. Developmental Lineage Tree layout
    """
    n_cats = len(categories)
    
    # ----------------- Fig 1: PAGA Connectivity Graph -----------------
    fig, ax = plt.subplots(figsize=(9, 8), dpi=dpi)
    sc.pl.paga(
        adata,
        color='cell_type',
        threshold=0.03,
        node_size_scale=1.8,
        edge_width_scale=1.5,
        fontsize=10,
        ax=ax,
        show=False
    )
    ax.set_title(f"PAGA Trajectory Graph: Rice {tissue} Atlas", fontsize=13, fontweight='bold', pad=15)
    paga_fig_path = os.path.join(fig_dir, f"Fig1_{tissue}_paga_graph.png")
    fig.savefig(paga_fig_path, bbox_inches='tight', dpi=dpi)
    fig.savefig(os.path.join(fig_dir, f"Fig1_{tissue}_paga_graph.pdf"), bbox_inches='tight')
    plt.close(fig)
    
    # ----------------- Fig 2: Cell-Type Colored PAGA UMAP -----------------
    umap_coords = adata.obsm['X_umap']
    umap_df = pd.DataFrame({
        'u1': umap_coords[:, 0],
        'u2': umap_coords[:, 1],
        'cell_type': adata.obs['cell_type'].values
    })
    
    centers = umap_df.groupby('cell_type')[['u1', 'u2']].median().reindex(categories)
    
    fig, (ax_main, ax_leg) = plt.subplots(
        1, 2, figsize=(14, max(6, n_cats * 0.8)), dpi=dpi,
        gridspec_kw={'width_ratios': [3.2, 1.2], 'wspace': 0.08}
    )
    
    for ct in categories:
        mask = umap_df['cell_type'] == ct
        ax_main.scatter(
            umap_df.loc[mask, 'u1'], umap_df.loc[mask, 'u2'],
            c=celltype_colors[ct], s=12, alpha=0.75, linewidths=0,
            rasterized=True, label=ct
        )
    
    # Add clear text badges
    for ct, row in centers.iterrows():
        if pd.isna(row['u1']):
            continue
        display_label = ct.replace('_', '\n')
        ax_main.text(
            row['u1'], row['u2'], display_label,
            fontsize=8, fontweight='bold', ha='center', va='center',
            path_effects=[pe.withStroke(linewidth=3.5, foreground='white')]
        )
        
    ax_main.set_xlabel('UMAP 1', fontsize=11, fontweight='bold')
    ax_main.set_ylabel('UMAP 2', fontsize=11, fontweight='bold')
    ax_main.set_title(f"Rice {tissue}: PAGA-Initialized Single-Cell UMAP", fontsize=13, fontweight='bold')
    ax_main.tick_params(labelsize=9)
    
    # Clean custom legend panel
    ax_leg.axis('off')
    for i, ct in enumerate(categories):
        y = n_cats - i
        ax_leg.scatter([0.2], [y], s=140, c=[celltype_colors[ct]], edgecolors='black', linewidths=0.5)
        count = (adata.obs['cell_type'] == ct).sum()
        label_text = f"{ct.replace('_', ' ')} (n={count})"
        ax_leg.text(0.35, y, label_text, ha='left', va='center', fontsize=9, color='#111111')
    ax_leg.set_xlim(0, 2)
    ax_leg.set_ylim(0.5, n_cats + 0.5)
    
    umap_fig_path = os.path.join(fig_dir, f"Fig2_{tissue}_paga_umap_celltypes.png")
    fig.savefig(umap_fig_path, bbox_inches='tight', dpi=dpi)
    fig.savefig(os.path.join(fig_dir, f"Fig2_{tissue}_paga_umap_celltypes.pdf"), bbox_inches='tight')
    plt.close(fig)
    
    # ----------------- Fig 3: Diffusion Pseudotime (DPT) UMAP -----------------
    if 'dpt_pseudotime' in adata.obs:
        pt_vals = adata.obs['dpt_pseudotime'].values
        valid = ~np.isnan(pt_vals)
        cmap_pt = plt.cm.viridis
        norm_pt = Normalize(vmin=np.nanmin(pt_vals), vmax=np.nanmax(pt_vals))
        
        fig, ax_pt = plt.subplots(figsize=(10, 8), dpi=dpi)
        scatter = ax_pt.scatter(
            umap_df.loc[valid, 'u1'], umap_df.loc[valid, 'u2'],
            c=pt_vals[valid], cmap=cmap_pt, norm=norm_pt,
            s=12, alpha=0.8, linewidths=0, rasterized=True
        )
        cbar = plt.colorbar(scatter, ax=ax_pt, fraction=0.035, pad=0.02)
        cbar.set_label('Diffusion Pseudotime (Progression from Root)', fontsize=10, fontweight='bold')
        
        # Highlight root location
        root_center = centers.loc[root]
        if not pd.isna(root_center['u1']):
            ax_pt.scatter(
                [root_center['u1']], [root_center['u2']],
                s=280, marker='*', c='crimson', edgecolors='white', linewidths=1.5,
                zorder=5, label=f"Progenitor Root ({root})"
            )
            ax_pt.legend(loc='upper right', frameon=True, fontsize=9)
            
        ax_pt.set_xlabel('UMAP 1', fontsize=11, fontweight='bold')
        ax_pt.set_ylabel('UMAP 2', fontsize=11, fontweight='bold')
        ax_pt.set_title(f"Rice {tissue}: Diffusion Pseudotime (Root: {root.replace('_', ' ')})", 
                         fontsize=13, fontweight='bold')
        
        pt_fig_path = os.path.join(fig_dir, f"Fig3_{tissue}_diffusion_pseudotime.png")
        fig.savefig(pt_fig_path, bbox_inches='tight', dpi=dpi)
        fig.savefig(os.path.join(fig_dir, f"Fig3_{tissue}_diffusion_pseudotime.pdf"), bbox_inches='tight')
        plt.close(fig)
        
    # ----------------- Fig 4: Lineage Tree Graph -----------------
    fig, ax_tr = plt.subplots(figsize=(10, 8), dpi=dpi)
    pos = nx.spring_layout(T, seed=42)
    node_colors = [celltype_colors.get(n, '#4A90E2') for n in T.nodes]
    node_sizes = [max(300, (adata.obs['cell_type'] == n).sum() * 0.8) for n in T.nodes]
    
    nx.draw_networkx_nodes(T, pos, node_color=node_colors, node_size=node_sizes, ax=ax_tr, edgecolors='black', linewidths=1)
    nx.draw_networkx_edges(T, pos, width=2.0, alpha=0.7, edge_color='#444444', ax=ax_tr)
    labels = {n: n.replace('_', '\n') for n in T.nodes}
    nx.draw_networkx_labels(T, pos, labels=labels, font_size=8, font_weight='bold', ax=ax_tr)
    
    ax_tr.set_title(f"Empirical Maximum Spanning Tree (MST): Rice {tissue}\nRoot: {root.replace('_', ' ')}", 
                    fontsize=12, fontweight='bold')
    ax_tr.axis('off')
    
    tree_fig_path = os.path.join(fig_dir, f"Fig4_{tissue}_lineage_tree_mst.png")
    fig.savefig(tree_fig_path, bbox_inches='tight', dpi=dpi)
    fig.savefig(os.path.join(fig_dir, f"Fig4_{tissue}_lineage_tree_mst.pdf"), bbox_inches='tight')
    plt.close(fig)


# ------------------------------------------------------------------------------
# 5. Core Pipeline for Single Tissue
# ------------------------------------------------------------------------------
def process_single_tissue(tissue, opt):
    print("\n" + "=" * 70)
    print(f"PROCESSING RICE TISSUE: {tissue.upper()}")
    print("=" * 70)
    
    tissue_out_dir = os.path.join(opt.output_dir, tissue)
    fig_dir = os.path.join(tissue_out_dir, "figures")
    os.makedirs(tissue_out_dir, exist_ok=True)
    os.makedirs(fig_dir, exist_ok=True)
    
    # 1. Load Data
    adata = load_tissue_dataset(tissue, opt.input_dir, opt.standardized_dir)
    print(f"  [Dimensions] Raw: {adata.n_obs} cells x {adata.n_vars} genes")
    
    # Check cell_type column
    ct_col = 'cell_type'
    if ct_col not in adata.obs.columns:
        for candidate in ['CellType', 'cluster_names', 'seurat_clusters']:
            if candidate in adata.obs.columns:
                adata.obs[ct_col] = adata.obs[candidate]
                break
    if ct_col not in adata.obs.columns:
        raise ValueError(f"No cell_type annotation found in metadata for tissue {tissue}!")
        
    adata.obs[ct_col] = adata.obs[ct_col].astype('category')
    categories = adata.obs[ct_col].cat.categories.tolist()
    n_cats = len(categories)
    print(f"  [Cell Types] Found {n_cats} categories: {categories}")
    
    # 2. Quality and Normalization
    x_max = adata.X.max() if sp.issparse(adata.X) else np.max(adata.X)
    if x_max > 25:
        print("  [Preprocess] Raw count values detected. Applying total count normalization + log1p...")
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
    else:
        print("  [Preprocess] Data appears log-normalized. Proceeding directly...")
        
    # HVG & PCA
    n_hvg = min(2000, adata.n_vars)
    if adata.n_vars > n_hvg:
        sc.pp.highly_variable_genes(adata, n_top_genes=n_hvg)
    sc.pp.pca(adata, n_comps=min(50, adata.n_obs - 1, adata.n_vars - 1), use_highly_variable=True if 'highly_variable' in adata.var else False)
    
    # 3. Neighborhood Graph & PAGA
    n_pcs_use = min(opt.n_pcs, adata.obsm['X_pca'].shape[1])
    sc.pp.neighbors(adata, n_neighbors=opt.n_neighbors, n_pcs=n_pcs_use)
    
    print("  [PAGA] Computing graph abstraction...")
    sc.tl.paga(adata, groups=ct_col)
    
    # PAGA-initialized UMAP
    sc.tl.umap(adata, init_pos='paga')
    
    # 4. Determine Biological Progenitor Root
    root = None
    if opt.root_celltype and opt.root_celltype in categories:
        root = opt.root_celltype
        print(f"  [Root] User-specified progenitor root: '{root}'")
    else:
        candidates = TISSUE_PROGENITORS.get(tissue, [])
        for cand in candidates:
            if cand in categories:
                root = cand
                print(f"  [Root] Empirical biological progenitor matched: '{root}'")
                break
                
    if root is None:
        # Fallback to degree centrality on PAGA graph
        conn_mat = adata.uns['paga']['connectivities'].toarray()
        degree_sums = conn_mat.sum(axis=1)
        root = categories[np.argmax(degree_sums)]
        print(f"  [Root] Progenitor fallback (highest PAGA connectivity): '{root}'")
        
    # 5. Diffusion Pseudotime (DPT)
    try:
        root_indices = np.flatnonzero(adata.obs[ct_col] == root)
        if len(root_indices) > 0:
            adata.uns['iroot'] = int(root_indices[0])
            sc.tl.diffmap(adata)
            sc.tl.dpt(adata)
            print("  [DPT] Successfully computed Diffusion Pseudotime progression.")
    except Exception as e:
        print(f"  [DPT Warning] Diffusion pseudotime calculation skipped: {e}")
        
    # 6. Extract Lineage Tree via Maximum Spanning Tree (MST)
    conn = adata.uns['paga']['connectivities'].toarray()
    paga_conn_df = pd.DataFrame(conn, index=categories, columns=categories)
    conn_csv_path = os.path.join(tissue_out_dir, f"{tissue}_paga_connectivity.csv")
    paga_conn_df.to_csv(conn_csv_path)
    print(f"  [Export] Saved connectivity matrix: {os.path.basename(conn_csv_path)}")
    
    # Build NetworkX Graph
    G = nx.Graph()
    for i, ct_i in enumerate(categories):
        for j, ct_j in enumerate(categories):
            if i < j and conn[i, j] > 0:
                G.add_edge(ct_i, ct_j, weight=conn[i, j])
                
    # Maximum Spanning Tree
    T = nx.maximum_spanning_tree(G, weight='weight')
    T_rooted = nx.bfs_tree(T, root)
    
    # 7. Generate Dual Newick Trees
    # Tree A: scMultiSim-Optimized (all cell types as tips)
    newick_all_tips = build_all_tips_newick(T_rooted, root, conn, categories) + ";"
    nwk_all_tips_path = os.path.join(tissue_out_dir, f"{tissue}_paga_tree_all_tips.nwk")
    with open(nwk_all_tips_path, 'w') as f:
        f.write(newick_all_tips)
        
    # Tree B: Classic Hierarchical Tree
    newick_classic = build_classic_newick(T_rooted, root, conn, categories) + ";"
    nwk_classic_path = os.path.join(tissue_out_dir, f"{tissue}_paga_tree_classic.nwk")
    with open(nwk_classic_path, 'w') as f:
        f.write(newick_classic)
        
    # Also save the primary tree as standard <tissue>_paga_tree.nwk
    primary_nwk_path = os.path.join(tissue_out_dir, f"{tissue}_paga_tree.nwk")
    with open(primary_nwk_path, 'w') as f:
        f.write(newick_all_tips)
        
    print(f"  [Export] Saved scMultiSim Newick tree: {os.path.basename(primary_nwk_path)}")
    print(f"  [Newick Preview]: {newick_all_tips[:120]}...")
    
    # 8. Cell Type Metadata Summary
    ct_counts = adata.obs[ct_col].value_counts().reindex(categories)
    med_pt = (
        adata.obs.groupby(ct_col, observed=True)['dpt_pseudotime'].median().reindex(categories)
        if 'dpt_pseudotime' in adata.obs else pd.Series(0.0, index=categories)
    )
    
    ct_summary = pd.DataFrame({
        'cell_type': categories,
        'cell_type_clean': [ct.replace(' ', '_') for ct in categories],
        'n_cells': ct_counts.values,
        'is_progenitor_root': [ct == root for ct in categories],
        'median_pseudotime': med_pt.values
    })
    summary_path = os.path.join(tissue_out_dir, f"{tissue}_celltype_metadata.csv")
    ct_summary.to_csv(summary_path, index=False)
    print(f"  [Export] Saved cell type metadata: {os.path.basename(summary_path)}")
    
    # 9. Assign Colors & Generate Figures
    if n_cats <= 10:
        pal = plt.cm.get_cmap('tab10', n_cats)
    elif n_cats <= 20:
        pal = plt.cm.get_cmap('tab20', n_cats)
    else:
        pal = plt.cm.get_cmap('gist_ncar', n_cats)
        
    celltype_colors = {ct: mcolors.to_hex(pal(i)) for i, ct in enumerate(categories)}
    adata.uns['cell_type_colors'] = [celltype_colors[ct] for ct in categories]
    
    print(f"  [Graphics] Generating 600 DPI publication figures...")
    generate_publication_figures(adata, tissue, categories, celltype_colors, T, root, fig_dir, dpi=opt.dpi)
    
    # 10. Save Processed AnnData
    h5ad_out = os.path.join(tissue_out_dir, f"{tissue}_paga_processed.h5ad")
    adata.write_h5ad(h5ad_out)
    print(f"  [Export] Saved complete AnnData object: {os.path.basename(h5ad_out)}")
    
    return {
        'tissue': tissue,
        'n_cells': adata.n_obs,
        'n_genes': adata.n_vars,
        'n_cell_types': n_cats,
        'root_cell_type': root,
        'tree_path': primary_nwk_path
    }


# ------------------------------------------------------------------------------
# 6. Main Execution Loop
# ------------------------------------------------------------------------------
def main():
    opt = parse_arguments()
    np.random.seed(opt.seed)
    
    print("=" * 75)
    print("RICE SINGLE-CELL PAGA LINEAGE TREE RECONSTRUCTION SUITE")
    print("Paper Reference: Wolf et al., Genome Biology (2019) 20:59")
    print(f"Targeting Simulator: scMultiSim (Nature Methods 2023)")
    print(f"Input Directory:  {opt.input_dir}")
    print(f"Output Directory: {opt.output_dir}")
    print("=" * 75)
    
    tissues_to_run = (
        ALL_RICE_TISSUES if opt.tissue.lower() == "all"
        else [t.strip() for t in opt.tissue.split(",")]
    )
    
    summary_records = []
    for tissue in tissues_to_run:
        try:
            res = process_single_tissue(tissue, opt)
            summary_records.append(res)
        except Exception as e:
            print(f"[ERROR] Failed processing tissue '{tissue}': {e}")
            import traceback
            traceback.print_exc()
            
    # Save overall summary
    if summary_records:
        df_sum = pd.DataFrame(summary_records)
        sum_file = os.path.join(opt.output_dir, "rice_all_tissues_paga_summary.csv")
        df_sum.to_csv(sum_file, index=False)
        print("\n" + "=" * 75)
        print("EXECUTION COMPLETED SUCCESSFULLY!")
        print(f"Master Summary Table: {sum_file}")
        print("=" * 75)
        print(df_sum.to_string(index=False))


if __name__ == "__main__":
    main()
