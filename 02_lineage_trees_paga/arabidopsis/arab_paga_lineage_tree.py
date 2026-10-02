#!/usr/bin/env python3
# ==============================================================================
# arab_paga_lineage_tree.py
# Publication-Grade PAGA Cell Lineage Tree Reconstruction for Arabidopsis Root
# Incorporates Curated Botanical Annotations from Farmer et al., Mol Plant (2021):
#   - 21 Developmental Clusters (1-21) -> 6 Biological Cell Types -> 4 Tissue Layers
#
# Scientific Framework:
#   - PAGA Algorithm: Wolf et al., Genome Biology (2019) 20:59
#   - Target Simulator: scMultiSim (Nature Methods 2023)
#   - Optimized for HP 280 G4 Workstation (i7-8700, 16 GB RAM, sc_env)
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
    print("Please activate your environment: C:\\Users\\kabil\\sc_env\\Scripts\\activate")
    sys.exit(1)

# Curated Biological Mapping (GSE155304; Farmer et al., 2021)
CLUSTER_TO_CELLTYPE = {
    'k0_11': 'Cortical cell',        'o1_15': 'Endodermal cell',  
    'e2_5':  'Atrichoblast',         'r3_18': 'Stele cell',  
    'g4_7':  'Atrichoblast',         'b5_2':  'Trichoblast',  
    'h6_8':  'Meristematic cell',    'q7_17': 'Stele cell',  
    'f8_6':  'Atrichoblast',         'i9_9':  'Meristematic cell',  
    'l10_12':'Cortical cell',        'm11_13':'Endodermal cell',  
    'd12_4': 'Atrichoblast',         'c13_3': 'Trichoblast',  
    'a14_1': 'Trichoblast',          'u15_21':'Stele cell',  
    'p16_16':'Endodermal cell',      'n17_14':'Endodermal cell',  
    's18_19':'Stele cell',           't19_20':'Stele cell',  
    'j20_10':'Meristematic cell'
}

CELLTYPE_TO_LAYER = {
    'Trichoblast':       'Epidermis',
    'Atrichoblast':      'Epidermis',
    'Meristematic cell': 'Meristem',
    'Cortical cell':     'Ground tissue',
    'Endodermal cell':   'Ground tissue',
    'Stele cell':        'Vasculature'
}

CELLTYPE_COLORS = {
    'Trichoblast':       '#7BA6A6',
    'Atrichoblast':      '#8BC34A',
    'Meristematic cell': '#757575',
    'Cortical cell':     '#B388FF',
    'Endodermal cell':   '#9C27B0',
    'Stele cell':        '#D8A15D'
}

CELLTYPE_ORDER = [
    'Trichoblast',
    'Atrichoblast',
    'Meristematic cell',
    'Cortical cell',
    'Endodermal cell',
    'Stele cell'
]

LAYER_GROUPS = [
    {'name': 'Epidermis',     'celltypes': ['Trichoblast', 'Atrichoblast'],      'color': '#7BA6A6'},
    {'name': 'Meristem',      'celltypes': ['Meristematic cell'],                'color': '#757575'},
    {'name': 'Ground tissue', 'celltypes': ['Cortical cell', 'Endodermal cell'], 'color': '#B388FF'},
    {'name': 'Vasculature',   'celltypes': ['Stele cell'],                       'color': '#D8A15D'}
]


def parse_arguments():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    drive_prefix = os.path.splitdrive(script_dir)[0]
    if not drive_prefix:
        drive_prefix = "G:"

    parser = argparse.ArgumentParser(
        description="Construct Informative PAGA Lineage Trees for Arabidopsis Root"
    )
    parser.add_argument(
        "--input_dir", type=str,
        default=f"{drive_prefix}/PhD/sc_datasets/arabidopsis_root/featured_datasets_anndata",
        help="Path containing exported tissue datasets (h5ad or matrix_market)"
    )
    parser.add_argument(
        "--standardized_dir", type=str,
        default=f"{drive_prefix}/PhD/sc_datasets/arabidopsis_root/featured_datasets_standardized",
        help="Fallback path containing standardized R RDS / CSV datasets"
    )
    parser.add_argument(
        "--output_dir", type=str,
        default=f"{drive_prefix}/PhD/data_simulation_comparison/scMultiSim/arab_root/arab_paga_results",
        help="Directory to save Newick trees, connectivity matrices, and publication figures"
    )
    parser.add_argument(
        "--n_neighbors", type=int, default=15,
        help="Number of nearest neighbors (default: 15)"
    )
    parser.add_argument(
        "--n_pcs", type=int, default=30,
        help="Number of principal components (default: 30)"
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


def load_dataset(input_dir, standardized_dir):
    candidates_h5 = [
        os.path.join(input_dir, "Root", "arab_processed_Root.h5ad"),
        os.path.join(input_dir, "arab_processed_Root.h5ad"),
        os.path.join(standardized_dir, "Root", "arab_processed_Root.h5ad"),
    ]
    for h5 in candidates_h5:
        if os.path.exists(h5):
            print(f"  [Loader] Loading AnnData from H5AD: {h5}")
            return sc.read_h5ad(h5)

    candidates_mm = [
        os.path.join(input_dir, "Root", "matrix_market"),
        os.path.join(input_dir, "matrix_market"),
        os.path.join(standardized_dir, "Root", "matrix_market")
    ]
    for mm_dir in candidates_mm:
        mtx = os.path.join(mm_dir, "matrix.mtx")
        barcodes_file = os.path.join(mm_dir, "barcodes.tsv")
        genes_file = os.path.join(mm_dir, "genes.tsv")
        meta_file = os.path.join(mm_dir, "metadata.csv")
        if os.path.exists(mtx) and os.path.exists(barcodes_file):
            print(f"  [Loader] Loading AnnData from Matrix Market: {mm_dir}")
            mat = sp.csr_matrix(sc.read_mtx(mtx).X.T)
            barcodes = pd.read_csv(barcodes_file, header=None)[0].astype(str).tolist()
            genes = pd.read_csv(genes_file, header=None)[0].astype(str).tolist()
            adata = ad.AnnData(X=mat, obs=pd.DataFrame(index=barcodes), var=pd.DataFrame(index=genes))
            if os.path.exists(meta_file):
                meta = pd.read_csv(meta_file, index_col=0)
                for c in meta.columns:
                    adata.obs[c] = meta.loc[adata.obs_names, c].values
            return adata

    raise FileNotFoundError(
        f"Could not locate dataset in:\n  1) {input_dir}\n  2) {standardized_dir}\n"
        "Please run 'export_arab_standardized_to_anndata.R' first."
    )


def build_two_level_newick(ct_tree, node, num_to_ct, conn, categories, use_informative_names=False):
    """
    Constructs the 2-level hierarchical tree:
    - Internal nodes: 6 biological cell types (Meristematic_cell, Stele_cell, etc.)
    - Leaves/Tips: 21 developmental clusters (cluster1..cluster21 or informative labels)
    Guarantees exactly 21 terminal tips in tree$tip.label for scMultiSim discrete.pop.size!
    """
    node_label = str(node).replace(' ', '_')
    children = list(ct_tree.successors(node))
    clusters_here = [c for c, ct in num_to_ct.items() if ct == node]
    clusters_here.sort()

    if use_informative_names:
        leaf_strings = [f"{node_label}_c{c}:0.1000" for c in clusters_here]
    else:
        leaf_strings = [f"cluster{c}:0.1000" for c in clusters_here]

    if not children:
        return f"({','.join(leaf_strings)}){node_label}"
    else:
        child_strings = []
        for c in children:
            idx_u = categories.index(node)
            idx_v = categories.index(c)
            w = conn[idx_u, idx_v]
            branch_len = round(max(0.05, min(2.0, 1.0 / w)), 4) if w > 0 else 0.5
            child_str = build_two_level_newick(ct_tree, c, num_to_ct, conn, categories, use_informative_names)
            child_strings.append(f"{child_str}:{branch_len}")

        all_children = child_strings + leaf_strings
        return f"({','.join(all_children)}){node_label}"


def generate_publication_figures(adata, fig_dir, dpi=600):
    os.makedirs(fig_dir, exist_ok=True)
    categories = CELLTYPE_ORDER

    # ----------------- Fig 1: PAGA Graph (6 Cell Types) -----------------
    fig, ax = plt.subplots(figsize=(9, 8), dpi=dpi)
    sc.pl.paga(
        adata, color='cell_type', threshold=0.03,
        node_size_scale=1.8, edge_width_scale=1.5, fontsize=10,
        ax=ax, show=False
    )
    ax.set_title("PAGA Trajectory: Arabidopsis Root (6 Cell Types)", fontsize=13, fontweight='bold', pad=15)
    fig.savefig(os.path.join(fig_dir, "Fig1_Arab_Root_paga_graph.png"), bbox_inches='tight', dpi=dpi)
    fig.savefig(os.path.join(fig_dir, "Fig1_Arab_Root_paga_graph.pdf"), bbox_inches='tight')
    plt.close(fig)

    # ----------------- Fig 2: UMAP Colored by Cell Type (with Layer Grouped Legend) -----------------
    umap_coords = adata.obsm['X_umap']
    umap_df = pd.DataFrame({
        'u1': umap_coords[:, 0], 'u2': umap_coords[:, 1],
        'cell_type': adata.obs['cell_type'].values,
        'cluster_num': adata.obs['cluster_num'].values
    })
    centers = umap_df.groupby('cell_type')[['u1', 'u2']].median().reindex(categories)

    fig, (ax_main, ax_leg) = plt.subplots(
        1, 2, figsize=(15, 8), dpi=dpi,
        gridspec_kw={'width_ratios': [3.2, 1.3], 'wspace': 0.08}
    )

    for ct in categories:
        mask = umap_df['cell_type'] == ct
        ax_main.scatter(
            umap_df.loc[mask, 'u1'], umap_df.loc[mask, 'u2'],
            c=CELLTYPE_COLORS[ct], s=8, alpha=0.75, linewidths=0, rasterized=True
        )

    for ct, row in centers.iterrows():
        if pd.isna(row['u1']): continue
        ax_main.text(
            row['u1'], row['u2'], ct.replace(' ', '\n'),
            fontsize=9, fontweight='bold', ha='center', va='center',
            path_effects=[pe.withStroke(linewidth=3.5, foreground='white')]
        )

    ax_main.set_xlabel('UMAP 1', fontsize=11, fontweight='bold')
    ax_main.set_ylabel('UMAP 2', fontsize=11, fontweight='bold')
    ax_main.set_title("Arabidopsis Root: PAGA Single-Cell UMAP (6 Cell Types)", fontsize=13, fontweight='bold')
    ax_main.tick_params(labelsize=9)

    # Multi-group legend panel
    ax_leg.axis('off')
    y = 12.0
    for grp in LAYER_GROUPS:
        ax_leg.text(0.1, y, grp['name'].upper(), fontsize=10, fontweight='bold', color=grp['color'])
        y -= 0.8
        for ct in grp['celltypes']:
            count = (adata.obs['cell_type'] == ct).sum()
            ax_leg.scatter([0.2], [y], s=120, c=[CELLTYPE_COLORS[ct]], edgecolors='black', linewidths=0.5)
            ax_leg.text(0.35, y, f"{ct} (n={count})", ha='left', va='center', fontsize=9, color='#111111')
            y -= 0.9
        y -= 0.4
    ax_leg.set_xlim(0, 2)
    ax_leg.set_ylim(0, 13)

    fig.savefig(os.path.join(fig_dir, "Fig2_Arab_Root_paga_umap_celltypes.png"), bbox_inches='tight', dpi=dpi)
    fig.savefig(os.path.join(fig_dir, "Fig2_Arab_Root_paga_umap_celltypes.pdf"), bbox_inches='tight')
    plt.close(fig)

    # ----------------- Fig 3: UMAP with 21 Cluster Numbers -----------------
    centers_21 = umap_df.groupby('cluster_num')[['u1', 'u2']].median()
    fig, ax_clust = plt.subplots(figsize=(10, 8), dpi=dpi)
    for ct in categories:
        mask = umap_df['cell_type'] == ct
        ax_clust.scatter(
            umap_df.loc[mask, 'u1'], umap_df.loc[mask, 'u2'],
            c=CELLTYPE_COLORS[ct], s=8, alpha=0.6, linewidths=0, rasterized=True
        )
    for c_num, row in centers_21.iterrows():
        ax_clust.text(
            row['u1'], row['u2'], str(c_num),
            fontsize=10, fontweight='bold', ha='center', va='center', color='black',
            path_effects=[pe.withStroke(linewidth=3.5, foreground='white')]
        )
    ax_clust.set_xlabel('UMAP 1', fontsize=11, fontweight='bold')
    ax_clust.set_ylabel('UMAP 2', fontsize=11, fontweight='bold')
    ax_clust.set_title("Arabidopsis Root: 21 Developmental Clusters Along Lineage Coordinates", fontsize=13, fontweight='bold')
    fig.savefig(os.path.join(fig_dir, "Fig3_Arab_Root_21clusters_umap.png"), bbox_inches='tight', dpi=dpi)
    fig.savefig(os.path.join(fig_dir, "Fig3_Arab_Root_21clusters_umap.pdf"), bbox_inches='tight')
    plt.close(fig)

    # ----------------- Fig 4: Diffusion Pseudotime UMAP -----------------
    if 'dpt_pseudotime' in adata.obs:
        pt_vals = adata.obs['dpt_pseudotime'].values
        valid = ~np.isnan(pt_vals)
        norm_pt = Normalize(vmin=np.nanmin(pt_vals), vmax=np.nanmax(pt_vals))

        fig, ax_pt = plt.subplots(figsize=(10, 8), dpi=dpi)
        scatter = ax_pt.scatter(
            umap_df.loc[valid, 'u1'], umap_df.loc[valid, 'u2'],
            c=pt_vals[valid], cmap=plt.cm.viridis, norm=norm_pt,
            s=8, alpha=0.8, linewidths=0, rasterized=True
        )
        cbar = plt.colorbar(scatter, ax=ax_pt, fraction=0.035, pad=0.02)
        cbar.set_label('Diffusion Pseudotime (Root: Meristematic cell)', fontsize=10, fontweight='bold')

        # Highlight Meristem center
        meristem_center = centers.loc['Meristematic cell']
        if not pd.isna(meristem_center['u1']):
            ax_pt.scatter(
                [meristem_center['u1']], [meristem_center['u2']],
                s=280, marker='*', c='crimson', edgecolors='white', linewidths=1.5,
                zorder=5, label='Stem Cell Progenitor Root'
            )
            ax_pt.legend(loc='upper right', frameon=True, fontsize=9)

        ax_pt.set_xlabel('UMAP 1', fontsize=11, fontweight='bold')
        ax_pt.set_ylabel('UMAP 2', fontsize=11, fontweight='bold')
        ax_pt.set_title("Arabidopsis Root: Diffusion Pseudotime Differentiation Trajectory", fontsize=13, fontweight='bold')
        fig.savefig(os.path.join(fig_dir, "Fig4_Arab_Root_diffusion_pseudotime.png"), bbox_inches='tight', dpi=dpi)
        fig.savefig(os.path.join(fig_dir, "Fig4_Arab_Root_diffusion_pseudotime.pdf"), bbox_inches='tight')
        plt.close(fig)


def main():
    opt = parse_arguments()
    np.random.seed(opt.seed)

    print("=" * 75)
    print("ARABIDOPSIS ROOT PAGA LINEAGE TREE RECONSTRUCTION (INFORMATIVE SUITE)")
    print("Biological Reference: Farmer et al., Molecular Plant (2021)")
    print("Target Simulator: scMultiSim (Nature Methods 2023)")
    print("=" * 75)

    os.makedirs(opt.output_dir, exist_ok=True)
    fig_dir = os.path.join(opt.output_dir, "figures")
    os.makedirs(fig_dir, exist_ok=True)

    adata = load_dataset(opt.input_dir, opt.standardized_dir)
    print(f"  [Dimensions] {adata.n_obs} cells x {adata.n_vars} genes")

    # Map curated annotations
    cluster_key = adata.obs['seurat_clusters_renamed'] if 'seurat_clusters_renamed' in adata.obs else adata.obs['cell_type']
    adata.obs['cell_type'] = cluster_key.map(CLUSTER_TO_CELLTYPE).astype('category')
    adata.obs['tissue_layer'] = adata.obs['cell_type'].map(CELLTYPE_TO_LAYER).astype('category')
    adata.obs['cluster_num'] = cluster_key.str.extract(r'_(\d+)$')[0].astype(int)
    adata.obs['cluster_tip'] = "cluster" + adata.obs['cluster_num'].astype(str)

    categories = CELLTYPE_ORDER
    adata.obs['cell_type'] = pd.Categorical(adata.obs['cell_type'], categories=categories)
    print(f"  [Biological Annotation] Successfully mapped 6 cell types across 4 tissue layers.")

    # Normalization & PCA
    x_max = adata.X.max() if sp.issparse(adata.X) else np.max(adata.X)
    if x_max > 25:
        print("  [Preprocess] Normalizing total counts + log1p...")
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)

    n_hvg = min(2000, adata.n_vars)
    if adata.n_vars > n_hvg:
        sc.pp.highly_variable_genes(adata, n_top_genes=n_hvg)
    sc.pp.pca(adata, n_comps=min(50, adata.n_obs - 1, adata.n_vars - 1), use_highly_variable=True if 'highly_variable' in adata.var else False)

    # Neighborhood Graph & PAGA on 6 Cell Types
    n_pcs_use = min(opt.n_pcs, adata.obsm['X_pca'].shape[1])
    sc.pp.neighbors(adata, n_neighbors=opt.n_neighbors, n_pcs=n_pcs_use)
    print("  [PAGA] Computing graph abstraction on 6 biological cell types...")
    sc.tl.paga(adata, groups='cell_type')
    sc.pl.paga(adata, show=False)
    sc.tl.umap(adata, init_pos='paga')

    # Diffusion Pseudotime rooted at Meristematic cell
    root_ct = 'Meristematic cell'
    meristem_indices = np.flatnonzero(adata.obs['cell_type'] == root_ct)
    adata.uns['iroot'] = int(meristem_indices[0])
    sc.tl.diffmap(adata)
    sc.tl.dpt(adata)
    print("  [DPT] Computed Diffusion Pseudotime rooted at Meristematic cell.")

    # Maximum Spanning Tree on Cell Types
    conn = adata.uns['paga']['connectivities'].toarray()
    paga_conn_df = pd.DataFrame(conn, index=categories, columns=categories)
    conn_csv_path = os.path.join(opt.output_dir, "arab_root_paga_connectivity.csv")
    paga_conn_df.to_csv(conn_csv_path)

    G = nx.Graph()
    for i, ct_i in enumerate(categories):
        for j, ct_j in enumerate(categories):
            if i < j and conn[i, j] > 0:
                G.add_edge(ct_i, ct_j, weight=conn[i, j])
    T_ct = nx.maximum_spanning_tree(G, weight='weight')
    T_ct_rooted = nx.bfs_tree(T_ct, root_ct)

    # Build 2-Level Trees (21 Tips for scMultiSim discrete.pop.size matching)
    num_to_ct = {}
    for seurat_id, ct in CLUSTER_TO_CELLTYPE.items():
        num = int(seurat_id.split('_')[1])
        num_to_ct[num] = ct

    newick_21tips = build_two_level_newick(T_ct_rooted, root_ct, num_to_ct, conn, categories, use_informative_names=False) + ":0.0;"
    newick_informative = build_two_level_newick(T_ct_rooted, root_ct, num_to_ct, conn, categories, use_informative_names=True) + ":0.0;"

    with open(os.path.join(opt.output_dir, "arab_root_paga_tree_21tips.nwk"), 'w') as f:
        f.write(newick_21tips)
    with open(os.path.join(opt.output_dir, "arab_root_paga_tree_informative.nwk"), 'w') as f:
        f.write(newick_informative)
    with open(os.path.join(opt.output_dir, "arab_root_paga_tree.nwk"), 'w') as f:
        f.write(newick_21tips)

    print(f"  [Export] Saved 21-tip Newick tree: arab_root_paga_tree_21tips.nwk")
    print(f"  [Preview]: {newick_21tips}")

    # Cluster Metadata Table (1-21)
    seurat_sorted = [k for k in sorted(CLUSTER_TO_CELLTYPE.keys(), key=lambda x: int(x.split('_')[1]))]
    med_pt_by_num = adata.obs.groupby('cluster_num')['dpt_pseudotime'].median()

    cluster_meta = pd.DataFrame({
        'tip_name':          [f'cluster{c}' for c in range(1, 22)],
        'cluster_num':       list(range(1, 22)),
        'seurat_id':         seurat_sorted,
        'cell_type':         [num_to_ct[c] for c in range(1, 22)],
        'tissue_layer':      [CELLTYPE_TO_LAYER[num_to_ct[c]] for c in range(1, 22)],
        'informative_label': [f"{num_to_ct[c].replace(' ', '_')}_c{c}" for c in range(1, 22)],
        'n_cells':           [int((adata.obs['cluster_num'] == c).sum()) for c in range(1, 22)],
        'median_pseudotime': [round(med_pt_by_num.get(c, 0.0), 4) for c in range(1, 22)]
    })
    meta_path = os.path.join(opt.output_dir, "cluster_metadata_21tips.csv")
    cluster_meta.to_csv(meta_path, index=False)
    print(f"  [Export] Saved cluster metadata: {os.path.basename(meta_path)}")

    # Visualizations
    print("  [Graphics] Generating 600 DPI publication figures with curated legends...")
    generate_publication_figures(adata, fig_dir, dpi=opt.dpi)

    h5_path = os.path.join(opt.output_dir, "arab_root_paga_processed.h5ad")
    adata.write_h5ad(h5_path)
    print(f"  [Export] Saved complete AnnData: {os.path.basename(h5_path)}")

    print("\n" + "=" * 75)
    print("ARABIDOPSIS ROOT PAGA LINEAGE TREE COMPLETED SUCCESSFULLY!")
    print(f"Output files located in: {opt.output_dir}")
    print("=" * 75)
    print(cluster_meta.to_string(index=False))


if __name__ == "__main__":
    main()
