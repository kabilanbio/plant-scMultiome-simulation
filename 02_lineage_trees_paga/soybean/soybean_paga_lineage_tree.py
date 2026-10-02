#!/usr/bin/env python3
# ==============================================================================
# soybean_paga_lineage_tree.py
# Publication-Grade PAGA Cell Lineage Tree Reconstruction for Soybean Datasets
# Single-Cell Multi-Omics Benchmark (GSE270192; Zhang et al., 2025)
#
# Scientific Framework:
#   - PAGA Algorithm: Wolf et al., Genome Biology (2019) 20:59
#   - Target Simulator: scMultiSim (Nature Methods 2023)
#   - System Optimization: HP 280 G4 Workstation (i7-8700, 16 GB RAM, sc_env)
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

SOYBEAN_PROGENITORS = {
    "Root": ["RAM", "Root_meristem", "Root_cap", "Meristem"],
    "Hypocotyl": ["Procambium", "Pith", "Vascular_parenchyma"],
    "Early_nodule": ["Nodule_meritem", "Nodule_meristem", "Procambium"],
    "Cotyledon_stage_seeds": ["Emb_proper_innitials", "Emb_proper_initials", "Endosperm"],
    "Globular_stage_seeds": ["Endosperm", "SC_inner_integument", "SC_endothelium"],
    "Heart_stage_seeds": ["SC_inner_integument", "Endosperm", "SC_endothelium"],
    "Early_maturation_stage_seeds": ["Emb_vasculature", "Emb_parenchyma", "Endosperm"]
}

ALL_SOY_TISSUES = [
    "Cotyledon_stage_seeds",
    "Early_maturation_stage_seeds",
    "Early_nodule",
    "Globular_stage_seeds",
    "Heart_stage_seeds",
    "Hypocotyl",
    "Root"
]


def parse_arguments():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    drive_prefix = os.path.splitdrive(script_dir)[0]
    if not drive_prefix:
        drive_prefix = "G:"

    parser = argparse.ArgumentParser(
        description="Construct PAGA Lineage Trees for Soybean Multiomics Datasets"
    )
    parser.add_argument(
        "--input_dir", type=str,
        default=f"{drive_prefix}/PhD/sc_datasets/soybean/featured_datasets_anndata",
        help="Path containing exported tissue datasets (h5ad or matrix_market)"
    )
    parser.add_argument(
        "--standardized_dir", type=str,
        default=f"{drive_prefix}/PhD/sc_datasets/soybean/featured_datasets_standardized",
        help="Fallback path containing standardized R RDS / CSV datasets"
    )
    parser.add_argument(
        "--output_dir", type=str,
        default=f"{drive_prefix}/PhD/data_simulation_comparison/scMultiSim/soybean_datasets/soybean_paga_results",
        help="Directory to save Newick trees, connectivity matrices, and publication figures"
    )
    parser.add_argument(
        "--tissue", type=str, default="all",
        help="Specific tissue name or 'all' to process all 7 tissues"
    )
    parser.add_argument(
        "--root_celltype", type=str, default=None,
        help="Custom root cell type name"
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


def load_soy_dataset(tissue, input_dir, standardized_dir):
    candidates_h5 = [
        os.path.join(input_dir, tissue, f"soy_processed_{tissue}.h5ad"),
        os.path.join(input_dir, f"soy_processed_{tissue}.h5ad"),
        os.path.join(standardized_dir, tissue, f"soy_processed_{tissue}.h5ad")
    ]
    for h5 in candidates_h5:
        if os.path.exists(h5):
            print(f"  [Loader] Loading AnnData from H5AD: {h5}")
            return sc.read_h5ad(h5)

    candidates_mm = [
        os.path.join(input_dir, tissue, "matrix_market"),
        os.path.join(standardized_dir, tissue, "matrix_market")
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
        f"Could not find Soybean dataset for '{tissue}' in:\n  1) {input_dir}\n  2) {standardized_dir}\n"
        "Please run 'export_soybean_standardized_to_anndata.R' first."
    )


def build_classic_newick(tree, node, conn, categories, parent=None):
    children = list(tree.successors(node))
    if parent is not None:
        idx_p = categories.index(parent)
        idx_c = categories.index(node)
        w = conn[idx_p, idx_c]
        branch_len = round(max(0.05, min(2.0, 1.0 / w)), 4) if w > 0 else 0.5
    else:
        branch_len = 0.0
    label = str(node).replace(' ', '_')
    if not children:
        return f"{label}:{branch_len}"
    child_str = ','.join(build_classic_newick(tree, c, conn, categories, node) for c in children)
    return f"({child_str}){label}:{branch_len}"


def build_all_tips_newick(tree, node, conn, categories, parent=None):
    children = list(tree.successors(node))
    if parent is not None:
        idx_p = categories.index(parent)
        idx_c = categories.index(node)
        w = conn[idx_p, idx_c]
        branch_len = round(max(0.05, min(2.0, 1.0 / w)), 4) if w > 0 else 0.5
    else:
        branch_len = 0.0
    label = str(node).replace(' ', '_')
    if not children:
        return f"{label}:{branch_len}"
    child_strings = [build_all_tips_newick(tree, c, conn, categories, node) for c in children]
    self_tip = f"{label}_state:0.05"
    child_strings.append(self_tip)
    child_str = ','.join(child_strings)
    return f"({child_str}){label}_node:{branch_len}"


def generate_publication_figures(adata, tissue, categories, celltype_colors, T, root, fig_dir, dpi=600):
    n_cats = len(categories)
    
    # 1. PAGA Graph
    fig, ax = plt.subplots(figsize=(9, 8), dpi=dpi)
    sc.pl.paga(
        adata, color='cell_type', threshold=0.03,
        node_size_scale=1.8, edge_width_scale=1.5, fontsize=9.5,
        ax=ax, show=False
    )
    ax.set_title(f"PAGA Trajectory: Soybean {tissue.replace('_', ' ')}", fontsize=13, fontweight='bold', pad=15)
    fig.savefig(os.path.join(fig_dir, f"Fig1_{tissue}_paga_graph.png"), bbox_inches='tight', dpi=dpi)
    fig.savefig(os.path.join(fig_dir, f"Fig1_{tissue}_paga_graph.pdf"), bbox_inches='tight')
    plt.close(fig)

    # 2. UMAP
    umap_coords = adata.obsm['X_umap']
    umap_df = pd.DataFrame({
        'u1': umap_coords[:, 0], 'u2': umap_coords[:, 1],
        'cell_type': adata.obs['cell_type'].values
    })
    centers = umap_df.groupby('cell_type')[['u1', 'u2']].median().reindex(categories)

    fig, (ax_main, ax_leg) = plt.subplots(
        1, 2, figsize=(14, max(6, n_cats * 0.7)), dpi=dpi,
        gridspec_kw={'width_ratios': [3.2, 1.2], 'wspace': 0.08}
    )
    for ct in categories:
        mask = umap_df['cell_type'] == ct
        ax_main.scatter(
            umap_df.loc[mask, 'u1'], umap_df.loc[mask, 'u2'],
            c=celltype_colors[ct], s=12, alpha=0.75, linewidths=0, rasterized=True
        )
    for ct, row in centers.iterrows():
        if pd.isna(row['u1']): continue
        ax_main.text(
            row['u1'], row['u2'], str(ct).replace('_', '\n'),
            fontsize=8, fontweight='bold', ha='center', va='center',
            path_effects=[pe.withStroke(linewidth=3.5, foreground='white')]
        )
    ax_main.set_xlabel('UMAP 1', fontsize=11, fontweight='bold')
    ax_main.set_ylabel('UMAP 2', fontsize=11, fontweight='bold')
    ax_main.set_title(f"Soybean {tissue.replace('_', ' ')}: PAGA-Initialized UMAP", fontsize=13, fontweight='bold')
    ax_main.tick_params(labelsize=9)

    ax_leg.axis('off')
    for i, ct in enumerate(categories):
        y = n_cats - i
        ax_leg.scatter([0.2], [y], s=120, c=[celltype_colors[ct]], edgecolors='black', linewidths=0.5)
        count = (adata.obs['cell_type'] == ct).sum()
        ax_leg.text(0.35, y, f"{ct.replace('_', ' ')} (n={count})", ha='left', va='center', fontsize=9, color='#111111')
    ax_leg.set_xlim(0, 2)
    ax_leg.set_ylim(0.5, n_cats + 0.5)
    fig.savefig(os.path.join(fig_dir, f"Fig2_{tissue}_paga_umap_celltypes.png"), bbox_inches='tight', dpi=dpi)
    fig.savefig(os.path.join(fig_dir, f"Fig2_{tissue}_paga_umap_celltypes.pdf"), bbox_inches='tight')
    plt.close(fig)

    # 3. Diffusion Pseudotime
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
        cbar.set_label('Diffusion Pseudotime (Root: ' + str(root) + ')', fontsize=10, fontweight='bold')
        
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
        ax_pt.set_title(f"Soybean {tissue}: Diffusion Pseudotime Progression", fontsize=13, fontweight='bold')
        fig.savefig(os.path.join(fig_dir, f"Fig3_{tissue}_diffusion_pseudotime.png"), bbox_inches='tight', dpi=dpi)
        fig.savefig(os.path.join(fig_dir, f"Fig3_{tissue}_diffusion_pseudotime.pdf"), bbox_inches='tight')
        plt.close(fig)

    # 4. MST Tree
    fig, ax_tr = plt.subplots(figsize=(10, 8), dpi=dpi)
    pos = nx.spring_layout(T, seed=42)
    node_colors = [celltype_colors.get(n, '#4A90E2') for n in T.nodes]
    node_sizes = [max(280, (adata.obs['cell_type'] == n).sum() * 0.6) for n in T.nodes]
    nx.draw_networkx_nodes(T, pos, node_color=node_colors, node_size=node_sizes, ax=ax_tr, edgecolors='black', linewidths=1)
    nx.draw_networkx_edges(T, pos, width=2.0, alpha=0.7, edge_color='#444444', ax=ax_tr)
    labels = {n: str(n).replace('_', '\n') for n in T.nodes}
    nx.draw_networkx_labels(T, pos, labels=labels, font_size=8, font_weight='bold', ax=ax_tr)
    ax_tr.set_title(f"Soybean {tissue} Maximum Spanning Tree (MST)\nRoot: {root}", fontsize=12, fontweight='bold')
    ax_tr.axis('off')
    fig.savefig(os.path.join(fig_dir, f"Fig4_{tissue}_lineage_tree_mst.png"), bbox_inches='tight', dpi=dpi)
    fig.savefig(os.path.join(fig_dir, f"Fig4_{tissue}_lineage_tree_mst.pdf"), bbox_inches='tight')
    plt.close(fig)


def process_single_tissue(tissue, opt):
    print("\n" + "=" * 70)
    print(f"PROCESSING SOYBEAN TISSUE: {tissue.upper()}")
    print("=" * 70)

    tissue_out_dir = os.path.join(opt.output_dir, tissue)
    fig_dir = os.path.join(tissue_out_dir, "figures")
    os.makedirs(tissue_out_dir, exist_ok=True)
    os.makedirs(fig_dir, exist_ok=True)

    adata = load_soy_dataset(tissue, opt.input_dir, opt.standardized_dir)
    print(f"  [Dimensions] {adata.n_obs} cells x {adata.n_vars} genes")

    ct_col = 'cell_type'
    if ct_col not in adata.obs.columns:
        for c in ['cluster', 'seurat_clusters', 'CellType']:
            if c in adata.obs.columns:
                adata.obs[ct_col] = adata.obs[c]
                break
    adata.obs[ct_col] = adata.obs[ct_col].astype('category')
    categories = adata.obs[ct_col].cat.categories.tolist()
    n_cats = len(categories)
    print(f"  [Cell Types] Found {n_cats} categories: {categories}")

    x_max = adata.X.max() if sp.issparse(adata.X) else np.max(adata.X)
    if x_max > 25:
        print("  [Preprocess] Normalizing total counts + log1p...")
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)

    n_hvg = min(2000, adata.n_vars)
    if adata.n_vars > n_hvg:
        sc.pp.highly_variable_genes(adata, n_top_genes=n_hvg)
    sc.pp.pca(adata, n_comps=min(50, adata.n_obs - 1, adata.n_vars - 1), use_highly_variable=True if 'highly_variable' in adata.var else False)

    n_pcs_use = min(opt.n_pcs, adata.obsm['X_pca'].shape[1])
    sc.pp.neighbors(adata, n_neighbors=opt.n_neighbors, n_pcs=n_pcs_use)
    sc.tl.paga(adata, groups=ct_col)
    sc.pl.paga(adata, show=False)
    sc.tl.umap(adata, init_pos='paga')

    # Root selection
    root = None
    if opt.root_celltype and opt.root_celltype in categories:
        root = opt.root_celltype
    else:
        for cand in SOYBEAN_PROGENITORS.get(tissue, []):
            if cand in categories:
                root = cand
                break
    if root is None:
        conn_mat = adata.uns['paga']['connectivities'].toarray()
        root = categories[np.argmax(conn_mat.sum(axis=1))]
    print(f"  [Root] Selected progenitor root: '{root}'")

    # DPT
    try:
        root_idx = np.flatnonzero(adata.obs[ct_col] == root)
        if len(root_idx) > 0:
            adata.uns['iroot'] = int(root_idx[0])
            sc.tl.diffmap(adata)
            sc.tl.dpt(adata)
    except Exception as e:
        print(f"  [DPT Notice] {e}")

    conn = adata.uns['paga']['connectivities'].toarray()
    paga_conn_df = pd.DataFrame(conn, index=categories, columns=categories)
    conn_csv_path = os.path.join(tissue_out_dir, f"{tissue}_paga_connectivity.csv")
    paga_conn_df.to_csv(conn_csv_path)

    G = nx.Graph()
    for i, ct_i in enumerate(categories):
        for j, ct_j in enumerate(categories):
            if i < j and conn[i, j] > 0:
                G.add_edge(ct_i, ct_j, weight=conn[i, j])
    T = nx.maximum_spanning_tree(G, weight='weight')
    T_rooted = nx.bfs_tree(T, root)

    nwk_all_tips = build_all_tips_newick(T_rooted, root, conn, categories) + ";"
    nwk_classic = build_classic_newick(T_rooted, root, conn, categories) + ";"

    with open(os.path.join(tissue_out_dir, f"{tissue}_paga_tree_all_tips.nwk"), 'w') as f:
        f.write(nwk_all_tips)
    with open(os.path.join(tissue_out_dir, f"{tissue}_paga_tree_classic.nwk"), 'w') as f:
        f.write(nwk_classic)
    with open(os.path.join(tissue_out_dir, f"{tissue}_paga_tree.nwk"), 'w') as f:
        f.write(nwk_all_tips)

    ct_counts = adata.obs[ct_col].value_counts().reindex(categories)
    med_pt = adata.obs.groupby(ct_col, observed=True)['dpt_pseudotime'].median().reindex(categories) if 'dpt_pseudotime' in adata.obs else pd.Series(0.0, index=categories)
    ct_summary = pd.DataFrame({
        'cell_type': categories,
        'n_cells': ct_counts.values,
        'is_progenitor_root': [ct == root for ct in categories],
        'median_pseudotime': med_pt.values
    })
    ct_summary.to_csv(os.path.join(tissue_out_dir, f"{tissue}_celltype_metadata.csv"), index=False)

    pal = plt.cm.get_cmap('tab20', n_cats) if n_cats <= 20 else plt.cm.get_cmap('gist_ncar', n_cats)
    celltype_colors = {ct: mcolors.to_hex(pal(i)) for i, ct in enumerate(categories)}
    adata.uns['cell_type_colors'] = [celltype_colors[ct] for ct in categories]

    print("  [Graphics] Generating 600 DPI publication figures...")
    generate_publication_figures(adata, tissue, categories, celltype_colors, T, root, fig_dir, dpi=opt.dpi)

    h5_path = os.path.join(tissue_out_dir, f"{tissue}_paga_processed.h5ad")
    adata.write_h5ad(h5_path)
    print(f"  [Export] Saved AnnData: {os.path.basename(h5_path)}")

    return {
        'tissue': tissue,
        'n_cells': adata.n_obs,
        'n_genes': adata.n_vars,
        'n_cell_types': n_cats,
        'root': root,
        'tree_path': os.path.join(tissue_out_dir, f"{tissue}_paga_tree.nwk")
    }


def main():
    opt = parse_arguments()
    np.random.seed(opt.seed)
    
    print("=" * 75)
    print("SOYBEAN PAGA LINEAGE TREE RECONSTRUCTION SUITE")
    print("Targeting Simulator: scMultiSim (Nature Methods 2023)")
    print("=" * 75)

    tissues_to_run = ALL_SOY_TISSUES if opt.tissue.lower() == "all" else [t.strip() for t in opt.tissue.split(",")]
    records = []
    for tissue in tissues_to_run:
        try:
            res = process_single_tissue(tissue, opt)
            records.append(res)
        except Exception as e:
            print(f"[ERROR] Failed processing tissue '{tissue}': {e}")
            import traceback
            traceback.print_exc()

    if records:
        df = pd.DataFrame(records)
        df.to_csv(os.path.join(opt.output_dir, "soybean_all_tissues_paga_summary.csv"), index=False)
        print("\n" + "=" * 75)
        print("ALL SOYBEAN TISSUES PROCESSED SUCCESSFULLY!")
        print("=" * 75)
        print(df.to_string(index=False))


if __name__ == "__main__":
    main()
