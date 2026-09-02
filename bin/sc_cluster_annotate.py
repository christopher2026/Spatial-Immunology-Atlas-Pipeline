#!/usr/bin/env python3
"""Analyze scRNA-seq QC output through clustering and annotation validation.

Normalization -> highly variable genes -> PCA -> neighbor graph -> Leiden clustering
-> UMAP visualization -> comparison with curated cell-type labels.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
import pandas as pd
import scanpy as sc

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def parse_args() -> argparse.Namespace:
    """Define command-line options for the clustering and validation workflow."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True, help="Filtered scRNA-seq .h5ad file")
    parser.add_argument("--output", type=Path, required=True, help="Annotated output .h5ad file")
    parser.add_argument(
        "--figures-dir",
        type=Path,
        default=None,
        help="Output directory for figures (default: <output parent>/figures)",
    )
    parser.add_argument(
        "--agreement",
        type=Path,
        default=None,
        help="Cluster-versus-curated-label TSV (default: <output stem>.agreement.tsv)",
    )
    parser.add_argument(
        "--celltype-key",
        default="Subset",
        help="Existing curated cell-type column in .obs (default: Subset)",
    )
    parser.add_argument("--n-hvg", type=int, default=2000)
    parser.add_argument("--n-pcs", type=int, default=50)
    parser.add_argument("--n-neighbors", type=int, default=15)
    parser.add_argument("--resolution", type=float, default=1.0)
    return parser.parse_args()


def preprocess_data(adata, n_hvg: int, n_pcs: int):
    """Normalize expression, select informative genes, and compute PCA."""
    sc.pp.normalize_total(adata,target_sum = 10_000)
    sc.pp.log1p(adata)
    
    # Preserve the full normalized dataset for later marker analysis.
    adata.raw = adata

    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=n_hvg,
        flavor="seurat",
    )

    sc.tl.pca(
        adata,
        n_comps=n_pcs,
        use_highly_variable=True,
        svd_solver="arpack",
        random_state=0
    )
    return adata

def cluster_cells(adata, n_neighbors: int, resolution: float):
    """Build a PCA-space neighbor graph and assign Leiden clusters."""
    n_pcs = adata.obsm["X_pca"].shape[1]

    sc.pp.neighbors(
        adata,
        n_neighbors=n_neighbors,
        n_pcs=n_pcs,
        random_state=0,
    )

    sc.tl.leiden(
        adata,
        resolution=resolution,
        key_added="leiden",
        flavor="leidenalg",
        random_state=0,
    )

    sc.tl.umap(
        adata,
        random_state=0,
    )
    
    print("connectivities:", adata.obsp["connectivities"].shape)
    print("number of graph edges:", adata.obsp["connectivities"].nnz)

    print("\nLeiden clusters:")
    print(adata.obs["leiden"].value_counts().sort_index())

    print("\nUMAP shape:", adata.obsm["X_umap"].shape)
    print("First 5 UMAP coordinates:")
    print(adata.obsm["X_umap"][:5])   

    return adata


def make_plots(adata, output_dir: Path, celltype_key: str) -> None:
    """Save UMAP and marker visualizations."""
    output_dir.mkdir(parents=True, exist_ok=True)    
    sc.pl.umap(adata, color="leiden",show=False)
    plt.savefig(output_dir / "sc_cluster_leiden_umap.png", dpi=150, bbox_inches="tight")
    plt.close()

    sc.pl.umap(adata, color=celltype_key, show=False)    
    plt.savefig(output_dir / "sc_cluster_subset_umap.png", dpi=150, bbox_inches="tight")
    plt.close()

    sc.tl.rank_genes_groups(
        adata,
        groupby="leiden",
    method="wilcoxon",
    )

    markers = {
    "B cells": ["MS4A1", "CD79A", "CD74", "HLA-DRA", "CD37"],
    "T cells": ["CD3D", "CD3E", "TRBC1", "TRBC2", "LTB"],
    "CD4 T / Tfh": ["IL7R", "CXCR5", "PDCD1", "ICOS", "BCL6"],
    "CD8 / cytotoxic": ["CD8A", "CD8B", "NKG7", "GNLY", "PRF1"],
    "NK cells": ["NKG7", "GNLY", "KLRD1", "PRF1", "GZMB"],
    "Plasma cells": ["MZB1", "JCHAIN", "XBP1", "SEC11C", "IGKC"],
    "Myeloid": ["LST1", "TYROBP", "FCER1G", "CTSS", "LYZ"],
    "Macrophages": ["C1QC", "APOC1", "CTSD", "CTSB", "LILRB1"],
    "Dendritic cells": ["FCER1A", "CLEC10A", "CST3", "GZMB", "IRF7"],
    "Mast cells": ["KIT", "TPSAB1", "TPSB2", "MS4A2"],
    "Endothelial": ["PECAM1", "VWF", "EMCN", "KDR"],
    "Stromal": ["COL1A1", "COL1A2", "DCN", "COL3A1", "PDGFRA"],
    "FDC / follicular": ["CXCL13", "CR1", "CR2", "FDCSP"],
    }
    
    #check if markers exist in reference 
    available_genes = set(adata.raw.var_names)

    for cell_type, genes in markers.items():
        missing = [gene for gene in genes if gene not in available_genes]

        if missing:
            print(f"{cell_type} markers not found: {missing}")

        markers[cell_type] = [
            gene for gene in genes if gene in available_genes
    ]
    
    sc.pl.dotplot(
    adata,
    var_names=markers,
    groupby="leiden",
    use_raw=True,
    show=False,   
    )
    plt.savefig(output_dir / "sc_cluster_gene_marker_dotplot.png", dpi=150, bbox_inches="tight")
    plt.close()


def validate_clusters(adata, celltype_key: str, output_path: Path) -> None:
    """Compare independently computed Leiden clusters with curated labels."""
    agreement = pd.crosstab(
        adata.obs["leiden"],
        adata.obs[celltype_key],
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    agreement.to_csv(output_path, sep="\t")


def main() -> int:
    """Run the complete clustering, visualization, and validation workflow."""
    args = parse_args()
    if not args.input.exists():
        raise FileNotFoundError(f"Input does not exist: {args.input}")
    
    args.output.parent.mkdir(parents=True, exist_ok=True)  

    figures_dir = (
        args.figures_dir
        or args.output.parent / "figures"
    )

    agreement_path = (
        args.agreement
        or args.output.with_suffix(".agreement.tsv")
    )

    adata = sc.read_h5ad(args.input)

    if args.celltype_key not in adata.obs.columns:
        raise KeyError(
            f"Missing cell-type column: {args.celltype_key}"
        )
    adata = preprocess_data(
    adata,
    args.n_hvg,
    args.n_pcs,
    )

    adata = cluster_cells(
        adata,
        args.n_neighbors,
        args.resolution,
    )

    make_plots(
        adata,
        figures_dir,
        args.celltype_key,
    )

    validate_clusters(
        adata,
        args.celltype_key,
        agreement_path,
    )
    
    adata.write_h5ad(args.output)

    print(f"Saved annotated data: {args.output}")
    print(f"Saved figures: {figures_dir}")
    print(f"Saved agreement table: {agreement_path}")

    return 0

if __name__ == "__main__":
    raise SystemExit(main())