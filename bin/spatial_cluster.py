#!/usr/bin/env python3
"""Cluster filtered Visium spots and visualize their spatial organization."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import scanpy as sc
import squidpy as sq


def parse_args() -> argparse.Namespace:
    """Parse command-line options for spatial clustering."""
    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--figures-dir", type=Path, default=None)
    parser.add_argument("--n-hvg", type=int, default=2000)
    parser.add_argument("--n-pcs", type=int, default=30)
    parser.add_argument("--n-neighbors", type=int, default=15)
    parser.add_argument("--resolution", type=float, default=1.0)

    return parser.parse_args()


def preprocess_spots(adata, n_hvg: int, n_pcs: int):
    """Normalize spots, select highly variable genes, and compute PCA."""
    sc.pp.normalize_total(adata, target_sum=10_000)
    sc.pp.log1p(adata)

    #Preserve the normalized full-gene matrix
    adata.raw = adata

    sc.pp.highly_variable_genes(
    adata,
    n_top_genes=n_hvg,
    flavor="seurat",
    )
    sc.tl.pca(
        adata,
        n_comps=n_pcs,
        mask_var="highly_variable",
        svd_solver="arpack",
        random_state=0,
    )
    return adata


def cluster_spots(adata, n_neighbors: int, resolution: float):
    """Build an expression-based neighbor graph and assign Leiden clusters."""
    n_pcs = adata.obsm["X_pca"].shape[1]

    sc.pp.neighbors(
        adata,
        n_neighbors=n_neighbors,
        n_pcs=n_pcs,
        key_added="expression_neighbors",
        random_state=0,
    )

    sc.tl.leiden(
        adata,
        resolution=resolution,
        neighbors_key="expression_neighbors",
        key_added="spatial_leiden",
        flavor="leidenalg",
        random_state=0,
    )

    sc.tl.umap(
        adata,
        neighbors_key="expression_neighbors",
        random_state=0,
    )

    #Seperate graph based on physical spot locations
    sq.gr.spatial_neighbors(adata)

    print("Number of spatial clusters:", adata.obs["spatial_leiden"].nunique())

    return adata


def make_spatial_plot(adata, output_dir: Path) -> None:
    """Save the spot clusters over the Visium tissue image."""
    
    output_dir.mkdir(parents=True, exist_ok=True)

    sq.pl.spatial_scatter(
        adata,
        color="spatial_leiden",
        img_res_key="hires",
        return_ax=True,
        legend_loc="right margin",
    )

    plt.savefig(output_dir / "spatial_clusters.png", dpi=150, bbox_inches="tight")
    plt.close()

    sc.pl.umap(adata, color="spatial_leiden", show=False)
    plt.savefig(output_dir / "spatial_clusters_umap.png", dpi=150, bbox_inches="tight")
    plt.close()


def main() -> int:
    """Run the complete Visium clustering workflow."""
    
    args = parse_args()

    if not args.input.exists():
        raise FileNotFoundError(
            f"Input file does not exist: {args.input}"
        )

    args.output.parent.mkdir(parents=True,exist_ok=True)

    figures_dir = (
        args.figures_dir 
        or args.output.parent / "figures"
        )

    adata = sc.read_h5ad(args.input)

    if "spatial" not in adata.obsm:
        raise ValueError(
            "Input is missing spatial coordinates in adata.obsm['spatial']"
        )

    adata = preprocess_spots(adata, args.n_hvg, args.n_pcs)
    adata = cluster_spots(adata, args.n_neighbors, args.resolution)
    
    make_spatial_plot(adata, figures_dir)
    adata.write_h5ad(args.output)

    print(f"Saved clustered data: {args.output}")
    print(f"Saved figures: {figures_dir}")
    print(
        f"Final data: "
        f"{adata.n_obs:,} spots x {adata.n_vars:,} genes"
    )
    return 0

if __name__ == "__main__":
    raise SystemExit(main())