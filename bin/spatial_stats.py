#!/usr/bin/env python3
"""Compute neighbourhood enrichment and Moran's I on deconvolved Visium spots.

1. load and validate the deconvolved AnnData object;
2. neighbourhood enrichment on the spatial graph;
3. Moran's I for spatially variable genes;
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import squidpy as sq


def parse_args():
    """Parse command-line options for the spatial statistics workflow."""
    parser = argparse.ArgumentParser(
        description="Run neighbourhood enrichment and Moran's I on deconvolved Visium spots."
    )
    parser.add_argument("--spatial-input", type=Path, required=True, help="Deconvolved spatial AnnData file.")
    parser.add_argument("--outdir", type=Path, required=True, help="Directory for spatial stats outputs.")
    parser.add_argument(
        "--cluster-key",
        default="spatial_leiden",
        help="obs column with categorical spot labels for neighbourhood enrichment.",
    )
    parser.add_argument(
        "--nhood-n-perms", type=int, default=1000, help="Permutations for neighbourhood enrichment."
    )
    parser.add_argument("--moran-n-perms", type=int, default=100, help="Permutations for Moran's I.")
    parser.add_argument("--n-top-svgs", type=int, default=4, help="Number of top Moran genes to plot.")
    return parser.parse_args()


def load_inputs(spatial_input: Path):
    """Load the deconvolved spatial AnnData object."""
    if not spatial_input.is_file():
        raise FileNotFoundError(f"Spatial input not found: {spatial_input}")

    return sc.read_h5ad(spatial_input)


def validate_inputs(adata, cluster_key):
    """Check labels, coordinates, and the spatial neighbour graph."""
    if adata.n_obs == 0 or adata.n_vars == 0:
        raise ValueError("The spatial AnnData object is empty.")

    if cluster_key not in adata.obs.columns:
        raise KeyError(
            f"Cluster column '{cluster_key}' not found. Available columns: {list(adata.obs.columns)}"
        )

    labels = adata.obs[cluster_key]
    if labels.isna().all():
        raise ValueError(f"All values in '{cluster_key}' are missing.")

    n_clusters = labels.nunique(dropna=True)
    if n_clusters < 2:
        raise ValueError(
            f"Neighbourhood enrichment needs at least two labels; '{cluster_key}' has {n_clusters}."
        )

    if "spatial" not in adata.obsm:
        raise KeyError("Spatial coordinates not found in adata.obsm['spatial'].")

    coordinates = np.asarray(adata.obsm["spatial"])
    if coordinates.shape[0] != adata.n_obs:
        raise ValueError(
            "The number of spatial coordinate rows does not match the number of spatial observations."
        )
    if not np.isfinite(coordinates).all():
        raise ValueError("Spatial coordinates contain NaN or infinite values.")

    if "spatial_connectivities" not in adata.obsp:
        raise KeyError(
            "Spatial neighbour graph not found in adata.obsp['spatial_connectivities']. "
            "Do not use expression_neighbors; rebuild spatial neighbours if needed."
        )


def run_nhood_enrichment(adata, cluster_key, n_perms, outdir):
    """Permute spatial-graph labels and write z-scores plus a heatmap."""
    outdir = Path(outdir)
    figdir = outdir / "figures"
    figdir.mkdir(parents=True, exist_ok=True)

    adata.obs[cluster_key] = adata.obs[cluster_key].astype("category")

    sq.gr.nhood_enrichment(
        adata,
        cluster_key=cluster_key,
        connectivity_key="spatial_connectivities",
        n_perms=n_perms,
        show_progress_bar=False,
    )

    result_key = f"{cluster_key}_nhood_enrichment"
    result = adata.uns[result_key]
    zscores = np.asarray(result["zscore"])
    labels = list(adata.obs[cluster_key].cat.categories)
    zscore_df = pd.DataFrame(zscores, index=labels, columns=labels)

    tsv_path = outdir / "nhood_enrichment.tsv"
    zscore_df.to_csv(tsv_path, sep="\t")

    sq.pl.nhood_enrichment(adata, cluster_key=cluster_key, return_ax=True)
    plt.savefig(figdir / "nhood_enrichment.png", dpi=150, bbox_inches="tight")
    plt.close()

    print("Neighbourhood enrichment key:", result_key)
    print("Z-score shape:", zscore_df.shape)
    print("Wrote", tsv_path)
    print("Wrote", figdir / "nhood_enrichment.png")
    return zscore_df


def get_moran_genes(adata):
    """Return highly variable gene names for Moran's I."""
    if "highly_variable" not in adata.var.columns:
        raise KeyError("adata.var['highly_variable'] is missing; cannot subset genes for Moran.")

    genes = adata.var_names[adata.var["highly_variable"].to_numpy()].tolist()
    if len(genes) == 0:
        raise ValueError("No highly variable genes are marked True; cannot run Moran.")

    print("Moran gene count (HVGs):", len(genes))
    return genes


def run_moran(adata, genes, n_perms, outdir):
    """Compute Moran's I on the spatial graph and write a ranked table."""
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    sq.gr.spatial_autocorr(
        adata, connectivity_key="spatial_connectivities", genes=genes, mode="moran", n_perms=n_perms
    )

    if "moranI" not in adata.uns:
        raise KeyError(
            f"Squidpy did not write adata.uns['moranI']. Available uns keys: {list(adata.uns.keys())}"
        )

    moran_table = pd.DataFrame(adata.uns["moranI"]).copy()
    if "I" not in moran_table.columns:
        raise KeyError(f"moranI table has no 'I' column. Columns: {list(moran_table.columns)}")

    moran_table = moran_table.sort_values("I", ascending=False)
    tsv_path = outdir / "moran.tsv"
    moran_table.to_csv(tsv_path, sep="\t")

    top_gene = moran_table.index[0]
    print("Moran table shape:", moran_table.shape)
    print("Top gene:", top_gene, "I =", float(moran_table.iloc[0]["I"]))
    print("Wrote", tsv_path)
    return moran_table


def plot_top_svgs(adata, moran_table, outdir, n=4):
    """Plot the top spatially autocorrelated genes on the H&E image."""
    figdir = Path(outdir) / "figures"
    figdir.mkdir(parents=True, exist_ok=True)
    top_genes = [str(gene) for gene in moran_table.index[:n]]
    written = []

    for gene in top_genes:
        if gene not in adata.var_names:
            raise KeyError(f"Top Moran gene '{gene}' is not in adata.var_names.")

        sq.pl.spatial_scatter(adata, color=gene, img_res_key="hires", return_ax=True)
        safe_name = gene.replace("/", "_").replace("+", "pos")
        path = figdir / f"moran_{safe_name}.png"
        plt.savefig(path, dpi=150, bbox_inches="tight")
        plt.close()
        written.append(path)
        print("Wrote", path)

    return written


def main():
    """Load, validate, run neighbourhood enrichment, then Moran's I."""
    args = parse_args()
    adata = load_inputs(args.spatial_input)
    validate_inputs(adata, args.cluster_key)

    print("Shape:", adata.shape)
    print("Cluster key:", args.cluster_key)
    print("Number of clusters:", adata.obs[args.cluster_key].nunique(dropna=True))
    print("Spatial graph present:", "spatial_connectivities" in adata.obsp)
    print("Neighbourhood permutations:", args.nhood_n_perms)
    print("Moran's I permutations:", args.moran_n_perms)
    print(adata.obs[args.cluster_key].astype(str).value_counts())

    args.outdir.mkdir(parents=True, exist_ok=True)
    run_nhood_enrichment(adata, args.cluster_key, args.nhood_n_perms, args.outdir)

    genes = get_moran_genes(adata)
    moran_table = run_moran(adata, genes, args.moran_n_perms, args.outdir)
    plot_top_svgs(adata, moran_table, args.outdir, n=args.n_top_svgs)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
