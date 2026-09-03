#!/usr/bin/env python3
"""Perform quality control and filtering for Visium spatial transcriptomics data.

Loads a Visium sample, calculates spot-level QC metrics, removes background and
low-quality spots, generates diagnostic figures, and saves a filtered AnnData file.

Unlike scRNA-seq QC, Visium observations are spatial spots containing mixtures of
multiple cells rather than individual cells.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import scanpy as sc
import squidpy as sq


def parse_args() -> argparse.Namespace:
    """Parse command-line options for Visium spot-level quality control."""
    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument(
        "--input",
        type=Path,
        required=True,
        help="Visium sample directory containing matrix and spatial files",
    )
    parser.add_argument(
        "--output",
        type=Path,
        required=True,
        help="Output filtered Visium .h5ad file",
    )
    parser.add_argument(
        "--figures-dir",
        type=Path,
        default=None,
        help="Directory for QC figures; defaults to <output parent>/figures",
    )
    parser.add_argument(
        "--min-counts",
        type=int,
        default=500,
        help="Minimum total counts required per spot",
    )
    parser.add_argument(
        "--min-genes",
        type=int,
        default=200,
        help="Minimum detected genes required per spot",
    )
    parser.add_argument(
        "--max-pct-mt",
        type=float,
        default=20.0,
        help="Maximum mitochondrial percentage, if mitochondrial genes exist",
    )

    return parser.parse_args()

def load_visium(input_dir: Path):
    """Load a Visium sample directory into an AnnData object."""
    if not input_dir.exists():
        raise FileNotFoundError(f"Visium directory does not exist: {input_dir}")

    if not input_dir.is_dir():
        raise NotADirectoryError(f"Expected a directory: {input_dir}")

    adata = sq.read.visium(input_dir)
    adata.var_names_make_unique()

    if "spatial" not in adata.obsm:
        raise ValueError("Visium data is missing spatial coordinates in adata.obsm['spatial']")

    print(f"Loaded Visium data: {adata.n_obs:,} spots x {adata.n_vars:,} genes")
    return adata

def calculate_qc(adata):
    """Calculate spot-level count, gene, tissue, and mitochondrial QC metrics."""

    adata.var["mt"] = (
        adata.var_names.str.upper().str.startswith(("MT-","MT."))
    )

    has_mito = bool(adata.var["mt"].any())

    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars = ["mt"] if has_mito else [],
        inplace=True,
        log1p=False,
    )

    if not has_mito:
        adata.obs["pct_counts_mt"] = 0.0
        print("No mitochondrial genes found; mitochondrial filtering will be skipped.")
    else:
        print("Mitochondrial QC metrics calculated.")

    print(
    f"QC metrics calculated for "
    f"{adata.n_obs:,} spots and {adata.n_vars:,} genes."
    )

    return adata


def filter_spots(adata, min_counts: int, min_genes: int, max_pct_mt: float):
    """Remove low-quality or off-tissue Visium spots."""

    initial_spots = adata.n_obs

    tissue_mask = adata.obs["in_tissue"].astype(str).isin(
        ["1","True","true"]
    )

    adata = adata[tissue_mask].copy()
    spots_after_tissue = adata.n_obs
    
    sc.pp.filter_cells(adata, min_counts=min_counts)

    sc.pp.filter_cells(adata, min_genes=min_genes)

    if adata.var["mt"].any():
        adata = adata[
            adata.obs["pct_counts_mt"] <= max_pct_mt
        ].copy()

    print(f"Initial spots: {initial_spots:,}")
    print(f"On-tissue spots: {spots_after_tissue:,}")
    print(f"Final spots: {adata.n_obs:,}")
    print(f"Spots removed: {initial_spots - adata.n_obs:,}")
    
    return adata

def save_qc_figures(adata, output_dir: Path) -> None:
    """Save diagnostic plots showing Visium spot quality."""

    output_dir.mkdir(parents=True, exist_ok=True)

    sc.pl.violin(
        adata,
        ["total_counts","n_genes_by_counts","pct_counts_mt"],
        jitter=0.2,
        multi_panel=True,
        show=False,
    )
    plt.savefig(output_dir / "spatial_qc_violin.png", dpi=150, bbox_inches="tight")
    plt.close()

    sq.pl.spatial_scatter(
        adata,
        color="total_counts",
        img_res_key="hires",
        return_ax=True,
    )
    plt.savefig(output_dir / "spatial_qc_counts.png", dpi=150, bbox_inches="tight")
    plt.close()

def main() -> int:
    """Run the complete Visium QC workflow."""

    args = parse_args()

    if not args.input.exists():
        raise FileNotFoundError(
            f"Visium directory does not exist: {args.input}"
        )
    
    args.output.parent.mkdir(parents=True, exist_ok=True)

    figures_dir = (
        args.figures_dir
        or args.output.parent / "figures"
    )

    adata = load_visium(args.input)
    adata = calculate_qc(adata)

    adata = filter_spots(adata, args.min_counts, args.min_genes, args.max_pct_mt)

    save_qc_figures(
        adata,
        figures_dir,
    )

    adata.write_h5ad(args.output)
    
    print(f"Saved filtered data: {args.output}")
    print(f"Saved figures: {figures_dir}")
    print(
        f"Final data: "
        f"{adata.n_obs:,} spots x {adata.n_vars:,} genes"
    )
    return 0

if __name__ == "__main__":
    raise SystemExit(main())