#!/usr/bin/env python3
"""Quality-control and filtering for the single-cell reference.

The reference already contains published QC columns, but this script deliberately recomputes QC
from the expression matrix. That makes the result reproducible and prevents us from treating prior
analysis as ground truth.

Example:
    python bin/sc_qc.py --input data/reference/sc.h5ad --output results/sc.filtered.h5ad
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import scanpy as sc


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True, help="Raw/reference AnnData .h5ad file")
    parser.add_argument("--output", type=Path, required=True, help="Output filtered .h5ad file")
    parser.add_argument(
        "--figures-dir",
        type=Path,
        default=None,
        help="QC figure directory (default: <output parent>/figures)",
    )
    parser.add_argument(
        "--summary",
        type=Path,
        default=None,
        help="QC summary TSV (default: <output stem>.qc_summary.tsv)",
    )
    parser.add_argument("--min-genes", type=int, default=200)
    parser.add_argument("--min-cells", type=int, default=3)
    parser.add_argument("--max-pct-mt", type=float, default=20.0)
    parser.add_argument(
        "--run-scrublet",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Run Scrublet (default: enabled; use --no-run-scrublet for tiny test data)",
    )
    return parser.parse_args()


def calculate_qc(adata) -> bool:
    """Add QC metrics and return whether mitochondrial genes were found."""
    adata.var["mt"] = adata.var_names.str.upper().str.startswith(("MT-", "MT."))
    has_mito = bool(adata.var["mt"].any())

    # calculate_qc_metrics accepts an empty mitochondrial set in current scanpy, but making the
    # no-mito branch explicit avoids version-dependent behaviour and documents the scientific choice.
    # This Scanpy release expects qc_vars to be iterable; None causes
    # `TypeError: 'NoneType' object is not iterable` when no MT genes exist.
    qc_vars = ["mt"] if has_mito else []
    sc.pp.calculate_qc_metrics(adata, qc_vars=qc_vars, inplace=True, log1p=False)
    if not has_mito:
        adata.obs["pct_counts_mt"] = 0.0
    return has_mito


def run_scrublet(adata) -> None:
    """Run Scrublet, retaining an explicit result even if a dataset is too small to model."""
    print(f"Running Scrublet on {adata.n_obs:,} cells; this may take several minutes...", flush=True)
    started = time.monotonic()
    try:
        sc.pp.scrublet(adata, random_state=0)
    except Exception as exc:  # third-party failures vary across scanpy/scrublet versions
        print(
            f"WARNING: Scrublet could not run ({type(exc).__name__}: {exc}). "
            "Continuing with an explicit unavailable result.",
            file=sys.stderr,
        )
        adata.obs["doublet_score"] = float("nan")
        adata.obs["predicted_doublet"] = False
        adata.uns["scrublet_status"] = f"unavailable: {type(exc).__name__}: {exc}"
    else:
        adata.uns["scrublet_status"] = "completed"
        print(f"Scrublet finished in {(time.monotonic() - started) / 60:.1f} minutes.", flush=True)


def save_qc_figures(adata, figures_dir: Path) -> None:
    figures_dir.mkdir(parents=True, exist_ok=True)

    # A violin plot answers: did filtering remove low-quality cells, and are there suspiciously
    # broad distributions? The mitochondrial panel remains informative as an explicit all-zero
    # panel when the reference has no MT genes.
    violin_keys = ["n_genes_by_counts", "total_counts", "pct_counts_mt"]
    if "predicted_doublet" in adata.obs:
        violin_keys.append("doublet_score")
    sc.pl.violin(adata, violin_keys, jitter=0.2, multi_panel=True, show=False)
    plt.savefig(figures_dir / "sc_qc_violin.png", dpi=150, bbox_inches="tight")
    plt.close()

    # Counts vs detected genes is the basic library-complexity diagnostic. A healthy sample usually
    # has a positive relationship; isolated high-count/low-gene cells are common QC outliers.
    sc.pl.scatter(
        adata,
        x="total_counts",
        y="n_genes_by_counts",
        color="pct_counts_mt",
        show=False,
    )
    plt.savefig(figures_dir / "sc_qc_counts_vs_genes.png", dpi=150, bbox_inches="tight")
    plt.close()


def main() -> int:
    args = parse_args()
    if not args.input.exists():
        raise FileNotFoundError(f"Input does not exist: {args.input}")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    figures_dir = args.figures_dir or args.output.parent / "figures"
    summary_path = args.summary or args.output.with_suffix(".qc_summary.tsv")

    print(f"Reading {args.input}")
    adata = sc.read_h5ad(args.input)
    adata.var_names_make_unique()
    input_cells, input_genes = adata.n_obs, adata.n_vars

    has_mito = calculate_qc(adata)
    print(f"Input: {input_cells:,} cells x {input_genes:,} genes")
    print(f"Mitochondrial genes found: {int(has_mito)}")

    sc.pp.filter_cells(adata, min_genes=args.min_genes)
    cells_after_min_genes = adata.n_obs
    sc.pp.filter_genes(adata, min_cells=args.min_cells)
    genes_after_min_cells = adata.n_vars

    if has_mito:
        adata = adata[adata.obs["pct_counts_mt"] <= args.max_pct_mt].copy()
    else:
        print("No mitochondrial genes present; mitochondrial filtering skipped.")

    if args.run_scrublet:
        run_scrublet(adata)
    else:
        adata.obs["doublet_score"] = float("nan")
        adata.obs["predicted_doublet"] = False
        adata.uns["scrublet_status"] = "disabled"

    save_qc_figures(adata, figures_dir)
    print("Saving filtered dataset and QC figures...", flush=True)
    adata.write_h5ad(args.output)

    summary = [
        ("input_cells", input_cells),
        ("cells_after_min_genes", cells_after_min_genes),
        ("output_cells", adata.n_obs),
        ("cells_removed", input_cells - adata.n_obs),
        ("input_genes", input_genes),
        ("genes_after_min_cells", genes_after_min_cells),
        ("output_genes", adata.n_vars),
        ("genes_removed", input_genes - adata.n_vars),
        ("mitochondrial_genes_found", int(has_mito)),
        ("max_pct_mt_applied", args.max_pct_mt if has_mito else "skipped"),
        ("scrublet_status", adata.uns["scrublet_status"]),
    ]
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    summary_path.write_text(
        "metric\tvalue\n" + "\n".join(f"{metric}\t{value}" for metric, value in summary) + "\n",
        encoding="utf-8",
    )

    print(f"Output: {adata.n_obs:,} cells x {adata.n_vars:,} genes")
    print(f"Removed: {input_cells - adata.n_obs:,} cells, {input_genes - adata.n_vars:,} genes")
    print(f"Figures: {figures_dir}")
    print(f"Summary: {summary_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
