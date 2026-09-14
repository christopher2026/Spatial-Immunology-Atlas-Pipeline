"""Map an annotated scRNA-seq reference onto Visium spatial spots with Tangram.

This command-line script is being built in stages:
1. validate the two AnnData inputs;
2. audit and prepare shared genes;
3. run Tangram's cell-to-spot optimization;
4. aggregate mappings into cell-type proportions and spatial figures.

Keeping these steps in a standalone CLI makes the analysis independently
testable before it is wrapped in a Nextflow process.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
import numpy as np
import scanpy as sc
import tangram as tg

matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Exact Subset labels from the reference; used for the first anatomical figures.
KEY_CELL_TYPES = (
    "B_GC_DZ",
    "B_GC_LZ",
    "FDC",
    "T_CD4+",
    "T_CD8+_cytotoxic",
)


def parse_args():
    """Parse command-line options for the Tangram mapping run."""
    parser = argparse.ArgumentParser(description="Map scRNA-seq cell types onto Visium spots with Tangram.")

    parser.add_argument("--spatial-input", type=Path, required=True, help="Clustered spatial AnnData file.")
    parser.add_argument("--sc-input", type=Path, required=True, help="Annotated scRNA-seq AnnData file.")
    parser.add_argument("--outdir", type=Path, required=True, help="Directory for Tangram outputs.")
    parser.add_argument(
        "--celltype-key", default="Subset", help="Reference obs column containing cell-type labels."
    )
    parser.add_argument(
        "--min-shared-genes", type=int, default=1000, help="Fail if fewer shared genes remain."
    )
    parser.add_argument(
        "--exclude-prefix",
        action="append",
        default=[],
        help="Gene prefix to exclude; may be supplied repeatedly.",
    )
    parser.add_argument("--mode", default="cells", choices=["cells", "clusters"])
    parser.add_argument("--num-epochs", type=int, default=1000)
    parser.add_argument("--device", default="cpu")
    parser.add_argument("--density-prior", default="rna_count_based")
    parser.add_argument(
        "--map-input",
        type=Path,
        default=None,
        help="Optional saved Tangram mapping h5ad. If set, skip training.",
    )
    parser.add_argument(
        "--plot-cell-type",
        action="append",
        default=[],
        help="Subset label to plot; may be supplied repeatedly.",
    )

    return parser.parse_args()


def load_inputs(sc_path, spatial_path):
    """Load the annotated reference and clustered spatial AnnData objects."""
    if not sc_path.is_file():
        raise FileNotFoundError(f"scRNA-seq input not found: {sc_path}")

    if not spatial_path.is_file():
        raise FileNotFoundError(f"Spatial input not found: {spatial_path}")

    adata_sc = sc.read_h5ad(sc_path)
    adata_sp = sc.read_h5ad(spatial_path)

    return adata_sc, adata_sp


def validate_inputs(adata_sc, adata_sp, celltype_key):
    """Validate metadata, gene names, and coordinates required by Tangram."""
    if adata_sc.n_obs == 0 or adata_sc.n_vars == 0:
        raise ValueError("The scRNA-seq AnnData object is empty.")

    if adata_sp.n_obs == 0 or adata_sp.n_vars == 0:
        raise ValueError("The spatial AnnData object is empty")

    if celltype_key not in adata_sc.obs.columns:
        raise KeyError(
            f"Cell-type column '{celltype_key}' not found. Available columns: {list(adata_sc.obs.columns)}"
        )

    if "spatial" not in adata_sp.obsm:
        raise KeyError("Spatial coordinates not found in adata_sp.obsm['spatial']")

    if not adata_sc.var_names.is_unique:
        raise ValueError("Reference gene names are not unique.")

    if not adata_sp.var_names.is_unique:
        raise ValueError("Spatial gene names are not unique.")

    labels = adata_sc.obs[celltype_key]

    if labels.isna().all():
        raise ValueError(f"All values in reference column '{celltype_key}' are missing")

    if labels.astype(str).str.strip().eq("").all():
        raise ValueError(f"Reference column '{celltype_key}' contains no usable labels")

    coordinates = np.asarray(adata_sp.obsm["spatial"])

    if coordinates.shape[0] != adata_sp.n_obs:
        raise ValueError(
            "The number of spatial coordinate rows does not match the number of spatial observations."
        )

    if coordinates.ndim != 2 or coordinates.shape[1] < 2:
        raise ValueError("Spatial coordinates must be a two-dimensional array with at least two columns.")

    if not np.isfinite(coordinates).all():
        raise ValueError("Spatial coordinates contain NaN or infinite values")


def get_shared_genes(adata_sc, adata_sp, exclude_prefixes):
    """Return sorted shared genes after optional prefix filtering."""
    shared_genes = np.intersect1d(adata_sc.var_names, adata_sp.var_names)

    normalized_prefixes = tuple(prefix.upper() for prefix in exclude_prefixes)

    if normalized_prefixes:
        shared_genes = np.array(
            [gene for gene in shared_genes if not gene.upper().startswith(normalized_prefixes)]
        )

    return shared_genes


def write_shared_genes(shared_genes, outdir):
    """Write the shared-gene audit table to the requested output directory."""
    outdir.mkdir(parents=True, exist_ok=True)

    output_path = outdir / "shared_genes.tsv"

    with output_path.open("w", encoding="utf-8") as handle:
        handle.write("gene\n")
        for gene in shared_genes:
            handle.write(f"{gene}\n")

    return output_path


def run_mapping(adata_sc, adata_sp, celltype_key, mode, num_epochs, device, density_prior):
    """Learn a Tangram mapping from reference cells onto spatial spots."""
    return tg.map_cells_to_space(
        adata_sc,
        adata_sp,
        mode=mode,
        cluster_label=celltype_key if mode == "clusters" else None,
        num_epochs=num_epochs,
        device=device,
        density_prior=density_prior,
    )


def plot_celltype_maps(adata_sp, cell_types, outdir):
    """Write one spatial proportion map per selected cell type."""
    figdir = Path(outdir) / "figures"
    figdir.mkdir(parents=True, exist_ok=True)

    props = adata_sp.obsm["tangram_ct_pred"]

    for cell_type in cell_types:
        if cell_type not in props.columns:
            raise KeyError(f"{cell_type} not in proportion table: {list(props.columns)}")

        adata_sp.obs[cell_type] = props[cell_type].to_numpy()

        sc.pl.spatial(adata_sp, color=cell_type, img_key="hires", show=False)
        safe_name = cell_type.replace("/", "_").replace("+", "pos")
        plt.savefig(figdir / f"{safe_name}.png", dpi=150, bbox_inches="tight")
        plt.close()
        print("Wrote", figdir / f"{safe_name}.png")


def main():
    """Run input loading and validation for the Tangram workflow."""
    args = parse_args()
    adata_sc, adata_sp = load_inputs(args.sc_input, args.spatial_input)

    validate_inputs(adata_sc, adata_sp, args.celltype_key)

    print("Reference shape:", adata_sc.shape)
    print("Spatial shape:", adata_sp.shape)

    shared_genes = get_shared_genes(adata_sc, adata_sp, args.exclude_prefix)

    if len(shared_genes) < args.min_shared_genes:
        raise ValueError(
            f"Only {len(shared_genes)} shared genes remain; minimum required is {args.min_shared_genes}."
        )

    shared_genes_path = write_shared_genes(shared_genes, args.outdir)

    print("Exact shared genes:", len(shared_genes))
    print("Shared-gene audit:", shared_genes_path)

    # from collections import defaultdict
    # nonzero_counts = {}

    # for gene in shared_genes:
    #     sc_index = adata_sc.var_names.get_loc(gene)
    #     sp_index = adata_sp.var_names.get_loc(gene)

    #     if sparse.issparse(adata_sc.X):
    #         sc_nonzero = adata_sc.X[:, sc_index].getnnz()
    #     else:
    #         sc_nonzero = np.count_nonzero(adata_sc.X[:, sc_index])
    #     if sparse.issparse(adata_sp.X):
    #         sp_nonzero = adata_sp.X[:, sp_index].getnnz()
    #     else:
    #         sp_nonzero = np.count_nonzero(adata_sp.X[:, sp_index])

    #     nonzero_counts[gene] = (sc_nonzero, sp_nonzero)

    # Tanagram Preprocessing
    tg.pp_adatas(adata_sc, adata_sp, genes=shared_genes, gene_to_lowercase=False)

    print("Tanagram preprocessing complete")

    training_genes = np.asarray(adata_sc.uns["training_genes"])

    print("Tangram training genes:", len(training_genes))
    print("Shared genes excluded by Tangram", len(shared_genes) - len(training_genes))

    if args.map_input is not None:
        if not args.map_input.is_file():
            raise FileNotFoundError(f"Saved mapping not found: {args.map_input}")
        ad_map = sc.read_h5ad(args.map_input)
        print("Loaded saved mapping:", args.map_input)
    else:
        ad_map = run_mapping(
            adata_sc,
            adata_sp,
            args.celltype_key,
            args.mode,
            args.num_epochs,
            args.device,
            args.density_prior,
        )
        ad_map.write_h5ad(args.outdir / "tangram_map.h5ad")
        print("Wrote", args.outdir / "tangram_map.h5ad")

    print("Mapping shape:", ad_map.shape)

    tg.project_cell_annotations(ad_map, adata_sp, annotation=args.celltype_key)

    props = adata_sp.obsm["tangram_ct_pred"].copy()
    row_sums = props.sum(axis=1)
    props = props.div(row_sums.replace(0, np.nan), axis=0)
    adata_sp.obsm["tangram_ct_pred"] = props
    props.to_csv(args.outdir / "proportions.tsv", sep="\t")

    print("Proportions shape:", props.shape)
    print("Cell-type columns:", props.shape[1])
    print("Proportion columns:", list(props.columns))
    print("Spots with zero total score:", int((row_sums == 0).sum()))
    finite_sums = props.sum(axis=1).dropna()
    print("Finite row-sum min/max:", float(finite_sums.min()), float(finite_sums.max()))

    cell_types = args.plot_cell_type or list(KEY_CELL_TYPES)
    plot_celltype_maps(adata_sp, cell_types, args.outdir)
    adata_sp.write_h5ad(args.outdir / "spatial.deconvolved.h5ad")


if __name__ == "__main__":
    raise SystemExit(main())
