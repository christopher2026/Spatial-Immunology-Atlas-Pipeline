#!/usr/bin/env python3
"""Build the committed CI fixtures under assets/test_data/.

This is not a random downsample. The test profile has to exercise the same DAG as production
without becoming either (a) a popularity sample that drops rare Subset labels or (b) a scatter of
isolated Visium spots with no spatial neighbors.

Example (WSL, after fetch_data.py):

    python bin/make_test_data.py \\
        --sc-input "$STPIPE_DATA/reference/sc.h5ad" \\
        --visium-input "$STPIPE_DATA/visium/V1_Human_Lymph_Node" \\
        --outdir assets/test_data
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
from scipy import sparse

REPO_ROOT = Path(__file__).resolve().parents[1]


def default_data_dir() -> Path:
    return Path(os.environ.get("STPIPE_DATA") or "./data").expanduser().resolve()


def human(n: float) -> str:
    value = float(n)
    for unit in ("B", "KB", "MB", "GB"):
        if abs(value) < 1024.0:
            return f"{value:,.1f} {unit}"
        value /= 1024.0
    return f"{value:,.1f} TB"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    data_dir = default_data_dir()
    parser.add_argument(
        "--sc-input",
        type=Path,
        default=data_dir / "reference" / "sc.h5ad",
        help="Full reference AnnData (default: $STPIPE_DATA/reference/sc.h5ad)",
    )
    parser.add_argument(
        "--visium-input",
        type=Path,
        default=data_dir / "visium" / "V1_Human_Lymph_Node",
        help="Full Visium sample directory",
    )
    parser.add_argument(
        "--outdir",
        type=Path,
        default=REPO_ROOT / "assets" / "test_data",
        help="Destination directory (default: assets/test_data)",
    )
    parser.add_argument("--n-cells", type=int, default=500, help="Target number of reference cells")
    parser.add_argument("--n-spots", type=int, default=200, help="Target number of in-tissue Visium spots")
    parser.add_argument("--celltype-key", default="Subset", help="obs column to stratify on")
    parser.add_argument("--seed", type=int, default=0, help="RNG seed for cell sampling")
    parser.add_argument(
        "--force", action="store_true", help="Overwrite existing sc_subsampled.h5ad / visium_subsampled"
    )
    return parser.parse_args()


def stratified_positions(labels: pd.Series, n_cells: int, rng: np.random.Generator) -> np.ndarray:
    """Equal-as-possible sample that keeps every category.

    Rare types with fewer cells than the per-type quota are taken in full; leftover slots go to
    types that still have unused cells. Uniform sc.pp.subsample is the thing this replaces.
    """
    values = labels.astype(str).to_numpy()
    types = sorted(pd.unique(values))
    n_types = len(types)
    if n_cells < n_types:
        raise ValueError(
            f"--n-cells={n_cells} is smaller than the {n_types} {labels.name} categories; "
            "every type must appear at least once."
        )

    shuffled = {t: rng.permutation(np.flatnonzero(values == t)) for t in types}
    quota, extra = divmod(n_cells, n_types)
    chosen: list[np.ndarray] = []
    leftover = 0
    unused: list[tuple[str, np.ndarray]] = []

    # Protect rare types first so their shortfall becomes leftover for abundant types.
    order = sorted(types, key=lambda t: len(shuffled[t]))
    for i, t in enumerate(order):
        want = quota + (1 if i < extra else 0)
        pool = shuffled[t]
        take = min(want, len(pool))
        chosen.append(pool[:take])
        leftover += want - take
        if take < len(pool):
            unused.append((t, pool[take:]))

    unused.sort(key=lambda item: len(item[1]), reverse=True)
    while leftover > 0 and unused:
        t, pool = unused.pop(0)
        chosen.append(pool[:1])
        leftover -= 1
        if len(pool) > 1:
            unused.append((t, pool[1:]))
            unused.sort(key=lambda item: len(item[1]), reverse=True)

    idx = np.concatenate(chosen)
    if idx.size < n_cells:
        print(
            f"WARNING: only {idx.size} cells available after protecting every {labels.name}; requested {n_cells}."
        )
    return np.sort(idx)


def subsample_reference(path: Path, celltype_key: str, n_cells: int, seed: int) -> tuple[ad.AnnData, int]:
    print(f"Loading reference (backed) from {path}")
    # backed='r' avoids pulling 73k cells into RAM just to choose indices.
    backed = ad.read_h5ad(path, backed="r")
    if celltype_key not in backed.obs.columns:
        backed.file.close()
        raise KeyError(f"{path} has no obs[{celltype_key!r}]. Columns: {list(backed.obs.columns)}")

    labels = backed.obs[celltype_key]
    n_types_in = int(labels.nunique(dropna=True))
    rng = np.random.default_rng(seed)
    idx = stratified_positions(labels, n_cells, rng)
    barcodes = backed.obs_names[idx].tolist()
    print(f"Materialising {len(barcodes)} stratified cells x {backed.n_vars} genes")
    view = backed[barcodes].to_memory()
    backed.file.close()
    # Drop embeddings / graphs from the 73k-cell object; sc_qc recomputes everything from X.
    slim = ad.AnnData(X=view.X.copy(), obs=view.obs.copy(), var=view.var.copy())
    slim.obs_names = view.obs_names
    slim.var_names = view.var_names
    slim.uns["make_test_data"] = {
        "source": str(path),
        "n_cells": int(slim.n_obs),
        "celltype_key": celltype_key,
        "n_types_in": n_types_in,
        "seed": seed,
    }
    del view
    return slim, n_types_in


def read_visium(path: Path):
    try:
        import squidpy as sq

        return sq.read.visium(path)
    except Exception:
        import scanpy as sc

        return sc.read_visium(path)


def in_tissue_mask(adata: ad.AnnData) -> np.ndarray:
    if "in_tissue" not in adata.obs.columns:
        return np.ones(adata.n_obs, dtype=bool)
    return adata.obs["in_tissue"].astype(str).isin(["1", "True", "true"]).to_numpy()


def contiguous_spot_indices(adata: ad.AnnData, n_spots: int) -> np.ndarray:
    """Keep a compact lattice around the in-tissue centroid, not a random barcode sample."""
    tissue = np.flatnonzero(in_tissue_mask(adata))
    if tissue.size < n_spots:
        raise ValueError(f"Only {tissue.size} in-tissue spots; cannot select {n_spots}.")
    for col in ("array_row", "array_col"):
        if col not in adata.obs.columns:
            raise KeyError(f"Visium obs is missing {col}; cannot crop a contiguous region.")

    grid = adata.obs.iloc[tissue][["array_row", "array_col"]].to_numpy(dtype=float)
    center = np.median(grid, axis=0)
    dist = np.linalg.norm(grid - center, axis=1)
    order = np.argsort(dist, kind="stable")
    return np.sort(tissue[order[:n_spots]])


def find_positions_file(spatial_dir: Path) -> Path:
    for name in ("tissue_positions_list.csv", "tissue_positions.csv"):
        candidate = spatial_dir / name
        if candidate.exists():
            return candidate
    raise FileNotFoundError(f"No tissue_positions*.csv under {spatial_dir}")


def read_positions(path: Path) -> tuple[pd.DataFrame, bool]:
    with path.open() as handle:
        first_cell = handle.readline().split(",")[0].strip().strip('"')
    has_header = first_cell.lower() == "barcode"
    coords = pd.read_csv(path, header=0 if has_header else None, index_col=0)
    coords.index = coords.index.astype(str)
    coords.columns = ["in_tissue", "array_row", "array_col", "pxl_row_in_fullres", "pxl_col_in_fullres"]
    return coords, has_header


def write_positions(dest: Path, coords: pd.DataFrame) -> None:
    """Classic Space Ranger v1 list format (no header), as named in assets/test_data/README.md."""
    out = coords.copy()
    out.index.name = None
    out.to_csv(dest, header=False)


def _bytes_array(values) -> np.ndarray:
    return np.asarray(values).astype(str).astype("S")


def visium_library_id(adata: ad.AnnData, default: str = "V1_Human_Lymph_Node") -> str:
    spatial = adata.uns.get("spatial") or {}
    return str(next(iter(spatial))) if spatial else default


def copy_h5_root_attrs(source_h5: Path | None, dest_handle: h5py.File, library_id: str) -> None:
    """Space Ranger stores library_ids on the file root; squidpy/scanpy pop that attr on read."""
    if source_h5 is not None and source_h5.exists():
        with h5py.File(source_h5, "r") as src:
            for key, value in src.attrs.items():
                dest_handle.attrs[key] = value
    if "library_ids" not in dest_handle.attrs:
        dest_handle.attrs["library_ids"] = np.array([library_id], dtype="S")


def write_10x_h5(
    adata: ad.AnnData,
    dest: Path,
    *,
    source_h5: Path | None = None,
    library_id: str = "V1_Human_Lymph_Node",
) -> None:
    """Write Cell Ranger v3 feature-barcode HDF5 (genes x barcodes CSC == AnnData CSR transposed)."""
    matrix = adata.X
    matrix = sparse.csr_matrix(matrix) if not sparse.issparse(matrix) else matrix.tocsr()
    data = (
        np.rint(matrix.data).astype(np.int32)
        if matrix.data.size and np.allclose(matrix.data, np.round(matrix.data))
        else matrix.data.astype(np.float32)
    )

    gene_ids = adata.var["gene_ids"] if "gene_ids" in adata.var.columns else adata.var_names
    feature_types = (
        adata.var["feature_types"]
        if "feature_types" in adata.var.columns
        else np.repeat("Gene Expression", adata.n_vars)
    )
    genome = adata.var["genome"] if "genome" in adata.var.columns else np.repeat("unknown", adata.n_vars)

    dest.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(dest, "w") as handle:
        copy_h5_root_attrs(source_h5, handle, library_id)
        grp = handle.create_group("matrix")
        gzip = {"compression": "gzip", "compression_opts": 4}
        grp.create_dataset("barcodes", data=_bytes_array(adata.obs_names), **gzip)
        grp.create_dataset("data", data=data, **gzip)
        grp.create_dataset("indices", data=matrix.indices.astype(np.int32), **gzip)
        grp.create_dataset("indptr", data=matrix.indptr.astype(np.int32), **gzip)
        grp.create_dataset("shape", data=np.array([adata.n_vars, adata.n_obs], dtype=np.int32))
        features = grp.create_group("features")
        features.create_dataset("id", data=_bytes_array(gene_ids), **gzip)
        features.create_dataset("name", data=_bytes_array(adata.var_names), **gzip)
        features.create_dataset("feature_type", data=_bytes_array(feature_types), **gzip)
        features.create_dataset("genome", data=_bytes_array(genome), **gzip)
        features.create_dataset("_all_tag_keys", data=_bytes_array(["genome"]))


def write_visium_dir(source_dir: Path, adata: ad.AnnData, dest_dir: Path) -> None:
    spatial_src = source_dir / "spatial"
    spatial_dest = dest_dir / "spatial"
    if spatial_dest.exists():
        shutil.rmtree(spatial_dest)
    spatial_dest.mkdir(parents=True, exist_ok=True)

    source_h5 = source_dir / "filtered_feature_bc_matrix.h5"
    write_10x_h5(
        adata,
        dest_dir / "filtered_feature_bc_matrix.h5",
        source_h5=source_h5 if source_h5.exists() else None,
        library_id=visium_library_id(adata),
    )

    coords, _has_header = read_positions(find_positions_file(spatial_src))
    selected = adata.obs_names.astype(str)
    missing = [b for b in selected if b not in coords.index]
    if missing:
        raise KeyError(
            f"{len(missing)} selected barcodes are absent from tissue positions, e.g. {missing[:3]}"
        )
    write_positions(spatial_dest / "tissue_positions_list.csv", coords.loc[list(selected)])

    scalef_src = spatial_src / "scalefactors_json.json"
    scalef = json.loads(scalef_src.read_text())
    lowres_src = spatial_src / "tissue_lowres_image.png"
    hires_src = spatial_src / "tissue_hires_image.png"
    if not lowres_src.exists():
        raise FileNotFoundError(f"Missing {lowres_src}")
    shutil.copy2(lowres_src, spatial_dest / "tissue_lowres_image.png")

    # Production modules plot img_res_key="hires". Shipping the full hires PNG often blows the git
    # size budget; aliasing lowres and matching the hires scale factor keeps overlays aligned.
    hires_limit = 2 * 1024 * 1024
    if hires_src.exists() and hires_src.stat().st_size <= hires_limit:
        shutil.copy2(hires_src, spatial_dest / "tissue_hires_image.png")
    else:
        shutil.copy2(lowres_src, spatial_dest / "tissue_hires_image.png")
        if "tissue_lowres_scalef" in scalef:
            scalef["tissue_hires_scalef"] = scalef["tissue_lowres_scalef"]
            print(
                "Copied lowres H&E as hires and set tissue_hires_scalef = tissue_lowres_scalef (size budget)."
            )
    (spatial_dest / "scalefactors_json.json").write_text(json.dumps(scalef, indent=2) + "\n")


def print_tree_sizes(root: Path) -> None:
    total = 0
    print(f"\nWritten under {root}:")
    for path in sorted(p for p in root.rglob("*") if p.is_file()):
        size = path.stat().st_size
        total += size
        print(f"  {path.relative_to(root)}  ({human(size)})")
    print(f"  total  {human(total)}")
    if total > 50 * 1024 * 1024:
        raise SystemExit(f"Test data is {human(total)}; subsample harder instead of using Git LFS.")


def print_sc_summary(adata: ad.AnnData, celltype_key: str, n_types_in: int) -> None:
    counts = adata.obs[celltype_key].astype(str).value_counts()
    print("\n--- scRNA-seq subsample ---")
    print(f"shape            : {adata.n_obs:,} cells x {adata.n_vars:,} genes")
    print(f"{celltype_key} categories : {counts.size} (source had {n_types_in})")
    print(
        f"cells per type   : min={int(counts.min())}  median={float(counts.median()):.1f}  max={int(counts.max())}"
    )
    if counts.size != n_types_in:
        raise SystemExit(f"Stratification dropped types: in={n_types_in} out={counts.size}")


def print_spot_summary(adata: ad.AnnData) -> None:
    rows = adata.obs["array_row"].astype(int)
    cols = adata.obs["array_col"].astype(int)
    print("\n--- Visium subsample ---")
    print(f"shape            : {adata.n_obs:,} spots x {adata.n_vars:,} genes")
    print(f"array_row range  : {int(rows.min())}–{int(rows.max())}  (span {int(rows.max() - rows.min())})")
    print(f"array_col range  : {int(cols.min())}–{int(cols.max())}  (span {int(cols.max() - cols.min())})")
    print(f"in_tissue        : {int(in_tissue_mask(adata).sum())}/{adata.n_obs}")


def main() -> int:
    args = parse_args()
    sc_input = args.sc_input.expanduser().resolve()
    visium_input = args.visium_input.expanduser().resolve()
    outdir = args.outdir.expanduser().resolve()
    sc_out = outdir / "sc_subsampled.h5ad"
    vis_out = outdir / "visium_subsampled"

    if not sc_input.exists():
        raise FileNotFoundError(f"Reference not found: {sc_input}")
    if not visium_input.is_dir():
        raise FileNotFoundError(f"Visium directory not found: {visium_input}")
    if (sc_out.exists() or vis_out.exists()) and not args.force:
        raise SystemExit(f"{outdir} already has outputs. Pass --force to overwrite.")

    outdir.mkdir(parents=True, exist_ok=True)

    slim, n_types_in = subsample_reference(sc_input, args.celltype_key, args.n_cells, args.seed)
    slim.write_h5ad(sc_out, compression="gzip")
    print_sc_summary(slim, args.celltype_key, n_types_in)
    print(f"Wrote {sc_out} ({human(sc_out.stat().st_size)})")

    print(f"\nLoading Visium from {visium_input}")
    vis = read_visium(visium_input)
    vis.var_names_make_unique()
    spot_idx = contiguous_spot_indices(vis, args.n_spots)
    vis_sub = vis[spot_idx].copy()
    print_spot_summary(vis_sub)

    if vis_out.exists():
        shutil.rmtree(vis_out)
    vis_out.mkdir(parents=True, exist_ok=True)
    write_visium_dir(visium_input, vis_sub, vis_out)

    reloaded = read_visium(vis_out)
    if reloaded.n_obs != vis_sub.n_obs:
        raise SystemExit(f"Reload n_obs={reloaded.n_obs} != written {vis_sub.n_obs}")
    if "spatial" not in reloaded.obsm:
        raise SystemExit("Reloaded Visium is missing obsm['spatial']")
    print(
        f"Reload check     : {reloaded.n_obs} spots, obsm['spatial'] present, "
        f"images={list(reloaded.uns.get('spatial', {}).keys())}"
    )

    print_tree_sizes(outdir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
