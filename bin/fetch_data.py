#!/usr/bin/env python3
"""Download the two public datasets this pipeline runs on, and describe what arrived.

    python bin/fetch_data.py --all
    python bin/fetch_data.py --reference          # just the scRNA-seq reference
    python bin/fetch_data.py --visium --inspect   # just Visium, then print a structure summary

Datasets
--------
reference : Integrated secondary lymphoid organ scRNA-seq atlas from Kleshchevnikov et al.,
            Nat Biotechnol 2022 (the cell2location paper). 73,260 cells, 10,237 genes, 34 curated
            immune and stromal cell types, pooled from lymph node, spleen and tonsil. The pooling
            matters: it captures germinal-centre populations that a single lymph node sample is too
            shallow to resolve, and germinal centres are the most interesting thing in this tissue.
            ~1-2 GB download.

visium    : 10x Genomics V1_Human_Lymph_Node. One Visium slide, ~4,000 spots under tissue, with the
            paired H&E image. Spots are ~55 um and therefore mixtures of cells, which is the entire
            reason the deconvolution step exists.

Both are used together in the canonical cell2location and Tangram tutorials, so published figures on
this exact pairing exist to validate your results against. That is a deliberate choice: it converts
"my output looks plausible" into "my output matches a peer-reviewed result".

Output location
---------------
Defaults to $STPIPE_DATA (set by setup/bootstrap_wsl.sh to a native Linux path) and falls back to
./data. Downloads are skipped if the target already exists, so this is safe to re-run; pass --force
to redownload.
"""

from __future__ import annotations

import argparse
import os
import shutil
import sys
import time
import urllib.request
from pathlib import Path

REFERENCE_URL = (
    "https://cell2location.cog.sanger.ac.uk/paper/integrated_lymphoid_organ_scrna/"
    "RegressionNBV4Torch_57covariates_73260cells_10237genes/sc.h5ad"
)
VISIUM_SAMPLE_ID = "V1_Human_Lymph_Node"


# --------------------------------------------------------------------------------------- helpers


def default_data_dir() -> Path:
    return Path(os.environ.get("STPIPE_DATA") or "./data").expanduser().resolve()


def human(n: float) -> str:
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if abs(n) < 1024.0:
            return f"{n:,.1f} {unit}"
        n /= 1024.0
    return f"{n:,.1f} PB"


def download(url: str, dest: Path, force: bool = False) -> Path:
    """Stream a URL to dest, atomically, with progress on stderr."""
    if dest.exists() and not force:
        print(f"[skip] {dest} already exists ({human(dest.stat().st_size)}). Use --force to redownload.")
        return dest

    dest.parent.mkdir(parents=True, exist_ok=True)
    # Download to a .part file and rename only on success. Without this, an interrupted download
    # leaves a truncated file that looks complete on the next run and fails much later with a
    # baffling HDF5 error.
    tmp = dest.with_suffix(dest.suffix + ".part")
    print(f"[get ] {url}\n       -> {dest}")

    start = time.monotonic()

    def progress(block_num: int, block_size: int, total_size: int) -> None:
        got = block_num * block_size
        elapsed = max(time.monotonic() - start, 1e-6)
        rate = got / elapsed
        if total_size > 0:
            pct = min(100.0, got * 100.0 / total_size)
            msg = f"       {pct:5.1f}%  {human(got)} / {human(total_size)}  at {human(rate)}/s"
        else:
            msg = f"       {human(got)}  at {human(rate)}/s"
        sys.stderr.write("\r" + msg + " " * 8)
        sys.stderr.flush()

    try:
        urllib.request.urlretrieve(url, tmp, reporthook=progress)
    except BaseException:
        tmp.unlink(missing_ok=True)
        raise
    sys.stderr.write("\n")

    tmp.replace(dest)
    print(f"[done] {dest} ({human(dest.stat().st_size)}) in {time.monotonic() - start:.0f}s")
    return dest


# --------------------------------------------------------------------------------------- fetchers


def fetch_reference(data_dir: Path, force: bool = False) -> Path:
    dest = data_dir / "reference" / "sc.h5ad"
    return download(REFERENCE_URL, dest, force=force)


def fetch_visium(data_dir: Path, force: bool = False) -> Path:
    """Fetch the Visium sample via scanpy, tolerating the loader rename across versions.

    scanpy renamed `datasets.visium_sge` to `datasets.visium` in 1.10 and dropped the old name later,
    and squidpy has its own loader. Rather than pinning ourselves to one API, try each in turn. Small
    shims like this are the difference between a repo that still runs next year and one that does not.
    """
    import scanpy as sc

    target = data_dir / "visium" / VISIUM_SAMPLE_ID
    if target.exists() and any(target.iterdir()) and not force:
        print(f"[skip] {target} already populated. Use --force to redownload.")
        return target

    if force and target.exists():
        shutil.rmtree(target)

    base_dir = data_dir / "visium"
    base_dir.mkdir(parents=True, exist_ok=True)

    attempts = []

    loader = getattr(sc.datasets, "visium", None)
    if loader is not None:
        attempts.append(("scanpy.datasets.visium", lambda: loader(VISIUM_SAMPLE_ID, base_dir=base_dir)))

    legacy = getattr(sc.datasets, "visium_sge", None)
    if legacy is not None:
        attempts.append(
            ("scanpy.datasets.visium_sge", lambda: legacy(sample_id=VISIUM_SAMPLE_ID, base_dir=base_dir))
        )

    def _squidpy():
        import squidpy as sq

        return sq.datasets.visium(VISIUM_SAMPLE_ID, base_dir=base_dir)

    attempts.append(("squidpy.datasets.visium", _squidpy))

    last_error: Exception | None = None
    for name, fn in attempts:
        try:
            print(f"[get ] Visium {VISIUM_SAMPLE_ID} via {name}")
            adata = fn()
            print(f"[done] {adata.n_obs:,} spots x {adata.n_vars:,} genes")
            break
        except TypeError as exc:
            # Signature drift between versions, e.g. base_dir not accepted. Try the next loader.
            print(f"[warn] {name} signature mismatch: {exc}")
            last_error = exc
        except Exception as exc:
            print(f"[warn] {name} failed: {type(exc).__name__}: {exc}")
            last_error = exc
    else:
        raise RuntimeError(f"All Visium loaders failed. Last error: {last_error!r}")

    if not target.exists():
        # Some loader versions place the sample directly in base_dir rather than base_dir/sample_id.
        candidates = [p for p in base_dir.iterdir() if p.is_dir() and p.name != VISIUM_SAMPLE_ID]
        raise RuntimeError(
            f"Expected {target} to exist after download. Found instead: "
            f"{[str(p) for p in candidates]}. Move the sample into place or adjust --visium-dir."
        )

    return target


# --------------------------------------------------------------------------------------- inspect


def inspect_reference(path: Path) -> None:
    """Print the structure of the reference. Read this output carefully before Phase 1.

    Two things you need from it: which .obs column holds the curated cell type labels (params
    sc_celltype_key in nextflow.config assumes 'Subset' -- verify), and whether mitochondrial genes
    are present at all. Many processed public references have already had them removed, in which case
    a mito-fraction QC filter is vacuous and you should say so rather than silently filtering nothing.
    """
    import anndata as ad

    print(f"\n{'=' * 78}\nREFERENCE: {path}\n{'=' * 78}")
    # backed='r' keeps the matrix on disk; loading 73k cells densely is how you OOM a laptop.
    adata = ad.read_h5ad(path, backed="r")
    print(f"shape            : {adata.n_obs:,} cells x {adata.n_vars:,} genes")
    print(f"X dtype          : {adata.X.dtype if adata.X is not None else 'None'}")
    print(f"layers           : {list(adata.layers.keys())}")
    print(f"obsm             : {list(adata.obsm.keys())}")

    print(f"\n.obs columns ({len(adata.obs.columns)}):")
    for col in adata.obs.columns:
        series = adata.obs[col]
        n_unique = int(series.nunique(dropna=True))
        is_categorical = str(series.dtype) in ("object", "category", "string")
        flag = "   <-- candidate cell type label" if is_categorical and 5 <= n_unique <= 60 else ""
        print(f"  {col:<32} n_unique={n_unique:<6} dtype={str(series.dtype):<10}{flag}")

    print(f"\n.var columns     : {list(adata.var.columns)}")
    print(f"first 8 var_names: {list(adata.var_names[:8])}")

    looks_ensembl = str(adata.var_names[0]).startswith("ENSG")
    print(f"gene ID style    : {'Ensembl IDs' if looks_ensembl else 'gene symbols'}")

    mito = [g for g in adata.var_names if str(g).upper().startswith(("MT-", "MT."))]
    print(f"mitochondrial genes present: {len(mito)}  {mito[:10]}")
    if not mito:
        print("  NOTE: no MT- genes found. A percent-mito QC filter would do nothing here.")
        print("        Document this in docs/methods.md rather than filtering on an empty set.")

    adata.file.close()


def read_visium_any(path: Path):
    """Read a Visium sample directory, tolerating scanpy's deprecation of read_visium.

    scanpy.read_visium was deprecated in favour of squidpy.read.visium. Try squidpy first since it is
    the maintained path, then fall back so this works across a range of installed versions.
    """
    try:
        import squidpy as sq

        return sq.read.visium(path)
    except Exception:
        import scanpy as sc

        return sc.read_visium(path)


def inspect_visium(path: Path) -> None:
    print(f"\n{'=' * 78}\nVISIUM: {path}\n{'=' * 78}")
    print("directory contents:")
    for p in sorted(path.rglob("*")):
        if p.is_file():
            print(f"  {p.relative_to(path)}  ({human(p.stat().st_size)})")

    adata = read_visium_any(path)
    adata.var_names_make_unique()
    print(f"\nshape            : {adata.n_obs:,} spots x {adata.n_vars:,} genes")
    print(f".obs columns     : {list(adata.obs.columns)}")
    print(f".obsm            : {list(adata.obsm.keys())}  (spatial coordinates live in obsm['spatial'])")
    print(f".uns['spatial']  : {list(adata.uns.get('spatial', {}).keys())}")
    print(f"first 8 var_names: {list(adata.var_names[:8])}")
    looks_ensembl = str(adata.var_names[0]).startswith("ENSG")
    print(f"gene ID style    : {'Ensembl IDs' if looks_ensembl else 'gene symbols'}")
    print(
        "\nCOMPARE the gene ID style above against the reference. If they differ, that mismatch is\n"
        "the single most likely source of a silent failure in Phase 4 deconvolution."
    )


# --------------------------------------------------------------------------------------- cli


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--all", action="store_true", help="fetch both datasets")
    parser.add_argument("--reference", action="store_true", help="fetch the scRNA-seq reference h5ad")
    parser.add_argument("--visium", action="store_true", help="fetch the Visium lymph node sample")
    parser.add_argument(
        "--outdir",
        type=Path,
        default=default_data_dir(),
        help="destination root (default: $STPIPE_DATA, else ./data)",
    )
    parser.add_argument("--force", action="store_true", help="redownload even if files exist")
    parser.add_argument(
        "--inspect",
        action="store_true",
        help="print a structure summary of whatever was fetched (recommended on first run)",
    )
    args = parser.parse_args()

    if not (args.all or args.reference or args.visium):
        parser.error("choose at least one of --all, --reference, --visium")

    want_ref = args.all or args.reference
    want_vis = args.all or args.visium
    data_dir = args.outdir.expanduser().resolve()
    print(f"data dir: {data_dir}")

    ref_path = fetch_reference(data_dir, force=args.force) if want_ref else None
    vis_path = fetch_visium(data_dir, force=args.force) if want_vis else None

    if args.inspect:
        if ref_path:
            inspect_reference(ref_path)
        if vis_path:
            inspect_visium(vis_path)

    print("\nNext step: PROJECT_PLAN.md Phase 1 -- write bin/sc_qc.py.")
    if ref_path:
        print(f"  --sc_reference {ref_path}")
    if vis_path:
        print(f"  --visium_dir   {vis_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
