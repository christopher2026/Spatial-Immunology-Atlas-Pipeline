#!/usr/bin/env python3
"""Assemble pipeline figures and tables into a single standalone HTML report."""

from __future__ import annotations

import argparse
import base64
from pathlib import Path

import pandas as pd
from jinja2 import Environment, FileSystemLoader, select_autoescape


def parse_args():
    """Parse command-line options for report assembly."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-dir", type=Path, required=True, help="Pipeline results directory.")
    parser.add_argument("--output", type=Path, required=True, help="Output HTML path.")
    parser.add_argument("--template", type=Path, default=None, help="Jinja2 template path.")
    parser.add_argument("--sample-id", default="V1_Human_Lymph_Node")
    return parser.parse_args()


def png_to_data_uri(path: Path) -> str:
    """Encode a PNG as a data URI so the HTML has no relative image paths."""
    encoded = base64.b64encode(path.read_bytes()).decode("ascii")
    return f"data:image/png;base64,{encoded}"


def collect_pngs(directory: Path) -> list[Path]:
    """Return PNG paths in a directory, if it exists."""
    if not directory.is_dir():
        return []
    return sorted(directory.glob("*.png"))


def figure_entries(paths: list[Path]) -> list[dict]:
    """Build template figure dicts from PNG paths."""
    return [{"caption": path.stem.replace("_", " "), "b64": png_to_data_uri(path)} for path in paths]


def tsv_to_html(path: Path, max_rows: int = 12) -> str | None:
    """Render a small TSV as an HTML table."""
    if not path.is_file():
        return None
    table = pd.read_csv(path, sep="\t", index_col=0)
    return table.head(max_rows).to_html(classes="data", border=0)


def collect_versions(results_dir: Path) -> dict[str, str]:
    """Merge any versions.yml files published under results/."""
    versions: dict[str, str] = {}
    for path in sorted(results_dir.glob("**/versions.yml")):
        current_tool = path.parent.name
        for line in path.read_text(encoding="utf-8").splitlines():
            stripped = line.strip()
            if stripped.endswith(":") and not stripped.startswith(" "):
                current_tool = stripped[:-1]
            elif ":" in stripped:
                tool, value = stripped.split(":", 1)
                versions[f"{current_tool}.{tool.strip()}"] = value.strip()
    return versions


def build_context(results_dir: Path, sample_id: str) -> dict:
    """Collect summary fields, sections, and versions for the template."""
    sections = [
        {
            "title": "Single-cell QC and clustering",
            "description": "Reference QC and independent Leiden clusters compared with curated Subset labels.",
            "figures": figure_entries(collect_pngs(results_dir / "sc_qc" / "figures"))
            + figure_entries(collect_pngs(results_dir / "sc_cluster" / "figures"))
            + figure_entries(collect_pngs(results_dir / "preprocess_figures")),
            "tables": [],
        },
        {
            "title": "Spatial QC and clustering",
            "description": "Spot-level QC and expression-based Leiden clusters over the H&E image.",
            "figures": figure_entries(collect_pngs(results_dir / "spatial_qc" / "figures"))
            + figure_entries(collect_pngs(results_dir / "spatial_cluster" / "figures")),
            "tables": [],
        },
        {
            "title": "Tangram deconvolution",
            "description": "Mapped cell-type proportions. Germinal-centre B and FDC scores concentrate in follicles; T-cell scores in the paracortex.",
            "figures": figure_entries(collect_pngs(results_dir / "tangram" / "figures"))
            + figure_entries(collect_pngs(results_dir / "tangram_figures")),
            "tables": [],
        },
        {
            "title": "Spatial statistics",
            "description": "Neighbourhood enrichment of spatial Leiden clusters, and Moran's I for highly variable genes. FDCSP marks follicle-associated structure; CCL21 marks the complementary T-zone/stroma.",
            "figures": figure_entries(collect_pngs(results_dir / "spatial_stats" / "figures"))
            + figure_entries(collect_pngs(results_dir / "stats_figures")),
            "tables": [],
        },
    ]

    moran_path = results_dir / "spatial_stats" / "moran.tsv"
    if not moran_path.is_file():
        moran_path = results_dir / "moran.tsv"
    moran_html = tsv_to_html(moran_path)
    if moran_html:
        sections[-1]["tables"].append(
            {"caption": "Top Moran's I genes (first rows of moran.tsv).", "html": moran_html}
        )

    agreement = results_dir / "sc_cluster" / "V1_Human_Lymph_Node.agreement.tsv"
    if not agreement.is_file():
        matches = list((results_dir / "sc_cluster").glob("*agreement.tsv"))
        agreement = matches[0] if matches else agreement
    agreement_html = tsv_to_html(agreement, max_rows=8)
    if agreement_html:
        sections[0]["tables"].append(
            {"caption": "Leiden vs curated Subset (excerpt).", "html": agreement_html}
        )

    n_pngs = sum(len(section["figures"]) for section in sections)
    sample = {
        "id": sample_id,
        "n_figure_files": n_pngs,
        "spatial_stats_tables": "moran.tsv + nhood_enrichment.tsv",
    }

    return {
        "sample": sample,
        "sections": sections,
        "versions": collect_versions(results_dir),
    }


def main() -> int:
    """Render report.html from published pipeline artifacts."""
    args = parse_args()
    if not args.results_dir.is_dir():
        raise FileNotFoundError(f"Results directory not found: {args.results_dir}")

    template_path = args.template
    if template_path is None:
        template_path = Path(__file__).resolve().parents[1] / "assets" / "report_template.html.j2"
    if not template_path.is_file():
        raise FileNotFoundError(f"Template not found: {template_path}")

    env = Environment(
        loader=FileSystemLoader(str(template_path.parent)),
        autoescape=select_autoescape(["html", "j2"]),
    )
    template = env.get_template(template_path.name)
    html = template.render(**build_context(args.results_dir, args.sample_id))

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(html, encoding="utf-8")
    print("Wrote", args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
