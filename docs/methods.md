# Methods

> Write this like a journal methods section: precise, past tense, every tool with its version and
> every non-default parameter with its value. Fill each section in **as you complete its phase**, not
> at the end — you will not remember which resolution you settled on or why by week three.
>
> This is one of the highest-value documents in the repo. A PI or hiring manager will actually read
> it, and it is the fastest way to demonstrate that you understand what your own pipeline did.

*Status: updated through Phase 5 (spatial statistics and HTML report).*

---

## Datasets

The spatial dataset was the 10x Genomics `V1_Human_Lymph_Node` Visium sample. It contained 4,035
spots and 36,601 genes, with tissue flags, array coordinates, spatial pixel coordinates, and a paired
H&E image. The single-cell reference was the integrated secondary lymphoid organ atlas from
Kleshchevnikov et al. (2022), downloaded from the public cell2location resource. It contained 73,260
cells and 10,237 genes, with curated cell-type labels in `Subset` (34 categories). Both datasets used
gene symbols, which makes their later gene matching straightforward. The data were downloaded during
Phase 1 and cached under the configured data directory.

## Single-cell RNA-seq quality control

QC metrics were recalculated from the reference expression matrix. Cells with fewer than 200 detected
genes and genes detected in fewer than 3 cells were filtered. The reference contained no mitochondrial
genes, so mitochondrial filtering was explicitly skipped rather than treated as zero evidence of
stress. Scrublet was attempted but was too slow on the full 73,260-cell reference; the successful
pipeline run therefore disabled it. The final reference retained all 73,260 cells and 10,237 genes,
which is consistent with this already-processed public reference and the conservative thresholds.

## Single-cell normalisation, clustering and annotation

The reference was library-size normalised to a target sum of 10,000 counts per cell and transformed
with `log1p`. The full normalised matrix was preserved in `adata.raw`; 2,000 highly variable genes
were selected using the Seurat method. PCA used 50 components, followed by a 15-neighbour graph and
Leiden clustering at resolution 1.0 using the `leidenalg` backend. UMAP used a fixed random seed of
0 for reproducibility. We computed clusters independently of the published labels, then compared
them with the 34-category curated `Subset` column using a Leiden-by-Subset crosstab. Canonical marker
genes were ranked with the Wilcoxon method and visualised in a dotplot.

## Spatial transcriptomics quality control and clustering

Visium spots were filtered using `in_tissue`, a minimum of 500 total counts, a minimum of 200 detected
genes, and a maximum mitochondrial fraction of 20% when mitochondrial genes were available. All
4,035 spots were flagged as on-tissue; 4,025 remained after quality filtering. A Visium spot is a
multi-cell measurement, not an individual cell, so doublet detection was not applied. Filtered spots
were normalised to 10,000 counts, log-transformed, reduced to 2,000 highly variable genes and 30 PCs,
then clustered using 15 expression neighbours at Leiden resolution 1.0. Squidpy also constructed a
spatial neighbour graph for later spatial analysis.

## Deconvolution

Annotated reference cells (`Subset`) were mapped onto QC-filtered Visium spots with Tangram
(`tangram-sc`) in PyTorch. Gene identifiers were compared as exact, case-sensitive symbols. The
intersection contained 10,141 genes; `tg.pp_adatas` retained 10,103 training genes (`gene_to_lowercase=False`).
The 38 excluded genes were present by name in both objects but had zero counts in the spatial matrix.
A minimum of 1,000 shared genes was required or the script exited.

A diagnostic run used cell-level mapping (`mode=cells`) for 20 epochs on CPU to confirm the API,
memory, and output shapes (73,260 × 4,025 mapping matrix; 4,025 × 34 proportions). The production
Nextflow process used cluster-level mapping (`mode=clusters`, `cluster_label=Subset`) for 200 epochs
on CPU with an RNA-count-based density prior, because cell-level mapping of the full atlas was too
slow for routine pipeline runs. Mapping weights were aggregated to cell-type scores with
`project_cell_annotations`. Each spot was normalised to a compositional vector (row sums of zero
were set to missing). Spatial maps were plotted with `scanpy.pl.spatial` on the hires H&E image for
`B_GC_DZ`, `B_GC_LZ`, `FDC`, `T_CD4+`, and `T_CD8+_cytotoxic`. The deconvolution process ran in the
`stpipe-tangram:0.1.0` image and joined the PREPROCESS single-cell and spatial outputs on sample
metadata.

## Spatial statistics

Neighbourhood enrichment and Moran’s I were computed on the deconvolved Visium AnnData object
(4,025 spots) with Squidpy, using the **spatial** neighbour graph already stored in
`obsp['spatial_connectivities']` from Phase 3. The expression neighbour graph was not used for
co-occurrence tests.

Neighbourhood enrichment (`squidpy.gr.nhood_enrichment`) treated `spatial_leiden` as the categorical
spot label (10 observed clusters, labels 0–9). Cluster labels were permuted on the fixed spatial graph
1,000 times in the intended production setting (20 permutations were used only as a smoke test).
Z-scores were written to `nhood_enrichment.tsv` and plotted as a heatmap. A large positive z-score
means two labels occupy neighbouring spots more often than the permutation null; a large negative
z-score means they occupy neighbouring spots less often than chance.

Moran’s I (`squidpy.gr.spatial_autocorr`, mode Moran) was restricted to the 2,000 highly variable
genes to keep permutation cost tractable. Production runs used 100 permutations per gene. Genes were
ranked by I descending (`moran.tsv`). The top genes included `IGHG1`, `CCL21`, `FDCSP`, and `IGHG2`.
Those four were plotted on the hires H&E image with `squidpy.pl.spatial_scatter`. Ligand–receptor
analysis (`squidpy.gr.ligrec`) was not run.

The spatial-statistics process used the `stpipe-scanpy:0.1.0` image and took the deconvolved `.h5ad`
from the Tangram process.

## Report

A Jinja2 template (`assets/report_template.html.j2`) plus `bin/make_report.py` assembled published
PNGs and selected TSVs into a single `report.html`. Figures were embedded as base64 data URIs so the
file has no relative image paths and no CDN. The Nextflow `REPORT` process ran in `stpipe-report:0.1.0`
(Jinja2 and pandas only), staging figure/TSV artifacts rather than the large `.h5ad` files. The
process was gated on `SPATIAL_STATS` so the HTML is produced after enrichment and Moran’s I.

## Reproducibility

QC, clustering, and spatial statistics used `stpipe-scanpy:0.1.0`; Tangram used `stpipe-tangram:0.1.0`;
the HTML report used `stpipe-report:0.1.0`. The full DAG through the report is reproducible with:

```bash
nextflow run main.nf -profile docker -resume \
  --sc_run_scrublet false \
  --tangram_mode clusters \
  --tangram_num_epochs 200 \
  --tangram_device cpu \
  --nhood_n_perms 1000 \
  --moran_n_perms 100
```

Lower `--nhood_n_perms` / `--moran_n_perms` (e.g. 20) were used for Docker smoke tests. Changing
Tangram `--tangram_num_epochs` or `--tangram_mode` invalidates the deconvolution cache under
`-resume`; changing only the permutation flags reruns spatial statistics on the same mapped object.

Nextflow stores task work directories separately from published results under `results/`. Runtime
provenance is written to `results/pipeline_info/` when those reports are enabled.

## Software versions

*This table will be generated from collected `versions.yml` files after the remaining modules emit them.*

| Tool | Version |
|---|---|
| Nextflow | |
| scanpy | |
| squidpy | |
| anndata | |
| Tangram | |
| PyTorch | |
| scrublet | |
| leidenalg | |
