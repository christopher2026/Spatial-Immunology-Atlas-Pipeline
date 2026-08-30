# Methods

> Write this like a journal methods section: precise, past tense, every tool with its version and
> every non-default parameter with its value. Fill each section in **as you complete its phase**, not
> at the end — you will not remember which resolution you settled on or why by week three.
>
> This is one of the highest-value documents in the repo. A PI or hiring manager will actually read
> it, and it is the fastest way to demonstrate that you understand what your own pipeline did.

*Status: skeleton. Sections are filled in per phase — see [`../PROJECT_PLAN.md`](../PROJECT_PLAN.md).*

---

## Datasets

<!-- Phase 1. Include: accessions/URLs, n cells, n genes, n spots, cell type label provenance,
     gene identifier convention for each dataset, and the date you downloaded them. -->

## Single-cell RNA-seq quality control

<!-- Phase 1. Include: filtering thresholds and the rationale for each; whether mitochondrial genes
     were present in the reference at all; doublet detection method, expected doublet rate, and the
     number of cells called as doublets; final n cells and n genes retained. -->

## Single-cell normalisation, clustering and annotation

<!-- Phase 2. Include: normalisation target sum, log transform, HVG method and count, n PCs and how
     you chose it, n neighbors, Leiden resolution and why, UMAP parameters. For annotation: state
     clearly that the reference carries curated labels, that clusters were computed independently,
     and how agreement was quantified. -->

## Spatial transcriptomics quality control and clustering

<!-- Phase 3. Include: per-spot thresholds; explicitly note that spots are multi-cell mixtures and
     that no doublet detection was applied, with the reason; clustering parameters. -->

## Deconvolution

<!-- Phase 4. Include: Tangram version, mode (cells vs clusters), number of shared genes after ID
     harmonisation and how mismatched identifiers were resolved, density prior, number of epochs,
     device, and how per-spot cell type proportions were derived from the mapping matrix. -->

## Spatial statistics

<!-- Phase 5. Include: spatial neighbor graph construction (Delaunay vs grid, n_neighs); permutation
     counts for neighborhood enrichment and Moran's I; multiple-testing correction. -->

## Reproducibility

<!-- Phase 6. Include: container images and tags with digests, Nextflow version, how to reproduce
     (`nextflow run main.nf -profile test,docker`), and where run provenance is recorded
     (results/pipeline_info/). -->

## Software versions

<!-- Phase 6. Generate this table from the collected versions.yml files rather than typing it by
     hand — hand-maintained version tables drift and an interviewer may spot it. -->

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
