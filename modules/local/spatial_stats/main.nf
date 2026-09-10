/*
 * MODULE: SPATIAL_STATS
 * PHASE:  5  (Days 10-11)
 *
 * I/O CONTRACT
 *
 *   INPUT:
 *     tuple val(meta), path(h5ad)         deconvolved spots from DECONVOLUTION
 *
 *   OUTPUT:
 *     tuple val(meta), path("figures/ *.png"),   emit: figures  (nhood enrichment heatmap,
 *                                                                spatial plots of top SVGs)
 *     tuple val(meta), path("*.nhood.tsv"),      emit: nhood    (z-scores)
 *     tuple val(meta), path("*.moran.tsv"),      emit: moran    (ranked spatially variable genes)
 *     tuple val(meta), path("*.ligrec.tsv"),     emit: ligrec, optional: true
 *     path "versions.yml",                       emit: versions
 *
 *   SCRIPT: bin/spatial_stats.py
 *   LABELS: process_medium, container_scanpy
 *
 * WHAT THE TWO CORE STATISTICS MEAN:
 *
 *   Neighborhood enrichment (squidpy.gr.nhood_enrichment) builds a spatial neighbor graph over
 *   spots, then counts how often each pair of cluster labels appears as neighbors, and compares that
 *   count against a null built by permuting the cluster labels while holding the graph fixed. A high
 *   z-score means two populations co-occur more than tissue geometry alone would produce; a strongly
 *   negative one means they exclude each other. In lymph node you expect B-follicle and T-zone
 *   clusters to be mutually depleted, with a shared boundary - a built-in sanity check.
 *
 *   Moran's I (squidpy.gr.spatial_autocorr) measures spatial autocorrelation of a gene's expression:
 *   do nearby spots have similar values? Null is spatial randomness. Genes with high Moran's I are
 *   spatially structured, and here they should recover follicle/T-zone marker programmes without ever
 *   being told where the anatomy is.
 *
 * TODO Phase 5: implement.
 */
process SPATIAL_STATS {
    tag "$meta.id"
    label 'process_medium'
    label 'container_scanpy'
    publishDir "${params.outdir}/spatial_stats", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(h5ad)

    output:
    tuple val(meta), path("figures/*.png"), emit: figures
    tuple val(meta), path("*nhood_enrichment.tsv"), emit: nhood
    tuple val(meta), path("moran.tsv"), emit: moran
    path "versions.yml", emit: versions

    script:
    """
    python "${projectDir}/bin/spatial_stats.py" \
        --spatial-input "$h5ad" \
        --outdir "." \
        --cluster-key "spatial_leiden" \
        --nhood-n-perms ${params.nhood_n_perms} \
        --moran-n-perms ${params.moran_n_perms}

    printf '%s\n' \
        'spatial_stats:' \
        "  scanpy: \$(python -c 'import scanpy; print(scanpy.__version__)')" \
        "  squidpy: \$(python -c 'import squidpy; print(squidpy.__version__)')" \
        > versions.yml
    """
}