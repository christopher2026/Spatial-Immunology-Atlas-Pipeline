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
