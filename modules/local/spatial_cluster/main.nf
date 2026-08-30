/*
 * MODULE: SPATIAL_CLUSTER
 * PHASE:  3  (Days 5-6)
 *
 * I/O CONTRACT
 *
 *   INPUT:
 *     tuple val(meta), path(h5ad)         filtered Visium h5ad from SPATIAL_QC
 *
 *   OUTPUT:
 *     tuple val(meta), path("*.sp_clustered.h5ad"), emit: h5ad
 *     tuple val(meta), path("figures/ *.png"),      emit: figures  (clusters overlaid on H&E via
 *                                                                   squidpy.pl.spatial_scatter,
 *                                                                   UMAP of spots)
 *     path "versions.yml",                          emit: versions
 *
 *   SCRIPT: bin/spatial_cluster.py
 *   LABELS: process_medium, container_scanpy
 *
 * DESIGN NOTE: clustering here is on expression only, so spatial coherence in the resulting map is
 * an *emergent* property, not something we imposed. If your clusters form contiguous, follicle-shaped
 * regions without ever having seen the coordinates, that is real biological signal and is a strong
 * thing to point at. (Spatially-aware clustering methods exist and would be a fair follow-up
 * question - know that you chose the simpler baseline deliberately.)
 *
 * TODO Phase 3: implement.
 */
