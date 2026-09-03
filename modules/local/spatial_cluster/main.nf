/*
 * MODULE: SPATIAL_CLUSTER
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
 */

process SPATIAL_CLUSTER {
    tag "$meta.id"
    label 'process_medium'
    label 'container_scanpy'
    publishDir "${params.outdir}/spatial_cluster", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(h5ad)

    output:
    tuple val(meta), path("*.sp_clustered.h5ad"), emit: h5ad
    tuple val(meta), path("figures/*.png"), emit: figures

    script:
    def output_h5ad = "${meta.id}.sp_clustered.h5ad"

    """
    python "${projectDir}/bin/spatial_cluster.py" \
        --input "$h5ad" \
        --output "$output_h5ad" \
        --figures-dir figures \
        --n-hvg ${params.sp_n_hvg} \
        --n-pcs ${params.sp_n_pcs} \
        --n-neighbors ${params.sp_n_neighbors} \
        --resolution ${params.sp_leiden_resolution}
    """
}