/*
 * SUBWORKFLOW: PREPROCESS
 * PHASE: 3  (Days 5-6)
 *
 * I/O CONTRACT
 *
 *   TAKE:
 *     ch_sc_reference   tuple val(meta), path(h5ad)
 *     ch_visium         tuple val(meta), path(visium_dir)
 *
 *   EMIT:
 *     sc_h5ad           tuple val(meta), path(h5ad)     filtered + annotated reference
 *     sp_h5ad           tuple val(meta), path(h5ad)     filtered + clustered spots
 *     figures           tuple val(meta), path(png)      all QC and clustering figures, merged
 *     versions          path(versions.yml)              collected
 *
 * WHY A SUBWORKFLOW AT ALL: the two input branches (single-cell and spatial) are independent until
 * deconvolution joins them. Grouping them here means main.nf reads as three steps - preprocess,
 * deconvolve, report - instead of seven, and the branches become independently testable and
 * independently reusable. This is composition, and it is the main practical argument for DSL2 over
 * DSL1: workflows are values you can import, not just top-level scripts.
 * TODO Phase 3: implement once each process passes its own Definition of Done.
 */

include { SC_QC }               from '../../modules/local/sc_qc/main.nf'
include { SC_CLUSTER_ANNOTATE } from '../../modules/local/sc_cluster_annotate/main.nf'
include { SPATIAL_QC }          from '../../modules/local/spatial_qc/main.nf'
include { SPATIAL_CLUSTER }     from '../../modules/local/spatial_cluster/main.nf'

workflow PREPROCESS {
    take:
    ch_sc_reference
    ch_visium

    main:
    SC_QC(ch_sc_reference)
    SC_CLUSTER_ANNOTATE(SC_QC.out.h5ad)

    SPATIAL_QC(ch_visium)
    SPATIAL_CLUSTER(SPATIAL_QC.out.h5ad)

    emit:
    sc_h5ad = SC_CLUSTER_ANNOTATE.out.h5ad
    sp_h5ad = SPATIAL_CLUSTER.out.h5ad
    figures = SC_QC.out.figures
        .mix(SC_CLUSTER_ANNOTATE.out.figures)
        .mix(SPATIAL_QC.out.figures)
        .mix(SPATIAL_CLUSTER.out.figures)
}
