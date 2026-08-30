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
 *
 * TODO Phase 3: implement, once SC_QC, SC_CLUSTER_ANNOTATE, SPATIAL_QC and SPATIAL_CLUSTER all pass
 * their own Definitions of Done individually.
 */
