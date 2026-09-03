/*
 * MODULE: SPATIAL_QC
 *
 * I/O CONTRACT
 *
 *   INPUT:
 *     tuple val(meta), path(visium_dir)   10x Visium sample dir (matrix + spatial/ with H&E image)
 *
 *   OUTPUT:
 *     tuple val(meta), path("*.sp_filtered.h5ad"), emit: h5ad
 *     tuple val(meta), path("figures/ *.png"),     emit: figures  (violin of per-spot metrics,
 *                                                                  counts overlaid on tissue)
 *     tuple val(meta), path("*.sp_qc_summary.tsv"), emit: summary
 *     path "versions.yml",                          emit: versions
 *
 *   SCRIPT: bin/spatial_qc.py
 *   LABELS: process_medium, container_scanpy
 *
 * CRITICAL DISTINCTION: a Visium spot is ~55 um across and contains roughly 1-10 cells. Spots are
 * NOT cells. Consequences:
 *   - "doublet detection" is meaningless here; every spot is a mixture by construction
 *   - a high per-spot count may mean dense tissue, not a technical artifact
 *   - this mixing is precisely why deconvolution (Phase 4) is necessary at all
 * Conflating spots and cells is the most common beginner error in spatial transcriptomics. Say this
 * out loud in the report and in interviews - it signals rigor cheaply.
 */

 process SPATIAL_QC {
    tag "$meta.id"
    label 'process_medium'
    label 'container_scanpy'
    publishDir "${params.outdir}/spatial_qc", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(visium_dir)

    output:
    tuple val(meta), path("*.filtered.h5ad"), emit: h5ad
    tuple val(meta), path("figures/*.png"), emit: figures

    script:
    def output_h5ad = "${meta.id}.filtered.h5ad"

    """
    python "${projectDir}/bin/spatial_qc.py" \
        --input "$visium_dir" \
        --output "$output_h5ad" \
        --figures-dir figures \
        --min-counts ${params.sp_min_counts} \
        --min-genes ${params.sp_min_genes} \
        --max-pct-mt ${params.sp_max_pct_mt}
    """
 }
