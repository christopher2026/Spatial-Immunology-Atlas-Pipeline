/*
 * MODULE: SC_CLUSTER_ANNOTATE
 * PHASE:  2  (Days 3-4)
 *
 * I/O CONTRACT
 *
 *   INPUT:
 *     tuple val(meta), path(h5ad)         filtered scRNA-seq reference from SC_QC
 *
 *   OUTPUT:
 *     tuple val(meta), path("*.annotated.h5ad"),  emit: h5ad
 *     tuple val(meta), path("figures/ *.png"),    emit: figures   (UMAP by celltype, UMAP by
 *                                                                  leiden, marker dotplot)
 *     tuple val(meta), path("*.agreement.tsv"),   emit: agreement (leiden x curated label crosstab)
 *     tuple val(meta), path("*.markers.tsv"),     emit: markers
 *     path "versions.yml",                        emit: versions
 *
 *   SCRIPT: bin/sc_cluster_annotate.py
 *   LABELS: process_high_memory, container_scanpy
 *
 * DESIGN NOTE: this reference already ships curated cell type labels. We do NOT re-annotate from
 * scratch. We reproduce clusters with our own pipeline and then *validate* against the curated
 * labels via the agreement table. That is both more honest and more defensible: it demonstrates the
 * annotation workflow while making the result checkable. Disagreement between your Leiden clusters
 * and the curated labels is a finding to explain, not a bug to hide.
 *
 * TODO Phase 2: implement.
 */

process SC_CLUSTER_ANNOTATE {
    tag "$meta.id"
    label 'process_high_memory'
    label 'container_scanpy'
    publishDir "${params.outdir}/sc_cluster", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(h5ad)

    output:
    tuple val(meta), path("*.annotated.h5ad"), emit: h5ad
    tuple val(meta), path("figures/*.png"), emit: figures
    tuple val(meta), path("*.agreement.tsv"), emit: agreement
    path "versions.yml", emit: versions

    /*
     * The module handles staging, parameters and artifact collection. The scientific workflow
     * remains in bin/sc_cluster_annotate.py so it can run standalone.
     */
    script:
    def output_h5ad = "${meta.id}.annotated.h5ad"
    def agreement_path = "${meta.id}.agreement.tsv"
    """
    python "${projectDir}/bin/sc_cluster_annotate.py" \
        --input "$h5ad" \
        --output "$output_h5ad" \
        --figures-dir figures \
        --agreement "$agreement_path" \
        --celltype-key "${params.sc_celltype_key}" \
        --n-hvg ${params.sc_n_hvg} \
        --n-pcs ${params.sc_n_pcs} \
        --n-neighbors ${params.sc_n_neighbors} \
        --resolution ${params.sc_leiden_resolution}

    printf '%s\n' \
        'sc_cluster_annotate:' \
        "  scanpy: \$(python -c 'import scanpy; print(scanpy.__version__)')" \
        "  anndata: \$(python -c 'import anndata; print(anndata.__version__)')" \
        "  leidenalg: \$(python -c 'import leidenalg; print(leidenalg.__version__)' 2>/dev/null || echo installed)" \
        > versions.yml
    """
}
