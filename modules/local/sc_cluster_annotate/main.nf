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
