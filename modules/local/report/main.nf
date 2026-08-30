/*
 * MODULE: REPORT
 * PHASE:  5  (Days 10-11)
 *
 * I/O CONTRACT
 *
 *   INPUT:
 *     tuple val(meta), path(figures, stageAs: 'figures/ *')
 *     tuple val(meta), path(tables,  stageAs: 'tables/ *')
 *     path versions                      collected versions.yml from every upstream process
 *
 *   OUTPUT:
 *     path "report.html",  emit: report
 *     path "versions.yml", emit: versions
 *
 *   SCRIPT: bin/make_report.py   (template: assets/report_template.html.j2)
 *   LABELS: process_low, container_report
 *
 * DESIGN NOTE: figures are embedded as base64 data URIs, so report.html is a single file with no
 * relative-path dependencies. You can email it, attach it to an application, or open it from a USB
 * stick and it still renders. A report that breaks when moved is not a deliverable.
 *
 * The report must tell the whole story without the reader opening a single PNG: sample summary table
 * (n cells, n spots, n clusters), the QC figures, the annotated UMAP, the spatial cluster map, the
 * per-cell-type deconvolution maps, the neighborhood enrichment heatmap, the top spatially variable
 * genes, and the tool versions that produced all of it.
 *
 * TODO Phase 5: implement.
 */
