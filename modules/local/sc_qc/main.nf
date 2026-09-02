/*
 * MODULE: SC_QC
 * PHASE:  1  (Days 1-2)
 *
 * I/O CONTRACT  -- write this before the process body, always.
 *
 *   INPUT:
 *     tuple val(meta), path(h5ad)      raw scRNA-seq reference (73,260 cells x 10,237 genes)
 *
 *   OUTPUT:
 *     tuple val(meta), path("*.filtered.h5ad"),  emit: h5ad
 *     tuple val(meta), path("figures/ *.png"),   emit: figures    (>= 2: violin + scatter)
 *     tuple val(meta), path("*.qc_summary.tsv"), emit: summary    (cells/genes in vs out)
 *     path "versions.yml",                       emit: versions
 *
 *   SCRIPT: bin/sc_qc.py
 *   LABELS: process_high_memory, container_scanpy
 *
 * WHY A CONTRACT FIRST: a Nextflow process is a pure function over channels. If you cannot state
 * its inputs and outputs in two lines, the step is doing too much and should be split. Writing the
 * contract first also forces you to decide what the *artifact* is - and every artifact here has to
 * be something you can put on screen in an interview and explain.
 *
 * WHY `versions.yml`: nf-core convention. Every process emits the version of the tool it ran, and
 * they are collected into the report. This is how you answer "what exactly produced this figure?"
 * six months later.
 *
 * TODO Phase 1: implement. See PROJECT_PLAN.md.
 */

process SC_QC {
    tag "$meta.id"
    label 'process_high_memory'
    label 'container_scanpy'
    publishDir "${params.outdir}/sc_qc", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(h5ad)

    output:
    tuple val(meta), path("*.filtered.h5ad"), emit: h5ad
    tuple val(meta), path("figures/*.png"), emit: figures
    tuple val(meta), path("*.qc_summary.tsv"), emit: summary
    path "versions.yml", emit: versions

    /*
     * Keep the Nextflow block thin: it stages files, passes parameters, and publishes artifacts.
     * The scientific logic stays in bin/sc_qc.py so it can be run and debugged without Nextflow.
     */
    script:
    def output_h5ad = "${meta.id}.filtered.h5ad"
    def run_scrublet = params.sc_run_scrublet.toString().toBoolean()
    """
    python "${projectDir}/bin/sc_qc.py" \
        --input "$h5ad" \
        --output "$output_h5ad" \
        --min-genes ${params.sc_min_genes} \
        --min-cells ${params.sc_min_cells} \
        --max-pct-mt ${params.sc_max_pct_mt} \
        ${run_scrublet ? '--run-scrublet' : '--no-run-scrublet'}
    printf '%s\n' \
        'sc_qc:' \
        "  scanpy: \$(python -c 'import scanpy; print(scanpy.__version__)')" \
        "  anndata: \$(python -c 'import anndata; print(anndata.__version__)')" \
        "  scrublet: \$(python -c 'import scrublet; print(scrublet.__version__)' 2>/dev/null || echo unavailable)" \
        > versions.yml
    """
}
