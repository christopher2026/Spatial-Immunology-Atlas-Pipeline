/*
 * MODULE: REPORT
 * PHASE:  5  (Days 10-11)
 *
 * Stages PNGs and TSVs from upstream processes and renders a standalone HTML file.
 * Large h5ad objects are not copied into this task.
 */

process REPORT {
    tag "$meta.id"
    label 'process_low'
    label 'container_report'
    publishDir "${params.outdir}", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(moran)
    path nhood
    path stats_pngs,     stageAs: 'stats_figures/*'
    path tangram_pngs,   stageAs: 'tangram_figures/*'
    path preprocess_pngs, stageAs: 'preprocess_figures/*'
    path template

    output:
    path "report.html", emit: report
    path "versions.yml", emit: versions

    script:
    """
    mkdir -p spatial_stats
    cp -L "${moran}" spatial_stats/moran.tsv
    cp -L "${nhood}" spatial_stats/nhood_enrichment.tsv

    python "${projectDir}/bin/make_report.py" \
        --results-dir . \
        --output report.html \
        --template "${template}" \
        --sample-id "${meta.id}"

    printf '%s\n' \
        'report:' \
        "  jinja2: \$(python -c 'import jinja2; print(jinja2.__version__)')" \
        "  pandas: \$(python -c 'import pandas; print(pandas.__version__)')" \
        > versions.yml
    """
}
