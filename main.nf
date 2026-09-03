#!/usr/bin/env nextflow
/*
========================================================================================
    spatial-immune-atlas
========================================================================================
    Integrates single-cell RNA-seq and Visium spatial transcriptomics from human lymph
    node to spatially map immune cell types, identify tissue niches, and characterise
    cell-cell communication.

    Github: https://github.com/christopher2026/Spatial-Immune-Atlas-Pipeline
----------------------------------------------------------------------------------------

    THE DAG

        sc_reference (h5ad)              visium_dir
              |                              |
        [ SC_QC ]                      [ SPATIAL_QC ]
              |                              |
        [ SC_CLUSTER_ANNOTATE ]        [ SPATIAL_CLUSTER ]
              |                              |
              +-------------+----------------+
                            |
                    [ DECONVOLUTION ]   (Tangram: sc cell types -> spatial spots)
                            |
                    [ SPATIAL_STATS ]   (nhood enrichment, Moran's I, optional ligrec)
                            |
                       [ REPORT ]       (single standalone results.html)

    Build order (see PROJECT_PLAN.md): modules are wired in one at a time, each verified in
    isolation before the next is added. Do not uncomment a block until its module passes its
    own Definition of Done.
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
========================================================================================
    IMPORTS
========================================================================================
*/

// Phases 1-3
include { PREPROCESS }          from './subworkflows/local/preprocess.nf'
// Phase 4
// include { DECONVOLUTION }       from './modules/local/deconvolution/main.nf'
// Phase 5
// include { SPATIAL_STATS }       from './modules/local/spatial_stats/main.nf'
// include { REPORT }              from './modules/local/report/main.nf'

/*
========================================================================================
    HELP
========================================================================================
*/

def helpMessage() {
    log.info """
    ${workflow.manifest.name} v${workflow.manifest.version}

    Usage:
        nextflow run main.nf -profile docker
        nextflow run main.nf -profile test,docker      # fast subsampled smoke test

    Required inputs (defaults read from \$STPIPE_DATA, populated by bin/fetch_data.py):
        --sc_reference    Path to the scRNA-seq reference h5ad
        --visium_dir      Path to the Visium sample directory
        --sample_id       Sample identifier used in output filenames and the report

    Common options:
        --outdir          Output directory                        [${params.outdir}]
        --max_cpus        Cap on CPUs any single task may request [${params.max_cpus}]
        --max_memory      Cap on memory any single task may request [${params.max_memory}]

    Full parameter list with defaults and rationale: nextflow.config
    """.stripIndent()
}

/*
========================================================================================
    WORKFLOW
========================================================================================
*/

workflow {

    if (params.help) {
        helpMessage()
        return
    }

    log.info """
    ${'='*72}
    ${workflow.manifest.name} v${workflow.manifest.version}
    ${'='*72}
    sample_id     : ${params.sample_id}
    sc_reference  : ${params.sc_reference}
    visium_dir    : ${params.visium_dir}
    outdir        : ${params.outdir}
    profile       : ${workflow.profile}
    work dir      : ${workflow.workDir}
    ${'='*72}
    """.stripIndent()

    // ------------------------------------------------------------------------------------
    // Input channels.
    //
    // `checkIfExists: true` is not optional politeness - without it a typo'd path produces an
    // empty channel, every downstream process is silently skipped, and the pipeline reports
    // success having done nothing. That failure mode has cost people days.
    // ------------------------------------------------------------------------------------

    ch_sc_reference = Channel.fromPath(params.sc_reference, checkIfExists: true)
                             .map { file -> tuple([ id: params.sample_id ], file) }

    ch_visium       = Channel.fromPath(params.visium_dir, type: 'dir', checkIfExists: true)
                             .map { dir -> tuple([ id: params.sample_id ], dir) }

    // Every channel element is a `tuple(meta, file)` where meta is a map of sample metadata.
    // This is the nf-core convention: metadata travels alongside the data instead of being
    // encoded in filenames, which is what makes multi-sample support (stretch goal 4) a config
    // change rather than a rewrite.

    // ------------------------------------------------------------------------------------
    // PHASES 1-3 - preprocessing subworkflow
    // ------------------------------------------------------------------------------------
    PREPROCESS(ch_sc_reference, ch_visium)

    // ------------------------------------------------------------------------------------
    // PHASE 4 - deconvolution. First process that joins the two branches of the DAG.
    // ------------------------------------------------------------------------------------
    // ch_deconv_in = SC_CLUSTER_ANNOTATE.out.h5ad.join( SPATIAL_CLUSTER.out.h5ad )
    // DECONVOLUTION( ch_deconv_in )

    // ------------------------------------------------------------------------------------
    // PHASE 5 - spatial statistics and report
    // ------------------------------------------------------------------------------------
    // SPATIAL_STATS( DECONVOLUTION.out.h5ad )
    // REPORT( ... collected figures and tables ... )

    log.info "Phases 1-3 preprocessing are wired; deconvolution is the next stage."
}

