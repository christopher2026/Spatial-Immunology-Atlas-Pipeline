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
