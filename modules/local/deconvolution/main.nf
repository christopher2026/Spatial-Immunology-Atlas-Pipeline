/*
 * MODULE: DECONVOLUTION  (Tangram)
 * PHASE:  4  (Days 7-9)  -- the centerpiece; budget the most time here
 *
 * I/O CONTRACT
 *
 *   INPUT:
 *     tuple val(meta), path(sc_h5ad), path(sp_h5ad)
 *         annotated reference from SC_CLUSTER_ANNOTATE + clustered spots from SPATIAL_CLUSTER
 *
 *   OUTPUT:
 *     tuple val(meta), path("*.deconvolved.h5ad"),  emit: h5ad     (spots + per-type proportions
 *                                                                   in .obs / .obsm)
 *     tuple val(meta), path("*.proportions.tsv"),   emit: props    (spots x cell types)
 *     tuple val(meta), path("figures/ *.png"),      emit: figures  (per-type spatial maps for GC B
 *                                                                   cells, FDCs, CD4 T, CD8 T)
 *     tuple val(meta), path("*.shared_genes.tsv"),  emit: genes    (audit trail: which genes matched)
 *     path "versions.yml",                          emit: versions
 *
 *   SCRIPT: bin/run_tangram.py
 *   LABELS: process_high_memory, container_tangram
 *
 * WHAT TANGRAM DOES: it learns a soft assignment matrix M (cells x spots) by gradient descent,
 * maximising the similarity between each spot's measured expression and the expression predicted by
 * the mixture of cells assigned to it, over the genes shared by both datasets. It is optimisation
 * over a matching objective rather than a generative probabilistic model - which is exactly the
 * contrast with cell2location (Bayesian hierarchical, gives posterior uncertainty, slower). Have an
 * opinion on that tradeoff.
 *
 * THE LANDMINE: gene identifier conventions. The reference and the Visium data may use gene symbols
 * vs Ensembl IDs, or differ in case, or carry version suffixes (ENSG00000141510.14). If the shared
 * gene set silently collapses to a handful of genes, Tangram will still run happily and produce a
 * confident-looking map made of noise. Hence the `shared_genes.tsv` output and a hard assertion on
 * the shared gene count in the script. This is the most likely candidate for your "hardest bug"
 * interview story - write down what actually happens in docs/learning-log.md.
 *
 * TODO Phase 4: implement. Get it working standalone on REAL data before containerising.
 */
