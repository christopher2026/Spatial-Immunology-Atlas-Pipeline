# Results

> **This is the single most differentiating document in the project.** Plenty of student portfolios
> contain working code. Very few contain evidence that the author understood what the output meant.
>
> Write for a reader who knows immunology but not your code. Every claim should point at a specific
> figure or table, and every figure should answer a question someone would actually ask.

*Status: skeleton. Filled in across Phases 3–5 and finalised in Phase 7.*

---

## Summary

<!-- 3-5 sentences. What did the pipeline show? Aim for something like: "Deconvolution recovered the
     canonical lymph node architecture: germinal-centre B cells and follicular dendritic cells
     localised to discrete round regions consistent with follicles, while CD4 and CD8 T cell
     populations occupied the surrounding paracortex..." Then one sentence on what the spatial
     statistics added quantitatively, and one on the limitations. -->

## Single-cell reference

<!-- Phase 2. The annotated UMAP, the marker dotplot, and the agreement table. Key question to
     answer: did your independent clustering recover the curated cell types, and where did it not?
     Name the specific populations that merged or split, and give a reason (transcriptional
     similarity, low cell numbers, resolution choice). -->

## Spatial structure of the tissue

<!-- Phase 3. The Leiden cluster map over H&E. Key point: clustering used expression only and never
     saw the coordinates, so spatial coherence in the result is emergent signal, not an artifact of
     the method. Describe which clusters look like follicles. -->

## Spatially mapped cell types

<!-- Phase 4. The per-cell-type deconvolution maps. This is the centrepiece figure. Address:
       - Do B-lineage types concentrate in follicle-shaped regions?
       - Do T-cell types occupy the surrounding paracortex?
       - Where do follicular dendritic cells sit relative to germinal-centre B cells?
       - Anything unexpected, and your best hypothesis for it. -->

### Validation strategy

<!-- Phase 4. Be explicit that this is validation, not just a result. Two independent checks:
       1. Anatomical plausibility against known lymph node architecture.
       2. Visual comparison against the published cell2location result on this same dataset
          (Kleshchevnikov et al. 2022, Fig. 4).
     State what agreement you saw and, honestly, what did not agree. "It matched perfectly" is less
     believable than a specific discrepancy with an explanation. -->

## Spatial statistics

<!-- Phase 5. Neighborhood enrichment heatmap: which populations co-occur, which exclude each other,
     and does the B-follicle / T-zone mutual depletion appear as expected? Moran's I: name the top
     spatially variable genes and say why each makes biological sense in this tissue. -->

## Limitations

<!-- Phase 7. Be specific and honest. Candidates:
       - One Visium slide, one donor: no biological replication, so nothing here is a population claim.
       - Visium spots are multi-cell; proportions are inferred, not measured.
       - Tangram gives a point estimate with no uncertainty quantification (contrast: cell2location's
         posterior).
       - Reference and spatial data come from different donors and tissues (spleen/tonsil included),
         so composition priors may not transfer cleanly.
     Interviewers weight honest limitations heavily. Overclaiming is the fastest way to lose a
     technical reader. -->

## What I would do with more time or compute

<!-- Phase 7. 2-3 concrete answers, drawn from the stretch goals. Prefer specific over aspirational:
     "run cell2location on the same data to compare a Bayesian posterior against Tangram's point
     estimate" beats "explore more methods". -->
