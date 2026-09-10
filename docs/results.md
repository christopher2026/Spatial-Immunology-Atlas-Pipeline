# Results

> **This is the single most differentiating document in the project.** Plenty of student portfolios
> contain working code. Very few contain evidence that the author understood what the output meant.
>
> Write for a reader who knows immunology but not your code. Every claim should point at a specific
> figure or table, and every figure should answer a question someone would actually ask.

*Status: updated through Phase 4 (Tangram). Spatial statistics are not yet complete.*

---

## Summary

The pipeline ran end to end through Tangram in Nextflow and Docker. The scRNA-seq reference contained
73,260 cells and 34 curated `Subset` labels; Visium QC retained 4,025 of 4,035 spots. Expression-based
clustering produced 34 scRNA-seq Leiden clusters and 11 spot clusters. Tangram then mapped those
reference cell types onto spots. Germinal-centre B cells and follicular dendritic cells concentrated
in follicle-shaped regions on the H&E image, while CD4 and cytotoxic CD8 T cells were higher in the
surrounding paracortex, consistent with lymph-node anatomy and with the published cell2location maps
on this same 10x sample.

## Single-cell reference

The independently generated analysis produced 34 scRNA-seq Leiden clusters. UMAPs and the marker
dotplot showed clear expression structure across the reference. Myeloid and stromal populations are
present but relatively rare, and some have weaker marker evidence; their annotations should therefore
be interpreted with more caution. The agreement table compares expression-based clusters with the
curated `Subset` labels rather than assuming the two definitions are identical.

## Spatial structure of the tissue

The filtered Visium data contained 4,025 spots. Expression-based clustering produced 11 spot clusters,
which were visualised both in UMAP space and over the H&E image. The overlay lets us assess whether
expression-defined groups are spatially coherent. These are spot-level expression domains, not
confirmed cell types or anatomical regions; that interpretation is reserved for deconvolution.

## Spatially mapped cell types

Tangram assigned reference cells (or `Subset` averages, in cluster mode) to Visium spots using 10,103
training genes. Of 10,141 genes that matched by exact symbol, 38 were dropped because they had zero
counts on this slide; name matching was not the failure mode. Spot-by-type proportions were obtained
by aggregating mapping weights by `Subset` and normalising each spot to sum to 1.

Five spatial maps were examined: dark-zone and light-zone germinal-centre B cells (`B_GC_DZ`,
`B_GC_LZ`), follicular dendritic cells (`FDC`), helper T cells (`T_CD4+`), and cytotoxic T cells
(`T_CD8+_cytotoxic`). GC B and FDC scores formed round hotspots that coincide with follicles on the
H&E image. CD4 and CD8 scores were higher in the tissue between those circles (paracortex) and lower
on follicle cores. That opposing pattern is the main biological result: B-lineage programmes in
follicles, T-lineage programmes in the T zone.

These values are relative mixture weights, not counted cells. Early cell-level maps (20 epochs) had
narrow colour scales, so localisation was interpreted from *where* scores peaked rather than from
absolute percentages. The production Nextflow run used cluster-level mapping so the job could finish
on CPU; both views support the same follicle versus paracortex contrast.

### Validation strategy

Preprocessing was checked with the Leiden-versus-`Subset` agreement table, marker dotplots, and
whether Visium clusters were spatially coherent on H&E. Deconvolution was checked in two further
ways: (1) whether mapped types occupy known lymph-node compartments, and (2) whether that layout
agrees in broad terms with cell2location Fig. 4 on `V1_Human_Lymph_Node`. A map that is pretty but
puts B cells in the paracortex would be rejected. CD4 inside follicles would not automatically mean
failure: GC T follicular helpers are expected there, whereas CD8 inside follicles would be more
suspicious.

## Spatial statistics

*IGHG1/2 are IgG constant region transcripts (plasma cell/medully programmes) and are spatially clumped rather than a follicle ring. CCL21 (*T-zone chemokine from stromal/FRC cells) *is in the surrounding paracortex and is complementary to the FDCSP gene, which is Follicular dendretic cell products and should sit in follicle/GC patches.* 

## Limitations

The reference is already processed and contains no mitochondrial genes, so mitochondrial filtering
could not be independently applied to its cells. Scrublet was disabled for the full reference because
it was computationally impractical in the available environment. Visium spots are multi-cell
mixtures, so spatial clusters cannot be interpreted as individual cell types. Tangram proportions are
model-based and have no posterior intervals. Cluster-level mapping averages cells within each
`Subset`, so within-type heterogeneity is not used. The 20-epoch cell-level run was a smoke test, not
a fully trained map. Cell2location comparison uses the same tissue class and a related reference, so
it is a sanity check rather than a fully independent assay. Neighbourhood enrichment and ligand–
receptor analysis have not been run.

## What I would do with more time or compute

