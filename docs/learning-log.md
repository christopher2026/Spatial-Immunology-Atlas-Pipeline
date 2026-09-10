# Learning log

Running log of decisions, surprises and bugs. Two purposes:

1. **Interview ammunition.** Handoff §11 asks "what was the hardest bug?" and "why did you choose X
  over Y?". Those need specific, honest answers with a debugging narrative. You will not remember the
   details three weeks from now unless you write them down the day they happen.
2. **Understanding.** Writing down why something broke is the step where you actually learn it.

Keep entries short. One entry per non-obvious thing. Symptom → cause → fix → what you learned.

---

## Template

```
### YYYY-MM-DD — one-line title
**Phase:** N
**Symptom:** what you observed, verbatim error if there was one
**Cause:** what was actually wrong
**Fix:** what you changed
**Learned:** the generalisable lesson — this is the part you say out loud in an interview
```

---



## Decisions



### 2026-08-30 — Repo source on Windows, Nextflow work dir on native Linux ext4

**Phase:** 0
**Decision:** pipeline source stays at `C:\Users\Chris\Documents\stpipe` (visible in WSL at
`/mnt/c/...`), but `NXF_WORK` and the data directory point at `~/nxf-work` and `~/nxf-data` on the
WSL ext4 filesystem.
**Why:** Nextflow stages each task's inputs as symlinks into a per-task working directory. The
`/mnt/c` drvfs bridge supports symlinks poorly and is roughly an order of magnitude slower for the
many-small-files access pattern Nextflow generates. Keeping source on Windows keeps the editor
workflow simple; keeping work and data on Linux keeps execution correct and fast.
**Learned:** Nextflow's execution model is symlink-based, which is exactly why it has no native
Windows port — it is not an oversight.

### 2026-08-30 — Three containers, one per tool, rather than one image

**Phase:** 0
**Decision:** separate `scanpy_squidpy`, `tangram` and `report` images.
**Why:** torch is the largest and most version-fragile dependency in the project. Isolating it means
a torch upgrade cannot break the QC modules, the QC image stays small enough to pull quickly in CI,
and the report image needs no scientific stack at all. This is standard nf-core practice.
**Learned:** container boundaries are a dependency-isolation decision, not just a packaging one.

### 2026-08-30 — Python 3.11 rather than the 3.13 already installed

**Phase:** 0
**Decision:** pin the analysis environment to Python 3.11.
**Why:** the scanpy / squidpy / torch / tangram stack has the widest binary wheel coverage on 3.11.
A three-week deconvolution project is not the place to spend a day compiling igraph from source.
**Learned:** in scientific Python, the newest interpreter is usually the wrong choice; pick the one
the ecosystem has settled on.

---



## Bugs and surprises



### 2026-08-31 — Package versions belonged to different Python channels

**Phase:** 0–1  
**Symptom:** micromamba could not solve `anndata==0.13.3.post0` and `scanpy==1.12.4`.  
**Cause:** those releases required Python 3.12 and were not available from the conda-forge pins we
started with, while the environment used Python 3.11 for Tangram/PyTorch compatibility.  
**Fix:** use Python 3.11-compatible AnnData 0.12.19, Scanpy 1.11.5, and Squidpy 1.8.2, installed
through pip where appropriate.  
**Learned:** a version pin is not enough; the interpreter version and package channel are part of the
reproducibility contract.

### 2026-08-31 — Scrublet was impractical on the full reference

**Phase:** 1  
**Symptom:** QC printed its metrics and then ran overnight without completing.  
**Cause:** Scrublet is computationally expensive on 73,260 cells, and the initial script had no progress
message around that step.  
**Fix:** disabled Scrublet for the successful QC run and added progress reporting.  
**Learned:** a pipeline needs observable stages and a small test path; an algorithm can be valid but
still be inappropriate for the available compute budget.

### 2026-09-01 — Docker image lacked a C++ compiler

**Phase:** 1  
**Symptom:** the Scanpy image failed while building the `annoy` wheel with `g++: No such file or directory`.  
**Cause:** PyPI had no matching wheel, so pip attempted a source build inside the image.  
**Fix:** added `build-essential` to the Scanpy/Squidpy and Tangram images.  
**Learned:** containerizing Python does not remove native build dependencies; compiled scientific
packages still need a compiler when wheels are unavailable.

### 2026-09-02 — Nextflow task could not find the Python script

**Phase:** 2  
**Symptom:** `python: can't open file 'sc_qc.py'`.  
**Cause:** Nextflow staged the input file in the task directory, but the Python script was not there;
`python sc_qc.py` does not search the project `bin/` directory.  
**Fix:** called the script explicitly through `${projectDir}/bin/sc_qc.py`.  
**Learned:** Nextflow task directories are isolated; input files and project scripts must be referenced
through their declared/staged paths.

### 2026-09-03 — Subworkflow imports were inside a comment

**Phase:** 3  
**Symptom:** Nextflow reported that `SC_QC`, `SPATIAL_QC`, and the other processes were undefined.  
**Cause:** the `include` statements had been placed before the closing `*/` of the documentation block.  
**Fix:** moved all imports outside the comment.  
**Learned:** a visually correct block of code is still invisible to the parser if it is inside a comment;
when every symbol is undefined, check whether the import was actually parsed.

### 2026-09-03 — Squidpy plotting API differed from Scanpy

**Phase:** 3  
**Symptom:** `spatial_scatter` forwarded `img_key` and then `show` to Matplotlib, causing
`PatchCollection` errors.  
**Cause:** this Squidpy version uses `img_res_key` and `return_ax`, not those Scanpy-style plotting
arguments.  
**Fix:** switched to `sq.read.visium`, `img_res_key="hires"`, and `return_ax=True`.  
**Learned:** similarly named plotting APIs are not interchangeable; inspect the installed function
signature when a keyword unexpectedly reaches Matplotlib.

### 2026-09-04 — Tangram dropped 38 named genes that were all-zero on Visium

**Phase:** 4  
**Symptom:** exact symbol intersection was 10,141 genes; `tg.pp_adatas` kept 10,103 training genes.  
**Cause:** the 38 genes existed in both `.var_names` but had zero counts on this slide (including
Y-chromosome genes and rare receptors). Tangram subsets the objects to training genes, so looking
those names up *after* `pp_adatas` raised `KeyError`. Counting nonzeros on all 10,141 genes first
took about an hour because of repeated sparse column slices on `/mnt/c`.  
**Fix:** treat the drop as expected; do not force the 38 genes back; document them; skip the slow
audit in production.  
**Learned:** name overlap is not expression overlap. A silent empty intersection is the dangerous
case; here the intersection was large and the 38 genes were spatially blank.

### 2026-09-06 — Cell-level Tangram was too slow for 1,000 epochs on CPU

**Phase:** 4  
**Decision:** diagnostic `mode=cells`, 20 epochs (~35 min); production Nextflow `mode=clusters`,
200 epochs, `cluster_label=Subset`, CPU.  
**Why:** 73,260 × 4,025 cell-level mapping on CPU would be many hours at 1,000 epochs. Cluster mode
maps 34 type-average profiles. Follicle vs paracortex still matched.  
**Learned:** Tangram’s library default of 1,000 cell-level epochs is not a biological requirement.
Reproducibility is image + params + `-resume`, not copying the tutorial epoch count.

### 2026-09-06 — Nextflow `join` must use PREPROCESS outputs

**Phase:** 4  
**Symptom:** `main.nf` cannot see `SC_CLUSTER_ANNOTATE.out`.  
**Cause:** those processes live inside the PREPROCESS subworkflow.  
**Fix:** `PREPROCESS.out.sc_h5ad.join(PREPROCESS.out.sp_h5ad)` then `DECONVOLUTION`.  
**Learned:** `join` pairs on `meta.id`, not filenames. Subworkflow emits are the only public contract.

### 2026-09-07 — A trailing backslash glued `printf` onto the Python command

**Phase:** 4–5  
**Symptom:** `versions.yml` was written in the work directory but missing from `results/spatial_stats/`,
or Python received extra arguments.  
**Cause:** a `\` after the last argparse flag continues the line, so `printf` became part of
`run_tangram.py` / `spatial_stats.py`. Separately, Nextflow only publishes files listed in `output:`.  
**Fix:** drop the final `\`; emit `path "versions.yml"`.  
**Learned:** in a Nextflow `script:` block, `\` is bash line continuation. The last flag must not have
one if the next line is a new command.

### 2026-09-07 — Neighbourhood enrichment must not use the expression graph

**Phase:** 5  
**Decision:** `nhood_enrichment` and Moran’s I use `obsp['spatial_connectivities']`, not
`expression_neighbors`. Labels for enrichment are `spatial_leiden` (10 clusters), not `Subset`.  
**Why:** expression neighbours connect transcriptionally similar spots that may be far apart on the
slide. Spatial co-occurrence is a question about **tissue adjacency**. `Subset` is a cell-level
reference column, not a spot label.  
**Learned:** two graphs can live on the same AnnData; using the wrong one still “runs” and produces
a heatmap that does not mean proximity.

### 2026-09-08 — Scanpy `pl.spatial` vs Squidpy `spatial_scatter`

**Phase:** 5  
**Symptom:** `FutureWarning: Use squidpy.pl.spatial_scatter instead.`  
**Cause:** Scanpy’s Visium plot is deprecated in this stack. Squidpy wants `img_res_key` / `return_ax`.  
**Fix:** Moran plots in `spatial_stats.py` use `sq.pl.spatial_scatter`. Tangram plots stayed on
`sc.pl.spatial` because `stpipe-tangram` does not install Squidpy.  
**Learned:** the same figure can be produced by two APIs; the container contents decide which one is
legal in that process.

### 2026-09-09 — Report process must not stage the full `results/` tree

**Phase:** 5  
**Decision:** `REPORT` waits on `SPATIAL_STATS`, stages PNGs and Moran/nhood TSVs plus the Jinja2
template, and embeds images as base64. It does not copy deconvolved `.h5ad` files.  
**Why:** staging `params.outdir` would pull multi-gigabyte AnnData into the report task.  
**Learned:** the last process should depend on upstream *channels* for correctness, and copy only
the artifacts the HTML actually needs.

### 2026-09-09 — `-resume` hashes include interpolated parameters

**Phase:** 5  
**Symptom:** changing `--tangram_num_epochs` made Tangram run again; changing only `--nhood_n_perms`
did not.  
**Cause:** Nextflow’s cache key includes the process script after `${params...}` substitution.  
**Fix:** pass the same Tangram CLI flags when you only want to rerun stats or the report.  
**Learned:** `-resume` is not “skip everything that ever succeeded”; it is “skip tasks whose hash
is unchanged.”