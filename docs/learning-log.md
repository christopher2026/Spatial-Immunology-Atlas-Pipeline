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

<!-- Add entries here as they happen. Do not skip the small ones; the "hardest bug" question is
     usually best answered by a bug that was hard to *find*, not hard to fix. -->
