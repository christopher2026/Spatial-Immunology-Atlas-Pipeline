# Project Plan: Spatial Immune Atlas Pipeline

Tracker for the build described in [`PROJECT_HANDOFF.md`](PROJECT_HANDOFF.md).
Update the checkboxes as you go. Each phase has a **Definition of Done (DoD)** — do not
advance until it is met, because the DoD checks are what make this defensible in an interview.

- **Start:** Mon 31 Aug 2026
- **Target finish:** Fri 18 Sep 2026 (15 working days / 3 calendar weeks)
- **Pace assumed:** 4–6 focused hours per working day

**Legend:** `[ ]` todo · `[~]` in progress · `[x]` done · `[-]` skipped (note why)

---

## Progress at a glance

| Phase | Focus | Days | Dates | Status |
|---|---|---|---|---|
| 0 | Environment & repo scaffold | 0 | Sun 30 Aug | `[x]` |
| 1 | Data fetch + scRNA-seq QC + first Nextflow module | 1–2 | Mon 31 Aug – Tue 1 Sep | `[x]` |
| 2 | scRNA-seq clustering & annotation | 3–4 | Wed 2 – Thu 3 Sep | `[x]` |
| 3 | Spatial QC & spatial clustering | 5–6 | Fri 4 Sep, Mon 7 Sep | `[~]` |
| 4 | Deconvolution (Tangram) + full end-to-end DAG | 7–9 | Tue 8 – Thu 10 Sep | `[ ]` |
| 5 | Spatial statistics + HTML report | 10–11 | Fri 11 Sep, Mon 14 Sep | `[ ]` |
| 6 | Test dataset, CI, container polish | 12–13 | Tue 15 – Wed 16 Sep | `[ ]` |
| 7 | Documentation & interview prep | 14–15 | Thu 17 – Fri 18 Sep | `[ ]` |
| S | Stretch goals (only if core is done) | — | — | `[ ]` |

---

## Phase 0 — Environment & repo scaffold  ·  Sun 30 Aug

**Goal:** a repo skeleton and a working Linux toolchain, so that no phase later is blocked on installs.

### Your tasks (need admin rights / reboot)
- [x] Install WSL2 + Ubuntu (`wsl --install -d Ubuntu` in an **Administrator** PowerShell, then reboot)
- [x] Install Docker Desktop and enable **Settings → Resources → WSL Integration → Ubuntu**
- [x] Run `bash setup/bootstrap_wsl.sh` inside Ubuntu (installs Java, Nextflow, micromamba)
- [x] Verify: `nextflow -version`, `docker run hello-world`, `micromamba --version` all succeed inside WSL

Full step-by-step instructions: [`setup/WSL_SETUP.md`](setup/WSL_SETUP.md)

### Agent tasks
- [x] Repo directory skeleton per handoff §4
- [x] `nextflow.config` + `conf/{base,docker,test}.config`
- [x] Three `docker/*.Dockerfile` with pinned versions
- [x] `env/environment.yml` (analysis env) and `env/environment-tangram.yml`
- [x] `.gitignore`, `README.md` skeleton, `docs/` skeletons
- [x] `bin/fetch_data.py` — downloads both public datasets
- [x] Git repo initialised with first commit

**DoD:** `nextflow -version` and `docker run hello-world` both work inside WSL; `git log` shows the scaffold commit.

**What you should be able to explain:** why the pipeline source lives on the Windows filesystem but the
Nextflow `work/` directory and `data/` live on the native Linux filesystem (symlink semantics + IO speed).

---

## Phase 1 — Data fetch + scRNA-seq QC  ·  Days 1–2 (Mon 31 Aug – Tue 1 Sep)

**Goal:** one Nextflow process, fully working in Docker. This is the pattern every later module reuses,
so getting it genuinely clean now saves days later.

- [x] Run `python bin/fetch_data.py --all` — reference `sc.h5ad` (73,260 cells) + Visium lymph node
- [x] Explore the reference interactively: `.obs` columns, existing cell type labels, gene ID convention, is there a mito gene set?
- [x] Write `bin/sc_qc.py` with `argparse` (`--input`, `--outdir`, QC thresholds as flags)
  - [x] `sc.pp.calculate_qc_metrics`
  - [x] cell/gene filtering by min counts + min genes
  - [x] mitochondrial fraction — **check whether mito genes are even present**, document rather than assume
  - [ ] doublet detection (`sc.pp.scrublet`) — deferred: too slow on the full 73k-cell reference
  - [x] ≥2 QC figures (violin of counts/genes/pct_mt, scatter counts-vs-genes)
  - [x] printed text summary of cells/genes removed
- [x] Run it standalone outside Nextflow and eyeball the figures
- [x] Build the `scanpy_squidpy` Docker image
- [x] Write `modules/local/sc_qc/main.nf` — **write the I/O contract comment before the process body**
- [x] Write a minimal `main.nf` that runs only `sc_qc`
- [x] `nextflow run main.nf -profile docker` succeeds

**DoD:** `python bin/sc_qc.py --input ... --outdir ...` produces a filtered h5ad + ≥2 QC figures + a text
summary; the same script runs through Nextflow in Docker and lands outputs in `results/`.

**Teaching checkpoints (be able to answer these before moving on):**
- What does a Nextflow *channel* do, and why is it a dataflow model rather than a dependency graph?
- Why does the Python live in `bin/` and get called by CLI, instead of being embedded in the `script:` block?
- Which QC thresholds did you pick and *why* — what does a very high count / very low gene count spot mean biologically?
- What is a doublet, and why does it matter for downstream deconvolution specifically?

---

## Phase 2 — scRNA-seq clustering & annotation  ·  Days 3–4 (Wed 2 – Thu 3 Sep)

**Goal:** reproduce the reference's cell types with your own pipeline, then *validate* against the curated labels.

- [x] `bin/sc_cluster_annotate.py`: normalize_total → log1p → HVGs → PCA → neighbors → Leiden → UMAP
- [x] Marker gene ranking (`sc.tl.rank_genes_groups`) + dotplot for canonical lineage markers
- [x] Agreement table: your Leiden clusters × provided cell type labels (crosstab / normalised confusion matrix)
- [x] Outputs: annotated h5ad, UMAP by cell type, UMAP by Leiden cluster, marker dotplot, agreement TSV
- [x] `modules/local/sc_cluster_annotate/main.nf`, wired downstream of `sc_qc`
- [ ] Note mismatches in `docs/learning-log.md` — disagreement is a real result, not a failure

**DoD:** UMAP visually separates B / T / myeloid / stromal lineages; agreement table exists and mismatches
are documented with a hypothesis for each.

**Teaching checkpoints:**
- Why log-normalize before HVG selection, and what does `log1p` actually fix?
- Why PCA before the neighbor graph — what happens if you skip it?
- What does Leiden `resolution` control, and how would you defend the value you chose?
- Why is UMAP distance between distant clusters *not* meaningful?

---

## Phase 3 — Spatial QC & clustering  ·  Days 5–6 (Fri 4 Sep, Mon 7 Sep)

**Goal:** the spatial half of the input, plus your first tissue-overlay figure.

- [x] `bin/spatial_qc.py` — spot-level QC (**spots ≠ cells**; state this explicitly in the report)
- [x] `bin/spatial_cluster.py` — normalize, HVG, PCA, neighbors, Leiden on spots
- [x] `squidpy.pl.spatial_scatter` cluster overlay on the H&E image
- [x] `modules/local/spatial_qc/main.nf` and `modules/local/spatial_cluster/main.nf`
- [x] `subworkflows/local/preprocess.nf` wiring `sc_qc` + `spatial_qc`
- [x] Compare your cluster map against known lymph node anatomy — do you see follicle-shaped regions?

**DoD:** spatial cluster plot overlaid on H&E, with clusters that correspond to plausible anatomical regions.

**Teaching checkpoints:**
- Why is a Visium spot (~55 µm) a *mixture* of cells, and why does that make deconvolution necessary?
- What QC metric is meaningful for spots but meaningless for cells (and vice versa)?
- What does a subworkflow buy you over just putting both processes in `main.nf`?

---

## Phase 4 — Deconvolution with Tangram  ·  Days 7–9 (Tue 8 – Thu 10 Sep)

**Goal:** the ML centerpiece, plus the first full end-to-end DAG run. Budget the most time here.

- [ ] Get Tangram working **standalone on real (not subsampled) data first** — confirm the biology before containerising
- [ ] Subset to shared genes between reference and spatial; handle gene-ID convention mismatch (symbols vs Ensembl)
- [ ] `tangram.mapping_utils.map_cells_to_space` (CPU is fine at this scale)
- [ ] Project cell type probabilities onto spots → per-spot cell type proportion matrix (TSV)
- [ ] Spatial plots for key types: germinal centre B cells, follicular dendritic cells, CD4 T, CD8 T
- [ ] Build `tangram.Dockerfile`; write `modules/local/deconvolution/main.nf`
- [ ] Wire the **full** `main.nf` DAG: sc_qc → sc_cluster_annotate ↘ deconvolution ↙ spatial_qc → spatial_cluster
- [ ] Validate visually against published cell2location Fig. 4 on the same dataset; write the comparison into `docs/results.md`

**DoD:** full pipeline runs start-to-finish on real data; B-cell types land in follicle-shaped regions and
T-cell types in the surrounding paracortex.

**Teaching checkpoints:**
- What objective is Tangram actually optimizing, and in what sense is it "deep-learning-adjacent"?
- Why does the shared-gene set matter so much, and what fails silently if gene IDs mismatch?
- How would you know if your mapping were garbage? (Name at least two independent checks.)
- Optimization-based (Tangram) vs Bayesian hierarchical (cell2location) — when would you pick each?

---

## Phase 5 — Spatial statistics + report  ·  Days 10–11 (Fri 11 Sep, Mon 14 Sep)

**Goal:** turn the map into quantitative claims, then package everything into one shareable file.

- [ ] `bin/spatial_stats.py`
  - [ ] `squidpy.gr.spatial_neighbors` + `squidpy.gr.nhood_enrichment` → heatmap
  - [ ] `squidpy.gr.spatial_autocorr` (Moran's I) → ranked spatially variable gene table
  - [ ] *(stretch)* `squidpy.gr.ligrec` ligand–receptor
- [ ] `bin/make_report.py` + `assets/report_template.html.j2` — figures + auto-generated summary table
- [ ] `modules/local/spatial_stats/main.nf`, `modules/local/report/main.nf`, `docker/report.Dockerfile`

**DoD:** `results/report.html` opens standalone and tells the whole pipeline story without opening any PNG.

**Teaching checkpoints:**
- What is Moran's I measuring, and what's the null hypothesis?
- Neighborhood enrichment gives a z-score against a permutation null — permutation of *what*?
- Name one spatially variable gene you found and explain why it makes biological sense here.

---

## Phase 6 — Test dataset, CI, containers  ·  Days 12–13 (Tue 15 – Wed 16 Sep)

**Goal:** the green badge and the one-command verification that make the repo credible at a glance.

- [ ] `bin/make_test_data.py` — subsample to ~500 cells / ~200 spots, save to `assets/test_data/`
- [ ] `conf/test.config` pointing at it with reduced Tangram epochs and low resource limits
- [ ] `.github/workflows/ci.yml` running `nextflow run main.nf -profile test,docker` on push/PR
- [ ] Push all three images to GHCR; reference them by digest-pinned tag in `nextflow.config`
- [ ] *(if time)* `nf-test` unit tests: assert QC output has fewer cells than input and expected `.obs` columns
- [ ] CI badge in README

**DoD:** green CI badge; `nextflow run main.nf -profile test,docker` completes in under ~5 minutes from a clean clone.

**Teaching checkpoints:**
- Why one container per tool rather than one big container for everything?
- What does version pinning actually protect you from, and where does pinning still leave you exposed?

---

## Phase 7 — Documentation & interview prep  ·  Days 14–15 (Thu 17 – Fri 18 Sep)

- [ ] `README.md`: elevator pitch, Mermaid DAG, quickstart, key results figure inline, badges
- [ ] `docs/methods.md`: paper-style methods with tools, versions, key parameters
- [ ] `docs/results.md`: 3–5 sentence biological interpretation — **the single most differentiating deliverable**
- [ ] **Clean-clone test**: clone to a fresh directory and run `-profile test,docker` from scratch
- [ ] Write out answers to handoff §11 talking points in your own words
- [ ] Rehearse the 60-second and 5-minute pipeline walkthroughs out loud

**DoD:** a stranger can clone the repo, run one command, and see results; you can talk for 5 minutes about
the biology without looking at notes.

---

## Stretch goals (ranked by interview payoff per hour)

- [ ] cell2location as a second deconvolution method + comparison figure
- [ ] Ligand–receptor analysis (`squidpy.gr.ligrec`) for cell–cell communication hypotheses
- [ ] `nf-core pipelines lint` compliance + parameter JSON schema
- [ ] Multi-sample support via a samplesheet CSV (nf-core convention)
- [ ] Supervised cell type classifier on PCA latent space with cross-validated accuracy
- [ ] Interactive report (Quarto — already installed on this machine — or Streamlit)

---

## Risk log

| Risk | Likelihood | Mitigation |
|---|---|---|
| WSL2/Docker install friction (admin, virtualisation disabled in BIOS) | Medium | Do it Phase 0, before anything depends on it |
| Gene-ID mismatch between reference (symbols) and Visium (Ensembl) | **High** | Explicitly normalise IDs in Phase 4; this is likely your "hardest bug" interview story |
| Tangram memory blowup on 73k cells × full gene set | Medium | Subset to shared HVGs; reduce cells per type; CPU mode |
| Reference h5ad download is ~2 GB and slow | Medium | Fetch on Day 1, cache in `data/`, never re-download |
| `sc.datasets.visium_sge` renamed in scanpy ≥1.10 | High | `fetch_data.py` already handles both APIs |
| Scope creep into stretch goals before core works | Medium | Phases 0–7 are locked; stretch only after Phase 7 DoD |

---

## Learning log

Keep a running log in [`docs/learning-log.md`](docs/learning-log.md) of every non-obvious bug and decision.
Handoff §11.6 asks "what was the hardest bug?" — you want a real, specific answer with a debugging
narrative, and you will not remember it three weeks from now unless you write it down the day it happens.
