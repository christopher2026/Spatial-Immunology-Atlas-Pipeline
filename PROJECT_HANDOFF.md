# Project Handoff: Spatial Immune Atlas Pipeline
**A Nextflow-orchestrated, Dockerized single-cell + spatial transcriptomics pipeline for human lymph node immune profiling**

> This document is written for an AI coding agent (Claude, Copilot, etc.) working alongside the project owner. It defines scope, architecture, a day-by-day build plan, and acceptance criteria for each stage. Follow it sequentially. Do not skip the "Definition of Done" checks — they are what makes this project defensible in an interview.

---

## 0. Context for the agent

The user is a computational biology internship candidate with working Python/numpy/pandas/matplotlib skills, some scanpy/squidpy exposure, light ML background, and comfort with AnnData. They know Nextflow only lightly and are building this specifically to be interview-ready in **~2–3 focused weeks**. Prioritize:

1. **Correctness and a working end-to-end pipeline over exhaustive scope.** A lean pipeline that runs cleanly beats a sprawling one that's half-broken.
2. **Explaining biology, not just code.** Every module should produce an artifact (figure/table) the user can talk about in an interview — what it shows and why it's biologically meaningful.
3. **Teaching as you go.** When introducing a new tool (DSL2 syntax, a deconvolution algorithm, a container pattern), briefly explain *why* it's built that way before generating code, so the user can defend design choices in an interview, not just recite that "the agent did it."
4. **Nextflow/nf-core idioms**, since the user explicitly chose Nextflow for its industry relevance over Snakemake.

---

## 1. Elevator pitch (for the README and interviews)

> "I built a reproducible Nextflow pipeline that integrates single-cell RNA-seq and spatial transcriptomics data from human lymph node to spatially map immune cell types, identify tissue niches (T-cell zones, B-cell follicles, germinal centers), and characterize cell-cell communication — containerized end-to-end with Docker and tested in CI."

This sentence should be true by the end of week 2, and demonstrably true (with figures) by the end of week 3.

---

## 2. Skills this project is designed to build (map these explicitly onto your resume/interview prep)

| Skill area | What in this project demonstrates it |
|---|---|
| Single-cell analysis | QC, doublet detection, normalization, HVG selection, Leiden clustering, marker-based cell type annotation |
| Spatial transcriptomics | Visium QC, spatial clustering, neighborhood enrichment, spatial statistics (Moran's I) |
| Data integration / ML | Reference-based cell type deconvolution onto spatial spots using Tangram (deep-learning-based mapping); optional cell2location (Bayesian) comparison |
| Pipeline engineering | Nextflow DSL2 modules, channels, subworkflows, parameter schemas, resource profiles |
| Software engineering | Docker containers per tool, reproducible environments, version pinning |
| Testing / CI | nf-test module tests, GitHub Actions running the pipeline on a subsampled test dataset |
| Reproducibility & docs | README with DAG diagram, methods write-up, parameterized config, `-profile test` for one-command verification |
| Scientific communication | A results report translating deconvolution + spatial stats into a biological narrative about lymph node architecture |

---

## 3. Scientific background (read this before building anything)

**Tissue**: human lymph node — a well-studied, spatially organized immune organ with clear anatomical zones: **B-cell follicles/germinal centers**, **T-cell zones (paracortex)**, and **medullary/stromal regions**. This structure gives you a strong "ground truth" to sanity-check your spatial results against — if your pipeline is right, the deconvolved B cells should land in follicles and T cells in the paracortex.

**Datasets** (both public, both directly loadable via Python — no manual account/EULA friction):

- **Spatial (Visium)**: 10x Genomics `V1_Human_Lymph_Node`, loadable via:
  ```python
  import scanpy as sc
  adata_vis = sc.datasets.visium_sge(sample_id="V1_Human_Lymph_Node")
  ```
  (squidpy has an equivalent loader: `squidpy.datasets.visium("V1_Human_Lymph_Node")`.)

- **scRNA-seq reference**: the integrated secondary lymphoid organ atlas from Kleshchevnikov et al. (cell2location paper, *Nat Biotechnol* 2022) — 73,260 cells, 34 annotated immune/stromal cell types, combining lymph node/spleen/tonsil datasets specifically to capture germinal-center biology that single lymph node datasets miss:
  ```python
  adata_ref = sc.read(
      "sc.h5ad",
      backup_url="https://cell2location.cog.sanger.ac.uk/paper/integrated_lymphoid_organ_scrna/"
                 "RegressionNBV4Torch_57covariates_73260cells_10237genes/sc.h5ad",
  )
  ```
  This is the same reference used in the canonical cell2location/Tangram tutorials, so your results are benchmarkable against published figures — useful for sanity-checking and for confidently discussing expected results in interviews.

**Why this pairing over something else**: it's real, peer-reviewed, widely used in tutorials (so you have places to unstick yourself), has clean cell type labels already curated (you don't have to solve annotation from scratch, though you should still touch it), and has strong spatial structure that makes your visualizations interpretable at a glance.

---

## 4. Repository structure

```
spatial-immune-atlas/
├── main.nf                      # entry workflow
├── nextflow.config               # profiles: standard, docker, test
├── modules/
│   ├── local/
│   │   ├── sc_qc/                # scRNA-seq QC + filtering
│   │   ├── sc_cluster_annotate/  # clustering + cell type annotation
│   │   ├── spatial_qc/           # Visium QC + filtering
│   │   ├── spatial_cluster/      # spatial Leiden clustering
│   │   ├── deconvolution/        # Tangram mapping
│   │   ├── spatial_stats/        # neighborhood enrichment, Moran's I, ligand-receptor
│   │   └── report/               # HTML report generation
├── subworkflows/
│   └── local/
│       └── preprocess.nf         # wires sc_qc + spatial_qc
├── bin/                          # python scripts called by processes
│   ├── sc_qc.py
│   ├── sc_cluster_annotate.py
│   ├── spatial_qc.py
│   ├── spatial_cluster.py
│   ├── run_tangram.py
│   ├── spatial_stats.py
│   └── make_report.py
├── conf/
│   ├── test.config               # tiny subsampled data, low resources
│   └── docker.config
├── docker/
│   ├── scanpy_squidpy.Dockerfile
│   ├── tangram.Dockerfile
│   └── report.Dockerfile
├── assets/
│   └── report_template.html.j2
├── tests/
│   └── main.nf.test              # nf-test
├── .github/workflows/ci.yml
├── data/                         # gitignored; populated by a fetch script
├── results/                      # gitignored; pipeline output
├── docs/
│   └── methods.md                # your written biological interpretation
└── README.md
```

Keep Python analysis logic in `bin/*.py` scripts callable from the command line with argparse — this is the nf-core convention and it means every step is independently testable/runnable outside Nextflow too, which matters for debugging and for showing interviewers you can run steps standalone.

---

## 5. Pipeline design (the DAG)

```
                 ┌────────────┐
                 │  sc_qc     │  (raw reference h5ad -> filtered h5ad)
                 └─────┬──────┘
                        │
                 ┌─────▼───────────────┐
                 │ sc_cluster_annotate │  (Leiden clusters, marker genes,
                 └─────┬───────────────┘   cell type labels -> annotated h5ad)
                        │
       ┌────────────┐   │
       │ spatial_qc │   │
       └─────┬──────┘   │
              │           │
       ┌─────▼──────┐    │
       │spatial_    │    │
       │cluster     │    │
       └─────┬──────┘    │
              │           │
              └─────┬─────┘
                     │
              ┌─────▼──────────┐
              │ deconvolution   │  (Tangram: map sc cell types -> spatial spots)
              │  (Tangram)      │
              └─────┬───────────┘
                     │
              ┌─────▼──────────┐
              │ spatial_stats   │  (neighborhood enrichment, Moran's I,
              │                 │   optional ligand-receptor)
              └─────┬───────────┘
                     │
              ┌─────▼──────────┐
              │ report          │  (HTML report: figures + summary tables)
              └─────────────────┘
```

Each Nextflow process should have a **clearly typed input/output channel contract** — write this as a comment at the top of every module before writing the process body, e.g.:

```
// INPUT:  path to raw scRNA-seq h5ad
// OUTPUT: path to filtered h5ad, path to QC report figures (png)
```

This habit (I/O contracts before code) is standard in production bioinformatics pipelines and is worth explicitly mentioning in interviews as your process, not just something you did once.

---

## 6. Module-by-module spec

### 6.1 `sc_qc`
- Input: raw reference `sc.h5ad`
- Steps: calculate QC metrics (`sc.pp.calculate_qc_metrics`), filter cells/genes by min counts/genes, mitochondrial fraction (may be sparse/absent in some public references — check and note this rather than assuming), doublet detection via `scrublet` or `sc.pp.scrublet`.
- Output: filtered h5ad + QC violin/scatter plots.
- **Definition of done**: script runs standalone via `python bin/sc_qc.py --input ... --output ...`, produces a filtered h5ad and at least 2 QC figures, and prints a short text summary of cells/genes removed.

### 6.2 `sc_cluster_annotate`
- Steps: normalize (`sc.pp.normalize_total` + `log1p`), HVGs, PCA, neighbors, Leiden clustering, UMAP. Annotate: since this reference *already has curated cell type labels* (use them as ground truth), your job here is to (a) reproduce clusters via your own pipeline and (b) validate them against the provided labels via a marker-gene dot plot or confusion-matrix-style crosstab. This is more defensible in an interview than inventing annotation from scratch, while still proving you understand the annotation workflow.
- Output: annotated h5ad, UMAP colored by cell type, marker gene dotplot, clustering-vs-reference-label agreement table.
- **Definition of done**: UMAP visually separates major lineages (B cells, T cells, myeloid, stromal); agreement table exists even if not perfect (document mismatches — this is a legitimate and interesting result, not a failure).

### 6.3 `spatial_qc` and `spatial_cluster`
- Same QC pattern applied to Visium spots (spot-level, not cell-level — note this distinction explicitly in your report, since conflating spots and cells is a common beginner mistake and calling it out shows rigor).
- Spatial Leiden clustering + `squidpy.pl.spatial_scatter` visualization of clusters over the tissue image.
- **Definition of done**: spatial cluster plot overlaid on the H&E image; clusters visually correspond to plausible anatomical regions.

### 6.4 `deconvolution` (Tangram)
- This is the ML/integration centerpiece. Tangram trains a mapping (via optimization of a matching objective, deep-learning-adjacent) that aligns single-cell profiles to spatial spots using shared genes.
- Steps: subset to shared genes between reference and spatial data, run `tangram.mapping_utils.map_cells_to_space`, project cell type probabilities onto spatial coordinates, threshold/normalize to get per-spot cell type composition.
- Output: per-spot cell type proportion matrix, spatial plots of key cell types (e.g., germinal center B cells, follicular dendritic cells, CD4/CD8 T cells) overlaid on tissue.
- **Definition of done**: spatial plots show B-cell-associated types concentrated in follicle-shaped regions and T-cell types in surrounding paracortex — compare visually against the published cell2location Fig. 4 result to sanity-check (cite this comparison explicitly in your report as a validation strategy, not just a result).
- **Stretch**: repeat with `cell2location` (heavier install, GPU-friendly, Bayesian rather than optimization-based) and compare the two methods' spatial maps — this comparison is a great interview talking point ("I evaluated two deconvolution paradigms and discussed their tradeoffs").

### 6.5 `spatial_stats`
- Neighborhood enrichment (`squidpy.gr.nhood_enrichment`) on the deconvolved/clustered spots — do certain cell types spatially co-occur or avoid each other?
- Spatial autocorrelation (`squidpy.gr.spatial_autocorr`, Moran's I) to find spatially variable genes.
- **Stretch**: ligand-receptor analysis via `squidpy.gr.ligrec` (CellPhoneDB-style) to hypothesize cell-cell communication within niches — strong biological storytelling material.
- **Definition of done**: at least the neighborhood enrichment heatmap and a ranked list of spatially variable genes.

### 6.6 `report`
- A single Python script (Jinja2 templated HTML, or a Quarto doc if you want to stretch) that assembles all figures + a short auto-generated summary table (sample name, n cells, n spots, n clusters, top spatially variable genes) into one shareable HTML report per run.
- **Definition of done**: `results/report.html` opens standalone and tells the whole pipeline's story without needing to open individual PNGs.

---

## 7. Containers

Don't over-engineer this — 2–3 containers is enough and is realistic nf-core practice:

1. `scanpy_squidpy.Dockerfile` — scanpy, squidpy, scrublet, leidenalg, python-igraph. Used by QC/clustering/spatial_stats modules.
2. `tangram.Dockerfile` — tangram-sc + torch (CPU is fine at this data scale) + scanpy. Used by the deconvolution module.
3. `report.Dockerfile` — lightweight: jinja2, pandas (or quarto if you go that route).

Pin exact versions in each Dockerfile (`scanpy==1.10.x`, etc.) — reproducibility is the entire point of doing this, and an interviewer may ask "how do you know this will still run in a year?"

In `nextflow.config`, reference containers per-process with `container` directives and use the `docker` profile to enable `docker.enabled = true`. Build and push to a free registry (GitHub Container Registry, `ghcr.io`) so CI and anyone cloning the repo can pull rather than rebuild.

---

## 8. Testing & CI

- Create a **tiny test dataset**: subsample the reference to ~500 cells and the Visium data to ~200 spots (`sc.pp.subsample` / manual `.obs` slicing), save under `assets/test_data/` or fetch via a small script — this keeps CI runs under a few minutes.
- `conf/test.config` should set `params` to point at this subsampled data and reduce any resource-heavy parameters (fewer Tangram epochs, etc.).
- Add `.github/workflows/ci.yml` that runs `nextflow run main.nf -profile test,docker` on push/PR. A green CI badge in your README is one of the highest-signal, lowest-effort things you can add for recruiters skimming GitHub.
- If time allows, add a couple of `nf-test` unit tests for individual modules (e.g., assert the QC module output h5ad has fewer cells than input and all expected `.obs` columns present).

---

## 9. Documentation deliverables (do not skip — this is what gets read in 90 seconds by a recruiter)

- **README.md**: elevator pitch, pipeline DAG image (Mermaid diagram renders natively on GitHub — use it), quickstart (`nextflow run main.nf -profile test,docker`), key results figure inline, tech stack badges.
- **docs/methods.md**: a short, precise methods section written like a paper's methods (tools + versions + key parameters) — this is exactly the kind of writing a PI or hiring manager will actually read.
- **A results narrative** (can live in the README or a `docs/results.md`): 3–5 sentences interpreting what the spatial map + neighborhood enrichment actually show biologically. This is the single most differentiating piece of the whole project — most student portfolios show code but not scientific reasoning.

---

## 10. Day-by-day plan (focused 2–3 weeks / ~15 working days)

**Days 1–2 — Setup & scRNA-seq QC**
- Install Nextflow, Docker, set up repo skeleton above.
- Write `bin/sc_qc.py`, test it standalone on the reference h5ad.
- Write and run the first Nextflow module (`sc_qc`) locally with `-profile docker` on a subsampled file. Get this one process fully working before touching anything else — it teaches you the DSL2/Docker pattern you'll reuse for every other module.

**Days 3–4 — scRNA-seq clustering & annotation**
- `bin/sc_cluster_annotate.py` + module. Validate against provided reference labels.
- Checkpoint: annotated UMAP + agreement table exist.

**Days 5–6 — Spatial QC & clustering**
- `bin/spatial_qc.py`, `bin/spatial_cluster.py` + modules.
- Checkpoint: spatial cluster plot over tissue image.

**Days 7–9 — Deconvolution (Tangram)**
- This is the hardest and most important stage — budget the most days here.
- Get Tangram running standalone first (outside Nextflow) on the real (not subsampled) data to confirm it produces sensible biology before wrapping it in a container/process.
- Build the container, wire the module, wire the full pipeline DAG end-to-end for the first time.
- Checkpoint: full `main.nf` runs start to finish on real data; deconvolved cell type maps look anatomically plausible.

**Days 10–11 — Spatial statistics & report**
- `spatial_stats` module (neighborhood enrichment + Moran's I at minimum).
- `report` module assembling everything.
- Checkpoint: `results/report.html` tells the full story in one file.

**Days 12–13 — Testing, CI, containers polish**
- Build subsampled test dataset, `conf/test.config`, GitHub Actions CI.
- Pin all container versions, push images to GHCR.
- Checkpoint: green CI badge, `-profile test,docker` runs in under ~5 minutes.

**Days 14–15 — Documentation & interview prep**
- Write README, methods.md, results narrative, Mermaid DAG diagram.
- Do a full clean-clone-and-run test (delete your local repo, re-clone, run — this catches "works on my machine" bugs before an interviewer finds them).
- Prepare 3–4 talking points (Section 11).

If you finish early, use remaining time on the stretch goals in Section 12 rather than compressing the schedule — a clean core pipeline beats a rushed feature-complete one.

---

## 11. Interview talking points to prepare in advance

Be ready to answer, in your own words, not memorized:

1. **"Walk me through your pipeline."** — Practice a 60-second version and a 5-minute version.
2. **"Why Nextflow over Snakemake?"** — Have a real opinion: DSL2 modularity, container-per-process isolation, industry adoption (nf-core, pharma pipelines), channel-based dataflow vs. Snakemake's rule/wildcard model.
3. **"Why Tangram over cell2location (or vice versa if you did both)?"** — Optimization-based mapping vs. Bayesian hierarchical model; speed/compute tradeoffs; when you'd choose one over the other in practice.
4. **"How do you know your deconvolution results are correct?"** — Anatomical plausibility check against known lymph node structure + comparison to the published cell2location result on the same data.
5. **"What would you do differently with more time/compute?"** — Have 2–3 honest answers ready from Section 12.
6. **"What was the hardest bug?"** — Pick a real one (there will be one, likely gene-ID mismatches between reference and spatial data, or a Docker/dependency conflict) and narrate the debugging process, not just the fix.

---

## 12. Stretch goals (only after the core pipeline in Sections 5–9 is fully working and documented)

Ranked roughly by effort-to-payoff for internship marketability:

1. **cell2location as a second deconvolution method** + a comparison figure/section — high payoff, moderate effort.
2. **Ligand-receptor analysis** (`squidpy.gr.ligrec`) for cell-cell communication hypotheses within niches — high biological storytelling payoff.
3. **nf-core schema/lint compliance** (`nf-core pipelines lint`) — signals you know the ecosystem's conventions, not just Nextflow syntax.
4. **Multi-sample support** (parameterize for >1 Visium sample via a samplesheet CSV, nf-core convention) — demonstrates you designed for scale, not just one dataset.
5. **A small supervised ML add-on**: train a cell type classifier (e.g., logistic regression or small MLP on PCA/scVI latent space) on the reference and evaluate cross-validated accuracy — a compact, legitimate way to show ML fundamentals beyond just using someone else's deep learning tool.
6. **Interactive report** (Quarto or a small Streamlit app) instead of static HTML — nice-to-have polish, not essential.

---

## 13. Key references / places to get unstuck

- Tangram: https://github.com/broadinstitute/Tangram (tutorials in `docs/` cover exactly this reference-mapping-to-space workflow)
- cell2location tutorial (same dataset pairing, useful for validating your Tangram results even if you don't implement cell2location): https://cell2location.readthedocs.io/en/latest/notebooks/cell2location_tutorial.html
- squidpy docs (spatial stats, neighborhood enrichment, plotting): https://squidpy.readthedocs.io
- scanpy docs: https://scanpy.readthedocs.io
- nf-core developer guidelines (module/subworkflow conventions worth mimicking even if you don't join nf-core formally): https://nf-co.re/docs/contributing/guidelines
- nf-test docs: https://www.nf-test.com

---

## 14. First message to send the agent to kick this off

> "Let's start Day 1–2: set up the repo skeleton from Section 4, then help me write and test `bin/sc_qc.py` standalone on the reference h5ad before we touch Nextflow at all."
