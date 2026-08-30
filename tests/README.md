# Tests

## nf-test (Phase 6, if time allows)

`nf-test` runs individual Nextflow processes in isolation against fixture inputs and asserts on their
outputs. The distinction from the `test` profile matters:

- **`-profile test`** is an integration test: does the whole DAG execute end to end?
- **`nf-test`** is a unit test: given this input, does *this one process* produce the expected output?

The integration test is the higher priority — it is what CI runs and what a stranger runs to verify
the repo. Unit tests are the polish that shows you think about testability at the module level.

Worthwhile assertions, in priority order:

1. `sc_qc` output h5ad has strictly fewer cells than its input (filtering actually filtered).
2. `sc_qc` output `.obs` contains every expected QC column (`n_genes_by_counts`, `total_counts`,
   `pct_counts_mt`, `predicted_doublet`).
3. `sc_cluster_annotate` output has a `leiden` column in `.obs` and an `X_umap` in `.obsm`.
4. `deconvolution` proportions matrix has one row per spot and sums to ~1 per row.

Assertion 4 is the most valuable one in the whole suite: it is the check that catches a silently
broken mapping, which is the failure mode most likely to produce a confident-looking wrong answer.

```bash
# install
curl -fsSL https://get.nf-test.com | bash && mv nf-test ~/.local/bin/

# run
nf-test test tests/
```
