# AnnData, Scanpy, and RAPIDS handoff

This document describes the intended boundary between `perturb-data-lab` and
standard single-cell analysis tools.

`perturb-data-lab` owns corpus materialization, canonical metadata, sparse count
storage, feature alignment for loaders, durable per-dataset HVG rankings, and
training/runtime access patterns. After a dataset or subset has been exported as
`AnnData`, Scanpy or RAPIDS should own the usual analysis preprocessing steps.

## Tool Boundary

Use `corpus.to_anndata(...)` when one dataset or one deterministic subset will
fit in memory and you want Scanpy or RAPIDS to own the next analysis steps.

Use Scanpy or RAPIDS for:

- normalization and log transforms for analysis objects
- PCA and alternative dimensionality reduction
- neighbors, UMAP, clustering, and plotting
- exploratory marker ranking or differential expression
- batch correction and other analysis-side transforms

Use corpus-native `perturb_data_lab.pp` helpers for:

- durable or recalculated per-dataset HVG rankings
- streamed summary statistics over the on-disk corpus
- quick QA/debug checks without constructing `AnnData`
- bounded-memory fallback workflows when Scanpy/RAPIDS handoff is not practical

Important current limits:

- `to_anndata(...)` is per-dataset only and eager; it is not a backed or fully on-disk Scanpy workflow.
- Corpus extraction is counts-only: it builds CSR `adata.X` and does not add normalized or log-transformed layers.
- Scanpy, RAPIDS, and `rapids-singlecell` are user-managed optional dependencies.
- Internal `pp` should not grow into a full replacement for Scanpy/RAPIDS.

## Imports

```python
from perturb_data_lab.loaders import load_corpus, select_obs_indices
from perturb_data_lab.pp import calculate_hvgs, run_pca
```

If you want centered `method="incremental_pca"`, install the optional PCA
dependency first:

```bash
pip install ".[pca]"
```

## Dry-run before eager AnnData construction

Load extra canonical metadata columns up front if you want to stratify on them
or include them in `adata.obs`:

```python
corpus = load_corpus(
    "/path/to/corpus",
    extra_metadata_columns=["donor_id", "batch_id"],
)

estimate = corpus.to_anndata(
    dataset_id="replogle_k562",
    obs_columns=["perturb_label", "donor_id"],
    var_columns=["gene_id"],
    dry_run=True,
    max_memory_bytes=8_000_000_000,
    on_exceed="warn",
)

print(estimate["n_obs"], estimate["n_vars"], estimate["nnz"])
print(estimate["csr_memory_bytes"])
print(estimate["selected_row_index_summary"])
```

- `dry_run=True` returns shape, `nnz`, CSR memory estimates, metadata footprint estimates, and memory-guard status without materializing `adata.X`.
- Set `on_exceed="raise"` when you want the same request to fail before eager construction if `max_memory_bytes` is exceeded.

## Deterministic observation selection

### Random subset

```python
random_selection = select_obs_indices(
    corpus,
    dataset_id="replogle_k562",
    strategy="random",
    max_cells=20_000,
    seed=17,
)
```

### Stratified subset

```python
stratified_selection = corpus.select_obs_indices(
    dataset_id="replogle_k562",
    strategy="stratified",
    max_cells=20_000,
    stratify_by=["perturb_label"],
    seed=17,
)
```

### Balanced subset with provenance

```python
balanced_selection = corpus.select_obs_indices(
    dataset_id="replogle_k562",
    strategy="balanced",
    max_cells=20_000,
    stratify_by=["perturb_label", "donor_id"],
    max_per_group=500,
    drop_null_groups=True,
    seed=17,
)

balanced_selection.write_provenance("./artifacts/selections/replogle-balanced")
```

- Returned `row_indices` are corpus-global row indices and can be passed directly into `to_anndata(...)` or `run_pca(...)`.
- Provenance outputs should go to repo-local real directories such as `./artifacts/...`, never to `data/`, `pertTF/`, or `perturb/`.

## Eager AnnData handoff

`to_anndata(...)` exports raw counts only. Any normalization, log transform,
dimensionality reduction, clustering, or plotting below is intentionally owned by
the Scanpy/RAPIDS-side workflow, not by the materialized corpus.

```python
adata = corpus.to_anndata(
    dataset_id="replogle_k562",
    row_indices=balanced_selection.row_indices,
    obs_columns=["perturb_label", "donor_id"],
)
```

`adata.obs` always includes stable provenance fields such as `dataset_id`,
`dataset_index`, `global_row_index`, and `local_row_index`; requested
`obs_columns` are added alongside them.

`adata.var` includes dataset-local feature identifiers plus global feature IDs.

## Scanpy example

```python
import scanpy as sc

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2000)
sc.tl.pca(adata)
sc.pp.neighbors(adata)
sc.tl.umap(adata)
sc.tl.leiden(adata)
```

Use this route for ordinary single-dataset or subset analysis when `adata.X`
fits in CPU memory.

## RAPIDS example

```python
import rapids_singlecell as rsc

rsc.get.anndata_to_GPU(adata)
rsc.pp.normalize_total(adata, target_sum=1e4)
rsc.pp.log1p(adata)
rsc.pp.highly_variable_genes(adata, n_top_genes=2000)
rsc.tl.pca(adata)
rsc.pp.neighbors(adata)
rsc.tl.umap(adata)
rsc.tl.leiden(adata)
```

Use this route when the exported subset fits GPU memory and the needed
`rapids-singlecell` functions support the matrix representation you are using.
For very large matrices, first reduce the cells or genes exported to `AnnData`.

## Corpus-native HVG and fallback PCA

The repo still keeps small streamed helpers because they work directly from the
on-disk corpus. The most important durable artifact is the per-dataset
`hvg.parquet` ranking table created during materialization or by `recalc-hvg`.

```python
hvg_frame = calculate_hvgs(
    corpus,
    dataset_id="replogle_k562",
    batch_size=1024,
    n_hvg=2000,
)
```

Use streamed PCA only when you specifically need a bounded-memory fallback over
the corpus instead of a normal Scanpy/RAPIDS analysis object:

```python
fit_selection = corpus.select_obs_indices(
    dataset_id="replogle_k562",
    strategy="balanced",
    max_cells=20_000,
    stratify_by=["perturb_label"],
    max_per_group=500,
    seed=17,
)

transform_selection = corpus.select_obs_indices(
    dataset_id="replogle_k562",
    strategy="random",
    max_cells=80_000,
    seed=23,
)

result = run_pca(
    corpus,
    dataset_id="replogle_k562",
    method="incremental_pca",
    batch_size=1024,
    n_components=50,
    hvg_frame=hvg_frame,
    fit_row_indices=fit_selection.row_indices,
    transform_row_indices=transform_selection.row_indices,
    max_dense_batch_bytes=2_000_000_000,
    output_dir="./artifacts/pp/replogle-ipca",
    overwrite=True,
)
```

- `fit_row_indices` lets you fit on a bounded deterministic subset while `transform_row_indices` controls which rows receive embeddings.
- Omit `transform_row_indices` to transform all rows in the requested dataset.
- `max_dense_batch_bytes` guards the dense `batch_size x selected_features` working set before outputs are written.
- Treat `run_pca(...)` as a fallback/debug path, not the preferred full analysis stack.

## Future large-scale handoff

The preferred large-scale direction is an AnnData-Zarr or Dask-backed CSR bridge
from the existing Zarr matrix artifacts. The current Zarr layout already has the
same CSR pieces expected by AnnData:

```text
current corpus Zarr      AnnData CSR group
counts              ->   X/data
indices             ->   X/indices
row_offsets         ->   X/indptr
```

This path should support per-dataset lazy handoff before a multi-dataset handoff.

Multi-dataset AnnData export needs an explicit column coordinate system. Current
materialized sparse rows store dataset-local feature indices, while a combined
AnnData object needs one shared feature axis. That means a multi-dataset export
must remap local feature indices to corpus-global feature IDs before writing or
presenting `X`.

Lance should remain optimized for training/runtime row reads for now. Modifying
AnnData backed internals to read Lance directly is not the first target because
it is more fragile than exporting or viewing the Zarr CSR arrays through standard
AnnData/Dask conventions.
