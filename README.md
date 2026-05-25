# perturb-data-lab

`perturb-data-lab` is a corpus-first preprocessing and loading workspace for
large perturb-seq datasets. It turns raw `.h5ad` files into sparse on-disk
corpora, adds reviewed canonical metadata, and exposes a common runtime API for
model training or downstream analysis.

The current repo is intentionally lightweight. It is not a polished Python
package with a large compatibility surface; it is a practical framework for
inspecting real datasets, materializing count matrices safely, and loading the
resulting corpora.

## Current Workflow

```text
raw h5ad
  -> inspect
  -> materialize counts + raw sidecars
  -> draft/review final-schema.yaml
  -> canonicalize
  -> load_corpus()
  -> optional pertTF loader, AnnData handoff, or streamed pp helpers
```

The main public CLI commands are:

- `inspect`: profile one or more `.h5ad` files and write `dataset-summary.yaml`
- `materialize`: write sparse count data, raw metadata sidecars, manifests, and corpus registration records
- `draft-schema`: create a starting `draft-schema.yaml` from materialized raw sidecars
- `canonicalize`: apply reviewed `final-schema.yaml` files and write canonical obs/var parquet files
- `corpus-validate`: check corpus index and manifest consistency
- `corpus-compose`: build a new federated corpus by symlinking or copying whole dataset directories
- `recalc-hvg`: recalculate per-dataset `hvg.parquet` ranking tables
- `corpus-gc`: remove unregistered per-dataset directories after failed corpus writes

## Supported Corpus Routes

The maintained storage backends are:

- `lance`
- `zarr`

The maintained topologies are:

- `aggregate`: one shared matrix root with per-dataset metadata under `meta/<dataset_id>/`
- `federated`: one self-contained directory per dataset under `<dataset_id>/`

Recommended default: aggregate Lance.

Use Zarr when chunked array artifacts or node-local staging are useful. Use
federated topology when whole-dataset movement, recomposition, or isolation is
more important than a single aggregate matrix object.

See `docs/backend_note.md` for the short backend policy note.

## Repo Layout

```text
perturb-data-lab/
├── README.md
├── DESIGN.md
├── docs/
│   ├── inspect_materialize.md
│   ├── canonicalization_handbook.md
│   ├── perttf_loader.md
│   ├── anndata_scanpy_handoff.md
│   └── backend_note.md
├── examples/
├── scripts/
│   └── optional smoke and validation utilities
├── src/perturb_data_lab/
│   ├── cli.py
│   ├── inspectors/
│   ├── materializers/
│   ├── canonical/
│   ├── loaders/
│   └── pp/
└── tests/
```

Important code areas:

- `src/perturb_data_lab/cli.py`: public command-line entry points
- `src/perturb_data_lab/inspectors/`: backed `.h5ad` metadata and count-source inspection
- `src/perturb_data_lab/materializers/`: count-matrix streaming, backend writers, manifests, and corpus registration
- `src/perturb_data_lab/canonical/`: schema drafting, schema contracts, transforms, and canonicalization runner
- `src/perturb_data_lab/loaders/`: `load_corpus()`, expression readers, metadata index, feature registry, samplers, and pertTF adapter
- `src/perturb_data_lab/pp/`: streamed stats, HVG, PCA, and differential-expression helpers

## Documentation Map

- Architecture and data flow: `DESIGN.md`
- Inspection and materialization how-to: `docs/inspect_materialize.md`
- Canonical schema review and canonicalization details: `docs/canonicalization_handbook.md`
- pertTF-specific paired loader usage: `docs/perttf_loader.md`
- AnnData, Scanpy, deterministic subsampling, and streamed PCA: `docs/anndata_scanpy_handoff.md`
- Backend policy note: `docs/backend_note.md`

## Minimal CLI Example

Inspect one dataset:

```bash
PYTHONPATH=src python -m perturb_data_lab.cli inspect \
  --source /path/to/dataset.h5ad \
  --dataset-id my_dataset \
  --output-dir ./artifacts/review
```

Materialize it into a new aggregate Lance corpus:

```bash
PYTHONPATH=src python -m perturb_data_lab.cli materialize \
  --mode create \
  --source /path/to/dataset.h5ad \
  --dataset-id my_dataset \
  --inspection-summary ./artifacts/review/my_dataset/dataset-summary.yaml \
  --output-corpus ./artifacts/corpus \
  --backend lance \
  --topology aggregate
```

Draft and review canonical metadata:

```bash
PYTHONPATH=src python -m perturb_data_lab.cli draft-schema \
  --corpus ./artifacts/corpus
```

Edit `draft-schema.yaml` into `final-schema.yaml`, then canonicalize:

```bash
PYTHONPATH=src python -m perturb_data_lab.cli canonicalize \
  --corpus ./artifacts/corpus
```

Load the corpus:

```python
from perturb_data_lab.loaders import load_corpus

corpus = load_corpus("./artifacts/corpus")
expr = corpus.read_expression([0, 1, 2])
meta = corpus.take_metadata([0, 1, 2], columns=["dataset_id", "perturb_label"])
```

## Runtime API Sketch

```python
from perturb_data_lab.loaders import load_corpus

corpus = load_corpus("/path/to/corpus")

corpus.set_sampler(batch_size=128, seed=0)

for batch in corpus.loader(seq_len=1024, processing="gpu", num_workers=4):
    train_step(batch)
```

Useful runtime methods:

- `corpus.read_expression(global_row_indices)`: read sparse expression rows
- `corpus.take_metadata(global_row_indices, columns=[...])`: read canonical metadata columns
- `corpus.inspect_batch(global_row_indices, metadata_columns=[...])`: inspect expression plus metadata together
- `corpus.loader(...)`: build an iterable sparse batch loader
- `corpus.to_anndata(...)`: eager counts-only AnnData export for one dataset or subset

By default, `load_corpus()` loads core canonical metadata. Pass
`extra_metadata_columns=[...]` when a downstream workflow needs additional
canonical obs columns.

## Safety Notes

- Treat raw `.h5ad` sources as read-only.
- For large `.h5ad` files, run inspection and materialization on Slurm CPU using the project environment.
- Do not write outputs into protected symlink roots such as `data/`, `pertTF/`, or `perturb/`.
- Materialization does not apply obs filtering. If only a subset of cells should be materialized, create a pre-filtered `.h5ad` first.
- Canonicalization is a reviewed metadata step. Do not blindly promote every raw field into a canonical field.
