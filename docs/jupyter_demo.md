# Jupyter Materialization Demo

This page walks through the full demo workflow as a series of Python cells you
can run in a Jupyter notebook, IPython session, or plain script. Each cell is
self-contained and builds on the previous one.

## Prerequisites

- [Installation](installation.md) completed and `pip install -e ".[demo]"` succeeded.
- A working Python environment with `perturb_data_lab`, `anndata`, and `scanpy` importable.
- Demo data downloaded (see Cell 1).

---

### Cell 1 — Download demo data

```python
# Run the bundled download script (or use huggingface_hub / wget directly)
import subprocess
import sys

result = subprocess.run(
    [sys.executable, "scripts/download_demo_data.py", "--output-dir", "./demo_data"],
    capture_output=False,
)
```

After download, verify the files exist:

```python
from pathlib import Path

demo_root = Path("./demo_data")
for rel in ["h5ad/demo_marson_d2_rest.h5ad", "h5ad/demo_xorion_hct116_dual_guide.h5ad"]:
    f = demo_root / rel
    size_mb = f.stat().st_size / 1e6 if f.exists() else 0
    print(f"  {rel}: {size_mb:.1f} MB {'✓' if f.exists() else '✗ missing'}")
```

---

### Cell 2 — Inspect raw h5ad metadata

Inspection reads metadata profiles and samples matrix candidates without loading
full count matrices. The output `dataset-summary.yaml` tells you whether the
file is ready for materialization.

```python
from pathlib import Path
from perturb_data_lab.inspectors import inspect_target
from perturb_data_lab.inspectors.models import InspectionTarget

review_dir = Path("./artifacts/review")
review_dir.mkdir(parents=True, exist_ok=True)

datasets = [
    ("marson_d2_rest", "./demo_data/h5ad/demo_marson_d2_rest.h5ad"),
    ("xorion_hct116_dual_guide", "./demo_data/h5ad/demo_xorion_hct116_dual_guide.h5ad"),
]

for ds_id, source_path in datasets:
    result = inspect_target(
        InspectionTarget(dataset_id=ds_id, source_path=source_path, source_release=ds_id),
        review_dir,
    )
    print(f"\n{ds_id}:")
    print(f"  cells: {result.n_obs}, features: {result.n_vars}")
    print(f"  count source: {result.count_source}")
    print(f"  readiness: {result.materialization_readiness}")
    print(f"  obs fields ({len(result.obs_summary)}): {sorted(result.obs_summary.keys())[:10]}...")
```

---

### Cell 3 — Quick peek at obs metadata

Let's look at a few rows of raw metadata to understand what each dataset
contains before canonicalization:

```python
import anndata

for ds_id, source_path in datasets:
    adata = anndata.read_h5ad(source_path, backed="r")
    print(f"\n--- {ds_id} ---")
    print(f"Shape: {adata.shape}")
    print(f"obs columns: {list(adata.obs.columns)}")

    # Show a few key columns
    cols = [c for c in ["guide_id", "guide_type", "perturbed_gene_name",
                        "guide_target", "gene_target", "sample", "pass_guide_filter"]
            if c in adata.obs.columns]
    if cols:
        display_cols = cols[:6]
        print(adata.obs[display_cols].head(5).to_string())
```

Expected observations:

- **Marson**: has `guide_id`, `guide_type`, `perturbed_gene_name`, `lane_id`, QC columns
- **Xorion**: has `guide_target`, `gene_target`, `sample`, `pass_guide_filter`, QC columns

---

### Cell 4 — Materialize the corpus

Create an aggregate Lance corpus and stream sparse counts from both `.h5ad` files:

```python
from perturb_data_lab.materializers import DatasetMaterializer
from perturb_data_lab.materializers.models import OutputRoots
from perturb_data_lab.materializers.paths import resolve_corpus_paths
import yaml

corpus_root = Path("./artifacts/demo_corpus")
corpus_root.mkdir(parents=True, exist_ok=True)

# Reuse the summaries from the CLI inspection for simplicity
inspect_summary_dir = Path("./artifacts/review")

def load_summary(dataset_id):
    path = inspect_summary_dir / dataset_id / "dataset-summary.yaml"
    return yaml.safe_load(path.read_text())

# Step 4a — Create corpus with Marson
ds_id = "marson_d2_rest"
paths = resolve_corpus_paths("aggregate", corpus_root, ds_id)
summary = load_summary(ds_id)

m = DatasetMaterializer(
    source_path="./demo_data/h5ad/demo_marson_d2_rest.h5ad",
    inspection_summary_path=str(inspect_summary_dir / ds_id / "dataset-summary.yaml"),
    output_roots=OutputRoots(
        metadata_root=str(paths.meta_root),
        matrix_root=str(paths.matrix_root),
    ),
    dataset_id=ds_id,
    backend="lance",
    topology="aggregate",
    corpus_index_path=str(corpus_root / "corpus-index.yaml"),
    corpus_id="demo_corpus",
    register=True,
    mode="create",
    dataset_index=0,
    global_row_start=0,
)
manifest = m.materialize()
print(f"Created {ds_id}: {manifest.cell_count} cells, {manifest.feature_count} features")
```

```python
# Step 4b — Append Xorion
ds_id = "xorion_hct116_dual_guide"
paths = resolve_corpus_paths("aggregate", corpus_root, ds_id)
summary = load_summary(ds_id)

# Read the current corpus index to learn global_row_start
index_doc = yaml.safe_load((corpus_root / "corpus-index.yaml").read_text())
first_entry = index_doc["datasets"][0]
next_global_start = first_entry["global_end"]
next_dataset_index = first_entry["dataset_index"] + 1

m2 = DatasetMaterializer(
    source_path="./demo_data/h5ad/demo_xorion_hct116_dual_guide.h5ad",
    inspection_summary_path=str(inspect_summary_dir / ds_id / "dataset-summary.yaml"),
    output_roots=OutputRoots(
        metadata_root=str(paths.meta_root),
        matrix_root=str(paths.matrix_root),
    ),
    dataset_id=ds_id,
    backend="lance",
    topology="aggregate",
    corpus_index_path=str(corpus_root / "corpus-index.yaml"),
    corpus_id="demo_corpus",
    register=True,
    mode="append",
    dataset_index=next_dataset_index,
    global_row_start=next_global_start,
)
manifest2 = m2.materialize()
print(f"Appended {ds_id}: {manifest2.cell_count} cells, {manifest2.feature_count} features")
```

??? tip "Using the CLI instead"
    The CLI handles append bookkeeping (`dataset_index`, `global_row_start`)
    automatically. For production use, the [Bash demo](bash_demo.md) CLI commands
    are simpler. The Python API is shown here for notebook-style transparency.

---

### Cell 5 — Install reviewed schemas

Copy the reviewed demo final schemas into the corpus. These contain the
biological decisions that turn raw metadata columns into canonical labels:

```python
# Install schemas using the bundled helper
import subprocess
subprocess.run(
    [sys.executable, "scripts/install_demo_schemas.py", "--corpus", str(corpus_root)],
    check=True,
)
```

```python
# Verify they landed
for ds_id in ["marson_d2_rest", "xorion_hct116_dual_guide"]:
    schema_path = corpus_root / "meta" / ds_id / "final-schema.yaml"
    print(f"  {ds_id}: {'✓' if schema_path.exists() else '✗ missing'}")
```

Read the quick decision sheets that explain what each schema does:

```python
examples_root = Path("examples/demo_canonicalization")
for ds_id in ["marson_d2_rest", "xorion_hct116_dual_guide"]:
    hints = yaml.safe_load((examples_root / f"{ds_id}.schema-hints.yaml").read_text())
    print(f"\n--- {ds_id} hints ---")
    for field, info in hints.items():
        print(f"  {field}: {info}")
```

For a deeper explanation, see [Canonicalization](demo_canonicalization.md).

---

### Cell 6 — Canonicalize

Apply the reviewed schemas to produce canonical obs/var metadata:

```python
# Dry-run first
subprocess.run(
    [sys.executable, "-m", "perturb_data_lab.cli", "canonicalize",
     "--corpus", str(corpus_root), "--dry-run"],
    check=True,
)
```

```python
# Real canonicalization
subprocess.run(
    [sys.executable, "-m", "perturb_data_lab.cli", "canonicalize",
     "--corpus", str(corpus_root)],
    check=True,
)
print("Canonicalization complete.")
```

---

### Cell 7 — Validate and load the corpus

```python
# Validate
subprocess.run(
    [sys.executable, "-m", "perturb_data_lab.cli", "corpus-validate",
     str(corpus_root / "corpus-index.yaml")],
    check=True,
)

# Load
from perturb_data_lab.loaders import load_corpus

corpus = load_corpus(str(corpus_root))
print(f"Datasets: {corpus.dataset_ids}")
print(f"Total cells: {len(corpus.metadata_index)}")
print(f"Global vocab size: {corpus.feature_registry.global_vocab_size}")
```

---

### Cell 8 — Inspect canonical metadata

Look at a few rows to confirm canonical labels are sensible:

```python
import polars as pl

# Grab first 10 rows from each dataset
marson_rows = corpus.take_metadata(
    list(range(0, 10)),
    columns=["dataset_id", "perturb_label", "condition", "perturb_type",
             "cell_context", "batch_id", "gene_id"],
)
xorion_rows = corpus.take_metadata(
    list(range(2720, 2730)),
    columns=["dataset_id", "perturb_label", "condition", "perturb_type",
             "cell_context", "batch_id", "gene_id"],
)

print("--- Marson ---")
print(marson_rows)
print("\n--- Xorion ---")
print(xorion_rows)
```

```python
# Quick stats
meta = corpus.metadata_index
all_perturb_labels = meta.get_column("perturb_label")

# Count control vs treated
ctrl_count = (pl.Series(all_perturb_labels) == "ctrl").sum()
total = len(all_perturb_labels)
print(f"Control rows: {ctrl_count} / {total} ({100*ctrl_count/total:.0f}%)")

# Unique conditions
conditions = meta.get_column("condition")
unique = set(conditions.to_list())
print(f"Unique conditions: {len(unique)}")
print(f"Sample conditions: {sorted(unique)[:10]}...")
```

---

### Cell 9 — PertTF loader preview

Produce one paired batch to confirm the loader works end-to-end:

```python
from perturb_data_lab.loaders import PertTFAdapterConfig, PertTFPairedBatchLoader

config = PertTFAdapterConfig(
    label_fields={
        "perturb_label": "perturbation",
        "cell_context": "celltype",
        "batch_id": "batch",
        "dataset_id": "dataset",
    },
    perturbation_label="perturbation",
    control_labels=("ctrl",),
    pairing_group_labels=("dataset", "celltype"),
    mask_ratio=0.0,
)

loader = PertTFPairedBatchLoader(
    corpus,
    batch_size=4,
    seq_len=64,
    config=config,
    sampling_mode="hvg",
    hvg_top_k=2000,
    num_workers=0,
)

batch = next(iter(loader))
print("Batch keys:", sorted(batch.keys()))
print(f"  gene_ids shape: {batch['gene_ids'].shape}")
print(f"  values shape: {batch['values'].shape}")
print(f"  target_values shape: {batch['target_values'].shape}")
print(f"  target_values_next shape: {batch['target_values_next'].shape}")
print(f"  index: {batch['index'].tolist()}")
print(f"  next_index: {batch['next_index'].tolist()}")

# Decode perturbation labels
src_labels = corpus.take_metadata(
    batch["index"].tolist(),
    columns=["perturb_label", "dataset_id"],
)
tgt_labels = corpus.take_metadata(
    batch["next_index"].tolist(),
    columns=["perturb_label", "dataset_id"],
)
print("\nSource pairs:")
for i in range(len(batch["index"])):
    src = src_labels.row(i)
    tgt = tgt_labels.row(i)
    print(f"  {src[1]} {src[0]} -> {tgt[1]} {tgt[0]}")
```

For more details on loader configuration, see [pertTF Loading](perttf_loader.md).

---

### Cell 10 — AnnData handoff (Dask-backed)

Export the combined corpus as a Dask-backed AnnData with inner-join on features:

```python
adata = corpus.to_anndata_lazy(
    dataset_id=list(corpus.dataset_ids),
    obs_columns=["perturb_label", "condition", "cell_context", "batch_id"],
    chunk_rows=1024,
    var_join="inner",
)

print(f"adata shape: {adata.shape}")
print(f"adata.X type: {type(adata.X)}")
print(f"obs columns: {list(adata.obs.columns)}")
print(f"var columns: {list(adata.var.columns)}")

# Show the intersection gene axis (first few)
print(f"\nIntersection features: {adata.n_vars}")
print(adata.var.head(5).to_string())
```

```python
# Quick Scanpy smoke on the Dask-backed AnnData
import scanpy as sc

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2000)

print(f"HVG selected: {adata.var['highly_variable'].sum()}")

# Force a small computation to confirm Dask works
subset = adata[:100, :100].to_memory()
print(f"Subset shape: {subset.shape}")
print(f"Subset X nnz: {subset.X.nnz}")
```

For the full Scanpy/RAPIDS story, see [Scanpy & RAPIDS](scanpy_rapids.md).

---

## Re-running from scratch

Delete the artifacts directory and re-run from Cell 1:

```python
import shutil
shutil.rmtree("./artifacts", ignore_errors=True)
```

The demo data under `./demo_data/` does not need to be re-downloaded.

## Next steps

- **[Bash Demo](bash_demo.md)** — the CLI equivalent for scripting
- **[Canonicalization](demo_canonicalization.md)** — understand the two schema decisions
- **[pertTF Loading](perttf_loader.md)** — full loader configuration and batch fields
- **[Scanpy & RAPIDS](scanpy_rapids.md)** — Dask-backed AnnData, CPU Scanpy, GPU RAPIDS
