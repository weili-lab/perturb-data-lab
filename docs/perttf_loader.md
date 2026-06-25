# pertTF Loading Demo

This page shows the smallest useful pertTF-style loader demo on the canonicalized
Marson + Xorion corpus from the earlier walkthroughs.

Goal: prove that `perturb-data-lab` can produce one valid
`PertTFPairedBatchLoader` batch directly from the demo corpus, without touching an
external `pertTF` checkout.

## Before you start

You should already have a canonicalized corpus from either:

- [Bash Materialization](bash_demo.md)
- [Jupyter Materialization](jupyter_demo.md)

This page assumes the corpus lives at:

```text
./artifacts/demo_corpus
```

## One-command smoke test

If you just want proof that one pertTF-style batch can be built, run:

```bash
python scripts/perttf_demo_smoke.py ./artifacts/demo_corpus \
  --output-json ./artifacts/perttf_smoke.json
```

Expected outcome:

- the corpus loads successfully
- one batch is produced
- the script prints batch keys and tensor shapes
- `./artifacts/perttf_smoke.json` records the exact source/target pairs used in the smoke test

## Minimal Python example

```python
from perturb_data_lab.loaders import (
    PertTFAdapterConfig,
    PertTFPairedBatchLoader,
    load_corpus,
)

corpus = load_corpus("./artifacts/demo_corpus")

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
    seed=0,
)

batch = next(iter(loader))
print(sorted(batch.keys()))
print(batch["gene_ids"].shape)
print(batch["values"].shape)
print(batch["target_values_next"].shape)
print(batch["index"].tolist())
print(batch["next_index"].tolist())
```

Why these settings?

- `control_labels=("ctrl",)` matches the demo canonicalization output
- `dataset_id -> dataset` keeps pairing within each dataset
- `cell_context -> celltype` keeps pairing within the same cell context
- `mask_ratio=0.0` makes the smoke batch easier to inspect by keeping `values`
  equal to the unmasked source counts
- `seq_len=64` and `batch_size=4` keep the demo fast

## What comes back in one batch

The most important fields are:

- `gene_ids`: sampled tokenizer IDs for the source cell
- `values`: source expression values after masking
- `target_values`: unmasked source expression values
- `target_values_next`: paired target-cell values aligned to the sampled source genes
- `sf`, `sf_next`: source and target size factors
- `index`, `next_index`: corpus-global row indices for the source and paired target cells
- `perturbation_labels`, `celltype_labels`, `batch_labels`, `dataset_labels`: integer label IDs for the configured metadata labels

You will also see:

- `next_gene_ids`: currently aligned to `gene_ids` in this paired format
- `ps`, `ps_next`: placeholder pertTF side tensors retained for compatibility

## Inspect the source/target pairs

The easiest way to make the batch human-readable is to decode the metadata for
`index` and `next_index`.

```python
source_meta = corpus.take_metadata(
    batch["index"].tolist(),
    columns=["dataset_id", "perturb_label", "cell_context", "batch_id"],
)
target_meta = corpus.take_metadata(
    batch["next_index"].tolist(),
    columns=["dataset_id", "perturb_label", "cell_context", "batch_id"],
)

for i in range(len(batch["index"])):
    print(
        f"{source_meta['dataset_id'][i]} {source_meta['perturb_label'][i]} "
        f"-> {target_meta['dataset_id'][i]} {target_meta['perturb_label'][i]}"
    )
```

In the demo corpus, this usually gives same-dataset pairs such as:

- treated Marson row -> treated Marson row or matched control row
- control Xorion row -> treated Xorion row

The exact rows depend on the loader seed and pairing policy.

## Why native gene axes still work

Marson and Xorion keep their own native feature axes in the demo corpus.
`PertTFPairedBatchLoader` does **not** require a single shared raw gene order.

Instead:

- `load_corpus()` reconstructs a feature registry across datasets
- genes are converted to runtime token IDs
- the batch is built on those token IDs, not on a hard-coded shared raw matrix layout

That is why the same demo corpus can support:

- per-dataset pertTF-style loading here
- cross-dataset AnnData export later with `var_join="inner"`

## Common gotchas

- Use `dataset_id`, not `dataset_index`, in `label_fields` when you want a dataset label for pairing or reporting. The loader label vocabulary is string-based.
- If loader construction says there are no valid control rows, your `control_labels` do not match canonical `perturb_label` values.
- If HVG-based sampling fails, re-run canonicalization/HVG generation or switch to a non-HVG sampling mode.
- If too many rows are dropped, inspect nulls in the metadata columns named in `label_fields`.

## Next step

Once this batch smoke passes, continue to [Scanpy & RAPIDS](scanpy_rapids.md) for
the AnnData handoff path.
