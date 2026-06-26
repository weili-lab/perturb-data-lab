# perturb-data-lab

A demo-first walkthrough for perturb-seq data preprocessing and loading.

`perturb-data-lab` turns raw `.h5ad` files into sparse on-disk corpora, adds
reviewed canonical metadata, and exposes a common runtime API for model training
or downstream analysis.

<div class="grid cards" markdown>

-   :material-download:{ .lg .middle } __Installation__

    ---

    Set up your environment with one command.

    [:octicons-arrow-right-24: Get started](installation.md)

-   :material-console:{ .lg .middle } __Bash Demo__

    ---

    Follow the terminal-only demo from raw files to a loaded corpus.

    [:octicons-arrow-right-24: Bash materialization](bash_demo.md)

-   :material-notebook:{ .lg .middle } __Jupyter Demo__

    ---

    Run the notebook version with metadata inspection cells.

    [:octicons-arrow-right-24: Jupyter materialization](jupyter_demo.md)

-   :material-tag:{ .lg .middle } __Canonicalization__

    ---

    Understand how raw metadata becomes canonical labels.

    [:octicons-arrow-right-24: Schema decisions](demo_canonicalization.md)

-   :material-dna:{ .lg .middle } __pertTF Loader__

    ---

    Feed canonicalized corpora into pertTF-style paired batches.

    [:octicons-arrow-right-24: pertTF loading](perttf_loader.md)

-   :material-chart-scatter-plot:{ .lg .middle } __Scanpy & RAPIDS__

    ---

    Export AnnData, run Scanpy preprocessing, and explore GPU acceleration.

    [:octicons-arrow-right-24: Scanpy/RAPIDS](scanpy_rapids.md)

</div>

## Quick demo overview

The demo starts from two small HuggingFace-hosted `.h5ad` subsets (one Marson
CRISPRi, one Xorion dual-guide) and ends with a canonicalized corpus that can
feed `PertTFPairedBatchLoader` and export Dask-backed AnnData for Scanpy.

```text
HuggingFace download
  -> inspect
  -> materialize (create + append)
  -> copy reviewed schemas
  -> canonicalize
  -> load_corpus()
  -> pertTF paired batches
  -> AnnData + Scanpy
```

Start with [Installation](installation.md) and then follow the
[Bash demo](bash_demo.md) or [Jupyter demo](jupyter_demo.md).

## Reference docs

- [Inspection & Materialization](inspect_materialize.md) — how raw `.h5ad`
  files become materialized corpora
- [Canonicalization Handbook](canonicalization_handbook.md) — the full schema
  contract and canonicalization reference
- [Backend Notes](backend_note.md) — storage backend policy and selection
- [AnnData Handoff API](anndata_scanpy_handoff.md) — corpus-to-AnnData export
  and Scanpy/RAPIDS boundary
