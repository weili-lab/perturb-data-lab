# Installation

This page gets you from zero to a working `perturb-data-lab` environment ready
to run the [Bash](bash_demo.md) or [Jupyter](jupyter_demo.md) demo.

## 1. Clone the repository

```bash
git clone https://github.com/weililab/perturb-data-lab.git
cd perturb-data-lab
git checkout docs/github-pages-demo
```

## 2. Set up a Python environment

`perturb-data-lab` requires **Python ≥ 3.10**. We recommend using conda or mamba
to create an isolated environment:

```bash
# Using conda/mamba (recommended)
conda create -n pdl-demo python=3.11 -y
conda activate pdl-demo

# Or using venv
python3 -m venv venv
source venv/bin/activate
```

## 3. Install the package

Install in editable mode with core dependencies:

```bash
pip install -e .
```

## 4. Verify the installation

Run a quick import check to confirm everything is wired correctly:

```python
from perturb_data_lab.loaders import load_corpus
from perturb_data_lab.materializers.paths import resolve_corpus_paths

print("perturb-data-lab is ready.")
```

## 5. Download the demo data

The demo uses two small `.h5ad` subsets hosted on HuggingFace:

- **Repository**: [weililab/perturb-data-lab-demo](https://huggingface.co/datasets/weililab/perturb-data-lab-demo)
- **Datasets**: Marson D2 Rest (CRISPRi, 2.7K cells) and Xorion HCT116 (dual-guide, 2.7K cells)

### Option A — using the bundled download script (recommended)

```bash
# Uses huggingface_hub by default
python scripts/download_demo_data.py --output-dir ./demo_data

# Or with wget (no extra packages)
python scripts/download_demo_data.py --output-dir ./demo_data --method wget
```

The script downloads `demo_marson_d2_rest.h5ad`,
`demo_xorion_hct116_dual_guide.h5ad`, and `checksums.txt`, then verifies file
integrity.

### Option B — using huggingface-cli

```bash
pip install huggingface_hub
huggingface-cli download weililab/perturb-data-lab-demo \
  --local-dir ./demo_data \
  --repo-type dataset
```

### Option C — direct download

```bash
mkdir -p ./demo_data/h5ad
wget -O ./demo_data/h5ad/demo_marson_d2_rest.h5ad \
  https://huggingface.co/datasets/weililab/perturb-data-lab-demo/resolve/main/h5ad/demo_marson_d2_rest.h5ad
wget -O ./demo_data/h5ad/demo_xorion_hct116_dual_guide.h5ad \
  https://huggingface.co/datasets/weililab/perturb-data-lab-demo/resolve/main/h5ad/demo_xorion_hct116_dual_guide.h5ad
```

## 6. Install optional dependencies

The demo workflow benefits from a few extra packages:

```bash
# For Scanpy, Dask, and download helpers
pip install -e ".[demo]"

# For building the documentation site locally
pip install -e ".[docs]"
```

The `demo` group includes `scanpy`, `dask[dataframe]`, `requests`, and `tqdm`.
The `docs` group includes `mkdocs-material` and related extensions.

## Environment check

After installation, the following should all succeed:

```python
import anndata
import dask
import scanpy
import perturb_data_lab
from perturb_data_lab.loaders import load_corpus, PertTFAdapterConfig, PertTFPairedBatchLoader
```

## Next steps

- **[Bash Demo](bash_demo.md)** — copy-paste CLI walkthrough
- **[Jupyter Demo](jupyter_demo.md)** — interactive Python walkthrough
- **[Canonicalization](demo_canonicalization.md)** — learn the two demo schema decisions
