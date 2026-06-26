from __future__ import annotations

import argparse
import shutil
from pathlib import Path


DEFAULT_DATASETS = (
    "marson_d2_rest",
    "xorion_hct116_dual_guide",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Copy the reviewed demo final schemas into a materialized corpus.",
    )
    parser.add_argument("--corpus", required=True, help="Path to the materialized corpus root.")
    parser.add_argument(
        "--schema-root",
        default=str(Path(__file__).resolve().parent.parent / "examples" / "demo_canonicalization"),
        help="Directory containing <dataset>.final-schema.yaml files.",
    )
    parser.add_argument(
        "--dataset-id",
        dest="dataset_ids",
        action="append",
        help="Optional dataset id to install. Repeat to limit to a subset.",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite an existing final-schema.yaml in the corpus.",
    )
    return parser.parse_args()


def resolve_meta_dir(corpus_root: Path, dataset_id: str) -> Path:
    candidates = (
        corpus_root / "meta" / dataset_id,
        corpus_root / dataset_id / "meta",
    )
    for path in candidates:
        if path.exists():
            return path
    raise FileNotFoundError(f"Could not find metadata directory for dataset '{dataset_id}' under {corpus_root}")


def install_one(schema_root: Path, corpus_root: Path, dataset_id: str, overwrite: bool) -> Path:
    source = schema_root / f"{dataset_id}.final-schema.yaml"
    if not source.exists():
        raise FileNotFoundError(f"Missing demo schema: {source}")
    meta_dir = resolve_meta_dir(corpus_root, dataset_id)
    target = meta_dir / "final-schema.yaml"
    if target.exists() and not overwrite:
        raise FileExistsError(f"Refusing to overwrite existing schema: {target}")
    shutil.copyfile(source, target)
    return target


def main() -> None:
    args = parse_args()
    corpus_root = Path(args.corpus).resolve()
    schema_root = Path(args.schema_root).resolve()
    dataset_ids = tuple(args.dataset_ids or DEFAULT_DATASETS)

    if not corpus_root.exists():
        raise FileNotFoundError(f"Corpus root does not exist: {corpus_root}")
    if not schema_root.exists():
        raise FileNotFoundError(f"Schema root does not exist: {schema_root}")

    for dataset_id in dataset_ids:
        target = install_one(schema_root, corpus_root, dataset_id, args.overwrite)
        print(f"installed {dataset_id} -> {target}")


if __name__ == "__main__":
    main()
