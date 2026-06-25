#!/usr/bin/env python
from __future__ import annotations

import argparse
import json
from pathlib import Path

from perturb_data_lab.loaders import (
    PertTFAdapterConfig,
    PertTFPairedBatchLoader,
    load_corpus,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Smoke-test one pertTF-style batch from a canonicalized demo corpus."
    )
    parser.add_argument("corpus_root", help="Path to a canonicalized corpus root")
    parser.add_argument("--batch-size", type=int, default=4)
    parser.add_argument("--seq-len", type=int, default=64)
    parser.add_argument("--hvg-top-k", type=int, default=2000)
    parser.add_argument("--output-json", type=Path, default=None)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    corpus_root = Path(args.corpus_root).resolve()
    corpus = load_corpus(corpus_root)

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
        batch_size=args.batch_size,
        seq_len=args.seq_len,
        config=config,
        sampling_mode="hvg",
        hvg_top_k=args.hvg_top_k,
        num_workers=0,
        seed=0,
    )
    batch = next(iter(loader))

    pair_columns = ["dataset_id", "perturb_label", "cell_context", "batch_id"]
    source_meta = corpus.take_metadata(batch["index"].tolist(), columns=pair_columns)
    target_meta = corpus.take_metadata(batch["next_index"].tolist(), columns=pair_columns)

    summary = {
        "corpus_root": str(corpus_root),
        "dataset_ids": list(corpus.dataset_ids),
        "n_cells": len(corpus.metadata_index),
        "loader_settings": {
            "batch_size": args.batch_size,
            "seq_len": args.seq_len,
            "sampling_mode": "hvg",
            "hvg_top_k": args.hvg_top_k,
            "num_workers": 0,
            "mask_ratio": 0.0,
            "control_labels": ["ctrl"],
            "pairing_group_labels": ["dataset", "celltype"],
        },
        "effective_rows": {
            "label_rows": int(len(loader.effective_label_row_indices)),
            "source_rows": int(len(loader.effective_source_indices)),
            "target_candidate_rows": int(len(loader.effective_target_candidate_indices)),
        },
        "batch_keys": sorted(batch.keys()),
        "batch_shapes": {
            key: list(batch[key].shape)
            for key in [
                "gene_ids",
                "next_gene_ids",
                "values",
                "target_values",
                "target_values_next",
                "sf",
                "sf_next",
                "index",
                "next_index",
                "perturbation_labels",
                "perturbation_labels_next",
                "celltype_labels",
                "celltype_labels_next",
                "batch_labels",
                "batch_labels_next",
                "dataset_labels",
                "dataset_labels_next",
            ]
        },
        "pairs": [
            {
                "source_global_row": int(batch["index"][i]),
                "target_global_row": int(batch["next_index"][i]),
                "source_dataset_id": source_meta["dataset_id"][i],
                "source_perturb_label": source_meta["perturb_label"][i],
                "target_dataset_id": target_meta["dataset_id"][i],
                "target_perturb_label": target_meta["perturb_label"][i],
                "cell_context": source_meta["cell_context"][i],
                "batch_id": source_meta["batch_id"][i],
            }
            for i in range(len(batch["index"]))
        ],
    }

    print(f"Loaded corpus: {summary['n_cells']} cells across {', '.join(summary['dataset_ids'])}")
    print(
        "Effective row pool: "
        f"labels={summary['effective_rows']['label_rows']} "
        f"source={summary['effective_rows']['source_rows']} "
        f"target={summary['effective_rows']['target_candidate_rows']}"
    )
    print(f"Batch keys: {', '.join(summary['batch_keys'])}")
    print(
        "Shapes: "
        f"gene_ids={summary['batch_shapes']['gene_ids']} "
        f"values={summary['batch_shapes']['values']} "
        f"target_values={summary['batch_shapes']['target_values']} "
        f"target_values_next={summary['batch_shapes']['target_values_next']}"
    )
    print("Pairs:")
    for pair in summary["pairs"]:
        print(
            "  "
            f"{pair['source_dataset_id']} {pair['source_perturb_label']} "
            f"-> {pair['target_dataset_id']} {pair['target_perturb_label']} "
            f"[{pair['cell_context']} | {pair['batch_id']}]"
        )

    if args.output_json is not None:
        args.output_json.write_text(json.dumps(summary, indent=2), encoding="utf-8")
        print(f"Wrote summary: {args.output_json}")


if __name__ == "__main__":
    main()
