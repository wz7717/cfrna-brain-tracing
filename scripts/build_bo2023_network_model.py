#!/usr/bin/env python
from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.run_bo2023_loso_validation import read_vsd_matrix  # noqa: E402
from scripts.run_bo2023_network_correlation_validation import (  # noqa: E402
    build_group_reference,
    select_group_discriminative_genes,
)
from scripts.run_bo2023_v2_loso_validation import (  # noqa: E402
    DEFAULT_GENE_MAP,
    DEFAULT_MATRIX,
    DEFAULT_SAMPLE_INFO,
    map_matrix_to_symbols,
)


DEFAULT_OUT = ROOT / "data" / "models" / "bo2023_saleem_network_top200_model.npz"


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the production Bo2023 SaleemNetworks correlation model.")
    parser.add_argument("--matrix", type=Path, default=DEFAULT_MATRIX)
    parser.add_argument("--sample-info", type=Path, default=DEFAULT_SAMPLE_INFO)
    parser.add_argument("--gene-map", type=Path, default=DEFAULT_GENE_MAP)
    parser.add_argument("--sample-sheet", default="mfas5_819samples_phenSet4")
    parser.add_argument("--endpoint", default="SaleemNetworks")
    parser.add_argument("--top-n-genes", type=int, default=200)
    parser.add_argument("--full-loso-summary", type=Path)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()

    raw = read_vsd_matrix(args.matrix)
    matrix = map_matrix_to_symbols(raw, args.gene_map)
    ann = pd.read_excel(args.sample_info, sheet_name=args.sample_sheet)
    ann["sample_id"] = ann["No."].astype(str).str.strip()
    ann["endpoint_label"] = ann[args.endpoint].fillna("NA").astype(str).str.strip()
    ann = ann[ann["sample_id"].isin(set(matrix.columns))].copy()
    labels = ann.set_index("sample_id").reindex(matrix.columns)["endpoint_label"].to_numpy(dtype=str)
    groups = sorted(set(labels))
    values = matrix.to_numpy(dtype=np.float32)
    reference, training = build_group_reference(values, labels, groups, heldout_idx=-1)
    selected, audit = select_group_discriminative_genes(values, groups, training, args.top_n_genes)
    genes = matrix.index.astype(str).to_numpy(dtype=str)[selected]
    selected_reference = reference[selected, :]
    args.out.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        args.out,
        genes=genes,
        networks=np.asarray(groups, dtype=str),
        reference=selected_reference.astype(np.float32),
        fisher_scores=audit["fisher_score"].to_numpy(dtype=np.float32),
    )
    metadata: dict[str, object] = {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "dataset": "Bo2023_WangLab_VSD_region",
        "endpoint": args.endpoint,
        "n_training_samples": int(matrix.shape[1]),
        "n_networks": int(len(groups)),
        "n_genes": int(len(selected)),
        "method": "fold-validated discriminative-gene Pearson correlation",
        "input_requirement": "sample values must be on the same Bo2023 VSD-normalized scale",
    }
    if args.full_loso_summary and args.full_loso_summary.exists():
        validation = json.loads(args.full_loso_summary.read_text(encoding="utf-8"))
        metadata["full_loso_validation"] = validation.get("routes", {}).get(
            "network_discriminative_correlation_top200", {}
        )
    args.out.with_suffix(".json").write_text(json.dumps(metadata, ensure_ascii=False, indent=2), encoding="utf-8")
    audit.assign(gene_symbol=genes).to_csv(
        args.out.with_name(args.out.stem + "_genes.csv"), index=False, encoding="utf-8-sig"
    )
    print(json.dumps(metadata, ensure_ascii=False, indent=2))
    print(f"Model written to: {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
