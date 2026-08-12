#!/usr/bin/env python
from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser(description="Prepare immutable Bo2023 metadata for the P0-3 DESeq2 audit.")
    parser.add_argument("--metadata-xlsx", type=Path, required=True)
    parser.add_argument("--counts", type=Path, required=True)
    parser.add_argument("--gene-map", type=Path, required=True)
    parser.add_argument("--locked-panel", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    metadata = pd.read_excel(args.metadata_xlsx, sheet_name="mfas5_819samples_phenSet4", dtype=str)
    required = ["No.", "MonkeyID", "Side", "Batch", "Region", "SaleemNetworks", "SaleemNetworksAB"]
    missing = sorted(set(required) - set(metadata.columns))
    if missing:
        raise ValueError(f"metadata missing columns: {missing}")
    metadata = metadata[required].rename(
        columns={"No.": "sample_id", "MonkeyID": "donor_id", "SaleemNetworks": "network"}
    )
    for column in metadata.columns:
        metadata[column] = metadata[column].fillna("").astype(str).str.strip()
    metadata.loc[metadata["Region"].eq("10m"), "network"] = "Orbitomedial Prefrontal Cortex (OMPFC)"
    metadata.loc[metadata["Region"].eq("V2"), "network"] = "Occipital/Temporal"

    with args.counts.open("r", encoding="utf-8") as handle:
        count_samples = handle.readline().rstrip("\r\n").split("\t")
    if len(count_samples) != 819:
        raise ValueError(f"expected 819 count columns, found {len(count_samples)}")
    if metadata["sample_id"].duplicated().any():
        raise ValueError("metadata sample IDs are not unique")
    meta_samples = set(metadata["sample_id"])
    if set(count_samples) != meta_samples:
        raise ValueError(
            f"sample mismatch: counts_only={sorted(set(count_samples)-meta_samples)[:10]}, "
            f"metadata_only={sorted(meta_samples-set(count_samples))[:10]}"
        )
    metadata = metadata.set_index("sample_id").loc[count_samples].reset_index()
    if metadata[["donor_id", "network", "Region"]].eq("").any().any():
        raise ValueError("blank donor/network/region values remain")
    if metadata["network"].nunique() != 10 or metadata["donor_id"].nunique() != 9:
        raise ValueError(
            f"unexpected design: networks={metadata['network'].nunique()}, donors={metadata['donor_id'].nunique()}"
        )

    cross = pd.crosstab(metadata["donor_id"], metadata["network"])
    if (cross.sum(axis=0) == 0).any() or (cross.sum(axis=1) == 0).any():
        raise ValueError("empty donor or Network level")
    metadata.to_csv(args.output_dir / "bo2023_819_deseq2_metadata.tsv", sep="\t", index=False)
    cross.to_csv(args.output_dir / "donor_by_network_sample_counts.tsv", sep="\t")

    gene_map = pd.read_csv(args.gene_map, dtype=str)
    gene_required = ["Gene.stable.ID", "Gene.name", "Gene.type_ensembl"]
    if not set(gene_required).issubset(gene_map.columns):
        raise ValueError("gene map lacks required columns")
    gene_map = gene_map[gene_required].rename(
        columns={"Gene.stable.ID": "gene_id", "Gene.name": "gene_symbol", "Gene.type_ensembl": "gene_type"}
    )
    gene_map = gene_map.drop_duplicates("gene_id", keep="first")
    gene_map.to_csv(args.output_dir / "bo2023_gene_annotation.tsv", sep="\t", index=False)

    locked = pd.read_csv(args.locked_panel)
    locked[["gene_symbol", "fisher_score", "gene_index"]].to_csv(
        args.output_dir / "locked_network_top200.tsv", sep="\t", index=False
    )
    design_rank = int(pd.get_dummies(metadata[["donor_id", "network"]], drop_first=True).assign(intercept=1).astype(float).pipe(lambda x: __import__('numpy').linalg.matrix_rank(x.to_numpy())))
    expected_rank = 1 + (metadata["donor_id"].nunique() - 1) + (metadata["network"].nunique() - 1)
    if design_rank != expected_rank:
        raise ValueError(f"design is rank deficient: rank={design_rank}, expected={expected_rank}")

    manifest = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "analysis": "P0-3 donor-by-Network pseudobulk DESeq2 Network marker audit",
        "design": "sum raw counts within observed donor-Network cells; full ~ donor_id + network; reduced ~ donor_id; LRT for any Network effect",
        "effect_followup": "DESeq2-normalized pseudobulk log2 counts; highest Network versus other-Network mean; descriptive only",
        "prefilter": "pseudobulk raw count >=10 in at least 3 donor-Network units",
        "canonicalization": {"10m": "Orbitomedial Prefrontal Cortex (OMPFC)", "V2": "Occipital/Temporal"},
        "n_samples": len(metadata),
        "n_donors": int(metadata["donor_id"].nunique()),
        "n_networks": int(metadata["network"].nunique()),
        "design_rank": design_rank,
        "expected_design_rank": expected_rank,
        "inputs": {
            "metadata_xlsx": {"source_file": args.metadata_xlsx.name, "sha256": sha256(args.metadata_xlsx)},
            "counts": {"source_file": args.counts.name, "sha256": sha256(args.counts)},
            "gene_map": {"source_file": args.gene_map.name, "sha256": sha256(args.gene_map)},
            "locked_panel": {"source_file": args.locked_panel.name, "sha256": sha256(args.locked_panel)},
        },
    }
    (args.output_dir / "input_manifest.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    print(json.dumps({"status": "ok", "n_samples": len(metadata), "design_rank": design_rank, "network_counts": metadata["network"].value_counts().to_dict()}, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
