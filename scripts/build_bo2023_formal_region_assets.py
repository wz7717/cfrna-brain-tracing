#!/usr/bin/env python
from __future__ import annotations

import json
import sys
from datetime import datetime, timezone
from itertools import combinations
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from core.reference_projection import compute_logcpm, map_index_to_symbols, read_bo2023_gene_matrix, read_gene_map
from core.bo2023_metadata import assert_unique_region_network_mapping
from scripts.build_bo2023_reference_projector import DEFAULT_COUNTS, DEFAULT_SAMPLE_INFO, DEFAULT_VSD
from scripts.run_ahba_projected_vsd_external_validation import read_bo_metadata
from scripts.run_bo2023_network_correlation_validation import select_group_discriminative_genes
from scripts.run_bo2023_projected_vsd_exact_region import DEFAULT_CLEANED_GENE_MAP
from scripts.run_bo2023_resolution_tier_validation import (
    build_resolution_groups,
    normalize_resolution_annotations,
)


OUT_MATRIX = ROOT / "data/models/bo2023_formal_region_logcpm_reference_matrix.npz"
OUT_BEAMS = ROOT / "data/models/bo2023_formal_region_beam_gene_panels.json"


def main() -> int:
    gene_map = read_gene_map(DEFAULT_CLEANED_GENE_MAP)
    counts, _ = map_index_to_symbols(read_bo2023_gene_matrix(DEFAULT_COUNTS, dtype="float32"), gene_map)
    vsd, _ = map_index_to_symbols(read_bo2023_gene_matrix(DEFAULT_VSD, dtype="float32"), gene_map)
    genes = sorted(set(counts.index.astype(str)) & set(vsd.index.astype(str)))
    metadata = read_bo_metadata(DEFAULT_SAMPLE_INFO, "mfas5_819samples_phenSet4", "Region", "SaleemNetworks")
    samples = [sample for sample in counts.columns.astype(str) if sample in set(metadata["sample_id"])]
    metadata = metadata[metadata["sample_id"].isin(samples)].copy()
    assert_unique_region_network_mapping(metadata)
    metadata["region_key"] = metadata["network_id"].astype(str) + "::" + metadata["region_id"].astype(str)
    logcpm = compute_logcpm(counts.loc[genes, samples]).astype("float32")
    values = logcpm.to_numpy(dtype=np.float32)
    sample_pos = {sample: idx for idx, sample in enumerate(samples)}
    training = {
        str(region): np.asarray([sample_pos[sample] for sample in rows["sample_id"].astype(str)], dtype=int)
        for region, rows in metadata.groupby("region_key")
    }
    regions = sorted(training)
    networks = sorted(metadata["network_id"].astype(str).unique())
    region_network = metadata.drop_duplicates("region_key").set_index("region_key")["network_id"].astype(str).to_dict()
    region_display = metadata.drop_duplicates("region_key").set_index("region_key")["region_id"].astype(str).to_dict()
    reference = np.column_stack([values[:, training[region]].mean(axis=1) for region in regions]).astype("float32")

    np.savez_compressed(
        OUT_MATRIX,
        genes=np.asarray(genes),
        regions=np.asarray(regions),
        display_regions=np.asarray([region_display[region] for region in regions]),
        networks=np.asarray([region_network[region] for region in regions]),
        matrix=reference,
    )

    panels: dict[str, dict[str, object]] = {}
    for beam in combinations(networks, 3):
        candidates = sorted(region for region in regions if region_network[region] in set(beam))
        candidate_training = {region: training[region] for region in candidates}
        rows, _ = select_group_discriminative_genes(values, candidates, candidate_training, 200)
        assignments = {region: region_network[region] for region in candidates}
        annotations, _ = build_resolution_groups(
            values,
            candidates,
            candidate_training,
            assignments,
            rows,
            min_resolution_samples=8,
            min_merge_samples=3,
            min_pair_errors=2,
            min_confusion_rate=0.15,
            similarity_threshold=0.95,
            merge_similarity_threshold=0.90,
            max_group_size=8,
        )
        annotations = normalize_resolution_annotations(annotations)
        key = "||".join(sorted(beam))
        panels[key] = {
            "networks": list(sorted(beam)),
            "candidate_regions": candidates,
            "genes": [genes[int(index)] for index in rows],
            "annotations": annotations,
        }
    payload = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "route": "formal full-fit region assets",
        "region_identity": "network_id::region_id",
        "region_ontology": "110 canonical Bo2023 region IDs with one parent Network per region",
        "gene_selection": "Top200 Fisher-like between-region/within-region score",
        "n_training_samples": len(samples),
        "n_genes": len(genes),
        "n_regions": len(regions),
        "n_networks": len(networks),
        "n_beams": len(panels),
        "beams": panels,
    }
    OUT_BEAMS.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps({key: payload[key] for key in payload if key != "beams"}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
