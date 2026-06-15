#!/usr/bin/env python
from __future__ import annotations

import argparse
import json
import sys
import zipfile
from pathlib import Path
from typing import Any

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from core.bo2023_region_tracing import trace_bo2023_secondary_regions  # noqa: E402
from core.network_tracing import trace_network_expression  # noqa: E402
from core.region_resolution import annotate_region_candidates  # noqa: E402


DEFAULT_DATA_DIR = ROOT / "data" / "ahba_human_rnaseq"
DEFAULT_ZIP_DIR = DEFAULT_DATA_DIR / "raw_zips"
DEFAULT_MAPPING_OUT = ROOT / "docs" / "deliverables" / "public_dataset_anatomical_label_mapping.csv"
DEFAULT_BO2023_MAPPING = ROOT / "docs" / "deliverables" / "bo2023_public_anatomical_label_mapping.csv"
DEFAULT_OUTDIR = ROOT / "results" / "ahba_human_rnaseq_external_validation_20260603"
DEFAULT_DB = ROOT / "cfrna_source_tracing.db"


AHBA_MAIN_STRUCTURE_MAPPING = {
    "FL": {
        "public_major_anatomy": "frontal lobe",
        "allowed_bo2023_lobes": ["Frontal"],
        "allowed_bo2023_networks": [
            "Lateral Prefrontal Cortex",
            "Orbitomedial Prefrontal Cortex (OMPFC)",
            "Frontal (agranular frontal motor areas)",
            "Operculum/Insula",
        ],
        "supported": True,
        "confidence": "medium",
    },
    "PL": {
        "public_major_anatomy": "parietal lobe",
        "allowed_bo2023_lobes": ["Parietal"],
        "allowed_bo2023_networks": ["Parietal, and Parieto-occipital region", "Operculum/Insula"],
        "supported": True,
        "confidence": "medium",
    },
    "TL": {
        "public_major_anatomy": "temporal lobe",
        "allowed_bo2023_lobes": ["Temporal"],
        "allowed_bo2023_networks": ["Temporal", "Occipital/Temporal"],
        "supported": True,
        "confidence": "medium",
    },
    "OL": {
        "public_major_anatomy": "occipital lobe",
        "allowed_bo2023_lobes": ["Occipital"],
        "allowed_bo2023_networks": ["Occipital/Temporal"],
        "supported": True,
        "confidence": "medium",
    },
    "CgG": {
        "public_major_anatomy": "cingulate cortex",
        "allowed_bo2023_lobes": ["Cingulate", "Frontal"],
        "allowed_bo2023_networks": ["Cingulate gyrus", "Orbitomedial Prefrontal Cortex (OMPFC)"],
        "supported": True,
        "confidence": "medium",
    },
    "Ins": {
        "public_major_anatomy": "insula",
        "allowed_bo2023_lobes": ["Insula"],
        "allowed_bo2023_networks": ["Operculum/Insula"],
        "supported": True,
        "confidence": "high",
    },
    "PHG": {
        "public_major_anatomy": "parahippocampal cortex",
        "allowed_bo2023_lobes": ["Temporal"],
        "allowed_bo2023_networks": ["Temporal"],
        "supported": True,
        "confidence": "high",
        "allowed_bo2023_regions": ["TF", "TFO"],
    },
    "Str": {
        "public_major_anatomy": "striatum",
        "allowed_bo2023_lobes": ["Subcortical"],
        "allowed_bo2023_networks": ["Subcortical"],
        "supported": True,
        "confidence": "high",
    },
    "GP": {
        "public_major_anatomy": "globus pallidus",
        "allowed_bo2023_lobes": ["Subcortical"],
        "allowed_bo2023_networks": ["Subcortical"],
        "allowed_bo2023_regions": ["GPeGPi"],
        "supported": True,
        "confidence": "high",
    },
    "CbCx": {
        "public_major_anatomy": "cerebellar cortex",
        "allowed_bo2023_lobes": [],
        "allowed_bo2023_networks": [],
        "supported": False,
        "confidence": "unsupported",
    },
}

AHBA_SUBSTRUCTURE_OVERRIDES = {
    "Caudate": {"allowed_bo2023_regions": ["cd"], "public_fine_anatomy_label": "caudate"},
    "Putamen": {"allowed_bo2023_regions": ["pu"], "public_fine_anatomy_label": "putamen"},
    "GP": {"allowed_bo2023_regions": ["GPeGPi"], "public_fine_anatomy_label": "globus pallidus"},
    "str_V1": {"allowed_bo2023_regions": ["V1"], "public_fine_anatomy_label": "primary visual cortex"},
    "pest_V2": {"allowed_bo2023_regions": ["V2"], "public_fine_anatomy_label": "secondary visual cortex"},
    "PHG": {"allowed_bo2023_regions": ["TF", "TFO"], "public_fine_anatomy_label": "parahippocampal cortex"},
    "Insula": {"allowed_bo2023_regions": ["Ia", "Ial", "Iapl", "Iapm", "Id", "Ig"], "public_fine_anatomy_label": "insula"},
    "CgG": {"allowed_bo2023_regions": ["23a", "23b", "23c", "24a", "24a'", "24b", "24b'", "24c", "24c'", "25", "31", "32"], "public_fine_anatomy_label": "cingulate cortex"},
    "PrG": {"allowed_bo2023_regions": ["F1", "F2", "F3", "F4", "F5", "F6", "F7"], "public_fine_anatomy_label": "precentral/motor cortex"},
    "PoG-l": {"allowed_bo2023_regions": ["3a/b", "5", "SII"], "public_fine_anatomy_label": "postcentral/somatosensory cortex"},
    "PoG-cs": {"allowed_bo2023_regions": ["3a/b", "5", "SII"], "public_fine_anatomy_label": "postcentral/somatosensory cortex"},
}


def split_pipe(value: Any) -> list[str]:
    if value is None or pd.isna(value):
        return []
    return [item.strip() for item in str(value).split("|") if item.strip()]


def load_zip_payloads(zip_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    metadata_rows: list[pd.DataFrame] = []
    matrix_parts: list[pd.DataFrame] = []
    for zip_path in sorted(zip_dir.glob("*.zip")):
        donor = zip_path.stem.replace("_rnaseq", "")
        with zipfile.ZipFile(zip_path) as zf:
            annot = pd.read_csv(zf.open("SampleAnnot.csv"))
            sample_ids = [f"{donor}|{sample}" for sample in annot["RNAseq_sample_name"].astype(str)]
            annot = annot.copy()
            annot["ahba_donor"] = donor
            annot["ahba_sample_uid"] = sample_ids
            tpm = pd.read_csv(zf.open("RNAseqTPM.csv"), header=None)
            if tpm.shape[1] - 1 != len(sample_ids):
                raise ValueError(f"{zip_path.name} TPM/sample count mismatch: {tpm.shape[1] - 1} vs {len(sample_ids)}")
            tpm = tpm.rename(columns={0: "gene_symbol"})
            tpm.columns = ["gene_symbol", *sample_ids]
            matrix_parts.append(tpm.set_index("gene_symbol"))
            metadata_rows.append(annot)
    metadata = pd.concat(metadata_rows, ignore_index=True)
    matrix = pd.concat(matrix_parts, axis=1)
    matrix.index = matrix.index.astype(str)
    matrix = matrix.groupby(matrix.index).mean()
    return metadata, matrix


def mapping_for_ahba_label(main_structure: str, sub_structure: str) -> dict[str, Any]:
    base = dict(AHBA_MAIN_STRUCTURE_MAPPING.get(main_structure, {}))
    if not base:
        base = {
            "public_major_anatomy": str(main_structure),
            "allowed_bo2023_lobes": [],
            "allowed_bo2023_networks": [],
            "supported": False,
            "confidence": "unsupported",
        }
    override = AHBA_SUBSTRUCTURE_OVERRIDES.get(sub_structure, {})
    allowed_regions = override.get("allowed_bo2023_regions", base.get("allowed_bo2023_regions", []))
    fine = override.get("public_fine_anatomy_label", str(sub_structure))
    return {
        "public_dataset": "AHBA_RNAseq",
        "public_label_type": "AHBA main_structure/sub_structure",
        "public_label": f"{main_structure}::{sub_structure}",
        "public_main_structure": main_structure,
        "public_sub_structure": sub_structure,
        "public_major_anatomy": base.get("public_major_anatomy", ""),
        "public_fine_anatomy_label": fine,
        "allowed_bo2023_lobes": " | ".join(base.get("allowed_bo2023_lobes", [])),
        "allowed_bo2023_networks": " | ".join(base.get("allowed_bo2023_networks", [])),
        "allowed_bo2023_regions": " | ".join(allowed_regions),
        "mapping_confidence": base.get("confidence", "unsupported"),
        "supported_for_accuracy": bool(base.get("supported", False)),
        "accuracy_level": "exact_region" if allowed_regions else ("coarse_anatomy" if base.get("supported", False) else "unsupported"),
        "mapping_rule": "AHBA normal human anatomical label harmonized to Bo2023 lobe/network; exact region only for stable shared labels.",
        "caution_note": "Human-to-macaque external validation is approximate; report coarse anatomy/network accuracy, not strict Bo2023 exact-region accuracy.",
    }


def write_public_mapping(metadata: pd.DataFrame, out_path: Path) -> pd.DataFrame:
    rows = []
    for row in metadata[["main_structure", "sub_structure"]].drop_duplicates().itertuples(index=False):
        rows.append(mapping_for_ahba_label(str(row.main_structure), str(row.sub_structure)))
    mapping = pd.DataFrame(rows).sort_values(["public_main_structure", "public_sub_structure"])
    out_path.parent.mkdir(parents=True, exist_ok=True)
    mapping.to_csv(out_path, index=False, encoding="utf-8-sig")
    return mapping


def expression_frame(matrix: pd.DataFrame, sample_id: str) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "gene_symbol": matrix.index.astype(str),
            "tpm_value": pd.to_numeric(matrix[sample_id], errors="coerce").fillna(0.0).to_numpy(),
        }
    )


def top_values(rows: list[dict[str, Any]], key: str, k: int) -> list[str]:
    return [str(row.get(key, "")) for row in sorted(rows, key=lambda r: int(r.get("rank", 999)))[:k]]


def hit_any(predicted: list[str], allowed: list[str]) -> bool | None:
    if not allowed:
        return None
    return any(item in set(allowed) for item in predicted)


def summarize_bool(series: pd.Series) -> float:
    valid = series.dropna()
    if valid.empty:
        return float("nan")
    return float(valid.astype(bool).mean())


def main() -> int:
    parser = argparse.ArgumentParser(description="Build public anatomical mapping and run AHBA RNA-seq external validation.")
    parser.add_argument("--zip-dir", type=Path, default=DEFAULT_ZIP_DIR)
    parser.add_argument("--data-dir", type=Path, default=DEFAULT_DATA_DIR)
    parser.add_argument("--mapping-out", type=Path, default=DEFAULT_MAPPING_OUT)
    parser.add_argument("--bo2023-mapping", type=Path, default=DEFAULT_BO2023_MAPPING)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument("--db", type=Path, default=DEFAULT_DB)
    parser.add_argument("--atlas-id", type=int, default=4)
    parser.add_argument("--topk-region", type=int, default=15)
    parser.add_argument("--prepare-only", action="store_true")
    args = parser.parse_args()

    args.data_dir.mkdir(parents=True, exist_ok=True)
    args.outdir.mkdir(parents=True, exist_ok=True)
    metadata, matrix = load_zip_payloads(args.zip_dir)
    metadata_path = args.data_dir / "ahba_human_rnaseq_sample_metadata_242.csv"
    matrix_path = args.data_dir / "ahba_human_rnaseq_tpm_gene_symbol_matrix.tsv"
    metadata.to_csv(metadata_path, index=False, encoding="utf-8-sig")
    matrix.to_csv(matrix_path, sep="\t", encoding="utf-8")
    public_mapping = write_public_mapping(metadata, args.mapping_out)
    if args.prepare_only:
        print(json.dumps({"n_samples": len(metadata), "matrix_shape": matrix.shape, "mapping_rows": len(public_mapping)}, indent=2))
        return 0

    bo2023_mapping = pd.read_csv(args.bo2023_mapping)
    bo2023_mapping["bo2023_region_id"] = bo2023_mapping["bo2023_region_id"].astype(str)
    region_to_lobe = bo2023_mapping.set_index("bo2023_region_id")["bo2023_lobe"].to_dict()
    region_to_public_major = bo2023_mapping.set_index("bo2023_region_id")["public_major_anatomy"].to_dict()
    mapping_lookup = {
        row.public_label: row._asdict()
        for row in public_mapping.itertuples(index=False)
    }

    network_rows: list[dict[str, Any]] = []
    region_rows: list[dict[str, Any]] = []
    sample_rows: list[dict[str, Any]] = []
    trace_meta: dict[str, Any] = {}
    for i, row in enumerate(metadata.itertuples(index=False), start=1):
        sample_id = str(row.ahba_sample_uid)
        public_label = f"{row.main_structure}::{row.sub_structure}"
        label_map = mapping_lookup[public_label]
        expr = expression_frame(matrix, sample_id)
        network_out = trace_network_expression(expr, min_overlap_fraction=0.20)
        region_out = trace_bo2023_secondary_regions(expr, network_out, str(args.db), int(args.atlas_id), topk=int(args.topk_region))
        region_out = annotate_region_candidates(region_out, network_out)
        trace_meta[sample_id] = {"network": network_out.get("meta", {}), "region": region_out.get("meta", {})}

        common = {
            "sample_id": sample_id,
            "ahba_donor": row.ahba_donor,
            "sample_name": row.sample_name,
            "main_structure": row.main_structure,
            "sub_structure": row.sub_structure,
            "public_label": public_label,
            "public_major_anatomy": label_map["public_major_anatomy"],
            "public_fine_anatomy_label": label_map["public_fine_anatomy_label"],
            "mapping_confidence": label_map["mapping_confidence"],
            "supported_for_accuracy": bool(label_map["supported_for_accuracy"]),
            "accuracy_level": label_map["accuracy_level"],
        }
        for out_row in network_out.get("results", [])[:10]:
            network_rows.append({**common, **out_row})
        for out_row in region_out.get("results", [])[: int(args.topk_region)]:
            region_rows.append({**common, **out_row})

        network_top3 = top_values(network_out.get("results", []), "network_id", 3)
        region_top3 = top_values(region_out.get("results", []), "region_id", 3)
        group_top3 = top_values(region_out.get("results", []), "resolution_group", 3)
        region_top1 = region_top3[0] if region_top3 else ""
        predicted_lobe_top1 = region_to_lobe.get(region_top1, "")
        predicted_public_major_top1 = region_to_public_major.get(region_top1, "")
        allowed_lobes = split_pipe(label_map["allowed_bo2023_lobes"])
        allowed_networks = split_pipe(label_map["allowed_bo2023_networks"])
        allowed_regions = split_pipe(label_map["allowed_bo2023_regions"])
        supported = bool(label_map["supported_for_accuracy"])
        top_region = (region_out.get("results", []) or [{}])[0]
        sample_rows.append(
            {
                **common,
                "network_top1": network_top3[0] if network_top3 else "",
                "network_top3": " | ".join(network_top3),
                "region_top1": region_top1,
                "region_top3": " | ".join(region_top3),
                "resolution_group_top1": group_top3[0] if group_top3 else "",
                "resolution_group_top3": " | ".join(group_top3),
                "predicted_lobe_top1": predicted_lobe_top1,
                "predicted_public_major_top1": predicted_public_major_top1,
                "allowed_bo2023_lobes": label_map["allowed_bo2023_lobes"],
                "allowed_bo2023_networks": label_map["allowed_bo2023_networks"],
                "allowed_bo2023_regions": label_map["allowed_bo2023_regions"],
                "network_top1_hit": hit_any(network_top3[:1], allowed_networks) if supported else None,
                "network_top3_hit": hit_any(network_top3, allowed_networks) if supported else None,
                "region_lobe_top1_hit": hit_any([predicted_lobe_top1], allowed_lobes) if supported else None,
                "region_top1_exact_hit": hit_any(region_top3[:1], allowed_regions) if supported and allowed_regions else None,
                "region_top3_exact_hit": hit_any(region_top3, allowed_regions) if supported and allowed_regions else None,
                "top1_resolution_tier": str(top_region.get("resolution_tier", "")),
                "top1_manual_review_recommended": bool(top_region.get("manual_review_recommended", False)),
                "network_overlap_genes": network_out.get("meta", {}).get("n_overlap_genes"),
                "network_overlap_fraction": network_out.get("meta", {}).get("overlap_fraction"),
                "region_overlap_genes": region_out.get("meta", {}).get("n_overlap_genes"),
            }
        )
        print(f"[{i:03d}/{len(metadata)}] {sample_id} {public_label}")

    network_df = pd.DataFrame(network_rows)
    region_df = pd.DataFrame(region_rows)
    sample_df = pd.DataFrame(sample_rows)
    network_df.to_csv(args.outdir / "ahba_rnaseq_network_tracing_per_sample_top10.csv", index=False, encoding="utf-8-sig")
    region_df.to_csv(args.outdir / "ahba_rnaseq_region_tracing_per_sample_top15.csv", index=False, encoding="utf-8-sig")
    sample_df.to_csv(args.outdir / "ahba_rnaseq_external_validation_sample_summary.csv", index=False, encoding="utf-8-sig")
    (args.outdir / "ahba_rnaseq_external_validation_trace_meta.json").write_text(json.dumps(trace_meta, ensure_ascii=False, indent=2), encoding="utf-8")

    supported = sample_df[sample_df["supported_for_accuracy"].astype(bool)].copy()
    exact_evaluable = supported[supported["allowed_bo2023_regions"].fillna("").astype(str).str.len() > 0].copy()
    metrics = {
        "n_samples_total": int(len(sample_df)),
        "n_samples_supported_for_accuracy": int(len(supported)),
        "n_samples_unsupported": int(len(sample_df) - len(supported)),
        "unsupported_labels": sorted(sample_df.loc[~sample_df["supported_for_accuracy"].astype(bool), "public_label"].unique().tolist()),
        "network_top1_accuracy_coarse": summarize_bool(supported["network_top1_hit"]),
        "network_top3_accuracy_coarse": summarize_bool(supported["network_top3_hit"]),
        "region_lobe_top1_accuracy_coarse": summarize_bool(supported["region_lobe_top1_hit"]),
        "exact_region_evaluable_samples": int(len(exact_evaluable)),
        "region_top1_exact_accuracy_on_exact_mapped_labels": summarize_bool(exact_evaluable["region_top1_exact_hit"]),
        "region_top3_exact_accuracy_on_exact_mapped_labels": summarize_bool(exact_evaluable["region_top3_exact_hit"]),
        "low_resolution_top1_fraction": float((sample_df["top1_resolution_tier"].astype(str) == "low_resolution").mean()),
        "manual_review_top1_fraction": float(sample_df["top1_manual_review_recommended"].astype(bool).mean()),
        "interpretation": "AHBA RNA-seq is human normal brain external validation. Metrics are label-harmonized coarse anatomy/network accuracy, not strict macaque exact-region accuracy.",
    }
    by_label = (
        sample_df.groupby(["main_structure", "sub_structure", "public_label", "supported_for_accuracy", "accuracy_level"])
        .agg(
            n_samples=("sample_id", "count"),
            network_top1_accuracy=("network_top1_hit", summarize_bool),
            network_top3_accuracy=("network_top3_hit", summarize_bool),
            region_lobe_top1_accuracy=("region_lobe_top1_hit", summarize_bool),
            exact_region_top1_accuracy=("region_top1_exact_hit", summarize_bool),
            exact_region_top3_accuracy=("region_top3_exact_hit", summarize_bool),
            dominant_network_top1=("network_top1", lambda s: s.value_counts().index[0] if len(s) else ""),
            dominant_region_top1=("region_top1", lambda s: s.value_counts().index[0] if len(s) else ""),
            dominant_lobe_top1=("predicted_lobe_top1", lambda s: s.value_counts().index[0] if len(s) else ""),
        )
        .reset_index()
    )
    by_label.to_csv(args.outdir / "ahba_rnaseq_external_validation_by_label.csv", index=False, encoding="utf-8-sig")
    (args.outdir / "ahba_rnaseq_external_validation_metrics.json").write_text(json.dumps(metrics, ensure_ascii=False, indent=2), encoding="utf-8")

    report = [
        "# AHBA RNA-seq external validation",
        "",
        "Metrics are computed after harmonizing AHBA human anatomical labels to Bo2023 lobe/network labels.",
        "Cerebellar cortex is unsupported because the current Bo2023 reference does not include cerebellum.",
        "",
        json.dumps(metrics, ensure_ascii=False, indent=2),
    ]
    (args.outdir / "ahba_rnaseq_external_validation_report.md").write_text("\n".join(report), encoding="utf-8")
    print(json.dumps(metrics, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
