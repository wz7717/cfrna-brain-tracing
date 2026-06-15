#!/usr/bin/env python
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from core.bo2023_region_tracing import trace_bo2023_secondary_regions  # noqa: E402
from core.network_tracing import trace_network_expression  # noqa: E402
from core.region_resolution import annotate_region_candidates  # noqa: E402


DEFAULT_MATRIX = ROOT / "data" / "tcga_brain_tumor_expression" / "tcga_gbm_lgg_primary_tumor_tpm_unstranded_sample_mean.tsv"
DEFAULT_MANIFEST = ROOT / "data" / "tcga_brain_tumor_expression" / "tcga_gbm_lgg_gdc_star_counts_manifest.csv"
DEFAULT_LABELS = ROOT / "data" / "tcia_tcga_glioma_mri" / "tcia_tcga_glioma_mri_derived_labels.csv"
DEFAULT_OUTDIR = ROOT / "results" / "tcga_gbm_lgg_sample_mri_label_tracing_20260605"
DEFAULT_DB = ROOT / "cfrna_source_tracing.db"


NETWORK_TO_BROAD = {
    "Cingulate gyrus": "cingulate",
    "Frontal (agranular frontal motor areas)": "frontal",
    "Hippocampal formation": "medial_temporal",
    "Lateral Prefrontal Cortex": "frontal",
    "Occipital/Temporal": "occipital_temporal",
    "Operculum/Insula": "insula",
    "Orbitomedial Prefrontal Cortex (OMPFC)": "frontal",
    "Parietal, and Parieto-occipital region": "parietal_occipital",
    "Subcortical": "subcortical",
    "Temporal": "temporal",
    "Dorsolateral Prefrontal Cortex (DLPFC)": "frontal",
    "Ventrolateral Prefrontal Cortex (VLPFC)": "frontal",
    "Premotor Cortex": "frontal",
    "Primary Motor Cortex": "frontal",
    "Somatosensory Cortex": "parietal",
    "Posterior Parietal Cortex": "parietal",
    "Insula": "insula",
    "Temporal Auditory Cortex": "temporal",
    "Inferotemporal Cortex": "temporal",
    "Medial Temporal Cortex": "temporal",
    "Temporal Visual Cortex": "temporal",
    "Visual Cortex": "occipital",
}


def expression_frame(matrix: pd.DataFrame, sample_id: str) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "gene_symbol": matrix.index.astype(str),
            "tpm_value": pd.to_numeric(matrix[sample_id], errors="coerce").fillna(0.0).to_numpy(),
        }
    )


def sample_to_patient(sample_id: str) -> str:
    parts = str(sample_id).upper().split("-")
    return "-".join(parts[:3]) if len(parts) >= 3 else str(sample_id).upper()


def normalize_label(value: Any) -> str:
    return str(value or "").strip().lower().replace(" ", "_").replace("-", "_")


def split_truths(value: Any) -> set[str]:
    text = str(value or "").strip()
    if not text or text.lower() == "nan":
        return set()
    parts = []
    for chunk in text.replace(";", "|").replace(",", "|").split("|"):
        label = normalize_label(chunk)
        if label:
            parts.append(label)
    return set(parts)


def top_values(df: pd.DataFrame, column: str, sample_id: str, k: int) -> list[str]:
    sub = df[df["sample_id"].eq(sample_id)].sort_values("rank").head(k)
    return sub[column].astype(str).tolist() if len(sub) else []


def any_match(predicted: list[str], truth: set[str]) -> bool | None:
    if not truth:
        return None
    return any(normalize_label(item) in truth for item in predicted)


def broad_candidates(networks: list[str]) -> list[str]:
    out: list[str] = []
    for network in networks:
        broad = NETWORK_TO_BROAD.get(str(network), str(network))
        out.append(broad)
        if broad == "medial_temporal":
            out.append("temporal")
    return out


def load_labels(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame(columns=["patient_barcode"])
    labels = pd.read_csv(path)
    if "patient_barcode" not in labels.columns:
        raise ValueError(f"{path} must contain patient_barcode")
    labels["patient_barcode"] = labels["patient_barcode"].astype(str).str.upper().str.slice(0, 12)
    return labels.drop_duplicates("patient_barcode")


def metric(rows: list[dict[str, Any]], column: str) -> dict[str, Any]:
    vals = [row[column] for row in rows if row.get(column) is not None]
    n = len(vals)
    return {"n_evaluable": n, "n_correct": int(sum(bool(x) for x in vals)), "accuracy": float(sum(bool(x) for x in vals) / n) if n else None}


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Run sample-level TCGA-GBM/LGG tracing and evaluate against patient-level TCIA/BraTS MRI labels."
    )
    parser.add_argument("--matrix", type=Path, default=DEFAULT_MATRIX)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--labels", type=Path, default=DEFAULT_LABELS)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument("--db", type=Path, default=DEFAULT_DB)
    parser.add_argument("--atlas-id", type=int, default=4)
    parser.add_argument("--max-samples", type=int, default=None)
    parser.add_argument("--min-overlap-fraction", type=float, default=0.20)
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    matrix = pd.read_csv(args.matrix, sep="\t", index_col=0)
    matrix.index = matrix.index.astype(str)
    sample_ids = matrix.columns.astype(str).tolist()
    if args.max_samples is not None:
        sample_ids = sample_ids[: max(1, int(args.max_samples))]

    manifest = pd.read_csv(args.manifest)
    sample_project = (
        manifest[["sample_submitter_id", "project_id", "case_submitter_id"]]
        .drop_duplicates("sample_submitter_id")
        .set_index("sample_submitter_id")
        .to_dict("index")
    )
    labels = load_labels(args.labels)

    network_rows: list[dict[str, Any]] = []
    region_rows: list[dict[str, Any]] = []
    meta: dict[str, Any] = {}
    for i, sample_id in enumerate(sample_ids, start=1):
        expr = expression_frame(matrix, sample_id)
        network_out = trace_network_expression(expr, min_overlap_fraction=args.min_overlap_fraction)
        region_out = trace_bo2023_secondary_regions(expr, network_out, str(args.db), int(args.atlas_id), topk=15)
        region_out = annotate_region_candidates(region_out, network_out)
        meta[sample_id] = {"network": network_out.get("meta", {}), "region": region_out.get("meta", {})}
        for row in network_out.get("results", [])[:10]:
            network_rows.append({"sample_id": sample_id, "patient_barcode": sample_to_patient(sample_id), **row})
        for row in region_out.get("results", [])[:15]:
            region_rows.append({"sample_id": sample_id, "patient_barcode": sample_to_patient(sample_id), **row})
        if i % 50 == 0:
            print(f"scored {i}/{len(sample_ids)} samples", flush=True)

    network_df = pd.DataFrame(network_rows)
    region_df = pd.DataFrame(region_rows)
    network_df.to_csv(args.outdir / "tcga_gbm_lgg_sample_network_tracing.csv", index=False, encoding="utf-8-sig")
    region_df.to_csv(args.outdir / "tcga_gbm_lgg_sample_region_tracing.csv", index=False, encoding="utf-8-sig")
    (args.outdir / "tcga_gbm_lgg_sample_tracing_meta.json").write_text(
        json.dumps(meta, ensure_ascii=False, indent=2), encoding="utf-8"
    )

    summary_rows: list[dict[str, Any]] = []
    for sample_id in sample_ids:
        patient = sample_to_patient(sample_id)
        label_row = labels[labels["patient_barcode"].eq(patient)]
        label = label_row.iloc[0].to_dict() if len(label_row) else {}
        networks_top3 = top_values(network_df, "network_id", sample_id, 3)
        regions_top5 = top_values(region_df, "region_id", sample_id, 5)
        groups_top3 = top_values(region_df, "resolution_group", sample_id, 3) if "resolution_group" in region_df.columns else []
        predicted_broad_top3 = broad_candidates(networks_top3)
        truth_network = split_truths(label.get("tumor_network", label.get("true_network", "")))
        truth_broad = split_truths(label.get("tumor_broad_anatomy", label.get("true_broad_anatomy", "")))
        truth_lobe = split_truths(label.get("tumor_lobe", label.get("true_lobe", "")))
        truth_region = split_truths(label.get("tumor_region", label.get("true_region", "")))
        project_info = sample_project.get(sample_id, {})
        summary_rows.append(
            {
                "sample_id": sample_id,
                "patient_barcode": patient,
                "project_id": project_info.get("project_id", ""),
                "has_mri_label": bool(label),
                "tumor_lobe": label.get("tumor_lobe", ""),
                "tumor_broad_anatomy": label.get("tumor_broad_anatomy", ""),
                "tumor_network": label.get("tumor_network", ""),
                "network_top1": networks_top3[0] if networks_top3 else "",
                "network_top3": " | ".join(networks_top3),
                "predicted_broad_top1": predicted_broad_top3[0] if predicted_broad_top3 else "",
                "predicted_broad_top3": " | ".join(predicted_broad_top3[:3]),
                "region_top1": regions_top5[0] if regions_top5 else "",
                "region_top5": " | ".join(regions_top5),
                "region_group_top1": groups_top3[0] if groups_top3 else "",
                "region_group_top3": " | ".join(groups_top3),
                "network_top1_match": any_match(networks_top3[:1], truth_network),
                "network_top3_match": any_match(networks_top3[:3], truth_network),
                "broad_top1_match": any_match(predicted_broad_top3[:1], truth_broad or truth_lobe),
                "broad_top3_match": any_match(predicted_broad_top3[:3], truth_broad or truth_lobe),
                "exact_region_top1_match": any_match(regions_top5[:1], truth_region),
                "exact_region_top5_match": any_match(regions_top5[:5], truth_region),
            }
        )

    summary_df = pd.DataFrame(summary_rows)
    summary_df.to_csv(args.outdir / "tcga_gbm_lgg_sample_mri_label_tracing_summary.csv", index=False, encoding="utf-8-sig")
    labeled = [row for row in summary_rows if row["has_mri_label"]]
    metrics = {
        "n_expression_samples": int(len(sample_ids)),
        "n_unique_patients_expression": int(summary_df["patient_barcode"].nunique()),
        "n_labeled_samples": int(len(labeled)),
        "n_labeled_patients": int(summary_df.loc[summary_df["has_mri_label"], "patient_barcode"].nunique()) if len(summary_df) else 0,
        "network_top1": metric(labeled, "network_top1_match"),
        "network_top3": metric(labeled, "network_top3_match"),
        "broad_top1": metric(labeled, "broad_top1_match"),
        "broad_top3": metric(labeled, "broad_top3_match"),
        "exact_region_top1": metric(labeled, "exact_region_top1_match"),
        "exact_region_top5": metric(labeled, "exact_region_top5_match"),
        "missing_label_patients": sorted(set(summary_df.loc[~summary_df["has_mri_label"], "patient_barcode"].astype(str))),
        "label_file_used": str(args.labels) if args.labels.exists() else "",
    }
    (args.outdir / "tcga_gbm_lgg_sample_mri_label_evaluation_metrics.json").write_text(
        json.dumps(metrics, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(json.dumps(metrics, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
