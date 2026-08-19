#!/usr/bin/env python
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
WORKSPACE = ROOT.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

GENE_TYPES = (
    WORKSPACE
    / "raw_datasets_v0.1.9_20260728"
    / "01_Bo2023"
    / "auxiliary_buildkit"
    / "04_expressed_genes_neocortex_plus_subcortical.cleaned_symbols.csv"
)
ORTHOLOGY = (
    WORKSPACE
    / "raw_datasets_v0.1.9_20260728"
    / "08_Ensembl_Orthology"
    / "ensembl_mfascicularis_hsapiens_homology.tsv"
)
TOP200 = ROOT / "data" / "models" / "bo2023_saleem_network_top200_model_genes.csv"
NETWORK_MODEL = ROOT / "data" / "models" / "bo2023_saleem_network_top200_model.npz"
TCGA = ROOT / "reproducibility" / "v4_p0_12_tcga_brats_ci_summary.csv"
OUT = ROOT / "reproducibility" / "p1_bio1_4"

MISSING = [
    "ENSMFAG00000042392",
    "ENSMFAG00000046456",
    "ENSMFAG00000031052",
    "ENSMFAG00000015786",
    "ENSMFAG00000046600",
    "ENSMFAG00000030665",
    "ENSMFAG00000018905",
    "ENSMFAG00000017288",
    "ENSMFAG00000021768",
    "ENSMFAG00000038166",
    "ENSMFAG00000033848",
    "ENSMFAG00000018969",
]


def write_json(name: str, payload: object) -> None:
    (OUT / name).write_text(
        json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )


def analyze_missing_top200() -> None:
    genes = pd.read_csv(TOP200)
    subset = genes[genes["gene_symbol"].isin(MISSING)].copy()
    with np.load(NETWORK_MODEL, allow_pickle=False) as payload:
        model_genes = payload["genes"].astype(str)
        networks = payload["networks"].astype(str)
        reference = payload["reference"].astype(float)
    gene_to_row = {gene: i for i, gene in enumerate(model_genes)}

    types = pd.read_csv(GENE_TYPES).set_index("Gene.stable.ID")
    orth = pd.read_csv(ORTHOLOGY, sep="\t").set_index("Gene stable ID")
    rows = []
    for record in subset.itertuples(index=False):
        gene = str(record.gene_symbol)
        values = reference[gene_to_row[gene]]
        order = np.argsort(values)[::-1]
        ortho = orth.loc[gene] if gene in orth.index else None
        human_id = "" if ortho is None or pd.isna(ortho["Human gene stable ID"]) else str(ortho["Human gene stable ID"])
        human_name = "" if ortho is None or pd.isna(ortho["Human gene name"]) else str(ortho["Human gene name"])
        confidence = "" if ortho is None or pd.isna(ortho["Human orthology confidence [0 low, 1 high]"]) else float(ortho["Human orthology confidence [0 low, 1 high]"])
        rows.append(
            {
                "macaque_gene_id": gene,
                "top200_rank": int(genes.index[genes["gene_symbol"].eq(gene)][0]) + 1,
                "fisher_score": float(record.fisher_score),
                "gene_type": str(types.loc[gene, "Gene.type_ensembl"]),
                "highest_centroid_network": networks[order[0]],
                "second_centroid_network": networks[order[1]],
                "centroid_difference_top1_minus_top2": float(values[order[0]] - values[order[1]]),
                "centroid_range": float(values.max() - values.min()),
                "low_confidence_human_id": human_id,
                "low_confidence_human_name": human_name,
                "orthology_confidence": confidence,
            }
        )
    result = pd.DataFrame(rows).sort_values("top200_rank")
    result.to_csv(OUT / "bio2_missing_top200_network_and_annotation.csv", index=False)
    write_json(
        "bio2_missing_top200_summary.json",
        {
            "n_missing": len(result),
            "gene_type_counts": result["gene_type"].value_counts().to_dict(),
            "highest_centroid_network_counts": result["highest_centroid_network"].value_counts().to_dict(),
            "n_with_any_local_human_orthology_record": int(result["low_confidence_human_id"].ne("").sum()),
            "interpretation": (
                "Highest-centroid Network is a reference-expression attribution, not a functional "
                "annotation. With only one low-confidence human orthology record, pathway enrichment "
                "would be underidentified and was not performed."
            ),
        },
    )


def analyze_edema_endpoint() -> None:
    data = pd.read_csv(TCGA, comment="#")
    data = data[
        data["region_type"].isin(["center", "core", "edema", "whole_tumor"])
        & data["level"].isin(["network", "broad"])
        & data["top_k"].eq("top3")
        & data["variant"].eq("strict")
    ].copy()
    data["accuracy"] = pd.to_numeric(data["accuracy"], errors="coerce")
    data = data.dropna(subset=["accuracy"])
    data.to_csv(OUT / "bio4_truth_basis_top3_sensitivity.csv", index=False)
    pivot = data.pivot(index="region_type", columns="level", values="accuracy")
    patients_by_basis = data.groupby("region_type")["n_patients"].agg(
        lambda values: sorted({int(value) for value in values.dropna()})
    )
    if any(len(values) != 1 for values in patients_by_basis):
        raise ValueError("TCGA/BraTS source reports inconsistent denominators within a truth basis")
    non_edema_counts = {
        patients_by_basis[basis][0]
        for basis in ("center", "core", "whole_tumor")
    }
    if len(non_edema_counts) != 1:
        raise ValueError("TCGA/BraTS source reports inconsistent non-edema case counts")
    write_json(
        "bio4_edema_protocol_and_sensitivity.json",
        {
            "segmentation_semantics": {
                "whole_tumor": "segmentation label > 0",
                "core": "segmentation labels 1 or 4",
                "edema": "segmentation label 2",
                "center": "geometric/label center summary from the whole-tumour mask",
            },
            "n_patients": {
                "center_core_whole_tumor": non_edema_counts.pop(),
                "edema": patients_by_basis["edema"][0],
            },
            "strict_top3_accuracy": {
                basis: {level: float(pivot.loc[basis, level]) for level in pivot.columns}
                for basis in pivot.index
            },
            "range_across_truth_bases_percentage_points": {
                level: float((pivot[level].max() - pivot[level].min()) * 100)
                for level in pivot.columns
            },
            "sampling_boundary": (
                "TCGA RNA-seq is bulk primary-tumour material. No spatial record links the aliquot "
                "to BraTS label-2 edema voxels, tumour core, enhancing tumour or a matched MRI "
                "coordinate. Edema is an MRI-derived anatomical truth basis, not a transcriptomic "
                "compartment sampled for RNA-seq."
            ),
        },
    )


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    analyze_missing_top200()
    analyze_edema_endpoint()
    print(OUT)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
