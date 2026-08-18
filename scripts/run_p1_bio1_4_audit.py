#!/usr/bin/env python
from __future__ import annotations

import json
import sys
from functools import lru_cache
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
WORKSPACE = ROOT.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

HUANG = (
    WORKSPACE
    / "raw_datasets_v0.1.9_20260728"
    / "06_Huang2025_PMC12041490"
    / "41698_2025_909_MOESM2_ESM.csv"
)
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
PROJECTOR = (
    ROOT
    / "reproducibility"
    / "p0_bio3_projector"
    / "bo2023_reference_projector_linear_full.npz"
)
TCGA = ROOT / "reproducibility" / "v4_p0_12_tcga_brats_ci_summary.csv"
OUT = ROOT / "reproducibility" / "p1_bio1_4"

OMPFC = "Orbitomedial Prefrontal Cortex (OMPFC)"
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
PLATELET = ["PF4", "PPBP", "RGS18", "GP9", "ITGA2B", "TUBB1", "SELP", "NRGN"]
EV = ["CD9", "CD63", "CD81", "TSG101", "PDCD6IP", "SDCBP", "FLOT1", "FLOT2"]
FRACTIONS = [0.0, 0.10, 0.25, 0.50, 0.75, 0.90, 1.0]


def write_json(name: str, payload: object) -> None:
    (OUT / name).write_text(
        json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )


@lru_cache(maxsize=1)
def load_projection_arrays():
    with np.load(PROJECTOR, allow_pickle=False) as payload:
        return tuple(
            payload[name].copy()
            for name in ("genes", "slope", "intercept", "clip_low", "clip_high")
        )


@lru_cache(maxsize=1)
def load_network_arrays():
    with np.load(NETWORK_MODEL, allow_pickle=False) as payload:
        return tuple(payload[name].copy() for name in ("genes", "networks", "reference"))


def trace(values: pd.Series) -> dict[str, object]:
    pgenes, slope, intercept, clip_low, clip_high = load_projection_arrays()
    pgenes = pgenes.astype(str)
    projected = slope.astype(float) * values.reindex(pgenes).fillna(0).to_numpy(float) + intercept.astype(float)
    projected = np.clip(projected, clip_low.astype(float), clip_high.astype(float))
    projected_series = pd.Series(projected, index=pgenes)
    genes, networks, reference = load_network_arrays()
    genes, networks, reference = genes.astype(str), networks.astype(str), reference.astype(float)
    vector = projected_series.reindex(genes).fillna(0).to_numpy(float)
    vector_centered = vector - vector.mean()
    ref_centered = reference - reference.mean(axis=0, keepdims=True)
    denominator = np.sqrt(np.square(ref_centered).sum(axis=0) * np.square(vector_centered).sum())
    score_array = np.divide(
        (ref_centered * vector_centered[:, None]).sum(axis=0),
        denominator,
        out=np.zeros(reference.shape[1], dtype=float),
        where=denominator > 0,
    )
    order = np.argsort(score_array)[::-1]
    ids = networks[order].tolist()
    scores = {str(networks[i]): float(score_array[i]) for i in range(len(networks))}
    return {
        "top1": ids[0],
        "top3": " | ".join(ids[:3]),
        "ompfc_rank": ids.index(OMPFC) + 1,
        "ompfc_score": scores[OMPFC],
        "top1_margin": float(score_array[order[0]] - score_array[order[1]]),
    }


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


def z_composite(expression: pd.DataFrame, genes: list[str]) -> tuple[pd.Series, list[str]]:
    present = [gene for gene in genes if gene in expression.index]
    matrix = expression.loc[present].T.astype(float)
    sd = matrix.std(axis=0, ddof=0).replace(0, np.nan)
    z = (matrix - matrix.mean(axis=0)) / sd
    return z.mean(axis=1, skipna=True), present


def read_huang_expression() -> tuple[pd.DataFrame, pd.DataFrame]:
    header = pd.read_csv(HUANG, nrows=0).columns.astype(str).tolist()
    sample_ids = header[1:]
    rows = []
    for sample_id in sample_ids:
        prefix, specimen_number = sample_id.split("_", 1)
        specimen = "CSF" if specimen_number.startswith("CSF") else "plasma"
        number = int(specimen_number.replace("CSF", "").replace("plasma", ""))
        disease = {"GLI": "glioma", "MEN": "meningioma", "NOR": "control"}[prefix]
        rows.append(
            {
                "sample_id": sample_id,
                "disease": disease,
                "tumor_status": "control" if disease == "control" else "tumor",
                "specimen": specimen,
                "patient_key": f"{prefix}_{number:02d}",
            }
        )
    with np.load(PROJECTOR, allow_pickle=False) as payload:
        required = set(payload["genes"].astype(str))
    with np.load(NETWORK_MODEL, allow_pickle=False) as payload:
        required.update(payload["genes"].astype(str))
    required.update(PLATELET + EV)
    raw = pd.read_csv(HUANG, low_memory=False).rename(columns={header[0]: "gene_symbol"})
    selected = raw[raw["gene_symbol"].astype(str).isin(required)].copy()
    selected["gene_symbol"] = selected["gene_symbol"].astype(str)
    values = selected[sample_ids].to_numpy(dtype=float) * np.log(2.0)
    selected.loc[:, sample_ids] = values
    return selected.set_index("gene_symbol").astype(np.float32), pd.DataFrame(rows)


def analyze_huang_and_admixture() -> None:
    expression, metadata = read_huang_expression()
    meta = metadata.set_index("sample_id")
    platelet, platelet_present = z_composite(expression, PLATELET)
    ev, ev_present = z_composite(expression, EV)

    pairs = []
    for row in metadata.itertuples(index=False):
        if row.specimen != "CSF":
            continue
        plasma_id = row.sample_id.replace("_CSF", "_plasma")
        if plasma_id not in expression.columns:
            continue
        csf_linear = np.expm1(expression[row.sample_id].astype(float))
        plasma_linear = np.expm1(expression[plasma_id].astype(float))
        baseline_top1 = None
        for fraction in FRACTIONS:
            mixed = np.log1p((1 - fraction) * csf_linear + fraction * plasma_linear)
            prediction = trace(mixed)
            if fraction == 0:
                baseline_top1 = prediction["top1"]
            pairs.append(
                {
                    "patient_key": row.patient_key,
                    "disease": row.disease,
                    "plasma_fraction": fraction,
                    "baseline_csf_top1": baseline_top1,
                    **prediction,
                }
            )
    detail = pd.DataFrame(pairs)
    detail.to_csv(OUT / "bio3_matched_csf_plasma_admixture_detail.csv", index=False)
    summary_rows = []
    for fraction, group in detail.groupby("plasma_fraction", sort=True):
        summary_rows.append(
            {
                "plasma_fraction": fraction,
                "n_pairs": len(group),
                "ompfc_top1_fraction": float(group["top1"].eq(OMPFC).mean()),
                "ompfc_top3_fraction": float(group["top3"].str.contains(OMPFC, regex=False).mean()),
                "baseline_top1_retained_fraction": float(group["top1"].eq(group["baseline_csf_top1"]).mean()),
                "median_ompfc_rank": float(group["ompfc_rank"].median()),
                "median_top1_margin": float(group["top1_margin"].median()),
            }
        )
    summary = pd.DataFrame(summary_rows)
    summary.to_csv(OUT / "bio3_matched_csf_plasma_admixture_summary.csv", index=False)

    baseline = detail[detail["plasma_fraction"].isin([0.0, 1.0])].copy()
    baseline["sample_id"] = baseline.apply(
        lambda x: f"{x.patient_key.split('_')[0]}_{'CSF' if x.plasma_fraction == 0 else 'plasma'}{int(x.patient_key.split('_')[1])}",
        axis=1,
    )
    baseline = baseline.set_index("sample_id")
    marker_rows = []
    for specimen in ["CSF", "plasma", "all"]:
        ids = meta.index if specimen == "all" else meta.index[meta["specimen"].eq(specimen)]
        ids = [sample for sample in ids if sample in baseline.index]
        for name, composite in [("platelet", platelet), ("EV", ev)]:
            x = composite.reindex(ids).astype(float)
            y = baseline.reindex(ids)["ompfc_score"].astype(float)
            valid = x.notna() & y.notna()
            rho = x[valid].rank(method="average").corr(y[valid].rank(method="average"))
            marker_rows.append(
                {
                    "specimen": specimen,
                    "marker_composite": name,
                    "n": int(np.isfinite(x).sum()),
                    "spearman_rho_vs_ompfc_score": float(rho),
                    "p_value": "",
                }
            )
    pd.DataFrame(marker_rows).to_csv(OUT / "bio1_marker_ompfc_association.csv", index=False)
    write_json(
        "bio1_bio3_huang_summary.json",
        {
            "platelet_markers_requested": PLATELET,
            "platelet_markers_present_in_projector_input": platelet_present,
            "ev_markers_requested": EV,
            "ev_markers_present_in_projector_input": ev_present,
            "n_matched_csf_plasma_pairs": int(detail["patient_key"].nunique()),
            "mixture_definition": (
                "Linear-RPM mixture of each patient's CSF and plasma profile, followed by log1p "
                "and the frozen Network route. Plasma is a systemic/blood-background proxy, not "
                "a purified non-brain tissue reference."
            ),
            "limitations": (
                "Marker correlations and plasma admixture are diagnostic stress tests. They cannot "
                "separate platelet contamination, EV enrichment, class-prior effects and other "
                "biofluid domain shifts, and they do not validate anatomical localization."
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
    analyze_huang_and_admixture()
    analyze_edema_endpoint()
    print(OUT)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
