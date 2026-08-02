#!/usr/bin/env python3
"""Reproducible descriptive audits for round-4 P1-NEURO1--4."""
from __future__ import annotations

import csv
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
groups_path = ROOT / "code/supplementary/p1_resolution_groups.csv"
macro_path = ROOT / "code/reproducibility/macro_f1_class_data.json"
out_dir = ROOT / "manuscript/calculations"
out_dir.mkdir(parents=True, exist_ok=True)

qualified = {
    "Lateral Prefrontal Cortex::10o + 46d + 46v": (
        "qualified", "frontopolar, dorsolateral and ventrolateral prefrontal fields; transcriptomic grouping, not one classical subdivision"
    ),
    "Lateral Prefrontal Cortex::12l + 12r + 45 + 8A": (
        "qualified", "orbital/ventrolateral and dorsal area-8 fields span distinct connectivity systems"
    ),
    "Operculum/Insula::G + Id + Ig + SII": (
        "qualified", "granular/dysgranular insula plus secondary somatosensory operculum; operculo-insular affinity is plausible but not identity"
    ),
    "Orbitomedial Prefrontal Cortex (OMPFC)::13a + 13b + 25": (
        "qualified", "orbitofrontal area 13 and subgenual area 25 require a broad limbic-orbitomedial interpretation"
    ),
    "Parietal, and Parieto-occipital region::5 + LIPv + VIP": (
        "qualified", "area 5 somatosensory association plus intraparietal visuospatial fields; cross-modal convergence rather than one cytoarchitectonic unit"
    ),
    "Temporal::CL + ML + Tpt": (
        "qualified", "auditory belt fields CL/ML plus temporoparietal association area Tpt; broad auditory-association grouping"
    ),
}

rows = []
with groups_path.open(encoding="utf-8-sig", newline="") as f:
    for row in csv.DictReader(f):
        tier, rationale = qualified.get(
            row["group_id"],
            ("clear", "members share a recognized local field family or functional-anatomical system"),
        )
        rows.append({**row, "anatomical_audit_tier": tier, "audit_rationale": rationale})

with (out_dir / "P1_NEURO1_resolution_group_anatomy_audit.csv").open(
    "w", encoding="utf-8-sig", newline=""
) as f:
    w = csv.DictWriter(f, fieldnames=rows[0].keys())
    w.writeheader()
    w.writerows(rows)

macro = json.loads(macro_path.read_text(encoding="utf-8"))
metrics = {(r["endpoint"], r["class"]): r for r in macro["data"]}

def metric(endpoint, cls):
    r = metrics[(endpoint, cls)]
    return {k: float(r[k]) if k != "class" and k != "endpoint" else r[k]
            for k in r}

temporal_loso = metric("LOSO_Network", "Temporal")
parietal_loso = metric("LOSO_Network", "Parietal, and Parieto-occipital region")
temporal_tp = round(temporal_loso["n"] * temporal_loso["recall"])
temporal_errors = int(temporal_loso["n"] - temporal_tp)
temporal_to_oi = 47
temporal_to_parietal = 46
temporal_to_visual = 13

parietal_tp = round(parietal_loso["n"] * parietal_loso["recall"])
parietal_predicted = round(parietal_tp / parietal_loso["precision"])
parietal_fp = parietal_predicted - parietal_tp
known_temporal_visual_fp = temporal_to_parietal + 13  # Visual/dorsal STS -> Parietal

payload = {
    "neuro1_resolution_groups": {
        "total": len(rows),
        "clear": sum(r["anatomical_audit_tier"] == "clear" for r in rows),
        "qualified": sum(r["anatomical_audit_tier"] == "qualified" for r in rows),
        "interpretation": "qualified groups are operational transcriptomic resolution groups, not classical anatomical units",
    },
    "neuro2_hippocampal": {
        "network_regions": 1,
        "training_samples": 8,
        "tiers": "Network -> resolution-group=exact-region",
        "independent_fine_tier_evidence": False,
        "tool_annotation": "single-region Network; group=exact",
    },
    "neuro3_temporal_loso": {
        "true_samples": int(temporal_loso["n"]),
        "true_positive": temporal_tp,
        "false_negative": temporal_errors,
        "to_operculum_insula": temporal_to_oi,
        "to_parietal": temporal_to_parietal,
        "to_visual_dorsal_sts": temporal_to_visual,
        "oi_plus_parietal_share_of_temporal_errors": (temporal_to_oi + temporal_to_parietal) / temporal_errors,
        "three_destinations_share_of_temporal_errors": (temporal_to_oi + temporal_to_parietal + temporal_to_visual) / temporal_errors,
    },
    "neuro4_parietal_loso": {
        "true_samples": int(parietal_loso["n"]),
        "recall": parietal_loso["recall"],
        "precision": parietal_loso["precision"],
        "true_positive": parietal_tp,
        "predicted_positive_approx": parietal_predicted,
        "false_positive_approx": parietal_fp,
        "known_temporal_plus_visual_false_positives": known_temporal_visual_fp,
        "known_share_of_false_positives": known_temporal_visual_fp / parietal_fp,
        "interpretation": "capture of temporal polymodal and dorsal-visual samples is consistent with a broad convergence profile, but no connectivity model was tested",
    },
}
(out_dir / "P1_NEURO1-4_audit.json").write_text(
    json.dumps(payload, indent=2), encoding="utf-8"
)
print(json.dumps(payload, indent=2))
