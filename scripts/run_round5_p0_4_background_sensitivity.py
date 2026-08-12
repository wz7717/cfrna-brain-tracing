#!/usr/bin/env python3
"""Round-5 P0-4: same-version GO:BP/KEGG background sensitivity audit."""
from __future__ import annotations

import hashlib
import json
import math
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import requests
from scipy.stats import spearmanr


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "reproducibility/round5_analysis/p0_4_background_sensitivity"
PANEL_PATH = ROOT / "data/models/bo2023_saleem_network_top200_model_genes.csv"
PROJECTOR_PATH = ROOT / "data/models/bo2023_reference_projector_linear_full.npz"
DESEQ2_PATH = ROOT / (
    "reproducibility/round5_analysis/p0_3_deseq2_marker_audit/outputs/"
    "primary_pseudobulk/pseudobulk_deseq2_network_lrt_all_genes.csv"
)
FROZEN_MANIFEST = ROOT / "reproducibility/independent_enrichment/independent_enrichment_manifest.json"
API = "https://biit.cs.ut.ee/gprofiler/api/gost/profile/"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def symbols(values) -> set[str]:
    return {
        str(value).strip().upper()
        for value in values
        if pd.notna(value) and str(value).strip()
    }


def query_gprofiler(panel: set[str], background: set[str] | None, label: str) -> dict:
    payload = {
        "organism": "hsapiens",
        "query": sorted(panel),
        "sources": ["GO:BP", "KEGG"],
        "user_threshold": 0.05,
        "significance_threshold_method": "fdr",
        "no_evidences": False,
        "domain_scope": "annotated" if background is None else "custom",
        "all_results": True,
    }
    if background is not None:
        payload["background"] = sorted(background)
    response = requests.post(API, json=payload, timeout=300)
    response.raise_for_status()
    result = response.json()
    if "result" not in result or "meta" not in result:
        raise RuntimeError(f"Malformed g:Profiler response for {label}")
    return {"label": label, "request": payload, "response": result}


def version_fields(payload: dict) -> dict:
    meta = payload["response"].get("meta", {})
    return {
        "version": meta.get("version"),
        "timestamp": meta.get("timestamp"),
    }


def result_frame(payload: dict) -> pd.DataFrame:
    frame = pd.DataFrame(payload["response"].get("result", []))
    if frame.empty:
        raise RuntimeError(f"No g:Profiler results for {payload['label']}")
    frame.insert(0, "background", payload["label"])
    frame["significant_q05"] = frame["p_value"].astype(float) < 0.05
    return frame


def pairwise_metrics(a: pd.DataFrame, b: pd.DataFrame, a_name: str, b_name: str) -> list[dict]:
    rows = []
    for source in ["GO:BP", "KEGG"]:
        aa = a.loc[a.source.eq(source)].set_index("native")
        bb = b.loc[b.source.eq(source)].set_index("native")
        sa = set(aa.index[aa.significant_q05])
        sb = set(bb.index[bb.significant_q05])
        common = sorted(set(aa.index) & set(bb.index))
        if len(common) >= 2:
            x = -np.log10(np.clip(aa.loc[common, "p_value"].astype(float), 1e-300, 1))
            y = -np.log10(np.clip(bb.loc[common, "p_value"].astype(float), 1e-300, 1))
            rho = float(spearmanr(x, y).statistic)
        else:
            rho = math.nan
        rows.append({
            "source": source,
            "background_a": a_name,
            "background_b": b_name,
            "significant_a_n": len(sa),
            "significant_b_n": len(sb),
            "significant_intersection_n": len(sa & sb),
            "significant_union_n": len(sa | sb),
            "significant_jaccard": len(sa & sb) / len(sa | sb) if sa | sb else math.nan,
            "retention_of_a": len(sa & sb) / len(sa) if sa else math.nan,
            "retention_of_b": len(sa & sb) / len(sb) if sb else math.nan,
            "common_tested_terms_n": len(common),
            "spearman_neglog10_q": rho,
        })
    return rows


def representative_themes(frames: dict[str, pd.DataFrame]) -> pd.DataFrame:
    # Prespecified from the representative biological themes already discussed in the manuscript.
    patterns = {
        "cell-cell signaling": "cell-cell signaling",
        "nervous system development": "nervous system development",
        "generation of neurons": "generation of neurons",
        "potassium ion transport": "potassium ion transport",
        "glutamatergic synaptic transmission": "synaptic transmission, glutamatergic",
        "neuroactive ligand-receptor interaction": "neuroactive ligand-receptor interaction",
        "calcium signaling pathway": "calcium signaling pathway",
        "axon guidance": "axon guidance",
        "cholinergic synapse": "cholinergic synapse",
        "glutamatergic synapse": "glutamatergic synapse",
    }
    rows = []
    for theme, exact_name in patterns.items():
        for label, frame in frames.items():
            hit = frame.loc[frame["name"].astype(str).str.casefold().eq(exact_name.casefold())]
            rows.append({
                "theme": theme,
                "background": label,
                "native": None if hit.empty else hit.iloc[0]["native"],
                "source": None if hit.empty else hit.iloc[0]["source"],
                "q_value": math.nan if hit.empty else float(hit.iloc[0]["p_value"]),
                "significant_q05": False if hit.empty else bool(hit.iloc[0]["significant_q05"]),
            })
    return pd.DataFrame(rows)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    outputs = OUT / "outputs"
    outputs.mkdir(exist_ok=True)

    panel = symbols(pd.read_csv(PANEL_PATH)["gene_symbol"])
    model = symbols(np.load(PROJECTOR_PATH, allow_pickle=True)["genes"])
    de = pd.read_csv(DESEQ2_PATH)
    expression = symbols(de["gene_symbol"])
    if len(panel) != 200:
        raise RuntimeError(f"Expected 200 locked genes, found {len(panel)}")
    if not panel <= model or not panel <= expression:
        raise RuntimeError(
            f"Panel missing from background: model={sorted(panel-model)}, expression={sorted(panel-expression)}"
        )

    jobs = [
        query_gprofiler(panel, model, "model_space_21668"),
        query_gprofiler(panel, None, "gprofiler_annotated_domain"),
        query_gprofiler(panel, expression, "bo2023_pseudobulk_tested_expression"),
    ]
    versions = [version_fields(job) for job in jobs]
    if len({v["version"] for v in versions}) != 1:
        raise RuntimeError(f"g:Profiler version changed within run: {versions}")

    frames: dict[str, pd.DataFrame] = {}
    for job in jobs:
        label = job["label"]
        (outputs / f"gprofiler_{label}.json").write_text(
            json.dumps(job, indent=2), encoding="utf-8"
        )
        frames[label] = result_frame(job)
        frames[label].to_csv(outputs / f"gprofiler_{label}_all_terms.csv", index=False)

    labels = list(frames)
    pairwise = []
    for i in range(len(labels)):
        for j in range(i + 1, len(labels)):
            pairwise.extend(pairwise_metrics(frames[labels[i]], frames[labels[j]], labels[i], labels[j]))
    pd.DataFrame(pairwise).to_csv(outputs / "pairwise_background_stability.csv", index=False)
    representative_themes(frames).to_csv(outputs / "prespecified_theme_stability.csv", index=False)

    counts = []
    for label, frame in frames.items():
        for source in ["GO:BP", "KEGG"]:
            sub = frame.loc[frame.source.eq(source)]
            counts.append({
                "background": label,
                "source": source,
                "tested_terms_n": len(sub),
                "significant_q05_n": int(sub.significant_q05.sum()),
            })
    pd.DataFrame(counts).to_csv(outputs / "significant_term_counts.csv", index=False)

    sig_sets = {
        source: [
            set(frame.loc[frame.source.eq(source) & frame.significant_q05, "native"])
            for frame in frames.values()
        ]
        for source in ["GO:BP", "KEGG"]
    }
    three_way = []
    for source, sets_ in sig_sets.items():
        intersection = set.intersection(*sets_)
        union = set.union(*sets_)
        three_way.append({
            "source": source,
            "three_way_intersection_n": len(intersection),
            "three_way_union_n": len(union),
            "three_way_jaccard": len(intersection) / len(union) if union else math.nan,
            "model_terms_retained_in_both_sensitivities": len(intersection) / len(sets_[0]) if sets_[0] else math.nan,
        })
    pd.DataFrame(three_way).to_csv(outputs / "three_way_background_stability.csv", index=False)

    frozen = json.loads(FROZEN_MANIFEST.read_text(encoding="utf-8"))
    old_version = frozen.get("gprofiler_meta_primary", {}).get("version")
    current = {
        "version": versions[0]["version"],
        "request_timestamps": [v["timestamp"] for v in versions],
    }
    manifest = {
        "analysis": "Round-5 P0-4 GO:BP/KEGG background sensitivity",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "panel_n": len(panel),
        "backgrounds": {
            "model_space_21668": len(model),
            "gprofiler_annotated_domain": None,
            "bo2023_pseudobulk_tested_expression": len(expression),
        },
        "expression_background_rule": "unique nonblank uppercase gene_symbol among all 23,331 DESeq2-tested primary donor-by-Network pseudobulk rows; no significance filter",
        "gprofiler_current": current,
        "gprofiler_frozen_version_detected": old_version,
        "same_version_within_run": True,
        "request_rule": "all terms requested; q<0.05 applied to g:Profiler FDR-adjusted p_value; GO:BP and KEGG only",
        "input_sha256": {str(p.relative_to(ROOT)): sha256(p) for p in [PANEL_PATH, PROJECTOR_PATH, DESEQ2_PATH, FROZEN_MANIFEST]},
    }
    (OUT / "p0_4_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")

    files = sorted(p for p in OUT.rglob("*") if p.is_file())
    with (OUT / "SHA256SUMS.txt").open("w", encoding="utf-8", newline="\n") as handle:
        for path in files:
            if path.name != "SHA256SUMS.txt":
                handle.write(f"{sha256(path)}  {path.relative_to(OUT).as_posix()}\n")


if __name__ == "__main__":
    main()
