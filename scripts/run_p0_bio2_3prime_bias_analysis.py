#!/usr/bin/env python
"""P0-BIO2: length/detection analysis and a transparent 3'-capture proxy.

The public GSE189919 matrices are gene-level, so they cannot identify read
position along a transcript.  The simulation therefore thins each observed
gene count with p=min(1, W/L), where W is a terminal capture window and L is
the Ensembl-canonical transcript length.  This is a sensitivity proxy, not a
mechanistic estimate of degradation.
"""
from __future__ import annotations

import argparse
import hashlib
import io
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import requests
from scipy.stats import spearmanr


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def get_panel(model_path: Path, orthology_path: Path) -> pd.DataFrame:
    with np.load(model_path, allow_pickle=False) as model:
        genes = pd.Series(model["genes"].astype(str), name="model_gene")
    orthology = pd.read_csv(orthology_path, sep="\t")
    mapped = (
        orthology.dropna(subset=["Human gene name"])
        .drop_duplicates("Gene stable ID")
        .set_index("Gene stable ID")["Human gene name"]
    )
    frame = pd.DataFrame({"model_gene": genes})
    frame["human_gene"] = frame["model_gene"].map(lambda value: mapped.get(value, value))
    return frame


def fetch_lengths(genes: list[str], cache_path: Path) -> pd.DataFrame:
    if cache_path.exists():
        return pd.read_csv(cache_path, sep="\t")
    query = f"""<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query virtualSchemaName="default" formatter="TSV" header="1" uniqueRows="1"
 count="" datasetConfigVersion="0.6">
 <Dataset name="hsapiens_gene_ensembl" interface="default">
  <Filter name="hgnc_symbol" value="{','.join(genes)}"/>
  <Attribute name="hgnc_symbol"/>
  <Attribute name="ensembl_transcript_id"/>
  <Attribute name="transcript_length"/>
  <Attribute name="transcript_is_canonical"/>
 </Dataset>
</Query>"""
    response = requests.get(
        "https://www.ensembl.org/biomart/martservice",
        params={"query": query},
        timeout=120,
    )
    response.raise_for_status()
    frame = pd.read_csv(io.StringIO(response.text), sep="\t")
    frame.to_csv(cache_path, sep="\t", index=False)
    return frame


def summarize(values: pd.Series) -> dict[str, float]:
    return {
        "median": float(values.median()),
        "q1": float(values.quantile(0.25)),
        "q3": float(values.quantile(0.75)),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--model", type=Path, required=True)
    parser.add_argument("--orthology", type=Path, required=True)
    parser.add_argument("--counts", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--draws", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=20260729)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    panel = get_panel(args.model, args.orthology)
    raw_lengths = fetch_lengths(
        panel["human_gene"].drop_duplicates().tolist(),
        args.outdir / "ensembl_biomart_panel_transcripts.tsv",
    )
    raw_lengths.columns = ["human_gene", "ensembl_transcript_id", "length_bp", "canonical"]
    canonical = (
        raw_lengths[raw_lengths["canonical"] == 1]
        .sort_values(["human_gene", "ensembl_transcript_id"])
        .drop_duplicates("human_gene")
    )
    panel = panel.merge(
        canonical[["human_gene", "ensembl_transcript_id", "length_bp"]],
        on="human_gene",
        how="left",
    )

    counts = pd.read_csv(args.counts, index_col="Geneid")
    counts.index = counts.index.astype(str)
    counts = counts.apply(pd.to_numeric, errors="raise")
    panel_counts = counts.reindex(panel["human_gene"]).set_axis(panel.index)
    panel["present_in_matrix"] = panel_counts.notna().any(axis=1)
    panel_counts = panel_counts.fillna(0).astype(np.int64)
    panel["detection_rate"] = (panel_counts > 0).mean(axis=1).to_numpy()
    analysable = panel["length_bp"].notna() & panel["present_in_matrix"]
    rho, p_value = spearmanr(
        panel.loc[analysable, "length_bp"],
        panel.loc[analysable, "detection_rate"],
    )
    panel["length_quartile"] = pd.NA
    panel.loc[panel["length_bp"].notna(), "length_quartile"] = pd.qcut(
        panel.loc[panel["length_bp"].notna(), "length_bp"],
        4,
        labels=["Q1 shortest", "Q2", "Q3", "Q4 longest"],
        duplicates="drop",
    ).astype(str)
    panel.to_csv(args.outdir / "panel_length_detection.csv", index=False)

    strata = (
        panel.loc[analysable]
        .groupby("length_quartile", observed=True)
        .agg(
            n_genes=("human_gene", "size"),
            median_length_bp=("length_bp", "median"),
            median_detection_rate=("detection_rate", "median"),
            mean_detection_rate=("detection_rate", "mean"),
        )
        .reset_index()
    )
    strata.to_csv(args.outdir / "length_quartile_summary.csv", index=False)

    rng = np.random.default_rng(args.seed)
    annotated = panel["length_bp"].notna()
    base = panel_counts.loc[annotated].to_numpy(dtype=np.int64)
    lengths = panel.loc[annotated, "length_bp"].to_numpy(dtype=float)
    rows: list[dict[str, float | int]] = []
    for window in (500, 1000, 2000):
        probability = np.minimum(1.0, window / lengths)[:, None]
        for draw in range(args.draws):
            thinned = rng.binomial(base, probability)
            per_sample = (thinned > 0).sum(axis=0)
            rows.append(
                {
                    "window_bp": window,
                    "draw": draw + 1,
                    "median_detected_genes": float(np.median(per_sample)),
                    "minimum_detected_genes": int(per_sample.min()),
                    "mean_detected_genes": float(per_sample.mean()),
                    "median_detection_fraction": float(np.median(per_sample / len(lengths))),
                }
            )
    simulation = pd.DataFrame(rows)
    simulation.to_csv(args.outdir / "terminal_window_simulation_draws.csv", index=False)
    sim_summary = (
        simulation.groupby("window_bp")
        .agg(
            median_detected_genes=("median_detected_genes", "median"),
            q025_detected_genes=("median_detected_genes", lambda x: x.quantile(0.025)),
            q975_detected_genes=("median_detected_genes", lambda x: x.quantile(0.975)),
            median_detection_fraction=("median_detection_fraction", "median"),
            minimum_detected_genes_across_draws=("minimum_detected_genes", "min"),
        )
        .reset_index()
    )
    sim_summary.to_csv(args.outdir / "terminal_window_simulation_summary.csv", index=False)

    observed_per_sample = (base > 0).sum(axis=0)
    length_summary = summarize(panel.loc[panel["length_bp"].notna(), "length_bp"])
    result = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "analysis_boundary": (
            "Gene-level GSE189919 matrices contain no within-transcript coverage. "
            "The terminal-window binomial thinning is a length-dependent 3'-capture "
            "proxy and is not a mechanistic degradation model."
        ),
        "panel_genes": int(len(panel)),
        "human_orthology_substitutions": int((panel.model_gene != panel.human_gene).sum()),
        "genes_with_canonical_length": int(panel.length_bp.notna().sum()),
        "genes_present_and_length_annotated": int(analysable.sum()),
        "canonical_length_bp": length_summary,
        "spearman_length_vs_detection": {
            "rho": float(rho),
            "p_value_two_sided": float(p_value),
            "n_genes": int(analysable.sum()),
        },
        "observed_annotated_panel_detection_per_sample": {
            **summarize(pd.Series(observed_per_sample)),
            "denominator": int(len(lengths)),
        },
        "simulation_draws_per_window": int(args.draws),
        "seed": int(args.seed),
        "windows_bp": [500, 1000, 2000],
        "input_sha256": {
            "model": sha256(args.model),
            "orthology": sha256(args.orthology),
            "counts": sha256(args.counts),
        },
    }
    (args.outdir / "summary.json").write_text(
        json.dumps(result, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(result, indent=2, ensure_ascii=False))
    print(sim_summary.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
