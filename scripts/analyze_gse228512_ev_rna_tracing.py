#!/usr/bin/env python
from __future__ import annotations

import argparse
import gzip
import json
import re
import sys
from collections import Counter
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import chi2_contingency, mannwhitneyu


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from core.network_tracing import trace_network_expression  # noqa: E402

DATA_DIR = ROOT / "data" / "external_validation" / "GSE228512"
OUTDIR = ROOT / "results" / "gse228512_ev_rna_tracing_20260613"
BASELINE_MODEL = ROOT / "data" / "models" / "bo2023_saleem_network_top200_model.npz"
ADAPTED_MODEL = ROOT / "data" / "models" / "bo2023_tcga_tumor_adapted_network_model.npz"
ORTHOLOGY = ROOT / "data" / "orthology" / "ensembl_mfascicularis_hsapiens_homology.tsv"

NETWORK_TO_PRIMARY_BROAD = {
    "Cingulate gyrus": "cingulate",
    "Frontal (agranular frontal motor areas)": "frontal",
    "Hippocampal formation": "medial_temporal",
    "Lateral Prefrontal Cortex": "frontal",
    "Occipital/Temporal": "occipital_temporal",
    "Operculum/Insula": "insula_operculum",
    "Orbitomedial Prefrontal Cortex (OMPFC)": "frontal",
    "Parietal, and Parieto-occipital region": "parietal_occipital",
    "Subcortical": "subcortical",
    "Temporal": "temporal",
}


def read_soft_metadata(path: Path) -> pd.DataFrame:
    raw = gzip.open(path, "rt", encoding="utf-8", errors="replace").read()
    rows = []
    for block in raw.split("^SAMPLE = ")[1:]:
        accession = re.search(r"!Sample_geo_accession = (.*)", block)
        title = re.search(r"!Sample_title = (.*)", block)
        platform = re.search(r"!Sample_platform_id = (.*)", block)
        characteristics = re.findall(r"!Sample_characteristics_ch1 = (.*)", block)
        parsed = {}
        for item in characteristics:
            if ":" in item:
                key, value = item.split(":", 1)
                parsed[key.strip().lower()] = value.strip()
        rows.append(
            {
                "geo_accession": accession.group(1) if accession else "",
                "geo_title": title.group(1) if title else "",
                "platform_id": platform.group(1) if platform else "",
                "group": "GBM" if parsed.get("disease") == "Glioblastoma" else "Healthy",
            }
        )
    return pd.DataFrame(rows)


def sample_key(value: str) -> str:
    text = re.sub(r"[^a-z0-9]", "", str(value).lower())
    return text.removeprefix("control")


def load_counts_and_metadata(data_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    matrices = []
    sample_platform = {}
    for name, platform in (
        ("GSE228512_hiseq_counts.txt.gz", "HiSeq"),
        ("GSE228512_novaseq_counts.txt.gz", "NovaSeq"),
    ):
        frame = pd.read_csv(data_dir / name, sep="\t", index_col=0)
        frame.index = frame.index.astype(str).str.replace(r"\.\d+$", "", regex=True)
        frame = frame.apply(pd.to_numeric, errors="coerce").fillna(0.0)
        matrices.append(frame)
        sample_platform.update({str(column): platform for column in frame.columns})
    counts = pd.concat(matrices, axis=1, verify_integrity=True)

    soft = read_soft_metadata(data_dir / "GSE228512_family.soft.gz")
    soft_by_key = {sample_key(row.geo_title): row for row in soft.itertuples(index=False)}
    rows = []
    for sample in counts.columns.astype(str):
        match = soft_by_key.get(sample_key(sample))
        if match is None:
            raise ValueError(f"No SOFT metadata match for matrix sample {sample}")
        rows.append(
            {
                "sample_id": sample,
                "geo_accession": match.geo_accession,
                "geo_title": match.geo_title,
                "group": match.group,
                "platform": sample_platform[sample],
                "platform_id": match.platform_id,
            }
        )
    metadata = pd.DataFrame(rows)
    return counts, metadata


def ensembl_symbol_map(path: Path) -> pd.Series:
    table = pd.read_csv(path, sep="\t", dtype=str)
    table = table[["Human gene stable ID", "Human gene name"]].dropna()
    table["Human gene stable ID"] = table["Human gene stable ID"].str.replace(
        r"\.\d+$", "", regex=True
    )
    table["Human gene name"] = table["Human gene name"].str.strip()
    table = table[
        table["Human gene stable ID"].str.startswith("ENSG")
        & table["Human gene name"].str.match(r"^[A-Za-z][A-Za-z0-9.-]+$")
    ]
    return table.drop_duplicates("Human gene stable ID").set_index("Human gene stable ID")[
        "Human gene name"
    ]


def counts_to_symbol_cpm(counts: pd.DataFrame, mapping: pd.Series) -> pd.DataFrame:
    symbols = counts.index.to_series().map(mapping)
    keep = symbols.notna()
    symbol_counts = counts.loc[keep].groupby(symbols.loc[keep], sort=True).sum()
    library_sizes = symbol_counts.sum(axis=0)
    return symbol_counts.divide(library_sizes.where(library_sizes > 0), axis=1) * 1_000_000.0


def corr_scores(reference: np.ndarray, samples: np.ndarray) -> np.ndarray:
    ref0 = reference - reference.mean(axis=0, keepdims=True)
    x0 = samples - samples.mean(axis=0, keepdims=True)
    numerator = ref0.T @ x0
    denominator = np.sqrt(
        np.square(ref0).sum(axis=0)[:, None] * np.square(x0).sum(axis=0)[None, :] + 1e-12
    )
    return np.nan_to_num(numerator / denominator)


def softmax(scores: np.ndarray) -> np.ndarray:
    shifted = scores - scores.max(axis=0, keepdims=True)
    values = np.exp(shifted)
    return values / values.sum(axis=0, keepdims=True)


def quantile_map(values: np.ndarray, target: np.ndarray) -> np.ndarray:
    mapped = np.empty_like(values, dtype=float)
    target = np.sort(np.asarray(target, dtype=float))
    for column in range(values.shape[1]):
        order = np.argsort(values[:, column], kind="mergesort")
        mapped[order, column] = target
    return mapped


def model_scores(
    expression: pd.DataFrame, model_path: Path, transform: str
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    model = np.load(model_path, allow_pickle=False)
    genes = model["genes"].astype(str)
    networks = model["networks"].astype(str)
    present = np.asarray([gene in expression.index for gene in genes], dtype=bool)
    values = expression.reindex(genes).fillna(0.0).to_numpy(dtype=float)
    if transform == "logcpm":
        values = np.log1p(values)
    elif transform == "harmonized":
        values = quantile_map(np.log1p(values), model["target_quantiles"].astype(float))
    elif transform != "cpm":
        raise ValueError(transform)
    scores = corr_scores(model["reference"].astype(float), values)
    return scores, networks, genes, present


def add_route(
    output: pd.DataFrame,
    route: str,
    scores: np.ndarray,
    networks: np.ndarray,
    ood_threshold: float | None = None,
) -> pd.DataFrame:
    probabilities = softmax(scores)
    order = np.argsort(scores, axis=0)[::-1]
    n_classes = scores.shape[0]
    for index, network in enumerate(networks):
        output[f"{route}_score_{index + 1:02d}"] = scores[index]
        output[f"{route}_prob_{index + 1:02d}"] = probabilities[index]
    output[f"{route}_top1_network"] = [networks[int(order[0, i])] for i in range(len(output))]
    output[f"{route}_top3_network"] = [
        " | ".join(networks[int(j)] for j in order[:3, i]) for i in range(len(output))
    ]
    output[f"{route}_top1_broad"] = output[f"{route}_top1_network"].map(
        NETWORK_TO_PRIMARY_BROAD
    )
    output[f"{route}_top3_broad"] = [
        " | ".join(
            dict.fromkeys(NETWORK_TO_PRIMARY_BROAD[networks[int(j)]] for j in order[:3, i])
        )
        for i in range(len(output))
    ]
    sorted_scores = np.sort(scores, axis=0)
    sorted_probabilities = np.sort(probabilities, axis=0)
    output[f"{route}_max_score"] = scores.max(axis=0)
    output[f"{route}_raw_margin"] = sorted_scores[-1] - sorted_scores[-2]
    output[f"{route}_top1_probability"] = probabilities.max(axis=0)
    output[f"{route}_probability_margin"] = (
        sorted_probabilities[-1] - sorted_probabilities[-2]
    )
    entropy = -(probabilities * np.log(probabilities + 1e-12)).sum(axis=0)
    output[f"{route}_normalized_entropy"] = entropy / np.log(n_classes)
    if ood_threshold is not None:
        output[f"{route}_ood_threshold"] = ood_threshold
        output[f"{route}_ood_accepted"] = output[f"{route}_max_score"] >= ood_threshold
    return output


def stratified_permutation_p(
    values: np.ndarray,
    groups: np.ndarray,
    strata: np.ndarray,
    statistic,
    iterations: int,
    seed: int,
) -> float:
    observed = abs(float(statistic(values, groups)))
    rng = np.random.default_rng(seed)
    exceed = 0
    for _ in range(iterations):
        permuted = groups.copy()
        for stratum in np.unique(strata):
            idx = np.flatnonzero(strata == stratum)
            permuted[idx] = rng.permutation(permuted[idx])
        exceed += abs(float(statistic(values, permuted))) >= observed
    return (exceed + 1) / (iterations + 1)


def bh_adjust(values: pd.Series) -> pd.Series:
    p = pd.to_numeric(values, errors="coerce").to_numpy(dtype=float)
    order = np.argsort(p)
    adjusted = np.empty_like(p)
    running = 1.0
    n = len(p)
    for rank_index in range(n - 1, -1, -1):
        original_index = order[rank_index]
        rank = rank_index + 1
        running = min(running, p[original_index] * n / rank)
        adjusted[original_index] = min(running, 1.0)
    return pd.Series(adjusted, index=values.index)


def summarize_groups(
    output: pd.DataFrame, routes: list[str], iterations: int
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    continuous_rows = []
    distribution_rows = []
    collapse_rows = []
    group_values = output["group"].to_numpy()
    strata = output["platform"].to_numpy()
    for route in routes:
        for metric in (
            "max_score",
            "raw_margin",
            "top1_probability",
            "probability_margin",
            "normalized_entropy",
        ):
            column = f"{route}_{metric}"
            gbm = output.loc[output.group.eq("GBM"), column].to_numpy(dtype=float)
            healthy = output.loc[output.group.eq("Healthy"), column].to_numpy(dtype=float)
            _, p = mannwhitneyu(gbm, healthy, alternative="two-sided")

            def difference(values: np.ndarray, labels: np.ndarray) -> float:
                return float(values[labels == "GBM"].mean() - values[labels == "Healthy"].mean())

            permutation_p = stratified_permutation_p(
                output[column].to_numpy(dtype=float),
                group_values.copy(),
                strata,
                difference,
                iterations,
                seed=1731,
            )
            continuous_rows.append(
                {
                    "route": route,
                    "metric": metric,
                    "gbm_median": float(np.median(gbm)),
                    "healthy_median": float(np.median(healthy)),
                    "gbm_mean": float(np.mean(gbm)),
                    "healthy_mean": float(np.mean(healthy)),
                    "mean_difference_gbm_minus_healthy": float(np.mean(gbm) - np.mean(healthy)),
                    "mann_whitney_p": float(p),
                    "platform_stratified_permutation_p": float(permutation_p),
                }
            )

        prediction = f"{route}_top1_network"
        table = pd.crosstab(output["group"], output[prediction])
        chi2, p, _, _ = chi2_contingency(table)

        def distribution_stat(_: np.ndarray, labels: np.ndarray) -> float:
            perm_table = pd.crosstab(labels, output[prediction])
            return float(chi2_contingency(perm_table)[0])

        perm_p = stratified_permutation_p(
            np.arange(len(output)),
            group_values.copy(),
            strata,
            distribution_stat,
            iterations,
            seed=8821,
        )
        distribution_rows.append(
            {
                "route": route,
                "chi_square": float(chi2),
                "degrees_of_freedom": int((table.shape[0] - 1) * (table.shape[1] - 1)),
                "chi_square_p": float(p),
                "platform_stratified_permutation_p": float(perm_p),
                "active_networks": int(output[prediction].nunique()),
            }
        )
        for group, sub in output.groupby("group"):
            counts = sub[prediction].value_counts()
            proportions = counts / counts.sum()
            hhi = float(np.square(proportions).sum())
            collapse_rows.append(
                {
                    "route": route,
                    "group": group,
                    "n": int(len(sub)),
                    "active_networks": int(len(counts)),
                    "dominant_network": str(counts.index[0]),
                    "dominant_count": int(counts.iloc[0]),
                    "dominant_fraction": float(proportions.iloc[0]),
                    "hhi": hhi,
                    "effective_network_count": float(1.0 / hhi),
                    "collapse_flag_dominant_fraction_ge_0p70": bool(
                        proportions.iloc[0] >= 0.70
                    ),
                }
            )
    return (
        pd.DataFrame(continuous_rows),
        pd.DataFrame(distribution_rows),
        pd.DataFrame(collapse_rows),
    )


def summarize_prediction_only(
    output: pd.DataFrame, route: str, iterations: int
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    prediction = f"{route}_top1_network"
    groups = output["group"].to_numpy()
    strata = output["platform"].to_numpy()
    table = pd.crosstab(output["group"], output[prediction])
    chi2, p, _, _ = chi2_contingency(table)

    def distribution_stat(_: np.ndarray, labels: np.ndarray) -> float:
        return float(chi2_contingency(pd.crosstab(labels, output[prediction]))[0])

    permutation_p = stratified_permutation_p(
        np.arange(len(output)),
        groups.copy(),
        strata,
        distribution_stat,
        iterations,
        seed=9917,
    )
    distribution = {
        "route": route,
        "chi_square": float(chi2),
        "degrees_of_freedom": int((table.shape[0] - 1) * (table.shape[1] - 1)),
        "chi_square_p": float(p),
        "platform_stratified_permutation_p": float(permutation_p),
        "active_networks": int(output[prediction].nunique()),
    }
    collapse = []
    for group, sub in output.groupby("group"):
        counts = sub[prediction].value_counts()
        proportions = counts / counts.sum()
        hhi = float(np.square(proportions).sum())
        collapse.append(
            {
                "route": route,
                "group": group,
                "n": int(len(sub)),
                "active_networks": int(len(counts)),
                "dominant_network": str(counts.index[0]),
                "dominant_count": int(counts.iloc[0]),
                "dominant_fraction": float(proportions.iloc[0]),
                "hhi": hhi,
                "effective_network_count": float(1.0 / hhi),
                "collapse_flag_dominant_fraction_ge_0p70": bool(
                    proportions.iloc[0] >= 0.70
                ),
            }
        )
    return distribution, collapse


def plot_distributions(output: pd.DataFrame, routes: list[str], outdir: Path) -> None:
    for route in routes:
        counts = (
            output.groupby(["group", f"{route}_top1_network"])
            .size()
            .rename("n")
            .reset_index()
        )
        counts["fraction"] = counts["n"] / counts.groupby("group")["n"].transform("sum")
        pivot = counts.pivot(
            index=f"{route}_top1_network", columns="group", values="fraction"
        ).fillna(0.0)
        pivot.plot(kind="bar", figsize=(11, 5), color=["#C44E52", "#4C72B0"])
        plt.ylabel("Top1 fraction")
        plt.xlabel("Network")
        plt.title(f"GSE228512 {route}: GBM vs healthy")
        plt.xticks(rotation=35, ha="right")
        plt.tight_layout()
        plt.savefig(outdir / f"{route}_network_top1_distribution.png", dpi=180)
        plt.close()

        fig, axes = plt.subplots(1, 3, figsize=(12, 4))
        for ax, metric, title in zip(
            axes,
            ("normalized_entropy", "raw_margin", "max_score"),
            ("Normalized entropy", "Top1/Top2 raw margin", "Maximum correlation"),
        ):
            data = [
                output.loc[output.group.eq(group), f"{route}_{metric}"].to_numpy()
                for group in ("GBM", "Healthy")
            ]
            ax.boxplot(data, tick_labels=["GBM", "Healthy"], showfliers=False)
            ax.set_title(title)
        fig.suptitle(f"GSE228512 {route}: confidence and OOD")
        fig.tight_layout()
        fig.savefig(outdir / f"{route}_confidence_ood.png", dpi=180)
        plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", type=Path, default=DATA_DIR)
    parser.add_argument("--outdir", type=Path, default=OUTDIR)
    parser.add_argument("--permutations", type=int, default=10000)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    counts, metadata = load_counts_and_metadata(args.data_dir)
    mapping = ensembl_symbol_map(ORTHOLOGY)
    cpm = counts_to_symbol_cpm(counts, mapping)
    output = metadata.copy()

    baseline_scores, baseline_networks, baseline_genes, baseline_present = model_scores(
        cpm, BASELINE_MODEL, "logcpm"
    )
    output = add_route(output, "baseline_logcpm", baseline_scores, baseline_networks)
    production_top1 = []
    production_top3 = []
    production_switched = []
    logcpm = np.log1p(cpm)
    for sample in logcpm.columns:
        expression = pd.DataFrame(
            {"gene_symbol": logcpm.index, "tpm_value": logcpm[sample].to_numpy()}
        )
        traced = trace_network_expression(expression, min_overlap_fraction=0.50)
        ranked = traced.get("results", [])
        production_top1.append(ranked[0]["network_id"])
        production_top3.append(" | ".join(row["network_id"] for row in ranked[:3]))
        production_switched.append(
            bool(traced.get("meta", {}).get("pairwise_rescue", {}).get("switched", False))
        )
    output["baseline_production_top1_network"] = production_top1
    output["baseline_production_top3_network"] = production_top3
    output["baseline_production_top1_broad"] = output[
        "baseline_production_top1_network"
    ].map(NETWORK_TO_PRIMARY_BROAD)
    output["baseline_production_top3_broad"] = output[
        "baseline_production_top3_network"
    ].map(
        lambda value: " | ".join(
            dict.fromkeys(
                NETWORK_TO_PRIMARY_BROAD[item.strip()]
                for item in str(value).split("|")
                if item.strip()
            )
        )
    )
    output["baseline_production_pairwise_switched"] = production_switched

    adapted_model = np.load(ADAPTED_MODEL, allow_pickle=False)
    threshold = float(adapted_model["ood_max_correlation_threshold"][0])
    route_specs = (
        ("adapted_cpm", "cpm"),
        ("adapted_logcpm", "logcpm"),
        ("adapted_harmonized", "harmonized"),
    )
    adapted_present = None
    for route, transform in route_specs:
        scores, networks, genes, present = model_scores(cpm, ADAPTED_MODEL, transform)
        adapted_present = present
        output = add_route(output, route, scores, networks, threshold)

    routes = ["baseline_logcpm", *[item[0] for item in route_specs]]
    continuous, distribution, collapse = summarize_groups(
        output, routes, args.permutations
    )
    production_distribution, production_collapse = summarize_prediction_only(
        output, "baseline_production", args.permutations
    )
    distribution = pd.concat(
        [distribution, pd.DataFrame([production_distribution])], ignore_index=True
    )
    collapse = pd.concat(
        [collapse, pd.DataFrame(production_collapse)], ignore_index=True
    )
    continuous["mann_whitney_fdr"] = bh_adjust(continuous["mann_whitney_p"])
    continuous["platform_stratified_permutation_fdr"] = bh_adjust(
        continuous["platform_stratified_permutation_p"]
    )
    distribution["chi_square_fdr"] = bh_adjust(distribution["chi_square_p"])
    distribution["platform_stratified_permutation_fdr"] = bh_adjust(
        distribution["platform_stratified_permutation_p"]
    )
    prediction_routes = ["baseline_production", *routes]
    network_long = []
    for route in prediction_routes:
        for group, sub in output.groupby("group"):
            for network, count in sub[f"{route}_top1_network"].value_counts().items():
                network_long.append(
                    {
                        "route": route,
                        "group": group,
                        "network": network,
                        "count": int(count),
                        "fraction": float(count / len(sub)),
                    }
                )
    broad_long = []
    for route in prediction_routes:
        for group, sub in output.groupby("group"):
            for broad, count in sub[f"{route}_top1_broad"].value_counts().items():
                broad_long.append(
                    {
                        "route": route,
                        "group": group,
                        "broad_anatomy": broad,
                        "count": int(count),
                        "fraction": float(count / len(sub)),
                    }
                )

    output.to_csv(args.outdir / "gse228512_sample_predictions.csv", index=False, encoding="utf-8-sig")
    pd.DataFrame(network_long).to_csv(
        args.outdir / "network_prediction_distribution.csv", index=False, encoding="utf-8-sig"
    )
    pd.DataFrame(broad_long).to_csv(
        args.outdir / "broad_anatomy_prediction_distribution.csv",
        index=False,
        encoding="utf-8-sig",
    )
    top3_network_rows = []
    top3_broad_rows = []
    for route in prediction_routes:
        for group, sub in output.groupby("group"):
            network_counts = Counter(
                item.strip()
                for value in sub[f"{route}_top3_network"]
                for item in str(value).split("|")
                if item.strip()
            )
            broad_counts = Counter(
                item.strip()
                for value in sub[f"{route}_top3_broad"]
                for item in str(value).split("|")
                if item.strip()
            )
            for label, count in network_counts.items():
                top3_network_rows.append(
                    {
                        "route": route,
                        "group": group,
                        "network": label,
                        "samples_with_label_in_top3": int(count),
                        "fraction": float(count / len(sub)),
                    }
                )
            for label, count in broad_counts.items():
                top3_broad_rows.append(
                    {
                        "route": route,
                        "group": group,
                        "broad_anatomy": label,
                        "samples_with_label_in_top3": int(count),
                        "fraction": float(count / len(sub)),
                    }
                )
    pd.DataFrame(top3_network_rows).to_csv(
        args.outdir / "network_top3_occurrence_distribution.csv",
        index=False,
        encoding="utf-8-sig",
    )
    pd.DataFrame(top3_broad_rows).to_csv(
        args.outdir / "broad_anatomy_top3_occurrence_distribution.csv",
        index=False,
        encoding="utf-8-sig",
    )
    continuous.to_csv(args.outdir / "group_continuous_metric_tests.csv", index=False, encoding="utf-8-sig")
    distribution.to_csv(args.outdir / "group_network_distribution_tests.csv", index=False, encoding="utf-8-sig")
    collapse.to_csv(args.outdir / "prediction_collapse_metrics.csv", index=False, encoding="utf-8-sig")
    score_rows = []
    for route, networks, scores in (
        ("baseline_logcpm", baseline_networks, baseline_scores),
        *[
            (
                route,
                np.load(ADAPTED_MODEL, allow_pickle=False)["networks"].astype(str),
                output[
                    [f"{route}_score_{index + 1:02d}" for index in range(10)]
                ].to_numpy().T,
            )
            for route, _ in route_specs
        ],
    ):
        for group in ("GBM", "Healthy"):
            mask = output["group"].eq(group).to_numpy()
            for index, network in enumerate(networks):
                score_rows.append(
                    {
                        "route": route,
                        "group": group,
                        "network": network,
                        "mean_score": float(scores[index, mask].mean()),
                        "median_score": float(np.median(scores[index, mask])),
                    }
                )
    pd.DataFrame(score_rows).to_csv(
        args.outdir / "network_score_bias_diagnostic.csv",
        index=False,
        encoding="utf-8-sig",
    )
    adapted_gene_table = pd.read_csv(
        ROOT
        / "results"
        / "tcga_lgg_65_tumor_adapted_network_20260612"
        / "tumor_adapted_network_genes.csv"
    )
    detection_rows = []
    for row in adapted_gene_table.itertuples(index=False):
        gene = str(row.gene_symbol)
        values = cpm.loc[gene] if gene in cpm.index else pd.Series(0.0, index=cpm.columns)
        item = {"gene_symbol": gene, "balanced_owner_network": row.balanced_owner_network}
        for group in ("GBM", "Healthy"):
            samples = output.loc[output.group.eq(group), "sample_id"]
            group_values = values.reindex(samples).fillna(0.0)
            item[f"{group.lower()}_detection_fraction"] = float((group_values > 0).mean())
            item[f"{group.lower()}_median_cpm"] = float(group_values.median())
        detection_rows.append(item)
    pd.DataFrame(detection_rows).to_csv(
        args.outdir / "adapted_marker_detection_diagnostic.csv",
        index=False,
        encoding="utf-8-sig",
    )
    plot_distributions(output, routes, args.outdir)

    summary = {
        "dataset": "GSE228512",
        "n_samples": int(len(output)),
        "group_counts": output["group"].value_counts().to_dict(),
        "platform_counts": output["platform"].value_counts().to_dict(),
        "platform_by_group": pd.crosstab(output["platform"], output["group"]).to_dict(),
        "n_ensembl_rows": int(len(counts)),
        "n_mapped_gene_symbols": int(len(cpm)),
        "model_gene_coverage": {
            "baseline": {
                "present": int(baseline_present.sum()),
                "total": int(len(baseline_genes)),
                "fraction": float(baseline_present.mean()),
            },
            "adapted": {
                "present": int(adapted_present.sum()),
                "total": int(len(adapted_present)),
                "fraction": float(adapted_present.mean()),
            },
        },
        "ood_threshold": threshold,
        "ood_acceptance": {
            route: {
                group: {
                    "accepted": int(
                        output.loc[output.group.eq(group), f"{route}_ood_accepted"].sum()
                    ),
                    "n": int(output.group.eq(group).sum()),
                    "coverage": float(
                        output.loc[output.group.eq(group), f"{route}_ood_accepted"].mean()
                    ),
                }
                for group in ("GBM", "Healthy")
            }
            for route, _ in route_specs
        },
        "baseline_production_pairwise_rescue": {
            "n_switched": int(np.sum(production_switched)),
            "fraction_switched": float(np.mean(production_switched)),
        },
        "no_localization_accuracy": (
            "GSE228512 has no patient-level tumor lobe, coordinate, or MRI identifier; "
            "all outputs are distribution/domain-shift analyses."
        ),
    }
    (args.outdir / "analysis_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
