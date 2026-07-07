#!/usr/bin/env python
from __future__ import annotations

import json
import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
OUTDIR = ROOT / "reports" / "p0_hard_evidence_20260629"
RESULT_ROOT = ROOT / "results" / "bo2023_reference_projection_20260616_cleaned_symbols"
FORMAL_LOSO_DIR = OUTDIR / "formal_loso_rerun_20260629"
if not (FORMAL_LOSO_DIR / "hybrid_formal_loso_network_detail.csv").exists():
    FORMAL_LOSO_DIR = RESULT_ROOT / "formal_three_tier_loso_hybrid"
HYBRID = "hybrid_projected_network_logcpm_exact"
LABEL_CURATION_MAP = ROOT / "reports" / "bo2023_publication_label_curation_map_20260704.csv"
RNG_SEED = 20260629
N_WEIGHTED_RANDOM_DRAWS = 10000


NETWORK_ANATOMY = {
    "Cingulate gyrus": "Medial cingulate cortex; limbic/medial association cortex.",
    "Frontal (agranular frontal motor areas)": "Agranular frontal motor and premotor cortex.",
    "Hippocampal formation": "Hippocampal/parahippocampal formation.",
    "Lateral Prefrontal Cortex": "Dorsolateral and ventrolateral prefrontal association cortex.",
    "Occipital/Temporal": "Visual occipital cortex and occipito-temporal visual association cortex.",
    "Operculum/Insula": "Opercular and insular cortex.",
    "Orbitomedial Prefrontal Cortex (OMPFC)": "Orbitofrontal and medial prefrontal cortex.",
    "Parietal, and Parieto-occipital region": "Posterior parietal and parieto-occipital association cortex.",
    "Subcortical": "Basal ganglia, amygdala, thalamic and related subcortical labels in this reference.",
    "Temporal": "Temporal association and auditory-related cortex.",
}


@dataclass
class Endpoint:
    name: str
    path: Path
    y_true: str
    y_pred: str
    hit1: str
    hit3: str
    candidate_count: str | None = None
    route_column: str | None = None
    route_value: str | None = None
    supported_column: str | None = None
    supported_value: object | None = None
    top_cols: tuple[str, str, str] | None = None
    single_label_truth: bool = True


def read_endpoint(ep: Endpoint) -> pd.DataFrame:
    df = pd.read_csv(ep.path)
    if ep.route_column and ep.route_value:
        df = df[df[ep.route_column].astype(str).eq(ep.route_value)].copy()
    if ep.supported_column is not None:
        df = df[df[ep.supported_column].eq(ep.supported_value)].copy()
    return df.reset_index(drop=True)


def read_label_curation_map(path: Path = LABEL_CURATION_MAP) -> dict[str, str]:
    if not path.exists():
        return {}
    curation = pd.read_csv(path, dtype=str).fillna("")
    return {
        str(row.old_region_id).strip(): str(row.new_region_id).strip()
        for row in curation.itertuples(index=False)
        if str(row.old_region_id).strip()
        and str(row.new_region_id).strip()
        and str(row.old_region_id).strip() != str(row.new_region_id).strip()
    }


def remap_label(value: object, mapping: dict[str, str]) -> object:
    if pd.isna(value):
        return value
    text = str(value)
    if text in mapping:
        return mapping[text]
    if "|" in text:
        parts = [part.strip() for part in text.split("|")]
        remapped = [mapping.get(part, part) for part in parts]
        if remapped != parts:
            return " | ".join(remapped)
    remapped_text = text
    for old, new in mapping.items():
        remapped_text = re.sub(
            rf"(?<![A-Za-z0-9_]){re.escape(old)}(?![A-Za-z0-9_])",
            new,
            remapped_text,
        )
    return remapped_text if remapped_text != text else value


def wilson_ci(k: int, n: int, z: float = 1.959963984540054) -> tuple[float, float]:
    if n == 0:
        return (float("nan"), float("nan"))
    p = k / n
    denom = 1 + z * z / n
    center = (p + z * z / (2 * n)) / denom
    half = z * math.sqrt((p * (1 - p) + z * z / (4 * n)) / n) / denom
    return center - half, center + half


def binom_sf(k: int, n: int, p: float) -> float:
    if n <= 0 or not np.isfinite(p):
        return float("nan")
    p = min(max(float(p), 0.0), 1.0)
    try:
        from scipy.stats import binomtest

        return float(binomtest(k, n, p, alternative="greater").pvalue)
    except Exception:
        # Normal approximation fallback; used only if scipy is unavailable.
        mean = n * p
        var = n * p * (1 - p)
        if var == 0:
            return 0.0 if k > mean else 1.0
        z = ((k - 0.5) - mean) / math.sqrt(var)
        return 0.5 * math.erfc(z / math.sqrt(2))


def mcnemar_exact(x: Iterable[int], y: Iterable[int]) -> dict[str, float | int]:
    a = np.asarray(list(x)).astype(int)
    b = np.asarray(list(y)).astype(int)
    only_x = int(((a == 1) & (b == 0)).sum())
    only_y = int(((a == 0) & (b == 1)).sum())
    discordant = only_x + only_y
    if discordant == 0:
        p = 1.0
    else:
        try:
            from scipy.stats import binomtest

            p = float(binomtest(min(only_x, only_y), discordant, 0.5, alternative="two-sided").pvalue)
        except Exception:
            # Conservative exact two-sided binomial fallback.
            tail = sum(math.comb(discordant, i) for i in range(0, min(only_x, only_y) + 1))
            p = min(1.0, 2 * tail / (2**discordant))
    return {"only_loso_hit": only_x, "only_lomo_hit": only_y, "discordant": discordant, "mcnemar_exact_p": p}


def class_f1_table(df: pd.DataFrame, y_true: str, y_pred: str) -> pd.DataFrame:
    labels = sorted(set(df[y_true].dropna().astype(str)) | set(df[y_pred].dropna().astype(str)))
    rows = []
    truth = df[y_true].astype(str)
    pred = df[y_pred].astype(str)
    for label in labels:
        tp = int(((truth == label) & (pred == label)).sum())
        fp = int(((truth != label) & (pred == label)).sum())
        fn = int(((truth == label) & (pred != label)).sum())
        precision = tp / (tp + fp) if tp + fp else float("nan")
        recall = tp / (tp + fn) if tp + fn else float("nan")
        f1 = 2 * precision * recall / (precision + recall) if precision + recall and np.isfinite(precision + recall) else float("nan")
        rows.append({"label": label, "support": int((truth == label).sum()), "tp": tp, "fp": fp, "fn": fn, "precision": precision, "recall": recall, "f1": f1})
    return pd.DataFrame(rows)


def confusion_matrix(df: pd.DataFrame, y_true: str, y_pred: str) -> pd.DataFrame:
    return pd.crosstab(df[y_true].astype(str), df[y_pred].astype(str), rownames=["truth"], colnames=["pred"], dropna=False)


def weighted_topk_baseline(df: pd.DataFrame, y_true: str, top_k: int, candidate_count: str | None) -> tuple[float, str]:
    df = df[df[y_true].notna()].copy()
    truth = df[y_true].astype(str)
    labels = sorted(truth.unique())
    if not labels:
        return float("nan"), "no_truth_labels"
    weights = truth.value_counts(normalize=True).reindex(labels).fillna(0).to_numpy()
    label_to_idx = {label: i for i, label in enumerate(labels)}
    truth_idx = df[y_true].astype(str).map(label_to_idx).to_numpy()
    if candidate_count and candidate_count in df.columns:
        ks = np.minimum(top_k, df[candidate_count].fillna(len(labels)).astype(int).clip(lower=1).to_numpy())
        # Candidate identities are not fully serialized for region endpoints; the weighted baseline is
        # therefore a global-prior approximation over the endpoint label space.
        note = "global_truth_prior_weighted_without_replacement; candidate_count_limits_k_only"
    else:
        ks = np.full(len(df), min(top_k, len(labels)), dtype=int)
        note = "global_truth_prior_weighted_without_replacement"
    rng = np.random.default_rng(RNG_SEED + top_k + len(df))
    inclusion_by_k: dict[int, np.ndarray] = {}
    for k in sorted(set(int(x) for x in ks)):
        kk = min(k, len(labels))
        counts = np.zeros(len(labels), dtype=float)
        for _ in range(N_WEIGHTED_RANDOM_DRAWS):
            draw = rng.choice(len(labels), size=kk, replace=False, p=weights)
            counts[draw] += 1.0
        inclusion_by_k[k] = counts / N_WEIGHTED_RANDOM_DRAWS
    probs = np.array([inclusion_by_k[int(k)][int(idx)] for k, idx in zip(ks, truth_idx)])
    return float(probs.mean()), note


def uniform_baseline(df: pd.DataFrame, y_true: str, top_k: int, candidate_count: str | None) -> tuple[float, str]:
    df = df[df[y_true].notna()].copy()
    if candidate_count and candidate_count in df.columns:
        n = df[candidate_count].fillna(0).astype(float).clip(lower=1)
        return float(np.minimum(top_k / n, 1.0).mean()), f"sample_specific_uniform_top{top_k}_over_{candidate_count}"
    n_classes = max(1, df[y_true].dropna().astype(str).nunique())
    return min(top_k / n_classes, 1.0), f"uniform_top{top_k}_over_{n_classes}_truth_classes"


def metric_rows(endpoints: list[Endpoint]) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    metrics = []
    baselines = []
    tests = []
    for ep in endpoints:
        df = read_endpoint(ep)
        for hit_col, top_k in [(ep.hit1, 1), (ep.hit3, 3)]:
            if hit_col not in df.columns:
                continue
            hits = int(pd.to_numeric(df[hit_col], errors="coerce").fillna(0).sum())
            n = int(pd.to_numeric(df[hit_col], errors="coerce").notna().sum())
            ci_low, ci_high = wilson_ci(hits, n)
            acc = hits / n if n else float("nan")
            metrics.append({
                "endpoint": ep.name,
                "metric": f"Top{top_k}",
                "n": n,
                "hits": hits,
                "accuracy": acc,
                "wilson95_low": ci_low,
                "wilson95_high": ci_high,
                "source_file": str(ep.path.relative_to(ROOT)),
            })
            uniform_p, uniform_note = uniform_baseline(df, ep.y_true, top_k, ep.candidate_count)
            weighted_p, weighted_note = weighted_topk_baseline(df, ep.y_true, top_k, ep.candidate_count)
            baselines.extend([
                {"endpoint": ep.name, "metric": f"Top{top_k}", "baseline": "uniform_random", "expected_accuracy": uniform_p, "method_note": uniform_note},
                {"endpoint": ep.name, "metric": f"Top{top_k}", "baseline": "weighted_random", "expected_accuracy": weighted_p, "method_note": weighted_note},
            ])
            tests.extend([
                {"endpoint": ep.name, "metric": f"Top{top_k}", "baseline": "uniform_random", "observed_hits": hits, "n": n, "p0": uniform_p, "one_sided_binomial_p": binom_sf(hits, n, uniform_p)},
                {"endpoint": ep.name, "metric": f"Top{top_k}", "baseline": "weighted_random", "observed_hits": hits, "n": n, "p0": weighted_p, "one_sided_binomial_p": binom_sf(hits, n, weighted_p)},
            ])
    return pd.DataFrame(metrics), pd.DataFrame(baselines), pd.DataFrame(tests)


def write_confusion_and_f1(endpoints: list[Endpoint]) -> None:
    (OUTDIR / "confusion_matrices").mkdir(parents=True, exist_ok=True)
    (OUTDIR / "class_f1").mkdir(parents=True, exist_ok=True)
    for ep in endpoints:
        if not ep.single_label_truth:
            continue
        df = read_endpoint(ep)
        safe = ep.name.replace(" ", "_").replace("/", "_")
        confusion_matrix(df, ep.y_true, ep.y_pred).to_csv(OUTDIR / "confusion_matrices" / f"{safe}_top1_confusion.csv")
        class_f1_table(df, ep.y_true, ep.y_pred).to_csv(OUTDIR / "class_f1" / f"{safe}_top1_class_f1.csv", index=False)


def write_paired_tests() -> pd.DataFrame:
    pairs = [
        (
            "formal_internal_network_top3",
            FORMAL_LOSO_DIR / "hybrid_formal_loso_network_detail.csv",
            RESULT_ROOT / "formal_three_tier_lomo_hybrid" / "formal_lomo_network_detail.csv",
            "hit3",
            None,
        ),
        (
            "formal_internal_resolution_group_top3",
            FORMAL_LOSO_DIR / "hybrid_formal_loso_resolution_group_detail.csv",
            RESULT_ROOT / "formal_three_tier_lomo_hybrid" / "formal_lomo_resolution_group_detail.csv",
            "group_hit3",
            HYBRID,
        ),
        (
            "formal_internal_exact_region_top3",
            FORMAL_LOSO_DIR / "hybrid_formal_loso_exact_region_detail.csv",
            RESULT_ROOT / "formal_three_tier_lomo_hybrid" / "formal_lomo_exact_region_detail.csv",
            "hit3",
            HYBRID,
        ),
    ]
    rows = []
    for name, loso_path, lomo_path, hit_col, lomo_route in pairs:
        loso = pd.read_csv(loso_path)
        lomo = pd.read_csv(lomo_path)
        if "route_family" in loso.columns:
            loso = loso[loso["route_family"].eq(HYBRID)]
        if lomo_route and "route_family" in lomo.columns:
            lomo = lomo[lomo["route_family"].eq(lomo_route)]
        elif "route_family" in lomo.columns:
            lomo = lomo[lomo["route_family"].eq(HYBRID)]
        merged = loso[["sample_id", hit_col]].merge(lomo[["sample_id", hit_col]], on="sample_id", suffixes=("_loso", "_lomo"))
        stats = mcnemar_exact(merged[f"{hit_col}_loso"], merged[f"{hit_col}_lomo"])
        rows.append({"comparison": name, "paired_n": len(merged), **stats, "test": "exact McNemar/binomial on discordant pairs"})
    out = pd.DataFrame(rows)
    out.to_csv(OUTDIR / "loso_lomo_paired_tests.csv", index=False)
    return out


def write_hierarchy_tables() -> None:
    label_map = read_label_curation_map()
    entries = json.loads((ROOT / "data" / "models" / "bo2023_region_resolution_groups.json").read_text(encoding="utf-8"))["entries"]
    rows = []
    for item in entries.values():
        region = remap_label(item.get("region_id"), label_map)
        resolution_group = remap_label(item.get("resolution_group"), label_map)
        group_members = [remap_label(member, label_map) for member in item.get("group_members", [])]
        rows.append({
            "network": item.get("network_id"),
            "region": region,
            "resolution_group": resolution_group,
            "group_members": " | ".join(group_members),
            "resolution_tier": item.get("resolution_tier"),
            "resolution_reasons": ";".join(item.get("resolution_reasons", [])),
            "training_samples": item.get("training_samples"),
            "nearest_centroid_corr": item.get("nearest_centroid_corr"),
            "group_plausibility_tier": item.get("group_plausibility_tier"),
        })
    hierarchy = pd.DataFrame(rows)
    hierarchy.to_csv(OUTDIR / "resolution_group_hierarchy.csv", index=False)
    net = (
        hierarchy.groupby("network")
        .agg(
            n_regions=("region", "nunique"),
            n_resolution_groups=("resolution_group", "nunique"),
            n_low_resolution_regions=("resolution_tier", lambda s: int((s == "low_resolution").sum())),
            n_training_samples=("training_samples", "sum"),
        )
        .reset_index()
    )
    net["classical_neuroanatomy_mapping"] = net["network"].map(NETWORK_ANATOMY).fillna("Manual review required.")
    net.to_csv(OUTDIR / "network_anatomy_table.csv", index=False)


def write_denominator_audit(endpoints: list[Endpoint], metrics: pd.DataFrame) -> pd.DataFrame:
    expected = {
        "projected_vsd_network_loso": 819,
        "projected_vsd_network_lomo": 819,
        "formal_internal_network_loso": 819,
        "formal_internal_network_lomo": 819,
        "formal_internal_resolution_group_loso": 814,
        "formal_internal_resolution_group_lomo": 812,
        "formal_internal_exact_region_loso": 814,
        "formal_internal_exact_region_lomo": 812,
        "ahba_hybrid_network_supported": 233,
        "ahba_hybrid_resolution_group_exact_mapped": 91,
        "ahba_hybrid_exact_region_exact_mapped": 91,
        "tcga_brats_hybrid_network": 65,
        "tcga_brats_hybrid_broad_anatomy": 65,
    }
    rows = []
    for ep in endpoints:
        df = read_endpoint(ep)
        rows.append({
            "endpoint": ep.name,
            "observed_detail_rows": len(df),
            "expected_n_from_current_manuscript_or_review": expected.get(ep.name),
            "matches_expected": len(df) == expected.get(ep.name) if ep.name in expected else "",
            "source_file": str(ep.path.relative_to(ROOT)),
        })
    out = pd.DataFrame(rows)
    out.to_csv(OUTDIR / "denominator_audit.csv", index=False)
    return out


def write_mapping_degradation() -> None:
    ahba = pd.read_csv(RESULT_ROOT / "ahba_external_formal_three_tier" / "ahba_formal_three_tier_sample_detail.csv")
    ahba = ahba[ahba["route"].eq(HYBRID)]
    rows = ahba["accuracy_level"].value_counts(dropna=False).rename_axis("mapping_level").reset_index(name="n_samples")
    rows["fraction"] = rows["n_samples"] / len(ahba)
    rows.to_csv(OUTDIR / "ahba_mapping_granularity_summary.csv", index=False)


def write_methods_and_checklist(denom: pd.DataFrame) -> None:
    methods = f"""# P0 statistical evidence methods

Generated by `scripts/generate_p0_hard_evidence.py`.

## Accuracy intervals

All Top1/Top3 proportions are reported with two-sided 95% Wilson score confidence intervals. Wilson intervals are used instead of normal Wald intervals because several endpoints have moderate sample sizes or imbalanced class distributions.

## Random baselines

Uniform random baselines are endpoint-aware:

- Network endpoints use uniform Top-k recovery over the observed Network truth label space.
- Internal exact-region and resolution-group endpoints use sample-specific candidate counts (`n_candidate_regions` / `n_candidate_groups`) when available.
- External AHBA/TCGA endpoints with multi-label truth are reported for accuracy and CI, but class-level confusion is restricted to single-label internal endpoints.

Weighted random baselines use the empirical truth-label prior and fixed-seed weighted sampling without replacement (`seed={RNG_SEED}`, draws={N_WEIGHTED_RANDOM_DRAWS}). For region endpoints where full candidate identity lists are not serialized, this is a global-prior approximation and the output table states that limitation.

## Significance tests

Each observed hit count is tested against the corresponding random baseline with a one-sided binomial test. LOSO-vs-LOMO comparisons use paired exact McNemar/binomial tests on intersecting `sample_id` values.

## Confusion and F1

Top1 confusion matrices and per-class precision/recall/F1 are generated for single-label internal endpoints. AHBA and TCGA use allowed-label sets rather than single truth labels, so confusion matrices would be misleading and are not generated for those external mapped-label/coarse-consistency analyses.
"""
    (OUTDIR / "STATISTICAL_METHODS.md").write_text(methods, encoding="utf-8")

    loso_network = denom[denom["endpoint"].eq("formal_internal_network_loso")]
    loso_network_done = bool(len(loso_network) and bool(loso_network["matches_expected"].iloc[0]))
    checklist = [
        ("95% CI for main accuracy metrics", True, "metric_summary_with_ci.csv"),
        ("Uniform and weighted random baselines", True, "random_baselines.csv"),
        ("Binomial tests vs random baselines", True, "binomial_tests_vs_random.csv"),
        ("LOSO vs LOMO paired significance tests", True, "loso_lomo_paired_tests.csv"),
        ("Confusion matrices for single-label internal endpoints", True, "confusion_matrices/*.csv"),
        ("Class-level precision/recall/F1 for single-label internal endpoints", True, "class_f1/*.csv"),
        ("10 Network table and classical neuroanatomy correspondence", True, "network_anatomy_table.csv"),
        ("Resolution-group hierarchy and reasons", True, "resolution_group_hierarchy.csv"),
        ("AHBA mapping granularity/degradation summary", True, "ahba_mapping_granularity_summary.csv"),
        ("Formal LOSO Network denominator 819 fully traceable in current detail file", loso_network_done, "denominator_audit.csv"),
        (
            "Marker selection criteria, counts, cross-species mapping policy and cfRNA degradation boundary",
            (OUTDIR / "marker_methodology_audit.csv").exists() and (OUTDIR / "MARKER_METHODOLOGY_REPORT.md").exists(),
            "marker_methodology_audit.csv; MARKER_METHODOLOGY_REPORT.md",
        ),
    ]
    lines = ["# P0 hard-evidence completion checklist", ""]
    for item, done, evidence in checklist:
        lines.append(f"- [{'x' if done else ' '}] {item} - {evidence}")
    lines.append("")
    if not loso_network_done:
        lines.append("## Open audit item")
        lines.append("")
        lines.append("The current formal LOSO Network detail file contains 814 rows, while the manuscript/review expectation is 819 Network-evaluable samples. Region-level LOSO denominators are correctly 814. This package therefore completes the statistical machinery but flags the LOSO Network denominator as requiring rerun or replacement with a traceable 819-row formal Network detail file before claiming full P0 closure.")
    (OUTDIR / "P0_COMPLETION_CHECKLIST.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    endpoints = [
        Endpoint("projected_vsd_network_loso", RESULT_ROOT / "bo2023_projected_vsd_loso_detail.csv", "label", "pred_top1", "hit1", "hit3", route_column="route", route_value="projected_vsd"),
        Endpoint("projected_vsd_network_lomo", RESULT_ROOT / "bo2023_projected_vsd_lomo_detail.csv", "label", "pred_top1", "hit1", "hit3", route_column="route", route_value="projected_vsd"),
        Endpoint("formal_internal_network_loso", FORMAL_LOSO_DIR / "hybrid_formal_loso_network_detail.csv", "label", "pred_top1", "hit1", "hit3"),
        Endpoint("formal_internal_network_lomo", RESULT_ROOT / "formal_three_tier_lomo_hybrid" / "formal_lomo_network_detail.csv", "label", "pred_top1", "hit1", "hit3", route_column="route_family", route_value=HYBRID),
        Endpoint("formal_internal_resolution_group_loso", FORMAL_LOSO_DIR / "hybrid_formal_loso_resolution_group_detail.csv", "true_resolution_group", "pred_group_top1", "group_hit1", "group_hit3", candidate_count="n_candidate_groups"),
        Endpoint("formal_internal_resolution_group_lomo", RESULT_ROOT / "formal_three_tier_lomo_hybrid" / "formal_lomo_resolution_group_detail.csv", "true_resolution_group", "pred_group_top1", "group_hit1", "group_hit3", candidate_count="n_candidate_groups", route_column="route_family", route_value=HYBRID),
        Endpoint("formal_internal_exact_region_loso", FORMAL_LOSO_DIR / "hybrid_formal_loso_exact_region_detail.csv", "label", "pred_top1", "hit1", "hit3", candidate_count="n_candidate_regions"),
        Endpoint("formal_internal_exact_region_lomo", RESULT_ROOT / "formal_three_tier_lomo_hybrid" / "formal_lomo_exact_region_detail.csv", "label", "pred_top1", "hit1", "hit3", candidate_count="n_candidate_regions", route_column="route_family", route_value=HYBRID),
        Endpoint("ahba_hybrid_network_supported", RESULT_ROOT / "ahba_external_formal_three_tier" / "ahba_formal_three_tier_sample_detail.csv", "allowed_bo2023_networks", "network_top1", "network_top1_hit", "network_top3_hit", route_column="route", route_value=HYBRID, supported_column="supported_for_accuracy", supported_value=True, single_label_truth=False),
        Endpoint("ahba_hybrid_resolution_group_exact_mapped", RESULT_ROOT / "ahba_external_formal_three_tier" / "ahba_formal_three_tier_sample_detail.csv", "allowed_resolution_groups", "group_top1", "group_top1_hit", "group_top3_hit", route_column="route", route_value=HYBRID, supported_column="accuracy_level", supported_value="exact_region", single_label_truth=False),
        Endpoint("ahba_hybrid_exact_region_exact_mapped", RESULT_ROOT / "ahba_external_formal_three_tier" / "ahba_formal_three_tier_sample_detail.csv", "allowed_bo2023_regions", "region_top1", "region_top1_exact_hit", "region_top3_exact_hit", route_column="route", route_value=HYBRID, supported_column="accuracy_level", supported_value="exact_region", single_label_truth=False),
        Endpoint("tcga_brats_hybrid_network", RESULT_ROOT / "tcga_labeled_hybrid_formal_external" / "tcga_labeled_hybrid_formal_sample_detail.csv", "truth_network_candidates", "network_top1", "network_top1_hit", "network_top3_hit", route_column="route", route_value=HYBRID, single_label_truth=False),
        Endpoint("tcga_brats_hybrid_broad_anatomy", RESULT_ROOT / "tcga_labeled_hybrid_formal_external" / "tcga_labeled_hybrid_formal_sample_detail.csv", "truth_broad_candidates", "pred_broad_top1", "broad_top1_hit", "broad_top3_hit", route_column="route", route_value=HYBRID, single_label_truth=False),
    ]
    metrics, baselines, tests = metric_rows(endpoints)
    metrics.to_csv(OUTDIR / "metric_summary_with_ci.csv", index=False)
    baselines.to_csv(OUTDIR / "random_baselines.csv", index=False)
    tests.to_csv(OUTDIR / "binomial_tests_vs_random.csv", index=False)
    write_confusion_and_f1(endpoints)
    write_paired_tests()
    write_hierarchy_tables()
    write_mapping_degradation()
    denom = write_denominator_audit(endpoints, metrics)
    write_methods_and_checklist(denom)
    manifest = {
        "created_by": "scripts/generate_p0_hard_evidence.py",
        "output_dir": str(OUTDIR.relative_to(ROOT)),
        "files": sorted(str(p.relative_to(OUTDIR)) for p in OUTDIR.rglob("*") if p.is_file()),
    }
    (OUTDIR / "manifest.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"Wrote P0 evidence package to {OUTDIR}")


if __name__ == "__main__":
    main()
