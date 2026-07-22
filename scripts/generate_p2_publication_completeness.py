#!/usr/bin/env python
from __future__ import annotations

import json
import subprocess
import sys
import warnings
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy.stats import f_oneway
from sklearn.base import clone
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.metrics import confusion_matrix, precision_recall_fscore_support
from sklearn.neighbors import KNeighborsClassifier, NearestCentroid
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler


ROOT = Path(__file__).resolve().parents[1]
OUTDIR = ROOT / "reports" / "p2_publication_completeness_20260629"
ML_DIR = OUTDIR / "ml_baselines"
ENG_DIR = OUTDIR / "engineering_reproducibility"
TOOL_DIR = OUTDIR / "tool_comparison"
SYNC_DIR = OUTDIR / "sync_text"
SEED = 20260629


@dataclass
class ModelSpec:
    name: str
    model: Any
    score_mode: str
    k_features: int = 500


def read_bo2023_gene_matrix(path: Path) -> pd.DataFrame:
    with path.open("rt", encoding="utf-8") as handle:
        samples = handle.readline().rstrip("\n\r").split("\t")
    names = ["gene_id", *samples]
    return pd.read_csv(path, sep="\t", header=None, names=names, skiprows=1).set_index("gene_id")


def read_metadata() -> pd.DataFrame:
    path = ROOT / "bo2023 data" / "Information of sequenced samples_update_full878_filter819.xlsx"
    info = pd.read_excel(path, sheet_name="mfas5_819samples_phenSet4", usecols=["No.", "Region", "SaleemNetworks", "MonkeyID"])
    info["sample_id"] = info["No."].astype(str).str.strip()
    info["network"] = info["SaleemNetworks"].astype(str).str.strip()
    info["region"] = info["Region"].astype(str).str.strip()
    info["monkey_id"] = info["MonkeyID"].astype(str).str.strip()
    return info.drop_duplicates("sample_id").set_index("sample_id")


def wilson_ci(k: int, n: int, z: float = 1.959963984540054) -> tuple[float, float]:
    if n == 0:
        return float("nan"), float("nan")
    p = k / n
    denom = 1 + z * z / n
    center = (p + z * z / (2 * n)) / denom
    half = z * np.sqrt((p * (1 - p) + z * z / (4 * n)) / n) / denom
    return float(center - half), float(center + half)


def topk_from_scores(scores: np.ndarray, classes: np.ndarray, all_classes: list[str], k: int = 3) -> list[str]:
    score_map = {str(c): float(s) for c, s in zip(classes, scores)}
    full = np.asarray([score_map.get(c, -np.inf) for c in all_classes], dtype=float)
    order = np.argsort(full)[::-1][:k]
    return [all_classes[int(i)] for i in order]


def nearest_centroid_scores(pipe: Pipeline, x: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    selector = pipe.named_steps["select"]
    scaler = pipe.named_steps["scale"]
    model = pipe.named_steps["model"]
    xt = scaler.transform(selector.transform(x))
    centroids = model.centroids_
    distances = np.linalg.norm(xt[:, None, :] - centroids[None, :, :], axis=2)
    return -distances[0], model.classes_


def evaluate_lomo_baselines() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    vsd_path = ROOT / "bo2023 data" / "mfas5_819samples_23605genes_vsd4_rmbatch.xls"
    matrix = read_bo2023_gene_matrix(vsd_path)
    metadata = read_metadata()
    samples = [sample for sample in matrix.columns if sample in metadata.index]
    matrix = matrix.loc[:, samples]
    metadata = metadata.loc[samples]
    x = matrix.T.to_numpy(dtype=np.float32)
    y = metadata["network"].to_numpy(dtype=str)
    monkeys = metadata["monkey_id"].to_numpy(dtype=str)
    all_classes = sorted(pd.Series(y).unique())

    specs = [
        ModelSpec(
            "nearest_centroid_vsd_lomo",
            Pipeline(
                [
                    ("select", SelectKBest(f_classif, k=500)),
                    ("scale", StandardScaler()),
                    ("model", NearestCentroid()),
                ]
            ),
            "nearest_centroid",
        ),
        ModelSpec(
            "knn5_cosine_vsd_lomo",
            Pipeline(
                [
                    ("select", SelectKBest(f_classif, k=500)),
                    ("scale", StandardScaler()),
                    ("model", KNeighborsClassifier(n_neighbors=5, weights="distance", metric="cosine", algorithm="brute")),
                ]
            ),
            "predict_proba",
        ),
        ModelSpec(
            "random_forest_vsd_lomo",
            Pipeline(
                [
                    ("select", SelectKBest(f_classif, k=500)),
                    (
                        "model",
                        RandomForestClassifier(
                            n_estimators=300,
                            random_state=SEED,
                            class_weight="balanced_subsample",
                            n_jobs=-1,
                            min_samples_leaf=2,
                        ),
                    ),
                ]
            ),
            "predict_proba",
        ),
    ]
    detail_rows: list[dict[str, Any]] = []
    fold_rows: list[dict[str, Any]] = []
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for heldout_monkey in sorted(set(monkeys)):
            train = monkeys != heldout_monkey
            test = monkeys == heldout_monkey
            fold_rows.append(
                {
                    "heldout_monkey_id": heldout_monkey,
                    "n_train": int(train.sum()),
                    "n_test": int(test.sum()),
                    "train_networks": int(pd.Series(y[train]).nunique()),
                    "test_networks": int(pd.Series(y[test]).nunique()),
                }
            )
            for spec in specs:
                pipe = clone(spec.model)
                pipe.fit(x[train], y[train])
                for sample_id, truth, row in zip(np.asarray(samples)[test], y[test], x[test]):
                    xx = row.reshape(1, -1)
                    if spec.score_mode == "nearest_centroid":
                        scores, classes = nearest_centroid_scores(pipe, xx)
                    else:
                        scores = pipe.predict_proba(xx)[0]
                        classes = pipe.named_steps["model"].classes_
                    top = topk_from_scores(scores, classes, all_classes, 3)
                    detail_rows.append(
                        {
                            "model": spec.name,
                            "sample_id": sample_id,
                            "heldout_monkey_id": heldout_monkey,
                            "truth_network": truth,
                            "pred_top1": top[0],
                            "pred_top2": top[1],
                            "pred_top3": top[2],
                            "hit1": int(top[0] == truth),
                            "hit3": int(truth in top),
                            "n_selected_features": spec.k_features,
                            "validation_design": "leave-one-monkey-out; feature selection and model fitting use training monkeys only",
                        }
                    )
    detail = pd.DataFrame(detail_rows)
    metrics = []
    for model, frame in detail.groupby("model", sort=True):
        for metric, col in [("Top1", "hit1"), ("Top3", "hit3")]:
            hits = int(frame[col].sum())
            n = int(len(frame))
            lo, hi = wilson_ci(hits, n)
            metrics.append({"model": model, "metric": metric, "n": n, "hits": hits, "accuracy": hits / n, "wilson95_low": lo, "wilson95_high": hi})
    metrics_df = pd.DataFrame(metrics)
    f1_rows = []
    for model, frame in detail.groupby("model", sort=True):
        labels = all_classes
        precision, recall, f1, support = precision_recall_fscore_support(frame["truth_network"], frame["pred_top1"], labels=labels, zero_division=0)
        for label, p, r, f, s in zip(labels, precision, recall, f1, support):
            f1_rows.append({"model": model, "network": label, "support": int(s), "precision": float(p), "recall": float(r), "f1": float(f)})
        cm = pd.DataFrame(confusion_matrix(frame["truth_network"], frame["pred_top1"], labels=labels), index=labels, columns=labels)
        cm.to_csv(ML_DIR / f"{model}_top1_confusion.csv")
    return detail, metrics_df, pd.DataFrame(f1_rows).merge(pd.DataFrame(fold_rows), how="cross").iloc[:0] if False else pd.DataFrame(f1_rows), pd.DataFrame(fold_rows)


def write_tool_comparison() -> None:
    rows = [
        {
            "tool_or_family": "CIBERSORTx / cell-type deconvolution tools",
            "primary_goal": "Cell-type fraction inference from bulk expression",
            "truth_space": "Cell types",
            "comparison_to_cfRNA_BrainTrace": "Related expression deconvolution concept, but not a brain-region hierarchical source-ranking tool.",
        },
        {
            "tool_or_family": "TissueEnrich / tissue-expression enrichment tools",
            "primary_goal": "Tissue enrichment from gene sets or expression signatures",
            "truth_space": "Broad tissues",
            "comparison_to_cfRNA_BrainTrace": "Broad tissue-level inference; does not provide macaque Network -> resolution-group -> exact-region hierarchy.",
        },
        {
            "tool_or_family": "Cell-free transcriptome tissue-of-origin studies",
            "primary_goal": "Infer broad tissue or cell-type contribution to cfRNA",
            "truth_space": "Tissues/cell types",
            "comparison_to_cfRNA_BrainTrace": "Important biological context but generally not brain-region atlas candidate ranking with app/CLI output.",
        },
        {
            "tool_or_family": "Allen/AHBA atlas query workflows",
            "primary_goal": "Map genes or samples to human brain atlas annotations",
            "truth_space": "Human anatomical labels",
            "comparison_to_cfRNA_BrainTrace": "Atlas reference context; BrainTrace adds packaged hierarchical candidate ranking and validation boundary audits.",
        },
        {
            "tool_or_family": "BrainTrace",
            "primary_goal": "Hierarchical brain-origin candidate ranking and resolution-limit auditing",
            "truth_space": "Macaque atlas Network, resolution group, exploratory exact region",
            "comparison_to_cfRNA_BrainTrace": "Niche: explicit hierarchy, Top3 candidate beam, resolution-limit diagnostics, Streamlit and CLI from one scoring core.",
        },
    ]
    pd.DataFrame(rows).to_csv(TOOL_DIR / "conceptual_tool_comparison.csv", index=False)
    text = """# Same-field tool comparison

No strict like-for-like public tool was identified that performs the same three-level brain-region candidate-ranking task from RNA expression profiles. The closest methods are broad tissue-of-origin cfRNA models, cell-type deconvolution tools, tissue-enrichment workflows and atlas query pipelines.

BrainTrace's niche is narrower and more explicit: it returns hierarchical brain-origin candidates, audits the resolution limit of each output, and exposes the same scoring core through a Streamlit app and CLI. The tool should therefore be positioned as a resolution-aware candidate-ranking and audit tool, not as a generic tissue deconvolution or deterministic localization model.
"""
    (TOOL_DIR / "TOOL_COMPARISON_SUMMARY.md").write_text(text, encoding="utf-8")


def write_reproducibility_files() -> None:
    freeze = subprocess.run([sys.executable, "-m", "pip", "freeze"], check=False, capture_output=True, text=True)
    (ENG_DIR / "requirements-lock.txt").write_text(freeze.stdout, encoding="utf-8")
    packages = [line.strip() for line in freeze.stdout.splitlines() if line.strip() and " @ " not in line]
    env = ["name: cfrna-braintrace-v016", "channels:", "  - conda-forge", "dependencies:", f"  - python={sys.version_info.major}.{sys.version_info.minor}", "  - pip", "  - pip:"]
    env.extend([f"      - {pkg}" for pkg in packages])
    (ENG_DIR / "environment.yml").write_text("\n".join(env) + "\n", encoding="utf-8")
    seeds = {
        "p2_report_seed": SEED,
        "random_forest_random_state": SEED,
        "p0_weighted_random_seed": 20260629,
        "legacy_marker_validation_seed": 20260528,
        "network_confirmation_seed": 20260531,
        "notes": "LOMO folds are deterministic by MonkeyID; model fitting uses training folds only.",
    }
    (ENG_DIR / "random_seed_registry.json").write_text(json.dumps(seeds, indent=2), encoding="utf-8")
    sync = """# README / Zenodo / testing sync text

## README validation note

For publication completeness, the repository includes P0/P2 evidence packages under `reports/`. The P2 package contains leakage-controlled leave-one-monkey-out ML baselines, conceptual same-field tool comparison, pinned dependency exports, random seed registry and coverage artifacts where available.

## Zenodo release note

Archive the `reports/p0_hard_evidence_20260629/` and `reports/p2_publication_completeness_20260629/` directories with the manuscript-associated release. These directories document statistical inference, ML baselines, marker methodology, dependency locking and reproducibility checks.

## Testing note

Run `python -m pytest tests` for the standard test suite. If `coverage` is available, run `python -m coverage run -m pytest tests` followed by `python -m coverage report` and `python -m coverage xml -o reports/p2_publication_completeness_20260629/engineering_reproducibility/coverage.xml`.
"""
    (SYNC_DIR / "README_ZENODO_TESTING_SYNC_TEXT.md").write_text(sync, encoding="utf-8")


def write_summary(detail: pd.DataFrame, metrics: pd.DataFrame, f1: pd.DataFrame, folds: pd.DataFrame) -> None:
    detail.to_csv(ML_DIR / "simple_ml_lomo_detail.csv", index=False)
    metrics.to_csv(ML_DIR / "simple_ml_lomo_metrics.csv", index=False)
    f1.to_csv(ML_DIR / "simple_ml_lomo_class_f1.csv", index=False)
    folds.to_csv(ML_DIR / "simple_ml_lomo_fold_audit.csv", index=False)
    top3 = metrics[metrics["metric"].eq("Top3")].sort_values("accuracy", ascending=False)
    lines = [
        "# P2 publication completeness report",
        "",
        "## Simple ML baselines",
        "",
        "Design: leave-one-monkey-out Network classification on Bo2023 VSD expression. Each fold holds out one monkey. Feature selection (`SelectKBest(f_classif, k=500)`) and model fitting are performed only on training monkeys.",
        "",
        "| Model | Top1 | Top3 |",
        "|---|---:|---:|",
    ]
    for model in sorted(metrics["model"].unique()):
        m1 = metrics[(metrics["model"].eq(model)) & (metrics["metric"].eq("Top1"))].iloc[0]
        m3 = metrics[(metrics["model"].eq(model)) & (metrics["metric"].eq("Top3"))].iloc[0]
        lines.append(f"| {model} | {float(m1.accuracy):.2%} ({int(m1.hits)}/{int(m1.n)}) | {float(m3.accuracy):.2%} ({int(m3.hits)}/{int(m3.n)}) |")
    lines.extend(
        [
            "",
            f"Best Top3 baseline: `{top3.iloc[0].model}` at {float(top3.iloc[0].accuracy):.2%}. These baselines are comparison points only; they do not replace the hierarchical route because they do not return resolution groups, exact-region candidate rankings or resolution-limit diagnostics.",
            "",
            "## Same-field comparison",
            "",
            "No strict like-for-like public tool was identified. The comparison package therefore provides a conceptual table against tissue-of-origin, tissue-enrichment, cell-type deconvolution and atlas-query workflows.",
            "",
            "## Engineering reproducibility",
            "",
            "- `requirements-lock.txt` captures `pip freeze` from the current environment.",
            "- `environment.yml` provides a conda-style environment wrapper with pinned pip packages.",
            "- `random_seed_registry.json` records seeds used by P0/P2 scripts and legacy validation scripts.",
            "- Coverage artifacts should be generated with `coverage` when available; see `README_ZENODO_TESTING_SYNC_TEXT.md`.",
        ]
    )
    (OUTDIR / "P2_PUBLICATION_COMPLETENESS_REPORT.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    coverage_done = (ENG_DIR / "coverage.xml").exists() and (ENG_DIR / "coverage_report.txt").exists()
    checklist = f"""# P2 completion checklist

- [x] Simple ML baselines: k-NN, random forest, nearest-centroid under no-leak LOMO design.
- [x] Same-field / conceptual tool comparison.
- [x] Tool niche statement: hierarchical candidate ranking, resolution-limit auditing, shared app/CLI scoring core.
- [x] Pinned `requirements-lock.txt`.
- [x] `environment.yml`.
- [x] Random seed registry.
- [x] README / Zenodo / testing sync text.
- [{'x' if coverage_done else ' '}] Coverage report generated. See `engineering_reproducibility/coverage_report.txt` and `coverage.xml`.
"""
    (OUTDIR / "P2_COMPLETION_CHECKLIST.md").write_text(checklist, encoding="utf-8")


def main() -> int:
    for d in [OUTDIR, ML_DIR, ENG_DIR, TOOL_DIR, SYNC_DIR]:
        d.mkdir(parents=True, exist_ok=True)
    detail, metrics, f1, folds = evaluate_lomo_baselines()
    write_tool_comparison()
    write_reproducibility_files()
    write_summary(detail, metrics, f1, folds)
    manifest = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "script": "scripts/generate_p2_publication_completeness.py",
        "files": sorted(str(p.relative_to(OUTDIR)) for p in OUTDIR.rglob("*") if p.is_file()),
    }
    (OUTDIR / "manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(f"Wrote P2 publication-completeness package to {OUTDIR}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
