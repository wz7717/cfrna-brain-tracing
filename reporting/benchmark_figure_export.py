from __future__ import annotations

import io
import json
import math
import os
import tempfile
import zipfile
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np
import pandas as pd

from core.methods import method_label


def _import_matplotlib():
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        from matplotlib.backends.backend_pdf import PdfPages
        return matplotlib, plt, PdfPages
    except Exception as e:
        raise RuntimeError(
            "matplotlib is required for benchmark figure/PDF export. "
            "Please install dependencies with `pip install -r requirements.txt` or `pip install matplotlib`."
        ) from e


def _clean_float(x: Any, default: float = 0.0) -> float:
    try:
        v = float(x)
        return v if math.isfinite(v) else default
    except Exception:
        return default


def _safe_text(x: Any) -> str:
    if x is None:
        return ''
    if isinstance(x, float) and not math.isfinite(x):
        return ''
    return str(x)


def _json_safe(obj: Any) -> Any:
    if isinstance(obj, dict):
        return {str(k): _json_safe(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_json_safe(v) for v in obj]
    try:
        if pd.isna(obj):
            return None
    except Exception:
        pass
    if isinstance(obj, (np.floating, float)):
        return None if not math.isfinite(float(obj)) else float(obj)
    if isinstance(obj, (np.integer, int)):
        return int(obj)
    return obj


def _metric_value(metrics_df: pd.DataFrame, metric: str, default: Any = None) -> Any:
    if metrics_df is None or metrics_df.empty or "metric" not in metrics_df.columns:
        return default
    hit = metrics_df.loc[metrics_df["metric"] == metric, "value"]
    return default if hit.empty else hit.iloc[0]


def build_benchmark_narrative(suite: Dict[str, pd.DataFrame], metadata: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
    metadata = metadata or {}
    detail_df = suite.get("detail_df", pd.DataFrame())
    metrics_df = suite.get("metrics_df", pd.DataFrame())
    confusion_long = suite.get("confusion_long_df", pd.DataFrame())
    summary = detail_df.attrs.get("summary", {}) if hasattr(detail_df, "attrs") else {}
    wrong = pd.DataFrame()
    if confusion_long is not None and not confusion_long.empty:
        wrong = confusion_long[confusion_long["truth_region"].astype(str) != confusion_long["pred_region"].astype(str)].copy()
        wrong = wrong.sort_values(["count", "row_fraction"], ascending=False).head(5)
    top1 = summary.get("top1_acc", _metric_value(metrics_df, "Top1_acc_valid"))
    auc = summary.get("auc", _metric_value(metrics_df, "MacroAUC_ovr_valid"))
    abstain = summary.get("abstain_rate", _metric_value(metrics_df, "Abstain_rate"))
    method = metadata.get("method", "ensemble")
    recommendation = "Use Multi-signal ensemble as the primary result; report NNLS simplex as fraction/CI support when needed."
    if method == "nnls_simplex":
        recommendation = "Use NNLS simplex for fraction interpretation; add ensemble as the primary classifier when publishing accuracy claims."
    elif method == "correlation":
        recommendation = "Use correlation for screening; add ensemble or NNLS simplex for final reporting."
    return {
        "methods_text": (
            "Samples with ground-truth labels were evaluated with top-k accuracy, row-normalized confusion matrices, "
            "one-vs-rest ROC/AUC, prediction-space separability, and bootstrap stability summaries. "
            "Low-overlap or insufficient-evidence samples may be abstained before metric calculation."
        ),
        "best_parameter_summary": {
            "method": method,
            "method_label": method_label(method),
            "atlas_id": metadata.get("atlas_id"),
            "sigset_id": metadata.get("sigset_id"),
            "use_value": metadata.get("use_value"),
            "l2": metadata.get("l2"),
            "ensemble_alpha": metadata.get("ensemble_alpha"),
            "bootstrap_n": metadata.get("bootstrap_n"),
            "bootstrap_gene_frac": metadata.get("bootstrap_gene_frac"),
        },
        "performance_summary": {
            "top1_acc": _json_safe(top1),
            "macro_auc": _json_safe(auc),
            "abstain_rate": _json_safe(abstain),
        },
        "confusable_regions": wrong.to_dict(orient="records") if not wrong.empty else [],
        "recommendation": recommendation,
    }


def _save_figure(fig, path_base: Path, dpi: int = 300) -> List[str]:
    _, plt, _ = _import_matplotlib()
    out = []
    png_path = f"{path_base}.png"
    pdf_path = f"{path_base}.pdf"
    fig.savefig(png_path, dpi=dpi, bbox_inches='tight')
    fig.savefig(pdf_path, bbox_inches='tight')
    plt.close(fig)
    out.extend([png_path, pdf_path])
    return out


def _plot_table(ax, df: pd.DataFrame, title: str, nrows: int = 8):
    ax.axis('off')
    ax.set_title(title)
    if df is None or df.empty:
        ax.text(0.5, 0.5, 'No data', ha='center', va='center')
        return
    d = df.head(nrows).copy()
    for c in d.columns:
        if pd.api.types.is_float_dtype(d[c]):
            d[c] = d[c].map(lambda x: '' if pd.isna(x) else f'{float(x):.3f}')
    table = ax.table(cellText=d.values, colLabels=d.columns, loc='center')
    table.auto_set_font_size(False)
    table.set_fontsize(8)
    table.scale(1, 1.2)


def _plot_heatmap(ax, df: pd.DataFrame, title: str, fmt: str = '.2f', cmap: str = 'viridis'):
    _, plt, _ = _import_matplotlib()
    ax.set_title(title)
    if df is None or df.empty:
        ax.text(0.5, 0.5, 'No data', ha='center', va='center')
        ax.axis('off')
        return
    arr = df.to_numpy(dtype=float)
    im = ax.imshow(arr, aspect='auto', cmap=cmap)
    ax.set_xticks(range(df.shape[1]))
    ax.set_xticklabels(df.columns, rotation=45, ha='right', fontsize=8)
    ax.set_yticks(range(df.shape[0]))
    ax.set_yticklabels(df.index, fontsize=8)
    for i in range(df.shape[0]):
        for j in range(df.shape[1]):
            val = arr[i, j]
            txt = format(val, fmt) if math.isfinite(float(val)) else ''
            ax.text(j, i, txt, ha='center', va='center', fontsize=7)
    plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    ax.set_xlabel('Predicted region')
    ax.set_ylabel('True region')


def _metric_map(metrics_df: pd.DataFrame) -> Dict[str, Any]:
    if metrics_df is None or metrics_df.empty or not {'metric', 'value'}.issubset(metrics_df.columns):
        return {}
    return {str(k): v for k, v in zip(metrics_df['metric'], metrics_df['value'])}


def _plot_metrics_overview(metrics_df: pd.DataFrame, summary: Dict[str, Any], path_base: Path) -> List[str]:
    _, plt, _ = _import_matplotlib()
    fig = plt.figure(figsize=(12, 8))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.3])
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[1, :])

    m = _metric_map(metrics_df)
    keys = ['Top1_acc_valid', 'Balanced_acc_valid', 'MacroAUC_ovr_valid', 'Abstain_rate']
    labels = ['Top1', 'Balanced', 'MacroAUC', 'Abstain']
    vals = [_clean_float(m.get(k, np.nan), 0.0) for k in keys]
    ax1.bar(labels, vals)
    ax1.set_ylim(0, 1)
    ax1.set_title('Figure 1A. Benchmark performance overview')
    for i, v in enumerate(vals):
        ax1.text(i, min(1.0, v + 0.03), f'{v:.3f}', ha='center', fontsize=9)

    ax2.axis('off')
    txt_lines = [
        'Figure 1B. Run summary',
        f"N total: {_safe_text(summary.get('n_total', m.get('N_samples_total', '')))}",
        f"N valid: {_safe_text(summary.get('n_valid', m.get('N_samples_valid', '')))}",
        f"N classes: {_safe_text(summary.get('n_classes', m.get('N_classes', '')))}",
        f"Mean top1 confidence: {_clean_float(summary.get('mean_top1_confidence', np.nan), np.nan):.3f}" if summary.get('mean_top1_confidence') is not None else 'Mean top1 confidence: ',
        f"Mean decision margin: {_clean_float(summary.get('mean_decision_margin', np.nan), np.nan):.3f}" if summary.get('mean_decision_margin') is not None else 'Mean decision margin: ',
        f"Mean top1 stability: {_clean_float(summary.get('mean_top1_stability', np.nan), np.nan):.3f}" if summary.get('mean_top1_stability') is not None else 'Mean top1 stability: ',
    ]
    strat = summary.get('strat', {}) if isinstance(summary, dict) else {}
    if isinstance(strat, dict) and strat:
        txt_lines.append('Stability-stratified accuracy:')
        for k, v in strat.items():
            if v is None:
                txt_lines.append(f'  {k}: NA')
            else:
                vv = _clean_float(v, np.nan)
                txt_lines.append(f'  {k}: {vv:.3f}' if math.isfinite(vv) else f'  {k}: NA')
    ax2.text(0.0, 1.0, '\n'.join(txt_lines), va='top', ha='left', fontsize=11)
    _plot_table(ax3, metrics_df if metrics_df is not None else pd.DataFrame(), 'Figure 1C. Core metric table', nrows=14)
    fig.suptitle('Paper-grade Benchmark Figure 1', fontsize=14)
    return _save_figure(fig, path_base)


def _plot_confusion(conf_raw: pd.DataFrame, conf_norm: pd.DataFrame, path_base: Path) -> List[str]:
    _, plt, _ = _import_matplotlib()
    fig, axes = plt.subplots(1, 2, figsize=(15, 6))
    _plot_heatmap(axes[0], conf_norm, 'Figure 2A. Row-normalized confusion matrix', fmt='.2f', cmap='viridis')
    _plot_heatmap(axes[1], conf_raw, 'Figure 2B. Raw confusion counts', fmt='.0f', cmap='magma')
    fig.suptitle('Paper-grade Benchmark Figure 2', fontsize=14)
    return _save_figure(fig, path_base)


def _plot_roc(roc_curve_df: pd.DataFrame, roc_summary_df: pd.DataFrame, path_base: Path) -> List[str]:
    _, plt, _ = _import_matplotlib()
    fig = plt.figure(figsize=(13, 6))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.5, 1])
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax1.plot([0, 1], [0, 1], linestyle='--', linewidth=1)
    if roc_curve_df is not None and not roc_curve_df.empty:
        for region, g in roc_curve_df.groupby('region_id'):
            auc_val = None
            if 'auc' in g.columns and not g['auc'].isna().all():
                auc_val = _clean_float(g['auc'].dropna().iloc[0], np.nan)
            label = f'{region} (AUC={auc_val:.3f})' if auc_val is not None and math.isfinite(auc_val) else str(region)
            ax1.plot(g['fpr'], g['tpr'], label=label, linewidth=1.8)
        ax1.legend(fontsize=8, loc='lower right')
    else:
        ax1.text(0.5, 0.5, 'No ROC data', ha='center', va='center')
    ax1.set_title('Figure 3A. One-vs-rest ROC curves')
    ax1.set_xlabel('False positive rate')
    ax1.set_ylabel('True positive rate')
    ax1.set_xlim(0, 1)
    ax1.set_ylim(0, 1)

    if roc_summary_df is not None and not roc_summary_df.empty:
        rs = roc_summary_df.copy().sort_values('auc', ascending=False, na_position='last')
        ax2.barh(rs['region_id'], rs['auc'].fillna(0.0))
        ax2.invert_yaxis()
        ax2.set_xlim(0, 1)
        ax2.set_title('Figure 3B. Region-wise AUC summary')
        ax2.set_xlabel('AUC')
    else:
        ax2.axis('off')
        ax2.text(0.5, 0.5, 'No AUC summary', ha='center', va='center')
    fig.suptitle('Paper-grade Benchmark Figure 3', fontsize=14)
    return _save_figure(fig, path_base)


def _plot_separability_and_stability(separability_df: pd.DataFrame, stability_region_df: pd.DataFrame, stability_bin_df: pd.DataFrame, path_base: Path) -> List[str]:
    _, plt, _ = _import_matplotlib()
    fig = plt.figure(figsize=(14, 10))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.3, 1])
    ax1 = fig.add_subplot(gs[0, :])
    ax2 = fig.add_subplot(gs[1, 0])
    ax3 = fig.add_subplot(gs[1, 1])

    if separability_df is not None and not separability_df.empty:
        labels = separability_df['label'].astype(str).fillna('NA')
        uniq = list(dict.fromkeys(labels.tolist()))
        cmap = plt.get_cmap('tab10')
        color_map = {lab: cmap(i % 10) for i, lab in enumerate(uniq)}
        for lab, g in separability_df.groupby('label'):
            ax1.scatter(g['pc1'], g['pc2'], label=str(lab), alpha=0.8, s=40, c=[color_map.get(str(lab), 'grey')])
        ax1.legend(fontsize=8, ncol=2)
        pc1_var = _clean_float(separability_df.get('pc1_var', pd.Series([np.nan])).iloc[0], np.nan)
        pc2_var = _clean_float(separability_df.get('pc2_var', pd.Series([np.nan])).iloc[0], np.nan)
        ax1.set_xlabel(f"PC1 ({pc1_var:.1%})" if math.isfinite(pc1_var) else 'PC1')
        ax1.set_ylabel(f"PC2 ({pc2_var:.1%})" if math.isfinite(pc2_var) else 'PC2')
        ax1.set_title('Figure 4A. Brain region separability in prediction space')
    else:
        ax1.axis('off')
        ax1.text(0.5, 0.5, 'No separability data', ha='center', va='center')

    if stability_region_df is not None and not stability_region_df.empty:
        s = stability_region_df.copy().sort_values('top1_acc_valid', ascending=False, na_position='last')
        ax2.barh(s['region_id'], s['top1_acc_valid'].fillna(0.0))
        ax2.invert_yaxis()
        ax2.set_xlim(0, 1)
        ax2.set_title('Figure 4B. Region-wise Top1 accuracy')
        ax2.set_xlabel('Top1 accuracy')
    else:
        ax2.axis('off')
        ax2.text(0.5, 0.5, 'No region-wise stability summary', ha='center', va='center')

    if stability_bin_df is not None and not stability_bin_df.empty:
        b = stability_bin_df.copy()
        ax3.bar(b['stability_bin'], b['top1_acc'].fillna(0.0))
        ax3.set_ylim(0, 1)
        ax3.set_title('Figure 4C. Accuracy by stability bin')
        ax3.set_ylabel('Top1 accuracy')
        ax3.tick_params(axis='x', rotation=15)
        for i, row in b.reset_index(drop=True).iterrows():
            n = row.get('n_samples', '')
            val = _clean_float(row.get('top1_acc', np.nan), np.nan)
            if math.isfinite(val):
                ax3.text(i, min(1.0, val + 0.03), f'n={int(n) if pd.notna(n) else 0}', ha='center', fontsize=8)
    else:
        ax3.axis('off')
        ax3.text(0.5, 0.5, 'No stability bin summary', ha='center', va='center')
    fig.suptitle('Paper-grade Benchmark Figure 4', fontsize=14)
    return _save_figure(fig, path_base)


def _plot_region_reliability(stability_region_df: pd.DataFrame, path_base: Path) -> List[str]:
    _, plt, _ = _import_matplotlib()
    fig = plt.figure(figsize=(13, 6))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.2, 1])
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    if stability_region_df is not None and not stability_region_df.empty:
        s = stability_region_df.copy().sort_values('mean_top1_confidence', ascending=False, na_position='last')
        ax1.barh(s['region_id'], s['mean_top1_confidence'].fillna(0.0), label='Mean confidence')
        if 'mean_decision_margin' in s.columns:
            ax1.scatter(s['mean_decision_margin'].fillna(0.0), np.arange(len(s)), marker='o', label='Mean margin')
        ax1.invert_yaxis()
        ax1.set_xlim(0, 1)
        ax1.set_title('Figure 5A. Region-wise confidence and margin')
        ax1.set_xlabel('Score')
        ax1.legend(fontsize=8)

        s2 = stability_region_df.copy().sort_values('abstain_rate', ascending=False, na_position='last')
        ax2.barh(s2['region_id'], s2['abstain_rate'].fillna(0.0))
        ax2.invert_yaxis()
        ax2.set_xlim(0, 1)
        ax2.set_title('Figure 5B. Region-wise abstain rate')
        ax2.set_xlabel('Abstain rate')
    else:
        for ax in [ax1, ax2]:
            ax.axis('off')
            ax.text(0.5, 0.5, 'No reliability data', ha='center', va='center')
    fig.suptitle('Paper-grade Benchmark Figure 5', fontsize=14)
    return _save_figure(fig, path_base)


def _plot_overlap_failure_modes(detail_df: pd.DataFrame, path_base: Path) -> List[str]:
    _, plt, _ = _import_matplotlib()
    fig = plt.figure(figsize=(14, 10))
    gs = fig.add_gridspec(2, 2)
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[1, :])

    if detail_df is not None and not detail_df.empty:
        d = detail_df.copy()
        if 'hit1' not in d.columns:
            d['hit1'] = np.nan
        if 'abstained' not in d.columns:
            d['abstained'] = 0
        if 'overlap_genes' not in d.columns:
            d['overlap_genes'] = np.nan
        if 'top1_confidence' not in d.columns:
            d['top1_confidence'] = np.nan
        state = np.where(d['abstained'] == 1, 'abstained', np.where(d['hit1'] == 1, 'correct', 'incorrect'))
        d['prediction_state'] = state
        order = ['correct', 'incorrect', 'abstained']
        box_data = [d.loc[d['prediction_state'] == s, 'overlap_genes'].dropna().to_numpy() for s in order]
        ax1.boxplot(box_data, tick_labels=order, showfliers=False)
        ax1.set_title('Figure 6A. Overlap genes by prediction state')
        ax1.set_ylabel('Overlap genes')

        counts = d['traceability'].fillna('unknown').astype(str).value_counts()
        ax2.bar(counts.index.tolist(), counts.values.tolist())
        ax2.set_title('Figure 6B. Traceability/failure mode counts')
        ax2.tick_params(axis='x', rotation=20)

        for state_name, g in d.groupby('prediction_state'):
            ax3.scatter(g['overlap_genes'], g['top1_confidence'], label=state_name, alpha=0.7, s=35)
        ax3.set_title('Figure 6C. Overlap vs confidence')
        ax3.set_xlabel('Overlap genes')
        ax3.set_ylabel('Top1 confidence')
        ax3.legend(fontsize=8)
    else:
        for ax in [ax1, ax2, ax3]:
            ax.axis('off')
            ax.text(0.5, 0.5, 'No detail data', ha='center', va='center')
    fig.suptitle('Paper-grade Benchmark Figure 6', fontsize=14)
    return _save_figure(fig, path_base)


def _suite_tables(suite: Dict[str, pd.DataFrame]) -> Dict[str, pd.DataFrame]:
    return {
        'benchmark_detail.csv': suite.get('detail_df', pd.DataFrame()),
        'benchmark_metrics.csv': suite.get('metrics_df', pd.DataFrame()),
        'benchmark_probability_matrix.csv': suite.get('probability_df', pd.DataFrame()),
        'benchmark_confusion_raw.csv': suite.get('confusion_raw_df', pd.DataFrame()),
        'benchmark_confusion_norm.csv': suite.get('confusion_norm_df', pd.DataFrame()),
        'benchmark_confusion_long.csv': suite.get('confusion_long_df', pd.DataFrame()),
        'benchmark_roc_curve.csv': suite.get('roc_curve_df', pd.DataFrame()),
        'benchmark_roc_summary.csv': suite.get('roc_summary_df', pd.DataFrame()),
        'benchmark_separability.csv': suite.get('separability_df', pd.DataFrame()),
        'benchmark_centroids.csv': suite.get('centroid_df', pd.DataFrame()),
        'benchmark_centroid_distance.csv': suite.get('centroid_distance_df', pd.DataFrame()),
        'benchmark_stability_region.csv': suite.get('stability_region_df', pd.DataFrame()),
        'benchmark_stability_bin.csv': suite.get('stability_bin_df', pd.DataFrame()),
    }


def export_benchmark_paper_figures(output_dir: str | os.PathLike[str], suite: Dict[str, pd.DataFrame], metadata: Optional[Dict[str, Any]] = None, prefix: str = 'benchmark_papergrade') -> Dict[str, Any]:
    outdir = Path(output_dir)
    figures_dir = outdir / 'figures'
    tables_dir = outdir / 'tables'
    figures_dir.mkdir(parents=True, exist_ok=True)
    tables_dir.mkdir(parents=True, exist_ok=True)

    detail_df = suite.get('detail_df', pd.DataFrame())
    metrics_df = suite.get('metrics_df', pd.DataFrame())
    conf_raw = suite.get('confusion_raw_df', pd.DataFrame())
    conf_norm = suite.get('confusion_norm_df', pd.DataFrame())
    roc_curve_df = suite.get('roc_curve_df', pd.DataFrame())
    roc_summary_df = suite.get('roc_summary_df', pd.DataFrame())
    separability_df = suite.get('separability_df', pd.DataFrame())
    stability_region_df = suite.get('stability_region_df', pd.DataFrame())
    stability_bin_df = suite.get('stability_bin_df', pd.DataFrame())
    summary = detail_df.attrs.get('summary', {}) if hasattr(detail_df, 'attrs') else {}

    figure_files: List[str] = []
    figure_files += _plot_metrics_overview(metrics_df, summary, figures_dir / f'{prefix}_Figure1_overview')
    figure_files += _plot_confusion(conf_raw, conf_norm, figures_dir / f'{prefix}_Figure2_confusion')
    figure_files += _plot_roc(roc_curve_df, roc_summary_df, figures_dir / f'{prefix}_Figure3_roc')
    figure_files += _plot_separability_and_stability(separability_df, stability_region_df, stability_bin_df, figures_dir / f'{prefix}_Figure4_separability_stability')
    figure_files += _plot_region_reliability(stability_region_df, figures_dir / f'{prefix}_Figure5_region_reliability')
    figure_files += _plot_overlap_failure_modes(detail_df, figures_dir / f'{prefix}_Figure6_failure_modes')

    table_files: List[str] = []
    for name, df in _suite_tables(suite).items():
        if df is None or df.empty:
            continue
        path = tables_dir / name
        if isinstance(df.index, pd.Index) and df.index.name is not None:
            df.to_csv(path, encoding='utf-8-sig')
        else:
            df.to_csv(path, index=False, encoding='utf-8-sig')
        table_files.append(str(path))

    narrative = build_benchmark_narrative(suite, metadata)
    manifest = {
        'created_at': datetime.utcnow().isoformat() + 'Z',
        'prefix': prefix,
        'metadata': _json_safe(metadata or {}),
        'summary': _json_safe(summary),
        'narrative': _json_safe(narrative),
        'figure_files': [str(Path(p).relative_to(outdir)) for p in figure_files],
        'table_files': [str(Path(p).relative_to(outdir)) for p in table_files],
    }
    with open(outdir / 'manifest.json', 'w', encoding='utf-8') as f:
        json.dump(manifest, f, ensure_ascii=False, indent=2)
    readme = (
        'Paper-grade benchmark figure export\n\n'
        'This bundle contains automatically generated Figure1-Figure6 style outputs for the benchmark run.\n\n'
        'Figures\n'
        '- Figure1: overall metric overview and summary table\n'
        '- Figure2: confusion matrix (normalized and raw)\n'
        '- Figure3: one-vs-rest ROC curves and AUC summary\n'
        '- Figure4: brain region separability and stability analyses\n'
        '- Figure5: region-wise reliability (confidence, margin, abstain)\n'
        '- Figure6: failure modes and overlap-quality analyses\n'
        '\nMethods\n'
        f"{narrative['methods_text']}\n\n"
        'Recommendation\n'
        f"{narrative['recommendation']}\n"
    )
    with open(outdir / 'README_figures.txt', 'w', encoding='utf-8') as f:
        f.write(readme)
    return manifest


def build_benchmark_figure_bundle_bytes(suite: Dict[str, pd.DataFrame], metadata: Optional[Dict[str, Any]] = None, prefix: str = 'benchmark_papergrade') -> bytes:
    with tempfile.TemporaryDirectory(prefix='bench_fig_export_') as tmpdir:
        export_benchmark_paper_figures(tmpdir, suite=suite, metadata=metadata, prefix=prefix)
        bio = io.BytesIO()
        with zipfile.ZipFile(bio, 'w', compression=zipfile.ZIP_DEFLATED) as zf:
            for root, _, files in os.walk(tmpdir):
                for fn in files:
                    full = os.path.join(root, fn)
                    arc = os.path.relpath(full, tmpdir)
                    zf.write(full, arcname=arc)
        return bio.getvalue()


def export_benchmark_report_pdf(output_pdf_path: str | os.PathLike[str], suite: Dict[str, pd.DataFrame], metadata: Optional[Dict[str, Any]] = None, title: str = 'cfRNA Source Tracing Benchmark Report') -> str:
    _, plt, PdfPages = _import_matplotlib()
    output_pdf_path = str(output_pdf_path)
    os.makedirs(os.path.dirname(os.path.abspath(output_pdf_path)), exist_ok=True)

    detail_df = suite.get('detail_df', pd.DataFrame())
    metrics_df = suite.get('metrics_df', pd.DataFrame())
    conf_raw = suite.get('confusion_raw_df', pd.DataFrame())
    conf_norm = suite.get('confusion_norm_df', pd.DataFrame())
    roc_curve_df = suite.get('roc_curve_df', pd.DataFrame())
    roc_summary_df = suite.get('roc_summary_df', pd.DataFrame())
    separability_df = suite.get('separability_df', pd.DataFrame())
    stability_region_df = suite.get('stability_region_df', pd.DataFrame())
    stability_bin_df = suite.get('stability_bin_df', pd.DataFrame())
    summary = detail_df.attrs.get('summary', {}) if hasattr(detail_df, 'attrs') else {}
    meta = _json_safe(metadata or {})
    narrative = build_benchmark_narrative(suite, metadata)

    with PdfPages(output_pdf_path) as pdf:
        fig = plt.figure(figsize=(11.69, 8.27))
        ax = fig.add_subplot(111)
        ax.axis('off')
        lines = [title, '', 'Metadata']
        for k, v in meta.items():
            lines.append(f'- {k}: {v}')
        lines += ['', 'Summary']
        for k, v in _json_safe(summary).items():
            lines.append(f'- {k}: {v}')
        ax.text(0.03, 0.97, '\n'.join(lines), va='top', ha='left', fontsize=12)
        pdf.savefig(fig, bbox_inches='tight')
        plt.close(fig)

        methods_page = plt.figure(figsize=(11.69, 8.27))
        axm = methods_page.add_subplot(111)
        axm.axis('off')
        method_lines = [
            'Methods and interpretation',
            '',
            narrative['methods_text'],
            '',
            'Best parameter combination used in this run:',
        ]
        for k, v in narrative['best_parameter_summary'].items():
            method_lines.append(f'- {k}: {v}')
        method_lines += [
            '',
            'Recommendation:',
            narrative['recommendation'],
        ]
        if narrative.get('confusable_regions'):
            method_lines += ['', 'Most confusable region pairs:']
            for row in narrative['confusable_regions'][:5]:
                method_lines.append(
                    f"- {row.get('truth_region')} -> {row.get('pred_region')}: "
                    f"count={row.get('count')}, row_fraction={row.get('row_fraction')}"
                )
        axm.text(0.04, 0.96, '\n'.join(method_lines), va='top', ha='left', fontsize=11, wrap=True)
        pdf.savefig(methods_page, bbox_inches='tight')
        plt.close(methods_page)

        for builder in [
            lambda: _plot_metrics_overview(metrics_df, summary, Path(tempfile.gettempdir()) / 'ignore')[0:0],
        ]:
            pass

        # Re-create each figure as a live figure and append it to the PDF report.
        # Figure 1
        fig1 = plt.figure(figsize=(12, 8))
        gs = fig1.add_gridspec(2, 2, height_ratios=[1.0, 1.3])
        ax1 = fig1.add_subplot(gs[0, 0]); ax2 = fig1.add_subplot(gs[0, 1]); ax3 = fig1.add_subplot(gs[1, :])
        m = _metric_map(metrics_df)
        keys = ['Top1_acc_valid', 'Balanced_acc_valid', 'MacroAUC_ovr_valid', 'Abstain_rate']
        labels = ['Top1', 'Balanced', 'MacroAUC', 'Abstain']
        vals = [_clean_float(m.get(k, np.nan), 0.0) for k in keys]
        ax1.bar(labels, vals); ax1.set_ylim(0, 1); ax1.set_title('Figure 1A. Benchmark performance overview')
        for i, v in enumerate(vals):
            ax1.text(i, min(1.0, v + 0.03), f'{v:.3f}', ha='center', fontsize=9)
        ax2.axis('off')
        lines = ['Figure 1B. Run summary'] + [f'{k}: {v}' for k, v in _json_safe(summary).items()]
        ax2.text(0.0, 1.0, '\n'.join(lines[:18]), va='top', ha='left', fontsize=10)
        _plot_table(ax3, metrics_df if metrics_df is not None else pd.DataFrame(), 'Figure 1C. Core metric table', nrows=14)
        fig1.suptitle('Paper-grade Benchmark Figure 1', fontsize=14)
        pdf.savefig(fig1, bbox_inches='tight'); plt.close(fig1)

        # Figure 2
        fig2, axes = plt.subplots(1, 2, figsize=(15, 6))
        _plot_heatmap(axes[0], conf_norm, 'Figure 2A. Row-normalized confusion matrix', fmt='.2f', cmap='viridis')
        _plot_heatmap(axes[1], conf_raw, 'Figure 2B. Raw confusion counts', fmt='.0f', cmap='magma')
        fig2.suptitle('Paper-grade Benchmark Figure 2', fontsize=14)
        pdf.savefig(fig2, bbox_inches='tight'); plt.close(fig2)

        # Figure 3
        fig3 = plt.figure(figsize=(13, 6))
        gs3 = fig3.add_gridspec(1, 2, width_ratios=[1.5, 1])
        ax31 = fig3.add_subplot(gs3[0, 0]); ax32 = fig3.add_subplot(gs3[0, 1])
        ax31.plot([0, 1], [0, 1], linestyle='--', linewidth=1)
        if roc_curve_df is not None and not roc_curve_df.empty:
            for region, g in roc_curve_df.groupby('region_id'):
                auc_val = None
                if 'auc' in g.columns and not g['auc'].isna().all():
                    auc_val = _clean_float(g['auc'].dropna().iloc[0], np.nan)
                label = f'{region} (AUC={auc_val:.3f})' if auc_val is not None and math.isfinite(auc_val) else str(region)
                ax31.plot(g['fpr'], g['tpr'], label=label, linewidth=1.8)
            ax31.legend(fontsize=7, loc='lower right')
        else:
            ax31.text(0.5, 0.5, 'No ROC data', ha='center', va='center')
        ax31.set_title('Figure 3A. One-vs-rest ROC curves'); ax31.set_xlim(0, 1); ax31.set_ylim(0, 1)
        ax31.set_xlabel('False positive rate'); ax31.set_ylabel('True positive rate')
        if roc_summary_df is not None and not roc_summary_df.empty:
            rs = roc_summary_df.copy().sort_values('auc', ascending=False, na_position='last')
            ax32.barh(rs['region_id'], rs['auc'].fillna(0.0)); ax32.invert_yaxis(); ax32.set_xlim(0, 1)
            ax32.set_title('Figure 3B. Region-wise AUC summary'); ax32.set_xlabel('AUC')
        else:
            ax32.axis('off'); ax32.text(0.5, 0.5, 'No AUC summary', ha='center', va='center')
        fig3.suptitle('Paper-grade Benchmark Figure 3', fontsize=14)
        pdf.savefig(fig3, bbox_inches='tight'); plt.close(fig3)

        # Figure 4
        fig4 = plt.figure(figsize=(14, 10))
        gs4 = fig4.add_gridspec(2, 2, height_ratios=[1.3, 1])
        ax41 = fig4.add_subplot(gs4[0, :]); ax42 = fig4.add_subplot(gs4[1, 0]); ax43 = fig4.add_subplot(gs4[1, 1])
        if separability_df is not None and not separability_df.empty:
            labels4 = separability_df['label'].astype(str).fillna('NA')
            uniq = list(dict.fromkeys(labels4.tolist())); cmap = plt.get_cmap('tab10'); color_map = {lab: cmap(i % 10) for i, lab in enumerate(uniq)}
            for lab, g in separability_df.groupby('label'):
                ax41.scatter(g['pc1'], g['pc2'], label=str(lab), alpha=0.8, s=40, c=[color_map.get(str(lab), 'grey')])
            ax41.legend(fontsize=8, ncol=2)
            pc1_var = _clean_float(separability_df.get('pc1_var', pd.Series([np.nan])).iloc[0], np.nan)
            pc2_var = _clean_float(separability_df.get('pc2_var', pd.Series([np.nan])).iloc[0], np.nan)
            ax41.set_xlabel(f"PC1 ({pc1_var:.1%})" if math.isfinite(pc1_var) else 'PC1')
            ax41.set_ylabel(f"PC2 ({pc2_var:.1%})" if math.isfinite(pc2_var) else 'PC2')
            ax41.set_title('Figure 4A. Brain region separability in prediction space')
        else:
            ax41.axis('off'); ax41.text(0.5, 0.5, 'No separability data', ha='center', va='center')
        if stability_region_df is not None and not stability_region_df.empty:
            s = stability_region_df.copy().sort_values('top1_acc_valid', ascending=False, na_position='last')
            ax42.barh(s['region_id'], s['top1_acc_valid'].fillna(0.0)); ax42.invert_yaxis(); ax42.set_xlim(0, 1)
            ax42.set_title('Figure 4B. Region-wise Top1 accuracy'); ax42.set_xlabel('Top1 accuracy')
        else:
            ax42.axis('off'); ax42.text(0.5, 0.5, 'No region-wise stability summary', ha='center', va='center')
        if stability_bin_df is not None and not stability_bin_df.empty:
            b = stability_bin_df.copy(); ax43.bar(b['stability_bin'], b['top1_acc'].fillna(0.0)); ax43.set_ylim(0, 1)
            ax43.set_title('Figure 4C. Accuracy by stability bin'); ax43.set_ylabel('Top1 accuracy'); ax43.tick_params(axis='x', rotation=15)
        else:
            ax43.axis('off'); ax43.text(0.5, 0.5, 'No stability bin summary', ha='center', va='center')
        fig4.suptitle('Paper-grade Benchmark Figure 4', fontsize=14)
        pdf.savefig(fig4, bbox_inches='tight'); plt.close(fig4)

        # Figure 5 and 6 via reusable builders
        for fig_func in [
            lambda pb: _plot_region_reliability(stability_region_df, pb),
            lambda pb: _plot_overlap_failure_modes(detail_df, pb),
        ]:
            with tempfile.TemporaryDirectory(prefix='bench_pdf_') as td:
                paths = fig_func(Path(td) / 'tmpfig')
                # The builder already saved and closed; reopen PNG on a page.
                img_path = next((p for p in paths if p.endswith('.png')), None)
                page = plt.figure(figsize=(11.69, 8.27))
                ax = page.add_subplot(111)
                ax.axis('off')
                if img_path and os.path.exists(img_path):
                    img = plt.imread(img_path)
                    ax.imshow(img)
                pdf.savefig(page, bbox_inches='tight')
                plt.close(page)

        # Appendix table page
        appendix = plt.figure(figsize=(11.69, 8.27))
        ax = appendix.add_subplot(111)
        _plot_table(ax, metrics_df if metrics_df is not None else pd.DataFrame(), 'Appendix. Core metrics', nrows=25)
        pdf.savefig(appendix, bbox_inches='tight'); plt.close(appendix)
    return output_pdf_path


def build_benchmark_report_bundle_bytes(suite: Dict[str, pd.DataFrame], metadata: Optional[Dict[str, Any]] = None, prefix: str = 'benchmark_papergrade') -> bytes:
    with tempfile.TemporaryDirectory(prefix='bench_report_export_') as tmpdir:
        manifest = export_benchmark_paper_figures(tmpdir, suite=suite, metadata=metadata, prefix=prefix)
        pdf_path = os.path.join(tmpdir, f'{prefix}_benchmark_report.pdf')
        export_benchmark_report_pdf(pdf_path, suite=suite, metadata=metadata, title=f'{prefix} benchmark report')
        manifest['report_pdf'] = os.path.relpath(pdf_path, tmpdir)
        with open(os.path.join(tmpdir, 'manifest.json'), 'w', encoding='utf-8') as f:
            json.dump(_json_safe(manifest), f, ensure_ascii=False, indent=2)
        bio = io.BytesIO()
        with zipfile.ZipFile(bio, 'w', compression=zipfile.ZIP_DEFLATED) as zf:
            for root, _, files in os.walk(tmpdir):
                for fn in files:
                    full = os.path.join(root, fn)
                    arc = os.path.relpath(full, tmpdir)
                    zf.write(full, arcname=arc)
        return bio.getvalue()
