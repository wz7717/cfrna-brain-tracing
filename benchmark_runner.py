from __future__ import annotations

import sqlite3
from typing import Callable, Dict, List, Optional, Tuple, Any

import numpy as np
import pandas as pd

from benchmark.label_utils import default_label_extractor, is_mixed_label, safe_soft_label as _safe_soft_label
from benchmark.metrics import compute_multiclass_auc as _compute_multiclass_auc
from core.methods import canonical_method
from source_tracing_v2 import SourceTracingEngineV2


def _json_safe_value(x):
    try:
        if pd.isna(x):
            return None
    except Exception:
        pass
    if isinstance(x, (np.floating, float)):
        return None if not np.isfinite(float(x)) else float(x)
    if isinstance(x, (np.integer, int)):
        return int(x)
    return x


def _attach_json_safe_summary(df: pd.DataFrame, summary: Dict[str, Any]) -> None:
    safe = {}
    for k, v in summary.items():
        if isinstance(v, dict):
            safe[k] = {kk: _json_safe_value(vv) for kk, vv in v.items()}
        else:
            safe[k] = _json_safe_value(v)
    df.attrs['summary'] = safe


def _load_labeled_samples(conn, label_extractor, limit=None, allow_mixed=False):
    samples = pd.read_sql_query("SELECT * FROM cfrna_samples", conn)
    labeled = []
    for _, r in samples.iterrows():
        label = _safe_soft_label(label_extractor(r.to_dict()))
        if not label:
            continue
        if (not allow_mixed) and is_mixed_label(label):
            continue
        labeled.append((r['sample_id'], label))
        if limit is not None and len(labeled) >= int(limit):
            break
    return labeled


def _normalize_method(method: str) -> str:
    return canonical_method(method)


def _safe_div(a: float, b: float) -> float:
    return float(a) / float(b) if float(b) != 0 else float('nan')


def _compute_binary_roc(y_true_bin: np.ndarray, y_score: np.ndarray):
    try:
        from sklearn.metrics import roc_curve, auc
    except Exception:
        return None, None, float('nan')
    if len(np.unique(y_true_bin)) < 2:
        return None, None, float('nan')
    fpr, tpr, _ = roc_curve(y_true_bin, y_score)
    return fpr, tpr, float(auc(fpr, tpr))


def _prepare_probability_frame(rows: List[Dict[str, Any]], class_list: Optional[List[str]]) -> pd.DataFrame:
    if not rows or not class_list:
        return pd.DataFrame()
    df = pd.DataFrame(rows)
    for c in class_list:
        if c not in df.columns:
            df[c] = 0.0
    return df[['sample_id', 'label'] + class_list]


def _compute_confusion_tables(valid_df: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    if valid_df.empty or 'label' not in valid_df or 'pred_top1' not in valid_df:
        return pd.DataFrame(), pd.DataFrame(), pd.DataFrame()
    labels = sorted(set(valid_df['label'].dropna().astype(str)) | set(valid_df['pred_top1'].dropna().astype(str)))
    if not labels:
        return pd.DataFrame(), pd.DataFrame(), pd.DataFrame()
    raw = pd.crosstab(
        valid_df['label'].astype(str),
        valid_df['pred_top1'].astype(str),
        rownames=['truth'],
        colnames=['pred'],
        dropna=False,
    ).reindex(index=labels, columns=labels, fill_value=0)
    norm = raw.div(raw.sum(axis=1).replace(0, np.nan), axis=0).fillna(0.0)
    long_rows = []
    for truth in raw.index:
        total_truth = int(raw.loc[truth].sum())
        for pred in raw.columns:
            count = int(raw.loc[truth, pred])
            long_rows.append({
                'truth_region': truth,
                'pred_region': pred,
                'count': count,
                'row_fraction': float(norm.loc[truth, pred]),
                'n_truth_samples': total_truth,
            })
    long_df = pd.DataFrame(long_rows)
    return raw, norm, long_df


def _compute_roc_tables(prob_df: pd.DataFrame, class_list: Optional[List[str]]) -> Tuple[pd.DataFrame, pd.DataFrame]:
    if prob_df.empty or not class_list or 'label' not in prob_df.columns:
        return pd.DataFrame(), pd.DataFrame()
    curve_rows = []
    summary_rows = []
    for region in class_list:
        y_true_bin = (prob_df['label'].astype(str).to_numpy() == str(region)).astype(int)
        y_score = pd.to_numeric(prob_df[region], errors='coerce').fillna(0.0).to_numpy()
        fpr, tpr, auc_val = _compute_binary_roc(y_true_bin, y_score)
        n_pos = int(y_true_bin.sum())
        n_neg = int((1 - y_true_bin).sum())
        summary_rows.append({
            'region_id': region,
            'auc': _json_safe_value(auc_val),
            'n_pos': n_pos,
            'n_neg': n_neg,
        })
        if fpr is not None and tpr is not None:
            for fp, tp in zip(fpr, tpr):
                curve_rows.append({'region_id': region, 'fpr': float(fp), 'tpr': float(tp), 'auc': _json_safe_value(auc_val)})
    return pd.DataFrame(curve_rows), pd.DataFrame(summary_rows)


def _compute_separability_tables(prob_df: pd.DataFrame, detail_df: pd.DataFrame, class_list: Optional[List[str]]) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    if prob_df.empty or not class_list:
        return pd.DataFrame(), pd.DataFrame(), pd.DataFrame()
    X = prob_df[class_list].to_numpy(dtype=float)
    if X.shape[0] < 2 or X.shape[1] < 2:
        return pd.DataFrame(), pd.DataFrame(), pd.DataFrame()
    try:
        from sklearn.decomposition import PCA
        pca = PCA(n_components=2)
        coords = pca.fit_transform(X)
        evr = pca.explained_variance_ratio_
    except Exception:
        return pd.DataFrame(), pd.DataFrame(), pd.DataFrame()
    sep_df = prob_df[['sample_id', 'label']].copy()
    sep_df['pc1'] = coords[:, 0]
    sep_df['pc2'] = coords[:, 1]
    sep_df['pc1_var'] = float(evr[0]) if len(evr) > 0 else float('nan')
    sep_df['pc2_var'] = float(evr[1]) if len(evr) > 1 else float('nan')
    if not detail_df.empty:
        merge_cols = [c for c in ['sample_id', 'pred_top1', 'hit1', 'top1_confidence', 'top1_stability', 'decision_margin'] if c in detail_df.columns]
        sep_df = sep_df.merge(detail_df[merge_cols], on='sample_id', how='left')
    centroid_rows = []
    for label, g in sep_df.groupby('label'):
        centroid_rows.append({
            'region_id': label,
            'centroid_pc1': float(g['pc1'].mean()),
            'centroid_pc2': float(g['pc2'].mean()),
            'n_samples': int(len(g)),
        })
    centroid_df = pd.DataFrame(centroid_rows)
    dist_rows = []
    if not centroid_df.empty:
        for _, r1 in centroid_df.iterrows():
            for _, r2 in centroid_df.iterrows():
                d = float(np.sqrt((r1['centroid_pc1'] - r2['centroid_pc1']) ** 2 + (r1['centroid_pc2'] - r2['centroid_pc2']) ** 2))
                dist_rows.append({'region_a': r1['region_id'], 'region_b': r2['region_id'], 'centroid_distance': d})
    dist_df = pd.DataFrame(dist_rows)
    return sep_df, centroid_df, dist_df


def _compute_stability_tables(detail_df: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
    if detail_df.empty:
        return pd.DataFrame(), pd.DataFrame()
    df = detail_df.copy()
    if 'top1_stability' not in df.columns:
        df['top1_stability'] = np.nan
    if 'decision_margin' not in df.columns:
        df['decision_margin'] = np.nan
    if 'top1_confidence' not in df.columns:
        df['top1_confidence'] = np.nan
    region_summary = []
    for label, g in df.groupby('label'):
        valid = g[g['abstained'] == 0] if 'abstained' in g.columns else g
        region_summary.append({
            'region_id': label,
            'n_samples': int(len(g)),
            'n_valid': int(len(valid)),
            'top1_acc_valid': float(valid['hit1'].mean()) if len(valid) and 'hit1' in valid.columns else float('nan'),
            'mean_top1_confidence': float(valid['top1_confidence'].mean()) if len(valid) else float('nan'),
            'median_top1_confidence': float(valid['top1_confidence'].median()) if len(valid) else float('nan'),
            'mean_top1_stability': float(valid['top1_stability'].mean()) if len(valid) else float('nan'),
            'median_top1_stability': float(valid['top1_stability'].median()) if len(valid) else float('nan'),
            'mean_decision_margin': float(valid['decision_margin'].mean()) if len(valid) else float('nan'),
            'abstain_rate': float(g['abstained'].mean()) if 'abstained' in g.columns else 0.0,
        })
    bins = [0.0, 0.5, 0.8, 1.0]
    labels = ['low(<0.5)', 'mid(0.5-0.8)', 'high(>=0.8)']
    stab_detail = df.copy()
    stab_detail['stability_bin'] = pd.cut(stab_detail['top1_stability'], bins=bins, labels=labels, include_lowest=True, right=True)
    stab_bins = []
    for b in labels:
        sub = stab_detail[stab_detail['stability_bin'] == b]
        stab_bins.append({
            'stability_bin': b,
            'n_samples': int(len(sub)),
            'top1_acc': float(sub['hit1'].mean()) if len(sub) and 'hit1' in sub.columns else float('nan'),
            'mean_confidence': float(sub['top1_confidence'].mean()) if len(sub) else float('nan'),
            'mean_margin': float(sub['decision_margin'].mean()) if len(sub) else float('nan'),
        })
    return pd.DataFrame(region_summary), pd.DataFrame(stab_bins)


def run_topk_benchmark(db_path: str, method: str = 'ensemble', k: int = 3, label_extractor: Callable[[Dict], Optional[str]] = default_label_extractor, limit: Optional[int] = None, allow_mixed: bool = False) -> pd.DataFrame:
    conn = sqlite3.connect(db_path)
    try:
        labeled = _load_labeled_samples(conn, label_extractor, limit=limit, allow_mixed=allow_mixed)
    finally:
        conn.close()
    engine = SourceTracingEngineV2(db_path)
    rows = []
    method = _normalize_method(method)
    for sample_id, label in labeled:
        out = engine.trace(sample_id=sample_id, method=method, return_all=True, bootstrap_n=0, persist=False, abstain_on_low_overlap=True, min_overlap_genes=5, min_overlap_fraction=0.003)
        res = out.get('results', [])
        preds = [x['region_id'] for x in res]
        rows.append({'sample_id': sample_id, 'label': label, 'pred_top1': preds[0] if preds else None, 'hit1': int(bool(preds) and preds[0] == label), f'hit{k}': int(label in preds[:k]), 'abstained': int(len(preds) == 0), 'traceability': out.get('meta', {}).get('traceability', 'unknown'), 'overlap_genes': out.get('meta', {}).get('sample_overlap', {}).get('n_overlap_genes')})
    df = pd.DataFrame(rows)
    if not df.empty:
        valid = df[df['abstained'] == 0]
        _attach_json_safe_summary(df, {'n': int(len(df)), 'n_valid': int(len(valid)), 'abstain_rate': float(df['abstained'].mean()), 'top1_acc_valid': float(valid['hit1'].mean()) if len(valid) else float('nan'), f'top{k}_acc_valid': float(valid[f'hit{k}'].mean()) if len(valid) else float('nan')})
    return df


def auto_tune_ensemble_weights(db_path: str, atlas_id: int = 1, sigset_id: Optional[int] = None, use_value: str = 'log1p', l2_grid: Optional[List[float]] = None, alpha_grid: Optional[List[float]] = None, label_extractor: Callable[[Dict], Optional[str]] = default_label_extractor, limit: Optional[int] = None, optimize_metric: str = 'top1_acc', allow_mixed: bool = False, n_splits: int = 3) -> Tuple[Dict[str, Any], pd.DataFrame]:
    """在标注样本上用网格搜索 + K-fold 交叉验证选出最优 (l2, ensemble_alpha) 组合。

    Args:
        n_splits: K-fold 折数。样本数不足时自动退化为 leave-one-out（n_splits=n_samples）。
                  设为 1 时退化为原有的不分折模式（不推荐，仅用于调试）。
    """
    if l2_grid is None:
        l2_grid = [0.0, 1e-4, 1e-3, 1e-2]
    if alpha_grid is None:
        alpha_grid = [i / 10 for i in range(0, 11)]

    conn = sqlite3.connect(db_path)
    try:
        labeled = _load_labeled_samples(conn, label_extractor, limit=limit, allow_mixed=allow_mixed)
    finally:
        conn.close()
    if not labeled:
        return {'error': 'No labeled samples'}, pd.DataFrame()

    n = len(labeled)
    # 样本不足时自动调整折数
    effective_splits = min(int(n_splits), n) if n >= 2 else 1
    if effective_splits < 2:
        # fallback: 不做交叉验证，沿用老逻辑（会有过拟合风险）
        import warnings as _warnings
        _warnings.warn(
            f"auto_tune: only {n} labeled samples, cannot do {n_splits}-fold CV. "
            "Falling back to in-sample evaluation (results may be overfit).",
            UserWarning,
            stacklevel=2,
        )
        fold_indices = [(list(range(n)), list(range(n)))]  # train=test=all
    else:
        # 构造 stratified K-fold（标签列表）
        labels_list = [lbl for _, lbl in labeled]
        try:
            from sklearn.model_selection import StratifiedKFold
            skf = StratifiedKFold(n_splits=effective_splits, shuffle=True, random_state=42)
            fold_indices = list(skf.split(range(n), labels_list))
        except Exception:
            # sklearn 不可用时退化为均匀 K-fold
            fold_size = n // effective_splits
            fold_indices = []
            for k in range(effective_splits):
                val_idx = list(range(k * fold_size, min((k + 1) * fold_size, n)))
                train_idx = [i for i in range(n) if i not in set(val_idx)]
                fold_indices.append((train_idx, val_idx))

    engine = SourceTracingEngineV2(db_path)
    rows = []
    best: Optional[Dict[str, Any]] = None
    best_val = -1e18

    for l2 in l2_grid:
        for alpha in alpha_grid:
            fold_accs = []
            fold_aucs = []
            for _train_idx, val_idx in fold_indices:
                val_samples = [labeled[i] for i in val_idx]
                y_true_fold, y_pred_fold = [], []
                prob_rows_fold, class_list_fold = [], None
                for sample_id, label in val_samples:
                    out = engine.trace(
                        sample_id=sample_id,
                        method='ensemble',
                        sigset_id=sigset_id,
                        atlas_id=atlas_id,
                        use_value=use_value,
                        l2=float(l2),
                        ensemble_alpha=float(alpha),
                        bootstrap_n=0,
                        persist=False,
                        return_all=True,
                        abstain_on_low_overlap=True,
                        min_overlap_genes=5,
                        min_overlap_fraction=0.003,
                    )
                    res = out.get('results', [])
                    if not res:
                        continue
                    if class_list_fold is None:
                        class_list_fold = [x['region_id'] for x in res]
                    prob = {x['region_id']: float(x.get('confidence', 0.0)) for x in res}
                    prob_rows_fold.append([prob.get(c, 0.0) for c in class_list_fold])
                    y_true_fold.append(label)
                    y_pred_fold.append(res[0]['region_id'])

                if len(y_true_fold) < 2:
                    continue
                fold_acc = float(np.mean([1.0 if p == t else 0.0 for p, t in zip(y_pred_fold, y_true_fold)]))
                fold_accs.append(fold_acc)
                fold_auc = float('nan')
                if class_list_fold is not None:
                    fold_auc = _compute_multiclass_auc(y_true_fold, pd.DataFrame(prob_rows_fold, columns=class_list_fold))
                fold_aucs.append(fold_auc)

            if not fold_accs:
                continue

            mean_acc = float(np.mean(fold_accs))
            valid_aucs = [a for a in fold_aucs if not np.isnan(a)]
            mean_auc = float(np.mean(valid_aucs)) if valid_aucs else float('nan')
            val = mean_acc if optimize_metric == 'top1_acc' else (mean_auc if not np.isnan(mean_auc) else -1e9)
            rows.append({
                'l2': float(l2),
                'ensemble_alpha': float(alpha),
                'cv_top1_acc': mean_acc,
                'cv_auc': mean_auc,
                'n_folds': len(fold_accs),
            })
            if val > best_val:
                best_val = val
                best = {
                    'atlas_id': atlas_id,
                    'sigset_id': sigset_id,
                    'use_value': use_value,
                    'l2': float(l2),
                    'ensemble_alpha': float(alpha),
                    'optimize_metric': optimize_metric,
                    'cv_top1_acc': mean_acc,
                    'cv_auc': mean_auc,
                    'n_splits_used': effective_splits,
                }

    grid_df = pd.DataFrame(rows)
    if not grid_df.empty:
        sort_col = 'cv_top1_acc' if optimize_metric == 'top1_acc' else 'cv_auc'
        grid_df = grid_df.sort_values(by=sort_col, ascending=False)
    if best is None:
        best = {'error': 'Tuning failed'}
    return best, grid_df


def run_publish_grade_benchmark(db_path: str, method: str = 'ensemble', k: int = 3, atlas_id: int = 1, sigset_id: Optional[int] = None, use_value: str = 'log1p', l2: float = 1e-4, ensemble_alpha: float = 0.5, bootstrap_n: int = 50, bootstrap_gene_frac: float = 0.7, label_extractor: Callable[[Dict], Optional[str]] = default_label_extractor, limit: Optional[int] = None, allow_mixed: bool = False) -> Tuple[pd.DataFrame, pd.DataFrame]:
    suite = run_paper_grade_benchmark_suite(
        db_path=db_path,
        method=method,
        k=k,
        atlas_id=atlas_id,
        sigset_id=sigset_id,
        use_value=use_value,
        l2=l2,
        ensemble_alpha=ensemble_alpha,
        bootstrap_n=bootstrap_n,
        bootstrap_gene_frac=bootstrap_gene_frac,
        label_extractor=label_extractor,
        limit=limit,
        allow_mixed=allow_mixed,
    )
    return suite['detail_df'], suite['metrics_df']


def run_paper_grade_benchmark_suite(db_path: str, method: str = 'ensemble', k: int = 3, atlas_id: int = 1, sigset_id: Optional[int] = None, use_value: str = 'log1p', l2: float = 1e-4, ensemble_alpha: float = 0.5, bootstrap_n: int = 50, bootstrap_gene_frac: float = 0.7, label_extractor: Callable[[Dict], Optional[str]] = default_label_extractor, limit: Optional[int] = None, allow_mixed: bool = False) -> Dict[str, pd.DataFrame]:
    conn = sqlite3.connect(db_path)
    try:
        labeled = _load_labeled_samples(conn, label_extractor, limit=limit, allow_mixed=allow_mixed)
    finally:
        conn.close()
    empty = {
        'detail_df': pd.DataFrame(),
        'metrics_df': pd.DataFrame(),
        'probability_df': pd.DataFrame(),
        'confusion_raw_df': pd.DataFrame(),
        'confusion_norm_df': pd.DataFrame(),
        'confusion_long_df': pd.DataFrame(),
        'roc_curve_df': pd.DataFrame(),
        'roc_summary_df': pd.DataFrame(),
        'separability_df': pd.DataFrame(),
        'centroid_df': pd.DataFrame(),
        'centroid_distance_df': pd.DataFrame(),
        'stability_region_df': pd.DataFrame(),
        'stability_bin_df': pd.DataFrame(),
    }
    if not labeled:
        return empty

    method = _normalize_method(method)
    engine = SourceTracingEngineV2(db_path)
    rows = []
    prob_rows = []
    y_true = []
    class_list = None
    for sample_id, label in labeled:
        out = engine.trace(
            sample_id=sample_id,
            method=method,
            sigset_id=sigset_id,
            atlas_id=atlas_id,
            use_value=use_value,
            l2=float(l2),
            ensemble_alpha=float(ensemble_alpha),
            bootstrap_n=int(bootstrap_n),
            bootstrap_gene_frac=float(bootstrap_gene_frac),
            persist=False,
            return_all=True,
            abstain_on_low_overlap=True,
            min_overlap_genes=5,
            min_overlap_fraction=0.003,
        )
        res = out.get('results', [])
        meta = out.get('meta', {})
        if not res:
            rows.append({
                'sample_id': sample_id,
                'label': label,
                'pred_top1': None,
                'pred_top2': None,
                'hit1': 0,
                f'hit{k}': 0,
                'top1_confidence': np.nan,
                'top2_confidence': np.nan,
                'decision_margin': np.nan,
                'top1_stability': np.nan,
                'abstained': 1,
                'traceability': meta.get('traceability'),
                'overlap_genes': meta.get('sample_overlap', {}).get('n_overlap_genes'),
            })
            continue
        preds = [x['region_id'] for x in res]
        pred1 = preds[0]
        pred2 = preds[1] if len(preds) > 1 else None
        top1_conf = float(res[0].get('confidence', 0.0))
        top2_conf = float(res[1].get('confidence', 0.0)) if len(res) > 1 else 0.0
        hit1 = int(pred1 == label)
        hitk = int(label in preds[:k])
        top1_stab = float(res[0].get('stability')) if 'stability' in res[0] and res[0].get('stability') is not None else np.nan
        decision_margin = float(meta.get('decision_margin')) if meta.get('decision_margin') is not None else float(top1_conf - top2_conf)
        rows.append({
            'sample_id': sample_id,
            'label': label,
            'pred_top1': pred1,
            'pred_top2': pred2,
            'hit1': hit1,
            f'hit{k}': hitk,
            'top1_confidence': top1_conf,
            'top2_confidence': top2_conf,
            'decision_margin': decision_margin,
            'top1_stability': top1_stab,
            'abstained': 0,
            'traceability': meta.get('traceability'),
            'overlap_genes': meta.get('sample_overlap', {}).get('n_overlap_genes'),
        })
        if class_list is None:
            class_list = [x['region_id'] for x in res]
        prob = {x['region_id']: float(x.get('confidence', 0.0)) for x in res}
        prob_row = {'sample_id': sample_id, 'label': label}
        for c in class_list:
            prob_row[c] = prob.get(c, 0.0)
        prob_rows.append(prob_row)
        y_true.append(label)

    detail_df = pd.DataFrame(rows)
    if detail_df.empty:
        return empty
    valid = detail_df[detail_df['abstained'] == 0].copy()
    probability_df = _prepare_probability_frame(prob_rows, class_list)
    top1_acc = float(valid['hit1'].mean()) if len(valid) else float('nan')
    topk_acc = float(valid[f'hit{k}'].mean()) if len(valid) else float('nan')
    auc = float('nan')
    if not probability_df.empty and class_list is not None and len(probability_df) >= 2:
        auc = _compute_multiclass_auc(list(probability_df['label'].astype(str)), probability_df[class_list])
    strat = {'high(>=0.8)': float('nan'), 'mid(0.5-0.8)': float('nan'), 'low(<0.5)': float('nan')}
    if valid['top1_stability'].notna().any():
        tmp = valid.dropna(subset=['top1_stability']).copy()
        if not tmp.empty:
            tmp['stab_bin'] = tmp['top1_stability'].apply(lambda s: 'high(>=0.8)' if s >= 0.8 else 'mid(0.5-0.8)' if s >= 0.5 else 'low(<0.5)')
            for bname, g in tmp.groupby('stab_bin'):
                strat[bname] = float(g['hit1'].mean())

    n_classes = int(len(class_list)) if class_list else int(valid['label'].nunique())
    balanced_acc = float(valid.groupby('label')['hit1'].mean().mean()) if len(valid) else float('nan')
    metrics_df = pd.DataFrame([
        {'metric': 'N_samples_total', 'value': int(len(detail_df))},
        {'metric': 'N_samples_valid', 'value': int(len(valid))},
        {'metric': 'N_classes', 'value': n_classes},
        {'metric': 'Abstain_rate', 'value': float(detail_df['abstained'].mean())},
        {'metric': 'Top1_acc_valid', 'value': top1_acc},
        {'metric': f'Top{k}_acc_valid', 'value': topk_acc},
        {'metric': 'Balanced_acc_valid', 'value': balanced_acc},
        {'metric': 'MacroAUC_ovr_valid', 'value': auc},
        {'metric': 'Mean_top1_confidence_valid', 'value': float(valid['top1_confidence'].mean()) if len(valid) else float('nan')},
        {'metric': 'Mean_decision_margin_valid', 'value': float(valid['decision_margin'].mean()) if len(valid) else float('nan')},
        {'metric': 'Mean_top1_stability_valid', 'value': float(valid['top1_stability'].mean()) if len(valid) else float('nan')},
        {'metric': 'Acc_stab_high', 'value': strat['high(>=0.8)']},
        {'metric': 'Acc_stab_mid', 'value': strat['mid(0.5-0.8)']},
        {'metric': 'Acc_stab_low', 'value': strat['low(<0.5)']},
    ])

    confusion_raw_df, confusion_norm_df, confusion_long_df = _compute_confusion_tables(valid)
    roc_curve_df, roc_summary_df = _compute_roc_tables(probability_df, class_list)
    separability_df, centroid_df, centroid_distance_df = _compute_separability_tables(probability_df, detail_df, class_list)
    stability_region_df, stability_bin_df = _compute_stability_tables(detail_df)

    safe_strat = {k: _json_safe_value(v) for k, v in strat.items()}
    _attach_json_safe_summary(detail_df, {
        'n_total': int(len(detail_df)),
        'n_valid': int(len(valid)),
        'n_classes': n_classes,
        'abstain_rate': float(detail_df['abstained'].mean()),
        'top1_acc': top1_acc,
        f'top{k}_acc': topk_acc,
        'balanced_acc': balanced_acc,
        'auc': auc,
        'mean_top1_confidence': float(valid['top1_confidence'].mean()) if len(valid) else float('nan'),
        'mean_decision_margin': float(valid['decision_margin'].mean()) if len(valid) else float('nan'),
        'mean_top1_stability': float(valid['top1_stability'].mean()) if len(valid) else float('nan'),
        'strat': safe_strat,
    })

    return {
        'detail_df': detail_df,
        'metrics_df': metrics_df,
        'probability_df': probability_df,
        'confusion_raw_df': confusion_raw_df,
        'confusion_norm_df': confusion_norm_df,
        'confusion_long_df': confusion_long_df,
        'roc_curve_df': roc_curve_df,
        'roc_summary_df': roc_summary_df,
        'separability_df': separability_df,
        'centroid_df': centroid_df,
        'centroid_distance_df': centroid_distance_df,
        'stability_region_df': stability_region_df,
        'stability_bin_df': stability_bin_df,
    }
