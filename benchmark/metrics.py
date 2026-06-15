from __future__ import annotations

from typing import List
import numpy as np
import pandas as pd


def compute_multiclass_auc(y_true: List[str], prob_df: pd.DataFrame) -> float:
    try:
        from sklearn.preprocessing import label_binarize
        from sklearn.metrics import roc_auc_score
    except Exception:
        return float('nan')
    classes = list(prob_df.columns)
    if len(set(y_true)) < 2 or len(classes) < 2:
        return float('nan')
    y_bin = label_binarize(y_true, classes=classes)
    valid_cols = []
    for j, c in enumerate(classes):
        col = y_bin[:, j]
        if col.sum() > 0 and col.sum() < len(col):
            valid_cols.append(j)
    if not valid_cols:
        return float('nan')
    y_bin = y_bin[:, valid_cols]
    y_score = prob_df.iloc[:, valid_cols].to_numpy()
    return float(roc_auc_score(y_bin, y_score, average='macro'))


def summarize_publish_metrics(detail_df: pd.DataFrame, k: int = 3) -> pd.DataFrame:
    if detail_df.empty:
        return pd.DataFrame()
    out = {
        'metric': ['n', 'Top1 acc', f'Top{k} acc', 'Macro AUC', 'Acc_stab_high', 'Acc_stab_mid', 'Acc_stab_low'],
        'value': [
            int(len(detail_df)),
            float(detail_df['hit1'].mean()) if 'hit1' in detail_df else float('nan'),
            float(detail_df[f'hit{k}'].mean()) if f'hit{k}' in detail_df else float('nan'),
            float(detail_df.attrs.get('summary', {}).get('auc', float('nan'))),
            float(detail_df[detail_df['stability_bin']=='high']['hit1'].mean()) if 'stability_bin' in detail_df and (detail_df['stability_bin']=='high').any() else float('nan'),
            float(detail_df[detail_df['stability_bin']=='mid']['hit1'].mean()) if 'stability_bin' in detail_df and (detail_df['stability_bin']=='mid').any() else float('nan'),
            float(detail_df[detail_df['stability_bin']=='low']['hit1'].mean()) if 'stability_bin' in detail_df and (detail_df['stability_bin']=='low').any() else float('nan'),
        ]
    }
    return pd.DataFrame(out)
