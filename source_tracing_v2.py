from __future__ import annotations

import json
import logging
import sqlite3
import uuid
from typing import Any, Dict, List, Optional, Tuple

logger = logging.getLogger(__name__)

import numpy as np
import pandas as pd
from scipy.stats import rankdata

from core.bootstrap import bootstrap_nnls
from core.gene_utils import guess_gene_id_type
from core.methods import canonical_method
from core.models import apply_value_transform, softmax_confidence, trace_corr, trace_nnls_simplex, zscore_safe
from core.params import TraceParams
from core.reference_loader import (
    _col_exists as rl_col_exists,
    _table_exists as rl_table_exists,
    get_latest_sigset_id as rl_get_latest_sigset_id,
    load_marker_signature_genes as rl_load_marker_signature_genes,
    load_signature_genes as rl_load_signature_genes,
)


class SourceTracingEngineV2:
    """发布级 cfRNA 脑区溯源引擎。

    这一版将原先主要依赖 TPM/logTPM 的集成算法升级为多信号联合判定：
    1. 加权表达相关性
    2. Simplex-NNLS 混合比例
    3. 全局表达秩相关
    4. 区域 marker 富集强度
    5. marker 检出覆盖度
    6. 丰度/reads 支持度

    前三项强调“整体表达形状”，后三项强调“区域特异 marker 是否真的被看到”。
    """

    def __init__(self, db_path: str):
        self.db_path = db_path
        self._cache_ref: Dict[Tuple[Any, ...], Tuple[np.ndarray, np.ndarray, np.ndarray, List[str], Dict[str, Any]]] = {}

    def _connect(self) -> sqlite3.Connection:
        return sqlite3.connect(self.db_path)

    def _table_exists(self, conn: sqlite3.Connection, name: str) -> bool:
        return rl_table_exists(conn, name)

    def _col_exists(self, conn: sqlite3.Connection, table: str, col: str) -> bool:
        return rl_col_exists(conn, table, col)

    def _parse_params(self, params: Dict[str, Any]) -> Tuple[TraceParams, Dict[str, Any]]:
        tp = TraceParams()
        for k, v in params.items():
            if hasattr(tp, k):
                setattr(tp, k, v)
        tp.bootstrap_gene_frac = float(np.clip(tp.bootstrap_gene_frac, 0.1, 1.0))
        tp.bootstrap_n = int(max(0, tp.bootstrap_n))
        tp.topk = int(max(1, tp.topk))
        tp.l2 = float(max(0.0, tp.l2))
        tp.atlas_id = int(tp.atlas_id)
        tp.ensemble_alpha = float(np.clip(tp.ensemble_alpha, 0.0, 1.0))
        tp.return_all = bool(tp.return_all)
        extra = {
            'enable_weighting': bool(params.get('enable_weighting', True)),
            'specificity_weight': float(max(0.0, params.get('specificity_weight', 1.25))),
            'marker_weight': float(max(0.0, params.get('marker_weight', 0.75))),
            'marker_bonus': float(max(0.0, params.get('marker_bonus', 0.75))),
            'min_weight': float(max(0.01, params.get('min_weight', 0.25))),
            'max_weight': float(max(1.0, params.get('max_weight', 6.0))),
            'normalize_region_weight': bool(params.get('normalize_region_weight', True)),
            'use_read_count_hint': bool(params.get('use_read_count_hint', True)),
            'read_count_weight': float(max(0.0, params.get('read_count_weight', 0.15))),
            'fallback_marker_topk': int(max(20, params.get('fallback_marker_topk', 200))),
            'min_overlap_genes': int(max(1, params.get('min_overlap_genes', 20))),
            'min_overlap_fraction': float(max(0.0, params.get('min_overlap_fraction', 0.01))),
            'abstain_on_low_overlap': bool(params.get('abstain_on_low_overlap', False)),
            'marker_panel_topk': int(max(10, params.get('marker_panel_topk', 40))),
            'corr_component_weight': float(max(0.0, params.get('corr_component_weight', 0.28))),
            'nnls_component_weight': float(max(0.0, params.get('nnls_component_weight', 0.24))),
            'rank_component_weight': float(max(0.0, params.get('rank_component_weight', 0.16))),
            'marker_component_weight': float(max(0.0, params.get('marker_component_weight', 0.16))),
            'detect_component_weight': float(max(0.0, params.get('detect_component_weight', 0.10))),
            'support_component_weight': float(max(0.0, params.get('support_component_weight', 0.06))),
            'abstain_on_ambiguous': bool(params.get('abstain_on_ambiguous', False)),
            'min_top1_confidence': float(np.clip(params.get('min_top1_confidence', 0.0), 0.0, 1.0)),
            'min_decision_margin': float(max(0.0, params.get('min_decision_margin', 0.0))),
        }
        return tp, extra

    def _load_signature_genes(self, sigset_id: int) -> Optional[List[str]]:
        return rl_load_signature_genes(self.db_path, sigset_id)

    def _get_latest_sigset_id(self, atlas_id: int = 1) -> Optional[int]:
        return rl_get_latest_sigset_id(self.db_path, atlas_id)

    def _get_atlas_meta(self, atlas_id: int) -> Dict[str, Any]:
        conn = self._connect()
        try:
            if not self._table_exists(conn, 'atlas_versions'):
                return {}
            row = conn.execute(
                """
                SELECT atlas_id, atlas_name, species, level, build_version,
                       gene_id_type, normalization, notes
                FROM atlas_versions
                WHERE atlas_id = ?
                """,
                (int(atlas_id),),
            ).fetchone()
            if row is None:
                return {}
            keys = ['atlas_id', 'atlas_name', 'species', 'level', 'build_version', 'gene_id_type', 'normalization', 'notes']
            return dict(zip(keys, row))
        finally:
            conn.close()

    @staticmethod
    def _is_vsd_reference(atlas_meta: Dict[str, Any]) -> bool:
        text = " ".join(
            str(atlas_meta.get(k, "") or "").lower()
            for k in ["atlas_name", "build_version", "normalization", "notes"]
        )
        return "vsd" in text or "batch_removed" in text or "batch removed" in text

    def _load_marker_signature_genes(self, topk_per_region: int = 200) -> Optional[List[str]]:
        return rl_load_marker_signature_genes(self.db_path, topk_per_region=topk_per_region)

    @staticmethod
    def _safe_numeric(s: pd.Series, fill: float = 0.0) -> pd.Series:
        return pd.to_numeric(s, errors='coerce').fillna(fill)

    @staticmethod
    def _clip_array(x: np.ndarray, lo: float, hi: float) -> np.ndarray:
        return np.clip(np.asarray(x, dtype=float), lo, hi)

    @staticmethod
    def _normalize_col_mean_one(M: np.ndarray) -> np.ndarray:
        M = np.asarray(M, dtype=float)
        col_mean = M.mean(axis=0, keepdims=True)
        col_mean[col_mean <= 0] = 1.0
        return M / col_mean

    @staticmethod
    def _zscore_safe(v: np.ndarray) -> np.ndarray:
        return zscore_safe(v)

    @staticmethod
    def _softmax_confidence(scores: np.ndarray) -> np.ndarray:
        return softmax_confidence(scores)

    @staticmethod
    def _apply_value_transform(x: np.ndarray, use_value: str) -> np.ndarray:
        return apply_value_transform(x, use_value)

    @staticmethod
    def _trace_corr(A: np.ndarray, b: np.ndarray) -> np.ndarray:
        return trace_corr(A, b)

    @staticmethod
    def _trace_nnls_simplex(A: np.ndarray, b: np.ndarray, l2: float) -> Tuple[np.ndarray, float]:
        return trace_nnls_simplex(A, b, l2)

    @staticmethod
    def _bootstrap_nnls(A: np.ndarray, b: np.ndarray, regions: List[str], n: int, gene_frac: float, l2: float, seed: int):
        return bootstrap_nnls(A, b, regions, n, gene_frac, l2, seed)

    def _load_reference_matrix(self, atlas_id: int, sigset_id: Optional[int], use_value: str):
        genes, A_base, _W, regions, _meta = self._load_reference_bundle(
            atlas_id, sigset_id, use_value, False, 0.0, 0.0, 0.0, 0.25, 6.0, False, 200
        )
        return genes, A_base, regions

    def _load_sample_vector(self, sample_id: str, genes: np.ndarray, use_value: str) -> np.ndarray:
        b, _w, _m = self._load_sample_bundle(sample_id, genes, use_value, np.ones((len(genes), 1)), False, 0.0, 1.0, 1.0)
        return b

    def _load_reference_bundle(
        self,
        atlas_id: int,
        sigset_id: Optional[int],
        use_value: str,
        enable_weighting: bool,
        specificity_weight: float,
        marker_weight: float,
        marker_bonus: float,
        min_weight: float,
        max_weight: float,
        normalize_region_weight: bool,
        fallback_marker_topk: int = 200,
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray, List[str], Dict[str, Any]]:
        cache_key = (
            int(atlas_id),
            int(sigset_id) if sigset_id is not None else None,
            str(use_value),
            bool(enable_weighting),
            float(specificity_weight),
            float(marker_weight),
            float(marker_bonus),
            float(min_weight),
            float(max_weight),
            bool(normalize_region_weight),
            int(fallback_marker_topk),
        )
        if cache_key in self._cache_ref:
            return self._cache_ref[cache_key]

        genes_filter = self._load_signature_genes(sigset_id) if sigset_id is not None else None
        if genes_filter is None:
            genes_filter = self._load_marker_signature_genes(topk_per_region=fallback_marker_topk)

        conn = self._connect()
        try:
            select_sql = (
                'SELECT region_id, gene_symbol, avg_tpm, atlas_id FROM reference_expression'
                if self._col_exists(conn, 'reference_expression', 'atlas_id')
                else 'SELECT region_id, gene_symbol, avg_tpm FROM reference_expression'
            )
            ref = pd.read_sql_query(select_sql, conn)
            if ref.empty:
                empty = (
                    np.array([], dtype=object),
                    np.zeros((0, 0), dtype=float),
                    np.zeros((0, 0), dtype=float),
                    [],
                    {'weighting_active': False, 'n_weighted_pairs': 0},
                )
                self._cache_ref[cache_key] = empty
                return empty

            ref = ref.dropna(subset=['region_id', 'gene_symbol', 'avg_tpm']).copy()
            ref['gene_symbol'] = ref['gene_symbol'].astype(str)
            ref['avg_tpm'] = self._safe_numeric(ref['avg_tpm'], fill=0.0)
            if 'atlas_id' in ref.columns:
                ref = ref[(ref['atlas_id'].isna()) | (ref['atlas_id'].astype(float).astype(int) == int(atlas_id))]
            if genes_filter is not None:
                ref = ref[ref['gene_symbol'].isin(set(map(str, genes_filter)))]

            mat = ref.pivot_table(index='gene_symbol', columns='region_id', values='avg_tpm', aggfunc='mean').fillna(0.0)
            mat = mat.loc[mat.sum(axis=1) > 0]
            genes = mat.index.values.astype(object)
            regions = list(mat.columns)
            # 同时保留原始 TPM 矩阵，供 trace() 中秩相关等信号直接使用，避免重复 IO
            A_tpm_raw = mat.values.astype(float)
            A_base = self._apply_value_transform(A_tpm_raw, use_value)
            W_ref = np.ones_like(A_base, dtype=float)
            meta: Dict[str, Any] = {
                'weighting_active': False,
                'n_weighted_pairs': 0,
                'reference_gene_id_type': guess_gene_id_type(list(map(str, genes[: min(len(genes), 200)]))),
                '_A_tpm_cache': A_tpm_raw,   # 内部缓存，供 trace() 使用
            }

            if enable_weighting and self._table_exists(conn, 'region_gene_signature'):
                sig_sql = (
                    'SELECT region_id, gene_symbol, marker_score, specificity_score, is_marker, atlas_id FROM region_gene_signature'
                    if self._col_exists(conn, 'region_gene_signature', 'atlas_id')
                    else 'SELECT region_id, gene_symbol, marker_score, specificity_score, is_marker FROM region_gene_signature'
                )
                sig = pd.read_sql_query(sig_sql, conn)
                if not sig.empty:
                    sig = sig.dropna(subset=['region_id', 'gene_symbol']).copy()
                    sig['gene_symbol'] = sig['gene_symbol'].astype(str)
                    if 'atlas_id' in sig.columns:
                        sig = sig[(sig['atlas_id'].isna()) | (sig['atlas_id'].astype(float).astype(int) == int(atlas_id))]
                    if genes_filter is not None:
                        sig = sig[sig['gene_symbol'].isin(set(map(str, genes_filter)))]
                    sig['marker_score'] = self._safe_numeric(sig.get('marker_score', 0.0), fill=0.0)
                    sig['specificity_score'] = self._safe_numeric(sig.get('specificity_score', 0.0), fill=0.0)
                    sig['is_marker'] = self._safe_numeric(sig.get('is_marker', 0.0), fill=0.0)
                    if not sig.empty:
                        ms = (
                            sig.pivot_table(index='gene_symbol', columns='region_id', values='marker_score', aggfunc='mean')
                            .reindex(index=mat.index, columns=mat.columns)
                            .fillna(0.0)
                            .values.astype(float)
                        )
                        sp = (
                            sig.pivot_table(index='gene_symbol', columns='region_id', values='specificity_score', aggfunc='mean')
                            .reindex(index=mat.index, columns=mat.columns)
                            .fillna(0.0)
                            .values.astype(float)
                        )
                        mk = (
                            sig.pivot_table(index='gene_symbol', columns='region_id', values='is_marker', aggfunc='max')
                            .reindex(index=mat.index, columns=mat.columns)
                            .fillna(0.0)
                            .values.astype(float)
                            > 0
                        ).astype(float)
                        W_ref = 1.0 + float(specificity_weight) * np.clip(sp, 0.0, None) + float(marker_weight) * np.clip(ms, 0.0, None) + float(marker_bonus) * mk
                        W_ref = self._clip_array(W_ref, min_weight, max_weight)
                        if normalize_region_weight:
                            W_ref = self._normalize_col_mean_one(W_ref)
                        meta = {
                            'weighting_active': True,
                            'n_weighted_pairs': int(np.sum((sp > 0) | (ms > 0) | (mk > 0))),
                            'used_marker_score': bool(np.any(ms > 0)),
                            'used_specificity_score': bool(np.any(sp > 0)),
                            'used_is_marker': bool(np.any(mk > 0)),
                            'mean_weight': float(np.mean(W_ref)),
                            'max_weight': float(np.max(W_ref)),
                            'reference_gene_id_type': guess_gene_id_type(list(map(str, genes[: min(len(genes), 200)]))),
                        }
            out = (genes, A_base, W_ref, regions, meta)
            self._cache_ref[cache_key] = out
            return out
        finally:
            conn.close()

    def _load_sample_features(self, sample_id: str, genes: np.ndarray) -> pd.DataFrame:
        if genes.size == 0:
            return pd.DataFrame(index=pd.Index([], name='gene_symbol'))
        conn = self._connect()
        try:
            cols = ['gene_symbol']
            available = {r[1] for r in conn.execute('PRAGMA table_info(cfrna_expression)').fetchall()}
            for col in ['tpm_value', 'log_tpm', 'zscore_tpm', 'read_count', 'detected', 'gene_id_type']:
                if col in available:
                    cols.append(col)
            df = pd.read_sql_query(
                f"SELECT {', '.join(cols)} FROM cfrna_expression WHERE sample_id = ?",
                conn,
                params=[sample_id],
            )
            if df.empty:
                return pd.DataFrame(index=pd.Index(genes, name='gene_symbol'))
            df = df.dropna(subset=['gene_symbol']).copy()
            df['gene_symbol'] = df['gene_symbol'].astype(str)
            grouped = df.groupby('gene_symbol', as_index=True).agg({
                c: ('max' if c == 'detected' else 'mean') for c in df.columns if c != 'gene_symbol' and c != 'gene_id_type'
            })
            feat = grouped.reindex(genes).copy()
            if 'tpm_value' not in feat.columns:
                feat['tpm_value'] = 0.0
            feat['tpm_value'] = pd.to_numeric(feat['tpm_value'], errors='coerce').fillna(0.0)
            if 'log_tpm' not in feat.columns:
                feat['log_tpm'] = np.log1p(feat['tpm_value'].clip(lower=0.0))
            else:
                feat['log_tpm'] = pd.to_numeric(feat['log_tpm'], errors='coerce').fillna(np.log1p(feat['tpm_value'].clip(lower=0.0)))
            if 'zscore_tpm' not in feat.columns:
                y = feat['log_tpm'].to_numpy(dtype=float)
                feat['zscore_tpm'] = 0.0 if len(y) == 0 else (y - y.mean()) / (y.std() + 1e-8)
            else:
                fill_z = ((feat['log_tpm'] - feat['log_tpm'].mean()) / (feat['log_tpm'].std() + 1e-8)).fillna(0.0)
                feat['zscore_tpm'] = pd.to_numeric(feat['zscore_tpm'], errors='coerce').fillna(fill_z)
            if 'detected' not in feat.columns:
                feat['detected'] = (feat['tpm_value'] >= 1.0).astype(int)
            else:
                feat['detected'] = pd.to_numeric(feat['detected'], errors='coerce').fillna((feat['tpm_value'] >= 1.0).astype(int)).astype(float)
            if 'read_count' not in feat.columns:
                feat['read_count'] = np.nan
            else:
                feat['read_count'] = pd.to_numeric(feat['read_count'], errors='coerce')
            feat['rank_pct'] = feat['tpm_value'].rank(pct=True, method='average').fillna(0.0)
            if feat['read_count'].notna().any() and float(feat['read_count'].fillna(0).max()) > 0:
                rc = np.log1p(feat['read_count'].fillna(0.0).to_numpy(dtype=float))
                feat['read_support'] = rc / (float(rc.max()) + 1e-8)
            else:
                feat['read_support'] = np.nan
            gene_ids = df['gene_id_type'].dropna().astype(str).tolist() if 'gene_id_type' in df.columns else []
            feat.attrs['sample_gene_id_type'] = gene_ids[0] if gene_ids else guess_gene_id_type(df['gene_symbol'].tolist())
            return feat
        finally:
            conn.close()

    def _load_sample_bundle(
        self,
        sample_id: str,
        genes: np.ndarray,
        use_value: str,
        W_ref: np.ndarray,
        use_read_count_hint: bool,
        read_count_weight: float,
        min_weight: float,
        max_weight: float,
    ) -> Tuple[np.ndarray, np.ndarray, Dict[str, Any]]:
        feat = self._load_sample_features(sample_id, genes)
        if feat.empty:
            return np.array([], dtype=float), np.array([], dtype=float), {'error': 'sample not found'}
        vals = feat['tpm_value'].to_numpy(dtype=float)
        if use_value in {'tpm', 'vsd'}:
            b_base = vals
            expr_col = 'tpm_value'
        elif use_value == 'zscore':
            b_base = feat['zscore_tpm'].to_numpy(dtype=float)
            expr_col = 'zscore_tpm'
        else:
            b_base = feat['log_tpm'].to_numpy(dtype=float)
            expr_col = 'log_tpm'

        W_sample = np.max(W_ref, axis=1).astype(float) if W_ref.size > 0 else np.ones_like(b_base, dtype=float)
        meta: Dict[str, Any] = {
            'used_expr_col': expr_col,
            'used_read_count_hint': False,
            'used_detected_hint': bool(np.any(feat['detected'].to_numpy(dtype=float) > 0)),
            'sample_gene_id_type': feat.attrs.get('sample_gene_id_type', 'unknown'),
            'n_input_rows': int(len(feat)),
            'has_read_count': bool(feat['read_count'].notna().any()),
        }
        if use_read_count_hint and read_count_weight > 0 and feat['read_count'].notna().any():
            rc = np.log1p(feat['read_count'].fillna(0.0).to_numpy(dtype=float))
            denom = float(rc.max())
            if denom > 0:
                rc_norm = rc / denom
                W_sample = self._clip_array(W_sample * (1.0 + float(read_count_weight) * rc_norm), min_weight, max_weight)
                meta['used_read_count_hint'] = True

        overlap_n = int(np.sum(feat['tpm_value'].to_numpy(dtype=float) > 0))
        meta['n_overlap_genes'] = overlap_n
        meta['overlap_fraction'] = float(overlap_n / max(len(genes), 1))
        return b_base.astype(float), W_sample.astype(float), meta

    @staticmethod
    def _apply_weights(A_base: np.ndarray, b_base: np.ndarray, W_ref: np.ndarray, W_sample: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        return np.asarray(A_base, dtype=float) * np.asarray(W_ref, dtype=float), np.asarray(b_base, dtype=float) * np.asarray(W_sample, dtype=float)

    @staticmethod
    def _rank_correlation_scores(A_tpm: np.ndarray, b_tpm: np.ndarray) -> np.ndarray:
        if A_tpm.size == 0 or b_tpm.size == 0:
            return np.array([], dtype=float)
        b_rank = rankdata(np.asarray(b_tpm, dtype=float), method='average')
        scores = []
        for j in range(A_tpm.shape[1]):
            a_rank = rankdata(np.asarray(A_tpm[:, j], dtype=float), method='average')
            a0 = a_rank - a_rank.mean()
            b0 = b_rank - b_rank.mean()
            den = np.sqrt((a0 @ a0) * (b0 @ b0) + 1e-12)
            scores.append(float((a0 @ b0) / den) if den > 0 else 0.0)
        return np.nan_to_num(np.asarray(scores, dtype=float), nan=0.0, posinf=0.0, neginf=0.0)

    @staticmethod
    def _normalize_signal_weights(weight_map: Dict[str, float], availability: Dict[str, bool]) -> Dict[str, float]:
        kept = {k: float(v) for k, v in weight_map.items() if float(v) > 0 and availability.get(k, True)}
        total = float(sum(kept.values()))
        if total <= 0:
            kept = {k: 1.0 for k in weight_map if availability.get(k, True)}
            total = float(sum(kept.values()))
        return {k: v / total for k, v in kept.items()}

    @staticmethod
    def _blend_ensemble_weights(weight_map: Dict[str, float], availability: Dict[str, bool], ensemble_alpha: float) -> Dict[str, float]:
        """Blend corr/nnls according to ensemble_alpha while preserving auxiliary evidence."""
        alpha = float(np.clip(ensemble_alpha, 0.0, 1.0))
        base = {k: float(v) for k, v in weight_map.items() if float(v) > 0 and availability.get(k, True)}
        if not base:
            return {}
        corr_base = base.pop('corr', 0.0)
        nnls_base = base.pop('nnls', 0.0)
        primary_total = corr_base + nnls_base
        blended = dict(base)
        if primary_total > 0:
            if availability.get('corr', False):
                blended['corr'] = primary_total * alpha
            if availability.get('nnls', False):
                blended['nnls'] = primary_total * (1.0 - alpha)
            if 'corr' not in blended and availability.get('nnls', False):
                blended['nnls'] = primary_total
            if 'nnls' not in blended and availability.get('corr', False):
                blended['corr'] = primary_total
        return SourceTracingEngineV2._normalize_signal_weights(blended, availability)

    @staticmethod
    def _safe_margin(scores: np.ndarray) -> float:
        s = np.asarray(scores, dtype=float)
        if s.size < 2:
            return 1.0
        order = np.sort(s)[::-1]
        return float(order[0] - order[1])

    def _marker_panel_scores(
        self,
        sample_feat: pd.DataFrame,
        W_ref: np.ndarray,
        regions: List[str],
        topk: int,
    ) -> Dict[str, np.ndarray]:
        n_regions = len(regions)
        if sample_feat.empty or W_ref.size == 0 or n_regions == 0:
            zeros = np.zeros(n_regions, dtype=float)
            return {'marker': zeros, 'detect': zeros, 'support': zeros}

        z = sample_feat['zscore_tpm'].to_numpy(dtype=float)
        log_tpm = sample_feat['log_tpm'].to_numpy(dtype=float)
        detected = sample_feat['detected'].to_numpy(dtype=float)
        rank_pct = sample_feat['rank_pct'].to_numpy(dtype=float)
        read_support = sample_feat['read_support'].to_numpy(dtype=float) if 'read_support' in sample_feat.columns else np.full(len(sample_feat), np.nan)
        log_scaled = log_tpm / (float(np.nanmax(log_tpm)) + 1e-8) if len(log_tpm) else log_tpm

        marker_scores = np.zeros(n_regions, dtype=float)
        detect_scores = np.zeros(n_regions, dtype=float)
        support_scores = np.zeros(n_regions, dtype=float)

        for j in range(n_regions):
            region_w = np.asarray(W_ref[:, j], dtype=float)
            if region_w.size == 0:
                continue
            baseline = float(np.nanmedian(region_w))
            priority = np.maximum(region_w - baseline, 0.0)
            if np.all(priority <= 0):
                priority = region_w.copy()
            sel = np.argsort(priority)[::-1][: int(topk)]
            w = np.maximum(priority[sel], 1e-8)
            w = w / (float(w.sum()) + 1e-12)
            marker_signal = np.average(np.maximum(z[sel], 0.0) + 0.20 * log_scaled[sel], weights=w)
            detect_signal = np.average(np.clip(detected[sel], 0.0, 1.0), weights=w)
            if np.isfinite(read_support).any():
                support_mix = 0.70 * rank_pct[sel] + 0.30 * np.nan_to_num(read_support[sel], nan=0.0)
            else:
                support_mix = 0.75 * rank_pct[sel] + 0.25 * log_scaled[sel]
            support_signal = np.average(np.clip(support_mix, 0.0, 1.0), weights=w)
            marker_scores[j] = float(marker_signal)
            detect_scores[j] = float(detect_signal)
            support_scores[j] = float(support_signal)

        return {'marker': marker_scores, 'detect': detect_scores, 'support': support_scores}

    def _persist(
        self,
        run_id: str,
        sample_id: str,
        atlas_id: int,
        sigset_id: Optional[int],
        method: str,
        params_json: str,
        rows: List[Dict[str, Any]],
    ) -> bool:
        """将分析结果持久化到数据库。

        Returns:
            True 表示写入成功，False 表示表不存在或写入失败（已记录日志）。
        """
        conn = self._connect()
        try:
            if not self._table_exists(conn, 'analysis_runs'):
                logger.warning(
                    "analysis_runs table not found; skipping persist. "
                    "run_id=%s sample_id=%s — 请先初始化数据库 schema。",
                    run_id, sample_id,
                )
                return False
            if not self._table_exists(conn, 'analysis_results'):
                logger.warning(
                    "analysis_results table not found; skipping persist. "
                    "run_id=%s sample_id=%s",
                    run_id, sample_id,
                )
                return False
            cur = conn.cursor()
            cur.execute(
                "INSERT INTO analysis_runs"
                "(run_id, sample_id, atlas_id, sigset_id, method, params_json)"
                " VALUES (?, ?, ?, ?, ?, ?)",
                (
                    run_id,
                    sample_id,
                    int(atlas_id),
                    int(sigset_id) if sigset_id is not None else None,
                    method,
                    params_json,
                ),
            )
            cur.executemany(
                "INSERT OR REPLACE INTO analysis_results"
                "(run_id, region_id, score, fraction, ci_low, ci_high,"
                " stability, reconstruction_error, rank)"
                " VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)",
                [
                    (
                        run_id,
                        r['region_id'],
                        float(r['score']),
                        r.get('fraction'),
                        r.get('ci_low'),
                        r.get('ci_high'),
                        r.get('stability'),
                        r.get('reconstruction_error'),
                        int(r['rank']),
                    )
                    for r in rows
                ],
            )
            conn.commit()
            return True
        except Exception as exc:
            logger.error(
                "Failed to persist analysis results. run_id=%s sample_id=%s error=%r",
                run_id, sample_id, exc,
            )
            return False
        finally:
            conn.close()

    def trace(self, sample_id: str, method: str = 'correlation', sigset_id: Optional[int] = None, **params: Any) -> Dict[str, Any]:
        tp, ep = self._parse_params(params)
        persist_results = bool(params.get('persist', True))
        atlas_id = tp.atlas_id
        atlas_meta = self._get_atlas_meta(atlas_id)
        is_vsd_reference = self._is_vsd_reference(atlas_meta)
        vsd_compatible = bool(params.get('vsd_compatible', True)) and is_vsd_reference
        if vsd_compatible:
            if tp.use_value in ('tpm', 'log1p'):
                tp.use_value = 'zscore'
            ep.update(
                {
                    'corr_component_weight': float(params.get('corr_component_weight', 0.34)),
                    'nnls_component_weight': float(params.get('nnls_component_weight', 0.06)),
                    'rank_component_weight': float(params.get('rank_component_weight', 0.24)),
                    'marker_component_weight': float(params.get('marker_component_weight', 0.18)),
                    'detect_component_weight': float(params.get('detect_component_weight', 0.10)),
                    'support_component_weight': float(params.get('support_component_weight', 0.08)),
                    'use_read_count_hint': bool(params.get('use_read_count_hint', True)),
                }
            )
        resolved_sigset_id = sigset_id if sigset_id is not None else self._get_latest_sigset_id(atlas_id)

        genes, A_base, W_ref, regions, ref_meta = self._load_reference_bundle(
            atlas_id,
            resolved_sigset_id,
            tp.use_value,
            ep['enable_weighting'],
            ep['specificity_weight'],
            ep['marker_weight'],
            ep['marker_bonus'],
            ep['min_weight'],
            ep['max_weight'],
            ep['normalize_region_weight'],
            ep['fallback_marker_topk'],
        )
        b_base, W_sample, sample_meta = self._load_sample_bundle(
            sample_id,
            genes,
            tp.use_value,
            W_ref,
            ep['use_read_count_hint'],
            ep['read_count_weight'],
            ep['min_weight'],
            ep['max_weight'],
        )
        if b_base.size == 0 or A_base.size == 0:
            return {'sample_id': sample_id, 'method': method, 'run_id': None, 'results': [], 'meta': {'error': 'Empty sample/reference after alignment'}}

        overlap_n = int(sample_meta.get('n_overlap_genes', int(np.sum(b_base > 0))))
        overlap_fraction = float(sample_meta.get('overlap_fraction', overlap_n / max(len(genes), 1)))
        traceability = 'high'
        warnings: List[str] = []
        ref_gene_id_type = str(ref_meta.get('reference_gene_id_type', 'unknown'))
        sample_gene_id_type = str(sample_meta.get('sample_gene_id_type', 'unknown'))
        if vsd_compatible:
            warnings.append(
                'VSD-compatible tracing enabled: reference values are interpreted as normalized expression patterns, not TPM abundance.'
            )
        if ref_gene_id_type != 'unknown' and sample_gene_id_type != 'unknown' and ref_gene_id_type != sample_gene_id_type:
            warnings.append(f'Gene ID type mismatch: sample={sample_gene_id_type}, reference={ref_gene_id_type}')
            if overlap_n < max(ep['min_overlap_genes'] * 2, 40):
                traceability = 'low'
        if overlap_n < ep['min_overlap_genes'] or overlap_fraction < ep['min_overlap_fraction']:
            traceability = 'low'
            warnings.append(f'Low sample/reference overlap: overlap_genes={overlap_n}, overlap_fraction={overlap_fraction:.4f}')
        if overlap_n < max(3, ep['min_overlap_genes'] // 2):
            traceability = 'insufficient'
            warnings.append('Too few overlapping genes to support reliable tracing')
        if ep['abstain_on_low_overlap'] and traceability == 'insufficient':
            return {
                'sample_id': sample_id,
                'method': method,
                'run_id': None,
                'results': [],
                'meta': {
                    'atlas_id': atlas_id,
                    'sigset_id': resolved_sigset_id,
                    'traceability': traceability,
                    'sample_overlap': {'n_overlap_genes': overlap_n, 'overlap_fraction': overlap_fraction},
                    'warnings': warnings,
                },
            }

        A, b = self._apply_weights(A_base, b_base, W_ref, W_sample)
        keep = (np.sum(np.abs(A), axis=1) > 0) & np.isfinite(b)
        if np.sum(keep) == 0:
            return {'sample_id': sample_id, 'method': method, 'run_id': None, 'results': [], 'meta': {'error': 'No informative genes after weighting'}}

        genes = genes[keep]
        A = A[keep, :]
        b = b[keep]
        W_ref_k = W_ref[keep, :] if W_ref.size else W_ref
        regions = list(regions)
        sample_feat = self._load_sample_features(sample_id, genes)

        # 从已缓存的 ref_meta 取原始 TPM 矩阵，避免重复数据库 IO
        # ref_meta['_A_tpm_cache'] 在 _load_reference_bundle 中与 A_base 同步构建
        _A_tpm_cached = ref_meta.get('_A_tpm_cache')
        if _A_tpm_cached is not None and len(_A_tpm_cached):
            # 需要对齐到 keep 过滤后的 gene 顺序
            # genes_before_keep 在上方 keep 操作前已记录（通过 A_base 原始行索引）
            # 此处通过 genes（已 keep 过滤）直接切片即可，因为 keep 是对 genes 的 bool mask
            A_tpm = np.asarray(_A_tpm_cached, dtype=float)[keep, :]
        else:
            # fallback：仍走原有路径，保持向后兼容
            genes_raw, A_tpm_full, _regions_chk = self._load_reference_matrix(
                atlas_id=atlas_id, sigset_id=resolved_sigset_id, use_value='tpm'
            )
            if len(genes_raw) and not np.array_equal(genes_raw.astype(str), genes.astype(str)):
                ref_raw_df = pd.DataFrame(A_tpm_full, index=genes_raw.astype(str), columns=_regions_chk)
                A_tpm = ref_raw_df.reindex(index=genes.astype(str), columns=regions).fillna(0.0).to_numpy(dtype=float)
            else:
                A_tpm = np.asarray(A_tpm_full, dtype=float)
            if A_tpm.shape[0] != len(genes):
                ref_raw_df = pd.DataFrame(A_tpm_full, index=np.asarray(genes_raw).astype(str), columns=_regions_chk)
                A_tpm = ref_raw_df.reindex(index=genes.astype(str), columns=regions).fillna(0.0).to_numpy(dtype=float)

        m = canonical_method(method)

        ci = stability = None
        recon_err = None
        fractions = None
        components: Dict[str, np.ndarray] = {}
        component_weights: Dict[str, float] = {}

        corr_scores = self._trace_corr(A, b)
        if m == 'correlation':
            scores = corr_scores
            confidence = self._softmax_confidence(scores)
            components = {'corr': corr_scores}
        elif m == 'nnls_simplex':
            fractions, recon_err = self._trace_nnls_simplex(A, b, l2=tp.l2)
            scores = fractions.copy()
            confidence = self._softmax_confidence(scores)
            components = {'nnls': scores}
            if vsd_compatible:
                warnings.append(
                    'NNLS/simplex on a VSD reference is reported as a normalized-expression-space fitting weight, not a biological tissue fraction.'
                )
            if tp.bootstrap_n > 0:
                ci, stability = self._bootstrap_nnls(A, b, regions, tp.bootstrap_n, tp.bootstrap_gene_frac, tp.l2, tp.random_seed)
        elif m == 'ensemble':
            fractions, recon_err = self._trace_nnls_simplex(A, b, l2=tp.l2)
            rank_scores = self._rank_correlation_scores(A_tpm, sample_feat['tpm_value'].to_numpy(dtype=float))
            panel_scores = self._marker_panel_scores(sample_feat, W_ref_k, regions, topk=ep['marker_panel_topk'])
            has_region_weighting = bool(ref_meta.get('weighting_active', False))
            availability = {
                'corr': True,
                'nnls': True,
                'rank': bool(np.isfinite(rank_scores).any()),
                'marker': has_region_weighting and bool(np.isfinite(panel_scores['marker']).any()),
                'detect': has_region_weighting and bool(sample_meta.get('used_detected_hint')) and bool(np.isfinite(panel_scores['detect']).any()),
                'support': has_region_weighting and bool(np.isfinite(panel_scores['support']).any()),
            }
            weight_map = {
                'corr': ep['corr_component_weight'],
                'nnls': ep['nnls_component_weight'],
                'rank': ep['rank_component_weight'],
                'marker': ep['marker_component_weight'],
                'detect': ep['detect_component_weight'],
                'support': ep['support_component_weight'],
            }
            if vsd_compatible and not has_region_weighting:
                warnings.append(
                    'VSD atlas has no region-specific marker weighting; marker/detection/support components are unavailable. Prefer correlation for primary region ranking until weighted marker evidence is validated.'
                )
            component_weights = self._blend_ensemble_weights(weight_map, availability, tp.ensemble_alpha)
            z_components = {
                'corr': self._zscore_safe(corr_scores),
                'nnls': self._zscore_safe(fractions),
                'rank': self._zscore_safe(rank_scores),
                'marker': self._zscore_safe(panel_scores['marker']),
                'detect': self._zscore_safe(panel_scores['detect']),
                'support': self._zscore_safe(panel_scores['support']),
            }
            scores = np.zeros(len(regions), dtype=float)
            for key, wt in component_weights.items():
                scores = scores + float(wt) * z_components[key]
            confidence = self._softmax_confidence(scores)
            components = {
                'corr': corr_scores,
                'nnls': fractions,
                'rank': rank_scores,
                'marker': panel_scores['marker'],
                'detect': panel_scores['detect'],
                'support': panel_scores['support'],
            }
            if tp.bootstrap_n > 0:
                ci, stability = self._bootstrap_nnls(A, b, regions, tp.bootstrap_n, tp.bootstrap_gene_frac, tp.l2, tp.random_seed)
        else:
            return {'sample_id': sample_id, 'method': method, 'run_id': None, 'results': [], 'meta': {'error': f'Unsupported method: {method}'}}

        order = np.argsort(scores)[::-1]
        decision_margin = self._safe_margin(scores)
        effective_topk = len(order) if tp.return_all else tp.topk
        out_rows: List[Dict[str, Any]] = []
        for rank, j in enumerate(order[:effective_topk], start=1):
            idx = int(j)
            region_id = regions[idx]
            row: Dict[str, Any] = {
                'region_id': region_id,
                'rank': int(rank),
                'score': float(scores[idx]),
                'confidence': float(confidence[idx]),
            }
            if fractions is not None:
                row['fraction'] = float(fractions[idx])
            if recon_err is not None:
                row['reconstruction_error'] = float(recon_err)
            if ci is not None and region_id in ci:
                row['ci_low'] = float(ci[region_id][0])
                row['ci_high'] = float(ci[region_id][1])
            if stability is not None and region_id in stability:
                row['stability'] = float(stability[region_id])
            if components:
                for key, vec in components.items():
                    if len(vec) == len(regions):
                        row[f'{key}_component'] = float(vec[idx])
            out_rows.append(row)

        if ep['abstain_on_ambiguous'] and len(order) > 0:
            top1_conf = float(confidence[int(order[0])])
            if top1_conf < ep['min_top1_confidence'] or decision_margin < ep['min_decision_margin']:
                warnings.append(
                    f'Ambiguous prediction: top1_confidence={top1_conf:.4f}, decision_margin={decision_margin:.4f}'
                )
                return {
                    'sample_id': sample_id,
                    'method': m,
                    'run_id': None,
                    'results': [],
                    'meta': {
                        'atlas_id': atlas_id,
                        'sigset_id': resolved_sigset_id,
                        'traceability': 'ambiguous',
                        'sample_overlap': {'n_overlap_genes': overlap_n, 'overlap_fraction': overlap_fraction},
                        'warnings': warnings,
                        'top1_confidence': top1_conf,
                        'decision_margin': decision_margin,
                        'abstained_reason': 'ambiguous_prediction',
                    },
                }

        ref_meta_public = {k: v for k, v in ref_meta.items() if not str(k).startswith('_')}
        meta = {
            'atlas_id': atlas_id,
            'atlas_name': atlas_meta.get('atlas_name'),
            'atlas_normalization': atlas_meta.get('normalization'),
            'atlas_build_version': atlas_meta.get('build_version'),
            'vsd_compatible_mode': bool(vsd_compatible),
            'result_interpretation': (
                'Region-level normalized-expression pattern similarity; NNLS/simplex weights are fitting weights, not absolute RNA contribution fractions.'
                if vsd_compatible
                else 'TPM-compatible source tracing; NNLS/simplex fractions may be interpreted as mixture-style fitting fractions when sample/reference units are aligned.'
            ),
            'recommended_interpretation': (
                'Interpret Top regions as the brain areas whose reference expression fingerprints best match the cfRNA sample.'
                if vsd_compatible
                else 'Interpret Top regions using the selected scoring method and sample/reference normalization.'
            ),
            'sigset_id': resolved_sigset_id,
            'use_value': tp.use_value,
            'bootstrap_n': tp.bootstrap_n,
            'bootstrap_gene_frac': tp.bootstrap_gene_frac,
            'l2': tp.l2,
            'ensemble_alpha': tp.ensemble_alpha,
            'n_genes': int(len(genes)),
            'n_regions': int(len(regions)),
            'enable_weighting': ep['enable_weighting'],
            'specificity_weight': ep['specificity_weight'],
            'marker_weight': ep['marker_weight'],
            'marker_bonus': ep['marker_bonus'],
            'use_read_count_hint': ep['use_read_count_hint'],
            'read_count_weight': ep['read_count_weight'],
            'reference_weighting': ref_meta_public,
            'sample_weighting': sample_meta,
            'sample_overlap': {'n_overlap_genes': overlap_n, 'overlap_fraction': overlap_fraction},
            'traceability': traceability,
            'warnings': warnings,
            'signal_strategy': (
                ' + '.join(component_weights.keys())
                if m == 'ensemble'
                else m
            ),
            'signal_component_weights': component_weights,
            'decision_margin': decision_margin,
            'top1_confidence': float(confidence[int(order[0])]) if len(order) else None,
            'reference_gene_id_type': ref_gene_id_type,
            'sample_gene_id_type': sample_gene_id_type,
            'signal_inputs': {
                'uses_tpm_value': tp.use_value == 'tpm',
                'uses_vsd_value': tp.use_value == 'vsd',
                'uses_normalized_reference_expression': bool(vsd_compatible),
                'uses_log_tpm': True,
                'uses_zscore_tpm': True,
                'uses_detected': bool(sample_feat['detected'].to_numpy(dtype=float).sum() > 0),
                'uses_read_count': bool(sample_feat['read_count'].notna().any()),
            },
        }
        run_id = str(uuid.uuid4()) if persist_results else None
        persisted = False
        if persist_results and run_id is not None:
            persisted = self._persist(run_id, sample_id, atlas_id, resolved_sigset_id, m, json.dumps(meta, ensure_ascii=False), out_rows)
        meta['persisted'] = persisted
        meta['persist_requested'] = persist_results
        return {'sample_id': sample_id, 'method': m, 'run_id': run_id, 'results': out_rows, 'meta': meta}
