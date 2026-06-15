"""
血浆cfRNA溯源分析兼容层（3.1）

本模块在 3.1 中只承担两类职责：
1. 对旧版前端/API 保持兼容；
2. 将核心算法统一收口到 SourceTracingEngineV2。

说明：
- 真正的 tracing 数学逻辑以 source_tracing_v2.py 为唯一事实来源；
- 本文件仍保留 legacy 风格返回结构，以免旧页面和旧脚本失效。
"""

from __future__ import annotations

import json
import sqlite3
from typing import Dict, Optional

import numpy as np
import pandas as pd

from source_tracing_v2 import SourceTracingEngineV2


class CSFRNASourceTracer:
    """血浆cfRNA溯源分析器（兼容包装器）"""

    def __init__(self, db_path: str = "cfrna_source_tracing.db"):
        self.db_path = db_path
        self.reference_expression: Optional[pd.DataFrame] = None
        self.region_signatures: Optional[pd.DataFrame] = None
        self.engine_v2 = SourceTracingEngineV2(db_path=db_path)

    def _connect(self) -> sqlite3.Connection:
        return sqlite3.connect(self.db_path)

    def load_reference_data(self):
        """加载参考数据到内存，供 legacy 页面与 marker 方法使用。"""
        with self._connect() as conn:
            self.reference_expression = pd.read_sql_query(
                "SELECT region_id, region_name, gene_symbol, avg_tpm, std_tpm, expression_class FROM reference_expression",
                conn,
            )
            self.region_signatures = pd.read_sql_query(
                "SELECT region_id, region_name, gene_symbol, marker_score, specificity_score, expression_level, is_marker FROM region_gene_signature",
                conn,
            )

    # ---------- New unified public API ----------
    def trace_source(self, sample_id: str, method: str = "ensemble", **params) -> Dict:
        """唯一推荐入口：直接调用 v2 引擎。"""
        return self.engine_v2.trace(sample_id=sample_id, method=method, **params)

    # ---------- Legacy result adapters ----------
    def _vector_from_df(self, cfrna_data: pd.DataFrame, genes: np.ndarray, use_value: str = "log1p") -> np.ndarray:
        c = cfrna_data.groupby("gene_symbol")["tpm_value"].mean()
        b = c.reindex(genes).fillna(0.0).values.astype(float)
        if use_value == "tpm":
            return b
        if use_value == "zscore":
            y = np.log1p(np.clip(b, 0, None))
            return (y - y.mean()) / (y.std() + 1e-8)
        return np.log1p(np.clip(b, 0, None))

    def correlation_based_tracing(self, cfrna_data: pd.DataFrame, top_regions: int = 5, atlas_id: int = 1, sigset_id: Optional[int] = None) -> Dict:
        genes, A, regions = self.engine_v2._load_reference_matrix(atlas_id=atlas_id, sigset_id=sigset_id or self.engine_v2._get_latest_sigset_id(atlas_id), use_value="log1p")
        results = {
            "method": "Correlation-Based Tracing (v2)",
            "top_regions": [],
            "region_scores": {},
            "marker_genes": {},
            "confidence": 0.0,
        }
        if len(genes) == 0:
            return results
        b = self._vector_from_df(cfrna_data, genes, use_value="log1p")
        corr = self.engine_v2._trace_corr(A, b)
        for j, region in enumerate(regions):
            results["region_scores"][region] = {
                "correlation": float(corr[j]),
                "p_value": None,
                "n_genes": int(len(genes)),
                "confidence": float(abs(corr[j])),
            }
        sorted_regions = sorted(results["region_scores"].items(), key=lambda x: x[1]["correlation"], reverse=True)
        results["top_regions"] = sorted_regions[:top_regions]
        if results["top_regions"]:
            results["confidence"] = float(np.mean([r[1]["confidence"] for r in results["top_regions"]]))
        return results

    def deconvolution_based_tracing(self, cfrna_data: pd.DataFrame, method: str = "NMF", n_components: Optional[int] = None, atlas_id: int = 1, sigset_id: Optional[int] = None) -> Dict:
        """旧接口保留：NMF/LS 统一收敛到 v2 的 simplex-NNLS。"""
        resolved_method = "nnls_simplex" if method in ("NMF", "LS") else method
        genes, A, regions = self.engine_v2._load_reference_matrix(atlas_id=atlas_id, sigset_id=sigset_id or self.engine_v2._get_latest_sigset_id(atlas_id), use_value="log1p")
        results = {
            "method": f"{method}-Based Deconvolution (v2)",
            "components": {},
            "contributions": {},
            "reconstruction_error": 0.0,
        }
        if len(genes) == 0:
            return results
        b = self._vector_from_df(cfrna_data, genes, use_value="log1p")
        w, rmse = self.engine_v2._trace_nnls_simplex(A, b, l2=1e-4)
        for i, region in enumerate(regions):
            results["contributions"][region] = float(w[i])
        results["reconstruction_error"] = float(rmse)
        sorted_contributions = sorted(results["contributions"].items(), key=lambda x: x[1], reverse=True)
        results["components"] = dict(sorted_contributions[:10])
        return results

    def marker_gene_based_tracing(self, cfrna_data: pd.DataFrame, n_markers: int = 20) -> Dict:
        """保留 marker-only legacy 方法，用于对比/诊断。"""
        if self.region_signatures is None:
            self.load_reference_data()
        results = {
            "method": "Marker Gene-Based Tracing",
            "region_scores": {},
            "top_markers": {},
            "confidence": 0.0,
        }
        top_markers = self.region_signatures[self.region_signatures["is_marker"] == 1].sort_values(
            "specificity_score", ascending=False
        ).groupby("region_id").head(n_markers).reset_index(drop=True)
        regions = top_markers["region_id"].unique()
        for region in regions:
            region_markers = top_markers[top_markers["region_id"] == region]
            merged = cfrna_data.merge(region_markers[["gene_symbol", "specificity_score", "marker_score"]], on="gene_symbol", how="inner")
            if len(merged) == 0:
                continue
            weighted_score = np.sum(merged["tpm_value"].values * merged["specificity_score"].values * merged["marker_score"].values)
            normalized_score = weighted_score / max(len(merged), 1)
            results["region_scores"][region] = {
                "score": float(normalized_score),
                "n_markers_found": int(len(merged)),
                "top_marker_genes": merged.nlargest(5, "tpm_value")["gene_symbol"].tolist(),
            }
        sorted_regions = sorted(results["region_scores"].items(), key=lambda x: x[1]["score"], reverse=True)
        if sorted_regions:
            max_score = sorted_regions[0][1]["score"] or 1.0
            for region in results["region_scores"]:
                results["region_scores"][region]["normalized_score"] = results["region_scores"][region]["score"] / max_score
            results["confidence"] = sorted_regions[0][1]["normalized_score"]
        return results

    def integrated_tracing(self, cfrna_data: pd.DataFrame, atlas_id: int = 1, sigset_id: Optional[int] = None, ensemble_alpha: float = 0.5) -> Dict:
        """旧版“集成分析”外观保留，但核心统一为 v2 ensemble。"""
        genes, A, regions = self.engine_v2._load_reference_matrix(atlas_id=atlas_id, sigset_id=sigset_id or self.engine_v2._get_latest_sigset_id(atlas_id), use_value="log1p")
        if len(genes) == 0:
            return {
                "method": "Integrated Multi-Method Tracing (v2)",
                "final_ranking": [],
                "integrated_scores": {},
                "method_results": {},
                "top_source": None,
                "overall_confidence": 0.0,
            }
        b = self._vector_from_df(cfrna_data, genes, use_value="log1p")
        corr = self.engine_v2._trace_corr(A, b)
        w, rmse = self.engine_v2._trace_nnls_simplex(A, b, l2=1e-4)
        alpha = float(np.clip(ensemble_alpha, 0.0, 1.0))
        scores = alpha * self.engine_v2._zscore_safe(corr) + (1.0 - alpha) * self.engine_v2._zscore_safe(w)
        integrated_scores = {regions[i]: float(scores[i]) for i in range(len(regions))}
        max_score = max(integrated_scores.values()) if integrated_scores else 1.0
        if max_score != 0:
            for k in list(integrated_scores.keys()):
                integrated_scores[k] = integrated_scores[k] / max_score
        final_ranking = sorted(integrated_scores.items(), key=lambda x: x[1], reverse=True)
        return {
            "method": "Integrated Multi-Method Tracing (v2)",
            "final_ranking": final_ranking[:10],
            "integrated_scores": integrated_scores,
            "method_results": {
                "correlation_vectorized": {"corr": {regions[i]: float(corr[i]) for i in range(len(regions))}},
                "nnls_simplex": {"contributions": {regions[i]: float(w[i]) for i in range(len(regions))}, "reconstruction_error": float(rmse)},
            },
            "top_source": final_ranking[0] if final_ranking else None,
            "overall_confidence": np.mean([r[1] for r in final_ranking[:3]]) if final_ranking else 0,
        }

    def save_tracing_results(self, sample_id: str, results: Dict):
        """保存 legacy 结果到旧表，保证历史页面可用。"""
        with self._connect() as conn:
            cursor = conn.cursor()
            if "final_ranking" in results:
                top_regions = [r[0] for r in results["final_ranking"]]
                top_markers = top_regions[:5]
                confidence = results.get("overall_confidence", 0)
                contributions = results.get("integrated_scores", {})
            else:
                if "region_scores" in results:
                    top_regions = list(results["region_scores"].keys())[:5]
                    top_markers = top_regions
                    confidence = results.get("confidence", 0)
                    contributions = {k: v.get("correlation", v.get("score", 0)) for k, v in results["region_scores"].items()}
                elif "contributions" in results:
                    top_regions = sorted(results["contributions"].keys(), key=lambda x: results["contributions"][x], reverse=True)[:5]
                    top_markers = top_regions
                    confidence = 0
                    contributions = results["contributions"]
                else:
                    top_regions, top_markers, confidence, contributions = [], [], 0, {}
            cursor.execute(
                """
                INSERT INTO source_tracing_results
                (sample_id, analysis_date, analysis_method, top_source_regions,
                 region_contributions, confidence_score, marker_genes_used,
                 cross_validation_score, results_json)
                VALUES (?, datetime('now'), ?, ?, ?, ?, ?, ?, ?)
                """,
                (
                    sample_id,
                    results["method"],
                    json.dumps(top_regions),
                    json.dumps(contributions),
                    confidence,
                    json.dumps(top_markers),
                    0.0,
                    json.dumps(results, ensure_ascii=False, indent=2),
                ),
            )
            conn.commit()

    def get_sample_results(self, sample_id: str) -> Optional[Dict]:
        with self._connect() as conn:
            df = pd.read_sql_query(
                "SELECT sample_id, analysis_date, analysis_method, top_source_regions, region_contributions, confidence_score, results_json FROM source_tracing_results WHERE sample_id = ?",
                conn,
                params=[sample_id],
            )
        if len(df) == 0:
            return None
        result = df.iloc[0].to_dict()
        result["top_source_regions"] = json.loads(result["top_source_regions"])
        result["region_contributions"] = json.loads(result["region_contributions"])
        result["results_json"] = json.loads(result["results_json"])
        return result

    def close(self):
        """短连接模式下无需显式关闭，保留兼容接口。"""
        return None


if __name__ == "__main__":
    tracer = CSFRNASourceTracer()
    tracer.load_reference_data()
    print("3.1 compatibility wrapper ready.")
