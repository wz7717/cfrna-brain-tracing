"""
数据上传与预处理模块。
支持 cfRNA 数据导入、清洗、样本元数据抽取、表达标准化以及基础 QC。
"""

from __future__ import annotations

import io
import json
import sqlite3
from typing import Any, Dict, List, Tuple

import numpy as np
import pandas as pd

from core.gene_utils import guess_gene_id_type
from data.qc import compute_cohort_qc, compute_sample_qc as qc_compute_sample_qc, grade_sample_qc


class DataProcessor:
    FILE_METADATA_ALIASES = {
        "sample_id": ["sample_id", "sample", "sampleid"],
        "subject_id": ["subject_id", "subject", "subjectid", "animal_id"],
        "species": ["species"],
        "age_years": ["age_years", "age"],
        "sex": ["sex", "gender"],
        "diagnosis": ["diagnosis", "group"],
        "sample_type": ["sample_type"],
        "ground_truth_region": ["ground_truth_region", "source_region", "injury_region", "label_region", "true_source"],
        "ground_truth_region_name": ["ground_truth_region_name"],
        "source_type": ["source_type"],
        "surgery_region": ["surgery_region"],
        "surgery_side": ["surgery_side"],
        "post_op_day": ["post_op_day", "postop_day"],
        "collection_date": ["collection_date"],
        "plasma_volume_ml": ["plasma_volume_ml", "plasma_volume"],
        "extraction_method": ["extraction_method"],
        "rna_concentration_ng_ul": ["rna_concentration_ng_ul"],
        "rin_value": ["rin_value", "rin"],
        "library_preparation": ["library_preparation"],
        "sequencing_platform": ["sequencing_platform"],
        "total_reads": ["total_reads"],
        "mapped_reads": ["mapped_reads"],
        "mapping_rate": ["mapping_rate"],
        "gene_id_type": ["gene_id_type"],
        "brain_traceability": ["brain_traceability"],
    }

    def __init__(self, db_path: str = "cfrna_source_tracing.db"):
        self.db_path = db_path

    def _get_conn(self):
        return sqlite3.connect(self.db_path)

    def parse_expression_file(self, file_content: Any, file_format: str = "csv") -> pd.DataFrame:
        try:
            if file_format == "csv":
                payload = file_content if isinstance(file_content, str) else file_content.decode("utf-8-sig")
                df = pd.read_csv(io.StringIO(payload))
            elif file_format in ["tsv", "txt"]:
                payload = file_content if isinstance(file_content, str) else file_content.decode("utf-8-sig")
                df = pd.read_csv(io.StringIO(payload), sep="\t")
            elif file_format in ["excel", "xlsx"]:
                if isinstance(file_content, (bytes, bytearray)):
                    df = pd.read_excel(io.BytesIO(file_content))
                else:
                    raise ValueError("Excel 格式需要传入二进制内容。")
            else:
                raise ValueError(f"不支持的文件格式: {file_format}")
            return self._normalize_columns(df)
        except Exception as e:
            raise ValueError(f"文件解析失败: {e}")

    def _normalize_columns(self, df: pd.DataFrame) -> pd.DataFrame:
        mapping = {
            "gene": "gene_symbol",
            "Gene": "gene_symbol",
            "GENE": "gene_symbol",
            "gene_id": "gene_symbol",
            "gene_name": "gene_symbol",
            "Gene_Name": "gene_symbol",
            "tpm": "tpm_value",
            "TPM": "tpm_value",
            "fpkm": "tpm_value",
            "FPKM": "tpm_value",
            "reads": "read_count",
            "count": "read_count",
            "readcount": "read_count",
        }
        for old, new in mapping.items():
            if old in df.columns and new not in df.columns:
                df = df.rename(columns={old: new})
        return df

    def extract_embedded_metadata(self, df: pd.DataFrame) -> Dict[str, Any]:
        df = self._normalize_columns(df.copy())
        out: Dict[str, Any] = {}
        for target, candidates in self.FILE_METADATA_ALIASES.items():
            found = next((c for c in candidates if c in df.columns), None)
            if found is None:
                continue
            vals = df[found].dropna().astype(str).str.strip()
            vals = vals[vals != ""]
            if vals.empty:
                continue
            out[target] = vals.iloc[0]

        for num_key in [
            "age_years",
            "post_op_day",
            "plasma_volume_ml",
            "rna_concentration_ng_ul",
            "rin_value",
            "total_reads",
            "mapped_reads",
            "mapping_rate",
        ]:
            if num_key in out:
                try:
                    out[num_key] = float(out[num_key]) if "." in str(out[num_key]) else int(out[num_key])
                except Exception:
                    pass

        if "gene_id_type" not in out and "gene_symbol" in df.columns:
            out["gene_id_type"] = guess_gene_id_type(df["gene_symbol"].astype(str).tolist())
        if "sample_type" not in out:
            out["sample_type"] = "plasma"
        return out

    def validate_expression_data(self, df: pd.DataFrame) -> Tuple[bool, List[str]]:
        df = self._normalize_columns(df.copy())
        errors = []
        for col in ["gene_symbol", "tpm_value"]:
            if col not in df.columns:
                errors.append(f"缺少必需列: {col}")
        if errors:
            return False, errors

        if df["gene_symbol"].astype(str).str.strip().eq("").any() or df["gene_symbol"].isnull().any():
            errors.append("存在空的基因符号")

        df["tpm_value"] = pd.to_numeric(df["tpm_value"], errors="coerce")
        if df["tpm_value"].isnull().any():
            errors.append("存在空的 TPM 值或不可解析数值")
        if (df["tpm_value"].fillna(0) < 0).any():
            errors.append("TPM 值不能为负数")
        if len(df) < 10:
            errors.append(f"数据量过少，至少需要 10 个基因（当前: {len(df)}）")
        return len(errors) == 0, errors

    def preprocess_expression_data(self, df: pd.DataFrame, min_tpm: float = 0.1, log_transform: bool = True) -> pd.DataFrame:
        df = self._normalize_columns(df.copy())
        df["gene_symbol"] = df["gene_symbol"].astype(str).str.strip()
        df["tpm_value"] = pd.to_numeric(df["tpm_value"], errors="coerce").fillna(0.0)
        if "read_count" in df.columns:
            df["read_count"] = pd.to_numeric(df["read_count"], errors="coerce").fillna(0).astype(int)

        df = df[df["tpm_value"] >= float(min_tpm)].copy()
        df = df.groupby("gene_symbol", as_index=False).agg({
            "tpm_value": "mean",
            **({"read_count": "max"} if "read_count" in df.columns else {}),
        })

        df["log_tpm"] = np.log1p(df["tpm_value"].clip(lower=0)) if log_transform else df["tpm_value"]
        std = float(df["log_tpm"].std()) if len(df) else 0.0
        df["zscore_tpm"] = 0.0 if std == 0 else (df["log_tpm"] - df["log_tpm"].mean()) / std
        df["detected"] = (df["tpm_value"] >= 1.0).astype(int)
        df["gene_id_type"] = guess_gene_id_type(df["gene_symbol"].tolist())
        df["expression_unit"] = "TPM"
        return df

    def compute_sample_qc(self, df: pd.DataFrame) -> Dict[str, float]:
        return qc_compute_sample_qc(df)

    def save_sample_qc(self, sample_id: str, qc: Dict[str, float]) -> None:
        with self._get_conn() as conn:
            cur = conn.cursor()
            cur.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='sample_qc'")
            if cur.fetchone() is None:
                return
            cur.execute(
                "INSERT OR REPLACE INTO sample_qc(sample_id, hemolysis_hbb_hba_ratio, immune_ptprc, albumin_alb, brain_signal_score) VALUES (?, ?, ?, ?, ?)",
                (
                    sample_id,
                    None if pd.isna(qc.get("hemolysis_hbb_hba_ratio", np.nan)) else float(qc.get("hemolysis_hbb_hba_ratio")),
                    None if pd.isna(qc.get("immune_ptprc", np.nan)) else float(qc.get("immune_ptprc")),
                    None if pd.isna(qc.get("albumin_alb", np.nan)) else float(qc.get("albumin_alb")),
                    None if pd.isna(qc.get("brain_signal_score", np.nan)) else float(qc.get("brain_signal_score")),
                ),
            )
            conn.commit()

    def save_sample_metadata(self, metadata: Dict) -> str:
        with self._get_conn() as conn:
            cursor = conn.cursor()
            cursor.execute("PRAGMA table_info(cfrna_samples)")
            cols = {row[1] for row in cursor.fetchall()}
            payload = {
                "sample_id": metadata.get("sample_id"),
                "subject_id": metadata.get("subject_id"),
                "species": metadata.get("species"),
                "age_years": metadata.get("age_years"),
                "sex": metadata.get("sex"),
                "diagnosis": metadata.get("diagnosis"),
                "csf_volume_ml": metadata.get("csf_volume_ml"),
                "collection_date": metadata.get("collection_date"),
                "extraction_method": metadata.get("extraction_method"),
                "rna_concentration_ng_ul": metadata.get("rna_concentration_ng_ul"),
                "rin_value": metadata.get("rin_value"),
                "library_preparation": metadata.get("library_preparation"),
                "sequencing_platform": metadata.get("sequencing_platform"),
                "total_reads": metadata.get("total_reads"),
                "mapped_reads": metadata.get("mapped_reads"),
                "mapping_rate": metadata.get("mapping_rate"),
                "qc_status": metadata.get("qc_status", "Pending"),
                "metadata": json.dumps(metadata, ensure_ascii=False),
            }
            extras = {
                "plasma_volume_ml": metadata.get("plasma_volume_ml"),
                "sample_type": metadata.get("sample_type", "plasma"),
                "gene_id_type": metadata.get("gene_id_type"),
                "brain_traceability": metadata.get("brain_traceability"),
                "post_op_day": metadata.get("post_op_day"),
                "surgery_region": metadata.get("surgery_region"),
                "surgery_side": metadata.get("surgery_side"),
            }
            for k, v in extras.items():
                if k in cols:
                    payload[k] = v
            cols_to_insert = list(payload.keys())
            values = [payload[c] for c in cols_to_insert]
            cursor.execute(
                f"INSERT OR REPLACE INTO cfrna_samples ({', '.join(cols_to_insert)}) VALUES ({','.join(['?'] * len(cols_to_insert))})",
                values,
            )
            conn.commit()
        return metadata.get("sample_id")

    def save_expression_data(self, sample_id: str, df: pd.DataFrame):
        df = df.copy()
        df["sample_id"] = sample_id
        with self._get_conn() as conn:
            cur = conn.cursor()
            cur.execute("PRAGMA table_info(cfrna_expression)")
            cols = {row[1] for row in cur.fetchall()}

        columns_to_save = ["sample_id", "gene_symbol", "tpm_value", "detected"]
        for opt in ["read_count", "log_tpm", "zscore_tpm", "gene_id_type", "expression_unit"]:
            if opt in cols and opt in df.columns:
                columns_to_save.append(opt)

        with self._get_conn() as conn:
            conn.execute("DELETE FROM cfrna_expression WHERE sample_id = ?", (sample_id,))
            df[columns_to_save].to_sql("cfrna_expression", conn, if_exists="append", index=False)
            conn.commit()

        try:
            qc = self.compute_sample_qc(df)
            self.save_sample_qc(sample_id, qc)
            grade = grade_sample_qc(qc)
            with self._get_conn() as conn2:
                gid = str(df["gene_id_type"].iloc[0]) if "gene_id_type" in df.columns and len(df) else None
                conn2.execute(
                    "UPDATE cfrna_samples SET qc_status = ?, gene_id_type = COALESCE(gene_id_type, ?), brain_traceability = COALESCE(brain_traceability, ?) WHERE sample_id = ?",
                    (grade, gid, grade, sample_id),
                )
                conn2.commit()
        except Exception:
            pass

    def get_sample_expression(self, sample_id: str) -> pd.DataFrame:
        with self._get_conn() as conn:
            return pd.read_sql_query(
                "SELECT gene_symbol, tpm_value, detected FROM cfrna_expression WHERE sample_id = ?",
                conn,
                params=[sample_id],
            )

    def get_sample_info(self, sample_id: str) -> Dict:
        with self._get_conn() as conn:
            df = pd.read_sql_query("SELECT * FROM cfrna_samples WHERE sample_id = ?", conn, params=[sample_id])
        if len(df) == 0:
            return None
        info = df.iloc[0].to_dict()
        if "metadata" in info and isinstance(info["metadata"], str):
            try:
                info["metadata"] = json.loads(info["metadata"])
            except Exception:
                pass
        return info

    def get_all_samples(self) -> pd.DataFrame:
        with self._get_conn() as conn:
            return pd.read_sql_query(
                "SELECT sample_id, subject_id, species, diagnosis, collection_date, qc_status FROM cfrna_samples ORDER BY collection_date DESC",
                conn,
            )

    def compute_database_cohort_qc(self) -> pd.DataFrame:
        samples_df = self.get_all_samples()
        if samples_df.empty:
            return pd.DataFrame()

        sample_map = {}
        for sample_id in samples_df["sample_id"].astype(str).tolist():
            expr_df = self.get_sample_expression(sample_id)
            if expr_df is not None and not expr_df.empty:
                sample_map[sample_id] = expr_df

        if not sample_map:
            return pd.DataFrame()

        cohort_qc = compute_cohort_qc(sample_map)
        rows = []
        for sample_id, qc in cohort_qc.items():
            rows.append(
                {
                    "sample_id": sample_id,
                    "overall_risk": qc.get("overall_risk"),
                    "gene_id_type": qc.get("gene_id_type"),
                    "rbc_score": qc.get("rbc_mrna_score"),
                    "rbc_percentile": qc.get("rbc_mrna_percentile"),
                    "rbc_risk": qc.get("rbc_mrna_risk"),
                    "immune_score": qc.get("immune_mrna_score"),
                    "immune_percentile": qc.get("immune_mrna_percentile"),
                    "immune_risk": qc.get("immune_mrna_risk"),
                    "brain_score": qc.get("brain_marker_score"),
                    "brain_percentile": qc.get("brain_marker_percentile"),
                    "brain_risk": qc.get("brain_marker_risk"),
                    "hemolysis_mirna_risk": qc.get("hemolysis_mirna_risk"),
                    "interpretation": qc.get("interpretation"),
                }
            )

        qc_df = pd.DataFrame(rows)
        if qc_df.empty:
            return qc_df
        return samples_df.merge(qc_df, on="sample_id", how="left")

    def delete_sample(self, sample_id: str):
        with self._get_conn() as conn:
            cursor = conn.cursor()
            cursor.execute("DELETE FROM cfrna_expression WHERE sample_id = ?", [sample_id])
            cursor.execute("DELETE FROM source_tracing_results WHERE sample_id = ?", [sample_id])
            tables = {r[0] for r in cursor.execute("SELECT name FROM sqlite_master WHERE type='table'").fetchall()}
            if "analysis_runs" in tables:
                run_ids = [r[0] for r in cursor.execute("SELECT run_id FROM analysis_runs WHERE sample_id=?", [sample_id]).fetchall()]
                if run_ids:
                    cursor.executemany("DELETE FROM analysis_results WHERE run_id = ?", [(rid,) for rid in run_ids])
                cursor.execute("DELETE FROM analysis_runs WHERE sample_id = ?", [sample_id])
            if "sample_qc" in tables:
                cursor.execute("DELETE FROM sample_qc WHERE sample_id = ?", [sample_id])
            cursor.execute("DELETE FROM cfrna_samples WHERE sample_id = ?", [sample_id])
            conn.commit()

    def generate_qc_report(self, sample_id: str) -> Dict:
        df = self.get_sample_expression(sample_id)
        report = {"sample_id": sample_id, "basic_stats": {}, "warnings": [], "status": "Pass"}
        if df is None or len(df) == 0:
            report["status"] = "Fail"
            report["warnings"].append("未找到表达数据。")
            return report

        report["basic_stats"] = {
            "total_genes": len(df),
            "detected_genes": int(df["detected"].sum()),
            "detection_rate": float(df["detected"].mean() * 100),
            "mean_tpm": float(df["tpm_value"].mean()),
            "median_tpm": float(df["tpm_value"].median()),
        }

        qc = self.compute_sample_qc(df)
        report["status"] = grade_sample_qc(qc)
        if not int(qc.get("qc_applicable", 0)):
            report["warnings"].append("当前基因 ID 不是 symbol-like，无法进行基于基因符号的 QC 风险评估。")

        overall_interp = str(qc.get("interpretation", "")).strip()
        if overall_interp:
            report["warnings"].append(overall_interp)

        panel_messages = [
            ("hemolysis_mirna_risk", qc.get("mir451a_mir23a_ratio_interpretation")),
            ("rbc_mrna_risk", qc.get("rbc_mrna_interpretation")),
            ("immune_mrna_risk", qc.get("immune_mrna_interpretation")),
            ("brain_marker_risk", qc.get("brain_marker_interpretation")),
        ]
        for risk_key, message in panel_messages:
            risk_value = qc.get(risk_key)
            message = str(message or "").strip()
            if risk_value in {"Moderate risk", "High risk", "Uncalibrated"} and message:
                if message not in report["warnings"]:
                    report["warnings"].append(message)
        return report
