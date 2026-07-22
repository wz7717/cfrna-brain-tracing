"""
数据上传与预处理模块。
支持 cfRNA 数据导入、清洗、样本元数据抽取、表达标准化以及基础 QC。
"""

from __future__ import annotations

import io
import json
import logging
import sqlite3
from typing import Any, Dict, List, Tuple

import numpy as np
import pandas as pd

from core.gene_utils import guess_gene_id_type
from data.qc import compute_cohort_qc, compute_sample_qc as qc_compute_sample_qc, grade_sample_qc


logger = logging.getLogger(__name__)


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

    def __init__(self, db_path: str = "braintrace_source_tracing.db"):
        self.db_path = db_path

    def _get_conn(self):
        conn = sqlite3.connect(self.db_path)
        conn.execute("PRAGMA foreign_keys = ON")
        return conn

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
            "logtpm": "log_tpm_input",
            "log_tpm": "log_tpm_input",
            "log1p_tpm": "log_tpm_input",
            "logcpm": "logcpm_value",
            "log_cpm": "logcpm_value",
            "log2cpm": "logcpm_value",
            "log2_cpm": "logcpm_value",
            "raw_counts": "read_count",
            "raw_count": "read_count",
            "counts": "read_count",
            "reads": "read_count",
            "count": "read_count",
            "readcount": "read_count",
        }
        normalized_targets = set(df.columns)
        renames = {}
        for column in df.columns:
            normalized = mapping.get(str(column).strip().lower())
            if normalized and column != normalized and normalized not in normalized_targets:
                renames[column] = normalized
                normalized_targets.add(normalized)
        return df.rename(columns=renames)

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
        if "gene_symbol" not in df.columns:
            errors.append("缺少必需列: gene_symbol")
        value_cols = [col for col in ["read_count", "logcpm_value", "log_tpm_input", "tpm_value"] if col in df.columns]
        if not value_cols:
            errors.append("缺少表达值列: 推荐 raw_counts/count/read_count 或 logCPM；TPM/logTPM 仅作为 fallback")
        if errors:
            return False, errors

        if df["gene_symbol"].astype(str).str.strip().eq("").any() or df["gene_symbol"].isnull().any():
            errors.append("存在空的基因符号")

        for col in value_cols:
            values = pd.to_numeric(df[col], errors="coerce")
            if values.isnull().any():
                errors.append(f"存在空的 {col} 值或不可解析数值")
            if (values.fillna(0) < 0).any():
                errors.append(f"{col} 值不能为负数")
        if len(df) < 10:
            errors.append(f"数据量过少，至少需要 10 个基因（当前: {len(df)}）")
        return len(errors) == 0, errors

    def preprocess_expression_data(self, df: pd.DataFrame, min_tpm: float = 0.1, log_transform: bool = True) -> pd.DataFrame:
        df = self._normalize_columns(df.copy())
        df["gene_symbol"] = df["gene_symbol"].astype(str).str.strip()
        if "read_count" in df.columns:
            df["read_count"] = pd.to_numeric(df["read_count"], errors="coerce").fillna(0.0)

        expression_unit = "TPM_fallback"
        if "read_count" in df.columns and "tpm_value" not in df.columns and "logcpm_value" not in df.columns:
            total = float(df["read_count"].clip(lower=0).sum())
            if total <= 0:
                raise ValueError("raw counts/read_count 总和必须大于 0")
            cpm = df["read_count"].clip(lower=0) / total * 1_000_000.0
            df["tpm_value"] = np.log2(cpm + 1.0)
            expression_unit = "logCPM_from_raw_counts"
        elif "logcpm_value" in df.columns:
            df["tpm_value"] = pd.to_numeric(df["logcpm_value"], errors="coerce").fillna(0.0)
            expression_unit = "logCPM"
        elif "log_tpm_input" in df.columns and "tpm_value" not in df.columns:
            df["tpm_value"] = pd.to_numeric(df["log_tpm_input"], errors="coerce").fillna(0.0)
            expression_unit = "logTPM_fallback"
        else:
            df["tpm_value"] = pd.to_numeric(df["tpm_value"], errors="coerce").fillna(0.0)

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
        df["expression_unit"] = expression_unit
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
            cursor.execute("PRAGMA table_info(braintrace_samples)")
            cols = {row[1] for row in cursor.fetchall()}
            payload = self._sample_metadata_payload(metadata, cols)
            self._upsert_sample_metadata(cursor, payload)
            conn.commit()
        return metadata.get("sample_id")

    @staticmethod
    def _sample_metadata_payload(metadata: Dict, available_columns: set[str]) -> Dict:
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
        for key, value in extras.items():
            if key in available_columns:
                payload[key] = value
        return {key: value for key, value in payload.items() if key in available_columns}

    @staticmethod
    def _upsert_sample_metadata(conn: sqlite3.Connection | sqlite3.Cursor, payload: Dict) -> None:
        columns = list(payload)
        if "sample_id" not in columns or not payload.get("sample_id"):
            raise ValueError("sample_id is required")
        update_columns = [column for column in columns if column != "sample_id"]
        update_sql = ", ".join(f"{column}=excluded.{column}" for column in update_columns)
        sql = (
            f"INSERT INTO braintrace_samples ({', '.join(columns)}) "
            f"VALUES ({', '.join(['?'] * len(columns))}) "
            f"ON CONFLICT(sample_id) DO UPDATE SET {update_sql}"
        )
        conn.execute(sql, [payload[column] for column in columns])

    def save_sample_with_expression(self, metadata: Dict, df: pd.DataFrame) -> str:
        """Atomically replace one sample's metadata and expression matrix."""
        sample_id = str(metadata.get("sample_id") or "").strip()
        if not sample_id:
            raise ValueError("sample_id is required")

        expression = df.copy()
        required = {"gene_symbol", "tpm_value", "detected"}
        missing = sorted(required.difference(expression.columns))
        if missing:
            raise ValueError(f"Expression data is missing required columns: {', '.join(missing)}")
        expression["sample_id"] = sample_id

        conn = self._get_conn()
        try:
            conn.execute("BEGIN IMMEDIATE")
            sample_columns = {row[1] for row in conn.execute("PRAGMA table_info(braintrace_samples)").fetchall()}
            expression_columns = {row[1] for row in conn.execute("PRAGMA table_info(braintrace_expression)").fetchall()}

            payload = self._sample_metadata_payload(metadata, sample_columns)
            self._upsert_sample_metadata(conn, payload)

            columns_to_save = ["sample_id", "gene_symbol", "tpm_value", "detected"]
            for optional in ["read_count", "log_tpm", "zscore_tpm", "gene_id_type", "expression_unit"]:
                if optional in expression_columns and optional in expression.columns:
                    columns_to_save.append(optional)

            conn.execute("DELETE FROM braintrace_expression WHERE sample_id = ?", (sample_id,))
            placeholders = ", ".join(["?"] * len(columns_to_save))
            insert_sql = (
                f"INSERT INTO braintrace_expression ({', '.join(columns_to_save)}) "
                f"VALUES ({placeholders})"
            )

            def sqlite_value(value):
                if pd.isna(value):
                    return None
                return value.item() if isinstance(value, np.generic) else value

            rows = (
                tuple(sqlite_value(value) for value in row)
                for row in expression[columns_to_save].itertuples(index=False, name=None)
            )
            conn.executemany(insert_sql, rows)
            conn.commit()
        except Exception:
            conn.rollback()
            raise
        finally:
            conn.close()
        return sample_id

    def save_expression_data(self, sample_id: str, df: pd.DataFrame, *, run_qc: bool = True):
        df = df.copy()
        df["sample_id"] = sample_id
        with self._get_conn() as conn:
            cur = conn.cursor()
            cur.execute("PRAGMA table_info(braintrace_expression)")
            cols = {row[1] for row in cur.fetchall()}

        columns_to_save = ["sample_id", "gene_symbol", "tpm_value", "detected"]
        for opt in ["read_count", "log_tpm", "zscore_tpm", "gene_id_type", "expression_unit"]:
            if opt in cols and opt in df.columns:
                columns_to_save.append(opt)

        with self._get_conn() as conn:
            conn.execute("DELETE FROM braintrace_expression WHERE sample_id = ?", (sample_id,))
            df[columns_to_save].to_sql("braintrace_expression", conn, if_exists="append", index=False)
            conn.commit()

        # Legacy marker-panel QC remains available to existing callers, but can be
        # disabled by the current upload workflow because it is not part of the
        # manuscript tracing route and is not calibrated for a single new sample.
        if not run_qc:
            return
        try:
            qc = self.compute_sample_qc(df)
            self.save_sample_qc(sample_id, qc)
            grade = grade_sample_qc(qc)
            with self._get_conn() as conn2:
                gid = str(df["gene_id_type"].iloc[0]) if "gene_id_type" in df.columns and len(df) else None
                conn2.execute(
                    "UPDATE braintrace_samples SET qc_status = ?, gene_id_type = COALESCE(gene_id_type, ?), brain_traceability = COALESCE(brain_traceability, ?) WHERE sample_id = ?",
                    (grade, gid, grade, sample_id),
                )
                conn2.commit()
        except Exception:
            logger.exception("Legacy sample QC failed after expression storage for sample_id=%s", sample_id)

    def get_sample_expression(self, sample_id: str) -> pd.DataFrame:
        with self._get_conn() as conn:
            cursor = conn.cursor()
            cursor.execute("PRAGMA table_info(braintrace_expression)")
            available_cols = {row[1] for row in cursor.fetchall()}
            select_cols = ["gene_symbol", "tpm_value", "detected"]
            for optional_col in ["read_count", "log_tpm", "zscore_tpm", "expression_unit"]:
                if optional_col in available_cols:
                    select_cols.append(optional_col)
            return pd.read_sql_query(
                f"SELECT {', '.join(select_cols)} FROM braintrace_expression WHERE sample_id = ?",
                conn,
                params=[sample_id],
            )

    def get_sample_info(self, sample_id: str) -> Dict:
        with self._get_conn() as conn:
            df = pd.read_sql_query("SELECT * FROM braintrace_samples WHERE sample_id = ?", conn, params=[sample_id])
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
        conn = self._get_conn()
        try:
            columns = ["sample_id", "subject_id", "species", "diagnosis", "collection_date", "qc_status"]
            has_samples_table = conn.execute(
                "SELECT 1 FROM sqlite_master WHERE type = 'table' AND name = 'braintrace_samples'"
            ).fetchone()
            if not has_samples_table:
                return pd.DataFrame(columns=columns)
            return pd.read_sql_query(
                f"SELECT {', '.join(columns)} FROM braintrace_samples ORDER BY collection_date DESC",
                conn,
            )
        finally:
            conn.close()

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
            cursor.execute("DELETE FROM braintrace_expression WHERE sample_id = ?", [sample_id])
            cursor.execute("DELETE FROM source_tracing_results WHERE sample_id = ?", [sample_id])
            tables = {r[0] for r in cursor.execute("SELECT name FROM sqlite_master WHERE type='table'").fetchall()}
            if "analysis_runs" in tables:
                run_ids = [r[0] for r in cursor.execute("SELECT run_id FROM analysis_runs WHERE sample_id=?", [sample_id]).fetchall()]
                if run_ids:
                    cursor.executemany("DELETE FROM analysis_results WHERE run_id = ?", [(rid,) for rid in run_ids])
                cursor.execute("DELETE FROM analysis_runs WHERE sample_id = ?", [sample_id])
            if "sample_qc" in tables:
                cursor.execute("DELETE FROM sample_qc WHERE sample_id = ?", [sample_id])
            cursor.execute("DELETE FROM braintrace_samples WHERE sample_id = ?", [sample_id])
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
