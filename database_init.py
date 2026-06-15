"""
基于猕猴全脑转录组图谱的血浆cfRNA溯源数据库
数据库初始化脚本
"""

import sqlite3
import pandas as pd
import numpy as np
from pathlib import Path
from typing import Dict, List, Optional
import json
from data.migrations import run_migrations

class CSFRNASourceDatabase:
    """血浆cfRNA溯源数据库"""

    def __init__(self, db_path: str = "cfrna_source_tracing.db"):
        self.db_path = db_path
        self.conn = None

    def connect(self):
        """连接数据库"""
        self.conn = sqlite3.connect(self.db_path)
        return self.conn

    def create_database_schema(self):
        """创建数据库表结构"""
        cursor = self.conn.cursor()

        # 1. 猕猴全脑转录组图谱表
        cursor.execute("""
            CREATE TABLE IF NOT EXISTS macaque_brain_atlas (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                region_id TEXT NOT NULL,
                region_name TEXT NOT NULL,
                region_acronym TEXT,
                parent_region_id TEXT,
                hemi TEXT,
                layer TEXT,
                atlas_version TEXT,
                coordinates TEXT,
                UNIQUE(region_id, layer, hemi)
            )
        """)

        # 2. 参考基因表达谱表
        cursor.execute("""
            CREATE TABLE IF NOT EXISTS reference_expression (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                gene_symbol TEXT NOT NULL,
                gene_name TEXT,
                ensembl_id TEXT,
                ncbi_id TEXT,
                region_id TEXT NOT NULL,
                region_name TEXT NOT NULL,
                avg_tpm REAL,
                std_tpm REAL,
                median_tpm REAL,
                sample_count INTEGER,
                expression_class TEXT,
                cell_type_marker TEXT,
                FOREIGN KEY (region_id) REFERENCES macaque_brain_atlas(region_id)
            )
        """)

        # 3. 血浆cfRNA样本表
        cursor.execute("""
            CREATE TABLE IF NOT EXISTS cfrna_samples (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                sample_id TEXT NOT NULL UNIQUE,
                subject_id TEXT,
                species TEXT,
                age_years REAL,
                sex TEXT,
                diagnosis TEXT,
                csf_volume_ml REAL,
                collection_date TEXT,
                extraction_method TEXT,
                rna_concentration_ng_ul REAL,
                rin_value REAL,
                library_preparation TEXT,
                sequencing_platform TEXT,
                total_reads INTEGER,
                mapped_reads INTEGER,
                mapping_rate REAL,
                qc_status TEXT,
                metadata TEXT
            )
        """)

        # 4. cfRNA基因表达表
        cursor.execute("""
            CREATE TABLE IF NOT EXISTS cfrna_expression (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                sample_id TEXT NOT NULL,
                gene_symbol TEXT NOT NULL,
                tpm_value REAL,
                read_count INTEGER,
                detected INTEGER,
                FOREIGN KEY (sample_id) REFERENCES cfrna_samples(sample_id)
            )
        """)

        # 5. 溯源分析结果表
        cursor.execute("""
            CREATE TABLE IF NOT EXISTS source_tracing_results (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                sample_id TEXT NOT NULL,
                analysis_date TEXT,
                analysis_method TEXT,
                top_source_regions TEXT,
                region_contributions TEXT,
                confidence_score REAL,
                marker_genes_used TEXT,
                cross_validation_score REAL,
                results_json TEXT,
                FOREIGN KEY (sample_id) REFERENCES cfrna_samples(sample_id)
            )
        """)

        # 6. 脑区-基因关联表（用于溯源）
        cursor.execute("""
            CREATE TABLE IF NOT EXISTS region_gene_signature (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                region_id TEXT NOT NULL,
                region_name TEXT NOT NULL,
                gene_symbol TEXT NOT NULL,
                marker_score REAL,
                specificity_score REAL,
                expression_level TEXT,
                is_marker INTEGER,
                FOREIGN KEY (region_id) REFERENCES macaque_brain_atlas(region_id)
            )
        """)

        # 7. 疾病关联表
        cursor.execute("""
            CREATE TABLE IF NOT EXISTS disease_associations (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                disease_name TEXT NOT NULL,
                source_regions TEXT,
                dysregulated_genes TEXT,
                enrichment_score REAL,
                evidence_level TEXT
            )
        """)

        # 8. 分析历史表
        cursor.execute("""
            CREATE TABLE IF NOT EXISTS analysis_history (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                user_id TEXT,
                analysis_type TEXT,
                parameters TEXT,
                results_summary TEXT,
                timestamp TEXT
            )
        """)

        # =============================
        # Publish-grade extensions (backward compatible)
        # =============================
        # 9. Atlas versions (optional, for reproducibility)
        cursor.execute("""
            CREATE TABLE IF NOT EXISTS atlas_versions (
                atlas_id INTEGER PRIMARY KEY,
                atlas_name TEXT NOT NULL,
                species TEXT,
                level TEXT,
                build_version TEXT,
                gene_id_type TEXT,
                normalization TEXT,
                created_at TEXT DEFAULT (datetime('now')),
                notes TEXT
            )
        """)

        # 10. Signature sets
        cursor.execute("""
            CREATE TABLE IF NOT EXISTS signature_sets (
                sigset_id INTEGER PRIMARY KEY AUTOINCREMENT,
                atlas_id INTEGER,
                method TEXT NOT NULL,
                topk_per_region INTEGER NOT NULL,
                remove_housekeeping INTEGER DEFAULT 1,
                remove_blood_background INTEGER DEFAULT 1,
                created_at TEXT DEFAULT (datetime('now')),
                params_json TEXT
            )
        """)

        # 11. Signature genes
        cursor.execute("""
            CREATE TABLE IF NOT EXISTS signature_genes (
                sigset_id INTEGER NOT NULL,
                region_id TEXT NOT NULL,
                gene_symbol TEXT NOT NULL,
                weight REAL DEFAULT 1.0,
                PRIMARY KEY(sigset_id, region_id, gene_symbol)
            )
        """)

        # 12. Analysis runs/results (v2)
        cursor.execute("""
            CREATE TABLE IF NOT EXISTS analysis_runs (
                run_id TEXT PRIMARY KEY,
                sample_id TEXT NOT NULL,
                atlas_id INTEGER,
                sigset_id INTEGER,
                method TEXT NOT NULL,
                params_json TEXT,
                created_at TEXT DEFAULT (datetime('now'))
            )
        """)

        cursor.execute("""
            CREATE TABLE IF NOT EXISTS analysis_results (
                run_id TEXT NOT NULL,
                region_id TEXT NOT NULL,
                score REAL NOT NULL,
                fraction REAL,
                ci_low REAL,
                ci_high REAL,
                stability REAL,
                reconstruction_error REAL,
                rank INTEGER NOT NULL,
                PRIMARY KEY(run_id, region_id)
            )
        """)

        # 13. Sample QC
        cursor.execute("""
            CREATE TABLE IF NOT EXISTS sample_qc (
                sample_id TEXT PRIMARY KEY,
                hemolysis_hbb_hba_ratio REAL,
                immune_ptprc REAL,
                albumin_alb REAL,
                brain_signal_score REAL,
                created_at TEXT DEFAULT (datetime('now'))
            )
        """)

        # 创建索引
        cursor.execute("CREATE INDEX IF NOT EXISTS idx_ref_exp_gene ON reference_expression(gene_symbol)")
        cursor.execute("CREATE INDEX IF NOT EXISTS idx_ref_exp_region ON reference_expression(region_id)")
        cursor.execute("CREATE INDEX IF NOT EXISTS idx_cfrna_sample ON cfrna_expression(sample_id)")
        cursor.execute("CREATE INDEX IF NOT EXISTS idx_cfrna_gene ON cfrna_expression(gene_symbol)")
        cursor.execute("CREATE INDEX IF NOT EXISTS idx_region_sig_region ON region_gene_signature(region_id)")
        cursor.execute("CREATE INDEX IF NOT EXISTS idx_region_sig_gene ON region_gene_signature(gene_symbol)")
        cursor.execute("CREATE INDEX IF NOT EXISTS idx_tracing_sample ON source_tracing_results(sample_id)")

        # New indexes
        cursor.execute("CREATE INDEX IF NOT EXISTS idx_sig_genes ON signature_genes(sigset_id, region_id, gene_symbol)")
        cursor.execute("CREATE INDEX IF NOT EXISTS idx_runs_sample ON analysis_runs(sample_id)")
        cursor.execute("CREATE INDEX IF NOT EXISTS idx_results_run_rank ON analysis_results(run_id, rank)")

        self.conn.commit()

        # Backward-compatible schema updates (ALTER TABLE) if needed
        self._ensure_schema_updates()
        try:
            run_migrations(self.db_path)
        except Exception:
            pass

    def _ensure_schema_updates(self):
        """Perform lightweight, backward-compatible ALTER TABLE updates."""
        cursor = self.conn.cursor()

        # Add optional normalized columns to cfrna_expression if missing
        cursor.execute("PRAGMA table_info(cfrna_expression)")
        cols = {row[1] for row in cursor.fetchall()}
        if 'log_tpm' not in cols:
            cursor.execute("ALTER TABLE cfrna_expression ADD COLUMN log_tpm REAL")
        if 'zscore_tpm' not in cols:
            cursor.execute("ALTER TABLE cfrna_expression ADD COLUMN zscore_tpm REAL")

        # Do not auto-seed the removed legacy 8-region rhesus atlas.
        # The current rhesus reference is Bo2023_WangLab_VSD_region.

        self.conn.commit()

    def load_macaque_atlas_data(self, data_file: str = None):
        """加载猕猴全脑转录组图谱数据

        参数:
            data_file: CSV文件路径，包含以下列:
                - region_id: 脑区ID (必需)
                - region_name: 脑区名称 (必需)
                - region_acronym: 脑区缩写
                - parent_region_id: 父脑区ID
                - hemi: 半脑 (L/R/None)
                - layer: 脑层
                - atlas_version: 图谱版本
                - coordinates: 坐标信息
        """
        if data_file is None:
            # 生成演示数据
            self._generate_demo_atlas_data()
        else:
            # 从文件加载
            df = pd.read_csv(data_file)

            # 验证必需列
            required_cols = ['region_id', 'region_name']
            missing_cols = [col for col in required_cols if col not in df.columns]

            if missing_cols:
                raise ValueError(f"数据文件缺少必需列: {missing_cols}")

            # 填充可选列的默认值
            if 'region_acronym' not in df.columns:
                df['region_acronym'] = df['region_id']
            if 'parent_region_id' not in df.columns:
                df['parent_region_id'] = None
            if 'hemi' not in df.columns:
                df['hemi'] = None
            if 'layer' not in df.columns:
                df['layer'] = None
            if 'atlas_version' not in df.columns:
                df['atlas_version'] = 'Custom'
            if 'coordinates' not in df.columns:
                df['coordinates'] = None

            # 选择需要的列并导入
            columns_to_import = ['region_id', 'region_name', 'region_acronym',
                              'parent_region_id', 'hemi', 'layer',
                              'atlas_version', 'coordinates']

            df[columns_to_import].to_sql('macaque_brain_atlas', self.conn,
                                       if_exists='append', index=False)
            print(f"已从 {data_file} 导入 {len(df)} 个脑区数据")

    def _generate_demo_atlas_data(self):
        """生成演示用的脑区数据。

        层级结构:
          BRAIN
          ├── CTX (大脑皮层, hemi=None 代表双侧)
          │   ├── V1 / A1 / M1 / PFC (皮层亚区, hemi=None)
          │   │   └── 各 Layer 亚区
          ├── HIP / THA / STR / BG (皮层下, parent=None)
          ├── CBX (小脑)
          └── BS (脑干)
              ├── MB / PON / MED
        """
        demo_regions = [
            # (region_id, region_name, region_acronym, parent_region_id, hemi, layer, atlas_version, coordinates)
            # — 顶层：大脑皮层（双侧，无 parent）
            ("CTX",     "大脑皮层",             "CTX",  None,   None, None,   "Paxinos2020", None),
            # — 皮层亚区（parent=CTX）
            ("V1",      "初级视觉皮层",          "V1",   "CTX",  None, None,   "Paxinos2020", None),
            ("V1.L1",   "初级视觉皮层-Layer1",   "V1",   "V1",   None, "L1",   "Paxinos2020", None),
            ("V1.L2/3", "初级视觉皮层-Layer2/3", "V1",   "V1",   None, "L2/3", "Paxinos2020", None),
            ("V1.L4",   "初级视觉皮层-Layer4",   "V1",   "V1",   None, "L4",   "Paxinos2020", None),
            ("V1.L5",   "初级视觉皮层-Layer5",   "V1",   "V1",   None, "L5",   "Paxinos2020", None),
            ("V1.L6",   "初级视觉皮层-Layer6",   "V1",   "V1",   None, "L6",   "Paxinos2020", None),
            ("A1",      "初级听觉皮层",          "A1",   "CTX",  None, None,   "Paxinos2020", None),
            ("A1.L2/3", "初级听觉皮层-Layer2/3", "A1",   "A1",   None, "L2/3", "Paxinos2020", None),
            ("A1.L4",   "初级听觉皮层-Layer4",   "A1",   "A1",   None, "L4",   "Paxinos2020", None),
            ("M1",      "初级运动皮层",          "M1",   "CTX",  None, None,   "Paxinos2020", None),
            ("PFC",     "前额叶皮层",            "PFC",  "CTX",  None, None,   "Paxinos2020", None),
            # — 皮层下结构（parent=None，双侧）
            ("HIP",     "海马体",               "HIP",  None,   None, None,   "Paxinos2020", None),
            ("THA",     "丘脑",                 "THA",  None,   None, None,   "Paxinos2020", None),
            ("STR",     "纹状体",               "STR",  None,   None, None,   "Paxinos2020", None),
            ("BG",      "基底节",               "BG",   None,   None, None,   "Paxinos2020", None),
            # — 小脑
            ("CBX",     "小脑",                 "CBX",  None,   None, None,   "Paxinos2020", None),
            # — 脑干（统一归入 BS）
            ("BS",      "脑干",                 "BS",   None,   None, None,   "Paxinos2020", None),
            ("MB",      "中脑",                 "MB",   "BS",   None, None,   "Paxinos2020", None),
            ("PON",     "脑桥",                 "PON",  "BS",   None, None,   "Paxinos2020", None),
            ("MED",     "延髓",                 "MED",  "BS",   None, None,   "Paxinos2020", None),
        ]

        cursor = self.conn.cursor()
        cursor.executemany("""
            INSERT OR IGNORE INTO macaque_brain_atlas
            (region_id, region_name, region_acronym, parent_region_id, hemi, layer, atlas_version, coordinates)
            VALUES (?, ?, ?, ?, ?, ?, ?, ?)
        """, demo_regions)

        self.conn.commit()
        print(f"已加载 {len(demo_regions)} 个脑区数据")

    def load_reference_expression_data(self, expression_file: str = None, n_genes=1000):
        """加载参考表达谱数据

        参数:
            expression_file: CSV文件路径，包含以下列:
                - gene_symbol: 基因符号 (必需)
                - gene_name: 基因名称
                - ensembl_id: Ensembl ID
                - ncbi_id: NCBI ID
                - region_id: 脑区ID (必需)
                - region_name: 脑区名称 (必需)
                - avg_tpm: 平均TPM值 (必需)
                - std_tpm: TPM标准差
                - median_tpm: TPM中位数
                - sample_count: 样本数
                - expression_class: 表达等级 (High/Medium/Low)
                - cell_type_marker: 细胞类型标记
            n_genes: 如果不提供文件，生成模拟数据的基因数
        """
        if expression_file is None:
            # 生成演示数据
            self._generate_reference_expression_data(n_genes)
        else:
            # 从文件加载
            df = pd.read_csv(expression_file)

            # 验证必需列
            required_cols = ['gene_symbol', 'region_id', 'region_name', 'avg_tpm']
            missing_cols = [col for col in required_cols if col not in df.columns]

            if missing_cols:
                raise ValueError(f"数据文件缺少必需列: {missing_cols}")

            # 填充可选列的默认值
            if 'gene_name' not in df.columns:
                df['gene_name'] = df['gene_symbol']
            if 'ensembl_id' not in df.columns:
                df['ensembl_id'] = None
            if 'ncbi_id' not in df.columns:
                df['ncbi_id'] = None
            if 'std_tpm' not in df.columns:
                df['std_tpm'] = df['avg_tpm'] * 0.3
            if 'median_tpm' not in df.columns:
                df['median_tpm'] = df['avg_tpm'] * 0.9
            if 'sample_count' not in df.columns:
                df['sample_count'] = 1
            if 'expression_class' not in df.columns:
                df['expression_class'] = np.where(df['avg_tpm'] > 50, 'High',
                                               np.where(df['avg_tpm'] > 10, 'Medium', 'Low'))
            if 'cell_type_marker' not in df.columns:
                df['cell_type_marker'] = 'Unknown'

            # 选择需要的列并导入
            columns_to_import = ['gene_symbol', 'gene_name', 'ensembl_id', 'ncbi_id',
                              'region_id', 'region_name', 'avg_tpm', 'std_tpm',
                              'median_tpm', 'sample_count', 'expression_class',
                              'cell_type_marker']

            df[columns_to_import].to_sql('reference_expression', self.conn,
                                       if_exists='append', index=False)
            print(f"已从 {expression_file} 导入 {len(df)} 条表达谱数据")

    def _generate_reference_expression_data(self, n_genes=1000):
        """生成参考表达谱数据（模拟数据）"""
        # 常用标记基因
        marker_genes = {
            "V1": ["SLC17A7", "RORB", "CUX2", "FEZF2", "BCL11B"],
            "A1": ["SLC17A7", "GAD1", "SST", "PVALB"],
            "M1": ["FEZF2", "BCL11B", "CTIP2"],
            "PFC": ["SATB2", "CUX1", "CUX2"],
            "HIP": ["PROX1", "CALB1", "RELN"],
            "THA": ["TCF4", "GRIK4"],
            "STR": ["DRD1", "DRD2", "GAD1"],
            "CBX": ["PCP2", "CALB1"],
            "BG": ["GAD1", "GAD2", "SLC32A1"],
        }

        cursor = self.conn.cursor()
        count = 0

        # 为每个脑区生成表达数据
        for region_id in marker_genes.keys():
            for gene_symbol in marker_genes[region_id]:
                # 基础表达值
                base_tpm = np.random.lognormal(4, 1)

                # 如果是标记基因，表达更高
                if gene_symbol in marker_genes[region_id]:
                    base_tpm *= np.random.uniform(2, 10)

                # 添加噪声
                tpm = max(0.1, base_tpm + np.random.normal(0, base_tpm * 0.2))

                cursor.execute("""
                    INSERT OR IGNORE INTO reference_expression
                    (gene_symbol, gene_name, ensembl_id, ncbi_id, region_id, region_name,
                     avg_tpm, std_tpm, median_tpm, sample_count, expression_class, cell_type_marker)
                    VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                """, (
                    gene_symbol,
                    f"Gene {gene_symbol}",
                    f"ENSMG{np.random.randint(100000000, 999999999)}",
                    str(np.random.randint(1000, 99999)),
                    region_id,
                    region_id,
                    tpm,
                    tpm * 0.3,
                    tpm * 0.9,
                    np.random.randint(5, 20),
                    "High" if tpm > 50 else "Medium" if tpm > 10 else "Low",
                    "Marker" if gene_symbol in marker_genes[region_id] else "Non-marker"
                ))
                count += 1

        # 添加更多随机基因
        for _ in range(n_genes - count):
            gene_symbol = f"Gene_{np.random.randint(1000, 9999)}"
            region_id = np.random.choice(list(marker_genes.keys()))
            tpm = max(0.1, np.random.lognormal(2, 1))

            cursor.execute("""
                INSERT OR IGNORE INTO reference_expression
                (gene_symbol, gene_name, ensembl_id, ncbi_id, region_id, region_name,
                 avg_tpm, std_tpm, median_tpm, sample_count, expression_class, cell_type_marker)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            """, (
                gene_symbol,
                f"Gene {gene_symbol}",
                f"ENSMG{np.random.randint(100000000, 999999999)}",
                str(np.random.randint(1000, 99999)),
                region_id,
                region_id,
                tpm,
                tpm * 0.3,
                tpm * 0.9,
                np.random.randint(3, 15),
                "Low",
                "Non-marker"
            ))

        self.conn.commit()
        print(f"已生成 {count + n_genes - count} 条参考表达谱数据")

    def load_region_signatures(self, signature_file: str = None):
        """加载脑区基因签名

        参数:
            signature_file: CSV文件路径，包含以下列:
                - region_id: 脑区ID (必需)
                - region_name: 脑区名称 (必需)
                - gene_symbol: 基因符号 (必需)
                - marker_score: 标记分数
                - specificity_score: 特异性分数
                - expression_level: 表达水平 (High/Medium/Low)
                - is_marker: 是否为标记基因 (1/0)
        """
        if signature_file is None:
            # 生成签名数据
            self._generate_region_signatures()
        else:
            # 从文件加载
            df = pd.read_csv(signature_file)

            # 验证必需列
            required_cols = ['region_id', 'region_name', 'gene_symbol']
            missing_cols = [col for col in required_cols if col not in df.columns]

            if missing_cols:
                raise ValueError(f"数据文件缺少必需列: {missing_cols}")

            # 填充可选列的默认值
            if 'marker_score' not in df.columns:
                df['marker_score'] = 1.0
            if 'specificity_score' not in df.columns:
                df['specificity_score'] = 1.0
            if 'expression_level' not in df.columns:
                df['expression_level'] = 'Medium'
            if 'is_marker' not in df.columns:
                df['is_marker'] = 1

            # 选择需要的列并导入
            columns_to_import = ['region_id', 'region_name', 'gene_symbol',
                              'marker_score', 'specificity_score',
                              'expression_level', 'is_marker']

            df[columns_to_import].to_sql('region_gene_signature', self.conn,
                                       if_exists='append', index=False)
            print(f"已从 {signature_file} 导入 {len(df)} 条基因签名数据")

    def _generate_region_signatures(self):
        """生成脑区基因签名（模拟数据）"""
        cursor = self.conn.cursor()

        # 从参考表达谱中计算签名
        cursor.execute("""
            SELECT region_id, region_name, gene_symbol, avg_tpm
            FROM reference_expression
            ORDER BY region_id, avg_tpm DESC
        """)

        results = cursor.fetchall()

        # 为每个区域选择top基因作为签名
        from collections import defaultdict
        region_genes = defaultdict(list)

        for region_id, region_name, gene_symbol, avg_tpm in results:
            region_genes[region_id].append((gene_symbol, avg_tpm))

        for region_id, genes in region_genes.items():
            # 排序并取前50个
            genes_sorted = sorted(genes, key=lambda x: x[1], reverse=True)[:50]

            for idx, (gene_symbol, avg_tpm) in enumerate(genes_sorted):
                marker_score = 1.0 - (idx / 50)  # 降序评分
                specificity_score = np.random.uniform(0.5, 1.0)

                cursor.execute("""
                    INSERT OR REPLACE INTO region_gene_signature
                    (region_id, region_name, gene_symbol, marker_score, specificity_score,
                     expression_level, is_marker)
                    VALUES (?, ?, ?, ?, ?, ?, ?)
                """, (
                    region_id,
                    region_id,
                    gene_symbol,
                    marker_score,
                    specificity_score,
                    "High" if idx < 10 else "Medium" if idx < 30 else "Low",
                    1 if idx < 20 else 0
                ))

        self.conn.commit()
        print("已生成脑区基因签名")

    def initialize_database(self, atlas_file=None, expression_file=None, signature_file=None):
        """初始化数据库

        参数:
            atlas_file: 脑区图谱文件路径 (可选)
            expression_file: 表达谱数据文件路径 (可选)
            signature_file: 基因签名文件路径 (可选)
            如果不提供这些参数，将生成模拟演示数据
        """
        print("正在初始化数据库...")
        self.connect()
        self.create_database_schema()
        self.load_macaque_atlas_data(atlas_file)
        self.load_reference_expression_data(expression_file, n_genes=500)
        self.load_region_signatures(signature_file)
        print("数据库初始化完成！")

    def close(self):
        """关闭数据库连接"""
        if self.conn:
            self.conn.close()

if __name__ == "__main__":
    # 初始化数据库
    db = CSFRNASourceDatabase()
    db.initialize_database()

    # 查看数据统计
    print("\n数据库统计:")
    cursor = db.conn.cursor()

    tables = ['macaque_brain_atlas', 'reference_expression',
              'cfrna_samples', 'region_gene_signature']

    for table in tables:
        cursor.execute(f"SELECT COUNT(*) FROM {table}")
        count = cursor.fetchone()[0]
        print(f"  {table}: {count} 条记录")

    db.close()
