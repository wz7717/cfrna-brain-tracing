from __future__ import annotations

import gzip
import json
import sqlite3
from pathlib import Path
from typing import Optional

import pandas as pd


REQUIRED_REGION_COL = 'brain_region'


def _read_table(path: str | Path) -> pd.DataFrame:
    p = Path(path)
    if p.suffix.lower() in {'.gz', '.bgz'}:
        with gzip.open(p, 'rt', encoding='utf-8') as fh:
            return pd.read_csv(fh, sep='\t')
    if p.suffix.lower() in {'.xlsx', '.xls'}:
        return pd.read_excel(p)
    sep = '\t' if p.suffix.lower() in {'.tsv', '.txt'} else ','
    return pd.read_csv(p, sep=sep)


def _next_atlas_id(conn: sqlite3.Connection) -> int:
    cur = conn.cursor()
    cur.execute('SELECT COALESCE(MAX(atlas_id), 0) + 1 FROM atlas_versions')
    return int(cur.fetchone()[0])


def _normalize_annotation(annotation_df: pd.DataFrame) -> pd.DataFrame:
    df = annotation_df.copy()
    if REQUIRED_REGION_COL not in df.columns:
        raise ValueError(f'annotation file 缺少必需列: {REQUIRED_REGION_COL}')
    rename_map = {
        'region_full_name': 'region_name',
        'brain_region': 'region_id',
        'regional_map': 'parent_region_id',
        'meta_lobe': 'lobe',
    }
    for old, new in rename_map.items():
        if old in df.columns and new not in df.columns:
            df = df.rename(columns={old: new})

    keep_cols = [c for c in [
        'region_id', 'region_name', 'lobe', 'parent_region_id',
        'neocortex_flag', 'roi173'
    ] if c in df.columns]
    df = df[keep_cols].drop_duplicates(subset=['region_id'])
    if 'region_name' not in df.columns:
        df['region_name'] = df['region_id']
    if 'lobe' not in df.columns:
        df['lobe'] = ''
    if 'parent_region_id' not in df.columns:
        df['parent_region_id'] = ''
    if 'neocortex_flag' not in df.columns:
        df['neocortex_flag'] = ''
    if 'roi173' not in df.columns:
        df['roi173'] = ''
    return df


def import_region_matrix(
    db_path: str,
    matrix_path: str,
    annotation_path: str,
    atlas_name: str = 'WangLab Bo2023 reconstructed bulk atlas',
    build_version: str = 'reconstructed_from_PRJNA905082',
    gene_id_type: str = 'gene_symbol',
    normalization: str = 'TPM',
    notes: Optional[str] = None,
) -> dict:
    """导入 gene×region TPM 矩阵到现有 5.2 SQLite 系统。

    matrix_path 预期格式：TSV/TSV.GZ，至少包含前两列 gene_id/gene_name，后续列为脑区。
    annotation_path 使用 sample_annotation_master_auto_brain_region.tsv 即可。
    """
    matrix = _read_table(matrix_path)
    ann_raw = _read_table(annotation_path)
    ann = _normalize_annotation(ann_raw)

    if matrix.shape[1] < 3:
        raise ValueError('matrix 文件列数过少，至少应包含 gene_id/gene_name + 一个脑区列。')

    # 兼容不同首列命名
    if 'gene_symbol' not in matrix.columns:
        if 'gene_name' in matrix.columns:
            matrix = matrix.rename(columns={'gene_name': 'gene_symbol'})
        elif 'gene_id' in matrix.columns:
            matrix['gene_symbol'] = matrix['gene_id'].astype(str)
        else:
            matrix = matrix.rename(columns={matrix.columns[0]: 'gene_symbol'})
    if 'gene_name' not in matrix.columns:
        matrix['gene_name'] = matrix['gene_symbol'].astype(str)

    fixed_cols = {'gene_id', 'gene_symbol', 'gene_name'}
    region_cols = [c for c in matrix.columns if c not in fixed_cols]
    region_cols = [c for c in region_cols if pd.notna(c) and str(c).strip()]
    if not region_cols:
        raise ValueError('没有识别到脑区列。')

    missing_regions = sorted(set(region_cols) - set(ann['region_id'].astype(str)))

    conn = sqlite3.connect(db_path)
    try:
        atlas_id = _next_atlas_id(conn)
        cur = conn.cursor()
        cur.execute(
            """
            INSERT INTO atlas_versions
            (atlas_id, atlas_name, species, level, build_version, gene_id_type, normalization, notes)
            VALUES (?, ?, ?, ?, ?, ?, ?, ?)
            """,
            (
                atlas_id,
                atlas_name,
                'Macaca mulatta',
                'region',
                build_version,
                gene_id_type,
                normalization,
                notes or json.dumps({
                    'matrix_path': str(Path(matrix_path).name),
                    'annotation_path': str(Path(annotation_path).name),
                    'missing_regions_without_annotation': missing_regions,
                }, ensure_ascii=False),
            ),
        )

        atlas_version_text = f'{atlas_name} | {build_version}'
        ann_use = ann[ann['region_id'].astype(str).isin(region_cols)].copy()
        atlas_rows = []
        for _, r in ann_use.iterrows():
            coords = json.dumps({
                'lobe': r.get('lobe', ''),
                'neocortex_flag': r.get('neocortex_flag', ''),
                'roi173': r.get('roi173', ''),
            }, ensure_ascii=False)
            atlas_rows.append((
                str(r['region_id']),
                str(r.get('region_name', r['region_id'])),
                str(r['region_id']),
                str(r.get('parent_region_id', '') or ''),
                None,
                None,
                atlas_version_text,
                coords,
                atlas_id,
            ))
        cur.executemany(
            """
            INSERT INTO macaque_brain_atlas
            (region_id, region_name, region_acronym, parent_region_id, hemi, layer, atlas_version, coordinates, atlas_id)
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
            """,
            atlas_rows,
        )

        # 长表写入参考表达矩阵
        region_name_map = {
            str(r['region_id']): str(r.get('region_name', r['region_id']))
            for _, r in ann_use.iterrows()
        }
        melted = matrix[['gene_symbol', 'gene_name'] + region_cols].melt(
            id_vars=['gene_symbol', 'gene_name'],
            value_vars=region_cols,
            var_name='region_id',
            value_name='avg_tpm',
        )
        melted = melted.dropna(subset=['gene_symbol', 'region_id', 'avg_tpm'])
        melted['region_id'] = melted['region_id'].astype(str)
        melted = melted[melted['region_id'].isin(region_name_map)]
        melted['region_name'] = melted['region_id'].map(region_name_map)
        melted['avg_tpm'] = pd.to_numeric(melted['avg_tpm'], errors='coerce').fillna(0.0)
        melted['atlas_id'] = atlas_id
        melted['expression_class'] = pd.cut(
            melted['avg_tpm'],
            bins=[-1, 0.1, 1, 10, float('inf')],
            labels=['silent', 'low', 'medium', 'high'],
        ).astype(str)
        melted['sample_count'] = None
        melted['std_tpm'] = None
        melted['median_tpm'] = melted['avg_tpm']
        melted['cell_type_marker'] = None
        rows = [
            (
                str(r.gene_symbol),
                str(r.gene_name),
                None,
                None,
                str(r.region_id),
                str(r.region_name),
                float(r.avg_tpm),
                None,
                float(r.median_tpm),
                None,
                str(r.expression_class),
                None,
                int(r.atlas_id),
            )
            for r in melted.itertuples(index=False)
        ]
        cur.executemany(
            """
            INSERT INTO reference_expression
            (gene_symbol, gene_name, ensembl_id, ncbi_id, region_id, region_name, avg_tpm, std_tpm, median_tpm, sample_count, expression_class, cell_type_marker, atlas_id)
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            """,
            rows,
        )
        conn.commit()
        return {
            'atlas_id': atlas_id,
            'regions_imported': int(len(atlas_rows)),
            'genes_imported': int(matrix['gene_symbol'].astype(str).nunique()),
            'expression_rows_imported': int(len(rows)),
            'missing_regions_without_annotation': missing_regions,
        }
    finally:
        conn.close()
