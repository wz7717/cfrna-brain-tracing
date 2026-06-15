# Bo2023 重建表达矩阵接入说明

## 现在这版项目已经接入了什么
1. `references/bo2023_reconstruction/sample_annotation_master_auto_brain_region.tsv`
   - 已根据 `SraRunInfo.csv` 自动补好 `run_id -> brain_region` 映射。
2. `references/bo2023_reconstruction/pipeline_pkg/`
   - 包含从 PRJNA905082 原始数据构建 counts / TPM / gene×region 矩阵的中文脚本包。
3. `cli.py`
   - 新增命令：

```bash
python cli.py import-bo2023-region-matrix \
  --db cfrna_source_tracing.db \
  --matrix results/bo2023_gene_by_region_mean_tpm.tsv.gz \
  --annotation references/bo2023_reconstruction/sample_annotation_master_auto_brain_region.tsv
```

4. `scripts/import_bo2023_reconstructed_atlas.sh`
   - 上述命令的简化封装。

## 这一步为什么还不能直接替你产出最终数值矩阵
因为上传给我的文件里没有真正的 `featureCounts` 结果，也没有 PRJNA905082 的 FASTQ/BAM 原始数据。
所以当前环境里无法生成真实的 `bo2023_gene_by_region_mean_tpm.tsv.gz` 数值文件。

## 你拿到原始数据后，最短路径
### 1. 进入 pipeline 包
```bash
cd references/bo2023_reconstruction/pipeline_pkg/prjna905082_bulk_matrix_pipeline_cn
```

### 2. 用 RunInfo 生成下载列表
```bash
python scripts/01_prepare_runinfo.py --runinfo /你的路径/SraRunInfo.csv --outdir meta
```

### 3. 用已自动补好的 annotation 替换模板
```bash
cp ../../sample_annotation_master_auto_brain_region.tsv meta/sample_annotation_master.tsv
```

### 4. 跑下载、比对、计数、合并
最后会生成：
- `bo2023_gene_by_sample_counts.tsv.gz`
- `bo2023_gene_by_sample_tpm.tsv.gz`
- `bo2023_gene_by_region_mean_tpm.tsv.gz`

### 5. 导入 5.2 溯源系统
```bash
bash scripts/import_bo2023_reconstructed_atlas.sh \
  cfrna_source_tracing.db \
  /你的results/bo2023_gene_by_region_mean_tpm.tsv.gz \
  references/bo2023_reconstruction/sample_annotation_master_auto_brain_region.tsv
```

导入后前端“溯源分析”里的 Atlas 下拉框会自动出现新的 Bo2023 atlas。

## 建议
导入后再运行一次 signature 构建，这样新 atlas 的区分度会明显更好。
