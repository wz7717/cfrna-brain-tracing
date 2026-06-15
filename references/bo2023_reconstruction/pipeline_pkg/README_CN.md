# PRJNA905082 批量下载 + 基因表达矩阵构建包（中文）

这套脚本是按 **PRJNA905082 原始数据 + 王征组公开代码思路**整理的实用版重建流程。
目标输出三张核心矩阵：

- `bo2023_gene_by_sample_counts.tsv.gz`
- `bo2023_gene_by_sample_tpm.tsv.gz`
- `bo2023_gene_by_region_mean_tpm.tsv.gz`

## 你现在需要自己准备的东西
1. `RunInfo.csv`
   - 推荐直接从 NCBI 网页导出
   - 或用链接：`https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/runinfo?acc=PRJNA905082`
2. 恒河猴/食蟹猴参考基因组 FASTA
3. 对应 GTF 注释文件
4. 足够的硬盘空间（非常重要）

## 推荐目录结构
```bash
project/
├── RunInfo.csv
├── ref/
│   ├── genome.fa
│   └── genes.gtf
├── work/
└── results/
```

## 最短使用流程
### 第 1 步：把 RunInfo 变成下载列表和样本注释模板
```bash
python scripts/01_prepare_runinfo.py \
  --runinfo RunInfo.csv \
  --outdir meta
```

输出：
- `meta/srr_list.txt`
- `meta/sample_annotation_master.tsv`

### 第 2 步：在 `meta/sample_annotation_master.tsv` 里补齐脑区信息
至少补这些列：
- `brain_region`
- `region_group`
- `monkey_id`
- `include`

如果你暂时没有完整脑区映射，也可以先只填 `brain_region = NA`，先生成 `gene × sample` 矩阵。

### 第 3 步：批量下载并转 FASTQ
```bash
bash scripts/02_prefetch_and_fasterq_batch.sh \
  meta/srr_list.txt \
  work/sra \
  work/fastq \
  8
```

### 第 4 步：建立 STAR 索引
```bash
bash scripts/03_build_star_index.sh \
  ref/genome.fa \
  ref/genes.gtf \
  work/star_index \
  16
```

### 第 5 步：批量比对 + featureCounts
```bash
bash scripts/04_align_and_count_batch.sh \
  meta/sample_annotation_master.tsv \
  ref/genome.fa \
  ref/genes.gtf \
  work/star_index \
  work/fastq \
  work/bam \
  work/counts \
  16
```

### 第 6 步：合并 counts 并计算 TPM / region matrix
```bash
python scripts/05_merge_featurecounts.py \
  --counts-dir work/counts \
  --sample-sheet meta/sample_annotation_master.tsv \
  --gtf ref/genes.gtf \
  --outdir results
```

如果 `sample_annotation_master.tsv` 里已经填了 `brain_region`，会额外输出：
- `results/bo2023_gene_by_region_mean_tpm.tsv.gz`

## 为什么这套流程比较接近论文/仓库
- 仓库里的 R 脚本明确依赖 featureCounts 产生的输入对象
- 仓库后续分析用的是已经整理好的 sample / region 级矩阵
- 所以这里优先走 `prefetch → fasterq-dump → STAR → featureCounts → TPM/region aggregation`

## 注意
1. 这套脚本能稳定生成你自己的重建版表达矩阵。
2. 但要和论文最终 819 samples / 110 brain regions 完全一致，还需要你把样本筛选和脑区映射补齐。
3. 如果中途下载断了，`prefetch` 支持续传。
