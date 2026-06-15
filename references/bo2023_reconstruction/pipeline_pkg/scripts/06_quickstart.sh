#!/usr/bin/env bash
set -euo pipefail

cat <<'EOF'
这个 quickstart 只是提醒执行顺序，不会自动下载参考基因组。

1) 准备 RunInfo.csv
   python scripts/01_prepare_runinfo.py --runinfo RunInfo.csv --outdir meta

2) 手工补齐 meta/sample_annotation_master.tsv 中的 brain_region / monkey_id

3) 下载并转 FASTQ
   bash scripts/02_prefetch_and_fasterq_batch.sh meta/srr_list.txt work/sra work/fastq 8

4) 建 STAR 索引
   bash scripts/03_build_star_index.sh ref/genome.fa ref/genes.gtf work/star_index 16

5) 比对 + featureCounts
   bash scripts/04_align_and_count_batch.sh meta/sample_annotation_master.tsv ref/genome.fa ref/genes.gtf work/star_index work/fastq work/bam work/counts 16

6) 生成最终矩阵
   python scripts/05_merge_featurecounts.py --counts-dir work/counts --sample-sheet meta/sample_annotation_master.tsv --gtf ref/genes.gtf --outdir results
EOF
