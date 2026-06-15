#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 4 ]]; then
  echo "用法: bash scripts/02_prefetch_and_fasterq_batch.sh srr_list.txt sra_dir fastq_dir threads"
  exit 1
fi

SRR_LIST="$1"
SRA_DIR="$2"
FASTQ_DIR="$3"
THREADS="$4"

mkdir -p "$SRA_DIR" "$FASTQ_DIR"

while read -r SRR; do
  [[ -z "$SRR" ]] && continue
  echo "==== 下载 $SRR ===="
  prefetch --max-size 200G --output-directory "$SRA_DIR" "$SRR"

  echo "==== 转 FASTQ $SRR ===="
  # fasterq-dump 默认会生成 _1/_2.fastq
  fasterq-dump \
    --threads "$THREADS" \
    --outdir "$FASTQ_DIR" \
    --temp "$FASTQ_DIR/tmp_${SRR}" \
    "$SRA_DIR/$SRR"

  gzip -f "$FASTQ_DIR/${SRR}"*.fastq || true
done < "$SRR_LIST"
