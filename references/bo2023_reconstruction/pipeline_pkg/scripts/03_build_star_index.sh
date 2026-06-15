#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 4 ]]; then
  echo "用法: bash scripts/03_build_star_index.sh genome.fa genes.gtf star_index_dir threads"
  exit 1
fi

GENOME_FA="$1"
GTF="$2"
INDEX_DIR="$3"
THREADS="$4"

mkdir -p "$INDEX_DIR"

STAR \
  --runThreadN "$THREADS" \
  --runMode genomeGenerate \
  --genomeDir "$INDEX_DIR" \
  --genomeFastaFiles "$GENOME_FA" \
  --sjdbGTFfile "$GTF" \
  --sjdbOverhang 149
