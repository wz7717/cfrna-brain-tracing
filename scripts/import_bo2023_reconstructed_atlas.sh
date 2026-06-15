#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 3 ]]; then
  echo "用法: bash scripts/import_bo2023_reconstructed_atlas.sh 数据库路径 gene_by_region_tpm.tsv.gz annotation.tsv"
  exit 1
fi

DB_PATH="$1"
MATRIX_PATH="$2"
ANNOTATION_PATH="$3"

python cli.py import-bo2023-region-matrix \
  --db "$DB_PATH" \
  --matrix "$MATRIX_PATH" \
  --annotation "$ANNOTATION_PATH" \
  --atlas-name "WangLab Bo2023 reconstructed bulk atlas" \
  --build-version "reconstructed_from_PRJNA905082"
