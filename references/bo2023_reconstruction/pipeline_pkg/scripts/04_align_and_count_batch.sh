#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 8 ]]; then
  echo "用法: bash scripts/04_align_and_count_batch.sh sample_annotation_master.tsv genome.fa genes.gtf star_index fastq_dir bam_dir counts_dir threads"
  exit 1
fi

SAMPLE_SHEET="$1"
GENOME_FA="$2"
GTF="$3"
STAR_INDEX="$4"
FASTQ_DIR="$5"
BAM_DIR="$6"
COUNTS_DIR="$7"
THREADS="$8"

mkdir -p "$BAM_DIR" "$COUNTS_DIR"

tail -n +2 "$SAMPLE_SHEET" | while IFS=$'\t' read -r RUN_ID SAMPLE_NAME LIBRARY_NAME SCI_NAME SPOTS BASES PLATFORM MONKEY_ID BRAIN_REGION REGION_GROUP INCLUDE NOTES; do
  [[ "$INCLUDE" == "0" ]] && continue

  R1="$FASTQ_DIR/${RUN_ID}_1.fastq.gz"
  R2="$FASTQ_DIR/${RUN_ID}_2.fastq.gz"

  if [[ ! -f "$R1" || ! -f "$R2" ]]; then
    echo "缺 FASTQ: $RUN_ID"
    continue
  fi

  OUT_PREFIX="$BAM_DIR/${RUN_ID}."
  BAM_OUT="$BAM_DIR/${RUN_ID}.Aligned.sortedByCoord.out.bam"

  if [[ ! -f "$BAM_OUT" ]]; then
    STAR \
      --runThreadN "$THREADS" \
      --genomeDir "$STAR_INDEX" \
      --readFilesIn "$R1" "$R2" \
      --readFilesCommand zcat \
      --outSAMtype BAM SortedByCoordinate \
      --outFileNamePrefix "$OUT_PREFIX"
  fi

  COUNTS_OUT="$COUNTS_DIR/${RUN_ID}.featureCounts.txt"
  if [[ ! -f "$COUNTS_OUT" ]]; then
    featureCounts \
      -T "$THREADS" \
      -p \
      -B \
      -C \
      -a "$GTF" \
      -o "$COUNTS_OUT" \
      "$BAM_OUT"
  fi
done
