#!/usr/bin/env python
from __future__ import annotations

import argparse
import csv
import gzip
import re
from collections import defaultdict
from pathlib import Path


GENE_ID_RE = re.compile(r"ENSMFAG\d+")


def open_text(path: Path):
    if path.suffix.lower() == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("r", encoding="utf-8", errors="replace")


def parse_attributes(attr: str) -> dict[str, str]:
    attr = attr.strip()
    out: dict[str, str] = {}
    if not attr:
        return out

    # GTF style: key "value"; key2 "value2";
    if ";" in attr and '"' in attr:
        for chunk in attr.split(";"):
            chunk = chunk.strip()
            if not chunk:
                continue
            parts = chunk.split(None, 1)
            if len(parts) != 2:
                continue
            key, value = parts
            out[key] = value.strip().strip('"')
        return out

    # GFF3 style: key=value;key2=value2
    for chunk in attr.split(";"):
        chunk = chunk.strip()
        if not chunk or "=" not in chunk:
            continue
        key, value = chunk.split("=", 1)
        out[key.strip()] = value.strip()
    return out


def split_ids(value: str | None) -> list[str]:
    if not value:
        return []
    ids: list[str] = []
    for item in str(value).replace(",", ";").split(";"):
        item = item.strip()
        if not item:
            continue
        # Common GFF3 values can look like transcript:ENSMFAT... or gene:ENSMFAG...
        item = item.split(":", 1)[-1]
        ids.append(item)
    return ids


def first_ensmfag(*values: str | None) -> str | None:
    for value in values:
        if not value:
            continue
        match = GENE_ID_RE.search(str(value))
        if match:
            return match.group(0)
    return None


def merge_intervals(intervals: list[tuple[int, int]]) -> tuple[int, int]:
    if not intervals:
        return 0, 0
    intervals = sorted(intervals)
    merged: list[list[int]] = []
    for start, end in intervals:
        if start > end:
            start, end = end, start
        if not merged or start > merged[-1][1] + 1:
            merged.append([start, end])
        else:
            merged[-1][1] = max(merged[-1][1], end)
    length = sum(end - start + 1 for start, end in merged)
    return length, len(merged)


def collect_transcript_to_gene(path: Path) -> dict[str, str]:
    transcript_to_gene: dict[str, str] = {}
    transcript_features = {
        "mrna",
        "transcript",
        "ncrna",
        "lncrna",
        "mirna",
        "rrna",
        "trna",
        "snorna",
        "snrna",
        "pseudogenic_transcript",
    }
    with open_text(path) as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            feature = parts[2].lower()
            if feature not in transcript_features:
                continue
            attrs = parse_attributes(parts[8])
            transcript_ids = split_ids(attrs.get("ID") or attrs.get("transcript_id"))
            gene_id = first_ensmfag(attrs.get("gene_id"), attrs.get("Parent"), attrs.get("gene"), attrs.get("gene_name"))
            if not gene_id:
                continue
            for transcript_id in transcript_ids:
                transcript_to_gene[transcript_id] = gene_id
    return transcript_to_gene


def exon_gene_id(attrs: dict[str, str], transcript_to_gene: dict[str, str]) -> str | None:
    gene_id = first_ensmfag(
        attrs.get("gene_id"),
        attrs.get("gene"),
        attrs.get("gene_name"),
        attrs.get("locus_tag"),
        attrs.get("ID"),
    )
    if gene_id:
        return gene_id
    for parent_id in split_ids(attrs.get("Parent") or attrs.get("transcript_id")):
        if parent_id in transcript_to_gene:
            return transcript_to_gene[parent_id]
        gene_id = first_ensmfag(parent_id)
        if gene_id:
            return gene_id
    return None


def compute_gene_lengths(path: Path, feature: str) -> tuple[list[dict[str, object]], dict[str, int]]:
    transcript_to_gene = collect_transcript_to_gene(path)
    intervals_by_gene: dict[str, list[tuple[str, int, int]]] = defaultdict(list)
    raw_exon_counts: dict[str, int] = defaultdict(int)
    feature = feature.lower()
    stats = {
        "parsed_feature_rows": 0,
        "feature_rows_without_gene_id": 0,
        "transcript_parent_links": len(transcript_to_gene),
    }

    with open_text(path) as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            if parts[2].lower() != feature:
                continue
            stats["parsed_feature_rows"] += 1
            try:
                start = int(parts[3])
                end = int(parts[4])
            except ValueError:
                continue
            attrs = parse_attributes(parts[8])
            gene_id = exon_gene_id(attrs, transcript_to_gene)
            if not gene_id:
                stats["feature_rows_without_gene_id"] += 1
                continue
            seqid = parts[0]
            intervals_by_gene[gene_id].append((seqid, start, end))
            raw_exon_counts[gene_id] += 1

    rows: list[dict[str, object]] = []
    for gene_id, seq_intervals in sorted(intervals_by_gene.items()):
        by_seqid: dict[str, list[tuple[int, int]]] = defaultdict(list)
        for seqid, start, end in seq_intervals:
            by_seqid[seqid].append((start, end))
        length = 0
        merged_count = 0
        for intervals in by_seqid.values():
            seq_length, seq_merged_count = merge_intervals(intervals)
            length += seq_length
            merged_count += seq_merged_count
        rows.append(
            {
                "Geneid": gene_id,
                "Length": int(length),
                "n_exon_intervals": int(raw_exon_counts[gene_id]),
                "n_merged_exon_intervals": int(merged_count),
                "seqids": ",".join(sorted(by_seqid)),
            }
        )

    stats["n_genes_with_length"] = len(rows)
    return rows, stats


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Compute gene-level exon union length from a mfas5 GTF/GFF/GFF3 annotation. "
            "The output can be passed to build_bo2023_tpm_atlas_from_counts.py via --gene-lengths."
        )
    )
    parser.add_argument("--annotation", type=Path, required=True, help="mfas5 GTF/GFF/GFF3 annotation, optionally .gz")
    parser.add_argument("--output", type=Path, default=Path("mfas5_gene_exon_union_lengths.csv"))
    parser.add_argument("--feature", default="exon", help="Feature type used for length calculation. Default: exon")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if not args.annotation.exists():
        raise SystemExit(f"ERROR: annotation file does not exist: {args.annotation}")

    rows, stats = compute_gene_lengths(args.annotation, args.feature)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["Geneid", "Length", "n_exon_intervals", "n_merged_exon_intervals", "seqids"],
        )
        writer.writeheader()
        writer.writerows(rows)

    print(f"Wrote {len(rows)} gene lengths to {args.output}")
    print("Summary:")
    for key, value in stats.items():
        print(f"- {key}: {value}")
    if not rows:
        raise SystemExit(
            "ERROR: no gene lengths were produced. Check whether the annotation uses exon features and ENSMFAG gene IDs."
        )


if __name__ == "__main__":
    main()
