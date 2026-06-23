#!/usr/bin/env python
from __future__ import annotations

import argparse
import shutil
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
ALLOWED_ROOT = Path("/storage/wangzhen")
DEFAULT_DEST = PROJECT_ROOT / "tests" / "controlled_data"


FILES = {
    "bo2023/mfas5_819samples_28415genes_featurecounts_counts.txt": [
        "bo2023 data/mfas5_819samples_28415genes_featurecounts_counts.txt",
    ],
    "bo2023/mfas5_819samples_23605genes_vsd4_rmbatch.xls": [
        "bo2023 data/mfas5_819samples_23605genes_vsd4_rmbatch.xls",
    ],
    "bo2023/Information of sequenced samples_update_full878_filter819.xlsx": [
        "bo2023 data/Information of sequenced samples_update_full878_filter819.xlsx",
    ],
    "bo2023/04_expressed_genes_neocortex_plus_subcortical.cleaned_symbols.csv": [
        "bo2023_bulk_atlas_buildkit/04_expressed_genes_neocortex_plus_subcortical.cleaned_symbols.csv",
    ],
    "bo2023/01_region_level_meta.csv": [
        "bo2023_bulk_atlas_buildkit/01_region_level_meta.csv",
    ],
    "models/bo2023_reference_projector_linear_full.npz": [
        "data/models/bo2023_reference_projector_linear_full.npz",
        "results/bo2023_reference_projection_20260616_cleaned_symbols/bo2023_reference_projector_linear_full.npz",
        "results/bo2023_reference_projection_20260616/bo2023_reference_projector_linear_full.npz",
    ],
    "ahba/raw_zips/H0351_2001_rnaseq.zip": [
        "data/ahba_human_rnaseq/raw_zips/H0351_2001_rnaseq.zip",
    ],
    "ahba/raw_zips/H0351_2002_rnaseq.zip": [
        "data/ahba_human_rnaseq/raw_zips/H0351_2002_rnaseq.zip",
    ],
    "tcga/tcga_gbm_lgg_primary_tumor_unstranded_counts_sample_sum.tsv": [
        "data/tcga_brain_tumor_expression/tcga_gbm_lgg_primary_tumor_unstranded_counts_sample_sum.tsv",
    ],
    "tcga/tcga_gbm_lgg_gdc_star_counts_manifest.csv": [
        "data/tcga_brain_tumor_expression/tcga_gbm_lgg_gdc_star_counts_manifest.csv",
    ],
    "tcga/corrected_direct_overlap_mri_truth.csv": [
        "results/brats_tcga_lgg_65_mri_truth_corrected_20260612/corrected_direct_overlap_mri_truth.csv",
    ],
    "gse189919/GSE189919_tpm_count.csv.gz": [
        "data/external_validation/GSE189919/GSE189919_tpm_count.csv.gz",
    ],
    "gse189919/GSE189919_count.csv.gz": [
        "data/external_validation/GSE189919/GSE189919_count.csv.gz",
    ],
    "gse189919/GSE189919_family.soft.gz": [
        "data/external_validation/GSE189919/GSE189919_family.soft.gz",
    ],
}


def within_allowed(path: Path, allowed_root: Path) -> bool:
    try:
        path.resolve().relative_to(allowed_root.resolve())
        return True
    except ValueError:
        return False


def candidate_roots(source_root: Path) -> list[Path]:
    roots = [
        PROJECT_ROOT,
        source_root,
        source_root / "cfrna-brain-tracing-0.1.6",
        source_root / "github_projects" / "cfrna-brain-tracing-0.1.6",
        source_root / "github_projects" / "cfrna-brain-tracing-streamlit-cloud-ready",
    ]
    seen: set[Path] = set()
    unique: list[Path] = []
    for root in roots:
        resolved = root.resolve()
        if resolved in seen:
            continue
        seen.add(resolved)
        unique.append(resolved)
    return unique


def find_source(relative_candidates: list[str], roots: list[Path], allowed_root: Path) -> Path | None:
    for root in roots:
        if not within_allowed(root, allowed_root):
            continue
        for rel in relative_candidates:
            path = root / rel
            if path.exists() and within_allowed(path, allowed_root):
                return path
    return None


def copy_file(src: Path, dest: Path, dry_run: bool = False) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dry_run:
        print(f"DRY-RUN copy: {src} -> {dest}")
        return
    shutil.copy2(src, dest)
    print(f"copied: {src} -> {dest}")


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Stage controlled validation data under tests/controlled_data using only /storage/wangzhen paths."
    )
    parser.add_argument("--source-root", type=Path, default=ALLOWED_ROOT)
    parser.add_argument("--dest", type=Path, default=DEFAULT_DEST)
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    allowed_root = ALLOWED_ROOT.resolve()
    source_root = args.source_root.resolve()
    dest = args.dest.resolve()
    if not within_allowed(source_root, allowed_root):
        raise SystemExit(f"Refusing source outside {allowed_root}: {source_root}")
    if not within_allowed(dest, allowed_root):
        raise SystemExit(f"Refusing destination outside {allowed_root}: {dest}")

    roots = candidate_roots(source_root)
    missing: list[str] = []
    for dest_rel, rel_candidates in FILES.items():
        source = find_source(rel_candidates, roots, allowed_root)
        if source is None:
            missing.append(dest_rel)
            print(f"missing source for: {dest_rel}")
            continue
        copy_file(source, dest / dest_rel, dry_run=args.dry_run)

    if missing:
        print("\nMissing files:")
        for rel in missing:
            print(" -", rel)
        return 2

    print("\nAll controlled validation data staged under:", dest)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
