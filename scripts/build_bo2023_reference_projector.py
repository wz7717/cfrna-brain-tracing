#!/usr/bin/env python
from __future__ import annotations

import argparse
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from core.reference_projection import (  # noqa: E402
    align_matrices,
    apply_projector,
    compute_logcpm,
    fit_linear_projector,
    map_index_to_symbols,
    missing_items,
    read_bo2023_gene_matrix,
    read_gene_map,
    save_projector_npz,
    summarize_fit,
    write_json,
)


DEFAULT_COUNTS = ROOT / "bo2023 data" / "mfas5_819samples_28415genes_featurecounts_counts.txt"
DEFAULT_VSD = ROOT / "bo2023 data" / "mfas5_819samples_23605genes_vsd4_rmbatch.xls"
DEFAULT_SAMPLE_INFO = ROOT / "bo2023 data" / "Information of sequenced samples_update_full878_filter819.xlsx"
DEFAULT_GENE_MAP = ROOT / "bo2023_bulk_atlas_buildkit" / "04_expressed_genes_neocortex_plus_subcortical.csv"
DEFAULT_MODEL_GENES = ROOT / "data" / "models" / "bo2023_saleem_network_top200_model_genes.csv"


def default_outdir() -> Path:
    return ROOT / "results" / f"bo2023_reference_projection_{datetime.now().strftime('%Y%m%d')}"


def read_sample_info(path: Path, sheet: str) -> pd.DataFrame:
    info = pd.read_excel(path, sheet_name=sheet)
    if "No." not in info.columns:
        raise ValueError("sample metadata must contain a 'No.' column")
    info = info.copy()
    info["sample_id"] = info["No."].astype(str).str.strip()
    return info.drop_duplicates("sample_id")


def read_locked_model_genes(path: Path) -> list[str]:
    if not path.exists():
        return []
    df = pd.read_csv(path)
    if "gene_symbol" not in df.columns:
        return []
    return sorted(df["gene_symbol"].dropna().astype(str).str.strip().unique().tolist())


def main() -> int:
    parser = argparse.ArgumentParser(description="Audit Bo2023 data and train a full linear reference projector.")
    parser.add_argument("--counts", type=Path, default=DEFAULT_COUNTS)
    parser.add_argument("--vsd", type=Path, default=DEFAULT_VSD)
    parser.add_argument("--sample-info", type=Path, default=DEFAULT_SAMPLE_INFO)
    parser.add_argument("--sample-sheet", default="mfas5_819samples_phenSet4")
    parser.add_argument("--gene-map", type=Path, default=DEFAULT_GENE_MAP)
    parser.add_argument("--locked-model-genes", type=Path, default=DEFAULT_MODEL_GENES)
    parser.add_argument("--outdir", type=Path, default=default_outdir())
    parser.add_argument("--audit-only", action="store_true")
    parser.add_argument("--save-projected-training-matrix", action="store_true")
    parser.add_argument("--min-nonzero-samples", type=int, default=10)
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    counts_raw = read_bo2023_gene_matrix(args.counts, dtype="float32")
    vsd_raw = read_bo2023_gene_matrix(args.vsd, dtype="float32")
    sample_info = read_sample_info(args.sample_info, args.sample_sheet)
    gene_map = read_gene_map(args.gene_map)

    counts, count_map_audit = map_index_to_symbols(counts_raw, gene_map)
    vsd, vsd_map_audit = map_index_to_symbols(vsd_raw, gene_map)
    locked_genes = read_locked_model_genes(args.locked_model_genes)
    _, _, common_genes, common_samples = align_matrices(counts, vsd)
    metadata_samples = sample_info["sample_id"].astype(str).tolist()

    common_panel = pd.DataFrame({"gene_symbol": common_genes})
    common_panel["in_locked_network_model"] = common_panel["gene_symbol"].isin(set(locked_genes))
    common_panel.to_csv(args.outdir / "common_gene_panel.csv", index=False)
    count_map_audit.to_csv(args.outdir / "count_gene_symbol_mapping_audit.csv", index=False)
    vsd_map_audit.to_csv(args.outdir / "vsd_gene_symbol_mapping_audit.csv", index=False)

    audit = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "counts_path": str(args.counts),
        "vsd_path": str(args.vsd),
        "sample_info_path": str(args.sample_info),
        "gene_map_path": str(args.gene_map),
        "n_count_genes_raw": int(counts_raw.shape[0]),
        "n_vsd_genes_raw": int(vsd_raw.shape[0]),
        "n_count_gene_symbols": int(counts.shape[0]),
        "n_vsd_gene_symbols": int(vsd.shape[0]),
        "n_count_samples": int(counts.shape[1]),
        "n_vsd_samples": int(vsd.shape[1]),
        "n_metadata_samples": int(sample_info["sample_id"].nunique()),
        "n_common_samples": int(len(common_samples)),
        "n_common_genes": int(len(common_genes)),
        "missing_count_samples_from_vsd": missing_items(counts.columns, vsd.columns),
        "missing_vsd_samples_from_counts": missing_items(vsd.columns, counts.columns),
        "missing_common_samples_from_metadata": missing_items(common_samples, metadata_samples),
        "n_locked_model_genes": int(len(locked_genes)),
        "n_locked_model_genes_in_common_panel": int(len(set(locked_genes) & set(common_genes))),
        "missing_locked_model_genes": missing_items(locked_genes, common_genes),
        "existing_vsd_reconstruction_artifacts": {
            "frozen_vst_reference": str(ROOT / "results" / "bo2023_vsd_reconstruction" / "bo2023_frozen_vst_reference.rds"),
            "best_reconstructed_vsd": str(ROOT / "results" / "bo2023_vsd_reconstruction" / "best_reconstructed_vsd.tsv.gz"),
        },
    }
    write_json(args.outdir / "data_audit_summary.json", audit)

    if args.audit_only:
        print(f"Wrote audit outputs to {args.outdir}")
        return 0

    counts_aligned, vsd_aligned, _, _ = align_matrices(counts, vsd)
    logcpm = compute_logcpm(counts_aligned)
    fit, params = fit_linear_projector(
        logcpm,
        vsd_aligned,
        min_nonzero_samples=args.min_nonzero_samples,
    )
    projected = apply_projector(fit, logcpm)
    fit_summary = summarize_fit(projected, vsd_aligned, params)
    params.to_csv(args.outdir / "projector_gene_parameters.csv", index=False)
    save_projector_npz(
        args.outdir / "bo2023_reference_projector_linear_full.npz",
        fit,
        {
            "created_at_utc": datetime.now(timezone.utc).isoformat(),
            "training_samples": common_samples,
            "training_gene_panel": "Bo2023 count/VSD common gene-symbol panel",
            "input_space": "Bo2023 raw count-derived logCPM",
            "target_space": "Bo2023 VSD batch-removed author matrix",
            "clip_quantiles": [0.005, 0.995],
            "min_nonzero_samples": args.min_nonzero_samples,
        },
    )
    qc = {**audit, "training_fit": fit_summary}
    write_json(args.outdir / "projector_qc_summary.json", qc)

    if args.save_projected_training_matrix:
        projected.to_csv(args.outdir / "bo2023_projected_vsd_training_full.tsv.gz", sep="\t", compression="gzip")

    note = f"""# Reference Projection Method Note

This exploratory branch fits a per-gene empirical projector from paired Bo2023 raw count-derived logCPM to the released Bo2023 VSD batch-removed matrix.

- Input transform: `log1p(count / sample_library_size * 1,000,000)`.
- Projector: one ordinary least-squares linear model per gene, `VSD_g = a_g * logCPM_g + b_g`.
- Fallback: genes with fewer than {args.min_nonzero_samples} nonzero training samples or near-zero logCPM SD use the training VSD mean.
- Clipping: projected values are clipped to the Bo2023 training VSD 0.5%-99.5% quantile range per gene.
- Boundary: projected values are Bo2023-like projected values, not native VSD and not a replacement for the locked production route before fold-local validation.

Full-training fit summary:

```json
{pd.Series(fit_summary).to_json(indent=2)}
```
"""
    (args.outdir / "method_note_reference_projection.md").write_text(note, encoding="utf-8")
    print(f"Wrote projector outputs to {args.outdir}")
    print(f"n_common_samples={len(common_samples)} n_common_genes={len(common_genes)}")
    print(f"median_sample_pearson={fit_summary['median_sample_pearson']:.6f} mae={fit_summary['mae']:.6f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
