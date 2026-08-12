"""Build the frozen 65-case TCGA TPM mixture matrix for CIBERSORTx.

This prepares an upload-ready, non-log, non-negative mixture file. It does not
run CIBERSORTx and deliberately does not redistribute LM22 or invent a brain
cell-type signature matrix.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workspace", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    detail = args.workspace / "reports/model_validation/model_validation_stage2_20260716/tcga_brats_frozen_no_pairwise/tcga_labeled_hybrid_formal_sample_detail.csv"
    if not detail.exists():
        detail = args.workspace / "reports/validation_recheck_20260713/tcga_65/tcga_labeled_hybrid_formal_sample_detail.csv"
    tpm = args.workspace / "data/tcga_brain_tumor_expression/tcga_gbm_lgg_primary_tumor_tpm_unstranded_sample_mean.tsv"
    if not detail.exists() or not tpm.exists():
        raise FileNotFoundError("Frozen TCGA detail or TPM matrix is missing")

    detail_df = pd.read_csv(detail)
    sample_ids = sorted(detail_df["sample_id"].dropna().astype(str).unique())
    if len(sample_ids) != 65:
        raise ValueError(f"Expected 65 frozen TCGA samples, found {len(sample_ids)}")

    header = pd.read_csv(tpm, sep="\t", nrows=0).columns.tolist()
    missing = sorted(set(sample_ids) - set(header))
    if missing:
        raise ValueError(f"Frozen samples absent from TPM matrix: {missing}")
    matrix = pd.read_csv(tpm, sep="\t", usecols=["gene_symbol", *sample_ids])
    matrix = matrix.dropna(subset=["gene_symbol"])
    matrix["gene_symbol"] = matrix["gene_symbol"].astype(str)
    matrix = matrix[matrix["gene_symbol"].str.len().gt(0)]
    matrix = matrix.groupby("gene_symbol", as_index=False)[sample_ids].mean()
    if matrix[sample_ids].isna().any().any() or (matrix[sample_ids] < 0).any().any():
        raise ValueError("CIBERSORTx mixture must be non-negative and contain no missing values")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    output = args.output_dir / "cibersortx_tcga65_tpm_mixture.tsv"
    matrix.rename(columns={"gene_symbol": "GeneSymbol"}).to_csv(output, sep="\t", index=False)
    manifest = {
        "status": "upload-ready mixture; CIBERSORTx not yet run",
        "sample_n": len(sample_ids),
        "gene_n": len(matrix),
        "scale": "GDC STAR TPM unstranded, sample-level mean, non-log linear",
        "quantile_normalization": False,
        "intended_signature": "LM22 only for immune-composition sensitivity, not brain-region localization",
        "source_detail": str(detail),
        "source_tpm": str(tpm),
        "source_detail_sha256": sha256(detail),
        "source_tpm_sha256": sha256(tpm),
        "output_sha256": sha256(output),
    }
    (args.output_dir / "cibersortx_tcga65_input_manifest.json").write_text(
        json.dumps(manifest, indent=2), encoding="utf-8"
    )


if __name__ == "__main__":
    main()
