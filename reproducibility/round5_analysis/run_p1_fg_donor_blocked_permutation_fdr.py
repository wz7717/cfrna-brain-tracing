"""P1 E3: donor-blocked permutation FDR for the F_g ranking score.

This is the formal donor-aware sensitivity requested for P1 E3.  It does not
replace the donor-by-Network DESeq2 LRT.  Repeated dissections are first
collapsed within each observed donor x region block, so a donor-region block
is the permutation unit.  After logCPM transformation, each gene is centered
within donor.  Network labels are then permuted within donor, preserving each
donor's observed block count and Network composition.

The earlier 110-region label-permutation run is intentionally left untouched
in ``p1_fg_permutation_fdr`` as a historical, non-donor-aware sensitivity.
"""

from __future__ import annotations

import hashlib
import json
import platform
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd


WORKSPACE = Path(__file__).resolve().parents[3]
SYNC_ROOT = Path(__file__).resolve().parents[2]
DATA_ROOT = WORKSPACE / "code" / "reproduction_validation_workspace_20260802"
BO_ROOT = DATA_ROOT / "bo2023 data"
OUTDIR = SYNC_ROOT / "reproducibility" / "round5_analysis" / "p1_fg_donor_blocked_permutation_fdr"

SEED = 20260809
N_PERMUTATIONS = 5000
PANEL_SIZE = 200
EPSILON = 1e-8


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def bh_adjust(p_values: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg adjusted p values with the required monotonic step."""
    p_values = np.asarray(p_values, dtype=float)
    m = p_values.size
    order = np.argsort(p_values, kind="mergesort")
    ranked = p_values[order] * m / np.arange(1, m + 1, dtype=float)
    monotone = np.minimum.accumulate(ranked[::-1])[::-1]
    q_values = np.empty_like(monotone)
    q_values[order] = np.minimum(monotone, 1.0)
    return q_values


def load_counts(counts_path: Path) -> tuple[list[str], pd.DataFrame]:
    # The featureCounts export has sample IDs on the first line and no gene-ID
    # header; reading with the inferred header would silently drop one sample.
    with counts_path.open("r", encoding="utf-8") as handle:
        sample_ids = handle.readline().rstrip("\r\n").split("\t")
    if len(sample_ids) != 819 or len(set(sample_ids)) != 819:
        raise ValueError(f"expected 819 unique count columns, found {len(sample_ids)}")
    counts = pd.read_csv(
        counts_path,
        sep="\t",
        skiprows=1,
        header=None,
        names=["gene_id", *sample_ids],
        dtype={sample_id: "int64" for sample_id in sample_ids},
    )
    if counts.shape[1] != 820 or counts["gene_id"].duplicated().any():
        raise ValueError("unexpected raw count table shape or duplicate gene IDs")
    counts["gene_id"] = counts["gene_id"].astype(str)
    return sample_ids, counts


def load_metadata(metadata_path: Path, sample_ids: list[str]) -> pd.DataFrame:
    metadata = pd.read_excel(metadata_path, sheet_name="mfas5_819samples_phenSet4", dtype=str)
    required = ["No.", "MonkeyID", "Region", "SaleemNetworks"]
    missing = sorted(set(required) - set(metadata.columns))
    if missing:
        raise ValueError(f"metadata missing columns: {missing}")
    metadata = metadata[required].rename(
        columns={"No.": "sample_id", "MonkeyID": "donor_id", "SaleemNetworks": "network"}
    )
    for column in metadata.columns:
        metadata[column] = metadata[column].fillna("").astype(str).str.strip()
    # Match the canonicalization used by the formal P0-3 DESeq2 audit.
    metadata.loc[metadata["Region"].eq("10m"), "network"] = "Orbitomedial Prefrontal Cortex (OMPFC)"
    metadata.loc[metadata["Region"].eq("V2"), "network"] = "Occipital/Temporal"
    if metadata["sample_id"].duplicated().any():
        raise ValueError("metadata sample IDs are not unique")
    if set(sample_ids) != set(metadata["sample_id"]):
        raise ValueError("count-table and metadata sample IDs do not match")
    metadata = metadata.set_index("sample_id").loc[sample_ids].reset_index()
    if metadata[["donor_id", "Region", "network"]].eq("").any().any():
        raise ValueError("blank donor, region, or Network values remain")
    if metadata["donor_id"].nunique() != 9 or metadata["network"].nunique() != 10:
        raise ValueError(
            f"unexpected design: donors={metadata['donor_id'].nunique()}, "
            f"Networks={metadata['network'].nunique()}"
        )
    if metadata["Region"].nunique() != 110:
        raise ValueError(f"expected 110 regions, found {metadata['Region'].nunique()}")
    return metadata


def prepare_model_matrix(
    counts: pd.DataFrame,
    mapping_path: Path,
    model_path: Path,
) -> tuple[np.ndarray, np.ndarray, int, int]:
    mapping = pd.read_csv(mapping_path, dtype=str)
    required = {"gene_id", "gene_symbol", "mapped_to_symbol"}
    if not required.issubset(mapping.columns):
        raise ValueError(f"gene mapping lacks columns: {sorted(required - set(mapping.columns))}")
    mapping = mapping.drop_duplicates("gene_id", keep="first").set_index("gene_id")
    symbols = counts["gene_id"].map(mapping["gene_symbol"])
    # ``gene_symbol`` is populated with the original Ensembl ID for unmapped
    # rows in the audit file, so count true symbol mappings from the explicit
    # audit flag rather than from the non-null fallback symbol column.
    mapped = mapping.loc[counts["gene_id"], "mapped_to_symbol"].astype(str).str.lower().eq("true").to_numpy()
    symbols = symbols.fillna(counts["gene_id"])
    symbols = symbols.where(symbols.ne(""), counts["gene_id"])
    count_values = counts.drop(columns="gene_id")
    count_values["gene_symbol"] = symbols.to_numpy()
    symbol_counts = count_values.groupby("gene_symbol", sort=False)[count_values.columns[:-1].tolist()].sum()

    model = np.load(model_path, allow_pickle=True)
    model_genes = model["genes"].astype(str)
    if len(model_genes) != 21668 or len(set(model_genes)) != len(model_genes):
        raise ValueError("frozen model gene universe is not the expected unique 21,668 genes")
    missing = sorted(set(model_genes) - set(symbol_counts.index))
    if missing:
        # Missing model symbols are retained as zero-count rows so the tested
        # model space is exactly the frozen 21,668-gene universe.
        symbol_counts = symbol_counts.reindex(symbol_counts.index.union(missing), fill_value=0)
    model_counts = symbol_counts.reindex(model_genes, fill_value=0).to_numpy(dtype=np.float64)
    return model_counts, model_genes, len(mapping), int(mapped.sum())


def build_donor_region_blocks(counts: np.ndarray, metadata: pd.DataFrame) -> tuple[np.ndarray, pd.DataFrame]:
    block_keys = metadata["donor_id"] + "::" + metadata["Region"]
    unique_keys = pd.unique(block_keys)
    block_rows: list[np.ndarray] = []
    block_meta: list[dict[str, object]] = []
    for key in unique_keys:
        sample_idx = np.flatnonzero(block_keys.to_numpy() == key)
        first = metadata.iloc[int(sample_idx[0])]
        networks = metadata.iloc[sample_idx]["network"].unique().tolist()
        if len(networks) != 1:
            raise ValueError(f"donor-region block maps to multiple Networks: {key}")
        block_rows.append(counts[:, sample_idx].sum(axis=1))
        block_meta.append(
            {
                "block_id": str(key),
                "donor_id": str(first["donor_id"]),
                "region": str(first["Region"]),
                "network": str(networks[0]),
                "n_tissue_samples": int(len(sample_idx)),
            }
        )
    block_counts = np.column_stack(block_rows)
    block_metadata = pd.DataFrame(block_meta)
    if len(block_metadata) != 459:
        raise ValueError(f"expected 459 donor-region blocks, found {len(block_metadata)}")
    return block_counts, block_metadata


def logcpm_donor_center(block_counts: np.ndarray, block_metadata: pd.DataFrame) -> np.ndarray:
    library_sizes = block_counts.sum(axis=0)
    if np.any(library_sizes <= 0):
        raise ValueError("one or more donor-region blocks has zero library size")
    logcpm = np.log1p(block_counts / library_sizes[np.newaxis, :] * 1_000_000.0)
    centered = logcpm.copy()
    for donor in sorted(block_metadata["donor_id"].unique()):
        idx = np.flatnonzero(block_metadata["donor_id"].to_numpy() == donor)
        centered[:, idx] -= centered[:, idx].mean(axis=1, keepdims=True)
    return centered


def compute_fg(expr: np.ndarray, group_index: np.ndarray, n_groups: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Compute F_g = B_g/(W_g + epsilon) for all genes."""
    n_genes, n_obs = expr.shape
    one_hot = np.eye(n_groups, dtype=np.float64)[group_index]
    n_by_group = one_hot.sum(axis=0)
    if np.any(n_by_group < 2):
        raise ValueError("a permuted Network has fewer than two donor-region blocks")
    sums = expr @ one_hot
    sums_sq = (expr * expr) @ one_hot
    means = sums / n_by_group[np.newaxis, :]
    grand = expr.mean(axis=1)
    between = (n_by_group[np.newaxis, :] * (means - grand[:, np.newaxis]) ** 2).sum(axis=1) / (n_groups - 1)
    within_ss = (sums_sq - (sums * sums) / n_by_group[np.newaxis, :]).sum(axis=1)
    within = np.maximum(within_ss, 0.0) / (n_obs - n_groups)
    fg = between / (within + EPSILON)
    return fg, between, within


def compare_deseq2(model_genes: np.ndarray, fg_obs: np.ndarray, locked_mask: np.ndarray, top200_idx: np.ndarray) -> dict[str, object]:
    path = SYNC_ROOT / "reproducibility" / "round5_analysis" / "p0_3_deseq2_marker_audit" / "outputs" / "primary_pseudobulk" / "pseudobulk_deseq2_network_lrt_all_genes.csv"
    if not path.exists():
        return {"note": "DESeq2 output not found"}
    de = pd.read_csv(path)
    if not {"gene_symbol", "padj"}.issubset(de.columns):
        return {"note": "DESeq2 output lacks gene_symbol/padj columns"}
    de = de.dropna(subset=["gene_symbol", "padj"]).sort_values(["padj", "gene_symbol"])
    de = de.drop_duplicates("gene_symbol", keep="first")
    de_map = de.set_index("gene_symbol")["padj"].to_dict()
    padj = np.array([de_map.get(g, np.nan) for g in model_genes], dtype=float)
    valid = np.isfinite(padj) & np.isfinite(fg_obs)
    valid_nonzero = valid & (fg_obs > 0) & (padj > 0)
    # Keep the analysis runnable in the bundled runtime (which intentionally
    # does not ship SciPy): pandas' rank transform followed by Pearson
    # correlation is the Spearman rho definition and is tie-aware.
    rank_fg = pd.Series(fg_obs[valid_nonzero]).rank(method="average").to_numpy()
    rank_de = pd.Series(-np.log10(padj[valid_nonzero])).rank(method="average").to_numpy()
    rho = float(np.corrcoef(rank_fg, rank_de)[0, 1])
    de_sig = valid & (padj < 0.05)
    return {
        "source": str(path.relative_to(WORKSPACE)),
        "n_deseq2_rows": int(len(de)),
        "n_model_genes_with_deseq2_padj": int(valid.sum()),
        "n_model_genes_deseq2_bh_005": int(de_sig.sum()),
        "recomputed_top200_deseq2_bh_005": int(de_sig[top200_idx].sum()),
        "locked_top200_deseq2_bh_005": int(de_sig[locked_mask].sum()),
        "spearman_fg_vs_neglog10_padj": rho,
        "spearman_p": None,
        "spearman_p_note": "rho computed from tie-aware ranks; no SciPy p-value is used in this sensitivity comparison",
        "interpretation": "DESeq2 remains a separate formal donor-by-Network LRT; this comparison is descriptive only.",
    }


def main() -> int:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    counts_path = next(BO_ROOT.glob("*featurecounts_counts.txt"))
    metadata_path = next(BO_ROOT.glob("Information*.xlsx"))
    mapping_path = DATA_ROOT / "reproducibility" / "p0_bio3_projector" / "count_gene_symbol_mapping_audit.csv"
    model_path = DATA_ROOT / "data" / "models" / "bo2023_formal_region_logcpm_reference_matrix.npz"
    locked_path = DATA_ROOT / "data" / "models" / "bo2023_saleem_network_top200_model_genes.csv"

    print("Loading raw counts and metadata...")
    sample_ids, counts = load_counts(counts_path)
    metadata = load_metadata(metadata_path, sample_ids)
    model_counts, model_genes, n_mapping_rows, n_mapped_rows = prepare_model_matrix(counts, mapping_path, model_path)
    locked = pd.read_csv(locked_path, dtype=str)
    locked_genes = set(locked["gene_symbol"].astype(str))
    if len(locked_genes) != PANEL_SIZE or not locked_genes.issubset(set(model_genes)):
        raise ValueError("frozen Top200 panel is not exactly represented in the model gene universe")
    locked_mask = np.isin(model_genes, list(locked_genes))

    print("Collapsing repeated dissections to donor-region blocks...")
    block_counts, block_metadata = build_donor_region_blocks(model_counts, metadata)
    block_metadata.to_csv(OUTDIR / "donor_region_block_metadata.tsv", sep="\t", index=False)
    expr = logcpm_donor_center(block_counts, block_metadata)
    network_names = sorted(block_metadata["network"].unique())
    group_index = pd.Categorical(block_metadata["network"], categories=network_names).codes
    if len(network_names) != 10:
        raise ValueError(f"expected 10 Networks after blocking, found {len(network_names)}")

    print("Computing observed F_g...")
    fg_obs, between_obs, within_obs = compute_fg(expr, group_index, len(network_names))
    observed_order = np.argsort(-fg_obs, kind="mergesort")
    fg_rank = np.empty(len(model_genes), dtype=int)
    fg_rank[observed_order] = np.arange(1, len(model_genes) + 1)
    recomputed_top200_idx = observed_order[:PANEL_SIZE]
    recomputed_top200_genes = set(model_genes[recomputed_top200_idx])
    overlap = len(recomputed_top200_genes & locked_genes)

    print(f"Running {N_PERMUTATIONS} donor-blocked permutations (seed={SEED})...")
    rng = np.random.default_rng(SEED)
    donor_block_indices = [
        np.flatnonzero(block_metadata["donor_id"].to_numpy() == donor)
        for donor in sorted(block_metadata["donor_id"].unique())
    ]
    exceed_counts = np.zeros(len(model_genes), dtype=np.int64)
    for iteration in range(N_PERMUTATIONS):
        permuted_groups = group_index.copy()
        for indices in donor_block_indices:
            permuted_groups[indices] = rng.permutation(permuted_groups[indices])
        fg_perm, _, _ = compute_fg(expr, permuted_groups, len(network_names))
        exceed_counts += fg_perm >= fg_obs
        if iteration == 0 or (iteration + 1) % 500 == 0:
            print(f"  permutation {iteration + 1}/{N_PERMUTATIONS}")

    p_values = (exceed_counts + 1.0) / (N_PERMUTATIONS + 1.0)
    q_values = bh_adjust(p_values)
    sig_005 = q_values < 0.05
    sig_001 = q_values < 0.01
    top200_sig = int(sig_005[recomputed_top200_idx].sum())
    locked_sig = int(sig_005[locked_mask].sum())
    top200_max_q = float(q_values[recomputed_top200_idx].max())
    locked_max_q = float(q_values[locked_mask].max())

    de_comparison = compare_deseq2(model_genes, fg_obs, locked_mask, recomputed_top200_idx)
    result = pd.DataFrame(
        {
            "gene_symbol": model_genes,
            "fg_score": fg_obs,
            "fg_rank": fg_rank,
            "between_variance": between_obs,
            "within_variance": within_obs,
            "permutation_exceed_count": exceed_counts,
            "permutation_p_value": p_values,
            "bh_fdr": q_values,
            "bh_fdr_005": sig_005,
            "bh_fdr_001": sig_001,
            "in_locked_top200": locked_mask,
            "in_recomputed_top200": np.isin(np.arange(len(model_genes)), recomputed_top200_idx),
        }
    )
    result.to_csv(OUTDIR / "fg_donor_blocked_permutation_fdr_per_gene.csv", index=False)
    result[result["in_locked_top200"]].sort_values("fg_rank").to_csv(
        OUTDIR / "locked_top200_donor_blocked_permutation_fdr.csv", index=False
    )
    np.save(OUTDIR / "permutation_exceed_counts.npy", exceed_counts)
    np.save(OUTDIR / "fg_observed.npy", fg_obs)

    donor_counts = block_metadata.groupby("donor_id").size().astype(int).to_dict()
    summary = {
        "analysis": "Donor-blocked permutation FDR for F_g feature-ranking scores",
        "status": "formal P1 E3 donor-aware sensitivity completed",
        "n_tissue_samples": int(len(metadata)),
        "n_raw_genes": int(counts.shape[0]),
        "n_model_genes": int(len(model_genes)),
        "n_mapping_rows": int(n_mapping_rows),
        "n_mapped_gene_rows": int(n_mapped_rows),
        "n_donor_region_blocks": int(len(block_metadata)),
        "n_donors": int(block_metadata["donor_id"].nunique()),
        "n_regions": int(block_metadata["region"].nunique()),
        "n_networks": int(len(network_names)),
        "donor_region_block_counts": donor_counts,
        "n_permutations": N_PERMUTATIONS,
        "seed": SEED,
        "permutation_unit": "donor x observed region block; repeated dissections summed within block",
        "preprocessing": "raw counts summed by donor-region, logCPM, then each gene centered within donor",
        "permutation_scheme": "Network labels permuted within donor, preserving donor block counts and observed Network composition",
        "p_value_correction": "add-one empirical p=(exceed+1)/(N_perm+1); Benjamini-Hochberg across 21,668 model genes",
        "fg_definition": "between-Network variance divided by within-Network variance, with df_between=9 and df_within=449",
        "genes_bh_fdr_005": int(sig_005.sum()),
        "genes_bh_fdr_005_pct": round(float(sig_005.mean() * 100), 3),
        "genes_bh_fdr_001": int(sig_001.sum()),
        "recomputed_top200_bh_fdr_005": top200_sig,
        "locked_top200_bh_fdr_005": locked_sig,
        "locked_top200_max_bh_fdr": locked_max_q,
        "recomputed_top200_max_bh_fdr": top200_max_q,
        "recomputed_top200_vs_locked_overlap": int(overlap),
        "locked_panel_reselected": False,
        "deseq2_comparison": de_comparison,
        "interpretation": "This donor-blocked permutation is a sensitivity analysis for the descriptive F_g score; it is not a replacement for the formal donor-by-Network DESeq2 pseudobulk LRT.",
        "historical_non_donor_run": "reproducibility/round5_analysis/p1_fg_permutation_fdr/ (110-region, 1,000-permutation exploratory sensitivity; not used for the donor-aware conclusion)",
    }
    (OUTDIR / "fg_donor_blocked_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )

    protocol = f"""# P1 E3 donor-blocked F_g permutation FDR

This run is the formal donor-aware sensitivity for the descriptive F_g
Network-ranking score. It does not replace the donor-by-Network DESeq2 LRT.

- Input: {len(metadata)} Bo2023 tissue samples, {counts.shape[0]} raw genes.
- Blocking: repeated dissections were summed within each observed donor x
  region block, producing {len(block_metadata)} blocks from
  {block_metadata['donor_id'].nunique()} donors and {block_metadata['region'].nunique()} regions.
- Transformation: logCPM on the block-level raw counts, followed by gene-wise
  centering within donor.
- Statistic: F_g = between-Network variance / within-Network variance, with
  df_between=9 and df_within=449.
- Null: Network labels were permuted within donor, preserving each donor's
  block count and observed Network composition.
- Multiplicity: {N_PERMUTATIONS} permutations, seed {SEED}, empirical p-values
  with the add-one correction, then BH across the {len(model_genes):,}-gene
  frozen model universe.
- Production panel: the frozen Top200 was not reselected or reordered.

The earlier `p1_fg_permutation_fdr` directory is retained as a historical
110-region, non-donor-aware exploratory sensitivity and must not be described
as donor-aware.
"""
    (OUTDIR / "PROTOCOL.md").write_text(protocol, encoding="utf-8")

    manifest = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "script": str(Path(__file__).relative_to(WORKSPACE)),
        "python": sys.version,
        "platform": platform.platform(),
        "inputs": {
            "counts": {"path": str(counts_path.relative_to(WORKSPACE)), "sha256": sha256(counts_path)},
            "metadata": {"path": str(metadata_path.relative_to(WORKSPACE)), "sha256": sha256(metadata_path)},
            "gene_mapping": {"path": str(mapping_path.relative_to(WORKSPACE)), "sha256": sha256(mapping_path)},
            "frozen_model": {"path": str(model_path.relative_to(WORKSPACE)), "sha256": sha256(model_path)},
            "locked_top200": {"path": str(locked_path.relative_to(WORKSPACE)), "sha256": sha256(locked_path)},
        },
        "outputs": {
            "summary": "fg_donor_blocked_summary.json",
            "per_gene": "fg_donor_blocked_permutation_fdr_per_gene.csv",
            "locked_panel": "locked_top200_donor_blocked_permutation_fdr.csv",
            "block_metadata": "donor_region_block_metadata.tsv",
            "exceed_counts": "permutation_exceed_counts.npy",
        },
    }
    (OUTDIR / "run_manifest.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    checksum_files = sorted(p for p in OUTDIR.iterdir() if p.is_file() and p.name != "SHA256SUMS.txt")
    with (OUTDIR / "SHA256SUMS.txt").open("w", encoding="utf-8") as handle:
        for path in checksum_files:
            handle.write(f"{sha256(path)}  {path.name}\n")

    print(json.dumps(summary, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
