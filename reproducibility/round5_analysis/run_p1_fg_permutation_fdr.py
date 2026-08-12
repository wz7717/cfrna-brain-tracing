"""
P1 [E3]: Donor-aware permutation FDR for F_g feature-ranking scores.

Computes F_g on Bo2023 region-level logCPM expression (110 regions, 10 Networks),
then runs network-label permutation (1000 perm) to build null distribution.
Computes permutation p-values, BH-FDR, and compares with:
- Locked Top200 gene panel
- Existing DESeq2 donor-aware pseudobulk LRT (separate audit)

Output: round5_analysis/p1_fg_permutation_fdr/
"""

import numpy as np
import pandas as pd
from pathlib import Path
import json

ROOT = Path(__file__).resolve().parents[2]
OUTDIR = ROOT / "reproducibility" / "round5_analysis" / "p1_fg_permutation_fdr"
OUTDIR.mkdir(parents=True, exist_ok=True)

SEED = 20260805
N_PERM = 1000
rng = np.random.default_rng(SEED)

# ═══ 1. Load data ═══
print("1. Loading Bo2023 region-level expression...")
mat = np.load(str(ROOT / "data" / "models" / "bo2023_formal_region_logcpm_reference_matrix.npz"))
gene_names = mat["genes"]      # 21668
region_labels = mat["regions"] # 110
network_per_region = np.array(mat["networks"])  # 110
expr = mat["matrix"].astype(np.float64)  # 21668 × 110
n_genes, n_regions = expr.shape
unique_networks = sorted(set(network_per_region))
n_networks = len(unique_networks)
print(f"  Genes: {n_genes}, Regions: {n_regions}, Networks: {n_networks}")

# ═══ 2. Observed F_g ═══
print("\n2. Computing observed F_g...")

def compute_fg_vectorized(expr_matrix, groups):
    """Vectorized F_g = B_g/(W_g+1e-8) for all genes.
    expr_matrix: genes × obs, groups: obs group labels"""
    n_g, n_o = expr_matrix.shape
    uniq = sorted(set(groups))
    n_grp = len(uniq)

    fg = np.zeros(n_g)
    between_v = np.zeros(n_g)
    within_v = np.zeros(n_g)

    grand_mean = expr_matrix.mean(axis=1)

    for grp in uniq:
        mask = groups == grp
        n_k = mask.sum()
        grp_expr = expr_matrix[:, mask]
        grp_mean = grp_expr.mean(axis=1)
        between_v += n_k * (grp_mean - grand_mean)**2 / (n_grp - 1)
        if n_k > 1:
            within_v += ((grp_expr - grp_mean[:, np.newaxis])**2).sum(axis=1)

    df_within = n_o - n_grp
    within_v /= max(df_within, 1)
    fg = between_v / (within_v + 1e-8)
    return fg, between_v, within_v

fg_obs, bg_obs, wg_obs = compute_fg_vectorized(expr, network_per_region)
rank_obs = np.argsort(-fg_obs)
top200_idx = set(rank_obs[:200])

print(f"  Observed F_g: median={np.median(fg_obs):.4f}, Q1={np.percentile(fg_obs,25):.4f}, Q3={np.percentile(fg_obs,75):.4f}")

# Load locked Top200
locked = pd.read_csv(str(ROOT / "data" / "models" / "bo2023_saleem_network_top200_model_genes.csv"))
locked_genes = set(locked["gene_symbol"])
overlap = len(set(gene_names[i] for i in top200_idx) & locked_genes)
print(f"  Recomputed Top200 vs locked Top200 overlap: {overlap}/200")

# ═══ 3. Permutation ═══
print(f"\n3. Running {N_PERM} network-label permutations...")
perm_fg_max = np.zeros((N_PERM, n_genes))

for p in range(N_PERM):
    if p % 250 == 0:
        print(f"  Perm {p}/{N_PERM}...")
    perm_nets = network_per_region.copy()
    rng.shuffle(perm_nets)
    fg_p, _, _ = compute_fg_vectorized(expr, perm_nets)
    perm_fg_max[p] = fg_p

# ═══ 4. P-values and FDR ═══
print("\n4. Computing permutation p-values and BH-FDR...")
pvals = np.array([(np.sum(perm_fg_max[:, g] >= fg_obs[g]) + 1) / (N_PERM + 1) for g in range(n_genes)])

# BH-FDR
order = np.argsort(pvals)
ranks = np.zeros(n_genes, dtype=int)
ranks[order] = np.arange(1, n_genes + 1)
bh_fdr_vals = np.minimum(1.0, pvals * n_genes / ranks)

sig_005 = bh_fdr_vals < 0.05
sig_001 = bh_fdr_vals < 0.01
n_sig_005 = sig_005.sum()
n_sig_001 = sig_001.sum()

top200_sig_005 = sig_005[list(top200_idx)].sum()
top200_sig_001 = sig_001[list(top200_idx)].sum()
top200_max_p = pvals[list(top200_idx)].max()
top200_max_fdr = bh_fdr_vals[list(top200_idx)].max()

print(f"  BH-FDR < 0.05: {n_sig_005}/{n_genes} ({100*n_sig_005/n_genes:.1f}%)")
print(f"  BH-FDR < 0.01: {n_sig_001}/{n_genes} ({100*n_sig_001/n_genes:.1f}%)")
print(f"  Top200 with BH-FDR < 0.05: {top200_sig_005}/200")
print(f"  Top200 max perm p-value: {top200_max_p:.6e}")
print(f"  Top200 max BH-FDR: {top200_max_fdr:.6e}")

# ═══ 5. DESeq2 comparison ═══
print("\n5. Comparing with DESeq2 donor-aware pseudobulk LRT...")
deseq2_file = ROOT / "reproducibility" / "round5_analysis" / "p0_3_deseq2_marker_audit" / "outputs" / "primary_pseudobulk" / "pseudobulk_deseq2_network_lrt_all_genes.csv"

deseq2_comparison = {}
if deseq2_file.exists():
    from scipy.stats import spearmanr
    de = pd.read_csv(str(deseq2_file))
    print(f"  Loaded DESeq2: {de.shape}, cols: {list(de.columns)[:8]}")

    # Find gene and padj columns
    gene_col = None
    for c in ["gene", "gene_symbol", "Gene", "symbol"]:
        if c in de.columns:
            gene_col = c; break
    padj_col = None
    for c in ["padj", "p_adj", "adjusted_p", "BH_p", "padj_BH"]:
        if c in de.columns:
            padj_col = c; break

    if gene_col and padj_col:
        de_map = dict(zip(de[gene_col], de[padj_col]))
        de_padj = np.array([de_map.get(g, 1.0) for g in gene_names])
        valid = ~np.isnan(de_padj) & ~np.isinf(de_padj) & (fg_obs > 0) & (de_padj < 1.0)
        rho, p_spear = spearmanr(fg_obs[valid], -np.log10(np.maximum(de_padj[valid], 1e-300)))

        de_sig = de_padj < 0.05
        n_de_sig = de_sig.sum()
        top200_de_sig = de_sig[list(top200_idx)].sum()

        print(f"  DESeq2 sig genes: {n_de_sig}/{n_genes}")
        print(f"  Top200 DESeq2 sig: {top200_de_sig}/200")
        print(f"  Spearman(F_g, -log10 DESeq2 padj): {rho:.4f}")

        deseq2_comparison = {
            "n_genes_tested": int(n_genes),
            "n_genes_sig_deseq2_bh_005": int(n_de_sig),
            "top200_deseq2_sig": int(top200_de_sig),
            "spearman_fg_vs_neglog10_padj": float(rho),
            "spearman_p": float(p_spear),
        }
    else:
        deseq2_comparison = {"note": f"Missing gene_col={gene_col} or padj_col={padj_col}"}
else:
    deseq2_comparison = {"note": "DESeq2 file not found"}

# ═══ 6. Save ═══
print("\n6. Saving results...")

results = pd.DataFrame({
    "gene": gene_names,
    "fg_score": fg_obs,
    "fg_rank": ranks,
    "between_variance": bg_obs,
    "within_variance": wg_obs,
    "permutation_p_value": pvals,
    "bh_fdr": bh_fdr_vals,
    "bh_fdr_005": sig_005,
    "bh_fdr_001": sig_001,
    "in_locked_top200": [g in locked_genes for g in gene_names],
})
results.to_csv(OUTDIR / "fg_permutation_fdr_per_gene.csv", index=False)

top200_detail = results[results["in_locked_top200"]].sort_values("fg_rank")
top200_detail.to_csv(OUTDIR / "locked_top200_permutation_fdr.csv", index=False)

summary = {
    "analysis": "Permutation FDR for F_g feature-ranking scores (region-level, 110 regions, 10 Networks)",
    "n_genes": n_genes,
    "n_regions": n_regions,
    "n_networks": n_networks,
    "n_permutations": N_PERM,
    "seed": SEED,
    "fg_obs_median": float(np.median(fg_obs)),
    "fg_obs_q1": float(np.percentile(fg_obs, 25)),
    "fg_obs_q3": float(np.percentile(fg_obs, 75)),
    "genes_bh_fdr_005": int(n_sig_005),
    "genes_bh_fdr_005_pct": round(100*n_sig_005/n_genes, 2),
    "genes_bh_fdr_001": int(n_sig_001),
    "genes_bh_fdr_001_pct": round(100*n_sig_001/n_genes, 2),
    "top200_bh_fdr_005": int(top200_sig_005),
    "top200_bh_fdr_001": int(top200_sig_001),
    "top200_max_permutation_p": float(top200_max_p),
    "top200_max_bh_fdr": float(top200_max_fdr),
    "top200_recomputed_vs_locked_overlap": int(overlap),
    "deseq2_comparison": deseq2_comparison,
    "limitation": "Region-level analysis (110 regions); not weighted by per-region sample count or donor structure. Serves as sensitivity check alongside formal donor-aware DESeq2 pseudobulk LRT (P0-3 audit)."
}

with open(OUTDIR / "fg_permutation_summary.json", "w") as f:
    json.dump(summary, f, indent=2)

np.savez_compressed(OUTDIR / "permutation_fg_max.npz",
    fg_obs=fg_obs, gene_names=gene_names, perm_fg_max=perm_fg_max)

# SHA256 manifest
import hashlib
manifest = {}
for fpath in sorted(OUTDIR.glob("*")):
    if fpath.is_file() and fpath.suffix in [".csv", ".json", ".npz"]:
        h = hashlib.sha256(fpath.read_bytes()).hexdigest()
        manifest[fpath.name] = h
with open(OUTDIR / "SHA256SUMS.txt", "w") as f:
    for k, v in manifest.items():
        f.write(f"{v}  {k}\n")

print(f"\nSaved to {OUTDIR}")
print(f"\n{'='*60}")
print("KEY RESULTS")
print("=" * 60)
print(f"Permutation BH-FDR < 0.05: {n_sig_005}/{n_genes} ({100*n_sig_005/n_genes:.1f}%)")
print(f"Top200 with BH-FDR < 0.05: {top200_sig_005}/200")
print(f"Top200 max BH-FDR: {top200_max_fdr:.6e}")
if "spearman_fg_vs_neglog10_padj" in deseq2_comparison:
    print(f"Spearman(F_g, DESeq2 -log10 padj): {deseq2_comparison['spearman_fg_vs_neglog10_padj']:.4f}")
