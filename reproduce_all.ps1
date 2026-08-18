param(
    [switch]$Quick
)

$ErrorActionPreference = "Stop"
$Root = Split-Path -Parent $MyInvocation.MyCommand.Path
Set-Location -LiteralPath $Root
$Out = Join-Path $Root "reproduced_runs"
New-Item -ItemType Directory -Force -Path $Out | Out-Null

function Run-Python {
    param(
        [Parameter(Mandatory = $true)]
        [string]$Name,
        [Parameter(ValueFromRemainingArguments = $true)]
        [string[]]$Arguments
    )
    Write-Host "`n=== $Name ==="
    & python @Arguments
    if ($LASTEXITCODE -ne 0) {
        throw "$Name failed with exit code $LASTEXITCODE"
    }
}

Run-Python "Unit and contract tests" -m pytest tests -q

if ($Quick) {
    Write-Host "`nQuick verification complete."
    exit 0
}

Run-Python "Canonical three-tier LOSO" `
    scripts/run_bo2023_hybrid_formal_loso.py `
    --outdir (Join-Path $Out "loso")

Run-Python "Frozen no-pairwise three-tier LOMO" `
    scripts/run_bo2023_projected_vsd_formal_lomo.py `
    --outdir (Join-Path $Out "lomo")

Run-Python "AHBA mapped-label transfer" `
    scripts/run_ahba_projected_vsd_formal_three_tier_external.py `
    --outdir (Join-Path $Out "ahba")

Run-Python "GSE189919 frozen route" `
    scripts/run_gse189919_latest_main_route.py `
    --data-dir (Join-Path $Root "data/external_validation/GSE189919") `
    --outdir (Join-Path $Out "gse189919") `
    --db-path (Join-Path $Root "braintrace_source_tracing.db")

Run-Python "Huang 2025 provenance-remediated external domain audit" `
    scripts/run_huang2025_external_candidate.py `
    --input-csv (Join-Path $Root "external_inputs/huang2025_pmc12041490/41698_2025_909_MOESM2_ESM.csv") `
    --source-xlsb (Join-Path $Root "external_inputs/huang2025_pmc12041490/41698_2025_909_MOESM2_ESM.xlsb") `
    --outdir (Join-Path $Root "reproducibility/huang_2025") `
    --db-path (Join-Path $Root "braintrace_source_tracing.db")

Run-Python "Sparse/domain-shift 30-repeat analysis" `
    scripts/run_p0_4_sparse_domain_shift_sensitivity.py `
    --replicates 30 `
    --seed 20260711 `
    --bootstrap-seed 20260716 `
    --n-bootstrap 50000 `
    --outdir (Join-Path $Out "sparse_30_repeat")

Run-Python "Donor-aware inference" `
    validation_runs/r02_small_donor_20260716/analyze_r02_small_donor.py `
    --loso-root (Join-Path $Out "loso") `
    --lomo-root (Join-Path $Out "lomo") `
    --metadata (Join-Path $Root "bo2023 data/Information of sequenced samples_update_full878_filter819.xlsx") `
    --bootstrap-reps 50000 `
    --seed 20260716 `
    --outdir (Join-Path $Out "donor_aware")

Run-Python "R08 internal and AHBA audit" `
    validation_runs/r08_high_feasibility_20260717/analyze_internal_ahba.py `
    --metadata (Join-Path $Root "bo2023 data/Information of sequenced samples_update_full878_filter819.xlsx") `
    --loso-network (Join-Path $Out "loso/hybrid_formal_loso_network_detail.csv") `
    --loso-exact (Join-Path $Out "loso/hybrid_formal_loso_exact_region_detail.csv") `
    --lomo-network (Join-Path $Out "lomo/formal_lomo_network_detail.csv") `
    --lomo-exact (Join-Path $Out "lomo/formal_lomo_exact_region_detail.csv") `
    --ahba-detail (Join-Path $Out "ahba/ahba_formal_three_tier_sample_detail.csv") `
    --ahba-replicate-audit (Join-Path $Out "ahba/ahba_technical_replicate_collapse_audit.csv") `
    --group-model (Join-Path $Root "data/models/bo2023_region_resolution_groups.json") `
    --reference-matrix (Join-Path $Root "data/models/bo2023_formal_region_logcpm_reference_matrix.npz") `
    --bootstrap-reps 50000 `
    --seed 20260717 `
    --outdir (Join-Path $Out "r08_internal_ahba")

Run-Python "Humanization audit" `
    validation_runs/r08_high_feasibility_20260717/analyze_humanization.py `
    --db (Join-Path $Root "braintrace_source_tracing.db") `
    --orthology (Join-Path $Root "data/orthology/ensembl_mfascicularis_hsapiens_homology.tsv") `
    --metadata (Join-Path $Root "bo2023 data/Information of sequenced samples_update_full878_filter819.xlsx") `
    --curation-map (Join-Path $Root "reports/supporting_inputs/bo2023_publication_label_curation_map_20260704.csv") `
    --top200 (Join-Path $Root "data/models/bo2023_saleem_network_top200_model_genes.csv") `
    --projector (Join-Path $Root "data/models/bo2023_reference_projector_linear_full.npz") `
    --reference-matrix (Join-Path $Root "data/models/bo2023_formal_region_logcpm_reference_matrix.npz") `
    --outdir (Join-Path $Out "humanization")

Run-Python "TCGA/BraTS and archived RF audit" `
    validation_runs/r08_high_feasibility_20260717/audit_tcga_rf.py `
    --tcga-web-detail (Join-Path $Root "reports/validation_recheck_20260713/web_core_tcga/web_core_tcga_sample_detail.csv") `
    --tcga-formal-detail (Join-Path $Root "reports/validation_recheck_20260713/tcga_65/tcga_labeled_hybrid_formal_sample_detail.csv") `
    --tcga-truth (Join-Path $Root "results/brats_tcga_lgg_65_mri_truth_corrected_20260612/corrected_direct_overlap_mri_truth.csv") `
    --tcga-claim-script (Join-Path $Root "scripts/apply_tcga_denominator_correction.py") `
    --tcga-supplement-script (Join-Path $Root "scripts/update_supplementary_methods_results_p0.py") `
    --legacy-baseline-script (Join-Path $Root "scripts/generate_p0_hard_evidence.py") `
    --rf-detail (Join-Path $Root "reports/archived_rf_ml_baselines/simple_ml_lomo_detail.csv") `
    --rf-metrics (Join-Path $Root "reports/archived_rf_ml_baselines/simple_ml_lomo_metrics.csv") `
    --rf-fold-audit (Join-Path $Root "reports/archived_rf_ml_baselines/simple_ml_lomo_fold_audit.csv") `
    --rf-code (Join-Path $Root "scripts/generate_p2_publication_completeness.py") `
    --lomo-network (Join-Path $Out "lomo/formal_lomo_network_detail.csv") `
    --outdir (Join-Path $Out "tcga_rf_audit")

Run-Python "Same-fold RF comparator" `
    validation_runs/r08_rf_fair_comparator_20260717/run_rf_fair_comparator.py `
    --mode full `
    --authorize-full `
    --vsd (Join-Path $Root "bo2023 data/mfas5_819samples_23605genes_vsd4_rmbatch.xls") `
    --network-detail (Join-Path $Out "lomo/formal_lomo_network_detail.csv") `
    --exact-detail (Join-Path $Out "lomo/formal_lomo_exact_region_detail.csv") `
    --group-detail (Join-Path $Out "lomo/formal_lomo_resolution_group_detail.csv") `
    --reference-matrix (Join-Path $Root "data/models/bo2023_formal_region_logcpm_reference_matrix.npz") `
    --historical-rf-code (Join-Path $Root "scripts/generate_p2_publication_completeness.py") `
    --group-code (Join-Path $Root "scripts/run_bo2023_resolution_tier_validation.py") `
    --outdir (Join-Path $Out "rf_fair_comparator")

foreach ($Workload in 1, 8, 51) {
    Run-Python "Real-input performance workload $Workload" `
        scripts/benchmark_real_input_inference.py `
        --data-dir (Join-Path $Root "data/external_validation/GSE189919") `
        --outdir (Join-Path $Out "performance_$Workload") `
        --formal-workload $Workload `
        --authorize-formal
}

Write-Host "`nAll registered workflows completed. Outputs: $Out"

