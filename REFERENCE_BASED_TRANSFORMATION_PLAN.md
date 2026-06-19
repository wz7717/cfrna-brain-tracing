# Reference-Based Transformation Implementation Plan

## Purpose

This file is a handoff guide for implementing and evaluating a reference-based
transformation branch in this project.

The goal is to test whether external or count-derived samples can be projected
into a Bo2023-like VSD space before running the existing brain-origin tracing
models.

This is an exploratory validation branch. It must not replace the current
locked production route unless it passes strict internal validation.

Current locked production route:

```text
Bo2023 VSD reference
-> fold-selected 200 genes
-> Pearson correlation
-> Top-3 pairwise rescue
```

## Key Idea

Do not let each external dataset create its own unrelated VST/VSD scale.
Instead, learn a transformation from Bo2023 training data:

```text
Bo2023 raw count-derived expression
-> Bo2023 VSD + batch-removed expression
```

Then apply that learned transformation to held-out or external count-derived
samples:

```text
external raw count
-> count-derived expression, usually logCPM
-> Bo2023-trained projector
-> projected Bo2023-like VSD
-> existing tracing model
```

The end-to-end route to evaluate is:

```text
validation dataset raw count files
-> harmonize gene IDs / orthologs / gene symbols
-> compute logCPM
-> apply a VSD parameter projector fitted only on Bo2023 training samples
-> generate a projected Bo2023-like VSD matrix
-> run tracing against the existing Bo2023 VSD reference
-> compare with existing TPM/logCPM/rank routes
```

Validation dataset raw counts are transform-only inputs. They must not be used
to fit projector parameters. The generated matrix is a `projected Bo2023-like
VSD matrix`, not native external-dataset VSD and not a strict replay of DESeq2's
native VSD transformation.

The most practical first implementation is not a full DESeq2 parameter replay.
Use a fold-local empirical projector trained from paired Bo2023 counts and
Bo2023 VSD values.

## Existing Inputs

Bo2023 raw count:

```text
bo2023 data/mfas5_819samples_28415genes_featurecounts_counts.txt
```

Bo2023 VSD + batch removed:

```text
bo2023 data/mfas5_819samples_23605genes_vsd4_rmbatch.xls
```

Bo2023 sample metadata:

```text
bo2023 data/Information of sequenced samples_update_full878_filter819.xlsx
```

Existing reconstructed/frozen VST artifacts to inspect before reimplementing:

```text
results/bo2023_vsd_reconstruction/bo2023_frozen_vst_reference.rds
results/bo2023_vsd_reconstruction/best_reconstructed_vsd.tsv.gz
scripts/reconstruct_bo2023_vsd.R
scripts/prepare_bo2023_vsd_metadata.py
```

External count-capable datasets already present locally:

```text
data/external_validation/GSE106804/GSE106804_Gene_counts.txt.gz
data/external_validation/GSE189919/GSE189919_count.csv.gz
data/external_validation/GSE228512/GSE228512_hiseq_counts.txt.gz
data/external_validation/GSE228512/GSE228512_novaseq_counts.txt.gz
data/tcga_brain_tumor_expression/tcga_gbm_lgg_primary_tumor_unstranded_counts_sample_sum.tsv
data/ahba_human_rnaseq/raw_zips/H0351_2001_rnaseq.zip
data/ahba_human_rnaseq/raw_zips/H0351_2002_rnaseq.zip
```

The AHBA zip files contain `RNAseqCounts.csv`, but the current AHBA validation
script uses `RNAseqTPM.csv`.

Ivy GAP currently has FPKM/TPM files locally, not a strict raw count matrix:

```text
data/ivy_gap_anatomic_rnaseq/ivy_gap_anatomic_structure_fpkm_gene_symbol_matrix.tsv
data/ivy_gap_anatomic_rnaseq/ivy_gap_anatomic_structure_tpm_gene_symbol_matrix.tsv
```

Do not claim strict VST projection for Ivy GAP unless raw counts are obtained.

## Non-Goals

Do not claim that projected values are true Bo2023 VSD unless the method is
validated internally.

Do not use test samples to fit any transformation parameter.

Do not use external cohort labels, MRI labels, tumor location labels, diagnosis
labels, or outcome labels to tune the projector.

Do not make projected VSD the production default until it improves or preserves
the locked Bo2023 LOSO/LOMO behavior.

## Recommended First Projector

Start with a simple, auditable per-gene mapping:

```text
input:  Bo2023 logCPM per gene
target: Bo2023 VSD_batch_removed per gene
model:  robust or ordinary linear regression per gene
```

For gene `g`:

```text
VSD_g = a_g * logCPM_g + b_g
```

Apply clipping or fallback rules:

```text
if gene has too few nonzero training samples:
    use z-score/median fallback or exclude gene
if predicted value is extreme:
    clip to Bo2023 training VSD quantile range, e.g. 0.5%-99.5%
```

Use logCPM from counts:

```text
CPM = count / sample_library_size * 1,000,000
logCPM = log1p(CPM)
```

Keep this implementation separate from production inference until validated.

## Alternative Projectors To Compare

Implement after the first linear projector is working.

1. Per-gene robust linear projector
2. Per-gene quantile mapping
3. Per-gene z-score-to-reference mapping
4. Rank-percentile projector
5. Frozen DESeq2/VST replay, only if the existing R artifacts make it feasible

Rank-based output should be treated as a separate robust baseline, not as VSD.

## Required Fold-Local Discipline

For Bo2023 LOSO:

```text
for each held-out sample:
    train_samples = all Bo2023 samples except held-out sample
    fit projector only on train_samples
    build or use reference only from train_samples
    transform held-out sample counts with train-only projector
    run tracing
```

For Bo2023 LOMO:

```text
for each held-out monkey:
    train_samples = all Bo2023 samples from other monkeys
    fit projector only on train_samples
    build or use reference only from train_samples
    transform all held-out monkey samples with train-only projector
    run tracing
```

For external datasets:

```text
fit projector on Bo2023 training/reference samples only
transform external samples
run tracing
```

If external samples are used to estimate distribution parameters, clearly mark
that route as transductive sensitivity analysis, not locked validation.

## Suggested Implementation Files

Create new files rather than modifying the locked route first:

```text
scripts/build_bo2023_reference_projector.py
scripts/run_bo2023_projected_vsd_loso.py
scripts/run_bo2023_projected_vsd_lomo.py
scripts/apply_projected_vsd_to_external_counts.py
core/reference_projection.py
```

Suggested output directory:

```text
results/bo2023_reference_projection_YYYYMMDD/
```

Suggested model artifact names:

```text
bo2023_reference_projector_linear_fold_{fold_id}.npz
bo2023_reference_projector_linear_full.npz
```

Suggested result files:

```text
projector_gene_parameters.csv
projector_qc_summary.json
bo2023_projected_vsd_loso_detail.csv
bo2023_projected_vsd_loso_summary.json
bo2023_projected_vsd_lomo_detail.csv
bo2023_projected_vsd_lomo_summary.json
external_projected_vsd_<dataset>_detail.csv
external_projected_vsd_<dataset>_summary.json
method_note_reference_projection.md
```

## Phase 1: Data Audit

Tasks:

1. Read Bo2023 raw count matrix.
2. Read Bo2023 VSD/batch-removed matrix.
3. Read Bo2023 sample metadata.
4. Confirm sample IDs match between count and VSD matrices.
5. Confirm gene IDs and gene symbols.
6. Collapse duplicate genes if needed.
7. Intersect count genes and VSD genes.
8. Report missing samples and missing genes.

Deliverables:

```text
results/bo2023_reference_projection_YYYYMMDD/data_audit_summary.json
results/bo2023_reference_projection_YYYYMMDD/common_gene_panel.csv
```

Minimum acceptable audit:

```text
n_common_samples close to 819
n_common_genes large enough for current model genes
all locked model genes either present or explicitly reported missing
```

## Phase 2: Projector Training

Tasks:

1. Convert Bo2023 raw counts to logCPM.
2. Align logCPM and VSD matrices to common genes and samples.
3. Fit per-gene projector on training samples.
4. Save slope, intercept, residual SD, training quantiles, and QC flags.
5. Add fallback behavior for unstable genes.

Recommended QC fields per gene:

```text
gene_symbol
n_train_samples
n_nonzero_count_samples
logcpm_mean
logcpm_sd
vsd_mean
vsd_sd
slope
intercept
r2
spearman_r
residual_sd
fallback_reason
clip_low
clip_high
```

Validation inside training folds:

```text
correlation(projected_vsd, native_vsd)
MAE(projected_vsd, native_vsd)
per-gene R2 distribution
per-sample correlation distribution
```

## Phase 3: Internal Bo2023 LOSO Evaluation

Compare at least these routes:

```text
native_vsd:
    held-out native Bo2023 VSD sample -> fold-local Bo2023 VSD reference

projected_vsd:
    held-out Bo2023 count -> fold-local projector -> projected VSD -> fold-local Bo2023 VSD reference

logcpm_baseline:
    held-out Bo2023 logCPM -> compatible reference or direct correlation baseline

rank_baseline:
    held-out Bo2023 rank-percentile vector -> rank-transformed reference
```

Primary metrics:

```text
Network Top1
Network Top3
Exact Region Top1 if supported
Exact Region Top3 if supported
median true rank
abstain rate
n_overlap_genes
decision margin
```

Pass/fail guidance:

```text
Projected VSD should preserve most native VSD Network Top3 performance.
Projected VSD should outperform or match plain logCPM/rank baseline.
If it performs materially worse than logCPM, do not continue to external claims.
```

## Phase 4: Internal Bo2023 LOMO Evaluation

Repeat the projected-VSD evaluation with leave-one-monkey-out splits.

This is more important than LOSO for generalization because the held-out unit is
a biological individual rather than one sample.

Required metadata:

```text
heldout_monkey_id
n_train_samples
n_test_samples
n_train_regions
n_test_regions
```

Compare:

```text
native_vsd LOMO
projected_vsd LOMO
logcpm/rank baselines
```

## Phase 5: External Count Dataset Projection

Only run after Phase 3 passes.

Recommended order:

1. GSE189919 CSF RNA
2. GSE228512 EV RNA
3. GSE106804 tumor-specific EV RNA
4. AHBA RNA-seq counts
5. TCGA GBM/LGG counts

For each dataset:

1. Read raw counts.
2. Convert IDs to gene symbols and, when needed, human-to-macaque orthologs.
3. Collapse duplicate symbols.
4. Convert to logCPM.
5. Apply full Bo2023-trained projector.
6. Export a projected Bo2023-like VSD matrix.
7. Run existing network tracing.
8. Compare against existing TPM/logCPM/rank routes.
9. Export confidence, margin, top-k distribution, and overlap-gene QC.

External outputs must be described as cross-domain stress tests unless true
anatomical labels exist.

## Phase 6: Reporting

Create a concise method note:

```text
results/bo2023_reference_projection_YYYYMMDD/method_note_reference_projection.md
```

It must state:

```text
The projector was trained from paired Bo2023 raw count-derived logCPM and Bo2023 VSD_batch_removed values.
All internal validation was fold-local.
No held-out sample was used to fit projector parameters.
External projected results are not native Bo2023 VSD and are interpreted as cross-domain projected-space analyses.
```

## Time Estimate

Minimum useful version:

```text
3-5 working days
```

Scope:

```text
Bo2023 data audit
linear per-gene projector
Bo2023 LOSO validation
one external proof-of-concept dataset
```

Research-grade version:

```text
7-10 working days
```

Scope:

```text
LOSO and LOMO
multiple projector variants
GSE106804, GSE189919, GSE228512, AHBA, and TCGA
QC reports and method note
```

Supplementary-material grade:

```text
10-15 working days
```

Scope:

```text
strict fold-local implementation
ortholog audit
external cohort QC
sensitivity analyses
publication-ready tables and figures
```

## Implementation Risks

Main risks:

```text
Bo2023 count and VSD sample IDs may need careful normalization.
Count matrix has more genes than VSD matrix.
The VSD matrix is batch removed; a simple logCPM-to-VSD projector may not fully capture batch removal.
Human external datasets require ortholog mapping.
Tumor and liquid-biopsy expression distributions are biologically far from macaque brain tissue.
```

Interpretation risks:

```text
Projected values are approximate Bo2023-like VSD, not true native VSD.
External results may reflect domain shift, tumor biology, or biofluid composition rather than anatomical origin.
ComBat or transductive quantile mapping can leak test-set information if not fold-local.
```

## Initial Commands For The Next Thread

Start by inspecting these files:

```powershell
Get-ChildItem "bo2023 data"
Get-Content scripts\reconstruct_bo2023_vsd.R -TotalCount 220
Get-Content scripts\run_bo2023_loso_validation.py -TotalCount 260
Get-Content scripts\run_bo2023_leave_one_monkey_out_validation.py -TotalCount 260
Get-Content core\network_tracing.py -TotalCount 260
```

Then run a small data audit script before building any model.

## Decision Gate

Do not proceed to broad external validation until this gate is met:

```text
Bo2023 projected-VSD LOSO Network Top3 is close to native-VSD LOSO Network Top3,
and projected-VSD is not worse than the simple logCPM/rank baseline.
```

If this gate fails, record the result and keep projected VSD as a negative
sensitivity analysis rather than a model improvement.
