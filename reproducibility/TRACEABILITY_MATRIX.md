# BrainTrace v0.1.17 scientific-provenance patch candidate - Calculation Traceability Matrix

> **Scope**: This matrix traces the current, unreleased v0.1.17 scientific-provenance patch candidate. It does not alter the frozen v0.1.12 scientific model or the immutable v0.1.16 release. Candidate-derived AHBA, TCGA/BraTS, orthology, tier-cascade and sign-flip artifacts are identified explicitly below. References to v0.1.15 GitHub/Zenodo locations in unchanged historical rows identify legacy public mirrors only; they do not describe the current candidate or a new release. No v0.1.17 DOI or release is claimed.

## How to Read This Matrix

Each trace entry has:
- **ID**: Unique identifier (T001–T0xx)
- **Manuscript Location**: Where the number appears (section + paragraph)
- **Manuscript Value**: The exact number as written
- **Source CSV**: Which CSV file contains the underlying data
- **Availability**: Public release or submission package only
- **CSV Location**: Row/column in the CSV
- **Formula**: How the number was computed
- **Raw Data Input**: What raw count or data the computation starts from
- **Verification**: Whether the chain is verified ✓

---

## 1. Internal Validation (LOSO/LOMO)

| ID | Manuscript Location | Manuscript Value | Source CSV | Availability | CSV Location | Formula | Raw Data Input | ✓ |
|----|-------------------|-----------------|------------|--------------|-------------|---------|---------------|---|
| T001 | §4 Validation ¶1 | LOSO Network Top1 = 58.97% (483/819) | `reproducibility/v4_p0_9_triple_ci.csv` | Public - GitHub v0.1.15 release; Zenodo full archive 10.5281/zenodo.21970278 | Row "LOSO Network Top1" | 483/819 = 0.5897 | LOSO validation: 483 correct out of 819 samples | ✓ |
| T002 | §4 Validation ¶1 | LOSO Network Top3 = 91.94% (753/819) | `reproducibility/v4_p0_9_triple_ci.csv` | Public - GitHub v0.1.15 release; Zenodo full archive 10.5281/zenodo.21970278 | Row "LOSO Network Top3" | 753/819 = 0.9194 | LOSO validation: 753 correct out of 819 | ✓ |
| T003 | §4 Validation ¶1 | LOMO Network Top1 = 55.56% (455/819) | `reproducibility/v4_p0_9_triple_ci.csv` | Public - GitHub v0.1.15 release; Zenodo full archive 10.5281/zenodo.21970278 | Row "LOMO Network Top1" | 455/819 = 0.5556 | LOMO validation: 455 correct out of 819 | ✓ |
| T004 | §4 Validation ¶1 | LOMO Network Top3 = 91.58% (750/819) | `reproducibility/v4_p0_9_triple_ci.csv` | Public - GitHub v0.1.15 release; Zenodo full archive 10.5281/zenodo.21970278 | Row "LOMO Network Top3" | 750/819 = 0.9158 | LOMO validation: 750 correct out of 819 | ✓ |
| T005 | §4 Validation ¶1 | Resolution-group Top3 LOSO = 72.48% | `reproducibility/v4_p0_9_triple_ci.csv` | Public - GitHub v0.1.15 release; Zenodo full archive 10.5281/zenodo.21970278 | (derived from RAW_COUNTS_RESOLUTION) | 590/814 = 0.7248 | LOSO ResGroup: 590 correct out of 814 | ✓ |
| T006 | §4 Validation ¶1 | Resolution-group Top3 LOMO = 70.07% | `reproducibility/v4_p0_9_triple_ci.csv` | Public - GitHub v0.1.15 release; Zenodo full archive 10.5281/zenodo.21970278 | (derived from RAW_COUNTS_RESOLUTION) | 569/812 = 0.7007 | LOMO ResGroup: 569 correct out of 812 | ✓ |
| T007 | §4 Validation ¶1 | Exact-region Top3 LOSO = 45.21% | `reproducibility/v4_p0_9_triple_ci.csv` | Public - GitHub v0.1.15 release; Zenodo full archive 10.5281/zenodo.21970278 | Row "LOSO Exact Top3" | 368/814 = 0.4521 | LOSO Exact: 368 correct out of 814 | ✓ |
| T008 | §4 Validation ¶1 | Exact-region Top3 LOMO = 42.61% | `reproducibility/v4_p0_9_triple_ci.csv` | Public - GitHub v0.1.15 release; Zenodo full archive 10.5281/zenodo.21970278 | Row "LOMO Exact Top3" | 346/812 = 0.4261 | LOMO Exact: 346 correct out of 812 | ✓ |
| T009 | §4 Validation ¶1 | Donor-macro Top3: 89.40% vs 89.24% | `reproducibility/v4_p0_13_macro_f1.csv` | Public - GitHub v0.1.15 release; Zenodo full archive 10.5281/zenodo.21970278 | SUMMARY rows (LOSO/LOMO Network) | Donor-weighted mean of per-donor Top3 | Per-donor hit rates from confusion matrices | ✓ |
| T010 | §4 Validation ¶1 | Current four-test sign-flip family: Network Top1 raw/BH=0.031250/0.125000; Network Top3=0.375000/0.500000; resolution-group Top3=0.593750/0.593750; exact-region Top3=0.324219/0.500000; none significant | `reproducibility/sign_flip_current_family.csv` | Unreleased candidate working tree | Rows 1-4 | Exhaustive source raw P values with BH across the four named tests | Donor-level exhaustive 2^9 sign-flip source | ✓ |
| T011 | §4 Validation ¶1 | MDE ~6pp at 80% power | (computed) | Submission package only | — | Minimum detectable effect: 80% power, α=0.05, two-sided, n=9 donors | Donor-level Top3 SD from validation | ✓ |

---

## 2. ML Baselines & Comparators

| ID | Manuscript Location | Manuscript Value | Source CSV | Availability | CSV Location | Formula | Raw Data Input | ✓ |
|----|-------------------|-----------------|------------|--------------|-------------|---------|---------------|---|
| T012 | §4 Validation ¶1 | Historical fold-local k=500 RF comparator: Network Top1/Top3 = 59.95%/91.09% | `reproducibility/v4_p0_11_rf_comparator.csv` | Public - GitHub v0.1.15 release; Zenodo full archive 10.5281/zenodo.21970278 | Rows "Network Top1/Top3" | 491/819=0.5995, 746/819=0.9109 | Historical RF LOMO: 491/746 correct out of 819; not the frozen 200-gene comparator | ✓ |
| T013 | §4 Validation ¶1 | 5-NN cosine Top1 = 43.71% | `reproducibility/v4_p0_11_ml_baselines.csv` | Public - GitHub v0.1.15 release; Zenodo full archive 10.5281/zenodo.21970278 | Row "5-NN cosine" | 358/819=0.4371 | 5-NN LOMO: 358 correct out of 819 | ✓ |
| T014 | §4 Validation ¶1 | Nearest Centroid Top1 = 19.29% | `reproducibility/v4_p0_11_ml_baselines.csv` | Public - GitHub v0.1.15 release; Zenodo full archive 10.5281/zenodo.21970278 | Row "Nearest Centroid" | 158/819=0.1929 | NC LOMO: 158 correct out of 819 | ✓ |
| T015 | §4 Validation ¶1 | Formal route Top1 = 55.56% (n=819) | `reproducibility/v4_p0_11_ml_baselines.csv` | Public - GitHub v0.1.15 release; Zenodo full archive 10.5281/zenodo.21970278 | Row "Formal Route LOMO" | 455/819=0.5556 | LOMO: 455 correct out of 819 | ✓ |

---

## 3. Subcortical PPV/Recall

| ID | Manuscript Location | Manuscript Value | Source CSV | Availability | CSV Location | Formula | Raw Data Input | ✓ |
|----|-------------------|-----------------|------------|--------------|-------------|---------|---------------|---|
| T016 | §4 Validation ¶1 | Subcortical PPV = 1.00 (42/42) | `reproducibility/v4_p0_4_subcortical_subsampling.csv` | Public - GitHub v0.1.15 release; Zenodo full archive 10.5281/zenodo.21970278 | Column "mean_ppv" (all rows) | 42/42 = 1.00 | Subcortical confusion matrix: 42 TP, 0 FP | ✓ |
| T017 | §4 Validation ¶1 | Subcortical recall = 0.78 (42/54) | `reproducibility/v4_p0_4_subcortical_subsampling.csv` | Public - GitHub v0.1.15 release; Zenodo full archive 10.5281/zenodo.21970278 | Full-sample row | 42/54 = 0.7778 | Subcortical confusion matrix: 42 TP, 12 FN | ✓ |

---

## 4. AHBA External Validation

| ID | Manuscript Location | Manuscript Value | Source CSV | Availability | CSV Location | Formula | Raw Data Input | ✓ |
|----|-------------------|-----------------|------------|--------------|-------------|---------|---------------|---|
| T018 | §4 Validation ¶2 | AHBA Network Top1/Top3 = 73.99%/94.62% (165/211 of 223) | `reproducibility/ahba/ahba_endpoint_evaluability_ledger.csv` | Unreleased candidate working tree | Network-evaluable rows | 165/223=0.7399, 211/223=0.9462 | Canonical formal AHBA sample-detail route | ✓ |
| T019 | §4 Validation ¶2 | AHBA ResGroup Top1/Top3 = 42.05%/68.18% (37/60 of 88) | `reproducibility/ahba/ahba_endpoint_evaluability_ledger.csv` | Unreleased candidate working tree | Group-evaluable rows | 37/88=0.4205, 60/88=0.6818 | Canonical formal AHBA sample-detail route | ✓ |
| T020 | §4 Validation ¶2 | AHBA Exact Top1/Top3 = 27.27%/45.45% (24/40 of 88) | `reproducibility/ahba/ahba_endpoint_evaluability_ledger.csv` | Unreleased candidate working tree | Exact-evaluable rows | 24/88=0.2727, 40/88=0.4545 | Canonical formal AHBA sample-detail route | ✓ |
| T021 | §4 Validation ¶2 | 231 replicate-collapsed AHBA tissues | `reproducibility/ahba/ahba_endpoint_evaluability_ledger.csv` | Unreleased candidate working tree | All ledger rows | count(rows)=231 | Canonical formal AHBA sample-detail route | ✓ |
| T022 | §4 Validation ¶2 | Endpoint-specific evaluability: Network n=223; resolution-group/exact n=88 | `reproducibility/ahba/ahba_endpoint_evaluability_ledger.csv` | Unreleased candidate working tree | `network_evaluable`, `group_evaluable`, `exact_evaluable` | sum(endpoint_evaluable) | Endpoint-specific eligibility; not a unique sequential attrition pipeline | ✓ |
| T023 | §4 Validation ¶2 | Single-label sensitivity subsets: Network n=56; group/exact n=40 | `reproducibility/ahba/ahba_endpoint_evaluability_ledger.csv` | Unreleased candidate working tree | `*_single_label_subset` columns | sum(subset flags) | Allowed Network/exact truth-label counts in canonical route | ✓ |
| T024 | §4 Validation ¶2 | Historical AHBA traces are not canonical endpoint accounting | `reproducibility/v4_p0_5_ahba_trace.csv`; `reproducibility/v4_p0_5_ahba_trace_manuscript_aligned.csv` | Historical engineering traces | `trace_classification` column | N/A | **HISTORICAL ENGINEERING TRACE — NOT THE CANONICAL ENDPOINT-EVALUABILITY LEDGER** | ✓ |

---

## 5. TCGA/BraTS Evaluation

| ID | Manuscript Location | Manuscript Value | Source CSV | Availability | CSV Location | Formula | Raw Data Input | ✓ |
|----|-------------------|-----------------|------------|--------------|-------------|---------|---------------|---|
| T025 | §4 Validation ¶2 | 63 primary edema-comparator patients | `reproducibility/tcga_brats_truth_basis_top3_summary.json` | Unreleased candidate working tree | `primary_edema_comparator` | 65 total - TCGA-HT-7686 (no label-2 edema voxels) - TCGA-HT-7680 (cerebellar/posterior-fossa out of scope) = 63; center/core/whole_tumor retain n=65 | Canonical per-patient truth/prediction output | ✓ |
| T026 | §4 Validation ¶2 | Network Top3 = 23.81% | `reproducibility/tcga_brats_truth_basis_top3_summary.csv` | Unreleased candidate working tree | edema/network/top3/strict | 15/63 = 0.2381 | TCGA/BraTS: 15 correct out of 63 | ✓ |
| T027 | §4 Validation ¶2 | Broad Top3 = 82.54% | `reproducibility/tcga_brats_truth_basis_top3_summary.csv` | Unreleased candidate working tree | edema/broad/top3/strict | 52/63 = 0.8254 | TCGA/BraTS: 52 correct out of 63 | ✓ |
| T028 | §4 Validation ¶2 | p=0.8888 (exploratory one-sided exact-binomial vs 30%) | `reproducibility/v4_p0_12_tcga_brats_ci_summary.csv` | Public - GitHub v0.1.15 release; Zenodo full archive 10.5281/zenodo.21970278 | manuscript_summary row | binomtest(15, 63, 0.30, "greater").pvalue | 15/63 vs null_p=0.30 (Top3/10 uniform) | ✓ |
| T029 | §4 Validation ¶2 | Strict truth-basis Top3 ranges: Network 15.38–29.23%=13.85 pp; broad 49.23–82.54%=33.31 pp | `reproducibility/tcga_brats_truth_basis_top3_summary.csv` | Unreleased candidate working tree | all 8 strict Top3 rows | max(percent)-min(percent) within level | Network: 19/65, 12/65, 15/63, 10/65; broad: 32/65, 45/65, 52/63, 46/65 | ✓ |

Candidate-set sensitivity: Network Top3 any-hit = 36.51% (23/63), descriptive only; the patient-specific truth set contains Networks with at least 20% edema overlap and is not an anatomical adjacency rule.

Endpoint provenance: the candidate values are recomputed from the canonical per-patient truth/prediction output into `reproducibility/tcga_brats_truth_basis_top3_summary.csv`. The primary edema comparator additionally excludes the no-edema and cerebellar/out-of-scope cases listed in T025. The archived 2026-06-09 endpoint is historical and is not a source for T026–T029.

---

## 6. Lambda Sensitivity & Friedman Test

| ID | Manuscript Location | Manuscript Value | Source CSV | Availability | CSV Location | Formula | Raw Data Input | ✓ |
|----|-------------------|-----------------|------------|--------------|-------------|---------|---------------|---|
| T030 | §2 System ¶4 | Network Top3 remained within 91.94% (<1pp deviation) | `reproducibility/v4_p0_10_lambda_sensitivity.csv` | Public - GitHub v0.1.15 release; Zenodo full archive 10.5281/zenodo.21970278 | All rows, network_top3 | 0.9194 at λ=0.25, 0.50, 0.75 (identical) | Lambda only affects exact-region fusion, not Network | ✓ |
| T031 | §2 System ¶4 | Friedman chi2=0.54, df=2, p=0.764 | `reproducibility/v4_p0_10_lambda_friedman.csv` | Public - GitHub v0.1.15 release; Zenodo full archive 10.5281/zenodo.21970278 | Friedman_test row | friedmanchisquare(hit3_λ0.25, hit3_λ0.50, hit3_λ0.75) | 9-donor per-donor hit3 at 3 lambda values | ✓ |
| T032 | Supp SR4 ¶3 | Friedman hit1 chi2=3.0, p=0.2231 | `reproducibility/v4_p0_10_lambda_friedman.csv` | Public - GitHub v0.1.15 release; Zenodo full archive 10.5281/zenodo.21970278 | Friedman_test_hit1 row | friedmanchisquare(hit1_λ0.25, hit1_λ0.50, hit1_λ0.75) | 9-donor per-donor hit1 at 3 lambda values | ✓ |

---

## 7. Engineering Performance

| ID | Manuscript Location | Manuscript Value | Source CSV | Availability | CSV Location | Formula | Raw Data Input | ✓ |
|----|-------------------|-----------------|------------|--------------|-------------|---------|---------------|---|
| T033 | §3 Implementation ¶1 | Cold inference 0.3841 s/profile (19.5901 s / 51 profiles) | `reproducibility/formal_real_input_performance_manifest.json` | Unreleased candidate working tree | cold/timing | 19.5901258 s / 51 profiles | GSE189919: 51-profile cold frozen-route workload | ✓ |
| T034 | §3 Implementation ¶1 | Cold p50/p95 0.3769/0.4026 s/profile | `reproducibility/formal_real_input_performance_manifest.json` | Unreleased candidate working tree | cold/timing | Percentile of 51 cold profile timings | GSE189919 cold per-profile timing | ✓ |
| T035 | §3 Implementation ¶1 | Cold peak 216.7 MiB; warm maximum 222.0 MiB across 153 timed inference events | `reproducibility/formal_real_input_performance_manifest.json` | Unreleased candidate working tree | cold/memory; warm/repeats | bytes / 2^20; 51 profiles × 3 warm repeats = 153 events | GSE189919 benchmark memory profiling | ✓ |

---

## 8. Cross-Species / Humanization

| ID | Manuscript Location | Manuscript Value | Source CSV | Availability | CSV Location | Formula | Raw Data Input | ✓ |
|----|-------------------|-----------------|------------|--------------|-------------|---------|---------------|---|
| T036 | §2 System ¶3 | 5,324/8,800 humanized and 3,476/8,800 unmapped **gene-by-region row occurrences** | `reproducibility/orthology_humanization_summary.json` | Unreleased candidate working tree | `region_signature_rows` | 5324/8800=60.50%; 3476/8800=39.50% | Frozen humanization output | ✓ |
| T037 | §2 System ¶3 | Top200 orthology humanizable = 188/200; g:Profiler mapped = 179/200 | `reproducibility/orthology_humanization_summary.json` | Unreleased candidate working tree | `network_top200_orthology_humanizable`; `gprofiler_mapped` | Separate numerator/denominator/filtering universes | Frozen orthology output and g:Profiler result metadata | ✓ |
| T038 | §2 System ¶3 | Row-level humanization is not independent-gene humanization | `reproducibility/orthology_humanization_summary.json` | Unreleased candidate working tree | `region_signature_rows.unit` | Unit is a gene-by-region row occurrence | 8,800 is not a count of independent macaque genes | ✓ |

---

## 9. Sparse-Query Sensitivity

| ID | Manuscript Location | Manuscript Value | Source CSV | Availability | CSV Location | Formula | Raw Data Input | ✓ |
|----|-------------------|-----------------|------------|--------------|-------------|---------|---------------|---|
| T039 | §5 Use ¶3 | Network Top3 declined 91.94% → 58.54% | (archived sparse simulation) | Submission package only | — | 30-repeat mean at extreme sparsity (20% genes, 1% depth) | 30-repeat sparse simulation, seed 20260711 | ✓ |
| T040 | §2 System ¶5 | 30 repeats each; 99,099 archived seeds | `seed_registry.json` | Public - GitHub v0.1.15 release; Zenodo full archive 10.5281/zenodo.21970278 | — | Base seed 20260711 + index → 99,099 unique seeds | Sparse simulation seed registry | ✓ |
| T041 | §2 System ¶5 | 50,000 bootstrap draws (seed 20260716) | `seed_registry.json` | Public - GitHub v0.1.15 release; Zenodo full archive 10.5281/zenodo.21970278 | — | 50,000 donor-cluster bootstrap resamples | Seed registry: donor_cluster_bootstrap=20260716 | ✓ |

---

## 10. GSE189919 & Huang2025

| ID | Manuscript Location | Manuscript Value | Source CSV | Availability | CSV Location | Formula | Raw Data Input | ✓ |
|----|-------------------|-----------------|------------|--------------|-------------|---------|---------------|---|
| T042 | §4 Validation ¶2 | GSE189919: 51 samples, 72.10% model-space gene overlap | (archived benchmark) | Submission package only | — | 15,622/21,668 = 0.7210 | GSE189919 gene symbols vs frozen Bo2023 projector gene space | ✓ |
| T043 | §4 Validation ¶2; Huang remediation working manuscript/supplement | Huang2025 full published matrix: 159/159 traceable profiles (77 CSF, 82 plasma); minimum nominal BH-adjusted profile-level P=0.722052 | `huang_2025_canonical_summary.json`; `huang_2025_sample_ledger.csv`; `huang_2025_tumour_control_comparisons.csv` | Unreleased working remediation; scientific approval required before publication | `reproducibility/huang_2025/` | 77 + 82 = 159; minimum nominal BH-adjusted value across six exploratory two-sided profile-level Mann-Whitney U comparisons; patient pairing unavailable and patient-level dependence not assessable | Huang2025 published cfRNA expression matrix; locked BrainTrace route | ✓ |

---

## 11. Confidence Intervals (Wilson/CP/AC)

All CI values in the manuscript and supplementary tables are traced to `v4_p0_9_triple_ci.csv` and `v4_p0_12_tcga_brats_ci_summary.csv`.

| CI Type | Formula | Reference | Used In |
|---------|---------|-----------|---------|
| Wilson score | p̃ = (p + z²/2n) / (1 + z²/n) ± z√[(p(1-p) + z²/4n)/n] / (1 + z²/n) | Wilson (1927) | Table S2, S4a, S13 |
| Clopper-Pearson | Beta⁻¹(α/2; x, n-x+1) to Beta⁻¹(1-α/2; x+1, n-x) | Clopper & Pearson (1934) | Table S2, S4a, S13 |
| Agresti-Coull | p̃ = (x+z²/2)/(n+z²) ± z√[p̃(1-p̃)/(n+z²)] | Agresti & Coull (1998) | Table S2, S4a |
| Bootstrap (percentile) | 2.5th–97.5th percentile of 50,000 resampled proportions | Efron (1979) | Table S13, TCGA/BraTS |
| Binomial test | P(X ≥ k \| X ~ Bin(n, p₀)) | — | TCGA/BraTS p-values |

All formulas are implemented in `reproducibility/generate_all_csvs.py` with documentation.

---

## 12. Macro F1 Summary Statistics

| ID | Manuscript Location | Manuscript Value | Source CSV | Availability | CSV Location | Formula | Raw Data Input | ✓ |
|----|-------------------|-----------------|------------|--------------|-------------|---------|---------------|---|
| T044 | Supp Table S8/S13 | LOSO Exact macro F1 = 0.2102 | `reproducibility/v4_p0_13_macro_f1.csv` | Public - GitHub v0.1.15 release; Zenodo full archive 10.5281/zenodo.21970278 | SUMMARY, LOSO_Exact_macro | mean(F1_i) for i=1..105 | 105 class-level F1 from LOSO Exact confusion matrix | ✓ |
| T045 | Supp Table S8/S13 | LOMO Exact macro F1 = 0.2034 | `reproducibility/formal_lomo_exact_region_f1.csv`; `reproducibility/v4_p0_13_macro_f1.csv` | Unreleased candidate working tree | Current formal class rows; SUMMARY, LOMO_Exact_macro | mean(F1_i) over the 104-label truth universe | 812 current formal prediction rows; prediction-only labels remain false positives | ✓ |
| T046 | Supp Table S8 | LOSO Network macro F1 | `reproducibility/v4_p0_13_macro_f1.csv` | Public - GitHub v0.1.15 release; Zenodo full archive 10.5281/zenodo.21970278 | SUMMARY, LOSO_Network_macro | mean(F1_i) for i=1..10 | 10 Network class-level F1 from LOSO confusion matrix | ✓ |
| T047 | Supp Table S8 | LOMO Network weighted F1 | `reproducibility/v4_p0_13_macro_f1.csv` | Public - GitHub v0.1.15 release; Zenodo full archive 10.5281/zenodo.21970278 | SUMMARY, LOMO_Network_weighted | Σ(n_i·F1_i) / Σ(n_i) for i=1..10 | 10 Network classes with sample counts from LOMO | ✓ |

---

**T055 (updated 2026-08-19):** Supplementary S11 reports Macro-F1 mean±SD, median(IQR), weighted-F1 and micro-F1 for LOSO/LOMO Exact and Network endpoints. The current LOMO Exact rows are regenerated from `reproducibility/p2_publication_completeness/formal_lomo_exact_region_detail.csv` by `scripts/generate_lomo_exact_f1_evidence.py`; the 104-class macro denominator is the truth-label universe, and `sum(TP)=177`, `sum(support)=812`, and micro-F1=`177/812`. The derived summaries are in `reproducibility/p1_cross1_5/cross3_f1_distribution_summary.csv` and `reproducibility/v4_p0_13_macro_f1.csv`.

**T056 (added 2026-08-10):** Supplementary S6 now reports all ten Network mapping fractions, AHBA allowed-label denominators, any-allowed Top1 hits, strict same-Network Top1 hits and percentage-point effects versus the pooled 165/223 mapped-label Top1 rate. The public derived source is `reproducibility/p1_e4_ahba_network_mapping_impact.csv`; its input hashes and formula are in `p1_e4_ahba_network_mapping_impact_manifest.json`. The primary denominators are membership-weighted because AHBA truth is multi-label; all results are descriptive with two AHBA donors.

**T057 (updated 2026-08-13):** CI runs Coverage.py branch coverage for `core/` and uploads `coverage.xml`. The v0.1.15 release-candidate run is 110 passed, 2 skipped and 64% total coverage (1,117 statements; 298 branches); the text and XML reports are `reproducibility/coverage_report.txt` and `reproducibility/coverage.xml`.

**T058 (added 2026-08-10):** Bo2023 provenance now records publisher file-level Supplementary Data 1–15 and Source Data URLs, the SRA BioProject and the authors' Zenodo code archive. The exact SHA-matched processed count/VSD author-package files remain non-redistributed; no stable public direct URL for those two matrices was identified, and the README now states this explicitly.

**T059 (updated 2026-08-17):** Formal LOMO Network F1 is regenerated from the
819-row prediction-level source
`reproducibility/p2_publication_completeness/formal_lomo_network_detail.csv`,
filtered to route family `hybrid_projected_network_logcpm_exact` and route
`network_discriminative_correlation_top200`. The locked endpoint is 455/819
Top1 and 750/819 Top3; integer class TP counts sum to 455, so micro-F1 is
455/819 = 0.5555555556. `formal_lomo_network_f1.csv` contains the class-level
TP/FP/FN/P/R/F1 rows, and `lomo_network_f1_provenance.json` plus
`LOMO_NETWORK_F1_PROVENANCE.md` record the source digest and superseded-route
root cause.

**T060 (updated 2026-08-17):** The derived LOMO Network distribution summary
in `reproducibility/p1_cross1_5/cross3_f1_distribution_summary.csv` and the
LOMO rows in `reproducibility/v4_p0_13_macro_f1.csv` are regenerated from the
same prediction-level source by `scripts/generate_lomo_f1_evidence.py` and
are regression-tested against the integer-count class table.

**T062 (added 2026-08-19):** Table S8 current-formal resolution-group Top3
random baselines are regenerated from the staged LOSO and LOMO prediction
details by `core/resolution_group_baselines.py`: LOSO uniform/weighted =
22.6%/6.5%; LOMO uniform/weighted = 21.3%/4.1%. The input origin and staged
SHA-256 values, formula, RNG seed and 10,000 weighted draws are recorded in
`reproducibility/formal_resolution_group_random_baselines.json`.

## 13. Frozen 200-gene Comparator, Truth Normalization, and Independent Analyses

| ID | Manuscript Location | Manuscript Value | Source CSV / Manifest | Availability | CSV Location | Formula | Raw Data Input | ✓? |
|----|-------------------|-----------------|-----------------------|--------------|-------------|---------|---------------|---|
| T048 | Supplementary S18 | Frozen 200-gene RF LOMO Top1/Top3 = 47.50%/83.03% (389/680 of 819) | `reproducibility/p2_publication_completeness/P2_RF200_lomo_detail.csv` | Public - GitHub v0.1.15 release; Zenodo full archive 10.5281/zenodo.21970278 | `hit1`, `hit3`, all 819 rows | Σhit1/819 = 389/819; Σhit3/819 = 680/819 | Exact frozen 200-gene Network panel, leave-one-macaque-out predictions | ✓ |
| T049 | Supp S18 / truth-normalization note | Two prespecified frozen-truth corrections | Same CSV as T048 plus frozen mapping audit | Submission package only | Samples `19R470` and `xz-L-10m-10o_BRRL200004691-1A` | Apply mapping before scoring: Parietal → Occipital/Temporal; Lateral Prefrontal Cortex → OMPFC | Raw metadata, anatomical label crosswalk, and frozen scoring truth | ✓ |
| T050 | Supp 3′-bias sensitivity | 185 evaluable genes; Spearman ρ=0.646, p=2.95×10⁻²³; median detected 151/185 | `reproducibility/p0_bio2_3prime_bias/panel_length_detection.csv` and archived summary | Submission package only | Per-gene detection/length rows | Spearman correlation and per-sample median; terminal-window simulations at 500/1,000/2,000 nt yield medians 74/106/136 | Observed gene detection plus transcript-coordinate proxy simulation; not molecule-level 3′-end sequencing simulation | ✓ |
| T051 | Supp independent panel annotation | 179/200 mapped; significant GO:BP/KEGG terms = 446/11 (primary) and 410/8 (annotated-domain sensitivity) | `reproducibility/independent_enrichment/independent_enrichment_manifest.json` | Public - GitHub v0.1.15 release; Zenodo full archive 10.5281/zenodo.21970278 | `gprofiler_*_significant` and g:Profiler metadata | g:Profiler FDR <0.05 separately by source | Frozen 200-gene panel, 21,668-gene model-space background, g:Profiler e114_eg62_p19_27110d83 (2026-07-31) | ✓ |
| T052 | Supp independent cell-type annotation | Rhesus excitatory/inhibitory enrichment q=2.67×10⁻⁷/4.20×10⁻⁴; human sensitivity q=0.0477/5.89×10⁻⁵ | `reproducibility/independent_enrichment/independent_celltype_enrichment.csv`; `reproducibility/independent_enrichment/independent_enrichment_manifest.json` | Public - GitHub v0.1.15 release; Zenodo full archive 10.5281/zenodo.21970278 | Chiou2023 and Siletti2023 rows | One-sided hypergeometric test; BH correction across seven prespecified families | Chiou 2023 rhesus markers (primary) and Siletti 2023 human markers (sensitivity), intersected with 21,668-gene background | ✓ |
| T053 | Supplementary S18 | Complete 110-row Bo2023-to-Saleem anatomical crosswalk | `reproducibility/crosswalks/P2_Bo2023_Saleem_crosswalk.csv` | Public - GitHub v0.1.15 release; Zenodo full archive 10.5281/zenodo.21970278 | All 110 rows | Direct row-level crosswalk | Locked Bo2023 region IDs, Network labels, Saleem-style names and broad mappings | ✓ |
| T054 | Supplementary Tables S18-S19 | Portable sensitivity blocks and frozen output | `reproducibility/s18_s19/run_s18_s19_sensitivity.py`; `reproducibility/s18_s19/sensitivity_analysis_results.json`; `reproducibility/s18_s19/input_manifest.json` | Public - GitHub v0.1.15 release; Zenodo full archive 10.5281/zenodo.21970278 | Four named result blocks and input hashes | Deterministic runner using frozen public inputs | Frozen model artifacts and corrected lambda/AHBA CSV inputs | ✓ |

## 13.1 Current Candidate Tier-Cascade Provenance

| ID | Manuscript Location | Manuscript Value | Source CSV / Manifest | Availability | CSV Location | Formula | Raw Data Input | ✓? |
|----|-------------------|-----------------|-----------------------|--------------|-------------|---------|---------------|---|
| T061 | Supplementary Results R10 | Same-LOSO-exact-universe cascade: exact n=814, Exact Top3=368/814, Network candidate truth retained/missed=750/64, Network miss share of exact misses=64/446=14.35%, conditional Exact/Group Top3=368/750=49.07% and 590/750=78.67%, recovery after a Network candidate miss=0% | `reproducibility/tier_cascade_loso_summary.json` | Unreleased candidate working tree | top-level fields | All candidate-set quantities use the exact-evaluable 814-sample universe | Canonical LOSO exact and resolution-group sample-detail outputs | ✓ |

---

## Verification Summary

| Category | Trace IDs | Verified | Coverage |
|----------|-----------|----------|----------|
| Internal validation (LOSO/LOMO) | T001–T011 | 11/11 | 100% |
| ML baselines & comparators | T012–T015 | 4/4 | 100% |
| Subcortical PPV/recall | T016–T017 | 2/2 | 100% |
| AHBA external validation | T018–T024 | 7/7 | 100% |
| TCGA/BraTS evaluation | T025–T029 | 5/5 | 100% |
| Lambda sensitivity & Friedman | T030–T032 | 3/3 | 100% |
| Engineering performance | T033–T035 | 3/3 | 100% |
| Cross-species/humanization | T036–T038 | 3/3 | 100% |
| Sparse-query sensitivity | T039–T041 | 3/3 | 100% |
| GSE189919 & Huang2025 | T042–T043 | 2/2 | 100% |
| Macro F1 summary | T044–T047, T055 | 5/5 | 100% |
| Mapping-impact and CI coverage addenda | T056–T058 | 3/3 | 100% |
| Frozen comparator, truth normalization & independent analyses | T048-T054 | 7/7 | 100% |
| Current candidate tier-cascade provenance | T061 | 1/1 | 100% |
| Current-formal resolution-group baseline provenance | T062 | 1/1 | 100% |
| **Total** | **T001-T062** | **62/62** | **100%** |

> **Availability note:** “Unreleased candidate working tree” means the artifact is present only in this v0.1.17 scientific-provenance patch candidate; it is not a published release and has no new DOI. References to GitHub/Zenodo v0.1.15 in unchanged rows are explicitly historical mirrors. The current immutable public release remains v0.1.16; submission-package-only status does not imply public availability.
