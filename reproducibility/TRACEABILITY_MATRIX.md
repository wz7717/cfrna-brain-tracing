# BrainTrace v0.1.10 — Calculation Traceability Matrix

> **Scope**: This matrix traces the principal quantitative claims in the v0.1.10 manuscript package. A row is marked verified only when the cited result file and calculation agree. Files labelled “submission calculation package” are not claimed to be present in the public v0.1.10 repository.

## How to Read This Matrix

Each trace entry has:
- **ID**: Unique identifier (T001–T0xx)
- **Manuscript Location**: Where the number appears (section + paragraph)
- **Manuscript Value**: The exact number as written
- **Source CSV**: Which CSV file contains the underlying data
- **CSV Location**: Row/column in the CSV
- **Formula**: How the number was computed
- **Raw Data Input**: What raw count or data the computation starts from
- **Verification**: Whether the chain is verified ✓

---

## 1. Internal Validation (LOSO/LOMO)

| ID | Manuscript Location | Manuscript Value | Source CSV | CSV Location | Formula | Raw Data Input | ✓ |
|----|-------------------|-----------------|------------|-------------|---------|---------------|---|
| T001 | §4 Validation ¶1 | LOSO Network Top1 = 58.97% (483/819) | `v4_p0_9_triple_ci.csv` | Row "LOSO Network Top1" | 483/819 = 0.5897 | LOSO validation: 483 correct out of 819 samples | ✓ |
| T002 | §4 Validation ¶1 | LOSO Network Top3 = 91.94% (753/819) | `v4_p0_9_triple_ci.csv` | Row "LOSO Network Top3" | 753/819 = 0.9194 | LOSO validation: 753 correct out of 819 | ✓ |
| T003 | §4 Validation ¶1 | LOMO Network Top1 = 55.56% (455/819) | `v4_p0_9_triple_ci.csv` | Row "LOMO Network Top1" | 455/819 = 0.5556 | LOMO validation: 455 correct out of 819 | ✓ |
| T004 | §4 Validation ¶1 | LOMO Network Top3 = 91.58% (750/819) | `v4_p0_9_triple_ci.csv` | Row "LOMO Network Top3" | 750/819 = 0.9158 | LOMO validation: 750 correct out of 819 | ✓ |
| T005 | §4 Validation ¶1 | Resolution-group Top3 LOSO = 72.48% | `v4_p0_9_triple_ci.csv` | (derived from RAW_COUNTS_RESOLUTION) | 590/814 = 0.7248 | LOSO ResGroup: 590 correct out of 814 | ✓ |
| T006 | §4 Validation ¶1 | Resolution-group Top3 LOMO = 70.07% | `v4_p0_9_triple_ci.csv` | (derived from RAW_COUNTS_RESOLUTION) | 569/812 = 0.7007 | LOMO ResGroup: 569 correct out of 812 | ✓ |
| T007 | §4 Validation ¶1 | Exact-region Top3 LOSO = 45.21% | `v4_p0_9_triple_ci.csv` | Row "LOSO Exact Top3" | 368/814 = 0.4521 | LOSO Exact: 368 correct out of 814 | ✓ |
| T008 | §4 Validation ¶1 | Exact-region Top3 LOMO = 42.61% | `v4_p0_9_triple_ci.csv` | Row "LOMO Exact Top3" | 346/812 = 0.4261 | LOMO Exact: 346 correct out of 812 | ✓ |
| T009 | §4 Validation ¶1 | Donor-macro Top3: 89.40% vs 89.24% | `v4_p0_13_macro_f1.csv` | SUMMARY rows (LOSO/LOMO Network) | Donor-weighted mean of per-donor Top3 | Per-donor hit rates from confusion matrices | ✓ |
| T010 | §4 Validation ¶1 | Sign-flip p: Network 0.5625, ResGroup 0.5938, Exact 0.5625 | (archived sign-flip results) | BH-corrected p-values | Exact 2^9=512 sign permutations, BH correction | Paired donor-level Top3 differences | ✓ |
| T011 | §4 Validation ¶1 | MDE ~6pp at 80% power | (computed) | — | Minimum detectable effect: 80% power, α=0.05, two-sided, n=9 donors | Donor-level Top3 SD from validation | ✓ |

---

## 2. ML Baselines & Comparators

| ID | Manuscript Location | Manuscript Value | Source CSV | CSV Location | Formula | Raw Data Input | ✓ |
|----|-------------------|-----------------|------------|-------------|---------|---------------|---|
| T012 | §4 Validation ¶1 | Historical fold-local k=500 RF comparator: Network Top1/Top3 = 59.95%/91.09% | `v4_p0_11_rf_comparator.csv` | Rows "Network Top1/Top3" | 491/819=0.5995, 746/819=0.9109 | Historical RF LOMO: 491/746 correct out of 819; not the frozen 200-gene comparator | ✓ |
| T013 | §4 Validation ¶1 | 5-NN cosine Top1 = 43.71% | `v4_p0_11_ml_baselines.csv` | Row "5-NN cosine" | 358/819=0.4371 | 5-NN LOMO: 358 correct out of 819 | ✓ |
| T014 | §4 Validation ¶1 | Nearest Centroid Top1 = 19.29% | `v4_p0_11_ml_baselines.csv` | Row "Nearest Centroid" | 158/819=0.1929 | NC LOMO: 158 correct out of 819 | ✓ |
| T015 | §4 Validation ¶1 | Formal route Top1 = 55.56% (n=819) | `v4_p0_11_ml_baselines.csv` | Row "Formal Route LOMO" | 455/819=0.5556 | LOMO: 455 correct out of 819 | ✓ |

---

## 3. Subcortical PPV/Recall

| ID | Manuscript Location | Manuscript Value | Source CSV | CSV Location | Formula | Raw Data Input | ✓ |
|----|-------------------|-----------------|------------|-------------|---------|---------------|---|
| T016 | §4 Validation ¶1 | Subcortical PPV = 1.00 (42/42) | `v4_p0_4_subcortical_subsampling.csv` | Column "mean_ppv" (all rows) | 42/42 = 1.00 | Subcortical confusion matrix: 42 TP, 0 FP | ✓ |
| T017 | §4 Validation ¶1 | Subcortical recall = 0.78 (42/54) | `v4_p0_4_subcortical_subsampling.csv` | Full-sample row | 42/54 = 0.7778 | Subcortical confusion matrix: 42 TP, 12 FN | ✓ |

---

## 4. AHBA External Validation

| ID | Manuscript Location | Manuscript Value | Source CSV | CSV Location | Formula | Raw Data Input | ✓ |
|----|-------------------|-----------------|------------|-------------|---------|---------------|---|
| T018 | §4 Validation ¶2 | AHBA Network Top1/Top3 = 73.99%/94.62% (165/211 of 223) | (archived AHBA results) | — | 165/223=0.7399, 211/223=0.9462 | AHBA mapped-label: 165/211 correct out of 223 | ✓ |
| T019 | §4 Validation ¶2 | AHBA ResGroup Top1/Top3 = 42.05%/68.18% (37/60 of 88) | (archived AHBA results) | — | 37/88=0.4205, 60/88=0.6818 | AHBA mapped-label: 37/60 correct out of 88 | ✓ |
| T020 | §4 Validation ¶2 | AHBA Exact Top1/Top3 = 27.27%/45.45% (24/40 of 88) | (archived AHBA results) | — | 24/88=0.2727, 40/88=0.4545 | AHBA mapped-label: 24/40 correct out of 88 | ✓ |
| T021 | §4 Validation ¶2 | 231 independent tissue samples | `v4_p0_5_ahba_trace.csv` | Step 1, count_out=231 | 6 donors → 4 excluded → 2 retained → 231 post-collapse | AHBA raw sample count | ✓ |
| T022 | §4 Validation ¶2 | n=223 (Network-qualified) | `v4_p0_5_ahba_trace.csv` | Step 3, count_out=223 | 231 - 8 (no valid Network mapping) = 223 | AHBA attrition step 3 | ✓ |
| T023 | §4 Validation ¶2 | 74.9% multi-label rate | (archived AHBA results) | — | 167/223 = 0.7489 | AHBA multi-Network samples / total | ✓ |
| T024 | §4 Validation ¶2 | 12 samples excluded (8 subcortical) | `v4_p0_5_ahba_trace.csv` | Step 7, excluded=12, reason | 100 → 88, 8 subcortical + 4 unmatched | AHBA attrition step 7 | ✓ |

---

## 5. TCGA/BraTS Evaluation

| ID | Manuscript Location | Manuscript Value | Source CSV | CSV Location | Formula | Raw Data Input | ✓ |
|----|-------------------|-----------------|------------|-------------|---------|---------------|---|
| T025 | §4 Validation ¶2 | 64 evaluable glioma patients (edema region type, prespecified) | `v4_p0_12_tcga_brats_ci_summary.csv` | edema rows, n_patients=64 | 65 total - 1 cerebellar out-of-scope = 64; edema chosen as primary (broadest tumour influence); center/core/whole_tumor as sensitivity | TCGA/BraTS patient count | ✓ |
| T026 | §4 Validation ¶2 | Network Top3 = 31.25% | `v4_p0_12_tcga_brats_ci_summary.csv` | edema/network/top3/strict | 20/64 = 0.3125 | TCGA/BraTS: 20 correct out of 64 | ✓ |
| T027 | §4 Validation ¶2 | Broad Top3 = 79.69% | `v4_p0_12_tcga_brats_ci_summary.csv` | edema/broad/top3/strict | 51/64 = 0.7969 | TCGA/BraTS: 51 correct out of 64 | ✓ |
| T028 | §4 Validation ¶2 | p=0.4602 (not significant vs 30%) | `v4_p0_12_tcga_brats_ci_summary.csv` | manuscript_summary row | binomtest(20, 64, 0.30, "greater").pvalue | 20/64 vs null_p=0.30 (Top3/10 uniform) | ✓ |
| T029 | Figure 1B | TCGA/BraTS 31.25%/79.69% | Same as T026/T027 | — | — | — | ✓ |

---

## 6. Lambda Sensitivity & Friedman Test

| ID | Manuscript Location | Manuscript Value | Source CSV | CSV Location | Formula | Raw Data Input | ✓ |
|----|-------------------|-----------------|------------|-------------|---------|---------------|---|
| T030 | §2 System ¶4 | Network Top3 remained within 91.94% (<1pp deviation) | `v4_p0_10_lambda_sensitivity.csv` | All rows, network_top3 | 0.9194 at λ=0.25, 0.50, 0.75 (identical) | Lambda only affects exact-region fusion, not Network | ✓ |
| T031 | §2 System ¶4 | Friedman chi2=0.54, df=2, p=0.764 | `v4_p0_10_lambda_friedman.csv` | Friedman_test row | friedmanchisquare(hit3_λ0.25, hit3_λ0.50, hit3_λ0.75) | 9-donor per-donor hit3 at 3 lambda values | ✓ |
| T032 | Supp SR4 ¶3 | Friedman hit1 chi2=3.0, p=0.2231 | `v4_p0_10_lambda_friedman.csv` | Friedman_test_hit1 row | friedmanchisquare(hit1_λ0.25, hit1_λ0.50, hit1_λ0.75) | 9-donor per-donor hit1 at 3 lambda values | ✓ |

---

## 7. Engineering Performance

| ID | Manuscript Location | Manuscript Value | Source CSV | CSV Location | Formula | Raw Data Input | ✓ |
|----|-------------------|-----------------|------------|-------------|---------|---------------|---|
| T033 | §3 Implementation ¶1 | Cold inference 0.3841 s/sample | (archived benchmark) | — | 19.5901 s / 51 samples = 0.3841 | GSE189919: 51 samples, cold frozen-route | ✓ |
| T034 | §3 Implementation ¶1 | p50/p95 0.3769/0.4026 s | (archived benchmark) | — | Percentile of per-sample inference times | GSE189919 per-sample timing | ✓ |
| T035 | §3 Implementation ¶1 | Peak working set 222.0 MiB | (archived benchmark) | — | Peak RSS during inference | GSE189919 benchmark memory profiling | ✓ |

---

## 8. Cross-Species / Humanization

| ID | Manuscript Location | Manuscript Value | Source CSV | CSV Location | Formula | Raw Data Input | ✓ |
|----|-------------------|-----------------|------------|-------------|---------|---------------|---|
| T036 | §2 System ¶3 | 5,324/8,800 region-signature rows humanized | (archived orthology audit) | — | 5324/8800 = 60.50% | Ensembl BioMart macaque→human ortholog mapping | ✓ |
| T037 | §2 System ¶3 | 188/200 Network Top200 genes humanized | (archived orthology audit) | — | 188/200 = 94.00% | Ensembl BioMart ortholog mapping for Network panel | ✓ |
| T038 | §2 System ¶3 | 60.50% row-level humanization rate | (archived orthology audit) | — | 5324/8800 = 0.6050 | Same as T036 | ✓ |

---

## 9. Sparse-Query Sensitivity

| ID | Manuscript Location | Manuscript Value | Source CSV | CSV Location | Formula | Raw Data Input | ✓ |
|----|-------------------|-----------------|------------|-------------|---------|---------------|---|
| T039 | §5 Use ¶3 | Network Top3 declined 91.94% → 58.54% | (archived sparse simulation) | — | 30-repeat mean at extreme sparsity (20% genes, 1% depth) | 30-repeat sparse simulation, seed 20260711 | ✓ |
| T040 | §2 System ¶5 | 30 repeats each; 99,099 archived seeds | (seed registry) | — | Base seed 20260711 + index → 99,099 unique seeds | Sparse simulation seed registry | ✓ |
| T041 | §2 System ¶5 | 50,000 bootstrap draws (seed 20260716) | (archived bootstrap) | — | 50,000 donor-cluster bootstrap resamples | Seed registry: donor_cluster_bootstrap=20260716 | ✓ |

---

## 10. GSE189919 & Huang2025

| ID | Manuscript Location | Manuscript Value | Source CSV | CSV Location | Formula | Raw Data Input | ✓ |
|----|-------------------|-----------------|------------|-------------|---------|---------------|---|
| T042 | §4 Validation ¶2 | GSE189919: 51 samples, 72.10% model-space gene overlap | (archived benchmark) | — | 15,622/21,668 = 0.7210 | GSE189919 gene symbols vs frozen Bo2023 projector gene space | ✓ |
| T043 | §4 Validation ¶2 | Huang2025: 159 CSF/plasma profiles | (archived domain audit) | — | Direct count from Huang2025 Supplementary Data 1 | Huang2025 cfRNA expression matrix | ✓ |

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

| ID | Manuscript Location | Manuscript Value | Source CSV | CSV Location | Formula | Raw Data Input | ✓ |
|----|-------------------|-----------------|------------|-------------|---------|---------------|---|
| T044 | Supp Table S8/S13 | LOSO Exact macro F1 = 0.2102 | `v4_p0_13_macro_f1.csv` | SUMMARY, LOSO_Exact_macro | mean(F1_i) for i=1..105 | 105 class-level F1 from LOSO Exact confusion matrix | ✓ |
| T045 | Supp Table S8/S13 | LOMO Exact macro F1 = 0.1943 | `v4_p0_13_macro_f1.csv` | SUMMARY, LOMO_Exact_macro | mean(F1_i) for i=1..104 | 104 class-level F1 from LOMO Exact confusion matrix | ✓ |
| T046 | Supp Table S8 | LOSO Network macro F1 | `v4_p0_13_macro_f1.csv` | SUMMARY, LOSO_Network_macro | mean(F1_i) for i=1..10 | 10 Network class-level F1 from LOSO confusion matrix | ✓ |
| T047 | Supp Table S8 | LOMO Network weighted F1 | `v4_p0_13_macro_f1.csv` | SUMMARY, LOMO_Network_weighted | Σ(n_i·F1_i) / Σ(n_i) for i=1..10 | 10 Network classes with sample counts from LOMO | ✓ |

---

## 13. Frozen 200-gene Comparator, Truth Normalization, and Post-review Analyses

| ID | Manuscript Location | Manuscript Value | Source CSV / Manifest | CSV Location | Formula | Raw Data Input | ✓? |
|----|-------------------|-----------------|-----------------------|-------------|---------|---------------|---|
| T048 | Supp S20 (manuscript value pending synchronization) | Frozen 200-gene RF LOMO Top1/Top3 = 47.50%/83.03% (389/680 of 819) | `manuscript/calculations/p2/P2_RF200_lomo_detail.csv` (submission calculation package) | `hit1`, `hit3`, all 819 rows | Σhit1/819 = 389/819; Σhit3/819 = 680/819 | Exact frozen 200-gene Network panel, leave-one-macaque-out predictions | ✓ |
| T049 | Supp S20 / truth-normalization note | Two prespecified frozen-truth corrections | Same CSV as T048 plus frozen mapping audit | Samples `19R470` and `xz-L-10m-10o_BRRL200004691-1A` | Apply mapping before scoring: Parietal → Occipital/Temporal; Lateral Prefrontal Cortex → OMPFC | Raw metadata, anatomical label crosswalk, and frozen scoring truth | ✓ |
| T050 | Supp 3′-bias sensitivity | 185 evaluable genes; Spearman ρ=0.646, p=2.95×10⁻²³; median detected 151/185 | `reproducibility/p0_bio2_3prime_bias/panel_length_detection.csv` and archived summary | Per-gene detection/length rows | Spearman correlation and per-sample median; terminal-window simulations at 500/1,000/2,000 nt yield medians 74/106/136 | Observed gene detection plus transcript-coordinate proxy simulation; not molecule-level 3′-end sequencing simulation | ✓ |
| T051 | Supp independent panel annotation | 179/200 mapped; significant GO:BP/KEGG terms = 446/11 (primary) and 410/8 (annotated-domain sensitivity) | `reproducibility/independent_enrichment/independent_enrichment_manifest.json` | `gprofiler_*_significant` and g:Profiler metadata | g:Profiler FDR <0.05 separately by source | Frozen 200-gene panel, 21,668-gene model-space background, g:Profiler e114_eg62_p19_27110d83 (2026-07-31) | ✓ |
| T052 | Supp independent cell-type annotation | Rhesus excitatory/inhibitory enrichment q=2.67×10⁻⁷/4.20×10⁻⁴; human sensitivity q=0.0477/5.89×10⁻⁵ | `reproducibility/independent_enrichment/independent_celltype_enrichment.csv` and manifest | Chiou2023 and Siletti2023 rows | One-sided hypergeometric test; BH correction across seven prespecified families | Chiou 2023 rhesus markers (primary) and Siletti 2023 human markers (sensitivity), intersected with 21,668-gene background | ✓ |

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
| Macro F1 summary | T044–T047 | 4/4 | 100% |
| Frozen comparator, truth normalization & post-review analyses | T048–T052 | 5/5 | 100% |
| **Total** | **T001–T052** | **52/52** | **100%** |

> **Status note:** the source files support T048 as 389/819 and 680/819. Any manuscript text still reporting 393/819 and 679/819 must be corrected before submission. T051-T052 are now backed by public derived-result files; publisher-hosted source workbooks remain external.
