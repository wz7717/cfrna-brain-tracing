# BrainTrace — Raw Data Provenance Manifest

> **Principle**: Every raw data file is traceable to a public repository with a documented accession number, download URL, and content fingerprint.

## 1. Primary Atlas Data

### 1.1 Bo2023 Macaque Transcriptomic Atlas

| Field | Value |
|-------|-------|
| **Label** | `bo2023_vsd_matrix` |
| **File** | `bo2023 data/mfas5_819samples_23605genes_vsd4_rmbatch.xls` |
| **Source** | NCBI SRA / Bo et al. (2023) supplementary data |
| **Accession** | PRJNA905082 |
| **Citation** | Bo, T. et al. (2023) Brain-wide and cell-specific transcriptomic insights into MRI-derived cortical morphology in macaque monkeys. *Nat. Commun.*, 14, 1499. |
| **Content** | 819 samples × 23,605 genes, variance-stabilized (VSD4), batch-corrected |
| **Donors** | 9 adult *Macaca fascicularis* (cynomolgus macaques) |
| **Regions** | 110 brain regions (cerebral cortex + subcortical; cerebellum excluded) |
| **Role** | Primary reference atlas; source of all internal validation (LOSO/LOMO) |

| Field | Value |
|-------|-------|
| **Label** | `bo2023_counts_matrix` |
| **File** | `bo2023 data/mfas5_819samples_28415genes_featurecounts_counts.txt` |
| **Source** | PRJNA905082 (raw counts before VSD transformation) |
| **Content** | 819 samples × 28,415 genes, raw featureCounts |
| **Role** | Used for logCPM normalization and marker selection |

| Field | Value |
|-------|-------|
| **Label** | `bo2023_sample_info` |
| **File** | `bo2023 data/Information of sequenced samples_update_full878_filter819.xlsx` |
| **Source** | PRJNA905082 (sample metadata) |
| **Content** | Sample IDs, MonkeyID (donor mapping), brain region annotations |
| **Role** | Donor mapping for LOSO/LOMO splits; region label assignment |

| Field | Value |
|-------|-------|
| **Label** | `bo2023_gene_map` |
| **File** | `bo2023_bulk_atlas_buildkit/04_expressed_genes_neocortex_plus_subcortical.csv` |
| **Source** | Bo2023 build kit (curated gene filter) |
| **Content** | Expressed genes in neocortex + subcortical structures |
| **Role** | Gene filtering for atlas database construction |

---

## 2. External Validation Data

### 2.1 Allen Human Brain Atlas (AHBA)

| Field | Value |
|-------|-------|
| **Label** | `ahba_tpm_matrix` |
| **File** | `data/ahba_human_rnaseq/ahba_human_rnaseq_tpm_gene_symbol_matrix.tsv` |
| **Source** | Allen Brain Map portal (https://portal.brain-map.org) |
| **Citation** | Hawrylycz, M.J. et al. (2012) An anatomically comprehensive atlas of the adult human brain transcriptome. *Nature*, 489, 391–399. |
| **Content** | Human brain RNA-seq TPM matrix, gene-symbol indexed |
| **Donors** | 6 total; 2 with whole-brain coverage (4 excluded for incomplete sampling) |
| **Samples** | 231 independent tissue samples (post technical-replicate collapse) → 223 Network-qualified → 88 exact-region evaluable |
| **Role** | Mapped-label transfer validation (not anatomical truth) |

| Field | Value |
|-------|-------|
| **Label** | `ahba_metadata` |
| **File** | `data/ahba_human_rnaseq/ahba_human_rnaseq_sample_metadata_242.csv` |
| **Source** | Allen Brain Map portal |
| **Content** | 242 sample metadata entries (donor ID, structure name, hemisphere) |
| **Role** | Sample-to-region mapping for AHBA label harmonization |

### 2.2 TCGA/BraTS Glioma Dataset

| Field | Value |
|-------|-------|
| **Label** | `tcga_expression_file_level` |
| **File** | `data/tcga_brain_tumor_expression/tcga_gbm_lgg_primary_tumor_tpm_unstranded_file_level.tsv` |
| **Source** | NCI Genomic Data Commons (https://portal.gdc.cancer.gov) |
| **Citation** | TCGA Research Network (https://www.cancer.gov/tcga) |
| **Content** | TCGA GBM/LGG primary tumor TPM expression (file-level) |
| **Patients** | 65 glioma patients; 63 primary edema-comparator cases after excluding TCGA-HT-7686 (no edema voxels) and TCGA-HT-7680 (cerebellar/out-of-scope edema) |
| **Role** | Tumor-tissue expression for brain-origin tracing evaluation |

| Field | Value |
|-------|-------|
| **Label** | `tcga_expression_sample_mean` |
| **File** | `data/tcga_brain_tumor_expression/tcga_gbm_lgg_primary_tumor_tpm_unstranded_sample_mean.tsv` |
| **Source** | NCI GDC (sample-mean aggregated) |
| **Role** | Per-patient mean expression for Network/lobe/broad evaluation |

| Field | Value |
|-------|-------|
| **Label** | `brats_training_dir` |
| **File** | `data/brats_tcga_lgg_training_65/` |
| **Source** | The Cancer Imaging Archive (TCIA) |
| **Citation** | Bakas, S. et al. (2017) Advancing The Cancer Genome Atlas glioma MRI collections with expert segmentation labels and radiomic features. *Sci. Data*, 4, 170117. |
| **Content** | 65 patient directories with pre-operative TCGA-LGG NIfTI images and segmentations |
| **Nested archive** | `PKG - BraTS-TCGA-LGG/BraTS-TCGA-LGG/Pre-operative_TCGA_LGG_NIfTI_and_Segmentations.zip` (562 MB, 454 files) |
| **Extraction** | Auto-extracted by `reproduce_all.py` PRE-BUILD phase |
| **Role** | MRI-derived tumor location truth (center, core, edema, whole_tumor) for anatomical consistency evaluation |

### 2.3 GSE189919 Engineering Benchmark

| Field | Value |
|-------|-------|
| **Label** | N/A (accessed via GEO) |
| **File** | Accessed from NCBI GEO (not stored as raw file) |
| **Source** | NCBI Gene Expression Omnibus (https://www.ncbi.nlm.nih.gov/geo/) |
| **Accession** | GSE189919 |
| **Citation** | Lee, B. et al. (2022) Medulloblastoma cerebrospinal fluid reveals metabolites and lipids indicative of hypoxia and cancer-specific RNAs. *Acta Neuropathol. Commun.*, 10, 25. |
| **Content** | 51 CSF cfRNA samples (no anatomical truth) |
| **Role** | Engineering performance benchmark (cold/warm inference timing) + gene-overlap concept feasibility |

### 2.4 Huang2025 cfRNA Domain-Shift Audit

| Field | Value |
|-------|-------|
| **Label** | `huang2025_cfRNA` |
| **File** | `external_inputs/huang2025_pmc12041490/41698_2025_909_MOESM2_ESM.csv` |
| **Source** | Published article Supplementary Data 1 |
| **Citation** | Huang, J. et al. (2025) Diagnostic and prognostic potential of cell-free RNAs in cerebrospinal fluid and plasma for brain tumors. *npj Precis. Oncol.*, 9, 123. [doi:10.1038/s41698-025-00909-6](https://doi.org/10.1038/s41698-025-00909-6) |
| **Content** | 159 CSF/plasma cfRNA expression profiles |
| **Role** | Domain-transfer audit (empirical domain-shift stress test) |

---

## 3. Build Artifacts (Derived, Not Raw)

These files are **generated by the pipeline**, not raw data. They are listed here to clarify the data lineage.

| Artifact | Generated By | From Raw Data | Size |
|----------|-------------|---------------|------|
| `braintrace_source_tracing.db` | `build_bo2023_atlas_from_wang_matrix.py` (Step 1) | Bo2023 VSD + counts + metadata + gene map | ~730 MB |
| `bo2023_saleem_network_top200_model.npz` | `build_bo2023_network_model.py` (Step 2) | Bo2023 VSD + metadata | ~2 MB |
| `bo2023_reference_projector_linear_full.npz` | `build_bo2023_reference_projector.py` (Step 3) | VSD + counts + metadata + network model | ~5 MB |
| `bo2023_formal_region_logcpm_reference_matrix.npz` | `build_bo2023_reference_projector.py` (Step 3) | VSD + counts + metadata + network model | ~8 MB |
| `braintrace_reference_v1.sqlite` | `build_sqlite_reference_db_v1.py` (Step 4) | Reference atlas manifest | ~50 MB |

> **Key architectural note**: `braintrace_source_tracing.db` is NOT a pre-existing dependency. It is regenerated from raw data by Step 1 of the pipeline. The pre-built DB is provided as a convenience to skip the expensive (~30 min) build step via `--skip-build`.

---

## 4. Manuscript CSV Files (Dynamically Generated)

All 10 supplementary CSV files accompanying the manuscript are **dynamically generated** from raw counts by `reproducibility/generate_all_csvs.py`. They are NOT manually created.

| CSV File | Raw Count Source | Formulas Applied |
|----------|-----------------|------------------|
| `v4_p0_9_triple_ci.csv` | LOSO/LOMO validation counts | Wilson CI, Clopper-Pearson CI, Agresti-Coull CI |
| `v4_p0_4_subcortical_subsampling.csv` | Subcortical 42/54 recall, 42/42 PPV | Bootstrap subsampling (2000 reps, 9 fractions) |
| `v4_p0_5_ahba_trace.csv` | AHBA attrition steps (documented) | N/A (trace table) |
| `v4_p0_5_ahba_trace_manuscript_aligned.csv` | AHBA 5-step manuscript version | N/A (trace table) |
| `v4_p0_10_lambda_friedman.csv` | Per-donor lambda hit rates | Friedman chi-squared test |
| `v4_p0_10_lambda_sensitivity.csv` | Lambda × endpoint accuracy | Direct from validation output |
| `v4_p0_11_ml_baselines.csv` | ML baseline LOMO counts | Wilson CI |
| `v4_p0_11_rf_comparator.csv` | Random Forest LOMO counts | Wilson CI, Clopper-Pearson CI |
| `v4_p0_12_tcga_brats_ci_summary.csv` | TCGA/BraTS per-patient evaluation | Wilson CI, CP CI, AC CI, Bootstrap CI, Binomial test |
| `v4_p0_13_macro_f1.csv` | LOSO class data plus prediction-level formal LOMO Network source | precision=TP/(TP+FP), recall=TP/(TP+FN), F1=2PR/(P+R), macro/weighted summary |

The formal LOMO Network F1 chain is generated by
`scripts/generate_lomo_f1_evidence.py` from
`p2_publication_completeness/formal_lomo_network_detail.csv`, which is the
819-row prediction-level source for route
`network_discriminative_correlation_top200` and route family
`hybrid_projected_network_logcpm_exact`. The accompanying
`formal_lomo_network_f1.csv`, `lomo_network_f1_provenance.json` and
`LOMO_NETWORK_F1_PROVENANCE.md` retain integer TP/FP/FN accounting and the
root-cause explanation for the superseded rounded summary.

**Regeneration command**: `python reproducibility/generate_all_csvs.py`
For the formal LOMO evidence chain, run:
`python scripts/generate_lomo_f1_evidence.py`.

---

## 5. Data Integrity Verification

### SHA256 Checksums
All raw data files are SHA256-verified against `SHA256SUMS.txt` before use.
The verification is performed by `reproduce_all.py` Step 0 (`verify_raw_data()`).

### Seed Registry
Every stochastic operation uses a documented seed:

| Seed Name | Value | Used For |
|-----------|-------|----------|
| `donor_cluster_inference` | 20260716 | 50,000 donor-cluster bootstrap |
| `p0_hard_evidence` | 20260629 | P0 hard-evidence generation |
| `benchmark_stratified_kfold` | 42 | Stratified k-fold benchmark |
| `bootstrap_default` | 13 | Core/params.py default bootstrap |
| `sparse_simulation_base` | 20260711 | 30-repeat sparse simulation (99,099 archival seeds) |
| `rf_fair_comparator` | 20260717 | Random Forest comparator |
| `per_region_bootstrap` | 20260717 | 418 per-region bootstrap records |
| `sign_flip_test` | exact_2^9 | Exhaustive 512 sign patterns (deterministic) |

### Version Pinning
The complete Python environment is pinned in `requirements_reproducible.txt`:
- Python 3.11.9
- numpy 1.26.4, scipy 1.13.1, pandas 2.2.3
- scikit-learn 1.5.2, streamlit 1.55.0
- Full list: 24 packages with exact versions

---

## 6. Public Access Summary

| Dataset | Repository | Accession/DOI | Access |
|---------|-----------|---------------|--------|
| Bo2023 macaque atlas | NCBI SRA | PRJNA905082 | Public |
| Allen Human Brain Atlas | Allen Brain Map | https://portal.brain-map.org | Public |
| TCGA glioma expression | NCI GDC | https://portal.gdc.cancer.gov | Public |
| BraTS-TCGA-LGG MRI | TCIA | Bakas et al. (2017) | Public |
| GSE189919 CSF cfRNA | NCBI GEO | GSE189919 | Public |
| Huang2025 cfRNA | npj Precis. Oncol. | Supplementary Data 1 | Public |
| BrainTrace code | GitHub | https://github.com/wz7717/cfrna-brain-tracing | Public (MIT) |
| BrainTrace archive | Zenodo | v0.1.14 historical software/full-archive records remain at 10.5281/zenodo.21920261 and 10.5281/zenodo.21920697; v0.1.12 scientific release 10.5281/zenodo.21911532; software concept 10.5281/zenodo.20773674 | Public |

> **FAIR statement**: Public source datasets are identified by stable accessions or article records and use standard formats. The exact Bo2023 processed author-package matrices are an explicitly documented exception: no stable public file URL was identified, and exact reuse requires lawful acquisition from the original authors plus SHA-256 verification.
