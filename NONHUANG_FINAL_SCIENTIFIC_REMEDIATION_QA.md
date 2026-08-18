# NON-HUANG final scientific-provenance remediation QA

## 1. Base commit

- `79edb23738f7167e277acb1890fe38af06147530` (`79edb23 Merge pull request #28 from wz7717/agent/v0.1.16-release-state-finalize`)
- `origin/main` was fetched at the end of this audit and resolves to the same SHA.

## 2. Branch

- `fix/nonhuang-final-scientific-provenance`

## 3. Files changed

- Provenance and traceability: `README.md`, `DATA_PROVENANCE.md`, `reproducibility/DATA_PROVENANCE.md`, `reproducibility/TRACEABILITY_MATRIX.md`, `NONHUANG_SCIENTIFIC_CONFLICT_LEDGER.csv`.
- Canonical derived records: `reproducibility/ahba/`, `reproducibility/tcga_brats_truth_basis_top3_summary.csv`, `reproducibility/tcga_brats_truth_basis_top3_summary.json`, `reproducibility/orthology_humanization_summary.json`, `reproducibility/tier_cascade_loso_summary.json`, `reproducibility/sign_flip_current_family.csv`, `reproducibility/sign_flip_current_family.json`, `reproducibility/nonhuang_scientific_arithmetic_qa.json`, and corrected `reproducibility/p1_bio1_4/` outputs.
- Reproducibility code and tests: `reproducibility/generate_all_csvs.py`, `scripts/generate_nonhuang_scientific_provenance_artifacts.py`, `scripts/verify_nonhuang_scientific_provenance.py`, `scripts/run_p1_bio1_4_audit.py`, `scripts/run_p1_stat1_4_audit.py`, `scripts/run_p2_audit.py`, `scripts/build_p2_change_list.py`, `tests/test_nonhuang_scientific_provenance.py`.
- AHBA historical-trace classification: `reproducibility/v4_p0_5_ahba_trace.csv` and `reproducibility/v4_p0_5_ahba_trace_manuscript_aligned.csv`.
- The two revised submission DOCX files are intentionally untracked and are not part of the Git change set.

## 4. AHBA authoritative source

The Tier-1 source is `ahba_formal_three_tier_sample_detail.csv` (external source basename; SHA-256 `82eb3ead890058b18599f4bcf5a77f0c352626f8e3a3f18eecf7a14210654515`), route `hybrid_projected_network_logcpm_exact`.

## 5. AHBA canonical ledger

- `reproducibility/ahba/ahba_endpoint_evaluability_ledger.csv`
- `reproducibility/ahba/ahba_endpoint_evaluability_summary.json`
- `reproducibility/ahba/README.md`

The ledger reports only fields available in the sample-level source. In particular, the group truth-label count is intentionally blank because the source contains a candidate-beam count rather than a scientific group-truth cardinality.

## 6. AHBA counts and performance

- Replicate-collapsed tissues: 231.
- Endpoint-specific evaluability: Network 223; resolution-group 88; exact-region 88.
- Single-label sensitivity subsets: Network 56; resolution-group/exact-region 40.
- Network Top1/Top3: 165/223 = 73.99%; 211/223 = 94.62%.
- Group Top1/Top3: 37/88 = 42.05%; 60/88 = 68.18%.
- Exact Top1/Top3: 24/88 = 27.27%; 40/88 = 45.45%.
- Single-label Top1/Top3: Network 39/56 = 69.64% and 44/56 = 78.57%; group 19/40 = 47.50% and 29/40 = 72.50%; exact 11/40 = 27.50% and 20/40 = 50.00%.

## 7. Historical AHBA traces

`reproducibility/v4_p0_5_ahba_trace.csv` and `reproducibility/v4_p0_5_ahba_trace_manuscript_aligned.csv` are retained and explicitly labelled **HISTORICAL ENGINEERING TRACE — NOT THE CANONICAL ENDPOINT-EVALUABILITY LEDGER**. No current document treats them as a unique sequential sample-flow account.

## 8. TCGA/BraTS primary edema comparator

The 65-patient Tier-1 truth/prediction CSV (`brats_tcga_lgg_65_mri_truth_and_predictions.csv`) yields 63 eligible edema comparisons after excluding `TCGA-HT-7686` (no label-2 edema voxels) and `TCGA-HT-7680` (cerebellar/posterior-fossa outside the current label space). The corrected audit derives this from patient-level data; no underlying hits changed.

## 9. TCGA/BraTS Network strict Top3 range

Center 19/65 = 29.23%; core 12/65 = 18.46%; edema 15/63 = 23.81%; whole tumour 10/65 = 15.38%. The source-derived range is **13.85 percentage points**.

## 10. TCGA/BraTS broad strict Top3 range

Center 32/65 = 49.23%; core 45/65 = 69.23%; edema 52/63 = 82.54%; whole tumour 46/65 = 70.77%. The source-derived range is **33.31 percentage points**.

## 11. Orthology/humanization counts

The Tier-1 `signature_humanization_8800.csv` reports 5,324 humanized and 3,476 unmapped values, totaling 8,800. The frozen Top200 orthology universe is 188/200 humanizable; the separate g:Profiler mapping universe is 179/200.

## 12. Orthology denominator unit

All 8,800 values are **gene-by-region row occurrences**, not independent genes. Thus 5,324/8,800 = 60.50% and 3,476/8,800 = 39.50% are row-level rates under the frozen mapping rule. Current wording explicitly keeps this unit distinct from the two 200-gene mapping universes.

## 13. Tier-cascade exact-evaluable universe

From the Tier-1 strict-LOSO exact detail, 368 exact Top3 hits plus 446 misses equals 814 exact-evaluable samples. In that same sample universe, Network truth is retained for 750 samples and absent from the candidate set for 64.

## 14. Tier-cascade Network-miss denominator

The denominator-consistent result is **64/446 = 14.35%**. It replaces the incorrect full-Network numerator of 66. The submission text renders the same result as “64 of 446 exact Top3 misses (14.35%).”

## 15. Tier-cascade conditional exact result

Conditional exact-region Top3 is **368/750 = 49.07%**. Exact-region recovery after a Network candidate-set miss is 0/64 = 0%.

## 16. Tier-cascade group conditional status

The same sample-level recomputation supports group Top3 = 590/750 = 78.67% and 0/64 recovery. The historical 78.32% claim is not retained in current scientific prose or derived outputs.

## 17. Current sign-flip family

The four-test Tier-1 family is: Network Top1 raw/BH 0.031250/0.125000; Network Top3 0.375000/0.500000; resolution-group Top3 0.593750/0.593750; exact-region Top3 0.324219/0.500000. None is significant.

## 18. TRACEABILITY_MATRIX

`reproducibility/TRACEABILITY_MATRIX.md` now identifies this as an unreleased v0.1.17 scientific-provenance patch candidate, points current AHBA accounting to the new ledger, records TCGA edema n=63 and 33.31 pp, uses the denominator-consistent tier results, and gives the current four-test sign-flip family. No DOI was added or changed.

## 19. README and provenance synchronization

`README.md`, `DATA_PROVENANCE.md`, and `reproducibility/DATA_PROVENANCE.md` now use endpoint-specific AHBA wording, 63 edema cases, and row-occurrence orthology semantics. They do not claim a v0.1.17 release.

## 20. Full test suite

- `python -m py_compile` on every changed Python file: PASS.
- `python -m pytest -q`: **123 passed, 2 skipped, 0 failed**; pytest-reported runtime 6.99 s (wall-clock 8.19 s).
- `scripts/verify_nonhuang_scientific_provenance.py`: **34/34 PASS**; machine-readable result at `reproducibility/nonhuang_scientific_arithmetic_qa.json`.

## 21. Residual grep

The required terms were searched across text/source files, excluding `.git` and the untracked DOCX render directory. There are no current incorrect claims: no `"edema": 64`, `edema=64`, `231 -> 223 -> 200`, `66/446`, `14.8%`, `78.32%`, or legacy sign-flip value as a current validation claim.

Remaining literal matches were classified as: (a) the conflict ledger and generator’s historical-statement fields, which deliberately preserve the superseded value for auditability; (b) DOCX updater source markers used solely to find and replace historical source-document text; (c) unrelated decimal substrings in unchanged per-gene numerical CSVs; or (d) correct negations such as “not a unique sequential attrition pipeline.” The revised DOCX outputs contain none of the prohibited stale phrases.

## 22. Frozen-model safety

No `data/models/*.npz`, `data/models/canonical110_model_lock.json`, Network scoring code, ontology, or locked model output was modified. The locked endpoints remain unchanged: LOSO Network 483/819 and 753/819; LOMO Network 455/819 and 750/819; LOSO group 368/814 and 590/814; LOMO group 344/812 and 569/812; LOSO exact 182/814 and 368/814; LOMO exact 177/812 and 346/812.

## 23. Git-diff safety

`git diff --check` is clean. The final diff contains provenance, derived QA, documentation, tests, AHBA ledger, and corrected audit code only. It excludes model NPZ/model-lock files, raw datasets, databases, DOCX files, release-integrity JSON, and tags. No merge, tag, release, or change to `v0.1.16` was performed; the feature branch was subsequently pushed for independent final audit.

## 24. DOCX synchronization

Untracked revised copies were produced from the user-supplied current submission files:

- `manuscript_remediation/BrainTrace_Main_Manuscript_NonHuangRemediated.docx`
- `manuscript_remediation/BrainTrace_Supplementary_File_NonHuangRemediated.docx`

They update only non-Huang AHBA, TCGA/BraTS, orthology, and tier-cascade wording. No cover letter was supplied, so its status is not applicable.

## 25. DOCX render and structure QA

All rendered pages were reviewed: main manuscript 12 pages and supplement 81 pages. No clipping, overflow, broken table, spurious blank page, wrong orientation, orphan heading, caption problem, or S20 row split was found. Both outputs retain Times New Roman 12 pt Normal style, zero tracked-change elements, zero comments, and the source section orientations (main portrait; supplement portrait/landscape/portrait/landscape/portrait).

## 26. Remaining unresolved issues

None. `NONHUANG_SCIENTIFIC_CONFLICT_LEDGER.csv` has 5 `FIXED_STALE_VALUE` and 1 `DOCUMENTATION_SYNC` entries; it has no `UNRESOLVED` entries.

## Acceptance checklist

| Item | Status |
| --- | --- |
| AHBA_ENDPOINT_LEDGER | PASS |
| AHBA_NUMERATORS_DENOMINATORS | PASS |
| TCGA_EDema_N63 | PASS |
| TCGA_RANGE_33_31 | PASS |
| ORTHOLOGY_ROW_UNIT | PASS |
| TIER_CASCADE_64_446 | PASS |
| TIER_EXACT_49_07 | PASS |
| SIGN_FLIP_TRACEABILITY | PASS |
| TRACEABILITY_MATRIX | PASS |
| README_PROVENANCE | PASS |
| FULL_TESTS | PASS |
| MODEL_LOCK | PASS |
| GIT_SAFETY | PASS |
| DOCX_SYNC_IF_APPLICABLE | PASS |

**READY FOR INDEPENDENT FINAL AUDIT**
