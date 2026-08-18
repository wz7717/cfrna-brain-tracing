# Scientific remediation report

All corrected secondary statistics are regenerated from the current formal prediction/detail sources. The frozen model and primary endpoints are unchanged.

## Corrected source chains

- TCGA/BraTS broad strict Top3 range: max-min = `33.3089133089133` pp (`33.31` pp displayed) from `reproducibility/tcga_brats_truth_basis_top3_summary.json`.
- LOMO Exact: `177/812` Top1 = micro-F1 `0.217980295566502`; macro-F1 `0.203415333979799` across `104` truth-label classes.
- LOMO Exact and LOMO Group origin/staged/generator-input path+SHA pairs: `reproducibility/LOMO_INPUT_CHAIN_PROVENANCE.md` (2 endpoint chains; staged and generator-input pairs are identical by assertion).
- Benchmark: `51` profiles × `28415` genes; `153` warm inference events; cold peak `216.7383` MiB and warm maximum `222.0039` MiB.
- Resolution-group Top3 random baselines: LOSO uniform/weighted `0.225520062166910` / `0.065434152334152`; LOMO `0.213447535105209` / `0.040752709359606`.
- Friedman: χ²=`0.5385`, df=`2`, P=`0.764`; exact-enumeration status: `REMOVED`.

## Provenance boundary

The LOMO Exact macro denominator is the truth-label universe. Top1 predictions outside that universe remain false positives and are disclosed in the formal F1 provenance artifact; they do not alter the frozen model or route.

## QA status

- Current generated-evidence stale scan: `PASS`
- DOCX remediation: `PASS: final DOCX rendered and visually inspected (12 main + 82 supplementary pages)`
- Regression tests: `PASS: python -m pytest -q (143 passed, 2 skipped)`
