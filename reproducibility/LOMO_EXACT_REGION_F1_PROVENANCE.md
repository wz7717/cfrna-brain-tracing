# Formal LOMO Exact F1 provenance

The current LOMO Exact endpoint is regenerated from the frozen prediction-level formal route.
Its macro denominator is the 104-label truth universe; Top1 labels outside that universe
remain false positives and do not create extra macro classes.

## Origin / staged / generator-input pairing

| Role | Path | SHA-256 | Binding |
| --- | --- | --- | --- |
| Origin | `external_source::formal_lomo_exact_origin/formal_lomo_exact_region_detail.csv` | `EB5F10F01B122F68D09256EA6866DEAE2B439AABAD27E076181EC8760E7AAF36` | frozen external formal prediction detail |
| Staged | `reproducibility/p2_publication_completeness/formal_lomo_exact_region_detail.csv` | `401441CFD7FF9B66408377CD854CF8A4C31B869F16C6D843E35DFCA63BE401C1` | repository-staged canonical table |
| Generator input | `reproducibility/p2_publication_completeness/formal_lomo_exact_region_detail.csv` | `401441CFD7FF9B66408377CD854CF8A4C31B869F16C6D843E35DFCA63BE401C1` | `scripts/generate_lomo_exact_f1_evidence.py` via `core.lomo_exact_f1.CANONICAL_FORMAL_PATH` |

The staged and generator-input path/SHA pairs must be identical; only the origin may differ because staging normalizes the frozen route rows.

- Route: `top3_beam_local_top50_top100_zfusion_w0p25`
- Route family: `hybrid_projected_network_logcpm_exact`
- Origin SHA-256: `EB5F10F01B122F68D09256EA6866DEAE2B439AABAD27E076181EC8760E7AAF36`
- Staged / generator-input SHA-256: `401441CFD7FF9B66408377CD854CF8A4C31B869F16C6D843E35DFCA63BE401C1`
- Top1 / Top3: `177/812`; `346/812`
- Prediction-only Top1 labels: `36c, Cla, MT`

## Derived summary

- macro_f1: `0.2034153339797993`
- sd_class_f1: `0.21593290846513288`
- median_class_f1: `0.15384615384615385`
- iqr_class_f1: `0.2857142857142857`
- weighted_f1: `0.21352686247014366`
- micro_f1: `0.21798029556650247`
- n_zero_f1_classes: `30`
- conditional_macro_f1_nonzero: `0.2858810099175558`

## Integer accounting

`sum(TP) = 177`; `sum(support) = 812`.
