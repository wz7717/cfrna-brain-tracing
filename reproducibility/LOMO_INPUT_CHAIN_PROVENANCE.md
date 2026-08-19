# LOMO origin / staged / generator-input path-SHA pairing

For each endpoint, generator_input path/SHA-256 must exactly equal the repository-staged path/SHA-256; origin is retained as its distinct source pair.

| Endpoint | Role | Path | SHA-256 | Generator binding |
| --- | --- | --- | --- | --- |
| LOMO Exact | Origin | `external_source::formal_lomo_exact_origin/formal_lomo_exact_region_detail.csv` | `EB5F10F01B122F68D09256EA6866DEAE2B439AABAD27E076181EC8760E7AAF36` | frozen external formal detail |
| LOMO Exact | Staged | `reproducibility/p2_publication_completeness/formal_lomo_exact_region_detail.csv` | `401441CFD7FF9B66408377CD854CF8A4C31B869F16C6D843E35DFCA63BE401C1` | repository-staged canonical table |
| LOMO Exact | Generator input | `reproducibility/p2_publication_completeness/formal_lomo_exact_region_detail.csv` | `401441CFD7FF9B66408377CD854CF8A4C31B869F16C6D843E35DFCA63BE401C1` | `scripts/generate_lomo_exact_f1_evidence.py` via `core.lomo_exact_f1.CANONICAL_FORMAL_PATH` |

| LOMO Group | Origin | `external_source::historical_formal_lomo_resolution_group_detail/formal_lomo_resolution_group_detail.csv` | `B9A17D20BA434F52BD812FAE361E1A3F51C55B705AEFA745B8853544276390F1` | frozen external formal detail |
| LOMO Group | Staged | `reproducibility/p2_publication_completeness/formal_lomo_resolution_group_detail.csv` | `685DA8F954490C70AAAEDA477EFBC86C9C4C622A8916D9BDEBC484747E1D736F` | repository-staged canonical table |
| LOMO Group | Generator input | `reproducibility/p2_publication_completeness/formal_lomo_resolution_group_detail.csv` | `685DA8F954490C70AAAEDA477EFBC86C9C4C622A8916D9BDEBC484747E1D736F` | `scripts/generate_scientific_remediation_artifacts.py` via `core.resolution_group_baselines.CANONICAL_PATHS['LOMO']` |
