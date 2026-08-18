# LOMO origin / staged / generator-input path-SHA pairing

For each endpoint, generator_input path/SHA-256 must exactly equal the repository-staged path/SHA-256; origin is retained as its distinct source pair.

| Endpoint | Role | Path | SHA-256 | Generator binding |
| --- | --- | --- | --- | --- |
| LOMO Exact | Origin | `D:\Download\文章改稿\code\reproduction_validation_workspace_20260802\reproduced_runs\release_gate_20260811_ai44563\lomo\formal_lomo_exact_region_detail.csv` | `EB5F10F01B122F68D09256EA6866DEAE2B439AABAD27E076181EC8760E7AAF36` | frozen external formal detail |
| LOMO Exact | Staged | `reproducibility/p2_publication_completeness/formal_lomo_exact_region_detail.csv` | `E90BB2ACF6C325530467A69252CBDF3B28C027EC41555FD03716527F5E3C101A` | repository-staged canonical table |
| LOMO Exact | Generator input | `reproducibility/p2_publication_completeness/formal_lomo_exact_region_detail.csv` | `E90BB2ACF6C325530467A69252CBDF3B28C027EC41555FD03716527F5E3C101A` | `scripts/generate_lomo_exact_f1_evidence.py` via `core.lomo_exact_f1.CANONICAL_FORMAL_PATH` |

| LOMO Group | Origin | `D:\Download\文章改稿\code\reproduction_validation_workspace_20260802\reproduced_runs\release_gate_20260811_ai44563\lomo\formal_lomo_resolution_group_detail.csv` | `B9A17D20BA434F52BD812FAE361E1A3F51C55B705AEFA745B8853544276390F1` | frozen external formal detail |
| LOMO Group | Staged | `reproducibility/p2_publication_completeness/formal_lomo_resolution_group_detail.csv` | `E3FB53B7B14F135B5B9781603515F8344CD9DE54CFA69CE2E06EEB8EAB21938A` | repository-staged canonical table |
| LOMO Group | Generator input | `reproducibility/p2_publication_completeness/formal_lomo_resolution_group_detail.csv` | `E3FB53B7B14F135B5B9781603515F8344CD9DE54CFA69CE2E06EEB8EAB21938A` | `scripts/generate_scientific_remediation_artifacts.py` via `core.resolution_group_baselines.CANONICAL_PATHS['LOMO']` |
