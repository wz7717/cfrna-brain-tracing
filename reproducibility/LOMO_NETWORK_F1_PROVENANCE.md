# Formal LOMO Network F1 provenance

This report is generated from the prediction-level formal LOMO Network source.
It is a reporting/provenance correction; the frozen model, ontology, route and
prediction set are unchanged.

- Route: `network_discriminative_correlation_top200`
- Route family: `hybrid_projected_network_logcpm_exact`
- Canonical source: `reproducibility/p2_publication_completeness/formal_lomo_network_detail.csv`
- Canonical source SHA-256: `37751431501D8A334AE3A73609A60EC284531A285C576225FA3CE9BE6651F0DA`
- Prediction rows: `819`
- Top1: `455/819 = 0.5555555556`
- Top3: `750/819 = 0.9157509158`

## Derived metrics

| class | support | TP | FP | FN | precision | recall | F1 |
|---|---:|---:|---:|---:|---:|---:|---:|
| Temporal | 193 | 89 | 45 | 104 | 0.6641791045 | 0.4611398964 | 0.5443425076 |
| Lateral Prefrontal Cortex | 107 | 42 | 24 | 65 | 0.6363636364 | 0.3925233645 | 0.4855491329 |
| Orbitomedial Prefrontal Cortex (OMPFC) | 104 | 47 | 32 | 57 | 0.5949367089 | 0.4519230769 | 0.5136612022 |
| Parietal, and Parieto-occipital region | 98 | 71 | 97 | 27 | 0.4226190476 | 0.7244897959 | 0.5338345865 |
| Operculum/Insula | 81 | 49 | 90 | 32 | 0.3525179856 | 0.6049382716 | 0.4454545455 |
| Occipital/Temporal | 70 | 48 | 16 | 22 | 0.7500000000 | 0.6857142857 | 0.7164179104 |
| Frontal (agranular frontal motor areas) | 65 | 36 | 23 | 29 | 0.6101694915 | 0.5538461538 | 0.5806451613 |
| Subcortical | 54 | 41 | 0 | 13 | 1.0000000000 | 0.7592592593 | 0.8631578947 |
| Cingulate gyrus | 39 | 24 | 30 | 15 | 0.4444444444 | 0.6153846154 | 0.5161290323 |
| Hippocampal formation | 8 | 8 | 7 | 0 | 0.5333333333 | 1.0000000000 | 0.6956521739 |

- Macro-F1: `0.5894844147`
- Class-level sample SD: `0.1291730211`
- Median: `0.5390885471`
- IQR: `0.1526222611`
- Weighted-F1: `0.5604715495`
- Micro-F1: `0.5555555556` (= Top1 accuracy)
- Zero-F1 classes: `0` (`0.0000000000`)

## Root cause of the superseded values

The historical `0.61845` macro-F1 and `0.5812749695` weighted-F1 came from the stale rounded `macro_f1_class_data.json` LOMO rows. The historical `0.5802962149` micro-F1 was then estimated as `sum(recall*n)/819`; because the recalls were rounded, this was not an integer prediction-level TP count. The historical validation report also described a separate pairwise-rescue route (`network_pairwise_correlation_rescue_top3`), which is retained as historical evidence but is not the current formal endpoint.

Integer accounting check: `sum(TP) = 455` and `sum(support) = 819`.
