# Round-5 projection and ranking audit

This directory contains three prespecified diagnostic analyses. They do not
retrain or replace the locked BrainTrace route.

## 1. Network-layer direct projection ablation

The same 819 Bo2023 samples, normalized 10-class Network truth, locked 200-gene
panel, fold-local leave-one-sample-out centroids and Pearson scoring were used
for all three rows. Only the Network-layer representation changed:

1. raw-count-derived logCPM;
2. released native VSD;
3. reference-fitted projected VSD.

No resolution-group/exact-region reranking or pairwise rescue was used. Thus
the analysis isolates the direct Network candidate-filtering stage and is not
the complete production-route validation. Donor-paired differences were
audited with an exact sign-flip test over nine donor-level estimates; those
P values are exploratory and unadjusted.

## 2. Projector OLS fit-quality audit

Gene-wise slopes and intercepts were fitted by ordinary least squares to paired
Bo2023 logCPM and released VSD values. The audit reports the distribution of
per-gene R-squared, Spearman correlation, residual standard deviation and
slope for all 21,668 projector genes and for the locked 200-gene Network panel.
These are in-sample paired-reference engineering-fit diagnostics. They are not
an equivalence test for DESeq2 VST, an external calibration result, or evidence
that query-cohort distributions are interchangeable with Bo2023.

## 3. MRR and NDCG@3

MRR and NDCG@3 were calculated from the frozen LOSO and no-pairwise LOMO rank
outputs at Network, resolution-group and exact-region levels. For the two fine
levels, a truth label absent from the retained candidate set receives reciprocal
rank 0 and NDCG@3 0. NDCG@3 uses one binary relevant truth label, so ideal DCG
is 1. These metrics supplement rather than replace Top1/Top3 accuracy.

## 4. P1 E4 per-Network mapping impact

`run_p1_e4_ahba_mapping_impact.py` joins the prespecified row-occurrence
humanization audit to the frozen AHBA sample-detail output. For each of the ten
Networks, the denominator is the number of AHBA samples whose allowed
Network-label set contains that Network. The primary numerator is the
set-valued any-allowed Network Top1 hit; a strict same-Network Top1 numerator is
reported as a sensitivity. Multi-label truth means that denominators overlap.
The percentage-point effect is each primary per-Network hit rate minus the
pooled 165/223 mapped-label Top1 rate. With two AHBA donors, these are
descriptive transfer metrics and not Network-level inferential tests.

## Files

- `network_projection_ablation_summary.csv`: sample-weighted direct-ablation
  results.
- `network_projection_ablation_by_donor.csv`: donor-level estimates.
- `network_projection_ablation_signflip.csv`: exploratory paired sign-flip
  comparisons.
- `projector_ols_quality_summary.csv`: OLS distribution summaries.
- `projector_ols_quality_manifest.json`: audited input scope and fallback count.
- `ranking_metrics_mrr_ndcg.csv`: LOSO/LOMO ranking metrics.
- `../p1_e4_ahba_network_mapping_impact.csv`: per-Network mapping fraction,
  AHBA denominator, primary/strict hits and percentage-point effects.
- `../p1_e4_ahba_network_mapping_impact_manifest.json`: input hashes and
  metric definitions for the P1 E4 output.
- `run_bo2023_projected_vsd_loso.py`: direct-ablation runner.
- `analyze_round5_projection_and_ranking.py`: summary and ranking-metric runner.
- `run_p1_e4_ahba_mapping_impact.py`: P1 E4 mapping-impact runner.

The scripts expect the external Bo2023 processed matrices and frozen validation
detail files from the complete reproduction package. Their placement and
SHA-256 values are documented in the repository `README.md` and
`DATA_PROVENANCE.md`.
