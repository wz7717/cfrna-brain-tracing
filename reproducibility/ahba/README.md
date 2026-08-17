# AHBA endpoint-evaluability ledger

`ahba_endpoint_evaluability_ledger.csv` is the canonical accounting artifact for
the current AHBA external-validation endpoints. It is derived from the
`hybrid_projected_network_logcpm_exact` rows in the formal sample-detail output
by `scripts/generate_nonhuang_scientific_provenance_artifacts.py`.

The unit is one replicate-collapsed AHBA tissue sample. The canonical counts are
231 replicate-collapsed tissues, 223 Network-evaluable samples, and 88
resolution-group/exact-region-evaluable samples. These are endpoint-specific
eligibility counts, not stages of one sequential attrition pipeline.

The ledger retains source-supported truth-label counts for Network and exact
region. A candidate-independent resolution-group truth-label count is not
available in the formal sample-detail output, so `group_truth_label_count` is
intentionally blank. The companion
`source_group_truth_label_count_in_candidate_beam` field is retained only as
candidate-beam-dependent provenance; it must not be used as an endpoint
denominator.

The single-label sensitivity subsets are 56 Network samples and 40
group/exact-region samples. Their hit counts and denominators are recomputed
from ledger rows in `ahba_endpoint_evaluability_summary.json`.

`../v4_p0_5_ahba_trace.csv` and
`../v4_p0_5_ahba_trace_manuscript_aligned.csv` are retained as **HISTORICAL
ENGINEERING TRACE — NOT THE CANONICAL ENDPOINT-EVALUABILITY LEDGER**. They may
describe intermediate engineering filters and cannot be substituted for this
ledger.
