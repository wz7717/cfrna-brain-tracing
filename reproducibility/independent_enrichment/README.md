# Independent frozen-panel enrichment package

This directory contains independent annotation analyses for the exact locked
200-gene Network panel. The files are derived results; the publisher-hosted
Chiou 2023 and Siletti 2023 source workbooks are not redistributed.

## Contents

- `analysis_protocol_20260731.json`: protocol frozen before computation.
- `independent_enrichment_manifest.json`: query settings, source versions,
  input SHA-256 values, and compact result summaries.
- `gprofiler_model_background.*`: primary GO:BP/KEGG results using the frozen
  21,668-gene model space as a custom background.
- `gprofiler_annotated_background.*`: annotated-domain sensitivity results.
- `gprofiler_GO_BP_representative_components.csv`: redundancy-reduced GO:BP
  representatives; full significant results remain in the two result exports.
- `independent_celltype_enrichment.csv`: all seven prespecified cell families
  for the Chiou rhesus primary analysis and Siletti human sensitivity analysis.
- `independent_primate_marker_sets.csv` and
  `independent_human_marker_sets.csv`: derived broad-family marker sets.

The public runner is `scripts/run_independent_panel_enrichment.py`. Re-querying
g:Profiler may change database-version-dependent results; the archived result
files and manifest are carried forward unchanged from the frozen v0.1.12
scientific model into the v0.1.14 release.

These analyses characterize annotation bias only. They do not establish cell
of origin, cell abundance, mechanism, causality, or predictive validity.
