# S18-S19 portable reproduction

Run from the repository root:

```bash
python reproducibility/s18_s19/run_s18_s19_sensitivity.py
```

The runner reads the frozen model artifacts under `data/models` and the
corrected lambda/AHBA CSV inputs under `reproducibility`. It writes:

- `sensitivity_analysis_results.json`: the four distinct sensitivity blocks.
- `input_manifest.json`: portable input paths, byte sizes, and SHA-256 values.

The blocks must remain distinct:

- panel-size analysis is a 110-region centroid development comparison;
- grouping-threshold analysis compares precomputed structures;
- local-gene analysis is strict LOSO with a fixed Network beam;
- lambda analysis is formal strict LOSO with corrected donor-level inference.

They are not interchangeable retrainings and must not overwrite formal-route
results.
