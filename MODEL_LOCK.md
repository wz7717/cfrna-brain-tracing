# Canonical 110 production model lock

The production model is frozen under lock ID `canonical110-v0.1.9-20260727`.
No model artifact, learned value, gene panel, projector coefficient, reference
matrix, resolution-group definition, scoring weight, threshold, candidate-beam
contract, or production-route switch may change in place.

The authoritative machine-readable lock is
`data/models/canonical110_model_lock.json`. Production inference verifies the
declared file sizes and SHA-256 hashes before loading the model and fails closed
on any mismatch. The production parameter dictionary is independently fixed in
`core/production_route.py` and must exactly match the manifest before inference.

`resolution_group_mean_weight=0.10` is retained as a legacy/inactive manifest
field so the frozen lock remains byte-for-byte semantically stable. It is not
used by the production route. Resolution groups are ranked independently using
the locked local Top200 panel; the former best-score-plus-mean aggregation must
not be reintroduced.

The executable region route binds the locked Network Top3 size, local Top200
panel, Top50 and Top100 exact-region panel sizes, 0.25 fusion weight, and the
existing minimum 20-gene Region overlap threshold through implementation
constants. The 20-gene threshold is lock coverage added for an unchanged
production guard; it is not a new threshold or model change.

Any future scientific model change requires all of the following rather than an
edit to the existing lock:

1. a new model lock ID and software version;
2. a complete canonical-region, Network, beam, and artifact-integrity audit;
3. rerunning the locked internal and external validation suite;
4. updating manuscript-facing metrics if any result changes;
5. a new reviewed GitHub release and immutable Zenodo archive.

Pure presentation, documentation, accessibility, or deployment-maintenance
changes may proceed without changing this model lock only when they do not alter
the locked files or production parameters.
