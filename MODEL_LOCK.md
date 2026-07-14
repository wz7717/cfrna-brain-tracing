# Canonical 110 production model lock

The production model is frozen under lock ID `canonical110-v0.1.8-20260714`.
No model artifact, learned value, gene panel, projector coefficient, reference
matrix, resolution-group definition, scoring weight, threshold, candidate-beam
contract, or production-route switch may change in place.

The authoritative machine-readable lock is
`data/models/canonical110_model_lock.json`. Production inference verifies the
declared file sizes and SHA-256 hashes before loading the model and fails closed
on any mismatch. The production parameter dictionary is independently fixed in
`core/model_lock.py` and must exactly match the manifest.

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
