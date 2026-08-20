# BrainTrace v0.1.18 pre-publication checklist

This immutable checklist records the state immediately before tag creation.
Post-tag GitHub and Zenodo publication will be verified separately by the
final release audit.

- [x] Version metadata synchronized from `release_manifest.json`.
- [x] Exact reviewed release tree committed.
- [x] Scientific-freeze verification passed with zero drift.
- [x] FULL reproduction reuse approved from an unchanged FULL-impact fingerprint.
- [x] Portable reproduction passed.
- [x] Application and reproducibility Docker build/smoke gates passed.
- [x] Two checksum-identical deterministic archive builds were produced.
- [x] GitHub Actions passed for the exact reviewed pre-publication commit.
- [x] Private-path and manuscript-production leakage audits passed.
- [x] Required Git LFS payloads were materialized.
- [ ] Create and push the annotated `v0.1.18` tag at the approved final main commit.
- [ ] Create the GitHub Release from the exact checked archives.
- [ ] Upload checksum-matched archives to the two existing-series Zenodo drafts.
- [ ] Verify remote metadata, files, checksums and concept-series membership.
- [ ] Publish the two existing Zenodo drafts without modifying v0.1.17.
- [ ] Verify DOI resolution, duplicate count and final post-publication consistency.
