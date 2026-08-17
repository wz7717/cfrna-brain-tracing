# Zenodo release checklist — BrainTrace v0.1.16

Completed 2026-08-17. The immutable v0.1.16 software version DOI is
`10.5281/zenodo.21974954`; the full reproducibility archive DOI is
`10.5281/zenodo.21974991`.

After the GitHub v0.1.16 release is created:

- [x] Create a new Zenodo version from the v0.1.16 GitHub release.
- [x] Confirm title, creators, MIT license, repository URL, publication date
      and version `v0.1.16`.
- [x] Confirm the archived source contains `examples/`, Help & Tutorial,
      `VERSION_STRING_AUDIT_v0.1.16.txt` and
      `SUBMISSION_ONLINE_RESOURCE_QA_v0.1.16.md`.
- [x] Upload/materialize the full reproducibility archive, including the
      approximately 164 MB Git LFS payload, if the GitHub-generated Zenodo
      source archive does not materialize LFS objects.
- [x] Verify the full archive against its release checksum manifest.
- [x] Record the real Zenodo v0.1.16 version DOI; never reuse the v0.1.15 DOI
      as though it identified v0.1.16.
- [x] Add `RELEASE_INTEGRITY_v0.1.16.json` with the real Git tag, GitHub
      release URL, version DOI, full-archive DOI/asset and SHA-256 values.
- [x] Synchronize the real version DOI and archive wording in `README.md` and
      `DATA_PROVENANCE.md`.
- [x] Synchronize the DOI in the manuscript, Supplementary Material and Cover
      Letter where the current software release is cited.
- [x] Confirm the Zenodo concept DOI still resolves to the latest version.
- [x] Re-run `braintrace validate`, `python examples/verify_examples.py`, the
      full test suite and package checksum verification after metadata edits.

The immutable v0.1.15 records remain historical and must not be edited.
