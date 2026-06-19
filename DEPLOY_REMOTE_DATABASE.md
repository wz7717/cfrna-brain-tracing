# Remote reference database deployment

The Streamlit Cloud app does not store `cfrna_source_tracing.db` in GitHub.
The database is large and remains ignored by Git. To enable full atlas browsing
and database-backed tracing in the cloud, host the SQLite file externally and
configure the app to download it on first startup.

## Upload target

Recommended public hosting options:

- Hugging Face Dataset
- Zenodo
- Figshare
- S3-compatible object storage

Upload the local file:

```text
cfrna_source_tracing.db
```

Use a direct download URL. For Hugging Face Datasets, the URL usually looks like:

```text
https://huggingface.co/datasets/<user>/<dataset>/resolve/main/cfrna_source_tracing.db
```

## Streamlit Cloud secret

In Streamlit Cloud, open:

```text
Manage app -> Settings -> Secrets
```

Add:

```toml
CFRNA_DB_URL = "https://your-direct-download-url/cfrna_source_tracing.db"
```

Then reboot the app. On first startup, the app downloads the database to:

```text
cfrna_source_tracing.db
```

If `CFRNA_DB_URL` is not configured, the app creates an empty SQLite database
with only the schema, so atlas browsing pages will not have reference atlas rows.

## Notes

- Keep `cfrna_source_tracing.db` out of GitHub.
- Keep `packages.txt` and `environment.yml` absent for Streamlit Cloud pip deployment.
- If the database is private, use a storage URL that Streamlit Cloud can access.
