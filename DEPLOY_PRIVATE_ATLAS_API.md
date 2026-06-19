# Private atlas API deployment

Do not upload the full `cfrna_source_tracing.db` to public object storage if it
contains Bo2023-derived matrices or other restricted reference data.

The safer architecture is:

```text
Streamlit Cloud app -> HTTPS API -> private server with cfrna_source_tracing.db
```

The Streamlit app receives only the rows returned by approved read-only API
endpoints. The SQLite database stays on your own server.

## Server setup

On the private server:

```bash
python -m venv .venv
. .venv/bin/activate
pip install -r requirements_api.txt
```

Place the full database on the server, for example:

```text
/srv/cfrna-braintrace/cfrna_source_tracing.db
```

Start the API:

```bash
export CFRNA_DB_PATH=/srv/cfrna-braintrace/cfrna_source_tracing.db
export CFRNA_API_KEY='replace-with-a-long-random-token'
uvicorn api.atlas_api:app --host 0.0.0.0 --port 8000
```

Put this behind HTTPS with Nginx, Caddy, a cloud load balancer, or a tunnel.

## Streamlit Cloud secrets

In Streamlit Cloud:

```text
Manage app -> Settings -> Secrets
```

Add:

```toml
CFRNA_API_URL = "https://your-api-domain.example.com"
CFRNA_API_KEY = "replace-with-the-same-long-random-token"
```

Then reboot the app.

## Available read-only endpoints

- `GET /health`
- `GET /atlas_versions`
- `GET /atlases/{atlas_id}/regions`
- `GET /atlases/{atlas_id}/celltypes`
- `GET /atlases/{atlas_id}/region-ranking`
- `GET /atlases/{atlas_id}/gene-candidates`
- `GET /atlases/{atlas_id}/expression`
- `GET /atlases/{atlas_id}/marker-evidence`

The API only exposes fixed read-only queries. It does not expose arbitrary SQL
or direct database file downloads.

## Security notes

- Keep `cfrna_source_tracing.db` off GitHub and public file hosts.
- Use HTTPS.
- Set `CFRNA_API_KEY` on the API server and Streamlit Cloud.
- Restrict firewall access if possible.
- Return only aggregate/reference rows that are allowed for your data license.
