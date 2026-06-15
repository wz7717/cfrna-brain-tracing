#!/usr/bin/env python
from __future__ import annotations

import argparse
import json
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path


BASE = "https://faspex.cancerimagingarchive.net"
CLIENT_ID = "ff9aa63a-72e1-436f-82ef-5677eb1f7aee"
ROOT = Path(__file__).resolve().parents[1]


class NoRedirect(urllib.request.HTTPRedirectHandler):
    def redirect_request(self, req, fp, code, msg, headers, newurl):  # type: ignore[no-untyped-def]
        return None


def read_response(req: urllib.request.Request, opener: urllib.request.OpenerDirector | None = None) -> tuple[int, dict[str, str], str]:
    opener = opener or urllib.request.build_opener()
    try:
        with opener.open(req, timeout=120) as resp:
            return int(resp.status), dict(resp.headers), resp.read().decode("utf-8", errors="ignore")
    except urllib.error.HTTPError as exc:
        return int(exc.code), dict(exc.headers), exc.read().decode("utf-8", errors="ignore")


def get_token(context: str) -> str:
    state = urllib.parse.quote(context, safe="")
    auth_url = (
        f"{BASE}/aspera/faspex/auth/authorize_public_link?response_type=code"
        f"&client_id={CLIENT_ID}"
        f"&redirect_uri=https%3A%2F%2Ffaspex.cancerimagingarchive.net%2Faspera%2Ffaspex%2Ftoken"
        f"&state={state}"
    )
    status, headers, _ = read_response(urllib.request.Request(auth_url), urllib.request.build_opener(NoRedirect))
    location = headers.get("location") or headers.get("Location") or ""
    if status not in {301, 302, 303, 307, 308} or not location:
        raise RuntimeError(f"authorization failed: status={status} location={location}")
    parsed = urllib.parse.urlparse(location)
    query = urllib.parse.parse_qs(parsed.query)
    payload = {
        "code": urllib.parse.unquote(query["code"][0]),
        "state": query["state"][0],
        "grant_type": "authorization_code",
        "client_id": CLIENT_ID,
        "redirect_uri": f"{BASE}/aspera/faspex/token",
    }
    req = urllib.request.Request(
        f"{BASE}/aspera/faspex/auth/token",
        data=json.dumps(payload).encode("utf-8"),
        headers={"Content-Type": "application/json", "Accept": "application/json"},
        method="POST",
    )
    status, _, body = read_response(req)
    if status >= 300:
        raise RuntimeError(f"token exchange failed: {status} {body[:500]}")
    return json.loads(body)["access_token"]


def api_get(path: str, token: str) -> tuple[int, str]:
    req = urllib.request.Request(
        f"{BASE}{path}",
        headers={"Accept": "application/json", "Authorization": token},
    )
    status, _, body = read_response(req)
    return status, body


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--context", required=True)
    parser.add_argument("--package-id", required=True)
    parser.add_argument("--outdir", type=Path, default=ROOT / "data" / "tcia_tcga_glioma_mri" / "aspera")
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    token = get_token(args.context)
    paths = [
        f"/aspera/faspex/api/v5/packages/{args.package_id}",
        f"/aspera/faspex/api/v5/packages/{args.package_id}/files",
        f"/aspera/faspex/api/v5/packages/{args.package_id}/files?offset=0&limit=100",
        f"/aspera/faspex/api/v5/packages/{args.package_id}/files/package_recipient",
        f"/aspera/faspex/api/v5/packages/{args.package_id}/files/package_recipient?offset=0&limit=100",
        f"/aspera/faspex/api/v5/packages/{args.package_id}/files/package_sender",
        f"/aspera/faspex/api/v5/packages/{args.package_id}/files/package_sender?offset=0&limit=100",
        f"/aspera/faspex/api/v5/packages/{args.package_id}/files/inbox",
        f"/aspera/faspex/api/v5/packages/{args.package_id}/files/inbox?offset=0&limit=100",
        f"/aspera/faspex/api/v5/packages/{args.package_id}/download_details",
        f"/aspera/faspex/api/v5/packages/{args.package_id}/transfer_spec/download?transfer_type=http_gateway",
        f"/aspera/faspex/api/v5/packages/{args.package_id}/transfer_spec/download?transfer_type=http_gateway&type=package_recipient",
        f"/aspera/faspex/api/v5/packages/{args.package_id}/transfer_spec/download?transfer_type=http_gateway&type=package_sender",
        f"/aspera/faspex/api/v5/packages/{args.package_id}/transfer_spec/download?transfer_type=http_gateway&type=inbox",
        f"/aspera/faspex/api/v5/packages/{args.package_id}/transfer_spec/download?transfer_type=connect&type=inbox",
    ]
    report = []
    for i, path in enumerate(paths):
        status, body = api_get(path, token)
        item = {"path": path, "status": status, "body_preview": body[:500]}
        report.append(item)
        if status < 300:
            suffix = "package" if i == 0 else f"endpoint_{i}"
            (args.outdir / f"package_{args.package_id}_{suffix}.json").write_text(body, encoding="utf-8")
    (args.outdir / f"package_{args.package_id}_api_probe.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(json.dumps(report, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
