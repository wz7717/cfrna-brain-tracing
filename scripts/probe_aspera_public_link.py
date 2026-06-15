#!/usr/bin/env python
from __future__ import annotations

import argparse
import json
import urllib.error
import urllib.parse
import urllib.request


BASE = "https://faspex.cancerimagingarchive.net"
CLIENT_ID = "ff9aa63a-72e1-436f-82ef-5677eb1f7aee"


def request(url: str, headers: dict[str, str] | None = None, data: bytes | None = None) -> tuple[int, dict[str, str], str]:
    req = urllib.request.Request(url, headers=headers or {}, data=data)
    opener = urllib.request.build_opener(urllib.request.HTTPRedirectHandler)
    try:
        with opener.open(req, timeout=60) as resp:
            body = resp.read().decode("utf-8", errors="ignore")
            return int(resp.status), dict(resp.headers), body
    except urllib.error.HTTPError as exc:
        body = exc.read().decode("utf-8", errors="ignore")
        return int(exc.code), dict(exc.headers), body


def no_redirect(url: str) -> tuple[int, dict[str, str], str]:
    class NoRedirect(urllib.request.HTTPRedirectHandler):
        def redirect_request(self, req, fp, code, msg, headers, newurl):  # type: ignore[no-untyped-def]
            return None

    opener = urllib.request.build_opener(NoRedirect)
    req = urllib.request.Request(url)
    try:
        with opener.open(req, timeout=60) as resp:
            return int(resp.status), dict(resp.headers), resp.read().decode("utf-8", errors="ignore")
    except urllib.error.HTTPError as exc:
        return int(exc.code), dict(exc.headers), exc.read().decode("utf-8", errors="ignore")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("context")
    args = parser.parse_args()
    state = urllib.parse.quote(args.context, safe="")
    paths = [
        f"/aspera/faspex/auth/authorize_public_link?response_type=code&client_id={CLIENT_ID}&redirect_uri=/aspera/faspex/token&state={state}",
        f"/aspera/faspex/auth/authorize_public_link?response_type=code&client_id={CLIENT_ID}&redirect_uri=https%3A%2F%2Ffaspex.cancerimagingarchive.net%2Faspera%2Ffaspex%2Ftoken&state={state}",
    ]
    token: str | None = None
    for path in paths:
        url = BASE + path
        status, headers, body = no_redirect(url)
        print("URL", url)
        print("STATUS", status)
        print("LOCATION", headers.get("Location", ""))
        print("HEADERS", json.dumps(headers, ensure_ascii=False, indent=2)[:2000])
        print("BODY", body[:1000])
        location = headers.get("location") or headers.get("Location") or ""
        if not location:
            continue
        parsed = urllib.parse.urlparse(location)
        code = urllib.parse.parse_qs(parsed.query).get("code", [""])[0]
        if not code:
            continue
        payload = {
            "code": urllib.parse.unquote(code),
            "state": urllib.parse.parse_qs(parsed.query).get("state", [""])[0],
            "grant_type": "authorization_code",
            "client_id": CLIENT_ID,
            "redirect_uri": f"{BASE}/aspera/faspex/token",
        }
        data = json.dumps(payload).encode("utf-8")
        st, hdr, txt = request(
            BASE + "/aspera/faspex/auth/token",
            headers={"Content-Type": "application/json", "Accept": "application/json"},
            data=data,
        )
        print("TOKEN_STATUS", st)
        print("TOKEN_BODY", txt[:1200])
        if st < 300:
            parsed_token = json.loads(txt)
            token = parsed_token.get("access_token") or parsed_token.get("token")
            break
    if token:
        headers = {"Accept": "application/json", "Authorization": token}
        api_paths = [
            "/aspera/faspex/api/v5/packages/576",
            "/aspera/faspex/api/v5/packages/576/files/received",
            "/aspera/faspex/api/v5/packages/576/files/sent",
            "/aspera/faspex/api/v5/packages/576/transfer_spec/download?transfer_type=http_gateway&type=received",
            "/aspera/faspex/api/v5/packages/576/transfer_spec/download?transfer_type=http_gateway&type=sent",
        ]
        for api_path in api_paths:
            st, hdr, txt = request(BASE + api_path, headers=headers)
            print("\nAPI", api_path, "STATUS", st)
            print(txt[:3000])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
