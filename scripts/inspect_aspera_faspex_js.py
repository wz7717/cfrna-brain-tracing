#!/usr/bin/env python
from __future__ import annotations

import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_JS = ROOT / "data" / "tcia_tcga_glioma_mri" / "manifests" / "aspera_faspex_index.js"


def main() -> int:
    text = DEFAULT_JS.read_text(encoding="utf-8", errors="ignore")
    patterns = [
        r"/api/v5/[A-Za-z0-9_./?=&{}:${}\-\[\]]+",
        r"/public/[A-Za-z0-9_./?=&{}:${}\-\[\]]+",
        r"/packages/[A-Za-z0-9_./?=&{}:${}\-\[\]]+",
    ]
    seen: set[str] = set()
    for pattern in patterns:
        for match in re.findall(pattern, text):
            seen.add(match[:240])
    for item in sorted(seen):
        print(item)
    print("n", len(seen))
    for token in ["external_download_package", "package_id", "passcode", "transfer_spec", "download_setup"]:
        idx = text.find(token)
        print("\nTOKEN", token, idx)
        if idx >= 0:
            print(text[max(0, idx - 800): idx + 1400])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
