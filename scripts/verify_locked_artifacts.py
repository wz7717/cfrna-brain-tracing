from __future__ import annotations

import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from core.model_lock import verify_locked_model_bundle  # noqa: E402


def main() -> None:
    manifest = verify_locked_model_bundle()
    print(
        f"Verified {len(manifest['artifacts'])} locked artifacts "
        f"for {manifest['lock_id']}"
    )


if __name__ == "__main__":
    main()
