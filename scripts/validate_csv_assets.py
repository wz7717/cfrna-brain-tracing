from __future__ import annotations

import csv
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def tracked_csv_paths() -> list[Path]:
    result = subprocess.run(
        ["git", "ls-files", "*.csv"],
        cwd=ROOT,
        check=True,
        capture_output=True,
        text=True,
        encoding="utf-8",
    )
    return [ROOT / line for line in result.stdout.splitlines() if line.strip()]


def validate_csv(path: Path) -> int:
    if not path.is_file() or path.stat().st_size == 0:
        raise ValueError("file is missing or empty")

    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle)
        try:
            header = next(reader)
        except StopIteration as exc:
            raise ValueError("file has no header") from exc

        if not header or any(not column.strip() for column in header):
            raise ValueError("header contains a blank column name")
        if len(set(header)) != len(header):
            raise ValueError("header contains duplicate column names")

        width = len(header)
        row_count = 0
        for row_number, row in enumerate(reader, start=2):
            if len(row) != width:
                raise ValueError(
                    f"row {row_number} has {len(row)} columns; expected {width}"
                )
            if any("\x00" in value for value in row):
                raise ValueError(f"row {row_number} contains a NUL byte")
            row_count += 1

    if row_count == 0:
        raise ValueError("file contains a header but no data rows")
    return row_count


def main() -> None:
    paths = tracked_csv_paths()
    if not paths:
        raise SystemExit("No tracked CSV files found")

    failures: list[str] = []
    total_rows = 0
    for path in paths:
        try:
            total_rows += validate_csv(path)
        except (OSError, UnicodeError, csv.Error, ValueError) as exc:
            failures.append(f"{path.relative_to(ROOT)}: {exc}")

    if failures:
        raise SystemExit("CSV validation failed:\n" + "\n".join(failures))
    print(f"Validated {len(paths)} tracked CSV files ({total_rows} data rows)")


if __name__ == "__main__":
    main()
