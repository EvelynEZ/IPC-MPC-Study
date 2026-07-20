"""Count NIS discharges for Asian/Pacific Islander females under age 11.

NIS is discharge-level data, so this analysis does not represent unique patients.
Results are cached against a lightweight dataset fingerprint for fast reruns.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from datetime import datetime, timezone
from pathlib import Path

import duckdb
from dotenv import load_dotenv


ANALYSIS_VERSION = "1"
FEMALE_CODE = 1
ASIAN_PACIFIC_ISLANDER_RACE_CODE = 4
MAX_AGE_EXCLUSIVE = 11


def dataset_files() -> list[Path]:
    load_dotenv()
    configured_path = os.getenv("NIS_DATASET_PATH")
    if not configured_path:
        raise SystemExit("Set NIS_DATASET_PATH in .env (see .env.example).")

    dataset_path = Path(configured_path).expanduser().resolve()
    files = sorted(dataset_path.glob("*.parquet"))
    if not files:
        raise SystemExit(f"No Parquet files found in: {dataset_path}")
    return files


def fingerprint(files: list[Path]) -> str:
    """Hash file metadata, avoiding an expensive hash of the 2+ GiB dataset."""
    digest = hashlib.sha256(ANALYSIS_VERSION.encode())
    for path in files:
        stat = path.stat()
        digest.update(f"{path.name}:{stat.st_size}:{stat.st_mtime_ns}\n".encode())
    return digest.hexdigest()


def run_query(files: list[Path]) -> dict[str, object]:
    parquet_glob = str(files[0].parent / "*.parquet")
    connection = duckdb.connect(database=":memory:")
    result = connection.execute(
        """
        SELECT
            count(*)::BIGINT AS discharge_count,
            sum(DISCWT)::DOUBLE AS weighted_national_estimate
        FROM read_parquet(?)
        WHERE FEMALE = ?
          AND RACE = ?
          AND AGE < ?
        """,
        [
            parquet_glob,
            FEMALE_CODE,
            ASIAN_PACIFIC_ISLANDER_RACE_CODE,
            MAX_AGE_EXCLUSIVE,
        ],
    ).fetchone()
    connection.close()
    return {
        "discharge_count": result[0],
        "weighted_national_estimate": result[1],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--refresh", action="store_true", help="Ignore a valid cached result."
    )
    args = parser.parse_args()

    files = dataset_files()
    dataset_fingerprint = fingerprint(files)
    cache_path = Path("outputs/cache/young_asian_female_discharges.json")

    if cache_path.exists() and not args.refresh:
        cached = json.loads(cache_path.read_text())
        if cached.get("dataset_fingerprint") == dataset_fingerprint:
            cached["cache_used"] = True
            print(json.dumps(cached, indent=2))
            return

    counts = run_query(files)
    result = {
        "analysis": "Asian/Pacific Islander female discharges, ages 0 through 10",
        "unit": "NIS inpatient discharge (not a unique patient)",
        **counts,
        "filters": {
            "FEMALE": FEMALE_CODE,
            "RACE": ASIAN_PACIFIC_ISLANDER_RACE_CODE,
            "AGE": f"< {MAX_AGE_EXCLUSIVE}",
        },
        "dataset_file_count": len(files),
        "dataset_fingerprint": dataset_fingerprint,
        "computed_at_utc": datetime.now(timezone.utc).isoformat(),
        "cache_used": False,
    }
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    cache_path.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
