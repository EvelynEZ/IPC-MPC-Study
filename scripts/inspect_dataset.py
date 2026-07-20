"""Validate the local NIS Parquet dataset without displaying record-level data."""

from __future__ import annotations

import os
from pathlib import Path

import duckdb
from dotenv import load_dotenv


def main() -> None:
    load_dotenv()
    configured_path = os.getenv("NIS_DATASET_PATH")
    if not configured_path:
        raise SystemExit("Set NIS_DATASET_PATH in .env (see .env.example).")

    dataset_path = Path(configured_path).expanduser().resolve()
    if not dataset_path.is_dir():
        raise SystemExit(f"Dataset directory not found: {dataset_path}")

    parquet_files = sorted(dataset_path.glob("*.parquet"))
    if not parquet_files:
        raise SystemExit(f"No .parquet files found in: {dataset_path}")

    total_bytes = sum(path.stat().st_size for path in parquet_files)
    parquet_glob = str(dataset_path / "*.parquet")

    connection = duckdb.connect(database=":memory:")
    # Parameters keep the local path separate from SQL syntax.
    schema = connection.execute(
        "DESCRIBE SELECT * FROM read_parquet(?)", [parquet_glob]
    ).fetchall()
    row_count = connection.execute(
        "SELECT count(*) FROM read_parquet(?)", [parquet_glob]
    ).fetchone()[0]

    print(f"Dataset: {dataset_path}")
    print(f"Parquet files: {len(parquet_files):,}")
    print(f"On-disk size: {total_bytes / (1024 ** 3):.2f} GiB")
    print(f"Rows: {row_count:,}")
    print(f"Columns: {len(schema):,}\n")
    print("Schema:")
    for column_name, column_type, *_ in schema:
        print(f"  {column_name}: {column_type}")


if __name__ == "__main__":
    main()
