"""Phase 1 data audit and Phase 2 cached cohort construction."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import duckdb
from dotenv import load_dotenv


REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CONFIG = REPO_ROOT / "config/hm_phenotype_v0_1.json"
OUTPUT_DIR = REPO_ROOT / "outputs/phase_1_2"
DATABASE_PATH = OUTPUT_DIR / "hm_cohort.duckdb"
ANALYSIS_VERSION = "2-stratum-derived-hospital-characteristics"


def load_dataset_files() -> list[Path]:
    load_dotenv(REPO_ROOT / ".env")
    configured_path = os.getenv("NIS_DATASET_PATH")
    if not configured_path:
        raise RuntimeError("Set NIS_DATASET_PATH in .env (see .env.example).")
    files = sorted(Path(configured_path).expanduser().resolve().glob("*.parquet"))
    if not files:
        raise RuntimeError(f"No Parquet files found under {configured_path}")
    return files


def dataset_fingerprint(files: list[Path], config_bytes: bytes) -> str:
    digest = hashlib.sha256(ANALYSIS_VERSION.encode() + config_bytes)
    for path in files:
        stat = path.stat()
        digest.update(f"{path.name}:{stat.st_size}:{stat.st_mtime_ns}\n".encode())
    return digest.hexdigest()


def normalize_code_expression(column: str) -> str:
    return (
        f"NULLIF(upper(replace(replace(trim(CAST({column} AS VARCHAR)), '.', ''), "
        "' ', '')), '')"
    )


def sql_string(value: str) -> str:
    return "'" + value.replace("'", "''") + "'"


def code_match_sql(code: str, rule: dict[str, Any]) -> str:
    clauses = [
        f"starts_with({code}, {sql_string(prefix)})"
        for prefix in rule.get("include_prefixes", [])
    ]
    clauses.extend(
        f"{code} = {sql_string(exact)}" for exact in rule.get("include_exact", [])
    )
    if not clauses:
        raise ValueError(f"Rule has no included codes: {rule}")
    expression = "(" + " OR ".join(clauses) + ")"
    exclusions = rule.get("exclude_exact", [])
    if exclusions:
        exclusion_sql = ", ".join(sql_string(value) for value in exclusions)
        expression += f" AND {code} NOT IN ({exclusion_sql})"
    return f"({expression})"


def subtype_case_sql(config: dict[str, Any], code: str = "code") -> str:
    cases = [
        f"WHEN {code_match_sql(code, subtype)} THEN {sql_string(subtype['id'])}"
        for subtype in config["subtypes"]
    ]
    return "CASE " + " ".join(cases) + " ELSE NULL END"


def get_columns(connection: duckdb.DuckDBPyConnection, parquet_glob: str) -> list[str]:
    rows = connection.execute(
        "DESCRIBE SELECT * FROM read_parquet(?)", [parquet_glob]
    ).fetchall()
    return [row[0] for row in rows]


def audit_rows(columns: list[str]) -> list[dict[str, str]]:
    column_set = set(columns)
    dx_columns = sorted(
        (name for name in columns if name.startswith("I10_DX") and name[6:].isdigit()),
        key=lambda name: int(name[6:]),
    )
    def row(component: str, available: bool, detail: str) -> dict[str, str]:
        return {
            "component": component,
            "status": "READY" if available else "MISSING/BLOCKED",
            "detail": detail,
        }

    rows = [
        row("Diagnosis fields", bool(dx_columns), f"{len(dx_columns)} fields available"),
        row("Adult cohort", "AGE" in column_set, "Requires AGE"),
        row("Discharge weights", "DISCWT" in column_set, "Requires DISCWT"),
        row("Survey strata", "NIS_STRATUM" in column_set, "Requires NIS_STRATUM"),
        row("Mortality", "DIED" in column_set, "Requires DIED"),
        row("Length of stay", "LOS" in column_set, "Requires LOS"),
    ]
    rows.append(
        {
            "component": "Hospital characteristics",
            "status": "DERIVED/REVIEW" if "NIS_STRATUM" in column_set else "MISSING/BLOCKED",
            "detail": "Division, region, control, location/teaching, and bed size decoded from four-digit NIS_STRATUM",
        }
    )
    rows.append(
        {
            "component": "Mechanical ventilation",
            "status": "OUT OF SCOPE",
            "detail": "Removed from the requested analysis; procedure fields are not required",
        }
    )
    return rows


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        return
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def create_cohort(
    connection: duckdb.DuckDBPyConnection,
    parquet_glob: str,
    columns: list[str],
    config: dict[str, Any],
) -> None:
    dx_columns = sorted(
        (name for name in columns if name.startswith("I10_DX") and name[6:].isdigit()),
        key=lambda name: int(name[6:]),
    )
    normalized_list = ", ".join(normalize_code_expression(name) for name in dx_columns)
    retained = [name for name in columns if name not in dx_columns]
    retained_sql = ", ".join(f'"{name}"' for name in retained)
    subtype_case = subtype_case_sql(config)
    sepsis_match = code_match_sql("code", config["sepsis"])
    palliative_match = code_match_sql("code", config["palliative_care"])
    minimum_age = int(config["minimum_age"])

    connection.execute("DROP TABLE IF EXISTS hm_cohort")
    connection.execute(
        f"""
        CREATE TABLE hm_cohort AS
        WITH normalized AS (
            SELECT
                {retained_sql},
                [{normalized_list}] AS diagnosis_codes
            FROM read_parquet(?)
            WHERE AGE >= {minimum_age}
        ),
        classified AS (
            SELECT
                *,
                list_transform(diagnosis_codes, code -> {subtype_case}) AS subtype_by_position,
                list_contains(list_transform(diagnosis_codes, code -> {sepsis_match}), TRUE) AS sepsis,
                list_contains(list_transform(diagnosis_codes, code -> {palliative_match}), TRUE) AS palliative_care
            FROM normalized
        ),
        derived AS (
            SELECT
                * EXCLUDE (subtype_by_position),
                list_distinct(list_filter(subtype_by_position, item -> item IS NOT NULL)) AS hm_groups,
                list_first(list_filter(subtype_by_position, item -> item IS NOT NULL)) AS hm_subtype
            FROM classified
        )
        SELECT
            *,
            len(hm_groups)::UTINYINT AS hm_group_count,
            (len(hm_groups) > 1) AS multiple_hm_groups,
            CAST(floor(NIS_STRATUM / 1000) AS UTINYINT) AS hospital_division_code,
            CASE
                WHEN floor(NIS_STRATUM / 1000) IN (1, 2) THEN 'Northeast'
                WHEN floor(NIS_STRATUM / 1000) IN (3, 4) THEN 'Midwest'
                WHEN floor(NIS_STRATUM / 1000) IN (5, 6, 7) THEN 'South'
                WHEN floor(NIS_STRATUM / 1000) IN (8, 9) THEN 'West'
                ELSE 'Unknown'
            END AS hospital_region,
            CASE CAST(floor(NIS_STRATUM / 100) % 10 AS INTEGER)
                WHEN 0 THEN 'Government or private (collapsed)'
                WHEN 1 THEN 'Government, nonfederal'
                WHEN 2 THEN 'Private, not-for-profit'
                WHEN 3 THEN 'Private, investor-owned'
                WHEN 4 THEN 'Private, type collapsed'
                ELSE 'Unknown'
            END AS hospital_control,
            CASE CAST(floor(NIS_STRATUM / 10) % 10 AS INTEGER)
                WHEN 1 THEN 'Rural'
                WHEN 2 THEN 'Urban nonteaching'
                WHEN 3 THEN 'Urban teaching'
                ELSE 'Unknown'
            END AS hospital_location_teaching,
            CASE CAST(NIS_STRATUM % 10 AS INTEGER)
                WHEN 1 THEN 'Small'
                WHEN 2 THEN 'Medium'
                WHEN 3 THEN 'Large'
                ELSE 'Unknown'
            END AS hospital_bed_size
        FROM derived
        WHERE len(hm_groups) >= 1
        """,
        [parquet_glob],
    )
    connection.execute("CREATE INDEX hm_cohort_year_idx ON hm_cohort(YEAR)")
    connection.execute("CREATE INDEX hm_cohort_subtype_idx ON hm_cohort(hm_subtype)")


def query_dicts(
    connection: duckdb.DuckDBPyConnection, sql: str
) -> list[dict[str, Any]]:
    cursor = connection.execute(sql)
    names = [item[0] for item in cursor.description]
    return [dict(zip(names, row)) for row in cursor.fetchall()]


def export_summaries(
    connection: duckdb.DuckDBPyConnection, config: dict[str, Any]
) -> dict[str, Any]:
    overall = query_dicts(
        connection,
        """
        SELECT
            count(*)::BIGINT AS unweighted_hm_discharges,
            round(sum(DISCWT), 0)::DOUBLE AS weighted_hm_discharges_2016_2022,
            count(*) FILTER (WHERE sepsis)::BIGINT AS unweighted_with_sepsis,
            round(sum(DISCWT * sepsis::INTEGER), 0)::DOUBLE AS weighted_with_sepsis,
            round(100.0 * avg(sepsis::INTEGER), 2) AS unweighted_sepsis_percent,
            round(100.0 * sum(DISCWT * sepsis::INTEGER) / sum(DISCWT), 2) AS weighted_sepsis_percent,
            count(*) FILTER (WHERE palliative_care)::BIGINT AS unweighted_with_palliative_care,
            round(sum(DISCWT * palliative_care::INTEGER), 0)::DOUBLE AS weighted_with_palliative_care,
            round(100.0 * sum(DISCWT * palliative_care::INTEGER) / sum(DISCWT), 2) AS weighted_palliative_care_percent,
            count(*) FILTER (WHERE multiple_hm_groups)::BIGINT AS unweighted_with_multiple_hm_groups,
            round(sum(DISCWT * multiple_hm_groups::INTEGER), 0)::DOUBLE AS weighted_with_multiple_hm_groups,
            round(100.0 * sum(DISCWT * multiple_hm_groups::INTEGER) / sum(DISCWT), 2) AS weighted_multiple_hm_groups_percent
        FROM hm_cohort
        """,
    )[0]
    by_year = query_dicts(
        connection,
        """
        SELECT
            YEAR::INTEGER AS year,
            count(*)::BIGINT AS unweighted_hm_discharges,
            round(sum(DISCWT), 0)::DOUBLE AS weighted_hm_discharges,
            count(*) FILTER (WHERE sepsis)::BIGINT AS unweighted_with_sepsis,
            round(100.0 * sum(DISCWT * sepsis::INTEGER) / sum(DISCWT), 2) AS weighted_sepsis_percent,
            count(*) FILTER (WHERE palliative_care)::BIGINT AS unweighted_with_palliative_care,
            round(100.0 * sum(DISCWT * palliative_care::INTEGER) / sum(DISCWT), 2) AS weighted_palliative_care_percent
        FROM hm_cohort
        GROUP BY YEAR
        ORDER BY YEAR
        """,
    )
    by_subtype = query_dicts(
        connection,
        """
        SELECT
            hm_subtype,
            count(*)::BIGINT AS unweighted_hm_discharges,
            round(sum(DISCWT), 0)::DOUBLE AS weighted_hm_discharges_2016_2022,
            count(*) FILTER (WHERE sepsis)::BIGINT AS unweighted_with_sepsis,
            round(100.0 * sum(DISCWT * sepsis::INTEGER) / sum(DISCWT), 2) AS weighted_sepsis_percent,
            count(*) FILTER (WHERE palliative_care)::BIGINT AS unweighted_with_palliative_care,
            round(100.0 * sum(DISCWT * palliative_care::INTEGER) / sum(DISCWT), 2) AS weighted_palliative_care_percent
        FROM hm_cohort
        GROUP BY hm_subtype
        ORDER BY unweighted_hm_discharges DESC
        """,
    )
    overlap = query_dicts(
        connection,
        """
        SELECT
            CASE WHEN hm_group_count >= 4 THEN '4+' ELSE CAST(hm_group_count AS VARCHAR) END AS hm_group_count,
            count(*)::BIGINT AS unweighted_hm_discharges,
            round(sum(DISCWT), 0)::DOUBLE AS weighted_hm_discharges_2016_2022,
            round(100.0 * sum(DISCWT) / sum(sum(DISCWT)) OVER (), 2) AS weighted_percent
        FROM hm_cohort
        GROUP BY CASE WHEN hm_group_count >= 4 THEN '4+' ELSE CAST(hm_group_count AS VARCHAR) END
        ORDER BY min(hm_group_count)
        """,
    )
    hospital_dimensions = {
        "hospital_by_region.csv": "hospital_region",
        "hospital_by_division.csv": "hospital_division_code",
        "hospital_by_location_teaching.csv": "hospital_location_teaching",
        "hospital_by_bed_size.csv": "hospital_bed_size",
        "hospital_by_control.csv": "hospital_control",
    }
    for filename, variable in hospital_dimensions.items():
        rows = query_dicts(
            connection,
            f"""
            SELECT
                {variable},
                count(*)::BIGINT AS unweighted_hm_discharges,
                round(sum(DISCWT), 0)::DOUBLE AS weighted_hm_discharges_2016_2022,
                round(100.0 * sum(DISCWT) / sum(sum(DISCWT)) OVER (), 2) AS weighted_percent
            FROM hm_cohort
            GROUP BY {variable}
            ORDER BY weighted_hm_discharges_2016_2022 DESC
            """,
        )
        write_csv(OUTPUT_DIR / filename, rows)
    stratum_validation = query_dicts(
        connection,
        """
        SELECT
            count(*)::BIGINT AS unweighted_hm_discharges,
            count(*) FILTER (
                WHERE hospital_division_code NOT BETWEEN 1 AND 9
                   OR hospital_region = 'Unknown'
                   OR hospital_control = 'Unknown'
                   OR hospital_location_teaching = 'Unknown'
                   OR hospital_bed_size = 'Unknown'
            )::BIGINT AS invalid_or_unknown_decodes
        FROM hm_cohort
        """,
    )
    write_csv(OUTPUT_DIR / "cohort_by_year.csv", by_year)
    write_csv(OUTPUT_DIR / "cohort_by_subtype.csv", by_subtype)
    write_csv(OUTPUT_DIR / "cohort_overlap.csv", overlap)
    write_csv(OUTPUT_DIR / "stratum_decode_validation.csv", stratum_validation)
    review_overview = [
        {
            "cohort_measure": "All adult HM discharges",
            "unweighted_n": overall["unweighted_hm_discharges"],
            "weighted_n_2016_2022": overall["weighted_hm_discharges_2016_2022"],
            "weighted_percent": 100.0,
        },
        {
            "cohort_measure": "Documented sepsis (A41*)",
            "unweighted_n": overall["unweighted_with_sepsis"],
            "weighted_n_2016_2022": overall["weighted_with_sepsis"],
            "weighted_percent": overall["weighted_sepsis_percent"],
        },
        {
            "cohort_measure": "Documented inpatient palliative-care use (Z51.5)",
            "unweighted_n": overall["unweighted_with_palliative_care"],
            "weighted_n_2016_2022": overall["weighted_with_palliative_care"],
            "weighted_percent": overall["weighted_palliative_care_percent"],
        },
        {
            "cohort_measure": "Multiple HM groups",
            "unweighted_n": overall["unweighted_with_multiple_hm_groups"],
            "weighted_n_2016_2022": overall["weighted_with_multiple_hm_groups"],
            "weighted_percent": overall["weighted_multiple_hm_groups_percent"],
        },
    ]
    write_csv(OUTPUT_DIR / "review_table_1_cohort_overview.csv", review_overview)

    subtype_labels = {item["id"]: item["label"] for item in config["subtypes"]}
    review_subtypes = [
        {
            "hm_subtype": subtype_labels.get(row["hm_subtype"], row["hm_subtype"]),
            "unweighted_n": row["unweighted_hm_discharges"],
            "weighted_n_2016_2022": row["weighted_hm_discharges_2016_2022"],
            "sepsis_unweighted_n": row["unweighted_with_sepsis"],
            "sepsis_weighted_percent": row["weighted_sepsis_percent"],
            "palliative_care_unweighted_n": row[
                "unweighted_with_palliative_care"
            ],
            "palliative_care_weighted_percent": row[
                "weighted_palliative_care_percent"
            ],
        }
        for row in by_subtype
    ]
    write_csv(OUTPUT_DIR / "review_table_3_subtypes.csv", review_subtypes)
    (OUTPUT_DIR / "cohort_summary.json").write_text(
        json.dumps(overall, indent=2) + "\n"
    )
    (OUTPUT_DIR / "phenotype_review.json").write_text(
        json.dumps(
            {
                "version": config["version"],
                "review_status": config["review_status"],
                "decisions_required_before_freeze": config[
                    "decisions_required_before_freeze"
                ],
            },
            indent=2,
        )
        + "\n"
    )
    return overall


def main(argv: list[str] | None = None) -> dict[str, Any]:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--refresh", action="store_true", help="Rebuild the cohort")
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    args = parser.parse_args(argv)

    config_bytes = args.config.read_bytes()
    config = json.loads(config_bytes)
    files = load_dataset_files()
    fingerprint = dataset_fingerprint(files, config_bytes)
    parquet_glob = str(files[0].parent / "*.parquet")
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    connection = duckdb.connect(str(DATABASE_PATH))
    columns = get_columns(connection, parquet_glob)
    audit = audit_rows(columns)
    write_csv(OUTPUT_DIR / "data_readiness_audit.csv", audit)

    cached_fingerprint = None
    if "analysis_metadata" in {
        row[0] for row in connection.execute("SHOW TABLES").fetchall()
    }:
        cached_fingerprint = connection.execute(
            "SELECT value FROM analysis_metadata WHERE key = 'fingerprint'"
        ).fetchone()
        cached_fingerprint = cached_fingerprint[0] if cached_fingerprint else None

    cache_used = cached_fingerprint == fingerprint and not args.refresh
    if not cache_used:
        create_cohort(connection, parquet_glob, columns, config)
        connection.execute("DROP TABLE IF EXISTS analysis_metadata")
        connection.execute("CREATE TABLE analysis_metadata(key VARCHAR, value VARCHAR)")
        connection.execute(
            "INSERT INTO analysis_metadata VALUES ('fingerprint', ?), ('built_at_utc', ?)",
            [fingerprint, datetime.now(timezone.utc).isoformat()],
        )

    overall = export_summaries(connection, config)
    connection.close()
    result = {
        "phenotype_version": config["version"],
        "phenotype_status": config["review_status"],
        "dataset_files": len(files),
        "cache_used": cache_used,
        **overall,
        "output_directory": str(OUTPUT_DIR),
    }
    print(json.dumps(result, indent=2))
    return result


if __name__ == "__main__":
    main()
