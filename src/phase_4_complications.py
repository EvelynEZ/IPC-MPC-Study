"""Compare co-documented complications in HM discharges by sepsis status."""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any

import duckdb


REPO_ROOT = Path(__file__).resolve().parents[1]
COHORT_DATABASE = REPO_ROOT / "outputs/phase_1_2/hm_cohort.duckdb"
OUTPUT_DIR = REPO_ROOT / "outputs/phase_4"

COMPLICATIONS = [
    {
        "id": "septic_shock",
        "label": "Severe sepsis with septic shock",
        "prefixes": ["R6521"],
    },
    {
        "id": "severe_sepsis",
        "label": "Severe sepsis, with or without shock",
        "prefixes": ["R652"],
    },
    {
        "id": "acute_respiratory_failure",
        "label": "Acute respiratory failure",
        "prefixes": ["J960"],
    },
    {
        "id": "respiratory_failure_unspecified",
        "label": "Respiratory failure, unspecified",
        "prefixes": ["J969"],
    },
    {"id": "aki", "label": "Acute kidney injury", "prefixes": ["N17"]},
    {
        "id": "pneumonia",
        "label": "Pneumonia",
        "prefixes": ["J12", "J13", "J14", "J15", "J16", "J17", "J18"],
    },
    {
        "id": "pulmonary_embolism",
        "label": "Pulmonary embolism",
        "exact": ["I2602", "I2609", "I2692", "I2693", "I2694", "I2699"],
    },
    {
        "id": "acute_lower_extremity_dvt",
        "label": "Acute lower-extremity DVT",
        "exact": [
            "I82401", "I82402", "I82403", "I82409", "I82491", "I82492",
            "I82493", "I82499", "I824Y1", "I824Y2", "I824Y3", "I824Y9",
            "I824Z1", "I824Z2", "I824Z3", "I824Z9"
        ],
    },
    {"id": "tls", "label": "Tumor lysis syndrome", "exact": ["E883"]},
    {"id": "dic", "label": "Disseminated intravascular coagulation", "prefixes": ["D65"]},
    {
        "id": "paroxysmal_tachycardia",
        "label": "Paroxysmal tachycardia/SVT/VT",
        "prefixes": ["I47"],
    },
    {
        "id": "afib_flutter",
        "label": "Atrial fibrillation/flutter, excluding chronic AF",
        "prefixes": ["I48"],
        "exclude_prefixes": ["I482"],
    },
    {
        "id": "other_arrhythmia",
        "label": "Other cardiac arrhythmias",
        "prefixes": ["I49"],
    },
]


def sql_string(value: str) -> str:
    return "'" + value.replace("'", "''") + "'"


def code_condition(rule: dict[str, Any], code: str = "code") -> str:
    included = [
        f"starts_with({code}, {sql_string(prefix)})"
        for prefix in rule.get("prefixes", [])
    ]
    included.extend(
        f"{code} = {sql_string(exact)}" for exact in rule.get("exact", [])
    )
    condition = "(" + " OR ".join(included) + ")"
    excluded = rule.get("exclude_prefixes", [])
    if excluded:
        exclusion = " OR ".join(
            f"starts_with({code}, {sql_string(prefix)})" for prefix in excluded
        )
        condition += f" AND NOT ({exclusion})"
    return f"({condition})"


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main() -> dict[str, Any]:
    if not COHORT_DATABASE.exists():
        raise RuntimeError("Run Phase 1–2 first to create the cached HM cohort.")
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    connection = duckdb.connect(str(COHORT_DATABASE), read_only=True)

    flag_sql = []
    for rule in COMPLICATIONS:
        condition = code_condition(rule)
        flag_sql.append(
            f"list_contains(list_transform(diagnosis_codes, code -> {condition}), TRUE) "
            f"AS {rule['id']}"
        )
    connection.execute(
        "CREATE TEMP TABLE complication_flags AS "
        "SELECT sepsis, DISCWT, " + ", ".join(flag_sql) + " FROM hm_cohort"
    )
    connection.execute(
        """
        CREATE TEMP TABLE complication_flags_with_any AS
        SELECT *,
            (paroxysmal_tachycardia OR afib_flutter OR other_arrhythmia) AS any_arrhythmia
        FROM complication_flags
        """
    )

    rules = COMPLICATIONS + [
        {"id": "any_arrhythmia", "label": "Any specified cardiac arrhythmia"}
    ]
    totals = {
        bool(row[0]): {"unweighted": row[1], "weighted": row[2]}
        for row in connection.execute(
            "SELECT sepsis, count(*)::BIGINT, sum(DISCWT)::DOUBLE "
            "FROM complication_flags_with_any GROUP BY sepsis"
        ).fetchall()
    }

    table_rows: list[dict[str, Any]] = []
    for rule in rules:
        result = {
            bool(row[0]): {"unweighted": row[1], "weighted": row[2] or 0.0}
            for row in connection.execute(
                f"""
                SELECT sepsis,
                    count(*) FILTER (WHERE {rule['id']})::BIGINT,
                    sum(DISCWT) FILTER (WHERE {rule['id']})::DOUBLE
                FROM complication_flags_with_any
                GROUP BY sepsis
                """
            ).fetchall()
        }
        p0 = result[False]["weighted"] / totals[False]["weighted"]
        p1 = result[True]["weighted"] / totals[True]["weighted"]
        table_rows.append(
            {
                "complication": rule["label"],
                "no_sepsis_unweighted_n": result[False]["unweighted"],
                "no_sepsis_weighted_n": round(result[False]["weighted"]),
                "no_sepsis_weighted_percent": round(100 * p0, 2),
                "sepsis_unweighted_n": result[True]["unweighted"],
                "sepsis_weighted_n": round(result[True]["weighted"]),
                "sepsis_weighted_percent": round(100 * p1, 2),
                "absolute_difference_percentage_points": round(100 * (p1 - p0), 2),
                "prevalence_ratio": round(p1 / p0, 2) if p0 else None,
            }
        )
    connection.close()
    write_csv(OUTPUT_DIR / "complications_by_sepsis.csv", table_rows)

    summary = {
        "unit": "NIS inpatient discharge, not unique patient",
        "exposure": "Documented sepsis defined by A41* in any diagnosis position",
        "interpretation": "Complications are co-documented during the same hospitalization; temporal order is not available.",
        "no_sepsis_unweighted_n": totals[False]["unweighted"],
        "sepsis_unweighted_n": totals[True]["unweighted"],
        "largest_absolute_differences": sorted(
            table_rows,
            key=lambda row: abs(row["absolute_difference_percentage_points"]),
            reverse=True,
        )[:5],
    }
    (OUTPUT_DIR / "phase_4_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    return summary


if __name__ == "__main__":
    main()
