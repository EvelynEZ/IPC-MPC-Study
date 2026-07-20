"""Phase 3 weighted baseline characteristics by documented sepsis status."""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path
from typing import Any

import duckdb


REPO_ROOT = Path(__file__).resolve().parents[1]
COHORT_DATABASE = REPO_ROOT / "outputs/phase_1_2/hm_cohort.duckdb"
OUTPUT_DIR = REPO_ROOT / "outputs/phase_3"


CATEGORICAL_CHARACTERISTICS = [
    (
        "Age group",
        "CASE WHEN AGE < 60 THEN '18–59' WHEN AGE >= 60 THEN '60 or older' ELSE 'Missing' END",
        ["18–59", "60 or older", "Missing"],
    ),
    (
        "Sex",
        "CASE FEMALE WHEN 0 THEN 'Male' WHEN 1 THEN 'Female' ELSE 'Missing' END",
        ["Male", "Female", "Missing"],
    ),
    (
        "Race/ethnicity",
        """CASE RACE
            WHEN 1 THEN 'White' WHEN 2 THEN 'Black'
            WHEN 3 THEN 'Hispanic' WHEN 4 THEN 'Asian/Pacific Islander'
            WHEN 5 THEN 'Native American' WHEN 6 THEN 'Other'
            ELSE 'Missing' END""",
        [
            "White",
            "Black",
            "Hispanic",
            "Asian/Pacific Islander",
            "Native American",
            "Other",
            "Missing",
        ],
    ),
    (
        "Primary payer",
        """CASE PAY1
            WHEN 1 THEN 'Medicare' WHEN 2 THEN 'Medicaid'
            WHEN 3 THEN 'Private insurance' WHEN 4 THEN 'Self-pay'
            WHEN 5 THEN 'No charge' WHEN 6 THEN 'Other'
            ELSE 'Missing' END""",
        [
            "Medicare",
            "Medicaid",
            "Private insurance",
            "Self-pay",
            "No charge",
            "Other",
            "Missing",
        ],
    ),
    (
        "Median household-income quartile",
        """CASE ZIPINC_QRTL
            WHEN 1 THEN '0–25th percentile' WHEN 2 THEN '26th–50th percentile'
            WHEN 3 THEN '51st–75th percentile' WHEN 4 THEN '76th–100th percentile'
            ELSE 'Missing' END""",
        [
            "0–25th percentile",
            "26th–50th percentile",
            "51st–75th percentile",
            "76th–100th percentile",
            "Missing",
        ],
    ),
    (
        "HM subtype",
        """CASE hm_subtype
            WHEN 'lymphoma' THEN 'Lymphoma'
            WHEN 'aml' THEN 'Acute myeloid leukemia'
            WHEN 'cml' THEN 'Chronic myeloid leukemia'
            WHEN 'cll_chronic_leukemia' THEN 'CLL/chronic leukemia'
            WHEN 'all' THEN 'Acute lymphoblastic/unspecified acute leukemia'
            WHEN 'other_leukemia' THEN 'Other leukemia'
            WHEN 'myeloma_plasma_cell' THEN 'Multiple myeloma/plasma-cell neoplasm'
            WHEN 'mds' THEN 'Myelodysplastic disease'
            WHEN 'mpn' THEN 'Myeloproliferative neoplasm'
            ELSE 'Missing' END""",
        [
            "Lymphoma",
            "Acute myeloid leukemia",
            "Chronic myeloid leukemia",
            "CLL/chronic leukemia",
            "Acute lymphoblastic/unspecified acute leukemia",
            "Other leukemia",
            "Multiple myeloma/plasma-cell neoplasm",
            "Myelodysplastic disease",
            "Myeloproliferative neoplasm",
            "Missing",
        ],
    ),
    (
        "Hospital region (stratum-derived)",
        "hospital_region",
        ["Northeast", "Midwest", "South", "West", "Unknown"],
    ),
    (
        "Hospital location/teaching (stratum-derived)",
        "hospital_location_teaching",
        ["Rural", "Urban nonteaching", "Urban teaching", "Unknown"],
    ),
    (
        "Hospital bed size (stratum-derived)",
        "hospital_bed_size",
        ["Small", "Medium", "Large", "Unknown"],
    ),
    (
        "Hospital ownership/control (stratum-derived)",
        "hospital_control",
        [
            "Government or private (collapsed)",
            "Government, nonfederal",
            "Private, not-for-profit",
            "Private, investor-owned",
            "Private, type collapsed",
            "Unknown",
        ],
    ),
    (
        "Admission year",
        "CAST(YEAR AS INTEGER)::VARCHAR",
        [str(year) for year in range(2016, 2023)] + ["Missing"],
    ),
]


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        return
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def standardized_difference(proportion_1: float, proportion_0: float) -> float:
    pooled_variance = (
        proportion_1 * (1 - proportion_1)
        + proportion_0 * (1 - proportion_0)
    ) / 2
    if pooled_variance <= 0:
        return 0.0
    return (proportion_1 - proportion_0) / math.sqrt(pooled_variance)


def weighted_age_summary(
    connection: duckdb.DuckDBPyConnection,
) -> list[dict[str, Any]]:
    rows = connection.execute(
        """
        SELECT
            sepsis,
            count(AGE)::BIGINT AS unweighted_n,
            sum(DISCWT) FILTER (WHERE AGE IS NOT NULL) AS weighted_n,
            sum(DISCWT * AGE) / sum(DISCWT) FILTER (WHERE AGE IS NOT NULL) AS mean_age,
            sqrt(
                sum(DISCWT * pow(AGE, 2)) / sum(DISCWT) FILTER (WHERE AGE IS NOT NULL)
                - pow(sum(DISCWT * AGE) / sum(DISCWT) FILTER (WHERE AGE IS NOT NULL), 2)
            ) AS sd_age
        FROM hm_cohort
        GROUP BY sepsis
        ORDER BY sepsis
        """
    ).fetchall()
    return [
        {
            "sepsis": bool(row[0]),
            "unweighted_n": row[1],
            "weighted_n": round(row[2]),
            "mean_age": round(row[3], 2),
            "sd_age": round(row[4], 2),
        }
        for row in rows
    ]


def main() -> dict[str, Any]:
    if not COHORT_DATABASE.exists():
        raise RuntimeError("Run the Phase 1–2 notebook first to build the HM cohort.")
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    connection = duckdb.connect(str(COHORT_DATABASE), read_only=True)

    totals = {
        bool(row[0]): {"unweighted": row[1], "weighted": row[2]}
        for row in connection.execute(
            """
            SELECT sepsis, count(*)::BIGINT, sum(DISCWT)::DOUBLE
            FROM hm_cohort GROUP BY sepsis
            """
        ).fetchall()
    }
    table_rows: list[dict[str, Any]] = []
    missing_rows: list[dict[str, Any]] = []

    age = weighted_age_summary(connection)
    age_by_status = {row["sepsis"]: row for row in age}
    age_0, age_1 = age_by_status[False], age_by_status[True]
    pooled_age_sd = math.sqrt((age_0["sd_age"] ** 2 + age_1["sd_age"] ** 2) / 2)
    table_rows.append(
        {
            "characteristic": "Age, years",
            "level": "Weighted mean (SD)",
            "no_sepsis_unweighted_n": age_0["unweighted_n"],
            "no_sepsis_weighted_n": age_0["weighted_n"],
            "no_sepsis_weighted_percent_or_mean_sd": f"{age_0['mean_age']:.2f} ({age_0['sd_age']:.2f})",
            "sepsis_unweighted_n": age_1["unweighted_n"],
            "sepsis_weighted_n": age_1["weighted_n"],
            "sepsis_weighted_percent_or_mean_sd": f"{age_1['mean_age']:.2f} ({age_1['sd_age']:.2f})",
            "difference_percentage_points": "",
            "standardized_difference": round(
                (age_1["mean_age"] - age_0["mean_age"]) / pooled_age_sd, 3
            ),
        }
    )

    for characteristic, expression, level_order in CATEGORICAL_CHARACTERISTICS:
        results = connection.execute(
            f"""
            WITH categorized AS (
                SELECT sepsis, DISCWT, {expression} AS level FROM hm_cohort
            )
            SELECT
                sepsis, level, count(*)::BIGINT AS unweighted_n,
                sum(DISCWT)::DOUBLE AS weighted_n
            FROM categorized
            GROUP BY sepsis, level
            """
        ).fetchall()
        result_map = {(bool(row[0]), str(row[1])): row for row in results}
        for level in level_order:
            row_0 = result_map.get((False, level), (False, level, 0, 0.0))
            row_1 = result_map.get((True, level), (True, level, 0, 0.0))
            p0 = row_0[3] / totals[False]["weighted"]
            p1 = row_1[3] / totals[True]["weighted"]
            if row_0[2] == 0 and row_1[2] == 0:
                continue
            table_rows.append(
                {
                    "characteristic": characteristic,
                    "level": level,
                    "no_sepsis_unweighted_n": row_0[2],
                    "no_sepsis_weighted_n": round(row_0[3]),
                    "no_sepsis_weighted_percent_or_mean_sd": f"{100 * p0:.2f}%",
                    "sepsis_unweighted_n": row_1[2],
                    "sepsis_weighted_n": round(row_1[3]),
                    "sepsis_weighted_percent_or_mean_sd": f"{100 * p1:.2f}%",
                    "difference_percentage_points": round(100 * (p1 - p0), 2),
                    "standardized_difference": round(
                        standardized_difference(p1, p0), 3
                    ),
                }
            )
        missing_level = next(
            (level for level in ["Missing", "Unknown"] if level in level_order), None
        )
        if missing_level:
            missing_0 = result_map.get((False, missing_level), (False, missing_level, 0, 0.0))
            missing_1 = result_map.get((True, missing_level), (True, missing_level, 0, 0.0))
            missing_rows.append(
                {
                    "characteristic": characteristic,
                    "no_sepsis_missing_unweighted_n": missing_0[2],
                    "no_sepsis_missing_weighted_percent": round(
                        100 * missing_0[3] / totals[False]["weighted"], 2
                    ),
                    "sepsis_missing_unweighted_n": missing_1[2],
                    "sepsis_missing_weighted_percent": round(
                        100 * missing_1[3] / totals[True]["weighted"], 2
                    ),
                }
            )

    connection.close()
    write_csv(OUTPUT_DIR / "table_1_baseline_characteristics.csv", table_rows)
    write_csv(OUTPUT_DIR / "missingness_review.csv", missing_rows)
    write_csv(OUTPUT_DIR / "age_summary.csv", age)

    ranked = sorted(
        (
            row
            for row in table_rows
            if row["level"] != "Weighted mean (SD)"
            and row["level"] not in {"Missing", "Unknown"}
        ),
        key=lambda row: abs(float(row["standardized_difference"])),
        reverse=True,
    )[:8]
    summary = {
        "unit": "NIS inpatient discharge, not unique patient",
        "no_sepsis_unweighted_n": totals[False]["unweighted"],
        "no_sepsis_weighted_n": round(totals[False]["weighted"]),
        "sepsis_unweighted_n": totals[True]["unweighted"],
        "sepsis_weighted_n": round(totals[True]["weighted"]),
        "age_no_sepsis_mean_sd": f"{age_0['mean_age']:.2f} ({age_0['sd_age']:.2f})",
        "age_sepsis_mean_sd": f"{age_1['mean_age']:.2f} ({age_1['sd_age']:.2f})",
        "largest_absolute_standardized_differences": [
            {
                "characteristic": row["characteristic"],
                "level": row["level"],
                "standardized_difference": row["standardized_difference"],
            }
            for row in ranked
        ],
        "limitations": [
            "Descriptive weighted estimates; no adjusted causal interpretation.",
            "HM and A41-only sepsis phenotypes remain draft pending investigator review.",
            "Hospital characteristics are decoded from NIS_STRATUM.",
            "Comorbidity burden is deferred until a validated Charlson or Elixhauser method is selected.",
        ],
    }
    (OUTPUT_DIR / "phase_3_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    return summary


if __name__ == "__main__":
    main()
