"""Command 15: palliative-care utilization by mutually exclusive HM subtype."""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any

import duckdb
import numpy as np
import pandas as pd
from scipy.stats import chi2

from src.phase_6_palliative_care import format_p, weighted_prevalence


REPO_ROOT = Path(__file__).resolve().parents[1]
COHORT_DATABASE = REPO_ROOT / "outputs/phase_1_2/hm_cohort.duckdb"
OUTPUT_DIR = REPO_ROOT / "outputs/phase_7"
SUBTYPES = [
    ("lymphoma", "Lymphoma"),
    ("aml", "AML"),
    ("cml", "CML"),
    ("cll_chronic_leukemia", "CLL/chronic leukemia"),
    ("all", "ALL/unspecified acute leukemia"),
    ("other_leukemia", "Other leukemia"),
    ("myeloma_plasma_cell", "Myeloma/plasma-cell neoplasm"),
    ("mds", "MDS"),
    ("mpn", "MPN"),
]


def stratified_covariance(influence: np.ndarray, strata: np.ndarray) -> np.ndarray:
    """Taylor covariance across strata, treating discharges as variance units."""
    covariance = np.zeros((influence.shape[1], influence.shape[1]), dtype=float)
    frame = pd.DataFrame({"stratum": strata})
    for indices in frame.groupby("stratum", sort=False).indices.values():
        values = influence[np.asarray(indices)]
        count = len(values)
        if count > 1:
            centered = values - values.mean(axis=0)
            covariance += count / (count - 1) * centered.T @ centered
    return covariance


def overall_subtype_test(estimates: np.ndarray, influence: np.ndarray, strata: np.ndarray) -> dict[str, Any]:
    """Wald test that all subtype-specific weighted prevalences are equal."""
    categories = len(estimates)
    contrast = np.zeros((categories - 1, categories))
    for row in range(categories - 1):
        contrast[row, 0] = -1
        contrast[row, row + 1] = 1
    differences = contrast @ estimates
    covariance = contrast @ stratified_covariance(influence, strata) @ contrast.T
    rank = int(np.linalg.matrix_rank(covariance))
    statistic = float(differences.T @ np.linalg.pinv(covariance) @ differences)
    return {"wald_chi_square": statistic, "degrees_of_freedom": rank, "p_value": float(chi2.sf(statistic, rank))}


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def cohort_table(frame: pd.DataFrame, cohort_label: str) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    weights = frame["weight"].to_numpy(dtype=float)
    outcome = frame["palliative_care"].to_numpy(dtype=float)
    subtype = frame["hm_subtype"].to_numpy()
    strata = (frame["year"].astype(str) + ":" + frame["nis_stratum"].astype(str)).to_numpy()
    results = []
    influences = []
    for subtype_id, subtype_label in SUBTYPES:
        result = weighted_prevalence(weights, outcome, strata, subtype == subtype_id)
        results.append((subtype_label, result))
        influences.append(result["influence"])
    test = overall_subtype_test(
        np.array([result["estimate"] for _, result in results]),
        np.column_stack(influences), strata,
    )
    rows = []
    for index, (subtype_label, result) in enumerate(results):
        rows.append({
            "cohort": cohort_label,
            "hm_subtype": subtype_label,
            "unweighted_sample_n": result["unweighted_n"],
            "unweighted_palliative_care_n": result["unweighted_events"],
            "weighted_hospitalizations": round(result["weighted_n"]),
            "weighted_palliative_care_n": round(result["weighted_events"]),
            "weighted_palliative_care_percent": round(100 * result["estimate"], 2),
            "ci_95_lower_percent": round(100 * result["ci_lower"], 2),
            "ci_95_upper_percent": round(100 * result["ci_upper"], 2),
            "overall_p_value": format_p(test["p_value"]) if index == 0 else "",
        })
    return rows, {
        "cohort": cohort_label,
        "wald_chi_square": round(test["wald_chi_square"], 3),
        "degrees_of_freedom": test["degrees_of_freedom"],
        "overall_p_value": format_p(test["p_value"]),
    }


def main() -> dict[str, Any]:
    if not COHORT_DATABASE.exists():
        raise RuntimeError("Run Phase 1–2 first to create the cached HM cohort.")
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    connection = duckdb.connect(str(COHORT_DATABASE), read_only=True)
    frame = connection.execute("""
        SELECT YEAR::INTEGER AS year, NIS_STRATUM::INTEGER AS nis_stratum,
               DISCWT::DOUBLE AS weight, sepsis::BOOLEAN AS sepsis,
               palliative_care::INTEGER AS palliative_care, hm_subtype
        FROM hm_cohort
        WHERE DISCWT IS NOT NULL AND NIS_STRATUM IS NOT NULL AND hm_subtype IS NOT NULL
    """).fetchdf()
    total = connection.execute("SELECT count(*) FROM hm_cohort").fetchone()[0]
    connection.close()

    definitions = [
        ("All adult HM hospitalizations", frame),
        ("HM without documented sepsis", frame.loc[~frame["sepsis"]].reset_index(drop=True)),
        ("HM with documented sepsis", frame.loc[frame["sepsis"]].reset_index(drop=True)),
    ]
    all_rows: list[dict[str, Any]] = []
    tests: list[dict[str, Any]] = []
    for filename, (label, cohort) in zip(
        ["all_hm", "no_sepsis", "sepsis"], definitions
    ):
        rows, test = cohort_table(cohort, label)
        write_csv(OUTPUT_DIR / f"subtype_palliative_care_{filename}.csv", rows)
        all_rows.extend(rows)
        tests.append(test)
    write_csv(OUTPUT_DIR / "subtype_palliative_care_all_tables.csv", all_rows)
    summary = {
        "definition": "Documented inpatient palliative-care use: normalized Z51.5 (Z515) in any diagnosis position.",
        "subtype_assignment": "Mutually exclusive first-listed qualifying HM subtype; overlaps retained in multiple_hm_groups for sensitivity analysis.",
        "records_excluded_for_missing_weight_stratum_or_subtype": int(total - len(frame)),
        "variance_note": "Year-specific NIS_STRATUM used in Taylor linearization; sampled discharges are variance units because HOSP_NIS is unavailable by study decision. Not a full NIS design variance estimate.",
        "overall_tests": tests,
        "tables": all_rows,
    }
    (OUTPUT_DIR / "phase_7_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    return summary


if __name__ == "__main__":
    main()
