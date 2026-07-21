"""Task 14: primary adjusted palliative-care analysis."""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path
from typing import Any

import duckdb
import numpy as np
import pandas as pd
from scipy.special import expit
from scipy.stats import norm

from src.phase_5_cci_los_mortality import CCI_COMPONENTS, CCI_WEIGHTS, component_condition
from src.phase_6_palliative_care import format_p


REPO_ROOT = Path(__file__).resolve().parents[1]
COHORT_DATABASE = REPO_ROOT / "outputs/phase_1_2/hm_cohort.duckdb"
OUTPUT_DIR = REPO_ROOT / "outputs/phase_8"

CATEGORY_SPECS = {
    "sex": ["Male", "Female", "Missing"],
    "race": ["White", "Black", "Hispanic", "Asian/Pacific Islander", "Native American", "Other", "Missing"],
    "payer": ["Medicare", "Medicaid", "Private insurance", "Self-pay", "No charge", "Other", "Missing"],
    "income": ["0–25th percentile", "26th–50th percentile", "51st–75th percentile", "76th–100th percentile", "Missing"],
    "cci_category": ["0", "1–2", "≥3"],
    "region": ["Northeast", "Midwest", "South", "West", "Unknown"],
    "location_teaching": ["Rural", "Urban nonteaching", "Urban teaching", "Unknown"],
    "bed_size": ["Small", "Medium", "Large", "Unknown"],
    "year_category": [str(year) for year in range(2016, 2023)],
    "hm_subtype_label": ["Lymphoma", "AML", "CML", "CLL/chronic leukemia", "ALL/unspecified acute leukemia", "Other leukemia", "Myeloma/plasma-cell neoplasm", "MDS", "MPN"],
}


def design_matrix(frame: pd.DataFrame, force_sepsis: int | None = None) -> tuple[np.ndarray, list[str]]:
    columns = [np.ones(len(frame)), (frame["age"].to_numpy(float) - 65.0)]
    names = ["Intercept", "Age, per year (centered at 65)"]
    sepsis = frame["sepsis"].to_numpy(float) if force_sepsis is None else np.full(len(frame), force_sepsis, dtype=float)
    columns.insert(1, sepsis)
    names.insert(1, "Documented sepsis")
    for variable, levels in CATEGORY_SPECS.items():
        values = frame[variable].astype(str).to_numpy()
        for level in levels[1:]:
            columns.append((values == level).astype(float))
            names.append(f"{variable}: {level}")
    return np.column_stack(columns), names


def fit_weighted_logistic(frame: pd.DataFrame, tolerance: float = 1e-9, max_iterations: int = 50) -> dict[str, Any]:
    x_full, all_names = design_matrix(frame)
    active_indices = np.array([0] + [index for index in range(1, x_full.shape[1]) if np.any(x_full[:, index] != x_full[0, index]) or np.any(x_full[:, index] != 0)])
    # Dummy columns that are identically zero represent categories absent from
    # the analytic cohort and must not enter the information matrix.
    active_indices = np.array([index for index in active_indices if index == 0 or np.any(x_full[:, index] != 0)])
    x = x_full[:, active_indices]
    names = [all_names[index] for index in active_indices]
    y = frame["palliative_care"].to_numpy(float)
    frequency = frame["frequency"].to_numpy(float)
    survey_weight = frame["weight"].to_numpy(float)
    combined_weight = frequency * survey_weight
    beta = np.zeros(x.shape[1])
    for iteration in range(1, max_iterations + 1):
        probability = expit(x @ beta)
        gradient = x.T @ (combined_weight * (y - probability))
        information = (x.T * (combined_weight * probability * (1 - probability))) @ x
        step = np.linalg.solve(information, gradient)
        beta += step
        if np.max(np.abs(step)) < tolerance:
            break
    else:
        raise RuntimeError("Weighted logistic regression failed to converge.")

    probability = expit(x @ beta)
    individual_score = (survey_weight * (y - probability))[:, None] * x
    strata = frame["stratum_id"].to_numpy()
    meat = np.zeros((x.shape[1], x.shape[1]))
    score_frame = pd.DataFrame({"stratum": strata, "row": np.arange(len(frame))})
    for indices in score_frame.groupby("stratum", sort=False).indices.values():
        idx = np.asarray(indices)
        scores = individual_score[idx]
        counts = frequency[idx]
        stratum_n = counts.sum()
        if stratum_n <= 1:
            continue
        score_sum = (scores * counts[:, None]).sum(axis=0)
        sum_outer = (scores.T * counts) @ scores
        meat += stratum_n / (stratum_n - 1) * (sum_outer - np.outer(score_sum, score_sum) / stratum_n)
    bread = np.linalg.inv(information)
    covariance = bread @ meat @ bread
    return {"beta": beta, "covariance": covariance, "names": names, "iterations": iteration, "active_indices": active_indices}


def standardized_probability(frame: pd.DataFrame, beta: np.ndarray, force_sepsis: int, active_indices: np.ndarray) -> tuple[float, np.ndarray]:
    x, _ = design_matrix(frame, force_sepsis=force_sepsis)
    x = x[:, active_indices]
    weights = frame["frequency"].to_numpy(float) * frame["weight"].to_numpy(float)
    probability = expit(x @ beta)
    total_weight = weights.sum()
    estimate = float(np.sum(weights * probability) / total_weight)
    gradient = np.sum((weights * probability * (1 - probability))[:, None] * x, axis=0) / total_weight
    return estimate, gradient


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def build_model_frame(connection: duckdb.DuckDBPyConnection) -> pd.DataFrame:
    flag_sql = []
    for component, prefixes in CCI_COMPONENTS.items():
        condition = component_condition(prefixes)
        flag_sql.append(f"list_contains(list_transform(diagnosis_codes, code -> {condition}), TRUE) AS {component}")
    score_terms = []
    for component, weight in CCI_WEIGHTS.items():
        present = "diab AND NOT diabwc" if component == "diab" else "mld AND NOT msld" if component == "mld" else component
        score_terms.append(f"CASE WHEN {present} THEN {weight} ELSE 0 END")
    connection.execute("CREATE TEMP TABLE model_flags AS SELECT *, " + ", ".join(flag_sql) + " FROM hm_cohort")
    connection.execute("CREATE TEMP TABLE model_scored AS SELECT *, " + " + ".join(score_terms) + " AS cci FROM model_flags")
    return connection.execute("""
        SELECT
            concat(CAST(YEAR AS INTEGER), ':', CAST(NIS_STRATUM AS INTEGER)) AS stratum_id,
            DISCWT::DOUBLE AS weight, sepsis::INTEGER AS sepsis,
            palliative_care::INTEGER AS palliative_care, AGE::DOUBLE AS age,
            CASE FEMALE WHEN 0 THEN 'Male' WHEN 1 THEN 'Female' ELSE 'Missing' END AS sex,
            CASE RACE WHEN 1 THEN 'White' WHEN 2 THEN 'Black' WHEN 3 THEN 'Hispanic'
                WHEN 4 THEN 'Asian/Pacific Islander' WHEN 5 THEN 'Native American'
                WHEN 6 THEN 'Other' ELSE 'Missing' END AS race,
            CASE PAY1 WHEN 1 THEN 'Medicare' WHEN 2 THEN 'Medicaid' WHEN 3 THEN 'Private insurance'
                WHEN 4 THEN 'Self-pay' WHEN 5 THEN 'No charge' WHEN 6 THEN 'Other' ELSE 'Missing' END AS payer,
            CASE ZIPINC_QRTL WHEN 1 THEN '0–25th percentile' WHEN 2 THEN '26th–50th percentile'
                WHEN 3 THEN '51st–75th percentile' WHEN 4 THEN '76th–100th percentile' ELSE 'Missing' END AS income,
            CASE WHEN cci = 0 THEN '0' WHEN cci BETWEEN 1 AND 2 THEN '1–2' ELSE '≥3' END AS cci_category,
            coalesce(hospital_region, 'Unknown') AS region,
            coalesce(hospital_location_teaching, 'Unknown') AS location_teaching,
            coalesce(hospital_bed_size, 'Unknown') AS bed_size,
            CAST(CAST(YEAR AS INTEGER) AS VARCHAR) AS year_category,
            CASE hm_subtype WHEN 'lymphoma' THEN 'Lymphoma' WHEN 'aml' THEN 'AML' WHEN 'cml' THEN 'CML'
                WHEN 'cll_chronic_leukemia' THEN 'CLL/chronic leukemia'
                WHEN 'all' THEN 'ALL/unspecified acute leukemia' WHEN 'other_leukemia' THEN 'Other leukemia'
                WHEN 'myeloma_plasma_cell' THEN 'Myeloma/plasma-cell neoplasm'
                WHEN 'mds' THEN 'MDS' WHEN 'mpn' THEN 'MPN' END AS hm_subtype_label,
            count(*)::BIGINT AS frequency
        FROM model_scored
        WHERE DISCWT IS NOT NULL AND NIS_STRATUM IS NOT NULL AND AGE IS NOT NULL
        GROUP BY ALL
    """).fetchdf()


def main() -> dict[str, Any]:
    if not COHORT_DATABASE.exists():
        raise RuntimeError("Run Phase 1–2 first to create the cached HM cohort.")
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    connection = duckdb.connect(str(COHORT_DATABASE), read_only=True)
    total_records = connection.execute("SELECT count(*) FROM hm_cohort").fetchone()[0]
    frame = build_model_frame(connection)
    connection.close()
    included_records = int(frame["frequency"].sum())
    model = fit_weighted_logistic(frame)
    beta = model["beta"]
    covariance = model["covariance"]
    sepsis_index = model["names"].index("Documented sepsis")
    coefficient = beta[sepsis_index]
    standard_error = math.sqrt(covariance[sepsis_index, sepsis_index])
    odds_ratio = math.exp(coefficient)
    odds_ratio_lower = math.exp(coefficient - 1.96 * standard_error)
    odds_ratio_upper = math.exp(coefficient + 1.96 * standard_error)
    odds_ratio_p = float(2 * norm.sf(abs(coefficient / standard_error)))

    probability_0, gradient_0 = standardized_probability(frame, beta, 0, model["active_indices"])
    probability_1, gradient_1 = standardized_probability(frame, beta, 1, model["active_indices"])
    difference = probability_1 - probability_0
    gradient_difference = gradient_1 - gradient_0
    difference_se = math.sqrt(float(gradient_difference @ covariance @ gradient_difference))
    difference_p = float(2 * norm.sf(abs(difference / difference_se)))
    probability_rows = []
    for label, estimate, gradient in [
        ("No documented sepsis", probability_0, gradient_0),
        ("Documented sepsis", probability_1, gradient_1),
    ]:
        se = math.sqrt(float(gradient @ covariance @ gradient))
        probability_rows.append({
            "sepsis_status": label,
            "adjusted_probability_percent": round(100 * estimate, 2),
            "ci_95_lower_percent": round(100 * max(0, estimate - 1.96 * se), 2),
            "ci_95_upper_percent": round(100 * min(1, estimate + 1.96 * se), 2),
        })
    probability_rows.append({
        "sepsis_status": "Total observed cohort",
        "adjusted_probability_percent": "—", "ci_95_lower_percent": "—", "ci_95_upper_percent": "—",
    })
    primary_rows = [
        {
            "measure": "Adjusted odds ratio for documented sepsis",
            "estimate": round(odds_ratio, 3),
            "ci_95": f"{odds_ratio_lower:.3f}–{odds_ratio_upper:.3f}",
            "p_value": format_p(odds_ratio_p),
        },
        {
            "measure": "Adjusted absolute probability difference, percentage points",
            "estimate": round(100 * difference, 2),
            "ci_95": f"{100*(difference-1.96*difference_se):.2f}–{100*(difference+1.96*difference_se):.2f}",
            "p_value": format_p(difference_p),
        },
        {"measure": "Total", "estimate": "—", "ci_95": "—", "p_value": "—"},
    ]
    write_csv(OUTPUT_DIR / "primary_adjusted_results.csv", primary_rows)
    write_csv(OUTPUT_DIR / "adjusted_probabilities.csv", probability_rows)
    coefficient_rows = []
    for index, name in enumerate(model["names"]):
        se = math.sqrt(covariance[index, index])
        coefficient_rows.append({
            "term": name, "adjusted_odds_ratio": round(math.exp(beta[index]), 3),
            "ci_95_lower": round(math.exp(beta[index] - 1.96 * se), 3),
            "ci_95_upper": round(math.exp(beta[index] + 1.96 * se), 3),
            "p_value": format_p(float(2 * norm.sf(abs(beta[index] / se)))),
        })
    coefficient_rows.append({"term": "Total", "adjusted_odds_ratio": "—", "ci_95_lower": "—", "ci_95_upper": "—", "p_value": "—"})
    write_csv(OUTPUT_DIR / "full_model_coefficients.csv", coefficient_rows)
    summary = {
        "outcome": "Documented inpatient palliative-care use (normalized Z515 in any diagnosis position)",
        "included_unweighted_records": included_records,
        "excluded_records": int(total_records - included_records),
        "iterations": model["iterations"],
        "reference_categories": {key: levels[0] for key, levels in CATEGORY_SPECS.items()},
        "variance_note": "DISCWT-weighted logistic model with year-specific NIS_STRATUM linearization and discharge-level variance units; not full NIS hospital-cluster-adjusted inference.",
        "primary_results": primary_rows,
        "adjusted_probabilities": probability_rows,
    }
    (OUTPUT_DIR / "phase_8_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    return summary


if __name__ == "__main__":
    main()
