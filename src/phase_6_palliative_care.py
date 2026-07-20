"""Task 12A–B: documented inpatient palliative-care use by sepsis status."""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path
from typing import Any

import duckdb
import numpy as np
import pandas as pd
from scipy.stats import norm


REPO_ROOT = Path(__file__).resolve().parents[1]
COHORT_DATABASE = REPO_ROOT / "outputs/phase_1_2/hm_cohort.duckdb"
OUTPUT_DIR = REPO_ROOT / "outputs/phase_6"


def stratified_variance(influence: np.ndarray, strata: np.ndarray) -> float:
    """Taylor variance across strata, treating discharges as variance units."""
    frame = pd.DataFrame({"influence": influence, "stratum": strata})
    result = 0.0
    for _, group in frame.groupby("stratum", sort=False):
        count = len(group)
        if count > 1:
            values = group["influence"].to_numpy(dtype=float)
            result += count / (count - 1) * np.square(values - values.mean()).sum()
    return float(result)


def weighted_prevalence(
    weights: np.ndarray,
    outcome: np.ndarray,
    strata: np.ndarray,
    domain: np.ndarray | None = None,
) -> dict[str, Any]:
    if domain is None:
        domain = np.ones(len(outcome), dtype=bool)
    domain_float = domain.astype(float)
    denominator = float(np.sum(weights * domain_float))
    numerator = float(np.sum(weights * domain_float * outcome))
    estimate = numerator / denominator
    influence = weights * domain_float * (outcome - estimate) / denominator
    standard_error = math.sqrt(stratified_variance(influence, strata))
    # Logit CI stays within [0, 1].
    if 0 < estimate < 1:
        logit = math.log(estimate / (1 - estimate))
        logit_se = standard_error / (estimate * (1 - estimate))
        lower = 1 / (1 + math.exp(-(logit - 1.96 * logit_se)))
        upper = 1 / (1 + math.exp(-(logit + 1.96 * logit_se)))
    else:
        lower = upper = estimate
    return {
        "unweighted_n": int(domain.sum()),
        "unweighted_events": int(np.sum(domain_float * outcome)),
        "weighted_n": denominator,
        "weighted_events": numerator,
        "estimate": estimate,
        "standard_error": standard_error,
        "ci_lower": lower,
        "ci_upper": upper,
        "influence": influence,
    }


def comparison(
    weights: np.ndarray,
    outcome: np.ndarray,
    strata: np.ndarray,
    exposure: np.ndarray,
) -> dict[str, float]:
    exposed = weighted_prevalence(weights, outcome, strata, exposure)
    unexposed = weighted_prevalence(weights, outcome, strata, ~exposure)
    difference = exposed["estimate"] - unexposed["estimate"]
    difference_influence = exposed["influence"] - unexposed["influence"]
    difference_se = math.sqrt(stratified_variance(difference_influence, strata))
    difference_z = difference / difference_se
    difference_p = 2 * norm.sf(abs(difference_z))

    odds_ratio = (
        exposed["estimate"] / (1 - exposed["estimate"])
        / (unexposed["estimate"] / (1 - unexposed["estimate"]))
    )
    log_or_influence = (
        exposed["influence"] / (exposed["estimate"] * (1 - exposed["estimate"]))
        - unexposed["influence"] / (unexposed["estimate"] * (1 - unexposed["estimate"]))
    )
    log_or_se = math.sqrt(stratified_variance(log_or_influence, strata))
    log_or = math.log(odds_ratio)
    return {
        "difference": difference,
        "difference_ci_lower": difference - 1.96 * difference_se,
        "difference_ci_upper": difference + 1.96 * difference_se,
        "difference_p": float(difference_p),
        "odds_ratio": odds_ratio,
        "odds_ratio_ci_lower": math.exp(log_or - 1.96 * log_or_se),
        "odds_ratio_ci_upper": math.exp(log_or + 1.96 * log_or_se),
        "odds_ratio_p": float(2 * norm.sf(abs(log_or / log_or_se))),
    }


def format_p(value: float) -> str:
    return "<0.001" if value < 0.001 else f"{value:.3f}"


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
    frame = connection.execute("""
        SELECT YEAR::INTEGER AS year, NIS_STRATUM::INTEGER AS nis_stratum,
               DISCWT::DOUBLE AS weight, sepsis::BOOLEAN AS sepsis,
               palliative_care::INTEGER AS palliative_care
        FROM hm_cohort
        WHERE DISCWT IS NOT NULL AND NIS_STRATUM IS NOT NULL
    """).fetchdf()
    total_records = connection.execute("SELECT count(*) FROM hm_cohort").fetchone()[0]
    normalized_check = connection.execute("""
        SELECT count(*) FILTER (WHERE palliative_care),
               count(*) FILTER (WHERE list_contains(diagnosis_codes, 'Z515')),
               count(*) FILTER (WHERE palliative_care != list_contains(diagnosis_codes, 'Z515'))
        FROM hm_cohort
    """).fetchone()
    connection.close()
    if normalized_check[2] != 0:
        raise RuntimeError("Cached palliative-care flag does not match normalized Z515 search.")

    weights_array = frame["weight"].to_numpy(dtype=float)
    outcome = frame["palliative_care"].to_numpy(dtype=float)
    exposure = frame["sepsis"].to_numpy(dtype=bool)
    # NIS strata are year-specific in the pooled 2016–2022 analysis.
    strata = (frame["year"].astype(str) + ":" + frame["nis_stratum"].astype(str)).to_numpy()

    overall = weighted_prevalence(weights_array, outcome, strata)
    no_sepsis = weighted_prevalence(weights_array, outcome, strata, ~exposure)
    sepsis = weighted_prevalence(weights_array, outcome, strata, exposure)
    contrast = comparison(weights_array, outcome, strata, exposure)

    rows = []
    for label, result in [
        ("All adult HM hospitalizations", overall),
        ("HM without documented sepsis", no_sepsis),
        ("HM with documented sepsis", sepsis),
    ]:
        rows.append({
            "cohort": label,
            "unweighted_sample_n": result["unweighted_n"],
            "unweighted_palliative_care_n": result["unweighted_events"],
            "weighted_hospitalizations": round(result["weighted_n"]),
            "weighted_palliative_care_n": round(result["weighted_events"]),
            "weighted_palliative_care_percent": round(100 * result["estimate"], 2),
            "ci_95_lower_percent": round(100 * result["ci_lower"], 2),
            "ci_95_upper_percent": round(100 * result["ci_upper"], 2),
        })
    write_csv(OUTPUT_DIR / "palliative_care_prevalence.csv", rows)

    annual_rows: list[dict[str, Any]] = []
    years = frame["year"].to_numpy(dtype=int)
    for year in sorted(frame["year"].unique()):
        for status, label in [
            (False, "No documented sepsis"),
            (True, "Documented sepsis"),
        ]:
            domain = (years == year) & (exposure == status)
            result = weighted_prevalence(weights_array, outcome, strata, domain)
            annual_rows.append({
                "year": int(year),
                "sepsis_status": label,
                "unweighted_sample_n": result["unweighted_n"],
                "unweighted_palliative_care_n": result["unweighted_events"],
                "weighted_hospitalizations": round(result["weighted_n"]),
                "weighted_palliative_care_n": round(result["weighted_events"]),
                "weighted_palliative_care_percent": round(100 * result["estimate"], 2),
                "ci_95_lower_percent": round(100 * result["ci_lower"], 2),
                "ci_95_upper_percent": round(100 * result["ci_upper"], 2),
            })
    write_csv(OUTPUT_DIR / "annual_palliative_care_by_sepsis.csv", annual_rows)
    summary = {
        "definition": "Documented inpatient palliative-care use: normalized Z51.5 (Z515) in any diagnosis position.",
        "normalization": "Uppercase; decimal points and spaces removed in the Phase 1–2 cohort build.",
        "records_excluded_for_missing_weight_or_stratum": int(total_records - len(frame)),
        "normalization_validation": {
            "flagged_records": int(normalized_check[0]),
            "records_containing_normalized_Z515": int(normalized_check[1]),
            "discordant_records": int(normalized_check[2]),
        },
        "variance_note": "NIS_STRATUM and year-specific strata used in Taylor linearization; sampled discharges are variance units because HOSP_NIS is unavailable by study decision. Not a full NIS design variance estimate.",
        "prevalence_table": rows,
        "annual_prevalence_table": annual_rows,
        "comparison": {
            "absolute_difference_percentage_points": round(100 * contrast["difference"], 2),
            "difference_ci_95_lower": round(100 * contrast["difference_ci_lower"], 2),
            "difference_ci_95_upper": round(100 * contrast["difference_ci_upper"], 2),
            "strata_adjusted_p_value": format_p(contrast["difference_p"]),
            "crude_weighted_odds_ratio": round(contrast["odds_ratio"], 3),
            "odds_ratio_ci_95_lower": round(contrast["odds_ratio_ci_lower"], 3),
            "odds_ratio_ci_95_upper": round(contrast["odds_ratio_ci_upper"], 3),
            "odds_ratio_p_value": format_p(contrast["odds_ratio_p"]),
        },
        "adjusted_probabilities_note": "Covariate-adjusted probabilities require the later multivariable model (Command 16B) and are not estimated in this unadjusted Task 12B analysis.",
    }
    (OUTPUT_DIR / "phase_6_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    return summary


if __name__ == "__main__":
    main()
