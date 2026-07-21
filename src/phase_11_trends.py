"""Commands 20A–20C: annual trends in documented palliative-care use."""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path
from typing import Any, Callable

import duckdb
import numpy as np
import pandas as pd
from scipy.stats import chi2, norm

from src.phase_6_palliative_care import comparison, format_p, weighted_prevalence
from src.phase_8_primary_adjusted import COHORT_DATABASE, build_model_frame
from src.phase_9_interaction import SUBTYPE_LEVELS


REPO_ROOT = Path(__file__).resolve().parents[1]
OUTPUT_DIR = REPO_ROOT / "outputs/phase_11"
SMALL_CELL_EVENT_THRESHOLD = 10


def fit_modified_poisson(frame: pd.DataFrame, builder: Callable[[pd.DataFrame], tuple[np.ndarray, list[str]]]) -> dict[str, Any]:
    x, names = builder(frame)
    y = frame["palliative_care"].to_numpy(float)
    frequency = frame["frequency"].to_numpy(float)
    survey_weight = frame["weight"].to_numpy(float)
    combined_weight = frequency * survey_weight
    beta = np.zeros(x.shape[1])
    beta[0] = math.log(np.average(y, weights=combined_weight))
    for iteration in range(1, 51):
        mu = np.exp(np.clip(x @ beta, -30, 20))
        gradient = x.T @ (combined_weight * (y - mu))
        information = (x.T * (combined_weight * mu)) @ x
        step = np.linalg.solve(information, gradient)
        beta += step
        if np.max(np.abs(step)) < 1e-10:
            break
    else:
        raise RuntimeError("Modified-Poisson trend model failed to converge.")
    mu = np.exp(np.clip(x @ beta, -30, 20))
    individual_score = (survey_weight * (y - mu))[:, None] * x
    meat = np.zeros((x.shape[1], x.shape[1]))
    grouped = pd.DataFrame({"stratum": frame["stratum_id"], "row": np.arange(len(frame))}).groupby("stratum", sort=False).indices
    for indices in grouped.values():
        idx = np.asarray(indices)
        scores, counts = individual_score[idx], frequency[idx]
        stratum_n = counts.sum()
        if stratum_n <= 1:
            continue
        score_sum = (scores * counts[:, None]).sum(axis=0)
        meat += stratum_n / (stratum_n - 1) * ((scores.T * counts) @ scores - np.outer(score_sum, score_sum) / stratum_n)
    bread = np.linalg.inv(information)
    return {"beta": beta, "covariance": bread @ meat @ bread, "names": names, "iterations": iteration}


def overall_design(frame: pd.DataFrame) -> tuple[np.ndarray, list[str]]:
    year = frame["year_category"].astype(int).to_numpy(float) - 2016
    return np.column_stack([np.ones(len(frame)), year]), ["Intercept", "Year"]


def sepsis_design(frame: pd.DataFrame) -> tuple[np.ndarray, list[str]]:
    year = frame["year_category"].astype(int).to_numpy(float) - 2016
    sepsis = frame["sepsis"].to_numpy(float)
    return np.column_stack([np.ones(len(frame)), year, sepsis, year * sepsis]), ["Intercept", "Year", "Sepsis", "Year × sepsis"]


def subtype_design(frame: pd.DataFrame) -> tuple[np.ndarray, list[str]]:
    year = frame["year_category"].astype(int).to_numpy(float) - 2016
    subtype = frame["hm_subtype_label"].astype(str).to_numpy()
    columns = [np.ones(len(frame)), year]
    names = ["Intercept", "Year"]
    dummy_columns = []
    for level in SUBTYPE_LEVELS[1:]:
        dummy = (subtype == level).astype(float)
        columns.append(dummy); dummy_columns.append(dummy); names.append(f"Subtype: {level}")
    for level, dummy in zip(SUBTYPE_LEVELS[1:], dummy_columns):
        columns.append(year * dummy); names.append(f"Year × {level}")
    return np.column_stack(columns), names


def eapc(beta: float, variance: float) -> dict[str, Any]:
    standard_error = math.sqrt(variance)
    return {
        "eapc_percent": round(100 * (math.exp(beta) - 1), 2),
        "ci_95_lower_percent": round(100 * (math.exp(beta - 1.96 * standard_error) - 1), 2),
        "ci_95_upper_percent": round(100 * (math.exp(beta + 1.96 * standard_error) - 1), 2),
        "p_value": format_p(float(2 * norm.sf(abs(beta / standard_error)))),
    }


def cell_prevalence(frame: pd.DataFrame) -> dict[str, Any]:
    expanded = frame.loc[frame.index.repeat(frame["frequency"])].reset_index(drop=True)
    weights = expanded["weight"].to_numpy(float)
    outcome = expanded["palliative_care"].to_numpy(float)
    strata = expanded["stratum_id"].to_numpy()
    return weighted_prevalence(weights, outcome, strata)


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)


def main() -> dict[str, Any]:
    if not COHORT_DATABASE.exists():
        raise RuntimeError("Run Phase 1–2 first to create the cached HM cohort.")
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    connection = duckdb.connect(str(COHORT_DATABASE), read_only=True)
    frame = build_model_frame(connection)
    connection.close()
    years = list(range(2016, 2023))

    annual_overall = []
    for year in years:
        result = cell_prevalence(frame.loc[frame["year_category"].astype(int) == year].reset_index(drop=True))
        annual_overall.append({
            "year": year, "unweighted_sample_n": result["unweighted_n"],
            "unweighted_palliative_care_n": result["unweighted_events"],
            "weighted_hospitalizations": round(result["weighted_n"]),
            "weighted_palliative_care_n": round(result["weighted_events"]),
            "weighted_palliative_care_percent": round(100 * result["estimate"], 2),
            "ci_95_lower_percent": round(100 * result["ci_lower"], 2),
            "ci_95_upper_percent": round(100 * result["ci_upper"], 2),
        })
    overall_model = fit_modified_poisson(frame, overall_design)
    overall_eapc = eapc(overall_model["beta"][1], overall_model["covariance"][1, 1])
    overall_total = cell_prevalence(frame)
    annual_overall.append({"year": "Total, 2016–2022", "unweighted_sample_n": overall_total["unweighted_n"], "unweighted_palliative_care_n": overall_total["unweighted_events"], "weighted_hospitalizations": round(overall_total["weighted_n"]), "weighted_palliative_care_n": round(overall_total["weighted_events"]), "weighted_palliative_care_percent": round(100 * overall_total["estimate"], 2), "ci_95_lower_percent": round(100 * overall_total["ci_lower"], 2), "ci_95_upper_percent": round(100 * overall_total["ci_upper"], 2)})
    write_csv(OUTPUT_DIR / "annual_overall_palliative_care.csv", annual_overall)

    annual_sepsis = []
    for year in years:
        year_frame = frame.loc[frame["year_category"].astype(int) == year].reset_index(drop=True)
        expanded_year = year_frame.loc[year_frame.index.repeat(year_frame["frequency"])].reset_index(drop=True)
        year_contrast = comparison(expanded_year["weight"].to_numpy(float), expanded_year["palliative_care"].to_numpy(float), expanded_year["stratum_id"].to_numpy(), expanded_year["sepsis"].to_numpy(bool))
        for status, label in [(0, "No documented sepsis"), (1, "Documented sepsis")]:
            result = weighted_prevalence(
                expanded_year["weight"].to_numpy(float),
                expanded_year["palliative_care"].to_numpy(float),
                expanded_year["stratum_id"].to_numpy(),
                expanded_year["sepsis"].to_numpy(int) == status,
            )
            annual_sepsis.append({"year": year, "sepsis_status": label, "unweighted_sample_n": result["unweighted_n"], "unweighted_palliative_care_n": result["unweighted_events"], "weighted_hospitalizations": round(result["weighted_n"]), "weighted_palliative_care_n": round(result["weighted_events"]), "weighted_palliative_care_percent": round(100 * result["estimate"], 2), "ci_95_lower_percent": round(100 * result["ci_lower"], 2), "ci_95_upper_percent": round(100 * result["ci_upper"], 2), "within_year_p_value": format_p(year_contrast["difference_p"])})
    sepsis_model = fit_modified_poisson(frame, sepsis_design)
    no_sepsis_eapc = eapc(sepsis_model["beta"][1], sepsis_model["covariance"][1, 1])
    sepsis_gradient = np.array([0, 1, 0, 1], dtype=float)
    sepsis_year_beta = float(sepsis_gradient @ sepsis_model["beta"])
    sepsis_year_variance = float(sepsis_gradient @ sepsis_model["covariance"] @ sepsis_gradient)
    sepsis_eapc = eapc(sepsis_year_beta, sepsis_year_variance)
    interaction_index = sepsis_model["names"].index("Year × sepsis")
    interaction_beta = sepsis_model["beta"][interaction_index]
    interaction_se = math.sqrt(sepsis_model["covariance"][interaction_index, interaction_index])
    trend_interaction_p = float(2 * norm.sf(abs(interaction_beta / interaction_se)))
    expanded_all = frame.loc[frame.index.repeat(frame["frequency"])].reset_index(drop=True)
    overall_sepsis_contrast = comparison(expanded_all["weight"].to_numpy(float), expanded_all["palliative_care"].to_numpy(float), expanded_all["stratum_id"].to_numpy(), expanded_all["sepsis"].to_numpy(bool))
    for status, label in [(0, "No documented sepsis"), (1, "Documented sepsis")]:
        total_result = weighted_prevalence(
            expanded_all["weight"].to_numpy(float), expanded_all["palliative_care"].to_numpy(float),
            expanded_all["stratum_id"].to_numpy(), expanded_all["sepsis"].to_numpy(int) == status,
        )
        annual_sepsis.append({"year": "Total, 2016–2022", "sepsis_status": label, "unweighted_sample_n": total_result["unweighted_n"], "unweighted_palliative_care_n": total_result["unweighted_events"], "weighted_hospitalizations": round(total_result["weighted_n"]), "weighted_palliative_care_n": round(total_result["weighted_events"]), "weighted_palliative_care_percent": round(100 * total_result["estimate"], 2), "ci_95_lower_percent": round(100 * total_result["ci_lower"], 2), "ci_95_upper_percent": round(100 * total_result["ci_upper"], 2), "within_year_p_value": format_p(overall_sepsis_contrast["difference_p"])})
    write_csv(OUTPUT_DIR / "annual_palliative_care_by_sepsis.csv", annual_sepsis)

    sepsis_only = frame.loc[frame["sepsis"] == 1].reset_index(drop=True)
    subtype_annual = []
    suppressed_cells = 0
    expanded_by_year = {
        year: frame.loc[frame["year_category"].astype(int) == year].reset_index(drop=True)
        for year in years
    }
    expanded_by_year = {
        year: year_frame.loc[year_frame.index.repeat(year_frame["frequency"])].reset_index(drop=True)
        for year, year_frame in expanded_by_year.items()
    }
    for subtype in SUBTYPE_LEVELS:
        for year in years:
            universe = expanded_by_year[year]
            domain = (universe["sepsis"].to_numpy(int) == 1) & (universe["hm_subtype_label"].astype(str).to_numpy() == subtype)
            result = weighted_prevalence(
                universe["weight"].to_numpy(float), universe["palliative_care"].to_numpy(float),
                universe["stratum_id"].to_numpy(), domain,
            )
            suppress = result["unweighted_events"] <= SMALL_CELL_EVENT_THRESHOLD
            suppressed_cells += int(suppress)
            subtype_annual.append({
                "hm_subtype": subtype, "year": year, "unweighted_sample_n": result["unweighted_n"],
                "unweighted_palliative_care_n": result["unweighted_events"],
                "weighted_hospitalizations": round(result["weighted_n"]),
                "weighted_palliative_care_n": "Suppressed" if suppress else round(result["weighted_events"]),
                "weighted_palliative_care_percent": "Suppressed" if suppress else round(100 * result["estimate"], 2),
                "ci_95_lower_percent": "Suppressed" if suppress else round(100 * result["ci_lower"], 2),
                "ci_95_upper_percent": "Suppressed" if suppress else round(100 * result["ci_upper"], 2),
            })
    subtype_model = fit_modified_poisson(sepsis_only, subtype_design)
    interaction_indices = [i for i, name in enumerate(subtype_model["names"]) if name.startswith("Year ×")]
    interaction_beta_vector = subtype_model["beta"][interaction_indices]
    interaction_cov = subtype_model["covariance"][np.ix_(interaction_indices, interaction_indices)]
    interaction_df = int(np.linalg.matrix_rank(interaction_cov))
    interaction_stat = float(interaction_beta_vector.T @ np.linalg.pinv(interaction_cov) @ interaction_beta_vector)
    subtype_interaction_p = float(chi2.sf(interaction_stat, interaction_df))
    subtype_trends = []
    for subtype in SUBTYPE_LEVELS:
        gradient = np.zeros(len(subtype_model["beta"])); gradient[subtype_model["names"].index("Year")] = 1
        if subtype != SUBTYPE_LEVELS[0]: gradient[subtype_model["names"].index(f"Year × {subtype}")] = 1
        estimate = float(gradient @ subtype_model["beta"]); variance = float(gradient @ subtype_model["covariance"] @ gradient)
        row = {"hm_subtype": subtype, **eapc(estimate, variance), "overall_year_by_subtype_p_value": format_p(subtype_interaction_p) if subtype == SUBTYPE_LEVELS[0] else ""}
        subtype_trends.append(row)
    subtype_trends.append({"hm_subtype": "Total", "eapc_percent": "—", "ci_95_lower_percent": "—", "ci_95_upper_percent": "—", "p_value": "—", "overall_year_by_subtype_p_value": format_p(subtype_interaction_p)})
    trend_p_by_subtype = {row["hm_subtype"]: row["p_value"] for row in subtype_trends if row["hm_subtype"] != "Total"}
    for row in subtype_annual:
        row["subtype_trend_p_value"] = trend_p_by_subtype[row["hm_subtype"]]
        row["overall_year_by_subtype_p_value"] = format_p(subtype_interaction_p) if row["hm_subtype"] == SUBTYPE_LEVELS[0] and row["year"] == years[0] else ""
    sepsis_total = weighted_prevalence(
        expanded_all["weight"].to_numpy(float), expanded_all["palliative_care"].to_numpy(float),
        expanded_all["stratum_id"].to_numpy(), expanded_all["sepsis"].to_numpy(int) == 1,
    )
    subtype_annual.append({"hm_subtype": "Total", "year": "2016–2022", "unweighted_sample_n": sepsis_total["unweighted_n"], "unweighted_palliative_care_n": sepsis_total["unweighted_events"], "weighted_hospitalizations": round(sepsis_total["weighted_n"]), "weighted_palliative_care_n": round(sepsis_total["weighted_events"]), "weighted_palliative_care_percent": round(100 * sepsis_total["estimate"], 2), "ci_95_lower_percent": round(100 * sepsis_total["ci_lower"], 2), "ci_95_upper_percent": round(100 * sepsis_total["ci_upper"], 2), "subtype_trend_p_value": "—", "overall_year_by_subtype_p_value": format_p(subtype_interaction_p)})
    write_csv(OUTPUT_DIR / "annual_sepsis_palliative_care_by_subtype.csv", subtype_annual)
    write_csv(OUTPUT_DIR / "sepsis_subtype_trend_tests.csv", subtype_trends)

    summary = {
        "command_20a_eapc": overall_eapc,
        "command_20b": {"no_sepsis_eapc": no_sepsis_eapc, "sepsis_eapc": sepsis_eapc, "year_by_sepsis_interaction_p_value": format_p(trend_interaction_p)},
        "command_20c": {"suppression_threshold": "10 or fewer unweighted palliative-care events", "suppressed_cells": suppressed_cells, "year_by_subtype_wald_chi_square": round(interaction_stat, 3), "degrees_of_freedom": interaction_df, "year_by_subtype_interaction_p_value": format_p(subtype_interaction_p), "subtype_trends": subtype_trends},
        "variance_note": "DISCWT-weighted modified-Poisson models with year-specific NIS_STRATUM linearization and discharge-level variance units; not full NIS hospital-cluster-adjusted inference.",
    }
    (OUTPUT_DIR / "phase_11_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    return summary


if __name__ == "__main__":
    main()
