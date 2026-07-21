"""Commands 17A and 17B: sepsis-by-HM-subtype interaction model."""

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
from scipy.stats import chi2, norm

from src.phase_6_palliative_care import format_p
from src.phase_8_primary_adjusted import (
    COHORT_DATABASE,
    CATEGORY_SPECS,
    build_model_frame,
    design_matrix,
)


REPO_ROOT = Path(__file__).resolve().parents[1]
OUTPUT_DIR = REPO_ROOT / "outputs/phase_9"
SUBTYPE_LEVELS = CATEGORY_SPECS["hm_subtype_label"]


def interaction_design_matrix(
    frame: pd.DataFrame,
    force_sepsis: int | None = None,
    force_subtype: str | None = None,
) -> tuple[np.ndarray, list[str]]:
    working = frame.copy() if force_subtype is not None else frame
    if force_subtype is not None:
        working = frame.copy()
        working["hm_subtype_label"] = force_subtype
    base, names = design_matrix(working, force_sepsis=force_sepsis)
    sepsis_index = names.index("Documented sepsis")
    subtype_indices = [names.index(f"hm_subtype_label: {level}") for level in SUBTYPE_LEVELS[1:]]
    interaction_columns = [base[:, sepsis_index] * base[:, index] for index in subtype_indices]
    interaction_names = [f"Sepsis × {level}" for level in SUBTYPE_LEVELS[1:]]
    return np.column_stack([base, *interaction_columns]), names + interaction_names


def fit_interaction_model(frame: pd.DataFrame, tolerance: float = 1e-9, max_iterations: int = 50) -> dict[str, Any]:
    x_full, all_names = interaction_design_matrix(frame)
    active = np.array([0] + [index for index in range(1, x_full.shape[1]) if np.any(x_full[:, index] != 0)])
    x = x_full[:, active]
    names = [all_names[index] for index in active]
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
        raise RuntimeError("Interaction model failed to converge.")

    probability = expit(x @ beta)
    individual_score = (survey_weight * (y - probability))[:, None] * x
    strata = frame["stratum_id"].to_numpy()
    meat = np.zeros((x.shape[1], x.shape[1]))
    grouped = pd.DataFrame({"stratum": strata, "row": np.arange(len(frame))}).groupby("stratum", sort=False).indices
    for indices in grouped.values():
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
    return {"beta": beta, "covariance": covariance, "names": names, "active": active, "iterations": iteration}


def standardized_probability(
    frame: pd.DataFrame,
    model: dict[str, Any],
    sepsis: int,
    subtype: str,
) -> tuple[float, np.ndarray]:
    x, _ = interaction_design_matrix(frame, force_sepsis=sepsis, force_subtype=subtype)
    x = x[:, model["active"]]
    weights = frame["frequency"].to_numpy(float) * frame["weight"].to_numpy(float)
    probability = expit(x @ model["beta"])
    total_weight = weights.sum()
    estimate = float(np.sum(weights * probability) / total_weight)
    gradient = np.sum((weights * probability * (1 - probability))[:, None] * x, axis=0) / total_weight
    return estimate, gradient


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
    total_records = connection.execute("SELECT count(*) FROM hm_cohort").fetchone()[0]
    frame = build_model_frame(connection)
    connection.close()
    model = fit_interaction_model(frame)
    covariance = model["covariance"]

    interaction_indices = [index for index, name in enumerate(model["names"]) if name.startswith("Sepsis ×")]
    interaction_beta = model["beta"][interaction_indices]
    interaction_covariance = covariance[np.ix_(interaction_indices, interaction_indices)]
    interaction_df = int(np.linalg.matrix_rank(interaction_covariance))
    interaction_statistic = float(interaction_beta.T @ np.linalg.pinv(interaction_covariance) @ interaction_beta)
    interaction_p = float(chi2.sf(interaction_statistic, interaction_df))

    rows: list[dict[str, Any]] = []
    figure_rows: list[dict[str, Any]] = []
    for subtype in SUBTYPE_LEVELS:
        probability_0, gradient_0 = standardized_probability(frame, model, 0, subtype)
        probability_1, gradient_1 = standardized_probability(frame, model, 1, subtype)
        difference = probability_1 - probability_0
        gradient_difference = gradient_1 - gradient_0
        difference_se = math.sqrt(float(gradient_difference @ covariance @ gradient_difference))
        p_value = float(2 * norm.sf(abs(difference / difference_se)))
        rows.append({
            "hm_subtype": subtype,
            "adjusted_probability_no_sepsis_percent": round(100 * probability_0, 2),
            "adjusted_probability_sepsis_percent": round(100 * probability_1, 2),
            "adjusted_absolute_difference_pp": round(100 * difference, 2),
            "difference_ci_95_lower_pp": round(100 * (difference - 1.96 * difference_se), 2),
            "difference_ci_95_upper_pp": round(100 * (difference + 1.96 * difference_se), 2),
            "subtype_p_value": format_p(p_value),
            "overall_interaction_p_value": format_p(interaction_p) if subtype == SUBTYPE_LEVELS[0] else "",
        })
        for status, estimate, gradient in [("No documented sepsis", probability_0, gradient_0), ("Documented sepsis", probability_1, gradient_1)]:
            se = math.sqrt(float(gradient @ covariance @ gradient))
            figure_rows.append({
                "hm_subtype": subtype, "sepsis_status": status,
                "adjusted_probability_percent": round(100 * estimate, 2),
                "ci_95_lower_percent": round(100 * max(0, estimate - 1.96 * se), 2),
                "ci_95_upper_percent": round(100 * min(1, estimate + 1.96 * se), 2),
            })
    rows.append({
        "hm_subtype": "Total", "adjusted_probability_no_sepsis_percent": "—",
        "adjusted_probability_sepsis_percent": "—", "adjusted_absolute_difference_pp": "—",
        "difference_ci_95_lower_pp": "—", "difference_ci_95_upper_pp": "—",
        "subtype_p_value": "—", "overall_interaction_p_value": format_p(interaction_p),
    })
    write_csv(OUTPUT_DIR / "subtype_adjusted_probabilities.csv", rows)
    write_csv(OUTPUT_DIR / "interaction_figure_data.csv", figure_rows)
    interaction_terms = []
    for index in interaction_indices:
        se = math.sqrt(covariance[index, index])
        interaction_terms.append({
            "term": model["names"][index], "interaction_odds_ratio": round(math.exp(model["beta"][index]), 3),
            "ci_95_lower": round(math.exp(model["beta"][index] - 1.96 * se), 3),
            "ci_95_upper": round(math.exp(model["beta"][index] + 1.96 * se), 3),
            "p_value": format_p(float(2 * norm.sf(abs(model["beta"][index] / se)))),
        })
    interaction_terms.append({"term": "Total", "interaction_odds_ratio": "—", "ci_95_lower": "—", "ci_95_upper": "—", "p_value": "—"})
    write_csv(OUTPUT_DIR / "interaction_terms_audit.csv", interaction_terms)
    summary = {
        "included_unweighted_records": int(frame["frequency"].sum()),
        "excluded_records": int(total_records - frame["frequency"].sum()),
        "iterations": model["iterations"],
        "joint_interaction_test": {
            "wald_chi_square": round(interaction_statistic, 3),
            "degrees_of_freedom": interaction_df,
            "overall_interaction_p_value": format_p(interaction_p),
        },
        "variance_note": "DISCWT-weighted interaction model with year-specific NIS_STRATUM linearization and discharge-level variance units; not full NIS hospital-cluster-adjusted inference.",
        "subtype_estimates": rows,
    }
    (OUTPUT_DIR / "phase_9_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    return summary


if __name__ == "__main__":
    main()
