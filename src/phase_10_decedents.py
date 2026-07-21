"""Commands 19A–19C: in-hospital decedent analyses."""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path
from typing import Any, Callable

import duckdb
import numpy as np
import pandas as pd
from scipy.special import expit
from scipy.stats import chi2, norm

from src.phase_6_palliative_care import comparison, format_p, weighted_prevalence
from src.phase_7_subtype_palliative_care import cohort_table
from src.phase_8_primary_adjusted import COHORT_DATABASE, build_model_frame, design_matrix
from src.phase_9_interaction import SUBTYPE_LEVELS, fit_interaction_model, interaction_design_matrix


REPO_ROOT = Path(__file__).resolve().parents[1]
OUTPUT_DIR = REPO_ROOT / "outputs/phase_10"


def subtype_only_design(frame: pd.DataFrame, force_subtype: str | None = None) -> tuple[np.ndarray, list[str]]:
    working = frame.copy() if force_subtype is not None else frame
    if force_subtype is not None:
        working["hm_subtype_label"] = force_subtype
    x, names = design_matrix(working, force_sepsis=0)
    sepsis_index = names.index("Documented sepsis")
    keep = [index for index in range(x.shape[1]) if index != sepsis_index]
    return x[:, keep], [names[index] for index in keep]


def fit_matrix_model(frame: pd.DataFrame, builder: Callable[[pd.DataFrame], tuple[np.ndarray, list[str]]]) -> dict[str, Any]:
    x_full, all_names = builder(frame)
    active = np.array([0] + [index for index in range(1, x_full.shape[1]) if np.any(x_full[:, index] != 0)])
    x = x_full[:, active]
    names = [all_names[index] for index in active]
    y = frame["palliative_care"].to_numpy(float)
    frequency = frame["frequency"].to_numpy(float)
    survey_weight = frame["weight"].to_numpy(float)
    combined_weight = frequency * survey_weight
    beta = np.zeros(x.shape[1])
    for iteration in range(1, 51):
        probability = expit(x @ beta)
        gradient = x.T @ (combined_weight * (y - probability))
        information = (x.T * (combined_weight * probability * (1 - probability))) @ x
        step = np.linalg.solve(information, gradient)
        beta += step
        if np.max(np.abs(step)) < 1e-9:
            break
    else:
        raise RuntimeError("Decedent subtype model failed to converge.")
    probability = expit(x @ beta)
    individual_score = (survey_weight * (y - probability))[:, None] * x
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
    return {"beta": beta, "covariance": bread @ meat @ bread, "names": names, "active": active, "iterations": iteration}


def margin(
    frame: pd.DataFrame, model: dict[str, Any], builder: Callable[..., tuple[np.ndarray, list[str]]], **overrides: Any
) -> tuple[float, np.ndarray]:
    x, _ = builder(frame, **overrides)
    x = x[:, model["active"]]
    weights = frame["frequency"].to_numpy(float) * frame["weight"].to_numpy(float)
    probability = expit(x @ model["beta"])
    estimate = float(np.sum(weights * probability) / weights.sum())
    gradient = np.sum((weights * probability * (1 - probability))[:, None] * x, axis=0) / weights.sum()
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
    died_missing = connection.execute("SELECT count(*) FROM hm_cohort WHERE DIED IS NULL").fetchone()[0]
    frame = build_model_frame(connection)
    connection.close()
    decedents = frame.loc[frame["died"] == 1].reset_index(drop=True)
    sepsis_decedents = decedents.loc[decedents["sepsis"] == 1].reset_index(drop=True)

    decedent_records = decedents.loc[
        decedents.index.repeat(decedents["frequency"])
    ].reset_index(drop=True)
    weights = decedent_records["weight"].to_numpy(float)
    outcome = decedent_records["palliative_care"].to_numpy(float)
    exposure = decedent_records["sepsis"].to_numpy(bool)
    strata = decedent_records["stratum_id"].to_numpy()
    contrast = comparison(weights, outcome, strata, exposure)
    prevalence_rows = []
    for label, domain in [("Decedents without documented sepsis", ~exposure), ("Decedents with documented sepsis", exposure)]:
        result = weighted_prevalence(weights, outcome, strata, domain)
        prevalence_rows.append({
            "cohort": label, "unweighted_sample_n": result["unweighted_n"],
            "unweighted_palliative_care_n": result["unweighted_events"],
            "weighted_hospitalizations": round(result["weighted_n"]),
            "weighted_palliative_care_n": round(result["weighted_events"]),
            "weighted_palliative_care_percent": round(100 * result["estimate"], 2),
            "ci_95_lower_percent": round(100 * result["ci_lower"], 2),
            "ci_95_upper_percent": round(100 * result["ci_upper"], 2),
            "sepsis_comparison_p_value": format_p(contrast["difference_p"]) if not prevalence_rows else "",
        })
    overall_decedent = weighted_prevalence(weights, outcome, strata)
    prevalence_rows.append({
        "cohort": "Total decedents", "unweighted_sample_n": overall_decedent["unweighted_n"],
        "unweighted_palliative_care_n": overall_decedent["unweighted_events"],
        "weighted_hospitalizations": round(overall_decedent["weighted_n"]),
        "weighted_palliative_care_n": round(overall_decedent["weighted_events"]),
        "weighted_palliative_care_percent": round(100 * overall_decedent["estimate"], 2),
        "ci_95_lower_percent": round(100 * overall_decedent["ci_lower"], 2),
        "ci_95_upper_percent": round(100 * overall_decedent["ci_upper"], 2),
        "sepsis_comparison_p_value": format_p(contrast["difference_p"]),
    })
    write_csv(OUTPUT_DIR / "decedent_palliative_care_by_sepsis.csv", prevalence_rows)

    interaction_model = fit_interaction_model(decedents)
    interaction_indices = [i for i, name in enumerate(interaction_model["names"]) if name.startswith("Sepsis ×")]
    interaction_beta = interaction_model["beta"][interaction_indices]
    interaction_covariance = interaction_model["covariance"][np.ix_(interaction_indices, interaction_indices)]
    interaction_df = int(np.linalg.matrix_rank(interaction_covariance))
    interaction_statistic = float(interaction_beta.T @ np.linalg.pinv(interaction_covariance) @ interaction_beta)
    interaction_p = float(chi2.sf(interaction_statistic, interaction_df))
    interaction_rows = []
    for subtype in SUBTYPE_LEVELS:
        p0, g0 = margin(decedents, interaction_model, interaction_design_matrix, force_sepsis=0, force_subtype=subtype)
        p1, g1 = margin(decedents, interaction_model, interaction_design_matrix, force_sepsis=1, force_subtype=subtype)
        diff, gradient = p1 - p0, g1 - g0
        se = math.sqrt(float(gradient @ interaction_model["covariance"] @ gradient))
        interaction_rows.append({
            "hm_subtype": subtype, "adjusted_probability_no_sepsis_percent": round(100 * p0, 2),
            "adjusted_probability_sepsis_percent": round(100 * p1, 2),
            "adjusted_difference_pp": round(100 * diff, 2),
            "difference_ci_95_lower_pp": round(100 * (diff - 1.96 * se), 2),
            "difference_ci_95_upper_pp": round(100 * (diff + 1.96 * se), 2),
            "subtype_p_value": format_p(float(2 * norm.sf(abs(diff / se)))),
            "overall_interaction_p_value": format_p(interaction_p) if subtype == SUBTYPE_LEVELS[0] else "",
        })
    interaction_rows.append({"hm_subtype": "Total", "adjusted_probability_no_sepsis_percent": "—", "adjusted_probability_sepsis_percent": "—", "adjusted_difference_pp": "—", "difference_ci_95_lower_pp": "—", "difference_ci_95_upper_pp": "—", "subtype_p_value": "—", "overall_interaction_p_value": format_p(interaction_p)})
    write_csv(OUTPUT_DIR / "adjusted_decedent_interaction.csv", interaction_rows)

    subtype_ids = {
        "Lymphoma": "lymphoma", "AML": "aml", "CML": "cml",
        "CLL/chronic leukemia": "cll_chronic_leukemia",
        "ALL/unspecified acute leukemia": "all", "Other leukemia": "other_leukemia",
        "Myeloma/plasma-cell neoplasm": "myeloma_plasma_cell", "MDS": "mds", "MPN": "mpn",
    }
    sepsis_decedents_for_table = sepsis_decedents.copy()
    sepsis_decedents_for_table["hm_subtype"] = sepsis_decedents_for_table["hm_subtype_label"].map(subtype_ids)
    sepsis_decedents_for_table = sepsis_decedents_for_table.loc[
        sepsis_decedents_for_table.index.repeat(sepsis_decedents_for_table["frequency"])
    ].reset_index(drop=True)
    stratum_parts = sepsis_decedents_for_table["stratum_id"].str.split(":", n=1, expand=True)
    sepsis_decedents_for_table["year"] = stratum_parts[0]
    sepsis_decedents_for_table["nis_stratum"] = stratum_parts[1]
    unadjusted_rows, unadjusted_test = cohort_table(sepsis_decedents_for_table, "Sepsis decedents")
    sepsis_total = prevalence_rows[1]
    unadjusted_rows.append({
        "cohort": "Sepsis decedents", "hm_subtype": "Total",
        "unweighted_sample_n": sepsis_total["unweighted_sample_n"],
        "unweighted_palliative_care_n": sepsis_total["unweighted_palliative_care_n"],
        "weighted_hospitalizations": sepsis_total["weighted_hospitalizations"],
        "weighted_palliative_care_n": sepsis_total["weighted_palliative_care_n"],
        "weighted_palliative_care_percent": sepsis_total["weighted_palliative_care_percent"],
        "ci_95_lower_percent": sepsis_total["ci_95_lower_percent"],
        "ci_95_upper_percent": sepsis_total["ci_95_upper_percent"],
        "overall_p_value": unadjusted_test["overall_p_value"],
    })
    write_csv(OUTPUT_DIR / "sepsis_decedent_subtype_unadjusted.csv", unadjusted_rows)
    subtype_model = fit_matrix_model(sepsis_decedents, subtype_only_design)
    subtype_indices = [i for i, name in enumerate(subtype_model["names"]) if name.startswith("hm_subtype_label:")]
    subtype_beta = subtype_model["beta"][subtype_indices]
    subtype_cov = subtype_model["covariance"][np.ix_(subtype_indices, subtype_indices)]
    subtype_df = int(np.linalg.matrix_rank(subtype_cov))
    subtype_statistic = float(subtype_beta.T @ np.linalg.pinv(subtype_cov) @ subtype_beta)
    subtype_p = float(chi2.sf(subtype_statistic, subtype_df))
    adjusted_sepsis_rows = []
    for subtype in SUBTYPE_LEVELS:
        estimate, gradient = margin(sepsis_decedents, subtype_model, subtype_only_design, force_subtype=subtype)
        se = math.sqrt(float(gradient @ subtype_model["covariance"] @ gradient))
        adjusted_sepsis_rows.append({
            "hm_subtype": subtype, "adjusted_probability_percent": round(100 * estimate, 2),
            "ci_95_lower_percent": round(100 * max(0, estimate - 1.96 * se), 2),
            "ci_95_upper_percent": round(100 * min(1, estimate + 1.96 * se), 2),
            "overall_subtype_p_value": format_p(subtype_p) if subtype == SUBTYPE_LEVELS[0] else "",
        })
    adjusted_sepsis_rows.append({"hm_subtype": "Total", "adjusted_probability_percent": "—", "ci_95_lower_percent": "—", "ci_95_upper_percent": "—", "overall_subtype_p_value": format_p(subtype_p)})
    write_csv(OUTPUT_DIR / "sepsis_decedent_subtype_adjusted.csv", adjusted_sepsis_rows)

    summary = {
        "died_missing_records": int(died_missing),
        "decedent_unweighted_n": int(decedents["frequency"].sum()),
        "sepsis_decedent_unweighted_n": int(sepsis_decedents["frequency"].sum()),
        "command_19a_p_value": format_p(contrast["difference_p"]),
        "command_19b_interaction_test": {"wald_chi_square": round(interaction_statistic, 3), "degrees_of_freedom": interaction_df, "p_value": format_p(interaction_p)},
        "command_19c_subtype_test": {"wald_chi_square": round(subtype_statistic, 3), "degrees_of_freedom": subtype_df, "p_value": format_p(subtype_p)},
        "variance_note": "DISCWT-weighted estimates with year-specific NIS_STRATUM linearization and discharge-level variance units; not full NIS hospital-cluster-adjusted inference.",
    }
    (OUTPUT_DIR / "phase_10_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    return summary


if __name__ == "__main__":
    main()
