"""Adjusted covariate effects on PC use within septic-shock strata."""

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
from src.phase_8_primary_adjusted import CATEGORY_SPECS, build_model_frame


ROOT = Path(__file__).resolve().parents[1]
COHORT_DATABASE = ROOT / "outputs/septic_shock/hm_cohort_septic_shock.duckdb"
OUTPUT_DIR = ROOT / "outputs/septic_shock/stratified_covariate_effects"

VARIABLE_LABELS = {
    "age": "Age, per year",
    "sex": "Sex",
    "race": "Race/ethnicity",
    "payer": "Primary payer",
    "income": "Income quartile",
    "cci_category": "Cancer-excluded CCI category",
    "region": "Hospital region",
    "location_teaching": "Hospital location/teaching status",
    "bed_size": "Hospital bed size",
    "year_category": "Admission year",
    "hm_subtype_label": "HM subtype",
}


def design_matrix(frame: pd.DataFrame) -> tuple[np.ndarray, list[str], dict[str, list[int]]]:
    columns = [np.ones(len(frame)), frame["age"].to_numpy(float) - 65.0]
    names = ["Intercept", "age"]
    blocks: dict[str, list[int]] = {"age": [1]}
    for variable, levels in CATEGORY_SPECS.items():
        values = frame[variable].astype(str).to_numpy()
        blocks[variable] = []
        for level in levels[1:]:
            blocks[variable].append(len(columns))
            columns.append((values == level).astype(float))
            names.append(f"{variable}: {level}")
    return np.column_stack(columns), names, blocks


def fit_model(frame: pd.DataFrame) -> dict[str, Any]:
    full, all_names, full_blocks = design_matrix(frame)
    active = np.array([0] + [i for i in range(1, full.shape[1]) if np.any(full[:, i] != 0)])
    x = full[:, active]
    names = [all_names[i] for i in active]
    old_to_new = {old: new for new, old in enumerate(active)}
    blocks = {key: [old_to_new[i] for i in indices if i in old_to_new] for key, indices in full_blocks.items()}
    y = frame["palliative_care"].to_numpy(float)
    frequency = frame["frequency"].to_numpy(float)
    survey_weight = frame["weight"].to_numpy(float)
    combined = frequency * survey_weight
    beta = np.zeros(x.shape[1])
    for iteration in range(1, 51):
        probability = expit(x @ beta)
        information = (x.T * (combined * probability * (1 - probability))) @ x
        step = np.linalg.solve(information, x.T @ (combined * (y - probability)))
        beta += step
        if np.max(np.abs(step)) < 1e-9:
            break
    else:
        raise RuntimeError("Stratified model failed to converge.")
    probability = expit(x @ beta)
    scores = (survey_weight * (y - probability))[:, None] * x
    meat = np.zeros((x.shape[1], x.shape[1]))
    groups = pd.DataFrame({"stratum": frame["stratum_id"], "row": np.arange(len(frame))}).groupby("stratum").indices
    for indices in groups.values():
        idx = np.asarray(indices)
        n = frequency[idx].sum()
        if n <= 1:
            continue
        score_sum = (scores[idx] * frequency[idx, None]).sum(axis=0)
        sum_outer = (scores[idx].T * frequency[idx]) @ scores[idx]
        meat += n / (n - 1) * (sum_outer - np.outer(score_sum, score_sum) / n)
    bread = np.linalg.inv(information)
    return {"beta": beta, "covariance": bread @ meat @ bread, "names": names,
            "blocks": blocks, "iterations": iteration}


def joint_test(model: dict[str, Any], indices: list[int]) -> tuple[float, int, float]:
    beta = model["beta"][indices]
    covariance = model["covariance"][np.ix_(indices, indices)]
    degrees = int(np.linalg.matrix_rank(covariance))
    statistic = float(beta.T @ np.linalg.pinv(covariance) @ beta)
    return statistic, degrees, float(chi2.sf(statistic, degrees))


def coefficient_result(model: dict[str, Any], index: int) -> tuple[float, float, float, float]:
    coefficient = model["beta"][index]
    se = math.sqrt(model["covariance"][index, index])
    return (math.exp(coefficient), math.exp(coefficient - 1.96 * se),
            math.exp(coefficient + 1.96 * se), float(2 * norm.sf(abs(coefficient / se))))


def result_rows(frame: pd.DataFrame, model: dict[str, Any], cohort: str) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    age_index = model["blocks"]["age"][0]
    estimate, low, high, p = coefficient_result(model, age_index)
    rows.append({"cohort": cohort, "variable": VARIABLE_LABELS["age"], "level": "Per 1-year increase",
                 "reference": "—", "adjusted_odds_ratio": round(estimate, 3),
                 "ci_95_lower": round(low, 3), "ci_95_upper": round(high, 3),
                 "level_p_value": format_p(p), "overall_joint_wald_chi_square": "—",
                 "joint_degrees_of_freedom": "—", "overall_joint_p_value": format_p(p)})
    for variable, levels in CATEGORY_SPECS.items():
        indices = model["blocks"][variable]
        statistic, degrees, overall_p = joint_test(model, indices)
        rows.append({"cohort": cohort, "variable": VARIABLE_LABELS[variable], "level": levels[0],
                     "reference": "Reference", "adjusted_odds_ratio": "1.000", "ci_95_lower": "—",
                     "ci_95_upper": "—", "level_p_value": "—",
                     "overall_joint_wald_chi_square": round(statistic, 3),
                     "joint_degrees_of_freedom": degrees, "overall_joint_p_value": format_p(overall_p)})
        index_by_name = {name: i for i, name in enumerate(model["names"])}
        for level in levels[1:]:
            name = f"{variable}: {level}"
            if name not in index_by_name:
                rows.append({"cohort": cohort, "variable": "", "level": level, "reference": "Absent",
                             "adjusted_odds_ratio": "—", "ci_95_lower": "—", "ci_95_upper": "—",
                             "level_p_value": "—", "overall_joint_wald_chi_square": "",
                             "joint_degrees_of_freedom": "", "overall_joint_p_value": ""})
                continue
            estimate, low, high, p = coefficient_result(model, index_by_name[name])
            rows.append({"cohort": cohort, "variable": "", "level": level, "reference": "",
                         "adjusted_odds_ratio": round(estimate, 3), "ci_95_lower": round(low, 3),
                         "ci_95_upper": round(high, 3), "level_p_value": format_p(p),
                         "overall_joint_wald_chi_square": "", "joint_degrees_of_freedom": "",
                         "overall_joint_p_value": ""})
    rows.append({"cohort": cohort, "variable": "Total", "level": "All included hospitalizations",
                 "reference": "—", "adjusted_odds_ratio": "—", "ci_95_lower": "—", "ci_95_upper": "—",
                 "level_p_value": "—", "overall_joint_wald_chi_square": "—",
                 "joint_degrees_of_freedom": "—", "overall_joint_p_value": "—"})
    return rows


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)


def table_markdown(rows: list[dict[str, Any]], cohort: str, n: int, weighted_n: int, events: int) -> list[str]:
    lines = [f"## {cohort}", "",
             f"The model included {n:,} sampled hospitalizations (weighted {weighted_n:,}) and {events:,} sampled documented inpatient palliative-care events.", "",
             "| Variable | Level | Adjusted OR | 95% CI | Level p-value | Overall joint Wald χ² (df) | Overall p-value |",
             "| --- | --- | ---: | ---: | ---: | ---: | ---: |"]
    for row in rows:
        ci = "—" if row["ci_95_lower"] == "—" else f'{float(row["ci_95_lower"]):.3f}–{float(row["ci_95_upper"]):.3f}'
        joint = "—" if row["overall_joint_wald_chi_square"] in ("", "—") else f'{float(row["overall_joint_wald_chi_square"]):.3f} ({row["joint_degrees_of_freedom"]})'
        lines.append(f'| {row["variable"]} | {row["level"]} | {row["adjusted_odds_ratio"]} | {ci} | {row["level_p_value"]} | {joint} | {row["overall_joint_p_value"]} |')
    return lines


def write_report(summary: dict[str, Any], rows_by_cohort: dict[str, list[dict[str, Any]]]) -> Path:
    lines = ["# Adjusted Covariate Associations With Documented Inpatient Palliative-Care Use by Septic-Shock Status", "",
             "**Outcome:** Documented inpatient palliative-care use (`Z51.5`) in any diagnosis position.", "",
             "Separate survey-weighted logistic-regression models were fitted among hospitalizations with and without documented septic shock (`R65.21`). Each model adjusted simultaneously for age, sex, race/ethnicity, primary payer, income quartile, cancer-excluded Charlson category, hospital region, hospital location/teaching status, hospital bed size, admission year, and mutually exclusive HM subtype.", "",
             "For categorical variables, the overall p-value is a joint Wald test of all displayed nonreference levels. Level-specific odds ratios compare each level with the stated reference category.", ""]
    for cohort in ["No septic shock", "Septic shock"]:
        info = summary["cohorts"][cohort]
        lines.extend(table_markdown(rows_by_cohort[cohort], cohort, info["unweighted_n"], info["weighted_n"], info["unweighted_pc_events"]))
        lines.append("")
    lines.extend(["## Interpretation and reporting notes", "",
                  "- Adjusted odds ratios describe conditional associations within each septic-shock stratum and should not be interpreted as causal effects.",
                  "- Differences in coefficient magnitude between the two separate models are descriptive; formal evidence that an association differs by septic-shock status requires a covariate-by-septic-shock interaction test in a pooled model.",
                  "- Estimates use `DISCWT` and year-specific `NIS_STRATUM`; `HOSP_NIS` is unavailable by study decision.",
                  "- The tables include admission year and HM subtype because both are covariates in the primary palliative-care model.", ""])
    path = OUTPUT_DIR / "stratified_covariate_effects_report.md"
    path.write_text("\n".join(lines), encoding="utf-8")
    return path


def main() -> dict[str, Any]:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    connection = duckdb.connect(str(COHORT_DATABASE), read_only=True)
    frame = build_model_frame(connection)
    connection.close()
    rows_by_cohort: dict[str, list[dict[str, Any]]] = {}
    summary: dict[str, Any] = {"cohorts": {}, "variance_note": "DISCWT-weighted models with year-specific NIS_STRATUM linearization and discharge-level variance units."}
    all_rows = []
    for shock, cohort in [(0, "No septic shock"), (1, "Septic shock")]:
        subset = frame[frame["sepsis"].eq(shock)].reset_index(drop=True)
        model = fit_model(subset)
        rows = result_rows(subset, model, cohort)
        rows_by_cohort[cohort] = rows; all_rows.extend(rows)
        summary["cohorts"][cohort] = {"unweighted_n": int(subset.frequency.sum()),
            "weighted_n": int(round((subset.frequency * subset.weight).sum())),
            "unweighted_pc_events": int((subset.frequency * subset.palliative_care).sum()),
            "iterations": model["iterations"]}
    write_csv(OUTPUT_DIR / "stratified_covariate_adjusted_odds_ratios.csv", all_rows)
    report = write_report(summary, rows_by_cohort)
    summary["report_path"] = str(report)
    (OUTPUT_DIR / "stratified_covariate_effects_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2)); return summary


if __name__ == "__main__":
    main()
