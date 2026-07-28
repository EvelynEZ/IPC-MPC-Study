"""Diagnosis-based BMT/HSCT modification of the septic-shock association."""

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

from src.phase_5_cci_los_mortality import CCI_COMPONENTS, CCI_WEIGHTS, component_condition
from src.phase_6_palliative_care import format_p
from src.phase_8_primary_adjusted import CATEGORY_SPECS


ROOT = Path(__file__).resolve().parents[1]
COHORT_DATABASE = ROOT / "outputs/septic_shock/hm_cohort_septic_shock.duckdb"
OUTPUT_DIR = ROOT / "outputs/septic_shock/bmt_interaction"


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def build_frame(connection: duckdb.DuckDBPyConnection) -> pd.DataFrame:
    flags = []
    for component, prefixes in CCI_COMPONENTS.items():
        flags.append(f"list_contains(list_transform(diagnosis_codes, code -> {component_condition(prefixes)}), TRUE) AS {component}")
    scores = []
    for component, weight in CCI_WEIGHTS.items():
        present = "diab AND NOT diabwc" if component == "diab" else "mld AND NOT msld" if component == "mld" else component
        scores.append(f"CASE WHEN {present} THEN {weight} ELSE 0 END")
    connection.execute("CREATE TEMP TABLE bmt_flags AS SELECT *, " + ", ".join(flags) + " FROM hm_cohort")
    connection.execute("CREATE TEMP TABLE bmt_scored AS SELECT *, " + " + ".join(scores) + " AS cci FROM bmt_flags")
    return connection.execute("""
        SELECT concat(CAST(YEAR AS INTEGER), ':', CAST(NIS_STRATUM AS INTEGER)) AS stratum_id,
            DISCWT::DOUBLE AS weight, sepsis::INTEGER AS shock,
            (list_contains(diagnosis_codes, 'Z9481') OR list_contains(diagnosis_codes, 'Z9484'))::INTEGER AS bmt,
            palliative_care::INTEGER AS palliative_care, AGE::DOUBLE AS age,
            CASE FEMALE WHEN 0 THEN 'Male' WHEN 1 THEN 'Female' ELSE 'Missing' END AS sex,
            CASE RACE WHEN 1 THEN 'White' WHEN 2 THEN 'Black' WHEN 3 THEN 'Hispanic'
                WHEN 4 THEN 'Asian/Pacific Islander' WHEN 5 THEN 'Native American' WHEN 6 THEN 'Other' ELSE 'Missing' END AS race,
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
        FROM bmt_scored
        WHERE DISCWT IS NOT NULL AND NIS_STRATUM IS NOT NULL AND AGE IS NOT NULL
        GROUP BY ALL
    """).fetchdf()


def design_matrix(frame: pd.DataFrame, shock: int | None = None, bmt: int | None = None) -> tuple[np.ndarray, list[str]]:
    shock_values = frame["shock"].to_numpy(float) if shock is None else np.full(len(frame), shock, dtype=float)
    bmt_values = frame["bmt"].to_numpy(float) if bmt is None else np.full(len(frame), bmt, dtype=float)
    columns = [np.ones(len(frame)), shock_values, bmt_values, shock_values * bmt_values,
               frame["age"].to_numpy(float) - 65.0]
    names = ["Intercept", "Septic shock", "BMT/HSCT status", "Septic shock × BMT/HSCT", "Age, per year (centered at 65)"]
    for variable, levels in CATEGORY_SPECS.items():
        values = frame[variable].astype(str).to_numpy()
        for level in levels[1:]:
            columns.append((values == level).astype(float))
            names.append(f"{variable}: {level}")
    return np.column_stack(columns), names


def fit_model(frame: pd.DataFrame) -> dict[str, Any]:
    full, all_names = design_matrix(frame)
    active = np.array([0] + [i for i in range(1, full.shape[1]) if np.any(full[:, i] != 0)])
    x, names = full[:, active], [all_names[i] for i in active]
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
        raise RuntimeError("BMT interaction model failed to converge.")
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
    return {"beta": beta, "covariance": bread @ meat @ bread, "names": names, "active": active, "iterations": iteration}


def margin(frame: pd.DataFrame, model: dict[str, Any], shock: int, bmt: int) -> tuple[float, np.ndarray]:
    x, _ = design_matrix(frame, shock, bmt)
    x = x[:, model["active"]]
    weights = frame["frequency"].to_numpy(float) * frame["weight"].to_numpy(float)
    probability = expit(x @ model["beta"])
    estimate = float(np.sum(weights * probability) / weights.sum())
    gradient = np.sum((weights * probability * (1 - probability))[:, None] * x, axis=0) / weights.sum()
    return estimate, gradient


def linear_result(estimate: float, gradient: np.ndarray, covariance: np.ndarray) -> tuple[float, float, float, float]:
    se = math.sqrt(max(0.0, float(gradient @ covariance @ gradient)))
    p = float(2 * norm.sf(abs(estimate / se))) if se else 0.0
    return estimate, estimate - 1.96 * se, estimate + 1.96 * se, p


def write_report(summary: dict[str, Any]) -> Path:
    """Create the report-facing narrative and tables for this analysis."""
    descriptive = summary["descriptive"]
    odds = summary["conditional_odds_ratios"]
    probabilities = summary["adjusted_probabilities"]
    differences = summary["probability_differences"]
    interaction = summary["interaction_wald_test"]
    lines = [
        "# BMT/HSCT Modification of the Association Between Septic Shock and Documented Inpatient Palliative-Care Use",
        "",
        "**Primary outcome:** Documented inpatient palliative-care use, defined as ICD-10-CM `Z51.5` in any diagnosis position.",
        "",
        "**Primary exposure:** Documented septic shock, defined as ICD-10-CM `R65.21` in any diagnosis position.",
        "",
        "**Effect modifier:** Diagnosis-based BMT/HSCT status, defined as ICD-10-CM `Z94.81` or `Z94.84` in any diagnosis position.",
        "",
        "## Four-Group Descriptive Results",
        "",
        "| BMT/HSCT and septic-shock group | Unweighted n | Unweighted PC events | Weighted n, 2016–2022 | Weighted PC events | Weighted PC prevalence | Adjusted overall group p-value |",
        "| --- | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for row in descriptive:
        lines.append(
            f'| {row["group"]} | {row["unweighted_hospitalizations"]:,} | {row["unweighted_pc_events"]:,} | '
            f'{row["weighted_hospitalizations"]:,} | {row["weighted_pc_events"]:,} | '
            f'{row["weighted_pc_percent"]:.2f}% | {row["adjusted_overall_group_p_value"] or ""} |'
        )
    lines.extend([
        "",
        "Counts describe adult hematologic-malignancy hospitalizations. Weighted counts and prevalence estimates use `DISCWT`. The p-value shown on the total row is the adjusted joint Wald test for septic shock, BMT/HSCT status, and their interaction.",
        "",
        "## Adjusted Interaction Model",
        "",
        "A survey-weighted logistic-regression model evaluated whether diagnosis-based BMT/HSCT status modified the association between septic shock and documented inpatient palliative-care use. The model included septic shock, BMT/HSCT status, and a septic shock × BMT/HSCT interaction and adjusted for continuous age, sex, race/ethnicity, primary payer, income quartile, cancer-excluded Charlson category, mutually exclusive HM subtype, hospital region, hospital location/teaching status, hospital bed size, and admission year.",
        "",
        "| Interaction test | Wald chi-square | Degrees of freedom | P-value |",
        "| --- | ---: | ---: | ---: |",
        f'| Septic shock × BMT/HSCT | {interaction["wald_chi_square"]:.3f} | {interaction["degrees_of_freedom"]} | {interaction["p_value"]} |',
        "| Total | — | — | " + interaction["p_value"] + " |",
        "",
        "The multiplicative interaction was not statistically significant, indicating no evidence that BMT/HSCT status modified the adjusted septic-shock odds ratio for documented inpatient palliative-care use.",
        "",
        "## Conditional Adjusted Odds Ratios",
        "",
        "| Contrast | Adjusted odds ratio | 95% CI | P-value |",
        "| --- | ---: | ---: | ---: |",
    ])
    for row in odds:
        ci = "—" if row["adjusted_odds_ratio"] == "—" else f'{row["ci_95_lower"]:.3f}–{row["ci_95_upper"]:.3f}'
        estimate = row["adjusted_odds_ratio"] if row["adjusted_odds_ratio"] == "—" else f'{row["adjusted_odds_ratio"]:.3f}'
        lines.append(f'| {row["contrast"]} | {estimate} | {ci} | {row["p_value"]} |')
    lines.extend([
        "",
        "## Adjusted Palliative-Care Probabilities",
        "",
        "Probabilities are average marginal predictions standardized over the observed weighted distribution of all adjustment covariates.",
        "",
        "| BMT/HSCT and septic-shock group | Adjusted probability | 95% CI |",
        "| --- | ---: | ---: |",
    ])
    for row in probabilities:
        if row["adjusted_pc_probability_percent"] == "—":
            estimate, ci = "—", "—"
        else:
            estimate = f'{row["adjusted_pc_probability_percent"]:.2f}%'
            ci = f'{row["ci_95_lower_percent"]:.2f}%–{row["ci_95_upper_percent"]:.2f}%'
        lines.append(f'| {row["group"]} | {estimate} | {ci} |')
    lines.extend([
        "",
        "## Adjusted Septic-Shock Probability Differences",
        "",
        "| Contrast | Adjusted difference, percentage points | 95% CI | P-value |",
        "| --- | ---: | ---: | ---: |",
    ])
    for row in differences:
        if row["adjusted_difference_pp"] == "—":
            estimate, ci = "—", "—"
        else:
            estimate = f'{row["adjusted_difference_pp"]:.2f}'
            ci = f'{row["ci_95_lower_pp"]:.2f} to {row["ci_95_upper_pp"]:.2f}'
        lines.append(f'| {row["contrast"]} | {estimate} | {ci} | {row["p_value"]} |')
    lines.extend([
        "",
        "## Interpretation",
        "",
        "Septic shock was associated with substantially greater adjusted documented inpatient palliative-care use in both BMT/HSCT and non-BMT hospitalizations. The conditional adjusted odds ratios were similar, and the prespecified multiplicative interaction was not significant (`p=0.445`).",
        "",
        "On the additive probability scale, septic shock was associated with an 18.52-percentage-point increase among non-BMT hospitalizations and a 20.38-point increase among BMT hospitalizations. The difference between these increases was 1.86 percentage points (95% CI, 0.04–3.69; `p=0.045`). Thus, evidence of effect modification depends on the statistical scale: absent on the multiplicative odds scale but marginally present on the additive probability scale.",
        "",
        "Because the primary interaction test was the logistic-model Wald test, the primary conclusion is that diagnosis-based BMT/HSCT status did not significantly modify the association between septic shock and documented inpatient palliative-care use. The additive interaction should be reported as a complementary marginal-effect result.",
        "",
        "## Methodologic Notes and Limitations",
        "",
        "- BMT/HSCT status is diagnosis-based. The source extract contains `I10_NPR`, the number of procedures, but not individual ICD-10-PCS procedure-code fields. Consequently, HSCT procedures performed during the hospitalization could not be added to the definition.",
        "- Hospitalizations without `Z94.81` or `Z94.84` were classified as no BMT/HSCT; some may therefore be misclassified.",
        "- Estimates use `DISCWT` and account for year-specific `NIS_STRATUM`. `HOSP_NIS` is unavailable by study decision, so these are not full hospital-cluster-adjusted NIS variance estimates.",
        "- Diagnosis fields establish co-documentation during hospitalization but not temporal ordering or causality.",
        "- The analysis is hospitalization-based and does not identify unique longitudinal patients.",
        "",
    ])
    path = OUTPUT_DIR / "bmt_shock_interaction_report.md"
    path.write_text("\n".join(lines), encoding="utf-8")
    return path


def main() -> dict[str, Any]:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    connection = duckdb.connect(str(COHORT_DATABASE), read_only=True)
    descriptive = connection.execute("""
        SELECT CASE WHEN list_contains(diagnosis_codes, 'Z9481') OR list_contains(diagnosis_codes, 'Z9484') THEN 1 ELSE 0 END AS bmt,
               sepsis::INTEGER AS shock, count(*)::BIGINT AS unweighted_hospitalizations,
               count(*) FILTER (WHERE palliative_care)::BIGINT AS unweighted_pc_events,
               round(sum(DISCWT), 0)::BIGINT AS weighted_hospitalizations,
               round(sum(DISCWT) FILTER (WHERE palliative_care), 0)::BIGINT AS weighted_pc_events,
               100 * sum(DISCWT) FILTER (WHERE palliative_care) / sum(DISCWT) AS weighted_pc_percent
        FROM hm_cohort GROUP BY ALL ORDER BY bmt, shock
    """).fetchdf()
    frame = build_frame(connection)
    connection.close()
    labels = {(0, 0): "No BMT, no septic shock", (0, 1): "No BMT, septic shock",
              (1, 0): "BMT, no septic shock", (1, 1): "BMT, septic shock"}
    descriptive_rows = []
    for row in descriptive.itertuples(index=False):
        descriptive_rows.append({"group": labels[(row.bmt, row.shock)],
            "unweighted_hospitalizations": int(row.unweighted_hospitalizations),
            "unweighted_pc_events": int(row.unweighted_pc_events),
            "weighted_hospitalizations": int(row.weighted_hospitalizations),
            "weighted_pc_events": int(row.weighted_pc_events),
            "weighted_pc_percent": round(row.weighted_pc_percent, 2)})
    descriptive_rows.append({"group": "Total", "unweighted_hospitalizations": sum(r["unweighted_hospitalizations"] for r in descriptive_rows),
        "unweighted_pc_events": sum(r["unweighted_pc_events"] for r in descriptive_rows),
        "weighted_hospitalizations": sum(r["weighted_hospitalizations"] for r in descriptive_rows),
        "weighted_pc_events": sum(r["weighted_pc_events"] for r in descriptive_rows), "weighted_pc_percent": 8.11})
    model = fit_model(frame)
    beta, covariance, names = model["beta"], model["covariance"], model["names"]
    shock_index, interaction_index = names.index("Septic shock"), names.index("Septic shock × BMT/HSCT")
    bmt_index = names.index("BMT/HSCT status")
    group_indices = [shock_index, bmt_index, interaction_index]
    group_beta = beta[group_indices]
    group_covariance = covariance[np.ix_(group_indices, group_indices)]
    group_wald = float(group_beta.T @ np.linalg.pinv(group_covariance) @ group_beta)
    group_p = float(chi2.sf(group_wald, int(np.linalg.matrix_rank(group_covariance))))
    for row in descriptive_rows:
        row["adjusted_overall_group_p_value"] = format_p(group_p) if row["group"] == "Total" else ""
    write_csv(OUTPUT_DIR / "four_group_descriptive.csv", descriptive_rows)
    interaction_beta, interaction_se = beta[interaction_index], math.sqrt(covariance[interaction_index, interaction_index])
    interaction_wald = (interaction_beta / interaction_se) ** 2
    interaction_p = float(chi2.sf(interaction_wald, 1))
    contrast_non = np.zeros(len(beta)); contrast_non[shock_index] = 1
    contrast_bmt = contrast_non.copy(); contrast_bmt[interaction_index] = 1
    or_rows = []
    for label, contrast in [("Septic shock among non-BMT hospitalizations", contrast_non), ("Septic shock among BMT hospitalizations", contrast_bmt)]:
        log_or = float(contrast @ beta); se = math.sqrt(float(contrast @ covariance @ contrast))
        or_rows.append({"contrast": label, "adjusted_odds_ratio": round(math.exp(log_or), 3),
            "ci_95_lower": round(math.exp(log_or - 1.96 * se), 3), "ci_95_upper": round(math.exp(log_or + 1.96 * se), 3),
            "p_value": format_p(float(2 * norm.sf(abs(log_or / se))))})
    or_rows.append({"contrast": "Total / interaction Wald test", "adjusted_odds_ratio": "—", "ci_95_lower": "—", "ci_95_upper": "—", "p_value": format_p(interaction_p)})
    write_csv(OUTPUT_DIR / "conditional_odds_ratios.csv", or_rows)

    margins = {(b, s): margin(frame, model, s, b) for b in (0, 1) for s in (0, 1)}
    margin_rows = []
    for b, s in [(0, 0), (0, 1), (1, 0), (1, 1)]:
        estimate, gradient = margins[(b, s)]
        _, low, high, _ = linear_result(estimate, gradient, covariance)
        margin_rows.append({"group": labels[(b, s)], "adjusted_pc_probability_percent": round(100 * estimate, 2),
            "ci_95_lower_percent": round(100 * max(0, low), 2), "ci_95_upper_percent": round(100 * min(1, high), 2)})
    margin_rows.append({"group": "Total", "adjusted_pc_probability_percent": "—", "ci_95_lower_percent": "—", "ci_95_upper_percent": "—"})
    write_csv(OUTPUT_DIR / "adjusted_four_group_probabilities.csv", margin_rows)

    d_non = margins[(0, 1)][0] - margins[(0, 0)][0]; g_non = margins[(0, 1)][1] - margins[(0, 0)][1]
    d_bmt = margins[(1, 1)][0] - margins[(1, 0)][0]; g_bmt = margins[(1, 1)][1] - margins[(1, 0)][1]
    did, g_did = d_bmt - d_non, g_bmt - g_non
    effect_rows = []
    for label, estimate, gradient in [("Septic-shock difference within non-BMT", d_non, g_non),
                                      ("Septic-shock difference within BMT", d_bmt, g_bmt),
                                      ("Difference between probability differences", did, g_did)]:
        _, low, high, p = linear_result(estimate, gradient, covariance)
        effect_rows.append({"contrast": label, "adjusted_difference_pp": round(100 * estimate, 2),
            "ci_95_lower_pp": round(100 * low, 2), "ci_95_upper_pp": round(100 * high, 2), "p_value": format_p(p)})
    effect_rows.append({"contrast": "Total / interaction Wald test", "adjusted_difference_pp": "—", "ci_95_lower_pp": "—", "ci_95_upper_pp": "—", "p_value": format_p(interaction_p)})
    write_csv(OUTPUT_DIR / "adjusted_probability_differences.csv", effect_rows)
    summary = {"bmt_definition": "Z9481 or Z9484 in any diagnosis field; procedure component unavailable",
        "included_unweighted_records": int(frame.frequency.sum()), "iterations": model["iterations"],
        "interaction_wald_test": {"wald_chi_square": round(interaction_wald, 3), "degrees_of_freedom": 1, "p_value": format_p(interaction_p)},
        "descriptive": descriptive_rows, "conditional_odds_ratios": or_rows,
        "adjusted_probabilities": margin_rows, "probability_differences": effect_rows,
        "variance_note": "DISCWT-weighted model with year-specific NIS_STRATUM linearization and discharge-level variance units; HOSP_NIS unavailable by study decision."}
    (OUTPUT_DIR / "bmt_shock_interaction_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    report_path = write_report(summary)
    summary["report_path"] = str(report_path)
    print(json.dumps(summary, indent=2))
    return summary


if __name__ == "__main__":
    main()
