"""Adjusted mortality model comparing cold with warm AIHA."""

from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
import statsmodels.api as sm


PROJECT_DIR = Path(__file__).resolve().parent
REPO_ROOT = PROJECT_DIR.parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from projects.warm_cold_aiha_nis_2020_22.baseline import build_frame


OUTPUT_DIR = PROJECT_DIR / "outputs"


FORMULA = (
    "died ~ C(aiha_type, Treatment(reference='Warm AIHA')) + age + "
    "C(sex, Treatment(reference='Male')) + "
    "C(race, Treatment(reference='White')) + associated_lymphoid_malignancy + "
    "C(region, Treatment(reference='Northeast')) + "
    "C(teaching_status, Treatment(reference='Nonteaching')) + "
    "C(bed_size, Treatment(reference='Small')) + cci"
)


def format_p(value: float) -> str:
    return "<0.001" if value < 0.001 else f"{value:.3f}"


def readable_term(term: str) -> str:
    if term == "Intercept":
        return "Intercept"
    if term == "age":
        return "Age, per year"
    if term == "cci":
        return "Charlson Comorbidity Index, per point"
    if term == "associated_lymphoid_malignancy":
        return "Associated lymphoid malignancy: yes vs no"
    if term.startswith("C(aiha_type"):
        return "Cold AIHA vs warm AIHA"
    for prefix, label in [("C(sex", "Sex"), ("C(race", "Race/ethnicity"),
                          ("C(region", "Hospital region"),
                          ("C(teaching_status", "Hospital teaching status"),
                          ("C(bed_size", "Hospital bed size")]:
        if term.startswith(prefix):
            level = term.split("[T.", 1)[1].rstrip("]")
            return f"{label}: {level}"
    return term


def main() -> dict:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    frame = build_frame().copy()
    frame = frame[frame.died.notna()].copy()
    frame["died"] = frame["died"].astype(int)
    frame["teaching_status"] = np.where(
        frame.location_teaching.eq("Urban teaching"), "Teaching", "Nonteaching"
    )
    # DISCWT equals 5 for every analytic record, so weighting changes neither
    # coefficients nor fitted probabilities. Unit weights avoid treating each
    # sampled discharge as five independent observations for model-based SEs.
    model = smf.glm(FORMULA, data=frame, family=sm.families.Binomial()).fit()
    confidence = model.conf_int()
    rows = []
    for term in model.params.index:
        rows.append({
            "variable": readable_term(term),
            "adjusted_odds_ratio": round(float(np.exp(model.params[term])), 3),
            "ci_lower": round(float(np.exp(confidence.loc[term, 0])), 3),
            "ci_upper": round(float(np.exp(confidence.loc[term, 1])), 3),
            "p_value": format_p(float(model.pvalues[term])),
        })
    exposure = next(row for row in rows if row["variable"] == "Cold AIHA vs warm AIHA")
    diagnostics = {
        "analysis_n": int(model.nobs),
        "deaths": int(frame.died.sum()),
        "warm_n": int(frame.aiha_type.eq("Warm AIHA").sum()),
        "warm_deaths": int(frame.loc[frame.aiha_type.eq("Warm AIHA"), "died"].sum()),
        "cold_n": int(frame.aiha_type.eq("Cold AIHA").sum()),
        "cold_deaths": int(frame.loc[frame.aiha_type.eq("Cold AIHA"), "died"].sum()),
        "parameters": int(len(model.params)),
        "converged": bool(model.converged),
        "aic": round(float(model.aic), 2),
    }
    csv_path = OUTPUT_DIR / "adjusted_mortality_logistic_regression.csv"
    with csv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)
    report = [
        "# Adjusted Association Between AIHA Subtype and In-Hospital Mortality", "",
        "A multivariable logistic-regression model used in-hospital mortality (`DIED=1`) as the outcome and AIHA subtype as the primary exposure. Warm AIHA was the reference group. The model adjusted for continuous age, sex, race/ethnicity, associated lymphoid malignancy, hospital region, hospital teaching status, hospital bed size, and continuous Charlson Comorbidity Index. Teaching status was coded as urban teaching versus nonteaching; the nonteaching category combines rural and urban nonteaching hospitals.", "",
        "## Primary exposure result", "",
        "| Exposure comparison | Adjusted odds ratio | 95% CI | P-value |",
        "| --- | ---: | ---: | ---: |",
        f'| Cold AIHA versus warm AIHA | {exposure["adjusted_odds_ratio"]:.3f} | {exposure["ci_lower"]:.3f}–{exposure["ci_upper"]:.3f} | {exposure["p_value"]} |', "",
        "## Model diagnostics", "",
        "| Measure | Value |", "| --- | ---: |",
        f'| Analysis hospitalizations | {diagnostics["analysis_n"]:,} |',
        f'| Deaths | {diagnostics["deaths"]:,} |',
        f'| Warm AIHA deaths / hospitalizations | {diagnostics["warm_deaths"]:,} / {diagnostics["warm_n"]:,} |',
        f'| Cold AIHA deaths / hospitalizations | {diagnostics["cold_deaths"]:,} / {diagnostics["cold_n"]:,} |',
        f'| Estimated model parameters | {diagnostics["parameters"]} |',
        f'| Model converged | {diagnostics["converged"]} |', "",
        "`DISCWT` is 5 for every record in this 2020–2022 analytic cohort. Therefore weighted and unweighted coefficient estimates are identical; the model uses each sampled hospitalization once so that model-based standard errors are not artificially reduced by treating the weight as replicated observations.", "",
        "Interpret cautiously: only 39 deaths were observed relative to the number of adjustment parameters, so the estimate may be imprecise and the model is vulnerable to sparse-data overfitting.", "",
        "## Full adjusted model", "",
        "| Covariate contrast | Adjusted odds ratio | 95% CI | P-value |",
        "| --- | ---: | ---: | ---: |",
    ]
    report.extend(f'| {row["variable"]} | {row["adjusted_odds_ratio"]:.3f} | {row["ci_lower"]:.3f}–{row["ci_upper"]:.3f} | {row["p_value"]} |' for row in rows if row["variable"] != "Intercept")
    report.append("")
    report_path = OUTPUT_DIR / "adjusted_mortality_logistic_regression_report.md"
    report_path.write_text("\n".join(report), encoding="utf-8")
    main_report_path = PROJECT_DIR / "report.md"
    marker = "\n<!-- ADJUSTED_MORTALITY_MODEL -->\n"
    main_report = main_report_path.read_text(encoding="utf-8").split(marker)[0].rstrip()
    main_report_path.write_text(main_report + marker + "\n".join(report), encoding="utf-8")
    summary = {"primary_exposure": exposure, "diagnostics": diagnostics, "all_coefficients": rows,
               "formula": FORMULA, "csv": str(csv_path), "report": str(report_path)}
    (OUTPUT_DIR / "adjusted_mortality_logistic_regression_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    return summary


if __name__ == "__main__":
    main()
