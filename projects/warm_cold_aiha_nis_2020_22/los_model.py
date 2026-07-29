"""Adjusted length-of-stay model comparing cold with warm AIHA."""

from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import numpy as np
import statsmodels.formula.api as smf


PROJECT_DIR = Path(__file__).resolve().parent
REPO_ROOT = PROJECT_DIR.parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from projects.warm_cold_aiha_nis_2020_22.baseline import build_frame
from projects.warm_cold_aiha_nis_2020_22.mortality_model import FORMULA as MORTALITY_FORMULA, readable_term


OUTPUT_DIR = PROJECT_DIR / "outputs"
FORMULA = MORTALITY_FORMULA.replace("died ~", "los ~")


def format_p(value: float) -> str:
    return "<0.001" if value < 0.001 else f"{value:.3f}"


def main() -> dict:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    frame = build_frame().copy()
    frame = frame[frame.los.notna()].copy()
    frame["teaching_status"] = np.where(
        frame.location_teaching.eq("Urban teaching"), "Teaching", "Nonteaching"
    )
    model = smf.ols(FORMULA, data=frame).fit(cov_type="HC3")
    confidence = model.conf_int()
    rows = []
    for term in model.params.index:
        rows.append({
            "variable": readable_term(term),
            "adjusted_coefficient": round(float(model.params[term]), 3),
            "ci_lower": round(float(confidence.loc[term, 0]), 3),
            "ci_upper": round(float(confidence.loc[term, 1]), 3),
            "p_value": format_p(float(model.pvalues[term])),
        })
    exposure = next(row for row in rows if row["variable"] == "Cold AIHA vs warm AIHA")
    diagnostics = {
        "analysis_n": int(model.nobs), "parameters": int(len(model.params)),
        "r_squared": round(float(model.rsquared), 4), "robust_covariance": "HC3",
        "warm_unadjusted_mean_los": round(float(frame.loc[frame.aiha_type.eq("Warm AIHA"), "los"].mean()), 2),
        "cold_unadjusted_mean_los": round(float(frame.loc[frame.aiha_type.eq("Cold AIHA"), "los"].mean()), 2),
    }
    csv_path = OUTPUT_DIR / "adjusted_los_linear_regression.csv"
    with csv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)
    report = [
        "# Adjusted Association Between AIHA Subtype and Length of Stay", "",
        "A multivariable linear-regression model used hospital length of stay in days as the continuous outcome and AIHA subtype as the primary exposure. Warm AIHA was the reference group. The model adjusted for continuous age, sex, race/ethnicity, associated lymphoid malignancy, hospital region, hospital teaching status, hospital bed size, and continuous Charlson Comorbidity Index. Teaching status was coded as urban teaching versus nonteaching; nonteaching combines rural and urban nonteaching hospitals.", "",
        "Heteroskedasticity-consistent HC3 standard errors were used because length of stay is right-skewed. The exposure coefficient is an adjusted mean difference in days, calculated as cold AIHA minus warm AIHA.", "",
        "## Primary exposure result", "",
        "| Exposure comparison | Adjusted mean difference, days | 95% CI | P-value |",
        "| --- | ---: | ---: | ---: |",
        f'| Cold AIHA versus warm AIHA | {exposure["adjusted_coefficient"]:.3f} | {exposure["ci_lower"]:.3f} to {exposure["ci_upper"]:.3f} | {exposure["p_value"]} |', "",
        "## Model diagnostics", "",
        "| Measure | Value |", "| --- | ---: |",
        f'| Analysis hospitalizations | {diagnostics["analysis_n"]:,} |',
        f'| Unadjusted mean LOS, warm AIHA | {diagnostics["warm_unadjusted_mean_los"]:.2f} days |',
        f'| Unadjusted mean LOS, cold AIHA | {diagnostics["cold_unadjusted_mean_los"]:.2f} days |',
        f'| R-squared | {diagnostics["r_squared"]:.4f} |',
        f'| Standard-error estimator | {diagnostics["robust_covariance"]} |', "",
        "`DISCWT` equals 5 for every record in this 2020–2022 cohort, so applying it does not change coefficients or fitted values. Each sampled hospitalization was used once for variance estimation.", "",
        "## Full adjusted model", "",
        "| Covariate contrast | Adjusted coefficient, days | 95% CI | P-value |",
        "| --- | ---: | ---: | ---: |",
    ]
    report.extend(f'| {row["variable"]} | {row["adjusted_coefficient"]:.3f} | {row["ci_lower"]:.3f} to {row["ci_upper"]:.3f} | {row["p_value"]} |' for row in rows if row["variable"] != "Intercept")
    report.append("")
    report_path = OUTPUT_DIR / "adjusted_los_linear_regression_report.md"
    report_path.write_text("\n".join(report), encoding="utf-8")
    main_report_path = PROJECT_DIR / "report.md"
    marker = "\n<!-- ADJUSTED_LOS_MODEL -->\n"
    main_report = main_report_path.read_text(encoding="utf-8").split(marker)[0].rstrip()
    main_report_path.write_text(main_report + marker + "\n".join(report), encoding="utf-8")
    summary = {"primary_exposure": exposure, "diagnostics": diagnostics, "all_coefficients": rows,
               "formula": FORMULA, "csv": str(csv_path), "report": str(report_path)}
    (OUTPUT_DIR / "adjusted_los_linear_regression_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    return summary


if __name__ == "__main__":
    main()
