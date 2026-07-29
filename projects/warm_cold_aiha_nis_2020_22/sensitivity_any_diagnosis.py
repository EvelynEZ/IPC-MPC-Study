"""Any-diagnosis-position AIHA sensitivity cohort and adjusted mortality model."""

from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import numpy as np
import statsmodels.api as sm
import statsmodels.formula.api as smf
import duckdb


PROJECT_DIR = Path(__file__).resolve().parent
REPO_ROOT = PROJECT_DIR.parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from projects.warm_cold_aiha_nis_2020_22.analysis import dataset_glob, normalize
from projects.warm_cold_aiha_nis_2020_22.baseline import build_frame
from projects.warm_cold_aiha_nis_2020_22.mortality_model import FORMULA, format_p


OUTPUT_DIR = PROJECT_DIR / "outputs"
DATABASE = OUTPUT_DIR / "warm_cold_aiha_cohort.duckdb"
TABLE = "aiha_sensitivity_cohort"


def main() -> dict:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    all_dx = [normalize(f"I10_DX{i}") for i in range(1, 41)]
    warm = " OR ".join(f"coalesce({code} = 'D5911', FALSE)" for code in all_dx)
    cold = " OR ".join(f"coalesce({code} = 'D5912', FALSE)" for code in all_dx)
    diagnosis_list = ", ".join(all_dx)
    connection = duckdb.connect(str(DATABASE))
    connection.execute(f"""
        CREATE TABLE IF NOT EXISTS {TABLE} AS
        SELECT *, [{diagnosis_list}] AS diagnosis_codes,
               CASE WHEN ({warm}) THEN 'Warm AIHA' ELSE 'Cold AIHA' END AS aiha_type
        FROM read_parquet(?)
        WHERE AGE >= 18 AND YEAR BETWEEN 2020 AND 2022
          AND (({warm}) OR ({cold}))
          AND NOT (({warm}) AND ({cold}))
    """, [dataset_glob()])
    cohort_counts = connection.execute(f"""
        SELECT aiha_type, count(*)::BIGINT, round(sum(DISCWT), 0)::BIGINT,
               count(*) FILTER (WHERE DIED = 1)::BIGINT
        FROM {TABLE} GROUP BY 1 ORDER BY 1
    """).fetchall()
    connection.close()
    counts = {row[0]: {"unweighted_n": int(row[1]), "weighted_n": int(row[2]), "deaths": int(row[3])}
              for row in cohort_counts}
    frame = build_frame(TABLE)
    frame = frame[frame.died.notna()].copy()
    frame["died"] = frame.died.astype(int)
    frame["teaching_status"] = np.where(
        frame.location_teaching.eq("Urban teaching"), "Teaching", "Nonteaching"
    )
    model = smf.glm(FORMULA, data=frame, family=sm.families.Binomial()).fit()
    term = next(term for term in model.params.index if term.startswith("C(aiha_type"))
    confidence = model.conf_int().loc[term]
    estimate = {
        "exposure": "Cold AIHA vs warm AIHA",
        "adjusted_odds_ratio": round(float(np.exp(model.params[term])), 3),
        "ci_lower": round(float(np.exp(confidence.iloc[0])), 3),
        "ci_upper": round(float(np.exp(confidence.iloc[1])), 3),
        "p_value": format_p(float(model.pvalues[term])),
    }
    total_n = sum(group["unweighted_n"] for group in counts.values())
    total_weighted = sum(group["weighted_n"] for group in counts.values())
    total_deaths = sum(group["deaths"] for group in counts.values())
    rows = [
        {"cohort": "Warm AIHA", **counts["Warm AIHA"]},
        {"cohort": "Cold AIHA", **counts["Cold AIHA"]},
        {"cohort": "Total sensitivity cohort", "unweighted_n": total_n,
         "weighted_n": total_weighted, "deaths": total_deaths},
    ]
    csv_path = OUTPUT_DIR / "any_diagnosis_sensitivity_cohort_counts.csv"
    with csv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)
    report = [
        "# Sensitivity Analysis: AIHA in Any Diagnosis Position", "",
        "The sensitivity cohort includes adult NIS hospitalizations from 2020–2022 with normalized `D59.11` or `D59.12` in any of DX1–DX40. Hospitalizations containing both warm and cold codes anywhere were excluded. This cohort is stored separately and does not replace the primary DX1–DX3 cohort.", "",
        "## Sensitivity cohort", "",
        "| Cohort | Unweighted hospitalizations | DISCWT-weighted hospitalizations | Unweighted deaths |",
        "| --- | ---: | ---: | ---: |",
    ]
    report.extend(f'| {row["cohort"]} | {row["unweighted_n"]:,} | {row["weighted_n"]:,} | {row["deaths"]:,} |' for row in rows)
    report.extend([
        "", "## Repeated primary adjusted mortality analysis", "",
        "The logistic-regression specification is identical to the primary analysis: mortality is the outcome; warm AIHA is the exposure reference; and adjustment includes continuous age, sex, race/ethnicity, associated lymphoid malignancy, hospital region, teaching status, bed size, and continuous CCI.", "",
        "| Exposure comparison | Adjusted odds ratio | 95% CI | P-value |",
        "| --- | ---: | ---: | ---: |",
        f'| Cold AIHA versus warm AIHA | {estimate["adjusted_odds_ratio"]:.3f} | {estimate["ci_lower"]:.3f}–{estimate["ci_upper"]:.3f} | {estimate["p_value"]} |', "",
        "As in the primary cohort, `DISCWT` is constant at 5, so coefficient estimates are unchanged by weighting. Each sampled discharge is used once for model-based variance estimation.", "",
    ])
    report_path = OUTPUT_DIR / "any_diagnosis_sensitivity_mortality_report.md"
    report_path.write_text("\n".join(report), encoding="utf-8")
    main_report_path = PROJECT_DIR / "report.md"
    marker = "\n<!-- ANY_DIAGNOSIS_SENSITIVITY -->\n"
    main_report = main_report_path.read_text(encoding="utf-8").split(marker)[0].rstrip()
    main_report_path.write_text(main_report + marker + "\n".join(report), encoding="utf-8")
    summary = {"definition": "D59.11 or D59.12 in DX1-DX40; exclude overlap anywhere",
               "counts": rows, "adjusted_mortality": estimate, "analysis_n": int(model.nobs),
               "parameters": int(len(model.params)), "converged": bool(model.converged),
               "csv": str(csv_path), "report": str(report_path)}
    (OUTPUT_DIR / "any_diagnosis_sensitivity_mortality_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    return summary


if __name__ == "__main__":
    main()
