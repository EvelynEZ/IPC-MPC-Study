"""One consolidated unweighted comparison of all calculated AIHA variables."""

from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import duckdb
import numpy as np
import pandas as pd
from scipy.stats import chi2_contingency, ttest_ind

PROJECT_DIR = Path(__file__).resolve().parent
REPO_ROOT = PROJECT_DIR.parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from projects.warm_cold_aiha_nis_2020_22.baseline import build_frame


OUTPUT_DIR = PROJECT_DIR / "outputs"
DATABASE = OUTPUT_DIR / "warm_cold_aiha_cohort.duckdb"


def format_p(value: float) -> str:
    return "<0.001" if value < 0.001 else f"{value:.3f}"


def binary_p(warm_event: int, cold_event: int, warm_n: int, cold_n: int) -> str:
    table = [[warm_event, warm_n - warm_event], [cold_event, cold_n - cold_event]]
    return format_p(float(chi2_contingency(table, correction=False)[1]))


def main() -> dict:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    frame = build_frame()
    connection = duckdb.connect(str(DATABASE), read_only=True)
    years = connection.execute("SELECT aiha_type, YEAR::INTEGER AS year FROM aiha_cohort").fetchdf()
    connection.close()
    frame = pd.concat([frame.reset_index(drop=True), years[["year"]].reset_index(drop=True)], axis=1)
    warm = frame[frame.aiha_type.eq("Warm AIHA")]
    cold = frame[frame.aiha_type.eq("Cold AIHA")]
    warm_n, cold_n = len(warm), len(cold)
    rows: list[dict] = []

    def add_continuous(section: str, characteristic: str, column: str) -> None:
        wx, cx = warm[column].dropna(), cold[column].dropna()
        p = float(ttest_ind(cx, wx, equal_var=False).pvalue)
        rows.append({"section": section, "characteristic": characteristic, "level": "Mean (SD)",
                     "warm_n_percent_or_summary": f"{wx.mean():.2f} ({wx.std(ddof=1):.2f})",
                     "cold_n_percent_or_summary": f"{cx.mean():.2f} ({cx.std(ddof=1):.2f})",
                     "difference_percentage_points": "—", "p_value": format_p(p), "test": "Welch t-test"})

    def add_categorical(section: str, characteristic: str, column: str, levels: list) -> None:
        nonmissing = [level for level in levels if level not in ("Missing", "Unknown")]
        table = np.asarray([[(warm[column] == level).sum(), (cold[column] == level).sum()] for level in nonmissing]).T
        table = table[:, table.sum(axis=0) > 0]
        p = float(chi2_contingency(table)[1]) if table.shape[1] > 1 else 1.0
        for index, level in enumerate(levels):
            wn, cn = int((warm[column] == level).sum()), int((cold[column] == level).sum())
            rows.append({"section": section, "characteristic": characteristic if index == 0 else "", "level": str(level),
                         "warm_n_percent_or_summary": f"{wn:,} ({100 * wn / warm_n:.2f}%)",
                         "cold_n_percent_or_summary": f"{cn:,} ({100 * cn / cold_n:.2f}%)",
                         "difference_percentage_points": f"{100 * (cn / cold_n - wn / warm_n):.2f}",
                         "p_value": format_p(p) if index == 0 else "", "test": "Pearson chi-square" if index == 0 else ""})

    add_continuous("Demographic and clinical characteristics", "Age, years", "age")
    add_categorical("Demographic and clinical characteristics", "Age group", "age_group", ["18–59", "≥60"])
    add_categorical("Demographic and clinical characteristics", "Sex", "sex", ["Male", "Female", "Missing"])
    add_categorical("Demographic and clinical characteristics", "Race/ethnicity", "race", ["White", "Black", "Hispanic", "Asian/Pacific Islander", "Native American", "Other", "Missing"])
    add_categorical("Demographic and clinical characteristics", "Primary payer", "payer", ["Medicare", "Medicaid", "Private insurance", "Self-pay", "No charge", "Other", "Missing"])
    add_categorical("Demographic and clinical characteristics", "ZIP-code income quartile", "income", ["0–25th percentile", "26th–50th percentile", "51st–75th percentile", "76th–100th percentile", "Missing"])
    add_categorical("Hospital characteristics", "Hospital region", "region", ["Northeast", "Midwest", "South", "West", "Unknown"])
    add_categorical("Hospital characteristics", "Hospital bed size", "bed_size", ["Small", "Medium", "Large", "Unknown"])
    add_categorical("Hospital characteristics", "Hospital location/teaching status", "location_teaching", ["Rural", "Urban nonteaching", "Urban teaching", "Unknown"])
    add_categorical("Admission year", "Admission year", "year", [2020, 2021, 2022])
    add_continuous("Clinical characteristics and outcomes", "Charlson Comorbidity Index", "cci")
    add_categorical("Clinical characteristics and outcomes", "Charlson category", "cci_category", ["0", "1–2", "≥3"])
    add_continuous("Clinical characteristics and outcomes", "Length of stay, days", "los")
    add_categorical("Clinical characteristics and outcomes", "In-hospital mortality", "mortality", ["Survived", "Died", "Missing"])

    phenotype_sources = [
        ("Lymphoid malignancy", "any_lymphoid_malignancy_by_aiha_type.csv", "malignancy_definition"),
        ("Acute thrombosis", "acute_thrombosis_by_aiha_type.csv", "diagnosis_family"),
        ("Selected complications", "selected_complications_by_aiha_type.csv", "complication"),
    ]
    for section, filename, label_column in phenotype_sources:
        source = pd.read_csv(OUTPUT_DIR / filename)
        source = source[~source[label_column].eq("Total cohort denominator")]
        for _, item in source.iterrows():
            wn, cn = int(item.warm_unweighted_n), int(item.cold_unweighted_n)
            rows.append({"section": section, "characteristic": str(item[label_column]), "level": "Present",
                         "warm_n_percent_or_summary": f"{wn:,} ({100 * wn / warm_n:.2f}%)",
                         "cold_n_percent_or_summary": f"{cn:,} ({100 * cn / cold_n:.2f}%)",
                         "difference_percentage_points": f"{100 * (cn / cold_n - wn / warm_n):.2f}",
                         "p_value": binary_p(wn, cn, warm_n, cold_n), "test": "Pearson chi-square"})

    rows.append({"section": "Total", "characteristic": "Total cohort denominator", "level": "All hospitalizations",
                 "warm_n_percent_or_summary": f"{warm_n:,} (100.00%)", "cold_n_percent_or_summary": f"{cold_n:,} (100.00%)",
                 "difference_percentage_points": "—", "p_value": "—", "test": "—"})
    csv_path = OUTPUT_DIR / "all_variables_unweighted_warm_vs_cold.csv"
    with csv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)

    report = ["# Unweighted Comparison of Warm Versus Cold AIHA Hospitalizations", "",
              "This table reports the actual sampled NIS hospitalization records without applying `DISCWT`. Warm AIHA is the reference group. Continuous variables are shown as unweighted mean (SD); categorical and binary variables are shown as unweighted n (%).", "",
              "P-values are unweighted Welch t-tests for continuous variables, overall Pearson chi-square tests for multilevel categorical variables, and Pearson chi-square tests for binary diagnoses. A p-value shown on the first level of a multilevel variable applies to the variable overall.", ""]
    current_section = None
    for row_index, row in enumerate(rows):
        if row["section"] != current_section:
            current_section = row["section"]
            report.extend([f"## {current_section}", "", "| Characteristic | Level | Warm AIHA | Cold AIHA | Difference, pp | P-value | Test |",
                           "| --- | --- | ---: | ---: | ---: | ---: | --- |"])
        report.append(f'| {row["characteristic"]} | {row["level"]} | {row["warm_n_percent_or_summary"]} | {row["cold_n_percent_or_summary"]} | {row["difference_percentage_points"]} | {row["p_value"]} | {row["test"]} |')
        if row_index == len(rows) - 1 or rows[row_index + 1]["section"] != current_section:
            report.append("")
    report.extend(["Lymphoid-malignancy component rows and thrombosis component rows may overlap. Composite rows count each hospitalization once. Diagnoses searched all 40 diagnosis positions.", ""])
    report_path = OUTPUT_DIR / "all_variables_unweighted_warm_vs_cold_report.md"
    report_path.write_text("\n".join(report), encoding="utf-8")
    main_report_path = PROJECT_DIR / "report.md"
    marker = "\n<!-- ALL_VARIABLES_UNWEIGHTED -->\n"
    main_report = main_report_path.read_text(encoding="utf-8").split(marker)[0].rstrip()
    main_report_path.write_text(main_report + marker + "\n".join(report), encoding="utf-8")
    summary = {"reference_group": "Warm AIHA", "warm_n": warm_n, "cold_n": cold_n,
               "rows": rows, "csv": str(csv_path), "report": str(report_path)}
    (OUTPUT_DIR / "all_variables_unweighted_warm_vs_cold_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps({"warm_n": warm_n, "cold_n": cold_n, "rows": len(rows), "report": str(report_path)}, indent=2))
    return summary


if __name__ == "__main__":
    main()
