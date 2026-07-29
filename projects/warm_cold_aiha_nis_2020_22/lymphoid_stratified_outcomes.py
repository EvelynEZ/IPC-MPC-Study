"""Unweighted outcomes by AIHA subtype and associated lymphoid malignancy."""

from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import duckdb
import numpy as np
from scipy.stats import chi2_contingency
from statsmodels.stats.oneway import anova_oneway


PROJECT_DIR = Path(__file__).resolve().parent
REPO_ROOT = PROJECT_DIR.parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from src.phase_1_2 import code_match_sql
from projects.warm_cold_aiha_nis_2020_22.thrombosis import FAMILIES, code_condition


OUTPUT_DIR = PROJECT_DIR / "outputs"
DATABASE = OUTPUT_DIR / "warm_cold_aiha_cohort.duckdb"
HM_CONFIG = REPO_ROOT / "config/hm_phenotype_v0_1.json"

GROUP_ORDER = [
    "Warm AIHA without lymphoid malignancy",
    "Warm AIHA with lymphoid malignancy",
    "Cold AIHA without lymphoid malignancy",
    "Cold AIHA with lymphoid malignancy",
]


def format_p(value: float) -> str:
    return "<0.001" if value < 0.001 else f"{value:.3f}"


def diagnosis_prefix(prefix: str) -> str:
    return f"list_contains(list_transform(diagnosis_codes, code -> starts_with(code, '{prefix}')), TRUE)"


def main() -> dict:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    config = json.loads(HM_CONFIG.read_text())
    rules = {rule["id"]: rule for rule in config["subtypes"]}
    lymphoid_parts = []
    for subtype in ("lymphoma", "cll_chronic_leukemia", "myeloma_plasma_cell"):
        condition = code_match_sql("code", rules[subtype])
        lymphoid_parts.append(f"list_contains(list_transform(diagnosis_codes, code -> {condition}), TRUE)")
    lymphoid = " OR ".join(f"({part})" for part in lymphoid_parts)
    family_expressions = {
        label: f"list_contains(list_transform(diagnosis_codes, code -> {code_condition(rule)}), TRUE)"
        for label, rule in FAMILIES.items()
    }
    any_vte = " OR ".join(f"({family_expressions[label]})" for label in list(FAMILIES)[:4])
    any_arterial = " OR ".join(f"({family_expressions[label]})" for label in list(FAMILIES)[4:])
    any_thrombosis = f"({any_vte}) OR ({any_arterial})"
    connection = duckdb.connect(str(DATABASE), read_only=True)
    frame = connection.execute(f"""
        SELECT CASE
                 WHEN aiha_type = 'Warm AIHA' AND NOT ({lymphoid}) THEN '{GROUP_ORDER[0]}'
                 WHEN aiha_type = 'Warm AIHA' AND ({lymphoid}) THEN '{GROUP_ORDER[1]}'
                 WHEN aiha_type = 'Cold AIHA' AND NOT ({lymphoid}) THEN '{GROUP_ORDER[2]}'
                 ELSE '{GROUP_ORDER[3]}'
               END AS analysis_group,
               ({any_thrombosis})::INTEGER AS any_acute_thrombosis,
               ({diagnosis_prefix('N17')})::INTEGER AS acute_kidney_injury,
               ({diagnosis_prefix('A41')})::INTEGER AS sepsis,
               ({diagnosis_prefix('J960')})::INTEGER AS acute_respiratory_failure,
               LOS::DOUBLE AS los, DIED::INTEGER AS mortality
        FROM aiha_cohort
    """).fetchdf()
    connection.close()
    denominators = {group: int((frame.analysis_group == group).sum()) for group in GROUP_ORDER}
    rows = []
    binary_outcomes = [
        ("Any acute thrombosis", "any_acute_thrombosis"),
        ("Acute kidney injury", "acute_kidney_injury"),
        ("Sepsis", "sepsis"),
        ("Acute respiratory failure", "acute_respiratory_failure"),
        ("In-hospital mortality", "mortality"),
    ]
    for label, column in binary_outcomes:
        counts = [int(frame.loc[frame.analysis_group == group, column].fillna(0).sum()) for group in GROUP_ORDER]
        available = [int(frame.loc[frame.analysis_group == group, column].notna().sum()) for group in GROUP_ORDER]
        contingency = np.asarray([[event, total - event] for event, total in zip(counts, available)])
        p_value = float(chi2_contingency(contingency, correction=False)[1])
        total_events, total_available = int(frame[column].fillna(0).sum()), int(frame[column].notna().sum())
        row = {"outcome": label, "summary_type": "n (%)"}
        for group, event, total in zip(GROUP_ORDER, counts, available):
            row[group] = f"{event:,} ({100 * event / total:.2f}%)"
        row["Total cohort"] = f"{total_events:,} ({100 * total_events / total_available:.2f}%)"
        row["overall_p_value"] = format_p(p_value)
        row["test"] = "Pearson chi-square"
        rows.append(row)
    los_samples = [frame.loc[frame.analysis_group == group, "los"].dropna().to_numpy() for group in GROUP_ORDER]
    los_test = anova_oneway(los_samples, use_var="unequal")
    los_row = {"outcome": "Length of stay, days", "summary_type": "Mean (SD)"}
    for group, values in zip(GROUP_ORDER, los_samples):
        los_row[group] = f"{values.mean():.2f} ({values.std(ddof=1):.2f})"
    all_los = frame.los.dropna().to_numpy()
    los_row["Total cohort"] = f"{all_los.mean():.2f} ({all_los.std(ddof=1):.2f})"
    los_row["overall_p_value"] = format_p(float(los_test.pvalue))
    los_row["test"] = "Welch ANOVA"
    # Keep the requested order with LOS before mortality.
    mortality_row = rows.pop()
    rows.extend([los_row, mortality_row])
    total_row = {"outcome": "Total cohort denominator", "summary_type": "n (100%)"}
    for group in GROUP_ORDER:
        total_row[group] = f"{denominators[group]:,} (100.00%)"
    total_row["Total cohort"] = f"{len(frame):,} (100.00%)"
    total_row["overall_p_value"] = "—"; total_row["test"] = "—"
    rows.append(total_row)
    csv_path = OUTPUT_DIR / "unweighted_outcomes_by_aiha_and_lymphoid_status.csv"
    with csv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)
    report = [
        "# Unweighted Outcomes by AIHA Subtype and Lymphoid-Malignancy Status", "",
        "This table uses sampled NIS hospitalization records without `DISCWT`. Associated lymphoid malignancy is the previously defined union of lymphoma, CLL/chronic leukemia, and plasma-cell neoplasms identified across all 40 diagnosis positions. Warm and cold AIHA remain mutually exclusive.", "",
        "Binary outcomes are unweighted n (%); LOS is unweighted mean (SD). Overall p-values compare all four groups using Pearson chi-square tests or Welch's unequal-variance ANOVA.", "",
        "| Outcome | Summary | Warm without lymphoid malignancy | Warm with lymphoid malignancy | Cold without lymphoid malignancy | Cold with lymphoid malignancy | Total cohort | Overall p-value | Test |",
        "| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |",
    ]
    report.extend(
        f'| {row["outcome"]} | {row["summary_type"]} | {row[GROUP_ORDER[0]]} | {row[GROUP_ORDER[1]]} | '
        f'{row[GROUP_ORDER[2]]} | {row[GROUP_ORDER[3]]} | {row["Total cohort"]} | {row["overall_p_value"]} | {row["test"]} |'
        for row in rows
    )
    report.extend(["", "Acute thrombosis is `ANY_ACUTE_THROMBOSIS`, the union of the previously defined acute venous and arterial thrombosis phenotypes. Acute respiratory failure uses `J96.0*` after code normalization. Diagnoses are co-documented and do not establish temporal ordering.", ""])
    report_path = OUTPUT_DIR / "unweighted_outcomes_by_aiha_and_lymphoid_status_report.md"
    report_path.write_text("\n".join(report), encoding="utf-8")
    main_report_path = PROJECT_DIR / "report.md"
    marker = "\n<!-- LYMPHOID_STRATIFIED_OUTCOMES -->\n"
    main_report = main_report_path.read_text(encoding="utf-8").split(marker)[0].rstrip()
    main_report_path.write_text(main_report + marker + "\n".join(report), encoding="utf-8")
    summary = {"group_denominators": denominators, "rows": rows, "csv": str(csv_path), "report": str(report_path)}
    (OUTPUT_DIR / "unweighted_outcomes_by_aiha_and_lymphoid_status_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    return summary


if __name__ == "__main__":
    main()
