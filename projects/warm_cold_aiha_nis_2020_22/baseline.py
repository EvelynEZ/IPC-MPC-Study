"""Warm-versus-cold AIHA baseline characteristics."""

from __future__ import annotations

import csv
import json
import math
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

from src.phase_5_cci_los_mortality import CCI_COMPONENTS, CCI_WEIGHTS, component_condition
from src.phase_1_2 import code_match_sql


OUTPUT_DIR = PROJECT_DIR / "outputs"
DATABASE = OUTPUT_DIR / "warm_cold_aiha_cohort.duckdb"
HM_CONFIG = REPO_ROOT / "config/hm_phenotype_v0_1.json"

CANCER_PREFIXES = ("C00", "C01", "C02", "C03", "C04", "C05", "C06", "C07", "C08", "C09", "C10", "C11", "C12", "C13", "C14", "C15", "C16", "C17", "C18", "C19", "C20", "C21", "C22", "C23", "C24", "C25", "C26", "C30", "C31", "C32", "C33", "C34", "C37", "C38", "C39", "C40", "C41", "C43", "C45", "C46", "C47", "C48", "C49", "C50", "C51", "C52", "C53", "C54", "C55", "C56", "C57", "C58", "C60", "C61", "C62", "C63", "C64", "C65", "C66", "C67", "C68", "C69", "C70", "C71", "C72", "C73", "C74", "C75", "C76", "C81", "C82", "C83", "C84", "C85", "C88", "C90", "C91", "C92", "C93", "C94", "C95", "C96", "C97")
METASTATIC_PREFIXES = ("C77", "C78", "C79", "C80")


def format_p(value: float) -> str:
    return "<0.001" if value < 0.001 else f"{value:.3f}"


def weighted_mean_sd(values: pd.Series, weights: pd.Series) -> tuple[float, float]:
    valid = values.notna() & weights.notna()
    x, w = values[valid].astype(float).to_numpy(), weights[valid].astype(float).to_numpy()
    mean = float(np.average(x, weights=w))
    sd = math.sqrt(max(0.0, float(np.average((x - mean) ** 2, weights=w))))
    return mean, sd


def proportion_smd(cold: float, warm: float) -> float:
    variance = (cold * (1 - cold) + warm * (1 - warm)) / 2
    return 0.0 if variance <= 0 else (cold - warm) / math.sqrt(variance)


def build_frame() -> pd.DataFrame:
    connection = duckdb.connect(str(DATABASE), read_only=True)
    subtype_rules = {rule["id"]: rule for rule in json.loads(HM_CONFIG.read_text())["subtypes"]}
    lymphoid_conditions = [code_match_sql("code", subtype_rules[subtype]) for subtype in
                           ("lymphoma", "cll_chronic_leukemia", "myeloma_plasma_cell")]
    lymphoid_expression = " OR ".join(
        f"list_contains(list_transform(diagnosis_codes, code -> {condition}), TRUE)"
        for condition in lymphoid_conditions
    )
    components = dict(CCI_COMPONENTS)
    components.update({"cancer": CANCER_PREFIXES, "metastatic": METASTATIC_PREFIXES})
    weights = dict(CCI_WEIGHTS); weights.update({"cancer": 2, "metastatic": 6})
    flags = [f"list_contains(list_transform(diagnosis_codes, code -> {component_condition(prefixes)}), TRUE) AS {name}"
             for name, prefixes in components.items()]
    connection.execute("CREATE TEMP TABLE flags AS SELECT *, " + ", ".join(flags) + " FROM aiha_cohort")
    terms = []
    for name, weight in weights.items():
        present = "diab AND NOT diabwc" if name == "diab" else "mld AND NOT msld" if name == "mld" else "cancer AND NOT metastatic" if name == "cancer" else name
        terms.append(f"CASE WHEN {present} THEN {weight} ELSE 0 END")
    frame = connection.execute("""
        SELECT aiha_type, DISCWT::DOUBLE AS weight, AGE::DOUBLE AS age, LOS::DOUBLE AS los,
               CASE WHEN (""" + lymphoid_expression + """) THEN 1 ELSE 0 END AS associated_lymphoid_malignancy,
               DIED::INTEGER AS died,
               CASE WHEN AGE < 60 THEN '18–59' ELSE '≥60' END AS age_group,
               CASE FEMALE WHEN 0 THEN 'Male' WHEN 1 THEN 'Female' ELSE 'Missing' END AS sex,
               CASE RACE WHEN 1 THEN 'White' WHEN 2 THEN 'Black' WHEN 3 THEN 'Hispanic'
                   WHEN 4 THEN 'Asian/Pacific Islander' WHEN 5 THEN 'Native American' WHEN 6 THEN 'Other' ELSE 'Missing' END AS race,
               CASE PAY1 WHEN 1 THEN 'Medicare' WHEN 2 THEN 'Medicaid' WHEN 3 THEN 'Private insurance'
                   WHEN 4 THEN 'Self-pay' WHEN 5 THEN 'No charge' WHEN 6 THEN 'Other' ELSE 'Missing' END AS payer,
               CASE ZIPINC_QRTL WHEN 1 THEN '0–25th percentile' WHEN 2 THEN '26th–50th percentile'
                   WHEN 3 THEN '51st–75th percentile' WHEN 4 THEN '76th–100th percentile' ELSE 'Missing' END AS income,
               CASE WHEN floor(NIS_STRATUM / 1000) IN (1,2) THEN 'Northeast'
                   WHEN floor(NIS_STRATUM / 1000) IN (3,4) THEN 'Midwest'
                   WHEN floor(NIS_STRATUM / 1000) IN (5,6,7) THEN 'South'
                   WHEN floor(NIS_STRATUM / 1000) IN (8,9) THEN 'West' ELSE 'Unknown' END AS region,
               CASE CAST(NIS_STRATUM % 10 AS INTEGER) WHEN 1 THEN 'Small' WHEN 2 THEN 'Medium' WHEN 3 THEN 'Large' ELSE 'Unknown' END AS bed_size,
               CASE CAST(floor(NIS_STRATUM / 10) % 10 AS INTEGER) WHEN 1 THEN 'Rural'
                   WHEN 2 THEN 'Urban nonteaching' WHEN 3 THEN 'Urban teaching' ELSE 'Unknown' END AS location_teaching,
               """ + " + ".join(terms) + """ AS cci
        FROM flags
    """).fetchdf()
    connection.close()
    frame["cci_category"] = pd.cut(frame.cci, bins=[-1, 0, 2, np.inf], labels=["0", "1–2", "≥3"]).astype(str)
    frame["mortality"] = frame.died.map({0: "Survived", 1: "Died"}).fillna("Missing")
    return frame


def main() -> dict:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    frame = build_frame()
    warm, cold = frame[frame.aiha_type.eq("Warm AIHA")], frame[frame.aiha_type.eq("Cold AIHA")]
    rows = []
    for label, column in [("Age, years", "age"), ("Charlson Comorbidity Index", "cci"), ("Length of stay, days", "los")]:
        wm, ws = weighted_mean_sd(warm[column], warm.weight); cm, cs = weighted_mean_sd(cold[column], cold.weight)
        p = float(ttest_ind(cold[column].dropna(), warm[column].dropna(), equal_var=False).pvalue)
        pooled = math.sqrt((ws * ws + cs * cs) / 2)
        rows.append({"characteristic": label, "level": "Mean (SD)", "warm_unweighted_n": int(warm[column].notna().sum()),
                     "warm_weighted_summary": f"{wm:.2f} ({ws:.2f})", "cold_unweighted_n": int(cold[column].notna().sum()),
                     "cold_weighted_summary": f"{cm:.2f} ({cs:.2f})", "difference_percentage_points": "—",
                     "standardized_mean_difference": round((cm - wm) / pooled, 3) if pooled else 0,
                     "p_value": format_p(p), "test": "Welch t-test"})
    specifications = [
        ("Age group", "age_group", ["18–59", "≥60"]),
        ("Sex", "sex", ["Male", "Female", "Missing"]),
        ("Race/ethnicity", "race", ["White", "Black", "Hispanic", "Asian/Pacific Islander", "Native American", "Other", "Missing"]),
        ("Primary payer", "payer", ["Medicare", "Medicaid", "Private insurance", "Self-pay", "No charge", "Other", "Missing"]),
        ("ZIP-code income quartile", "income", ["0–25th percentile", "26th–50th percentile", "51st–75th percentile", "76th–100th percentile", "Missing"]),
        ("Hospital region", "region", ["Northeast", "Midwest", "South", "West", "Unknown"]),
        ("Hospital bed size", "bed_size", ["Small", "Medium", "Large", "Unknown"]),
        ("Hospital location/teaching status", "location_teaching", ["Rural", "Urban nonteaching", "Urban teaching", "Unknown"]),
        ("Charlson category", "cci_category", ["0", "1–2", "≥3"]),
        ("In-hospital mortality", "mortality", ["Survived", "Died", "Missing"]),
    ]
    for characteristic, column, levels in specifications:
        valid_levels = [level for level in levels if level != "Missing"]
        contingency = [[int((warm[column] == level).sum()), int((cold[column] == level).sum())] for level in valid_levels]
        contingency = [row for row in contingency if sum(row) > 0]
        p = float(chi2_contingency(np.asarray(contingency).T)[1]) if len(contingency) > 1 else 1.0
        warm_weight_total, cold_weight_total = warm.weight.sum(), cold.weight.sum()
        for index, level in enumerate(levels):
            wn, cn = int((warm[column] == level).sum()), int((cold[column] == level).sum())
            ww = float(warm.loc[warm[column] == level, "weight"].sum()); cw = float(cold.loc[cold[column] == level, "weight"].sum())
            wp, cp = ww / warm_weight_total, cw / cold_weight_total
            rows.append({"characteristic": characteristic if index == 0 else "", "level": level,
                         "warm_unweighted_n": wn, "warm_weighted_summary": f"{round(ww):,} ({100*wp:.2f}%)",
                         "cold_unweighted_n": cn, "cold_weighted_summary": f"{round(cw):,} ({100*cp:.2f}%)",
                         "difference_percentage_points": round(100 * (cp - wp), 2),
                         "standardized_mean_difference": round(proportion_smd(cp, wp), 3),
                         "p_value": format_p(p) if index == 0 else "", "test": "Pearson chi-square" if index == 0 else ""})
    rows.append({"characteristic": "Total", "level": "All final-cohort hospitalizations",
                 "warm_unweighted_n": len(warm), "warm_weighted_summary": f"{round(warm.weight.sum()):,}",
                 "cold_unweighted_n": len(cold), "cold_weighted_summary": f"{round(cold.weight.sum()):,}",
                 "difference_percentage_points": "—", "standardized_mean_difference": "—", "p_value": "—", "test": "—"})
    csv_path = OUTPUT_DIR / "warm_vs_cold_baseline_characteristics.csv"
    with csv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)
    report = ["# Warm Versus Cold AIHA: Baseline Characteristics", "",
              "Warm AIHA is the reference group. Percentages, means, standard deviations, and displayed weighted counts use `DISCWT`. P-values are Welch tests for continuous variables and overall Pearson chi-square tests for categorical variables; missing levels are excluded from categorical tests.", "",
              "The Charlson Comorbidity Index uses all diagnosis positions, standard Quan ICD-10 components and weights, cancer/metastatic-cancer and diabetes/liver hierarchy, and no age points.", "",
              "| Characteristic | Level | Warm unweighted n | Warm weighted summary | Cold unweighted n | Cold weighted summary | Difference, pp | SMD | P-value | Test |",
              "| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |"]
    report.extend(f'| {r["characteristic"]} | {r["level"]} | {r["warm_unweighted_n"]:,} | {r["warm_weighted_summary"]} | {r["cold_unweighted_n"]:,} | {r["cold_weighted_summary"]} | {r["difference_percentage_points"]} | {r["standardized_mean_difference"]} | {r["p_value"]} | {r["test"]} |' for r in rows)
    report_path = OUTPUT_DIR / "warm_vs_cold_baseline_characteristics_report.md"
    report_path.write_text("\n".join(report) + "\n")
    main_report_path = PROJECT_DIR / "report.md"
    marker = "\n<!-- BASELINE_CHARACTERISTICS -->\n"
    main_report = main_report_path.read_text(encoding="utf-8").split(marker)[0].rstrip()
    main_report_path.write_text(main_report + marker + "\n".join(report) + "\n", encoding="utf-8")
    summary = {"reference_group": "Warm AIHA", "comparison_group": "Cold AIHA", "rows": rows,
               "csv": str(csv_path), "report": str(report_path)}
    (OUTPUT_DIR / "warm_vs_cold_baseline_characteristics_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps({"warm_n": len(warm), "cold_n": len(cold), "report": str(report_path)}, indent=2)); return summary


if __name__ == "__main__":
    main()
