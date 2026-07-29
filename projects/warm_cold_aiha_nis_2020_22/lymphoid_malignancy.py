"""Lymphoid-malignancy prevalence in warm versus cold AIHA."""

from __future__ import annotations

import csv
import json
import math
import sys
from pathlib import Path

import duckdb
import numpy as np
from scipy.stats import chi2_contingency


PROJECT_DIR = Path(__file__).resolve().parent
REPO_ROOT = PROJECT_DIR.parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from src.phase_1_2 import code_match_sql


OUTPUT_DIR = PROJECT_DIR / "outputs"
DATABASE = OUTPUT_DIR / "warm_cold_aiha_cohort.duckdb"
CONFIG = REPO_ROOT / "config/hm_phenotype_v0_1.json"


def format_p(value: float) -> str:
    return "<0.001" if value < 0.001 else f"{value:.3f}"


def smd(cold: float, warm: float) -> float:
    variance = (cold * (1 - cold) + warm * (1 - warm)) / 2
    return 0.0 if variance <= 0 else (cold - warm) / math.sqrt(variance)


def main() -> dict:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    config = json.loads(CONFIG.read_text())
    rules = {rule["id"]: rule for rule in config["subtypes"]}
    component_rules = {
        "Lymphoma spectrum (Hodgkin, NHL, WM/MALT)": rules["lymphoma"],
        "CLL/chronic leukemia group": rules["cll_chronic_leukemia"],
        "Plasma-cell neoplasms": rules["myeloma_plasma_cell"],
    }
    conditions = {label: code_match_sql("code", rule) for label, rule in component_rules.items()}
    flag_expressions = {label: f"list_contains(list_transform(diagnosis_codes, code -> {condition}), TRUE)"
                        for label, condition in conditions.items()}
    any_expression = " OR ".join(f"({expression})" for expression in flag_expressions.values())
    connection = duckdb.connect(str(DATABASE))
    connection.execute("ALTER TABLE aiha_cohort DROP COLUMN IF EXISTS any_lymphoid_malignancy")
    connection.execute(f"ALTER TABLE aiha_cohort ADD COLUMN any_lymphoid_malignancy BOOLEAN")
    connection.execute(f"UPDATE aiha_cohort SET any_lymphoid_malignancy = ({any_expression})")
    totals = {row[0]: {"n": int(row[1]), "weighted": float(row[2])}
              for row in connection.execute("SELECT aiha_type, count(*), sum(DISCWT) FROM aiha_cohort GROUP BY 1").fetchall()}
    rows = []
    all_flags = list(flag_expressions.items()) + [("ANY_LYMPHOID_MALIGNANCY", "any_lymphoid_malignancy")]
    for label, expression in all_flags:
        by_group = {row[0]: {"n": int(row[1]), "weighted": float(row[2] or 0)} for row in connection.execute(
            f"SELECT aiha_type, count(*) FILTER (WHERE {expression}), sum(DISCWT) FILTER (WHERE {expression}) FROM aiha_cohort GROUP BY 1"
        ).fetchall()}
        wn, cn = by_group["Warm AIHA"]["n"], by_group["Cold AIHA"]["n"]
        wp, cp = by_group["Warm AIHA"]["weighted"] / totals["Warm AIHA"]["weighted"], by_group["Cold AIHA"]["weighted"] / totals["Cold AIHA"]["weighted"]
        contingency = [[wn, totals["Warm AIHA"]["n"] - wn], [cn, totals["Cold AIHA"]["n"] - cn]]
        p = float(chi2_contingency(contingency, correction=False)[1])
        rows.append({"malignancy_definition": label, "warm_unweighted_n": wn,
                     "warm_weighted_n": round(by_group["Warm AIHA"]["weighted"]), "warm_weighted_percent": round(100 * wp, 2),
                     "cold_unweighted_n": cn, "cold_weighted_n": round(by_group["Cold AIHA"]["weighted"]), "cold_weighted_percent": round(100 * cp, 2),
                     "difference_percentage_points": round(100 * (cp - wp), 2), "standardized_mean_difference": round(smd(cp, wp), 3),
                     "p_value": format_p(p)})
    rows.append({"malignancy_definition": "Total cohort denominator", "warm_unweighted_n": totals["Warm AIHA"]["n"],
                 "warm_weighted_n": round(totals["Warm AIHA"]["weighted"]), "warm_weighted_percent": 100.0,
                 "cold_unweighted_n": totals["Cold AIHA"]["n"], "cold_weighted_n": round(totals["Cold AIHA"]["weighted"]), "cold_weighted_percent": 100.0,
                 "difference_percentage_points": "—", "standardized_mean_difference": "—", "p_value": "—"})
    connection.close()
    csv_path = OUTPUT_DIR / "any_lymphoid_malignancy_by_aiha_type.csv"
    with csv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)
    report = ["# Any Lymphoid Malignancy by AIHA Type", "",
              "Warm AIHA is the reference group. The phenotype searches all 40 diagnosis positions and reuses the prior HM project’s lymphoma, CLL/chronic-leukemia, and plasma-cell inclusion and exclusion rules.", "",
              "Component rows can overlap. `ANY_LYMPHOID_MALIGNANCY` is their union and counts each hospitalization once.", "",
              "| Malignancy definition | Warm unweighted n | Warm weighted n (%) | Cold unweighted n | Cold weighted n (%) | Difference, pp | SMD | P-value |",
              "| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |"]
    report.extend(f'| {r["malignancy_definition"]} | {r["warm_unweighted_n"]:,} | {r["warm_weighted_n"]:,} ({r["warm_weighted_percent"]:.2f}%) | {r["cold_unweighted_n"]:,} | {r["cold_weighted_n"]:,} ({r["cold_weighted_percent"]:.2f}%) | {r["difference_percentage_points"]} | {r["standardized_mean_difference"]} | {r["p_value"]} |' for r in rows)
    report_path = OUTPUT_DIR / "any_lymphoid_malignancy_report.md"
    report_path.write_text("\n".join(report) + "\n")
    main_report_path = PROJECT_DIR / "report.md"; marker = "\n<!-- LYMPHOID_MALIGNANCY -->\n"
    main_report = main_report_path.read_text().split(marker)[0].rstrip()
    main_report_path.write_text(main_report + marker + "\n".join(report) + "\n")
    summary = {"reference_group": "Warm AIHA", "phenotype_config": str(CONFIG), "rows": rows,
               "csv": str(csv_path), "report": str(report_path)}
    (OUTPUT_DIR / "any_lymphoid_malignancy_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2)); return summary


if __name__ == "__main__":
    main()
