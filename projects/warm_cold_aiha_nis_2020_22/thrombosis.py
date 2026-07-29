"""Acute venous and arterial thrombosis phenotypes by AIHA type."""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path

import duckdb
import numpy as np
from scipy.stats import chi2_contingency


PROJECT_DIR = Path(__file__).resolve().parent
OUTPUT_DIR = PROJECT_DIR / "outputs"
DATABASE = OUTPUT_DIR / "warm_cold_aiha_cohort.duckdb"


FAMILIES = {
    "Acute pulmonary embolism": {"exact": ["I2602", "I2609", "I2692", "I2693", "I2694", "I2699"]},
    "Acute deep-vein thrombosis": {"exact": [
        "I82401", "I82402", "I82403", "I82409", "I82411", "I82412", "I82413", "I82419",
        "I82421", "I82422", "I82423", "I82429", "I82431", "I82432", "I82433", "I82439",
        "I82441", "I82442", "I82443", "I82449", "I82451", "I82452", "I82453", "I82459",
        "I82461", "I82462", "I82463", "I82469", "I82491", "I82492", "I82493", "I82499",
        "I824Y1", "I824Y2", "I824Y3", "I824Y9", "I824Z1", "I824Z2", "I824Z3", "I824Z9",
        "I82601", "I82602", "I82603", "I82609", "I82621", "I82622", "I82623", "I82629",
        "I82A11", "I82A12", "I82A13", "I82A19", "I82B11", "I82B12", "I82B13", "I82B19",
        "I82C11", "I82C12", "I82C13", "I82C19",
    ]},
    "Splanchnic-vein thrombosis": {"exact": ["I81", "I820", "I823"]},
    "Other acute venous thrombosis": {"exact": [
        "I636", "I82210", "I82220", "I82290", "I82811", "I82812", "I82813", "I82819", "I82890",
    ]},
    "Acute ischemic stroke": {"prefix": ["I630", "I631", "I632", "I633", "I634", "I635", "I638", "I639"]},
    "Acute myocardial infarction": {"prefix": ["I210", "I211", "I212", "I213", "I214", "I219", "I220", "I221", "I222", "I228", "I229"]},
    "Other acute arterial embolism or thrombosis": {"prefix": ["I740", "I741", "I742", "I743", "I744", "I745", "I748", "I749"]},
}


def code_condition(rule: dict, code: str = "code") -> str:
    clauses = [f"{code} = '{value}'" for value in rule.get("exact", [])]
    clauses += [f"starts_with({code}, '{value}')" for value in rule.get("prefix", [])]
    return "(" + " OR ".join(clauses) + ")"


def format_p(value: float) -> str:
    return "<0.001" if value < 0.001 else f"{value:.3f}"


def smd(cold: float, warm: float) -> float:
    variance = (cold * (1 - cold) + warm * (1 - warm)) / 2
    return 0.0 if variance <= 0 else (cold - warm) / math.sqrt(variance)


def bh_adjust(values: list[float]) -> list[float]:
    order = np.argsort(values); adjusted = np.empty(len(values)); running = 1.0
    for reverse_rank, index in enumerate(order[::-1], start=1):
        rank = len(values) - reverse_rank + 1
        running = min(running, values[index] * len(values) / rank)
        adjusted[index] = running
    return adjusted.tolist()


def main() -> dict:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    expressions = {label: f"list_contains(list_transform(diagnosis_codes, code -> {code_condition(rule)}), TRUE)"
                   for label, rule in FAMILIES.items()}
    vte_labels = list(FAMILIES)[:4]; arterial_labels = list(FAMILIES)[4:]
    expressions["ANY_VTE"] = " OR ".join(f"({expressions[label]})" for label in vte_labels)
    expressions["ANY_ARTERIAL_THROMBOSIS"] = " OR ".join(f"({expressions[label]})" for label in arterial_labels)
    expressions["ANY_ACUTE_THROMBOSIS"] = f"({expressions['ANY_VTE']}) OR ({expressions['ANY_ARTERIAL_THROMBOSIS']})"
    connection = duckdb.connect(str(DATABASE))
    column_names = {label: label.lower().replace("-", "_").replace("/", "_").replace(" ", "_") for label in expressions}
    for label in ["ANY_VTE", "ANY_ARTERIAL_THROMBOSIS", "ANY_ACUTE_THROMBOSIS"]:
        column = column_names[label]
        connection.execute(f"ALTER TABLE aiha_cohort DROP COLUMN IF EXISTS {column}")
        connection.execute(f"ALTER TABLE aiha_cohort ADD COLUMN {column} BOOLEAN")
        connection.execute(f"UPDATE aiha_cohort SET {column} = ({expressions[label]})")
    totals = {row[0]: {"n": int(row[1]), "weighted": float(row[2])}
              for row in connection.execute("SELECT aiha_type, count(*), sum(DISCWT) FROM aiha_cohort GROUP BY 1").fetchall()}
    rows, p_values = [], []
    for label, expression in expressions.items():
        by_group = {row[0]: {"n": int(row[1]), "weighted": float(row[2] or 0)} for row in connection.execute(
            f"SELECT aiha_type, count(*) FILTER (WHERE {expression}), sum(DISCWT) FILTER (WHERE {expression}) FROM aiha_cohort GROUP BY 1"
        ).fetchall()}
        wn, cn = by_group["Warm AIHA"]["n"], by_group["Cold AIHA"]["n"]
        wp = by_group["Warm AIHA"]["weighted"] / totals["Warm AIHA"]["weighted"]
        cp = by_group["Cold AIHA"]["weighted"] / totals["Cold AIHA"]["weighted"]
        p = float(chi2_contingency([[wn, totals["Warm AIHA"]["n"] - wn], [cn, totals["Cold AIHA"]["n"] - cn]], correction=False)[1])
        p_values.append(p)
        rows.append({"diagnosis_family": label, "warm_unweighted_n": wn, "warm_weighted_n": round(by_group["Warm AIHA"]["weighted"]),
                     "warm_weighted_percent": round(100 * wp, 2), "cold_unweighted_n": cn,
                     "cold_weighted_n": round(by_group["Cold AIHA"]["weighted"]), "cold_weighted_percent": round(100 * cp, 2),
                     "difference_percentage_points": round(100 * (cp - wp), 2), "standardized_mean_difference": round(smd(cp, wp), 3),
                     "p_value": format_p(p), "fdr_adjusted_p_value": ""})
    for row, adjusted in zip(rows, bh_adjust(p_values)):
        row["fdr_adjusted_p_value"] = format_p(adjusted)
    rows.append({"diagnosis_family": "Total cohort denominator", "warm_unweighted_n": totals["Warm AIHA"]["n"],
                 "warm_weighted_n": round(totals["Warm AIHA"]["weighted"]), "warm_weighted_percent": 100.0,
                 "cold_unweighted_n": totals["Cold AIHA"]["n"], "cold_weighted_n": round(totals["Cold AIHA"]["weighted"]),
                 "cold_weighted_percent": 100.0, "difference_percentage_points": "—", "standardized_mean_difference": "—",
                 "p_value": "—", "fdr_adjusted_p_value": "—"})
    connection.close()
    csv_path = OUTPUT_DIR / "acute_thrombosis_by_aiha_type.csv"
    with csv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)
    config_path = OUTPUT_DIR / "acute_thrombosis_code_families.json"
    config_path.write_text(json.dumps(FAMILIES, indent=2) + "\n")
    report = ["# Acute Thrombosis by AIHA Type", "",
              "Warm AIHA is the reference group. Every phenotype searches all 40 diagnosis positions. Code-specific families can overlap; each composite counts a hospitalization once.", "",
              "`ANY_VTE` is the union of pulmonary embolism, acute DVT, splanchnic-vein thrombosis, and other acute venous thrombosis. `ANY_ARTERIAL_THROMBOSIS` is the union of acute ischemic stroke, acute MI, and other acute arterial embolism/thrombosis. `ANY_ACUTE_THROMBOSIS` is the union of those two composites.", "",
              "| Diagnosis family | Warm unweighted n | Warm weighted n (%) | Cold unweighted n | Cold weighted n (%) | Difference, pp | SMD | P-value | FDR-adjusted p |",
              "| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |"]
    report.extend(f'| {r["diagnosis_family"]} | {r["warm_unweighted_n"]:,} | {r["warm_weighted_n"]:,} ({r["warm_weighted_percent"]:.2f}%) | {r["cold_unweighted_n"]:,} | {r["cold_weighted_n"]:,} ({r["cold_weighted_percent"]:.2f}%) | {r["difference_percentage_points"]} | {r["standardized_mean_difference"]} | {r["p_value"]} | {r["fdr_adjusted_p_value"]} |' for r in rows)
    report_path = OUTPUT_DIR / "acute_thrombosis_report.md"; report_path.write_text("\n".join(report) + "\n")
    main_report_path = PROJECT_DIR / "report.md"; marker = "\n<!-- ACUTE_THROMBOSIS -->\n"
    main_report = main_report_path.read_text().split(marker)[0].rstrip(); main_report_path.write_text(main_report + marker + "\n".join(report) + "\n")
    summary = {"reference_group": "Warm AIHA", "rows": rows, "code_config": str(config_path), "csv": str(csv_path), "report": str(report_path)}
    (OUTPUT_DIR / "acute_thrombosis_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2)); return summary


if __name__ == "__main__":
    main()
