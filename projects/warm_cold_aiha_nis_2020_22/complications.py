"""Selected inpatient complications by warm versus cold AIHA type."""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path

import duckdb
from scipy.stats import chi2_contingency


PROJECT_DIR = Path(__file__).resolve().parent
OUTPUT_DIR = PROJECT_DIR / "outputs"
DATABASE = OUTPUT_DIR / "warm_cold_aiha_cohort.duckdb"

COMPLICATIONS = {
    "Acute kidney injury": {"prefix": "N17"},
    "Acute respiratory failure": {"prefix": "J960"},
    "Sepsis": {"prefix": "A41"},
}


def format_p(value: float) -> str:
    return "<0.001" if value < 0.001 else f"{value:.3f}"


def smd(cold: float, warm: float) -> float:
    variance = (cold * (1 - cold) + warm * (1 - warm)) / 2
    return 0.0 if variance <= 0 else (cold - warm) / math.sqrt(variance)


def main() -> dict:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    connection = duckdb.connect(str(DATABASE))
    totals = {
        row[0]: {"n": int(row[1]), "weighted": float(row[2])}
        for row in connection.execute(
            "SELECT aiha_type, count(*), sum(DISCWT) FROM aiha_cohort GROUP BY 1"
        ).fetchall()
    }
    rows = []
    for label, rule in COMPLICATIONS.items():
        prefix = rule["prefix"]
        expression = (
            "list_contains(list_transform(diagnosis_codes, "
            f"code -> starts_with(code, '{prefix}')), TRUE)"
        )
        by_group = {
            row[0]: {"n": int(row[1]), "weighted": float(row[2] or 0)}
            for row in connection.execute(
                f"""SELECT aiha_type,
                           count(*) FILTER (WHERE {expression}),
                           sum(DISCWT) FILTER (WHERE {expression})
                    FROM aiha_cohort GROUP BY 1"""
            ).fetchall()
        }
        wn, cn = by_group["Warm AIHA"]["n"], by_group["Cold AIHA"]["n"]
        wp = by_group["Warm AIHA"]["weighted"] / totals["Warm AIHA"]["weighted"]
        cp = by_group["Cold AIHA"]["weighted"] / totals["Cold AIHA"]["weighted"]
        p_value = float(chi2_contingency(
            [[wn, totals["Warm AIHA"]["n"] - wn],
             [cn, totals["Cold AIHA"]["n"] - cn]], correction=False
        )[1])
        rows.append({
            "complication": label,
            "icd_10_cm_definition": f"{prefix}*",
            "warm_unweighted_n": wn,
            "warm_weighted_n": round(by_group["Warm AIHA"]["weighted"]),
            "warm_weighted_percent": round(100 * wp, 2),
            "cold_unweighted_n": cn,
            "cold_weighted_n": round(by_group["Cold AIHA"]["weighted"]),
            "cold_weighted_percent": round(100 * cp, 2),
            "difference_percentage_points": round(100 * (cp - wp), 2),
            "standardized_mean_difference": round(smd(cp, wp), 3),
            "p_value": format_p(p_value),
        })
    connection.close()
    rows.append({
        "complication": "Total cohort denominator", "icd_10_cm_definition": "—",
        "warm_unweighted_n": totals["Warm AIHA"]["n"],
        "warm_weighted_n": round(totals["Warm AIHA"]["weighted"]),
        "warm_weighted_percent": 100.0,
        "cold_unweighted_n": totals["Cold AIHA"]["n"],
        "cold_weighted_n": round(totals["Cold AIHA"]["weighted"]),
        "cold_weighted_percent": 100.0,
        "difference_percentage_points": "—", "standardized_mean_difference": "—",
        "p_value": "—",
    })
    csv_path = OUTPUT_DIR / "selected_complications_by_aiha_type.csv"
    with csv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)
    report = [
        "# Selected Complications by AIHA Type", "",
        "Warm AIHA is the reference group. Diagnoses were identified in any of the 40 diagnosis positions after removing decimal points and spaces and converting codes to uppercase. Acute respiratory failure uses ICD-10-CM `J96.0*` (`J960*` normalized); `A960` was treated as a typographical error because it does not denote acute respiratory failure.", "",
        "| Complication | ICD-10-CM | Warm unweighted n | Warm weighted n (%) | Cold unweighted n | Cold weighted n (%) | Difference, pp | SMD | P-value |",
        "| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    report.extend(
        f'| {r["complication"]} | {r["icd_10_cm_definition"]} | {r["warm_unweighted_n"]:,} | '
        f'{r["warm_weighted_n"]:,} ({r["warm_weighted_percent"]:.2f}%) | {r["cold_unweighted_n"]:,} | '
        f'{r["cold_weighted_n"]:,} ({r["cold_weighted_percent"]:.2f}%) | {r["difference_percentage_points"]} | '
        f'{r["standardized_mean_difference"]} | {r["p_value"]} |' for r in rows
    )
    report.extend(["", "P-values are two-sided Pearson chi-square tests comparing warm and cold AIHA. Diagnoses are co-documented during the hospitalization and do not establish temporal ordering.", ""])
    report_path = OUTPUT_DIR / "selected_complications_report.md"
    report_path.write_text("\n".join(report), encoding="utf-8")
    main_report_path = PROJECT_DIR / "report.md"
    marker = "\n<!-- SELECTED_COMPLICATIONS -->\n"
    main_report = main_report_path.read_text().split(marker)[0].rstrip()
    main_report_path.write_text(main_report + marker + "\n".join(report), encoding="utf-8")
    summary = {"reference_group": "Warm AIHA", "rows": rows, "csv": str(csv_path), "report": str(report_path)}
    (OUTPUT_DIR / "selected_complications_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    return summary


if __name__ == "__main__":
    main()
