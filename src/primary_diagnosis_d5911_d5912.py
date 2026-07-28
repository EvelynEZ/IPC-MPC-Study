"""Count D59.11/D59.12 when recorded as the primary diagnosis."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import duckdb

from src.phase_1_2 import load_dataset_files, normalize_code_expression


ROOT = Path(__file__).resolve().parents[1]
OUTPUT_DIR = ROOT / "outputs/septic_shock/primary_diagnosis_d5911_d5912"


def main() -> dict:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    files = load_dataset_files()
    parquet_glob = str(files[0].parent / "*.parquet")
    connection = duckdb.connect()
    primary_code = normalize_code_expression("I10_DX1")
    rows = connection.execute("""
        SELECT {primary_code} AS normalized_primary_diagnosis,
               count(*)::BIGINT AS unweighted_hospitalizations,
               round(sum(DISCWT), 0)::BIGINT AS weighted_hospitalizations
        FROM read_parquet(?)
        WHERE AGE >= 18 AND {primary_code} IN ('D5911', 'D5912')
        GROUP BY 1 ORDER BY 1
    """.format(primary_code=primary_code), [parquet_glob]).fetchall()
    annual_raw = connection.execute("""
        SELECT CAST(YEAR AS INTEGER) AS admission_year,
               {primary_code} AS normalized_primary_diagnosis,
               count(*)::BIGINT AS unweighted_hospitalizations,
               round(sum(DISCWT), 0)::BIGINT AS weighted_hospitalizations
        FROM read_parquet(?)
        WHERE AGE >= 18 AND {primary_code} IN ('D5911', 'D5912')
        GROUP BY 1, 2 ORDER BY 1, 2
    """.format(primary_code=primary_code), [parquet_glob]).fetchall()
    connection.close()
    results = [{"primary_diagnosis": code[:3] + "." + code[3:],
                "unweighted_hospitalizations": unweighted,
                "weighted_hospitalizations": weighted} for code, unweighted, weighted in rows]
    results.append({"primary_diagnosis": "Total: D59.11 or D59.12",
                    "unweighted_hospitalizations": sum(row["unweighted_hospitalizations"] for row in results),
                    "weighted_hospitalizations": sum(row["weighted_hospitalizations"] for row in results)})
    annual = []
    lookup = {(year, code): (unweighted, weighted) for year, code, unweighted, weighted in annual_raw}
    for year in range(2016, 2023):
        d5911 = lookup.get((year, "D5911"), (0, 0))
        d5912 = lookup.get((year, "D5912"), (0, 0))
        annual.append({"admission_year": year, "d5911_unweighted": d5911[0], "d5911_weighted": d5911[1],
                       "d5912_unweighted": d5912[0], "d5912_weighted": d5912[1],
                       "combined_unweighted": d5911[0] + d5912[0], "combined_weighted": d5911[1] + d5912[1]})
    annual.append({"admission_year": "Total",
                   "d5911_unweighted": sum(row["d5911_unweighted"] for row in annual),
                   "d5911_weighted": sum(row["d5911_weighted"] for row in annual),
                   "d5912_unweighted": sum(row["d5912_unweighted"] for row in annual),
                   "d5912_weighted": sum(row["d5912_weighted"] for row in annual),
                   "combined_unweighted": sum(row["combined_unweighted"] for row in annual),
                   "combined_weighted": sum(row["combined_weighted"] for row in annual)})
    with (OUTPUT_DIR / "primary_diagnosis_counts.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(results[0])); writer.writeheader(); writer.writerows(results)
    with (OUTPUT_DIR / "primary_diagnosis_counts_by_year.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(annual[0])); writer.writeheader(); writer.writerows(annual)
    report = ["# Primary-Discharge-Diagnosis Counts for D59.11 and D59.12", "",
              "Counts include all adult NIS hospitalizations from 2016–2022, without restriction to hematologic malignancy. A code qualifies only when it is an exact normalized match in the first diagnosis field.", "",
              "| Primary discharge diagnosis | Unweighted hospitalizations | DISCWT-weighted hospitalizations |",
              "| --- | ---: | ---: |"]
    report.extend(f'| {row["primary_diagnosis"]} | {row["unweighted_hospitalizations"]:,} | {row["weighted_hospitalizations"]:,} |' for row in results)
    report.extend(["", "## Admission-Year Distribution", "",
                   "| Admission year | D59.11 unweighted | D59.11 weighted | D59.12 unweighted | D59.12 weighted | Combined unweighted | Combined weighted |",
                   "| --- | ---: | ---: | ---: | ---: | ---: | ---: |"])
    report.extend(f'| {row["admission_year"]} | {row["d5911_unweighted"]:,} | {row["d5911_weighted"]:,} | {row["d5912_unweighted"]:,} | {row["d5912_weighted"]:,} | {row["combined_unweighted"]:,} | {row["combined_weighted"]:,} |' for row in annual)
    report.extend(["", "No qualifying adult hospitalization was identified in 2016–2019. The unweighted count is the number of sampled NIS discharge records. The weighted count is the estimated national number of hospitalizations.", ""])
    report_path = OUTPUT_DIR / "primary_diagnosis_counts_report.md"
    report_path.write_text("\n".join(report), encoding="utf-8")
    summary = {"population": "All NIS hospitalizations with AGE >=18, 2016–2022; no HM restriction",
               "definition": "Exact normalized D5911 or D5912 in I10_DX1", "results": results,
               "annual_results": annual,
               "report_path": str(report_path)}
    (OUTPUT_DIR / "primary_diagnosis_counts_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2)); return summary


if __name__ == "__main__":
    main()
