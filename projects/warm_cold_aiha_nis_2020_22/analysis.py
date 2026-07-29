"""Preliminary cohort-yield analysis for warm and cold AIHA, NIS 2020–2022."""

from __future__ import annotations

import csv
import json
import os
from pathlib import Path

import duckdb
from dotenv import load_dotenv


PROJECT_DIR = Path(__file__).resolve().parent
REPO_ROOT = PROJECT_DIR.parents[1]
OUTPUT_DIR = PROJECT_DIR / "outputs"
REPORT_PATH = PROJECT_DIR / "report.md"


def normalize(column: str) -> str:
    return f"NULLIF(upper(replace(replace(trim(CAST({column} AS VARCHAR)), '.', ''), ' ', '')), '')"


def dataset_glob() -> str:
    load_dotenv(REPO_ROOT / ".env")
    configured = os.getenv("NIS_DATASET_PATH")
    if not configured:
        raise RuntimeError("Set NIS_DATASET_PATH in the repository .env file.")
    path = Path(configured).expanduser().resolve()
    if not list(path.glob("*.parquet")):
        raise RuntimeError(f"No Parquet files found under {path}")
    return str(path / "*.parquet")


def write_csv(name: str, rows: list[dict]) -> None:
    with (OUTPUT_DIR / name).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)


def main() -> dict:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    dx = [normalize(f"I10_DX{i}") for i in range(1, 4)]
    warm = " OR ".join(f"{code} = 'D5911'" for code in dx)
    cold = " OR ".join(f"{code} = 'D5912'" for code in dx)
    either = f"({warm}) OR ({cold})"
    connection = duckdb.connect()
    denominator_raw = connection.execute("""
        SELECT CAST(YEAR AS INTEGER), count(*)::BIGINT, round(sum(DISCWT), 0)::BIGINT
        FROM read_parquet(?) WHERE AGE >= 18 AND YEAR BETWEEN 2020 AND 2022
        GROUP BY 1 ORDER BY 1
    """, [dataset_glob()]).fetchall()
    primary_raw = connection.execute(f"""
        SELECT {dx[0]} AS code, count(*)::BIGINT, round(sum(DISCWT), 0)::BIGINT
        FROM read_parquet(?)
        WHERE AGE >= 18 AND YEAR BETWEEN 2020 AND 2022 AND {dx[0]} IN ('D5911', 'D5912')
        GROUP BY 1 ORDER BY 1
    """, [dataset_glob()]).fetchall()
    top_three = connection.execute(f"""
        SELECT count(*)::BIGINT, round(sum(DISCWT), 0)::BIGINT,
               count(*) FILTER (WHERE {warm})::BIGINT,
               round(sum(DISCWT) FILTER (WHERE {warm}), 0)::BIGINT,
               count(*) FILTER (WHERE {cold})::BIGINT,
               round(sum(DISCWT) FILTER (WHERE {cold}), 0)::BIGINT,
               count(*) FILTER (WHERE ({warm}) AND ({cold}))::BIGINT,
               round(sum(DISCWT) FILTER (WHERE ({warm}) AND ({cold})), 0)::BIGINT
        FROM read_parquet(?) WHERE AGE >= 18 AND YEAR BETWEEN 2020 AND 2022 AND ({either})
    """, [dataset_glob()]).fetchone()
    annual_raw = connection.execute(f"""
        SELECT CAST(YEAR AS INTEGER),
               count(*) FILTER (WHERE {warm})::BIGINT,
               round(sum(DISCWT) FILTER (WHERE {warm}), 0)::BIGINT,
               count(*) FILTER (WHERE {cold})::BIGINT,
               round(sum(DISCWT) FILTER (WHERE {cold}), 0)::BIGINT,
               count(*)::BIGINT, round(sum(DISCWT), 0)::BIGINT
        FROM read_parquet(?)
        WHERE AGE >= 18 AND YEAR BETWEEN 2020 AND 2022 AND ({either})
        GROUP BY 1 ORDER BY 1
    """, [dataset_glob()]).fetchall()
    connection.close()
    primary = [{"phenotype": "Warm AIHA (D59.11)" if code == "D5911" else "Cold AIHA (D59.12)",
                "unweighted_hospitalizations": n, "weighted_hospitalizations": weighted}
               for code, n, weighted in primary_raw]
    primary.append({"phenotype": "Total unique hospitalizations", "unweighted_hospitalizations": sum(r["unweighted_hospitalizations"] for r in primary),
                    "weighted_hospitalizations": sum(r["weighted_hospitalizations"] for r in primary)})
    expanded = [
        {"phenotype": "Warm AIHA (D59.11)", "unweighted_hospitalizations": top_three[2], "weighted_hospitalizations": top_three[3]},
        {"phenotype": "Cold AIHA (D59.12)", "unweighted_hospitalizations": top_three[4], "weighted_hospitalizations": top_three[5]},
        {"phenotype": "Both warm and cold codes", "unweighted_hospitalizations": top_three[6], "weighted_hospitalizations": top_three[7]},
        {"phenotype": "Total unique hospitalizations", "unweighted_hospitalizations": top_three[0], "weighted_hospitalizations": top_three[1]},
    ]
    annual = [{"admission_year": year, "warm_unweighted": wn, "warm_weighted": ww,
               "cold_unweighted": cn, "cold_weighted": cw, "unique_unweighted": n, "unique_weighted": weighted}
              for year, wn, ww, cn, cw, n, weighted in annual_raw]
    annual.append({"admission_year": "Total", "warm_unweighted": top_three[2], "warm_weighted": top_three[3],
                   "cold_unweighted": top_three[4], "cold_weighted": top_three[5],
                   "unique_unweighted": top_three[0], "unique_weighted": top_three[1]})
    denominators = [{"admission_year": year, "unweighted_all_adult_hospitalizations": n,
                     "weighted_all_adult_hospitalizations": weighted}
                    for year, n, weighted in denominator_raw]
    denominators.append({"admission_year": "Total",
                         "unweighted_all_adult_hospitalizations": sum(r["unweighted_all_adult_hospitalizations"] for r in denominators),
                         "weighted_all_adult_hospitalizations": sum(r["weighted_all_adult_hospitalizations"] for r in denominators)})
    write_csv("primary_diagnosis_cohort_yield.csv", primary)
    write_csv("top_three_diagnosis_cohort_yield.csv", expanded)
    write_csv("top_three_cohort_yield_by_year.csv", annual)
    write_csv("all_adult_hospitalizations_by_year.csv", denominators)
    report = ["# Warm and Cold AIHA — NIS 2020–2022", "",
              "## Preliminary cohort definition", "",
              "Population: all adult NIS hospitalizations (`AGE >= 18`) from 2020 through 2022. Warm AIHA is exact normalized `D59.11`; cold AIHA is exact normalized `D59.12`.", "",
              "## All Adult Hospitalization Denominators", "",
              "| Admission year | Unweighted adult hospitalizations | DISCWT-weighted adult hospitalizations |", "| --- | ---: | ---: |"]
    report.extend(f'| {r["admission_year"]} | {r["unweighted_all_adult_hospitalizations"]:,} | {r["weighted_all_adult_hospitalizations"]:,} |' for r in denominators)
    report.extend(["", "## Primary-diagnosis definition", "",
              "| Phenotype | Unweighted hospitalizations | DISCWT-weighted hospitalizations |", "| --- | ---: | ---: |"])
    report.extend(f'| {r["phenotype"]} | {r["unweighted_hospitalizations"]:,} | {r["weighted_hospitalizations"]:,} |' for r in primary)
    report.extend(["", "## Expanded definition: first three diagnosis positions", "",
                   "| Phenotype | Unweighted hospitalizations | DISCWT-weighted hospitalizations |", "| --- | ---: | ---: |"])
    report.extend(f'| {r["phenotype"]} | {r["unweighted_hospitalizations"]:,} | {r["weighted_hospitalizations"]:,} |' for r in expanded)
    report.extend(["", "Code-specific counts overlap when both codes are present; the total row counts unique hospitalizations.", "",
                   "## Annual distribution using the first three positions", "",
                   "| Year | Warm unweighted | Warm weighted | Cold unweighted | Cold weighted | Unique unweighted | Unique weighted |", "| --- | ---: | ---: | ---: | ---: | ---: | ---: |"])
    report.extend(f'| {r["admission_year"]} | {r["warm_unweighted"]:,} | {r["warm_weighted"]:,} | {r["cold_unweighted"]:,} | {r["cold_weighted"]:,} | {r["unique_unweighted"]:,} | {r["unique_weighted"]:,} |' for r in annual)
    report.extend(["", "These preliminary counts are hospitalization-based and do not identify unique patients. The definitive position rule and exclusion criteria remain to be specified in the study protocol.", ""])
    REPORT_PATH.write_text("\n".join(report), encoding="utf-8")
    summary = {"project": "Warm and Cold AIHA NIS 2020–2022", "all_adult_denominators": denominators,
               "primary": primary, "top_three": expanded,
               "annual_top_three": annual, "report": str(REPORT_PATH), "output_directory": str(OUTPUT_DIR)}
    (OUTPUT_DIR / "cohort_yield_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2)); return summary


if __name__ == "__main__":
    main()
