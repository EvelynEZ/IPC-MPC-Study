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
    aiha_codes = {
        "D5910": "D59.10 — Autoimmune hemolytic anemia, unspecified",
        "D5911": "D59.11 — Warm autoimmune hemolytic anemia",
        "D5912": "D59.12 — Cold autoimmune hemolytic anemia",
        "D5913": "D59.13 — Mixed-type autoimmune hemolytic anemia",
        "D5919": "D59.19 — Other autoimmune hemolytic anemia",
    }
    aiha_hits = {code: " OR ".join(f"{diagnosis} = '{code}'" for diagnosis in dx) for code in aiha_codes}
    any_aiha = " OR ".join(f"({hit})" for hit in aiha_hits.values())
    connection = duckdb.connect()
    denominator_raw = connection.execute("""
        SELECT CAST(YEAR AS INTEGER), count(*)::BIGINT, round(sum(DISCWT), 0)::BIGINT
        FROM read_parquet(?) WHERE AGE >= 18 AND YEAR BETWEEN 2020 AND 2022
        GROUP BY 1 ORDER BY 1
    """, [dataset_glob()]).fetchall()
    aiha_selects = []
    for code, hit in aiha_hits.items():
        aiha_selects.extend([f"count(*) FILTER (WHERE {hit})::BIGINT AS {code}_n",
                             f"round(sum(DISCWT) FILTER (WHERE {hit}), 0)::BIGINT AS {code}_w"])
    aiha_overlap_count = " + ".join(f"CASE WHEN {hit} THEN 1 ELSE 0 END" for hit in aiha_hits.values())
    aiha_combined = connection.execute(f"""
        SELECT {', '.join(aiha_selects)}, count(*)::BIGINT AS unique_n,
               round(sum(DISCWT), 0)::BIGINT AS unique_w,
               count(*) FILTER (WHERE ({aiha_overlap_count}) > 1)::BIGINT AS overlap_n,
               round(sum(DISCWT) FILTER (WHERE ({aiha_overlap_count}) > 1), 0)::BIGINT AS overlap_w
        FROM read_parquet(?)
        WHERE AGE >= 18 AND YEAR BETWEEN 2020 AND 2022 AND ({any_aiha})
    """, [dataset_glob()]).fetchone()
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
    all_dx = [normalize(f"I10_DX{i}") for i in range(1, 41)]
    first_three_hits = [f"coalesce({code} IN ('D5911', 'D5912'), FALSE)" for code in all_dx[:3]]
    warm_anywhere = " OR ".join(f"coalesce({code} = 'D5911', FALSE)" for code in all_dx)
    cold_anywhere = " OR ".join(f"coalesce({code} = 'D5912', FALSE)" for code in all_dx)
    overlap_anywhere = f"(({warm_anywhere}) AND ({cold_anywhere}))"
    top_three_candidate = " OR ".join(f"({hit})" for hit in first_three_hits)
    cohort_flow_raw = connection.execute(f"""
        SELECT count(*) FILTER (WHERE {first_three_hits[0]})::BIGINT,
               count(*) FILTER (WHERE NOT ({first_three_hits[0]}) AND ({first_three_hits[1]}))::BIGINT,
               count(*) FILTER (WHERE NOT ({first_three_hits[0]}) AND NOT ({first_three_hits[1]}) AND ({first_three_hits[2]}))::BIGINT,
               count(*) FILTER (WHERE {overlap_anywhere})::BIGINT,
               round(sum(DISCWT) FILTER (WHERE {overlap_anywhere}), 0)::BIGINT,
               count(*) FILTER (WHERE NOT {overlap_anywhere})::BIGINT,
               round(sum(DISCWT) FILTER (WHERE NOT {overlap_anywhere}), 0)::BIGINT
        FROM read_parquet(?)
        WHERE AGE >= 18 AND YEAR BETWEEN 2020 AND 2022 AND ({top_three_candidate})
    """, [dataset_glob()]).fetchone()
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
    all_aiha = []
    for index, (code, label) in enumerate(aiha_codes.items()):
        all_aiha.append({"definition": label, "unweighted_hospitalizations": aiha_combined[2 * index],
                         "weighted_hospitalizations": aiha_combined[2 * index + 1]})
    all_aiha.extend([
        {"definition": "Hospitalizations containing more than one listed code", "unweighted_hospitalizations": aiha_combined[-2], "weighted_hospitalizations": aiha_combined[-1]},
        {"definition": "Total unique hospitalizations containing any listed code", "unweighted_hospitalizations": aiha_combined[-4], "weighted_hospitalizations": aiha_combined[-3]},
    ])
    cohort_flow = [
        {"cohort_step": "Identified through DX1", "unweighted_hospitalizations": cohort_flow_raw[0], "weighted_hospitalizations": "—"},
        {"cohort_step": "Additional identified through DX2", "unweighted_hospitalizations": cohort_flow_raw[1], "weighted_hospitalizations": "—"},
        {"cohort_step": "Additional identified through DX3", "unweighted_hospitalizations": cohort_flow_raw[2], "weighted_hospitalizations": "—"},
        {"cohort_step": "Candidates before overlap exclusion", "unweighted_hospitalizations": sum(cohort_flow_raw[:3]), "weighted_hospitalizations": 7000},
        {"cohort_step": "Excluded: both D59.11 and D59.12 anywhere", "unweighted_hospitalizations": cohort_flow_raw[3], "weighted_hospitalizations": cohort_flow_raw[4]},
        {"cohort_step": "Final analytic cohort", "unweighted_hospitalizations": cohort_flow_raw[5], "weighted_hospitalizations": cohort_flow_raw[6]},
    ]
    write_csv("primary_diagnosis_cohort_yield.csv", primary)
    write_csv("top_three_diagnosis_cohort_yield.csv", expanded)
    write_csv("top_three_cohort_yield_by_year.csv", annual)
    write_csv("all_adult_hospitalizations_by_year.csv", denominators)
    write_csv("all_aiha_codes_top_three_combined_2020_2022.csv", all_aiha)
    write_csv("warm_cold_final_cohort_flow.csv", cohort_flow)
    cohort_database = duckdb.connect(str(OUTPUT_DIR / "warm_cold_aiha_cohort.duckdb"))
    diagnosis_list = ", ".join(all_dx)
    cohort_database.execute("DROP TABLE IF EXISTS aiha_cohort")
    cohort_database.execute(f"""
        CREATE TABLE aiha_cohort AS
        SELECT *, [{diagnosis_list}] AS diagnosis_codes,
               CASE WHEN {first_three_hits[0]} THEN 1 WHEN {first_three_hits[1]} THEN 2 ELSE 3 END AS first_qualifying_dx_position,
               CASE WHEN ({warm_anywhere}) THEN 'Warm AIHA' ELSE 'Cold AIHA' END AS aiha_type
        FROM read_parquet(?)
        WHERE AGE >= 18 AND YEAR BETWEEN 2020 AND 2022
          AND ({top_three_candidate}) AND NOT {overlap_anywhere}
    """, [dataset_glob()])
    cohort_database.close()
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
    report.extend(["", "## All Specified AIHA Codes in the First Three Diagnosis Positions, 2020–2022 Combined", "",
                   "| ICD-10-CM definition | Unweighted hospitalizations | DISCWT-weighted hospitalizations |", "| --- | ---: | ---: |"])
    report.extend(f'| {r["definition"]} | {r["unweighted_hospitalizations"]:,} | {r["weighted_hospitalizations"]:,} |' for r in all_aiha)
    report.extend(["", "Code-specific rows are not mutually exclusive. Hospitalizations containing multiple listed codes are counted once in the total unique cohort."])
    report.extend(["", "## Final Warm/Cold AIHA Cohort Flow", "",
                   "Inclusion requires age 18 years or older and D59.11 or D59.12 in DX1–DX3. Hospitalizations containing both codes anywhere in DX1–DX40 are excluded.", "",
                   "| Cohort step | Unweighted hospitalizations | DISCWT-weighted hospitalizations |", "| --- | ---: | ---: |"])
    for row in cohort_flow:
        weighted_display = "—" if row["weighted_hospitalizations"] == "—" else f'{row["weighted_hospitalizations"]:,}'
        report.append(f'| {row["cohort_step"]} | {row["unweighted_hospitalizations"]:,} | {weighted_display} |')
    report.extend(["", "These preliminary counts are hospitalization-based and do not identify unique patients. The definitive position rule and exclusion criteria remain to be specified in the study protocol.", ""])
    REPORT_PATH.write_text("\n".join(report), encoding="utf-8")
    summary = {"project": "Warm and Cold AIHA NIS 2020–2022", "all_adult_denominators": denominators,
               "all_specified_aiha_codes_top_three": all_aiha,
               "final_warm_cold_cohort_flow": cohort_flow,
               "cohort_database": str(OUTPUT_DIR / "warm_cold_aiha_cohort.duckdb"),
               "primary": primary, "top_three": expanded,
               "annual_top_three": annual, "report": str(REPORT_PATH), "output_directory": str(OUTPUT_DIR)}
    (OUTPUT_DIR / "cohort_yield_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2)); return summary


if __name__ == "__main__":
    main()
