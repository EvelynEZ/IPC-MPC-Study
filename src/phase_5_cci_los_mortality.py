"""Compare cancer-excluded Charlson score, LOS, and mortality by sepsis."""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path
from typing import Any

import duckdb
from scipy.stats import chi2_contingency, ttest_ind_from_stats


REPO_ROOT = Path(__file__).resolve().parents[1]
COHORT_DATABASE = REPO_ROOT / "outputs/phase_1_2/hm_cohort.duckdb"
OUTPUT_DIR = REPO_ROOT / "outputs/phase_5"
# Quan et al. ICD-10 adaptation of the Charlson conditions. The cancer and
# metastatic-cancer components are deliberately absent for this HM analysis.
CCI_COMPONENTS = {
    "ami": ("I21", "I22", "I252"),
    "chf": ("I099", "I110", "I130", "I132", "I255", "I420", "I425", "I426", "I427", "I428", "I429", "I43", "I50", "P290"),
    "pvd": ("I70", "I71", "I731", "I738", "I739", "I771", "I790", "I792", "K551", "K558", "K559", "Z958", "Z959"),
    "cevd": ("G45", "G46", "H340", "I60", "I61", "I62", "I63", "I64", "I65", "I66", "I67", "I68", "I69"),
    "dementia": ("F00", "F01", "F02", "F03", "F051", "G30", "G311"),
    "copd": ("I278", "I279", "J40", "J41", "J42", "J43", "J44", "J45", "J46", "J47", "J60", "J61", "J62", "J63", "J64", "J65", "J66", "J67", "J684", "J701", "J703"),
    "rheumd": ("M05", "M06", "M315", "M32", "M33", "M34", "M351", "M353", "M360"),
    "pud": ("K25", "K26", "K27", "K28"),
    "mld": ("B18", "K700", "K701", "K702", "K703", "K709", "K713", "K714", "K715", "K717", "K73", "K74", "K760", "K762", "K763", "K764", "K768", "K769", "Z944"),
    "diab": ("E100", "E101", "E106", "E108", "E109", "E110", "E111", "E116", "E118", "E119", "E120", "E121", "E126", "E128", "E129", "E130", "E131", "E136", "E138", "E139", "E140", "E141", "E146", "E148", "E149"),
    "diabwc": ("E102", "E103", "E104", "E105", "E107", "E112", "E113", "E114", "E115", "E117", "E122", "E123", "E124", "E125", "E127", "E132", "E133", "E134", "E135", "E137", "E142", "E143", "E144", "E145", "E147"),
    "hp": ("G041", "G114", "G801", "G802", "G81", "G82", "G830", "G831", "G832", "G833", "G834", "G839"),
    "rend": ("I120", "I131", "N032", "N033", "N034", "N035", "N036", "N037", "N052", "N053", "N054", "N055", "N056", "N057", "N18", "N19", "N250", "Z490", "Z491", "Z492", "Z940", "Z992"),
    "msld": ("I850", "I859", "I864", "I982", "K704", "K711", "K721", "K729", "K765", "K766", "K767"),
    "aids": ("B20", "B21", "B22", "B24"),
}
CCI_WEIGHTS = {
    "ami": 1, "chf": 1, "pvd": 1, "cevd": 1, "dementia": 1,
    "copd": 1, "rheumd": 1, "pud": 1, "mld": 1, "diab": 1,
    "diabwc": 2, "hp": 2, "rend": 2, "msld": 3, "aids": 6,
}


def sql_string(value: str) -> str:
    return "'" + value.replace("'", "''") + "'"


def component_condition(prefixes: tuple[str, ...], code: str = "code") -> str:
    return "(" + " OR ".join(
        f"starts_with({code}, {sql_string(prefix)})" for prefix in prefixes
    ) + ")"


def score_from_flags(flags: dict[str, bool]) -> int:
    """Original Charlson score, without cancer or age points, with hierarchy."""
    adjusted = {key: bool(flags.get(key, False)) for key in CCI_COMPONENTS}
    if adjusted["diabwc"]:
        adjusted["diab"] = False
    if adjusted["msld"]:
        adjusted["mld"] = False
    return sum(CCI_WEIGHTS[key] for key, present in adjusted.items() if present)


def weighted_mean_sd(rows: tuple[float, float, float, float]) -> tuple[float, float]:
    sum_w, sum_wx, sum_wx2, _ = rows
    mean = sum_wx / sum_w
    variance = max(sum_wx2 / sum_w - mean * mean, 0.0)
    return mean, math.sqrt(variance)


def format_p_value(value: float) -> str:
    return "<0.001" if value < 0.001 else f"{value:.3f}"


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main() -> dict[str, Any]:
    if not COHORT_DATABASE.exists():
        raise RuntimeError("Run Phase 1–2 first to create the cached HM cohort.")
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    connection = duckdb.connect(str(COHORT_DATABASE), read_only=True)

    flag_sql = []
    for component, prefixes in CCI_COMPONENTS.items():
        condition = component_condition(prefixes)
        flag_sql.append(
            "list_contains(list_transform(diagnosis_codes, "
            f"code -> {condition}), TRUE) AS {component}"
        )
    score_terms = []
    for component, weight in CCI_WEIGHTS.items():
        if component == "diab":
            present = "diab AND NOT diabwc"
        elif component == "mld":
            present = "mld AND NOT msld"
        else:
            present = component
        score_terms.append(f"CASE WHEN {present} THEN {weight} ELSE 0 END")
    connection.execute(
        "CREATE TEMP TABLE cci_flags AS SELECT sepsis, DISCWT, LOS, DIED, "
        + ", ".join(flag_sql)
        + " FROM hm_cohort"
    )
    connection.execute(
        "CREATE TEMP TABLE outcomes AS SELECT sepsis, DISCWT, LOS, DIED, "
        + " + ".join(score_terms)
        + " AS cci_excluding_cancer FROM cci_flags"
    )

    totals = {
        bool(row[0]): {"unweighted": int(row[1]), "weighted": float(row[2])}
        for row in connection.execute(
            "SELECT sepsis, count(*), sum(DISCWT) FROM outcomes GROUP BY sepsis"
        ).fetchall()
    }

    def continuous(variable: str) -> dict[bool, dict[str, float]]:
        output: dict[bool, dict[str, float]] = {}
        query = f"""
            SELECT sepsis, count({variable}), sum(DISCWT),
                   sum(DISCWT * {variable}), sum(DISCWT * {variable} * {variable}),
                   avg({variable}), stddev_samp({variable})
            FROM outcomes WHERE {variable} IS NOT NULL GROUP BY sepsis
        """
        for status, n, sum_w, sum_wx, sum_wx2, raw_mean, raw_sd in connection.execute(query).fetchall():
            mean, sd = weighted_mean_sd((sum_w, sum_wx, sum_wx2, n))
            output[bool(status)] = {
                "n": int(n), "mean": mean, "sd": sd,
                "raw_mean": float(raw_mean), "raw_sd": float(raw_sd),
            }
        return output

    cci_stats = continuous("cci_excluding_cancer")
    los_stats = continuous("LOS")

    def welch(stats: dict[bool, dict[str, float]]) -> float:
        return float(ttest_ind_from_stats(
            stats[True]["raw_mean"], stats[True]["raw_sd"], stats[True]["n"],
            stats[False]["raw_mean"], stats[False]["raw_sd"], stats[False]["n"],
            equal_var=False,
        ).pvalue)

    cci_p = welch(cci_stats)
    los_p = welch(los_stats)

    cci_category_rows: list[dict[str, Any]] = []
    contingency: list[list[int]] = []
    for label, predicate in [("0", "= 0"), ("1–2", "BETWEEN 1 AND 2"), ("≥3", ">= 3")]:
        by_status = {
            bool(row[0]): (int(row[1]), float(row[2] or 0))
            for row in connection.execute(
                f"SELECT sepsis, count(*) FILTER (WHERE cci_excluding_cancer {predicate}), "
                f"sum(DISCWT) FILTER (WHERE cci_excluding_cancer {predicate}) "
                "FROM outcomes GROUP BY sepsis"
            ).fetchall()
        }
        contingency.append([by_status[False][0], by_status[True][0]])
        cci_category_rows.append({
            "level": label,
            "no_sepsis_unweighted_n": by_status[False][0],
            "no_sepsis_weighted_n": round(by_status[False][1]),
            "no_sepsis_weighted_percent": round(100 * by_status[False][1] / totals[False]["weighted"], 2),
            "sepsis_unweighted_n": by_status[True][0],
            "sepsis_weighted_n": round(by_status[True][1]),
            "sepsis_weighted_percent": round(100 * by_status[True][1] / totals[True]["weighted"], 2),
        })
    cci_category_p = float(chi2_contingency(contingency)[1])

    mortality = {
        bool(row[0]): {
            "valid_n": int(row[1]), "events": int(row[2]),
            "valid_weighted": float(row[3]), "event_weighted": float(row[4] or 0),
        }
        for row in connection.execute("""
            SELECT sepsis, count(DIED), count(*) FILTER (WHERE DIED = 1),
                   sum(DISCWT) FILTER (WHERE DIED IS NOT NULL),
                   sum(DISCWT) FILTER (WHERE DIED = 1)
            FROM outcomes GROUP BY sepsis
        """).fetchall()
    }
    mortality_table = [
        [mortality[False]["events"], mortality[False]["valid_n"] - mortality[False]["events"]],
        [mortality[True]["events"], mortality[True]["valid_n"] - mortality[True]["events"]],
    ]
    mortality_p = float(chi2_contingency(mortality_table)[1])

    summary_rows = [
        {
            "outcome": "Charlson Comorbidity Index excluding cancer",
            "level": "Mean (SD)",
            "no_sepsis": f'{cci_stats[False]["mean"]:.2f} ({cci_stats[False]["sd"]:.2f})',
            "sepsis": f'{cci_stats[True]["mean"]:.2f} ({cci_stats[True]["sd"]:.2f})',
            "p_value": format_p_value(cci_p), "test": "Welch t-test",
        }
    ]
    for index, row in enumerate(cci_category_rows):
        summary_rows.append({
            "outcome": "Charlson category" if index == 0 else "",
            "level": row["level"],
            "no_sepsis": f'{row["no_sepsis_weighted_n"]:,} ({row["no_sepsis_weighted_percent"]:.2f}%)',
            "sepsis": f'{row["sepsis_weighted_n"]:,} ({row["sepsis_weighted_percent"]:.2f}%)',
            "p_value": format_p_value(cci_category_p) if index == 0 else "",
            "test": "Pearson chi-square" if index == 0 else "",
        })
    summary_rows.extend([
        {
            "outcome": "Length of stay, days", "level": "Mean (SD)",
            "no_sepsis": f'{los_stats[False]["mean"]:.2f} ({los_stats[False]["sd"]:.2f})',
            "sepsis": f'{los_stats[True]["mean"]:.2f} ({los_stats[True]["sd"]:.2f})',
            "p_value": format_p_value(los_p), "test": "Welch t-test",
        },
        {
            "outcome": "In-hospital mortality", "level": "Died",
            "no_sepsis": f'{round(mortality[False]["event_weighted"]):,} ({100*mortality[False]["event_weighted"]/mortality[False]["valid_weighted"]:.2f}%)',
            "sepsis": f'{round(mortality[True]["event_weighted"]):,} ({100*mortality[True]["event_weighted"]/mortality[True]["valid_weighted"]:.2f}%)',
            "p_value": format_p_value(mortality_p), "test": "Pearson chi-square",
        },
    ])
    unweighted_rows = [
        {
            "outcome": "Charlson Comorbidity Index excluding cancer",
            "level": "Mean (SD)",
            "no_sepsis": f'{cci_stats[False]["raw_mean"]:.2f} ({cci_stats[False]["raw_sd"]:.2f})',
            "sepsis": f'{cci_stats[True]["raw_mean"]:.2f} ({cci_stats[True]["raw_sd"]:.2f})',
            "p_value": format_p_value(cci_p), "test": "Welch t-test",
        }
    ]
    for index, row in enumerate(cci_category_rows):
        unweighted_rows.append({
            "outcome": "Charlson category" if index == 0 else "",
            "level": row["level"],
            "no_sepsis": f'{row["no_sepsis_unweighted_n"]:,} ({100 * row["no_sepsis_unweighted_n"] / totals[False]["unweighted"]:.2f}%)',
            "sepsis": f'{row["sepsis_unweighted_n"]:,} ({100 * row["sepsis_unweighted_n"] / totals[True]["unweighted"]:.2f}%)',
            "p_value": format_p_value(cci_category_p) if index == 0 else "",
            "test": "Pearson chi-square" if index == 0 else "",
        })
    unweighted_rows.extend([
        {
            "outcome": "Length of stay, days", "level": "Mean (SD)",
            "no_sepsis": f'{los_stats[False]["raw_mean"]:.2f} ({los_stats[False]["raw_sd"]:.2f})',
            "sepsis": f'{los_stats[True]["raw_mean"]:.2f} ({los_stats[True]["raw_sd"]:.2f})',
            "p_value": format_p_value(los_p), "test": "Welch t-test",
        },
        {
            "outcome": "In-hospital mortality", "level": "Died",
            "no_sepsis": f'{mortality[False]["events"]:,} ({100 * mortality[False]["events"] / mortality[False]["valid_n"]:.2f}%)',
            "sepsis": f'{mortality[True]["events"]:,} ({100 * mortality[True]["events"] / mortality[True]["valid_n"]:.2f}%)',
            "p_value": format_p_value(mortality_p), "test": "Pearson chi-square",
        },
    ])
    write_csv(OUTPUT_DIR / "cci_categories_by_sepsis.csv", cci_category_rows)
    write_csv(OUTPUT_DIR / "clinical_outcomes_by_sepsis.csv", summary_rows)
    write_csv(OUTPUT_DIR / "clinical_outcomes_by_sepsis_unweighted.csv", unweighted_rows)
    summary = {
        "unit": "NIS inpatient discharge, not unique patient",
        "cci_definition": "Quan ICD-10 mapping; original Charlson weights; cancer and metastatic cancer excluded; no age points; hierarchy applied.",
        "inference_note": "Descriptive values use DISCWT. P-values use sampled discharge counts and do not account for hospital clustering.",
        "no_sepsis_unweighted_n": totals[False]["unweighted"],
        "sepsis_unweighted_n": totals[True]["unweighted"],
        "los_missing": {"no_sepsis": totals[False]["unweighted"] - los_stats[False]["n"], "sepsis": totals[True]["unweighted"] - los_stats[True]["n"]},
        "died_missing": {"no_sepsis": totals[False]["unweighted"] - mortality[False]["valid_n"], "sepsis": totals[True]["unweighted"] - mortality[True]["valid_n"]},
        "table": summary_rows,
        "unweighted_table": unweighted_rows,
    }
    (OUTPUT_DIR / "phase_5_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    connection.close()
    print(json.dumps(summary, indent=2))
    return summary


if __name__ == "__main__":
    main()
