"""Rerun the validated study pipeline with septic shock (R65.21) as exposure."""

from __future__ import annotations

import json
import csv
from pathlib import Path
from typing import Any

import duckdb

from src import phase_3
from src import phase_4_complications
from src import phase_5_cci_los_mortality
from src import phase_6_palliative_care
from src import phase_7_subtype_palliative_care
from src import phase_8_primary_adjusted
from src import phase_9_interaction
from src import phase_10_decedents
from src import phase_11_trends
from src import bmt_shock_interaction
from src import stratified_covariate_effects


REPO_ROOT = Path(__file__).resolve().parents[1]
SOURCE_DATABASE = REPO_ROOT / "outputs/phase_1_2/hm_cohort.duckdb"
OUTPUT_ROOT = REPO_ROOT / "outputs/septic_shock"
EXPOSURE_DATABASE = OUTPUT_ROOT / "hm_cohort_septic_shock.duckdb"
PIPELINE_VERSION = "1.0.0"


def _report_label(value: object) -> str:
    """Replace legacy engine labels and escape Markdown table delimiters."""
    text = str(value)
    replacements = {
        "without documented sepsis": "without septic shock",
        "with documented sepsis": "with septic shock",
        "No documented sepsis": "No septic shock",
        "Documented sepsis": "Septic shock",
        "No sepsis": "No septic shock",
        "Sepsis": "Septic shock",
        "documented sepsis": "documented septic shock",
    }
    for old, new in replacements.items():
        text = text.replace(old, new)
    return text.replace("|", "\\|").replace("\n", " ")


def _csv_markdown(path: Path) -> str:
    """Render a generated CSV as a portable Markdown table."""
    with path.open(newline="") as handle:
        rows = list(csv.reader(handle))
    if not rows:
        return "_No rows were produced._"
    header, body = rows[0], rows[1:]
    friendly_header = []
    for value in header:
        value = value.replace("no_sepsis", "no_septic_shock").replace("sepsis", "septic_shock")
        friendly_header.append(value.replace("_", " ").capitalize())
    lines = [
        "| " + " | ".join(_report_label(value) for value in friendly_header) + " |",
        "| " + " | ".join("---" for _ in header) + " |",
    ]
    lines.extend("| " + " | ".join(_report_label(value) for value in row) + " |" for row in body)
    return "\n".join(lines)


def write_markdown_report(results: dict[str, Any]) -> Path:
    """Write a human-readable companion to the machine-readable JSON summary."""
    x = results["exposure"]
    sections = [
        ("1. Cohort definition", None),
        ("2. Baseline characteristics", OUTPUT_ROOT / "phase_3/table_1_baseline_characteristics.csv"),
        ("3. Clinical outcomes — weighted", OUTPUT_ROOT / "phase_5/clinical_outcomes_by_sepsis.csv"),
        ("3. Clinical outcomes — unweighted", OUTPUT_ROOT / "phase_5/clinical_outcomes_by_sepsis_unweighted.csv"),
        ("6. Complications", OUTPUT_ROOT / "phase_4/complications_by_sepsis.csv"),
        ("7. Overall documented inpatient palliative-care utilization", OUTPUT_ROOT / "phase_6/palliative_care_prevalence.csv"),
        ("8. Annual palliative-care utilization by septic-shock status", OUTPUT_ROOT / "phase_6/annual_palliative_care_by_sepsis.csv"),
        ("9. Palliative-care utilization by HM subtype — all HM", OUTPUT_ROOT / "phase_7/subtype_palliative_care_all_hm.csv"),
        ("10. Palliative-care utilization by HM subtype — no septic shock", OUTPUT_ROOT / "phase_7/subtype_palliative_care_no_sepsis.csv"),
        ("11. Palliative-care utilization by HM subtype — septic shock", OUTPUT_ROOT / "phase_7/subtype_palliative_care_sepsis.csv"),
        ("12. Primary adjusted association", OUTPUT_ROOT / "phase_8/primary_adjusted_results.csv"),
        ("12. Adjusted marginal probabilities", OUTPUT_ROOT / "phase_8/adjusted_probabilities.csv"),
        ("13. Decedent analysis", OUTPUT_ROOT / "phase_10/decedent_palliative_care_by_sepsis.csv"),
        ("14. Adjusted decedent interaction", OUTPUT_ROOT / "phase_10/adjusted_decedent_interaction.csv"),
        ("15. Septic-shock decedents by HM subtype — unadjusted", OUTPUT_ROOT / "phase_10/sepsis_decedent_subtype_unadjusted.csv"),
        ("15. Septic-shock decedents by HM subtype — adjusted", OUTPUT_ROOT / "phase_10/sepsis_decedent_subtype_adjusted.csv"),
        ("16. Overall annual trends", OUTPUT_ROOT / "phase_11/annual_overall_palliative_care.csv"),
        ("17. Trends by septic-shock status", OUTPUT_ROOT / "phase_11/annual_palliative_care_by_sepsis.csv"),
        ("18. Annual trends among septic-shock admissions by HM subtype", OUTPUT_ROOT / "phase_11/annual_sepsis_palliative_care_by_subtype.csv"),
        ("18. Subtype-specific EAPC", OUTPUT_ROOT / "phase_11/sepsis_subtype_trend_tests.csv"),
        ("20. Subtype-specific adjusted marginal effects", OUTPUT_ROOT / "phase_9/subtype_adjusted_probabilities.csv"),
        ("21. BMT/HSCT four-group descriptive analysis", OUTPUT_ROOT / "bmt_interaction/four_group_descriptive.csv"),
        ("21. BMT/HSCT conditional septic-shock odds ratios", OUTPUT_ROOT / "bmt_interaction/conditional_odds_ratios.csv"),
        ("21. BMT/HSCT adjusted four-group probabilities", OUTPUT_ROOT / "bmt_interaction/adjusted_four_group_probabilities.csv"),
        ("21. BMT/HSCT adjusted probability differences", OUTPUT_ROOT / "bmt_interaction/adjusted_probability_differences.csv"),
        ("22. Covariate effects stratified by septic-shock status", OUTPUT_ROOT / "stratified_covariate_effects/stratified_covariate_adjusted_odds_ratios.csv"),
    ]
    lines = [
        "# Septic Shock and Documented Inpatient Palliative-Care Use in Hematologic Malignancy Hospitalizations",
        "",
        "## Analysis definition",
        "",
        "Primary exposure: documented septic shock, defined as exact normalized ICD-10-CM code `R65.21` (`R6521`) in any diagnosis position. The outcome is documented inpatient palliative-care use (`Z51.5`). Estimates use `DISCWT`; variance estimation accounts for year-specific `NIS_STRATUM` without `HOSP_NIS`.",
        "",
        "This Markdown file is intended for human review. `septic_shock_pipeline_summary.json` remains the machine-readable audit record.",
        "",
    ]
    for heading, path in sections:
        lines.extend([f"## {heading}", ""])
        if path is None:
            lines.extend([
                "| Cohort | Unweighted hospitalizations | Weighted hospitalizations |",
                "| --- | ---: | ---: |",
                f'| No septic shock | {x["no_shock_unweighted"]:,} | {x["no_shock_weighted"]:,} |',
                f'| Septic shock (`R65.21`) | {x["shock_unweighted"]:,} | {x["shock_weighted"]:,} |',
                f'| Total | {x["total_unweighted"]:,} | {x["total_weighted"]:,} |',
            ])
        else:
            lines.append(_csv_markdown(path))
        lines.append("")
    report_path = OUTPUT_ROOT / "septic_shock_report.md"
    report_path.write_text("\n".join(lines), encoding="utf-8")
    return report_path


def build_exposure_database() -> dict[str, Any]:
    """Create a compact cached cohort where `sepsis` means exact R6521."""
    if not SOURCE_DATABASE.exists():
        raise RuntimeError("The Phase 1–2 HM cohort database is missing.")
    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
    source_mtime = SOURCE_DATABASE.stat().st_mtime_ns
    metadata_path = OUTPUT_ROOT / "exposure_cache_metadata.json"
    rebuild = True
    if EXPOSURE_DATABASE.exists() and metadata_path.exists():
        metadata = json.loads(metadata_path.read_text())
        rebuild = metadata.get("source_mtime_ns") != source_mtime or metadata.get("pipeline_version") != PIPELINE_VERSION
    if rebuild:
        if EXPOSURE_DATABASE.exists():
            EXPOSURE_DATABASE.unlink()
        connection = duckdb.connect(str(EXPOSURE_DATABASE))
        source_path_sql = str(SOURCE_DATABASE).replace("'", "''")
        connection.execute(f"ATTACH '{source_path_sql}' AS source_db (READ_ONLY)")
        connection.execute("""
            CREATE TABLE hm_cohort AS
            SELECT * EXCLUDE (sepsis),
                   sepsis AS original_documented_sepsis,
                   list_contains(diagnosis_codes, 'R6521') AS sepsis
            FROM source_db.hm_cohort
        """)
        connection.execute("DETACH source_db")
        connection.close()
        metadata_path.write_text(json.dumps({
            "pipeline_version": PIPELINE_VERSION,
            "source_mtime_ns": source_mtime,
            "exposure": "Documented septic shock",
            "normalized_exact_code": "R6521",
        }, indent=2) + "\n")
    connection = duckdb.connect(str(EXPOSURE_DATABASE), read_only=True)
    counts = connection.execute("""
        SELECT count(*)::BIGINT,
               count(*) FILTER (WHERE sepsis)::BIGINT,
               count(*) FILTER (WHERE NOT sepsis)::BIGINT,
               round(sum(DISCWT), 0)::BIGINT,
               round(sum(DISCWT) FILTER (WHERE sepsis), 0)::BIGINT,
               round(sum(DISCWT) FILTER (WHERE NOT sepsis), 0)::BIGINT
        FROM hm_cohort
    """).fetchone()
    validation = connection.execute("""
        SELECT count(*) FILTER (
            WHERE sepsis != list_contains(diagnosis_codes, 'R6521')
        ) FROM hm_cohort
    """).fetchone()[0]
    connection.close()
    return {
        "rebuilt": rebuild,
        "total_unweighted": counts[0],
        "shock_unweighted": counts[1],
        "no_shock_unweighted": counts[2],
        "total_weighted": counts[3],
        "shock_weighted": counts[4],
        "no_shock_weighted": counts[5],
        "validation_discordances": validation,
    }


def configure_module(module: Any, phase_name: str) -> None:
    module.COHORT_DATABASE = EXPOSURE_DATABASE
    module.OUTPUT_DIR = OUTPUT_ROOT / phase_name


def main() -> dict[str, Any]:
    exposure = build_exposure_database()
    modules = [
        (phase_3, "phase_3"),
        (phase_4_complications, "phase_4"),
        (phase_5_cci_los_mortality, "phase_5"),
        (phase_6_palliative_care, "phase_6"),
        (phase_7_subtype_palliative_care, "phase_7"),
        (phase_8_primary_adjusted, "phase_8"),
        (phase_9_interaction, "phase_9"),
        (phase_10_decedents, "phase_10"),
        (phase_11_trends, "phase_11"),
    ]
    for module, phase_name in modules:
        configure_module(module, phase_name)
    # Modules that imported COHORT_DATABASE by value need their globals patched.
    phase_9_interaction.COHORT_DATABASE = EXPOSURE_DATABASE
    phase_10_decedents.COHORT_DATABASE = EXPOSURE_DATABASE
    phase_11_trends.COHORT_DATABASE = EXPOSURE_DATABASE

    results = {
        "exposure": exposure,
        "phase_3": phase_3.main(),
        "phase_4": phase_4_complications.main(),
        "phase_5": phase_5_cci_los_mortality.main(),
        "phase_6": phase_6_palliative_care.main(),
        "phase_7": phase_7_subtype_palliative_care.main(),
        "phase_8": phase_8_primary_adjusted.main(),
        "phase_9": phase_9_interaction.main(),
        "phase_10": phase_10_decedents.main(),
        "phase_11": phase_11_trends.main(),
        "bmt_shock_interaction": bmt_shock_interaction.main(),
        "stratified_covariate_effects": stratified_covariate_effects.main(),
    }
    summary_path = OUTPUT_ROOT / "septic_shock_pipeline_summary.json"
    summary_path.write_text(json.dumps(results, indent=2) + "\n")
    report_path = write_markdown_report(results)
    print(json.dumps({
        "exposure": exposure,
        "primary_adjusted": results["phase_8"]["primary_results"],
        "interaction": results["phase_9"]["joint_interaction_test"],
        "decedent_interaction": results["phase_10"]["command_19b_interaction_test"],
        "overall_eapc": results["phase_11"]["command_20a_eapc"],
        "human_readable_report": str(report_path),
    }, indent=2))
    return results


if __name__ == "__main__":
    main()
