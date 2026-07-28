"""Rerun the validated study pipeline with septic shock (R65.21) as exposure."""

from __future__ import annotations

import json
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


REPO_ROOT = Path(__file__).resolve().parents[1]
SOURCE_DATABASE = REPO_ROOT / "outputs/phase_1_2/hm_cohort.duckdb"
OUTPUT_ROOT = REPO_ROOT / "outputs/septic_shock"
EXPOSURE_DATABASE = OUTPUT_ROOT / "hm_cohort_septic_shock.duckdb"
PIPELINE_VERSION = "1.0.0"


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
    }
    summary_path = OUTPUT_ROOT / "septic_shock_pipeline_summary.json"
    summary_path.write_text(json.dumps(results, indent=2) + "\n")
    print(json.dumps({
        "exposure": exposure,
        "primary_adjusted": results["phase_8"]["primary_results"],
        "interaction": results["phase_9"]["joint_interaction_test"],
        "decedent_interaction": results["phase_10"]["command_19b_interaction_test"],
        "overall_eapc": results["phase_11"]["command_20a_eapc"],
    }, indent=2))
    return results


if __name__ == "__main__":
    main()
