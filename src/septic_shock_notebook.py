"""Presentation helpers for the septic-shock master notebook."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
from IPython.display import Markdown, display


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs/septic_shock"


def _csv(phase: int, name: str) -> pd.DataFrame:
    return pd.read_csv(OUT / f"phase_{phase}" / name, dtype=str).fillna("")


def _rename_shock(frame: pd.DataFrame) -> pd.DataFrame:
    result = frame.copy()
    result.columns = [c.replace("no_sepsis", "no_shock").replace("sepsis", "septic_shock") for c in result.columns]
    for column in result.columns:
        if result[column].dtype == object:
            result[column] = (result[column].str.replace("without documented sepsis", "without septic shock", case=False)
                              .str.replace("with documented sepsis", "with septic shock", case=False)
                              .str.replace("No sepsis", "No septic shock", case=False)
                              .str.replace("Sepsis", "Septic shock", case=False))
    return result


def _show(frame: pd.DataFrame, note: str = "") -> None:
    def clean(value: object) -> str:
        return str(value).replace("|", "\\|").replace("\n", " ")
    header = "| " + " | ".join(clean(column) for column in frame.columns) + " |"
    rule = "| " + " | ".join("---" for _ in frame.columns) + " |"
    rows = ["| " + " | ".join(clean(value) for value in row) + " |" for row in frame.itertuples(index=False, name=None)]
    display(Markdown("\n".join([header, rule, *rows])))
    if note:
        display(Markdown(f"*{note}*"))


def _total_row(frame: pd.DataFrame, label_column: str, values: dict | None = None) -> pd.DataFrame:
    if frame.astype(str).apply(lambda col: col.str.lower().eq("total")).any().any():
        return frame
    row = {column: "—" for column in frame.columns}
    row[label_column] = "Total"
    row.update(values or {})
    return pd.concat([frame, pd.DataFrame([row])], ignore_index=True)


def rerun_all() -> dict:
    """Rerun and cache all analyses before presenting results."""
    from src.septic_shock_pipeline import main
    return main()


def show_step(step: int) -> None:
    """Display one copy/paste-friendly requested result section."""
    if step == 1:
        x = json.loads((OUT / "septic_shock_pipeline_summary.json").read_text())["exposure"]
        table = pd.DataFrame([
            ["No septic shock", x["no_shock_unweighted"], x["no_shock_weighted"]],
            ["Septic shock (R65.21)", x["shock_unweighted"], x["shock_weighted"]],
            ["Total", x["total_unweighted"], x["total_weighted"]],
        ], columns=["Cohort", "Unweighted hospitalizations", "Weighted hospitalizations"])
        _show(table, "Exposure is an exact normalized match to R6521 in any diagnosis position; one row represents one hospitalization.")
    elif step == 2:
        t = _rename_shock(_csv(3, "table_1_baseline_characteristics.csv"))
        keep = ["Age, years", "Age group", "Sex", "Race/ethnicity", "Primary payer", "Income quartile",
                "Hospital region", "Hospital location/teaching", "Hospital bed size", "Hospital ownership"]
        t = t[t.characteristic.isin(keep)]
        t = _total_row(t, "characteristic", {"level": "All adult HM hospitalizations",
            "no_shock_unweighted_n": "951721", "no_shock_weighted_n": "4758603",
            "septic_shock_unweighted_n": "43271", "septic_shock_weighted_n": "216355"})
        _show(t, "Percentages and means are DISCWT-weighted. P-values use NIS_STRATUM-based variance estimation; SMD describes balance independently of sample size.")
    elif step == 3:
        w = _rename_shock(_csv(5, "clinical_outcomes_by_sepsis.csv"))
        u = _rename_shock(_csv(5, "clinical_outcomes_by_sepsis_unweighted.csv"))
        w = w.rename(columns={"no_shock": "no_shock_weighted", "septic_shock": "septic_shock_weighted"})
        u = u.rename(columns={"no_shock": "no_shock_unweighted", "septic_shock": "septic_shock_unweighted"})
        t = w.merge(u[["outcome", "level", "no_shock_unweighted", "septic_shock_unweighted"]], on=["outcome", "level"])
        t = _total_row(t, "outcome", {"level": "Hospitalizations", "no_shock_weighted": "4,758,603",
            "septic_shock_weighted": "216,355", "no_shock_unweighted": "951,721", "septic_shock_unweighted": "43,271"})
        _show(t, "Means/SDs and percentages are weighted; event counts are shown in both weighted and sampled forms.")
    elif step == 4:
        t = _rename_shock(_csv(3, "table_1_baseline_characteristics.csv"))
        t = t[t.characteristic.eq("HM subtype")]
        t = _total_row(t, "characteristic", {"level": "All mutually exclusive subtypes", "no_shock_weighted_percent_or_mean_sd": "100.00%",
            "septic_shock_weighted_percent_or_mean_sd": "100.00%", "no_shock_unweighted_n": "951721",
            "septic_shock_unweighted_n": "43271", "no_shock_weighted_n": "4758603", "septic_shock_weighted_n": "216355"})
        _show(t, "Subtypes are mutually exclusive and assigned from the first qualifying HM diagnosis field; overlapping groups remain separately flagged for sensitivity analysis.")
    elif step == 5:
        t = _rename_shock(_csv(3, "table_1_baseline_characteristics.csv"))
        t = t[t.characteristic.eq("Admission year")]
        t = _total_row(t, "characteristic", {"level": "2016–2022", "no_shock_weighted_percent_or_mean_sd": "100.00%",
            "septic_shock_weighted_percent_or_mean_sd": "100.00%", "no_shock_unweighted_n": "951721",
            "septic_shock_unweighted_n": "43271", "no_shock_weighted_n": "4758603", "septic_shock_weighted_n": "216355"})
        _show(t)
    elif step == 6:
        t = _rename_shock(_csv(4, "complications_by_sepsis.csv"))
        wanted = ["Acute respiratory failure", "Acute kidney injury", "Pneumonia", "Pulmonary embolism",
                  "Acute lower-extremity DVT", "Tumor lysis syndrome", "Disseminated intravascular coagulation",
                  "Any specified cardiac arrhythmia"]
        t = t[t.complication.isin(wanted)]
        t = _total_row(t, "complication", {"no_shock_unweighted_n": "951721", "no_shock_weighted_n": "4758603",
            "septic_shock_unweighted_n": "43271", "septic_shock_weighted_n": "216355"})
        _show(t, "Arrhythmias are the union of I47, I48 excluding chronic I482*, and I49, so a hospitalization is counted once.")
    elif step == 7:
        t = _rename_shock(_csv(6, "palliative_care_prevalence.csv"))
        summary = json.loads((OUT / "phase_6/phase_6_summary.json").read_text())["comparison"]
        t["absolute_difference_pp"] = ""
        t["crude_odds_ratio_95_ci"] = ""
        t["p_value"] = ""
        t.loc[t.index[-1], "absolute_difference_pp"] = f'{summary["absolute_difference_percentage_points"]} ({summary["difference_ci_95_lower"]}–{summary["difference_ci_95_upper"]})'
        t.loc[t.index[-1], "crude_odds_ratio_95_ci"] = f'{summary["crude_weighted_odds_ratio"]} ({summary["odds_ratio_ci_95_lower"]}–{summary["odds_ratio_ci_95_upper"]})'
        t.loc[t.index[-1], "p_value"] = summary["strata_adjusted_p_value"]
        t = _total_row(t, "cohort", {"unweighted_sample_n": "994992", "weighted_hospitalizations": "4974958"})
        _show(t, "Outcome is documented inpatient palliative-care use (Z51.5), not necessarily a consultation.")
    elif step == 8:
        t = _rename_shock(_csv(6, "annual_palliative_care_by_sepsis.csv"))
        t = _total_row(t, "year", {"septic_shock_status": "All years; see Step 7 totals"})
        _show(t)
    elif step in (9, 10, 11):
        name = {9: "subtype_palliative_care_all_hm.csv", 10: "subtype_palliative_care_no_sepsis.csv", 11: "subtype_palliative_care_sepsis.csv"}[step]
        _show(_rename_shock(_csv(7, name)), "The final row contains cohort totals and the survey-adjusted overall subtype p-value.")
    elif step == 12:
        a = _rename_shock(_csv(8, "primary_adjusted_results.csv"))
        p = _rename_shock(_csv(8, "adjusted_probabilities.csv"))
        _show(a)
        _show(p, "Marginal probabilities average over the observed covariate distribution.")
    elif step == 13:
        _show(_rename_shock(_csv(10, "decedent_palliative_care_by_sepsis.csv")))
    elif step == 14:
        _show(_rename_shock(_csv(10, "adjusted_decedent_interaction.csv")), "The total row reports the joint Wald interaction p-value.")
    elif step == 15:
        _show(_rename_shock(_csv(10, "sepsis_decedent_subtype_unadjusted.csv")))
        _show(_rename_shock(_csv(10, "sepsis_decedent_subtype_adjusted.csv")), "Adjusted table follows the weighted descriptive table; each includes an overall subtype p-value.")
    elif step == 16:
        _show(_csv(11, "annual_overall_palliative_care.csv"))
        s = json.loads((OUT / "phase_11/phase_11_summary.json").read_text())["command_20a_eapc"]
        _show(pd.DataFrame([["Total / overall EAPC", s["eapc_percent"], f'{s["ci_95_lower_percent"]}–{s["ci_95_upper_percent"]}', s["p_value"]]], columns=["Measure", "EAPC, %", "95% CI", "p-value"]))
    elif step == 17:
        _show(_rename_shock(_csv(11, "annual_palliative_care_by_sepsis.csv")))
        s = json.loads((OUT / "phase_11/phase_11_summary.json").read_text())["command_20b"]
        t = pd.DataFrame([
            ["No septic shock", s["no_sepsis_eapc"]["eapc_percent"], f'{s["no_sepsis_eapc"]["ci_95_lower_percent"]}–{s["no_sepsis_eapc"]["ci_95_upper_percent"]}', s["no_sepsis_eapc"]["p_value"], ""],
            ["Septic shock", s["sepsis_eapc"]["eapc_percent"], f'{s["sepsis_eapc"]["ci_95_lower_percent"]}–{s["sepsis_eapc"]["ci_95_upper_percent"]}', s["sepsis_eapc"]["p_value"], ""],
            ["Total / year × shock test", "—", "—", "—", s["year_by_sepsis_interaction_p_value"]],
        ], columns=["Cohort", "EAPC, %", "95% CI", "trend p-value", "interaction p-value"])
        _show(t)
    elif step == 18:
        _show(_rename_shock(_csv(11, "annual_sepsis_palliative_care_by_subtype.csv")))
        _show(_rename_shock(_csv(11, "sepsis_subtype_trend_tests.csv")), "No subtype-year cell was suppressed; threshold was 10 or fewer unweighted palliative-care events.")
    elif step == 19:
        s = json.loads((OUT / "phase_9/phase_9_summary.json").read_text())["joint_interaction_test"]
        _show(pd.DataFrame([["Septic shock × HM subtype", s["wald_chi_square"], s["degrees_of_freedom"], s["overall_interaction_p_value"]], ["Total", "—", "—", s["overall_interaction_p_value"]]], columns=["Test", "Wald chi-square", "df", "p-value"]))
    elif step == 20:
        _show(_rename_shock(_csv(9, "subtype_adjusted_probabilities.csv")), "The total row reports the overall interaction p-value.")
    elif step == 21:
        directory = OUT / "bmt_interaction"
        def read(name: str) -> pd.DataFrame:
            return pd.read_csv(directory / name, dtype=str).fillna("")
        display(Markdown("### Four observed groups"))
        _show(read("four_group_descriptive.csv"))
        display(Markdown("### Conditional adjusted odds ratios"))
        _show(read("conditional_odds_ratios.csv"), "The final row reports the 1-df multiplicative interaction Wald p-value.")
        display(Markdown("### Adjusted marginal probabilities"))
        _show(read("adjusted_four_group_probabilities.csv"))
        display(Markdown("### Adjusted probability differences"))
        _show(read("adjusted_probability_differences.csv"), "The difference-in-differences is an additive-scale interaction; it can differ from the logistic model's multiplicative interaction test.")
        display(Markdown("*BMT/HSCT status is diagnosis-based (`Z94.81` or `Z94.84` in any diagnosis position). Actual inpatient HSCT procedures cannot be added because procedure-code fields are unavailable in the source extract.*"))
    else:
        raise ValueError("Step must be from 1 through 21.")
