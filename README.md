# IPC-MPC Study

Reproducible analysis workspace for the NIS 2016–2022 dataset. Source data stays outside this Git repository and is queried directly as a partitioned Parquet dataset.

## Complete study notebook

Use `notebooks/IPC_MPC_Septic_Shock_Master.ipynb` as the primary single-notebook reference. On macOS, double-click `launch_septic_shock_notebook.command`, then choose **Run → Run All Cells**. It runs the 20 requested steps using exact normalized ICD-10-CM code `R6521` as the primary exposure. Results are isolated under `outputs/septic_shock/`, preserving the earlier A41 sepsis analysis and `notebooks/IPC_MPC_Study_Master.ipynb` for auditability.

## Local setup

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
cp .env.example .env
```

Edit `.env` so `NIS_DATASET_PATH` points to the local dataset directory. Then validate access and print metadata only:

```bash
python scripts/inspect_dataset.py
```

The inspection script reports file count, total size, schema, and row count. It deliberately does not print patient-level rows.

## Smoke-test analysis

Count Asian/Pacific Islander female discharges for children ages 0–10:

```bash
python scripts/count_young_asian_female_discharges.py
```

The script reports the unweighted discharge count and the `DISCWT` survey-weighted national estimate. NIS identifies discharge records rather than unique longitudinal patients. A lightweight dataset fingerprint makes unchanged reruns use the cached aggregate in `outputs/cache/`; pass `--refresh` to force recomputation.

### Run without using the terminal

On macOS, double-click `launch_notebook.command` in Finder. It opens `notebooks/demographic_smoke_test.ipynb` in JupyterLab; click the notebook's **Run** button to execute the analysis. You can also open the notebook directly in an editor with Jupyter support and select the repository's `.venv` Python environment.

## Study phases 1–2

Open `notebooks/phase_1_2_cohort_review.ipynb` and choose **Run → Run All Cells**. It presents:

- the Phase 1 data-readiness audit;
- phenotype decisions still requiring review;
- the cached Phase 2 adult HM cohort;
- annual, subtype, sepsis, palliative-care, and overlap summaries;
- stratum-derived region, division, ownership, location/teaching, and bed-size summaries;
- copy/paste-friendly CSV output below every table.

The versioned phenotype is `config/hm_phenotype_v0_1.json` and remains explicitly marked as a draft until its ICD-10 rules are clinically reviewed. The record-level derived cohort and aggregate outputs remain local under `outputs/` and are not committed.

For a no-terminal workflow on macOS, double-click `launch_phase_1_2.command` in Finder.

## Study phase 3

Double-click `launch_phase_3.command` to open the executed baseline-characteristics notebook. It compares HM discharges with versus without documented sepsis and includes weighted Table 1, missingness review, standardized differences, and copy/paste report text.

## Complication analysis

Double-click `launch_phase_4.command` to open the executed complication-comparison notebook. It reports weighted counts and prevalence for the protocol-specified co-documented complications in HM discharges with versus without documented sepsis.

## CCI, length of stay, and mortality

Double-click `launch_phase_5.command` to open the executed clinical-outcomes notebook. It compares the cancer-excluded Charlson Comorbidity Index, length of stay, and in-hospital mortality between HM discharges with and without documented sepsis.

## Palliative-care utilization

Double-click `launch_phase_6.command` to open the executed Task 12A–B notebook. It reports documented inpatient palliative-care use overall and by sepsis status, including weighted prevalence, confidence intervals, an absolute difference, and a crude weighted odds ratio.

## Palliative-care utilization by HM subtype

Double-click `launch_phase_7.command` to open the executed Command 15 notebook. It reports documented inpatient palliative-care use across mutually exclusive HM subtypes overall and separately for hospitalizations with and without documented sepsis.

## Primary adjusted analysis

Double-click `launch_phase_8.command` to open the executed Commands 16 and 16B notebook. It reports the adjusted association between documented sepsis and documented inpatient palliative-care use, including the adjusted odds ratio, standardized probabilities, and adjusted absolute difference.

## Sepsis-by-HM-subtype interaction

Double-click `launch_phase_9.command` to open the executed Commands 17A and 17B notebook. It reports the joint interaction test and subtype-specific adjusted probabilities and absolute differences.

## In-hospital decedent analysis

Double-click `launch_phase_10.command` to open the executed Commands 19A–19C notebook. It reports documented inpatient palliative-care use among in-hospital decedents overall, by sepsis status, and by mutually exclusive HM subtype.

## Annual palliative-care trends

Double-click `launch_phase_11.command` to open the executed Commands 20A–20C notebook. It reports overall, sepsis-stratified, and sepsis-subtype annual trends in documented inpatient palliative-care use.

## Data handling

- Do not copy source medical data into this repository.
- Do not commit `.env`, Parquet files, DuckDB databases, or generated outputs.
- Keep derived results aggregate and de-identified; review small cell counts before sharing.
- Put reusable analysis code under `src/`, tests under `tests/`, and exploratory notebooks under `notebooks/`.

## Table-output convention

Every report-facing table ends with a clearly labeled total or cohort-denominator row. Only mathematically additive quantities are summed: mutually exclusive counts are totaled and exhaustive category percentages may show 100%. Means, SDs, confidence limits, differences, ratios, SMDs, and p-values are not summed and are displayed as `—`. Tables of overlapping diagnoses use the cohort denominator rather than a misleading sum of condition counts.
