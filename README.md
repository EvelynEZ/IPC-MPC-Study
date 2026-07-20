# IPC-MPC Study

Reproducible analysis workspace for the NIS 2016–2022 dataset. Source data stays outside this Git repository and is queried directly as a partitioned Parquet dataset.

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

## Data handling

- Do not copy source medical data into this repository.
- Do not commit `.env`, Parquet files, DuckDB databases, or generated outputs.
- Keep derived results aggregate and de-identified; review small cell counts before sharing.
- Put reusable analysis code under `src/`, tests under `tests/`, and exploratory notebooks under `notebooks/`.
