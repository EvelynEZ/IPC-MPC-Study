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

## Data handling

- Do not copy source medical data into this repository.
- Do not commit `.env`, Parquet files, DuckDB databases, or generated outputs.
- Keep derived results aggregate and de-identified; review small cell counts before sharing.
- Put reusable analysis code under `src/`, tests under `tests/`, and exploratory notebooks under `notebooks/`.
