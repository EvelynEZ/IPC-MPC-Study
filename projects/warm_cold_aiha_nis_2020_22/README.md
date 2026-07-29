# Warm and Cold AIHA — NIS 2020–2022

This project is isolated from the IPC-MPC septic-shock study while using the same external NIS Parquet source.

## Current preliminary cohort definitions

- Population: all adult (`AGE >= 18`) NIS hospitalizations, 2020–2022.
- Warm AIHA: exact normalized ICD-10-CM `D59.11` (`D5911`).
- Cold AIHA: exact normalized ICD-10-CM `D59.12` (`D5912`).
- Two definitions are reported: code in the primary diagnosis position and code in any of the first three diagnosis positions.
- Each hospitalization is counted once in combined estimates even if both codes appear.
- National estimates use `DISCWT`.

Double-click `launch_notebook.command`, then choose **Run → Run All Cells**. Generated CSV and JSON files are written only to this project's `outputs/` directory. The report-facing Markdown file is `report.md`.

The analytic plan and final cohort definition remain open for clinical review.
