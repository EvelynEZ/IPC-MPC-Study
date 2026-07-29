# Warm and Cold AIHA — NIS 2020–2022

## Preliminary cohort definition

Population: all adult NIS hospitalizations (`AGE >= 18`) from 2020 through 2022. Warm AIHA is exact normalized `D59.11`; cold AIHA is exact normalized `D59.12`.

**Regression reference category:** Warm AIHA. All future regression coefficients and odds ratios for AIHA type will compare cold AIHA with warm AIHA unless explicitly labeled otherwise.

## All Adult Hospitalization Denominators

| Admission year | Unweighted adult hospitalizations | DISCWT-weighted adult hospitalizations |
| --- | ---: | ---: |
| 2020 | 5,533,477 | 27,667,386 |
| 2021 | 5,688,355 | 28,441,739 |
| 2022 | 5,571,320 | 27,856,590 |
| Total | 16,793,152 | 83,965,715 |

## Primary-diagnosis definition

| Phenotype | Unweighted hospitalizations | DISCWT-weighted hospitalizations |
| --- | ---: | ---: |
| Warm AIHA (D59.11) | 571 | 2,855 |
| Cold AIHA (D59.12) | 232 | 1,160 |
| Total unique hospitalizations | 803 | 4,015 |

## Expanded definition: first three diagnosis positions

| Phenotype | Unweighted hospitalizations | DISCWT-weighted hospitalizations |
| --- | ---: | ---: |
| Warm AIHA (D59.11) | 827 | 4,135 |
| Cold AIHA (D59.12) | 576 | 2,880 |
| Both warm and cold codes | 3 | 15 |
| Total unique hospitalizations | 1,400 | 7,000 |

Code-specific counts overlap when both codes are present; the total row counts unique hospitalizations.

## Annual distribution using the first three positions

| Year | Warm unweighted | Warm weighted | Cold unweighted | Cold weighted | Unique unweighted | Unique weighted |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 2020 | 68 | 340 | 62 | 310 | 130 | 650 |
| 2021 | 356 | 1,780 | 249 | 1,245 | 604 | 3,020 |
| 2022 | 403 | 2,015 | 265 | 1,325 | 666 | 3,330 |
| Total | 827 | 4,135 | 576 | 2,880 | 1,400 | 7,000 |

## All Specified AIHA Codes in the First Three Diagnosis Positions, 2020–2022 Combined

| ICD-10-CM definition | Unweighted hospitalizations | DISCWT-weighted hospitalizations |
| --- | ---: | ---: |
| D59.10 — Autoimmune hemolytic anemia, unspecified | 1,758 | 8,790 |
| D59.11 — Warm autoimmune hemolytic anemia | 827 | 4,135 |
| D59.12 — Cold autoimmune hemolytic anemia | 576 | 2,880 |
| D59.13 — Mixed-type autoimmune hemolytic anemia | 46 | 230 |
| D59.19 — Other autoimmune hemolytic anemia | 237 | 1,185 |
| Hospitalizations containing more than one listed code | 13 | 65 |
| Total unique hospitalizations containing any listed code | 3,431 | 17,155 |

Code-specific rows are not mutually exclusive. Hospitalizations containing multiple listed codes are counted once in the total unique cohort.

## Final Warm/Cold AIHA Cohort Flow

Inclusion requires age 18 years or older and D59.11 or D59.12 in DX1–DX3. Hospitalizations containing both codes anywhere in DX1–DX40 are excluded.

| Cohort step | Unweighted hospitalizations | DISCWT-weighted hospitalizations |
| --- | ---: | ---: |
| Identified through DX1 | 803 | — |
| Additional identified through DX2 | 314 | — |
| Additional identified through DX3 | 283 | — |
| Candidates before overlap exclusion | 1,400 | 7,000 |
| Excluded: both D59.11 and D59.12 anywhere | 21 | 105 |
| Final analytic cohort | 1,379 | 6,895 |

These preliminary counts are hospitalization-based and do not identify unique patients. The definitive position rule and exclusion criteria remain to be specified in the study protocol.
