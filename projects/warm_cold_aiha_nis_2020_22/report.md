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
<!-- BASELINE_CHARACTERISTICS -->
# Warm Versus Cold AIHA: Baseline Characteristics

Warm AIHA is the reference group. Percentages, means, standard deviations, and displayed weighted counts use `DISCWT`. P-values are Welch tests for continuous variables and overall Pearson chi-square tests for categorical variables; missing levels are excluded from categorical tests.

The Charlson Comorbidity Index uses all diagnosis positions, standard Quan ICD-10 components and weights, cancer/metastatic-cancer and diabetes/liver hierarchy, and no age points.

| Characteristic | Level | Warm unweighted n | Warm weighted summary | Cold unweighted n | Cold weighted summary | Difference, pp | SMD | P-value | Test |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| Age, years | Mean (SD) | 811 | 60.15 (19.90) | 568 | 70.96 (15.69) | — | 0.603 | <0.001 | Welch t-test |
| Charlson Comorbidity Index | Mean (SD) | 811 | 2.09 (2.16) | 568 | 2.26 (2.04) | — | 0.084 | 0.122 | Welch t-test |
| Length of stay, days | Mean (SD) | 811 | 6.57 (7.64) | 568 | 5.93 (6.64) | — | -0.088 | 0.103 | Welch t-test |
| Age group | 18–59 | 337 | 1,685 (41.55%) | 101 | 505 (17.78%) | -23.77 | -0.539 | <0.001 | Pearson chi-square |
|  | ≥60 | 474 | 2,370 (58.45%) | 467 | 2,335 (82.22%) | 23.77 | 0.539 |  |  |
| Sex | Male | 392 | 1,960 (48.34%) | 210 | 1,050 (36.97%) | -11.36 | -0.231 | <0.001 | Pearson chi-square |
|  | Female | 419 | 2,095 (51.66%) | 358 | 1,790 (63.03%) | 11.36 | 0.231 |  |  |
|  | Missing | 0 | 0 (0.00%) | 0 | 0 (0.00%) | 0.0 | 0.0 |  |  |
| Race/ethnicity | White | 527 | 2,635 (64.98%) | 401 | 2,005 (70.60%) | 5.62 | 0.12 | <0.001 | Pearson chi-square |
|  | Black | 83 | 415 (10.23%) | 54 | 270 (9.51%) | -0.73 | -0.024 |  |  |
|  | Hispanic | 130 | 650 (16.03%) | 52 | 260 (9.15%) | -6.87 | -0.208 |  |  |
|  | Asian/Pacific Islander | 27 | 135 (3.33%) | 29 | 145 (5.11%) | 1.78 | 0.088 |  |  |
|  | Native American | 3 | 15 (0.37%) | 0 | 0 (0.00%) | -0.37 | -0.086 |  |  |
|  | Other | 18 | 90 (2.22%) | 25 | 125 (4.40%) | 2.18 | 0.122 |  |  |
|  | Missing | 23 | 115 (2.84%) | 7 | 35 (1.23%) | -1.6 | -0.114 |  |  |
| Primary payer | Medicare | 397 | 1,985 (48.95%) | 397 | 1,985 (69.89%) | 20.94 | 0.437 | <0.001 | Pearson chi-square |
|  | Medicaid | 118 | 590 (14.55%) | 33 | 165 (5.81%) | -8.74 | -0.292 |  |  |
|  | Private insurance | 247 | 1,235 (30.46%) | 120 | 600 (21.13%) | -9.33 | -0.214 |  |  |
|  | Self-pay | 34 | 170 (4.19%) | 11 | 55 (1.94%) | -2.26 | -0.131 |  |  |
|  | No charge | 1 | 5 (0.12%) | 1 | 5 (0.18%) | 0.05 | 0.014 |  |  |
|  | Other | 13 | 65 (1.60%) | 5 | 25 (0.88%) | -0.72 | -0.065 |  |  |
|  | Missing | 1 | 5 (0.12%) | 1 | 5 (0.18%) | 0.05 | 0.014 |  |  |
| ZIP-code income quartile | 0–25th percentile | 194 | 970 (23.92%) | 115 | 575 (20.25%) | -3.67 | -0.089 | 0.275 | Pearson chi-square |
|  | 26th–50th percentile | 187 | 935 (23.06%) | 150 | 750 (26.41%) | 3.35 | 0.078 |  |  |
|  | 51st–75th percentile | 222 | 1,110 (27.37%) | 149 | 745 (26.23%) | -1.14 | -0.026 |  |  |
|  | 76th–100th percentile | 198 | 990 (24.41%) | 145 | 725 (25.53%) | 1.11 | 0.026 |  |  |
|  | Missing | 10 | 50 (1.23%) | 9 | 45 (1.58%) | 0.35 | 0.03 |  |  |
| Hospital region | Northeast | 172 | 860 (21.21%) | 123 | 615 (21.65%) | 0.45 | 0.011 | 0.024 | Pearson chi-square |
|  | Midwest | 202 | 1,010 (24.91%) | 124 | 620 (21.83%) | -3.08 | -0.073 |  |  |
|  | South | 261 | 1,305 (32.18%) | 223 | 1,115 (39.26%) | 7.08 | 0.148 |  |  |
|  | West | 176 | 880 (21.70%) | 98 | 490 (17.25%) | -4.45 | -0.112 |  |  |
|  | Unknown | 0 | 0 (0.00%) | 0 | 0 (0.00%) | 0.0 | 0.0 |  |  |
| Hospital bed size | Small | 118 | 590 (14.55%) | 101 | 505 (17.78%) | 3.23 | 0.088 | 0.030 | Pearson chi-square |
|  | Medium | 227 | 1,135 (27.99%) | 181 | 905 (31.87%) | 3.88 | 0.085 |  |  |
|  | Large | 466 | 2,330 (57.46%) | 286 | 1,430 (50.35%) | -7.11 | -0.143 |  |  |
|  | Unknown | 0 | 0 (0.00%) | 0 | 0 (0.00%) | 0.0 | 0.0 |  |  |
| Hospital location/teaching status | Rural | 43 | 215 (5.30%) | 41 | 205 (7.22%) | 1.92 | 0.079 | 0.320 | Pearson chi-square |
|  | Urban nonteaching | 103 | 515 (12.70%) | 67 | 335 (11.80%) | -0.9 | -0.028 |  |  |
|  | Urban teaching | 665 | 3,325 (82.00%) | 460 | 2,300 (80.99%) | -1.01 | -0.026 |  |  |
|  | Unknown | 0 | 0 (0.00%) | 0 | 0 (0.00%) | 0.0 | 0.0 |  |  |
| Charlson category | 0 | 230 | 1,150 (28.36%) | 138 | 690 (24.30%) | -4.06 | -0.092 | 0.049 | Pearson chi-square |
|  | 1–2 | 326 | 1,630 (40.20%) | 217 | 1,085 (38.20%) | -1.99 | -0.041 |  |  |
|  | ≥3 | 255 | 1,275 (31.44%) | 213 | 1,065 (37.50%) | 6.06 | 0.128 |  |  |
| In-hospital mortality | Survived | 789 | 3,945 (97.29%) | 551 | 2,755 (97.01%) | -0.28 | -0.017 | 0.886 | Pearson chi-square |
|  | Died | 22 | 110 (2.71%) | 17 | 85 (2.99%) | 0.28 | 0.017 |  |  |
|  | Missing | 0 | 0 (0.00%) | 0 | 0 (0.00%) | 0.0 | 0.0 |  |  |
| Total | All final-cohort hospitalizations | 811 | 4,055 | 568 | 2,840 | — | — | — | — |
<!-- LYMPHOID_MALIGNANCY -->
# Any Lymphoid Malignancy by AIHA Type

Warm AIHA is the reference group. The phenotype searches all 40 diagnosis positions and reuses the prior HM project’s lymphoma, CLL/chronic-leukemia, and plasma-cell inclusion and exclusion rules.

Component rows can overlap. `ANY_LYMPHOID_MALIGNANCY` is their union and counts each hospitalization once.

| Malignancy definition | Warm unweighted n | Warm weighted n (%) | Cold unweighted n | Cold weighted n (%) | Difference, pp | SMD | P-value |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Lymphoma spectrum (Hodgkin, NHL, WM/MALT) | 42 | 210 (5.18%) | 89 | 445 (15.67%) | 10.49 | 0.348 | <0.001 |
| CLL/chronic leukemia group | 84 | 420 (10.36%) | 35 | 175 (6.16%) | -4.2 | -0.153 | 0.006 |
| Plasma-cell neoplasms | 9 | 45 (1.11%) | 5 | 25 (0.88%) | -0.23 | -0.023 | 0.676 |
| ANY_LYMPHOID_MALIGNANCY | 132 | 660 (16.28%) | 124 | 620 (21.83%) | 5.55 | 0.142 | 0.009 |
| Total cohort denominator | 811 | 4,055 (100.00%) | 568 | 2,840 (100.00%) | — | — | — |
<!-- ACUTE_THROMBOSIS -->
# Acute Thrombosis by AIHA Type

Warm AIHA is the reference group. Every phenotype searches all 40 diagnosis positions. Code-specific families can overlap; each composite counts a hospitalization once.

`ANY_VTE` is the union of pulmonary embolism, acute DVT, splanchnic-vein thrombosis, and other acute venous thrombosis. `ANY_ARTERIAL_THROMBOSIS` is the union of acute ischemic stroke, acute MI, and other acute arterial embolism/thrombosis. `ANY_ACUTE_THROMBOSIS` is the union of those two composites.

| Diagnosis family | Warm unweighted n | Warm weighted n (%) | Cold unweighted n | Cold weighted n (%) | Difference, pp | SMD | P-value | FDR-adjusted p |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Acute pulmonary embolism | 29 | 145 (3.58%) | 9 | 45 (1.58%) | -1.99 | -0.126 | 0.026 | 0.262 |
| Acute deep-vein thrombosis | 25 | 125 (3.08%) | 17 | 85 (2.99%) | -0.09 | -0.005 | 0.924 | 0.924 |
| Splanchnic-vein thrombosis | 5 | 25 (0.62%) | 3 | 15 (0.53%) | -0.09 | -0.012 | 0.832 | 0.924 |
| Other acute venous thrombosis | 6 | 30 (0.74%) | 1 | 5 (0.18%) | -0.56 | -0.084 | 0.147 | 0.490 |
| Acute ischemic stroke | 8 | 40 (0.99%) | 7 | 35 (1.23%) | 0.25 | 0.023 | 0.665 | 0.831 |
| Acute myocardial infarction | 9 | 45 (1.11%) | 11 | 55 (1.94%) | 0.83 | 0.068 | 0.206 | 0.515 |
| Other acute arterial embolism or thrombosis | 3 | 15 (0.37%) | 1 | 5 (0.18%) | -0.19 | -0.037 | 0.510 | 0.729 |
| ANY_VTE | 54 | 270 (6.66%) | 24 | 120 (4.23%) | -2.43 | -0.107 | 0.054 | 0.271 |
| ANY_ARTERIAL_THROMBOSIS | 20 | 100 (2.47%) | 19 | 95 (3.35%) | 0.88 | 0.052 | 0.333 | 0.608 |
| ANY_ACUTE_THROMBOSIS | 71 | 355 (8.75%) | 42 | 210 (7.39%) | -1.36 | -0.05 | 0.365 | 0.608 |
| Total cohort denominator | 811 | 4,055 (100.00%) | 568 | 2,840 (100.00%) | — | — | — | — |
<!-- SELECTED_COMPLICATIONS -->
# Selected Complications by AIHA Type

Warm AIHA is the reference group. Diagnoses were identified in any of the 40 diagnosis positions after removing decimal points and spaces and converting codes to uppercase. Acute respiratory failure uses ICD-10-CM `J96.0*` (`J960*` normalized); `A960` was treated as a typographical error because it does not denote acute respiratory failure.

| Complication | ICD-10-CM | Warm unweighted n | Warm weighted n (%) | Cold unweighted n | Cold weighted n (%) | Difference, pp | SMD | P-value |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Acute kidney injury | N17* | 169 | 845 (20.84%) | 92 | 460 (16.20%) | -4.64 | -0.12 | 0.030 |
| Acute respiratory failure | J960* | 47 | 235 (5.80%) | 27 | 135 (4.75%) | -1.04 | -0.047 | 0.398 |
| Sepsis | A41* | 29 | 145 (3.58%) | 30 | 150 (5.28%) | 1.71 | 0.083 | 0.123 |
| Total cohort denominator | — | 811 | 4,055 (100.00%) | 568 | 2,840 (100.00%) | — | — | — |

P-values are two-sided Pearson chi-square tests comparing warm and cold AIHA. Diagnoses are co-documented during the hospitalization and do not establish temporal ordering.
