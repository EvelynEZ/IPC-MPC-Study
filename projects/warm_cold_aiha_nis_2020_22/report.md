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
<!-- ALL_VARIABLES_UNWEIGHTED -->
# Unweighted Comparison of Warm Versus Cold AIHA Hospitalizations

This table reports the actual sampled NIS hospitalization records without applying `DISCWT`. Warm AIHA is the reference group. Continuous variables are shown as unweighted mean (SD); categorical and binary variables are shown as unweighted n (%).

P-values are unweighted Welch t-tests for continuous variables, overall Pearson chi-square tests for multilevel categorical variables, and Pearson chi-square tests for binary diagnoses. A p-value shown on the first level of a multilevel variable applies to the variable overall.

## Demographic and clinical characteristics

| Characteristic | Level | Warm AIHA | Cold AIHA | Difference, pp | P-value | Test |
| --- | --- | ---: | ---: | ---: | ---: | --- |
| Age, years | Mean (SD) | 60.15 (19.92) | 70.96 (15.70) | — | <0.001 | Welch t-test |
| Age group | 18–59 | 337 (41.55%) | 101 (17.78%) | -23.77 | <0.001 | Pearson chi-square |
|  | ≥60 | 474 (58.45%) | 467 (82.22%) | 23.77 |  |  |
| Sex | Male | 392 (48.34%) | 210 (36.97%) | -11.36 | <0.001 | Pearson chi-square |
|  | Female | 419 (51.66%) | 358 (63.03%) | 11.36 |  |  |
|  | Missing | 0 (0.00%) | 0 (0.00%) | 0.00 |  |  |
| Race/ethnicity | White | 527 (64.98%) | 401 (70.60%) | 5.62 | <0.001 | Pearson chi-square |
|  | Black | 83 (10.23%) | 54 (9.51%) | -0.73 |  |  |
|  | Hispanic | 130 (16.03%) | 52 (9.15%) | -6.87 |  |  |
|  | Asian/Pacific Islander | 27 (3.33%) | 29 (5.11%) | 1.78 |  |  |
|  | Native American | 3 (0.37%) | 0 (0.00%) | -0.37 |  |  |
|  | Other | 18 (2.22%) | 25 (4.40%) | 2.18 |  |  |
|  | Missing | 23 (2.84%) | 7 (1.23%) | -1.60 |  |  |
| Primary payer | Medicare | 397 (48.95%) | 397 (69.89%) | 20.94 | <0.001 | Pearson chi-square |
|  | Medicaid | 118 (14.55%) | 33 (5.81%) | -8.74 |  |  |
|  | Private insurance | 247 (30.46%) | 120 (21.13%) | -9.33 |  |  |
|  | Self-pay | 34 (4.19%) | 11 (1.94%) | -2.26 |  |  |
|  | No charge | 1 (0.12%) | 1 (0.18%) | 0.05 |  |  |
|  | Other | 13 (1.60%) | 5 (0.88%) | -0.72 |  |  |
|  | Missing | 1 (0.12%) | 1 (0.18%) | 0.05 |  |  |
| ZIP-code income quartile | 0–25th percentile | 194 (23.92%) | 115 (20.25%) | -3.67 | 0.275 | Pearson chi-square |
|  | 26th–50th percentile | 187 (23.06%) | 150 (26.41%) | 3.35 |  |  |
|  | 51st–75th percentile | 222 (27.37%) | 149 (26.23%) | -1.14 |  |  |
|  | 76th–100th percentile | 198 (24.41%) | 145 (25.53%) | 1.11 |  |  |
|  | Missing | 10 (1.23%) | 9 (1.58%) | 0.35 |  |  |

## Hospital characteristics

| Characteristic | Level | Warm AIHA | Cold AIHA | Difference, pp | P-value | Test |
| --- | --- | ---: | ---: | ---: | ---: | --- |
| Hospital region | Northeast | 172 (21.21%) | 123 (21.65%) | 0.45 | 0.024 | Pearson chi-square |
|  | Midwest | 202 (24.91%) | 124 (21.83%) | -3.08 |  |  |
|  | South | 261 (32.18%) | 223 (39.26%) | 7.08 |  |  |
|  | West | 176 (21.70%) | 98 (17.25%) | -4.45 |  |  |
|  | Unknown | 0 (0.00%) | 0 (0.00%) | 0.00 |  |  |
| Hospital bed size | Small | 118 (14.55%) | 101 (17.78%) | 3.23 | 0.030 | Pearson chi-square |
|  | Medium | 227 (27.99%) | 181 (31.87%) | 3.88 |  |  |
|  | Large | 466 (57.46%) | 286 (50.35%) | -7.11 |  |  |
|  | Unknown | 0 (0.00%) | 0 (0.00%) | 0.00 |  |  |
| Hospital location/teaching status | Rural | 43 (5.30%) | 41 (7.22%) | 1.92 | 0.320 | Pearson chi-square |
|  | Urban nonteaching | 103 (12.70%) | 67 (11.80%) | -0.90 |  |  |
|  | Urban teaching | 665 (82.00%) | 460 (80.99%) | -1.01 |  |  |
|  | Unknown | 0 (0.00%) | 0 (0.00%) | 0.00 |  |  |

## Admission year

| Characteristic | Level | Warm AIHA | Cold AIHA | Difference, pp | P-value | Test |
| --- | --- | ---: | ---: | ---: | ---: | --- |
| Admission year | 2020 | 68 (8.38%) | 62 (10.92%) | 2.53 | 0.243 | Pearson chi-square |
|  | 2021 | 349 (43.03%) | 246 (43.31%) | 0.28 |  |  |
|  | 2022 | 394 (48.58%) | 260 (45.77%) | -2.81 |  |  |

## Clinical characteristics and outcomes

| Characteristic | Level | Warm AIHA | Cold AIHA | Difference, pp | P-value | Test |
| --- | --- | ---: | ---: | ---: | ---: | --- |
| Charlson Comorbidity Index | Mean (SD) | 2.09 (2.16) | 2.26 (2.04) | — | 0.122 | Welch t-test |
| Charlson category | 0 | 230 (28.36%) | 138 (24.30%) | -4.06 | 0.049 | Pearson chi-square |
|  | 1–2 | 326 (40.20%) | 217 (38.20%) | -1.99 |  |  |
|  | ≥3 | 255 (31.44%) | 213 (37.50%) | 6.06 |  |  |
| Length of stay, days | Mean (SD) | 6.57 (7.65) | 5.93 (6.65) | — | 0.103 | Welch t-test |
| In-hospital mortality | Survived | 789 (97.29%) | 551 (97.01%) | -0.28 | 0.886 | Pearson chi-square |
|  | Died | 22 (2.71%) | 17 (2.99%) | 0.28 |  |  |
|  | Missing | 0 (0.00%) | 0 (0.00%) | 0.00 |  |  |

## Lymphoid malignancy

| Characteristic | Level | Warm AIHA | Cold AIHA | Difference, pp | P-value | Test |
| --- | --- | ---: | ---: | ---: | ---: | --- |
| Lymphoma spectrum (Hodgkin, NHL, WM/MALT) | Present | 42 (5.18%) | 89 (15.67%) | 10.49 | <0.001 | Pearson chi-square |
| CLL/chronic leukemia group | Present | 84 (10.36%) | 35 (6.16%) | -4.20 | 0.006 | Pearson chi-square |
| Plasma-cell neoplasms | Present | 9 (1.11%) | 5 (0.88%) | -0.23 | 0.676 | Pearson chi-square |
| ANY_LYMPHOID_MALIGNANCY | Present | 132 (16.28%) | 124 (21.83%) | 5.55 | 0.009 | Pearson chi-square |

## Acute thrombosis

| Characteristic | Level | Warm AIHA | Cold AIHA | Difference, pp | P-value | Test |
| --- | --- | ---: | ---: | ---: | ---: | --- |
| Acute pulmonary embolism | Present | 29 (3.58%) | 9 (1.58%) | -1.99 | 0.026 | Pearson chi-square |
| Acute deep-vein thrombosis | Present | 25 (3.08%) | 17 (2.99%) | -0.09 | 0.924 | Pearson chi-square |
| Splanchnic-vein thrombosis | Present | 5 (0.62%) | 3 (0.53%) | -0.09 | 0.832 | Pearson chi-square |
| Other acute venous thrombosis | Present | 6 (0.74%) | 1 (0.18%) | -0.56 | 0.147 | Pearson chi-square |
| Acute ischemic stroke | Present | 8 (0.99%) | 7 (1.23%) | 0.25 | 0.665 | Pearson chi-square |
| Acute myocardial infarction | Present | 9 (1.11%) | 11 (1.94%) | 0.83 | 0.206 | Pearson chi-square |
| Other acute arterial embolism or thrombosis | Present | 3 (0.37%) | 1 (0.18%) | -0.19 | 0.510 | Pearson chi-square |
| ANY_VTE | Present | 54 (6.66%) | 24 (4.23%) | -2.43 | 0.054 | Pearson chi-square |
| ANY_ARTERIAL_THROMBOSIS | Present | 20 (2.47%) | 19 (3.35%) | 0.88 | 0.333 | Pearson chi-square |
| ANY_ACUTE_THROMBOSIS | Present | 71 (8.75%) | 42 (7.39%) | -1.36 | 0.365 | Pearson chi-square |

## Selected complications

| Characteristic | Level | Warm AIHA | Cold AIHA | Difference, pp | P-value | Test |
| --- | --- | ---: | ---: | ---: | ---: | --- |
| Acute kidney injury | Present | 169 (20.84%) | 92 (16.20%) | -4.64 | 0.030 | Pearson chi-square |
| Acute respiratory failure | Present | 47 (5.80%) | 27 (4.75%) | -1.04 | 0.398 | Pearson chi-square |
| Sepsis | Present | 29 (3.58%) | 30 (5.28%) | 1.71 | 0.123 | Pearson chi-square |

## Total

| Characteristic | Level | Warm AIHA | Cold AIHA | Difference, pp | P-value | Test |
| --- | --- | ---: | ---: | ---: | ---: | --- |
| Total cohort denominator | All hospitalizations | 811 (100.00%) | 568 (100.00%) | — | — | — |

Lymphoid-malignancy component rows and thrombosis component rows may overlap. Composite rows count each hospitalization once. Diagnoses searched all 40 diagnosis positions.
<!-- ADJUSTED_MORTALITY_MODEL -->
# Adjusted Association Between AIHA Subtype and In-Hospital Mortality

A multivariable logistic-regression model used in-hospital mortality (`DIED=1`) as the outcome and AIHA subtype as the primary exposure. Warm AIHA was the reference group. The model adjusted for continuous age, sex, race/ethnicity, associated lymphoid malignancy, hospital region, hospital teaching status, hospital bed size, and continuous Charlson Comorbidity Index. Teaching status was coded as urban teaching versus nonteaching; the nonteaching category combines rural and urban nonteaching hospitals.

## Primary exposure result

| Exposure comparison | Adjusted odds ratio | 95% CI | P-value |
| --- | ---: | ---: | ---: |
| Cold AIHA versus warm AIHA | 1.199 | 0.598–2.401 | 0.609 |

## Model diagnostics

| Measure | Value |
| --- | ---: |
| Analysis hospitalizations | 1,379 |
| Deaths | 39 |
| Warm AIHA deaths / hospitalizations | 22 / 811 |
| Cold AIHA deaths / hospitalizations | 17 / 568 |
| Estimated model parameters | 18 |
| Model converged | True |

`DISCWT` is 5 for every record in this 2020–2022 analytic cohort. Therefore weighted and unweighted coefficient estimates are identical; the model uses each sampled hospitalization once so that model-based standard errors are not artificially reduced by treating the weight as replicated observations.

Interpret cautiously: only 39 deaths were observed relative to the number of adjustment parameters, so the estimate may be imprecise and the model is vulnerable to sparse-data overfitting.

## Full adjusted model

| Covariate contrast | Adjusted odds ratio | 95% CI | P-value |
| --- | ---: | ---: | ---: |
| Cold AIHA vs warm AIHA | 1.199 | 0.598–2.401 | 0.609 |
| Sex: Female | 0.420 | 0.207–0.853 | 0.016 |
| Race/ethnicity: Asian/Pacific Islander | 0.420 | 0.053–3.316 | 0.410 |
| Race/ethnicity: Black | 0.263 | 0.034–2.026 | 0.200 |
| Race/ethnicity: Hispanic | 0.827 | 0.268–2.553 | 0.741 |
| Race/ethnicity: Missing | 5.515 | 1.373–22.159 | 0.016 |
| Race/ethnicity: Native American | 22.792 | 0.982–529.200 | 0.051 |
| Race/ethnicity: Other | 1.454 | 0.316–6.691 | 0.631 |
| Hospital region: Midwest | 1.184 | 0.346–4.044 | 0.788 |
| Hospital region: South | 1.739 | 0.600–5.033 | 0.308 |
| Hospital region: West | 3.264 | 1.111–9.588 | 0.031 |
| Hospital teaching status: Teaching | 2.994 | 0.895–10.015 | 0.075 |
| Hospital bed size: Large | 1.165 | 0.418–3.248 | 0.770 |
| Hospital bed size: Medium | 1.478 | 0.506–4.322 | 0.475 |
| Age, per year | 1.017 | 0.993–1.042 | 0.171 |
| Associated lymphoid malignancy: yes vs no | 1.030 | 0.471–2.253 | 0.941 |
| Charlson Comorbidity Index, per point | 1.159 | 1.001–1.343 | 0.048 |
<!-- ADJUSTED_LOS_MODEL -->
# Adjusted Association Between AIHA Subtype and Length of Stay

A multivariable linear-regression model used hospital length of stay in days as the continuous outcome and AIHA subtype as the primary exposure. Warm AIHA was the reference group. The model adjusted for continuous age, sex, race/ethnicity, associated lymphoid malignancy, hospital region, hospital teaching status, hospital bed size, and continuous Charlson Comorbidity Index. Teaching status was coded as urban teaching versus nonteaching; nonteaching combines rural and urban nonteaching hospitals.

Heteroskedasticity-consistent HC3 standard errors were used because length of stay is right-skewed. The exposure coefficient is an adjusted mean difference in days, calculated as cold AIHA minus warm AIHA.

## Primary exposure result

| Exposure comparison | Adjusted mean difference, days | 95% CI | P-value |
| --- | ---: | ---: | ---: |
| Cold AIHA versus warm AIHA | -0.600 | -1.436 to 0.236 | 0.160 |

## Model diagnostics

| Measure | Value |
| --- | ---: |
| Analysis hospitalizations | 1,379 |
| Unadjusted mean LOS, warm AIHA | 6.57 days |
| Unadjusted mean LOS, cold AIHA | 5.93 days |
| R-squared | 0.0343 |
| Standard-error estimator | HC3 |

`DISCWT` equals 5 for every record in this 2020–2022 cohort, so applying it does not change coefficients or fitted values. Each sampled hospitalization was used once for variance estimation.

## Full adjusted model

| Covariate contrast | Adjusted coefficient, days | 95% CI | P-value |
| --- | ---: | ---: | ---: |
| Cold AIHA vs warm AIHA | -0.600 | -1.436 to 0.236 | 0.160 |
| Sex: Female | -0.399 | -1.159 to 0.361 | 0.303 |
| Race/ethnicity: Asian/Pacific Islander | -0.360 | -2.195 to 1.475 | 0.701 |
| Race/ethnicity: Black | 0.898 | -1.119 to 2.915 | 0.383 |
| Race/ethnicity: Hispanic | -0.200 | -1.413 to 1.012 | 0.746 |
| Race/ethnicity: Missing | -0.981 | -2.018 to 0.055 | 0.064 |
| Race/ethnicity: Native American | 2.685 | -9.139 to 14.509 | 0.656 |
| Race/ethnicity: Other | 1.114 | -1.712 to 3.939 | 0.440 |
| Hospital region: Midwest | 0.059 | -0.782 to 0.899 | 0.891 |
| Hospital region: South | 1.004 | 0.119 to 1.888 | 0.026 |
| Hospital region: West | 1.113 | -0.252 to 2.477 | 0.110 |
| Hospital teaching status: Teaching | 1.416 | 0.552 to 2.281 | 0.001 |
| Hospital bed size: Large | 1.317 | 0.395 to 2.239 | 0.005 |
| Hospital bed size: Medium | 0.199 | -0.694 to 1.092 | 0.662 |
| Age, per year | 0.001 | -0.022 to 0.023 | 0.935 |
| Associated lymphoid malignancy: yes vs no | -0.120 | -1.197 to 0.957 | 0.827 |
| Charlson Comorbidity Index, per point | 0.337 | 0.080 to 0.595 | 0.010 |
