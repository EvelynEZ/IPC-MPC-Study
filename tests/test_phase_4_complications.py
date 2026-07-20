from __future__ import annotations

import unittest

import duckdb

from src.phase_4_complications import (
    COMPLICATIONS,
    benjamini_hochberg,
    code_condition,
    two_proportion_p_value,
)


class ComplicationCodeTests(unittest.TestCase):
    def setUp(self) -> None:
        self.connection = duckdb.connect(database=":memory:")

    def tearDown(self) -> None:
        self.connection.close()

    def matches(self, rule_id: str, code: str) -> bool:
        rule = next(item for item in COMPLICATIONS if item["id"] == rule_id)
        expression = code_condition(rule)
        return self.connection.execute(
            f"SELECT {expression} FROM (SELECT ? AS code)", [code]
        ).fetchone()[0]

    def test_septic_shock_is_also_severe_sepsis(self) -> None:
        self.assertTrue(self.matches("septic_shock", "R6521"))
        self.assertTrue(self.matches("severe_sepsis", "R6521"))

    def test_chronic_af_is_excluded(self) -> None:
        self.assertTrue(self.matches("afib_flutter", "I480"))
        self.assertFalse(self.matches("afib_flutter", "I482"))
        self.assertFalse(self.matches("afib_flutter", "I4821"))

    def test_dvt_list_does_not_include_chronic_code(self) -> None:
        self.assertTrue(self.matches("acute_lower_extremity_dvt", "I82401"))
        self.assertFalse(self.matches("acute_lower_extremity_dvt", "I82501"))

    def test_two_proportion_test(self) -> None:
        self.assertEqual(two_proportion_p_value(10, 100, 10, 100), 1.0)
        self.assertLess(two_proportion_p_value(50, 100, 10, 100), 0.001)

    def test_benjamini_hochberg_is_monotonic_and_bounded(self) -> None:
        adjusted = benjamini_hochberg([0.001, 0.01, 0.5])
        self.assertTrue(all(0 <= value <= 1 for value in adjusted))
        self.assertLessEqual(adjusted[0], adjusted[1])
        self.assertLessEqual(adjusted[1], adjusted[2])


if __name__ == "__main__":
    unittest.main()
