from __future__ import annotations

import unittest

import duckdb

from src.phase_4_complications import COMPLICATIONS, code_condition


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


if __name__ == "__main__":
    unittest.main()
