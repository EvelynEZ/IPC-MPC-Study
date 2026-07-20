from __future__ import annotations

import json
import unittest

import duckdb

from src.phase_1_2 import DEFAULT_CONFIG, code_match_sql, subtype_case_sql


class PhenotypeSqlTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.config = json.loads(DEFAULT_CONFIG.read_text())
        cls.connection = duckdb.connect(database=":memory:")

    @classmethod
    def tearDownClass(cls) -> None:
        cls.connection.close()

    def classify(self, code: str) -> str | None:
        expression = subtype_case_sql(self.config, "code")
        return self.connection.execute(
            f"SELECT {expression} FROM (SELECT ? AS code)", [code]
        ).fetchone()[0]

    def test_expected_subtypes(self) -> None:
        examples = {
            "C8339": "lymphoma",
            "C9200": "aml",
            "C9210": "cml",
            "C9110": "cll_chronic_leukemia",
            "C9100": "all",
            "C9140": "other_leukemia",
            "C9000": "myeloma_plasma_cell",
            "D469": "mds",
            "D473": "mpn",
        }
        for code, expected in examples.items():
            with self.subTest(code=code):
                self.assertEqual(self.classify(code), expected)

    def test_remission_exclusions_from_proposal(self) -> None:
        for code in ["C9201", "C9211", "C9501", "C9001"]:
            with self.subTest(code=code):
                self.assertIsNone(self.classify(code))

    def test_sepsis_is_currently_a41_only(self) -> None:
        expression = code_match_sql("code", self.config["sepsis"])
        a41, a40 = self.connection.execute(
            f"SELECT {expression} FROM (SELECT 'A419' AS code) "
            f"UNION ALL SELECT {expression} FROM (SELECT 'A409' AS code)"
        ).fetchall()
        self.assertTrue(a41[0])
        self.assertFalse(a40[0])


if __name__ == "__main__":
    unittest.main()
