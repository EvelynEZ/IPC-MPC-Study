from __future__ import annotations

import unittest

import duckdb


class SepticShockDefinitionTests(unittest.TestCase):
    def test_exact_normalized_r6521_definition(self) -> None:
        connection = duckdb.connect(":memory:")
        for codes, expected in [(["R6521"], True), (["R652"], False), (["R6520"], False), (["A419"], False)]:
            actual = connection.execute("SELECT list_contains(?, 'R6521')", [codes]).fetchone()[0]
            self.assertEqual(actual, expected)
        connection.close()

    def test_diagnosis_based_bmt_definition(self) -> None:
        connection = duckdb.connect(":memory:")
        for codes, expected in [
            (["Z9481"], True),
            (["Z9484"], True),
            (["Z9481", "Z9484"], True),
            (["Z948"], False),
            (["Z9482"], False),
        ]:
            actual = connection.execute(
                "SELECT list_contains(?, 'Z9481') OR list_contains(?, 'Z9484')",
                [codes, codes],
            ).fetchone()[0]
            self.assertEqual(actual, expected)
        connection.close()


if __name__ == "__main__":
    unittest.main()
