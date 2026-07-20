from __future__ import annotations

import unittest

from src.phase_5_cci_los_mortality import CCI_COMPONENTS, score_from_flags


class CharlsonScoreTests(unittest.TestCase):
    def test_cancer_components_are_excluded(self) -> None:
        self.assertNotIn("canc", CCI_COMPONENTS)
        self.assertNotIn("metacanc", CCI_COMPONENTS)

    def test_empty_score_is_zero(self) -> None:
        self.assertEqual(score_from_flags({}), 0)

    def test_original_weights_are_added(self) -> None:
        self.assertEqual(score_from_flags({"chf": True, "rend": True}), 3)

    def test_diabetes_is_not_double_counted(self) -> None:
        self.assertEqual(score_from_flags({"diab": True, "diabwc": True}), 2)

    def test_liver_disease_is_not_double_counted(self) -> None:
        self.assertEqual(score_from_flags({"mld": True, "msld": True}), 3)


if __name__ == "__main__":
    unittest.main()
