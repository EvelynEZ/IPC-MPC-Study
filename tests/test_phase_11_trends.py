from __future__ import annotations

import unittest
import numpy as np
import pandas as pd

from src.phase_11_trends import eapc, sepsis_design, subtype_design


class TrendModelTests(unittest.TestCase):
    def test_eapc_zero_coefficient(self) -> None:
        result = eapc(0.0, 0.01)
        self.assertEqual(result["eapc_percent"], 0.0)

    def test_sepsis_interaction_column(self) -> None:
        frame = pd.DataFrame({"year_category": [2016, 2017], "sepsis": [0, 1]})
        x, names = sepsis_design(frame)
        self.assertEqual(names[-1], "Year × sepsis")
        self.assertEqual(x[:, -1].tolist(), [0.0, 1.0])

    def test_subtype_design_has_eight_interactions(self) -> None:
        frame = pd.DataFrame({"year_category": [2016], "hm_subtype_label": ["Lymphoma"]})
        _, names = subtype_design(frame)
        self.assertEqual(sum(name.startswith("Year ×") for name in names), 8)


if __name__ == "__main__":
    unittest.main()
