from __future__ import annotations

import unittest

import numpy as np
import pandas as pd

from src.phase_8_primary_adjusted import design_matrix


class PrimaryAdjustedModelTests(unittest.TestCase):
    def test_forced_sepsis_changes_only_sepsis_column(self) -> None:
        row = {"age": 65, "sepsis": 0}
        for variable, levels in __import__("src.phase_8_primary_adjusted", fromlist=["CATEGORY_SPECS"]).CATEGORY_SPECS.items():
            row[variable] = levels[0]
        frame = pd.DataFrame([row])
        x0, names = design_matrix(frame, 0)
        x1, _ = design_matrix(frame, 1)
        changed = np.flatnonzero((x1 - x0)[0])
        self.assertEqual(changed.tolist(), [names.index("Documented sepsis")])

    def test_reference_profile_has_zero_dummy_columns(self) -> None:
        from src.phase_8_primary_adjusted import CATEGORY_SPECS
        row = {"age": 65, "sepsis": 0, **{variable: levels[0] for variable, levels in CATEGORY_SPECS.items()}}
        x, _ = design_matrix(pd.DataFrame([row]))
        self.assertEqual(x[0, 0], 1)
        self.assertTrue(np.all(x[0, 1:] == 0))


if __name__ == "__main__":
    unittest.main()
