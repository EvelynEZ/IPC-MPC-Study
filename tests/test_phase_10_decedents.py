from __future__ import annotations

import unittest
import pandas as pd

from src.phase_8_primary_adjusted import CATEGORY_SPECS
from src.phase_10_decedents import subtype_only_design


class DecedentModelTests(unittest.TestCase):
    def test_subtype_only_design_excludes_sepsis(self) -> None:
        row = {"age": 65, "sepsis": 1, **{variable: levels[0] for variable, levels in CATEGORY_SPECS.items()}}
        _, names = subtype_only_design(pd.DataFrame([row]))
        self.assertNotIn("Documented sepsis", names)

    def test_forced_subtype_changes_design(self) -> None:
        row = {"age": 65, "sepsis": 1, **{variable: levels[0] for variable, levels in CATEGORY_SPECS.items()}}
        x0, _ = subtype_only_design(pd.DataFrame([row]), "Lymphoma")
        x1, _ = subtype_only_design(pd.DataFrame([row]), "AML")
        self.assertNotEqual(x0.tolist(), x1.tolist())


if __name__ == "__main__":
    unittest.main()
