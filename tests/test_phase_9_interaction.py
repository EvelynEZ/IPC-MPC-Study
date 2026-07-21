from __future__ import annotations

import unittest

import numpy as np
import pandas as pd

from src.phase_8_primary_adjusted import CATEGORY_SPECS
from src.phase_9_interaction import SUBTYPE_LEVELS, interaction_design_matrix


class InteractionDesignTests(unittest.TestCase):
    def reference_frame(self) -> pd.DataFrame:
        row = {"age": 65, "sepsis": 0, **{variable: levels[0] for variable, levels in CATEGORY_SPECS.items()}}
        return pd.DataFrame([row])

    def test_reference_lymphoma_has_zero_interactions(self) -> None:
        x, names = interaction_design_matrix(self.reference_frame(), force_sepsis=1, force_subtype="Lymphoma")
        indices = [i for i, name in enumerate(names) if name.startswith("Sepsis ×")]
        self.assertTrue(np.all(x[0, indices] == 0))

    def test_nonreference_subtype_activates_one_interaction(self) -> None:
        x, names = interaction_design_matrix(self.reference_frame(), force_sepsis=1, force_subtype=SUBTYPE_LEVELS[1])
        indices = [i for i, name in enumerate(names) if name.startswith("Sepsis ×")]
        self.assertEqual(int(x[0, indices].sum()), 1)

    def test_no_sepsis_has_zero_interactions(self) -> None:
        x, names = interaction_design_matrix(self.reference_frame(), force_sepsis=0, force_subtype=SUBTYPE_LEVELS[-1])
        indices = [i for i, name in enumerate(names) if name.startswith("Sepsis ×")]
        self.assertTrue(np.all(x[0, indices] == 0))


if __name__ == "__main__":
    unittest.main()
