from __future__ import annotations

import unittest

import numpy as np

from src.phase_6_palliative_care import comparison, weighted_prevalence


class PalliativeCareSurveyTests(unittest.TestCase):
    def test_equal_weights_reproduce_sample_prevalence(self) -> None:
        result = weighted_prevalence(
            np.ones(6), np.array([0, 0, 1, 0, 1, 1]), np.array([1, 1, 1, 2, 2, 2])
        )
        self.assertAlmostEqual(result["estimate"], 0.5)
        self.assertEqual(result["unweighted_events"], 3)

    def test_domain_counts_and_weighted_estimate(self) -> None:
        result = weighted_prevalence(
            np.array([1, 1, 2, 2]), np.array([0, 1, 0, 1]),
            np.array([1, 1, 2, 2]), np.array([False, True, True, True]),
        )
        self.assertEqual(result["unweighted_n"], 3)
        self.assertAlmostEqual(result["estimate"], 3 / 5)

    def test_comparison_odds_ratio(self) -> None:
        result = comparison(
            np.ones(8), np.array([0, 0, 0, 1, 0, 1, 1, 1]),
            np.array([1, 1, 2, 2, 1, 1, 2, 2]),
            np.array([False, False, False, False, True, True, True, True]),
        )
        self.assertAlmostEqual(result["odds_ratio"], 9.0)


if __name__ == "__main__":
    unittest.main()
