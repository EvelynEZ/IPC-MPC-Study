from __future__ import annotations

import unittest

import numpy as np

from src.phase_7_subtype_palliative_care import overall_subtype_test, stratified_covariance


class SubtypeComparisonTests(unittest.TestCase):
    def test_covariance_is_symmetric(self) -> None:
        influence = np.array([[1.0, 0.0], [0.0, 1.0], [2.0, 1.0], [1.0, 2.0]])
        covariance = stratified_covariance(influence, np.array([1, 1, 2, 2]))
        np.testing.assert_allclose(covariance, covariance.T)

    def test_equal_estimates_have_zero_wald_statistic(self) -> None:
        influence = np.array([
            [0.1, -0.1, 0.0], [-0.1, 0.1, 0.0],
            [0.0, 0.1, -0.1], [0.0, -0.1, 0.1],
        ])
        result = overall_subtype_test(
            np.array([0.2, 0.2, 0.2]), influence, np.array([1, 1, 2, 2])
        )
        self.assertAlmostEqual(result["wald_chi_square"], 0.0)
        self.assertAlmostEqual(result["p_value"], 1.0)


if __name__ == "__main__":
    unittest.main()
