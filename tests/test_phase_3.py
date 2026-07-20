from __future__ import annotations

import unittest

from src.phase_3 import format_p_value, standardized_difference


class StandardizedDifferenceTests(unittest.TestCase):
    def test_equal_proportions_have_no_difference(self) -> None:
        self.assertEqual(standardized_difference(0.5, 0.5), 0.0)

    def test_direction_is_sepsis_minus_nonsepsis(self) -> None:
        self.assertGreater(standardized_difference(0.6, 0.4), 0)
        self.assertLess(standardized_difference(0.4, 0.6), 0)

    def test_degenerate_equal_proportions_are_safe(self) -> None:
        self.assertEqual(standardized_difference(1.0, 1.0), 0.0)

    def test_p_value_format(self) -> None:
        self.assertEqual(format_p_value(0.0001), "<0.001")
        self.assertEqual(format_p_value(0.0126), "0.013")


if __name__ == "__main__":
    unittest.main()
