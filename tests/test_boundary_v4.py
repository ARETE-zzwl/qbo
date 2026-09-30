"""Frozen safety and continuity checks for the column-wise QBO mask."""
import unittest

import numpy as np

from qbo import boundary as MODULE


class SafePressureTests(unittest.TestCase):
    def test_uses_highest_of_all_three_tropopauses(self):
        pressure = [[100, 95], [100, 90], [80, 85]]
        np.testing.assert_array_equal(
            MODULE.safe_pressure(pressure, np.ones((3, 2))), [80, 85]
        )

    def test_missing_primary_suppresses_column_even_with_fallback(self):
        found = np.ones((3, 2))
        found[1, 0] = 0
        result = MODULE.safe_pressure(np.full((3, 2), 100.0), found)
        self.assertTrue(np.isnan(result[0]))
        self.assertEqual(result[1], 100)

    def test_masked_nonfinite_and_nonpositive_pressures_fail_closed(self):
        pressure = np.ma.array(np.full((3, 5), 100.0))
        pressure[0, 0] = np.nan
        pressure[1, 1] = np.inf
        pressure[2, 2] = -1
        pressure[0, 3] = 0
        pressure[1, 4] = np.ma.masked
        self.assertTrue(np.isnan(MODULE.safe_pressure(pressure, np.ones((3, 5)))).all())

    def test_nonfinite_found_flag_cannot_be_true(self):
        found = np.array([[1], [np.nan], [1]])
        self.assertTrue(np.isnan(MODULE.safe_pressure([[100], [100], [100]], found)[0]))

    def test_three_definitions_are_required(self):
        with self.assertRaises(ValueError):
            MODULE.safe_pressure(np.ones((2, 3)), np.ones((2, 3)))


class ColumnTaperTests(unittest.TestCase):
    def test_layer_bottom_not_midpoint_controls_safety(self):
        np.testing.assert_array_equal(
            MODULE.column_taper([80, 94.941, 111.693], 90), [MODULE.column_taper(80, 90), 0, 0]
        )

    def test_pressure_margin_and_equality_have_exact_zero(self):
        np.testing.assert_array_equal(
            MODULE.column_taper([100, 99.999, 110], 100), [0, 0, 0]
        )

    def test_half_transition_has_half_weight(self):
        bottom = 100 * np.exp(-MODULE.LOG_WIDTH / 2) - MODULE.MARGIN_HPA
        self.assertAlmostEqual(float(MODULE.column_taper(bottom, 100)), .5, places=13)

    def test_saturated_upper_region_is_unchanged(self):
        np.testing.assert_array_equal(MODULE.column_taper([10, 30, 50], 100), [1, 1, 1])

    def test_weights_are_finite_bounded_and_monotone(self):
        result = MODULE.column_taper(np.linspace(.1, 150, 10000), 100)
        self.assertTrue(np.isfinite(result).all())
        self.assertTrue(((result >= 0) & (result <= 1)).all())
        self.assertTrue((np.diff(result) <= 0).all())

    def test_invalid_inputs_have_zero_weight_without_warnings(self):
        with np.errstate(all="raise"):
            np.testing.assert_array_equal(
                MODULE.column_taper([np.nan, np.inf, 0, -1, 70], [100, 100, 100, 100, np.nan]),
                np.zeros(5),
            )

    def test_slope_tends_to_zero_at_both_transition_ends(self):
        delta = 1e-6
        bottom = 100 * np.exp(-MODULE.LOG_WIDTH * np.array([delta, 1 - delta])) - MODULE.MARGIN_HPA
        weight = MODULE.column_taper(bottom, 100)
        self.assertLess(float(weight[0] / delta), 1e-5)
        self.assertLess(float((1 - weight[1]) / delta), 1e-5)

    def test_missing_column_is_zero_at_every_level(self):
        result = MODULE.column_taper(np.array([20, 70, 100])[:, None], [np.nan, 90])
        np.testing.assert_array_equal(result[:, 0], 0)
        self.assertEqual(result[0, 1], 1)

    def test_explicit_scalar_step_is_contractive_and_never_amplified(self):
        multiplier = MODULE.column_taper(np.linspace(1, 120, 100), 100)
        coefficient = multiplier * (1800 / (10 * 86400))
        self.assertLessEqual(float(coefficient.max()), 1800 / (10 * 86400))
        for initial, target in ((-60, 30), (60, -30), (0, 0)):
            updated = initial + coefficient * (target - initial)
            self.assertTrue((updated >= min(initial, target)).all())
            self.assertTrue((updated <= max(initial, target)).all())


if __name__ == "__main__":
    unittest.main()
