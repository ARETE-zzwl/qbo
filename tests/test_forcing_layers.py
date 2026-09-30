"""Check the WACCM forcing support against hand-calculated small grids."""
import unittest

import numpy as np

from qbo import forcing as MODULE
forcing_layers = MODULE.forcing_layers


class ForcingLayerTests(unittest.TestCase):
    def test_includes_half_strength_layers_outside_target_range(self):
        rows = forcing_layers([5, 9, 11, 70, 95, 110, 150], [10, 100])
        self.assertEqual([row["index"] for row in rows], [2, 3, 4, 5, 6])
        self.assertEqual(
            [row["relative_strength"] for row in rows], [.5, 1, 1, 1, .5]
        )
        self.assertEqual(rows[-1]["pressure_hpa"], 110)

    def test_pressure_boundary_is_inclusive_for_layer_selection(self):
        rows = forcing_layers([5, 10, 70, 100, 150], [10, 100])
        self.assertEqual(
            [row["relative_strength"] for row in rows], [.5, 1, 1, 1, .5]
        )

    def test_buffer_is_clipped_at_model_top_and_bottom(self):
        rows = forcing_layers([10, 50, 100], [10, 100])
        self.assertEqual([row["index"] for row in rows], [1, 2, 3])
        self.assertTrue(all(row["relative_strength"] == 1 for row in rows))

    def test_rejects_reversed_or_empty_support(self):
        with self.assertRaises(ValueError):
            forcing_layers([100, 50, 10], [10, 100])
        with self.assertRaises(ValueError):
            forcing_layers([1, 2, 3], [10, 100])


class NativeValueTests(unittest.TestCase):
    def setUp(self):
        self.rows = np.array([
            [1, 1, 10901, 0, 10, 1.25],
            [1, 2, 10901, 0, 100, -2.5],
            [2, 1, 11001, 0, 10, 3.75],
            [2, 2, 11001, 0, 100, -4.5],
        ])
        self.arguments = (
            [10, 100], [10901, 11001], [0, 0],
            [[1.25, -2.5], [3.75, -4.5]],
        )

    def test_every_native_value_and_coordinate_matches(self):
        result = MODULE.check_native_values(self.rows, *self.arguments)
        self.assertEqual(result["values_checked"], 4)
        self.assertEqual(result["max_abs_error_mps"], 0)

    def test_wrong_native_value_is_rejected(self):
        self.rows[2, -1] += .25
        with self.assertRaisesRegex(ValueError, "wind"):
            MODULE.check_native_values(self.rows, *self.arguments)

    def test_wrong_time_coordinate_is_rejected(self):
        self.rows[2, 2] = 10901
        with self.assertRaisesRegex(ValueError, "coordinate"):
            MODULE.check_native_values(self.rows, *self.arguments)


class RuntimeDiagnosticTests(unittest.TestCase):
    def test_buffer_repeats_edge_wind_before_half_strength_weighting(self):
        values = MODULE.expected_qbo_diagnostic(
            [5, 9, 11, 70, 95, 110, 150], [0, 10, 30], [10, 100], [0, 9]
        )
        np.testing.assert_allclose(values[:, 0], [0, .05, .1, 6, 8.5, 4.25, 0])
        np.testing.assert_allclose(values[:, 2], 0)
        self.assertAlmostEqual(values[3, 1] / values[3, 0], .60653066, places=7)

    def test_targets_at_or_above_50_follow_waccm_skip_rule(self):
        values = MODULE.expected_qbo_diagnostic(
            [5, 20, 80, 120], [0], [10, 100], [50, 60]
        )
        np.testing.assert_array_equal(values, np.zeros((4, 1)))

    def test_tropopause_overlap_uses_area_and_reports_missing_coverage(self):
        trop = np.array([[80, 120], [106, 90]])
        valid = np.ones((2, 2), dtype=bool)
        active = np.ones((2, 2), dtype=bool)
        result = MODULE.tropopause_overlap(trop, valid, [1, 3], 103, 111, active)
        self.assertEqual(result["valid_area_fraction"], 1)
        self.assertEqual(result["midpoint_below_fraction"], .5)
        self.assertEqual(result["interface_below_fraction"], .875)
        valid[1, 0] = False
        result = MODULE.tropopause_overlap(trop, valid, [1, 3], 103, 111, active)
        self.assertEqual(result["valid_area_fraction"], .625)
        self.assertEqual(result["midpoint_below_fraction"], .8)

    def test_missing_tropopause_does_not_become_a_zero_overlap_claim(self):
        result = MODULE.tropopause_overlap(
            np.ones((2, 2)), np.zeros((2, 2), bool), [1, 3], 103, 111,
            np.ones((2, 2), bool),
        )
        self.assertEqual(result["valid_area_fraction"], 0)
        self.assertIsNone(result["midpoint_below_fraction"])

    def test_float32_pressure_roundoff_is_not_a_below_tropopause_signal(self):
        result = MODULE.tropopause_overlap(
            np.array([[102.999999]]), np.ones((1, 1), bool), [1], 103, 111,
            np.ones((1, 1), bool),
        )
        self.assertEqual(result["midpoint_below_fraction"], 0)


if __name__ == "__main__":
    unittest.main()
