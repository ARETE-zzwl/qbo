"""Analytic controls for the reviewed longitude-uniform envelope."""
import unittest

import numpy as np

from qbo import zonal as MODULE

from qbo.boundary import column_taper, safe_pressure


class ZonalEnvelopeTests(unittest.TestCase):
    def test_uniform_weight_is_unchanged(self):
        actual = MODULE.zonal_envelope([[.3, .3, .3], [1, 1, 1]])
        np.testing.assert_array_equal(actual, [[.3], [1]])

    def test_minimum_is_the_largest_pointwise_dominated_uniform_weight(self):
        columns = np.array([[.2, .8, .6]])
        weight = MODULE.zonal_envelope(columns)
        np.testing.assert_array_equal(weight, [[.2]])
        self.assertTrue(np.all(weight <= columns))
        self.assertTrue(np.any(weight + .001 > columns))

    def test_one_unsafe_column_suppresses_the_entire_ring(self):
        columns = np.array([[0, .8, 1]])
        self.assertEqual(MODULE.zonal_envelope(columns)[0, 0], 0)
        self.assertGreater(columns.mean(), 0)

    def test_missing_definition_cannot_be_averaged_away(self):
        pressures = np.full((3, 1, 3), 100.)
        pressures[0, 0, 1] = np.nan
        weights = column_taper(50, safe_pressure(pressures, np.ones_like(pressures)))
        np.testing.assert_array_equal(MODULE.zonal_envelope(weights), [[0]])

    def test_equality_margin_zero_propagates_to_ring(self):
        weights = column_taper([[80, 50]], [[80.001, 100]])
        np.testing.assert_array_equal(MODULE.zonal_envelope(weights), [[0]])

    def test_levels_and_latitudes_remain_independent(self):
        weights = np.array([[[.2, .3], [.8, .9]], [[.5, .4], [.7, .6]]])
        np.testing.assert_array_equal(
            MODULE.zonal_envelope(weights), [[[.2], [.8]], [[.4], [.6]]]
        )

    def test_longitude_order_does_not_change_the_envelope(self):
        weights = np.random.default_rng(20260914).uniform(size=(2, 3, 7))
        np.testing.assert_array_equal(
            MODULE.zonal_envelope(weights),
            MODULE.zonal_envelope(weights[..., [6, 0, 3, 2, 4, 1, 5]]),
        )

    def test_no_amplification_preserves_scalar_relaxation_bound(self):
        weights = np.random.default_rng(5).uniform(size=(4, 7))
        rate = MODULE.zonal_envelope(weights) / 864000
        self.assertTrue(np.all(1800 * rate <= 1800 / 864000))

    def test_invalid_weights_are_not_silently_ignored(self):
        for weights in ([[np.nan, 1]], [[np.inf, 1]], [[-.1, 1]], [[1.1, 1]]):
            with self.subTest(weights=weights):
                with self.assertRaises(ValueError):
                    MODULE.zonal_envelope(weights)

    def test_absent_longitude_axis_is_rejected(self):
        for weights in (np.zeros((2, 0)), np.array(1.)):
            with self.subTest(shape=weights.shape):
                with self.assertRaises(ValueError):
                    MODULE.zonal_envelope(weights)


if __name__ == "__main__":
    unittest.main()
