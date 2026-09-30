"""Hand-calculated area/rate weighting controls for the offline assessment."""
import unittest

import numpy as np

from qbo import metrics as MODULE


class RetentionMetricsTests(unittest.TestCase):
    def test_distinguishes_area_fraction_from_retained_relaxation_rate(self):
        mask = np.array([[1, 0], [.5, .5]])
        result = MODULE.retention_metrics(mask, mask > 0, [1, .5], [1, 3], 1)
        self.assertEqual(result["positive_area_fraction"], .875)
        self.assertEqual(result["retained_rate_fraction"], .5)
        self.assertEqual(result["hard_mask_upper_bound"], .8)
        self.assertEqual(result["equivalent_relaxation_days"], 32)
        self.assertAlmostEqual(result["longitude_modulation_rms_ratio"], np.sqrt(4 / 7))

    def test_missing_control_is_reported_as_zero_not_a_finite_timescale(self):
        result = MODULE.retention_metrics(
            np.zeros((2, 2)), np.zeros((2, 2), bool), [1, 1], [1, 1], .5
        )
        self.assertEqual(result["retained_rate_fraction"], 0)
        self.assertIsNone(result["equivalent_relaxation_days"])
        self.assertIsNone(result["longitude_modulation_rms_ratio"])

    def test_unchanged_half_strength_buffer_keeps_its_twenty_day_timescale(self):
        result = MODULE.retention_metrics(
            np.ones((1, 2)), np.ones((1, 2), bool), [1], [1], .5
        )
        self.assertEqual(result["retained_rate_fraction"], 1)
        self.assertEqual(result["equivalent_relaxation_days"], 20)
        self.assertEqual(result["longitude_modulation_rms_ratio"], 0)


if __name__ == "__main__":
    unittest.main()
