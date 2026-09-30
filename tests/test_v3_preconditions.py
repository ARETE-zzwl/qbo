"""Reject the observed forcing-off case and screen the proposed hard cutoff."""
import unittest

import numpy as np

from qbo import preconditions as MODULE


class NamelistGateTests(unittest.TestCase):
    def setUp(self):
        self.values = {
            "qbo_use_forcing": [".true."], "qbo_cyclic": [".false."],
            "qbo_forcing_file": ["'/inputs/qbo_target_W.nc'"],
            "nhtfrq": ["1"], "mfilt": ["2"], "avgflag_pertape": ["'I'"],
            "empty_htapes": [".true."],
            "fincl1": ["'" + name + "'" for name in MODULE.FIELDS],
        }

    def test_explicit_valid_diagnostic_configuration_passes(self):
        self.assertEqual(MODULE.namelist_failures(self.values, "/inputs/qbo_target_W.nc"), [])

    def test_completed_but_forcing_disabled_is_rejected(self):
        self.values["qbo_use_forcing"] = [".false."]
        self.assertIn("forcing_disabled", MODULE.namelist_failures(
            self.values, "/inputs/qbo_target_W.nc"
        ))

    def test_empty_namelist_fails_closed(self):
        self.assertIn("forcing_disabled", MODULE.namelist_failures(
            {}, "/inputs/qbo_target_W.nc"
        ))

    def test_wrong_phase_target_is_rejected(self):
        self.assertIn("wrong_target", MODULE.namelist_failures(
            self.values, "/inputs/qbo_target_E.nc"
        ))

    def test_default_history_frequency_is_rejected(self):
        self.values["nhtfrq"] = ["0", "-24", "-6"]
        self.assertIn("history_frequency", MODULE.namelist_failures(
            self.values, "/inputs/qbo_target_W.nc"
        ))

    def test_missing_tendency_diagnostic_is_rejected(self):
        self.values["fincl1"].remove("'QBOTEND'")
        self.assertIn("missing_diagnostics", MODULE.namelist_failures(
            self.values, "/inputs/qbo_target_W.nc"
        ))


class HardCutoffTests(unittest.TestCase):
    def test_removes_both_buffers_without_changing_interior(self):
        np.testing.assert_array_equal(
            MODULE.hard_cutoff_strength([5, 9, 11, 70, 95, 110, 150], [10, 100]),
            [0, 0, 1, 1, 1, 0, 0],
        )

    def test_zero_buffers_do_not_make_interior_troposphere_safe(self):
        result = MODULE.cutoff_overlap(
            [5, 50, 95, 110], [1, 10, 70, 100, 120],
            np.array([[80.0]]), np.ones((1, 1), bool), [1.0],
        )
        self.assertEqual(result["midpoint_below_fraction"], 1)
        self.assertEqual(result["interface_below_fraction"], 1)

    def test_missing_tropopause_cannot_be_a_safety_pass(self):
        result = MODULE.cutoff_overlap(
            [5, 50, 95, 110], [1, 10, 70, 100, 120],
            np.array([[np.nan]]), np.zeros((1, 1), bool), [1.0],
        )
        self.assertIsNone(result["midpoint_below_fraction"])
        self.assertEqual(result["valid_area_fraction"], 0)


if __name__ == "__main__":
    unittest.main()
