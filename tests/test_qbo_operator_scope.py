"""Analytic controls for the original-intervention scope audit."""
import importlib.util
from pathlib import Path
import sys
import tempfile
import unittest

import numpy as np

SCRIPTS = Path(__file__).resolve().parents[1] / "scripts/analysis"
sys.path.insert(0, str(SCRIPTS))
SPEC = importlib.util.spec_from_file_location(
    "audit_qbo_operator_scope", SCRIPTS / "audit_qbo_operator_scope.py"
)
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class OperatorProjectionTests(unittest.TestCase):
    def test_zonal_uniform_control_has_no_eddy_component(self):
        result = MODULE.operator_metrics([[1, 1], [2, 2]], [1, 3])
        self.assertEqual(result["area_mean_rate_per_s"], 1.75)
        self.assertEqual(result["eddy_rms_per_s"], 0)
        self.assertEqual(result["eddy_squared_amplitude_fraction"], 0)

    def test_column_mask_creates_expected_non_zonal_component(self):
        result = MODULE.operator_metrics([[0, 2]], [1])
        self.assertEqual(result["zonal_rms_per_s"], 1)
        self.assertEqual(result["eddy_rms_per_s"], 1)
        self.assertEqual(result["eddy_to_zonal_rms_ratio"], 1)
        self.assertEqual(result["eddy_squared_amplitude_fraction"], .5)

    def test_weights_act_on_squared_components_before_aggregation(self):
        result = MODULE.operator_metrics([[0, 2], [1, 1]], [1, 3])
        self.assertEqual(result["eddy_rms_per_s"], .5)
        self.assertEqual(result["eddy_to_zonal_rms_ratio"], .5)
        self.assertEqual(result["eddy_squared_amplitude_fraction"], .2)
        self.assertEqual(result["squared_norm_relative_error"], 0)

    def test_zero_operator_is_retained_with_undefined_ratios(self):
        result = MODULE.operator_metrics(np.zeros((2, 3)), [1, 2])
        self.assertEqual(result["eddy_rms_per_s"], 0)
        self.assertIsNone(result["eddy_to_zonal_rms_ratio"])
        self.assertIsNone(result["eddy_squared_amplitude_fraction"])

    def test_reconstruction_measures_roundoff_after_components_are_added(self):
        result = MODULE.operator_metrics([[1e-30, 1e-6]], [1])
        self.assertEqual(result["maximum_reconstruction_error_per_s"], 1e-30)

    def test_uniform_error_scaling_does_not_change_spatial_fraction(self):
        original = MODULE.operator_metrics([[0, 2]], [1])
        doubled = MODULE.operator_metrics([[0, 4]], [1])
        self.assertEqual(doubled["eddy_rms_per_s"], 2 * original["eddy_rms_per_s"])
        self.assertEqual(
            doubled["eddy_squared_amplitude_fraction"],
            original["eddy_squared_amplitude_fraction"],
        )

    def test_nonfinite_rate_is_not_silently_excluded(self):
        with self.assertRaisesRegex(ValueError, "finite"):
            MODULE.operator_metrics([[0, np.nan]], [1])

    def test_existing_assessment_is_never_overwritten(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "assessment.json"
            path.write_text("old evidence", encoding="ascii")
            with self.assertRaises(FileExistsError):
                MODULE.write_assessment(path)
            self.assertEqual(path.read_text(), "old evidence")


if __name__ == "__main__":
    unittest.main()
