"""Preparation and source-scope contracts for the synchronous CAM adapter."""
import hashlib
import importlib.util
from pathlib import Path
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "prepare_cam_adapter_v4", ROOT / "scripts/prepare_cam_adapter.py"
)
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)
BASELINE = ROOT / "vendor/cam/qbo.F90"
KERNEL = ROOT / "native/boundary/qbo_boundary_v4.f90"


class AdapterPreparationTests(unittest.TestCase):
    def test_changed_baseline_is_rejected_before_patching(self):
        with self.assertRaisesRegex(ValueError, "baseline"):
            MODULE.build_source(BASELINE.read_bytes() + b"\n")

    def test_inherited_parameters_and_public_interface_are_unchanged(self):
        source = MODULE.build_source(BASELINE.read_bytes())
        for text in (
            "tauz(ktop-1)    = 2._r8",
            "tauz(kbot+1)    = 2._r8",
            "tauz(ktop:kbot) = 1._r8",
            "tconst = 10._r8",
            "subroutine qbo_relax( state, pbuf, ptend )",
            "qbo_use_forcing  = .FALSE.",
        ):
            self.assertIn(text, source)
        original = BASELINE.read_text(encoding="utf-8")
        for name in ("qbo_readnl", "qbo_timestep_init", "taux"):
            self.assertEqual(
                MODULE.procedure(source, name), MODULE.procedure(original, name)
            )

    def test_adapter_runs_only_inside_enabled_forcing_branch(self):
        body = MODULE.procedure(MODULE.build_source(BASELINE.read_bytes()), "qbo_relax")
        self.assertLess(
            body.index("if( qbo_use_forcing ) then"),
            body.index("call qbo_boundary_weights"),
        )
        self.assertNotIn("tropopause_output", body)

    def test_tendency_and_diagnostics_share_the_frozen_weight(self):
        body = MODULE.procedure(MODULE.build_source(BASELINE.read_bytes()), "qbo_relax")
        self.assertIn("crelax = boundary_weight(i,k) / crelax", body)
        self.assertIn("if(u < 50.0_r8 .and. crelax > 0._r8)", body)
        self.assertIn("qbo_u0(i,k) = u/tauzz/tauxi(i)*tconst1 * boundary_weight(i,k)", body)
        self.assertIn("qbo_rate(i,k) = crelax", body)
        self.assertIn("qbo_mask(i,k) = boundary_weight(i,k)", body)

    def test_preparation_is_flat_and_preserves_the_original_kernel(self):
        before = hashlib.sha256(BASELINE.read_bytes()).hexdigest()
        with tempfile.TemporaryDirectory() as temporary:
            destination = Path(temporary) / "prepared"
            MODULE.prepare(destination)
            flat = destination / "SourceMods/src.cam"
            self.assertEqual(
                sorted(p.name for p in flat.iterdir()),
                ["qbo.F90", "qbo_boundary_cam_v4.F90", "qbo_boundary_v4.F90"],
            )
            self.assertEqual((flat / "qbo_boundary_v4.F90").read_bytes(), KERNEL.read_bytes())
            self.assertEqual(before, hashlib.sha256(BASELINE.read_bytes()).hexdigest())

    def test_existing_destination_is_not_overwritten(self):
        with tempfile.TemporaryDirectory() as temporary:
            destination = Path(temporary)
            sentinel = destination / "keep.txt"
            sentinel.write_text("old material", encoding="ascii")
            with self.assertRaises(FileExistsError):
                MODULE.prepare(destination)
            self.assertEqual(sentinel.read_text(), "old material")
            self.assertEqual(list(destination.iterdir()), [sentinel])

    def test_native_envelope_contains_the_exact_patched_procedure_bodies(self):
        with tempfile.TemporaryDirectory() as temporary:
            destination = Path(temporary) / "prepared"
            MODULE.prepare(destination)
            actual = (destination / "SourceMods/src.cam/qbo.F90").read_text(encoding="utf-8")
            harness = (destination / "native/qbo_test.F90").read_text(encoding="utf-8")
            for name in ("qbo_relax", "taux"):
                self.assertEqual(
                    MODULE.procedure(actual, name), MODULE.procedure(harness, name)
                )


if __name__ == "__main__":
    unittest.main()
