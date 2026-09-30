"""Prepare isolated, hash-pinned CAM source; never edit an existing case."""
import argparse
import hashlib
import json
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[1]
BASELINE = ROOT / "vendor/cam/qbo.F90"
BASELINE_SHA256 = "7e098beb9455d7c45a6a63ffc2f457d29a3d9e971cf9f1643b748a59a0212a64"
KERNEL = ROOT / "native/boundary/qbo_boundary_v4.f90"
KERNEL_SHA256 = "92c519139c0883414cc326dfd7493d34a2c8d9e92b1cfd302a085ddab8fd39aa"


def procedure(source, name):
    pattern = rf"(?mi)^[ \t]*(subroutine|function)[ \t]+{name}\b"
    starts = list(re.finditer(pattern, source))
    if len(starts) != 1:
        raise ValueError("Procedure is not unique: " + name)
    start = starts[0]
    end = re.search(
        rf"(?mi)^[ \t]*end[ \t]+{start[1]}[ \t]+{name}\b[^\r\n]*",
        source[start.end():],
    )
    if end is None:
        raise ValueError("Procedure end not found: " + name)
    return source[start.start():start.end() + end.end()]


def build_source(baseline):
    if hashlib.sha256(baseline).hexdigest() != BASELINE_SHA256:
        raise ValueError("Unexpected QBO baseline hash.")
    source = baseline.decode("utf-8")
    fields = [
        "    call addfld('QBOMASK', (/ 'lev' /), 'I', '1', 'Applied QBO boundary weight')",
        "    call addfld('QBO_RATE', (/ 'lev' /), 'I', '1/s', 'Applied QBO relaxation rate')",
    ]
    for name in ("QBO_TROP", "QBO_TRPP", "QBO_TRPF"):
        fields.append(f"    call addfld('{name}', horiz_only, 'I', 'Pa', 'Pre-QBO tropopause pressure')")
    for name in ("QBO_FD", "QBO_FDP", "QBO_FDF"):
        fields.append(f"    call addfld('{name}', horiz_only, 'I', '1', 'Pre-QBO tropopause found')")
    changes = [
        ("use cam_history,  only : addfld, add_default",
         "use cam_history,  only : addfld, add_default, horiz_only"),
        ("    call addfld ('QBO_U0', (/ 'lev' /), 'A','M/S','Specified wind used for QBO')",
         "    call addfld ('QBO_U0', (/ 'lev' /), 'A','M/S','Specified wind used for QBO')\n"
         + "\n".join(fields)),
        ("      use cam_history,    only: outfld",
         "      use cam_history,    only: outfld\n"
         "      use qbo_boundary_cam_v4, only: qbo_boundary_weights"),
        ("    real(r8), pointer                  :: uzm(:,:)",
         "    real(r8), pointer                  :: uzm(:,:)\n"
         "    real(r8) :: boundary_weight(pcols,pver), qbo_mask(pcols,pver), qbo_rate(pcols,pver)\n"
         "    real(r8) :: boundary_pressure(pcols,3), boundary_found(pcols,3)"),
        ("    qbo_u0(:,:) = 0._r8",
         "    call qbo_boundary_weights(state, boundary_weight, boundary_pressure, boundary_found)\n"
         "    qbo_u0(:,:) = 0._r8\n"
         "    qbo_mask(:,:) = 0._r8\n"
         "    qbo_rate(:,:) = 0._r8"),
        ("            crelax = 1._r8 / crelax",
         "            crelax = boundary_weight(i,k) / crelax"),
        ("         if(u < 50.0_r8) then",
         "         if(u < 50.0_r8 .and. crelax > 0._r8) then"),
        ("            qbo_u0(i,k) = u/tauzz/tauxi(i)*tconst1",
         "            qbo_u0(i,k) = u/tauzz/tauxi(i)*tconst1 * boundary_weight(i,k)\n"
         "            qbo_mask(i,k) = boundary_weight(i,k)\n"
         "            qbo_rate(i,k) = crelax"),
        ("    call outfld( 'QBO_U0', qbo_u0, pcols, lchnk )",
         "    call outfld( 'QBO_U0', qbo_u0, pcols, lchnk )\n"
         "    call outfld('QBOMASK', qbo_mask, pcols, lchnk)\n"
         "    call outfld('QBO_RATE', qbo_rate, pcols, lchnk)\n"
         + "\n".join(
             f"    call outfld('{name}', {array}(:,{index}), pcols, lchnk)"
             for array, names in (
                 ("boundary_pressure", ("QBO_TROP", "QBO_TRPP", "QBO_TRPF")),
                 ("boundary_found", ("QBO_FD", "QBO_FDP", "QBO_FDF")),
             )
             for index, name in enumerate(names, 1)
         )),
    ]
    for old, new in changes:
        if source.count(old) != 1:
            raise ValueError("QBO source anchor is not unique: " + old)
        source = source.replace(old, new, 1)
    return source


def prepare(directory):
    if directory.exists():
        raise FileExistsError("Use a new CAM adapter directory.")
    baseline = BASELINE.read_bytes()
    source = build_source(baseline)
    kernel = KERNEL.read_bytes()
    if hashlib.sha256(kernel).hexdigest() != KERNEL_SHA256:
        raise ValueError("The frozen v4 kernel changed.")
    envelope = """module qbo_test
  use shr_kind_mod, only: r8 => shr_kind_r8
  use ppgrid, only: pcols, pver
  use physics_types, only: physics_state, physics_ptend, physics_ptend_init
  use phys_grid, only: get_rlat_p
  implicit none
  integer :: ktop=2, kbot=pver-1, uzm_idx=1
  logical :: qbo_use_forcing=.true.
  real(r8) :: tauz(pver)=1._r8, u_tstep(pver)=0._r8
contains
"""
    baseline_text = baseline.decode("utf-8")
    negative = envelope + procedure(baseline_text, "qbo_relax") + "\n" + procedure(baseline_text, "taux")
    negative += "\nend module qbo_test\n"
    envelope += procedure(source, "qbo_relax") + "\n" + procedure(source, "taux")
    envelope += "\nend module qbo_test\n"
    files = {
        "SourceMods/src.cam/qbo.F90": source.encode("utf-8"),
        "SourceMods/src.cam/qbo_boundary_v4.F90": kernel,
        "SourceMods/src.cam/qbo_boundary_cam_v4.F90": (
            ROOT / "native/cam_adapter/qbo_boundary_cam_v4.F90"
        ).read_bytes(),
        "native/qbo_test.F90": envelope.encode("utf-8"),
        "native/qbo_baseline_test.F90": negative.encode("utf-8"),
    }
    support = ROOT / "native/cam_adapter"
    for name in ("cam_test_support.F90", "probe_cam_adapter_v4.F90"):
        files["native/" + name] = (support / name).read_bytes()
    directory.mkdir(parents=True)
    for relative, content in files.items():
        path = directory / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("xb") as handle:
            handle.write(content)
    manifest = {
        "scope": "source_adapter_and_instrumented_test_not_full_CESM",
        "baseline_sha256": BASELINE_SHA256,
        "frozen_kernel_sha256": KERNEL_SHA256,
        "files": {name: hashlib.sha256(content).hexdigest() for name, content in files.items()},
        "preparation_script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "full_model_run_approved": False,
    }
    files["manifest.json"] = (json.dumps(manifest, indent=2) + "\n").encode("utf-8")
    with (directory / "manifest.json").open("xb") as handle:
        handle.write(files["manifest.json"])
    with (directory / "SHA256SUMS").open("x", encoding="ascii", newline="\n") as handle:
        for name, content in sorted(files.items()):
            handle.write(f"{hashlib.sha256(content).hexdigest()}  {name}\n")
    return manifest


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("directory", type=Path)
    args = parser.parse_args()
    print(json.dumps(prepare(args.directory), indent=2))
