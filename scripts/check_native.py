"""Build and run the standalone boundary, zonal and MPI contract checks."""
import argparse
import json
from pathlib import Path
import shutil
import subprocess

import numpy as np

from qbo.native import prepare, verify

ROOT = Path(__file__).resolve().parents[1]


def run(command, directory, name, check=True):
    result = subprocess.run(
        [str(part) for part in command], cwd=directory,
        capture_output=True, text=True, timeout=120,
    )
    (directory / f"{name}.stdout").write_text(result.stdout, encoding="utf-8")
    (directory / f"{name}.stderr").write_text(result.stderr, encoding="utf-8")
    if check:
        result.check_returncode()
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    compilers = {name: shutil.which(name) for name in ("gfortran", "mpifort", "mpiexec")}
    missing = [name for name, path in compilers.items() if path is None]
    if missing:
        parser.error("Add these tools to PATH: " + ", ".join(missing))
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    flags = ["-O0", "-fcheck=all", "-ffree-line-length-none"]
    sources = ROOT / "native"

    vectors = output / "boundary-vectors"
    prepare(vectors)
    boundary = output / "boundary"
    run([compilers["gfortran"], *flags,
         sources / "boundary/qbo_boundary_v4.f90",
         sources / "boundary/probe_boundary_v4.f90", "-o", boundary],
        output, "compile-boundary")
    run([boundary, vectors / "vectors.txt"], output, "boundary")
    verify(vectors, output / "boundary.stdout", output / "boundary.json")

    zonal = output / "zonal"
    run([compilers["mpifort"], *flags,
         sources / "zonal/probe_zonal_v5_mpi.f90", "-o", zonal],
        output, "compile-zonal")
    weights = np.array([[1., .75, .25, .5], [0., 1., 1., 1.], [1., 1., 1., 1.]])
    fixture = output / "rings.txt"
    np.savetxt(fixture, weights, header="3 4", comments="", fmt="%.17e")
    layouts = [(1, "cyclic"), (2, "cyclic"), (2, "block"), (2, "single")]
    for ranks, layout in layouts:
        name = f"zonal-{ranks}-{layout}"
        run([compilers["mpiexec"], "-n", ranks, zonal, fixture, layout], output, name)
        actual = np.loadtxt(output / f"{name}.stdout")
        np.testing.assert_array_equal(actual[:, 0], np.arange(1, 4))
        np.testing.assert_allclose(actual[:, 1], weights.min(axis=1), atol=5e-14, rtol=0)
    rejected = run([compilers["mpiexec"], "-n", 2, zonal, fixture, "incomplete"],
                   output, "zonal-incomplete", check=False)
    if rejected.returncode == 0 or "PARTITION_COVERAGE_REJECTED got=3 expected=4" not in rejected.stdout:
        raise AssertionError("Incomplete longitude coverage was not rejected.")

    contract = output / "contract"
    run([compilers["mpifort"], *flags,
         sources / "contract/qbo_contract.f90", sources / "contract/test_contract.f90",
         "-o", contract], output, "compile-contract")
    result = run([compilers["mpiexec"], "-n", 2, contract], output, "contract")
    expected = json.loads((sources / "contract/cases.json").read_text())["cases"]
    rows = [line.split() for line in result.stdout.splitlines() if line.startswith("CASE ")]
    if len(rows) != len(expected) or "CONTRACT_SYNTHETIC_PASS" not in result.stdout:
        raise AssertionError("The twelve contract cases did not complete.")
    for row, (case, _, code) in zip(rows, expected):
        if [int(value) for value in row[1:6]] != [case, code, code, code, 1]:
            raise AssertionError("Unexpected contract case: " + " ".join(row))
    summary = {
        "boundary": json.loads((output / "boundary.json").read_text()),
        "zonal_layouts": len(layouts), "zonal_incomplete_coverage_rejected": True,
        "contract_cases": len(rows), "passed": True,
    }
    (output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
