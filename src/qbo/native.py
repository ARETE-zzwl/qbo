"""Prepare deterministic kernel vectors or verify native output against them."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from .boundary import column_taper, safe_pressure, LOG_WIDTH, MARGIN_HPA


def prepare(directory):
    directory.mkdir()
    rng = np.random.default_rng(20260914)
    vectors = np.ones((10000, 7))
    vectors[:, 0] = rng.uniform(.1, 180, 10000)
    vectors[:, 1:4] = rng.uniform(50, 150, (10000, 3))
    extra = []
    for x in (-1, -1e-12, 0, 1e-12, .1, .5, .9, 1 - 1e-12, 1, 2):
        extra.append([100 * np.exp(-LOG_WIDTH * x) - MARGIN_HPA, 100, 100, 100, 1, 1, 1])
    for column in range(7):
        for invalid in (np.nan, np.inf, -np.inf, 0, -1):
            row = np.array([70, 100, 100, 100, 1, 1, 1], dtype=float)
            row[column] = invalid
            extra.append(row)
    extra.extend([[1e-300, 100, 100, 100, 1, 1, 1],
                  [1e300, 2e300, 2e300, 2e300, 1, 1, 1]])
    vectors = np.vstack([vectors, extra])
    safe = safe_pressure(vectors[:, 1:4].T, vectors[:, 4:7].T)
    expected = column_taper(vectors[:, 0], safe)
    np.savetxt(directory / "vectors.txt", vectors, fmt="%.17e", header=str(len(vectors)), comments="")
    np.savetxt(directory / "expected_weights.txt", expected, fmt="%.17e")
    manifest = {
        "seed": 20260914, "vectors": len(vectors),
        "absolute_tolerance": 5e-14,
        "files": {name: hashlib.sha256((directory / name).read_bytes()).hexdigest()
                  for name in ("vectors.txt", "expected_weights.txt")},
        "reference_sha256": hashlib.sha256(
            Path(__file__).with_name("boundary.py").read_bytes()).hexdigest(),
    }
    (directory / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(manifest, indent=2))


def verify(directory, native, output):
    if output.exists():
        raise FileExistsError("Use a new native assessment path.")
    manifest = json.loads((directory / "manifest.json").read_text())
    for name, expected_hash in manifest["files"].items():
        if hashlib.sha256((directory / name).read_bytes()).hexdigest() != expected_hash:
            raise ValueError("Test vector/reference file changed: " + name)
    expected = np.loadtxt(directory / "expected_weights.txt")
    values = np.loadtxt(native)
    if values.shape != (manifest["vectors"], 2) or not np.array_equal(values[:, 0], np.arange(1, len(expected) + 1)):
        raise ValueError("Native vector IDs/count do not match.")
    weight = values[:, 1]
    error = float(np.max(np.abs(weight - expected)))
    passed = (np.isfinite(weight).all() and np.all((weight >= 0) & (weight <= 1))
              and error <= manifest["absolute_tolerance"] and np.all(weight[expected == 0] == 0))
    result = {
        "scope": "standalone_native_kernel_not_full_CESM",
        "vectors": len(expected), "maximum_absolute_error": error,
        "exact_zero_vectors": int(np.count_nonzero(expected == 0)),
        "native_kernel_passed": bool(passed),
        "native_output_sha256": hashlib.sha256(native.read_bytes()).hexdigest(),
    }
    with output.open("x", encoding="utf-8") as handle:
        json.dump(result, handle, indent=2, allow_nan=False)
        handle.write("\n")
    print(json.dumps(result, indent=2))
    if not passed:
        raise ValueError("Native kernel disagrees with the frozen reference.")


def main():
    parser = argparse.ArgumentParser()
    commands = parser.add_subparsers(dest="command", required=True)
    create = commands.add_parser("prepare")
    create.add_argument("directory", type=Path)
    compare = commands.add_parser("verify")
    compare.add_argument("directory", type=Path)
    compare.add_argument("native", type=Path)
    compare.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.command == "prepare":
        prepare(args.directory)
    else:
        verify(args.directory, args.native, args.output)


if __name__ == "__main__":
    main()
