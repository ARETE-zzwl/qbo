"""Validate target-file structure and the frozen same-QBO treatment invariants."""
from __future__ import annotations

import argparse
import json
import os
import shutil
import tempfile
from pathlib import Path

from .io import sha256

import numpy as np
from netCDF4 import Dataset


def read_nc(path: Path):
    fd, temp_name = tempfile.mkstemp(suffix=".nc")
    os.close(fd)
    Path(temp_name).unlink()
    try:
        shutil.copy2(path, temp_name)
        with Dataset(temp_name) as ds:
            if ds.variables["qbo"].dimensions != ("time", "level"):
                raise ValueError("WACCM requires on-disk qbo(time,level) for its Fortran reader")
            levels = np.asarray(ds.variables["level"][:], dtype=float)
            dates = np.asarray(ds.variables["date"][:], dtype=int)
            secs = np.asarray(ds.variables["secs"][:], dtype=int)
            qbo = np.asarray(ds.variables["qbo"][:], dtype=float).T
            return levels, dates, secs, qbo
    finally:
        Path(temp_name).unlink(missing_ok=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("directory", type=Path)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    w_path = args.directory / "qbo_target_W.nc"
    e_path = args.directory / "qbo_target_E.nc"
    wl, wd, ws, w = read_nc(w_path)
    el, ed, es, e = read_nc(e_path)
    expected_dates = np.array([10901, 11001, 11101, 11201, 20101, 20201, 20301, 20401, 20501, 20601, 20701, 20801, 20901])
    checks = {
        "same_shape": bool(w.shape == e.shape == (7, 13)),
        "levels_equal": bool(np.array_equal(wl, el)),
        "levels_expected_top_to_bottom": bool(np.array_equal(wl, np.array([10, 20, 30, 40, 50, 70, 100], dtype=float))),
        "dates_equal": bool(np.array_equal(wd, ed)),
        "dates_expected": bool(np.array_equal(wd, expected_dates)),
        "secs_zero": bool(np.all(ws == 0) and np.all(es == 0)),
        "finite": bool(np.isfinite(w).all() and np.isfinite(e).all()),
        "march_70_separation_positive": bool(w[5, 6] > 0 and e[5, 6] < 0 and w[5, 6] - e[5, 6] > 10),
        "max_abs_below_relax_guard": bool(np.max(np.abs(w)) < 50 and np.max(np.abs(e)) < 50),
    }
    payload = {
        "scope": "file_format_and_target_values_only",
        "model_fidelity_verified": False,
        "checks": checks,
        "passed": all(checks.values()),
        "qbo_target_W_sha256": sha256(w_path),
        "qbo_target_E_sha256": sha256(e_path),
        "march_70_W_mps": float(w[5, 6]),
        "march_70_E_mps": float(e[5, 6]),
        "march_70_delta_mps": float(w[5, 6] - e[5, 6]),
        "max_abs_W_mps": float(np.max(np.abs(w))),
        "max_abs_E_mps": float(np.max(np.abs(e))),
    }
    args.output.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(json.dumps(payload, indent=2))
    if not payload["passed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
