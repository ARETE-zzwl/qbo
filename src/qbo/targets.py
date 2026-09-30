"""Build deterministic full-profile WQBO/EQBO monthly targets for the CESM QBO module.

The only selection variable is March 70-hPa zonal-wind anomaly.  Targets are
12-month composites (September before the selected March through the following
September endpoint) from the fixed upper/lower eight-year tails in 1981--2024.
"""
from __future__ import annotations

import argparse
import csv
import json
import re
import shutil
import tempfile
from pathlib import Path

from .io import sha256

import numpy as np
from netCDF4 import Dataset


START_YEAR = 1981
END_YEAR = 2024
TAIL_N = 8
PROFILE_HPA = np.array(
    [10.0, 20.0, 30.0, 40.0, 50.0, 70.0, 100.0],
    dtype=np.float64,
)


def parse_merra(path: Path):
    lines = path.read_text(encoding="ascii", errors="strict").splitlines()
    pressure = None
    rows = {}
    for line in lines:
        if line.lstrip().startswith("P (hPa):"):
            pressure = np.array([float(x) for x in line.split(":", 1)[1].split()], dtype=float)
        elif re.match(r"^\s*\d{6}\s+", line):
            if pressure is None:
                raise ValueError("Data row encountered before pressure header")
            fields = line.split()
            stamp = fields[0]
            if len(fields) != 1 + len(pressure):
                raise ValueError(f"Unexpected field count for {stamp}: {len(fields)}")
            rows[(int(stamp[:4]), int(stamp[4:]))] = np.array(
                [float(x) for x in fields[1:]], dtype=np.float64
            )
    if pressure is None or len(rows) == 0:
        raise ValueError("MERRA-2 header or data rows not found")
    return pressure, rows


def date_sequence(year: int):
    # September of the year before March through September endpoint.
    return [(year - 1, month) for month in range(9, 13)] + [(year, month) for month in range(1, 10)]


def write_target(path: Path, values: np.ndarray, source_years: list[int], phase: str):
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(suffix=".nc", delete=False) as handle:
        temporary = Path(handle.name)
    with Dataset(temporary, "w", format="NETCDF4_CLASSIC") as ds:
        ds.createDimension("level", values.shape[0])
        ds.createDimension("time", values.shape[1])
        lev = ds.createVariable("level", "f4", ("level",))
        tim = ds.createVariable("time", "i4", ("time",))
        date = ds.createVariable("date", "i4", ("time",))
        secs = ds.createVariable("secs", "i4", ("time",))
        # NetCDF's C order is reversed by nf90_get_var into u_inp(level,time).
        qbo = ds.createVariable("qbo", "f8", ("time", "level"))
        lev.units = "hPa"
        tim.long_name = "monthly target index"
        date.long_name = "calendar date YYYYMMDD"
        secs.long_name = "seconds of day"
        qbo.units = "m s-1"
        lev[:] = PROFILE_HPA
        tim[:] = np.arange(values.shape[1], dtype=np.int32)
        dates = [(1, m) for m in range(9, 13)] + [(2, m) for m in range(1, 10)]
        date[:] = np.array([y * 10000 + m * 100 + 1 for y, m in dates], dtype=np.int32)
        secs[:] = 0
        qbo[:] = values.T
        ds.title = f"Deterministic full-profile {phase} QBO target"
        ds.source = "MERRA-2 GSFC 00N monthly zonal wind; fixed March 70-hPa tail selection"
        ds.selection_years = ",".join(str(y) for y in source_years)
        ds.selection_rule = f"{TAIL_N} highest/lowest March 70-hPa anomalies, {START_YEAR}-{END_YEAR}"
        ds.pressure_order = "top_to_bottom"
    shutil.copy2(temporary, path)
    temporary.unlink()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    pressure, rows = parse_merra(args.input)
    idx70 = int(np.argmin(np.abs(pressure - 70.0)))
    months = np.array([rows[(y, m)] for y in range(START_YEAR, END_YEAR + 1) for m in range(1, 13)])
    climatology = np.stack([months[np.arange(m - 1, len(months), 12)].mean(axis=0) for m in range(1, 13)])
    march = []
    for year in range(START_YEAR, END_YEAR + 1):
        value = rows[(year, 3)][idx70] - climatology[2, idx70]
        march.append((float(value), year))
    march.sort()
    east_years = [y for _, y in march[:TAIL_N]]
    west_years = [y for _, y in march[-TAIL_N:]][::-1]

    # Average each relative month, not a single 12-month mean.
    def profile(years):
        return np.stack([
            np.mean([rows[key] for y in years for key in [date_sequence(y)[i]]], axis=0)
            for i in range(13)
        ], axis=1)

    # Reorder from MERRA pressure descending to QBO module top-to-bottom.
    selected = np.array([np.argmin(np.abs(pressure - p)) for p in PROFILE_HPA])
    west = profile(west_years)[selected, :]
    east = profile(east_years)[selected, :]
    out = args.output
    write_target(out / "qbo_target_W.nc", west, west_years, "WQBO")
    write_target(out / "qbo_target_E.nc", east, east_years, "EQBO")
    anomalies = {"WQBO": west_years, "EQBO": east_years}
    with (out / "selection.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["phase", "year", "march_70_anomaly_mps", "rank"])
        for rank, (value, year) in enumerate(march, start=1):
            phase = "EQBO" if year in east_years else "WQBO" if year in west_years else "unused"
            writer.writerow([phase, year, f"{value:.8f}", rank])
    metadata = {
        "input": str(args.input),
        "input_sha256": sha256(args.input),
        "period": [START_YEAR, END_YEAR],
        "tail_n": TAIL_N,
        "selection_variable": "March 70-hPa MERRA-2 00N zonal-wind anomaly",
        "profile_hpa_top_to_bottom": PROFILE_HPA.tolist(),
        "qbo_months": ["Sep(-1)", "Oct(-1)", "Nov(-1)", "Dec(-1)", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep(+)"],
        "source_years": anomalies,
        "target_files": {"WQBO": sha256(out / "qbo_target_W.nc"), "EQBO": sha256(out / "qbo_target_E.nc")},
        "same_across_backgrounds": True,
        "selection_uses_outcomes": False,
    }
    (out / "target_metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
