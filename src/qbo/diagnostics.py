"""Evaluate the frozen boundary operator on one complete latitude–longitude snapshot."""
import argparse
import csv
import json
from pathlib import Path

from netCDF4 import Dataset
import numpy as np

from .boundary import column_taper, safe_pressure
from .io import sha256
from .zonal import zonal_envelope


def diagnose(path: Path, output: Path):
    with Dataset(path) as ds:
        def read(name, dimensions):
            var = ds[name]
            if var.dimensions not in dimensions:
                raise ValueError(f"Unexpected dimensions for {name}: {var.dimensions}")
            return np.ma.asarray(var[:], dtype=float).filled(np.nan)

        lat = read("lat", [("lat",)])
        lon = read("lon", [("lon",)])
        pressure = read("tropopause_hpa", [("diagnostic", "lat", "lon")])
        found = read("found", [("diagnostic", "lat", "lon")])
        bottom = read("layer_bottom_hpa", [("level",), ("level", "lat", "lon")])
        input_kind = str(getattr(ds, "input_kind", "unspecified"))
    if (not lat.size or not np.isfinite(lat).all() or np.any(np.abs(lat) >= 90)
            or len(np.unique(lat)) != len(lat)):
        raise ValueError("Provide unique finite latitudes between -90 and 90 degrees")
    spacing = np.diff(np.sort(lat))
    if spacing.size and not np.allclose(spacing, spacing[0], atol=1e-6, rtol=0):
        raise ValueError("Area summaries require equally spaced latitudes")
    sorted_lon = np.sort(np.mod(lon, 360))
    gaps = np.diff(np.r_[sorted_lon, sorted_lon[:1] + 360])
    if (lon.size < 2 or not np.isfinite(lon).all()
            or not np.allclose(gaps, 360 / lon.size, atol=1e-6, rtol=0)):
        raise ValueError("Provide a complete, unique, equally spaced longitude ring")
    if not bottom.size or not np.isfinite(bottom).all() or np.any(bottom <= 0):
        raise ValueError("Layer-bottom pressures must be finite and positive")
    if bottom.ndim == 1:
        bottom = bottom[:, None, None]
    safe = safe_pressure(pressure, found)
    column = column_taper(bottom, safe)
    zonal = np.broadcast_to(zonal_envelope(column), column.shape)
    area = np.broadcast_to(np.cos(np.deg2rad(lat))[:, None], safe.shape)
    rows = []
    for level in range(len(column)):
        rows.append({
            "level_index": level,
            "layer_bottom_min_hpa": float(bottom[level].min()),
            "layer_bottom_max_hpa": float(bottom[level].max()),
            "column_mean_weight": float(np.average(column[level], weights=area)),
            "zonal_mean_weight": float(np.average(zonal[level], weights=area)),
            "zonal_positive_area_fraction": float(np.average(zonal[level] > 0, weights=area)),
        })
    report = {"input_sha256": sha256(path), "input_kind": input_kind,
              "latitude_weighting": "cos(latitude), equally spaced latitude centers",
              "invalid_tropopause_columns": int((~np.isfinite(safe)).sum()), "levels": rows}
    output.mkdir(parents=True, exist_ok=False)
    with Dataset(output / "weights.nc", "w", format="NETCDF4_CLASSIC") as ds:
        for name, size in (("level", len(column)), ("lat", len(lat)), ("lon", len(lon))):
            ds.createDimension(name, size)
        for name, values, units in (("lat", lat, "degrees_north"), ("lon", lon, "degrees_east")):
            var = ds.createVariable(name, "f8", (name,))
            var[:] = values
            var.units = units
        for name, values in (("column_weight", column), ("zonal_weight", zonal),
                             ("layer_bottom_hpa", np.broadcast_to(bottom, column.shape))):
            var = ds.createVariable(name, "f8", ("level", "lat", "lon"))
            var[:] = values
            var.units = "hPa" if name == "layer_bottom_hpa" else "1"
        ds.input_sha256 = report["input_sha256"]
        ds.input_kind = input_kind
    with (output / "summary.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    (output / "report.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    diagnose(args.input, args.output)


if __name__ == "__main__":
    main()
