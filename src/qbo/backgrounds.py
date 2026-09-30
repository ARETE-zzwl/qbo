"""Build monthly HadISST SST and sea-ice climatologies on a baseline FV grid.

Early and late periods are 1981–1994 and 1995–2024. Missing source cells
use the corresponding baseline values after periodic nearest-neighbour mapping.
"""
from __future__ import annotations

import argparse
import gzip
import json
import os
import shutil
import tempfile
from pathlib import Path

from .io import sha256

import numpy as np
from netCDF4 import Dataset, num2date
from scipy.interpolate import RegularGridInterpolator


def unzip_to_temp(path: Path, suffix: str) -> Path:
    fd, name = tempfile.mkstemp(suffix=suffix)
    os.close(fd)
    Path(name).unlink()
    with gzip.open(path, "rb") as src, open(name, "wb") as dst:
        shutil.copyfileobj(src, dst)
    return Path(name)


def climatology(path: Path, variable: str, years: tuple[int, int]):
    with Dataset(path) as ds:
        time_var = ds.variables["time"]
        dates = num2date(time_var[:], time_var.units, getattr(time_var, "calendar", "standard"))
        idx = [i for i, date in enumerate(dates) if years[0] <= date.year <= years[1]]
        if len(idx) != (years[1] - years[0] + 1) * 12:
            raise ValueError(f"Incomplete HadISST period {years}: {len(idx)} months")
        lat = np.asarray(ds.variables["latitude"][:], dtype=float)
        lon = np.asarray(ds.variables["longitude"][:], dtype=float)
        values = np.asarray(ds.variables[variable][idx], dtype=np.float32)
    values[values < -100] = np.nan
    # HadISST is latitude-descending and longitude -179.5..179.5.
    lat_order = np.argsort(lat)
    lon360 = np.mod(lon, 360.0)
    lon_order = np.argsort(lon360)
    lat = lat[lat_order]
    lon360 = lon360[lon_order]
    values = values[:, lat_order][:, :, lon_order]
    monthly = np.stack([np.nanmean(values[m::12], axis=0) for m in range(12)], axis=0)
    return lat, lon360, monthly


def regrid_nearest(data, lat, lon, target_lat, target_lon, fallback):
    lon_ext = np.r_[lon[-1] - 360.0, lon, lon[0] + 360.0]
    # data is (month, latitude, longitude); close the longitude seam only.
    data_ext = np.concatenate([data[:, :, -1:], data, data[:, :, :1]], axis=2)
    points = np.column_stack([np.repeat(target_lat, target_lon.size), np.tile(np.mod(target_lon, 360.0), target_lat.size)])
    out = np.empty((12, target_lat.size, target_lon.size), dtype=np.float32)
    for month in range(12):
        interp = RegularGridInterpolator((lat, lon_ext), data_ext[month], method="nearest", bounds_error=False, fill_value=np.nan)
        values = interp(points).reshape(target_lat.size, target_lon.size)
        missing = ~np.isfinite(values)
        values[missing] = fallback[month][missing]
        out[month] = values
    return out


def write_background(path: Path, baseline: dict, sst: np.ndarray, ice: np.ndarray, period: tuple[int, int], source_hashes: dict):
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(suffix=".nc", delete=False) as handle:
        temp = Path(handle.name)
    try:
        with Dataset(temp, "w", format="NETCDF4_CLASSIC") as ds:
            ds.createDimension("lon", baseline["lon"].size)
            ds.createDimension("lat", baseline["lat"].size)
            ds.createDimension("time", 12)
            for name in ("SST_cpl", "SST_cpl_prediddle", "ice_cov", "ice_cov_prediddle"):
                var = ds.createVariable(name, "f4", ("time", "lat", "lon"))
                var.long_name = baseline["attrs"][name]["long_name"]
                var.units = baseline["attrs"][name]["units"]
            lat = ds.createVariable("lat", "f8", ("lat",))
            lon = ds.createVariable("lon", "f8", ("lon",))
            time = ds.createVariable("time", "f8", ("time",))
            date = ds.createVariable("date", "i4", ("time",))
            datesec = ds.createVariable("datesec", "i4", ("time",))
            lat.long_name, lat.units = "latitude", "degrees_north"
            lon.long_name, lon.units = "longitude", "degrees_east"
            time.units = baseline["time_units"]
            lat[:] = baseline["lat"]
            lon[:] = baseline["lon"]
            time[:] = baseline["time"]
            date[:] = baseline["date"]
            datesec[:] = baseline["datesec"]
            ds.variables["SST_cpl"][:] = sst
            ds.variables["SST_cpl_prediddle"][:] = sst
            ds.variables["ice_cov"][:] = np.clip(ice, 0.0, 1.0)
            ds.variables["ice_cov_prediddle"][:] = np.clip(ice, 0.0, 1.0)
            ds.title = f"HadISST {period[0]}-{period[1]} monthly climatological boundary"
            ds.source = "Met Office HadISST SST and sea ice; nearest-neighbour map to f09 FV grid"
            ds.period = f"{period[0]}-{period[1]}"
            ds.regrid_method = "nearest_neighbor"
            ds.source_sha256 = json.dumps(source_hashes, sort_keys=True)
        shutil.copy2(temp, path)
    finally:
        temp.unlink(missing_ok=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sst-gz", type=Path, required=True)
    ap.add_argument("--ice-gz", type=Path, required=True)
    ap.add_argument("--baseline", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    sst_nc = unzip_to_temp(args.sst_gz, ".sst.nc")
    ice_nc = unzip_to_temp(args.ice_gz, ".ice.nc")
    try:
        with tempfile.NamedTemporaryFile(suffix=".nc", delete=False) as handle:
            baseline_temp = Path(handle.name)
        baseline_temp.unlink()
        shutil.copy2(args.baseline, baseline_temp)
        with Dataset(baseline_temp) as ds:
            baseline = {
                "lat": np.asarray(ds.variables["lat"][:], dtype=float),
                "lon": np.asarray(ds.variables["lon"][:], dtype=float),
                "time": np.asarray(ds.variables["time"][:], dtype=float),
                "date": np.asarray(ds.variables["date"][:], dtype=np.int32),
                "datesec": np.asarray(ds.variables["datesec"][:], dtype=np.int32),
                "time_units": ds.variables["time"].units,
                "attrs": {name: {"long_name": ds.variables[name].long_name, "units": ds.variables[name].units} for name in ("SST_cpl", "ice_cov")},
            }
        baseline["attrs"]["SST_cpl_prediddle"] = baseline["attrs"]["SST_cpl"]
        baseline["attrs"]["ice_cov_prediddle"] = baseline["attrs"]["ice_cov"]
        with Dataset(baseline_temp) as ds:
            fallback_sst = np.asarray(ds.variables["SST_cpl"][:], dtype=np.float32)
            fallback_ice = np.asarray(ds.variables["ice_cov"][:], dtype=np.float32)
        source_hashes = {"HadISST_sst.nc.gz": sha256(args.sst_gz), "HadISST_ice.nc.gz": sha256(args.ice_gz), "baseline": sha256(args.baseline)}
        outputs = {}
        for label, period in (("early", (1981, 1994)), ("late", (1995, 2024))):
            slat, slon, sst_clim = climatology(sst_nc, "sst", period)
            ilat, ilon, ice_clim = climatology(ice_nc, "sic", period)
            if not (np.array_equal(slat, ilat) and np.array_equal(slon, ilon)):
                raise ValueError("HadISST SST and ice grids differ")
            sst = regrid_nearest(sst_clim, slat, slon, baseline["lat"], baseline["lon"], np.full_like(fallback_sst, np.nan))
            ice = regrid_nearest(ice_clim, ilat, ilon, baseline["lat"], baseline["lon"], np.full_like(fallback_ice, np.nan))
            sst_missing = ~np.isfinite(sst)
            ice_missing = ~np.isfinite(ice)
            sst[sst_missing] = fallback_sst[sst_missing]
            ice[ice_missing] = fallback_ice[ice_missing]
            path = args.output / f"sst_HadISST_{label}_0.9x1.25_clim.nc"
            write_background(path, baseline, sst, ice, period, source_hashes)
            outputs[label] = {
                "path": str(path),
                "sha256": sha256(path),
                "sst_range": [float(np.nanmin(sst)), float(np.nanmax(sst))],
                "ice_range": [float(np.nanmin(ice)), float(np.nanmax(ice))],
                "sst_missing_filled": int(sst_missing.sum()),
                "ice_missing_filled": int(ice_missing.sum()),
                "ice_values_clipped": int(np.count_nonzero((ice < 0) | (ice > 1))),
            }
        metadata = {"source_hashes": source_hashes, "grid": {"lat_count": int(baseline["lat"].size), "lon_count": int(baseline["lon"].size)}, "method": "nearest_neighbor", "outputs": outputs}
        (args.output / "background_metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    finally:
        sst_nc.unlink(missing_ok=True)
        ice_nc.unlink(missing_ok=True)
        baseline_temp.unlink(missing_ok=True)


if __name__ == "__main__":
    main()
