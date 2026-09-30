"""Exercise input-builder entry points using small synthetic source datasets."""
import gzip
import json
import subprocess
import sys

import netCDF4
import numpy as np


def command(module, *args):
    return subprocess.run(
        [sys.executable, "-m", module, *map(str, args)],
        capture_output=True, text=True, check=True, timeout=60,
    )


def test_qbo_build_and_validate(tmp_path):
    table = tmp_path / "wind.txt"
    lines = ["P (hPa): 100 70 50 40 30 20 10"]
    for year in range(1980, 2025):
        for month in range(1, 13):
            wind = (year - 2002.5) * .8
            lines.append(f"{year}{month:02d} " + " ".join([str(wind)] * 7))
    table.write_text("\n".join(lines), encoding="ascii")
    output = tmp_path / "targets"
    command("qbo.targets", "--input", table, "--output", output)
    command("qbo.validate", output, "--output", output / "qc.json")
    qc = json.loads((output / "qc.json").read_text())
    metadata = json.loads((output / "target_metadata.json").read_text())
    assert qc["passed"]
    assert metadata["source_years"]["WQBO"] == list(range(2024, 2016, -1))
    assert metadata["source_years"]["EQBO"] == list(range(1981, 1989))
    np.testing.assert_allclose(qc["march_70_W_mps"], 14.4, rtol=0, atol=1e-12)
    np.testing.assert_allclose(qc["march_70_E_mps"], -14.4, rtol=0, atol=1e-12)
    for phase in ("W", "E"):
        with netCDF4.Dataset(output / f"qbo_target_{phase}.nc") as data:
            assert data["qbo"].dimensions == ("time", "level")
            assert data["qbo"].shape == (13, 7)


def write_hadisst(path, variable):
    with netCDF4.Dataset(path, "w") as data:
        data.createDimension("time", 44 * 12)
        data.createDimension("latitude", 2)
        data.createDimension("longitude", 3)
        time = data.createVariable("time", "f8", ("time",))
        time.units = "days since 1981-01-01"
        time.calendar = "360_day"
        time[:] = np.arange(44 * 12) * 30
        data.createVariable("latitude", "f8", ("latitude",))[:] = [10, -10]
        data.createVariable("longitude", "f8", ("longitude",))[:] = [-179.5, -.5, .5]
        values = data.createVariable(variable, "f4", ("time", "latitude", "longitude"))
        if variable == "sst":
            year = np.repeat(np.arange(1, 45), 12) * .01
            month = np.tile(np.arange(1, 13), 44) * .1
            values[:] = np.broadcast_to((year + month)[:, None, None], (528, 2, 3))
        else:
            values[:] = .5
    compressed = path.with_suffix(".nc.gz")
    with gzip.open(compressed, "wb") as stream:
        stream.write(path.read_bytes())
    return compressed


def test_hadisst_background_workflow(tmp_path):
    sst = write_hadisst(tmp_path / "sst.nc", "sst")
    ice = write_hadisst(tmp_path / "ice.nc", "sic")
    baseline = tmp_path / "baseline.nc"
    with netCDF4.Dataset(baseline, "w") as data:
        for name, size in [("time", 12), ("lat", 1), ("lon", 3)]:
            data.createDimension(name, size)
        data.createVariable("lat", "f8", ("lat",))[:] = [0]
        data.createVariable("lon", "f8", ("lon",))[:] = [0, 180, 359.5]
        time = data.createVariable("time", "f8", ("time",))
        time.units = "days since 0001-01-01"
        time[:] = np.arange(12) * 30
        data.createVariable("date", "i4", ("time",))[:] = 10101 + np.arange(12) * 100
        data.createVariable("datesec", "i4", ("time",))[:] = 0
        for name, value, units in [("SST_cpl", 999., "deg_C"), ("ice_cov", 0., "1")]:
            field = data.createVariable(name, "f4", ("time", "lat", "lon"))
            field.long_name = name
            field.units = units
            field[:] = value
    output = tmp_path / "backgrounds"
    command("qbo.backgrounds", "--sst-gz", sst, "--ice-gz", ice,
            "--baseline", baseline, "--output", output)
    metadata = json.loads((output / "background_metadata.json").read_text())
    for label, year_mean in [("early", .075), ("late", .295)]:
        with netCDF4.Dataset(output / f"sst_HadISST_{label}_0.9x1.25_clim.nc") as data:
            expected = np.broadcast_to((year_mean + np.arange(1, 13) * .1)[:, None, None], (12, 1, 3))
            np.testing.assert_allclose(data["SST_cpl"][:], expected, rtol=1e-6)
            np.testing.assert_array_equal(data["ice_cov"][:], np.full((12, 1, 3), .5))
            assert data["SST_cpl"].dimensions == ("time", "lat", "lon")
        assert metadata["outputs"][label]["sst_missing_filled"] == 0
