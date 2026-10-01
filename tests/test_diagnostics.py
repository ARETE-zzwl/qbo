import json

from netCDF4 import Dataset
import numpy as np
import pytest

from qbo.diagnostics import diagnose


def snapshot(path, bottom=(50., 80.), invalid_column=False):
    with Dataset(path, "w") as ds:
        for name, size in (("diagnostic", 3), ("level", 2), ("lat", 2), ("lon", 4)):
            ds.createDimension(name, size)
        ds.input_kind = "synthetic"
        ds.createVariable("lat", "f8", ("lat",))[:] = [-5, 5]
        ds.createVariable("lon", "f8", ("lon",))[:] = [0, 90, 180, 270]
        ds.createVariable("layer_bottom_hpa", "f8", ("level",))[:] = bottom
        ds.createVariable("tropopause_hpa", "f8", ("diagnostic", "lat", "lon"))[:] = 100
        found = ds.createVariable("found", "i4", ("diagnostic", "lat", "lon"))
        found[:] = 1
        if invalid_column:
            found[0, 0, 1] = 0
    return path


def test_diagnostics_preserves_zero_support_and_uniform_ring(tmp_path):
    path = snapshot(tmp_path / "state.nc", invalid_column=True)
    out = tmp_path / "report"
    result = diagnose(path, out)
    with Dataset(out / "weights.nc") as ds:
        column = ds["column_weight"][:]
        zonal = ds["zonal_weight"][:]
        assert np.all(column[:, 0, 1] == 0)
        assert np.all(zonal[:, 0, :] == 0)
        np.testing.assert_allclose(zonal[:, 1, :], 1)
        assert np.all(zonal <= column)
        assert ds.input_kind == "synthetic"
    assert result["levels"][0]["zonal_mean_weight"] == pytest.approx(.5)
    assert result["input_kind"] == "synthetic"
    assert json.loads((out / "report.json").read_text())["input_sha256"]
    with pytest.raises(FileExistsError):
        diagnose(path, out)


def test_diagnostics_reports_no_support_at_deep_levels(tmp_path):
    path = snapshot(tmp_path / "state.nc", bottom=(100., 110.))
    result = diagnose(path, tmp_path / "report")
    assert all(level["zonal_mean_weight"] == 0 for level in result["levels"])


def test_diagnostics_accepts_column_dependent_layer_bottoms(tmp_path):
    path = snapshot(tmp_path / "state.nc")
    with Dataset(path, "a") as ds:
        ds.renameVariable("layer_bottom_hpa", "old_bottom")
        ds.createVariable("layer_bottom_hpa", "f8", ("level", "lat", "lon"))[:] = 110
    result = diagnose(path, tmp_path / "report")
    assert result["levels"][0]["column_mean_weight"] == 0


def test_diagnostics_rejects_ambiguous_longitude_grid(tmp_path):
    path = snapshot(tmp_path / "state.nc")
    with Dataset(path, "a") as ds:
        ds["lon"][:] = [0, 90, 180, 360]
    with pytest.raises(ValueError, match="longitude"):
        diagnose(path, tmp_path / "report")


def test_diagnostics_rejects_irregular_latitude_spacing(tmp_path):
    path = tmp_path / "state.nc"
    with Dataset(path, "w") as ds:
        for name, size in (("diagnostic", 3), ("level", 1), ("lat", 3), ("lon", 4)):
            ds.createDimension(name, size)
        ds.createVariable("lat", "f8", ("lat",))[:] = [-10, 0, 30]
        ds.createVariable("lon", "f8", ("lon",))[:] = [0, 90, 180, 270]
        ds.createVariable("layer_bottom_hpa", "f8", ("level",))[:] = 70
        ds.createVariable("tropopause_hpa", "f8", ("diagnostic", "lat", "lon"))[:] = 100
        ds.createVariable("found", "i4", ("diagnostic", "lat", "lon"))[:] = 1
    with pytest.raises(ValueError, match="latitude"):
        diagnose(path, tmp_path / "report")
