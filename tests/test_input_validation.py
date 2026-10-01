"""Reject ambiguous input records and preserve missing-data masks."""
import netCDF4
import numpy as np
import pytest

from qbo.backgrounds import climatology
from qbo.targets import parse_merra


def write_merra(path, rows, pressure=(100, 70, 50, 40, 30, 20, 10)):
    lines = ["P (hPa): " + " ".join(map(str, pressure))]
    lines.extend(stamp + " " + " ".join(map(str, winds)) for stamp, winds in rows)
    path.write_text("\n".join(lines), encoding="ascii")
    return path


def write_hadisst(path, month_indices, values, fill_value=9999.0):
    values = np.asarray(values, dtype=np.float32)
    with netCDF4.Dataset(path, "w") as data:
        data.createDimension("time", len(month_indices))
        data.createDimension("latitude", values.shape[1])
        data.createDimension("longitude", values.shape[2])
        time = data.createVariable("time", "f8", ("time",))
        time.units = "days since 2000-01-01"
        time.calendar = "360_day"
        time[:] = np.asarray(month_indices) * 30
        data.createVariable("latitude", "f8", ("latitude",))[:] = np.arange(values.shape[1])
        data.createVariable("longitude", "f8", ("longitude",))[:] = np.arange(values.shape[2])
        data.createVariable(
            "sst", "f4", ("time", "latitude", "longitude"), fill_value=fill_value,
        )[:] = values
    return path


def test_merra_rejects_duplicate_month(tmp_path):
    path = write_merra(tmp_path / "wind.txt", [("200001", [1] * 7), ("200001", [2] * 7)])
    with pytest.raises(ValueError, match="(?i)duplicate|month"):
        parse_merra(path)


@pytest.mark.parametrize("stamp", ["200000", "200013"])
def test_merra_rejects_invalid_month(tmp_path, stamp):
    path = write_merra(tmp_path / "wind.txt", [(stamp, [1] * 7)])
    with pytest.raises(ValueError, match="(?i)month"):
        parse_merra(path)


@pytest.mark.parametrize("invalid", [float("nan"), float("inf"), float("-inf")])
def test_merra_rejects_nonfinite_wind(tmp_path, invalid):
    path = write_merra(tmp_path / "wind.txt", [("200001", [invalid] + [1] * 6)])
    with pytest.raises(ValueError, match="(?i)finite"):
        parse_merra(path)


def test_merra_requires_exact_target_pressure_levels(tmp_path):
    path = write_merra(
        tmp_path / "wind.txt", [("200001", [1] * 7)],
        pressure=(100, 75, 50, 40, 30, 20, 10),
    )
    with pytest.raises(ValueError, match="(?i)pressure|level|70"):
        parse_merra(path)


def test_merra_rejects_duplicate_pressure_levels(tmp_path):
    path = write_merra(
        tmp_path / "wind.txt", [("200001", [1] * 8)],
        pressure=(100, 70, 70, 50, 40, 30, 20, 10),
    )
    with pytest.raises(ValueError, match="(?i)duplicate|pressure"):
        parse_merra(path)


def test_hadisst_rejects_duplicate_month_despite_correct_record_count(tmp_path):
    path = write_hadisst(
        tmp_path / "sst.nc", [0, 0, *range(2, 12)], np.ones((12, 1, 1)),
    )
    with pytest.raises(ValueError, match="(?i)duplicate|month"):
        climatology(path, "sst", (2000, 2000))


def test_hadisst_groups_shuffled_records_by_calendar_month(tmp_path):
    months = np.array([11, 4, 0, 8, 2, 7, 1, 6, 10, 3, 9, 5])
    path = write_hadisst(tmp_path / "sst.nc", months, (months + 1)[:, None, None])
    _, _, monthly = climatology(path, "sst", (2000, 2000))
    np.testing.assert_allclose(monthly[:, 0, 0], np.arange(1, 13))


def test_hadisst_excludes_positive_fill_value_from_monthly_mean(tmp_path):
    values = np.full((24, 1, 1), 2.0)
    values[0, 0, 0] = 9999.0
    values[12, 0, 0] = 14.0
    path = write_hadisst(tmp_path / "sst.nc", np.arange(24), values)
    _, _, monthly = climatology(path, "sst", (2000, 2001))
    np.testing.assert_allclose(monthly[:, 0, 0], [14] + [2] * 11)


def test_hadisst_keeps_all_missing_cells_nan_for_fallback(tmp_path):
    values = np.full((12, 1, 2), 10.0)
    values[:, :, 0] = 9999.0
    path = write_hadisst(tmp_path / "sst.nc", np.arange(12), values)
    _, _, monthly = climatology(path, "sst", (2000, 2000))
    assert np.isnan(monthly[:, 0, 0]).all()
    np.testing.assert_allclose(monthly[:, 0, 1], 10)
