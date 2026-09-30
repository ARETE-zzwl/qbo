"""Regression checks against the project's actual input builders."""
import ctypes
import sys
import tempfile
import unittest
from pathlib import Path

import netCDF4
import numpy as np


from qbo import targets, backgrounds, validate as validator


def native_fortran_shape_read(path, levels, times):
    """Use NetCDF-C with the reversed count used by nf90_get_var(level,time)."""
    library_dir = Path(netCDF4.__file__).resolve().parent.parent / "netcdf4.libs"
    library = ctypes.CDLL(str(next(library_dir.glob("netcdf-*.dll"))) if sys.platform == "win32" else netCDF4._netCDF4.__file__)
    integer_p = ctypes.POINTER(ctypes.c_int)
    size_p = ctypes.POINTER(ctypes.c_size_t)
    double_p = ctypes.POINTER(ctypes.c_double)
    library.nc_open.argtypes = [ctypes.c_char_p, ctypes.c_int, integer_p]
    library.nc_inq_varid.argtypes = [ctypes.c_int, ctypes.c_char_p, integer_p]
    library.nc_get_vara_double.argtypes = [
        ctypes.c_int, ctypes.c_int, size_p, size_p, double_p,
    ]
    library.nc_close.argtypes = [ctypes.c_int]
    ncid, varid = ctypes.c_int(), ctypes.c_int()
    status = library.nc_open(str(path).encode("ascii"), 0, ctypes.byref(ncid))
    if status != 0:
        raise RuntimeError(f"nc_open failed: {status}")
    try:
        status = library.nc_inq_varid(ncid, b"qbo", ctypes.byref(varid))
        if status != 0:
            raise RuntimeError(f"nc_inq_varid failed: {status}")
        start = (ctypes.c_size_t * 2)(0, 0)
        count = (ctypes.c_size_t * 2)(times, levels)
        values = np.empty((times, levels), dtype=np.float64)
        status = library.nc_get_vara_double(
            ncid, varid, start, count, values.ctypes.data_as(double_p)
        )
        return status, values
    finally:
        library.nc_close(ncid)


class QboInterfaceTests(unittest.TestCase):
    def test_validator_rejects_archived_incompatible_layout(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "legacy.nc"
            with netCDF4.Dataset(path, "w") as dataset:
                dataset.createDimension("level", 7)
                dataset.createDimension("time", 13)
                dataset.createVariable("qbo", "f8", ("level", "time"))
            with self.assertRaisesRegex(ValueError, "on-disk qbo"):
                validator.read_nc(path)

    def test_validator_returns_canonical_level_time_values(self):
        values = np.arange(91, dtype=float).reshape(7, 13) / 10
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "target.nc"
            targets.write_target(path, values, [1985], "test")
            _, _, _, actual = validator.read_nc(path)
            np.testing.assert_array_equal(actual, values)

    def test_written_dimensions_match_fortran_reader(self):
        values = np.arange(91, dtype=float).reshape(7, 13) / 10
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "target.nc"
            targets.write_target(path, values, [1985], "test")
            with netCDF4.Dataset(path) as dataset:
                self.assertEqual(dataset["qbo"].dimensions, ("time", "level"))
                np.testing.assert_array_equal(dataset["qbo"][:], values.T)

    def test_native_read_accepts_fortran_array_shape(self):
        values = np.arange(91, dtype=float).reshape(7, 13) / 10
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "target.nc"
            targets.write_target(path, values, [1985], "test")
            status, actual = native_fortran_shape_read(path, 7, 13)
            self.assertEqual(status, 0, "NetCDF shape read must not return NC_EEDGE")
            np.testing.assert_array_equal(actual, values.T)


class LongitudeTests(unittest.TestCase):
    def test_interior_missing_source_still_uses_fallback(self):
        data = np.broadcast_to([1.0, np.nan, 3.0], (12, 2, 3)).copy()
        actual = backgrounds.regrid_nearest(
            data, [-10, 10], np.array([0.5, 179.5, 359.5]),
            np.array([0.0]), np.array([0.5, 179.5, 359.5]),
            np.full((12, 1, 3), 999.0),
        )
        expected = np.broadcast_to([1.0, 999.0, 3.0], (12, 1, 3))
        np.testing.assert_array_equal(actual, expected)

    def test_greenwich_uses_source_not_fallback(self):
        data = np.broadcast_to([1.0, 2.0, 3.0], (12, 2, 3)).copy()
        actual = backgrounds.regrid_nearest(
            data, [-10, 10], np.array([0.5, 179.5, 359.5]),
            np.array([0.0]), np.array([0.0]), np.full((12, 1, 1), 999.0),
        )
        np.testing.assert_array_equal(actual, np.full((12, 1, 1), 3.0))

    def test_periodic_equivalent_longitudes_match(self):
        data = np.broadcast_to([1.0, 2.0, 3.0], (12, 2, 3)).copy()
        actual = backgrounds.regrid_nearest(
            data, [-10, 10], np.array([0.5, 179.5, 359.5]),
            np.array([0.0]), np.array([-0.25, 359.75, 719.75]),
            np.full((12, 1, 3), 999.0),
        )
        np.testing.assert_array_equal(actual, np.full((12, 1, 3), 3.0))


if __name__ == "__main__":
    unittest.main()
